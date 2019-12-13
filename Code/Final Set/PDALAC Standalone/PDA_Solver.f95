! This file is responsible for conducting the progressive damage analysis and postprocessing 

    module PDA

    use material_parameters

    use analysis_parameters

    implicit none

    integer :: J, step, IMAX

    ! Arrays and varaibles declaration
    real*8, parameter :: Pi = 3.14159265358979d+0
    real*8, dimension  (nblock,6) :: Stress, Strain
    real*8, dimension  (nblock,6) :: dmg, e, fflags, d_index
    real*8:: F1,F2,F3,F11,F22,F33,F44,F55,F66,F12,F13,F23,Xphi(6)
    real*8:: E11D, E22D, E33D, G12D, G13D, G23D
    real*8:: THETA,Rn,Rn1,Rnt,SFP,TN1,TNT,Pnt,Pn1,THETAMAX,IFF,A
    real*8, dimension (nblock) :: Phi

! +------------+---------------------------------------+
! | Parameter  |              Description              |
! +------------+---------------------------------------+
! |  Stress    | Stress array                          |
! |  Strain    | Strain array                          |
! |  dmg       | Damage coefficients array             |
! |  e         | Failure indicies array                |
! |  fflags    | Failure flags array                   |
! |  d_index   | damage indicies for continium damage  |
! |  Phi       | Tsai-Wu/Hoffman failure index         |
! |  Xphi,Fij  | Tsai-Wu/Hoffman polynomials           |
! |  Eij       | Degeared material properties          |
! |  THETA     | Fracture plane angle for Puck         |
! | Rn,Rn1,Rnt | Puck fracture strength parameters     |
! |SFP,TN1,TNT | Stresses on fracrure plane for Puck   |
! | Pnt,Pn1    | Puck slope parameters                 |
! +------------+---------------------------------------+

    contains

    subroutine PDA_Solver

    print"(/)"
    print *, '___________________________________'
    print"(/,a33)","Step 3: Results Post-Processing"
    print *, '___________________________________'

!###########          Start Analysis Loop         ###########
!#----------------------------------------------------------#

    do step = 1,nsteps

        print"(/,/,a10,i2)","Step No:",step
        print"(a)","----------------------"


    if (step.eq.1) then 

!###########         Initialize  Variables        ###########
!#----------------------------------------------------------#

! At Step 1 the strain tensors are initialized to zero

! At Step 1 the material is undamaged (linear elastic)

! At Step 1 the failure flags and indicies == 0

        Strain = 0.0d+0
        dmg = 1.0d+0
        e = 0.0d+0
        fflags = 0.0d+0
        IFF=0.0d+0
        THETAMAX=0.0d+0

!###########     Calculate Initial Stresses       ###########
!#----------------------------------------------------------#

! call ortho3d subroutine to: form constitutive tensor,calculate stresses
                            
        call ortho3d(E11,E22,E33,G12,G13,G23,ANU12,ANU13,ANU23,Ypsilon,nblock,dmg,Stress,Strain) 

    else

!###########         Update Total Strain          ###########
!#----------------------------------------------------------#

! Loop through all material points and accumulate strains            

        do I= 1,nblock 
            do j=1,6  

                Strain(I, J)= Strain(I,J)+ StrainInc(I,J) 

            end do
        end do

!###########           Update Stresses            ###########
!#----------------------------------------------------------#

        call ortho3d(E11,E22,E33,G12,G13,G23,ANU12,ANU13,ANU23,Ypsilon,nblock,dmg,Stress,Strain)


!###########         Failure Evaluation           ###########
!#----------------------------------------------------------#

        call failure_check(failure_id)


!###########         Damage Evaluation            ###########
!#----------------------------------------------------------#

! Apply damage if any of the failure indicies excceds one
! or previous damage presist

    if ((maxval(e).GT.1.0) .OR. (minval(dmg).LT.1.0)) then

        call damage_check(damage_id)
        
!###########           Update Stresses            ###########
!#----------------------------------------------------------#

        call ortho3d(E11,E22,E33,G12,G13,G23,ANU12,ANU13,ANU23,Ypsilon,nblock,dmg,Stress,Strain) 
      
! Exception handle for damage values

        do I= 1,nblock
            do J=1,6

            if ((dmg(I,J).GT.1.0).OR.(dmg(I,J).LT.0.0)) then

                print*,'Encountered Incorrect results for damage variables'

                print*,'Check your damage model variables'

                call exit (0)
            end if 
            end do 
        end do
        
    end if
    
    end if

!###########           Post-Processing            ###########
!#----------------------------------------------------------#    
        print"(/)"
        print"(a31,/)","Effective Stress Components: "
        print"(6a11)","S11","S22","S33","S12","S13","S23"   
        print"(6e11.3)",Stress(1,:)

        print"(/)"
        print"(a27,/)","Total Strain Components: "
        print"(6a11)","Eps11","Eps22","Eps33","Eps12","Eps13","Eps23"   
        print"(6e11.3)",Strain(1,:)
        
        print"(/)"
        print"(a20,/)","Failure Indicies: "
        print"(6a11)","e1","e2","e3","e4","e5","e6"   
        print"(6f11.3)",e(1,:)             

        print"(/)"
        print"(a22,/)","Damage Coefficents: "
        print"(6a11)","d11","d22","d33","d12","d13","d23"   
        print"(6f11.3)",dmg(1,:)
        
        ! Write output file for GNU Plot 
        
        open (unit=1, file='data.dat')
        write(1,*) Strain(1,:), Stress(1,:)


    end do
    
    close(unit=1)

    end subroutine PDA_Solver


!###########           End Analysis Loop          ###########
!#----------------------------------------------------------#


!###########    Utility Subroutine 1: ortho3d     ###########
!#----------------------------------------------------------#

! The following subroutine :

! 1. Build the effective constiutive matrix

! 2. Return the updated stresses

    subroutine ortho3D (E11,E22,E33,G12,G13,G23,ANU12,ANU13,ANU23,Ypsilon,nblock,dmg,Stress,Strain) 

    implicit none
    integer :: I, nblock
    real*8, dimension (nblock,6) :: dmg, Stress, Strain
    real*8, dimension (6,6) :: DDSDDE
    real*8 :: E11, E22, E33, G12, G13, G23, ANU12, ANU13, ANU23, Ypsilon

    select case (damage_id)

    case (5)  ! effective stress for continium damage mechanics

    do I =1,nblock   ! Loop through all material points 

        Ypsilon=(1)/(1-dmg(I,1)*dmg(I,2)*ANU12*ANU21&
        &-dmg(I,2)*dmg(I,3)*ANU23*ANU32&
        &-dmg(I,1)*dmg(I,3)*ANU13*ANU31&
        &-2*dmg(I,1)*dmg(I,2)*dmg(I,3)*ANU21*ANU32*ANU13)

! Form effective orthotropic material matrix including damage (Jacobian Matrix) 

        DDSDDE(1,1)=(E11*Ypsilon*(1.0-dmg(I,2)*dmg(I,3)*ANU23*ANU32))*dmg(I,1)

        DDSDDE(2,2)=(E22*Ypsilon*(1.0-dmg(I,1)*dmg(I,3)*ANU13*ANU31))*dmg(I,2)

        DDSDDE(3,3)=(E33*Ypsilon*(1.0-dmg(I,1)*dmg(I,2)*ANU12*ANU21))*dmg(I,3)

        DDSDDE(1,2)=(E11*Ypsilon*(ANU21+ANU31*ANU23*dmg(I,3)))*dmg(I,1)*dmg(I,2)

        DDSDDE(1,3)=(E11*Ypsilon*(ANU31+ANU21*ANU32*dmg(I,2)))*dmg(I,1)*dmg(I,3)
                                
        DDSDDE(2,3)=(E22*Ypsilon*(ANU32+dmg(I,1)*ANU12*ANU31))*dmg(I,3)*dmg(I,2)

        DDSDDE(2,1)=DDSDDE(1,2)

        DDSDDE(3,1)=DDSDDE(1,3) 

        DDSDDE(3,2)=DDSDDE(2,3)

        DDSDDE(4,4)=G12*dmg(I,4)

        DDSDDE(5,5)=G13*dmg(I,5)

        DDSDDE(6,6)=G23*dmg(I,6)

! Calculate Stress Components

        Stress(I,1)= (DDSDDE(1,1)*Strain(I,1)+DDSDDE(1,2)*Strain(I,2)+DDSDDE(1,3)*Strain(I,3))
        Stress(I,2)= (DDSDDE(2,1)*Strain(I,1)+DDSDDE(2,2)*Strain(I,2)+DDSDDE(2,3)*Strain(I,3))
        Stress(I,3)= (DDSDDE(3,1)*Strain(I,1)+DDSDDE(3,2)*Strain(I,2)+DDSDDE(3,3)*Strain(I,3)) 
        Stress(I,4)= (DDSDDE(4,4)*Strain(I,4))
        Stress(I,5)= (DDSDDE(5,5)*Strain(I,5))        
        Stress(I,6)= (DDSDDE(6,6)*Strain(I,6))

    end do 

    case default  ! effective stress for ply-discount damage

    do I = 1,nblock   ! Loop through all material points 

! Form effective orthotropic material matrix including damage (Jacobian Matrix) 

        DDSDDE(1,1)=(E11*Ypsilon*(1.0-ANU23*ANU32))*abs(dmg(I,1))

        DDSDDE(2,2)=(E22*Ypsilon*(1.0-ANU13*ANU31))*abs(dmg(I,2))

        DDSDDE(3,3)=(E33*Ypsilon*(1.0-ANU12*ANU21))*abs(dmg(I,3))       

        DDSDDE(1,2)=(E11*Ypsilon*(ANU21+ANU31*ANU23))*abs(dmg(I,1)*dmg(I,2))

        DDSDDE(1,3)=(E11*Ypsilon*(ANU31+ANU21*ANU32))*abs(dmg(I,1)*dmg(I,3))

        DDSDDE(2,3)=(E22*Ypsilon*(ANU32+ANU12*ANU31))*abs(dmg(I,3)*dmg(I,2))       

        DDSDDE(2,1)=DDSDDE(1,2)

        DDSDDE(3,1)=DDSDDE(1,3)

        DDSDDE(3,2)=DDSDDE(2,3)

        DDSDDE(4,4)=G12*abs(dmg(I,4))

        DDSDDE(5,5)=G13*abs(dmg(I,5))

        DDSDDE(6,6)=G23*abs(dmg(I,6))

! Calculate Stress Components

        Stress(I,1)= (DDSDDE(1,1)*Strain(I,1)+DDSDDE(1,2)*Strain(I,2)+DDSDDE(1,3)*Strain(I,3))
        Stress(I,2)= (DDSDDE(2,1)*Strain(I,1)+DDSDDE(2,2)*Strain(I,2)+DDSDDE(2,3)*Strain(I,3))
        Stress(I,3)= (DDSDDE(3,1)*Strain(I,1)+DDSDDE(3,2)*Strain(I,2)+DDSDDE(3,3)*Strain(I,3)) 
        Stress(I,4)= (DDSDDE(4,4)*Strain(I,4))
        Stress(I,5)= (DDSDDE(5,5)*Strain(I,5))        
        Stress(I,6)= (DDSDDE(6,6)*Strain(I,6))

    end do

    end select

    end subroutine ortho3d


!########### Utility Subroutine 2: failure_check  ###########
!#----------------------------------------------------------#

! The following subroutine :

! Evaluates the failure indicies [e] based on selected failure criteria


    subroutine failure_check (id)

    implicit none 

    integer :: id

    select case (id)



    case(1)
!#-------------------Max Stress Criteria--------------------#

    do I= 1,nblock   ! Loop through all material points and calculate failure indicies

! Failure along fiber (11 axis) (Tension/Compression)

        if(Stress(I,1).GE.0.0) then
            e(I,1)=Stress(I,1)/Xt
        else
            e(I,1)=-Stress(I,1)/Xc
        end if 

! Failure along matrix (22 axis) (Tension/Compression)

        if(Stress(I,2).GE.0.0) then 
            e(I,2)=Stress(I,2)/Yt 
        else
            e(I,2)=-Stress(I,2)/Yc 
        end if

! Interlaminar failure (33 axis) (Tension/Compression)                

        if(Stress(I,3).GE.0.0) then 
            e(I,3)=Stress(I,3)/Zt
        else         
            e(I,3)=-Stress(I,3)/Zc
        end if

! Shear Failure (1-2 plane) 

        e(I,4)=abs(Stress(I,4))/S12

! Shear Failure (1-3 plane)

        e(I,5)=abs(Stress(I,5))/S13

! Shear Failure (2-3 plane)

        e(I,6)=abs(Stress(I,6))/S23

    end do


    case(2) 
!#-------------------Max Strain Criteria--------------------#

    do I= 1,nblock   ! Loop through all material points and calculate failure indicies

! Failure along fiber (11 axis) (Tension/Compression)

    if(Strain(I,1).GE.0.0) then
        e(I,1)=Strain(I,1)/EpsXt
    else
        e(I,1)=-Strain(I,1)/EpsXc
    end if 

! Failure along matrix (22 axis) (Tension/Compression)

    if(Stress(I,2).GE.0.0) then 
        e(I,2)=Strain(I,2)/EpsYt 
    else
        e(I,2)=-Strain(I,2)/EpsYc  
    end if

! Interlaminar failure (33 axis) (Tension/Compression)                    

    if(Stress(I,3).GE.0.0) then 
        e(I,3)=Strain(I,3)/EpsZt
    else         
        e(I,3)=-Strain(I,3)/EpsZc
    end if

! Shear Failure (1-2 plane) 

        e(I,4)=abs(Strain(I,4))/GamS12

! Shear Failure (1-3 plane)

        e(I,5)=Strain(I,5)/GamS13

! Shear Failure (2-3 plane)

        e(I,6)=abs(Strain(I,6))/GamS23

    end do


    case(3) 
!#-------------------Tsai-Wu Criteria--------------------#

! Calculate Tsai-Wu Strength parameters 

    F1=(1.0d+0/Xt)-(1.0d+0/Xc)
    F2=(1.0d+0/Yt)-(1.0d+0/Yc)
    F3=(1.0d+0/Zt)-(1.0d+0/Zc)
    F11=(1.0d+0/(Xt*Xc))
    F22=(1.0d+0/(Yt*Yc))
    F33=(1.0d+0/(Zt*Zc))
    F44=(1.0d+0/(S13*S13))
    F55=(1.0d+0/(S23*S23))
    F66=(1.0d+0/(S12*S12))
    F12=-0.5d+0*(1.0d+0/SQRT(Xt*Xc*Yt*Yc))
    F13=-0.5d+0*(1.0d+0/SQRT(Xt*Xc*Zt*Zc))
    F23=-0.5d+0*(1.0d+00/SQRT(Yt*Yc*Zt*Zc))
    
! Compute Tsai-Wu polynomials 

    do I= 1,nblock   ! Loop through all material points

        Xphi(1)=(F1*Stress(I,1))+(F11*Stress(I,1)*Stress(I,1))&
        +(F12*Stress(I,1)*Stress(I,2))+(F13*Stress(I,1)*Stress(I,3))

        Xphi(2)=(F2*Stress(I,2))+(F22*Stress(I,2)*Stress(I,2))&
        +(F12*Stress(I,1)*Stress(I,2))+(F23*Stress(I,2)*Stress(I,3))

        Xphi(3)=(F3*Stress(I,3))+(F33*Stress(I,3)*Stress(I,3))&
        +(F13*Stress(I,1)*Stress(I,3))+(F23*Stress(I,2)*Stress(I,3))

        Xphi(4)=(F66*Stress(I,4)*Stress(I,4))

        Xphi(5)=(F44*Stress(I,5)*Stress(I,5))

        Xphi(6)=(F55*Stress(I,6)*Stress(I,6))

! Obtain and store the index of maximum contributing polynomial

        IMAX = maxloc(Xphi, dim =1)
        
        Phi(I) =Xphi(1)+Xphi(2)+Xphi(3)+Xphi(4)+Xphi(5)+Xphi(6)

        e(I,IMAX) = Phi(I)

    end do

    case(4)
!#-------------------Hoffman Criteria--------------------#

! Calculate Hoffman Strength parameters

    F1=(1.0d+0/Xt)-(1.0d+0/Xc)
    F2=(1.0d+0/Yt)-(1.0d+0/Yc)
    F3=(1.0d+0/Zt)-(1.0d+0/Zc)
    F11=(1.0d+0/(Xt*Xc))
    F22=(1.0d+0/(Yt*Yc))
    F33=(1.0d+0/(Zt*Zc))
    F44=(1.0d+0/(S13*S13))
    F55=(1.0d+0/(S23*S23))
    F66=(1.0d+0/(S12*S12))
    F12=-0.50d+0*((1.0d+0/(Xt*Xc))+(1.0d+0/(Yt*Yc))-(1.0d+0/(Zt*Zc)))
    F13=-0.50d+0*((1.0d+0/(Xt*Xc))+(1.0d+0/(Zt*Zc))-(1.0d+0/(Yt*Yc)))
    F23=-0.50d+0*((1.0d+0/(Zt*Zc))+(1.0d+0/(Yt*Yc))-(1.0d+0/(Xt*Xc)))

! Compute Hoffman polynomials 

    do I= 1,nblock   ! Loop through all material points

        Xphi(1)=(F1*Stress(I,1))+(F11*Stress(I,1)*Stress(I,1))&
        +(F12*Stress(I,1)*Stress(I,2))+(F13*Stress(I,1)*Stress(I,3))

        Xphi(2)=(F2*Stress(I,2))+(F22*Stress(I,2)*Stress(I,2))&
        +(F12*Stress(I,1)*Stress(I,2))+(F23*Stress(I,2)*Stress(I,3))

        Xphi(3)=(F3*Stress(I,3))+(F33*Stress(I,3)*Stress(I,3))&
        +(F13*Stress(I,1)*Stress(I,3))+(F23*Stress(I,2)*Stress(I,3))

        Xphi(4)=(F66*Stress(I,4)*Stress(I,4))

        Xphi(5)=(F44*Stress(I,5)*Stress(I,5))

        Xphi(6)=(F55*Stress(I,6)*Stress(I,6))

! Obtain and store the index of maximum contributing polynomial
        
        IMAX = maxloc(Xphi, dim=1)

        Phi(I) =Xphi(1)+Xphi(2)+Xphi(3)+Xphi(4)+Xphi(5)+Xphi(6)   

        e(I,IMAX) = Phi(I)      

    end do


    case(5)
!#-------------------Hashin Criteria--------------------#

    do I= 1,nblock   ! Loop through all material points and calculate failure indicies

! Tensile/Compressive Fiber Failure

    if(Stress(I,1).GE.0.0) then
        e(I,1)= SQRT(((Stress(I,1)/Xt)**2) + (((Stress(I,4)**2)+(Stress(I,5)**2))/(S12)))
    else
        e(I,1)= SQRT(((Stress(I,1)/Xc)**2))
    end if            

! Tensile/Compressive Matrix Failure

    if(Stress(I,2)+Stress(I,3).GE.0.0) then 

        e(I,2)= SQRT((((Stress(I,2)+Stress(I,3))**2)/(Yt)**2)+(((Stress(I,6)**2)&
        -(Stress(I,2)*Stress(I,3)))/(S23**2))+(((Stress(I,4)**2)+(Stress(I,5)**2))/(S12**2)))
    else                
        e(I,2)= SQRT((((Yc/(2*S23))**2)-1)*((Stress(I,2)+Stress(I,3))/(Yc))&
        +((Stress(I,2)+Stress(I,3))**2/(4*S23**2))&
        +((Stress(I,6)**2-Stress(I,2)*Stress(I,3))/(S23**2))&
        +((Stress(I,4)**2+Stress(I,5)**2)/(S12**2)))
    end if

! Interlaminar normal Tensile/Compressive failure      

    if(Stress(I,3).GE.0.0) then 
        e(I,3)=SQRT(((Stress(I,3))/(Zt))**2)
    else
        e(I,3)=SQRT(((Stress(I,3))/(Zc))**2)
    end if              

    end do

    case(6)
!#-------------------Hashin-Rotem Criteria--------------------# 

    do I= 1,nblock   ! Loop through all material points and calculate failure indicies

! Tensile/Compressive Fiber Failure

    if(Strain(I,1).GE.0.0) then
        e(I,1)= (Strain(I,1)/(EpsXt))**2
    else
        e(I,1)= (Strain(I,1)/(EpsXc))**2
    end if            

! In-Plane Matrix Tensile/Compressive Failure 

    if(Strain(I,2).GE.0.0) then 
        e(I,2)= (Strain(I,2)/(EpsYt))**2 +(Strain(I,6)/(GamS23))**2 +(Strain(I,4)/(GamS12))**2
    else
        e(I,2)= (Strain(I,2)/(EpsYc))**2 +(Strain(I,6)/(GamS23))**2 +(Strain(I,4)/(GamS12))**2
    end if             

! Transverse Matrix Tensile/Compressive Failure 

    if(Strain(I,3).GE.0.0) then 

        e(I,3)= (Strain(I,3)/(EpsZt))**2 +(Strain(I,6)/(GamS23))**2 +(Strain(I,5)/(GamS13))**2
    else
        e(I,3)= (Strain(I,3)/(EpsZc))**2 +(Strain(I,6)/(GamS23))**2 +(Strain(I,5)/(GamS13))**2
    end if

    end do


    case(7) 
!#-------------------Puck Criteria--------------------#

    do I= 1,nblock   ! Loop through all material points 

! Fiber fracture Tensile/Compressive
           
        A = Stress(I,1)-(ANU12-ANU12f*MGF*(E11/E11f))*(Stress(I,2)+Stress(I,3))
                
    if(A.GE.0.0) then
            
        e(I,1)=(1.0d+0/Xt)*A
                 
    else if(A.LT.0.0) then 
                 
        e(I,1)=(1.0d+0/Xc)*A
              
    end if
              
! Inter-Fiber fracture 

        ! Step-wise search for failure plane 
                            
        do J= 0,180 ! loop from -90 to 90 degrees
            
            THETA = -Pi/2.0d+0 + J*(Pi/180)  
            
            ! Stresses on fracture plane          
            
            SFP =(Stress(I,2)*cos(THETA)**2)&
            &+(Stress(I,3)*sin(THETA)**2)+(2*Stress(I,6)*sin(THETA)*cos(THETA))
            
            TNT =-Stress(I,2)*sin(THETA)*cos(THETA)+Stress(I,3)*sin(THETA)*cos(THETA)+&
            &Stress(I,6)*(cos(THETA)**2-sin(THETA)**2)
                            
            TN1 =Stress(I,4)*cos(THETA)+Stress(I,5)*sin(THETA)
            
            ! Puck Strength Parameters
                   
            Rn= Yt
            
            Rn1 = S12
            
            Rnt = Yc/(2.0d+0*tan(THETAF))
                
            ! Puck Slope Parameters                
            
            Pnt = -1.0d+0/(2.0d+0*tan(2*THETAF))
            
            Pn1 = Pnt*(Rn1/Rnt)        
           
        if(SFP.GE.0.0) then
            
            IFF = (SFP/Rn)**2 +(TN1/(Rn1-Pn1*SFP))**2 + (TNT/(Rnt-Pnt*SFP))**2                       
            
        else if(SFP.LT.0.0) then
            
            IFF = (TN1/(Rn1-Pn1*SFP))**2 +(TNT/(Rnt-Pnt*SFP))**2
                
        end if
        
                                        
        if(IFF.GT.e(I,2)) then
            
            e(I,2)=IFF
            
            THETAMAX=THETA  
     
        end if 
            
        end do
                          
    end do

    case default 

! Print error message and exit the program for invalid selection

    print*, "Invalid Failure ID"

    call exit (0)

    end select

    end subroutine failure_check


!###########  Utility Subroutine 3: damage_check  ###########
!#----------------------------------------------------------#

! The following subroutine :

! Evaluates the damage variables(dmg)based on selected damage model

    subroutine damage_check (id)

    implicit none

    integer :: id

    select case (id)

    case(1) 
!#----------------Instantaneous Degredation-----------------#


    do I= 1,nblock ! Loop through all material points

! Fiber Tensile/Compressive damage

    if((e(I,1).GT.1.0) .AND. (Stress(I,1).GE.0.0)) then  

        fflags(I,1)=step    

        dmg(I,1)= (0.99999999-beta_ft)  

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(I,4)= (0.99999999-beta_s)

        end if
               

    else if((e(I,1).GT.1.0) .AND. (Stress(I,1).LT.0.0)) then  

        fflags(I,1)=-step  

        dmg(I,1)= (0.99999999-beta_fc)  

    end if

! Matrix Tensile/Compressive damage

    if((e(I,2).GT.1.0) .AND. (Stress(I,2).GE.0.0)) then  

        fflags(I,2)=step    

        dmg(I,2)= (0.99999999-beta_mt)  

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(I,4)= (0.99999999-beta_s)
            dmg(I,6)= (0.99999999-beta_s)

        end if
        
        if (failure_id.EQ.7) then ! Additional induced damage for Puck

            dmg(I,3)= (0.99999999-beta_mt)
            dmg(I,4)= (0.99999999-beta_s)
            dmg(I,6)= (0.99999999-beta_s)

        end if        
        
        

    else if((e(I,2).GT.1.0) .AND. (Stress(I,2).LT.0.0)) then  

        fflags(I,2)=-step  

        dmg(I,2)= (0.99999999-beta_mc)  

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(I,4)= (0.99999999-beta_s)
            dmg(I,6)= (0.99999999-beta_s)

        end if
        
        if (failure_id.EQ.7) then ! Additional induced damage for Puck

            dmg(I,3)= (0.99999999-beta_mc)
            dmg(I,4)= (0.99999999-beta_s)
            dmg(I,6)= (0.99999999-beta_s)

        end if                  

    end if

! Interlaminar Tensile/Compressive damage

    if((e(I,3).GT.1.0) .AND. (Stress(I,3).GE.0.0)) then  

        fflags(I,3)=step    

        dmg(I,3)= (0.99999999-beta_mt)  

    else if((e(I,3).GT.1.0) .AND. (Stress(I,2).LT.0.0)) then  

        fflags(I,3)=-step  

        dmg(I,3)= (0.99999999-beta_mc)  

    end if
         
! Shear damage

    do J=4,6    

    if((e(I,J).GT.1.0)) THEN  

        fflags(I,J)=step  

        dmg(I,J)= (0.99999999-beta_s)  

    end if                  

    end do                  

    end do

    case(2)
!#----------------Recursive Degredation-----------------#

    do I= 1,nblock ! Loop through all material points

! Fiber Tensile/Compressive damage

    if(((e(I,1).GT.1.0).OR.(dmg(I,1).LT.1.0)) .AND. (Stress(I,1).GE.0.0)) then  

        fflags(I,1)=step    

        dmg(I,1)= dmg(I,1)*(0.99999999-beta_ft)  

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(I,4)= dmg(I,4)*(0.99999999-beta_s)

        end if        

    else if(((e(I,1).GT.1.0).OR.(dmg(I,1).LT.1.0)) .AND. (Stress(I,1).LT.0.0)) then  

        fflags(I,1)=-step  

        dmg(I,1)= dmg(I,1)*(0.99999999-beta_fc)  

    end if

! Matrix Tensile/Compressive damage

    if(((e(I,2).GT.1.0).OR.(dmg(I,2).LT.1.0)) .AND. (Stress(I,2).GE.0.0)) then 

        fflags(I,2)=step    

        dmg(I,2)= dmg(I,2)*(0.99999999-beta_mt)  

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(I,4)= dmg(I,4)*(0.99999999-beta_s)
            dmg(I,6)= dmg(I,6)*(0.99999999-beta_s)

        end if  
        
        if (failure_id.EQ.7) then ! Additional induced damage for Puck

            dmg(I,3)= dmg(I,3)*(0.99999999-beta_mt)
            dmg(I,4)= dmg(I,4)*(0.99999999-beta_s)
            dmg(I,6)= dmg(I,6)*(0.99999999-beta_s)

        end if

    else if(((e(I,2).GT.1.0).OR.(dmg(I,2).LT.1.0)) .AND. (Stress(I,2).LT.0.0)) then  

        fflags(I,2)=-step  

        dmg(I,2)= dmg(I,2)*(0.99999999-beta_mc) 

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(I,4)= dmg(I,4)*(0.99999999-beta_s)
            dmg(I,6)= dmg(I,6)*(0.99999999-beta_s)

        end if 
        
        if (failure_id.EQ.7) then ! Additional induced damage for Puck

            dmg(I,3)= dmg(I,3)*(0.99999999-beta_mc)
            dmg(I,4)= dmg(I,4)*(0.99999999-beta_s)
            dmg(I,6)= dmg(I,6)*(0.99999999-beta_s)

        end if              

    end if

! Interlaminar Tensile/Compressive damage

    if(((e(I,3).GT.1.0).OR.(dmg(I,3).LT.1.0)) .AND. (Stress(I,3).GE.0.0)) then  

        fflags(I,3)=step   

        dmg(I,3)= dmg(I,3)*(0.99999999-beta_mt) 

    else if(((e(I,3).GT.1.0).OR.(dmg(I,3).LT.1.0)) .AND. (Stress(I,3).LT.0.0)) then  

        fflags(I,3)=-step  

        dmg(I,3)= dmg(I,3)*(0.99999999-beta_mc)  

    end if     

    do J=4,6    

    if(((e(I,J).GT.1.0).OR.(dmg(I,J).LT.1.0))) THEN 

        fflags(I,J)=step  

        dmg(I,J)= dmg(I,J)*(0.99999999-beta_s) 

    end if                  

    end do                  

    end do


    case(3)
!#----------------Exponential Degredation-----------------#

    do I= 1,nblock ! Loop through all material points

! Fiber Tensile/Compressive damage

    if(((e(I,1).GT.1.0).OR.(dmg(I,1).LT.1.0)) .AND. (Strain(I,1).GE.0.0)) then 

        fflags(I,1)=step    

        dmg(I,1)= EXP(-a_ft*(Strain(I,1)- EpsXt)/(n_ft*EpsXt)) 

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(I,4)= EXP(-a_s*(Strain(I,4)- GamS12)/(n_s*GamS12))

        end if  

    else if(((e(I,1).GT.1.0).OR.(dmg(I,1).LT.1.0)) .AND. (Strain(I,1).LT.0.0)) then 

        fflags(I,1)=-step  

        dmg(I,1)= EXP(-a_fc*(abs(Strain(I,1))- EpsXc)/(n_fc*EpsXc))  

    end if

! Matrix Tensile/Compressive damage

    if(((e(I,2).GT.1.0).OR.(dmg(I,2).LT.1.0)) .AND. (Strain(I,2).GE.0.0)) then  

        fflags(I,2)=step    

        dmg(I,2)= EXP(-a_mt*(Strain(I,2)- EpsYt)/(n_mt*EpsYt))   

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(I,4)= EXP(-a_s*(Strain(I,4)- GamS12)/(n_s*GamS12))

            dmg(I,6)= EXP(-a_s*(Strain(I,6)- GamS23)/(n_s*GamS23))           

        end if
        
        if (failure_id.EQ.7) then ! Additional induced damage for Puck

            dmg(I,3)= EXP(-a_mt*(Strain(I,3)- EpsZt)/(n_mt*EpsZt))
            dmg(I,4)= EXP(-a_s*(Strain(I,4)- GamS12)/(n_s*GamS12))
            dmg(I,6)= EXP(-a_s*(Strain(I,6)- GamS23)/(n_s*GamS23))

        end if                 


    else if(((e(I,2).GT.1.0).OR.(dmg(I,2).LT.1.0)) .AND. (Strain(I,2).LT.0.0)) then  

        fflags(I,2)=-step  

        dmg(I,2)= EXP(-a_mc*(abs(Strain(I,2))- EpsYc)/(n_mc*EpsYc))  

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(I,4)= EXP(-a_s*(Strain(I,4)- GamS12)/(n_s*GamS12))

            dmg(I,6)= EXP(-a_s*(Strain(I,6)- GamS23)/(n_s*GamS23))           

        end if 
        
        if (failure_id.EQ.7) then ! Additional induced damage for Puck

            dmg(I,3)= EXP(-a_mc*abs((Strain(I,3))- EpsZc)/(n_mc*EpsZc)) 
            dmg(I,4)= EXP(-a_s*(Strain(I,4)- GamS12)/(n_s*GamS12))
            dmg(I,6)= EXP(-a_s*(Strain(I,6)- GamS23)/(n_s*GamS23))

        end if                

    end if

! Interlaminar Tensile/Compressive damage

    if(((e(I,3).GT.1.0).OR.(dmg(I,3).LT.1.0)) .AND. (Stress(I,3).GE.0.0)) then  

        fflags(I,3)=step    

        dmg(I,3)= EXP(-a_mt*(Strain(I,3)- EpsZt)/(n_mt*EpsZt)) 

    else if(((e(I,3).GT.1.0).OR.(dmg(I,3).LT.1.0)) .AND. (Stress(I,3).LT.0.0)) then  

        fflags(I,3)=-step  

        dmg(I,3)= EXP(-a_mc*abs((Strain(I,3))- EpsZc)/(n_mc*EpsZc))  

    end if     

! Shear damage


    if((e(I,4).GT.1.0).OR.(dmg(I,4).LT.1.0)) then  

        fflags(I,4)=step  

        dmg(I,4)= EXP(-a_s*abs((Strain(I,4))- GamS12)/(n_s*GamS12)) 

    end if                  

    if(((e(I,5).GT.1.0).OR.(dmg(I,5).LT.1.0))) then  

        fflags(I,5)=step  

        dmg(I,5)= EXP(-a_s*abs((Strain(I,5))- GamS13)/(n_s*GamS13))  

    end if

    if((e(I,6).GT.1.0).OR.(dmg(I,6).LT.1.0)) then  

        fflags(I,6)=step  

        dmg(I,6)= EXP(-a_s*abs((Strain(I,6))- GamS23)/(n_s*GamS23))  

    end if                  

    end do  


    case(4) 
!#--------------Constant Stress Degredation---------------#

    do I= 1,nblock ! Loop through all material points

! Fiber Tensile/Compressive damage

    if((e(I,1).GT.1.0) .AND. (Stress(I,1).GE.0.0)) then  

        fflags(I,1)=step    

        dmg(I,1)= dmg(I,1)*(1.0d+0/e(I,1))  

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(I,4)= dmg(I,4)*(1.0d+0/e(I,1)) 

        end if

    else if((e(I,1).GT.1.0) .AND. (Stress(I,1).LT.0.0)) then  

        fflags(I,1)=-step  

        dmg(I,1)= dmg(I,1)*(1.0d+0/e(I,1))  

    end if

! Matrix Tensile/Compressive damage

    if((e(I,2).GT.1.0) .AND. (Stress(I,2).GE.0.0)) then  

        fflags(I,2)=step    

        dmg(I,2)= dmg(I,2)*(1.0d+0/e(I,2))  

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(I,4)= dmg(I,4)*(1.0d+0/e(I,2))
            dmg(I,6)= dmg(I,6)*(1.0d+0/e(I,2))

        end if
        
        if (failure_id.EQ.7) then ! Additional induced damage for Puck

            dmg(I,3)=dmg(I,3)*(1.0d+0/e(I,2)) 
            dmg(I,4)=dmg(I,4)*(1.0d+0/e(I,2))
            dmg(I,6)=dmg(I,6)*(1.0d+0/e(I,2))

        end if         
        

    else if((e(I,2).GT.1.0) .AND. (Stress(I,2).LT.0.0)) then  

            fflags(I,2)=-step  

            dmg(I,2)= dmg(I,2)*(1.0d+0/e(I,2))  

            if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

                dmg(I,4)= dmg(I,4)*(1.0d+0/e(I,2))
                dmg(I,6)= dmg(I,6)*(1.0d+0/e(I,2))

            end if
 
         if (failure_id.EQ.7) then ! Additional induced damage for Puck

            dmg(I,3)=dmg(I,3)*(1.0d+0/e(I,2)) 
            dmg(I,4)=dmg(I,4)*(1.0d+0/e(I,2))
            dmg(I,6)=dmg(I,6)*(1.0d+0/e(I,2))

        end if           
            

    end if

! Interlaminar Tensile/Compressive damage

    if((e(I,3).GT.1.0) .AND. (Stress(I,3).GE.0.0)) then  

        fflags(I,3)=step    

        dmg(I,3)= dmg(I,3)*(1.0d+0/e(I,3))  

    else if((e(I,3).GT.1.0) .AND. (Stress(I,2).LT.0.0)) then  

        fflags(I,3)=-step  

        dmg(I,3)= dmg(I,3)*(1.0d+0/e(I,3))

    end if 
    
! Shear damage

    do J=4,6    

    if((e(I,J).GT.1.0)) then  

        fflags(I,J)=step  

        dmg(I,J)= dmg(I,J)*(1.0d+0/e(I,J))  

    end if                  

    end do                  
    
    end do

    case(5) 
!#------Continium damage mechanics: Crack Band Theory-----#

    do I= 1,nblock ! Loop through all material points

! Fiber Tensile/Compressive damage

    if((e(I,1).GT.1.0) .AND. (Strain(I,1).GE.0.0)) then 

        fflags(I,1)=step    

! Calculate Degraded E11 modulus

        E11D =((1.0/E11)+(Strain(I,1)-EpsXt)/(Xt*(1.0-(le*Xt*(Strain(I,1)- EpsXt))/(2*G_ft))))**(-1)
        d_index(I,1) = 1.0d+0 -(E11D/E11) 

    else if((e(I,1).GT.1.0) .AND. (Strain(I,1).LT.0.0)) then  

        fflags(I,1)=-step                     

! Calculate Degraded E11 modulus 

        E11D =((1.0/E11)+(abs(Strain(I,1))-EpsXc))/(Xc*(1.0-(le*Xc*((abs(Strain(I,1)-EpsXc)))/(2*G_fc))))**(-1)
        d_index(I,1) = 1.0d+0 - (E11D/E11) 

    end if  
            

! Matrix Tensile/Compressive damage

    if((e(I,2).GT.1.0) .AND. (Strain(I,2).GE.0.0)) then  

        fflags(I,2)=step   

! Calculate Degraded E22,G12,G23 modulI                    

        E22D =((1.0/E22)+(Strain(I,2)-EpsYt)/(Yt*(1.0-(le*Yt*(Strain(I,2)-EpsYt))/(2*G_mt))))**(-1)

        G12D =((1.0/G12)+(abs(Strain(I,4))-GamS12))/(2*S12*(1.0-(le*S12*((abs(Strain(I,4)-GamS12)))/(4*G_IIC))))**(-1)

        G23D =((1.0/G23)+(abs(Strain(I,6))-GamS23))/(2*S23*(1.0-(let*S23*((abs(Strain(I,6)-GamS23)))/(4*G_IIC))))**(-1)
                  
        d_index(I,2) = 1.0d+0 -(E22D/E22)
        d_index(I,4) = 1.0d+0 -(G12D/G12)                                         
        d_index(I,6) = 1.0d+0 -(G23D/G23)  


    else if((e(I,2).GT.1.0) .AND. (Strain(I,2).LT.0.0)) then  

        fflags(I,2)=-step                      

! Calculate Degraded E22,G12,G23 modulI                    

        E22D =((1.0/E22)+(abs(Strain(I,2))-EpsYc)/(Yc*(1.0-(le*Yc*((abs(Strain(I,2))-EpsYc)))/(2*G_mc))))**(-1)

        G12D =((1.0/G12)+(abs(Strain(I,4))-GamS12)/(2*S12*(1.0-(le*S12*((abs(Strain(I,4))-GamS12)))/(4*G_IIC))))**(-1)

        G23D =((1.0/G23)+(abs(Strain(I,6))-GamS23)/(2*S23*(1.0-(let*S23*((abs(Strain(I,6))-GamS23)))/(4*G_IIC))))**(-1)
                  
        d_index(I,2) = 1.0d+0 -(E22D/E22)
        d_index(I,4) = 1.0d+0 -(G12D/G12)                                         
        d_index(I,6) = 1.0d+0 -(G23D/G23)

    end if 

! Transverse Matrix Damage Modes

    if((e(I,3).GT.1.0) .AND. (Strain(I,3).GE.0.0)) then  

        fflags(I,3)=step                      

! Calculate Degraded (E33,G23,G13) moduli 

        E33D =((1.0/E33)+(Strain(I,3)-EpsZt)/(Zt*(1.0-(let*Zt*(Strain(I,3)-EpsZt))/(2*G_mt))))**(-1)

        G23D =((1.0/G23)+((abs(Strain(I,6))-GamS23))/(2*S23*(1.0-(let*S23*((abs(Strain(I,6))-GamS23)))/(4*G_IIC))))**(-1)

        G13D =((1.0/G13)+((abs(Strain(I,5))-GamS13))/(2*S13*(1.0-(let*S13*((abs(Strain(I,5))-GamS13)))/(4*G_IIC))))**(-1)                    

        d_index(I,3) = 1.0d+0 -(E33D/E33)
        d_index(I,5) = 1.0d+0 -(G13D/G13)                                         
        d_index(I,6) = 1.0d+0 -(G23D/G23)                  

    else if((e(I,3).GT.1.0) .AND. (Strain(I,3).LT.0.0)) then  

        fflags(I,3)=-step   

! Calculate Degraded (E33,G23,G13) moduli

        E33D =((1.0/E33)+(abs(Strain(I,3))-EpsZc)/(Zc*(1.0-(let*Zc*(Strain(I,3)-EpsZc))/(2*G_mc))))**(-1)

        G23D =((1.0/G23)+((abs(Strain(I,6))-GamS23))/(2*S23*(1.0-(let*S23*((abs(Strain(I,6))-GamS23)))/(4*G_IIC))))**(-1)

        G13D =((1.0/G13)+((abs(Strain(I,5))-GamS13))/(2*S13*(1.0-(let*S13*((abs(Strain(I,5))-GamS13)))/(4*G_IIC))))**(-1)                    

        d_index(I,3) = 1.0d+0 -(E33D/E33)
        d_index(I,5) = 1.0d+0 -(G13D/G13)                                      
        d_index(I,6) = 1.0d+0 -(G23D/G23)   

    end if 

! Assign Constiutive Matrix Damage Variables 

    dmg(I,1) = abs(0.999999-(d_index(I,1)))
    dmg(I,2) = abs(0.999999-(d_index(I,2)))
    dmg(I,3) = abs(0.999999-(d_index(I,3)))               
    dmg(I,4) = abs(0.999999-((1.0d+0-Dmg(I,1))*(1.0d+0-alpha*Dmg(I,2))*(1.0d+0-beta*d_Index(I,4))))
    dmg(I,5) = abs(0.999999-((1.0d+0-Dmg(I,1))*(1.0d+0-alpha*Dmg(I,3))*(1.0d+0-beta*d_Index(I,5))))
    dmg(I,6) = abs(0.999999-((1.0d+0-Dmg(I,1))*(1.0d+0-alpha*Dmg(I,2))*(1.0d+0-alpha*Dmg(I,3))*(1.0-beta*d_Index(I,5))))

    end do

    case default 

! Print error message and exit the program for invalid selection

    print*, "Invalid damage ID"

    call exit (0)

    end select

    end subroutine damage_check

    end module PDA
