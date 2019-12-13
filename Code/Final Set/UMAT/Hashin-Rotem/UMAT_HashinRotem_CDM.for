C*************************************************************************
! Technische Universität München
! Chair of Computational Mechanics
! Software Lab Project: Development of Failure Criteria for Composites
C*************************************************************************
! 3D ORTHOTROPIC ELATICITY WITH HASHIN-ROTEM CRITERION AND CONTINIUM DAMAGE MECHANICS
! DAMAGE MODEL: CRACK BAND THEORY
! (CAN NOT BE USED FOR 2D PROBLEMS)
! VERSION: 1.4
C-------------------------------------------------------------------------
C   GROUP: 11
C	WRITTEN BY AMMAR KHALLOUF
C   DATE: 06.12.2019
C   
C=========================================================================
!	HEADER OF THE SUBROUTINE

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
	 

      INCLUDE 'ABA_PARAM.INC'
	  
!     WARNING - the aba_param.inc file declares
!        Implicit real*8(a-h,o-z)
!     This means that, by default, any variables with
!     first letter between a-h or o-z are double precision.
!     The rest are integers.
!     Note that this also means that if you type a variable
!     name incorrectly, the compiler won't catch your typo.	  

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4),UPSTRAN(NTENS),dmg(6),e(6)
	 
	  INTEGER fflags(6), Del_Stat
      DOUBLE PRECISION le,let 
      
!DIR$ FREEFORM     
        
!###########           Read UMAT Input            ###########
!#----------------------------------------------------------#         

!                ****Orthotropic Material Constants****
! +-----------+-------------------------------------------------------------------+
! | Parameter |                            Description                            |
! +-----------+-------------------------------------------------------------------+
! | E11       | Longitudinal elastic Young's modulus(Fiber: 11 axis)              |
! | E22       | Transverse elastic Young's modulus   (Matrix: 22 axis)            |
! | E33       | Interlaminar elastic Young's modulus (Interlaminar: 33 axis)      |
! | G12       | Elastic Shear Modulus  (1-2 plane)                                |
! | G13       | Elastic Shear Modulus  (1-3 plane)                                |
! | G23       | Elastic Shear Modulus  (2-3 plane)                                |
! | ANU12     | Poisson ratio   (1-2 plane)                                       |
! | ANU13     | Poisson ratio   (1-3 plane)                                       |
! | ANU23     | Poisson ration  (2-3 plane)                                       |
! | Xt        | Laminate tensile strength (11 fiber direction)                    |
! | Xc        | Laminate compressive strength  (11 fiber direction)               |
! | Yt        | Laminate tensile strength  (22 matrix direction)                  |
! | Yc        | Laminate compressive strength (22 matrix direction)               |
! | Zt        | Laminate interlaminar tensile strength (33 direction)             |
! | Zc        | Laminate interlaminar compressive strength (33 direction)         |
! | S12       | Laminate Shear Strength  (1-2 plane)                              |
! | S13       | Laminate Shear Strength  (1-3 plane)                              |
! | S23       | Laminate Shear Strength  (2-3 plane)                              |
! | EpsXt     | Laminate allowable tensile strain (11 fiber direction)            |
! | EpsXc     | Laminate allowable compressive strain  (11 fiber direction)       |
! | EpsYt     | Laminate allowable tensile strain (Matrix: 22 axis)               |
! | EpsYc     | Laminate allowable compressive strain (Matrix: 22 axis)           |
! | EpsZt     | Laminate allowable interlaminar tensile strain (33 direction)     |
! | EpsZc     | Laminate allowable interlaminar compressive strain (33 direction) |
! | GamS12    | Laminate allowable shear Strain  (1-2 plane)                      |
! | GamS13    | Laminate allowable shear Strain  (1-3 plane)                      |
! | GamS23    | Laminate allowable shear Strain  (2-3 plane)                      |
! +-----------+-------------------------------------------------------------------+

! Important Remarks: 

! *Engineering shear strain convention is used 

! *In FEA, the material local axis need to be assigned with repspect to the   
!  above convention in order to establish the correct material properties orienatation

! *The elastic constants need to satisfy the stability limits for an orthotropic 
!  material law as detailed in the project documentation

    E11=PROPS(1)
    E22=PROPS(2)
    E33=PROPS(3)
    ANU12=PROPS(4)
    ANU13=PROPS(5)
    ANU23=PROPS(6)   
    G12=PROPS(7)
    G13=PROPS(8)
    G23=PROPS(9)
                  
    Xt=PROPS(10)
    Xc=PROPS(11)
    Yt=PROPS(12)
    Yc=PROPS(13)
    Zt=PROPS(14)
    Zc=PROPS(15)
    S12=PROPS(16)
    S13=PROPS(17)
    S23=PROPS(18)
      
    EpsXt=PROPS(19)
    EpsXc=PROPS(20)
    EpsYt=PROPS(21)
    EpsYc=PROPS(22)
    EpsZt=PROPS(23)      
    EpsZc=PROPS(24)  
    GamS12=PROPS(25)
    GamS13=PROPS(26)
    GamS23=PROPS(27)      
 
!        ****Continium Damage Mechanics Parameters****
! +-----------+-----------------------------------------+
! | Parameter |               Description               |
! +-----------+-----------------------------------------+
! | G_ft      | fracture energy (fiber tension)         |
! | G_fc      | fracture energy (fiber compression)     |
! | G_mt      | fracture energy (matrix tension)        |
! | G_mc      | fracture energy (matrix compression)    |
! | G_IIC     | fracture energy (Mode II shear failure) |
! | le        | characteristic element length           |
! | let       | characteristic element thickness        |
! | alpha,beta| Nonlinear Shear degradation factors     |
! +-----------+-----------------------------------------+

! Important Remarks: 

!*Degradation factors values are in the range of (0,1) 
!*Degredation factor == 0 implies no damage
!*Degredation factor == 1 implies full damage (i.e stiffness approx zero)

 
    G_ft=PROPS(28)
    G_fc=PROPS(29)
    G_mt=PROPS(30) 
    G_mc=PROPS(31)
    G_IC=PROPS(32)
    G_IIC=PROPS(33)
    le=PROPS(34)
    let=PROPS(35)
    alpha=PROPS(36)
    beta=PROPS(37)



!###########      Start Progressive Analysis      ###########
!#----------------------------------------------------------#             

 ! Assign state/solution dependent variables:

 ! includes: damage factors, failure flags, and element deletion status  

    dmg(1:6)=STATEV(1:6)
                
    fflags(1:6)=STATEV(7:12)
                
    Del_Stat=STATEV(13) ! Control element deletion status

! At Step 1 the material is undamaged (linear elastic)

! At Step 1 the failure flags and indicies == 0
 
if(KINC.EQ.1) then
  
    dmg=1.0   
    e =0.0   
    fflags=0      
    Del_Stat = 1
        
! Initial Strain increment
      
    do I=1,NTENS
    
        UPSTRAN(I)= STRAN(I)+DSTRAN(I)
       
    end do
   
!###########      Calculate Initial Stresses      ###########
!#----------------------------------------------------------#

    call ortho3d(E11,E22,E33,G12,G13,G23,ANU12,ANU13,ANU23,dmg,STRESS,UPSTRAN,DDSDDE)

!###########    Return Initial State Variables    ###########
!#----------------------------------------------------------# 

    STATEV(1:6)=dmg(1:6)  
    STATEV(7:12)=fflags(1:6) 
    STATEV(13)=Del_Stat

else
		  
!###########         Update Total UPSTRAN          ###########
!#----------------------------------------------------------#
      
! Accumulate Strains

    do I=1,NTENS   
        UPSTRAN(I)= STRAN(I)+DSTRAN(I)       
    end do
    
! Read State variables from previous increment

    dmg(1:6)=STATEV(1:6)
                
    fflags(1:6)=STATEV(7:12)
                
    Del_Stat=STATEV(13) ! Control element deletion status

!###########     Calculate Trial Stresses         ###########
!#----------------------------------------------------------#

    call ortho3d(E11,E22,E33,G12,G13,G23,ANU12,ANU13,ANU23,dmg,STRESS,UPSTRAN,DDSDDE)
       
!###########    Failure Evaluation: Hashin-Rotem  ###########
!#----------------------------------------------------------#       	  
   
    call failure_calc(STRESS,UPSTRAN,e,EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23)    
          
!###########        Damage Evaluation: CDM        ###########
!#----------------------------------------------------------#
    
    if ((maxval(e).GT.1.0) .OR. (minval(dmg).LT.1.0)) then

    call damage_calc (STRESS,UPSTRAN,dmg,e,fflags,KINC,G_ft,G_fc,G_mt,G_mc,G_IIC,le,let,alpha,beta&
    &,EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23,E11,E22,E33,G12,G13,G23,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23)   
    
    do I=1,6 

    if ((dmg(I).GT.1.0) .OR. (dmg(I).LT.0.0)) then

        write(*,*) "Encountered Negative or Greater than one Damage Multipliers"
    
        write(*,*) "Please Check User Material Input Variables"

        call XIT
    
    end if
    
    end do
    
    end if    
         
!###########           Update Stresses            ###########
!#----------------------------------------------------------#

    call ortho3d(E11,E22,E33,G12,G13,G23,ANU12,ANU13,ANU23,dmg,STRESS,UPSTRAN,DDSDDE)           
                        
!###########       Update State Variables         ###########
!#----------------------------------------------------------# 

! Delete Failed elements with excessive Stress/Strength ratio (default=2)

    if(maxval(e) .GT. 2.0) then    
        Del_Stat= 0
    end if
    
    STATEV(1:6)=dmg(1:6)  
    STATEV(7:12)=fflags(1:6) 
    STATEV(13)=Del_Stat

!###########      End Progressive Analysis        ###########
!#----------------------------------------------------------#         

end if	  	  
	        
      RETURN
      END
      

!###########    Utility Subroutine 1: ortho3d     ###########
!#----------------------------------------------------------#
        
    subroutine ortho3D (E11,E22,E33,G12,G13,G23,ANU12,ANU13,ANU23,dmg,STRESS,UPSTRAN,DDSDDE)  
       
    implicit none
    integer :: I,J
    real*8, dimension (6) :: dmg,STRESS,UPSTRAN
    real*8, dimension (6,6) :: DDSDDE
    real*8 :: E11,E22,E33,G12,G13,G23,ANU12,ANU13,ANU23,ANU21,ANU31,ANU32,Ypsilon
    
    !  NOTE: ABAQUS uses engineering shear strains,
    !  i.e. stran(ndi+1) = 2*e_12, etc...

    ANU21=(E22*ANU12)/E11
    ANU31=(E33*ANU13)/E11
    ANU32=(E33*ANU23)/E22
    
    Ypsilon=(1.0d+0)/(1-dmg(1)*dmg(2)*ANU12*ANU21&
    &-dmg(2)*dmg(3)*ANU23*ANU32&
    &-dmg(1)*dmg(3)*ANU13*ANU31&
    &-2.0d+0*dmg(1)*dmg(2)*dmg(3)*ANU21*ANU32*ANU13)
                 
    ! Form effective orthotropic material matrix including damage (Jacobian Matrix) 
    
    DDSDDE(1,1)=(E11*Ypsilon*(1.0-dmg(2)*dmg(3)*ANU23*ANU32))*dmg(1)
    
    DDSDDE(2,2)=(E22*Ypsilon*(1.0-dmg(1)*dmg(3)*ANU13*ANU31))*dmg(2)
    
    DDSDDE(3,3)=(E33*Ypsilon*(1.0-dmg(1)*dmg(2)*ANU12*ANU21))*dmg(3)
                    
    DDSDDE(1,2)=(E11*Ypsilon*(ANU21+ANU31*ANU23*dmg(3)))*dmg(1)*dmg(2)
                               
    DDSDDE(1,3)=(E11*Ypsilon*(ANU31+ANU21*ANU32*dmg(2)))*dmg(1)*dmg(3)
                                                                    
    DDSDDE(2,3)=(E22*Ypsilon*(ANU32+dmg(1)*ANU12*ANU31))*dmg(3)*dmg(2)
    
    DDSDDE(2,1)=DDSDDE(1,2)
        
    DDSDDE(3,1)=DDSDDE(1,3) 
        
    DDSDDE(3,2)=DDSDDE(2,3)
                                
    DDSDDE(4,4)=G12*dmg(4)
        
    DDSDDE(5,5)=G13*dmg(5)
         
    DDSDDE(6,6)=G23*dmg(6)

    ! Calculate Stress Components
    
    do I=1,6    
        STRESS(I)=0.0d+0
                       
        do J=1,6

        STRESS(I)=STRESS(I)+DDSDDE(I,J)*UPSTRAN(J)
        
    end do
        end do
   
    end subroutine ortho3d

!###########  Utility Subroutine 2: failure_calc  ###########
!#----------------------------------------------------------#

    subroutine failure_calc(STRESS,UPSTRAN,e,EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23)
    
    implicit none
    real*8, dimension  (6) :: e,STRESS,UPSTRAN
    real*8 :: EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23

!#-------------------Hashin-Rotem Criteria--------------------#

! Tensile/Compressive Fiber Failure

    if(UPSTRAN(1).GE.0.0) then
        e(1)=(UPSTRAN(1)/(EpsXt))**2
    else
        e(1)=(UPSTRAN(1)/(EpsXc))**2
    end if            

! In-Plane Matrix Tensile/Compressive Failure 

    if(UPSTRAN(2).GE.0.0) then 
        e(2)=(UPSTRAN(2)/(EpsYt))**2 +(UPSTRAN(6)/(GamS23))**2 +(UPSTRAN(4)/(GamS12))**2
    else
        e(2)=(UPSTRAN(2)/(EpsYc))**2 +(UPSTRAN(6)/(GamS23))**2 +(UPSTRAN(4)/(GamS12))**2
    end if             

! Transverse Matrix Tensile/Compressive Failure 

    if(UPSTRAN(3).GE.0.0) then 
        e(3)=(UPSTRAN(3)/(EpsZt))**2 +(UPSTRAN(6)/(GamS23))**2 +(UPSTRAN(5)/(GamS13))**2
    else
        e(3)=(UPSTRAN(3)/(EpsZc))**2 +(UPSTRAN(6)/(GamS23))**2 +(UPSTRAN(5)/(GamS13))**2
    end if
      
    end subroutine failure_calc    
    
!###########   Utility Subroutine 2: damage_calc  ###########
!#----------------------------------------------------------#
 
    subroutine damage_calc (STRESS,UPSTRAN,dmg,e,fflags,KINC,G_ft,G_fc,G_mt,G_mc,G_IIC,le,let,alpha,beta&
    &,EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23,E11,E22,E33,G12,G13,G23,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23)    
    
    implicit none
    real*8, dimension  (6) :: e,STRESS,UPSTRAN,dmg,d_index
    real*8 :: G_ft,G_fc,G_mt,G_mc,G_IIC,le,let,alpha,beta,E11D,E22D,E33D,G12D,G13D,G23D
    real*8 :: EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23
    real*8 :: E11,E22,E33,G12,G13,G23,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23    
    integer:: I,KINC,fflags(6)
    
!#------Continium damage mechanics: Crack Band Theory-----#

! Fiber Tensile/Compressive damage

    if((e(1).GT.1.0) .AND. (UPSTRAN(1).GE.0.0)) then  

        fflags(1)=KINC   

        ! Calculate Degraded E11 modulus 

        E11D =((1.0/E11)+(UPSTRAN(1)- EpsXt)/(Xt*(1.0-(le*Xt*(UPSTRAN(1)- EpsXt))/(2*G_ft))))**(-1)
        d_index(1) = 1.0d+0 -(E11D/E11) 

    else if((e(1).GT.1.0) .AND. (UPSTRAN(1).LT.0.0)) then  

        fflags(1)=-KINC                   

        ! Calculate Degraded E11 modulus as per Eq (12)

        E11D =((1.0/E11)+(abs(UPSTRAN(1))-EpsXc))/(Xc*(1.0-(le*Xc*((abs(UPSTRAN(1)-EpsXc)))/(2*G_fc))))**(-1)
        d_index(1) = 1.0d+0 - (E11D/E11) 

    end if          

! Matrix Tensile/Compressive damage

    if((e(2).GT.1.0) .AND. (UPSTRAN(2).GE.0.0)) then  

        fflags(2)=KINC    

        ! Calculate Degraded E22,G12,G23 modulI                   

        E22D =((1.0/E22)+(UPSTRAN(2)-EpsYt)/(Yt*(1.0-(le*Yt*(UPSTRAN(2)-EpsYt))/(2*G_mt))))**(-1)

        G12D =((1.0/G12)+(abs(UPSTRAN(4))-GamS12))/&
        &(2*S12*(1.0-(le*S12*((abs(UPSTRAN(4)-GamS12)))/(4*G_IIC))))**(-1)

        G23D =((1.0/G23)+(abs(UPSTRAN(6))-GamS23))/&
        &(2*S23*(1.0-(let*S23*((abs(UPSTRAN(6)-GamS23)))/(4*G_IIC))))**(-1)                   

        d_index(2) = 1.0d+0 -(E22D/E22)
        d_index(4) = 1.0d+0 -(G12D/G12)                                      
        d_index(6) = 1.0d+0 -(G23D/G23)  


    else if((e(2).GT.1.0) .AND. (UPSTRAN(2).LT.0.0)) then

        fflags(2)=-KINC                      

        ! Calculate Degraded E22,G12,G23 modulI                    

        E22D =((1.0/E22)+(abs(UPSTRAN(2))-EpsYc)/(Yc*(1.0-(le*Yc*((abs(UPSTRAN(2))-EpsYc)))/(2*G_mc))))**(-1)

        G12D =((1.0/G12)+(abs(UPSTRAN(4))-GamS12)/&
        &(2*S12*(1.0-(le*S12*((abs(UPSTRAN(4))-GamS12)))/(4*G_IIC))))**(-1)

        G23D =((1.0/G23)+(abs(UPSTRAN(6))-GamS23)/&
        &(2*S23*(1.0-(let*S23*((abs(UPSTRAN(6))-GamS23)))/(4*G_IIC))))**(-1)                   

        d_index(2) = 1.0d+0 -(E22D/E22)
        d_index(4) = 1.0d+0 -(G12D/G12)                                        
        d_index(6) = 1.0d+0 -(G23D/G23)

    end if 

! Interlaminar Tensile/Compressive damage

    if((e(3).GT.1.0) .AND. (UPSTRAN(3).GE.0.0)) then 

        fflags(3)=KINC                

        ! Calculate Degraded (E33,G23,G13) moduli

        E33D =((1.0/E33)+(UPSTRAN(3)-EpsZt)/(Zt*(1.0-(let*Zt*(UPSTRAN(3)-EpsZt))/(2*G_mt))))**(-1)

        G23D =((1.0/G23)+((abs(UPSTRAN(6))-GamS23))/&
        &(2*S23*(1.0-(let*S23*((abs(UPSTRAN(6))-GamS23)))/(4*G_IIC))))**(-1)

        G13D =((1.0/G13)+((abs(UPSTRAN(5))-GamS13))/&
        &(2*S13*(1.0-(let*S13*((abs(UPSTRAN(5))-GamS13)))/(4*G_IIC))))**(-1)                    

        d_index(3) = 1.0d+0 -(E33D/E33)
        d_index(5) = 1.0d+0 -(G13D/G13)                                        
        d_index(6) = 1.0d+0 -(G23D/G23)                 

    else if((e(3).GT.1.0) .AND. (UPSTRAN(3).LT.0.0)) then  

        fflags(3)=-KINC    

        ! Calculate Degraded (E33,G23,G13) moduli 

        E33D =((1.0/E33)+(abs(UPSTRAN(3))-EpsZc)/(Zc*(1.0-(let*Zc*(UPSTRAN(3)-EpsZc))/(2*G_mc))))**(-1)

        G23D =((1.0/G23)+((abs(UPSTRAN(6))-GamS23))/&
        &(2*S23*(1.0-(let*S23*((abs(UPSTRAN(6))-GamS23)))/(4*G_IIC))))**(-1)

        G13D =((1.0/G13)+((abs(UPSTRAN(5))-GamS13))/&
        &(2*S13*(1.0-(let*S13*((abs(UPSTRAN(5))-GamS13)))/(4*G_IIC))))**(-1)                    

        d_index(3) = 1.0d+0 -(E33D/E33)
        d_index(5) = 1.0d+0 -(G13D/G13)                                        
        d_index(6) = 1.0d+0 -(G23D/G23) 

    end if 

! Assign Constiutive Matrix Damage Variables 

    dmg(1) = abs(0.999999-(d_index(1)))
    dmg(2) = abs(0.999999-(d_index(2)))
    dmg(3) = abs(0.999999-(d_index(3)))               
    dmg(4) = abs(0.999999-((1.0d+0-dmg(1))*(1.0d+0-alpha*dmg(2))*(1.0d+0-beta*d_index(4))))
    dmg(5) = abs(0.999999-((1.0d+0-dmg(1))*(1.0d+0-alpha*dmg(3))*(1.0d+0-beta*d_index(5))))
    dmg(6) = abs(0.999999-((1.0d+0-dmg(1))*(1.0d+0-alpha*dmg(2))*(1.0d+0-alpha*dmg(3))*(1.0d+0-beta*d_index(5))))     

    end subroutine damage_calc  