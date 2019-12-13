C*************************************************************************
! Technische Universität München
! Chair of Computational Mechanics
! Software Lab Project: Development of Failure Criteria for Composites
C*************************************************************************
! PDALAC UMAT Subroutine (Progressive damage analysis of Laminated Composites)
! (Applicable only for 3D Elements)
! VERSION: 1.0
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
     
	  INTEGER fflags(6), Del_Stat,failure_id,damage_id
      DOUBLE PRECISION MGF,le,let
      
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
        
!               ****Progressive Analysis Variables****
! +-------------+-----------------------------------------------------+
! |  Parameter  |                     Description                     |
! +-------------+-----------------------------------------------------+
! | analysis_id | Index of analysis option set to be used             |
! | nblock      | Number of material integration points               |
! | nsteps      | Number of load steps                                |
! | StrainInc   | Strain increments vector                            |
! | THETAF      | Maximum fracture angle for Puck in radians          |
! | MGF         | Magnification factor for Puck                       |
! | E11F        | Fiber elastic modulus for Puck                      |
! | ANU12F      | Fiber poisson ratio for Puck                        |
! | failure_id  | Index for failure criteria selection:               |
! |             | --->1-Max Stress                                    |
! |             | --->2-Max Strain                                    |
! |             | --->3-Tsai-Wu                                       |
! |             | --->4-Hoffman                                       |
! |             | --->5-Hashin                                        |
! |             | --->6-Hashin-Rotem                                  |
! |             | --->7-Puck                                          |
! | damage_id   | Index for damage model selection:                   |
! |             | --->1-Ply Discount:Instantaneous Degradation        |
! |             | --->2-Ply Discount:Recursive Degradation            |
! |             | --->3-Ply Discount:Exponential Degradation          |
! |             | --->4-Ply Discount:Constant Stress Degradation      |
! |             | --->5-Continium damage mechanics: Crack Band Theory |
! +-------------+-----------------------------------------------------+

!                  ****Ply Discount Damage Parameters****
! +------------+----------------------------------------------+
! | Parameter  |               Description                    |
! +------------+----------------------------------------------+
! +-----------------------------------------------------------+
! |            Instant/Recursive Damage                       |
! +-----------------------------------------------------------+
! |  beta_ft   | degradation factor (fiber tension)           |
! |  beta_fc   | degradation factor (fiber compression)       |
! |  beta_mt   | degradation factor (matrix tension)          |
! |  beta_mc   | degradation factor (matrix compression)      |
! |  beta_s    | degradation factor (shear)                   |
! +-----------------------------------------------------------+
! |            Exponential Damage                             |
! +-----------------------------------------------------------+
! | a_ft, n_ft | Exp. degradation factors (fiber tension)     |
! | a_fc, n_fc | Exp. degradation factor (fiber compression)  |
! | a_mt, n_mt | Exp. degradation factor (matrix tension)     |
! | a_mc, n_mc | Exp. degradation factor (matrix compression) |
! | a_s,  n_s  | Exp. degradation factor (shear)              |
! +------------+----------------------------------------------+


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

    failure_id=PROPS(28)         
    damage_id= PROPS(29)
    beta_ft=PROPS(30) 
    beta_fc=PROPS(31)     
    beta_mt=PROPS(32)   
    beta_mc=PROPS(33)   
    beta_s=PROPS(34)    
    a_ft=PROPS(35)
    n_ft=PROPS(36)
    a_fc=PROPS(37)
    n_fc=PROPS(38)
    a_mt=PROPS(39)
    n_mt=(40)
    a_mc=(41)
    n_mc=PROPS(42)
    a_s=PROPS(43)
    n_s=PROPS(44)
    G_ft=PROPS(45)
    G_fc=PROPS(46)
    G_mt=PROPS(47) 
    G_mc=PROPS(48)
    G_IC=PROPS(49)
    G_IIC=PROPS(50)
    le=PROPS(51)
    let=PROPS(52)
    alpha=PROPS(53)
    beta=PROPS(54)
    THETAF=PROPS(55)
    MGF=PROPS(56)
    E11F=PROPS(57)
    ANU12F=PROPS(58)
    
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

    call ortho3D (E11,E22,E33,G12,G13,G23,ANU12,ANU13,ANU23,dmg,STRESS,UPSTRAN,DDSDDE,damage_id) 


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

    call ortho3D (E11,E22,E33,G12,G13,G23,ANU12,ANU13,ANU23,dmg,STRESS,UPSTRAN,DDSDDE,damage_id) 
 
!###########           Failure Evaluation         ###########
!#----------------------------------------------------------#

! call failure_calc subroutine to calculate failure indicies
 
    call failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23,&
    &EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23,THETAF,MGF,ANU12,ANU12f,E11,E11f)
    
!###########         Damage Evaluation            ###########
!#----------------------------------------------------------#

! Apply damage if any of the failure indicies excceds one
! or previous damage presist      

   if ((maxval(e).GT.1.0) .OR. (minval(dmg).LT.1.0)) then

        call damage_calc (failure_id,damage_id,STRESS,UPSTRAN,dmg,e,fflags,KINC,beta_ft,beta_fc,beta_mt,beta_mc,beta_s,&
    &EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23,a_ft,a_fc,a_mt,a_mc,a_s,n_ft,n_fc,n_mt,n_mc,n_s,&
    &G_ft,G_fc,G_mt,G_mc,G_IIC,le,let,alpha,beta,E11,E22,E33,G12,G13,G23,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23)  
    
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

    call ortho3D (E11,E22,E33,G12,G13,G23,ANU12,ANU13,ANU23,dmg,STRESS,UPSTRAN,DDSDDE,damage_id)     
    
    
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

endif
                    
    RETURN
    END
    
!###########  Utility Subroutines implemnetation  ###########
!#----------------------------------------------------------#

! The following section include the numerical implementation
! for the ortho3d,failure_calc,damage_calc subrotuines

! Ortho3d subroutine: form the effective orthotropic elasticity
! matrix based on selected damage model and calculated damage 
! it returns the updated stresses 

! failure_calc subroutine: returns the six failure indicies
! based on selected failure criteria

! damage_calc subroutine: returns the six damage indicies
! based on selected damage model

!###########    Utility Subroutine 1: ortho3d     ###########
!#----------------------------------------------------------#
        
    subroutine ortho3D (E11,E22,E33,G12,G13,G23,ANU12,ANU13,ANU23,dmg,STRESS,UPSTRAN,DDSDDE,damage_id) 
    implicit none
    integer :: I,J,damage_id
    real*8, dimension (6) :: dmg,STRESS,UPSTRAN
    real*8, dimension (6,6) :: DDSDDE
    real*8 :: E11,E22,E33,G12,G13,G23,ANU12,ANU13,ANU23,ANU21,ANU31,ANU32,Ypsilon
    
    !  NOTE: ABAQUS uses engineering shear strains,
    !  i.e. stran(ndi+1) = 2*e_12, etc...

    select case (damage_id)
    
    case (5)  ! effective stress for continium damage mechanics
    
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
                    
    case default  ! effective stress for ply-discount damage
    
    !  NOTE: ABAQUS uses engineering shear strains,
    !  i.e. stran(ndi+1) = 2*e_12, etc...    
    
    ANU21=(E22*ANU12)/E11
    ANU31=(E33*ANU13)/E11
    ANU32=(E33*ANU23)/E22
    Ypsilon=(1.0d+0)/(1.0d+0-ANU12*ANU21-ANU23*ANU32-ANU13*ANU31-2.0d+0*ANU21*ANU32*ANU13)
    
    ! Form effective orthotropic material matrix including damage (Jacobian Matrix) 
    
    DDSDDE(1,1)=(E11*Ypsilon*(1.0d+0-ANU23*ANU32))*abs(dmg(1))
    
    DDSDDE(2,2)=(E22*Ypsilon*(1.0d+0-ANU13*ANU31))*abs(dmg(2))
    
    DDSDDE(3,3)=(E33*Ypsilon*(1.0d+0-ANU12*ANU21))*abs(dmg(3))       
    
    DDSDDE(1,2)=(E11*Ypsilon*(ANU21+ANU31*ANU23))*abs(dmg(1)*dmg(2))
    
    DDSDDE(1,3)=(E11*Ypsilon*(ANU31+ANU21*ANU32))*abs(dmg(1)*dmg(3))
    
    DDSDDE(2,3)=(E22*Ypsilon*(ANU32+ANU12*ANU31))*abs(dmg(3)*dmg(2))       
    
    DDSDDE(2,1)=DDSDDE(1,2)
                
    DDSDDE(3,1)=DDSDDE(1,3)
    
    DDSDDE(3,2)=DDSDDE(2,3)
    
    DDSDDE(4,4)=G12*abs(dmg(4))
    
    DDSDDE(5,5)=G13*abs(dmg(5))
    
    DDSDDE(6,6)=G23*abs(dmg(6))
    
    ! Calculate Stress Components
    
    do I=1,6   
        STRESS(I)=0.0d+0
                       
        do J=1,6

        STRESS(I)=STRESS(I)+DDSDDE(I,J)*UPSTRAN(J)

        end do
    end do      
    
    end select

    end subroutine ortho3d


!###########  Utility Subroutine 2: failure_calc  ###########
!#----------------------------------------------------------#
      
    subroutine failure_calc(failure_id,STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23,&
    &EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23,THETAF,MGF,ANU12,ANU12f,E11,E11f)
    
    implicit none
    
    real*8, dimension  (6) :: e,STRESS,UPSTRAN,Xphi
    real*8 :: Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23,Phi,E11,THETAF,MGF,ANU12,ANU12f,E11f    
    real*8 :: EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23
    real*8 :: F1,F2,F3,F11,F22,F33,F44,F55,F66,F12,F13,F23
    real*8 :: THETA,Rn,Rn1,Rnt,SFP,TN1,TNT,Pnt,Pn1,THETAMAX,IFF,A
    real*8, parameter :: Pi = 3.14159265358979d+0
    integer :: failure_id,IMAX,I,J

    
    select case (failure_id)

    case(1)
!#-------------------Max Stress Criteria--------------------#

! Failure along fiber (11 axis) (Tension/Compression)

    if(STRESS(1).GE.0.0) then
        e(1)=STRESS(1)/Xt
    else
        e(1)=-STRESS(1)/Xc
    end if 

! Failure along matrix (22 axis) (Tension/Compression)

    if(STRESS(2).GE.0.0) then 
        e(2)=STRESS(2)/Yt 
    else
        e(2)=-STRESS(2)/Yc 
    end if

! Interlaminar failure (33 axis) (Tension/Compression)                   

    if(STRESS(3).GE.0.0) then 
        e(3)=STRESS(3)/Zt   
    else         
        e(3)=-STRESS(3)/Zc  
    end if

! Shear Failure (1-2 plane)
 
        e(4)=abs(STRESS(4))/S12

! Shear Failure (1-3 plane)

        e(5)=abs(STRESS(5))/S13

! Shear Failure (2-3 plane)

        e(6)=abs(STRESS(6))/S23 
        
    case(2)
!#-------------------Max Strain Criteria--------------------#
 
! Failure along fiber (11 axis) (Tension/Compression)

    if (UPSTRAN(1).GE.0.0) then
        e(1)=UPSTRAN(1)/EpsXt
    else
        e(1)=-UPSTRAN(1)/EpsXc
    end if 

! Failure along matrix (22 axis) (Tension/Compression)

    if (UPSTRAN(2).GE.0.0) then 
        e(2)=UPSTRAN(2)/EpsYt 
    else
        e(2)=-UPSTRAN(2)/EpsYc  
    end if

! Interlaminar failure (33 axis) (Tension/Compression)                 

    if (UPSTRAN(3).GE.0.0) then 
        e(3)=UPSTRAN(3)/EpsZt
    else         
        e(3)=-UPSTRAN(3)/EpsZc
    end if

! Shear Failure (1-2 plane) 

        e(4)=abs(UPSTRAN(4))/GamS12

! Shear Failure (1-3 plane)

        e(5)=abs(UPSTRAN(5))/GamS13

! Shear Failure (2-3 plane)

        e(6)=abs(UPSTRAN(6))/GamS23

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
    F12=-0.50d+0*(1.0d+0/SQRT(Xt*Xc*Yt*Yc))
    F13=-0.50d+0*(1.0d+0/SQRT(Xt*Xc*Zt*Zc))
    F23=-0.50d+0*(1.0d+0/SQRT(Yt*Yc*Zt*Zc))
    
! Compute Tsai-Wu polynomials 

    Xphi(1)=(F1*Stress(1))+(F11*Stress(1)*Stress(1))&
    +(F12*Stress(1)*Stress(2))+(F13*Stress(1)*Stress(3))

    Xphi(2)=(F2*Stress(2))+(F22*Stress(2)*Stress(2))&
    +(F12*Stress(1)*Stress(2))+(F23*Stress(2)*Stress(3))

    Xphi(3)=(F3*Stress(3))+(F33*Stress(3)*Stress(3))&
    +(F13*Stress(1)*Stress(3))+(F23*Stress(2)*Stress(3))

    Xphi(4)=(F66*Stress(4)*Stress(4))
    
    Xphi(5)=(F44*Stress(5)*Stress(5))
    
    Xphi(6)=(F55*Stress(6)*Stress(6))        

    Phi =Xphi(1)+Xphi(2)+Xphi(3)+Xphi(4)+Xphi(5)+Xphi(6)

! Obtain and store the index of maximum contributing polynomial

    IMAX = maxloc(Xphi, dim =1)

    e(IMAX) = Phi

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
    F12=-0.5d+0*((1.0d+0/(Xt*Xc))+(1.0d+0/(Yt*Yc))-(1.0d+0/(Zt*Zc)))
    F13=-0.50d+0*((1.0d+0/(Xt*Xc))+(1.0d+0/(Zt*Zc))-(1.0d+0/(Yt*Yc)))
    F23=-0.50d+0*((1.0d+0/(Zt*Zc))+(1.0d+0/(Yt*Yc))-(1.0d+0/(Xt*Xc)))

! Compute Hoffman polynomials 

    Xphi(1)=(F1*Stress(1))+(F11*Stress(1)*Stress(1))&
    +(F12*Stress(1)*Stress(2))+(F13*Stress(1)*Stress(3))
    
    Xphi(2)=(F2*Stress(2))+(F22*Stress(2)*Stress(2))&
    +(F12*Stress(1)*Stress(2))+(F23*Stress(2)*Stress(3))
    
    Xphi(3)=(F3*Stress(3))+(F33*Stress(3)*Stress(3))&
    +(F13*Stress(1)*Stress(3))+(F23*Stress(2)*Stress(3))
    
    Xphi(4)=(F66*Stress(4)*Stress(4))
    
    Xphi(5)=(F44*Stress(5)*Stress(5))
    
    Xphi(6)=(F55*Stress(6)*Stress(6))
    
    Phi =Xphi(1)+Xphi(2)+Xphi(3)+Xphi(4)+Xphi(5)+Xphi(6)       

! Obtain and store the index of maximum contributing polynomial

    IMAX = maxloc(Xphi, dim =1)

    e(IMAX) = Phi

    case(5)
!#-------------------Hashin Criteria--------------------#

! Tensile/Compressive Fiber Failure

    if(STRESS(1).GE.0.0) then
        e(1)= SQRT(((STRESS(1)/Xt)**2)+(((STRESS(4)**2)+(STRESS(5)**2))/(S12)))
    else
        e(1)=SQRT(((STRESS(1)/Xc)**2))
    end if            

! Tensile/Compressive Matrix Failure

    if(STRESS(2)+STRESS(3).GE.0.0) then 
        e(2)= SQRT((((STRESS(2)+STRESS(3))**2)/(Yt)**2)+(((STRESS(6)**2)&
        -(STRESS(2)*STRESS(3)))/(S23**2))+(((STRESS(4)**2)+(STRESS(5)**2))/(S12**2)))
    else                
        e(2)=SQRT((((Yc/(2*S23))**2)-1)*((STRESS(2)+STRESS(3))/(Yc))&
        +((STRESS(2)+STRESS(3))**2/(4*S23**2))&
        +((STRESS(6)**2-STRESS(2)*STRESS(3))/(S23**2))&
        +((STRESS(4)**2+STRESS(5)**2)/(S12**2)))
    end if

! Interlaminar normal Tensile/Compressive failure       

    if(STRESS(3).GE.0.0) then 
        e(3)=SQRT(((STRESS(3))/(Zt))**2)
    else
        e(3)=SQRT(((STRESS(3))/(Zc))**2)
    end if              

    case(6) 
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

    case(7) 
!#-------------------Puck Criteria--------------------#

! Fiber fracture Tensile/Compressive
           
    A = STRESS(1)-(ANU12-ANU12f*MGF*(E11/E11f))*(Stress(2)+Stress(3))
                
    if(A.GE.0.0) then
            
        e(1)=(1.0d+0/Xt)*A
                 
    else if(A.LT.0.0) then 
                 
        e(1)=(1.0d+0/Xc)*A
              
    end if
              
! Inter-Fiber fracture 

    ! Step-wise search for failure plane 
                            
    do J= 0,180 ! loop from -90 to 90 degrees
            
        THETA = -Pi/2.0d+0 + J*(Pi/180)  
        
        ! Stresses on fracture plane          
        
        SFP =(Stress(2)*cos(THETA)**2)&
        &+(Stress(3)*sin(THETA)**2)+(2*Stress(6)*sin(THETA)*cos(THETA))
        
        TNT =-Stress(2)*sin(THETA)*cos(THETA)+Stress(3)*sin(THETA)*cos(THETA)+&
        &Stress(6)*(cos(THETA)**2-sin(THETA)**2)
                        
        TN1 =Stress(4)*cos(THETA)+Stress(5)*sin(THETA)
        
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
                                        
    if(IFF.GT.e(2)) then
            
        e(2)=IFF            
        THETAMAX=THETA  
     
    end if 
            
    end do
                         
    case default 

! Print error message and exit the program for invalid selection

    write(*,*) "Invalid Failure ID"

    call XIT

    end select

    end subroutine failure_calc


!###########   Utility Subroutine 2: damage_calc  ###########
!#----------------------------------------------------------#
 
    subroutine damage_calc (failure_id,damage_id,STRESS,UPSTRAN,dmg,e,fflags,KINC,beta_ft,beta_fc,beta_mt,beta_mc,beta_s,&
    &EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23,a_ft,a_fc,a_mt,a_mc,a_s,n_ft,n_fc,n_mt,n_mc,n_s,&
    &G_ft,G_fc,G_mt,G_mc,G_IIC,le,let,alpha,beta,E11,E22,E33,G12,G13,G23,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23)    
    
    implicit none
    real*8, dimension  (6) :: e,STRESS,UPSTRAN,dmg,d_index
    real*8 :: beta_ft,beta_fc,beta_mc,beta_mt,beta_s
    real*8 :: a_ft,a_fc,a_mt,a_mc,a_s,n_ft,n_fc,n_mt,n_mc,n_s
    real*8 :: EpsXt,EpsXc,EpsYt,EpsYc,EpsZt,EpsZc,GamS12,GamS13,GamS23
    real*8 :: E11,E22,E33,G12,G13,G23,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23    
    real*8 :: G_ft,G_fc,G_mt,G_mc,G_IIC,le,let,alpha,beta,E11D,E22D,E33D,G12D,G13D,G23D    
    integer:: damage_id,failure_id,I,KINC,fflags(6)
    
    select case (damage_id)

    case(1) 
!#----------------Instantaneous Degredation-----------------#

! Fiber Tensile/Compressive damage

    if((e(1).GT.1.0) .AND. (STRESS(1).GE.0.0)) then
     
        fflags(1)=KINC   
        dmg(1)=(0.99999999-beta_ft)

        if(failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)=(0.99999999-beta_s)

        end if
                  
    else if((e(1).GT.1.0) .AND. (STRESS(1).LT.0.0)) then  

        fflags(1)=-KINC  
        dmg(1)=(0.99999999-beta_fc)  

    end if

! Matrix Tensile/Compressive damage

    if((e(2).GT.1.0) .AND. (STRESS(2).GE.0.0)) then  

        fflags(2)=KINC    
        dmg(2)= (0.99999999-beta_mt)  

        if(failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)=(0.99999999-beta_s)
            dmg(6)=(0.99999999-beta_s)

        end if 
        
        if(failure_id.EQ.7) then ! Additional induced shear damage for Puck

            dmg(3)= (0.99999999-beta_mt)
            dmg(4)= (0.99999999-beta_s)
            dmg(6)= (0.99999999-beta_s)
            
        end if           
        
    else if((e(2).GT.1.0) .AND. (STRESS(2).LT.0.0)) then  

        fflags(2)=-KINC  
        dmg(2)= (0.99999999-beta_mc)  

        if(failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)=(0.99999999-beta_s)
            dmg(6)=(0.99999999-beta_s)

        end if 

        if(failure_id.EQ.7) then ! Additional induced shear damage for Puck

            dmg(3)= (0.99999999-beta_mc)
            dmg(4)= (0.99999999-beta_s)
            dmg(6)= (0.99999999-beta_s)
            
        end if    

    end if

! Interlaminar Tensile/Compressive damage

    if((e(3).GT.1.0) .AND. (STRESS(3).GE.0.0)) then  

        fflags(3)=KINC    
        dmg(3)=(0.99999999-beta_mt) 

    else if((e(3).GT.1.0) .AND. (STRESS(2).LT.0.0)) then  

        fflags(3)=-KINC  
        dmg(3)=(0.99999999-beta_mc)  

    end if     

! Shear damage

    do I=4,6 

    if(e(I).GT.1.0) then 

        fflags(I)=KINC  
        dmg(I)= (0.99999999-beta_s)  

    end if
                      
    end do                  


    case(2)

!#----------------Recursive Degredation-----------------#

! Fiber Tensile/Compressive damage

    if(((e(1).GT.1.0).OR.(dmg(1).LT.1.0)) .AND. (STRESS(1).GE.0.0)) then 

        fflags(1)=KINC    
        dmg(1)= dmg(1)*(0.99999999d+0-beta_ft) 

        
        if(failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)= dmg(4)*(0.99999999-beta_s)

        end if   

    else if(((e(1).GT.1.0).OR.(dmg(1).LT.1.0)) .AND. (STRESS(1).LT.0.0)) then  

        fflags(1)=-KINC 
        dmg(1)= dmg(1)*(0.99999999d+0-beta_fc)  

    end if

! Matrix Tensile/Compressive damage

    if(((e(2).GT.1.0).OR.(dmg(2).LT.1.0)) .AND. (STRESS(2).GE.0.0)) then  
    
        fflags(2)=KINC   
        dmg(2)=dmg(2)*(0.99999999d+0-beta_mt)

        if(failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)=dmg(4)*(0.99999999-beta_s)
            dmg(6)=dmg(6)*(0.99999999-beta_s)

        end if 

        if(failure_id.EQ.7) then ! Additional induced shear damage for Puck

            dmg(3)= dmg(3)*(0.99999999-beta_mt)
            dmg(4)= dmg(4)*(0.99999999-beta_s)
            dmg(6)= dmg(6)*(0.99999999-beta_s)
            
        end if         

    else if(((e(2).GT.1.0).OR.(dmg(2).LT.1.0)) .AND. (STRESS(2).LT.0.0)) then 

        fflags(2)=-KINC  
        dmg(2)=dmg(2)*(0.99999999d+0-beta_mc)  
 
         if(failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)=dmg(4)*(0.99999999-beta_s)
            dmg(6)=dmg(6)*(0.99999999-beta_s)

        end if  

        if(failure_id.EQ.7) then ! Additional induced shear damage for Puck

            dmg(3)= dmg(3)*(0.99999999-beta_mc)
            dmg(4)= dmg(4)*(0.99999999-beta_s)
            dmg(6)= dmg(6)*(0.99999999-beta_s)
            
        end if
 
    end if

! Interlaminar Tensile/Compressive damage

    if(((e(3).GT.1.0).OR.(dmg(3).LT.1.0)) .AND. (STRESS(3).GE.0.0)) then  

        fflags(3)=KINC    
        dmg(3)= dmg(3)*(0.99999999d+0-beta_mt)  

    else if(((e(3).GT.1.0).OR.(dmg(3).LT.1.0)) .AND. (STRESS(3).LT.0.0)) then  

        fflags(3)=-KINC  
        dmg(3)= dmg(3)*(0.99999999d+0-beta_mc)  

    end if     

! Shear damage

    do I=4,6 

    if((e(I).GT.1.0).OR.(dmg(I).LT.1.0)) THEN

        fflags(I)=KINC  
        dmg(I)=dmg(I)*(0.99999999d+0-beta_s)  

    end if
                    
    end do

    case(3) 
    
!#----------------Exponential Degredation-----------------#

! Fiber Tensile/Compressive damage

    if(((e(1).GT.1.0).OR.(dmg(1).LT.1.0)) .AND. (UPSTRAN(1).GE.0.0)) then  

        fflags(1)=KINC    
        dmg(1)= EXP(-a_ft*(UPSTRAN(1)- EpsXt)/(n_ft*EpsXt))

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)= EXP(-a_ft*(UPSTRAN(1)- EpsXt)/(n_ft*EpsXt))

        end if

   
    else if(((e(1).GT.1.0).OR.(dmg(1).LT.1.0)) .AND. (UPSTRAN(1).LT.0.0)) then  

        fflags(1)=-KINC  
        dmg(1)= EXP(-a_fc*(abs(UPSTRAN(1))- EpsXc)/(n_fc*EpsXc))  


    end if

! Matrix Tensile/Compressive damage

    if(((e(2).GT.1.0).OR.(dmg(2).LT.1.0)) .AND. (UPSTRAN(2).GE.0.0)) then  

        fflags(2)=KINC    
        dmg(2)= EXP(-a_mt*(UPSTRAN(2)- EpsYt)/(n_mt*EpsYt))   

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)= EXP(-a_mt*(UPSTRAN(2)- EpsYt)/(n_mt*EpsYt))
            dmg(6)= EXP(-a_mt*(UPSTRAN(2)- EpsYt)/(n_mt*EpsYt))       

        end if 
 
        if(failure_id.EQ.7) then ! Additional induced shear damage for Puck

            dmg(3)= EXP(-a_mt*(UPSTRAN(2)- EpsYt)/(n_mt*EpsYt))
            dmg(4)= EXP(-a_mt*(UPSTRAN(2)- EpsYt)/(n_mt*EpsYt))
            dmg(6)= EXP(-a_mt*(UPSTRAN(2)- EpsYt)/(n_mt*EpsYt))
            
        end if  
 
    else if(((e(2).GT.1.0).OR.(dmg(2).LT.1.0)) .AND. (UPSTRAN(2).LT.0.0)) then  

        fflags(2)=-KINC 
        dmg(2)= EXP(-a_mc*(abs(UPSTRAN(2))- EpsYc)/(n_mc*EpsYc))  
          
        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)= EXP(-a_mc*(abs(UPSTRAN(2))- EpsYc)/(n_mc*EpsYc))
            dmg(6)= EXP(-a_mc*(abs(UPSTRAN(2))- EpsYc)/(n_mc*EpsYc))       

        end if 

        if(failure_id.EQ.7) then ! Additional induced shear damage for Puck

            dmg(3)= EXP(-a_mc*(abs(UPSTRAN(2))- EpsYc)/(n_mc*EpsYc))
            dmg(4)= EXP(-a_mc*(abs(UPSTRAN(2))- EpsYc)/(n_mc*EpsYc))
            dmg(6)= EXP(-a_mc*(abs(UPSTRAN(2))- EpsYc)/(n_mc*EpsYc))
            
        end if  

    end if

! Interlaminar Tensile/Compressive damage

    if(((e(3).GT.1.0).OR.(dmg(3).LT.1.0)) .AND. (UPSTRAN(3).GE.0.0)) then  

        fflags(3)=KINC    
        dmg(3)= EXP(-a_mt*(UPSTRAN(3)- EpsZt)/(n_mt*EpsZt)) 

    else if(((e(3).GT.1.0).OR.(dmg(3).LT.1.0)) .AND. (UPSTRAN(3).LT.0.0)) then  

        fflags(3)=-KINC  
        dmg(3)= EXP(-a_mc*abs((UPSTRAN(3))- EpsZc)/(n_mc*EpsZc))  

    end if     

! Shear damage

    if((e(4).GT.1.0).OR.(dmg(4).LT.1.0)) then  

        fflags(4)=KINC  
        dmg(4)= EXP(-a_s*abs((UPSTRAN(4))- GamS12)/(n_s*GamS12))  

    end if                  

    if((e(5).GT.1.0).OR.(dmg(5).LT.1.0)) then  

        fflags(5)=KINC  
        dmg(5)= EXP(-a_s*abs((UPSTRAN(5))- GamS13)/(n_s*GamS13))  

    end if

    if((e(6).GT.1.0).OR.(dmg(6).LT.1.0)) then 

        fflags(6)=KINC  
        dmg(6)= EXP(-a_s*abs((UPSTRAN(6))- GamS23)/(n_s*GamS23))  

    end if                      
    
    case(4) 
    
!#--------------Constant Stress Degredation---------------#

! Fiber Tensile/Compressive damage

    if((e(1).GT.1.0) .AND. (STRESS(1).GE.0.0)) then
     
        fflags(1)=KINC   
        dmg(1)=dmg(1)*(1.0d+0/e(1)) 

        if (failure_id.EQ.5) then ! Additional induced shear damage for Hashin
        
            dmg(4)=dmg(4)*(1.0d+0/e(1))
            
        end if         


    else if((e(1).GT.1.0) .AND. (STRESS(1).LT.0.0)) then  

        fflags(1)=-KINC  
        dmg(1)=dmg(1)*(1.0d+0/e(1))   

    end if

! Matrix Tensile/Compressive damage

    if((e(2).GT.1.0) .AND. (STRESS(2).GE.0.0)) then  

        fflags(2)=KINC    
        dmg(2)=dmg(2)*(1.0d+0/e(2))   
 
        if(failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)=dmg(4)*(1.0d+0/e(2))
            dmg(6)=dmg(4)*(1.0d+0/e(2))

        end if

        if(failure_id.EQ.7) then ! Additional induced shear damage for Puck

            dmg(3)= dmg(3)*(1.0d+0/e(2))
            dmg(4)= dmg(4)*(1.0d+0/e(2))
            dmg(6)= dmg(6)*(1.0d+0/e(2))
            
        end if 
 
 
    else if((e(2).GT.1.0) .AND. (STRESS(2).LT.0.0)) then  

        fflags(2)=-KINC  
        dmg(2)=dmg(2)*(1.0d+0/e(2)) 

        if(failure_id.EQ.5) then ! Additional induced shear damage for Hashin

            dmg(4)=dmg(4)*(1.0d+0/e(2))
            dmg(6)=dmg(4)*(1.0d+0/e(2))

        end if

        if(failure_id.EQ.7) then ! Additional induced shear damage for Puck

            dmg(3)= dmg(3)*(1.0d+0/e(2))
            dmg(4)= dmg(4)*(1.0d+0/e(2))
            dmg(6)= dmg(6)*(1.0d+0/e(2))
            
        end if 

    end if

! Interlaminar Tensile/Compressive damage

    if((e(3).GT.1.0) .AND. (STRESS(3).GE.0.0)) then  

        fflags(3)=KINC    
        dmg(3)=dmg(3)*(1.0d+0/e(3))  

    else if((e(3).GT.1.0) .AND. (STRESS(2).LT.0.0)) then  

        fflags(3)=-KINC  
        dmg(3)=dmg(3)*(1.0d+0/e(3))   

    end if     

! Shear damage

    do I=4,6 

    if(e(I).GT.1.0) then 

        fflags(I)=KINC  
        dmg(I)=dmg(I)*(1.0d+0/e(I)) 

    end if
                      
    end do      
            
    case(5) 
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


    case default 

! Print error message and exit the program for invalid selection

    write(*,*) "Invalid damage ID"
        
    call XIT
    
    end select 
    

    end subroutine damage_calc      

