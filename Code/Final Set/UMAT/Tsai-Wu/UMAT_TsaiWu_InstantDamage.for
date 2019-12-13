C*************************************************************************
! Technische Universität München
! Chair of Computational Mechanics
! Software Lab Project: Development of Failure Criteria for Composites
C*************************************************************************
! 3D ORTHOTROPIC ELATICITY WITH TSAIWU CRITERION AND PLY-DISCOUNT DAMAGE
! DAMAGE MODEL: INSTANT DEGREDATION
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
    
!            ****Ply Discount Damage Parameters****
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

! Important Remarks: 

!*Degradation factors values are in the range of (0,1) 
!*Degredation factor == 0 implies no damage
!*Degredation factor == 1 implies full damage (i.e stiffness approx zero)

    beta_ft=PROPS(19) 
    beta_fc=PROPS(20)     
    beta_mt=PROPS(21)   
    beta_mc=PROPS(22)   
    beta_s=PROPS(23)

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
       
!###########     Failure Evaluation: Tsai-Wu      ###########
!#----------------------------------------------------------#       	  
   
    call failure_calc(STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23)    
          
!###########      Damage Evaluation: Instant      ###########
!#----------------------------------------------------------#
    
    if ((maxval(e).GT.1.0) .OR. (minval(dmg).LT.1.0)) then

        call damage_calc (STRESS,UPSTRAN,dmg,e,fflags,KINC,beta_ft,beta_fc,beta_mt,beta_mc,beta_s) 
    
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
    
    end subroutine ortho3d

!###########  Utility Subroutine 2: failure_calc  ###########
!#----------------------------------------------------------#

    subroutine failure_calc(STRESS,UPSTRAN,e,Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23)
    
    implicit none
    real*8, dimension  (6) :: e,STRESS,UPSTRAN,Xphi
    real*8 :: Xt,Xc,Yt,Yc,Zt,Zc,S12,S13,S23,Phi
    real*8 :: F1,F2,F3,F11,F22,F33,F44,F55,F66,F12,F13,F23
    integer :: IMAX
    
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
    
    end subroutine failure_calc
    
!###########   Utility Subroutine 2: damage_calc  ###########
!#----------------------------------------------------------#
 
    subroutine damage_calc (STRESS,UPSTRAN,dmg,e,fflags,KINC,beta_ft,beta_fc,beta_mt,beta_mc,beta_s)    
    
    implicit none
    real*8, dimension  (6) :: e,STRESS,UPSTRAN,dmg
    real*8 :: beta_ft,beta_fc,beta_mc,beta_mt,beta_s    
    integer:: I,KINC,fflags(6)
    
!#----------------Instantaneous Degredation-----------------#

! Fiber Tensile/Compressive damage

    if((e(1).GT.1.0) .AND. (STRESS(1).GE.0.0)) then
     
        fflags(1)=KINC   
        dmg(1)=(0.99999999-beta_ft) 
                   
    else if((e(1).GT.1.0) .AND. (STRESS(1).LT.0.0)) then  

        fflags(1)=-KINC  
        dmg(1)=(0.99999999-beta_fc)  

    end if

! Matrix Tensile/Compressive damage

    if((e(2).GT.1.0) .AND. (STRESS(2).GE.0.0)) then  

        fflags(2)=KINC    
        dmg(2)= (0.99999999-beta_mt)  

        
    else if((e(2).GT.1.0) .AND. (STRESS(2).LT.0.0)) then  

        fflags(2)=-KINC  
        dmg(2)= (0.99999999-beta_mc)  

    end if

! Interlaminar Tensile/Compressive damage

    if((e(3).GT.1.0) .AND. (STRESS(3).GE.0.0)) then  

        fflags(3)=KINC    
        dmg(3)= (0.99999999-beta_mt) 

    else if((e(3).GT.1.0) .AND. (STRESS(2).LT.0.0)) then  

        fflags(3)=-KINC  
        dmg(3)= (0.99999999-beta_mc)  

    end if     

! Shear damage

    do I=4,6 

    if(e(I).GT.1.0) then 

        fflags(I)=KINC  
        dmg(I)= (0.99999999-beta_s) 
        
    end if                      
    end do

    end subroutine damage_calc   