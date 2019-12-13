C*************************************************************************
! Technische Universität München
! Chair of Computational Mechanics
! Software Lab Project: Development of Failure Criteria for Composites
C*************************************************************************
! 3D ORTHOTROPIC ELATICITY WITH HASHIN CRITERIA AND PLY-DISCOUNT DAMAGE
! DAMAGE MODEL: LINEAR SLOPE DEGREDATION
! (CAN NOT BE USED FOR 2D PROBLEMS)
! VERSION: 1.2
C-------------------------------------------------------------------------
C   GROUP: 11
C	WRITTEN BY AMMAR KHALLOUF
C   DATE: 19.7.2019
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
     4 JSTEP(4),UPSTRAN(NTENS),CD(NTENS,NTENS),Dgrd(3),dmg(6),Xphi(6),e(3)
	 
    
	 
C  	USER DEFINED VARAIBLES IN UMAT SUBROUTINE

      PARAMETER(ZERO = 0.0D0,ONE=1.0D0, TWO=2.0D0)
	  
	  
	  DOUBLE PRECISION emax
	  INTEGER fflags(3)
	  
C=========================================================================		  
  	  	 
C ELASTIC PROPERTIES

! Read Nine Indendent Constants of orthoptropic materials
! Material model can be viewed here: 
! https://www.sharcnet.ca/Software/Abaqus610/Documentation/docs/v6.10/books/usb/default.htm?startat=pt05ch19s02abm02.html#usb-mat-clinearelastic-engconst

      E11=PROPS(1)
      E22=PROPS(2)
      E33=PROPS(3)
      XNU12=PROPS(4)
      XNU13=PROPS(5)
      XNU23=PROPS(6)   
      G12=PROPS(7)
      G13=PROPS(8)
      G23=PROPS(9)
	  
	  XNU21=(E22*XNU12)/E11
      XNU31=(E33*XNU13)/E11
      XNU32=(E33*XNU23)/E22
	  Upsilon=(ONE)/(ONE-XNU12*XNU21-XNU23*XNU32-XNU13*XNU31-TWO*XNU21*XNU32*XNU13)
	  
C Failure Stresses (Tension,Compression,Shear)

! Read material strength components

	  Xt=PROPS(10)
	  Xc=PROPS(11)
	  Yt=PROPS(12)
	  Yc=PROPS(13)
	  Zt=PROPS(14)
	  Zc=PROPS(15)
	  S12=PROPS(16)
	  S13=PROPS(17)
	  S23=PROPS(18)

C=========================================================================	
	  
C Degredation Parameters 

! Read Damage Factors Array
! Tension, Compression, Shear

      Dgrd(1:3)=PROPS(19:21)
	  emax=props(22)
	  
C=========================================================================
 	  
C     UPDATE THE STRAIN AT THE END OF INCREMENT  
      
      DO I=1,NTENS
          UPSTRAN(I)=STRAN(I)+DSTRAN(I)
      END DO
C=========================================================================

C     GET THE STATE VARIABLES FROM PREVIOUS INCREMENT

! STATE VARIABLES INCLUDE DAMAGE FACTORS AND FAILURE FLAGS

        dmg(1:6)=STATEV(1:6)
		fflags(1:3)=STATEV(7:9)
  	    Ddelete=STATEV(10) ! Delete Element
		
! INITIALIZE TO ONE AT FIRST INCREMENT 
		
		IF(KINC.EQ.1) THEN
		
		DO I=1,6
		
		dmg(I)=1.0
		
		
		END DO
		  END IF
		
C=========================================================================

C 	BUILD CONSTITUTIVE TENSOR 

      
      DO I=1,NTENS
       DO J=1,NTENS
       DDSDDE(I,J)=0.0D0       
	   ENDDO
      ENDDO
	  
    
      DDSDDE(1,1)=(E11*Upsilon*(ONE-XNU23*XNU32))
      DDSDDE(1,2)=(E11*Upsilon*(XNU21+XNU31*XNU23))
      DDSDDE(1,3)=(E11*Upsilon*(XNU31+XNU21*XNU32))
      DDSDDE(2,1)=(E11*Upsilon*(XNU21+XNU31*XNU23))
      DDSDDE(2,2)=(E22*Upsilon*(ONE-XNU13*XNU31))
      DDSDDE(2,3)=(E22*Upsilon*(XNU32+XNU12*XNU31))
      DDSDDE(3,1)=(E11*Upsilon*(XNU31+XNU21*XNU32))
      DDSDDE(3,2)=(E22*Upsilon*(XNU32+XNU12*XNU31))
      DDSDDE(3,3)=(E33*Upsilon*(ONE-XNU12*XNU21))
      DDSDDE(4,4)=G12
      DDSDDE(5,5)=G13
      DDSDDE(6,6)=G23
	  
	

C=========================================================================	  
C     COMPUTE THE STRESSES

!     NOTE: ABAQUS uses engineering shear strains,
!     i.e. stran(ndi+1) = 2*e_12, etc...
   	
	   DO I=1,NTENS
          STRESS(I)=0.0D0	   
		  DO J=1, NTENS
	    
	 	   STRESS(I)=STRESS(I)+DDSDDE(I,J)*UPSTRAN(J)
   
          END DO
      END DO	  

 
C=========================================================================
! TEMPORARY VARAIBLES TO SAVE STRESSES

	  ST11=STRESS(1) 
      ST22=STRESS(2) 
      ST33=STRESS(3) 
      ST12=STRESS(4) 
      ST13=STRESS(5) 
      ST23=STRESS(6) 
C=========================================================================	
C     CALCULATE HASHIN FAILURE PARAMETERS

      IF(ST11.GE.0.0) THEN
	  
	  e(1)= SQRT(((ST11/Xt)**2) + (((ST12**2)+(ST13**2))/(S12)))  
	  
	  ELSE
	  
	  e(1)= SQRT(((ST11/Xc)**2))
	  
	  END IF
	  
	  IF((ST22+ST33).GE.0.0) THEN
	  
	  e(2)= SQRT((((ST22+ST33)**2)/(Yt)**2) + (((ST23**2)-(ST22*ST33))/(S23**2)) + (((ST12**2)+(ST13**2))/(S12**2)))
	  
	  ELSEIF((ST22+ST33).LT.0.0) THEN
	  
	  e(2)= SQRT((((Yc/(2*S23))**2)-1)*((ST22+ST33)/(Yc))+((ST22+ST33)**2/(4*S23**2))+((ST23**2-ST22*ST33)/(S23**2))+((ST12**2+ST13**2)/(S12**2)))

      END IF
	  
      IF(ST33.GE.0) THEN
	  
	  e(3)= SQRT(((ST33)/(Zt))**2)
	  
	  ELSE
	  
	  e(3)= SQRT(((ST33)/(Zc))**2)
	  
	  END IF
	  
	  
C=========================================================================
 
C CHECK FAILURE INITIATION AS PER HASHIN

	DO I=1,3
	
     	IF (e(i).LE.1.0D0) THEN

!!   IF FAILURE NOT DETECTED ALL SDV'S SET TO ONE

   			   dmg(I)=1.0
			   fflags(I)=0
			   dmg(4)=1.0
			   dmg(5)=1.0
			   dmg(6)=1.0 
			   
	    END IF  			   
		  END DO
			 
 ! CHECK FAILURE MODES INDEPENDENTLEY AND ASSIGN DAMAGE VARIABLES
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			           			   
        IF((e(1).GT.1.0) .AND. (ST11.GE.0.0)) THEN
		
	!  Tensile fiber failure
	
		fflags(1)=KINC
		dmg(1)= (0.99-Dgrd(1))*(emax/e(1))
		dmg(4)= (0.99-Dgrd(3))*(emax/e(1))	   
		
		ELSE IF((e(1).GT.1.0) .AND. (ST11.LT.0.0)) THEN
		
  ! Compressive fiber failure
		
		fflags(1)=-1*KINC
		dmg(1)= (0.99-Dgrd(2))*(emax/e(1))
		
		END IF
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        IF((e(2).GT.1.0) .AND. ((ST22+ST33).GE.0.0)) THEN
	
  !  Tensile Matrix failure	
  
		fflags(2)=KINC
		dmg(2)= (0.99-Dgrd(1))*(emax/e(2))
		dmg(4)= (0.99-Dgrd(3))*(emax/e(2))
		
		ELSE IF((e(2).GT.1.0) .AND. ((ST22+ST33).LT.0.0)) THEN
		
  ! Compressive Matrix failure
		
		fflags(2)=-KINC
		dmg(2)= (0.99-Dgrd(2))*(emax/e(2))
		dmg(4)= (0.99-Dgrd(3))*(emax/e(2))
		
		END IF		
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
	
	    IF((e(3).GT.1.0) .AND. (ST33.GE.0.0)) THEN
	
  !  Interlaminar normal tensile failure
  
		fflags(3)=KINC
		dmg(3)= (0.99-Dgrd(1))*(emax/e(3))
		
		
		ELSE IF((e(3).GT.1.0) .AND. (ST33.LT.0.0)) THEN
		
  ! Interlaminar normal compression failure
		
		fflags(3)=-KINC
		dmg(3)= (0.99-Dgrd(2))*(emax/e(3))
		
		END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
  ! NO DEGREDATION FACTORS PROVIDED BY HASHIN FOR SHEAR INDEPENDENTLEY
		
		dmg(5)=1.0
		dmg(6)=1.0				
		
C=========================================================================

C 	BUILD DEGRADED CONSTITUTIVE TENSOR IN CASE OF FAILURE DETECTION  
   
      DO I=1,NTENS
       DO J=1,NTENS
       DDSDDE(I,J)=0.0D0       
	   ENDDO
      ENDDO
	  
	  DDSDDE(1,1)=(E11*Upsilon*(ONE-XNU23*XNU32))*abs(dmg(1))
      DDSDDE(1,2)=(E11*Upsilon*(XNU21+XNU31*XNU23))*abs(dmg(1)*dmg(2))
      DDSDDE(1,3)=(E11*Upsilon*(XNU31+XNU21*XNU32))*abs(dmg(1)*dmg(3))
      DDSDDE(2,1)=(E11*Upsilon*(XNU21+XNU31*XNU23))*abs(dmg(1)*dmg(2))
      DDSDDE(2,2)=(E22*Upsilon*(ONE-XNU13*XNU31))*abs(dmg(2))
      DDSDDE(2,3)=(E22*Upsilon*(XNU32+XNU12*XNU31))*abs(dmg(3)*dmg(2))
      DDSDDE(3,1)=(E11*Upsilon*(XNU31+XNU21*XNU32))*abs(dmg(3)*dmg(1))
      DDSDDE(3,2)=(E22*Upsilon*(XNU32+XNU12*XNU31))*abs(dmg(3)*dmg(2))
      DDSDDE(3,3)=(E33*Upsilon*(ONE-XNU12*XNU21))*abs(dmg(3))
      DDSDDE(4,4)=G12*abs(dmg(4))
      DDSDDE(5,5)=G13*abs(dmg(5))
      DDSDDE(6,6)=G23*abs(dmg(6))
	  
C=========================================================================	  

C     RECOMPUTE THE STRESSES

!     NOTE: ABAQUS uses engineering shear strains,
!     i.e. stran(ndi+1) = 2*e_12, etc...

	   DO I=1,NTENS
          STRESS(I)=0.0D0	   
		  DO J=1, NTENS
		
		STRESS(I)=STRESS(I)+DDSDDE(I,J)*UPSTRAN(J)
          END DO
      END DO	  
	  
C=========================================================================	  	  
C     UPDATE STATE VARIABLES

		STATEV(1:6)=dmg(1:6)  
		STATEV(7:9)=fflags(1:3) 
		STATEV(10)=Ddelete


	  	  
	        ! IF (STATEV(1).GE. 3.5) THEN 
            ! Ddelete=0.0
              ! END IF 
			! IF (STATEV(2).GE. 3.5) THEN 
            ! Ddelete=0.0 
          	  ! END IF
		    ! IF (STATEV(3).GE. 3.5) THEN 
			! Ddelete=0.0
              ! END IF
	        ! IF (STATEV(4).GE. 3.5) THEN 
            ! Ddelete=0.0
              ! END IF 
			! IF (STATEV(5).GE. 3.5) THEN 
            ! Ddelete=0.0 
          	  ! END IF
		    ! IF (STATEV(6).GE. 3.5) THEN 
			! Ddelete=0.0
              ! END IF
		    ! IF (STATEV(7).GE. 3.5) THEN 
			! Ddelete=0.0
              ! END IF
		    ! IF (STATEV(8).GE. 3.5) THEN 
			! Ddelete=0.0
              ! END IF		
		    ! IF (STATEV(9).GE. 3.5) THEN 
			! Ddelete=0.0
              ! END IF				  	  	  	  


	  
	  
	        
      RETURN
      END
	  
      
