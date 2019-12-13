C*************************************************************************
! Technische Universität München
! Chair of Computational Mechanics
! Software Lab Project: Development of Failure Criteria for Composites
C*************************************************************************
! 3D ORTHOTROPIC ELATICITY WITH MAX STRESS CRITERIA AND LINEAR DAMAGE
! (CAN NOT BE USED FOR 2D PROBLEMS)
! VERSION: 1.1
C-------------------------------------------------------------------------
C   GROUP: 11
C	WRITTEN BY AMMAR KHALLOUF
C   DATE: 20.6.2019
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
     4 JSTEP(4),UPSTRAN(NTENS)
	 
C  	USER DEFINED VARAIBLES IN UMAT SUBROUTINE

      PARAMETER(ONE=1.0D0, TWO=2.0D0)
	  
C     GET THE STATE VARIABLES FROM PREVIOUS INCREMENT

! STATE VARIABLES INCLUDE DEMAND CAPACITY RATIOS AT EACH SOLUTION STEP

		Xe1=STATEV(1)  ! D/C ratio X-direction (Local Fiber Direction)
		Xe2=STATEV(2)  ! D/C ratio Y-direction 
		Xe3=STATEV(3)  ! D/C ratio Z-direction 
		Xe4=STATEV(4)  ! D/C ratio XY-direction
		Xe5=STATEV(5)  ! D/C ratio XZ-direction 
		Xe6=STATEV(6)  ! D/C ratio YZ-direction
	    DMAX=STATEV(7) ! Damage Threshold for element deletion
		Ddelete=STATEV(8)
	  	  	 
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

! Read material strength components as per Max Stress Criterion

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
	  
C Damage Initiatilization Stresses  

! Read ratios of Ultimate to residual stresses for progressive damage.

 
      RDXt=PROPS(19)
	  RDXc=PROPS(20)
	  RDYt=PROPS(21)
	  RDYt=PROPS(22)
	  RDYc=PROPS(23)
	  RDZt=PROPS(24)
	  RDS12=PROPS(25)
	  RDS13=PROPS(26)
	  RDS23=PROPS(27)
! Degredation parameters for components of constitutive matrix
	  
	  DMAX=0.98 !not exaclty one for numerical stability    
	  D1=0.D0
	  D2=0.D0
	  D3=0.D0
	  D4=0.D0
	  D5=0.D0
	  D6=0.D0
	  Dtemp=0.D0
	  Ddelete=1.0
C=========================================================================
 	  
C     UPDATE THE STRAIN AT THE END OF ITERATION  
      
      DO I=1,NTENS
          UPSTRAN(I)=STRAN(I)+DSTRAN(I)
      END DO
C=========================================================================

C 	COMPUTE STIFFNESS MATRIX WITH DEGREDATION ! Damage parameters set initially to zero
   
      DO I=1,NTENS
       DO J=1,NTENS
       DDSDDE(I,J)=0.0D0       
	   ENDDO
      ENDDO	  

      DDSDDE(1,1)=(E11*Upsilon*(ONE-XNU23*XNU32))*(1-D1)
      DDSDDE(1,2)=(E11*Upsilon*(XNU21+XNU31*XNU23))*(1-D1)*(1-D2)
      DDSDDE(1,3)=(E11*Upsilon*(XNU31+XNU21*XNU32))*(1-D1)*(1-D3)
      DDSDDE(2,1)=(E11*Upsilon*(XNU21+XNU31*XNU23))*(1-D1)*(1-D2)
      DDSDDE(2,2)=(E22*Upsilon*(ONE-XNU13*XNU31))*(1-D2)
      DDSDDE(2,3)=(E22*Upsilon*(XNU32+XNU12*XNU31))*(1-D2)*(1-D3)
      DDSDDE(3,1)=(E11*Upsilon*(XNU31+XNU21*XNU32))*(1-D1)*(1-D3)
      DDSDDE(3,2)=(E22*Upsilon*(XNU32+XNU12*XNU31))*(1-D2)*(1-D3)
      DDSDDE(3,3)=(E33*Upsilon*(ONE-XNU12*XNU21))*(1-D3)
      DDSDDE(4,4)=G12*(1-D4)
      DDSDDE(5,5)=G13*(1-D5)
      DDSDDE(6,6)=G23*(1-D6)
C=========================================================================

C     COMPUTE THE STRESSES

!     NOTE: ABAQUS uses engineering shear strains,
!     i.e. stran(ndi+1) = 2*e_12, etc...

	   DO I=1,NTENS
          STRESS(I)=0.D0	   
		  DO J=1, NTENS
		
		STRESS(I)=STRESS(I)+DDSDDE(I,J)*UPSTRAN(J)
          END DO
      END DO	  
 
C=========================================================================
 
C SEPERATE TENSILE AND COMPRESSIVE NORMAL STRESSES 

        IF (STRESS(1).GE. 0.0) THEN 
            X=Xt
        ELSE 
            X=-1.0*Xc
        END IF 
        IF (STRESS(2).GE. 0.0) THEN 
            Y=Yt 
        ELSE 
            Y=-1.0*Yc 
        END IF
		
        IF (STRESS(3).GE. 0.0) THEN 
            Z=Zt
        ELSE 
            Z=-1.0*Zc
        END IF
		
        IF (STRESS(4).GE. 0.0) THEN 
            S12=S12
        ELSE 
            S12=-1.0*S12
        END IF		
		
        IF (STRESS(5).GE. 0.0) THEN 
            S13=S13
        ELSE 
            S13=-1.0*S13
        END IF

        IF (STRESS(6).GE. 0.0) THEN 
            S23=S23
        ELSE 
            S23=-1.0*S23
        END IF				
		
C=========================================================================
		
C MAXIMUM STRESS CRITERIA EVALUATION

! Ratios of Apllied Stress components to material strength 
! in respective direction.

	 Xe1=STRESS(1)/X
	 Xe2=STRESS(2)/Y
	 Xe3=STRESS(3)/Z
	 Xe4=STRESS(4)/S12
	 Xe5=STRESS(5)/S13
	 Xe6=STRESS(6)/S23
	 

	 
C=========================================================================	
	
C  DAMAGE EVOLUTION

! Linear reduction of material properies
! if STRESS(I)<SD(I) DAMAGE =0
! If < STRESS(I)< X(I),Y(I),Z(I),S(I) DAMAGE BETWEEN 0 AND DMAX

! https://abaqus-docs.mit.edu/2017/English/SIMACAEMATRefMap/simamat-c-damageevolductile.htm#simamat-c-damageevolductile-t-TabularForm-sma-topic4__simamat-c-dmgevol-plastic-disp

       IF (STATEV(1).GE.RDXt) THEN 
          Dtemp= DMAX
		  
      ELSE IF ((STATEV(1).LT.RDXt) .AND. (STATEV(1).GT.1.0)) THEN
	          Dtemp= DMAX*((STATEV(1)-1.0)/(RDXt-1.0))
      ELSE 
		  Dtemp=0.0
      END IF
	  
       IF (Dtemp.GT.D1) THEN
		  D1= Dtemp
      END IF 
	   
      IF (STATEV(2).GE.RDYt) THEN 
          Dtemp= DMAX
      ELSE IF ((STATEV(2).LT.RDYt) .AND. (STATEV(2).GT.1.0)) THEN
	      Dtemp= DMAX*((STATEV(2)-1.0)/(RDYt-1.0))
      ELSE 
		  Dtemp=0.0
      END IF
	  
	   IF (Dtemp.GT.D2) THEN
		   D2= Dtemp
	  END IF 
	   
      IF (STATEV(3).GE.RDZt) THEN 
          Dtemp= DMAX
      ELSE IF ((STATEV(3).LT.RDZt) .AND. (STATEV(3).GT.1.0)) THEN
	      Dtemp= DMAX*((STATEV(3)-1.0)/(RDZt-1.0))
      ELSE 
		  Dtemp=0.0
      END IF
	  
       IF (Dtemp.GT.D3) THEN
	      D3= Dtemp
	   END IF 
	   
       IF (STATEV(4).GE.RDS12) THEN 
          Dtemp= DMAX
      ELSE IF ((STATEV(4).LT.RDS12) .AND. (STATEV(4).GT.1.0)) THEN
	      Dtemp= DMAX*((STATEV(4)-1.0)/(RDS12-1.0))
      ELSE 
          Dtemp=0.0
      END IF
	  
       IF (Dtemp.GT.D4) THEN
	     D4= Dtemp
      END IF 
	  
       IF (STATEV(5).GE.RDS13) THEN 
          Dtemp= DMAX
      ELSE IF ( (STATEV(5).LT.RDS13) .AND. (STATEV(5).GT.1.0) ) THEN
	      Dtemp= DMAX*((STATEV(5)-1.0)/(RDS13-1.0))
      ELSE 
          Dtemp=0.0
      END IF
	  
       IF (Dtemp.GT.D5) THEN
	     D5= Dtemp
      END IF
	  
       IF (STATEV(6).GE.RDS23) THEN 
          Dtemp= DMAX
      ELSE IF ( (STATEV(6).LT.RDS23) .AND. (STATEV(6).GT.1.0) ) THEN
	      Dtemp= DMAX*((STATEV(6)-1.0)/(RDS23-1.0))
      ELSE 
          Dtemp=0.0
      END IF
	  
      IF (Dtemp.GT.D6) THEN
	      D6= Dtemp
      END IF	

C=========================================================================	  
		
C     UPDATE STIFFNESS MATRIX 
	  
      DO I=1,NTENS
       DO J=1,NTENS
       DDSDDE(I,J)=0.0D0       
	   ENDDO
      ENDDO	  

      DDSDDE(1,1)=(E11*Upsilon*(ONE-XNU23*XNU32))*(1-D1)
      DDSDDE(1,2)=(E11*Upsilon*(XNU21+XNU31*XNU23))*(1-D1)*(1-D2)
      DDSDDE(1,3)=(E11*Upsilon*(XNU31+XNU21*XNU32))*(1-D1)*(1-D3)
      DDSDDE(2,1)=(E11*Upsilon*(XNU21+XNU31*XNU23))*(1-D1)*(1-D2)
      DDSDDE(2,2)=(E22*Upsilon*(ONE-XNU13*XNU31))*(1-D2)
      DDSDDE(2,3)=(E22*Upsilon*(XNU32+XNU12*XNU31))*(1-D2)*(1-D3)
      DDSDDE(3,1)=(E11*Upsilon*(XNU31+XNU21*XNU32))*(1-D1)*(1-D3)
      DDSDDE(3,2)=(E22*Upsilon*(XNU32+XNU12*XNU31))*(1-D2)*(1-D3)
      DDSDDE(3,3)=(E33*Upsilon*(ONE-XNU12*XNU21))*(1-D3)
      DDSDDE(4,4)=G12*(1-D4)
      DDSDDE(5,5)=G13*(1-D5)
      DDSDDE(6,6)=G23*(1-D6)
	  
C=========================================================================
   
C     RECOMPUTE STRESSES

	   DO I=1,NTENS
          STRESS(I)=0.D0	   
		  DO J=1,NTENS
		
		STRESS(I)=STRESS(I)+DDSDDE(I,J)*UPSTRAN(J)
          END DO
      END DO	

C=========================================================================	  
            
C     UPDATE STATE VARIABLES

      STATEV(1)=Xe1
      STATEV(2)=Xe2
	  STATEV(3)=Xe3
	  STATEV(4)=Xe4
	  STATEV(5)=Xe5
	  STATEV(6)=Xe6
	  STATEV(7)=DMAX
	  
	        IF (STATEV(1).GE. 2.5) THEN 
            Ddelete=0.0
              END IF 
			IF (STATEV(2).GE. 2.5) THEN 
            Ddelete=0.0 
          	  END IF
		    IF (STATEV(3).GE. 2.5) THEN 
			Ddelete=0.0
              END IF
	        IF (STATEV(4).GE. 2.5) THEN 
            Ddelete=0.0
              END IF 
			IF (STATEV(5).GE. 2.5) THEN 
            Ddelete=0.0 
          	  END IF
		    IF (STATEV(6).GE. 2.5) THEN 
			Ddelete=0.0
              END IF	
			  
	  STATEV(8)=Ddelete
	  
	        
      RETURN
      END
	  
      
