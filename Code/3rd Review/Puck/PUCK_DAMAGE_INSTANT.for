C*************************************************************************
! Technische Universität München
! Chair of Computational Mechanics
! Software Lab Project: Development of Failure Criteria for Composites
C*************************************************************************
! PUCK'S CRITERION WITH PLY-DISCOUNT DAMAGE
! PUCK THEORY FROM https://d-nb.info/1010526227/34
! DAMAGE MODEL: RECURSIVE DAMAGE
! AXIS 1 SHOULD BE IN FIBER DIRECTION
! VERSION: 2.0
C-------------------------------------------------------------------------
C   GROUP: 11
C	WRITTEN BY PANAGIOTIS GAVALLAS
C   DATE: 05.11.2019
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
      DIMENSION STRESS(NTENS),STATEV(NSTATV),FARRAY(91),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4),UPSTRAN(NTENS),e(6)
	 
C  	USER DEFINED VARAIBLES IN UMAT SUBROUTINE

      PARAMETER(ONE=1.0D0, TWO=2.0D0, pi = 3.14)
	  

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
	  
C=========================================================================		  
	  
C Failure Stresses (Tension,Compression,Shear)

! Read material strength components

	  Xt=PROPS(10)
	  Xc=PROPS(11)
	  Yt=PROPS(12)
	  Yc=PROPS(13)
	  S12=PROPS(14)
	  S23=PROPS(15)   
	  Dgrd=PROPS(16)   !CHOOSE DAMAGE FACTOR
  
C=========================================================================

C     UPDATE THE STRAIN AT THE END OF ITERATION  

      
      DO I=1,NTENS
          UPSTRAN(I)=STRAN(I)+DSTRAN(I)
      END DO
	
C=========================================================================
C     GET THE STATE VARIABLES FROM PREVIOUS INCREMENT


  	    dmg1=STATEV(1)
		dmg2=STATEV(2)
		dmg3=STATEV(3)
		
C INITIALIZE DEGRADATION TO ONE AT FIRST INCREMENT 
		
		IF(KINC.EQ.1) THEN
             dmg1= 1.0
			 dmg2= 1.0
			 dmg3= 1.0
		END IF

		

C 	BUILD CONSTITUTIVE TENSOR 

      
      DO I=1,NTENS
        DO J=1,NTENS
             DDSDDE(I,J)=0.0D0       
	     END DO
      END DO
	  
    
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
            STRESS(I)=0.D0	   
		    DO J=1, NTENS
		         STRESS(I)=STRESS(I)+DDSDDE(I,J)*UPSTRAN(J)
            END DO
        END DO
	  
       S1=STRESS(1) 
       S2=STRESS(2) 
       S3=STRESS(3) 
       T12=STRESS(4) 
       T23=STRESS(5) 
       T13=STRESS(6) 
 
C CRITERIA EVALUATION 

C=========================================================================
C   FIBER FAILURE
        IF (S1.GE. 0.0) THEN
          fFF = S1/Xt
        ELSE IF (S1.LE. 0.0) THEN
          fFF = abs(S1/Xc)
        END IF



C   INTERFIBER FAILURE

C   FIND ACTION PLANE MOST LIKELY TO FAIL BY CALCULATING MAXIMUM FAILURE INDEX fIFF
         FMAX = 0
         DO THETA =-pi/2, pi/2 ,pi/180   !LOOK THROUGH THETA, FIND PLANE WITH MAX FAILURE INDEX

             C = COS(THETA)
             S = SIN(THETA)
			 !CALCULATE ACTION PLANE STRESSES
             SN = S2 * C**2 + S3 * S**2 + 2 * T23 * C * S 
             TNT = - S2 * S * C + S3 * S * C + T23 * (C**2 - S**2) 
             TN1 = T13 * S + T21 * C

             RAPP = Yc*2 / 5

             Cpsi = (TNT**2)/( TNT**2 - TN1**2)
             Spsi = 1 - Cpsi

             A = (1/(4 * RAPP))* Cpsi + (1/(4 * S12)) * Spsi

             IF ( SN.GE. 0.0) THEN
                 F = SQRT( (( 1/Yt - A) *SN)**2 + (TNT/RAPP)**2 + (TN1/S12)**2) + A * SN 
             ELSE
                 F = SQRT( (TNT/RAPP)**2 + (TN1/S12)**2 + (A * SN)**2) + A* SN
	         END IF
			 !DIFFERENT EQUATIONS FOR FAILURE INDEX IN CASE OF TENSION/COMPR. IN ACTION PLANE
			 
             IF (F.GE. FMAX) THEN        !FIND MAX INDEX
			     FMAX = F
				 SNF = SN
			END IF
         END DO
        fIFF = FMAX
				 
	 	 	 	 	 
C=========================================================================	
c CHECK FAILURE INITIATION
! FIBER FAILURE
        IF ( fFF.GE. 1.0) THEN
             dmg1 = (0.999999 -Dgrd)

		ELSE
		     dmg1 = 1     !DECREASE STIFFNESS IN FIBER DIRECTION
		END IF
! INTERFIBER FAILURE
        IF (fIFF.GE. 1.0) THEN
		    dmg3 = (0.999999 -Dgrd)      !DECREASE SHEAR STIFNESS
			IF  (SNF.GE. 0.0) THEN
			    dmg2 = (0.999999 -Dgrd)   !DECREASE LATERAL STIFFNESS IF NORMAL STRESS ON ACTION PLANE IS POSITIVE
			ELSE
			    dmg2 = 1
			END IF
		ELSE 
		    dmg3 = 1
		END IF




C=========================================================================	  
		
C 	BUILD DEGRADED CONSTITUTIVE TENSOR IN CASE OF FAILURE DETECTION  
   
      DO I=1,NTENS
       DO J=1,NTENS
       DDSDDE(I,J)=0.0D0       
	   ENDDO
      ENDDO
	  
	  DDSDDE(1,1)=(E11*Upsilon*(ONE-XNU23*XNU32))* abs(dmg1)
      DDSDDE(1,2)=(E11*Upsilon*(XNU21+XNU31*XNU23))
      DDSDDE(1,3)=(E11*Upsilon*(XNU31+XNU21*XNU32))
      DDSDDE(2,1)=(E11*Upsilon*(XNU21+XNU31*XNU23))
      DDSDDE(2,2)=(E22*Upsilon*(ONE-XNU13*XNU31))*abs(dmg2)
      DDSDDE(2,3)=(E22*Upsilon*(XNU32+XNU12*XNU31))
      DDSDDE(3,1)=(E11*Upsilon*(XNU31+XNU21*XNU32))
      DDSDDE(3,2)=(E22*Upsilon*(XNU32+XNU12*XNU31))
      DDSDDE(3,3)=(E33*Upsilon*(ONE-XNU12*XNU21))*abs(dmg2)
      DDSDDE(4,4)=G12*abs(dmg3)
      DDSDDE(5,5)=G13*abs(dmg3)
      DDSDDE(6,6)=G23
	  
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

		STATEV(1) = dmg1
		STATEV(2) = dmg2
		STATEV(3) = dmg3
	    STATEV(4) = fFF ! FIBER FAILURE INDEX
		STATEV(5) = fIFF ! INTER-FIBER FAILURE INDEX
		STATEV(6) = DDSDDE(1,1) !STIFFNESS IN FIBER DIRECTION
		STATEV(7) = DDSDDE(2,2) !LATERAL STIFFNESS
		STATEV(8) = DDSDDE(4,4) !SHEAR STIFFFNESS
! SDVs 6-8 CAN BE USED TO CHECK WHETHER STIFFNESS IS DEGRADED

C=========================================================================	  	

      RETURN
      END
	  
      
