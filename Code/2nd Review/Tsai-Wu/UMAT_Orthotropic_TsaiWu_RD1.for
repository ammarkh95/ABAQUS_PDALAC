C*************************************************************************
! Technische Universität München
! Chair of Computational Mechanics
! Software Lab Project: Development of Failure Criteria for Composites
C*************************************************************************
! 3D ORTHOTROPIC ELATICITY WITH TSAI-WU CRITERIA AND PLY-DISCOUNT DAMAGE
! DAMAGE MODEL: CONSTANT DEGREDATION (RECURSIVE)
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
     4 JSTEP(4),UPSTRAN(NTENS),CD(NTENS,NTENS),Dgrd(3),dmg(6),Xphi(6)
	 
    
	 
C  	USER DEFINED VARAIBLES IN UMAT SUBROUTINE

      PARAMETER(ZERO = 0.0D0,ONE=1.0D0, TWO=2.0D0)
	  
	  
	  DOUBLE PRECISION F1,F2,F3,F11,F22,F33,F44,F55,F66,F12,F13,F23,Phi,Phimax,DMAX
	  INTEGER fflags(6),IMAX(6)
	  
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
	  
C=========================================================================
 	  
C     UPDATE THE STRAIN AT THE END OF INCREMENT  
      
      DO I=1,NTENS
          UPSTRAN(I)=STRAN(I)+DSTRAN(I)
      END DO
C=========================================================================

C     GET THE STATE VARIABLES FROM PREVIOUS INCREMENT

! STATE VARIABLES INCLUDE DAMAGE FACTORS AND FAILURE FLAGS

        dmg(1:6)=STATEV(1:6)
		fflags(1:6)=STATEV(7:12)
  	    Ddelete=STATEV(13) ! Delete Element
		
! INITIALIZE TO ONE AT FIRST INCREMENT 
		
		IF(KINC.EQ.1) THEN
		
		DO I=1,6
		
		dmg(I)=1.0
		fflags(I)=0
		
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

172	  ST11=STRESS(1) 
      ST22=STRESS(2) 
      ST33=STRESS(3) 
      ST12=STRESS(4) 
      ST13=STRESS(5) 
      ST23=STRESS(6) 
C=========================================================================	
C     CALCULATE TSAI-WU STRENGTH PARARMETERS

	  F1=(1.0D0/Xt)-(1.0D0/Xc)
	  F2=(1.0D0/Yt)-(1.0D0/Yc)
      F3=(1.0D0/Zt)-(1.0D0/Zc)
      F11=(1.0D0/(Xt*Xc))
      F22=(1.0D0/(Yt*Yc))
      F33=(1.0D0/(Zt*Zc))
      F44=(1.0D0/(S13*S13))
      F55=(1.0D0/(S23*S23))
      F66=(1.0D0/(S12*S12))
      F12=-0.5*(1.0D0/SQRT(Xt*Xc*Yt*Yc))
      F13=-0.5*(1.0D0/SQRT(Xt*Xc*Zt*Zc))
      F23=-0.5*(1.0D0/SQRT(Yt*Yc*Zt*Zc))
	  


C=========================================================================
 
C   COMPUTE TSAI-WU POLYNOMIALS

      Xphi(1)=(F1*ST11)+(F11*ST11*ST11)+(F12*ST11*ST22)+(F13*ST11*ST33)

      Xphi(2)=(F2*ST22)+(F22*ST22*ST22)+(F12*ST11*ST22)+(F23*ST22*ST33)

      Xphi(3)=(F3*ST33)+(F33*ST33*ST33)+(F13*ST11*ST33)+(F23*ST22*ST33)

      Xphi(4)=(F66*ST12*ST12)

      Xphi(5)=(F44*ST13*ST13)

      Xphi(6)=(F55*ST23*ST23)
        
      Phi=Xphi(1)+Xphi(2)+Xphi(3)+Xphi(4)+Xphi(5)+Xphi(6) 

 
	 
C=========================================================================	

C CHECK FAILURE INITIATION AS PER TSAI-WU CRITERIA


        IF (Phi.LE.1.0D0) THEN

! IF FAILURE NOT DETECTED ALL SDV'S SET TO ZERO

            DO I=1,6
			   dmg(I)=1.0
			   fflags(I)=0
			   IMAX(I)=0
			   
			 END DO
            
			Phimax=0.0
           			   
        ELSE 
		
! WHEN FAILURE IS DETECTED LOOP THROUGH TSAI-WU POLYNOMIAL AND FIND MAXIMUM CONTRIBUTION

            Phimax= maxval(Xphi)
		
			DO I=1,6
		
			 IF(Xphi(I).EQ.Phimax) THEN
			   
! FAILURE FLAGS SET TO INCREMENT NUMBER AND DOMINANT FAILURE MODE DETECTED
			   
			   fflags(I)=KINC
			   IMAX(I)=I
			   
			   ELSE
			   
			   fflags(I)=0
			   IMAX(I)=0	   
			   
		     ENDIF
			   END DO
			  
		ENDIF
		
C=========================================================================
			   
C   CALCULATE DAMAGE FACTORS AS PER DOMINANT FAILURE MODE

! HANDLE CHECK FOR PREVIOUS DAMAGE INCREMENTS

        DO I=1,6 
		    IF(dmg(I).NE.1.0) THEN
			
		    GOTO 318
			
			ELSE
			
			GOTO 281
			
			ENDIF
             END DO

	
! TENSION OR COMPRESSION FAILURE
			   
281		    DO I=1,3
			 
			    IF ((I.EQ.IMAX(I)).AND.(STRESS(I).GE.0.0)) THEN
			 
			     dmg(I)=(0.99-Dgrd(1))
				
			    ELSE IF((I.EQ.IMAX(I)).AND.(STRESS(I).LT.0.0)) THEN
				
				 dmg(I)=(0.99-Dgrd(2))
				
			    ELSE 
			
			    dmg(I)=1.0
			
			     END IF
			      END DO
			  		 

! SHEAR FAILURE			 
			
			    DO I=4,6
			
			     IF (I.EQ.IMAX(I)) THEN
			 
			         dmg(I)=(0.99-Dgrd(3))
			 
			     ELSE
			
			     dmg(I)=1.0
			 
			     END IF
			       END DO
				   
				
			
 ! RECURSIVE DAMAGE FOR NON ZERO DAMAGE VARIABLES		
		     
318		        DO I=1,3 
					 
		    IF ((dmg(I).NE.1.0) .AND. (STRESS(I).GE.0.0)) THEN
			
	                dmg(I)=(1.0-Dgrd(1))**(25)     !dmg(I)*Dgrd(1)
	 
	        ELSE IF((dmg(I).NE.1.0) .AND. (STRESS(I).LT.0.0)) THEN

                 dmg(I)=(1.0-Dgrd(2))**(25) !dmg(I)*Dgrd(2)

	          END IF
	             END DO
					
	        DO I=4,6
			
			IF (dmg(I).NE.1.0) THEN
									
	        dmg(I)=(1.0-Dgrd(3))**(25) !dmg(I)*Dgrd(3)
				  
				  END IF
			          END DO

            					  
   	
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
		STATEV(7:12)=fflags(1:6) 
		STATEV(13)=Ddelete
		
		
		DMAX= minval(dmg)
			
		IF ((DMAX.NE.1.0) .AND. (DMAX.GT.0.0001)) THEN
			
		GOTO 172
			
		ELSE
		
		GOTO 455
		
		END IF
		 


	  	  
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



	  ! write(6,*) ' INCREMENT '
	  ! write(6,*) KINC
	  ! write(6,*) ' ---------------------------- '
	  ! write(6,*) 'damage1'
	  ! write(6,*) dmg(1)
	  ! write(6,*) 'damage2'
	  ! write(6,*) dmg(2)
	  ! write(6,*) 'damage3'
	  ! write(6,*) dmg(3)
	  ! write(6,*) 'damage4'
	  ! write(6,*) dmg(4)
	  ! write(6,*) 'damage5'
	  ! write(6,*) dmg(5)
	  ! write(6,*) 'damage6'
	  ! write(6,*) dmg(6)
      ! write(6,*) ' ---------------------------- '
	  
	  
	  
	  
	  
	        
455     RETURN
      END
	  
      
