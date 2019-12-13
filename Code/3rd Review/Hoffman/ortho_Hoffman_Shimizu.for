!	Written by Yasuyuki Shimizu
!	Hoffman criteria
!	3D model with using immediate failure criteria

      SUBROUTINE UMAT(STRESS1, STATEV, DDSDDE1, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN1, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, 
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)
	 
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
	  
C
      DIMENSION STRESS1(NTENS), STATEV(NSTATV), DDSDDE1(NTENS, NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN1(NTENS),
     2 PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DROT(3, 3),
     3 DFGRD0(3, 3), DFGRD1(3, 3)
      
      DIMENSION EELAS(6), EPLAS(6), FLOW(6)
	  !elastic strain, plastic strain, flow strain
C
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0,
     1 ENUMAX=.4999D0, NEWTON=10, TOLER=1.0D-6)

	  *SECTION CONTROLS,ELEMENT DELETION=YES"
	  !*SECTION CONTROLS,ELEMENT DELETION=YES,MAX DEGRADATION=0.75,NAME=CONT"

C	<Elastic properties>
C	------------------------------------------------------------------
	!Basic Elastic properties(poisson's ratio, young's moduli, shear moduli)    
	  pr12=PROPS(1)
      pr13=PROPS(2)
      pr23=PROPS(3)
      e11=PROPS(4)
      e22=PROPS(5)
      e33=PROPS(6)   
      g12=PROPS(7)
      g13=PROPS(8)
      g23=PROPS(9)
	  
	!For poisson's ratio, This should be fulfilled(reciprocity relations):
      pr21=e22*pr12/e11
      pr31=e33*pr13/e11
      pr32=e33*pr23/e22

 	!Denominator of stiffness coefficient
      tri=(ONE-pr12*pr21-pr23*pr32-pr31*pr13-TWO*pr12*pr23*pr31)
C	------------------------------------------------------------------

C	<Stiffness matrix for orthotropic model>
	!This should be degraded afterwards
C	------------------------------------------------------------------
	DDSDDE1(1,1)=(1-pr23*pr32)/(tri)
	DDSDDE1(1,2)=(pr21+pr31*pr23)/(tri)
	DDSDDE1(1,3)=(pr31+pr21*pr32)/(tri)
	DDSDDE1(1,4)=0
	DDSDDE1(1,5)=0
	DDSDDE1(1,6)=0
    
	DDSDDE1(2,1)=(pr12+pr13*pr32)/(tri)
	DDSDDE1(2,2)=(1-pr31*pr13)/(tri)
	DDSDDE1(2,3)=(pr32+pr31*pr12)/(tri)
	DDSDDE1(2,4)=0
	DDSDDE1(2,5)=0
	DDSDDE1(2,6)=0
				  
	DDSDDE1(3,1)=(pr13+pr12*pr23)/(tri)
	DDSDDE1(3,2)=(pr23+pr13*pr21)/(tri)
	DDSDDE1(3,3)=(1-pr12*pr21)/(tri)
	DDSDDE1(3,4)=0
	DDSDDE1(3,5)=0
	DDSDDE1(3,6)=0
				  
	DDSDDE1(4,1)=0
	DDSDDE1(4,2)=0
	DDSDDE1(4,3)=0
	DDSDDE1(4,4)=ONE*g11
	DDSDDE1(4,5)=0
	DDSDDE1(4,6)=0
		  
	DDSDDE1(5,1)=0
	DDSDDE1(5,2)=0
	DDSDDE1(5,3)=0
	DDSDDE1(5,4)=0
	DDSDDE1(5,5)=ONE*g22
	DDSDDE1(5,6)=0
				  
	DDSDDE1(6,1)=0
	DDSDDE1(6,2)=0
	DDSDDE1(6,3)=0
	DDSDDE1(6,4)=0
	DDSDDE1(6,5)=0
	DDSDDE1(6,6)=ONE*g33
C	------------------------------------------------------------------

	<Set user defined properties>
C	------------------------------------------------------------------
	! Transverse normal strength
	  XT = PROPS(10)
	  XC = PROPS(11)
	  YT = PROPS(12)
	  YC = PROPS(13)
	  ZT = PROPS(14)
	  ZC = PROPS(15)
	! Shear strength  
	  S12 = PROPS(16)
	  S13 = PROPS(17)
	  S23 = PROPS(18)
	   
	! Degradation factors beta
	! ply-discounting approach (for composite material)
	  BetaT11 = PROPS(19)
	  BetaC11 = PROPS(20)
	  BetaT22 = PROPS(21)
	  BetaC22 = PROPS(22)
	  BetaT33 = PROPS(23)
	  BetaC33 = PROPS(24)
	  BetaS12 = PROPS(25)
	  BetaS13 = PROPS(26)
	  BetaS23 = PROPS(27)
	  
	! Recursive approach
	  DamageCoefficientT11 = PROPS(28)
	  DamageCoefficientC11 = PROPS(29)	  
	  DamageCoefficientT22 = PROPS(30)
	  DamageCoefficientC22 = PROPS(31)
	  DamageCoefficientT33 = PROPS(32)
	  DamageCoefficientC33 = PROPS(33)
	  DamageCoefficientS12 = PROPS(34)
	  DamageCoefficientS13 = PROPS(35)
	  DamageCoefficientS23 = PROPS(36)
	  
	! for material deletion
	  UltStress11 = PROPS(37)
	  UltStress22 = PROPS(38)
	  UltStress33 = PROPS(39)
	  UltStress12 = PROPS(40)
	  UltStress13 = PROPS(41)
	  UltStress23 = PROPS(42)
C	------------------------------------------------------------------  

C	<Compute strain>
C	------------------------------------------------------------------  
      DO K1=1, NTENS
          STRAN1(K1)=STRAN1(K1) + DSTRAN1(K1)
      END DO
C	------------------------------------------------------------------

C	<Hoffman criteria>
C   ! https://webthesis.biblio.polito.it/6883/1/tesi.pdf
	! https://web.fe.up.pt/~stpinho/teaching/feup/y0506/fcriteria.pdf
C	------------------------------------------------------------------
    ! failure indices
	! F12,F13,F23(?)
	  C1 = (1/(ZT*ZC)+1/(YT*YC)-1/(XT*XC))/2
	  C2 = (1/(ZT*ZC)+1/(XT*XC)-1/(YT*YC))/2
	  C3 = (1/(XT*XC)+1/(YT*YC)-1/(ZT*ZC))/2
	! F11,F22,F33(STRESS1(1),STRESS1(2),STRESS1(3))
	  C4 = 1/XT - 1/XC
	  C5 = 1/YT - 1/YC
	  C6 = 1/ZT - 1/ZC
	! S12,S13,S23(STRESS1(4),STRESS1(5),STRESS1(6))
	  C7 = 1/(S12^2)
	  C8 = 1/(S13^2)
	  C9 = 1/(S23^2)
	  
	  phi1 = C2(STRESS1(3)-STRESS1(1)) + C3(STRESS1(1)-STRESS1(2)) + C4*STRESS1(1) 
	  phi2 = C1(STRESS1(2)-STRESS1(3)) + C3(STRESS1(1)-STRESS1(2)) + C5*STRESS1(2)
	  phi3 = C1(STRESS1(2)-STRESS1(3)) + C2(STRESS1(3)-STRESS1(1)) + C6*STRESS1(3)
	  phi4 = C7*STRESS1(4)^2
	  phi5 = C8*STRESS1(5)^2
	  phi6 = C9*STRESS1(6)^2
	  
	  F = C1(STRESS1(2)-STRESS1(3)) + C2(STRESS1(3)-STRESS1(1)) + C3(STRESS1(1)-STRESS1(2))  + C4*STRESS1(1) +  C5*STRESS1(2) + C6*STRESS1(3) + C7*STRESS1(4)^2 + C8*STRESS1(5)^2 + C9*STRESS1(6)^2
	  
	! if the failure index F exceeds unity, we must get a failure mode
	
	  if(F.ge.1) then
		fm = max(phi1,phi2,phi3,phi4,phi5,phi6)
			  
	
	! if STRESS1(1) is main factor of failure 
	  if(fm == phi1) then
		  if(STRESS1(1).ge.ZERO) then
			Beta11 = BetaT11 !for ply-approach 
			DamageCoefficient11 = DamageCoefficientT11 !for Recursive approach
		  else if(STRESS1(1).lt.ZERO) then
			Beta11 = BetaC11
			DamageCoefficient11 = DamageCoefficientC11
		  !else
		  end if
	  end if
	
	 
	! if STRESS1(2) is main factor of failure 
	  if(fm == phi2) then
		  if(STRESS1(2).ge.ZERO) then
			Beta22 = BetaT22 !for ply-approach 
			DamageCoefficient22 = DamageCoefficientT22!for Recursive approach
		  else if(STRESS1(2).lt.ZERO) then
			Beta22 = BetaC22
			DamageCoefficient22 = DamageCoefficientC22
		  !else
		  end if
	  end if
	  
	! if STRESS1(3) is main factor of failure 
	  if(fm == phi3) then
		  if(STRESS1(3).ge.ZERO) then
			Beta33 = BetaT33 !for ply-approach 
			DamageCoefficient33 = DamageCoefficientT33 !for Recursive approach
		  else if(STRESS1(3).lt.ZERO) then
			Beta33 = BetaC33
			DamageCoefficient33 = DamageCoefficientC33
		  !else
		  end if
	  end if
	  
	! if STRESS1(4) is main factor of failure 
	  if(fm == phi4) then
		Beta12 = BetaS12 !for ply-approach 
		DamageCoefficient12 = DamageCoefficientS12 !for Recursive approach
	  end if	

	! if STRESS1(5) is main factor of failure 
	  if(fm == phi5) then
		Beta13 = BetaS13 !for ply-approach 
		DamageCoefficient13 = DamageCoefficientS13 !for Recursive approach
	  end if		  
	  
	! if STRESS1(6) is main factor of failure 
	  if(fm == phi6) then
		Beta23 = BetaS23 !for ply-approach 
		DamageCoefficient23 = DamageCoefficientS23 !for Recursive approach
	  end if		  	  

C	------------------------------------------------------------------


  
	<Degradation of stiffness matrix>
C	------------------------------------------------------------------
	flagForFailureApproach = PROPS(43)
	! choose the failure model
	! 1 - immediate brittle fracture
	! 2 - ply-discounting approach - two linear approximation (as normal material)  
	! 3 - recursive degradation (for composite material) 
	! 4 - multiple linear approximation  
	
	! analytical strain >= ultimate strain : Degradation factor is applied
	! ultimate strain >= analytical strain : Stiffness is kept original
	
	! 1 - for material which occurs immediate brittle fracture after reaching failure stress 
	  ElemDeletion = STATEV(7)
	  
	  if(flagForFailureApproach == 1) then
	    if(ElemDeletion == 0)
		  if(et1.GE.1.0 .OR. ec1.GE.1.0 .OR. et2.GE.1.0 .OR. ec2.GE.1.0 .OR. et3.GE.1.0 .OR. ec3.GE.1.0 .OR.) then
		    DO I = 1, NTENS
			  STRESS1(I) = 0
			END DO
			ElemDeletion = 1
		  endif
		  
	! 2 - for material which has another (lower) constant tangent stiffness 
		if(F.ge.1) then
			  
	
	    ! if STRESS1(1) is main factor of failure 
	      if(fm == phi1) then  
			DDSDDE1(1,1) = Beta11 * DDSDDE1(1,1)
			DDSDDE1(1,2) = Beta11 * DDSDDE1(1,2)
			DDSDDE1(1,3) = Beta11 * DDSDDE1(1,3)
		  endif

		  if(fm == phi2) then
			DDSDDE1(2,2) = Beta22 * DDSDDE1(2,2)
			DDSDDE1(2,1) = Beta22 * DDSDDE1(2,1)
			DDSDDE1(2,3) = Beta22 * DDSDDE1(2,3)
		  endif

		  if(fm == phi3) then
			DDSDDE1(3,3) = Beta33 * DDSDDE1(3,3)
			DDSDDE1(3,1) = Beta33 * DDSDDE1(3,1)
			DDSDDE1(3,2) = Beta33 * DDSDDE1(3,2)
		  endif
		  
		  if(fm == phi4) then
			DDSDDE1(4,4) = Beta12 * DDSDDE1(4,4)
		  !else if((ratioStran12.LT.1.0) .and.(ratioStran12.GT.))
		  !else
		  endif
		  
		  if(fm == phi5) then
			DDSDDE1(5,5) = Beta13 * DDSDDE1(5,5)
		  !else if((ratioStran13.LT.1.0) .and.(ratioStran13.GT.))
		  !else
		  endif
		  
		  if(fm == phi6) then
			DDSDDE1(6,6) = Beta23 * DDSDDE1(6,6)
		  !else if((ratioStran23.LT.1.0) .and.(ratioStran23.GT.))
		  !else
		  endif
		endif
		
		
	! 3 - Recursive Degradation
		if(flagForFailureApproach == 3) then
		  
		  Gamma11 = Gamma11*DamageCoefficient11
		  Gamma22 = Gamma22*DamageCoefficient22
		  Gamma33 = Gamma33*DamageCoefficient33
		  Gamma12 = Gamma12*DamageCoefficient12
		  Gamma13 = Gamma13*DamageCoefficient13
		  Gamma23 = Gamma23*DamageCoefficient23
		  
		  if(fm == phi1) then
			DDSDDE1(1,1) = Gamma11 * DDSDDE1(1,1)
			DDSDDE1(1,2) = Gamma11 * DDSDDE1(1,2)
			DDSDDE1(1,3) = Gamma11 * DDSDDE1(1,3)
		  endif

		  if(fm == phi2) then
			DDSDDE1(2,2) = Gamma22 * DDSDDE1(2,2)
			DDSDDE1(2,1) = Gamma22 * DDSDDE1(2,1)
			DDSDDE1(2,3) = Gamma22 * DDSDDE1(2,3)
		  endif

		  if(fm == phi3) then
			DDSDDE1(3,3) = Gamma33 * DDSDDE1(3,3)
			DDSDDE1(3,1) = Gamma33 * DDSDDE1(3,1)
			DDSDDE1(3,2) = Gamma33 * DDSDDE1(3,2)
		  endif
		  
		  if(fm == phi4) then
			DDSDDE1(4,4) = Gamma12 * DDSDDE1(4,4)
		  !else if((ratioStran12.LT.1.0) .and.(ratioStran12.GT.))
		  !else
		  endif
		  
		  if(fm == phi5) then
			DDSDDE1(5,5) = Gamma13 * DDSDDE1(5,5)
		  !else if((ratioStran13.LT.1.0) .and.(ratioStran13.GT.))
		  !else
		  endif
		  
		  if(fm == phi6) then
			DDSDDE1(6,6) = Gamma23 * DDSDDE1(6,6)
		  !else if((ratioStran23.LT.1.0) .and.(ratioStran23.GT.))
		  !else
		  endif
		  
		endif
	
C	------------------------------------------------------------------
	  
C	<Compute stress>
C	------------------------------------------------------------------  
      DO K1=1, NTENS
          DO K2=1, NTENS
              STRESS1(K2)=STRESS1(K2)+DDSDDE1(K2, K1)*DSTRAN1(K1)
          END DO
      END DO
C	------------------------------------------------------------------

	<Material Deletion>
C	------------------------------------------------------------------
	! With using a STATUS
	! If max_stress, elemdeletion = 0 - Element Deleted
	! Can be used in Abaqus/Standard analysis
	  IF (STRESS1(1).GE.UltStress11.OR.STRESS1(2).GE.UltStress22.OR.STRESS1(3).GE.UltStress33.OR.STRESS1(4).GE.UltStress12.OR.STRESS1(5).GE.UltStress13.OR.STRESS1(6).GE.UltStress23.OR.) THEN
		ElemDeletion = 0
	  ENDIF
C	------------------------------------------------------------------

	<state variables>
C	------------------------------------------------------------------
	! for Recursive-approach, recent coefficeint of degradation should be preserved
	  STATEV(1) = Gamma11
	  STATEV(2) = Gamma22
	  STATEV(3) = Gamma33
	  STATEV(4) = Gamma12
	  STATEV(5) = Gamma13
	  STATEV(6) = Gamma23
	! for material deletion 
	  STATEV(7) = ElemDeletion
	! Flag for failure
	  !STATEV(8)
C	------------------------------------------------------------------
	  

      RETURN
      END