C*************************************************************************
! Technische Universität München
! Chair of Computational Mechanics
! Software Lab Project: Development of Failure Criteria for Composites
C*************************************************************************
! 3D ORTHOTROPIC ELATICITY WITH MAXIMUM STRESS CRITERIA AND PLY-DISCOUNT DAMAGE
! DAMAGE MODEL: INSTANT DEGREDATION
! (CAN NOT BE USED FOR 2D PROBLEMS)
! VERSION: 1.0
C-------------------------------------------------------------------------
C   GROUP: 11
C	WRITTEN BY AMMAR KHALLOUF
C   DATE: 25.10.2019
C   
C=========================================================================
!	HEADER OF THE SUBROUTINE

      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
C 


C  	USER DEFINED VARAIBLES IN VUMAT SUBROUTINE

  
 	  INTEGER fflags(nblock,6)
	  
	  DOUBLE PRECISION Dgrd(3), dmg(nblock,6), e(nblock,6), DDSDDE(6,6), Strain(nblock,6) 
	   

C=========================================================================		  

C ELASTIC PROPERTIES

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
	  
	  
C Initialize Damage indicies and Failure Flags at begining of the Analysis
C=========================================================================  
		
		IF(totalTime.EQ.0) THEN   
		
		DO I = 1,nblock   ! Loop through all material points
		
		! Assume purely elastic undamaged material at the begining of analysis
		
				dmg(I,1)= 1.0
                dmg(I,2)= 1.0
                dmg(I,3)= 1.0
                dmg(I,4)= 1.0
                dmg(I,5)= 1.0
                dmg(I,6)= 1.0
                
		! Initialize all failure flags
            
                fflags(I,1)= 0
                fflags(I,2)= 0
                fflags(I,3)= 0
                fflags(I,4)= 0
                fflags(I,5)= 0
                fflags(I,6)= 0
				
		! Initialize strain vector with strain increments
		
			    Strain(I,1)= strainInc(I,1) 
				Strain(I,2)= strainInc(I,2) 
                Strain(I,3)= strainInc(I,3) 
                Strain(I,4)= strainInc(I,4) 
                Strain(I,5)= strainInc(I,5) 
				Strain(I,6)= strainInc(I,6) 
								
		END DO
			  



C Calculate Stresses
C==================================================================

CALL Ortho3D (E11,E22,E33,G12,G13,G23,XNU12,XNU13,XNU23,XNU21,XNU31,XNU32,Upsilon,nblock,DDSDDE,Strain,dmg,stressNew)

C Update Total Strains
C=========================================================================

	ELSE
		
           DO I= 1,nblock   ! Loop through all material points

				DO J= 1,6       ! Loop through all strain components
			
     			Strain(I,J)=  Strain(I,J)+ strainInc(I,J) ! Accumulate Strain

				END DO
				
					END DO
					
C Calculate Stresses
C=========================================================================


					
CALL Ortho3D (E11,E22,E33,G12,G13,G23,XNU12,XNU13,XNU23,XNU21,XNU31,XNU32,Upsilon,nblock,DDSDDE,Strain,dmg,stressNew)

C Failure Evaluation as per Maximum Stress Criteria
C=========================================================================

		 ! Failure Evaluation as per Maximum Stress Criteria

           DO I= 1,nblock   ! Loop through all material points and calculate failure indicies

			  ! Failure along 1-direction (Tension/Compression)
				
              IF (stressNew(I,1).GE. 0.0) THEN
                
            	 e(I,1)=stressNew(I,1)/Xt
                 
        	  ELSE
               
            	e(I,1)=-stressNew(I,1)/Xc
                
        	  END IF

            ! Failure along 2-direction  (Tension/Compression)
	
        	 IF (stressNew(I,2).GE. 0.0) THEN 
            
                e(I,2)=stressNew(I,2)/Yt
				
        	 ELSE
             
            	e(I,2)=-stressNew(I,2)/Yc 
                
        	 END IF

            ! Failure along 3-direction  (Tension/Compression)       
		
        	 IF (stressNew(I,3).GE. 0.0) THEN 
            
            	e(I,3)=stressNew(I,3)/Zt
                
        	 ELSE
            
            e(I,3)=-stressNew(I,3)/Zc
            
        	 END IF

            ! Shear Failure along 1-2 direction
		
        	 IF (stressNew(I,4).GE. 0.0) THEN 
            
            	e(I,4)=stressNew(I,4)/S12
                
        	 ELSE 
            
            	e(I,4)=-stressNew(I,4)/S12
                
        	 END IF	
            	
            ! Shear Failure along 1-3 direction
		
        	 IF (stressNew(I,5).GE. 0.0) THEN
             
            	e(I,5)=stressNew(I,5)/S13
                
        	 ELSE
            
            	e(I,5)=-stressNew(I,5)/S13
                
        	 END IF

            ! Shear Failure along 2-3 direction

        	 IF (stressNew(I,6).GE. 0.0) THEN
             
           	   e(I,6)=stressNew(I,6)/S23
               
        	 ELSE
            
           	  e(I,6)=-stressNew(I,6)/S23
              
       	     END IF	
           
		   END DO
		      
		   
C Damage Application based on Failure Mode
C=========================================================================

		DO I= 1,nblock ! Loop through all material points

        	DO J=1,3   ! Loop through normal stress components 

        		IF((e(I,J).GT.1.0) .AND. (stressNew(I,J).GE.0.0)) THEN  ! Check for tensile failure

    				fflags(I,J)=totalTime    ! Assign Failure Flag to step number as a positive number

					dmg(I,J)= (0.999999-Dgrd(1))  ! Apply tensile damage 


    			ELSE IF((e(I,J).GT.1.0) .AND. (stressNew(I,J).LT.0.0)) THEN  ! Check for compressive failure

    				fflags(I,J)=-totalTime    ! Assign Failure Flag to step number as a positive number

					dmg(I,J)= (0.999999-Dgrd(2))  ! Apply compressive damage


                   END IF

					END DO
					

        	DO J=4,6   ! Loop through shear stress components 
            

				IF((e(I,J).GT.1.0)) THEN  ! Check for shear failure


    				fflags(I,J)=totalTime    ! Assign Failure Flag to step number as a positive number

					dmg(I,J)= (0.999999-Dgrd(3))  ! Apply shear damage

                   END IF

                     END DO                  	

						END DO       
		
C Calculate Stresses
C=========================================================================		
	   

           
CALL Ortho3D (E11,E22,E33,G12,G13,G23,XNU12,XNU13,XNU23,XNU21,XNU31,XNU32,Upsilon,nblock,DDSDDE,Strain,dmg,stressNew)
		
			
C Update State Variables
C=========================================================================	
 
		DO I= 1,nblock ! Loop through all material points
        	
		   DO J= 1,6   ! Loop through all damage indicies 
			
        	stateNew(I,J)= dmg(I,J)   ! Assign updated state variables related to damage indicies

			END DO


            ! Assign updated state variables related to failure flags 
			
			 stateNew(I,7)= fflags(I,1)
             
			 stateNew(I,8)= fflags(I,2)

			 stateNew(I,9)= fflags(I,3)

			 stateNew(I,10)= fflags(I,4)             
             
			 stateNew(I,11)= fflags(I,5)

			 stateNew(I,12)= fflags(I,6)          

            	END DO

			END IF

		RETURN
        END

C================================================================================================	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Utility Subroutine: Ortho3D 

! The following subroutine : 1. Build the effective constiutive matrix based on orthtropic elasticity

!							 2. Return the updated stresses

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  		
 	   SUBROUTINE Ortho3D (E11,E22,E33,G12,G13,G23,XNU12,XNU13,XNU23,XNU21,XNU31,XNU32,Upsilon,nblock,DDSDDE,Strain,dmg,Stress)
		 	 
		  INTEGER I, nblock
   		  DOUBLE PRECISION  Strain(nblock,6), DDSDDE(6,6)
   		  DOUBLE PRECISION  dmg(nblock,6), Stress(nblock,6)

        
		    Do I = 1,nblock   ! Loop through all material points 

			
			! Form effective orthotropic material matrix including damage (Jacobian Matrix) 

      		DDSDDE(1,1)=(E11*Upsilon*(1.0-XNU23*XNU32))*abs(dmg(I,1))
      		DDSDDE(1,2)=(E11*Upsilon*(XNU21+XNU31*XNU23))*abs(dmg(I,1)*dmg(I,2))
      		DDSDDE(1,3)=(E11*Upsilon*(XNU31+XNU21*XNU32))*abs(dmg(I,1)*dmg(I,3))
      		DDSDDE(2,1)=(E11*Upsilon*(XNU21+XNU31*XNU23))*abs(dmg(I,1)*dmg(I,2))
      		DDSDDE(2,2)=(E22*Upsilon*(1.0-XNU13*XNU31))*abs(dmg(I,2))
      		DDSDDE(2,3)=(E22*Upsilon*(XNU32+XNU12*XNU31))*abs(dmg(I,3)*dmg(I,2))
      		DDSDDE(3,1)=(E11*Upsilon*(XNU31+XNU21*XNU32))*abs(dmg(I,3)*dmg(I,1))
      		DDSDDE(3,2)=(E22*Upsilon*(XNU32+XNU12*XNU31))*abs(dmg(I,3)*dmg(I,2))
      		DDSDDE(3,3)=(E33*Upsilon*(1.0-XNU12*XNU21))*abs(dmg(I,3))
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


		     END DO   

			
		END SUBROUTINE Ortho3D 