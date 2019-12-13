!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
																										!	
!  Chair of Computational Mecahnics - Technical University of Munich                                    !
                                                                                                        !
!  Software Lab Project 2019: "Development of Failure Criteria For Composites"							!
																										!	
!  ##### Standalone Test Program for ABAQUS UMAT/VUMAT SUBROUTINE ##### 					     		!
																										!
!  Can be used to debug and verify the execution of constitutive law prior to applying it to ABAQUS     !
																										!
!  Failure Criterion: (Max Strain)  Damage Criteria: (Recursive Stiffness reduction)                    !                                                                                            		!
																										!	
!  Date: 25.10.2019																						!
																										!	
!  Written By: Ammar Khallouf  																			!
																									    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM UserMaterialTest


! 1.(Define Data Variables, Type and Arrays Dimension)
!=====================================================

   integer, PARAMETER :: nblock=5, ndir=3, nshr=3, nstatev=12 !(Number of Material Points, Stress Components and State Variables)
   integer :: fflags(nblock,6), STEP, NSTEPS
   real :: StrainInc(nblock,ndir+nshr), Strain(nblock, ndir+nshr),DDSDDE(ndir+nshr,ndir+nshr)
   real :: StressNew(nblock, ndir+nshr), StateNew(nblock,nstatev), e(nblock,6), Dmg(nblock,6)
   real :: E11, E22, E33, G12, G13, G23, XNU12, XNU13, XNU23, Upsilon, XNU21, XNU31, XNU32
   real :: Xt, Xc, Yt, Yc, Zt, Zc, S12, S13, S23, XbetaT ,XbetaC, XbetaS
   

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

! 2.(Input data for the User Material Subroutine)
!=========================================================================================================================

!2.1 (Material Elastic Properties)
!====================================================

! The following input is required to form the constiutive law matrix (Material Jacobain Matrix)

    E11=210 ! Elastic Modulus along 1-direction 
    E22=210   ! Elastic Moduls  along 2-direction
    E33=210   ! Elastic Modulus along 3-direction (thickness direction)
    XNU12=0.30 ! Poisson ratio along 1-2
    XNU13=0.30 ! Poisson ratio along 1-3
    XNU23=0.30 ! Poisson ratio along 2-3
    G12=80.769   ! Elastic Shear Modulus along 1-2 
    G13=80.769   ! Elastic Shear Modulus along 1-3 
    G23=80.769   ! Elastic Shear Modulus along 2-3

 ! The rest of the porperies are calculated based on previous input using reciprocal rule of consitutive material law
 
    XNU21=(E22*XNU12)/E11
    XNU31=(E33*XNU13)/E11
    XNU32=(E33*XNU23)/E22
	Upsilon=(1)/(1-XNU12*XNU21-XNU23*XNU32-XNU13*XNU31-2*XNU21*XNU32*XNU13)

!2.2 (Material Strength Properties)
!==================================================

	Xt=100.0   ! Tension Strength in normal 1-direction
    Xc=100.0  ! Compressive Strength in normal 1-direction
    Yt=0.015  ! Tension Strength in normal 2-direction
    Yc=100.0  ! Compressive Strength in normal 2-direction
    Zt=100.0  ! Tension Strength in normal 3-direction
    Zc=100.0  ! Compressive Strength in normal 3-direction
    S12=100.0 ! Shear Strength in  1-2 plane
    S13=100.0 ! Shear Strength in  1-3 plane
    S23=100.0 ! Shear Strength in  2-3 plane
    

!2.3 (Degredation Factors for Tension, Compression and Shear Failure)
!===============================================================================
    
    XbetaT=0.0001 ! Tension Stiffness Degredation Factor
    XbetaC=0.0001 ! Compression Stiffness Degredation Factor
    XbetaS=0.0001 ! Shear Stiffness Degredation Factor

            
!2.4 (Define Strain Increment and number of analysis steps)
!===============================================================================

! Assign the strain increment to be applied to the material points 
! In real FEA this increment is obtained from ABAQUS based on the calculated displacements of the solved system

!  For simplicity here we consider a uniform strain increment applied to all material points in the normal 1-direction
!  You can also assign different strain increments in different directions


	DO I = 1,nblock ! Loop through all material points 

  
	    ! Assign strain increments from ABAQUS

  		StrainInc(I,1)= -0.00000285718
		StrainInc(I,2)= 0.00000952369   
		StrainInc(I,3)= -0.00000285704
        StrainInc(I,4)= -0.000000000428766
        StrainInc(I,5)= -0.00000000000943274
        StrainInc(I,6)= 0.00000000000984941

    END DO


 ! Define number of analysis steps, in real ABAQUS FEA this represents the pseudo-time of the analysis


  	NSTEPS = 11 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 


! 3.(Analysis Phase of VUMAT Subroutine)
!=========================================================================================================================
    
  DO STEP = 1,NSTEPS  ! Analysis steps counter


!3.1 (At Begining of analysis initialize strain vector and assume elastic material with no damage)
!====================================================================================================   	
				
	    IF (STEP.EQ.1) THEN

            DO I= 1,nblock   ! Loop through all material points

              ! Initialize strain vector 

     			Strain(I,1)= 0 
				Strain(I,2)= 0
                Strain(I,3)= 0
                Strain(I,4)= 0
                Strain(I,5)= 0
				Strain(I,6)= 0
                
			    ! Initialize stress vector 
                
				StressNew(I,1)=0
				StressNew(I,2)=0
                StressNew(I,3)=0
                StressNew(I,4)=0
                StressNew(I,5)=0
 				StressNew(I,6)=0
                
			! Assume purely elastic undamaged material at the begining of analysis
						
				Dmg(I,1)= 1.0
                Dmg(I,2)= 1.0
                Dmg(I,3)= 1.0
                Dmg(I,4)= 1.0
                Dmg(I,5)= 1.0
                Dmg(I,6)= 1.0
                
			! Initialize all failure flags
            
                fflags(I,1)= 0
                fflags(I,2)= 0
                fflags(I,3)= 0
                fflags(I,4)= 0
                fflags(I,5)= 0
                fflags(I,6)= 0

             ! Initialize all stress/strength ratios
            
                e(I,1)= 0
                e(I,2)= 0
                e(I,3)= 0
                e(I,4)= 0
                e(I,5)= 0
                e(I,6)= 0


             END DO
             
!3.2 (Calculate Inital Stresses using Ortho	3D Utility Subroutine)
!==================================================================

	  CALL Ortho3D (E11,E22,E33,G12,G13,G23,XNU12,XNU13,XNU23,XNU21,XNU31,XNU32,Upsilon,DDSDDE,Strain,Dmg,StressNew)


!3.3 (Update  Total strains)
!==================================================================
			
		ELSE

           DO I= 1,nblock   ! Loop through all material points and accumulate strains
     		
     			Strain(I, 1)=  Strain(I,1)+ StrainInc(I,1) 
                Strain(I, 2)=  Strain(I,2)+ StrainInc(I,2) 
                Strain(I, 3)=  Strain(I,3)+ StrainInc(I,3) 
                Strain(I, 4)=  Strain(I,4)+ StrainInc(I,4) 
                Strain(I, 5)=  Strain(I,5)+ StrainInc(I,5) 
 				Strain(I, 6)=  Strain(I,6)+ StrainInc(I,6)

				END DO

!3.4 (Update stresses)
!==================================================================

          ! Calculate Updated Stresses
           
       	   CALL Ortho3D (E11,E22,E33,G12,G13,G23,XNU12,XNU13,XNU23,XNU21,XNU31,XNU32,Upsilon,DDSDDE,Strain,Dmg,StressNew)
        
!3.5 (Failure Evaluation)
!==================================================================

		 ! Failure Evaluation as per Maximum Strain Criteria

           DO I= 1,nblock   ! Loop through all material points and calculate failure indicies

			  ! Failure along 1-direction (Tension/Compression)
				
              IF (Strain(I,1).GE. 0.0) THEN
                
            	 e(I,1)=Strain(I,1)/(Xt/E11)
                 
        	  ELSE
               
            	e(I,1)=-Strain(I,1)/(Xc/E11)
                
        	  END IF

            ! Failure along 2-direction  (Tension/Compression)
	
        	 IF (Strain(I,2).GE. 0.0) THEN 
            
                e(I,2)=Strain(I,2)/(Yt/E22) 
        	 ELSE
             
            	e(I,2)=-Strain(I,2)/(Yc/E22) 
                
        	 END IF

            ! Failure along 3-direction  (Tension/Compression)       
		
        	 IF (Strain(I,3).GE. 0.0) THEN 
            
            	e(I,3)=Strain(I,3)/(Zt/E33)
                
        	 ELSE
            
            e(I,3)=-Strain(I,3)/(Zc/E33)
            
        	 END IF

            ! Shear Failure along 1-2 direction
		
        	 IF (Strain(I,4).GE. 0.0) THEN 
            
            	e(I,4)=Strain(I,4)/(S12/G12)
                
        	 ELSE 
            
            	e(I,4)=-Strain(I,4)/(S12/G12)
                
        	 END IF	
            	
            ! Shear Failure along 1-3 direction
		
        	 IF (Strain(I,5).GE. 0.0) THEN
             
            	e(I,5)=Strain(I,5)/(S13/G13)
                
        	 ELSE
            
            	e(I,5)=-Strain(I,5)/(S13/G13)
                
        	 END IF

            ! Shear Failure along 2-3 direction

        	 IF (Strain(I,6).GE. 0.0) THEN
             
           	   e(I,6)=Strain(I,6)/(S23/G23)
               
        	 ELSE
            
           	  e(I,6)=-Strain(I,6)/(S23/G23)
              
       	     END IF	
           
		   END DO
	  	 
!3.6 (Damage Application based on failure mode)
!==================================================================

		DO I= 1,nblock ! Loop through all material points

        	DO J=1,3   ! Loop through normal stress components 

        		IF(((e(I,J).GT.1.0).OR.(Dmg(I,J).LT.1.0)) .AND. (Strain(I,J).GE.0.0)) THEN  ! Check for tensile failure

    				fflags(I,J)=STEP    ! Assign Failure Flag to step number as a positive number

					Dmg(I,J)= Dmg(I,J)*(0.999999-XbetaT)  ! Apply tensile damage 


    			ELSE IF(((e(I,J).GT.1.0).OR.(Dmg(I,J).LT.1.0)) .AND. (Strain(I,J).LT.0.0)) THEN  ! Check for compressive failure

    				fflags(I,J)=-STEP    ! Assign Failure Flag to step number as a negative number

					Dmg(I,J)= Dmg(I,J)*(0.999999-XbetaC)  ! Apply compressive damage


                   END IF

                     END DO

        	DO J=4,6   ! Loop through shear stress components 
            

				IF(((e(I,J).GT.1.0).OR.(Dmg(I,J).LT.1.0))) THEN  ! Check for shear failure


    				fflags(I,J)=STEP    ! Assign Failure Flag to step number as a positive number

					Dmg(I,J)= Dmg(I,J)*(0.999999-XbetaS)  ! Apply shear damage

                   END IF

                     END DO                  	

						END DO

!3.7 (Update Stresses)
!==================================================================        

          ! Calculate Updated Stresses
           
       	   CALL Ortho3D (E11,E22,E33,G12,G13,G23,XNU12,XNU13,XNU23,XNU21,XNU31,XNU32,Upsilon,DDSDDE,Strain,Dmg,StressNew)
		
!3.8 (Update State Variables)
!================================================================== 
 
		DO I= 1,nblock ! Loop through all material points
        	
		   DO J= 1,6   ! Loop through all damage indicies 
			
        	StateNew(I,J)= Dmg(I,J)   ! Assign updated state variables related to damage indicies

			END DO


            ! Assign updated state variables related to failure flags 
			
			 StateNew(I,7)= fflags(I,1)
             
			 StateNew(I,8)= fflags(I,2)

			 StateNew(I,9)= fflags(I,3)

			 StateNew(I,10)= fflags(I,4)             
             
			 StateNew(I,11)= fflags(I,5)

			 StateNew(I,12)= fflags(I,6)
             


            	END DO



                	END IF


!3.9 (Output Results)
!================================================================== 

 Print *, 'Step No:', STEP  ! Analysis step number
 Print *, '****************************************************************'
 
 WRITE(*,*)  ! print blank line

 print '(" Failure ratio in normal 2-direction: ",f6.3)',e(1,2) ! print stress/strength ratio (2) for material point number 1
 Print *, '________________________________________________________________'
 
 Print *, 'Damage Factor in normal 2-direction: ', Dmg(1,2) ! print damage (2) factor for material point number 1
 Print *, '________________________________________________________________'
 
 print '(" Stress in normal 2-direction: ",f6.3)',StressNew(1,2) ! print Stress (2) for material point number 1
 Print *, '________________________________________________________________'
 
 Print *, 'Strain in normal 2-direction: ', Strain(1,2)     ! print Stain (2) for material point number 1
 Print *, '________________________________________________________________'
 
 Print *, 'Strain Increment in normal 2-direction: ', StrainInc(1,2)  ! print Stain increment (2) for material point number 1
 Print *, '________________________________________________________________'
 

 WRITE(*,*)
 WRITE(*,*)


                    	END DO  ! END ANALYSIS STEPS COUNTER (END OF ANALYSIS)

                        
               
 END PROGRAM UserMaterialTest
  
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Utility Subroutine: Ortho3D 

! The following subroutine : 1. Build the effective constiutive matrix based on orthtropic elasticity

!							 2. Return the updated stresses

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  		
 		SUBROUTINE Ortho3D (E11,E22,E33,G12,G13,G23,XNU12,XNU13,XNU23,XNU21,XNU31,XNU32,Upsilon,DDSDDE,Strain,Dmg,Stress)

   		integer, PARAMETER :: nblock=5, ndir=3, nshr=3
   		integer :: I
   		real :: Strain(nblock, ndir+nshr),DDSDDE(ndir+nshr,ndir+nshr)
   		real :: Stress(nblock, ndir+nshr), Dmg(nblock,6)

        
		Do I = 1,nblock   ! Loop through all material points 

			
			! Form effective orthotropic material matrix including damage (Jacobian Matrix) 

      		DDSDDE(1,1)=(E11*Upsilon*(1.0-XNU23*XNU32))*abs(Dmg(I,1))
      		DDSDDE(1,2)=(E11*Upsilon*(XNU21+XNU31*XNU23))*abs(Dmg(I,1)*Dmg(I,2))
      		DDSDDE(1,3)=(E11*Upsilon*(XNU31+XNU21*XNU32))*abs(Dmg(I,1)*Dmg(I,3))
      		DDSDDE(2,1)=(E11*Upsilon*(XNU21+XNU31*XNU23))*abs(Dmg(I,1)*Dmg(I,2))
      		DDSDDE(2,2)=(E22*Upsilon*(1.0-XNU13*XNU31))*abs(Dmg(I,2))
      		DDSDDE(2,3)=(E22*Upsilon*(XNU32+XNU12*XNU31))*abs(Dmg(I,3)*Dmg(I,2))
      		DDSDDE(3,1)=(E11*Upsilon*(XNU31+XNU21*XNU32))*abs(Dmg(I,3)*Dmg(I,1))
      		DDSDDE(3,2)=(E22*Upsilon*(XNU32+XNU12*XNU31))*abs(Dmg(I,3)*Dmg(I,2))
      		DDSDDE(3,3)=(E33*Upsilon*(1.0-XNU12*XNU21))*abs(Dmg(I,3))
      		DDSDDE(4,4)=G12*abs(Dmg(I,4))
      		DDSDDE(5,5)=G13*abs(Dmg(I,5))
      		DDSDDE(6,6)=G23*abs(Dmg(I,6))

			
			! Calculate Stress Components
            
      		Stress(I,1)= (DDSDDE(1,1)*Strain(I,1)+DDSDDE(1,2)*Strain(I,2)+DDSDDE(1,3)*Strain(I,3))
	  		Stress(I,2)= (DDSDDE(2,1)*Strain(I,1)+DDSDDE(2,2)*Strain(I,2)+DDSDDE(2,3)*Strain(I,3))
	  		Stress(I,3)= (DDSDDE(3,1)*Strain(I,1)+DDSDDE(3,2)*Strain(I,2)+DDSDDE(3,3)*Strain(I,3)) 
	  		Stress(I,4)= (DDSDDE(4,4)*Strain(I,4))
	  		Stress(I,5)= (DDSDDE(5,5)*Strain(I,5))   	  
	  		Stress(I,6)= (DDSDDE(6,6)*Strain(I,6)) 


            		END DO   
            

      		END SUBROUTINE Ortho3D   