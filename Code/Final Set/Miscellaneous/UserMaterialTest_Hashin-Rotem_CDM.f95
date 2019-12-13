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
!  Failure Criterion: (Hashin)  Damage Criteria: (Continium Damage Mechanics- Crack Band Theory)        !                                                                                            		!
																										!	
!  Date: 25.10.2019																						!
																										!	
!  Written By: Ammar Khallouf  																			!
																									    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM UserMaterialTest


! 1.(Define Data Variables, Type and Arrays Dimension)
!=====================================================

   integer, PARAMETER :: nblock=5, ndir=3, nshr=3, nstatev=9 !(Number of Material Points, Stress Components and State Variables)
   integer :: fflags(nblock,3), STEP, NSTEPS
   real :: StrainInc(nblock,ndir+nshr), Strain(nblock, ndir+nshr),DDSDDE(ndir+nshr,ndir+nshr)
   real :: StressNew(nblock, ndir+nshr), StateNew(nblock,nstatev), e(nblock,3), d_index(nblock,6), Dmg(nblock,6)
   real :: E11, E22, E33, G12, G13, G23, XNU12, XNU13, XNU23, Upsilon, XNU21, XNU31, XNU32
   real :: Xt, Xc, Yt, Yc, Zt, Zc, S12, S13, S23, le, let, SD(6), E11D, E22D, E33D, G12D, G13D, G23D ! Strength properites and degraded moduli
   real :: G_ft, G_mt, G_fc, G_mc, G_IIC ! fiber and matrix fracture toughness, Mode II critical energy release for shear

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

! 2.(Input data for the User Material Subroutine)
!=========================================================================================================================

!2.1 (Material Elastic Properties)
!====================================================

! The following input is required to form the constiutive law matrix (Material Jacobain Matrix)

    E11=164.0 ! Elastic Modulus along 1-direction 
    E22=8.98   ! Elastic Moduls  along 2-direction
    E33=8.98   ! Elastic Modulus along 3-direction (thickness direction)
    XNU12=0.32 ! Poisson ratio along 1-2
    XNU13=0.32 ! Poisson ratio along 1-3
    XNU23=0.496 ! Poisson ratio along 2-3
    G12=5.02   ! Elastic Shear Modulus along 1-2 
    G13=5.02   ! Elastic Shear Modulus along 1-3 
    G23=3.00   ! Elastic Shear Modulus along 2-3

 ! The rest of the porperies are calculated based on previous input using reciprocal rule of consitutive material law
 
    XNU21=(E22*XNU12)/E11
    XNU31=(E33*XNU13)/E11
    XNU32=(E33*XNU23)/E22
	Upsilon=(1)/(1-XNU12*XNU21-XNU23*XNU32-XNU13*XNU31-2*XNU21*XNU32*XNU13)

!2.2 (Material Strength Properties)
!==================================================

	Xt=2.90   ! Tension Strength in normal 1-direction
    Xc=1.68   ! Compressive Strength in normal 1-direction
    Yt=0.1    ! Tension Strength in normal 2-direction
    Yc=0.247  ! Compressive Strength in normal 2-direction
    Zt=0.1    ! Tension Strength in normal 3-direction
    Zc=0.247  ! Compressive Strength in normal 3-direction
    S12=0.08  ! Shear Strength in  1-2 plane
    S13=0.08 ! Shear Strength in  1-3 plane
    S23=0.08 ! Shear Strength in  2-3 plane
    

!2.3 (Continium Damage Mechanics Input variables)
!===============================================================================
    
	G_ft = 0.081534      ! Fiber Tension fracture toughness
	
	G_mt = 0.000256     ! Matrix Tension fracture toughness

    G_fc = 0.024533     ! Fiber Compression fracture toughness
	
	G_mc = 0.001156      		! Matrix Compression fracture toughness
	
	G_IIC = 0.001156     ! Mode II critical energy release rate for shear

    le = 4         ! Length scale for ply in-plane strains

    let = 0.127        ! Length scale for ply transverse strains


    ! alpha and beta  are correction factors for nonlinear shear stiffness degredation

    ! They take a value between (0-1)

    ! When alpha and beta == 0 no effect on shear stiffness degredation
    
	! When alpha and beta == 1 maximum effect on shear stiffness degredation

    
    alpha = 0.0     

    beta  = 0.0

            
!2.4 (Define Strain Increment and number of analysis steps)
!===============================================================================

! Assign the strain increment to be applied to the material points 
! In real FEA this increment is obtained from ABAQUS based on the calculated displacements of the solved system

!  For simplicity here we consider a uniform strain increment applied to all material points in the normal 1-direction
!  You can also assign different strain increments in different directions


	DO I = 1,nblock ! Loop through all material points 

  
	    ! Assign strain increments from ABAQUS

  		StrainInc(I,1)= 0!0.001d+0
		StrainInc(I,2)= 0.0007d+0 
		StrainInc(I,3)= 0.0d+0
        StrainInc(I,4)= 0.0d+0
        StrainInc(I,5)= 0.0d+0
        StrainInc(I,6)= 0.0d+0

    END DO


 ! Define number of analysis steps, in real ABAQUS FEA this represents the pseudo-time of the analysis


  	NSTEPS = 20

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

                d_index(I,1)=0.0
                d_index(I,2)=0.0
				d_index(I,3)=0.0
                d_index(I,4)=0.0
				d_index(I,5)=0.0
				d_index(I,6)=0.0


                
			! Initialize all failure flags
            
                fflags(I,1)= 0
                fflags(I,2)= 0
                fflags(I,3)= 0


             ! Initialize all stress/strength ratios
            
                e(I,1)= 0
                e(I,2)= 0
                e(I,3)= 0



             END DO
             
!3.2 (Calculate Inital Stresses using Ortho	3D Utility Subroutine)
!==================================================================

	  CALL Ortho3D (E11,E22,E33,G12,G13,G23,XNU12,XNU13,XNU23,XNU21,XNU31,XNU32,Upsilon,DDSDDE,Strain,Dmg,StressNew)


!3.3 (Update  Total strains)
!==================================================================
			
		ELSE

           DO I= 1,nblock   ! Loop through all material points and accumulate strains
     		
     			Strain(I,1)=  Strain(I,1)+ StrainInc(I,1) 
                Strain(I,2)=  Strain(I,2)+ StrainInc(I,2) 
                Strain(I,3)=  Strain(I,3)+ StrainInc(I,3) 
                Strain(I,4)=  Strain(I,4)+ StrainInc(I,4) 
                Strain(I,5)=  Strain(I,5)+ StrainInc(I,5) 
 				Strain(I,6)=  Strain(I,6)+ StrainInc(I,6)

				END DO

!3.4 (Update stresses)
!==================================================================

          ! Calculate Updated Stresses
           
       	   CALL Ortho3D (E11,E22,E33,G12,G13,G23,XNU12,XNU13,XNU23,XNU21,XNU31,XNU32,Upsilon,DDSDDE,Strain,Dmg,StressNew)
        
!3.5 (Failure Evaluation Hashin-Rotem)
!==================================================================

		 ! Failure Evaluation as per Hashin-Rotem Strain Based Criteria
         

           DO I= 1,nblock   ! Loop through all material points and calculate failure indicies

			  ! Fiber Failure (Similar in Tension/Compression)
				
              IF (Strain(I,1).GE. 0.0) THEN
                
            	 e(I,1)= (Strain(I,1)/(Xt/E11))**2
                 
        	  ELSE
                          
            	e(I,1)= (Strain(I,1)/(Xc/E11))**2
                
        	  END IF
             
            !  In-Plane Matrix Failure (Similar in Tension/Compression)
	
        	 IF (Strain(I,2).GE. 0.0) THEN 
            
                e(I,2)= (Strain(I,2)/(Yt/E22))**2 + (Strain(I,6)/(S23/G23))**2 + (Strain(I,4)/(S12/G12))**2

                
        	 ELSE
           
                e(I,2)= (Strain(I,2)/(Yc/E22))**2 + (Strain(I,6)/(S23/G23))**2 + (Strain(I,4)/(S12/G12))**2
                
        	 END IF
             

            ! Transverse Matrix Failure (Similar in Tension/Compression)
		
        	 IF (Strain(I,3).GE. 0.0) THEN 
            
            	e(I,3)= (Strain(I,3)/(Zt/E33))**2 + (Strain(I,6)/(S23/G23))**2 +(Strain(I,5)/(S13/G13))**2
                
        	 ELSE

            
            	e(I,3)= (Strain(I,3)/(Zc/E33))**2 + (Strain(I,6)/(S23/G23))**2 +(Strain(I,5)/(S13/G13))**2
            
        	 END IF
	
           
		   END DO
	  	 
!3.6 (Damage Application based on Continium Damage Mechanics- Crack Band Theory)
!================================================================================

	IF (maxval(e) .GT. 1.0) THEN  !if any of the failure inidcies excced one initiate damage

		DO I= 1,nblock ! Loop through all material points check for fiber/matrix failure and calculate damage indicies
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		   ! Fiber Damage Modes

           IF (e(I,1).GT.1.0 .AND. (Strain(I,1).GE.0.0)) THEN  ! Check for fiber Tension failure

    				fflags(I,1)=STEP    ! Assign Failure Flag to step number (Positive)

					! Calculate differnce between incremented strain and damage initiation strain component

                    SD(1) = Strain(I,1)-(Xt/E11)
            
                    ! Calculate Degraded E11 modulus as per Eq (12)
					
					E11D = ((1.0/E11)+ SD(1)/(Xt*(1.0-(le*Xt*SD(1))/(2*G_ft))))**(-1)
					
					d_index(I,1) = 1.0 - (E11D/E11) ! Calculate damage index as per Eq (18)

           ELSE IF (e(I,1).GT.1.0 .AND. (Strain(I,1).LT.0.0)) THEN  ! Check for fiber Compression failure

    				fflags(I,1)=-STEP    ! Assign Failure Flag to step number (Negative)
                   
					! Calculate differnce between incremented strain and damage initiation strain component

                    SD(1) = abs(Strain(I,1))-(Xc/E11)
                    
                    ! Calculate Degraded E11 modulus as per Eq (12)
					
					E11D = ((1.0/E11)+ SD(1)/(Xc*(1.0-(le*Xc*SD(1))/(2*G_fc))))**(-1)
					
					d_index(I,1) = 1.0 - (E11D/E11) ! Calculate damage index as per Eq (18)

                   

           END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           

		   ! In-Plane Matrix Damage Modes

           IF(e(I,2).GT.1.0 .AND. (Strain(I,2).GE.0.0)) THEN  ! Check for matrix tensile failure 

    				fflags(I,2)=STEP    ! Assign Failure Flag to step number (Positive)

					! Calculate differnce between incremented strain and damage initiation strain component

                    SD(2) = Strain(I,2)-(Yt/E22)

                    SD(4) = abs(Strain(I,4))-(S12/G12)

					SD(6) = abs(Strain(I,6))-(S23/G23)                    

                    ! Calculate Degraded E22,G12,G23 modulI as per Eqs (13,15,17)                    
                    
					E22D = ((1.0/E22)+ SD(2)/(Yt*(1.0-(le*Yt*SD(2))/(2*G_mt))))**(-1)
                      
                    G12D = ((1.0/G12)+ SD(4)/(2*S12*(1.0-(le*S12*SD(4))/(4*G_IIC))))**(-1)
                                      
                    G23D = ((1.0/G23)+ SD(6)/(2*S23*(1.0-(let*S23*SD(6))/(4*G_IIC))))**(-1)
                                        
				    ! Calculate damage indicies as per Eqs (19,21,23)                   

                    d_index(I,2) = 1.0 - (E22D/E22)

                    d_index(I,4) = 1.0 - (G12D/G12)  					                   

                    d_index(I,6) = 1.0 - (G23D/G23)  

                   
           ELSE IF(e(I,2).GT.1.0 .AND. (Strain(I,2).LT.0.0)) THEN  ! Check for matrix compressive failure 

    				fflags(I,2)=-STEP    ! Assign Failure Flag to step number (Negative)

					! Calculate differnce between incremented strain and damage initiation strain component

                    SD(2) = abs(Strain(I,2))-(Yc/E22)

                    SD(4) = abs(Strain(I,4))-(S12/G12)

					SD(6) = abs(Strain(I,6))-(S23/G23)                    

                    ! Calculate Degraded E22,G12,G23 modulI as per Eqs (13,15,17)                    
                    
					E22D = ((1.0/E22)+ SD(2)/(Yc*(1.0-(le*Yc*SD(2))/(2*G_mc))))**(-1)
                      
                    G12D = ((1.0/G12)+ SD(4)/(2*S12*(1.0-(le*S12*SD(4))/(4*G_IIC))))**(-1)
                                      
                    G23D = ((1.0/G23)+ SD(6)/(2*S23*(1.0-(let*S23*SD(6))/(4*G_IIC))))**(-1)
					
				    ! Calculate damage indicies as per Eqs (19,21,23)                   

                    d_index(I,2) = 1.0 - (E22D/E22)

                    d_index(I,4) = 1.0 - (G12D/G12)  					                   

                    d_index(I,6) = 1.0 - (G23D/G23)


           END IF



           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           

		   ! Transverse Matrix Damage Modes

           IF((e(I,3).GT.1.0) .AND. (Strain(I,3).GE.0.0)) THEN  ! Check for matrix tensile failure

    				fflags(I,3)=STEP    ! Assign Failure Flag to step number (Positive)

					! Calculate differnce between incremented strain and damage initiation strain component

                    SD(3) = Strain(I,3)-(Zt/E33)

                    SD(5) = abs(Strain(I,5))-(S12/G13)

					SD(6) = abs(Strain(I,6))-(S23/G23)                     

                   ! Calculate Degraded (E33,G23,G13) moduli as per Eqs (14,16,17)
                    
					E33D = ((1.0/E33)+ SD(3)/(Zt*(1.0-(let*Zt*SD(3))/(2*G_mt))))**(-1)
                    
					G23D = ((1.0/G23)+ SD(6)/(2*S23*(1.0-(let*S23*SD(6))/(4*G_IIC))))**(-1)
                    
					G13D = ((1.0/G13)+ SD(5)/(2*S13*(1.0-(let*S13*SD(6))/(4*G_IIC))))**(-1)                    
				
				   ! Calculate damage indicies as per Eqs (20,22,23) 
                   
                    d_index(I,3) = 1.0 - (E33D/E33)

                    d_index(I,5) = 1.0 - (G13D/G13)  					                   

                    d_index(I,6) = 1.0 - (G23D/G23)					
					
    	   ELSE IF((e(I,3).GT.1.0) .AND. (Strain(I,3).LT.0.0)) THEN  ! Check for matrix compressive failure
           
    				fflags(I,3)=-STEP    ! Assign Failure Flag to step number (Negative)
                    
					! Calculate differnce between incremented strain and damage initiation strain component

                    SD(3) = abs(Strain(I,3))-(Zc/E33)

                    SD(5) = abs(Strain(I,5))-(S12/G13)

					SD(6) = abs(Strain(I,6))-(S23/G23)                     

                   ! Calculate Degraded (E33,G23,G13) moduli as per Eqs (14,16,17)
                    
					E33D = ((1.0/E33)+ SD(3)/(Zc*(1.0-(let*Zc*SD(3))/(2*G_mc))))**(-1)
                    
					G23D = ((1.0/G23)+ SD(6)/(2*S23*(1.0-(let*S23*SD(6))/(4*G_IIC))))**(-1)
                    
					G13D = ((1.0/G13)+ SD(5)/(2*S13*(1.0-(let*S13*SD(5))/(4*G_IIC))))**(-1)                  
				
				   ! Calculate damage indicies as per Eqs (20,22,23) 
                   
                    d_index(I,3) = 1.0 - (E33D/E33)

                    d_index(I,5) = 1.0 - (G13D/G13)  					                   

                    d_index(I,6) = 1.0 - (G23D/G23)	

           END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           

	   ! Assign Constiutive Matrix Damage Variables as per Eqs (28-33)

       ! Remark: Damage Indicies grows from 0 (No damage) to 1 (Full Damage)
       
       !       : Damage Variables degrade form 1 (No Degredation) to 0 (Full Degredation)
				

				Dmg(I,1) = abs(0.9999 - (d_index(I,1)))
        		Dmg(I,2) = abs(0.9999 - (d_index(I,2)))
        		Dmg(I,3) = abs(0.9999 - (d_index(I,3)))               
           		Dmg(I,4) = abs(0.9999 - ((1.0 - Dmg(I,1))*(1.0-alpha*Dmg(I,2))*(1.0-beta*d_Index(I,4))))
                Dmg(I,5) = abs(0.9999 - ((1.0 - Dmg(I,1))*(1.0-alpha*Dmg(I,3))*(1.0-beta*d_Index(I,5))))
                Dmg(I,6) = abs(0.9999 - ((1.0 - Dmg(I,1))*(1.0-alpha*Dmg(I,2))*(1.0-alpha*Dmg(I,3))*(1.0-beta*d_Index(I,5))))
                

		 END DO

	END IF

    
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

            	END DO

                	END IF


!3.9 (Output Results)
!================================================================== 

 Print *, 'Step No:', STEP  ! Analysis step number
 Print *, '****************************************************************'
 
 WRITE(*,*)  ! print blank line

 print '(" Failure ratio in normal 2-direction: ",f6.3)',e(1,2) ! print failure ratio for material point number 1
 Print *, '________________________________________________________________'
 
 Print *, 'Damage Factor in normal 2-direction: ', Dmg(1,2) ! print damage (2) factor for material point number 1
 Print *, '________________________________________________________________'
 
 print '(" Stress in normal 2-direction: ",f6.3)',StressNew(1,2) ! print Stress (2) for material point number 1
 Print *, '________________________________________________________________'
 
 Print *, 'Strain in normal 1-direction: ', Strain(1,2)     ! print Stain (2) for material point number 1
 Print *, '________________________________________________________________'
 
 Print *, 'Strain Increment in normal 1-direction: ', StrainInc(1,1)  ! print Stain increment (2) for material point number 1
 Print *, '________________________________________________________________'
 

 WRITE(*,*)
 WRITE(*,*)


                    	END DO  ! END ANALYSIS STEPS COUNTER (END OF ANALYSIS)

                        
               
 END PROGRAM UserMaterialTest
  
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Utility Subroutine: Ortho3D 

! The following subroutine : 1. Build the effective constiutive matrix based on orthtoropic elasticity and calculated damage 
!        						In case of no damage (i.e Dmg(nblock,6)==1)the subroutine returns elastic constitutive tensor

!							 2. Return the updated stresses

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  		
 		SUBROUTINE Ortho3D (E11,E22,E33,G12,G13,G23,XNU12,XNU13,XNU23,XNU21,XNU31,XNU32,Upsilon,DDSDDE,Strain,Dmg,Stress)

   		integer, PARAMETER :: nblock=5, ndir=3, nshr=3
   		integer :: I
   		real :: Strain(nblock, ndir+nshr),DDSDDE(ndir+nshr,ndir+nshr)
   		real :: Stress(nblock, ndir+nshr), Dmg(nblock,6)        

        
		Do I = 1,nblock   ! Loop through all material points 

			Upsilon=(1)/(1-Dmg(I,1)*Dmg(I,2)*XNU12*XNU21&
                    &-Dmg(I,2)*Dmg(I,3)*XNU23*XNU32&
                    &-Dmg(I,1)*Dmg(I,3)*XNU13*XNU31&
                    &-2*Dmg(I,1)*Dmg(I,2)*Dmg(I,3)*XNU21*XNU32*XNU13)

                    
			
			! Form effective orthotropic material matrix including damage (Jacobian Matrix) 

      		DDSDDE(1,1)=(E11*Upsilon*(1.0-Dmg(I,2)*Dmg(I,3)*XNU23*XNU32))*Dmg(I,1)

            DDSDDE(2,2)=(E22*Upsilon*(1.0-Dmg(I,1)*Dmg(I,3)*XNU13*XNU31))*Dmg(I,2)

            DDSDDE(3,3)=(E33*Upsilon*(1.0-Dmg(I,1)*Dmg(I,2)*XNU12*XNU21))*Dmg(I,3)
                        
      		DDSDDE(1,2)=(E11*Upsilon*(XNU21+XNU31*XNU23*Dmg(I,3)))*Dmg(I,1)*Dmg(I,2)
                                   
      		DDSDDE(1,3)=(E11*Upsilon*(XNU31+XNU21*XNU32*Dmg(I,2)))*Dmg(I,1)*Dmg(I,3)
                 		     		                  		            
      		DDSDDE(2,3)=(E22*Upsilon*(XNU32+Dmg(I,1)*XNU12*XNU31))*Dmg(I,3)*Dmg(I,2)

            DDSDDE(2,1)=DDSDDE(1,2)
            
            DDSDDE(3,1)=DDSDDE(1,3) 
            
      		DDSDDE(3,2)=DDSDDE(2,3)
                 		            
      		DDSDDE(4,4)=G12*Dmg(I,4)
            
      		DDSDDE(5,5)=G13*Dmg(I,5)
            
      		DDSDDE(6,6)=G23*Dmg(I,6)

			
			! Calculate Stress Components
            
      		Stress(I,1)= (DDSDDE(1,1)*Strain(I,1)+DDSDDE(1,2)*Strain(I,2)+DDSDDE(1,3)*Strain(I,3))
	  		Stress(I,2)= (DDSDDE(2,1)*Strain(I,1)+DDSDDE(2,2)*Strain(I,2)+DDSDDE(2,3)*Strain(I,3))
	  		Stress(I,3)= (DDSDDE(3,1)*Strain(I,1)+DDSDDE(3,2)*Strain(I,2)+DDSDDE(3,3)*Strain(I,3)) 
	  		Stress(I,4)= (DDSDDE(4,4)*Strain(I,4))
	  		Stress(I,5)= (DDSDDE(5,5)*Strain(I,5))   	  
	  		Stress(I,6)= (DDSDDE(6,6)*Strain(I,6)) 


            		END DO   
            

      		END SUBROUTINE Ortho3D   