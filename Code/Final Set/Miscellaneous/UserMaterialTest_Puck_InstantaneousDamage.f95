!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                                        !   
!  Chair of Computational Mecahnics - Technical University of Munich                                    !
                                                                                                        !
!  Software Lab Project 2019: "Development of Failure Criteria For Composites"                          !
                                                                                                        !   
!  ##### Standalone Test Program for ABAQUS UMAT/VUMAT SUBROUTINE #####                                 !
                                                                                                        !
!  Can be used to debug and verify the execution of constitutive law prior to applying it to ABAQUS     !
                                                                                                        !
!  Failure Criterion: (Hashin)  Damage Criteria: (Instantaneous Stiffness reduction)                        !                                                                                                   !
                                                                                                        !   
!  Date: 25.10.2019                                                                                     !
                                                                                                        !   
!  Written By: Ammar Khallouf                                                                           !
                                                                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM UserMaterialTest


! 1.(Define Data Variables, Type and Arrays Dimension)
!=====================================================

   integer, PARAMETER :: nblock=5, ndir=3, nshr=3, nstatev=9 !(Number of Material Points, Stress Components and State Variables)
   integer :: fflags(nblock,3), STEP, NSTEPS
   real :: StrainInc(nblock,ndir+nshr), Strain(nblock, ndir+nshr),DDSDDE(ndir+nshr,ndir+nshr)
   real :: StressNew(nblock, ndir+nshr), StateNew(nblock,nstatev), e(nblock,3), Dmg(nblock,6)
   real :: E11, E22, E33, G12, G13, G23, XNU12, XNU13, XNU23, Upsilon, XNU21, XNU31, XNU32
   real :: Xt, Xc, Yt, Yc, Zt, Zc, S12, S13, S23, XbetaT ,XbetaC, XbetaS
   real :: FF, XNU12f, E11f, A, MGF, Pn1, Pnt, THETA, THETAMAX, THETAF, IFF
   real :: SFP, TNT, TN1, Rn, Rnt, Rn1, FMAX
   real, PARAMETER :: Pi = 3.14159265358979d+0
   

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

! 2.(Input data for the User Material Subroutine)
!=========================================================================================================================

!2.1 (Material Elastic Properties)
!====================================================

! The following input is required to form the constiutive law matrix (Material Jacobain Matrix)

    E11=117.212 ! Elastic Modulus along 1-direction 
    E22=9.653   ! Elastic Moduls  along 2-direction
    E33=9.653  ! Elastic Modulus along 3-direction (thickness direction)
    XNU12=0.34 ! Poisson ratio along 1-2
    XNU13=0.34 ! Poisson ratio along 1-3
    XNU23=0.34 ! Poisson ratio along 2-3
    G12=3.103   ! Elastic Shear Modulus along 1-2 
    G13=3.103   ! Elastic Shear Modulus along 1-3 
    G23=3.103   ! Elastic Shear Modulus along 2-3

 ! The rest of the porperies are calculated based on previous input using reciprocal rule of consitutive material law
 
    XNU21=(E22*XNU12)/E11
    XNU31=(E33*XNU13)/E11
    XNU32=(E33*XNU23)/E22
    Upsilon=(1)/(1-XNU12*XNU21-XNU23*XNU32-XNU13*XNU31-2*XNU21*XNU32*XNU13)

!2.2 (Material Strength Properties)
!==================================================

    Xt=1378.96/1000   ! Tension Strength in normal 1-direction
    Xc=689.48/1000  ! Compressive Strength in normal 1-direction
    Yt=37.5  ! Tension Strength in normal 2-direction
    Yc=130.3 ! Compressive Strength in normal 2-direction
    Zt=100.0  ! Tension Strength in normal 3-direction
    Zc=100.0  ! Compressive Strength in normal 3-direction
    S12=66.5 ! Shear Strength in  1-2 plane
    S13=66.5 ! Shear Strength in  1-3 plane
    S23=66.5 ! Shear Strength in  2-3 plane
    THETAF = 0.92502450356  ! Angle of fracture plane for puck criterion (in radians)
    E11f= 117.212
    XNU12f=0.2

!2.3 (Degredation Factors for Tension, Compression and Shear Failure)
!===============================================================================
    
    XbetaT=0.5 ! Tension Stiffness Degredation Factor
    XbetaC=0.7 ! Compression Stiffness Degredation Factor
    XbetaS=0.5 ! Shear Stiffness Degredation Factor
    MGF= 1.2

            
!2.4 (Define Strain Increment and number of analysis steps)
!===============================================================================

! Assign the strain increment to be applied to the material points 
! In real FEA this increment is obtained from ABAQUS based on the calculated displacements of the solved system

!  For simplicity here we consider a uniform strain increment applied to all material points in the normal 1-direction
!  You can also assign different strain increments in different directions


    DO I = 1,nblock ! Loop through all material points 

  
        ! Assign strain increments from ABAQUS

        StrainInc(I,1)= 0
        StrainInc(I,2)= -0.008  
        StrainInc(I,3)= 0
        StrainInc(I,4)= 0
        StrainInc(I,5)= 0
        StrainInc(I,6)= 0

    END DO


 ! Define number of analysis steps, in real ABAQUS FEA this represents the pseudo-time of the analysis


    NSTEPS =2

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


             ! Initialize all stress/strength ratios
            
                e(I,1)= 0
                e(I,2)= 0
                e(I,3)= 0

                A=0.0
                THETA= 0.0
                THETAMAX=0.0
                IFF =0.0
                SFP=0.0
                TNT=0.0
                TN1=0.0
                Pn1=0.0
                Pnt=0.0
                Rn=0.0
                Rnt= 0.0
                Rn1=0.0
                FMAX =0.0              

             END DO
             
!3.2 (Calculate Inital Stresses using Ortho 3D Utility Subroutine)
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

        
!3.5 (Failure Evaluation Hashin)
!==================================================================

         ! Failure Evaluation as per Puck 3D Criteria
         

           DO I= 1,nblock   ! Loop through all material points and calculate failure indicies

           ! Fiber fracture Tensile/Compressive
           
        A = StressNew(I,1)-(XNU12-XNU12f*MGF*(E11/E11f))*(StressNew(I,2)+StressNew(I,3))
                
              IF(A.GE. 0.0) THEN        
                FF = (1.0/Xt)*A
                e(I,1) = FF
                 
              ELSEIF(A.LT.0.0) THEN 
                 
                FF = (1.0/Xc)*A
                e(I,1) = FF
              
              END IF
              
            ! Inter-Fiber fracture by finidng fracture plane angle
                            
            DO J= 0,180 !loop from -90 to 90 degress
            
                THETA = -Pi/2.0 + J*(Pi/180)
            	
				StressNew(I,2)= 37.5
                StressNew(I,3)= 0.0
                StressNew(I,5)= 0.0
                StressNew(I,4)= 66.5  
                StressNew(I,6)= 0.0
            
				
                SFP = (StressNew(I,2)*cos(THETA)**2)&
                &+(StressNew(I,3)*sin(THETA)**2)+(2*StressNew(I,6)*sin(THETA)*cos(THETA))
            
                TNT = -StressNew(I,2)*sin(THETA)*cos(THETA)+StressNew(I,3)*sin(THETA)*cos(THETA)+&
                &StressNew(I,6)*(cos(THETA)**2-sin(THETA)**2)
                            
                TN1 = StressNew(I,4)*cos(THETA) + StressNew(I,5)*sin(THETA)
                   
                Rn= Yt
            
                Rn1 = S12
            
                Rnt = Yc/(2.0*tan(THETAF))
            
                Pnt = -1.0/(2*tan(2*THETAF))
            
                Pn1 = Pnt*(Rn1/Rnt)
            
                IF (SFP.GE.0.0) THEN
            
                    IFF = (SFP/Rn)**2 +(TN1/(Rn1-Pn1*SFP))**2 + (TNT/(Rnt-Pnt*SFP))**2
                          open (unit=1, file='puck.dat')
                          write(1,*) SFP,TN1
            
                ELSE IF (SFP.LT.0.0) THEN
            
                    IFF = (TN1/(Rn1-Pn1*SFP))**2 +(TNT/(Rnt-Pnt*SFP))**2
            
                end if
                
                 ! Store global maximum failure index and facture plane angle
                                        
                if(IFF.GT.FMAX) then
            
                    FMAX =IFF
                    THETAMAX=THETA
     
                END IF 
            
                END DO

               close(unit=1)
            
                    e(I,2) = FMAX
            
           
           END DO
         
!3.6 (Damage Application based on failure mode)
!==================================================================

       DO I= 1,nblock ! Loop through all material points
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Fiber Damage Modes

           IF((e(I,1).GT.1.0) .AND. (A.GE.0.0)) THEN  ! Check for fiber tensile failure

                    fflags(I,1)=STEP    ! Assign Failure Flag to step number as a positive number

                    Dmg(I,1)= (0.999999-XbetaT)  ! Apply fiber tensile damage 
  

           ELSE IF((e(I,1).GT.1.0) .AND. (A.LT.0.0)) THEN  ! Check for fiber compressive failure

                    fflags(I,1)=-STEP    ! Assign Failure Flag to step number as a negarive number

                    Dmg(I,1)= (0.999999-XbetaC)  ! Apply fiber compressive damage
           END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           

           ! Matrix Damage Modes

           IF((e(I,2).GT.1.0) .AND. (SFP.GE.0.0)) THEN  ! Check for matrix tensile failure

                    fflags(I,2)=STEP    ! Assign Failure Flag to step number as a positive number

                    Dmg(I,2)= (0.999999-XbetaT)  ! Apply matrix tensile damage 
                    Dmg(I,3)= (0.999999-XbetaT)
                    Dmg(I,4)= (0.999999-XbetaS)  ! Apply shear damage (1-2)
                    Dmg(I,6)= (0.999999-XbetaS)  ! Apply shear damage (2-3)

           ELSE IF((e(I,2).GT.1.0) .AND. (SFP.LT.0.0)) THEN  ! Check for matrix compressive failure

                    fflags(I,2)=-STEP    ! Assign Failure Flag to step number as a negarive number

                    Dmg(I,2)= (0.999999-XbetaC)  ! Apply matrix compressive damage
                    Dmg(I,3)= (0.999999-XbetaC)
                    Dmg(I,4)= (0.999999-XbetaS)  ! Apply shear damage (1-2)
                    Dmg(I,6)= (0.999999-XbetaS)  ! Apply shear damage (2-3)                    
           END IF                                              

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

                END DO

                    END IF


!3.9 (Output Results)
!================================================================== 

!$$$$$$  Print *, 'Step No:', STEP  ! Analysis step number
!$$$$$$  Print *, '****************************************************************'
!$$$$$$  
!$$$$$$  WRITE(*,*)  ! print blank line
!$$$$$$ 
!$$$$$$  print '(" Failure ratio in normal 2-direction: ",f6.3)',e(1,2) ! print stress/strength ratio (2) for material point number 1
!$$$$$$  Print *, '________________________________________________________________'
!$$$$$$  
!$$$$$$  Print *, 'Damage Factor in normal 2-direction: ', Dmg(1,2) ! print damage (2) factor for material point number 1
!$$$$$$  Print *, '________________________________________________________________'
!$$$$$$  
!$$$$$$  print '(" Stress in normal 2-direction: ",f6.3)',StressNew(1,2) ! print Stress (2) for material point number 1
!$$$$$$  Print *, '________________________________________________________________'
!$$$$$$  
!$$$$$$  Print *, 'Strain in normal 2-direction: ', Strain(1,2)     ! print Stain (2) for material point number 1
!$$$$$$  Print *, '________________________________________________________________'
!$$$$$$  
!$$$$$$  Print *, 'Strain Increment in normal 2-direction: ', StrainInc(1,2)  ! print Stain increment (2) for material point number 1
!$$$$$$  Print *, '________________________________________________________________'
!$$$$$$ 
!$$$$$$  Print *, 'Fracture angle: ', THETAMAX  ! print Stain increment (2) for material point number 1
!$$$$$$  Print *, '________________________________________________________________'
!$$$$$$  
!$$$$$$ 
!$$$$$$  WRITE(*,*)
!$$$$$$  WRITE(*,*)


                        END DO  ! END ANALYSIS STEPS COUNTER (END OF ANALYSIS)

                        
               
 END PROGRAM UserMaterialTest
  
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Utility Subroutine: Ortho3D 

! The following subroutine : 1. Build the effective constiutive matrix based on orthtoropic elasticity

!                            2. Return the updated stresses

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
