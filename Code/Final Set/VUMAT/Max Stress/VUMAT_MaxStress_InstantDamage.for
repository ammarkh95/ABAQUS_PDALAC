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

  
 	  INTEGER fflags(6)
	  
	  DOUBLE PRECISION Dgrd(3),dmg(6),e(6),DDSDDE(6,6)

!DIR$ FREEFORM		   

!C=========================================================================		  

!C ELASTIC PROPERTIES

    E11=PROPS(1)
    E22=PROPS(2)
    E33=PROPS(3)
    XNU12=PROPS(4)
    XNU13=PROPS(5)
    XNU23=PROPS(6)   
    G12=PROPS(7)
    G13=PROPS(8)
    G23=PROPS(9)
	  
!C=========================================================================		  

!C Failure Stresses (Tension,Compression,Shear)

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
      
!C=========================================================================	
	  
!C Degredation Parameters 

! Read Damage Factors Array
! Tension, Compression, Shear

    Dgrd(1:3)=PROPS(19:21)     
          	  	  
!C Process element block
!C=========================================================================  

        if(stepTime.EQ.0) then

 
 ! Purely Elastic material at begining of analysis 
 
            dmg = 1.0d+0
            fflags=0

            do I = 1,nblock   ! Loop through all material points     

!C Calculate Stresses
!C==================================================================

        XNU21=(E22*XNU12)/E11
        XNU31=(E33*XNU13)/E11
        XNU32=(E33*XNU23)/E22
        Upsilon=(1.0d+0)/(1.0d+0-XNU12*XNU21-XNU23*XNU32-XNU13*XNU31-2.0d+0*XNU21*XNU32*XNU13)        
        		
        ! Form effective orthotropic material matrix including damage (Jacobian Matrix) 

        DDSDDE(1,1)=(E11*Upsilon*(1.0-XNU23*XNU32))*abs(dmg(1))
        DDSDDE(1,2)=(E11*Upsilon*(XNU21+XNU31*XNU23))*abs(dmg(1)*dmg(2))
        DDSDDE(1,3)=(E11*Upsilon*(XNU31+XNU21*XNU32))*abs(dmg(1)*dmg(3))
        DDSDDE(2,1)=(E11*Upsilon*(XNU21+XNU31*XNU23))*abs(dmg(1)*dmg(2))
        DDSDDE(2,2)=(E22*Upsilon*(1.0-XNU13*XNU31))*abs(dmg(2))
        DDSDDE(2,3)=(E22*Upsilon*(XNU32+XNU12*XNU31))*abs(dmg(3)*dmg(2))
        DDSDDE(3,1)=(E11*Upsilon*(XNU31+XNU21*XNU32))*abs(dmg(3)*dmg(1))
        DDSDDE(3,2)=(E22*Upsilon*(XNU32+XNU12*XNU31))*abs(dmg(3)*dmg(2))
        DDSDDE(3,3)=(E33*Upsilon*(1.0-XNU12*XNU21))*abs(dmg(3))
        DDSDDE(4,4)=2.0d+0*G12*abs(dmg(4))
        DDSDDE(5,5)=2.0d+0*G23*abs(dmg(6))
        DDSDDE(6,6)=2.0d+0*G13*abs(dmg(5))

			
	    ! Calculate Stress Components
            
            stressNew(I,1)= stressOld(I,1)+(DDSDDE(1,1)*strainInc(I,1)+DDSDDE(1,2)*strainInc(I,2)+DDSDDE(1,3)*strainInc(I,3))
            stressNew(I,2)= stressOld(I,2)+(DDSDDE(2,1)*strainInc(I,1)+DDSDDE(2,2)*strainInc(I,2)+DDSDDE(2,3)*strainInc(I,3))
            stressNew(I,3)= stressOld(I,3)+(DDSDDE(3,1)*strainInc(I,1)+DDSDDE(3,2)*strainInc(I,2)+DDSDDE(3,3)*strainInc(I,3)) 
            stressNew(I,4)= stressOld(I,4)+(DDSDDE(4,4)*strainInc(I,4))
            stressNew(I,5)= stressOld(I,5)+(DDSDDE(5,5)*strainInc(I,5))   	  
            stressNew(I,6)= stressOld(I,6)+(DDSDDE(6,6)*strainInc(I,6))

        	stateNew(I,1)= 1.0d+0-dmg(1)
        	stateNew(I,2)= 1.0d+0-dmg(2)            
        	stateNew(I,3)= 1.0d+0-dmg(3)             
        	stateNew(I,4)= 1.0d+0-dmg(4) 
        	stateNew(I,5)= 1.0d+0-dmg(5) 
        	stateNew(I,6)= 1.0d+0-dmg(6) 
            
        ! Assign updated state variables related to failure flags 
			
            stateNew(I,7)= fflags(1)
            stateNew(I,8)= fflags(2)
            stateNew(I,9)= fflags(3)
            stateNew(I,10)= fflags(4)                       
            stateNew(I,11)= fflags(5)
            stateNew(I,12)= fflags(6)
            
           end do  
              
     else
     
     do i=1,nblock
        
        dmg(1)=1.0d+0-stateOld(I,1)
        dmg(2)=1.0d+0-stateOld(I,2)
        dmg(3)=1.0d+0-stateOld(I,3)        
        dmg(4)=1.0d+0-stateOld(I,4)
        dmg(5)=1.0d+0-stateOld(I,5)
        dmg(6)=1.0d+0-stateOld(I,6)
        
        fflags(1)=stateOld(I,7)
        fflags(2)=stateOld(I,8)
        fflags(3)=stateOld(I,9) 
        fflags(4)=stateOld(I,10)                       
        fflags(5)=stateOld(I,11)
        fflags(6)=stateOld(I,12)        

                   
!C Update Stresses
!C=========================================================================

        XNU21=(E22*XNU12)/E11
        XNU31=(E33*XNU13)/E11
        XNU32=(E33*XNU23)/E22
        Upsilon=(1.0d+0)/(1.0d+0-XNU12*XNU21-XNU23*XNU32-XNU13*XNU31-2.0d+0*XNU21*XNU32*XNU13)        
        		
        ! Form effective orthotropic material matrix including damage (Jacobian Matrix) 

        DDSDDE(1,1)=(E11*Upsilon*(1.0-XNU23*XNU32))*abs(dmg(1))
        DDSDDE(1,2)=(E11*Upsilon*(XNU21+XNU31*XNU23))*abs(dmg(1)*dmg(2))
        DDSDDE(1,3)=(E11*Upsilon*(XNU31+XNU21*XNU32))*abs(dmg(1)*dmg(3))
        DDSDDE(2,1)=(E11*Upsilon*(XNU21+XNU31*XNU23))*abs(dmg(1)*dmg(2))
        DDSDDE(2,2)=(E22*Upsilon*(1.0-XNU13*XNU31))*abs(dmg(2))
        DDSDDE(2,3)=(E22*Upsilon*(XNU32+XNU12*XNU31))*abs(dmg(3)*dmg(2))
        DDSDDE(3,1)=(E11*Upsilon*(XNU31+XNU21*XNU32))*abs(dmg(3)*dmg(1))
        DDSDDE(3,2)=(E22*Upsilon*(XNU32+XNU12*XNU31))*abs(dmg(3)*dmg(2))
        DDSDDE(3,3)=(E33*Upsilon*(1.0-XNU12*XNU21))*abs(dmg(3))
        DDSDDE(4,4)=2.0d+0*G12*abs(dmg(4))
        DDSDDE(5,5)=2.0d+0*G23*abs(dmg(6))
        DDSDDE(6,6)=2.0d+0*G13*abs(dmg(5))


				
            stressNew(I,1)= stressOld(I,1)+(DDSDDE(1,1)*strainInc(I,1)+DDSDDE(1,2)*strainInc(I,2)+DDSDDE(1,3)*strainInc(I,3))
            stressNew(I,2)= stressOld(I,2)+(DDSDDE(2,1)*strainInc(I,1)+DDSDDE(2,2)*strainInc(I,2)+DDSDDE(2,3)*strainInc(I,3))
            stressNew(I,3)= stressOld(I,3)+(DDSDDE(3,1)*strainInc(I,1)+DDSDDE(3,2)*strainInc(I,2)+DDSDDE(3,3)*strainInc(I,3)) 
            stressNew(I,4)= stressOld(I,4)+(DDSDDE(4,4)*strainInc(I,4))
            stressNew(I,5)= stressOld(I,5)+(DDSDDE(5,5)*strainInc(I,5))   	  
            stressNew(I,6)= stressOld(I,6)+(DDSDDE(6,6)*strainInc(I,6))   

! !C Failure Evaluation as per Maximum Stress Criteria
! !C=========================================================================

		 ! Failure Evaluation as per Maximum Stress Criteria

			  ! Failure along 1-direction (Tension/Compression)
				
              IF (stressNew(I,1).GE. 0.0) THEN
                
            	 e(1)=stressNew(I,1)/Xt
                 
        	  ELSE
               
            	e(1)=-stressNew(I,1)/Xc
                
        	  END IF

            ! Failure along 2-direction  (Tension/Compression)
	
        	 IF (stressNew(I,2).GE. 0.0) THEN 
            
                e(2)=stressNew(I,2)/Yt
				
        	 ELSE
             
            	e(2)=-stressNew(I,2)/Yc 
                
        	 END IF

            ! Failure along 3-direction  (Tension/Compression)       
		
        	 IF (stressNew(I,3).GE. 0.0) THEN 
            
            	e(3)=stressNew(I,3)/Zt
                
        	 ELSE
            
                e(3)=-stressNew(I,3)/Zc
            
        	 END IF

            ! Shear Failure along 1-2 direction

            	e(4)=abs(stressNew(I,4))/S12              
	
            ! Shear Failure along 2-3 direction
             
            	e(5)=abs(stressNew(I,5))/S13
                
            ! Shear Failure along 1-3 direction

           	    e(6)=abs(stressNew(I,6))/S23
                   
! !C Damage Application based on Failure Mode
! !C=========================================================================

        	DO J=1,3   ! Loop through normal stress components 

        		IF((e(J).GT.1.0) .AND. (stressNew(I,J).GE.0.0)) THEN  ! Check for tensile failure

    				fflags(J)=stepTime    ! Assign Failure Flag to step number as a positive number

					dmg(J)= (0.999999-Dgrd(1))  ! Apply tensile damage 


    			ELSE IF((e(J).GT.1.0) .AND. (stressNew(I,J).LT.0.0)) THEN  ! Check for compressive failure or previous damage

    				fflags(J)=-stepTime    ! Assign Failure Flag to step number as a positive number

					dmg(J)= (0.999999-Dgrd(2))  ! Apply compressive damage


                   END IF

					END DO

        	DO J=4,6   ! Loop through shear stress components 
            

				IF((e(J).GT.1.0)) THEN  ! Check for shear failure or previos damage


    				fflags(J)=stepTime    ! Assign Failure Flag to step number as a positive number

					dmg(J)= (0.999999-Dgrd(3))  ! Apply shear damage

                   END IF

                     END DO                  	  
		
!C Calculate Stresses
!C=========================================================================		
	   
        XNU21=(E22*XNU12)/E11
        XNU31=(E33*XNU13)/E11
        XNU32=(E33*XNU23)/E22
        Upsilon=(1.0d+0)/(1.0d+0-XNU12*XNU21-XNU23*XNU32-XNU13*XNU31-2.0d+0*XNU21*XNU32*XNU13)        
        		
        ! Form effective orthotropic material matrix including damage (Jacobian Matrix) 

        DDSDDE(1,1)=(E11*Upsilon*(1.0-XNU23*XNU32))*abs(dmg(1))
        DDSDDE(1,2)=(E11*Upsilon*(XNU21+XNU31*XNU23))*abs(dmg(1)*dmg(2))
        DDSDDE(1,3)=(E11*Upsilon*(XNU31+XNU21*XNU32))*abs(dmg(1)*dmg(3))
        DDSDDE(2,1)=(E11*Upsilon*(XNU21+XNU31*XNU23))*abs(dmg(1)*dmg(2))
        DDSDDE(2,2)=(E22*Upsilon*(1.0-XNU13*XNU31))*abs(dmg(2))
        DDSDDE(2,3)=(E22*Upsilon*(XNU32+XNU12*XNU31))*abs(dmg(3)*dmg(2))
        DDSDDE(3,1)=(E11*Upsilon*(XNU31+XNU21*XNU32))*abs(dmg(3)*dmg(1))
        DDSDDE(3,2)=(E22*Upsilon*(XNU32+XNU12*XNU31))*abs(dmg(3)*dmg(2))
        DDSDDE(3,3)=(E33*Upsilon*(1.0-XNU12*XNU21))*abs(dmg(3))
        DDSDDE(4,4)=2.0d+0*G12*abs(dmg(4))
        DDSDDE(5,5)=2.0d+0*G23*abs(dmg(6))
        DDSDDE(6,6)=2.0d+0*G13*abs(dmg(5))

			
	    ! Calculate Stress Components
            
            stressNew(I,1)= stressOld(I,1)+(DDSDDE(1,1)*strainInc(I,1)+DDSDDE(1,2)*strainInc(I,2)+DDSDDE(1,3)*strainInc(I,3))
            stressNew(I,2)= stressOld(I,2)+(DDSDDE(2,1)*strainInc(I,1)+DDSDDE(2,2)*strainInc(I,2)+DDSDDE(2,3)*strainInc(I,3))
            stressNew(I,3)= stressOld(I,3)+(DDSDDE(3,1)*strainInc(I,1)+DDSDDE(3,2)*strainInc(I,2)+DDSDDE(3,3)*strainInc(I,3)) 
            stressNew(I,4)= stressOld(I,4)+(DDSDDE(4,4)*strainInc(I,4))
            stressNew(I,5)= stressOld(I,5)+(DDSDDE(5,5)*strainInc(I,5))   	  
            stressNew(I,6)= stressOld(I,6)+(DDSDDE(6,6)*strainInc(I,6))  
		
			
!C Update State Variables
!C=========================================================================	
        		
        	stateNew(I,1)= 1.0d+0-dmg(1)
        	stateNew(I,2)= 1.0d+0-dmg(2)            
        	stateNew(I,3)= 1.0d+0-dmg(3)             
        	stateNew(I,4)= 1.0d+0-dmg(4) 
        	stateNew(I,5)= 1.0d+0-dmg(5) 
        	stateNew(I,6)= 1.0d+0-dmg(6)       
            
        ! Assign updated state variables related to failure flags 
			
            stateNew(I,7)= fflags(1)
            stateNew(I,8)= fflags(2)
            stateNew(I,9)= fflags(3)
            stateNew(I,10)= fflags(4)                       
            stateNew(I,11)= fflags(5)
            stateNew(I,12)= fflags(6)
            
           END DO
           
           end if

		RETURN
        END

