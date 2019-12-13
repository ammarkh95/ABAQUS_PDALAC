! This is the main input interface for the program

! Here you can define: Material propoerties 
!                      Analysis variables
!                      Output results options   

    module material_parameters

    implicit none

! Declare 3D Material Elastic Properties 
    real*8 :: E11, E22, E33, G12, G13, G23, ANU12, ANU13, ANU23, Ypsilon, ANU21, ANU31, ANU32

! Declare 3D Material Stress Strength Parameters
    real*8 :: Xt, Xc, Yt, Yc, Zt, Zc, S12, S13, S23

! Declare 3D Material Strain Strength Parameters
    real*8 :: EpsXt, EpsXc, EpsYt, EpsYc, EpsZt, EpsZc, GamS12, GamS13, GamS23

! In this section you can define multiple sets of composite laminate properties 

! referred to as  "Material Cards", each material card has a unique "mat_id"

! material cards can be added by creating additional cases in the selection list

! The following legend provides description for the material input variables

!=================================================================================
!                        ****Orthotropic Material Constants****
! +-----------+-------------------------------------------------------------------+
! | Parameter |                            Description                            |
! +-----------+-------------------------------------------------------------------+
! | E11       | Longitudinal elastic Young's modulus(Fiber: 11 axis)              |
! | E22       | Transverse elastic Young's modulus   (Matrix: 22 axis)            |
! | E33       | Interlaminar elastic Young's modulus (Interlaminar: 33 axis)      |
! | G12       | Elastic Shear Modulus  (1-2 plane)                                |
! | G13       | Elastic Shear Modulus  (1-3 plane)                                |
! | G23       | Elastic Shear Modulus  (2-3 plane)                                |
! | ANU12     | Poisson ratio   (1-2 plane)                                       |
! | ANU13     | Poisson ratio   (1-3 plane)                                       |
! | ANU23     | Poisson ration  (2-3 plane)                                       |
! | Xt        | Laminate tensile strength (11 fiber direction)                    |
! | Xc        | Laminate compressive strength  (11 fiber direction)               |
! | Yt        | Laminate tensile strength  (22 matrix direction)                  |
! | Yc        | Laminate compressive strength (22 matrix direction)               |
! | Zt        | Laminate interlaminar tensile strength (33 direction)             |
! | Zc        | Laminate interlaminar compressive strength (33 direction)         |
! | S12       | Laminate Shear Strength  (1-2 plane)                              |
! | S13       | Laminate Shear Strength  (1-3 plane)                              |
! | S23       | Laminate Shear Strength  (2-3 plane)                              |
! | EpsXt     | Laminate allowable tensile strain (11 fiber direction)            |
! | EpsXc     | Laminate allowable compressive strain  (11 fiber direction)       |
! | EpsYt     | Laminate allowable tensile strain (Matrix: 22 axis)               |
! | EpsYc     | Laminate allowable compressive strain (Matrix: 22 axis)           |
! | EpsZt     | Laminate allowable interlaminar tensile strain (33 direction)     |
! | EpsZc     | Laminate allowable interlaminar compressive strain (33 direction) |
! | GamS12    | Laminate allowable shear Strain  (1-2 plane)                      |
! | GamS13    | Laminate allowable shear Strain  (1-3 plane)                      |
! | GamS23    | Laminate allowable shear Strain  (2-3 plane)                      |
! +-----------+-------------------------------------------------------------------+

! Important Remarks: 

! *Engineering shear strain convention is used 

! *In FEA, the material local axis need to be assigned with repspect to the   
!  above convention in order to establish the correct material properties orienatation

! *The elastic constants need to satisfy the stability limits for an orthotropic 
!  material law as detailed in the project documentation

!=================================================================================

    contains

    subroutine mat_props(id)

    implicit none 

    integer, intent(in):: id

    select case (id)  

!###########         Material cards List          ###########
!#----------------------------------------------------------#

    case(1) ! Material Card No.1 
    
    !T300/5208 Graphite/epoxy
    !laminate stacking sequence [(±45)2/90/(±45)/̅0]s.

    E11=117.212d+0 
    E22=9.653d+0   
    E33=9.653d+0   
    ANU12=0.33d+0
    ANU13=0.33d+0
    ANU23=0.34d+0 
    G12=3.103d+0   
    G13=3.103d+0   
    G23=3.103d+0 
    
    ! Complementary elastic properties 
       
    ANU21=(E22*ANU12)/E11
    ANU31=(E33*ANU13)/E11
    ANU32=(E33*ANU23)/E22
    Ypsilon=(1.0d+0)/(1.0d+0-ANU12*ANU21-ANU23*ANU32-ANU13*ANU31-2.0d+0*ANU21*ANU32*ANU13)

    
    ! Allowable Stresses 
    
    Xt=1.378d+0  
    Xc=0.689d+0 
    Yt=0.103d+0 
    Yc=0.206d+0  
    Zt=0.103d+0  
    Zc=0.206d+0 
    S12=0.103d+0 
    S13=0.103d+0 
    S23=0.103d+0  


    case(2) ! Material Card No.2 
    
    ! CFRP IM7/977-3 Laminate [45/0/90/-45]2s 

    E11=164.0d+0
    E22=8.98d+0   
    E33=8.98d+0   
    ANU12=0.32d+0 
    ANU13=0.32d+0 
    ANU23=0.496d+0 
    G12=5.02d+0   
    G13=5.02d+0   
    G23=3.00d+0
    
    ! Complementary elastic properties
    
    ANU21=(E22*ANU12)/E11
    ANU31=(E33*ANU13)/E11
    ANU32=(E33*ANU23)/E22
    Ypsilon=(1.0d+0)/(1.0d+0-ANU12*ANU21-ANU23*ANU32-ANU13*ANU31-2.0d+0*ANU21*ANU32*ANU13)
    
    ! Allowable Stresses 
    
    Xt=2.9d+0  
    Xc=1.68d+0 
    Yt=0.10d+0 
    Yc=0.247d+0  
    Zt=0.10d+0 
    Zc=0.247d+0 
    S12=0.08d+0 
    S13=0.08d+0 
    S23=0.08d+0
    
    ! Allowable Strains
    
    
    EpsXt=Xt/E11
    
    EpsXc=Xc/E11
    
    EpsYt=Yt/E22
    
    EpsYc=Yc/E22
    
    EpsZt=Zt/E33
    
    EpsZc=Zc/E33
    
    GamS12=S12/G12
    
    GamS13=S13/G13

    GamS23=S23/G23

    case default

! Print error message and exit the program for invalid selection

    print*, "Invalid Material ID"

    call exit (0)

    end select

    end subroutine mat_props

    end module material_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module analysis_parameters

    implicit none

    integer :: I, nsteps, damage_id, failure_id
    integer, parameter :: nblock =5

! Declare Progressive Damage analysis parameters
    real*8 :: beta_ft, beta_fc, beta_mt, beta_mc, beta_s
    real*8 :: a_ft, a_fc, a_mt, a_mc, a_s, n_ft, n_fc, n_mt, n_mc, n_s     
    real*8 :: G_ft, G_fc, G_mt, G_mc, G_IIC, le, let, alpha, beta
    real*8 :: THETAF, MGF, ANU12f, E11f    

! Declare Strain increment array  
    real*8, dimension  (nblock,6) :: StrainInc
    
! Declare array of characters for failure ciriteria and damage models names

    character(LEN=16) :: failure_criteria(7),damage_criteria(5)
    
    data failure_criteria /'Max Stress','Max Strain','Tsai-Wu','Hoffman','Hashin','Hashin-Rotem','Puck'/
    
    data damage_criteria /'Instant','Recursive','Exponential','Constant Stress','CDM'/    

!=================================================================================
!               ****Progressive Analysis Variables****
! +-------------+-----------------------------------------------------+
! |  Parameter  |                     Description                     |
! +-------------+-----------------------------------------------------+
! | analysis_id | Index of analysis option set to be used             |
! | nblock      | Number of material integration points               |
! | nsteps      | Number of load steps                                |
! | StrainInc   | Strain increments vector                            |
! | THETAF      | Maximum fracture angle for Puck in radians          |
! | MGF         | Magnification factor for Puck                       |
! | E11F        | Fiber elastic modulus for Puck                      |
! | ANU12F      | Fiber poisson ratio for Puck                        |
! | failure_id  | Index for failure criteria selection:               |
! |             | --->1-Max Stress                                    |
! |             | --->2-Max Strain                                    |
! |             | --->3-Tsai-Wu                                       |
! |             | --->4-Hoffman                                       |
! |             | --->5-Hashin                                        |
! |             | --->6-Hashin-Rotem                                  |
! |             | --->7-Puck                                          |
! | damage_id   | Index for damage model selection:                   |
! |             | --->1-Ply Discount:Instantaneous Degradation        |
! |             | --->2-Ply Discount:Recursive Degradation            |
! |             | --->3-Ply Discount:Exponential Degradation          |
! |             | --->4-Ply Discount:Constant Stress Degradation      |
! |             | --->5-Continium damage mechanics: Crack Band Theory |
! +-------------+-----------------------------------------------------+

!           ****Ply Discount Damage Parameters****
! +------------+----------------------------------------------+
! | Parameter  |               Description                    |
! +------------+----------------------------------------------+
! +-----------------------------------------------------------+
! |            Instant/Recursive Damage                       |
! +-----------------------------------------------------------+
! |  beta_ft   | degradation factor (fiber tension)           |
! |  beta_fc   | degradation factor (fiber compression)       |
! |  beta_mt   | degradation factor (matrix tension)          |
! |  beta_mc   | degradation factor (matrix compression)      |
! |  beta_s    | degradation factor (shear)                   |
! +-----------------------------------------------------------+
! |            Exponential Damage                             |
! +-----------------------------------------------------------+
! | a_ft, n_ft | Exp. degradation factors (fiber tension)     |
! | a_fc, n_fc | Exp. degradation factor (fiber compression)  |
! | a_mt, n_mt | Exp. degradation factor (matrix tension)     |
! | a_mc, n_mc | Exp. degradation factor (matrix compression) |
! | a_s,  n_s  | Exp. degradation factor (shear)              |
! +------------+----------------------------------------------+


!       ****Continium Damage Mechanics Parameters****
! +-----------+-----------------------------------------+
! | Parameter |               Description               |
! +-----------+-----------------------------------------+
! | G_ft      | fracture energy (fiber tension)         |
! | G_fc      | fracture energy (fiber compression)     |
! | G_mt      | fracture energy (matrix tension)        |
! | G_mc      | fracture energy (matrix compression)    |
! | G_IIC     | fracture energy (Mode II shear failure) |
! | le        | characteristic element length           |
! | let       | characteristic element thickness        |
! | alpha,beta| Nonlinear Shear degradation factors     |
! +-----------+-----------------------------------------+

! Important Remarks: 

!*Degradation factors values are in the range of (0,1) 
!*Degredation factor == 0 implies no damage
!*Degredation factor == 1 implies full damage (i.e stiffness approx zero)

!=================================================================================

    contains

    subroutine analysis_set (id)

    implicit none 

    integer, intent(in):: id

    select case (id)

!###########        Analysis Options sets         ###########
!#----------------------------------------------------------#

    case(1) ! Analaysis options set.1 

    nsteps = 35

    do I = 1,nblock      
    StrainInc(I,1)= 0.0d+0
    StrainInc(I,2)= -0.001d+0 ! (Uni-axial compression)
    StrainInc(I,3)= 0.0d+0
    StrainInc(I,4)= 0.0d+0
    StrainInc(I,5)= 0.0d+0
    StrainInc(I,6)= 0.0d+0        
    end do

    failure_id = 7 !(Puck)
    damage_id =  2 !(Recursive damage)
    
    
    beta_ft = 0.8d+0
    beta_fc = 0.8d+0
    beta_mt = 0.5d+0
    beta_mc = 0.7d+0
    beta_s  = 0.5d+0
       
    ! Puck related parameters   
      
    E11f= 117.212d+0
    MGF= 1.2d+0
    THETAF= 0.89d+0
    ANU12f = 0.2d+0
    

    case(2) ! Analaysis options set.2 


    nsteps = 25

    do I = 1,nblock      
    StrainInc(I,1)= 0.001d+0 !(Uni-axial tension)
    StrainInc(I,2)= 0.0
    StrainInc(I,3)= 0.0
    StrainInc(I,4)= 0.0
    StrainInc(I,5)= 0.0
    StrainInc(I,6)= 0.0        
    end do

    failure_id = 6 !(Hashin-Rotem)
    damage_id =  5 !(Continium Damage Mechanics)
    
    ! Continium damage mechanics parameters
    
    G_ft=81.53/1000D+0
    G_fc=24.53/1000D+0
    G_mt=0.256/1000D+0
    G_mc=1.156/1000D+0
    G_IIC=1.156/1000D+0
    le=4d+0
    let=0.127d+0
    
    alpha=0.5d+0
    beta=0.5d+0
    
    

    case default 

! Print error message and exit the program for invalid selection

    print*, "Invalid Analysis ID"

    call exit (0)

    end select

    end subroutine analysis_set 
    
    end module analysis_parameters
