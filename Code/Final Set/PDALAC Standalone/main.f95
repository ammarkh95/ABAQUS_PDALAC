! This is main execution file for the program

    program main

    use material_parameters

    use analysis_parameters

    use PDA

    implicit none

    integer :: mat_id, analysis_id

    print*,'********************************************************************'
    print"(a32)","Technical University of Munich"                                                                                                               
    print"(/,a34)","Chair of Computational Mecahnics"                                                                                                                                          
    print"(/,a10,a47,/)","SL 2019:","Development of failure criteria for composites"
    print *,'********************************************************************'
    print"(a60)"," _____  _____          _               _____ "
    print"(a60)","|  __ \|  __ \   /\   | |        /\   / ____|"
    print"(a60)","| |__) | |  | | /  \  | |       /  \ | |     "
    print"(a60)","|  ___/| |  | |/ /\ \ | |      / /\ \| |     "
    print"(a60)","| |    | |__| / ____ \| |____ / ____ \ |____ "
    print"(a60)","|_|    |_____/_/    \_\______/_/    \_\_____|"
    print"(/,a65,/)","Progressive damage analysis of laminated composites"
    print"(/,a19,/)","V 1.0"
    print *,'********************************************************************'


!###########              User Input              ###########
!#----------------------------------------------------------#

! Select material card number from predfined list in parameters file
    mat_id = 2
    
! Select analysis option number from predfined list in parameters file
    analysis_id =2

!###########    Load Material Analysis Model      ###########
!#----------------------------------------------------------#

    print"(/,a28)","Step 1: Loading Model Data"
    print *, '___________________________________'

! Read material porperties

    call mat_props(mat_id)

    print"(/,a27,i2,/)","Loaded Material Card No: ", mat_id

!------------------------------------------------------------------------

! Read analysis settings

    call analysis_set(analysis_id)

    print"(a29,i2,/)","Loaded Analysis Option No: ", analysis_id   

! Confirm successful loading of data

    print"(a32,/,/)","Model data read successfully.."
    print *, '___________________________________'

!###########         Start Analysis Phase         ###########
!#----------------------------------------------------------#

    print"(/,a32)","Step 2: Analysis Phase Summary"
    print*, '___________________________________'
    print"(/,a20,a13)","Failure Criteria: ",failure_criteria(failure_id)
    print"(/,2a16)","Damage Model: ",damage_criteria(damage_id)
    print"(/,a16,i2)","No. of Steps: ",nsteps
    print"(/,a27,/)","Applied Strain Increments"
    print"(6a11)","Eps11","Eps22","Eps33","Gamm12","Gamm13","Gam23"
    print"(6e11.3)",StrainInc(1,:)

! Call the progressive damage analysis solver to start the analysis

    call PDA_Solver 

! Confirm Successful completion of analysis

    print"(/,/,a36)","Analysis Completed successfully.."
    print *, '___________________________________'



!###########          Post-Processing             ###########
!#----------------------------------------------------------#

! GNU Plot is used for plotting stress,strain graphs

! Refer to data_plot.plt for customizing data output

    call system('gnuplot -p data_plot.plt')
                                    
    end program main
