!--------------------------------------------------------------------------------------------------
! PROGRAM: cis_nto_prog
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date April, 2019
!
!> @brief Analyze CIS type wave functions in terms of their NTOs.
!--------------------------------------------------------------------------------------------------
program nto_prog
    ! General
    use global_defs
    ! Variables
    use nto_variables
    ! I/O
    use nto_input_mod
    use nto_output_mod
    ! Calculation
    use nto_cis_mod
    implicit none

    ! All program variables are defined in overlap_variables_mod


    ! Begin run.
    time00 = omp_get_wtime()
    call set_defaults()
    call command_line_interface()
    if (print_level >= 1) then
        write(stdout, '(5x,a)') '-------------------------------------------------------------'
        write(stdout, '(5x,a)') '                         nto program                         '
        write(stdout, '(5x,a)') '                         version 1.1                         '
        write(stdout, '(5x,a)') '                                                             '
        write(stdout, '(5x,a)') ' Program compiled on '//__DATE__//' '//__TIME__//'.          '
        write(stdout, '(5x,a)') '-------------------------------------------------------------'
        write(stdout, *) 
        write(stdout, '(1x,a,i0,a)') 'Using ', omp_get_max_threads(), ' threads.'
        write(stdout, *) 
        write(stdout, '(1x,a)') 'Reading orbitals/states: '
        write(stdout, '(5x,a,a)') 'Format: ', input_format
        write(stdout, '(5x,a,a)') 'Path: ', path
    end if

    call nto_cis()

    call output_nto()


end program nto_prog
