!--------------------------------------------------------------------------------------------------
! PROGRAM: cis_overlap_prog
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @brief Calculate overlap between two sets of CIS type wave functions.
!--------------------------------------------------------------------------------------------------
program cis_overlap_prog
    ! General
    use global_defs
    ! I/O
    use overlap_output_mod
    use overlap_input_mod
    ! Calculation
    use ao_overlap_mod
    use mo_overlap_mod
    use cis_overlap_mod
    ! Output
    use overlap_output_mod
    implicit none

    ! Results
    real(dp), allocatable :: s_ao(:, :) !< Atomic orbital overlaps.
    real(dp), allocatable :: s_mo_a(:, :) !< Molecular orbital overlaps alpha.
    real(dp), allocatable :: s_mo_b(:, :) !< Molecular orbital overlaps beta.
    real(dp), allocatable :: s_wf(:, :) !< Wave function overlaps.
    ! Option and timing variables are defined in overlap_input_mod


    ! Begin run.
    time00 = omp_get_wtime()
    call set_defaults()
    call command_line_interface()
    if (print_level >= 1) then
        write(stdout, '(5x,a)') '-------------------------------------------------------------'
        write(stdout, '(5x,a)') '                     cis_overlap program                     '
        write(stdout, '(5x,a)') '                         version 1.1                         '
        write(stdout, '(5x,a)') '                                                             '
        write(stdout, '(5x,a)') ' Program compiled on '//__DATE__//' '//__TIME__//'.          '
        write(stdout, '(5x,a)') '-------------------------------------------------------------'
        write(stdout, *) 
        write(stdout, '(1x,a,i0,a)') 'Using ', omp_get_max_threads(), ' threads.'
        write(stdout, *) 
        write(stdout, '(1x,a)') 'Input for bra orbitals/states: '
        write(stdout, '(5x,a,a)') 'Format: ', input_format_1
        write(stdout, '(5x,a,a)') 'Path: ', path1
        write(stdout, '(1x,a)') 'Input for ket orbitals/states: '
        write(stdout, '(5x,a,a)') 'Format: ', input_format_2
        write(stdout, '(5x,a,a)') 'Path: ', path2
    end if

    ! Calculate AO overlaps.
    time0 = omp_get_wtime()
    call ao_overlap(path1, input_format_1, path2, input_format_2, s_ao, center_atoms, center_pairs)
    call output_ao(outfile_ao, s_ao)
    time_ao = omp_get_wtime() - time0
    if (ao_stop) goto 999

    ! Calculate MO overlaps.
    time0 = omp_get_wtime()
    call mo_overlap(path1, input_format_1, path2, input_format_2, s_ao, s_mo_a, s_mo_b)
    call output_mo(outfile_mo, s_mo_a, s_mo_b)
    time_mo = omp_get_wtime() - time0
    if (mo_stop) goto 999

    ! Calculate CIS overlaps.
    time0 = omp_get_wtime()
    call cis_overlap(path1, input_format_1, path2, input_format_2, s_mo_a, s_mo_b, cis_algorithm,  &
    &                s_wf, norm_states, orth_states, wf_threshold, freeze_mo_norm, freeze_mo_norm_t)
    call output_orth(orth_overlap, s_wf)
    call output_phase(match_phase, s_wf)
    call output_wf(outfile_wf, s_wf)
    time_wf = omp_get_wtime() - time0

    ! Timings
    999 time_tot = omp_get_wtime() - time00
    if (print_level >= 2) then
        write(stdout, *) 
        write(stdout,'(1x,a)') 'Program time:'
        write(stdout, '(5x, a40, f14.4)') 'AO overlap                - time (sec):', time_ao
        write(stdout, '(5x, a40, f14.4)') 'MO overlap                - time (sec):', time_mo
        write(stdout, '(5x, a40, f14.4)') 'WF overlap                - time (sec):', time_wf
        write(stdout, '(5x, 40x, a14)') '--------------'
        write(stdout, '(5x, a40, f14.4)') 'Total                     - time (sec):', time_tot
    end if

    if (print_level >= 1) then
        write(stdout, *) 
        write(stdout, '(1x,a)') 'cis_overlap done                                             '
        write(stdout, '(1x,a)') '-------------------------------------------------------------'
    end if


end program cis_overlap_prog
