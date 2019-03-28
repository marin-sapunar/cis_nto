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
    use file_mod
    use matrix_mod
    ! Chem
    use ao_overlap_mod
    use mo_overlap_mod
    use cis_overlap_mod
    use ang_mom_defs
    ! Output
    use overlap_output_mod
    implicit none


    ! Options
    character(len=:), allocatable :: path1
    character(len=:), allocatable :: path2
    character(len=:), allocatable :: input_format_1
    character(len=:), allocatable :: input_format_2
    logical :: ao_stop
    logical :: mo_stop
    character(len=4) :: cis_algorithm
    real(dp) :: wf_threshold
    logical :: norm_states
    logical :: orth_states
    logical :: orth_overlap
    logical :: match_phase
    logical :: center1
    logical :: center_pairs
    logical :: freeze_mo_norm
    real(dp) :: freeze_mo_norm_t
    character(len=:), allocatable :: outfile_ao
    character(len=:), allocatable :: outfile_mo
    character(len=:), allocatable :: outfile_wf
    ! Results
    real(dp), allocatable :: s_ao(:, :) !< Atomic orbital overlaps.
    real(dp), allocatable :: s_mo_a(:, :) !< Molecular orbital overlaps alpha.
    real(dp), allocatable :: s_mo_b(:, :) !< Molecular orbital overlaps beta.
    real(dp), allocatable :: s_wf(:, :) !< Wave function overlaps.
    ! Help
    integer, external :: omp_get_max_threads
    real(dp), external :: omp_get_wtime
    real(dp) :: time00, time0
    real(dp) :: time_ao, time_mo, time_wf, time_tot

    call set_defaults()
    call amp%init(7)
    time00 = omp_get_wtime()

    ! Begin run.
    call command_line_interface()
    if (print_level >= 1) then
        write(stdout, '(5x,a)') '-------------------------------------------------------------'
        write(stdout, '(5x,a)') '                     cis_overlap program                     '
        write(stdout, '(5x,a)') '                         version 1.0                         '
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
    call ao_overlap(path1, input_format_1, path2, input_format_2, s_ao, center1, center_pairs)
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


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: print_help
    !> @brief Print help message and exit program.
    !----------------------------------------------------------------------------------------------
    subroutine print_help()
        write(stdout, '(a)') 'usage: cis_overlap.exe [optional arguments] path1 path2'
        write(stdout, '(a)')
        write(stdout, '(a)') 'Calculate overlaps between two sets CIS type wave functions.'
        write(stdout, '(a)')
        write(stdout, '(a)') 'positional arguments:'
        write(stdout, '(a)') '  path1                       directory containing calculation for bra states      '
        write(stdout, '(a)') '  path2                       directory containing calculation for ket states      '
        write(stdout, '(a)')
        write(stdout, '(a)') 'optional arguments:'
        write(stdout, '(a)') '  -h, --help                  show this help message and exit                      '
        write(stdout, '(a)') '  -in1, --input-format-1      format of input files for bra orbitals/states        '
        write(stdout, '(32x,a,l1)') 'default: ', input_format_1
        write(stdout, '(a)') '  -in2, --input-format-2      format of input files for ket orbitals/states        '
        write(stdout, '(32x,a,l1)') 'default: ', input_format_2
        write(stdout, '(a)') '  -ao, --ao-stop              stop program after AO overlap calculation            '
        write(stdout, '(a)') '  -mo, --mo-stop              stop program after MO overlap calculation            '
        write(stdout, '(a)') '  -alg, --algorithm ALG       algorithm to use for the overlap calculation         '
        write(stdout, '(a)') '                              Available options:                                   '
        write(stdout, '(a)') '                                CIS (very slow, not recommended)                   '
        write(stdout, '(a)') '                                L2M                                                '
        write(stdout, '(a)') '                                NTO                                                '
        write(stdout, '(32x,a,a)') 'default: ', cis_algorithm
        write(stdout, '(a)') '  -t, --wf-threshold t        truncate wave functions using given threshold        '
        write(stdout, '(a)') '  -ns, --(no-)norm-states     renormalize input states before calculation          '
        write(stdout, '(32x,a,l1)') 'default: ', norm_states
        write(stdout, '(a)') '  -os, --(no-)orth-states     reorthogonalize input states before calculation      '
        write(stdout, '(32x,a,l1)') 'default: ', orth_states
        write(stdout, '(a)') '  -oo, --(no-)orth-overlap    orthogonalize overlap matrix                         '
        write(stdout, '(32x,a,l1)') 'default: ', orth_overlap
        write(stdout, '(a)') '  -mp, --(no-)match-phase     match phase between assigned bra/ket states          '
        write(stdout, '(32x,a,l1)') 'default: ', match_phase
        write(stdout, '(a)') '  -cp, --(no-)center-pairs    attempt to remove effect of basis set translation by '
        write(stdout, '(a)') '                              recentering pairs of AOs in AO overlap calculation   '
        write(stdout, '(a)') '                              (untested, not recommended)                          '
        write(stdout, '(32x,a,l1)') 'default: ', center_pairs
        write(stdout, '(a)') '  -fmn, --freeze-mo-norm t    freeze occupied ket MOs when their norm in bra basis '
        write(stdout, '(a)') '                              is smaller than given threshold. Same number of bra  '
        write(stdout, '(a)') '                              MOs with smallest norms in ket basis is also frozen. '
        write(stdout, '(a)') '                              Used when geometry of a small part of a system is    '
        write(stdout, '(a)') '                              significantly different between bra and ket states.  '
        write(stdout, '(a)') '                              (untested, not recommended)                          '
        write(stdout, '(32x,a,l1)') 'default: ', freeze_mo_norm
        write(stdout, '(a)') '  -sao, --sao-output file     output final ao overlap matrix to file               '
        write(stdout, '(32x,a,a)') 'default: ',  outfile_ao
        write(stdout, '(a)') '  -smo, --smo-output file     output final mo overlap matrix to file               '
        write(stdout, '(32x,a,a)') 'default: ',  outfile_mo
        write(stdout, '(a)') '  -swf, --swf-output file     output final wf overlap matrix to file               '
        write(stdout, '(32x,a,a)') 'default: ',  outfile_wf
        write(stdout, '(a)') '  -p, --print-level p         control output level of program (0 = quiet)          '
        write(stdout, '(32x,a,i0)') 'default: ', print_level
        stop
    end subroutine print_help


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: set_defaults
    !> @brief Set default values for program options.
    !----------------------------------------------------------------------------------------------
    subroutine set_defaults()
        ! CLI options.
        input_format_1 = 'turbomole'
        input_format_2 = 'turbomole'
        ao_stop = .false.
        mo_stop = .false.
        cis_algorithm = 'NTO'
        wf_threshold = 1.0_dp
        norm_states = .true.
        orth_states = .false.
        orth_overlap = .false.
        match_phase = .false.
        center1 = .false. !unused
        center_pairs = .false.
        freeze_mo_norm = .false.
        freeze_mo_norm_t = -1.0_dp
        outfile_ao = 'None'
        outfile_mo = 'None'
        outfile_wf = 's_wf'
        print_level = 2
        ! Set timings to 0
        time_ao = 0.0_dp
        time_mo = 0.0_dp
        time_wf = 0.0_dp
        time_tot = 0.0_dp
    end subroutine set_defaults


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: command_line_interface
    !> @brief Read command line options.
    !----------------------------------------------------------------------------------------------
    subroutine command_line_interface()
        integer :: narg, i
        character(len=1000) :: temp
        narg = command_argument_count()
        i = 0
        if (narg == 0) call print_help()
        do while (i < narg)
            i = i + 1
            if (i == narg) call print_help()
            call get_command_argument(i, temp)
            if (i == narg - 1) exit
            select case(temp)
            case('--help', '-h')
                call print_help()
            case('--input-format-1', '-in1')
                i = i + 1
                call get_command_argument(i, temp)
                input_format_1 = trim(adjustl(temp))
            case('--input-format-2', '-in2')
                i = i + 1
                call get_command_argument(i, temp)
                input_format_2 = trim(adjustl(temp))
            case('--ao-stop', '-ao')
                ao_stop = .true.
                if (outfile_ao == 'None') outfile_ao = 's_ao'
            case('--mo-stop', '-mo')
                mo_stop = .true.
                if (outfile_mo == 'None') outfile_mo = 's_mo'
            case('--algorithm', '-alg')
                i = i + 1
                call get_command_argument(i, temp)
                cis_algorithm = temp(1:4)
            case('--wf-threshold', '-t')
                i = i + 1
                call get_command_argument(i, temp)
                read(temp, *) wf_threshold
            case('--norm-states', '-ns')
                norm_states = .true.
            case('--no-norm-states', '-nns')
                norm_states = .false.
            case('--orth-states', '-os')
                orth_states = .true.
            case('--no-orth-states', '-nos')
                orth_states = .false.
            case('--orth-overlap', '-oo')
                orth_overlap = .true.
            case('--no-orth-overlap', '-noo')
                orth_overlap = .false.
            case('--match-phase', '-mp')
                match_phase = .true.
            case('--no-match-phase', '-nmp')
                match_phase = .false.
            case('--center-pairs', '-cp')
                center_pairs = .true.
            case('--no-center-pairs', '-ncp')
                center_pairs = .false.
            case('--freeze-mo-norm', '-fmn')
                i = i + 1
                call get_command_argument(i, temp)
                read(temp, *) freeze_mo_norm_t
                if (freeze_mo_norm_t > 0.0_dp) freeze_mo_norm = .true.
            case('--sao-output', '-sao')
                i = i + 1
                call get_command_argument(i, temp)
                outfile_ao = trim(adjustl(temp))
            case('--smo-output', '-smo')
                i = i + 1
                call get_command_argument(i, temp)
                outfile_mo = trim(adjustl(temp))
            case('--swf-output', '-swf')
                i = i + 1
                call get_command_argument(i, temp)
                outfile_wf = trim(adjustl(temp))
            case('--print-level', '-p')
                i = i + 1
                call get_command_argument(i, temp)
                read(temp, *) print_level
            end select
        end do
        call get_command_argument(narg-1, temp)
        path1 = trim(adjustl(temp))
        call get_command_argument(narg, temp)
        path2 = trim(adjustl(temp))
    end subroutine command_line_interface


end program cis_overlap_prog
