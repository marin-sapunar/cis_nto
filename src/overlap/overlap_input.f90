!--------------------------------------------------------------------------------------------------
! MODULE: overlap_input_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @brief Hold subroutines for overlap program input handling.
!--------------------------------------------------------------------------------------------------
module overlap_input_mod
    use global_defs
    use overlap_variables
    implicit none


    private
    public :: print_help
    public :: set_defaults
    public :: command_line_interface


    character(len=4), parameter :: default_file_ao = 's_ao'
    character(len=4), parameter :: default_file_mo = 's_mo'
    character(len=4), parameter :: default_file_wf = 's_wf'


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
        write(stdout, '(a)') '  -in, --input-format         format of input files                                '
        write(stdout, '(32x,a,a)') 'default: ', input_format_1
        write(stdout, '(a)') '  -in1, --input-format-1      format of input files for bra orbitals/states        '
        write(stdout, '(32x,a,a)') 'default: ', input_format_1
        write(stdout, '(a)') '  -in2, --input-format-2      format of input files for ket orbitals/states        '
        write(stdout, '(32x,a,a)') 'default: ', input_format_2
        write(stdout, '(a)') '  -ao, --ao-stop              stop program after AO overlap calculation            '
        write(stdout, '(a)') '  -mo, --mo-stop              stop program after MO overlap calculation            '
        write(stdout, '(a)') '  -alg, --algorithm ALG       algorithm to use for the overlap calculation         '
        write(stdout, '(a)') '                              Available options:                                   '
        write(stdout, '(a)') '                                CIS (very slow, not recommended)                   '
        write(stdout, '(a)') '                                L2M                                                '
        write(stdout, '(a)') '                                NTO                                                '
        write(stdout, '(32x,a,a)') 'default: ', cis_algorithm
        write(stdout, '(a)') '  -dys, --dyson               calculate Dyson orbitals instead of overlaps         '
        write(stdout, '(32x,a,l1)') 'default: ', dyson_c
        write(stdout, '(a)') '  -t, --wf-threshold t        truncate wave functions using given threshold        '
        write(stdout, '(a)') '  -tnex, --truncate-nex n     truncate wave functions to n dominant excitations    '
        write(stdout, '(a)') '  -ns, --(no-)norm-states     renormalize input states before calculation          '
        write(stdout, '(32x,a,l1)') 'default: ', norm_states
        write(stdout, '(a)') '  -os, --(no-)orth-states     reorthogonalize input states before calculation      '
        write(stdout, '(32x,a,l1)') 'default: ', orth_states
        write(stdout, '(a)') '  -oo, --(no-)orth-overlap    orthogonalize overlap matrix                         '
        write(stdout, '(32x,a,l1)') 'default: ', orth_overlap
        write(stdout, '(a)') '  -mp, --(no-)match-phase     match phase between assigned bra/ket states          '
        write(stdout, '(32x,a,l1)') 'default: ', match_phase
        write(stdout, '(a)') '  -op, --operator op          operator whose expectation values are cacluclated.   '
        write(stdout, '(a)') '                              Currently, only overlaps can be calculated for wave  '
        write(stdout, '(a)') '                              functions.                                           '
        write(stdout, '(a)') '                              Available options:                                   '
        write(stdout, '(a)') '                                analytic expression (eg. "x * y^2 - z^2")          '
        write(stdout, '(a)') '                                  note: use quotation marks                        '
        write(stdout, '(a)') '                                  note: see read_operator_string subroutine docs   '
        write(stdout, '(a)') '                                overlap                                            '
        write(stdout, '(a)') '                                r^2 (equivalent to "x^2 + y^2 + z^2")              '
        write(stdout, '(32x,a,a)') 'default: ', operator_string
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
        use ang_mom_defs

        call amp%init(7)
        ! CLI options.
        input_format_1 = 'turbomole'
        input_format_2 = 'turbomole'
        ao_stop = .false.
        mo_stop = .false.
        cis_algorithm = 'NTO'
        dyson_c = .false.
        wf_threshold = 1.0_dp
        truncate_nex = 0
        norm_states = .true.
        orth_states = .false.
        orth_overlap = .false.
        match_phase = .false.
        center_atoms = .false. !unused
        center_pairs = .false.
        freeze_mo_norm = .false.
        freeze_mo_norm_t = -1.0_dp
        operator_string = 'overlap'
        outfile_ao = 'None'
        outfile_mo = 'None'
        outfile_wf = default_file_wf
        prefix_dyson = 'dys'
        print_level = 2
        ! Set timings to 0
        time_in = 0.0_dp
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
            if ((i == narg - 1) .and. (temp /= '--help') .and. (temp /= '-h')) exit
            select case(temp)
            case('--help', '-h')
                call print_help()
            case('--input-format', '-in')
                i = i + 1
                call get_command_argument(i, temp)
                input_format_1 = trim(adjustl(temp))
                input_format_2 = trim(adjustl(temp))
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
                if (outfile_ao == 'None') outfile_ao = default_file_ao
            case('--mo-stop', '-mo')
                mo_stop = .true.
                if (outfile_mo == 'None') outfile_mo = default_file_mo
            case('--algorithm', '-alg')
                i = i + 1
                call get_command_argument(i, temp)
                cis_algorithm = temp(1:4)
            case('--dyson', '-dys')
                dyson_c = .true.
            case('--no-dyson', '-ndys')
                dyson_c = .false.
            case('--wf-threshold', '-t')
                i = i + 1
                call get_command_argument(i, temp)
                read(temp, *) wf_threshold
            case('--truncate-nex', '-tnex')
                i = i + 1
                call get_command_argument(i, temp)
                read(temp, *) truncate_nex
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
            case('--operator', '-op')
                i = i + 1
                call get_command_argument(i, temp)
                operator_string = trim(adjustl(temp))
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

        if (operator_string /= 'overlap') then
            if (.not. (mo_stop .or. ao_stop)) then
                write(stderr, *) 'Warning. Only overlap operator is implemented for wave functions.'
                write(stderr, *) '  Setting mo_stop option to .true.'
                mo_stop = .true.
                if (outfile_mo == 'None') outfile_mo = default_file_mo
            end if
        end if
    end subroutine command_line_interface


end module overlap_input_mod
