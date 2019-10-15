!--------------------------------------------------------------------------------------------------
! MODULE: nto_input_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @brief Hold subroutines for overlap program input handling.
!--------------------------------------------------------------------------------------------------
module nto_input_mod
    use global_defs
    use nto_variables
    implicit none


    private
    public :: print_help
    public :: set_defaults
    public :: command_line_interface


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: print_help
    !> @brief Print help message and exit program.
    !----------------------------------------------------------------------------------------------
    subroutine print_help()
        write(stdout, '(a)') 'usage: cis_nto.exe [optional arguments] path'
        write(stdout, '(a)')
        write(stdout, '(a)') 'Analyze CIS wave function in terms of NTOs.'
        write(stdout, '(a)')
        write(stdout, '(a)') 'positional arguments:'
        write(stdout, '(a)') '  path                        directory containing calculation                     '
        write(stdout, '(a)')
        write(stdout, '(a)') 'optional arguments:'
        write(stdout, '(a)') '  -h, --help                  show this help message and exit                      '
        write(stdout, '(a)') '  -in, --input-format         format of input files for bra orbitals/states        '
        write(stdout, '(32x,a,a)') 'default: ', input_format
        write(stdout, '(a)') '  -t, --wf-threshold t        truncate wave functions using given threshold        '
        write(stdout, '(a)') '  -tnex, --truncate-nex n     truncate wave functions to n dominant excitations    '
        write(stdout, '(a)') '  -ns, --(no-)norm-states     renormalize input states before calculation          '
        write(stdout, '(32x,a,l1)') 'default: ', norm_states
        write(stdout, '(a)') '  -os, --(no-)orth-states     reorthogonalize input states before calculation      '
        write(stdout, '(32x,a,l1)') 'default: ', orth_states
        write(stdout, '(a)') '  -cp, --(no-)center-pairs    attempt to remove effect of basis set translation by '
        write(stdout, '(a)') '                              recentering pairs of AOs in AO overlap calculation   '
        write(stdout, '(a)') '                              (untested, not recommended)                          '
        write(stdout, '(32x,a,l1)') 'default: ', center_pairs
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
        input_format = 'turbomole'
        output_format = 'molden_cart'
        outfile_nto = 'nto.molden'
        wf_threshold = 1.0_dp
        truncate_nex = 0
        norm_states = .true.
        orth_states = .false.
        center_atoms = .false. !unused
        center_pairs = .false.
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
            call get_command_argument(i, temp)
            if ((i == narg) .and. (temp /= '--help') .and. (temp /= '-h')) exit
            select case(temp)
            case('--help', '-h')
                call print_help()
            case('--input-format', '-in')
                i = i + 1
                call get_command_argument(i, temp)
                input_format = trim(adjustl(temp))
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
            case('--center-pairs', '-cp')
                center_pairs = .true.
            case('--no-center-pairs', '-ncp')
                center_pairs = .false.
            case('--print-level', '-p')
                i = i + 1
                call get_command_argument(i, temp)
                read(temp, *) print_level
            end select
        end do
        call get_command_argument(narg, temp)
        path = trim(adjustl(temp))
    end subroutine command_line_interface


end module nto_input_mod
