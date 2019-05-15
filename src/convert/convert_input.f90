!--------------------------------------------------------------------------------------------------
! MODULE: convert_input_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @brief Hold subroutines for convert program input handling.
!--------------------------------------------------------------------------------------------------
module convert_input_mod
    use global_defs
    use convert_variables
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
        write(stdout, '(a)') 'usage: convert.exe [optional arguments] input_path output_path'
        write(stdout, '(a)')
        write(stdout, '(a)') 'Convert geometry/basis/MOs to different format.'
        write(stdout, '(a)')
        write(stdout, '(a)') 'positional arguments:'
        write(stdout, '(a)') '  input_path                  path for input file(s) in original format            '
        write(stdout, '(a)') '  output_path                 path for output file(s)                              '
        write(stdout, '(a)')
        write(stdout, '(a)') 'optional arguments:'
        write(stdout, '(a)') '  -h, --help                  show this help message and exit                      '
        write(stdout, '(a)') '  -in, --input-format         format of input file(s)                              '
        write(stdout, '(32x,a,l1)') 'default: ', input_format
        write(stdout, '(a)') '  -out, --output-format       output format                                        '
        write(stdout, '(32x,a,l1)') 'default: ', output_format
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
        print_level = 0
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
                input_format = trim(adjustl(temp))
            case('--output-format', '-out')
                i = i + 1
                call get_command_argument(i, temp)
                output_format = trim(adjustl(temp))
            case('--print-level', '-p')
                i = i + 1
                call get_command_argument(i, temp)
                read(temp, *) print_level
            end select
        end do
        call get_command_argument(narg-1, temp)
        input_path = trim(adjustl(temp))
        call get_command_argument(narg, temp)
        output_path = trim(adjustl(temp))
    end subroutine command_line_interface


end module convert_input_mod
