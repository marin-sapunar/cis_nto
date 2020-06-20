!--------------------------------------------------------------------------------------------------
! PROGRAM: cis_nto_prog
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date April, 2019
!
!> @brief Convert between different geometry/AO/MO formats.
!> @todo Only molden format write subroutines implemented ATM.
!> @todo Add checks/options to read/convert only parts of the data that are available in the
!!       given input files.
!--------------------------------------------------------------------------------------------------
program nto_prog
    ! General
    use global_defs
    ! Variables
    use convert_variables
    ! I/O
    use convert_input_mod
    use read_all_mod
    use write_all_mod
    implicit none

    ! All program variables are defined in convert_variables

    ! Begin run.
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
        write(stdout, '(1x,a)') 'Reading orbitals/states: '
        write(stdout, '(5x,a,a)') 'Format: ', input_format
        write(stdout, '(5x,a,a)') 'Path: ', input_path
    end if

    call read_geom(input_format, input_path, geom, atsym, atnum)
    call read_basis(input_format, input_path, bs)
    call read_mo(input_format, input_path, mos, bs)

    call bs_trans%init(bs, bs%source_format, output_format)
    call bs_trans%transform(mos%ca, 1)
    if (allocated(mos%cb)) call bs_trans%transform(mos%cb, 1)
    call write_all(output_format, output_path, geom=geom, atom_symbol=atsym, atom_number=atnum,&
    &              basis=bs, mos=mos)



end program nto_prog
