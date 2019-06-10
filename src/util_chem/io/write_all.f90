!----------------------------------------------------------------------------------------------
! MODULE: write_all
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date May, 2019
!
!> @brief Contains subroutines for writing data to a variety of quantum chemistry formats.
!----------------------------------------------------------------------------------------------
module write_all_mod
    use global_defs
    use molden_write_mod
    implicit none


    private
    public :: write_all


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: write_all
    !> @brief Write all available data to a given format.
    !< @todo make input variables optional
    !----------------------------------------------------------------------------------------------
    subroutine write_all(dformat, path, title, geom, atom_symbol, atom_number, basis, mos)
        use element_mod
        use basis_set_mod, only : basis_set
        use molecular_orbitals_mod, only : molecular_orbitals
        character(len=*), intent(in) :: dformat
        character(len=*), intent(in) :: path
        character(len=*), intent(in), optional :: title
        real(dp), intent(in) :: geom(:)
        character(len=2), intent(in) :: atom_symbol(:)
        integer, intent(in) :: atom_number(:)
        type(basis_set), intent(in) :: basis
        type(molecular_orbitals), intent(in) :: mos
        integer :: outunit

        select case(dformat)
        case('molden_cart')
            open(newunit=outunit, file=path, action='write')
            write(outunit, '(a)') '[Molden Format]'
            write(outunit, '(a)') '[Title]'
            if (present(title)) then
                write(outunit, '(a)') title
            else
                write(outunit, *) 'TheoChemUtil write'
            end if
            call molden_write_atoms(outunit, geom, atom_symbol, atom_number)
            call molden_write_gto(outunit, basis)
            call molden_write_mo(outunit, mos)
        case default
            write(stderr, *) 'Error in write_all subroutine.'
            write(stderr, *) 'Input format not implemented: '//trim(adjustl(dformat))
        end select
    end subroutine write_all


end module write_all_mod
