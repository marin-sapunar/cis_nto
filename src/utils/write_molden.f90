!----------------------------------------------------------------------------------------------
! MODULE: write_molden_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date June, 2017
!
! DESCRIPTION:
!> @brief Contains subroutines for writing a Molden format file.
!> @details
!! The Molden file format is organized into sections marked by keywords (ex. [Atoms]). Each 
!! subroutine in this module writes a section of the file. 
!----------------------------------------------------------------------------------------------
module write_molden_mod
    use global_defs
    implicit none

    private
    public :: write_molden_mo_single


contains

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: write_molden_mo_single
    !> @brief Write a single molecular orbital to the [MO] section.
    !----------------------------------------------------------------------------------------------
    subroutine write_molden_mo_single(outunit, n, sym, spin, ene, occup, coef)
        integer, intent(in) :: outunit !< Unit in which the file is open.
        integer, intent(in) :: n !< Number of the MO.
        character(len=*), intent(in) :: sym !< Symmetry group.
        integer, intent(in) :: spin !< Spin.
        real(dp), intent(in) :: ene !< Energy.
        real(dp), intent(in) :: occup !< Occupation number.
        real(dp), intent(in) :: coef(:) !< Coefficients.
        character(len=19), parameter :: fsym="(1x, 'Sym=', i6, a)"
        character(len=24), parameter :: fene="(1x, 'Ene=', 1x, e20.14)"
        character(len=20), parameter :: fspi="(1x, 'Spin=', 1x, a)"
        character(len=24), parameter :: focc="(1x, 'Occup=', 1x, f8.6)"
        character(len=16), parameter :: fcoe="(i6, 1x, e20.14)"
        integer :: i

        write(outunit, fsym) n, sym
        write(outunit, fene) ene
        if (spin == 1) write(outunit, fspi) 'Alpha'
        if (spin == 2) write(outunit, fspi) 'Beta'
        write(outunit, focc) occup
        do i = 1, size(coef)
            write(outunit, fcoe) i, coef(i)
        end do
    end subroutine write_molden_mo_single


end module write_molden_mod
