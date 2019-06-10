!----------------------------------------------------------------------------------------------
! MODULE: xyz_write_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date June, 2017
!
!> @brief Contains subroutines for writing an xyz file
!----------------------------------------------------------------------------------------------
module xyz_write_mod
    use global_defs
    implicit none

    private
    public :: xyz_write


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: xyz_write
    !> @brief Write a geometry to a xyz format file.
    !----------------------------------------------------------------------------------------------
    subroutine xyz_write(outunit, title, geom, at_symbol)
        integer, intent(in) :: outunit !< Output unit
        character(len=*), intent(in), optional :: title !< Title line to print
        real(dp), intent(in) :: geom(:) !< Geometry
        character(len=2), intent(in) :: at_symbol(:) !< Atom symbols
        character(len=*), parameter :: fatm = '(1x, a2, 3(1x,e22.16))'
        character(len=:), allocatable :: ttitle
        real(dp) :: lenconv
        integer :: i

        lenconv = 1.88972612456 !< @todo: Add unit conversion constants module.
        ttitle = ''
        if (present(title)) ttitle = title
        write(outunit, '(i0)') size(geom)/3
        write(outunit, '(a)') ttitle
        do i = 1, size(at_symbol)
            write(outunit, fatm) at_symbol(i), geom(3*i-2:3*i)*lenconv
        end do
    end subroutine xyz_write


end module xyz_write_mod
