!----------------------------------------------------------------------------------------------
! MODULE: xyz_read_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date June, 2017
!
!> @brief Contains subroutines for reading data from an xyz format file.
!> @todo Add option to read multiple concatenated geometries.
!----------------------------------------------------------------------------------------------
module xyz_read_mod
    use global_defs
    use file_mod, only : reader
    implicit none


    private
    public :: xyz_read_geom


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: xyz_read_geom
    !> @brief Read the geometry and atom symbols from an xyz file.
    !----------------------------------------------------------------------------------------------
    subroutine xyz_read_geom(fname, title, geom, at_symbol)
        character(len=*), intent(in) :: fname !< Name of output file
        character(len=:), allocatable, intent(out), optional :: title !< Title line
        real(dp), allocatable :: geom(:) !< Geometry
        character(len=2), allocatable :: at_symbol(:) !< Atom symbols
        type(reader) :: readf
        integer :: natom, i

        call readf%open(fname, skip_empty=.false.)
        call readf%next()
        read(readf%line, *) natom
        allocate(geom(natom*3))
        allocate(at_symbol(natom))
        call readf%next()
        if (present(title)) title = readf%line
        do i = 1, natom
            call readf%next()
            read(readf%line, *) at_symbol(i), geom(3*i-2:3*i)
        end do
        call readf%close()
    end subroutine xyz_read_geom


end module xyz_read_mod
