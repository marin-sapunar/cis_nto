!----------------------------------------------------------------------------------------------
! MODULE: molden_read_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date June, 2017
!
!> @brief Contains subroutines for reading data from a Molden format file.
!> @details
!! The Molden file format is organized into sections marked by keywords (ex. [Atoms]). Each 
!! subroutine in this module reads a section from the file. It should be possible to call the
!! subroutines in any order, as the subroutines will attempt to read any missing required
!! information by calling each other.
!> @todo Check different formats for writing the geometry in molden files (such as [fr-coord])
!----------------------------------------------------------------------------------------------
module molden_read_mod
    use global_defs
    use molden_sections_mod
    use molden_definitions_mod
    use file_mod, only : reader
    implicit none


    private
    public :: molden_read_natom
    public :: molden_read_geom
    public :: molden_read_basis
    public :: molden_read_mo


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: molden_read_natom
    !> @brief Read the number of atoms from a molden file.
    !----------------------------------------------------------------------------------------------
    subroutine molden_read_natom(fname, natom)
        character(len=*), intent(in) :: fname
        integer, intent(out) :: natom
        type(reader) :: readf
        logical :: check

        call readf%open(fname, skip_empty=.false.)
        call readf%go_to_keyword('[n_atoms]', found=check)
        if (check) then
            call readf%next()
            read(readf%line, *) natom
            goto 900
        end if
        call readf%rewind()
        call readf%go_to_keyword('[atoms]', found=check, case_sensitive=.False.)
        if (check) then
            call readf%go_to_keyword('[', count_lines=natom, found=check)
            goto 900
        else
            write(stderr, *) ' Error in molden_read_natom subroutine.'
            write(stderr, *) ' Failed to read number of atoms.'
            stop
        end if
900     call readf%close()
    end subroutine molden_read_natom


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: molden_read_geom
    !> @brief Read the geometry and atom numbers/symbols from a molden file.
    !----------------------------------------------------------------------------------------------
    subroutine molden_read_geom(fname, geom, at_symbol, at_number)
        character(len=*), intent(in) :: fname
        real(dp), allocatable :: geom(:) !< Geometry.
        character(len=2), allocatable :: at_symbol(:) !< Atom symbols.
        integer, allocatable :: at_number(:) !< Atomic numbers.

        type(reader) :: readf
        logical :: check
        integer :: natom

        call molden_read_natom(fname, natom)
        allocate(geom(natom*3))
        allocate(at_symbol(natom))
        allocate(at_number(natom))
        call readf%open(fname, skip_empty = .false.)
        call readf%go_to_keyword('[Atoms]', found=check, case_sensitive=.False.)
        if (.not. check) then
            write(stderr, *) ' Error in molden_read_geom subroutine.'
            write(stderr, *) ' [Atoms] section not found.'
            stop
        end if
        call section_read_atoms(readf, natom, at_symbol, at_number, geom)
        call readf%close()
    end subroutine molden_read_geom


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: molden_read_basis
    !> @brief Read basis set information from a molden file.
    !----------------------------------------------------------------------------------------------
    subroutine molden_read_basis(fname, bs)
        use basis_set_mod, only : basis_set
        use file_mod, only : check_string_in_file
        character(len=*), intent(in) :: fname
        type(basis_set), intent(out) :: bs
        type(reader) :: readf
        integer :: ncart, nsphe


        call molden_read_natom(fname, bs%n_center)
        bs%n_bs = bs%n_center
        allocate(bs%bs(bs%n_bs))
        allocate(bs%center_i_bs(bs%n_bs))
        bs%sphe_mo = check_string_in_file(fname, '[5D') !> @todo this assumes no mixed sphe/cart, 
                               !! should throw error in case of [5D10F] or similar (or handle it).
        if (bs%sphe_mo) then
            bs%source_format = 'molden_sphe'
        else
            bs%source_format = 'molden_cart'
        end if

        call readf%open(fname, skip_empty = .false.)
        call readf%go_to_keyword('[GTO]')
        call section_read_gto(readf, bs%n_bs, bs%bs, bs%center_i_bs)
        call readf%close()
    end subroutine molden_read_basis


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: molden_read_mo
    !> @brief Read molecular orbitals from a molden file.
    !----------------------------------------------------------------------------------------------
    subroutine molden_read_mo(fname, mos)
        use molecular_orbitals_mod, only : molecular_orbitals
        character(len=*), intent(in) :: fname
        type(molecular_orbitals) :: mos
        real(dp), allocatable :: tmoa_c(:, :)
        real(dp), allocatable :: tmoa_e(:)
        real(dp), allocatable :: tmoa_o(:)
        real(dp), allocatable :: tmob_c(:, :)
        real(dp), allocatable :: tmob_e(:)
        real(dp), allocatable :: tmob_o(:)
        type(reader) :: readf
        integer :: nmo_a, nmo_b, nbas

        call readf%open(fname, skip_empty = .false.)
        call readf%go_to_keyword('[MO]')
        call section_count_mo(readf, nmo_a, nmo_b, nbas)
        allocate(tmoa_c(nbas, nmo_a))
        allocate(tmoa_e(nmo_a))
        allocate(tmoa_o(nmo_a))
        allocate(tmob_c(nbas, nmo_b))
        allocate(tmob_e(nmo_b))
        allocate(tmob_o(nmo_b))
        call section_read_mo(readf, nbas, nmo_a, tmoa_e, tmoa_o, tmoa_c, nmo_b, tmob_e, tmob_o, tmob_c)
        call readf%close()
        call mos%init(ca=tmoa_c, ea=tmoa_e, oa=tmoa_o, cb=tmob_c, eb=tmob_e, ob=tmob_o)
    end subroutine molden_read_mo


end module molden_read_mod
