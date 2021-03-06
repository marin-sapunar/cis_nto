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
    subroutine molden_read_basis(fname, bs, fix_gto_norm)
        use basis_set_mod, only : basis_set
        use file_mod, only : check_string_in_file
        character(len=*), intent(in) :: fname
        type(basis_set), intent(out) :: bs
        integer, intent(in) :: fix_gto_norm
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
        select case(fix_gto_norm)
        case(3)
            call section_read_gto(readf, bs%n_bs, bs%bs, bs%center_i_bs, .true.)
        case default
            call section_read_gto(readf, bs%n_bs, bs%bs, bs%center_i_bs, .false.)
        end select
        call readf%close()
    end subroutine molden_read_basis


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: molden_read_mo
    !> @brief Read molecular orbitals from a molden file.
    !----------------------------------------------------------------------------------------------
    subroutine molden_read_mo(fname, mos, bs, fix_mo)
        use molecular_orbitals_mod, only : molecular_orbitals
        use basis_set_mod, only : basis_set
        character(len=*), intent(in) :: fname
        type(molecular_orbitals) :: mos
        type(basis_set), intent(in), optional :: bs
        integer, intent(in), optional :: fix_mo
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
        if (present(fix_mo)) then
            select case(fix_mo)
            case(0)
            case(1)
                call tm2molden_fix(bs, tmoa_c)
                if (nmo_b > 0) call tm2molden_fix(bs, tmob_c)
            case(2)
                call cfour_fix(bs, tmoa_c)
                if (nmo_b > 0) call cfour_fix(bs, tmob_c)
            case(3)
                call orca_fix(bs, tmoa_c)
                if (nmo_b > 0) call orca_fix(bs, tmob_c)
            end select
        end if
        call mos%init(ca=tmoa_c, ea=tmoa_e, oa=tmoa_o, cb=tmob_c, eb=tmob_e, ob=tmob_o)
    end subroutine molden_read_mo


    subroutine orca_fix(bs, moa)
        use basis_set_mod, only : basis_set
        use ang_mom_defs, only : amp
        use math_mod, only : factorial2
        type(basis_set), intent(in) :: bs
        real(dp), intent(inout) :: moa(:, :)
        integer :: i, j, k, pos, l, cbs
        integer, allocatable :: ls(:, :)
        real(dp) :: mult

        pos = 1
        do i = 1, bs%n_center
            cbs = bs%center_i_bs(i)
            do j = 1, bs%bs(cbs)%n_subshell
                l = bs%bs(cbs)%cg(j)%l
                ls = molden_sphe_order(l)
                do k = 1, amp%n_sphe(l)
                    if ((abs(ls(2, k)) > 2) .and. (abs(ls(2, k)) < 5)) then
                        moa(pos, :) = - moa(pos, :)
                    end if
                    pos = pos + 1
                end do
            end do
        end do
    end subroutine orca_fix


    subroutine tm2molden_fix(bs, moa)
        use basis_set_mod, only : basis_set
        use ang_mom_defs, only : amp
        use math_mod, only : factorial2
        type(basis_set), intent(in) :: bs
        real(dp), intent(inout) :: moa(:, :)
        integer :: i, j, pos, epos, cbs
        real(dp) :: mult

        pos = 1
        do i = 1, bs%n_center
            cbs = bs%center_i_bs(i)
            do j = 1, bs%bs(cbs)%n_subshell
                mult = sqrt(real(factorial2(2*bs%bs(cbs)%cg(j)%l-1), kind=dp))
                epos = pos+amp%n_cart(bs%bs(cbs)%cg(j)%l) - 1
                moa(pos:epos, :) = moa(pos:epos, :) * mult
                pos = epos + 1
            end do
        end do
    end subroutine tm2molden_fix


    subroutine cfour_fix(bs, moa)
        use basis_set_mod, only : basis_set
        use ang_mom_defs, only : amp
        use math_mod, only : factorial2
        type(basis_set), intent(in) :: bs
        real(dp), intent(inout) :: moa(:, :)
        integer :: i, j, k, pos, l, cbs
        integer, allocatable :: lc(:, :)
        real(dp) :: mult

        pos = 1
        do i = 1, bs%n_center
            cbs = bs%center_i_bs(i)
            do j = 1, bs%bs(cbs)%n_subshell
                l = bs%bs(cbs)%cg(j)%l
                lc = molden_cart_order(l)
                do k = 1, amp%n_cart(l)
                    mult = 1.0_dp
                    mult = mult * sqrt(real(factorial2(2*lc(1, k)-1), kind=dp))
                    mult = mult * sqrt(real(factorial2(2*lc(2, k)-1), kind=dp))
                    mult = mult * sqrt(real(factorial2(2*lc(3, k)-1), kind=dp))
                    moa(pos, :) = moa(pos, :) * mult
                    pos = pos + 1
                end do
            end do
        end do
    end subroutine cfour_fix




end module molden_read_mod
