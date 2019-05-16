!----------------------------------------------------------------------------------------------
! MODULE: read_molden_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date June, 2017
!
!> @brief Contains subroutines for reading data from a Molden format file.
!> @details
!! The Molden file format is organized into sections marked by keywords (ex. [Atoms]). Each 
!! subroutine in this module reads a section from the file. It should be possible to call the
!! subroutines in any order, as the subroutines will attempt to read any missing required
!! information by calling each other.
!----------------------------------------------------------------------------------------------
module read_molden_mod
    use global_defs
    use file_mod, only : reader
    implicit none

    private

    public :: read_natom_molden
    public :: read_geom_molden
    public :: read_basis_molden
    public :: read_mo_molden
    public :: molden_sphe_order
    public :: molden_cart_order

contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_natom_molden
    !> @brief Read the number of atoms from a molden file.
    !----------------------------------------------------------------------------------------------
    subroutine read_natom_molden(fname, natom)
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
            write(stderr, *) ' Error in read_natom_molden subroutine.'
            write(stderr, *) ' Failed to read number of atoms.'
            stop
        end if
900     call readf%close()
    end subroutine read_natom_molden


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_natom_molden
    !> @brief Read the geometry and atom numbers/symbols from a molden file.
    !----------------------------------------------------------------------------------------------
    subroutine read_geom_molden(fname, geom, at_symbol, at_number)
        character(len=*), intent(in) :: fname
        real(dp), allocatable :: geom(:) !< Geometry.
        character(len=2), allocatable :: at_symbol(:) !< Atom symbols.
        integer, allocatable :: at_number(:) !< Atomic numbers.

        type(reader) :: readf
        logical :: check
        integer :: natom

        call read_natom_molden(fname, natom)
        allocate(geom(natom*3))
        allocate(at_symbol(natom))
        allocate(at_number(natom))
        call readf%open(fname, skip_empty = .false.)
        call readf%go_to_keyword('[Atoms]', found=check, case_sensitive=.False.)
        if (.not. check) then
            write(stderr, *) ' Error in read_natom_molden subroutine.'
            write(stderr, *) ' [Atoms] section not found.'
            stop
        end if
        call section_read_atoms(readf, natom, at_symbol, at_number, geom)
        call readf%close()
    end subroutine read_geom_molden


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_basis_molden
    !> @brief Read basis set information from a molden file.
    !----------------------------------------------------------------------------------------------
    subroutine read_basis_molden(fname, bs)
        use basis_set_mod, only : basis_set
        use file_mod, only : check_string_in_file
        character(len=*), intent(in) :: fname
        type(basis_set), intent(out) :: bs
        type(reader) :: readf
        integer :: ncart, nsphe


        call read_natom_molden(fname, bs%n_center)
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
        call section_read_gto(readf, bs%n_bs, bs%bs, bs%center_i_bs, ncart, nsphe)
        call readf%close()
        bs%n_bf_cart = ncart
        bs%n_bf_sphe = nsphe
    end subroutine read_basis_molden


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_mo_molden
    !> @brief Read molecular orbitals from a molden file.
    !----------------------------------------------------------------------------------------------
    subroutine read_mo_molden(fname, mos)
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
        call molden_count_mo(readf, nmo_a, nmo_b, nbas)
        allocate(tmoa_c(nbas, nmo_a))
        allocate(tmoa_e(nmo_a))
        allocate(tmoa_o(nmo_a))
        allocate(tmob_c(nbas, nmo_b))
        allocate(tmob_e(nmo_b))
        allocate(tmob_o(nmo_b))
        call section_read_mo(readf, nbas, nmo_a, tmoa_e, tmoa_o, tmoa_c, nmo_b, tmob_e, tmob_o, tmob_c)
        call readf%close()
        call mos%init(ca=tmoa_c, ea=tmoa_e, oa=tmoa_o, cb=tmob_c, eb=tmob_e, ob=tmob_o)
    end subroutine read_mo_molden


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_atoms
    !
    !> @brief Read the molden [Atoms] section.
    !> @details
    !! Read the geometry (and atom symbols/numbers).
    !----------------------------------------------------------------------------------------------
    subroutine section_read_atoms(readf, natom, asym, anum, geom)
        type(reader), intent(inout) :: readf !< Reader.
        integer, intent(in) :: natom !< Number of atoms.
        character(len=2), intent(out) :: asym(natom) !< Atom symbols.
        integer, intent(out) :: anum(natom) !< Atomic numbers.
        real(dp), intent(out) :: geom(3*natom)
        integer :: i, idum
        real(dp) :: lenconv

        lenconv = 1.0_dp
        call readf%parseline(' ')
        if (readf%narg > 1) then
            select case(readf%args(2)%s)
            case('Angs')
                lenconv =  1.88972612456_dp !< @todo conversion constants module
            end select
        end if
        do i = 1, natom
            call readf%next()
            read(readf%line, *) asym(i), idum, anum(i), geom(3*i-2:3*i)
        end do
        geom = geom * lenconv
    end subroutine section_read_atoms


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_gto
    !
    !> @brief Read the molden [GTO] section.
    !> @details
    !! Read the full basis set of each atom into an atombasis array.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_gto(readf, nabas, abas, abas_index, tot_ncart, tot_nsphe)
        use basis_set_mod, only : basis_set_single
        use ang_mom_defs
        type(reader), intent(inout) :: readf
        integer, intent(in) :: nabas !< Number of atoms.
        integer, intent(out) :: abas_index(nabas)
        type(basis_set_single), intent(out) :: abas(nabas)
        integer, intent(out) :: tot_ncart
        integer, intent(out) :: tot_nsphe
        character, parameter :: orborder(6) = ['s', 'p', 'd', 'f', 'g', 'h']
        integer, parameter :: maxnfunc = 100
        integer, parameter :: maxnprim = 25
        integer :: nfunc
        integer :: ncart
        integer :: nsphe
        integer :: nprim(maxnfunc)
        real(dp) :: zeta(maxnfunc, maxnprim)
        real(dp) :: beta(maxnfunc, maxnprim)
        character :: orbtype(maxnfunc)
        integer :: i
        integer :: j
       
        tot_ncart = 0 
        tot_nsphe = 0
        do i = 1, nabas
            call readf%next()
            read(readf%line, *) abas_index(i)
            nfunc = 0
            ncart = 0
            nsphe = 0
            do
                call readf%next()
                if (readf%line == '') exit
                call readf%parseline(' ')
                nfunc = nfunc + 1
                read(readf%args(1)%s, *) orbtype(nfunc)
                read(readf%args(2)%s, *) nprim(nfunc)
                do j = 1, nprim(nfunc)
                    call readf%next()
                    read(readf%line, *) zeta(nfunc, j), beta(nfunc, j)
                end do
            end do
          
            abas(i)%n_subshell = nfunc
            if (allocated(abas(i)%cg)) deallocate(abas(i)%cg)
            allocate(abas(i)%cg(nfunc))
            do j = 1, nfunc
                allocate(abas(i)%cg(j)%z(nprim(j)))
                allocate(abas(i)%cg(j)%b(nprim(j)))
                abas(i)%cg(j)%typ = orbtype(j)
                abas(i)%cg(j)%l = amp%l(orbtype(j))
                abas(i)%cg(j)%n_prim = nprim(j)
                abas(i)%cg(j)%z = zeta(j, 1:nprim(j))
                abas(i)%cg(j)%b = beta(j, 1:nprim(j))
                call abas(i)%cg(j)%norm_b()
                nsphe = nsphe + amp%n_sphe(abas(i)%cg(j)%l)
                ncart = ncart + amp%n_cart(abas(i)%cg(j)%l)
            end do
            abas(i)%n_cart = ncart
            abas(i)%n_sphe = nsphe
            tot_ncart = tot_ncart + ncart
            tot_nsphe = tot_nsphe + nsphe
        end do
    end subroutine section_read_gto


    subroutine molden_count_mo(readf, nmo_a, nmo_b, nbas)
        type(reader), intent(inout) :: readf
        integer, intent(out) :: nmo_a
        integer, intent(out) :: nmo_b
        integer, intent(out) :: nbas
        integer :: iline, chk

        iline = readf%line_num
 outer: do
            call readf%next(abort_on_eof=.false.)
            call readf%parseline()
            if (readf%args(1)%s /= '1') cycle
            do
                call readf%next(abort_on_eof=.false.)
                read(readf%line, *, iostat=chk) nbas
                if (chk /= 0) exit outer
            end do
        end do outer
        call readf%go_to_line(iline)
        nmo_a = readf%count_keyword_appearances('Spin= Alpha')
        nmo_b = readf%count_keyword_appearances('Spin= Beta')
    end subroutine molden_count_mo


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_mo
    !
    !> @brief Read the molden [MO] section.
    !> @details
    !! Read the full set of molecular orbitals from a molden file.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_mo(readf, nbas, nmo_a, mo_e_a, mo_o_a, mo_c_a, nmo_b, mo_e_b, mo_o_b, &
    &                          mo_c_b)
        type(reader), intent(inout) :: readf
        integer, intent(in) :: nbas !< Number of basis functions.
        integer, intent(in) :: nmo_a !< Number of molecular orbitals alpha.
        real(dp), intent(out) :: mo_e_a(nmo_a) !< Molecular orbital energies alpha.
        real(dp), intent(out) :: mo_o_a(nmo_a) !< Molecular orbital energies alpha.
        real(dp), intent(out) :: mo_c_a(nbas, nmo_a) !< Molecular orbital coefficients alpha.
        integer, intent(in) :: nmo_b !< Number of molecular orbitals beta.
        real(dp), intent(out) :: mo_e_b(nmo_b) !< Molecular orbital energies beta.
        real(dp), intent(out) :: mo_o_b(nmo_b) !< Molecular orbital energies beta.
        real(dp), intent(out) :: mo_c_b(nbas, nmo_b) !< Molecular orbital coefficients beta.
        integer :: cs
        real(dp) :: ce, co
        integer :: cmo_a
        integer :: cmo_b
        integer :: i, j, ca, cb

        cmo_a = 0
        cmo_b = 0
        i = 0
        ca = 0
        cb = 0
        do while (i < nmo_a + nmo_b)
            call readf%next()
            call readf%parseline('= ')
            select case(readf%args(1)%s) 
            case('Spin')
                select case(readf%args(2)%s)
                case('Alpha')
                    cs = 1
                    cmo_a = cmo_a + 1
                case('Beta')
                    cs = 2
                    cmo_b = cmo_b + 1
                end select
            case('Ene')
                read(readf%args(2)%s, *) ce
            case('Occup')
                read(readf%args(2)%s, *) co
            case('1')
                i = i + 1
                if (cs == 1) then
                    ca = ca + 1
                    read(readf%args(2)%s, *) mo_c_a(1, cmo_a)
                    mo_e_a(ca) = ce
                    mo_o_a(ca) = co
                else
                    cb = cb + 1
                    read(readf%args(2)%s, *) mo_c_b(1, cmo_b)
                    mo_e_b(cb) = ce
                    mo_o_b(cb) = co
                end if
                do j = 2, nbas
                    call readf%next()
                    call readf%parseline('= ')
                    if (cs == 1) then
                        read(readf%args(2)%s, *) mo_c_a(j, cmo_a)
                    else
                        read(readf%args(2)%s, *) mo_c_b(j, cmo_b)
                    end if
                end do
            end select
        end do
    end subroutine section_read_mo


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: molden_cart_order
    !> @brief Return ordered ang. mom. numbers of cartesian bfs within subshell of given l.
    !----------------------------------------------------------------------------------------------
    function molden_cart_order(ll) result(lxyz)
        use ang_mom_defs, only : amp, cart_nums
        integer, intent(in) :: ll
        integer, allocatable :: lxyz(:, :)
        integer, parameter :: f_order(3,10) = reshape( [ &
        &      3, 0, 0, 1, 2, 2, 1, 0, 0, 1, &
        &      0, 3, 0, 2, 1, 0, 0, 1, 2, 1, &
        &      0, 0, 3, 0, 0, 1, 2, 2, 1, 1], &
        &     shape(f_order), order=[2,1])


        allocate(lxyz(3, amp%n_cart(ll)))
        select case(ll)
        case(3)
            lxyz = f_order
        case default
            lxyz = cart_nums(ll)
        end select
    end function molden_cart_order


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: molden_sphe_order
    !> @brief Return ordered ang. mom. numbers of spherical bfs within subshell of given l.
    !----------------------------------------------------------------------------------------------
    function molden_sphe_order(ll) result(lm)
        use ang_mom_defs, only : amp, sphe_nums
        integer, intent(in) :: ll
        integer, allocatable :: lm(:, :)
        allocate(lm(2, amp%n_sphe(ll)))
        lm = sphe_nums(ll)
    end function molden_sphe_order


end module read_molden_mod
