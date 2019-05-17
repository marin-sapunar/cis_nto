!----------------------------------------------------------------------------------------------
! MODULE: molpro_output_read_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date March, 2019
!
!> @brief Contains subroutines for reading data from a molpro output file.
!> @details
!! The Molden file format is organized into sections marked by keywords (ex. [Atoms]). Each 
!! subroutine in this module reads a section from the file. It should be possible to call the
!! subroutines in any order, as the subroutines will attempt to read any missing required
!! information by calling each other.
!----------------------------------------------------------------------------------------------
module molpro_output_read_mod
    use global_defs
    use file_mod, only : reader
    implicit none

    private

    public :: molpro_output_read_geom
    public :: molpro_output_read_basis
    public :: molpro_output_read_mo
    public :: molpro_sphe_order
   !public :: molden_cart_order


contains



    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: molpro_output_read_geom
    !> @brief Read the geometry and atom numbers/symbols from a molden file.
    !----------------------------------------------------------------------------------------------
    subroutine molpro_output_read_geom(fname, geom, at_symbol, at_number)
        use string_mod, only : tolower
        character(len=*), intent(in) :: fname
        real(dp), allocatable :: geom(:) !< Geometry.
        character(len=2), allocatable :: at_symbol(:) !< Atom symbols.
        integer, allocatable :: at_number(:) !< Atomic numbers.

        type(reader) :: readf
        integer :: natom, i
        real(dp) :: charge

        call readf%open(fname, skip_empty = .false.)
        call find_n_atom(readf, natom)
        allocate(geom(3*natom))
        allocate(at_symbol(natom))
        allocate(at_number(natom))
        do i = 1, natom
            call readf%next()
            call readf%parseline()
            at_symbol(i) = tolower(readf%args(2)%s)
            read(readf%args(3)%s, *) charge
            at_number(i) = int(charge)
            read(readf%args(4)%s, *) geom(3*i-2)
            read(readf%args(5)%s, *) geom(3*i-1)
            read(readf%args(6)%s, *) geom(3*i)
        end do
        call readf%close()
    end subroutine molpro_output_read_geom


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: molpro_output_read_basis
    !> @brief Read basis set information from a molden file.
    !----------------------------------------------------------------------------------------------
    subroutine molpro_output_read_basis(fname, bs)
        use basis_set_mod, only : basis_set
        use file_mod, only : check_string_in_file
        character(len=*), intent(in) :: fname
        type(basis_set), intent(out) :: bs
        type(reader) :: readf
        integer :: i

        call readf%open(fname, skip_empty = .false.)
        call find_n_atom(readf, bs%n_center)
        bs%n_bs = bs%n_center
        allocate(bs%bs(bs%n_bs))
        allocate(bs%center_i_bs(1:bs%n_bs), source=[ (i, i = 1, bs%n_bs ) ])
        bs%sphe_mo = .true.
        bs%source_format = 'molpro'
        call readf%rewind()
        call readf%go_to_keyword('BASIS DATA')
        call readf%next() ! Empty line
        call readf%next() ! Nr Sym  Nuc  Type         Exponents   Contraction coefficients
        call readf%next() ! Empty line
        call section_read_basis(readf, bs%n_bs, bs%bs)
        call readf%close()
    end subroutine molpro_output_read_basis


    subroutine section_read_basis(readf, nabas, abas)
        use basis_set_mod, only : basis_set_single
        use ang_mom_defs
        type(reader), intent(inout) :: readf
        integer, intent(in) :: nabas !< Number of atoms.
        type(basis_set_single), intent(out) :: abas(nabas)
        integer, parameter :: maxnfunc = 100
        integer, parameter :: maxnprim = 25
        integer :: nprim(maxnfunc)
        real(dp) :: zeta(maxnfunc, maxnprim)
        real(dp) :: beta(maxnfunc, maxnprim)
        character :: orbtype(maxnfunc)
        logical :: skip
        integer :: c_prim, c_block_len, c_block_start, c_block_end, bi, i, atom, start
        character(len=1) :: orb
        character(len=1) :: chkm

        bi = 1
        c_prim = 0
        c_block_len = 0
        c_block_end = 0
        do
            call readf%next()
            if (readf%line == '') exit
            call readf%parseline(' ')
            c_prim = c_prim + 1
            if (c_prim <= c_block_len) then
                start = 3
                if (skip) cycle
            else if (readf%args(2)%s /= 'A') then
                start = 1
                if (skip) cycle
            else 
                start = 5
                ! Check start of new atom block
                read(readf%args(3)%s, *) atom
                if (atom /= bi) then
                    if (bi /= 0) call abas(bi)%init(c_block_end, orbtype, nprim, zeta, beta)
                    bi = atom
                    c_block_len = 0
                    c_block_end = 0
                end if
                ! Check if function block should be skipped
                read(readf%args(start-1)%s, '(1x,1a)') orb
                if (orb /= 's') read(readf%args(start-1)%s, '(2x,a)') chkm
                skip = .false.
                select case(orb)
                case('s')
                    skip = .false.
                case('p')
                    if (chkm /= 'x') skip = .true.
                case default
                    if (chkm /= '0') skip = .true.
                end select
                if (skip) cycle

                c_block_len = readf%narg - start
                c_block_start = c_block_end + 1
                c_block_end = c_block_start + c_block_len - 1
                c_prim = 1
            end if

            nprim(c_block_start:c_block_end) = c_prim
            orbtype(c_block_start:c_block_end) = orb
            do i = 0, c_block_len - 1
                read(readf%args(start)%s, *) zeta(c_block_start+i, c_prim)
                read(readf%args(start+i+1)%s, *) beta(c_block_start+i, c_prim)
            end do
        end do
        call abas(bi)%init(c_block_end, orbtype, nprim, zeta, beta)
    end subroutine section_read_basis



    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: molpro_output_read_mo
    !> @brief Read molecular orbitals from a molpro output file.
    !> @details
    !! The input has to have the thrprint=-1 keyword to print the full molecular orbitals.
    !----------------------------------------------------------------------------------------------
    subroutine molpro_output_read_mo(fname, mos)
        use molecular_orbitals_mod, only : molecular_orbitals
        character(len=*), intent(in) :: fname
        type(molecular_orbitals) :: mos
        real(dp), allocatable :: moa_c(:, :)
        real(dp), allocatable :: mob_c(:, :)
        real(dp), allocatable :: moa_e(:)
        real(dp), allocatable :: mob_e(:)
        real(dp), allocatable :: moa_o(:)
        real(dp), allocatable :: mob_o(:)
        type(reader) :: readf
        integer :: nmo_a, nmo_b, nbas

        call readf%open(fname, skip_empty = .false.)
        call readf%go_to_keyword('MOLECULAR ORBITALS, SYMMETRY 1:')
        call readf%next() ! ========== line
        call count_mo(readf, nmo_a, nbas)

        allocate(moa_c(nbas, nmo_a))
        allocate(moa_e(nmo_a))
        allocate(moa_o(nmo_a))
        call section_read_mo(readf, nbas, nmo_a, moa_e, moa_c)
        call readf%close()
        call mos%init(ca=moa_c, ea=moa_e, oa=moa_o)
    end subroutine molpro_output_read_mo


    subroutine count_mo(readf, nmo, nao)
        type(reader), intent(inout) :: readf
        integer, intent(out) :: nmo
        integer, intent(out) :: nao
        integer :: iline

        iline = readf%line_num
        call readf%next() ! Empty line
        nmo = 0
        do
            call readf%next()
            call readf%parseline(' ')
            if (readf%line == '') exit
            if (readf%args(1)%s /= 'Orb.') exit
            nmo = nmo + readf%narg - 2
            call readf%next() ! Energies
            call readf%next() ! Empty line
            call readf%next() ! Nr  Atom  Typ     Orbital Coefficients:
            call readf%go_to_keyword('', count_lines=nao)
        end do
        call readf%rewind()
        call readf%go_to_line(iline)
    end subroutine count_mo



    subroutine section_read_mo(readf, nbas, nmo, mo_e, mo_c)
        type(reader), intent(inout) :: readf
        integer, intent(in) :: nbas !< Number of basis functions
        integer, intent(in) :: nmo !< Number of molecular orbitals
        real(dp), intent(out) :: mo_e(nmo) !< Molecular orbital energies
        real(dp), intent(out) :: mo_c(nbas, nmo) !< Molecular orbital coefficients
        integer :: i, j, k, n_col

        i = 0
        do
            call readf%next() ! empty
            call readf%next() ! nums
            call readf%next() ! energ
            call readf%parseline(' ')
            n_col = readf%narg - 2
            do k = 1, n_col
                read(readf%args(2+k)%s, *) mo_e(i+k)
            end do
            call readf%next() ! empty
            call readf%next() ! nr atom typ orb coeff
            do j = 1, nbas
                call readf%next()
                call readf%parseline(' ')
                do k = 1, n_col
                    read(readf%args(3+k)%s, *) mo_c(j, i+k)
                end do
            end do
            i = i + n_col
            if (i == nmo) exit
        end do
    end subroutine section_read_mo


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: molpro_sphe_order
    !> @brief Return ordered ang. mom. numbers of spherical bfs within subshell of given l.
    !----------------------------------------------------------------------------------------------
    function molpro_sphe_order(ll) result(lm)
        use ang_mom_defs, only : amp, sphe_nums
        integer, intent(in) :: ll
        integer, allocatable :: lm(:, :)
        allocate(lm(2, amp%n_sphe(ll)))
        lm(1, :) = ll
        select case(ll)
        case(0)
            lm(2, :) = [0]
        case(1)
            lm(2, :) = [1, -1, 0]
        case(2)
            lm(2, :) = [0, -2, 1, 2, -1]
        case(3)
            lm(2, :) = [1, -1, 0, 3, -2, -3, 2]
        case(4)
            lm(2, :) = [0, -2, 1, 4, -1, 2, -4, 3, -3]
        case default
            stop 'molpro_sphe_order not implemented'
        end select
    end function molpro_sphe_order


    subroutine find_n_atom(readf, n)
        type(reader), intent(inout) :: readf
        integer, intent(out) :: n
        logical :: check

        call readf%go_to_keyword('ATOMIC COORDINATES', found=check)
        if (.not. check) then
            write(stderr, *) ' Error reading molpro output.'
            write(stderr, *) '     "ATOMIC COORDINATES" not found.'
            stop
        end if
        call readf%next() ! Empty line
        call readf%next() ! NR  ATOM    CHARGE       X              Y              Z
        call readf%next() ! Empty line
        n = readf%count_lines_to_keyword('')
    end subroutine find_n_atom


end module molpro_output_read_mod
