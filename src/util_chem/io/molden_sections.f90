!----------------------------------------------------------------------------------------------
! MODULE: molden_sections_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date June, 2017
!
!> @brief Subroutines for reading sections of molden files.
!----------------------------------------------------------------------------------------------
module molden_sections_mod
    use global_defs
    use file_mod, only : reader
    implicit none


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_atoms
    !> @brief Read the molden [Atoms] section.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_atoms(readf, natom, asym, anum, geom)
        use string_mod, only : tolower
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
            select case(tolower(readf%args(2)%s))
            case('angs')
                lenconv = 1.88972612456_dp !< @todo conversion constants module
            case('au')
                lenconv = 1.0_dp
            !< @todo Other possible ways to set units in molden format?
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
    !> @brief Read the molden [GTO] section.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_gto(readf, nabas, abas, abas_index, fix_gto_norm)
        use basis_set_mod, only : basis_set_single
        use cgto_mod, only : gto_norms
        use ang_mom_defs
        type(reader), intent(inout) :: readf
        integer, intent(in) :: nabas !< Number of atoms.
        integer, intent(out) :: abas_index(nabas)
        type(basis_set_single), intent(out) :: abas(nabas)
        logical, intent(in) :: fix_gto_norm
        character, parameter :: orborder(6) = ['s', 'p', 'd', 'f', 'g', 'h']
        integer, parameter :: maxnfunc = 100
        integer, parameter :: maxnprim = 25
        integer :: nfunc
        integer :: nprim(maxnfunc)
        real(dp) :: zeta(maxnfunc, maxnprim)
        real(dp) :: beta(maxnfunc, maxnprim)
        character :: orbtype(maxnfunc)
        integer :: i, j
       
        do i = 1, nabas
            call readf%next()
            read(readf%line, *) abas_index(i)
            nfunc = 0
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
                if (fix_gto_norm) then
                    beta(nfunc, :) = beta(nfunc, :) / gto_norms(zeta(nfunc, :), &
                    &                                           amp%l(orbtype(nfunc)), 0, 0)
                end if
            end do
            call abas(abas_index(i))%init(nfunc, orbtype, nprim, zeta, beta)
        end do
    end subroutine section_read_gto


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_count_mo
    !> @brief Count number of AOs/MOs in a molden [MO] section.
    !----------------------------------------------------------------------------------------------
    subroutine section_count_mo(readf, nmo_a, nmo_b, nbas)
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
        nmo_a = readf%count_keyword_appearances('Alpha')
        nmo_b = readf%count_keyword_appearances('Beta')
    end subroutine section_count_mo


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_mo
    !> @brief Read the molden [MO] section.
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


end module molden_sections_mod
