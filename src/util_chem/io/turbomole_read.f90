!----------------------------------------------------------------------------------------------
! MODULE: turbomole_read_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @brief Read data turbomole input/output files.
!----------------------------------------------------------------------------------------------
module turbomole_read_mod
    use global_defs
    use file_mod
    use turbomole_sections_mod
    use turbomole_definitions_mod
    implicit none


    private
    public :: turbomole_read_geom
    public :: turbomole_read_basis
    public :: turbomole_read_mo
    public :: turbomole_read_cis


    ! Helper variables while searching through a turbomole directory.
    type(reader) :: tm_reader
    character(len=:), allocatable :: tm_dname
    character(len=:), allocatable :: tm_fname
    character(len=1000) :: tm_initdir


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: turbomole_read_geom
    !> @brief Read geometry (and atom symbols) from a turbomole format file.
    !----------------------------------------------------------------------------------------------
    subroutine turbomole_read_geom(path, geom, at_symbol)
        character(len=*), intent(in) :: path
        real(dp), allocatable :: geom(:) !< Geometry.
        character(len=2), allocatable, optional :: at_symbol(:) !< Atomic symbols.
        character(len=2), allocatable :: tsymbol(:) !< Temp. atomic symbols.
        integer :: natom

        call tm_enter_dir(path)
        call tm_find_section('$coord')
        natom = tm_reader%count_lines_to_keyword('$')
        allocate(geom(3*natom))
        allocate(tsymbol(natom))
        call section_read_coord(tm_reader, natom, tsymbol, geom)
        call tm_reader%close()
        call chdir(tm_initdir)
        if (present(at_symbol)) at_symbol=tsymbol
    end subroutine turbomole_read_geom


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: turbomole_read_basis
    !> @brief Read all information about the basis set from a turbomole calculation.
    !----------------------------------------------------------------------------------------------
    subroutine turbomole_read_basis(path, bs, read_atoms)
        use basis_set_mod, only : basis_set
        use ivec_mod
        use string_mod
        character(len=*), intent(in) :: path
        type(basis_set), intent(out) :: bs
        logical, intent(in), optional :: read_atoms
        integer :: nabas
        type(string), allocatable :: key(:) !< Basis set keys from atoms section.
        character(len=2), allocatable :: key_asym(:) !< Atomic symbols corresponding to keys.
        type(ivec), allocatable :: key_atom_index(:) !< Atom indexes corresponding to keys.
        logical :: tread_atoms
        integer :: i, n

        bs%sphe_mo = .true.
        bs%source_format = 'turbomole'
        tread_atoms = .true.
        if (present(read_atoms)) tread_atoms = read_atoms

        call tm_enter_dir(path)
        call tm_find_section('$basis')
        bs%n_bs = section_count_basis(tm_reader)
        allocate(bs%bs(bs%n_bs))
        call section_read_basis(tm_reader, bs%n_bs, bs%bs)
        call tm_reader%close()
        if (tread_atoms) then
            call tm_find_section('$atoms')
            nabas = tm_reader%count_lines_to_keyword('$')
            allocate(key_asym(nabas))
            allocate(key_atom_index(nabas))
            allocate(key(nabas))
            call section_read_atoms(tm_reader, nabas, key_asym, key_atom_index, key)
            bs%n_center = maxval(key_atom_index)
            if (allocated(bs%center_i_bs)) deallocate(bs%center_i_bs)
            allocate(bs%center_i_bs(bs%n_center))
            call tm_match_basis_with_atoms(bs%bs, key, key_atom_index, bs%center_i_bs)
            call tm_reader%close()
            do i = 1, bs%n_center
                bs%n_bf_sphe = bs%n_bf_sphe + bs%bs(bs%center_i_bs(i))%n_sphe
                bs%n_bf_cart = bs%n_bf_cart + bs%bs(bs%center_i_bs(i))%n_cart
            end do
        end if
        call chdir(tm_initdir)
    end subroutine turbomole_read_basis


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: turbomole_read_mo
    !> @brief Read all information about molecular orbitals from a turbomole calculation.
    !----------------------------------------------------------------------------------------------
    subroutine turbomole_read_mo(path, mos)
        use molecular_orbitals_mod, only : molecular_orbitals
        character(len=*), intent(in) :: path
        type(molecular_orbitals) :: mos 
        logical :: check
        integer :: nbas, nmo
        real(dp), allocatable :: tmoa_c(:, :)
        real(dp), allocatable :: tmoa_e(:)
        real(dp), allocatable :: tmob_c(:, :)
        real(dp), allocatable :: tmob_e(:)

        call tm_enter_dir(path)
        call tm_find_section('$scfmo', check)
        if (.not. check) then
            call tm_reader%close()
            call tm_find_section('$uhfmo_alpha', check)
        end if
        if (check) then
            call section_count_mo(tm_reader, nbas, nmo)
            allocate(tmoa_c(nbas, nmo))
            allocate(tmoa_e(nmo))
            call section_read_mo(tm_reader, nmo, nbas, tmoa_c, tmoa_e)
        end if
        call tm_reader%close()
        call tm_find_section('$uhfmo_beta', check)
        if (check) then
            call section_count_mo(tm_reader, nbas, nmo)
            allocate(tmob_c(nbas, nmo))
            allocate(tmob_e(nmo))
            call section_read_mo(tm_reader, nmo, nbas, tmob_c, tmob_e)
        end if
        call tm_reader%close()
        call chdir(tm_initdir)

        call mos%init(ca=tmoa_c, ea=tmoa_e, cb=tmob_c, eb=tmob_e)
    end subroutine turbomole_read_mo


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: turbomole_read_cis
    !> @brief Read all information about CIS vectors from escf or ricc2 calculation.
    !----------------------------------------------------------------------------------------------
    subroutine turbomole_read_cis(path, wfa, wfb, occ_mo, act_mo)
        use occupation_mod
        use basis_set_mod, only : basis_set
        character(len=*), intent(in) :: path
        real(dp), allocatable :: wfa(:, :)
        real(dp), allocatable :: wfb(:, :)
        logical, allocatable :: occ_mo(:, :)
        logical, allocatable :: act_mo(:, :)
        type(occupation_numbers) :: onum
        character(len=:), allocatable :: method
        type(basis_set) :: bs
        integer :: natom
        integer :: ncart
        integer :: rhf
        integer :: nex
        logical :: check

        call tm_enter_dir(path)
        call tm_find_section('$rundimensions', check)
        natom = 0
        ncart = 0
        onum%n = 0
        rhf = 0
        if (check) then
            call section_read_rundim(tm_reader, natom, ncart, onum%n, rhf)
        end if
        if (rhf == 0) then
            call tm_find_section('$closed shells', check)
            if (check) then
                rhf = 1
            else
                call tm_find_section('$alpha shells', check)
                if (check) rhf = 2
            end if
        end if
        if (rhf == 0) then
            write(stderr, *) 'Error in turbomole_read module,&
                            & turbomole_read_cis subroutine.'
            write(stderr, *) ' No rhf information'
            write(stderr, *) ' Looking for: rhfshells/closed shells/alpha shells.'
            call abort()
        end if
        if ((ncart == 0) .or. (onum%n == 0) .or. (natom == 0)) then
            call turbomole_read_basis(path, bs, .true.)
            natom = bs%n_center
            ncart = bs%n_bf_cart
            onum%n = bs%n_bf_sphe
        end if

        call tm_reader%close()
        allocate(occ_mo(onum%n, 2), source=.false.)
        allocate(act_mo(onum%n, 2), source=.true.)
        if (rhf == 1) then
            call tm_find_section('$closed shells')
            call section_read_shells(tm_reader, onum%n, occ_mo(:, 1))
            occ_mo(:, 2) = occ_mo(:, 1)
            call tm_reader%close()
        else
            call tm_find_section('$alpha shells')
            call section_read_shells(tm_reader, onum%n, occ_mo(:, 1))
            call tm_reader%close()
            call tm_find_section('$beta shells')
            call section_read_shells(tm_reader, onum%n, occ_mo(:, 2))
            call tm_reader%close()
        end if
        call tm_find_section('$freeze', check)
        if (check) then
            call section_read_freeze(tm_reader, onum%n, act_mo(:, 1))
            act_mo(:, 2) = act_mo(:, 1)
        end if
        call tm_reader%close()

        call occ_mask_to_nums(rhf, occ_mo, act_mo, onum)
        call tm_find_section('ricc2', check)
        call tm_reader%close()
        if (check) then
            method = 'ricc2'
        else
            method = 'escf'
        end if
        nex = tm_count_n_ex_states(method)

        if (allocated(wfa)) deallocate(wfa)
        if (allocated(wfb)) deallocate(wfb)
        allocate(wfa(onum%ao(1)*onum%av(1), nex))
        if (rhf == 2) allocate(wfb(onum%ao(2)*onum%av(2), nex))

        select case(method)
        case('ricc2')
            call read_cis_ricc2(rhf, wfa, wfb)
        case('escf')
            call read_cis_escf(rhf, wfa, wfb)
        end select
        call chdir(tm_initdir)
    end subroutine turbomole_read_cis


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_cis_ricc2
    !
    !> @brief Read the CIS coefficients after an ricc2 calculation.
    !> @details
    !! The coefficients are read from unformmated CCRE0... files.
    !----------------------------------------------------------------------------------------------
    subroutine read_cis_ricc2(rhf, wfa, wfb)
        integer, intent(in) :: rhf
        real(dp) :: wfa(:, :)
        real(dp) :: wfb(:, :)
        integer, parameter :: lenstk = 4
        character(len=lenstk) :: statekey
        character(len=:), allocatable :: inputfile
        character(len=3) :: dummy
        integer :: i, inunit

        do i = 1, size(wfa, 2)
            statekey = '----'
            write(statekey(lenstk-int(log10(real(i))):lenstk), '(i0)') i
            inputfile = 'CCRE0-1--1'//statekey
            call need_file(inputfile, 'Error in read_cis_ricc2 subroutine.')
            open(newunit=inunit, file=inputfile, form='UNFORMATTED', action='READ')
                read(unit=inunit) dummy(1:3)
                read(unit=inunit) wfa(:, i)
                if (rhf == 2) read(unit=inunit) wfb(:, i)
            close(inunit)
        end do
    end subroutine read_cis_ricc2


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_cis_escf
    !
    !> @brief Read the CIS coefficients after an escf calculation.
    !> @details
    !! Read type of functional and type of calculation (regular or TDA) to determine dimensions of
    !! eigenvectors and read eigenvectors from appropriate file.
    !----------------------------------------------------------------------------------------------
    subroutine read_cis_escf(rhf, wfa, wfb)
        integer, intent(in) :: rhf
        real(dp) :: wfa(:, :)
        real(dp) :: wfb(:, :)
        real(dp), allocatable :: eigval(:)
        real(dp), allocatable :: eigvec(:, :)
        character(len=:), allocatable :: dft_func, instab
        integer :: dima, dimb, edim, nex, supdim
        integer :: i


        call tm_find_section('$dft')
        call section_read_dft(tm_reader, dft_func)
        call tm_reader%close()
        call tm_find_section('$scfinstab')
        call tm_reader%parseline(' ')
        instab = trim(adjustl(tm_reader%args(2)%s))
        call tm_reader%close()

        supdim = 2
        if ((instab == 'ucis') .or. (instab == 'ciss')) supdim = 1
        if (turbomole_functional_type(dft_func) < 3) supdim = 1
        dima = size(wfa, 1)
        dimb = 0
        if (rhf == 2) dimb = size(wfb, 1)
        nex = size(wfa, 2)
        edim = (dima + dimb) * supdim
        allocate(eigvec(edim, nex))
        allocate(eigval(nex))

        call tm_find_eigenpairs(instab)
        call section_read_eigenpairs(tm_reader, nex, edim, eigval, eigvec)
        call tm_reader%close()
        do i = 1, nex
            if (rhf == 2) then
                call tm_eigvec_to_wf(1, supdim, rhf, eigvec(:, i), wfa(:, i), wfb(:, i))
            else
                call tm_eigvec_to_wf(1, supdim, rhf, eigvec(:, i), wfa(:, i), wfa(:, i))
            end if
        end do
    end subroutine read_cis_escf


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: tm_match_basis_with atoms
    !
    !> @brief Match information from $atoms and $basis sections.
    !> @details
    !! Take list of indices of atoms corresponding to each basis set and create list with index of
    !! of basis set corresponding to each atom.
    !----------------------------------------------------------------------------------------------
    subroutine tm_match_basis_with_atoms(basis, key, key_atom_index, atom_basis_index)
        use basis_set_mod, only : basis_set_single
        use ivec_mod
        use string_mod
        type(basis_set_single), intent(in) :: basis(:)
        type(string), intent(in) :: key(:)
        type(ivec), intent(in) :: key_atom_index(:)
        integer, intent(out) :: atom_basis_index(:)
        integer :: i, j, k, chk_natom

        chk_natom = 0
        do i = 1, size(key)
            do j = 1, size(basis)
                if (key(i)%s == basis(j)%key) then
                    do k = 1, size(key_atom_index(i)%c)
                        atom_basis_index(key_atom_index(i)%c(k)) = j
                        chk_natom = chk_natom + 1
                    end do
                    exit
                else
                    if (j == size(basis)) then
                        write(stderr, *) 'Error in turbomole_read module,&
                                        & tm_match_basis_with_atoms subroutine.'
                        write(stderr, *) key(i)%s//' not found in basis section.'
                        call abort()
                    end if
                end if
            end do
        end do
        if (chk_natom /= size(atom_basis_index)) then
            write(stderr, *) 'Error in turbomole_read module, tm_match_basis_with_atoms subroutine.'
            write(stderr, '(a,i0,a)') ' Found ', chk_natom, ' atoms.'
            write(stderr, '(a,i0,a)') ' Expected ', size(atom_basis_index), 'atoms.'
        end if
    end subroutine tm_match_basis_with_atoms


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: tm_eigvec_to_wf
    !
    !> @brief Change format of TDDFT eigenvector coefficients and split alpha and beta parts.
    !> @details
    !! The format in which the output is given is determined by the outopt variable:
    !!    0 - X or X+Y        - Dimension: no * nv.
    !!    1 - sqrt(X^2 - Y^2) - Single electron excitation "coefficients". Dimension: no * nv.
    !!    2 - X+Y and X-Y     - Format of input files. Dimension: 2 * no * nv
    !! When dimension of super-tensorspace is 1 (LDA and GGA functionals or TDA approximation),
    !! output option 0 is always used.
    !----------------------------------------------------------------------------------------------
    subroutine tm_eigvec_to_wf(outopt, supdim, rhf, ev, wfa, wfb)
        integer, intent(in) :: outopt
        integer, intent(in) :: supdim
        integer, intent(in) :: rhf
        real(dp), intent(in) :: ev(:)
        real(dp), intent(out) :: wfa(:)
        real(dp), intent(out) :: wfb(:)
        integer :: da
        integer :: db
        integer :: dt
        dt = size(ev)
        da = size(wfa)
        db = 0
        if (rhf == 2) db = size(wfb)
        if (outopt == 2) then
            da = da / supdim
            db = db / supdim
        end if
        if (supdim == 1) then
            wfa(:) = ev(1:da)
            if (rhf == 2) wfb(:) = ev(da+1:dt)
        else
            select case(outopt)
            case(0)
                wfa(:) = ev(1:da)
                if (rhf == 2) wfb(:) = ev(da+1:da+db)
            case(1)
                wfa(:) = sqrt(abs(ev(1:da)*ev(da+db+1:2*da+db)))
                wfa(:) = sign(wfa, -ev(1:da))
                if (rhf == 2) then
                  wfb(:) = sqrt(abs(ev(da+1:da+db)*ev(2*da+db+1:dt)))
                  wfb(:) = sign(wfb, -ev(da+1:da+db))
                end if
            case(2)
                wfa(1:da) = ev(1:da)
                wfa(da+1:2*da) = ev(da+db+1:2*da+db)
                if (rhf == 2) then
                    wfb(1:db) = ev(da+1:da+db)
                    wfb(db+1:2*db) = ev(2*da+db+1:dt)
                end if
            end select
        end if
    end subroutine tm_eigvec_to_wf


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: tm_enter_dir
    !
    !> @brief Enter directory containing turbomole files.
    !> @details
    !! If given path is a directory enter the directory and set main file name to 'control'. If 
    !! not, work in current directory and start from the file given in path.
    !----------------------------------------------------------------------------------------------
    subroutine tm_enter_dir(path)
        character(len=*), intent(in) :: path

        if (check_is_dir(path//'/')) then
            tm_dname = path
            tm_fname = 'control'
        else
            tm_dname = '.'
            tm_fname = path
        end if
        call getcwd(tm_initdir)
        call chdir(tm_dname)
    end subroutine tm_enter_dir


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: tm_find_section
    !
    !> @brief Go to a turbomole section.
    !> @details
    !! Finds a section in the tm_fname file. If the section is redirected to a different file
    !! (as in: $section file=fname), open the new file and find the section in that file.
    !----------------------------------------------------------------------------------------------
    subroutine tm_find_section(section, found)
        character(len=*), intent(in) :: section
        logical, intent(out), optional :: found
        integer :: i
        logical :: check
        
        call tm_reader%open(tm_fname, comment='#', continuation='\')
        call tm_reader%go_to_keyword(section, found=check)
        if (.not. check) then
            if (.not. present(found)) then
                write(stderr, *) 'Error in tm_find_section.'
                write(stderr, *) 'Section "'//section//'" not found.'
                call abort()
            else
                found = .false.
                return
            end if
        else
            if (present(found)) found = .true.
        end if
        call tm_reader%parseline(' =')
        do i = 2, tm_reader%narg - 1
            if (tm_reader%args(i)%s == 'file') then
                call tm_reader%close()
                call tm_reader%open(tm_reader%args(i+1)%s)
                call tm_reader%go_to_keyword(tm_reader%args(1)%s)
            end if
        end do
    end subroutine tm_find_section


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: tm_find_eigenpairs
    !
    !> @brief Find eigenpairs section depending on type of calculation.
    !> @details
    !! Eigenpairs section is not mentioned in control file, the default name of the file containing
    !! the section depends on the scfinstab option in the control file.
    !----------------------------------------------------------------------------------------------
    subroutine tm_find_eigenpairs(instab)
        character(len=*), intent(in) :: instab
        logical :: check
        call tm_find_section('$eigenpairs', check)
        if (.not. check) then
            call tm_reader%close()
            select case(instab)
            case('rpas')
                call tm_reader%open('sing_a')
            case('ciss')
                call tm_reader%open('ciss_a')
            case('urpa')
                call tm_reader%open('unrs_a')
            case('ucis')
                call tm_reader%open('ucis_a')
            end select
            call tm_reader%go_to_keyword('$eigenpairs')
        end if
    end subroutine tm_find_eigenpairs


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: tm_count_n_ex_states
    !> @brief Find number of excited states included in calculation (depending on method).
    !----------------------------------------------------------------------------------------------
    function tm_count_n_ex_states(method) result(nex)
        character(len=*) :: method
        integer :: nex
        integer :: i
        select case(method)
        case('ricc2')
            call tm_find_section('$excitations')
            outer: do
                call tm_reader%next()
                if (index(tm_reader%line, '$') == 1) then
                    write(stderr, *) 'Error in turbomole_read_mod.'
                    write(stderr, *) 'irrep line not found in $excitations section.'
                    call abort()
                end if
                if (index(tm_reader%line, 'irrep') /= 0) then
                    call tm_reader%parseline(' =')
                    do i = 1, tm_reader%narg - 1
                        if (tm_reader%args(i)%s == 'nexc') then
                            read(tm_reader%args(i+1)%s, *) nex
                            exit outer
                        end if
                    end do
                end if
            end do outer
            call tm_reader%close()
        case('escf')
            call tm_find_section('$soes')
            call tm_reader%next()
            call tm_reader%parseline(' ')
            read(tm_reader%args(2)%s, *) nex
            call tm_reader%close()
        case default
            write(stderr, *) 'Error in turbomole_read_mod subroutine.'
            write(stderr, *) 'Unrecognized method '//method//' while checking number of states.'
        end select
    end function tm_count_n_ex_states


end module turbomole_read_mod
