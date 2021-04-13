!----------------------------------------------------------------------------------------------
! MODULE: turbomole_sections_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @brief Subroutines for reading data from sections present in turbomole input/output files.
!----------------------------------------------------------------------------------------------
module turbomole_sections_mod
    use global_defs
    use file_mod
    implicit none


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_rundim
    !
    !> @brief Read the turbomole $rundimensions section.
    !> @details
    !! Contains number of atoms and basis functions, and info on (un)restricted calculation.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_rundim(readf, natom, ncart, nsphe, rhf)
        type(reader), intent(inout) :: readf
        integer, intent(out) :: natom
        integer, intent(out) :: ncart
        integer, intent(out) :: nsphe
        integer, intent(out) :: rhf

        rhf = 0
        nsphe = 0
        ncart = 0
        natom = 0
        do
            call readf%next()
            if (index(readf%line, '$') /= 0) return
            call readf%parseline('=')
            select case (readf%args(1)%s)
                case('natoms')
                    read(readf%args(2)%s,*) natom
                case('nbf(CAO)')
                    read(readf%args(2)%s,*) ncart
                case('nbf(AO)')
                    read(readf%args(2)%s,*) nsphe
                case('rhfshells')
                    read(readf%args(2)%s,*) rhf
            end select
            cycle
        end do
    end subroutine section_read_rundim


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_coord
    !> @brief Read the geometry and atom symbol list from a turbomole $coord section.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_coord(readf, natom, asym, geom)
        type(reader), intent(inout) :: readf
        integer, intent(in) :: natom
        character(len=2), intent(out) :: asym(natom)
        real(dp), intent(out) :: geom(3*natom)
        integer :: i

        do i = 1, 3*natom, 3
            call readf%next()
            read(readf%line, *) geom(i:i+2), asym((i+2)/3)
        end do
    end subroutine section_read_coord


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_atoms
    !
    !> @brief Read the turbomole $atoms section.
    !> @details
    !! The $atoms section contains list of basis sets and indices of atoms corresponding to each
    !! type of basis set. Also contains aux. basis set information, but is not read here.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_atoms(readf, nabas, basis_asym, basis_index, basis_key)
        use ivec_mod
        use string_mod
        type(reader), intent(inout) :: readf
        integer, intent(in) :: nabas
        character(len=2), intent(out) :: basis_asym(nabas)
        type(ivec), intent(out) :: basis_index(nabas)
        type(string), intent(out) :: basis_key(nabas)
        integer :: i
        integer :: j
        character(len=:), allocatable :: indexstr
        
        do i = 1, nabas
           indexstr = ''
           call readf%next()
           call readf%parseline('=\ ')
           read(readf%args(1)%s, '(a)') basis_asym(i)
           do j = 2, readf%narg
               if (char_is_num(readf%args(j)%s(1:1))) then
                   indexstr = indexstr//' '//readf%args(j)%s
               else
                   exit
               end if
           end do
           call read_index_list(indexstr, basis_index(i)%c)
           do j = 3, readf%narg - 1
               if (readf%args(j)%s == 'basis') then
                   basis_key(i)%s = readf%args(j+1)%s//' '//readf%args(j+2)%s
               end if
           end do
        end do
    end subroutine section_read_atoms


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_shells
    !
    !> @brief Read the occupation numbers of molecular orbitals from a Turbomole 'control' file.
    !> @details
    !! The occupation numbers are written in the "$closed shells" section of the control file for 
    !! restricted calculations, and "$alpha shells" and "$beta shells" sections for unrestricted
    !! calculations.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_shells(readf, nmo, omask)
        use string_mod
        type(reader), intent(inout) :: readf
        integer, intent(in) :: nmo
        logical, intent(out) :: omask(nmo)
        character(len=2) :: sym
        integer, allocatable :: indexlist(:)
        integer :: ipos

        call readf%next()
        call readf%parseline('()')
        read(readf%args(1)%s, *) sym
        ipos = index(readf%args(1)%s, sym) + 2
        call read_index_list(readf%args(1)%s(ipos:), indexlist)
        omask(indexlist) = .true.
    end subroutine section_read_shells


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_freeze
    !
    !> @brief Read the turbomole $freeze section.
    !> @details
    !! Frozen orbitals can be defined as implicit (lowest core or highest virtual) or as a list
    !! of frozen orbitals.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_freeze(readf, nmo, amask)
        use string_mod
        type(reader), intent(inout) :: readf
        integer, intent(in) :: nmo
        logical, intent(out) :: amask(nmo)

        integer :: i
        integer :: num
        integer, allocatable :: indexlist(:)

        call readf%next()
        if (index(readf%line, 'implicit') /= 0) then
            call readf%parseline(' =')
            do i = 1, readf%narg - 1
                if (trim(readf%args(i)%s) == 'core') then
                    read(readf%args(i+1)%s, *) num
                    amask(1:num) = .false.
                else if (trim(readf%args(i)%s) == 'virt') then
                    read(readf%args(i+1)%s, *) num
                    num = nmo - num + 1
                    amask(num:) = .false.
                end if
            end do
        else
            readf%line = readf%line(index(readf%line, ' '):len(readf%line))
            call read_index_list(readf%line, indexlist)
            do i = 1, size(indexlist)
                amask(indexlist(i)) = .false.
            end do
        end if
    end subroutine section_read_freeze


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_basis
    !
    !> @brief Read all basis sets from a turbomole $basis section
    !> @details
    !! Coefficients of contracted Gaussian functions are read from the $basis section, sorted in
    !! order of increasing angular momentum and stored into an array of atombasis type.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_basis(readf, nabas, abas)
        use basis_set_mod, only : basis_set_single
        type(reader), intent(inout) :: readf
        integer, intent(in) :: nabas
        type(basis_set_single), intent(out) :: abas(nabas)
       
        integer :: i
        integer :: j
        integer :: l
        integer, parameter :: maxnfunc = 100
        integer, parameter :: maxnprim = 25
        integer :: nfunc
        integer :: ncart
        integer :: nsphe
        integer :: nprim(maxnfunc)
        real(dp) :: zeta(maxnfunc, maxnprim)
        real(dp) :: beta(maxnfunc, maxnprim)
        character :: orbtype(maxnfunc)
        integer :: cfunc
       
        call readf%next() ! *
        do i = 1, nabas
            nfunc = 0
            ncart = 0
            nsphe = 0
            cfunc = 0
            call readf%next() 
            abas(i)%key = trim(adjustl(readf%line))
            call readf%next() ! *
            do
                call readf%next()
                if (index(readf%line, '*') == 1) exit
                nfunc = nfunc + 1
                read(readf%line, *) nprim(nfunc), orbtype(nfunc)
                do j = 1, nprim(nfunc)
                    call readf%next()
                    read(readf%line, *) zeta(nfunc,j), beta(nfunc,j)
                end do
            end do
            call abas(i)%init(nfunc, orbtype, nprim, zeta, beta, .true.)
        end do
    end subroutine section_read_basis


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_mo
    !> @brief Read all MOs and MO energies from a $uhfmo or $rhfmo section.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_mo(readf, nmo, nbas, mo_c, mo_e)
        type(reader), intent(inout) :: readf
        integer, intent(in) :: nmo !< Number of molecular orbitals.
        integer, intent(in) :: nbas !< Number of basis functions.
        real(dp), intent(inout) :: mo_c(nbas, nmo) !< Molecular orbital coefficients.
        real(dp), intent(inout) :: mo_e(nmo) !< Molecular orbital energies.
        integer :: i
        integer :: j
        character(len=:), allocatable :: readformat
        integer :: ninline
        integer :: cmo


        if (index(readf%line, 'format') /= 0) then
            call readf%parseline('()') ! Extract the format string.
            readformat = '('//readf%args(2)%s//')'
            ! To get the number of orbitals per line, the string is split after the 
            ! first occurence of one of the fortran format specifiers for reals.
            readf%line = readf%args(2)%s
            call readf%parseline('def')
            read(readf%args(1)%s, *) ninline
        end if
        do i = 1, nmo
            call readf%next()
            if (index(readf%line, '$') == 1) then
                write(stderr, *) 'Error in turbomole_sections_mod.'
                write(stderr, *) '  Not all molecular orbitals found.'
                stop
            end if
            call readf%parseline(' =')
            read(readf%args(1)%s, *) cmo
            read(readf%args(4)%s, *) mo_e(cmo)
            do j = 1, nbas, ninline
               call readf%next()
               read(readf%line, readformat) mo_c(j:min(nbas,j+ninline-1), cmo)
            end do
        end do
    end subroutine section_read_mo


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_eigenpairs
    !> @brief Read the eigenvectors and eigenvalues from the $eigenpairs section.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_eigenpairs(readf, nstate, edim, eigval, eigvec)
        type(reader), intent(inout) :: readf
        integer, intent(in) :: nstate
        integer, intent(in) :: edim
        real(dp), intent(out) :: eigval(nstate)
        real(dp), intent(out) :: eigvec(edim, nstate)
        character(len=9), parameter :: linefmt = '(4d20.14)'
        integer :: i
        integer :: cstate

        do i = 1, nstate
            call readf%next()
            call readf%parseline(' =')
            read(readf%args(1)%s, *) cstate
            read(readf%args(3)%s, *) eigval(cstate)
            read(readf%unit, linefmt) eigvec(:, cstate)
        end do
    end subroutine section_read_eigenpairs


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_energy
    !> @brief Read the full turbomole $energy section
    !----------------------------------------------------------------------------------------------
    subroutine section_read_energy(readf, nener, scf, kin, pot, ene)
        type(reader), intent(inout) :: readf
        integer, intent(in) :: nener
        real(dp), intent(out) :: scf(nener)
        real(dp), intent(out) :: kin(nener)
        real(dp), intent(out) :: pot(nener)
        real(dp), intent(out) :: ene(nener)
        integer :: i
        integer :: itmp

        do i = 1, nener
            call readf%next()
            read(readf%line, *) itmp, scf(i), kin(i), pot(i), ene(i)
            !> @todo some calculations don't have final column.
        end do
    end subroutine section_read_energy


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_spectrum
    !> @brief Read the spectrum file.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_spectrum(readf, nstate, ene, osc)
        type(reader), intent(inout) ::readf
        integer, intent(in) :: nstate
        real(dp), intent(out) :: ene(nstate)
        real(dp), intent(out) :: osc(nstate)
        integer :: i

        do i = 1, nstate
            call readf%next()
            read(readf%line, *) ene(i), osc(i)
        end do
    end subroutine section_read_spectrum
 

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_grad_last
    !> @brief Read the last gradient from a turbomole $grad section.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_grad_last(readf, natom, grad)
        type(reader), intent(inout) :: readf
        integer, intent(in) :: natom
        real(dp), intent(out) :: grad(3, natom)
        integer :: i

        do
            call readf%next()
            if (index(readf%line,'$') == 1) return
            do i = 1, natom
                call readf%next()
            end do
            do i = 1, natom
                call readf%next()
                read(readf%line, *) grad(:, i)
            end do
        end do
    end subroutine section_read_grad_last


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_grad
    !> @brief Read the full turbomole $grad section.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_grad(readf, ncycle, dof, all_geom, all_grad, all_en)
        type(reader), intent(inout) :: readf
        integer, intent(in) :: ncycle
        integer, intent(in) :: dof
        real(dp), intent(out) :: all_geom(dof, ncycle)
        real(dp), intent(out) :: all_grad(dof, ncycle)
        real(dp), intent(out) :: all_en(ncycle) !< @todo read the energy.
        integer :: i, c

        do c = 1, ncycle
            call readf%next()
            do i = 1, dof, 3
                call readf%next()
                read(readf%line, *) all_geom(i:i+2, c)
            end do
            do i = 1, dof, 3
                call readf%next()
                read(readf%line, *) all_grad(i:i+2, c)
            end do
        end do
    end subroutine section_read_grad
 

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_hess
    !> @brief Read hessian.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_hess(readf, dof, hess)
        type(reader), intent(inout) :: readf
        integer, intent(in) :: dof
        real(dp), intent(inout) :: hess(dof, dof)
        integer, parameter :: ninline = 5
        integer :: i, j, du

        do i = 1, dof
            do j = 1, dof, ninline
                call readf%next()
                read(readf%line, *) du, du, hess(i, j:min(dof, j+ninline-1))
            end do
        end do
    end subroutine section_read_hess


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_read_dft 
    !> @brief Read selected DFT functional.
    !----------------------------------------------------------------------------------------------
    subroutine section_read_dft(readf, dft_functional)
        type(reader), intent(inout) :: readf
        character(len=:), allocatable, intent(out) :: dft_functional

        do
            call readf%next()
            if (index(readf%line, '$') == 1) exit
            call readf%parseline(' ')
            if (readf%args(1)%s == 'functional') then
                dft_functional = trim(adjustl(readf%args(2)%s))
            end if
        end do
    end subroutine section_read_dft


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: section_count_basis
    !> @brief Get number of basis sets in a $basis section.
    !----------------------------------------------------------------------------------------------
    function section_count_basis(readf) result (nabas)
        type(reader), intent(inout) :: readf
        integer :: nabas
        logical :: check
        integer :: iline
        iline = readf%line_num
        nabas = 0
        do
            call readf%go_to_keyword('*', found=check)
            if (.not. check) exit
            nabas = nabas + 1
        end do
        call readf%go_to_line(iline)
        nabas = (nabas - 1) / 2
    end function section_count_basis


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: section_count_mo
    !> @brief Get dimensions of MO array.
    !----------------------------------------------------------------------------------------------
    subroutine section_count_mo(readf, nbas, nmo)
        type(reader) :: readf
        integer, intent(out) :: nbas
        integer, intent(out) :: nmo
        logical :: check
        integer :: iline
        integer :: i
        iline = readf%line_num
        nmo = 0
        do 
            call readf%go_to_keyword('eigenvalue=', found=check)
            if (.not. check) exit
            nmo = nmo + 1
            if (nmo == 1) then
                call readf%parseline(' =')
                do i = 1, readf%narg
                    if (readf%args(i)%s == 'nsaos') then
                        read(readf%args(i+1)%s, *) nbas
                        exit
                    end if
                end do
            end if
        end do
        call readf%go_to_line(iline)
    end subroutine section_count_mo


end module turbomole_sections_mod
