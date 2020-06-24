!----------------------------------------------------------------------------------------------
! MODULE: read_all
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date June, 2017
!
!> @brief Contains subroutines for reading data from a variety of quantum chemistry formats.
!> @details
!----------------------------------------------------------------------------------------------
module read_all_mod
    use global_defs
    use molden_read_mod
    use turbomole_read_mod
    use molpro_output_read_mod
    implicit none


    private
    public :: read_geom
    public :: read_basis
    public :: read_mo
    public :: read_cis


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_geom
    !> @brief Read geometry and list of atoms.
    !----------------------------------------------------------------------------------------------
    subroutine read_geom(dformat, path, geom, atom_symbol, atom_number)
        use element_mod
        character(len=*), intent(in) :: dformat
        character(len=*), intent(in) :: path
        real(dp), allocatable :: geom(:)
        character(len=2), allocatable, optional :: atom_symbol(:)
        integer, allocatable, optional :: atom_number(:)
        character(len=2), allocatable :: tsymbol(:)
        integer, allocatable :: tnumber(:)

        select case(dformat(:6))
        case('turbom')
            call turbomole_read_geom(path, geom, tsymbol)
            tnumber = element_s2z(tsymbol)
        case('molden')
            call molden_read_geom(path, geom, tsymbol, tnumber)
        case('molpro')
            call molpro_output_read_geom(path, geom, tsymbol, tnumber)
        case default
            write(stderr, *) 'Error in read_geom subroutine.'
            write(stderr, *) 'Input format not implemented: '//trim(adjustl(dformat))
        end select

        if (present(atom_symbol)) atom_symbol = tsymbol
        if (present(atom_number)) atom_number = tnumber
    end subroutine read_geom


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_basis
    !> @brief Read basis set information.
    !----------------------------------------------------------------------------------------------
    subroutine read_basis(dformat, path, basis)
        use basis_set_mod, only : basis_set
        character(len=*), intent(in) :: dformat
        character(len=*), intent(in) :: path
        type(basis_set), intent(inout) :: basis
        
        select case(dformat(:6))
        case('turbom')
            call turbomole_read_basis(path, basis)
        case('molden') 
            select case(dformat)
            case('molden_orca')
                call molden_read_basis(path, basis, 3)
            case default
                call molden_read_basis(path, basis, 1)
            end select
        case('molpro')
            call molpro_output_read_basis(path, basis)
        case default
            write(stderr, *) 'Error in read_basis subroutine.'
            write(stderr, *) 'Input format not implemented: '//trim(adjustl(dformat))
            stop
        end select

        call basis%check_init()
    end subroutine read_basis


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_mo
    !> @brief Read molecular orbitals.
    !----------------------------------------------------------------------------------------------
    subroutine read_mo(dformat, path, mos, bs)
        use molecular_orbitals_mod, only : molecular_orbitals
        use basis_set_mod, only : basis_set
        character(len=*), intent(in) :: dformat
        character(len=*), intent(in) :: path
        type(molecular_orbitals) :: mos
        type(basis_set), intent(in), optional :: bs

        select case(dformat(:6))
        case('turbom')
            call turbomole_read_mo(path, mos)
        case('molden') 
            select case(dformat)
            case('molden_tm2molden')
                if (.not. present(bs)) then
                    write(stderr, *) 'Error, read_mo called for tm2molden without passing bs.'
                    stop
                end if
                call molden_read_mo(path, mos, bs, 1)
            case('molden_cfour')
                if (.not. present(bs)) then
                    write(stderr, *) 'Error, read_mo called for cfour without passing bs.'
                    stop
                end if
                call molden_read_mo(path, mos, bs, 2)
            case('molden_orca')
                call molden_read_mo(path, mos, bs, 3)
            case default
                call molden_read_mo(path, mos)
            end select
        case('molpro_output')
            call molpro_output_read_mo(path, mos)
        case default
            write(stderr, *) 'Error in read_mo subroutine.'
            write(stderr, *) 'Input format not implemented: '//trim(adjustl(dformat))
        end select

        !< @todo add check for proper mos initialization
    end subroutine read_mo


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_cis
    !> @brief Read CIS wave function coefficients.
    !> @todo Add CIS type to hold relevant arrays.
    !> @todo Code cleaning
    !> @todo Orthogonalization for unrestricted wave functions?
    !----------------------------------------------------------------------------------------------
    subroutine read_cis(dformat, path, wfa, wfb, occ_mo, act_mo, cisa, cisb, occ_num, norm, orthog)
        use occupation_mod
        use wf_norm_mod
        use orthog_mod
        character(len=*), intent(in) :: dformat
        character(len=*), intent(in) :: path
        real(dp), allocatable, optional :: wfa(:, :)
        real(dp), allocatable, optional :: wfb(:, :)
        logical, allocatable, optional :: occ_mo(:, :)
        logical, allocatable, optional :: act_mo(:, :)
        real(dp), allocatable, optional :: cisa(:, :, :)
        real(dp), allocatable, optional :: cisb(:, :, :)
        type(occupation_numbers), optional :: occ_num
        logical, optional :: norm
        logical, optional :: orthog
        real(dp), allocatable :: twfa(:, :)
        real(dp), allocatable :: twfb(:, :)
        logical, allocatable :: tocc_mo(:, :)
        logical, allocatable :: tact_mo(:, :)
        type(occupation_numbers) :: on
        integer :: nex, rhf
        logical :: tnorm
        logical :: torthog

        tnorm = .true.
        if (present(norm)) tnorm = norm
        torthog = .true.
        if (present(orthog)) torthog = orthog

        select case(dformat)
        case('turbomole')
            call turbomole_read_cis(path, twfa, twfb, tocc_mo, tact_mo)
        case default
            write(stderr, *) 'Error in read_cis subroutine.'
            write(stderr, *) 'Input format not implemented: '//trim(adjustl(dformat))
        end select

        if (tnorm) call wfab_normalize(allocated(twfb), twfa, twfb)
        !---- Bad for unrestricted
        if (torthog) then
            call orthog_lowdin(twfa)
            if (allocated(twfb)) call orthog_lowdin(twfa)
        end if
        !----
        rhf = 1
        if (allocated(twfb)) rhf = 2
        call occ_mask_to_nums(rhf, tocc_mo, tact_mo, on)
        nex = size(twfa, 2)

        if (present(wfa) .and. allocated(twfa)) wfa = twfa
        if (present(wfb) .and. allocated(twfb)) wfb = twfb
        if (present(occ_mo) .and. allocated(tocc_mo)) occ_mo = tocc_mo
        if (present(act_mo) .and. allocated(tact_mo)) act_mo = tact_mo
        if (present(cisa) .and. allocated(twfa)) cisa = reshape(twfa, [on%av(1), on%ao(1), nex])
        if (present(cisb) .and. allocated(twfb)) cisb = reshape(twfb, [on%av(2), on%ao(2), nex])
        if (present(occ_num)) occ_num = on
    end subroutine read_cis


end module read_all_mod
