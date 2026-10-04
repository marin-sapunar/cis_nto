!----------------------------------------------------------------------------------------------
! MODULE: cp2k_read_mod
!> @author Hamid Zabihi Hashjin
!> @date October, 2026
!
!> @brief Read data from CP2K .molden and .log files.
!> @details
!! Path convention: "molden_path:log_path" -- a single string with the two file paths
!! separated by a colon.
!!
!! The CP2K backend reads:
!!   - geometry, basis set, and molecular orbitals from a .molden file
!!     (produced by CP2K's &DFT &PRINT &MO_MOLDEN block)
!!   - excitation energies and CIS amplitudes from the &TDDFPT block of
!!     the CP2K log file.
!!
!! Notes:
!!  - occ_mo / act_mo are sized (nmo_scf + added_mos, 2). The two-column
!!    shape matches the occupation_numbers type, which declares
!!    o(2), v(2), ao(2), ... . Restricted (closed-shell) CP2K calculations
!!    mirror the single spin channel into both columns.
!!  - If the CP2K log's TDDFPT block references more MOs than the molden
!!    actually contains (CP2K <= 2024.1 does not write the &MO_MOLDEN
!!    ADDED_MOS virtuals to the molden by default), added_mos is capped
!!    at the molden's MO count and a warning is printed.
!!  - Uses cp2k_read_tddft_block_dynamic so the excitation arrays are
!!    sized exactly to the log content (no fixed 1000-amplitude limit).
!!  - Restricted TDDFPT only. Unrestricted support would require reading
!!    a beta-spin TDDFPT block; not implemented.
!----------------------------------------------------------------------------------------------
module cp2k_read_mod
    use global_defs
    use file_mod, only : reader
    use molden_read_mod
    use cp2k_sections_mod
    implicit none

    private
    public :: cp2k_read_geom
    public :: cp2k_read_basis
    public :: cp2k_read_mo
    public :: cp2k_read_cis

contains


    !------------------------------------------------------------------------------------------
    ! SUBROUTINE: cp2k_split_paths
    !> @brief Split "molden_path:log_path" into its two halves.
    !> @details
    !! Whitespace around the colon is stripped. If no colon is present, both
    !! halves equal the input (useful when only a molden file is needed).
    !------------------------------------------------------------------------------------------
    subroutine cp2k_split_paths(path, molden_path, log_path)
        character(len=*), intent(in)  :: path
        character(len=:), allocatable, intent(out) :: molden_path, log_path
        integer :: k
        character(len=:), allocatable :: work

        work = trim(adjustl(path))

        k = index(work, ':')
        if (k == 0) then
            molden_path = work
            log_path    = work
        else
            molden_path = trim(adjustl(work(:k-1)))
            log_path    = trim(adjustl(work(k+1:)))
        end if
    end subroutine cp2k_split_paths


    !------------------------------------------------------------------------------------------
    ! SUBROUTINE: cp2k_read_geom
    !> @brief Read geometry and atom list from the molden half of the path.
    !------------------------------------------------------------------------------------------
    subroutine cp2k_read_geom(path, geom, at_symbol, at_number)
        character(len=*), intent(in) :: path
        real(dp), allocatable :: geom(:)
        character(len=2), allocatable :: at_symbol(:)
        integer, allocatable :: at_number(:)
        character(len=:), allocatable :: molden_path, log_path

        call cp2k_split_paths(path, molden_path, log_path)
        call molden_read_geom(molden_path, geom, at_symbol, at_number)
    end subroutine cp2k_read_geom


    !------------------------------------------------------------------------------------------
    ! SUBROUTINE: cp2k_read_basis
    !> @brief Read basis set from the molden half of the path.
    !------------------------------------------------------------------------------------------
    subroutine cp2k_read_basis(path, bs)
        use basis_set_mod, only : basis_set
        character(len=*), intent(in) :: path
        type(basis_set), intent(out) :: bs
        character(len=:), allocatable :: molden_path, log_path

        call cp2k_split_paths(path, molden_path, log_path)
        call molden_read_basis(molden_path, bs, 1)
    end subroutine cp2k_read_basis


    !------------------------------------------------------------------------------------------
    ! SUBROUTINE: cp2k_read_mo
    !> @brief Read molecular orbitals from the molden half of the path.
    !------------------------------------------------------------------------------------------
    subroutine cp2k_read_mo(path, mos, bs, fix_mo)
        use molecular_orbitals_mod, only : molecular_orbitals
        use basis_set_mod, only : basis_set
        character(len=*), intent(in) :: path
        type(molecular_orbitals) :: mos
        type(basis_set), intent(in), optional :: bs
        integer, intent(in), optional :: fix_mo
        character(len=:), allocatable :: molden_path, log_path
        type(basis_set) :: bs_dummy
        integer :: fix

        call cp2k_split_paths(path, molden_path, log_path)

        fix = 0
        if (present(fix_mo)) fix = fix_mo

        if (present(bs)) then
            call molden_read_mo(molden_path, mos, bs, fix)
        else
            call molden_read_mo(molden_path, mos, bs_dummy, fix)
        end if
    end subroutine cp2k_read_mo


    !------------------------------------------------------------------------------------------
    ! SUBROUTINE: cp2k_count_molden_mos
    !> @brief Count the number of "Ene=" lines in a molden file.
    !> @details
    !! Used to detect the case where the CP2K log's TDDFPT block references
    !! more MOs than the molden actually contains. In CP2K <= 2024.1 the
    !! &MO_MOLDEN block prints only the SCF MOs, while the TDDFPT section
    !! may request ADDED_MOS virtuals beyond that manifold.
    !------------------------------------------------------------------------------------------
    subroutine cp2k_count_molden_mos(molden_path, nmo)
        character(len=*), intent(in) :: molden_path
        integer, intent(out) :: nmo
        integer :: unit, ios
        character(len=512) :: line

        nmo = 0
        open(newunit=unit, file=molden_path, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(stderr, *) 'Error in cp2k_count_molden_mos: cannot open '//trim(molden_path)
            return
        end if
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'Ene=') > 0) nmo = nmo + 1
        end do
        close(unit)
    end subroutine cp2k_count_molden_mos


    !------------------------------------------------------------------------------------------
    ! SUBROUTINE: cp2k_read_cis
    !> @brief Read CIS/TDDFT excitation data from the log half of the path.
    !> @details
    !! Produces wfa(nocc*nvir, nstates) in the ordering expected by cis_nto:
    !!   idx = (occ-1)*nvir + virt
    !!
    !! occ_mo and act_mo are sized (nmo_scf + added_mos, 2). The second
    !! spin column mirrors the first (restricted calculations only).
    !!
    !! If the molden file contains fewer MOs than nmo_scf + added_mos, the
    !! excess added virtuals cannot be represented in the AO basis, so
    !! added_mos is capped and a warning is printed.
    !------------------------------------------------------------------------------------------
    subroutine cp2k_read_cis(path, wfa, occ_mo, act_mo)
        character(len=*), intent(in) :: path
        real(dp), allocatable, intent(out) :: wfa(:,:)
        logical, allocatable, intent(out) :: occ_mo(:,:)
        logical, allocatable, intent(out) :: act_mo(:,:)
        character(len=:), allocatable :: molden_path, log_path

        type(reader) :: readf
        integer :: nocc, nmo_scf, added_mos, nstates, nvir, ntot
        integer :: nmo_molden, n_capped
        integer, allocatable :: n_excitations(:), occ_indices(:,:), virt_indices(:,:)
        real(dp), allocatable :: exc_energies(:), ci_amplitudes(:,:)
        integer :: istate, j, occ, virt, idx

        call cp2k_split_paths(path, molden_path, log_path)

        ! --- Parse nocc, nmo_scf, added_mos from the log --------------------------
        call readf%open(log_path, skip_empty=.true., abort_on_eof=.false.)
        call cp2k_parse_nocc_log(readf, nocc)
        call cp2k_parse_nmo_scf(readf, nmo_scf)
        call cp2k_parse_added_mos(readf, added_mos)
        call readf%close()

        if (nocc < 1) then
            write(stderr, *) 'Error in cp2k_read_cis: could not read nocc from log.'
            call abort()
        end if
        if (nmo_scf < 1) then
            write(stderr, *) 'Error in cp2k_read_cis: could not read nmo_scf from log.'
            call abort()
        end if
        if (added_mos < 0) added_mos = 0

        ! --- Cap added_mos at what the molden actually provides ------------------
        call cp2k_count_molden_mos(molden_path, nmo_molden)
        if (nmo_molden < 1) then
            write(stderr, *) 'Error in cp2k_read_cis: molden file contains no MOs.'
            write(stderr, *) '  path: '//molden_path
            call abort()
        end if
        if (nmo_scf + added_mos > nmo_molden) then
            n_capped = max(0, nmo_molden - nmo_scf)
            write(*, '(A)') 'Warning in cp2k_read_cis:'
            write(*, '(A,I0)') '  CIS log references ', nmo_scf + added_mos, ' MOs'
            write(*, '(A,I0)') '  molden contains ',      nmo_molden, ' MOs'
            write(*, '(A,I0)') '  capping added_mos to ', n_capped
            write(*, '(A)') '  (CP2K <= 2024.1 does not write &MO_MOLDEN added virtuals.'
            write(*, '(A)') '   Set ADDED_MOS inside &MO_MOLDEN to include them.)'
            added_mos = n_capped
        end if

        nvir = (nmo_scf + added_mos) - nocc
        ntot =  nmo_scf + added_mos

        if (nvir < 1) then
            write(stderr, *) 'Error in cp2k_read_cis: nvir <= 0. Check nocc/nmo.'
            call abort()
        end if

        ! --- Parse the TDDFT block ----------------------------------------------
        call readf%open(log_path, skip_empty=.true., abort_on_eof=.false.)
        call cp2k_read_tddft_block_dynamic(readf, nocc, nstates, &
                                           n_excitations, exc_energies, &
                                           occ_indices, virt_indices, ci_amplitudes)
        call readf%close()

        if (nstates < 1) then
            write(stderr, *) 'Error in cp2k_read_cis: no excited states found.'
            call abort()
        end if

        ! --- Allocate output arrays ---------------------------------------------
        allocate(wfa(nocc * nvir, nstates), source=0.0_dp)
        allocate(occ_mo(ntot, 2), source=.false.)
        allocate(act_mo(ntot, 2), source=.true.)

        ! Mark occupied orbitals in both spin channels.
        occ_mo(1:nocc, 1) = .true.
        occ_mo(1:nocc, 2) = .true.

        ! --- Populate wfa -------------------------------------------------------
        do istate = 1, nstates
            do j = 1, n_excitations(istate)
                occ  = occ_indices(j, istate)
                virt = virt_indices(j, istate) - nocc

                if (occ  < 1 .or. occ  > nocc) cycle
                if (virt < 1 .or. virt > nvir) cycle

                idx = (occ - 1) * nvir + virt
                wfa(idx, istate) = ci_amplitudes(j, istate)
            end do
        end do
    end subroutine cp2k_read_cis


end module cp2k_read_mod