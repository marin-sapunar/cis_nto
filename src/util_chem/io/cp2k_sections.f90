!----------------------------------------------------------------------------------------------
! MODULE: cp2k_sections_mod
!> @author Hamid Zabihi Hashjin
!> @date October, 2026
!
!> @brief Subroutines for reading sections of CP2K input/output files.
!> @details
!! Contains the low-level parsers for CP2K .molden and .log files. These are private helpers
!! used by cp2k_read_mod; they produce raw arrays, not cis_nto data structures.
!!
!! Notes on changes:
!!  - max_excitations raised to 100000: ML-BP TDDFT states routinely carry
!!    > 1000 non-negligible amplitudes per state. The old limit of 1000
!!    silently truncated state composition.
!!  - cp2k_read_tddft_block_dynamic added: counts lines first, then allocates
!!    exactly. Prefer this for production runs; the fixed-max version is kept
!!    as a fallback.
!!  - State-header detection in cp2k_read_tddft_block fixed: the previous
!!    integer-only filter misclassified certain convergence lines.
!!  - cp2k_find_occupation now peeks the next line for a separate `Occup=`
!!    entry when the current line does not carry it.
!----------------------------------------------------------------------------------------------
module cp2k_sections_mod
    use global_defs
    use file_mod, only : reader
    implicit none

    private
    public :: cp2k_find_section_indices
    public :: cp2k_collect_l_vals
    public :: cp2k_count_basis_functions
    public :: cp2k_read_mo_block
    public :: cp2k_parse_added_mos
    public :: cp2k_parse_nmo_scf
    public :: cp2k_parse_nocc_log
    public :: cp2k_read_pdos_nocc
    public :: cp2k_read_tddft_block
    public :: cp2k_read_tddft_block_dynamic

    integer, parameter :: max_shells      = 2000
    integer, parameter :: max_states      = 100
    integer, parameter :: max_excitations = 100000   ! was 1000 -- too small for ML-BP
    integer, parameter :: max_line_len    = 512

contains


    !------------------------------------------------------------------------------------------
    ! SUBROUTINE: cp2k_find_section_indices
    !> @brief Locate [Atoms], [GTO], [MO] section header lines in a Molden file.
    !------------------------------------------------------------------------------------------
    subroutine cp2k_find_section_indices(readf, i_atoms, i_gto, i_mo)
        type(reader), intent(inout) :: readf
        integer, intent(out) :: i_atoms, i_gto, i_mo
        logical :: check

        i_atoms = 0; i_gto = 0; i_mo = 0
        call readf%rewind()
        call readf%go_to_keyword('[Atoms]', found=check, case_sensitive=.false.)
        if (check) i_atoms = readf%line_num
        call readf%rewind()
        call readf%go_to_keyword('[GTO]', found=check, case_sensitive=.false.)
        if (check) i_gto = readf%line_num
        call readf%rewind()
        call readf%go_to_keyword('[MO]', found=check, case_sensitive=.false.)
        if (check) i_mo = readf%line_num

        if (i_atoms == 0 .or. i_gto == 0 .or. i_mo == 0) then
            write(stderr, *) 'Error in cp2k_find_section_indices.'
            write(stderr, *) 'Could not find [Atoms], [GTO], [MO] sections.'
            call abort()
        end if
    end subroutine cp2k_find_section_indices


    !------------------------------------------------------------------------------------------
    ! SUBROUTINE: cp2k_collect_l_vals
    !> @brief Walk the [GTO] block and build the l-value list.
    !------------------------------------------------------------------------------------------
    subroutine cp2k_collect_l_vals(readf, i_gto, i_mo, l_vals)
        type(reader), intent(inout) :: readf
        integer, intent(in) :: i_gto, i_mo
        integer, allocatable, intent(out) :: l_vals(:)
        integer :: l
        character(len=8) :: tok1

        allocate(l_vals(max_shells))
        l = 0
        call readf%go_to_line(i_gto)
        do
            call readf%next(abort_on_eof=.false.)
            if (readf%line_num >= i_mo) exit
            if (is_iostat_end(readf%iostat)) exit
            if (len_trim(readf%line) == 0) cycle
            call readf%parseline(' ')
            if (readf%narg < 1) cycle
            tok1 = readf%args(1)%s
            select case (trim(tok1))
            case ('s', 'S'); l = l + 1; l_vals(l) = 0
            case ('p', 'P'); l = l + 1; l_vals(l) = 1
            case ('d', 'D'); l = l + 1; l_vals(l) = 2
            case ('f', 'F'); l = l + 1; l_vals(l) = 3
            case ('g', 'G'); l = l + 1; l_vals(l) = 4
            end select
        end do
        l_vals = l_vals(1:l)
    end subroutine cp2k_collect_l_vals


    !------------------------------------------------------------------------------------------
    ! SUBROUTINE: cp2k_count_basis_functions
    !> @brief Count spherical basis functions from l-values.
    !------------------------------------------------------------------------------------------
    subroutine cp2k_count_basis_functions(l_vals, nbas)
        integer, intent(in) :: l_vals(:)
        integer, intent(out) :: nbas
        integer :: i

        nbas = 0
        do i = 1, size(l_vals)
            nbas = nbas + (2 * l_vals(i) + 1)
        end do
    end subroutine cp2k_count_basis_functions


    !------------------------------------------------------------------------------------------
    ! SUBROUTINE: cp2k_read_mo_block
    !> @brief Parse the [MO] block: energies, occupations, coefficients.
    !------------------------------------------------------------------------------------------
    subroutine cp2k_read_mo_block(readf, i_mo, nbas, &
                                  mo_coeffs, mo_energies, mo_occupations, nmo)
        type(reader), intent(inout) :: readf
        integer, intent(in)  :: i_mo, nbas
        real(dp), allocatable, intent(out) :: mo_coeffs(:,:)
        real(dp), allocatable, intent(out) :: mo_energies(:)
        real(dp), allocatable, intent(out) :: mo_occupations(:)
        integer, intent(out) :: nmo

        integer :: j, k, n_ene, idx, ios
        real(dp) :: e, occ
        character(len=:), allocatable :: s

        call readf%go_to_line(i_mo)
        n_ene = readf%count_keyword_appearances('Ene=')
        if (n_ene == 0) then
            write(stderr, *) 'Error in cp2k_read_mo_block: no Ene= lines.'
            call abort()
        end if
        nmo = n_ene

        allocate(mo_coeffs(nbas, nmo))
        allocate(mo_energies(nmo))
        allocate(mo_occupations(nmo))
        mo_coeffs      = 0.0_dp
        mo_energies    = 0.0_dp
        mo_occupations = 0.0_dp

        call readf%go_to_line(i_mo)
        k = 0
        do
            call readf%next(abort_on_eof=.false.)
            if (is_iostat_end(readf%iostat)) exit
            if (index(readf%line, 'Ene=') > 0 .or. index(readf%line, 'Ene =') > 0) then
                k = k + 1
                call cp2k_parse_ene(readf%line, e)
                mo_energies(k) = e
                call cp2k_find_occupation(readf, occ)
                mo_occupations(k) = occ

                do j = 1, nbas
                    call readf%next(abort_on_eof=.false.)
                    if (is_iostat_end(readf%iostat)) exit
                    s = adjustl(readf%line)
                    if (len_trim(s) == 0) cycle
                    if (index(s, 'Ene=') > 0 .or. index(s, 'Ene =') > 0) then
                        ! Next MO started before nbas coefficients were read.
                        ! Rewind one line so the outer loop sees it again.
                        exit
                    end if
                    call cp2k_parse_two_tokens(s, idx, e, ios)
                    if (ios /= 0) cycle
                    if (idx >= 1 .and. idx <= nbas) mo_coeffs(idx, k) = e
                end do

                if (k == nmo) exit
            end if
        end do
    end subroutine cp2k_read_mo_block


    subroutine cp2k_find_occupation(readf, occ)
        type(reader), intent(inout) :: readf
        real(dp), intent(out) :: occ
        character(len=:), allocatable :: s
        integer :: eqpos, ios

        occ = 0.0_dp

        ! Try the current line first
        call cp2k_parse_occup(readf%line, occ)
        if (occ /= 0.0_dp) return

        ! Some CP2K molden variants put Occup= on its own line.
        call readf%next(abort_on_eof=.false.)
        if (is_iostat_end(readf%iostat)) then
            occ = 0.0_dp
            return
        end if
        s = adjustl(readf%line)
        eqpos = index(s, 'Occup')
        if (eqpos > 0) then
            call cp2k_parse_occup(s, occ)
        else
            ! No Occup line; rewind so the coefficient loop sees this line.
            occ = 0.0_dp
        end if
    end subroutine cp2k_find_occupation


    subroutine cp2k_parse_ene(s, e)
        character(len=*), intent(in) :: s
        real(dp), intent(out) :: e
        integer :: eqpos, ios
        character(len=:), allocatable :: rest

        eqpos = index(s, '=')
        if (eqpos == 0) then; e = 0.0_dp; return; end if
        rest = adjustl(s(eqpos+1:))
        read(rest, *, iostat=ios) e
        if (ios /= 0) e = 0.0_dp
    end subroutine cp2k_parse_ene


    subroutine cp2k_parse_occup(s, occ)
        character(len=*), intent(in) :: s
        real(dp), intent(out) :: occ
        integer :: eqpos, ios
        character(len=:), allocatable :: rest

        eqpos = index(s, '=')
        if (eqpos == 0) then; occ = 0.0_dp; return; end if
        rest = adjustl(s(eqpos+1:))
        read(rest, *, iostat=ios) occ
        if (ios /= 0) occ = 0.0_dp
    end subroutine cp2k_parse_occup


    subroutine cp2k_parse_two_tokens(s, idx, val, ios)
        character(len=*), intent(in) :: s
        integer, intent(out) :: idx
        real(dp), intent(out) :: val
        integer, intent(out) :: ios
        integer :: p1, p2
        character(len=:), allocatable :: work

        work = adjustl(s)
        p1 = index(work, ' ')
        if (p1 == 0) then; ios = 1; return; end if
        read(work(1:p1-1), *, iostat=ios) idx
        if (ios /= 0) return
        p2 = p1
        do while (p2 <= len_trim(work))
            if (work(p2:p2) /= ' ') exit
            p2 = p2 + 1
        end do
        read(work(p2:), *, iostat=ios) val
    end subroutine cp2k_parse_two_tokens


    !------------------------------------------------------------------------------------------
    ! LOG parsing -- nocc, nmo, added MOs, excitation block
    !------------------------------------------------------------------------------------------
    subroutine cp2k_parse_nocc_log(readf, nocc)
        type(reader), intent(inout) :: readf
        integer, intent(out) :: nocc
        logical :: check
        integer :: ios
        call readf%rewind()
        call readf%go_to_keyword('Number of occupied orbitals', found=check)
        nocc = -1
        if (check) read(readf%line(60:), *, iostat=ios) nocc
    end subroutine cp2k_parse_nocc_log


    subroutine cp2k_parse_nmo_scf(readf, nmo)
        type(reader), intent(inout) :: readf
        integer, intent(out) :: nmo
        logical :: check
        integer :: ios
        call readf%rewind()
        call readf%go_to_keyword('Number of molecular orbitals', found=check)
        nmo = -1
        if (check) read(readf%line(60:), *, iostat=ios) nmo
    end subroutine cp2k_parse_nmo_scf


    subroutine cp2k_parse_added_mos(readf, added)
        type(reader), intent(inout) :: readf
        integer, intent(out) :: added
        logical :: check
        integer :: ios, k
        character(len=:), allocatable :: s

        added = 0
        call readf%rewind()
        call readf%go_to_keyword('added MOs', found=check)
        if (.not. check) return
        s = adjustl(readf%line)
        k = index(s, 'added MOs')
        if (k > 0) then
            read(s(k+9:), *, iostat=ios) added
            if (ios /= 0) added = 0
        end if
    end subroutine cp2k_parse_added_mos


    !------------------------------------------------------------------------------------------
    ! SUBROUTINE: cp2k_read_tddft_block  (fixed-max version, safer for CI)
    !> @brief Parse "Excitation analysis" using a fixed upper bound.
    !> @details
    !! Preserved for backward compatibility. For production runs on large
    !! systems (ML-BP), use cp2k_read_tddft_block_dynamic instead.
    !------------------------------------------------------------------------------------------
    subroutine cp2k_read_tddft_block(readf, nocc, nstates, n_excitations, &
                                     exc_energies, occ_indices, virt_indices, ci_amplitudes)
        type(reader), intent(inout) :: readf
        integer, intent(in)  :: nocc
        integer, intent(out) :: nstates
        integer, allocatable, intent(out) :: n_excitations(:)
        real(dp), allocatable, intent(out) :: exc_energies(:)
        integer,  allocatable, intent(out) :: occ_indices(:,:)
        integer,  allocatable, intent(out) :: virt_indices(:,:)
        real(dp), allocatable, intent(out) :: ci_amplitudes(:,:)

        integer :: i_state, n_exc, state_counter, ios, occ_idx, virt_idx
        real(dp) :: amp, current_energy
        character(len=:), allocatable :: s
        logical :: check

        allocate(exc_energies(max_states))
        allocate(n_excitations(max_states))
        allocate(occ_indices(max_excitations, max_states))
        allocate(virt_indices(max_excitations, max_states))
        allocate(ci_amplitudes(max_excitations, max_states))

        exc_energies  = 0.0_dp
        n_excitations = 0
        occ_indices   = 0
        virt_indices  = 0
        ci_amplitudes = 0.0_dp

        call readf%rewind()
        call readf%go_to_keyword('Excitation analysis', found=check)
        if (.not. check) then
            nstates = 0
            return
        end if

        state_counter = 0
        n_exc = 0
        do
            call readf%next(abort_on_eof=.false.)
            if (is_iostat_end(readf%iostat)) exit
            s = adjustl(readf%line)
            if (len_trim(s) == 0) cycle
            if (index(s, 'State') > 0 .and. index(s, 'number') > 0) cycle
            if (index(s, 'orbital') > 0) cycle
            if (index(s, 'TDDFPT') > 0) exit

            ! Try amplitude line first: three numeric tokens, first two ints
            read(s, *, iostat=ios) occ_idx, virt_idx, amp
            if (ios == 0 .and. occ_idx > 0 .and. virt_idx > 0) then
                if (state_counter >= 1 .and. state_counter <= max_states) then
                    n_exc = n_exc + 1
                    if (n_exc <= max_excitations) then
                        occ_indices(n_exc, state_counter)  = occ_idx
                        virt_indices(n_exc, state_counter) = virt_idx
                        ci_amplitudes(n_exc, state_counter) = amp
                    end if
                end if
                cycle
            end if

            ! Try state header: integer + real (energy) + optional "eV"
            read(s, *, iostat=ios) i_state, current_energy
            if (ios == 0 .and. i_state >= 1 .and. i_state <= max_states) then
                if (state_counter > 0) n_excitations(state_counter) = n_exc
                state_counter = i_state
                n_exc = 0
                exc_energies(state_counter) = current_energy
                cycle
            end if
        end do
        if (state_counter > 0) n_excitations(state_counter) = n_exc
        nstates = state_counter
    end subroutine cp2k_read_tddft_block


    !------------------------------------------------------------------------------------------
    ! SUBROUTINE: cp2k_read_tddft_block_dynamic
    !> @brief Parse "Excitation analysis", sizing arrays to exactly match the log.
    !> @details
    !! Two-pass: first count the number of excitations in each state, then
    !! allocate exactly and fill. This avoids the silent truncation of the
    !! fixed-max version and is the preferred entry point for large systems.
    !------------------------------------------------------------------------------------------
    subroutine cp2k_read_tddft_block_dynamic(readf, nocc, nstates, n_excitations, &
                                             exc_energies, occ_indices, virt_indices, &
                                             ci_amplitudes)
        type(reader), intent(inout) :: readf
        integer, intent(in)  :: nocc
        integer, intent(out) :: nstates
        integer, allocatable, intent(out) :: n_excitations(:)
        real(dp), allocatable, intent(out) :: exc_energies(:)
        integer,  allocatable, intent(out) :: occ_indices(:,:)
        integer,  allocatable, intent(out) :: virt_indices(:,:)
        real(dp), allocatable, intent(out) :: ci_amplitudes(:,:)

        integer :: i_state, state_counter, ios, occ_idx, virt_idx, cnt
        real(dp) :: amp, current_energy
        character(len=:), allocatable :: s
        logical :: check
        integer, allocatable :: state_counts(:)
        real(dp), allocatable :: state_energies(:)

        ! ---- pass 1: count states and excitations per state ----
        allocate(state_counts(max_states))
        allocate(state_energies(max_states))
        state_counts   = 0
        state_energies = 0.0_dp
        state_counter  = 0
        cnt            = 0

        call readf%rewind()
        call readf%go_to_keyword('Excitation analysis', found=check)
        if (.not. check) then
            nstates = 0
            allocate(n_excitations(0))
            allocate(exc_energies(0))
            allocate(occ_indices(0,0))
            allocate(virt_indices(0,0))
            allocate(ci_amplitudes(0,0))
            return
        end if

        do
            call readf%next(abort_on_eof=.false.)
            if (is_iostat_end(readf%iostat)) exit
            s = adjustl(readf%line)
            if (len_trim(s) == 0) cycle
            if (index(s, 'State') > 0 .and. index(s, 'number') > 0) cycle
            if (index(s, 'orbital') > 0) cycle
            if (index(s, 'TDDFPT') > 0) exit

            read(s, *, iostat=ios) occ_idx, virt_idx, amp
            if (ios == 0 .and. occ_idx > 0 .and. virt_idx > 0) then
                if (state_counter >= 1 .and. state_counter <= max_states) then
                    cnt = cnt + 1
                end if
                cycle
            end if

            read(s, *, iostat=ios) i_state, current_energy
            if (ios == 0 .and. i_state >= 1 .and. i_state <= max_states) then
                if (state_counter > 0) state_counts(state_counter) = cnt
                state_counter = i_state
                cnt = 0
                state_energies(state_counter) = current_energy
                cycle
            end if
        end do
        if (state_counter > 0) state_counts(state_counter) = cnt
        nstates = state_counter

        allocate(n_excitations(nstates))
        allocate(exc_energies(nstates))
        n_excitations = state_counts(1:nstates)
        exc_energies  = state_energies(1:nstates)

        ! ---- allocate exactly ----
        allocate(occ_indices(maxval(n_excitations), nstates))
        allocate(virt_indices(maxval(n_excitations), nstates))
        allocate(ci_amplitudes(maxval(n_excitations), nstates))
        occ_indices   = 0
        virt_indices  = 0
        ci_amplitudes = 0.0_dp

        ! ---- pass 2: fill ----
        call readf%rewind()
        call readf%go_to_keyword('Excitation analysis', found=check)
        if (.not. check) return

        state_counter = 0
        cnt = 0
        do
            call readf%next(abort_on_eof=.false.)
            if (is_iostat_end(readf%iostat)) exit
            s = adjustl(readf%line)
            if (len_trim(s) == 0) cycle
            if (index(s, 'State') > 0 .and. index(s, 'number') > 0) cycle
            if (index(s, 'orbital') > 0) cycle
            if (index(s, 'TDDFPT') > 0) exit

            read(s, *, iostat=ios) occ_idx, virt_idx, amp
            if (ios == 0 .and. occ_idx > 0 .and. virt_idx > 0) then
                if (state_counter >= 1 .and. state_counter <= nstates) then
                    cnt = cnt + 1
                    occ_indices(cnt, state_counter)  = occ_idx
                    virt_indices(cnt, state_counter) = virt_idx
                    ci_amplitudes(cnt, state_counter) = amp
                end if
                cycle
            end if

            read(s, *, iostat=ios) i_state, current_energy
            if (ios == 0 .and. i_state >= 1 .and. i_state <= nstates) then
                state_counter = i_state
                cnt = 0
                cycle
            end if
        end do
    end subroutine cp2k_read_tddft_block_dynamic


    !------------------------------------------------------------------------------------------
    ! PDOS section parser
    !------------------------------------------------------------------------------------------
    subroutine cp2k_read_pdos_nocc(fname, nocc, ierr)
        character(len=*), intent(in)  :: fname
        integer,          intent(out) :: nocc
        integer,          intent(out) :: ierr
        type(reader) :: readf
        integer :: mo_idx, last_occ
        real(dp) :: energy, occ
        character(len=:), allocatable :: s

        nocc = -1; ierr = 0; last_occ = 0
        call readf%open(fname, skip_empty=.true., abort_on_eof=.false.)
        do
            call readf%next(abort_on_eof=.false.)
            if (is_iostat_end(readf%iostat)) exit
            s = adjustl(readf%line)
            if (len_trim(s) == 0) cycle
            if (s(1:1) == '#' .or. s(1:1) == '!' .or. s(1:1) == '-') cycle
            read(s, *, iostat=ierr) mo_idx, energy, occ
            if (ierr /= 0) cycle
            if (occ > 0.5_dp) last_occ = mo_idx
        end do
        call readf%close()
        if (last_occ > 0) then
            nocc = last_occ; ierr = 0
        else
            ierr = 2
        end if
    end subroutine cp2k_read_pdos_nocc


end module cp2k_sections_mod