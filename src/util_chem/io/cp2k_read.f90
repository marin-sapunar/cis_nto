module cp2k_read_mod
    use global_defs
    use file_mod, only : reader, check_string_in_file
    use molden_read_mod
    use cp2k_sections_mod
    use cp2k_virtuals_mod, only : reconstruct_virtuals, reconstruction_info
    implicit none
    private
    public :: cp2k_read_geom, cp2k_read_basis, cp2k_read_mo, cp2k_read_cis
contains
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
    end subroutine
    subroutine cp2k_read_geom(path, geom, at_symbol, at_number)
        character(len=*), intent(in) :: path
        real(dp), allocatable :: geom(:)
        character(len=2), allocatable :: at_symbol(:)
        integer, allocatable :: at_number(:)
        character(len=:), allocatable :: molden_path, log_path
        call cp2k_split_paths(path, molden_path, log_path)
        call molden_read_geom(molden_path, geom, at_symbol, at_number)
    end subroutine
    subroutine cp2k_read_basis(path, bs)
        use basis_set_mod, only : basis_set
        character(len=*), intent(in) :: path
        type(basis_set), intent(out) :: bs
        character(len=:), allocatable :: molden_path, log_path
        call cp2k_split_paths(path, molden_path, log_path)
        call molden_read_basis(molden_path, bs, 1)
    end subroutine
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
        if (present(bs)) then
            call molden_read_mo(molden_path, mos, bs, fix)
        else
            call molden_read_mo(molden_path, mos, bs_dummy, fix)
        end if
    end subroutine
    subroutine cp2k_count_molden_mos(molden_path, nmo)
        character(len=*), intent(in) :: molden_path
        integer, intent(out) :: nmo
        integer :: unit, ios
        character(len=512) :: line
        nmo = 0
        open(newunit=unit, file=molden_path, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'Ene=') > 0) nmo = nmo + 1
        end do
        close(unit)
    end subroutine
    subroutine cp2k_read_cis(path, wfa, wfb, occ_mo, act_mo)
        character(len=*), intent(in) :: path
        real(dp), allocatable, intent(out) :: wfa(:,:), wfb(:,:)
        logical, allocatable, intent(out) :: occ_mo(:,:), act_mo(:,:)
        character(len=:), allocatable :: molden_path, log_path

        type(reader) :: readf
        integer :: nocc_a, nocc_b, nmo_scf, added_mos, nstates
        integer :: nvir_a, nvir_b, ntot, nmo_molden
        integer :: istate, occ, virt, i, io, nheader_dashes
        integer :: io_u, ios_u
        real(dp) :: amp
        character(len=16) :: s1, s2
        character(len=256) :: raw, uks_line
        logical :: is_uks, header_done, found_spin2
        type(reconstruction_info) :: recon_info

        call cp2k_split_paths(path, molden_path, log_path)

        ! Detect UKS directly via file scan
        is_uks = .false.
        open(newunit=io_u, file=log_path, status='old', action='read', iostat=ios_u)
        if (ios_u == 0) then
            do
                read(io_u, '(A)', iostat=ios_u) uks_line
                if (ios_u /= 0) exit
                if (index(uks_line, "U-TDDFPT") > 0) then
                    is_uks = .true.
                    exit
                end if
            end do
            close(io_u)
        end if

        ! --- Parse nocc_a via reader class ---
        call readf%open(log_path, skip_empty=.true., abort_on_eof=.false.)
        call cp2k_parse_nocc_log(readf, nocc_a)
        call readf%close()

        ! --- Parse nocc_b directly via I/O for UKS ---
        if (is_uks) then
            found_spin2 = .false.
            nocc_b = -1
            open(newunit=io_u, file=log_path, status='old', action='read', iostat=ios_u)
            do
                read(io_u, '(A)', iostat=ios_u) uks_line
                if (ios_u /= 0) exit
                if (.not. found_spin2) then
                    if (index(uks_line, "Spin 2") > 0) found_spin2 = .true.
                else
                    if (index(uks_line, "Number of occupied orbitals") > 0) then
                        read(uks_line(scan(uks_line, ':')+1:), *) nocc_b
                        exit
                    end if
                end if
            end do
            close(io_u)
            if (nocc_b < 1) nocc_b = nocc_a
        else
            nocc_b = nocc_a
        end if

        ! --- Parse nmo_scf and added_mos via reader ---
        call readf%open(log_path, skip_empty=.true., abort_on_eof=.false.)
        call cp2k_parse_nmo_scf(readf, nmo_scf)
        call cp2k_parse_added_mos(readf, added_mos)
        call readf%close()

        write(*, '(A,I0,A,I0,A,L1)') 'DEBUG: nocc_a=', nocc_a, ' nocc_b=', nocc_b, ' is_uks=', is_uks

        if (nocc_a < 1 .or. nocc_b < 1 .or. nmo_scf < 1) call abort()
        if (added_mos < 0) added_mos = 0

        call cp2k_count_molden_mos(molden_path, nmo_molden)
        if (nmo_molden < 1) call abort()
        if (nmo_scf + added_mos > nmo_molden) then
            call reconstruct_virtuals(molden_path, log_path, nmo_scf, nocc_a, added_mos, recon_info)
            if (recon_info%method == "truncated") added_mos = max(0, nmo_molden - nmo_scf)
        end if

        nvir_a = (nmo_scf + added_mos) - nocc_a
        nvir_b = (nmo_scf + added_mos) - nocc_b
        ntot =  nmo_scf + added_mos
        if (nvir_a < 1 .or. nvir_b < 1) call abort()

        ! --- Count states via direct I/O ---
        nstates = 0
        open(newunit=io, file=log_path, status='old', action='read')
        nheader_dashes = 0
        do
            read(io, '(A)', end=100) raw
            if (index(raw, "Excitation analysis") > 0) then
                do
                    read(io, '(A)', end=100) raw
                    if (index(raw, "-------") > 0) nheader_dashes = nheader_dashes + 1
                    if (nheader_dashes == 2) exit
                end do
                do
                    read(io, '(A)', end=100) raw
                    if (index(raw, " eV") > 0) nstates = nstates + 1
                    if (index(raw, "ENERGY|") > 0) exit
                    if (index(raw, "PROGRAM ENDED") > 0) exit
                    if (index(raw, "DBCSR STATISTICS") > 0) exit
                    if (index(raw, "Natural Orbital") > 0) exit
                end do
                exit
            end if
        end do
100     close(io)
        if (nstates < 1) then
            write(stderr, *) 'Error: no excited states found.'
            call abort()
        end if

        allocate(wfa(nocc_a * nvir_a, nstates), source=0.0_dp)
        if (is_uks) then
            allocate(wfb(nocc_b * nvir_b, nstates), source=0.0_dp)
        else
            allocate(wfb(nocc_a * nvir_a, nstates), source=0.0_dp)
        end if

        ! --- Parse amplitudes via direct I/O ---
        open(newunit=io, file=log_path, status='old', action='read')
        header_done = .false.
        nheader_dashes = 0
        istate = 0
        do
            read(io, '(A)', end=200) raw
            if (.not. header_done) then
                if (index(raw, "Excitation analysis") > 0) then
                    do
                        read(io, '(A)', end=200) raw
                        if (index(raw, "-------") > 0) nheader_dashes = nheader_dashes + 1
                        if (nheader_dashes == 2) exit
                    end do
                    header_done = .true.
                end if
                cycle
            end if
            if (index(raw, "ENERGY|") > 0) exit
            if (index(raw, "PROGRAM ENDED") > 0) exit
            if (index(raw, "Natural Orbital") > 0) exit
            if (index(raw, "DBCSR STATISTICS") > 0) exit
            if (index(raw, " eV") > 0) then
                istate = istate + 1
                cycle
            end if
            if (istate < 1 .or. istate > nstates) cycle
            if (len_trim(raw) < 5) cycle

            if (is_uks) then
                read(raw, *, iostat=i) occ, s1, virt, s2, amp
                if (i == 0) then
                    if (index(s1, "alp") > 0) wfa((occ-1)*nvir_a + (virt-nocc_a), istate) = amp
                    if (index(s1, "bet") > 0) wfb((occ-1)*nvir_b + (virt-nocc_b), istate) = amp
                end if
            else
                read(raw, *, iostat=i) occ, virt, amp
                if (i == 0) then
                    wfa((occ-1)*nvir_a + (virt-nocc_a), istate) = amp
                    wfb((occ-1)*nvir_a + (virt-nocc_a), istate) = amp
                end if
            end if
        end do
200     close(io)

        allocate(occ_mo(ntot, 2), source=.false.)
        allocate(act_mo(ntot, 2), source=.true.)
        occ_mo(1:nocc_a, 1) = .true.
        occ_mo(1:nocc_b, 2) = .true.
    end subroutine cp2k_read_cis
end module cp2k_read_mod
