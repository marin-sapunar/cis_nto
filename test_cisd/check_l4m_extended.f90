program check_l4m_extended
    use cis_overlap_cisd_mod,          only: cis_overlap_cisd_bruteforce, dp
    use cis_overlap_cisd_l4m_mod,      only: cis_overlap_cisd_l4m
    use cis_overlap_cisd_l4m_fast_mod, only: cis_overlap_cisd_l4m_fast
    implicit none

    write(*,'(A)') '=== Extended L4M-CISD validation ==='
    write(*,'(A45, 2X, A12, 1X, A12, 1X, A6, 2X, A12)') &
        'Test', '|brute-L4M|', '|brute-fast|', 'Status', 'S(1,1)'

    call test_case('identity-MO, same dets',             60,20,20, 2,'id',  'same')
    call test_case('identity-MO, singles (rank-1)',      60,20,20, 2,'id',  'sing')
    call test_case('identity-MO, doubles (rank-2)',      60,20,20, 2,'id',  'doub')
    call test_case('tiny-MO  (0.001), singles',          60,20,20, 2,'tiny','sing')
    call test_case('tiny-MO  (0.001), doubles',          60,20,20, 2,'tiny','doub')
    call test_case('medium-MO(0.05),  singles',          60,20,20, 2,'med', 'sing')
    call test_case('medium-MO(0.05),  doubles',          60,20,20, 2,'med', 'doub')
    call test_case('medium-MO, doubles, n_dets=40',      60,20,40, 2,'med', 'doub')
    call test_case('medium-MO, doubles, nstates=4',      60,20,20, 4,'med', 'doub')
    call test_case('medium-MO, doubles, nstates=8',      60,20,20, 8,'med', 'doub')
    call test_case('medium-MO, mixed sing+doub',         60,20,20, 2,'med', 'mix')
    call test_case('small n_occ=10, doubles',            30,10,20, 2,'med', 'doub')
    call test_case('large n_occ=40, doubles',           100,40,20, 2,'med', 'doub')
    call test_case('realistic MO, singles',              60,20,20, 2,'realmo','sing')
    call test_case('realistic MO, doubles',              60,20,20, 2,'realmo','doub')
    call test_case('realistic MO, doubles, n_dets=40',   60,20,40, 2,'realmo','doub')
    call test_case('near-singular A_ref, doubles',       60,20,20, 2,'near_sing','doub')
    call test_case('sign-flip MO, doubles',              60,20,20, 2,'signflip','doub')

    write(*,'(A)') '=== done ==='
contains

    subroutine test_case(label, n_mo, n_occ, n_dets, nstates, mo_kind, det_kind)
        character(len=*), intent(in) :: label, mo_kind, det_kind
        integer, intent(in) :: n_mo, n_occ, n_dets, nstates

        real(dp), allocatable :: s_mo_a(:,:), s_mo_b(:,:)
        integer,  allocatable :: dets_a_bra(:,:,:), dets_b_bra(:,:,:)
        integer,  allocatable :: dets_a_ket(:,:,:), dets_b_ket(:,:,:)
        real(dp), allocatable :: coeffs_bra(:,:), coeffs_ket(:,:)
        real(dp), allocatable :: s_brute(:,:), s_l4m(:,:), s_fast(:,:)
        real(dp) :: d_l4m, d_fast
        character(len=6) :: status

        allocate(s_mo_a(n_mo,n_mo), s_mo_b(n_mo,n_mo))
        allocate(dets_a_bra(n_occ,n_dets,nstates), dets_b_bra(n_occ,n_dets,nstates))
        allocate(dets_a_ket(n_occ,n_dets,nstates), dets_b_ket(n_occ,n_dets,nstates))
        allocate(coeffs_bra(n_dets,nstates), coeffs_ket(n_dets,nstates))
        allocate(s_brute(nstates,nstates), s_l4m(nstates,nstates), s_fast(nstates,nstates))

        call build_mo(s_mo_a, s_mo_b, n_mo, mo_kind)
        call build_dets(dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
                        n_occ, n_dets, nstates, det_kind)
        coeffs_bra = 0.1_dp
        coeffs_ket = 0.1_dp

        call cis_overlap_cisd_bruteforce(n_mo, n_occ, n_dets, nstates, &
            s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
            coeffs_bra, coeffs_ket, s_brute)
        call cis_overlap_cisd_l4m(n_mo, n_occ, n_dets, nstates, &
            s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
            coeffs_bra, coeffs_ket, s_l4m)
        call cis_overlap_cisd_l4m_fast(n_mo, n_occ, n_dets, nstates, &
            s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
            coeffs_bra, coeffs_ket, s_fast)

        d_l4m  = maxval(abs(s_brute - s_l4m))
        d_fast = maxval(abs(s_brute - s_fast))
        if (d_fast < 1.0e-12_dp .and. d_l4m < 1.0e-12_dp) then
            status = 'PASS'
        else
            status = 'FAIL'
        end if
        write(*,'(A45, 2X, ES12.4, 1X, ES12.4, 1X, A6, 2X, ES12.4)') &
            label, d_l4m, d_fast, status, s_brute(1,1)

        deallocate(s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
                   coeffs_bra, coeffs_ket, s_brute, s_l4m, s_fast)
    end subroutine

    subroutine build_mo(sa, sb, n, kind)
        real(dp), intent(out) :: sa(:,:), sb(:,:)
        integer, intent(in) :: n
        character(len=*), intent(in) :: kind
        real(dp) :: x
        integer :: i, j
        sa = 0.0_dp; sb = 0.0_dp
        do i = 1, n
            sa(i,i) = 1.0_dp; sb(i,i) = 1.0_dp
        end do
        select case (kind)
        case ('id');   x = 0.0_dp
        case ('tiny'); x = 0.001_dp
        case ('med');  x = 0.05_dp
        case ('realmo')
            ! realistic: decaying off-diagonals, mixed signs, cond ~ 20
            do i = 1, n
                do j = i+1, n
                    x = 0.15_dp * exp(-real(abs(i-j), dp)/8.0_dp) * &
                        cos(real(i*j, dp)*0.37_dp)
                    sa(i,j) = x; sa(j,i) = x
                    sb(i,j) = x; sb(j,i) = x
                end do
            end do
            return
        case ('near_sing')
            ! two nearly-degenerate orbitals (trivial-crossing geometry)
            do i = 1, n
                do j = i+1, n
                    x = 0.49_dp * exp(-real(abs(i-j), dp))
                    sa(i,j) = x; sa(j,i) = x
                    sb(i,j) = x; sb(j,i) = x
                end do
            end do
            return
        case ('signflip')
            ! realistic + one orbital phase-flipped in ket
            do i = 1, n
                do j = i+1, n
                    x = 0.05_dp
                    if (mod(i,3) == 0 .or. mod(j,3) == 0) x = -x
                    sa(i,j) = x; sa(j,i) = x
                    sb(i,j) = x; sb(j,i) = x
                end do
            end do
            do i = 1, n
                sb(i, 3) = -sb(i, 3)
            end do
            return
        case default;  x = 0.02_dp
        end select
        if (x > 0.0_dp) then
            do i = 1, n
                do j = i+1, n
                    sa(i,j) = x; sa(j,i) = x
                    sb(i,j) = x; sb(j,i) = x
                end do
            end do
        end if
    end subroutine

    subroutine build_dets(ab, bb, ak, bk, n_occ, n_dets, nstates, kind)
        integer, intent(out) :: ab(:,:,:), bb(:,:,:), ak(:,:,:), bk(:,:,:)
        integer, intent(in) :: n_occ, n_dets, nstates
        character(len=*), intent(in) :: kind
        integer :: d, i, k

        do k = 1, nstates
            do d = 1, n_dets
                do i = 1, n_occ
                    ab(i,d,k) = i; bb(i,d,k) = i
                    ak(i,d,k) = i; bk(i,d,k) = i
                end do
            end do
        end do

        do k = 1, nstates
            do d = 1, n_dets
                select case (kind)
                case ('same')
                    ! no change
                case ('sing')
                    if (d <= n_occ) then
                        ak(d,d,k) = n_occ + d; bk(d,d,k) = n_occ + d
                    end if
                case ('doub')
                    if (d <= n_occ - 1) then
                        ak(d,  d,k) = n_occ + d;   bk(d,  d,k) = n_occ + d
                        ak(d+1,d,k) = n_occ + d+1; bk(d+1,d,k) = n_occ + d+1
                    else if (d <= n_occ) then
                        ak(d,d,k) = n_occ + d; bk(d,d,k) = n_occ + d
                    end if
                case ('mix')
                    if (mod(d,3) == 1 .and. d <= n_occ) then
                        ak(d,d,k) = n_occ + d; bk(d,d,k) = n_occ + d
                    else if (mod(d,3) == 2 .and. d <= n_occ - 1) then
                        ak(d,  d,k) = n_occ + d;   bk(d,  d,k) = n_occ + d
                        ak(d+1,d,k) = n_occ + d+1; bk(d+1,d,k) = n_occ + d+1
                    end if
                end select
            end do
        end do
    end subroutine

end program
