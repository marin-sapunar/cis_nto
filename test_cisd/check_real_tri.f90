program check_real_tri
    use cis_overlap_cisd_mod,          only: cis_overlap_cisd_bruteforce, dp
    use cis_overlap_cisd_l4m_mod,      only: cis_overlap_cisd_l4m
    use cis_overlap_cisd_l4m_fast_mod, only: cis_overlap_cisd_l4m_fast
    implicit none

    integer :: n_mo, n_occ, n_dets, nstates
    real(dp), allocatable :: s_mo_a(:,:), s_mo_b(:,:)
    integer, allocatable :: dets_a_bra(:,:,:), dets_b_bra(:,:,:)
    integer, allocatable :: dets_a_ket(:,:,:), dets_b_ket(:,:,:)
    real(dp), allocatable :: coeffs_bra(:,:), coeffs_ket(:,:)
    real(dp) :: s_brute(1,1), s_l4m(1,1), s_fast(1,1)
    integer :: i, j, k
    real(dp) :: d_l4m, d_fast

    open(10, file='real_data.txt', status='old')
    read(10, *) n_mo, n_occ, n_dets, nstates

    allocate(s_mo_a(n_mo, n_mo), s_mo_b(n_mo, n_mo))
    allocate(dets_a_bra(n_occ, n_dets, 1), dets_b_bra(n_occ, n_dets, 1))
    allocate(dets_a_ket(n_occ, n_dets, 1), dets_b_ket(n_occ, n_dets, 1))
    allocate(coeffs_bra(n_dets, 1), coeffs_ket(n_dets, 1))

    do i = 1, n_mo
        read(10, *) (s_mo_a(i, j), j = 1, n_mo)
    end do
    s_mo_b = s_mo_a

    do k = 1, n_dets
        read(10, *) (dets_a_bra(i, k, 1), i = 1, n_occ)
        read(10, *) (dets_b_bra(i, k, 1), i = 1, n_occ)
        read(10, *) coeffs_bra(k, 1), coeffs_ket(k, 1)
        dets_a_ket(:, k, 1) = dets_a_bra(:, k, 1)
        dets_b_ket(:, k, 1) = dets_b_bra(:, k, 1)
    end do
    close(10)

    print *, "=== Three-way validation on real amyloid data ==="
    print *, "n_mo =", n_mo, " n_occ =", n_occ, " n_dets =", n_dets

    call cis_overlap_cisd_bruteforce(n_mo, n_occ, n_dets, 1, &
        s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
        coeffs_bra, coeffs_ket, s_brute)

    call cis_overlap_cisd_l4m(n_mo, n_occ, n_dets, 1, &
        s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
        coeffs_bra, coeffs_ket, s_l4m)

    call cis_overlap_cisd_l4m_fast(n_mo, n_occ, n_dets, 1, &
        s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
        coeffs_bra, coeffs_ket, s_fast)

    d_l4m  = abs(s_brute(1,1) - s_l4m(1,1))
    d_fast = abs(s_brute(1,1) - s_fast(1,1))

    write(*, '(A, F14.8)') ' S_brute = ', s_brute(1,1)
    write(*, '(A, F14.8)') ' S_L4M   = ', s_l4m(1,1)
    write(*, '(A, F14.8)') ' S_fast  = ', s_fast(1,1)
    write(*, '(A, ES14.4)') ' |brute - slow_L4M| = ', d_l4m
    write(*, '(A, ES14.4)') ' |brute - fast_L4M| = ', d_fast

    if (d_l4m < 1e-10_dp .and. d_fast < 1e-10_dp) then
        print *, "PASS (both L4M variants agree with brute)"
    else if (d_l4m < 1e-10_dp) then
        print *, "PARTIAL: slow L4M correct, fast L4M is broken"
    else if (d_fast < 1e-10_dp) then
        print *, "PARTIAL: fast L4M correct, slow L4M is broken"
    else
        print *, "FAIL (both L4M variants disagree with brute -> inputs are wrong)"
    end if
end program check_real_tri
