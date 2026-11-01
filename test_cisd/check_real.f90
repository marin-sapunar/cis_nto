program check_real
    use cis_overlap_cisd_mod,          only: cis_overlap_cisd_bruteforce, dp
    use cis_overlap_cisd_l4m_fast_mod, only: cis_overlap_cisd_l4m_fast
    implicit none

    integer :: n_mo, n_occ, n_dets
    real(dp), allocatable :: s_mo_a(:,:), s_mo_b(:,:)
    integer, allocatable :: dets_a_bra(:,:,:), dets_b_bra(:,:,:)
    integer, allocatable :: dets_a_ket(:,:,:), dets_b_ket(:,:,:)
    real(dp), allocatable :: coeffs_bra(:,:), coeffs_ket(:,:)
    real(dp) :: s_brute(1,1), s_fast(1,1)
    integer :: i, j, k
    real(dp) :: val

    open(10, file='real_data.txt', status='old')
    read(10, *) n_mo, n_occ, n_dets

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

    call cis_overlap_cisd_bruteforce(n_mo, n_occ, n_dets, 1, &
        s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
        coeffs_bra, coeffs_ket, s_brute)

    call cis_overlap_cisd_l4m_fast(n_mo, n_occ, n_dets, 1, &
        s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
        coeffs_bra, coeffs_ket, s_fast)

    val = abs(s_brute(1,1) - s_fast(1,1))
    write(*, '(A, ES16.6)') 'max|S_brute - S_fast| for Benzene NAMD step = ', val
    write(*, '(A, F10.6)') 'Overlap S(t1,t2) = ', s_brute(1,1)
    ! Use 1e-10 tolerance for 36x36 matrices
    if (val < 1.0e-10_dp) then
        write(*, *) 'PASS'
    else
        write(*, *) 'FAIL'
    end if

end program check_real
