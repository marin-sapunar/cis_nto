program check_l4m
    use cis_overlap_cisd_mod,          only: cis_overlap_cisd_bruteforce, dp
    use cis_overlap_cisd_l4m_mod,      only: cis_overlap_cisd_l4m
    use cis_overlap_cisd_l4m_fast_mod, only: cis_overlap_cisd_l4m_fast
    implicit none

    integer, parameter :: n_mo = 729, n_occ = 121, nstates = 2, n_dets = 60
    real(dp) :: s_mo_a(n_mo, n_mo), s_mo_b(n_mo, n_mo)
    integer :: dets_a_bra(n_occ, n_dets, nstates), dets_b_bra(n_occ, n_dets, nstates)
    integer :: dets_a_ket(n_occ, n_dets, nstates), dets_b_ket(n_occ, n_dets, nstates)
    real(dp) :: coeffs_bra(n_dets, nstates), coeffs_ket(n_dets, nstates)
    real(dp) :: s_brute(nstates, nstates), s_l4m(nstates, nstates), s_fast(nstates, nstates)
    integer  :: i, j, d, k

    s_mo_a = 0.0_dp; s_mo_b = 0.0_dp
    do i = 1, n_mo
        s_mo_a(i, i) = 1.0_dp; s_mo_b(i, i) = 1.0_dp
    end do
    do i = 1, n_mo
        do j = i+1, n_mo
            s_mo_a(i, j) = 0.02_dp; s_mo_a(j, i) = 0.02_dp
            s_mo_b(i, j) = 0.02_dp; s_mo_b(j, i) = 0.02_dp
        end do
    end do

    do k = 1, nstates
        do d = 1, n_dets
            do i = 1, n_occ
                dets_a_bra(i, d, k) = i
                dets_b_bra(i, d, k) = i
                dets_a_ket(i, d, k) = i
                dets_b_ket(i, d, k) = i
            end do
            if (d <= n_occ - 1) then
                dets_a_ket(d, d, k) = n_occ + d
                dets_b_ket(d, d, k) = n_occ + d
            end if
            coeffs_bra(d, k) = 0.1_dp
            coeffs_ket(d, k) = 0.1_dp
        end do
    end do

    call cis_overlap_cisd_bruteforce(n_mo, n_occ, n_dets, nstates, &
        s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
        coeffs_bra, coeffs_ket, s_brute)

    call cis_overlap_cisd_l4m(n_mo, n_occ, n_dets, nstates, &
        s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
        coeffs_bra, coeffs_ket, s_l4m)

    call cis_overlap_cisd_l4m_fast(n_mo, n_occ, n_dets, nstates, &
        s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
        coeffs_bra, coeffs_ket, s_fast)

    write(*, '(A)') 'brute S:'
    do i = 1, nstates
        write(*, '(2(F20.12, 1X))') (s_brute(i, j), j = 1, nstates)
    end do
    write(*, '(A)') 'L4M S:'
    do i = 1, nstates
        write(*, '(2(F20.12, 1X))') (s_l4m(i, j), j = 1, nstates)
    end do
    write(*, '(A)') 'fast S:'
    do i = 1, nstates
        write(*, '(2(F20.12, 1X))') (s_fast(i, j), j = 1, nstates)
    end do
    write(*, '(A, ES16.6)') 'max|S_brute - S_l4m|  = ', maxval(abs(s_brute - s_l4m))
    write(*, '(A, ES16.6)') 'max|S_brute - S_fast| = ', maxval(abs(s_brute - s_fast))
    if (maxval(abs(s_brute - s_fast)) < 1.0e-10_dp) then
        write(*, '(A)') 'PASS'
    else
        write(*, '(A)') 'FAIL'
    end if
end program
