program test_fortran_l4m
    ! Verify cis_overlap_cisd_l4m matches cis_overlap_cisd_bruteforce.
    use cis_overlap_cisd_mod, only: cis_overlap_cisd_bruteforce, dp
    use cis_overlap_cisd_l4m_mod, only: cis_overlap_cisd_l4m
    implicit none

    integer, parameter :: n_mo = 3, n_occ = 1, n_dets = 4, nstates = 2
    real(dp) :: s_mo_a(n_mo, n_mo), s_mo_b(n_mo, n_mo)
    integer :: dets_a_bra(n_occ, n_dets, nstates), dets_b_bra(n_occ, n_dets, nstates)
    integer :: dets_a_ket(n_occ, n_dets, nstates), dets_b_ket(n_occ, n_dets, nstates)
    real(dp) :: coeffs_bra(n_dets, nstates), coeffs_ket(n_dets, nstates)
    real(dp) :: s_wf_brute(nstates, nstates), s_wf_l4m(nstates, nstates)
    real(dp) :: coeffs(n_dets)
    integer :: i, j
    real(dp) :: max_diff

    ! MO overlap
    s_mo_a = 0.0_dp; s_mo_b = 0.0_dp
    do i = 1, n_mo
        s_mo_a(i, i) = 1.0_dp; s_mo_b(i, i) = 1.0_dp
    end do
    s_mo_a(1, 2) = 0.05_dp; s_mo_a(2, 1) = 0.05_dp
    s_mo_b = s_mo_a

    ! Determinants
    dets_a_bra(:, 1, 1) = [1]; dets_b_bra(:, 1, 1) = [1]
    dets_a_ket(:, 1, 1) = [1]; dets_b_ket(:, 1, 1) = [1]
    dets_a_bra(:, 2, 1) = [2]; dets_b_bra(:, 2, 1) = [1]
    dets_a_ket(:, 2, 1) = [2]; dets_b_ket(:, 2, 1) = [1]
    dets_a_bra(:, 3, 1) = [1]; dets_b_bra(:, 3, 1) = [2]
    dets_a_ket(:, 3, 1) = [1]; dets_b_ket(:, 3, 1) = [2]
    dets_a_bra(:, 4, 1) = [2]; dets_b_bra(:, 4, 1) = [2]
    dets_a_ket(:, 4, 1) = [2]; dets_b_ket(:, 4, 1) = [2]
    do i = 1, n_dets
        dets_a_bra(:, i, 2) = dets_a_bra(:, i, 1)
        dets_b_bra(:, i, 2) = dets_b_bra(:, i, 1)
        dets_a_ket(:, i, 2) = dets_a_ket(:, i, 1)
        dets_b_ket(:, i, 2) = dets_b_ket(:, i, 1)
    end do

    coeffs = [0.9_dp, 0.3_dp, 0.2_dp, 0.1_dp]
    coeffs = coeffs / sqrt(sum(coeffs**2))
    coeffs_bra(:, 1) = coeffs; coeffs_ket(:, 1) = coeffs
    coeffs_bra(:, 2) = coeffs; coeffs_ket(:, 2) = coeffs

    ! Brute-force
    call cis_overlap_cisd_bruteforce(n_mo, n_occ, n_dets, nstates, &
                                      s_mo_a, s_mo_b, &
                                      dets_a_bra, dets_b_bra, &
                                      dets_a_ket, dets_b_ket, &
                                      coeffs_bra, coeffs_ket, s_wf_brute)

    ! L4M
    call cis_overlap_cisd_l4m(n_mo, n_occ, n_dets, nstates, &
                               s_mo_a, s_mo_b, &
                               dets_a_bra, dets_b_bra, &
                               dets_a_ket, dets_b_ket, &
                               coeffs_bra, coeffs_ket, s_wf_l4m)

    max_diff = 0.0_dp
    do i = 1, nstates
        do j = 1, nstates
            max_diff = max(max_diff, abs(s_wf_brute(i, j) - s_wf_l4m(i, j)))
        end do
    end do

    print '(A)', '=== Fortran L4M vs Brute-force ==='
    print '(A, ES20.12)', 'S_11 (brute) = ', s_wf_brute(1, 1)
    print '(A, ES20.12)', 'S_11 (L4M)   = ', s_wf_l4m(1, 1)
    print '(A, ES20.12)', 'S_22 (brute) = ', s_wf_brute(2, 2)
    print '(A, ES20.12)', 'S_22 (L4M)   = ', s_wf_l4m(2, 2)
    print '(A, ES20.12)', 'Max |brute - L4M| = ', max_diff
    write(*, *)
    if (max_diff < 1.0d-13) then
        print '(A)', '✅ MATCH: Fortran L4M verified against brute-force'
    else
        print '(A)', '❌ MISMATCH'
    end if
end program
