program test_fortran_cisd
    ! Verifies cis_overlap_cisd_bruteforce matches Python brute-force.
    ! Minimal case: n_mo=3, n_occ=1 -> 4 determinants (ref + 2 singles + 1 double)
    use cis_overlap_cisd_mod
    implicit none

    integer, parameter :: n_mo = 3, n_occ = 1, n_dets = 4, nstates = 2
    real(dp) :: s_mo_a(n_mo, n_mo), s_mo_b(n_mo, n_mo)
    integer :: dets_a_bra(n_occ, n_dets, nstates), dets_b_bra(n_occ, n_dets, nstates)
    integer :: dets_a_ket(n_occ, n_dets, nstates), dets_b_ket(n_occ, n_dets, nstates)
    real(dp) :: coeffs_bra(n_dets, nstates), coeffs_ket(n_dets, nstates)
    real(dp) :: s_wf(nstates, nstates)
    real(dp) :: coeffs(n_dets)
    integer :: i

    ! MO overlap: identity + small off-diagonal (like a 0.5 fs step)
    s_mo_a = 0.0_dp; s_mo_b = 0.0_dp
    do i = 1, n_mo
        s_mo_a(i, i) = 1.0_dp; s_mo_b(i, i) = 1.0_dp
    end do
    s_mo_a(1, 2) = 0.05_dp; s_mo_a(2, 1) = 0.05_dp
    s_mo_b = s_mo_a

    ! Determinants (1-based Fortran indices; orbital 1 = occupied, 2,3 = virtual)
    ! Det 1: reference (a1, b1)
    dets_a_bra(:, 1, 1) = [1]; dets_b_bra(:, 1, 1) = [1]
    dets_a_ket(:, 1, 1) = [1]; dets_b_ket(:, 1, 1) = [1]
    ! Det 2: alpha single (a2, b1)
    dets_a_bra(:, 2, 1) = [2]; dets_b_bra(:, 2, 1) = [1]
    dets_a_ket(:, 2, 1) = [2]; dets_b_ket(:, 2, 1) = [1]
    ! Det 3: beta single (a1, b2)
    dets_a_bra(:, 3, 1) = [1]; dets_b_bra(:, 3, 1) = [2]
    dets_a_ket(:, 3, 1) = [1]; dets_b_ket(:, 3, 1) = [2]
    ! Det 4: alpha-beta double (a2, b2)
    dets_a_bra(:, 4, 1) = [2]; dets_b_bra(:, 4, 1) = [2]
    dets_a_ket(:, 4, 1) = [2]; dets_b_ket(:, 4, 1) = [2]

    ! Duplicate state 1 as state 2 (for a simple 2-state test)
    do i = 1, n_dets
        dets_a_bra(:, i, 2) = dets_a_bra(:, i, 1)
        dets_b_bra(:, i, 2) = dets_b_bra(:, i, 1)
        dets_a_ket(:, i, 2) = dets_a_ket(:, i, 1)
        dets_b_ket(:, i, 2) = dets_b_ket(:, i, 1)
    end do

    ! CI coefficients (matched between bra/ket — self-overlap case)
    coeffs = [0.9_dp, 0.3_dp, 0.2_dp, 0.1_dp]
    coeffs = coeffs / sqrt(sum(coeffs**2))
    coeffs_bra(:, 1) = coeffs; coeffs_ket(:, 1) = coeffs
    coeffs_bra(:, 2) = coeffs; coeffs_ket(:, 2) = coeffs

    call cis_overlap_cisd_bruteforce(n_mo, n_occ, n_dets, nstates, &
                                     s_mo_a, s_mo_b, &
                                     dets_a_bra, dets_b_bra, &
                                     dets_a_ket, dets_b_ket, &
                                     coeffs_bra, coeffs_ket, s_wf)

    print '(A)', '=== Fortran CISD overlap ==='
    print '(A,F18.14)', 'S_11 = ', s_wf(1, 1)
    print '(A,F18.14)', 'S_22 = ', s_wf(2, 2)
    print '(A,F18.14)', 'S_12 = ', s_wf(1, 2)
    write(*, *)
    print '(A)', 'Expected: Fortran S_11 matches Python brute-force to machine precision.'
end program
