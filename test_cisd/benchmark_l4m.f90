program benchmark_l4m
    use cis_overlap_cisd_mod, only: cis_overlap_cisd_bruteforce, dp
    use cis_overlap_cisd_l4m_mod, only: cis_overlap_cisd_l4m
    use cis_overlap_cisd_l4m_fast_mod, only: cis_overlap_cisd_l4m_fast
    implicit none

    integer, parameter :: n_mo = 729, n_occ = 121, nstates = 16
    integer, parameter :: max_dets = 200
    real(dp) :: s_mo_a(n_mo, n_mo), s_mo_b(n_mo, n_mo)
    integer :: dets_a_bra(n_occ, max_dets, nstates), dets_b_bra(n_occ, max_dets, nstates)
    integer :: dets_a_ket(n_occ, max_dets, nstates), dets_b_ket(n_occ, max_dets, nstates)
    real(dp) :: coeffs_bra(max_dets, nstates), coeffs_ket(max_dets, nstates)
    real(dp) :: s_wf(nstates, nstates), coeffs(max_dets)
    real(dp) :: t0, t1
    real(dp) :: t_brute, t_l4m, t_fast
    integer :: n_dets_list(2) = [60, 100]
    integer :: it, i, j, k, d, trial, n_dets, n_reps, v

    ! MO overlap
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

    print '(A)', '========================================================'
    print '(A)', ' Phase 3 Benchmark: brute vs L4M vs L4M-fast'
    print '(A)', '========================================================'
    print '(A,I0,A,I0,A,I0)', '  n_mo=', n_mo, ', n_occ=', n_occ, ', nstates=', nstates
    print '(A)', ''
    print '(A)', '  n_dets |  brute(ms) |  L4M(ms)   |  L4M-fast(ms) |  Speedup'
    print '(A)', '  -------+------------+------------+---------------+---------'

    do it = 1, size(n_dets_list)
        n_dets = n_dets_list(it)

        ! Zero all
        dets_a_bra = 0; dets_b_bra = 0; dets_a_ket = 0; dets_b_ket = 0
        d = 0

        ! Reference (occupations 1..n_occ)
        d = d + 1
        do i = 1, n_occ
            dets_a_bra(i, d, 1) = i
            dets_b_bra(i, d, 1) = i
            dets_a_ket(i, d, 1) = i
            dets_b_ket(i, d, 1) = i
        end do

        ! Singles (alpha only): cycle through occ/virt pairs
        do i = 1, n_occ
            do v = n_occ+1, n_mo
                if (d >= n_dets) exit
                d = d + 1
                do k = 1, n_occ
                    dets_a_bra(k, d, 1) = k
                    dets_b_bra(k, d, 1) = k
                    dets_a_ket(k, d, 1) = k
                    dets_b_ket(k, d, 1) = k
                end do
                dets_a_bra(i, d, 1) = v
                dets_a_ket(i, d, 1) = v
            end do
            if (d >= n_dets) exit
        end do

        ! Doubles (alpha+beta): cycle through a few combos
        do i = 1, n_occ
            do v = n_occ+1, n_mo
                if (d >= n_dets) exit
                d = d + 1
                do k = 1, n_occ
                    dets_a_bra(k, d, 1) = k
                    dets_b_bra(k, d, 1) = k
                    dets_a_ket(k, d, 1) = k
                    dets_b_ket(k, d, 1) = k
                end do
                dets_a_bra(i, d, 1) = v
                dets_a_ket(i, d, 1) = v
                dets_b_bra(mod(i, n_occ) + 1, d, 1) = n_mo - (v - n_occ - 1)
                dets_b_ket(mod(i, n_occ) + 1, d, 1) = n_mo - (v - n_occ - 1)
            end do
            if (d >= n_dets) exit
        end do

        ! Fill remaining by cycling
        do while (d < n_dets)
            d = d + 1
            do k = 1, n_occ
                dets_a_bra(k, d, 1) = dets_a_bra(k, d-1, 1)
                dets_b_bra(k, d, 1) = dets_b_bra(k, d-1, 1)
                dets_a_ket(k, d, 1) = dets_a_ket(k, d-1, 1)
                dets_b_ket(k, d, 1) = dets_b_ket(k, d-1, 1)
            end do
        end do

        ! State 2 = state 1
        do k = 1, n_dets
            dets_a_bra(:, k, 2) = dets_a_bra(:, k, 1)
            dets_b_bra(:, k, 2) = dets_b_bra(:, k, 1)
            dets_a_ket(:, k, 2) = dets_a_ket(:, k, 1)
            dets_b_ket(:, k, 2) = dets_b_ket(:, k, 1)
        end do

        ! Coefficients
        do i = 1, n_dets
            coeffs(i) = 1.0_dp / real(n_dets, dp)
        end do
        do k = 1, nstates
            coeffs_bra(:, k) = coeffs
            coeffs_ket(:, k) = coeffs
        end do

        ! Warm-up
        call cis_overlap_cisd_bruteforce(n_mo, n_occ, n_dets, nstates, &
                                          s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, &
                                          dets_a_ket, dets_b_ket, coeffs_bra, coeffs_ket, s_wf)

        ! Determine reps so total work is reasonable
        n_reps = 1

        ! Brute-force timing
        call cpu_time(t0)
        do trial = 1, n_reps
            call cis_overlap_cisd_bruteforce(n_mo, n_occ, n_dets, nstates, &
                                              s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, &
                                              dets_a_ket, dets_b_ket, coeffs_bra, coeffs_ket, s_wf)
        end do
        call cpu_time(t1)
        t_brute = (t1 - t0) / real(n_reps, dp) * 1000.0_dp

        ! L4M direct
        call cpu_time(t0)
        do trial = 1, n_reps
            call cis_overlap_cisd_l4m(n_mo, n_occ, n_dets, nstates, &
                                       s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, &
                                       dets_a_ket, dets_b_ket, coeffs_bra, coeffs_ket, s_wf)
        end do
        call cpu_time(t1)
        t_l4m = (t1 - t0) / real(n_reps, dp) * 1000.0_dp

        ! L4M fast
        call cpu_time(t0)
        do trial = 1, n_reps
            call cis_overlap_cisd_l4m_fast(n_mo, n_occ, n_dets, nstates, &
                                            s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, &
                                            dets_a_ket, dets_b_ket, coeffs_bra, coeffs_ket, s_wf)
        end do
        call cpu_time(t1)
        t_fast = (t1 - t0) / real(n_reps, dp) * 1000.0_dp

        print '(I8, A, F10.4, A, F10.4, A, F13.4, A, F8.2, A)', &
            n_dets, ' | ', t_brute, ' | ', t_l4m, ' | ', t_fast, ' | ', &
            t_brute / max(t_fast, 1.0e-6_dp), 'x'
    end do

    print '(A)', '  -------+------------+------------+---------------+---------'
end program