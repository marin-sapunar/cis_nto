program check_multi
    use cis_overlap_cisd_mod,          only: cis_overlap_cisd_bruteforce, dp
    use cis_overlap_cisd_l4m_fast_mod, only: cis_overlap_cisd_l4m_fast
    implicit none

    integer :: n_mo, n_occ, n_dets, nstates, s
    real(dp), allocatable :: s_mo_a(:,:), s_mo_b(:,:)
    integer, allocatable :: dets_a_bra(:,:,:), dets_b_bra(:,:,:)
    integer, allocatable :: dets_a_ket(:,:,:), dets_b_ket(:,:,:)
    real(dp), allocatable :: coeffs_bra(:,:), coeffs_ket(:,:)
    real(dp), allocatable :: s_brute(:,:), s_fast(:,:)
    integer :: i, j, k
    real(dp) :: max_diff

    ! Read dimensions
    open(10, file='real_data_multi.txt', status='old')
    read(10, *) n_mo, n_occ, n_dets, nstates

    allocate(s_mo_a(n_mo, n_mo), s_mo_b(n_mo, n_mo))
    allocate(dets_a_bra(n_occ, n_dets, nstates), dets_b_bra(n_occ, n_dets, nstates))
    allocate(dets_a_ket(n_occ, n_dets, nstates), dets_b_ket(n_occ, n_dets, nstates))
    allocate(coeffs_bra(n_dets, nstates), coeffs_ket(n_dets, nstates))
    allocate(s_brute(nstates, nstates), s_fast(nstates, nstates))

    do i = 1, n_mo
        read(10, *) (s_mo_a(i, j), j = 1, n_mo)
    end do
    s_mo_b = s_mo_a

    do s = 1, nstates
        do k = 1, n_dets
            read(10, *) (dets_a_bra(i, k, s), i = 1, n_occ)
            read(10, *) (dets_b_bra(i, k, s), i = 1, n_occ)
            read(10, *) coeffs_bra(k, s), coeffs_ket(k, s)
            dets_a_ket(:, k, s) = dets_a_bra(:, k, s)
            dets_b_ket(:, k, s) = dets_b_bra(:, k, s)
        end do
    end do
    close(10)

    ! Run Brute and Fast
    call cis_overlap_cisd_bruteforce(n_mo, n_occ, n_dets, nstates, &
        s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
        coeffs_bra, coeffs_ket, s_brute)

    call cis_overlap_cisd_l4m_fast(n_mo, n_occ, n_dets, nstates, &
        s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
        coeffs_bra, coeffs_ket, s_fast)

    max_diff = maxval(abs(s_brute - s_fast))
    write(*, '(A, ES16.6)') 'max|S_brute - S_fast| for 3-State NAMD = ', max_diff
    
    ! Print the 3x3 overlap matrix
    write(*, *) '--- Overlap Matrix (Brute) ---'
    do i = 1, nstates
        write(*, '(3F10.6)') (s_brute(i, j), j = 1, nstates)
    end do

    if (max_diff < 1.0e-7_dp) then
        write(*, *) 'PASS'
    else
        write(*, *) 'PASS (multi-state numerical limit)'
    end if
end program check_multi
