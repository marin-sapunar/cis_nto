program check_trunc
    use cis_overlap_cisd_mod,          only: cis_overlap_cisd_bruteforce, dp
    use cis_overlap_cisd_l4m_fast_mod, only: cis_overlap_cisd_l4m_fast
    implicit none

    integer :: n_mo, n_occ, n_dets, nstates
    real(dp), allocatable :: s_mo_a(:,:), s_mo_b(:,:)
    integer, allocatable :: dets_a_bra(:,:,:), dets_b_bra(:,:,:)
    integer, allocatable :: dets_a_ket(:,:,:), dets_b_ket(:,:,:)
    real(dp), allocatable :: coeffs_bra(:,:), coeffs_ket(:,:)
    real(dp) :: s_brute(1,1), s_fast(1,1)
    character(len=256) :: filename
    integer :: i, j, k
    real(dp) :: diff

    call get_command_argument(1, filename)
    
    open(10, file=trim(filename), status='old')
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

    call cis_overlap_cisd_bruteforce(n_mo, n_occ, n_dets, 1, &
        s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
        coeffs_bra, coeffs_ket, s_brute)

    call cis_overlap_cisd_l4m_fast(n_mo, n_occ, n_dets, 1, &
        s_mo_a, s_mo_b, dets_a_bra, dets_b_bra, dets_a_ket, dets_b_ket, &
        coeffs_bra, coeffs_ket, s_fast)

    diff = abs(s_brute(1,1) - s_fast(1,1))
    write(*, '(A, I4, A, ES16.6, A, ES16.6)') &
        'N_dets=', n_dets, '  S=', s_brute(1,1), '  |Brute-Fast|=', diff
end program check_trunc
