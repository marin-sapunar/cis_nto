module cis_overlap_cisd_mod
    ! Brute-force CISD overlap calculator (reference implementation).
    ! This is the Fortran counterpart of test_cisd/bruteforce_cisd.py.
    ! Slower than the L2M-based CIS code, but handles doubles exactly.
    implicit none
    private
    public :: cis_overlap_cisd_bruteforce
    public :: dp

    integer, parameter :: dp = kind(1.0d0)

contains

    subroutine lu_det(A, n, det)
        integer, intent(in) :: n
        real(dp), intent(inout) :: A(n, n)
        real(dp), intent(out) :: det
        integer :: i, j, k, pivot
        real(dp) :: sign, temp, max_val
        det = 1.0_dp
        sign = 1.0_dp
        do k = 1, n
            pivot = k
            max_val = abs(A(k, k))
            do i = k+1, n
                if (abs(A(i, k)) > max_val) then
                    max_val = abs(A(i, k)); pivot = i
                end if
            end do
            if (max_val < 1.0d-14) then
                det = 0.0_dp; return
            end if
            if (pivot /= k) then
                do j = 1, n
                    temp = A(k, j); A(k, j) = A(pivot, j); A(pivot, j) = temp
                end do
                sign = -sign
            end if
            det = det * A(k, k)
            do i = k+1, n
                temp = A(i, k) / A(k, k)
                do j = k+1, n
                    A(i, j) = A(i, j) - temp * A(k, j)
                end do
            end do
        end do
        det = det * sign
    end subroutine

    subroutine det_overlap(n_occ, a1, b1, a2, b2, s_mo_a, s_mo_b, val)
        ! Overlap of two Slater determinants:
        ! <det(a1,b1)|det(a2,b2)> = det(S_MO_a[a1,a2]) * det(S_MO_b[b1,b2])
        integer, intent(in) :: n_occ
        integer, intent(in) :: a1(n_occ), b1(n_occ), a2(n_occ), b2(n_occ)
        real(dp), intent(in) :: s_mo_a(:, :), s_mo_b(:, :)
        real(dp), intent(out) :: val
        real(dp) :: det_a, det_b
        real(dp), allocatable :: sub_a(:, :), sub_b(:, :)
        integer :: i, j
        allocate(sub_a(n_occ, n_occ), sub_b(n_occ, n_occ))
        do i = 1, n_occ
            do j = 1, n_occ
                sub_a(i, j) = s_mo_a(a1(i), a2(j))
                sub_b(i, j) = s_mo_b(b1(i), b2(j))
            end do
        end do
        call lu_det(sub_a, n_occ, det_a)
        call lu_det(sub_b, n_occ, det_b)
        val = det_a * det_b
    end subroutine

    subroutine cis_overlap_cisd_bruteforce(n_mo, n_occ, n_dets, nstates, &
                                          s_mo_a, s_mo_b, &
                                          dets_a_bra, dets_b_bra, &
                                          dets_a_ket, dets_b_ket, &
                                          coeffs_bra, coeffs_ket, s_wf)
        ! Inputs:
        !   s_mo_a, s_mo_b:  (n_mo, n_mo) MO overlaps (alpha, beta)
        !   dets_a_bra, dets_a_ket: (n_occ, n_dets, nstates) occupied indices
        !   dets_b_bra, dets_b_ket: same for beta
        !   coeffs_bra, coeffs_ket: (n_dets, nstates) CI coefficients
        ! Output:
        !   s_wf: (nstates, nstates) overlap matrix
        integer, intent(in) :: n_mo, n_occ, n_dets, nstates
        real(dp), intent(in) :: s_mo_a(n_mo, n_mo), s_mo_b(n_mo, n_mo)
        integer, intent(in) :: dets_a_bra(n_occ, n_dets, nstates)
        integer, intent(in) :: dets_b_bra(n_occ, n_dets, nstates)
        integer, intent(in) :: dets_a_ket(n_occ, n_dets, nstates)
        integer, intent(in) :: dets_b_ket(n_occ, n_dets, nstates)
        real(dp), intent(in) :: coeffs_bra(n_dets, nstates)
        real(dp), intent(in) :: coeffs_ket(n_dets, nstates)
        real(dp), intent(out) :: s_wf(nstates, nstates)
        integer :: I, J, da, db
        integer :: a1(n_occ), b1(n_occ), a2(n_occ), b2(n_occ)
        real(dp) :: val, ov
        do I = 1, nstates
            do J = 1, nstates
                val = 0.0_dp
                do da = 1, n_dets
                    a1 = dets_a_bra(:, da, I)
                    b1 = dets_b_bra(:, da, I)
                    do db = 1, n_dets
                        a2 = dets_a_ket(:, db, J)
                        b2 = dets_b_ket(:, db, J)
                        call det_overlap(n_occ, a1, b1, a2, b2, s_mo_a, s_mo_b, ov)
                        val = val + coeffs_bra(da, I) * coeffs_ket(db, J) * ov
                    end do
                end do
                s_wf(I, J) = val
            end do
        end do
    end subroutine

end module cis_overlap_cisd_mod
