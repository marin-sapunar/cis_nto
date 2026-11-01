module cis_overlap_cisd_l4m_fast_mod
    implicit none
    private
    public :: cis_overlap_cisd_l4m_fast, dp
    integer, parameter :: dp = kind(1.0d0)
contains

    subroutine lu_det(A, n, det)
        integer, intent(in) :: n
        real(dp), intent(inout) :: A(n, n)
        real(dp), intent(out) :: det
        integer :: i, j, k, pivot
        real(dp) :: sign, temp, max_val
        det = 1.0_dp; sign = 1.0_dp
        do k = 1, n
            pivot = k; max_val = abs(A(k, k))
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

    subroutine inv_and_det(A_in, n, A_inv, det, ok)
        ! Hand-rolled LU inverse (cache-friendly for small n on this system)
        integer, intent(in) :: n
        real(dp), intent(in) :: A_in(n, n)
        real(dp), intent(out) :: A_inv(n, n), det
        logical, intent(out) :: ok
        real(dp) :: M(n, n), L(n, n), U(n, n), Y(n, n)
        integer :: i, j, k, pivot
        real(dp) :: factor, max_val, temp, sign, s
        M = A_in; L = 0.0_dp
        do i = 1, n
            L(i, i) = 1.0_dp
        end do
        sign = 1.0_dp; ok = .true.
        do k = 1, n
            pivot = k; max_val = abs(M(k, k))
            do i = k+1, n
                if (abs(M(i, k)) > max_val) then
                    max_val = abs(M(i, k)); pivot = i
                end if
            end do
            if (max_val < 1.0d-14) then
                ok = .false.; return
            end if
            if (pivot /= k) then
                do j = 1, n
                    temp = M(k, j); M(k, j) = M(pivot, j); M(pivot, j) = temp
                    temp = L(k, j); L(k, j) = L(pivot, j); L(pivot, j) = temp
                end do
                sign = -sign
            end if
            do i = k+1, n
                factor = M(i, k) / M(k, k)
                L(i, k) = factor
                M(i, k) = 0.0_dp
                do j = k+1, n
                    M(i, j) = M(i, j) - factor * M(k, j)
                end do
            end do
        end do
        U = M
        det = sign
        do k = 1, n
            det = det * U(k, k)
        end do
        do i = 1, n
            do j = 1, n
                s = 0.0_dp
                if (i == j) s = 1.0_dp
                do k = 1, i-1
                    s = s - L(i, k) * Y(k, j)
                end do
                Y(i, j) = s
            end do
        end do
        do i = n, 1, -1
            do j = 1, n
                s = Y(i, j)
                do k = i+1, n
                    s = s - U(i, k) * A_inv(k, j)
                end do
                A_inv(i, j) = s / U(i, i)
            end do
        end do
    end subroutine

    function det_from_occ_fast(bra_occ, ket_occ, s_mo, n_occ_ref, A_inv, det_A, Ainv_smo) result(det)
        integer, intent(in) :: n_occ_ref
        integer, intent(in) :: bra_occ(n_occ_ref), ket_occ(n_occ_ref)
        real(dp), intent(in) :: s_mo(:, :), A_inv(n_occ_ref, n_occ_ref), det_A, Ainv_smo(:,:)
        real(dp) :: det
        real(dp), allocatable :: M(:,:), sub(:,:)
        integer :: k_bra, k_ket, r, c, i, exc_r, exc_c
        real(dp) :: delta
        k_bra = 0; k_ket = 0; exc_r = 0; exc_c = 0
        do i = 1, n_occ_ref
            if (bra_occ(i) /= i) then
                k_bra = k_bra + 1; exc_r = i
            end if
            if (ket_occ(i) /= i) then
                k_ket = k_ket + 1; exc_c = i
            end if
        end do
        if (k_bra == 0 .and. k_ket == 0) then
            det = det_A
            return
        end if
        if (k_bra == 0 .and. k_ket == 1) then
            delta = 0.0_dp
            do r = 1, n_occ_ref
                delta = delta + A_inv(exc_c, r) * (s_mo(r, ket_occ(exc_c)) - s_mo(r, exc_c))
            end do
            det = det_A * (1.0_dp + delta)
            return
        end if
        if (k_bra == 1 .and. k_ket == 0) then
            delta = 0.0_dp
            do c = 1, n_occ_ref
                delta = delta + A_inv(c, exc_r) * (s_mo(bra_occ(exc_r), c) - s_mo(exc_r, c))
            end do
            det = det_A * (1.0_dp + delta)
            return
        end if

        ! (1,1): rank-2 Schur via Jacobi's formula
        ! det = det_A * det(I_2 + W), W[i,j] = v_i^T A_inv u_j
        ! where: u_1 = e_{r0}, v_1 = M_new[r0,:] - A_ref[r0,:]  (full row delta)
        !        u_2 = (0 except column c0), v_2 = e_{c0}
        if (k_bra == 1 .and. k_ket == 1) then
            block
                real(dp) :: v1(n_occ_ref), u2(n_occ_ref), Ainv_u2(n_occ_ref)
                real(dp) :: W11, W12, W21, W22, T(2,2), detT
                integer :: r, c
                ! v1 = M_new[r0, :] - A_ref[r0, :]  (row difference, full row)
                do c = 1, n_occ_ref
                    if (c == exc_c) then
                        v1(c) = s_mo(bra_occ(exc_r), ket_occ(exc_c)) - s_mo(exc_r, exc_c)
                    else
                        v1(c) = s_mo(bra_occ(exc_r), c) - s_mo(exc_r, c)
                    end if
                end do
                ! u2[r] = 0 at r0, = M_new[r, c0] - A_ref[r, c0] for r != r0
                u2 = 0.0_dp
                do r = 1, n_occ_ref
                    if (r == exc_r) cycle
                    u2(r) = s_mo(r, ket_occ(exc_c)) - s_mo(r, exc_c)
                end do
                ! Ainv_u2 = A_inv @ u2, with u2(exc_r) = 0
                ! = Ainv_smo(:, ket_c) - Ainv_smo(:, exc_c) - A_inv(:, exc_r) * u2(exc_r)
                do r = 1, n_occ_ref
                    Ainv_u2(r) = Ainv_smo(r, ket_occ(exc_c)) - Ainv_smo(r, exc_c) &
                               - A_inv(r, exc_r) * (s_mo(exc_r, ket_occ(exc_c)) - s_mo(exc_r, exc_c))
                end do
                ! W[1,1] = v1^T A_inv e_r0 = (A_inv[:,r0]) · v1 = sum over c of A_inv(c, r0)*v1(c)
                W11 = 0.0_dp
                do c = 1, n_occ_ref
                    W11 = W11 + A_inv(c, exc_r) * v1(c)
                end do
                ! W[1,2] = v1^T A_inv u2 = v1 · Ainv_u2
                W12 = 0.0_dp
                do r = 1, n_occ_ref
                    W12 = W12 + v1(r) * Ainv_u2(r)
                end do
                ! W[2,1] = e_c0^T A_inv e_r0 = A_inv(c0, r0)
                W21 = A_inv(exc_c, exc_r)
                ! W[2,2] = e_c0^T A_inv u2 = Ainv_u2(c0)
                W22 = Ainv_u2(exc_c)
                T(1,1) = 1.0_dp + W11
                T(1,2) = W12
                T(2,1) = W21
                T(2,2) = 1.0_dp + W22
                detT = T(1,1)*T(2,2) - T(1,2)*T(2,1)
                det = det_A * detT
            end block
            return
        end if
        allocate(M(n_occ_ref, n_occ_ref), sub(n_occ_ref, n_occ_ref))
        do r = 1, n_occ_ref
            do c = 1, n_occ_ref
                M(r, c) = s_mo(bra_occ(r), ket_occ(c))
            end do
        end do
        sub = M
        call lu_det(sub, n_occ_ref, det)
        deallocate(M, sub)
    end function

    subroutine cis_overlap_cisd_l4m_fast(n_mo, n_occ_ref, n_dets, nstates, &
                                          s_mo_a, s_mo_b, &
                                          dets_a_bra, dets_b_bra, &
                                          dets_a_ket, dets_b_ket, &
                                          coeffs_bra, coeffs_ket, s_wf)
        integer, intent(in) :: n_mo, n_occ_ref, n_dets, nstates
        real(dp), intent(in) :: s_mo_a(n_mo, n_mo), s_mo_b(n_mo, n_mo)
        integer, intent(in) :: dets_a_bra(n_occ_ref, n_dets, nstates)
        integer, intent(in) :: dets_b_bra(n_occ_ref, n_dets, nstates)
        integer, intent(in) :: dets_a_ket(n_occ_ref, n_dets, nstates)
        integer, intent(in) :: dets_b_ket(n_occ_ref, n_dets, nstates)
        real(dp), intent(in) :: coeffs_bra(n_dets, nstates), coeffs_ket(n_dets, nstates)
        real(dp), intent(out) :: s_wf(nstates, nstates)
        real(dp) :: A_a(n_occ_ref, n_occ_ref), A_b(n_occ_ref, n_occ_ref)
        real(dp) :: Ainv_smo_a(200, 800), Ainv_smo_b(200, 800)
        real(dp) :: A_inv_a(n_occ_ref, n_occ_ref), A_inv_b(n_occ_ref, n_occ_ref)
        real(dp) :: det_A_a, det_A_b, det_a, det_b, val
        logical :: ok
        integer :: sI, sJ, da, db, m1
        integer :: a1(n_occ_ref), b1(n_occ_ref), a2(n_occ_ref), b2(n_occ_ref)
        do m1 = 1, n_occ_ref
            A_a(m1, :) = s_mo_a(m1, 1:n_occ_ref)
            A_b(m1, :) = s_mo_b(m1, 1:n_occ_ref)
        end do
        call inv_and_det(A_a, n_occ_ref, A_inv_a, det_A_a, ok)
        if (.not. ok) return
        call inv_and_det(A_b, n_occ_ref, A_inv_b, det_A_b, ok)
        if (.not. ok) return

        ! Precompute A_inv * s_mo over the FULL MO space (n_mo columns)
        Ainv_smo_a(1:n_occ_ref, 1:n_mo) = matmul(A_inv_a, s_mo_a(1:n_occ_ref, 1:n_mo))
        Ainv_smo_b(1:n_occ_ref, 1:n_mo) = matmul(A_inv_b, s_mo_b(1:n_occ_ref, 1:n_mo))

        s_wf = 0.0_dp
        do sI = 1, nstates
            do sJ = 1, nstates
                val = 0.0_dp
                do da = 1, n_dets
                    if (coeffs_bra(da, sI) == 0.0_dp) cycle
                    a1 = dets_a_bra(:, da, sI)
                    b1 = dets_b_bra(:, da, sI)
                    do db = 1, n_dets
                        if (coeffs_ket(db, sJ) == 0.0_dp) cycle
                        a2 = dets_a_ket(:, db, sJ)
                        b2 = dets_b_ket(:, db, sJ)
                        det_a = det_from_occ_fast(a1, a2, s_mo_a, n_occ_ref, A_inv_a, det_A_a, Ainv_smo_a(1:n_occ_ref, 1:n_mo))
                        det_b = det_from_occ_fast(b1, b2, s_mo_b, n_occ_ref, A_inv_b, det_A_b, Ainv_smo_b(1:n_occ_ref, 1:n_mo))
                        val = val + coeffs_bra(da, sI) * coeffs_ket(db, sJ) * det_a * det_b
                    end do
                end do
                s_wf(sI, sJ) = val
            end do
        end do
    end subroutine

end module cis_overlap_cisd_l4m_fast_mod