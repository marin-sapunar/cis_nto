module cis_overlap_mod
    use global_defs
    use cis_nto_mod

    implicit none

contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_overlap
    !
    !> @brief Calculate overlap matrix between two sets of CIS wave functions.
    !> @details
    !! Transform the wave functions into expansion in terms of excitations between NTOs before
    !! calculating the overlaps. If a threshold trunc between 0 and 1 is given, the wave function 
    !! expansions are truncated to include only the dominant contributions until the square of the
    !! norm for each state is higher than the threshold.
    !
    !> @note The first row/column of the output matrix s_wf are the overlaps between the reference
    !! and the CIS states.
    !> @note The rr, sr, rs and ss arrays, hold the results of the sums for each pair of states.
    !! The RR block is the same for all pairs of states, and the SR and RS blocks are the same for
    !! each ket and bra state, respectively. However, the sign of the determinants can get fliped 
    !! by the SVD procedure so we calculate and store them for each pair of states so we don't need
    !! to keep track of the sign.
    !----------------------------------------------------------------------------------------------
    subroutine cis_overlap(trunc, s_mo, wf_a1, wf_a2, wf_b1, wf_b2, s_wf)
        use blas95, only : gemm
        use matrix_mod, only : mat_ge_det
        real(dp), intent(in) :: trunc !< Threshold for truncating the wave functions.
        real(dp), intent(in) :: s_mo(:, :, :) !< Overlaps of alpha/beta molecular orbitals.
        real(dp), intent(in) :: wf_a1(:, :, :) !< Bra alpha wave function coefficients.
        real(dp), intent(in) :: wf_a2(:, :, :) !< Ket alpha wave function coefficients.
        real(dp), intent(in) :: wf_b1(:, :, :) !< Bra beta wave function coefficients.
        real(dp), intent(in) :: wf_b2(:, :, :) !< Ket beta wave function coefficients.
        real(dp), allocatable, intent(out) :: s_wf(:, :) !< Wave function overlaps.
        logical :: beta !< Restricted/unrestricted calculation.
        integer :: n_2 !< Number of bra orbitals.
        integer :: no_a !< Number of occupied alpha orbitals.
        integer :: no_b !< Number of occupied beta orbitals.
        integer :: nwf_1 !< Number of bra states.
        integer :: nwf_2 !< Number of ket states.
        real(dp), allocatable :: c_a1(:, :) !< Bra NTO coefficients alpha.
        real(dp), allocatable :: c_a2(:, :) !< Ket NTO coefficients alpha.
        real(dp), allocatable :: c_b1(:, :) !< Bra NTO coefficients beta.
        real(dp), allocatable :: c_b2(:, :) !< Ket NTO coefficients beta.
        real(dp), allocatable :: mo_a1(:, :, :) !< Bra NTOs in MO basis alpha.
        real(dp), allocatable :: mo_a2(:, :, :) !< Ket NTOs in MO basis alpha.
        real(dp), allocatable :: mo_b1(:, :, :) !< Bra NTOs in MO basis beta.
        real(dp), allocatable :: mo_b2(:, :, :) !< Ket NTOs in MO basis beta.
        real(dp), allocatable :: s_nto_a(:, :) !< NTO overlaps alpha.
        real(dp), allocatable :: s_nto_b(:, :) !< NTO overlaps beta.
        real(dp), allocatable :: wrk_a(:, :) !< Work array for alpha overlaps.
        real(dp), allocatable :: wrk_b(:, :) !< Work array for beta overlaps.
        real(dp), allocatable :: rr_a(:, :) !< RR block alpha.
        real(dp), allocatable :: rr_b(:, :) !< RR block beta.
        real(dp), allocatable :: sr_a(:, :) !< SR block alpha.
        real(dp), allocatable :: sr_b(:, :) !< SR block beta.
        real(dp), allocatable :: rs_a(:, :) !< RS block alpha.
        real(dp), allocatable :: rs_b(:, :) !< RS block beta.
        real(dp), allocatable :: ss_a(:, :) !< SS block alpha.
        real(dp), allocatable :: ss_b(:, :) !< SS block beta.
        integer, allocatable :: na_a1(:) !< Number of active bra orbitals alpha.
        integer, allocatable :: na_a2(:) !< Number of active ket orbitals alpha.
        integer, allocatable :: na_b1(:) !< Number of active bra orbitals beta.
        integer, allocatable :: na_b2(:) !< Number of active ket orbitals beta.
        integer :: i, j

        ! Get dimensions.
        beta = .false.
        if (size(s_mo, 3) == 2) beta = .true.
        n_2 = size(s_mo, 2)
        no_a = size(wf_a1, 2)
        if (beta) no_b = size(wf_b1, 2)
        nwf_1 = size(wf_a1, 3)
        nwf_2 = size(wf_a2, 3)
        ! Allocate work arrays.
        allocate(s_nto_a(no_a*2, no_a*2))
        allocate(rr_a(nwf_1, nwf_2))
        allocate(sr_a(nwf_1, nwf_2))
        allocate(rs_a(nwf_1, nwf_2))
        allocate(ss_a(nwf_1, nwf_2))
        allocate(wrk_a(no_a*2, n_2))
        if (beta) then
            allocate(s_nto_b(no_b*2, no_b*2))
            allocate(rr_b(nwf_1, nwf_2))
            allocate(sr_b(nwf_1, nwf_2))
            allocate(rs_b(nwf_1, nwf_2))
            allocate(ss_b(nwf_1, nwf_2))
            allocate(wrk_b(no_b*2, n_2))
        end if

        ! Calculate NTOs.
        call cis_nto(wf_a1, c_a1, mo_a1)
        call cis_nto(wf_a2, c_a2, mo_a2)
        if (beta) then
            call cis_nto(wf_b1, c_b1, mo_b1)
            call cis_nto(wf_b2, c_b2, mo_b2)
        end if

        ! Truncate wave functions:
        call cis_nto_truncate(beta, trunc, c_a1, c_b1, na_a1, na_b1)
        call cis_nto_truncate(beta, trunc, c_a2, c_b2, na_a2, na_b2)
        if (trunc < 1.0_dp) then
            write(stdout, '(a,f0.8)') 'Truncating wave functions based on threshold ', trunc
            write(stdout, '(a)') 'Number of remaining determinants for bra states:'
            write(stdout, '(3x,1000(i0,1x))') na_a1
            if (beta) write(stdout, '(3x,1000(i0,1x))') na_b1
            write(stdout, '(a)') 'Number of remaining determinants for ket states:'
            write(stdout, '(3x,1000(i0,1x))') na_a2
            if (beta) write(stdout, '(3x,1000(i0,1x))') na_b2
        end if

        ! Calculate alpha determinant blocks.
        do i = 1, nwf_1
            call gemm(mo_a1(:, :, i), s_mo(:, :, 1), wrk_a, transa='T')
            do j = 1, nwf_2
                call gemm(wrk_a, mo_a2(:, :, j), s_nto_a)
                rr_a(i, j) = mat_ge_det(s_nto_a(1:no_a, 1:no_a))
                call nto_rs(no_a, na_a1(i), s_nto_a, c_a1(:, i), sr_a(i, j), .true.)
                call nto_rs(no_a, na_a2(j), s_nto_a, c_a2(:, j), rs_a(i, j), .false.)
                call nto_ss(no_a, na_a1(i), na_a2(j), s_nto_a, c_a1(:, i), c_a2(:, j), ss_a(i, j))
            end do
        end do

        ! Calculate beta determinant blocks.
        if (beta) then
            do i = 1, nwf_1
                call gemm(mo_b1(:, :, i), s_mo(:, :, 2), wrk_b, transa='T')
                do j = 1, nwf_2
                    call gemm(wrk_b, mo_b2(:, :, j), s_nto_b)
                    rr_b(i, j) = mat_ge_det(s_nto_b(1:no_b, 1:no_b))
                    call nto_rs(no_b, na_b1(i), s_nto_b, c_b1(:, i), sr_b(i, j), .true.)
                    call nto_rs(no_b, na_b2(j), s_nto_b, c_b2(:, j), rs_b(i, j), .false.)
                    call nto_ss(no_b, na_b1(i), na_b2(j), s_nto_b, c_b1(:, i), c_b2(:, j), ss_b(i, j))
                end do
            end do
        end if

        if (allocated(s_wf)) deallocate(s_wf)
        allocate(s_wf(nwf_1+1, nwf_2+1))

        if (beta) then
            ! Ground state overlap:
            s_wf(1, 1) = rr_a(1, 1) * rr_b(1, 1)
            ! Ground-excited state overlap:
            s_wf(2:, 1) = rr_a(:, 1) * sr_b(:, 1) + rr_b(:, 1) * sr_a(:, 1)
            ! Excited-ground state overlap:
            s_wf(1, 2:) = rr_b(1, :) * rs_a(1, :) + rr_a(1, :) * rs_b(1, :)
            ! Excited-excited state overlap:
            s_wf(2:, 2:) = rr_b * ss_a + sr_a * rs_b + rr_a * ss_b + sr_b * rs_a
        else
            ! Ground state overlap:
            s_wf(1, 1) = rr_a(1, 1) * rr_a(1, 1)
            ! Ground-excited state overlap:
            s_wf(2:, 1) = rr_a(:, 1) * sr_a(:, 1)
            ! Excited-ground state overlap:
            s_wf(1, 2:) = rr_a(1, :) * rs_a(1, :)
            ! Excited-excited state overlap:
            s_wf(2:, 2:) = rr_a * ss_a + sr_a * rs_a
        end if
    end subroutine cis_overlap


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: nto_rs
    !> @brief Calculate the RS or SR block between two states in a cis overlap calculation.
    !> @note Row is .false. for RS block calculations and .true. for SR block calculations.
    !----------------------------------------------------------------------------------------------
    subroutine nto_rs(n, na, s_nto, c, rs, row)
        use matrix_mod, only : mat_ge_det
        use blas95, only : dot
        integer, intent(in) :: n !< Number of occupied orbitals.
        integer, intent(in) :: na !< Number of active excitations.
        real(dp), intent(in) :: s_nto(:, :) !< Overlap matrix between NTOs.
        real(dp), intent(in) :: c(n) !< Coefficients of excitations.
        real(dp), intent(out) :: rs !< RS/SR block sum.
        logical, intent(in) :: row !< RS or SR block.
        real(dp) :: dets(na) !< Determinants.
        real(dp) :: wrk(n, n) !< Work array.
        integer :: i

        rs = 0.0_dp
        if (na == 0) return

        !$omp parallel
        !$omp do private(wrk)
        do i = 1, na
            wrk = s_nto(1:n, 1:n)
            if (row) then
                wrk(i, :) = s_nto(n+i, 1:n)
            else
                wrk(:, i) = s_nto(1:n, n+i)
            end if
            dets(i) = mat_ge_det(wrk)
        end do
        !$omp end do
        !$omp end parallel

        rs = dot(c(1:na), dets)
    end subroutine nto_rs


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: nto_ss
    !> @brief Calculate the SS block between two states in a cis overlap calculation.
    !----------------------------------------------------------------------------------------------
    subroutine nto_ss(n, na1, na2, s_nto, c1, c2, ss)
        use matrix_mod, only : mat_ge_det
        use blas95, only : gemv, dot
        integer, intent(in) :: n !< Number of occupied orbitals.
        integer, intent(in) :: na1 !< Number of active bra excitations.
        integer, intent(in) :: na2 !< Number of active ket excitations.
        real(dp), intent(in) :: s_nto(2*n, 2*n) !< Overlap matrix between NTOs.
        real(dp), intent(in) :: c1(n) !< Bra excitation coefficients.
        real(dp), intent(in) :: c2(n) !< Ket excitation coefficients.
        real(dp), intent(out) :: ss !< SS block sum.
        real(dp) :: dets(na1, na2) !< Determinants.
        real(dp) :: wrk(n, n) !< Work array.
        real(dp) :: vec(na1) !< Work vector.
        integer :: i, j

        ss = 0.0_dp
        if ((na1 == 0) .or. (na2 == 0)) return

        !$omp parallel default(shared)
        !$omp do private(wrk, j) schedule(dynamic)
        do i = 1, na1
            do j = 1, na2
                wrk = s_nto(1:n, 1:n)
                wrk(i, :) = s_nto(n+i, 1:n)
                wrk(:, j) = s_nto(1:n, n+j)
                wrk(i, j) = s_nto(n+i, n+j)
                dets(i, j) = mat_ge_det(wrk)
            end do
        end do
        !$omp end do
        !$omp end parallel

        call gemv(dets, c2(1:na2), vec)
        ss = dot(c1(1:na1), vec)
    end subroutine nto_ss


end module cis_overlap_mod
