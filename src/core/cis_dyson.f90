module cis_dyson_mod
    use global_defs
    use cis_nto_mod
    use cis_overlap_mod

    implicit none

contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_dyson
    !
    !> @brief Calculate Dyson orbitals between two sets of CIS wave functions.
    !> @details
    !
    !----------------------------------------------------------------------------------------------
    subroutine cis_dyson(trunc, s_mo, wf_a1, wf_a2, wf_b1, wf_b2, dys_mo)
        use blas95, only : gemv, gemm
        use matrix_mod, only : mat_ge_det
        real(dp), intent(in) :: trunc !< Threshold for truncating the wave functions.
        real(dp), intent(in) :: s_mo(:, :, :) !< Overlaps of alpha/beta molecular orbitals.
        real(dp), intent(in) :: wf_a1(:, :, :) !< Bra alpha wave function coefficients.
        real(dp), intent(in) :: wf_a2(:, :, :) !< Ket alpha wave function coefficients.
        real(dp), intent(in) :: wf_b1(:, :, :) !< Bra beta wave function coefficients.
        real(dp), intent(in) :: wf_b2(:, :, :) !< Ket beta wave function coefficients.
        real(dp), allocatable, intent(out) :: dys_mo(:, :, :) !< Dyson orbital in MO basis.
        integer :: n_1 !< Number of bra orbitals.
        integer :: n_2 !< Number of bra orbitals.
        integer :: no_a !< Number of occupied alpha orbitals.
        integer :: no_b1 !< Number of occupied beta orbitals.
        integer :: no_b2 !< Number of occupied beta orbitals.
        integer :: nwf_1 !< Number of bra states.
        integer :: nwf_2 !< Number of ket states.
        integer, allocatable :: na_a1(:) !< Number of active bra orbitals alpha.
        integer, allocatable :: na_a2(:) !< Number of active ket orbitals alpha.
        integer, allocatable :: na_b1(:) !< Number of active bra orbitals beta.
        integer, allocatable :: na_b2(:) !< Number of active ket orbitals beta.
        real(dp), allocatable :: c_a1(:, :) !< Bra NTO coefficients alpha.
        real(dp), allocatable :: c_a2(:, :) !< Ket NTO coefficients alpha.
        real(dp), allocatable :: c_b1(:, :) !< Bra NTO coefficients beta.
        real(dp), allocatable :: c_b2(:, :) !< Ket NTO coefficients beta.
        real(dp), allocatable :: mo_a1(:, :, :) !< Bra NTOs in MO basis alpha.
        real(dp), allocatable :: mo_a2(:, :, :) !< Ket NTOs in MO basis alpha.
        real(dp), allocatable :: mo_b1(:, :, :) !< Bra NTOs in MO basis beta.
        real(dp), allocatable :: mo_b2(:, :, :) !< Ket NTOs in MO basis beta.
        real(dp), allocatable :: wrk_a(:, :) !< Work array for alpha overlaps.
        real(dp), allocatable :: wrk_b(:, :) !< Work array for beta overlaps.
        real(dp), allocatable :: wrk0_a(:, :) !< Work array for alpha overlaps.
        real(dp), allocatable :: wrk0_b(:, :) !< Work array for beta overlaps.
        real(dp) :: rr_a !< RR block alpha.
        real(dp), allocatable :: sr_a(:) !< SR block alpha.
        real(dp), allocatable :: rs_a(:) !< RS block alpha.
        real(dp), allocatable :: ss_a(:, :) !< SS block alpha.
        real(dp), allocatable :: rr_b(:) !< RR Dyson block beta.
        real(dp), allocatable :: sr_b(:, :) !< SR Dyson block beta.
        real(dp), allocatable :: rs_b(:, :) !< RS Dyson block beta.
        real(dp), allocatable :: ss_b(:, :, :) !< SS Dyson block beta.
        real(dp), allocatable :: dys_nto(:, :, :) !< Dyson orbitals in NTO basis.
        integer :: i, j
 
        ! Get dimensions.
        n_1 = size(s_mo, 1)
        n_2 = size(s_mo, 2)
        no_a = size(wf_a2, 2)
        no_b1 = size(wf_b1, 2)
        no_b2 = size(wf_b2, 2)
        nwf_1 = size(wf_a1, 3)
        nwf_2 = size(wf_a2, 3)
        ! Allocate work arrays.
        allocate(sr_a(nwf_1), source = 0.0_dp)
        allocate(rs_a(nwf_2), source = 0.0_dp)
        allocate(ss_a(nwf_1, nwf_2), source = 0.0_dp)
        allocate(rr_b(2*no_b2), source = 0.0_dp)
        allocate(sr_b(2*no_b2, nwf_1), source = 0.0_dp)
        allocate(rs_b(2*no_b2, nwf_2), source = 0.0_dp)
        allocate(ss_b(2*no_b2, nwf_1, nwf_2), source = 0.0_dp)


        ! Calculate NTOs.
        call cis_nto(wf_a1, c_a1, mo_a1)
        call cis_nto(wf_a2, c_a2, mo_a2)
        call cis_nto(wf_b1, c_b1, mo_b1)
        call cis_nto(wf_b2, c_b2, mo_b2)

        ! Truncate wave functions.
        call cis_nto_truncate(.true., trunc, c_a1, c_b1, na_a1, na_b1)
        call cis_nto_truncate(.true., trunc, c_a2, c_b2, na_a2, na_b2)
        if (trunc < 1.0_dp) then
            write(stdout, '(a,f0.8)') 'Truncating wave functions based on threshold ', trunc
            write(stdout, '(a)') 'Number of remaining determinants for bra states:'
            write(stdout, '(3x,1000(i0,1x))') na_a1
            write(stdout, '(3x,1000(i0,1x))') na_b1
            write(stdout, '(a)') 'Number of remaining determinants for ket states:'
            write(stdout, '(3x,1000(i0,1x))') na_a2
            write(stdout, '(3x,1000(i0,1x))') na_b2
        end if

        ! Ref - Ref.
        rr_a = mat_ge_det(s_mo(1:no_a, 1:no_a, 1))
        call nto_rr_dys(no_b2, s_mo(:, :, 2), rr_b(:))

        allocate(wrk0_a(no_a*2, n_2))
        allocate(wrk0_b(no_b1*2, n_2))
        allocate(wrk_a(no_a*2, no_a*2))
        allocate(wrk_b(no_b1*2, no_b2*2))
        do i = 1, nwf_1
            ! CIS - Ref.
            call gemm(mo_a1(:, :, i), s_mo(:, :, 1), wrk0_a, transa='T')
            call gemm(mo_b1(:, :, i), s_mo(:, :, 2), wrk0_b, transa='T')
            call nto_rs(no_a, na_a1(i), wrk0_a, c_a1(:, i), sr_a(i), .true.)
            call nto_sr_dys(no_b2, na_b1(i), wrk0_b, c_b1(:, i), sr_b(:, i))
            ! CIS - CIS.
            do j = 1, nwf_2
                call gemm(wrk0_a, mo_a2(:, :, j), wrk_a)
                call gemm(wrk0_b, mo_b2(:, :, j), wrk_b)
                if (i == 1) then
                    call nto_rs(no_a, na_a2(j), wrk_a, c_a2(:, j), rs_a(j), .false.)
                    call nto_rs_dys(no_b2, na_b2(j), wrk_b, c_b2(:, j), rr_b(:), rs_b(:, j))
                end if
                call nto_ss(no_a, na_a1(i), na_a2(j), wrk_a, c_a1(:, i), c_a2(:, j), ss_a(i, j))
                call nto_ss_dys(no_b2, na_b1(i), na_b2(j), wrk_b, c_b1(:, i), c_b2(:, j), ss_b(:, i, j))
            end do
        end do
        deallocate(wrk_a)
        deallocate(wrk_b)
        deallocate(wrk0_a)
        deallocate(wrk0_b)

        allocate(dys_nto(2*no_b2, 0:nwf_1, 0:nwf_2))
        dys_nto(:, 0, 0) = rr_b * rr_a
        do i = 1, nwf_1
            dys_nto(:, i, 0) = rr_b * sr_a(i) + sr_b(:, i) * rr_a
        end do
        do j = 1, nwf_2
            dys_nto(:, 0, j) = rr_b * rs_a(j) + rs_b(:, j) * rr_a
        end do
        do i = 1, nwf_1
            do j = 1, nwf_2
                dys_nto(:, i, j) = rr_b * ss_a(i, j) + rs_b(:, j) * sr_a(i) + &
                                   ss_b(:, i, j) * rr_a + sr_b(:, i) * rs_a(j)
            end do
        end do

        if (allocated(dys_mo)) deallocate(dys_mo)
        allocate(dys_mo(n_2, nwf_1+1, nwf_2+1), source = 0.0_dp)
        dys_mo(1:no_b2, :, 1) = dys_nto(1:no_b2, :, 1)
        do j = 1, nwf_2 
            call gemm(mo_b2(:, :, j), dys_nto(:, :, j+1), dys_mo(:, :, j+1), alpha=sqrt(2.0_dp))
        end do
    end subroutine cis_dyson


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: nto_rr_dys
    !----------------------------------------------------------------------------------------------
    subroutine nto_rr_dys(nn, s_nto, rr)
        use matrix_mod, only : mat_ge_det
        integer, intent(in) :: nn !< Number of occupied orbitals in neutral.
        real(dp), intent(in) :: s_nto(:, :) !< Overlap matrix between NTOs.
        real(dp), intent(out) :: rr(2*nn) !< RR block coefficients.
        real(dp) :: wrk(nn-1, nn-1) !< Work array.
        integer :: seq(nn) !< All columns.
        integer :: cols(nn-1) !< Currently active columns.
        integer :: sgn !< Sign of the current minor.
        integer :: nc !< Number of occupied orbitals in cation.
        integer :: i

        seq = [ (i, i=1, nn) ]
        nc = nn - 1

        !$omp parallel default(shared)
        !$omp do private(wrk, cols, sgn) schedule(dynamic)
        do i = 1, nn
            sgn = (-1)**(i+1)
            cols = pack(seq, (seq /= i))
            wrk = s_nto(1:nc, cols)
            rr(i) = sgn * mat_ge_det(wrk)
        end do
        !$omp end do
        !$omp end parallel
    end subroutine nto_rr_dys


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: nto_sr_dys
    !----------------------------------------------------------------------------------------------
    subroutine nto_sr_dys(nn, na, s_nto, c, sr)
        use matrix_mod, only : mat_ge_det
        integer, intent(in) :: nn !< Number of occupied orbitals in neutral.
        integer, intent(in) :: na !< Number of active excitations.
        real(dp), intent(in) :: s_nto(:, :) !< Overlap matrix between NTOs.
        real(dp), intent(in) :: c(nn-1) !< Coefficients of excitations.
        real(dp), intent(out) :: sr(2*nn) !< RS block coefficients.
        real(dp) :: wrk(nn-1, nn-1) !< Work array.
        integer :: seq(nn) !< All columns.
        integer :: cols(nn-1) !< Currently active columns.
        integer :: sgn !< Sign of the current minor.
        integer :: nc !< Number of occupied orbitals in cation.
        integer :: i, j

        if (na == 0) return
        nc = nn - 1
        seq = [ (i, i=1, nn+1) ]

        !$omp parallel default(shared)
        !$omp do private(wrk, j, cols, sgn) schedule(dynamic)
        do i = 1, nn
            sgn = (-1)**(i+1)
            cols = pack(seq, (seq /= i))
            do j = 1, na
                wrk = s_nto(1:nc, cols)
                wrk(j, :) = s_nto(nc+j, cols)
                sr(i) = sr(i) + sgn * c(j) * mat_ge_det(wrk)
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine nto_sr_dys


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: nto_rs_dys
    !----------------------------------------------------------------------------------------------
    subroutine nto_rs_dys(nn, na, s_nto, c, rr, rs)
        use matrix_mod, only : mat_ge_det
        integer, intent(in) :: nn !< Number of occupied orbitals in neutral.
        integer, intent(in) :: na !< Number of active excitations.
        real(dp), intent(in) :: s_nto(:, :) !< Overlap matrix between NTOs.
        real(dp), intent(in) :: c(nn) !< Coefficients of excitations.
        real(dp), intent(in) :: rr(2*nn) !< RR block coefficients.
        real(dp), intent(out) :: rs(2*nn) !< RS block coefficients.
        real(dp) :: wrk(nn-1, nn-1) !< Work array.
        integer :: seq(nn) !< All columns.
        integer :: cols(nn-1) !< Currently active columns.
        integer :: sgn !< Sign of the current minor.
        integer :: nc !< Number of occupied orbitals in cation.
        integer :: i, j

        if (na == 0) return
        nc = nn - 1
        seq = [ (i, i=1, nn+1) ]

        !$omp parallel default(shared)
        !$omp do private(wrk, j, cols, sgn) schedule(dynamic)
        do i = 1, nn
            sgn = (-1)**(i+1)
            cols = pack(seq, (seq /= i))
            rs(nn+i) = rs(nn+i) + c(i) * rr(i)
            do j = 1, na
                if (j == i) cycle
                wrk = s_nto(1:nc, cols)
                if (j < i) then
                    wrk(:, j) = s_nto(1:nc, nn+j)
                else
                    wrk(:, j-1) = s_nto(1:nc, nn+j)
                end if
                rs(i) = rs(i) + sgn * c(j) * mat_ge_det(wrk)
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine nto_rs_dys


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: nto_ss_dys
    !----------------------------------------------------------------------------------------------
    subroutine nto_ss_dys(nn, na1, na2, s_nto, c1, c2, ss)
        use matrix_mod, only : mat_ge_det
        use blas95, only : gemv, dot
        integer, intent(in) :: nn !< Number of occupied orbitals in neutral.
        integer, intent(in) :: na1 !< Number of active bra excitations.
        integer, intent(in) :: na2 !< Number of active ket excitations.
        real(dp), intent(in) :: s_nto(2*nn-2, 2*nn) !< Overlap matrix between NTOs.
        real(dp), intent(in) :: c1(nn-1) !< Bra excitation coefficients.
        real(dp), intent(in) :: c2(nn) !< Ket excitation coefficients.
        real(dp), intent(out) :: ss(2*nn) !< SS block sum.
        real(dp) :: wrk(nn-1, nn-1) !< Work array.
        real(dp) :: wrk_c(nn-1, nn-1) !< Work array.
        integer :: seq(nn) !< All columns.
        integer :: cols(nn-1) !< Currently active columns.
        integer :: sgn !< Sign of the current minor.
        integer :: nc !< Number of occupied orbitals in cation.
        integer :: i, j, k

        if ((na1 == 0) .or. (na2 == 0)) return
        nc = nn - 1
        seq = [ (i, i=1, nn+1) ]

        !$omp parallel default(shared)
        !$omp do private(wrk, wrk_c, j, k, cols, sgn) schedule(dynamic)
        do i = 1, nn
            sgn = (-1)**(i+1)
            cols = pack(seq, (seq /= i))
            do j = 1, na1
                wrk_c = s_nto(1:nc, cols)
                wrk_c(j, :) = s_nto(nc+j, cols)
                ss(nn+i) = ss(nn+i) + sgn * c1(j) * c2(i) * mat_ge_det(wrk_c)
                do k = 1, na2
                    if (k == i) cycle
                    wrk = wrk_c
                    if (k < i) then
                        wrk(:, k) = s_nto(1:nc, nn+k)
                        wrk(j, k) = s_nto(nc+j, nn+k)
                    else
                        wrk(:, k-1) = s_nto(1:nc, nn+k)
                        wrk(j, k-1) = s_nto(nc+j, nn+k)
                    end if
                    ss(i) = ss(i) + sgn * c2(k) * c1(j) * mat_ge_det(wrk)
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine nto_ss_dys


end module cis_dyson_mod
