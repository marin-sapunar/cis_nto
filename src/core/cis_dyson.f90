!----------------------------------------------------------------------------------------------
! MODULE: cis_dyson_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date October, 2018
!
!> @brief Subroutines for calculating Dyson orbitals from CIS wave functions.
!----------------------------------------------------------------------------------------------
module cis_dyson_mod
    use global_defs
    use cis_nto_mod
    use cis_overlap_nto_mod
    implicit none


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_dyson
    !
    !> @brief Calculate Dyson orbitals between two sets of CIS wave functions.
    !> @details
    !! Natural transition orbitals are calculated for the CIS wave functions. Determinants are
    !! calculated in this basis and converted back to the canonical MO basis.
    !
    !> @note The Dyson orbitals between the N-1 el. reference and N el. CIS states are returned as
    !! dys_mo(:, 0, :). The Dyson orbitals between the N-1 el. CIS states and N el. reference are
    !! returned as dys_mo(:, :, 0).
    !
    !> @note The rr, sr, rs and ss blocks for alpha spin are just the corresponding overlap blocks
    !! as the number of alpha electrons is the same.
    !----------------------------------------------------------------------------------------------
    subroutine cis_dyson(trunc, s_mo_a, s_mo_b, wf_a1, wf_a2, wf_b1, wf_b2, dys_mo)
        use blas95, only : gemv, gemm
        use matrix_mod, only : mat_ge_det
        real(dp), intent(in) :: trunc !< Threshold for truncating the wave functions.
        real(dp), intent(in) :: s_mo_a(:, :) !< Overlaps of alpha molecular orbitals.
        real(dp), intent(in) :: s_mo_b(:, :) !< Overlaps of beta molecular orbitals.
        real(dp), intent(in) :: wf_a1(:, :, :) !< Bra alpha wave function coefficients.
        real(dp), intent(in) :: wf_a2(:, :, :) !< Ket alpha wave function coefficients.
        real(dp), intent(in) :: wf_b1(:, :, :) !< Bra beta wave function coefficients.
        real(dp), intent(in) :: wf_b2(:, :, :) !< Ket beta wave function coefficients.
        real(dp), allocatable, intent(out) :: dys_mo(:, :, :) !< Dyson orbital in MO basis.
        integer :: n_a2 !< Number of bra alpha orbitals.
        integer :: n_b2 !< Number of bra beta orbitals.
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
        real(dp), allocatable :: wrk_nto(:) !< Work array for blocks in NTO basis.
        real(dp) :: rr_a !< RR block alpha.
        real(dp), allocatable :: sr_a(:) !< SR block alpha.
        real(dp), allocatable :: rs_a(:) !< RS block alpha.
        real(dp), allocatable :: ss_a(:, :) !< SS block alpha.
        real(dp), allocatable :: rr_b(:) !< RR Dyson block beta.
        real(dp), allocatable :: sr_b(:, :) !< SR Dyson block beta.
        real(dp), allocatable :: rs_b(:, :) !< RS Dyson block beta.
        real(dp), allocatable :: ss_b(:, :, :) !< SS Dyson block beta.
        integer :: i, j
        real(dp), external :: omp_get_wtime
        real(dp) :: time00, time0
        real(dp) :: time_nto, time_det, time_tot
 
        time00 = omp_get_wtime()
        if (print_level >= 2) then
            write(stdout, *)
            write(stdout, '(5x,a)') '---- start cis_dyson subroutine ----'
        end if

        ! Get dimensions.
        n_a2 = size(s_mo_a, 2)
        n_b2 = size(s_mo_b, 2)
        no_a = size(wf_a2, 2)
        no_b1 = size(wf_b1, 2)
        no_b2 = size(wf_b2, 2)
        nwf_1 = size(wf_a1, 3)
        nwf_2 = size(wf_a2, 3)
        ! Allocate work arrays.
        allocate(sr_a(nwf_1), source = 0.0_dp)
        allocate(rs_a(nwf_2), source = 0.0_dp)
        allocate(ss_a(nwf_1, nwf_2), source = 0.0_dp)
        allocate(rr_b(n_b2), source = 0.0_dp)
        allocate(sr_b(n_b2, nwf_1), source = 0.0_dp)
        allocate(rs_b(n_b2, nwf_2), source = 0.0_dp)
        allocate(ss_b(n_b2, nwf_1, nwf_2), source = 0.0_dp)

        ! Calculate NTOs.
        time0 = omp_get_wtime()
        if (print_level >= 2) then
            write(stdout, *)
            write(stdout, '(5x,a)') 'Generating NTOs...'
        end if
        call cis_nto(wf_a1, c_a1, mo_a1)
        call cis_nto(wf_a2, c_a2, mo_a2)
        call cis_nto(wf_b1, c_b1, mo_b1)
        call cis_nto(wf_b2, c_b2, mo_b2)
        time_nto = omp_get_wtime() - time0

        ! Truncate wave functions.
        call cis_nto_truncate(.true., trunc, c_a1, c_b1, na_a1, na_b1)
        call cis_nto_truncate(.true., trunc, c_a2, c_b2, na_a2, na_b2)
        if (trunc < 1.0_dp) then
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(5x,a,f0.8)') 'Truncating wave functions based on threshold ', trunc
                write(stdout, '(5x,a)') 'Number of remaining determinants for N-1 el. states:'
                write(stdout, '(9x,1000(i0,1x))') na_a1
                write(stdout, '(9x,1000(i0,1x))') na_b1
                write(stdout, '(5x,a)') 'Number of remaining determinants for N el. states:'
                write(stdout, '(9x,1000(i0,1x))') na_a2
                write(stdout, '(9x,1000(i0,1x))') na_b2
            end if
        end if

        ! Ref - Ref.
        time0 = omp_get_wtime()
        if (print_level >= 2) then
            write(stdout, *)
            write(stdout, '(5x,a)') 'Computing determinant blocks...'
            write(stdout, '(5x,a)') 'Status:'
        end if
        rr_a = mat_ge_det(s_mo_a(1:no_a, 1:no_a))
        call nto_rr_dys(no_b2, s_mo_b, rr_b(1:no_b2))

        allocate(wrk0_a(no_a*2, n_a2))
        allocate(wrk0_b(no_b1*2, n_b2))
        allocate(wrk_a(no_a*2, no_a*2))
        allocate(wrk_b(no_b1*2, no_b2*2))
        allocate(wrk_nto(2*no_b2), source = 0.0_dp)
        do i = 1, nwf_1
            if (print_level >= 2) write(stdout, '(9x,a,i0,a,i0)') 'N-1 el. state ', i, '/', nwf_1
            ! CIS - Ref
            call gemm(mo_a1(:, :, i), s_mo_a, wrk0_a, transa='T')
            call gemm(mo_b1(:, :, i), s_mo_b, wrk0_b, transa='T')
            do j = 1, nwf_2
                if (print_level >= 2) write(stdout, '(13x,a,i0,a,i0)') 'N el. state ', j, '/', nwf_2
                call gemm(wrk0_a, mo_a2(:, :, j), wrk_a)
                call gemm(wrk0_b, mo_b2(:, :, j), wrk_b)
                ! Ref - CIS
                if (i == 1) then
                    call nto_rs(no_a, na_a2(j), wrk_a, c_a2(:, j), rs_a(j), .false.)
                    call nto_rs_dys(no_b2, na_b2(j), wrk_b, c_b2(:, j), wrk_nto)
                    call gemv(mo_b2(:, :, j), wrk_nto, rs_b(:, j))
                end if
                ! CIS - Ref
                if (j == 1) then
                    call nto_rs(no_a, na_a1(i), wrk_a, c_a1(:, i), sr_a(i), .true.)
                    call nto_sr_dys(no_b2, na_b1(i), wrk_b, c_b1(:, i), wrk_nto(1:no_b2))
                    call gemv(mo_b2(:, 1:no_b2, i), wrk_nto(1:no_b2), sr_b(1:no_b2, i))
                end if
                ! CIS - CIS
                call nto_ss(no_a, na_a1(i), na_a2(j), wrk_a, c_a1(:, i), c_a2(:, j), ss_a(i, j))
                call nto_ss_dys(no_b2, na_b1(i), na_b2(j), wrk_b, c_b1(:, i), c_b2(:, j), wrk_nto)
                call gemv(mo_b2(:, :, j), wrk_nto, ss_b(:, i, j))
            end do
        end do
        deallocate(wrk_a)
        deallocate(wrk_b)
        deallocate(wrk0_a)
        deallocate(wrk0_b)
        time_det = omp_get_wtime() - time0

        allocate(dys_mo(n_b2, 0:nwf_1, 0:nwf_2))
        dys_mo(:, 0, 0) = rr_b * rr_a
        do i = 1, nwf_1
            dys_mo(:, i, 0) = rr_b * sr_a(i) + sr_b(:, i) * rr_a
        end do
        do j = 1, nwf_2
            dys_mo(:, 0, j) = rr_b * rs_a(j) + rs_b(:, j) * rr_a
        end do
        do i = 1, nwf_1
            do j = 1, nwf_2
                dys_mo(:, i, j) = rr_b * ss_a(i, j) + rs_b(:, j) * sr_a(i) + &
                                  ss_b(:, i, j) * rr_a + sr_b(:, i) * rs_a(j)
            end do
        end do

        if (print_level >= 2) then
            write(stdout, *)
            write(stdout,'(5x, a)') 'cis_dyson time:'
            time_tot = omp_get_wtime() - time00
            write(stdout, '(9x, a40, f14.4)') 'NTO generation            - time (sec):', time_nto
            write(stdout, '(9x, a40, f14.4)') 'Determinant blocks        - time (sec):', time_det
            write(stdout, '(9x, 40x, a14)') '--------------'
            write(stdout, '(9x, a40, f14.4)') 'Total                     - time (sec):', time_tot
        end if
        if (print_level >= 2) then
            write(stdout, *)
            write(stdout, '(5x,a)') '---- end cis_dyson subroutine ----'
        end if
    end subroutine cis_dyson


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: nto_rr_dys
    !> @brief Calculate the RR block of determinants for a Dyson orbital calculation.
    !----------------------------------------------------------------------------------------------
    subroutine nto_rr_dys(nn, s_nto, rr)
        use matrix_mod, only : mat_ge_det
        integer, intent(in) :: nn !< Number of occupied orbitals in neutral.
        real(dp), intent(in) :: s_nto(:, :) !< Overlap matrix between NTOs.
        real(dp), intent(out) :: rr(nn) !< RR block coefficients.
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
    !> @brief Calculate the SR block of determinants for a Dyson orbital calculation.
    !----------------------------------------------------------------------------------------------
    subroutine nto_sr_dys(nn, na, s_nto, c, sr)
        use matrix_mod, only : mat_ge_det
        integer, intent(in) :: nn !< Number of occupied orbitals in neutral.
        integer, intent(in) :: na !< Number of active excitations.
        real(dp), intent(in) :: s_nto(:, :) !< Overlap matrix between NTOs.
        real(dp), intent(in) :: c(nn-1) !< Coefficients of excitations.
        real(dp), intent(out) :: sr(nn) !< RS block coefficients.
        real(dp) :: wrk(nn-1, nn-1) !< Work array.
        integer :: seq(nn) !< All columns.
        integer :: cols(nn-1) !< Currently active columns.
        integer :: sgn !< Sign of the current minor.
        integer :: nc !< Number of occupied orbitals in cation.
        integer :: i, j

        if (na == 0) return
        nc = nn - 1
        seq = [ (i, i=1, nn) ]
        sr = 0.0_dp

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
    !> @brief Calculate the RS block of determinants for a Dyson orbital calculation.
    !----------------------------------------------------------------------------------------------
    subroutine nto_rs_dys(nn, na, s_nto, c, rs)
        use matrix_mod, only : mat_ge_det
        integer, intent(in) :: nn !< Number of occupied orbitals in neutral.
        integer, intent(in) :: na !< Number of active excitations.
        real(dp), intent(in) :: s_nto(:, :) !< Overlap matrix between NTOs.
        real(dp), intent(in) :: c(nn) !< Coefficients of excitations.
        real(dp), intent(out) :: rs(2*nn) !< RS block coefficients.
        real(dp) :: wrk(nn-1, nn-1) !< Work array.
        integer :: seq(nn) !< All columns.
        integer :: cols(nn-1) !< Currently active columns.
        integer :: sgn !< Sign of the current minor.
        integer :: nc !< Number of occupied orbitals in cation.
        integer :: i, j

        if (na == 0) return
        nc = nn - 1
        seq = [ (i, i=1, nn) ]
        rs = 0.0_dp

        !$omp parallel default(shared)
        !$omp do private(wrk, j, cols, sgn) schedule(dynamic)
        do i = 1, nn
            sgn = (-1)**(i+1)
            cols = pack(seq, (seq /= i))
            wrk = s_nto(1:nc, cols)
            rs(nn+i) = rs(nn+i) + sgn * c(i) * mat_ge_det(wrk)
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
    !> @brief Calculate the SS block of determinants for a Dyson orbital calculation.
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
        seq = [ (i, i=1, nn) ]
        ss = 0.0_dp

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
                    else if (k > i) then
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
