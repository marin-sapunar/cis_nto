!----------------------------------------------------------------------------------------------
! MODULE: cis_overlap_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date October, 2018
!
!> @brief Subroutines for calculating overlaps between CIS wave functions.
!----------------------------------------------------------------------------------------------
module cis_overlap_nto_mod
    use global_defs
    use cis_nto_mod
    implicit none


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_overlap_nto
    !
    !> @brief Calculate overlap matrix between two sets of CIS wave functions.
    !> @details
    !! Transform the wave functions into expansion in terms of excitations between NTOs before
    !! calculating the overlaps. If a threshold trunc between 0 and 1 is given, the wave function 
    !! expansions are truncated to include only the dominant contributions until the square of the
    !! norm for each state is higher than the threshold.
    !
    !> @note The overlaps between the bra(ket) reference and ket(bra) CIS states are returned as
    !! s_wf(0, :) and s_wf(:, 0), respectively.
    !----------------------------------------------------------------------------------------------
    subroutine cis_overlap_nto(rhf, trunc_norm, trunc_nex, s_mo_a, s_mo_b, wf_a1, wf_a2, wf_b1,    &
    &                          wf_b2, s_wf)
        use linalg_wrapper_mod, only : gemm
        use matrix_mod, only : mat_ge_det, vec_outer
        integer, intent(in) :: rhf !< Restricted (1) or unrestricted (2) calculation
        real(dp), intent(in) :: trunc_norm !< Threshold for norm based truncation of wave functions
        integer, intent(in) :: trunc_nex !< Truncation to fixed number of dominant excitations
        real(dp), intent(in) :: s_mo_a(:, :) !< Overlaps of alpha molecular orbitals
        real(dp), intent(in) :: s_mo_b(:, :) !< Overlaps of beta molecular orbitals
        real(dp), intent(in) :: wf_a1(:, :, :) !< Bra alpha wave function coefficients
        real(dp), intent(in) :: wf_a2(:, :, :) !< Ket alpha wave function coefficients
        real(dp), intent(in) :: wf_b1(:, :, :) !< Bra beta wave function coefficients
        real(dp), intent(in) :: wf_b2(:, :, :) !< Ket beta wave function coefficients
        real(dp), allocatable, intent(out) :: s_wf(:, :) !< Wave function overlaps
        logical :: beta !< Restricted/unrestricted calculation
        integer :: n_a2 !< Number of ket alpha orbitals
        integer :: n_b2 !< Number of ket beta orbitals
        integer :: no_a !< Number of occupied alpha orbitals
        integer :: no_b !< Number of occupied beta orbitals
        integer :: nwf_1 !< Number of bra states
        integer :: nwf_2 !< Number of ket states
        real(dp), allocatable :: c_a1(:, :) !< Bra NTO coefficients alpha
        real(dp), allocatable :: c_a2(:, :) !< Ket NTO coefficients alpha
        real(dp), allocatable :: c_b1(:, :) !< Bra NTO coefficients beta
        real(dp), allocatable :: c_b2(:, :) !< Ket NTO coefficients beta
        real(dp), allocatable :: mo_a1(:, :, :) !< Bra NTOs in MO basis alpha
        real(dp), allocatable :: mo_a2(:, :, :) !< Ket NTOs in MO basis alpha
        real(dp), allocatable :: mo_b1(:, :, :) !< Bra NTOs in MO basis beta
        real(dp), allocatable :: mo_b2(:, :, :) !< Ket NTOs in MO basis beta
        real(dp), allocatable :: s_nto_a(:, :) !< NTO overlaps alpha
        real(dp), allocatable :: s_nto_b(:, :) !< NTO overlaps beta
        real(dp), allocatable :: wrk(:, :) !< Work array
        real(dp) :: rr_a !< RR block alpha
        real(dp) :: rr_b !< RR block beta
        real(dp), allocatable :: sr_a(:) !< SR block alpha
        real(dp), allocatable :: sr_b(:) !< SR block beta
        real(dp), allocatable :: rs_a(:) !< RS block alpha
        real(dp), allocatable :: rs_b(:) !< RS block beta
        real(dp), allocatable :: ss_a(:, :) !< SS block alpha
        real(dp), allocatable :: ss_b(:, :) !< SS block beta
        integer, allocatable :: na_a1(:) !< Number of active bra orbitals alpha
        integer, allocatable :: na_a2(:) !< Number of active ket orbitals alpha
        integer, allocatable :: na_b1(:) !< Number of active bra orbitals beta
        integer, allocatable :: na_b2(:) !< Number of active ket orbitals beta
        integer :: i, j
        real(dp), external :: omp_get_wtime
        real(dp) :: time00, time0
        real(dp) :: time_nto, time_det, time_tot

        time00 = omp_get_wtime()
        if (print_level >= 2) then
            write(stdout, *) 
            write(stdout, '(5x,a)') '---- start cis_overlap_nto subroutine ----'
        end if

        ! Get dimensions:
        beta = .false.
        if (rhf == 2) beta = .true.
        n_a2 = size(s_mo_a, 2)
        n_b2 = size(s_mo_b, 2)
        no_a = size(wf_a1, 2)
        if (beta) no_b = size(wf_b1, 2)
        nwf_1 = size(wf_a1, 3)
        nwf_2 = size(wf_a2, 3)
        ! Allocate work arrays:
        allocate(s_nto_a(no_a*2, no_a*2))
        allocate(sr_a(nwf_1))
        allocate(rs_a(nwf_2))
        allocate(ss_a(nwf_1, nwf_2))
        if (beta) then
            allocate(s_nto_b(no_b*2, no_b*2))
            allocate(sr_b(nwf_1))
            allocate(rs_b(nwf_2))
            allocate(ss_b(nwf_1, nwf_2))
        end if

        ! Calculate NTOs:
        time0 = omp_get_wtime()
        if (print_level >= 2) then
            write(stdout, *) 
            write(stdout, '(5x,a)') 'Generating NTOs...'
        end if
        call cis_nto(wf_a1, c_a1, mo_a1)
        call cis_nto(wf_a2, c_a2, mo_a2)
        if (beta) then
            call cis_nto(wf_b1, c_b1, mo_b1)
            call cis_nto(wf_b2, c_b2, mo_b2)
        end if
        time_nto = omp_get_wtime() - time0

        ! Truncate wave functions:
        call cis_nto_truncate(beta, trunc_norm, trunc_nex, c_a1, c_b1, na_a1, na_b1)
        call cis_nto_truncate(beta, trunc_norm, trunc_nex, c_a2, c_b2, na_a2, na_b2)
        if (print_level >= 1) then
            if (trunc_norm < 1.0_dp) then
                write(stdout, *)
                write(stdout, '(5x,a,f0.8)') 'Truncating wfs based on norm threshold: ', trunc_norm
                write(stdout, '(5x,a)') 'Number of remaining determinants for bra states:'
                write(stdout, '(9x,1000(i0,1x))') na_a1
                if (beta) write(stdout, '(9x,1000(i0,1x))') na_b1
                write(stdout, '(5x,a)') 'Number of remaining determinants for ket states:'
                write(stdout, '(9x,1000(i0,1x))') na_a2
                if (beta) write(stdout, '(9x,1000(i0,1x))') na_b2
            end if
            if (trunc_nex > 0) then
                write(stdout, *)
                write(stdout, '(5x,a,i0,a)') 'Using ', trunc_nex, ' dominant excitations per state.'
                write(stdout, '(5x,a)') 'Norms of states after truncation:'
                write(stdout,'(5x,8es11.3)') cis_norms_nto(beta, c_a1, c_b1, na_a1, na_b1)
                write(stdout, '(5x,a)') 'Norms of states after truncation:'
                write(stdout,'(5x,8es11.3)') cis_norms_nto(beta, c_a2, c_b2, na_a2, na_b2)
            end if
            write(stdout, *)
            if (beta) then
                i = sum(vec_outer(na_a1, na_a2))
                write(stdout, '(5x,a,i0)') 'Total number of alpha overlap determinants: ', i
                i = sum(vec_outer(na_b1, na_b2))
                write(stdout, '(5x,a,i0)') 'Total number of beta overlap determinants: ', i
            else
                i = sum(vec_outer(na_a1, na_a2))
                write(stdout, '(5x,a,i0)') 'Total number of overlap determinants: ', i
            end if
        end if

        ! Calculate determinant blocks:
        time0 = omp_get_wtime()
        if (print_level >= 2) then
            write(stdout, *)
            if (beta) then
                write(stdout, '(5x,a)') 'Computing alpha determinant blocks...'
            else
                write(stdout, '(5x,a)') 'Computing determinant blocks...'
            end if
            write(stdout, '(5x,a)') 'Status:'
        end if
        rr_a = mat_ge_det(s_mo_a(1:no_a, 1:no_a))
        allocate(wrk(no_a*2, n_a2))
        do i = 1, nwf_1
            if (print_level >= 2) write(stdout, '(9x,a,i0,a,i0)') 'bra state ', i, '/', nwf_1
            call gemm(mo_a1(:, :, i), s_mo_a, wrk, transa='T')
            do j = 1, nwf_2
                if (print_level >= 2) write(stdout, '(13x,a,i0,a,i0)') 'ket state ', j, '/', nwf_2
                call gemm(wrk, mo_a2(:, :, j), s_nto_a)
                if (j == 1) call nto_rs(no_a, na_a1(i), s_nto_a, c_a1(:, i), sr_a(i), .true.)
                if (i == 1) call nto_rs(no_a, na_a2(j), s_nto_a, c_a2(:, j), rs_a(j), .false.)
                call nto_ss(no_a, na_a1(i), na_a2(j), s_nto_a, c_a1(:, i), c_a2(:, j), ss_a(i, j))
            end do
        end do
        deallocate(wrk)

        if (beta) then
            if (print_level >= 2) then
                write(stdout, '(5x,a)') 'Computing beta determinant blocks...'
                write(stdout, '(5x,a)') 'Status:'
            end if
            rr_b = mat_ge_det(s_mo_b(1:no_b, 1:no_b))
            allocate(wrk(no_b*2, n_b2))
            do i = 1, nwf_1
                if (print_level >= 2) write(stdout, '(9x,a,i0,a,i0)') 'bra state ', i, '/', nwf_1
                call gemm(mo_b1(:, :, i), s_mo_b, wrk, transa='T')
                do j = 1, nwf_2
                    if (print_level >= 2) write(stdout, '(13x,a,i0,a,i0)') 'ket state ', j, '/', nwf_2
                    call gemm(wrk, mo_b2(:, :, j), s_nto_b)
                    if (j == 1) call nto_rs(no_b, na_b1(i), s_nto_b, c_b1(:, i), sr_b(i), .true.)
                    if (i == 1) call nto_rs(no_b, na_b2(j), s_nto_b, c_b2(:, j), rs_b(j), .false.)
                    call nto_ss(no_b, na_b1(i), na_b2(j), s_nto_b, c_b1(:, i), c_b2(:, j), ss_b(i, j))
                end do
            end do
            deallocate(wrk)
        end if
        time_det = omp_get_wtime() - time0

        if (allocated(s_wf)) deallocate(s_wf)
        allocate(s_wf(0:nwf_1, 0:nwf_2))
        if (beta) then
            ! Ground state overlap:
            s_wf(0, 0) = rr_a*rr_b
            ! Ground-excited state overlap:
            s_wf(1:, 0) = rr_a*sr_b + rr_b*sr_a
            ! Excited-ground state overlap:
            s_wf(0, 1:) = rr_b*rs_a + rr_a*rs_b
            ! Excited-excited state overlap:
            do i = 1, nwf_1
                s_wf(i, 1:) = rr_b*ss_a(i, :) + sr_a(i)*rs_b + rr_a*ss_b(i, :) + sr_b(i)*rs_a
            end do
        else
            ! Ground state overlap:
            s_wf(0, 0) = rr_a*rr_a
            ! Ground-excited state overlap:
            s_wf(1:, 0) = rr_a*sr_a
            ! Excited-ground state overlap:
            s_wf(0, 1:) = rr_a*rs_a
            ! Excited-excited state overlap:
            do i = 1, nwf_1
                s_wf(i, 1:) = rr_a*ss_a(i, :) + sr_a(i)*rs_a
            end do
        end if

        if (print_level >= 2) then
            write(stdout, *) 
            write(stdout,'(5x, a)') 'cis_overlap_nto time:'
            time_tot = omp_get_wtime() - time00
            write(stdout, '(9x, a40, f14.4)') 'NTO generation            - time (sec):', time_nto
            write(stdout, '(9x, a40, f14.4)') 'Determinant blocks        - time (sec):', time_det
            write(stdout, '(9x, 40x, a14)') '--------------'
            write(stdout, '(9x, a40, f14.4)') 'Total                     - time (sec):', time_tot
        end if
        if (print_level >= 2) then
            write(stdout, *) 
            write(stdout, '(5x,a)') '---- end cis_overlap_nto subroutine ----'
        end if
    end subroutine cis_overlap_nto


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: nto_rs
    !> @brief Calculate the RS or SR block between two states in a cis overlap calculation.
    !> @note Row is .false. for RS block calculations and .true. for SR block calculations.
    !----------------------------------------------------------------------------------------------
    subroutine nto_rs(n, na, s_nto, c, rs, row)
        use matrix_mod, only : mat_ge_det
        use linalg_wrapper_mod, only : dot
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
        use linalg_wrapper_mod, only : gemv, dot
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


end module cis_overlap_nto_mod
