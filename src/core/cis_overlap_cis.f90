!----------------------------------------------------------------------------------------------
! MODULE: cis_overlap_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date November, 2018
!
! DESCRIPTION:
!> @brief Subroutines for calculating overlaps between CIS wave functions.
!----------------------------------------------------------------------------------------------
module cis_overlap_cis_mod
    use global_defs
    implicit none


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_overlap_cis
    !
    !> @brief Calculate overlap matrix between two sets of CIS wave functions.
    !> @details
    !! If a threshold trunc between 0 and 1 is given, the wave function expansions are truncated 
    !! to include only the dominant contributions until the square of the norm for each state is
    !! higher than the threshold.
    !
    !> @note The overlaps between the bra(ket) reference and ket(bra) CIS states are returned as
    !! s_wf(0, :) and s_wf(:, 0), respectively.
    !----------------------------------------------------------------------------------------------
    subroutine cis_overlap_cis(rhf, trunc, s_mo_a, s_mo_b, cis_a1, cis_a2, cis_b1, cis_b2, s_wf)
        use blas95, only : gemm
        use matrix_mod, only : mat_ge_det
        use truncate_wf_mod
        integer, intent(in) :: rhf !< Restricted (1) or unrestricted (2) calculation.
        real(dp), intent(in) :: trunc !< Threshold for truncating the wave functions.
        real(dp), intent(in) :: s_mo_a(:, :) !< Overlaps of alpha molecular orbitals.
        real(dp), intent(in) :: s_mo_b(:, :) !< Overlaps of beta molecular orbitals.
        real(dp), intent(in) :: cis_a1(:, :, :) !< CIS matrix alpha 1.
        real(dp), intent(in) :: cis_a2(:, :, :) !< CIS matrix alpha 2.
        real(dp), intent(in) :: cis_b1(:, :, :) !< CIS matrix beta 1.
        real(dp), intent(in) :: cis_b2(:, :, :) !< CIS matrix beta 2.
        real(dp), allocatable, intent(out) :: s_wf(:, :) !< Wave function overlaps.
        logical :: beta !< Restricted/unrestricted calculation.
        integer :: n_1 !< Number of bra orbitals.
        integer :: n_2 !< Number of ket orbitals.
        integer :: no_a !< Number of occupied alpha orbitals.
        integer :: no_b !< Number of occupied beta orbitals.
        integer :: nv_a1 !< Number of virtual alpha orbitals.
        integer :: nv_b1 !< Number of virtual beta orbitals.
        integer :: nv_a2 !< Number of virtual alpha orbitals.
        integer :: nv_b2 !< Number of virtual beta orbitals.
        integer :: nwf_1 !< Number of bra states.
        integer :: nwf_2 !< Number of ket states.
        real(dp) :: rr_a !< RR block alpha.
        real(dp) :: rr_b !< RR block beta.
        real(dp), allocatable :: wf_a1(:, :) !< Bra alpha wave function coefficients.
        real(dp), allocatable :: wf_a2(:, :) !< Ket alpha wave function coefficients.
        real(dp), allocatable :: wf_b1(:, :) !< Bra beta wave function coefficients.
        real(dp), allocatable :: wf_b2(:, :) !< Ket beta wave function coefficients.
        real(dp), allocatable :: sr_a(:) !< SR block alpha.
        real(dp), allocatable :: sr_b(:) !< SR block beta.
        real(dp), allocatable :: rs_a(:) !< RS block alpha.
        real(dp), allocatable :: rs_b(:) !< RS block beta.
        real(dp), allocatable :: ss_a(:, :) !< SS block alpha.
        real(dp), allocatable :: ss_b(:, :) !< SS block beta.
        logical, allocatable :: a_a1(:, :) !< Bra active alpha excitations.
        logical, allocatable :: a_a2(:, :) !< Ket active alpha excitations.
        logical, allocatable :: a_b1(:, :) !< Bra active beta excitations.
        logical, allocatable :: a_b2(:, :) !< Ket active beta excitations.
        integer :: i
        real(dp), external :: omp_get_wtime
        real(dp) :: time00, time0
        real(dp) :: time_rs_sr, time_ss, time_tot

        time00 = omp_get_wtime()
        if (print_level >= 2) then
            write(stdout, *)
            write(stdout, '(5x,a)') '---- start cis_overlap_nto subroutine ----'
        end if

        ! Get dimensions.
        beta = .false.
        if (rhf == 2) beta = .true.
        no_a = size(cis_a1, 2)
        no_b = size(cis_b1, 2)
        n_1 = size(s_mo_a, 1)
        n_2 = size(s_mo_a, 2)
        nv_a1 = n_1 - no_a
        nv_a2 = n_2 - no_a
        if (beta) then
            nv_b1 = n_1 - no_b
            nv_b2 = n_2 - no_b
        end if
        nwf_1 = size(cis_a1, 3)
        nwf_2 = size(cis_a2, 3)
        ! Allocate work arrays.
        allocate(wf_a1(nv_a1*no_a, nwf_1), source=reshape(cis_a1, [nv_a1*no_a, nwf_1]))
        allocate(wf_a2(nv_a2*no_a, nwf_2), source=reshape(cis_a2, [nv_a2*no_a, nwf_2]))
        allocate(sr_a(nwf_1))
        allocate(rs_a(nwf_2))
        allocate(ss_a(nwf_1, nwf_2))
        if (beta) then
            allocate(wf_b1(nv_b1*no_b, nwf_1), source=reshape(cis_b1, [nv_b1*no_b, nwf_1]))
            allocate(wf_b2(nv_b2*no_b, nwf_2), source=reshape(cis_b2, [nv_b2*no_b, nwf_2]))
            allocate(sr_b(nwf_1))
            allocate(rs_b(nwf_2))
            allocate(ss_b(nwf_1, nwf_2))
        end if

        ! Truncate wave functions:
        call truncate_wf(trunc, beta, wf_a1, wf_b1, a_a1, a_b1)
        call truncate_wf(trunc, beta, wf_a2, wf_b2, a_a2, a_b2)
        if (trunc < 1.0_dp) then
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(5x,a,f0.8)') 'Truncating wave functions based on threshold ', trunc
                write(stdout, '(5x,a)') 'Number of remaining determinants for bra states:'
                write(stdout, '(9x,1000(i0,1x))') count(a_a1, 1)
                if (beta) write(stdout, '(9x,1000(i0,1x))') count(a_b1, 1)
                write(stdout, '(5x,a)') 'Number of remaining determinants for ket states:'
                write(stdout, '(9x,1000(i0,1x))') count(a_a2, 1)
                if (beta) write(stdout, '(9x,1000(i0,1x))') count(a_b2, 1)
            end if
        end if


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
        call cis_rs(no_a, nv_a1, s_mo_a, wf_a1, a_a1, sr_a, row = .true.)
        if (print_level >= 2) write(stdout, '(9x,a)') 'RS block done.'
        call cis_rs(no_a, nv_a2, s_mo_a, wf_a2, a_a2, rs_a, row = .false.)
        if (print_level >= 2) write(stdout, '(9x,a)') 'SR block done.'
        time_rs_sr = omp_get_wtime() - time0

        time0 = omp_get_wtime()
        call cis_ss(no_a, nv_a1, nv_a2, s_mo_a, wf_a1, wf_a2, a_a1, a_a2, ss_a)
        if (print_level >= 2) write(stdout, '(9x,a)') 'SS block done.'
        time_ss = omp_get_wtime() - time0
        
        if (beta) then
            time0 = omp_get_wtime()
            if (print_level >= 2) then
                write(stdout, '(5x,a)') 'Computing beta determinant blocks...'
                write(stdout, '(5x,a)') 'Status:'
            end if
            rr_b = mat_ge_det(s_mo_b(1:no_b, 1:no_b))
            call cis_rs(no_b, nv_b1, s_mo_b, wf_b1, a_b1, sr_b, row = .true.)
            if (print_level >= 2) write(stdout, '(9x,a)') 'RS block done.'
            call cis_rs(no_b, nv_b2, s_mo_b, wf_b2, a_b2, rs_b, row = .false.)
            if (print_level >= 2) write(stdout, '(9x,a)') 'SR block done.'
            time_rs_sr = time_rs_sr + omp_get_wtime() - time0
            
            time0 = omp_get_wtime()
            call cis_ss(no_b, nv_b1, nv_b2, s_mo_b, wf_b1, wf_b2, a_b1, a_b2, ss_b)
            if (print_level >= 2) write(stdout, '(9x,a)') 'SS block done.'
            time_ss = time_ss + omp_get_wtime() - time0
        end if

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
            write(stdout,'(5x, a)') 'cis_overlap_cis time:'
            time_tot = omp_get_wtime() - time00
            write(stdout, '(9x, a40, f14.4)') 'RS and SR blocks          - time (sec):', time_rs_sr
            write(stdout, '(9x, a40, f14.4)') 'SS block                  - time (sec):', time_ss
            write(stdout, '(9x, 40x, a14)') '--------------'
            write(stdout, '(9x, a40, f14.4)') 'Total                     - time (sec):', time_tot
        end if
        if (print_level >= 2) then
            write(stdout, *)
            write(stdout, '(5x,a)') '---- end cis_overlap_cis subroutine ----'
        end if
    end subroutine cis_overlap_cis


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_rs
    !> @brief Calculate the RS or SR block between two states in a cis overlap calculation.
    !> @note Row is .false. for RS block calculations and .true. for SR block calculations.
    !----------------------------------------------------------------------------------------------
    subroutine cis_rs(no, nv, s_mo, c, m, rs, row)
        use matrix_mod, only : mat_ge_det
        integer, intent(in) :: no !< Number of occupied orbitals.
        integer, intent(in) :: nv !< Number of bra virtual orbitals.
        real(dp), intent(in) :: s_mo(:, :) !< Overlap matrix between NTOs.
        real(dp), intent(in) :: c(:, :) !< Excitation coefficients.
        logical, intent(in) :: m(:, :) !< Active excitations mask.
        real(dp), intent(out) :: rs(:) !< RS/SR block sum.
        logical, intent(in) :: row !< RS or SR block.
        real(dp) :: wrk(no, no) !< Work array.
        integer :: i, a, nwf, o, v
        real(dp) :: cdet

        nwf = size(c, 2)
        rs = 0.0_dp

        !$omp parallel default(shared)
        !$omp do private(wrk, a, o, v, cdet) schedule(dynamic) reduction(+:rs)
        do i = 1, no*nv
            if (.not. any(m(i, :))) cycle
            o = int((i-1)/nv) + 1
            v = mod((i-1),nv) + 1
            wrk = s_mo(1:no, 1:no)
            if (row) then
                wrk(o, :) = s_mo(no+v, 1:no)
            else
                wrk(:, o) = s_mo(1:no, no+v)
            end if
            cdet = mat_ge_det(wrk)
            do a = 1, nwf
                if (.not. m(i, a)) cycle
                rs(a) = rs(a) + cdet * c(i, a)
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine cis_rs


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_ss
    !> @brief Calculate the SS block between two states in a cis overlap calculation.
    !----------------------------------------------------------------------------------------------
    subroutine cis_ss(no, nv1, nv2, s_mo, c1, c2, m1, m2, ss)
        use matrix_mod, only : mat_ge_det
        integer, intent(in) :: no !< Number of occupied orbitals.
        integer, intent(in) :: nv1 !< Number of bra virtual orbitals.
        integer, intent(in) :: nv2 !< Number of ket virtual orbitals.
        real(dp), intent(in) :: s_mo(:, :) !< Overlap matrix between NTOs.
        real(dp), intent(in) :: c1(:, :) !< Bra excitation coefficients.
        real(dp), intent(in) :: c2(:, :) !< Ket excitation coefficients.
        logical, intent(in) :: m1(:, :) !< Bra active excitations mask.
        logical, intent(in) :: m2(:, :) !< Ket active excitations mask.
        real(dp), intent(out) :: ss(:, :) !< SS block sum.
        real(dp) :: wrk(no, no) !< Work array.
        integer :: i, j, a, b, nwf1, nwf2, o1, v1, o2, v2
        real(dp) :: cdet

        nwf1 = size(c1, 2)
        nwf2 = size(c2, 2)
        ss = 0.0_dp

        !$omp parallel default(shared)
        !$omp do private(wrk, j, a, b, o1, v1, o2, v2, cdet) schedule(dynamic) reduction(+:ss)
        do i = 1, no*nv1
            if (.not. any(m1(i, :))) cycle
            o1 = int((i-1)/nv1) + 1
            v1 = mod((i-1),nv1) + 1
            do j = 1, no*nv2
                if (.not. any(m2(j, :))) cycle
                o2 = int((j-1)/nv2) + 1
                v2 = mod((j-1),nv2) + 1
                wrk = s_mo(1:no, 1:no)
                wrk(o1, :) = s_mo(no+v1, 1:no)
                wrk(:, o2) = s_mo(1:no, no+v2)
                wrk(o1, o2) = s_mo(no+v1, no+v2)
                cdet = mat_ge_det(wrk)
                do a = 1, nwf1
                    if (.not. m1(i, a)) cycle
                    do b = 1, nwf2
                        if (.not. m2(j, b)) cycle
                        ss(a, b) = ss(a, b) + cdet * c1(i, a) * c2(j, b)
                    end do
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine cis_ss


end module cis_overlap_cis_mod
