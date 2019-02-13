module cis_overlap_l2m_mod
    use global_defs
    use matrix_mod, only : mat_ge_det

    implicit none
    
    private
    public :: cis_overlap_l2m

contains


    subroutine cis_overlap_l2m(rhf, s_mo_a, s_mo_b, cis_a1, cis_a2, cis_b1, cis_b2, s_wf)
     !  real(dp), intent(in) :: csc(:, :, :) ! MO overlap matrix.
     !  integer, intent(in) :: no(:) !< Number of occupied orbitals.
     !  integer, intent(in) :: nv1(:) !< Number of virtual orbitals (bra).
     !  integer, intent(in) :: nv2(:) !< Number of virtual orbitals (ket).
     !  type(rmat), intent(inout) :: wf1(:) !< Ket wave function coefficients.
     !  type(rmat), intent(inout) :: wf2(:) !< Bra wave function coefficients.
     !  real(dp), allocatable, intent(out) :: omat(:, :) !< Overlap matrix.

        integer, intent(in) :: rhf !< Restricted (1) or unrestricted (2) calculation.
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
        real(dp), allocatable :: wf_a1(:, :) !< Bra alpha wave function coefficients.
        real(dp), allocatable :: wf_a2(:, :) !< Ket alpha wave function coefficients.
        real(dp), allocatable :: wf_b1(:, :) !< Bra beta wave function coefficients.
        real(dp), allocatable :: wf_b2(:, :) !< Ket beta wave function coefficients.

        real(dp) :: rr_a !< RR determinant alpha.
        real(dp) :: rr_b !< RR determinant beta.
        real(dp), allocatable :: sr_a(:) !< SR block alpha.
        real(dp), allocatable :: sr_b(:) !< SR block beta.
        real(dp), allocatable :: rs_a(:) !< RS block alpha.
        real(dp), allocatable :: rs_b(:) !< RS block beta.
        real(dp), allocatable :: ss_a(:, :) !< SS block alpha.
        real(dp), allocatable :: ss_b(:, :) !< SS block beta.

        real(dp), allocatable :: ref(:, :)
        real(dp), allocatable :: l1cminor(:, :)
        real(dp), allocatable :: l1rminor(:, :)
        ! Timing variables
        real(dp) :: start, finish
        real(dp) :: starttot, finishtot, total
        real(dp) :: wfbt, rrat, l1minort, l2minort, rsblkt, srblkt, ssblkt, omatt 

        ! External routines
        external :: ssblock
        integer, external :: mkl_get_max_threads
        integer :: num_threads

        integer :: i
        real(dp), external :: omp_get_wtime
        real(dp) :: time00, time0
        real(dp) :: time_l1m, time_rs_sr, time_ss, time_tot

        time00 = omp_get_wtime()
        if (print_level >= 2) then
            write(stdout, *)
            write(stdout, '(5x,a)') '---- start cis_overlap_l2m subroutine ----'
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

        time0 = omp_get_wtime()
       !num_threads = mkl_get_max_threads()
       !call mkl_set_num_threads(1)
        if (print_level >= 2) then
            write(stdout, *)
            if (beta) then
                write(stdout, '(5x,a)') 'Computing alpha level 1 minors...'
            else
                write(stdout, '(5x,a)') 'Computing level 1 minors...'
            end if
        end if
        allocate(l1cminor(no_a, no_a))
        allocate(l1rminor(no_a, no_a))
        allocate(ref(no_a, no_a))
        ref = s_mo_a(1:no_a, 1:no_a)
        call getlvl1minors(no_a, ref, l1rminor, 'r')
        call getlvl1minors(no_a, ref, l1cminor, 'c')
        time_l1m = omp_get_wtime() - time0

        if (print_level >= 2) then
            write(stdout, *)
            if (beta) then
                write(stdout, '(5x,a)') 'Computing alpha determinant blocks...'
            else
                write(stdout, '(5x,a)') 'Computing determinant blocks...'
            end if
        end if
        time0 = omp_get_wtime()
        rr_a = mat_ge_det(ref)
        call rsblock(s_mo_a, no_a, nv_a2, nwf_2, wf_a2, l1cminor, rs_a)
        if (print_level >= 2) write(stdout, '(9x,a)') 'RS block done.'
        call srblock(s_mo_a, no_a, nv_a1, nwf_1, wf_a1, l1rminor, sr_a)
        if (print_level >= 2) write(stdout, '(9x,a)') 'SR block done.'
        time_rs_sr = omp_get_wtime() - time0
       !call mkl_set_num_threads(num_threads)
        time0 = omp_get_wtime()
        call ssblock(s_mo_a, no_a, nv_a1, nv_a2, nwf_1, nwf_2, wf_a1, wf_a2, l1rminor, ss_a)
        time_ss = omp_get_wtime() - time0
        deallocate(l1cminor)
        deallocate(l1rminor)
        deallocate(ref)

        if (beta) then
            time0 = omp_get_wtime()
           !num_threads = mkl_get_max_threads()
           !call mkl_set_num_threads(1)
            if (print_level >= 2) then
                write(stdout, *)
                write(stdout, '(5x,a)') 'Computing beta level 1 minors...'
            end if
            allocate(l1cminor(no_b, no_a))
            allocate(l1rminor(no_b, no_b))
            allocate(ref(no_b, no_b))
            ref = s_mo_b(1:no_b, 1:no_b)
            call getlvl1minors(no_b, ref, l1rminor, 'r')
            call getlvl1minors(no_b, ref, l1cminor, 'c')
            time_l1m = time_l1m + omp_get_wtime() - time0
         
            if (print_level >= 2) then
                write(stdout, *)
                write(stdout, '(5x,a)') 'Computing beta determinant blocks...'
            end if
            time0 = omp_get_wtime()
            rr_b = mat_ge_det(ref)
            if (print_level >= 2) write(stdout, '(9x,a)') 'RS block done.'
            call rsblock(s_mo_b, no_b, nv_b2, nwf_2, wf_b2, l1cminor, rs_b)
            if (print_level >= 2) write(stdout, '(9x,a)') 'SR block done.'
            call srblock(s_mo_b, no_b, nv_b1, nwf_1, wf_b1, l1rminor, sr_b)
            time_rs_sr = time_rs_sr + omp_get_wtime() - time0
           !call mkl_set_num_threads(num_threads)
            time0 = omp_get_wtime()
            call ssblock(s_mo_b, no_b, nv_b1, nv_b2, nwf_1, nwf_2, wf_b1, wf_b2, l1rminor, ss_b)
            time_ss = time_ss + omp_get_wtime() - time0
            deallocate(l1cminor)
            deallocate(l1rminor)
            deallocate(ref)
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
            write(stdout,'(5x, a)') 'cis_overlap_l2m time:'
            time_tot = omp_get_wtime() - time00
            write(stdout, '(9x, a40, f14.4)') 'L1 minors                 - time (sec):', time_l1m
            write(stdout, '(9x, a40, f14.4)') 'RS and SR blocks          - time (sec):', time_rs_sr
            write(stdout, '(9x, a40, f14.4)') 'L2 minors and SS block    - time (sec):', time_ss
            write(stdout, '(9x, 40x, a14)') '--------------'
            write(stdout, '(9x, a40, f14.4)') 'Total                     - time (sec):', time_tot
        end if
        if (print_level >= 2) then
            write(stdout, *)
            write(stdout, '(5x,a)') '---- end cis_overlap_l2m subroutine ----'
        end if
    end subroutine cis_overlap_l2m


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: SRBlock
    !----------------------------------------------------------------------------------------------
    subroutine srblock(csc, no, nv1, ns1, wf1, l1minor, sr)
        real(dp), intent(in) :: csc(:, :) !< Molecular orbital overlap matrix.
        integer, intent(in) :: no !< Number of occupied orbitals.
        integer, intent(in) :: nv1 !< Number of bra virtual orbitals.
        integer, intent(in) :: ns1 !< Number of bra states.
        real(dp), intent(in) :: wf1(no*nv1, ns1) !< Bra wf coefficients.
        real(dp), intent(in) :: l1minor(no, no) !< Minors of occupied part of csc matrix.
        real(dp), intent(out) :: sr(ns1)

        integer :: o
        integer :: v
        integer :: i
        integer :: st
        real(dp) :: msum

        sr = 0.0_dp
        !$omp parallel shared (sr)
        !$omp do private(msum) schedule(dynamic) reduction(+:sr)
        do o = 1, no
            do v = 1, nv1
                msum = 0
                do i = 1, no
                    msum = msum + (-1)**(o + i) * csc(no + v, i) * l1minor(o, i)
                end do
                do st = 1, ns1
                    sr(st) = sr(st) + msum * wf1((o-1)*nv1 + v, st)
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine srblock


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: RSBlock
    !----------------------------------------------------------------------------------------------
    subroutine rsblock(csc, no, nv2, ns2, wf2, l1minor, rs)
        real(dp), intent(in) :: csc(:, :) !< Molecular orbital overlap matrix.
        integer, intent(in) :: no !< Number of occupied orbitals.
        integer, intent(in) :: nv2 !< Number of ket virtual orbitals.
        integer, intent(in) :: ns2 !< Number of ket states.
        real(dp), intent(in) :: wf2(no*nv2, ns2) !< Ket wf coefficients.
        real(dp), intent(in) :: l1minor(no, no) !< Minors of occupied part of csc matrix.
        real(dp), intent(out) :: rs(ns2)

        integer :: o
        integer :: v
        integer :: i
        integer :: st
        real(dp) :: msum

        
        rs = 0.0_dp
        !$omp parallel shared(rs)
        !$omp do private(msum) schedule(dynamic) reduction(+:rs)
        do o = 1, no
            do v = 1, nv2
                msum = 0
                do i = 1, no
                    msum = msum + (-1)**(o + i) * csc(i, no + v) * l1minor(i, o)
                end do
                do st = 1, ns2
                    rs(st) = rs(st) + msum * wf2((o-1)*nv2 + v, st)
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine rsblock

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: GetLvl1Minors
    ! The subroutines computes the minors of the matrix per columns or rows depending
    ! on the provided value direction = {r, c}.
    !----------------------------------------------------------------------------------------------
    subroutine getlvl1minors(n, mat, dets, direction)
        integer, intent(in) :: n
        real(dp), intent(in) :: mat(n,n)
        character, intent(in) :: direction
        real(dp), intent(out) :: dets(n,n)
        
        integer :: r
        integer :: c
        integer :: seq(n)
        integer :: cmask(n-1)
        integer :: rmask(n-1)
        real(dp) :: minor(n-1, n-1)

        do r = 1, n
            seq(r) = r
        end do
        do r = 1, n
            rmask = pack(seq, seq/=r)
            do c = 1, n
                cmask = pack(seq, seq/=c)
                if ( direction == 'r') then
                   minor = mat(rmask, cmask)
                   dets(r,c) = mat_ge_det(minor)
                else
                   minor = mat(cmask, rmask)
                   dets(c,r) = mat_ge_det(minor)
                end if
            end do
        end do
    end subroutine getlvl1minors


end module cis_overlap_l2m_mod
