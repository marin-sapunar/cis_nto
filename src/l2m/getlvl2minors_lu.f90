subroutine getlvl2minors_lu(n, nv2, csc, ld_csc, l2minor, ld_l2minor)

        use global_defs
        use matrix_mod

        implicit none

        integer, intent(in) :: n
        integer, intent(in) :: nv2
        integer, intent(in) :: ld_csc, ld_l2minor
        real(dp), intent(in) :: csc(ld_csc, n+nv2)
        real(dp), intent(out) :: l2minor(ld_l2minor*n)
        
        integer :: r1, r2
        integer :: c1, c2
        integer :: tmp, c11
        integer :: n2, n3
        integer :: r1_2, r1_3
        integer :: r2_2, r2_3
        integer :: c1_1, c2_1
        integer :: c1_sign
        integer :: c2p
        integer :: pos(4)
        integer :: pos_e(3)
        integer :: seq(n), seq_cols(n)
        integer :: cmask(n-2)
        integer :: rmask(n-2)
        real(dp) :: minor(n-2, n-2)
        real(dp) :: ref_minor(n-2, n-2)
        real(dp) :: ref(n,n)
        real(dp) :: det_minor
        real(dp) :: det_U1, det_U2
        integer  :: det_P1, det_P2

        ! Auxiliary variables
        integer :: num_thr
        integer, external :: omp_get_num_threads
        real(dp), external :: omp_get_wtime
        integer :: i
        integer :: num_cols, nc_minor
        integer :: ipiv1(n-2), ipiv2(n-2)
        integer :: info

        ! Timing variables
        real(dp) :: startc1, endc1
        real(dp) :: startr1, endr1
        real(dp) :: totalc1, totalr1

        ! External functions
        integer, external :: mkl_get_max_threads

        ! External LAPACK function
        external DGETRF

! ************************************************************************************
!
!       Execution part
!
! ************************************************************************************

        l2minor = 0.0_dp

        ref = csc(1:n,1:n)

        ! Set the number of BLAS threads to 1
        num_thr = mkl_get_max_threads()
        call mkl_set_num_threads(1)

        if (n <= 2) then
            ! Matrices of size <= 2 don't have lvl 2 minors.
            write(*,*) '[getlvl2minors_b] Matrix size <= 2! Returning 1...'
            l2minor = 1.0_dp
            return
        end if
        do r1 = 1, n
            seq(r1) = r1
        end do

        n2 = n * n
        n3 = n2 * n

        write(*,*) 'Getlvl2minors - computing minors using updated LU factorization'
        write(*,*) '____________________________'
        write(*,101) 'Outer loop pass', 'time(sec)'
        write(*,*) '____________________________'
 101    format(A17, A11)

        do r1 = n, 1, -1

                startr1 = omp_get_wtime()

                ! Remove row r1

                ! Store pre-computed factors (decrease no. operations)
                r1_3 = (r1-1)*n3
                r1_2 = (r1-1)*n2

                ! Position in the output matrix
                pos(1) = r1_2

                do r2 = r1-1, 1, -1

                        ! Mask rows be removed
                        rmask = pack(seq, (seq/=r1 .and. seq/=r2))

                        ! Remove row r2

                        ! Store pre-computed factors (decrease no. operations)
                        r2_3 = (r2-1)*n3
                        r2_2 = (r2-1)*n2

                        ! Position in the output matrix
                        pos(2) = pos(1) + r2_3
                        !pos(2) = pos(1)

                        !$omp parallel default(shared) num_threads(num_thr)
                        !$omp do private(det_minor, minor, c2, pos_e, c2_1, c1_1, c1_sign, &
                        !$omp& nc_minor, det_P2, det_P1, det_U1, det_U2, c2p, ipiv1, ipiv2, &
                        !$omp& cmask, seq_cols, info, i, ref_minor, num_cols) &
                        !$omp& firstprivate(pos) schedule(dynamic)
                        do c1 = n, 1, -1

!                                startc1 = omp_get_wtime()
                                ! Move last n-c1 columns to the beginning of the ref matrix
                                ! Update only the seq_cols variable holding the final columns permutation
                                seq_cols(1:n-c1) = (/ (i, i=c1+1,n ) /)
                                seq_cols(n-c1+1:n) = (/ (i, i=1, c1) /)

                                ! Store sign after permutation (the determinant sign is now changed)
                                num_cols = n-c1
                                c1_sign = (-1)**((n-num_cols)*num_cols)

                                ! ALTERNATIVE
                                !c1_sign = 1
                                !if ( MOD((n-num_cols)*num_cols, 2) /= 0 ) then
                                !       c1_sign = -1
                                !endif

                                ! Remove column c1 (effectively, remove the last column)

                                ! Store pre-computed factors (decrease no. operations)
                                c1_1 = (c1-1)*n

                                ! Position in the output matrix
                                pos(3) = pos(2) + c1_1

                                ! l2 minor number of columns (in practice it should be = num_cols-2)
                                nc_minor = n - 2

                                do c2 = c1-1, 1, -1

                                        ! Now, column c2 is on position no-c1+c2
                                        ! (no-c1 columns are inserted at the beginning)
                                        c2p = n - c1 + c2

                                        ! Mask columns that will be removed
                                        !cmask = pack(seq_cols, (seq_cols/=n .and. seq_cols/=c2p))
                                        cmask = pack(seq_cols, (seq_cols/=c1 .and. seq_cols/=c2))
                                        
                                        ! Create minor
                                        minor = ref(rmask, cmask) 
        
                                        ! Compute determinant of the minor.
                                        ! Two cases:
                                        ! In the first case (c2 == c1-1) c2 is the last column
                                        ! Compute the LU of the base case. It is a refernce LU factorization
                                        if ( c2 == c1 - 1 ) then

                                                ! Compute LU factorization
                                                call dgetrf( nc_minor, nc_minor, minor, nc_minor, ipiv1, info )

                                                ! Get determinant of the minor
                                                det_P1 = 1
                                                det_U1 = 1.0D+0
                                                do i = 1, nc_minor
                                                        if (ipiv1(i) /= i) det_P1 = -det_P1
                                                        det_U1 = det_U1 * minor(i, i)
                                                end do
                                                det_minor = det_P1 * det_U1
                                                ref_minor = minor

                                        else
                                                ! Compute for other cases, i.e. when c2 < c1-1. Update with L1

                                                ! Update the last no - c2 columns with L1
                                                ! First apply row interchanges to the last no-c2 columns
                                                call dlaswp( nc_minor-c2p+1, minor(1,c2p), nc_minor, 1, nc_minor, ipiv1, 1 )

                                                call dtrsm( 'L', 'L', 'N', 'U', nc_minor, nc_minor-c2p+1, 1.0D+0, ref_minor, nc_minor, minor(1,c2p), nc_minor )

                                                ! Compute LU of last nc_minor - c2p columns
                                                call dgetrf( nc_minor-c2p+1, nc_minor-c2p+1, minor(c2p, c2p), nc_minor, ipiv2(c2p), info )

                                                ! Get determinant of the minor
                                                det_P2 = 1
                                                det_U2 = 1.0D+0
                                                do i = c2p, nc_minor
                                                        if(ipiv2(i) /= i-c2p+1) det_P2 = -det_P2
                                                        det_U2 = det_U2 * minor(i, i)
                                                end do

                                                ! Do not compute det of U1 from scratch, as the c2 moves towards the beginning the
                                                ! U1 shrinks one columns at the time. Therefore, from the diagonal product only
                                                ! remove the factor of the last column (the one removed in this step)
                                                if ( c2p > 1 ) then
                                                        det_U1 = det_U1 / ref_minor(c2p,c2p)
                                                else
                                                        det_U1 = 1
                                                endif

                                                det_minor = det_P2 * det_P1 * det_U1 * det_U2
                                        endif
                                        
                                        ! Apply sign generated after column permutation at the beginning of c1 loop (moving c1
                                        ! column to the end)

                                        ! Final determinant of the minor
                                        det_minor = c1_sign * det_minor
                                        
                                        c2_1 = (c2-1)*n
                                        pos(4) = pos(3) + c2

                                        ! Store determinant of the current minor
                                        l2minor(pos(4)) = (-1)**(c1+r2-1) * det_minor

                                        ! Get the other 3 locations where the same determinant resides in the output l2minors matrix
                                        pos_e(1) = r1_2 + c2_1 + c1 + r2_3 ! (r1, c2, r2, c1)
                                        pos_e(2) = r2_2 + c1_1 + c2 + r1_3 ! (r2, c1, r1, c2)
                                        pos_e(3) = r2_2 + c2_1 + c1 + r1_3 ! (r2, c2, r1, c1)

                                        ! Copy the determinant (with a valid sign) at the computed locations
                                        l2minor(pos_e(1)) = (-1)**(c2+r2) * det_minor
                                        l2minor(pos_e(2)) = (-1)**(c1+r1) * det_minor
                                        l2minor(pos_e(3)) = (-1)**(c2+r1-1) * det_minor
                                end do
!                                endc1 = omp_get_wtime()
!                                totalc1 = endc1 - startc1
!                                write(*,101) '[getlvl2minors] c1 ', c1, 'time (sec): ', totalc1                                                       
                        end do
                        !$omp end do
                        !$omp end parallel
                end do
                endr1 = omp_get_wtime()
                totalr1 = endr1 - startr1
                !write(*,102) '[getlvl2minors] r1: ', r1, ' time (sec): ', totalr1
                write(*,102) r1, ' ', totalr1
        end do

 102    FORMAT(I10, A9, F7.4)
        write(*,*) '____________________________'
        write(*,*) ''

        ! Restore the original number of threads
        call mkl_set_num_threads(num_thr)

end subroutine getlvl2minors_lu

