subroutine print_mat(m, n, mat, lda)

        use global_defs

        integer, intent(in) :: m
        integer, intent(in) :: n
        integer, intent(in) :: lda
	real(dp), intent(in) :: mat(lda, n)
        integer :: k

        integer i,j

        do i = 1, m
                do j = 1, n
                        write(*,101, advance="no") mat(i,j), ' ' 
                enddo
                write(*,*)''
        enddo

 101    format(F8.5, A2)

 end subroutine print_mat


! Computing the determinants of the level-2 minors.
! Structure-aware algorithm described in the paper
! Davidovic, Quintana-Orti: "Structure-Aware Calculation of Many-Electron Wave Function Overlaps on Multicore Processors"

subroutine getlvl2minors_lu_csa(n, nv2, csc, ld_csc, l2minor, ld_l2minor)

        use global_defs
!       use mathmod

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
        real(dp), allocatable :: ref(:,:)   ! dimension (n,n)
        real(dp), allocatable :: minor(:,:) ! dimension (n-2,n)
        real(dp), allocatable :: U(:,:)     ! dimension (n-2,n)
        real(dp), allocatable :: Ucopy(:,:) ! dimension (n-2,n-2)
        real(dp) :: det_minor
        real(dp) :: det_U1, det_U2, det_U3
        integer  :: det_P1, det_P2, det_P3

        ! Auxiliary variables
        integer :: num_thr
        integer, external :: omp_get_num_threads
        real(dp), external :: omp_get_wtime
        integer :: i
        integer :: ipiv1(n-2), ipiv2(n-2)
        integer :: info
        real(dp) :: cs, sn, rr

        ! Timing variables
        real(dp) :: startc1, endc1
        real(dp) :: startr1, endr1
        real(dp) :: totalc1, totalr1

        ! External functions
#ifdef OPENBLAS
        external openblas_set_num_threads
        integer, external :: openblas_get_num_threads 
#elif MKL
        integer, external :: mkl_get_max_threads
#endif

        ! External LAPACK function
        external DGETRF

! ************************************************************************************
!
!       Execution part
!
! ************************************************************************************

        l2minor = 0.0_dp

        ! Allocate arrays 
        allocate(ref(n,n))
        allocate(minor(n-2,n))
        allocate(U(n-2,n))
        allocate(Ucopy(n-2,n-2))

        call dlacpy('A', n, n, csc, ld_csc, ref, n)

        ! Set the number of BLAS threads to 1
#ifdef OPENBLAS
        num_thr = openblas_get_num_threads()
        call openblas_set_num_threads(1)
#elif MKL
        num_thr = mkl_get_max_threads()
        call mkl_set_num_threads(1)
#endif

        if (n <= 2) then
            ! Matrices of size <= 2 don't have lvl 2 minors.
            write(*,*) '[getlvl2minors_b] Matrix size <= 2! Returning 1...'
            l2minor = 1.0_dp
            return
        end if

        n2 = n * n
        n3 = n2 * n

        do r1 = n, 1, -1

                startr1 = omp_get_wtime()

                ! Store pre-computed factors (decrease no. operations)
                r1_3 = (r1-1)*n3
                r1_2 = (r1-1)*n2

                ! Position in the output matrix
                pos(1) = r1_2

                do r2 = r1-1, 1, -1

                        ! Store pre-computed factors (decrease no. operations)
                        r2_3 = (r2-1)*n3
                        r2_2 = (r2-1)*n2

                        ! Position in the output matrix
                        pos(2) = pos(1) + r2_3

                        ! Copy from rows 1:r2-1
                        call dlacpy('A', r2-1, n, ref(1,1), n, minor(1,1), n-2)

                        ! Copy rows from r2+1:r1-1
                        call dlacpy('A', r1-r2-1, n, ref(r2+1,1), n, minor(r2,1), n-2)

                        ! Copy rows from r1+1:n
                        call dlacpy('A', n-r1, n, ref(r1+1,1), n, minor(r1-1,1), n-2)

                        ! Compute LU factorization
                        call dgetrf( n-2, n, minor, n-2, ipiv1, info )

                        ! Determinant of the permutation matrix
                        det_P1 = 1;
                        do i = 1, n-2;
                                if (ipiv1(i) /= i) then
                                        det_P1 = -det_P1
                                endif
                        enddo

                        ! Copy U factor to a separate matrix
                        call dlacpy('U', n-2, n, minor, n-2, U, n-2)

                        !$omp parallel default(shared) num_threads(num_thr)
                        !$omp do private(det_minor, Ucopy, c2, pos_e, c2_1, c1_1, &
                        !$omp& det_P2, det_P3, det_U1, det_U3, det_U2, ipiv2, &
                        !$omp& info, i, cs, sn, rr) &
                        !$omp& firstprivate(pos) schedule(dynamic)
                        do c1 = n, 1, -1

                                ! Store pre-computed factors (decrease no. operations)
                                c1_1 = (c1-1)*n

                                ! Position in the output matrix
                                pos(3) = pos(2) + c1_1

                                do c2 = c1-1, 1, -1

                                        ! Copy columns c2+1:c1-1
                                        call dlacpy('A', n-c2-1, c1-c2-1, U(c2,c2+1), n-2, Ucopy(c2,c2), n-2)

                                        ! Copy columns c1+1:n
                                        if(c1 < n) then
                                                call dlacpy('A', n-c2-1, n-c1, U(c2, c1+1), n-2, Ucopy(c2,c1-1), n-2)
                                        endif

                                        ! Product of the first c2-1 diagonal elements of U
                                        det_U1 = 1.0D+0
                                        do i = 1, c2-1
                                                det_U1 = det_U1 * U(i,i)
                                        enddo

                                        ! Annihilate first subdiagonal of the columns between c1 and c2
                                        det_P2 = 1
                                        det_U2 = 1.0D+0

                                        ! Compute and apply givens rotations to columns c2:c1
                                        if(c1-c2-1>0) then
                                        
                                                do i = c2, min(n-3,c1-2)

                                                        ! Givens rotation of the Ucopy(i,i) and Ucopy(i+1,i) element
                                                        call dlartg(Ucopy(i,i), Ucopy(i+1,i), cs, sn, rr)
                                                        call dlasr('L', 'V', 'F', 2, n-i-1, cs, sn, Ucopy(i,i), n-2)
                                                        det_U2 = det_U2 * Ucopy(i,i)
                                                enddo
                                                if(c1==n) det_U2 = det_U2 * Ucopy(i,i)
                                        endif

                                        det_P3 = 1
                                        det_U3 = 1.0D+0

                                        ! Eliminates the two subdiagonals in the columns c1:end
                                        ! Compute and apply givens rotations to columns c1:end
                                        do i = c1-1, n-3

                                                if (i+2 <= n-2 ) then
                                                        ! Remove the second subdiagonal element
                                                        call dlartg(Ucopy(i+1,i), Ucopy(i+2,i), cs, sn, rr)

                                                        ! Update rows i+1 and i+2
                                                        call dlasr('L', 'V', 'F', 2, n-i-1, cs, sn, Ucopy(i+1,i), n-2)
                                                endif

                                                ! Remove the first subdiagonal element
                                                call dlartg(Ucopy(i,i), Ucopy(i+1,i), cs, sn, rr)

                                                ! Update rows i and i+2
                                                call dlasr('L', 'T', 'F', 2, n-i-1, cs, sn, Ucopy(i,i), n-2)

                                                det_U3 = det_U3 * Ucopy(i,i)
                                        enddo
                                        if(c1 < n) det_U3 = det_u3 * Ucopy(n-2,n-2)

                                        det_minor = det_P1 * det_U1 * det_P2 * det_U2 * det_P3 * det_U3
                                        
                                        ! Position in the output matrix holding the determinants
                                        c2_1 = (c2-1)*n
                                        pos(4) = pos(3) + c2

                                        ! Store determinant of the current l2 minor
                                        l2minor(pos(4)) = (-1)**(c1+r2-1) * det_minor

                                        ! Get the other 3 locations where the same determinant resides in the output l2minors matrix
                                        ! There are 4 possible combinations of rows/columns to be remove (r1,r2,c1,c2) resulting the 
                                        ! same l2 minor and same determinants (expect of the sign). Therefor only one combination 
                                        ! is computed. 

                                        ! Get the positions in the output matrix for the other 3 cominations
                                        pos_e(1) = r1_2 + c2_1 + c1 + r2_3 ! (r1, c2, r2, c1)
                                        pos_e(2) = r2_2 + c1_1 + c2 + r1_3 ! (r2, c1, r1, c2)
                                        pos_e(3) = r2_2 + c2_1 + c1 + r1_3 ! (r2, c2, r1, c1)

                                        ! Copy the determinant (with a valid sign) at the computed locations
                                        l2minor(pos_e(1)) = (-1)**(c2+r2) * det_minor
                                        l2minor(pos_e(2)) = (-1)**(c1+r1) * det_minor
                                        l2minor(pos_e(3)) = (-1)**(c2+r1-1) * det_minor
                                end do
                        end do
                        !$omp end do
                        !$omp end parallel
                end do
                endr1 = omp_get_wtime()
                totalr1 = endr1 - startr1
        end do

        ! Restore the original number of threads
#ifdef OPENBLAS
        call openblas_set_num_threads(num_thr)
#elif MKL
        call mkl_set_num_threads(num_thr)
#endif

        ! Deallocate arrays
        deallocate(U)
        deallocate(ref)
        deallocate(minor)
        deallocate(Ucopy)

end subroutine getlvl2minors_lu_csa

