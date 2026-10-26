subroutine getlvl2minors_blk(r1, c1, n, dimmat, mat, l2minors, nthrs)

        use global_defs
        use matrix_mod

        implicit none

        integer, intent(in) :: r1, c1
        integer, intent(in) :: n
        integer, intent(in) :: dimmat
        integer, intent(in) :: nthrs
        real(dp), intent(in) :: mat(dimmat,dimmat)
        real(dp), intent(out) :: l2minors(n*n)
        
        integer :: r2
        integer :: c2
        integer :: r2_new, c1_new
        integer :: pos
        integer :: seq(n)
        integer :: cmask(n-2)
        integer :: rmask(n-2)
        real(dp) :: minor(n-2, n-2)
	real(dp) :: ref(n,n)

! *********************************************************************
!       Executable part       
! *********************************************************************

        ! Set l2minors array to zeros
        l2minors = 0.0_dp

        ! Extract reference matrix from csc (input) matrix
        ref = mat(1:n,1:n)

        if (n <= 2) then
            ! Matrices of size <= 2 don't have lvl 2 minors.
            write(*,*) '[getlvl2minors_b] Matrix size <= 2! Returning 1...'
            l2minors = 1.0_dp
            return
        end if

        ! Auxiliary array to store row/column permutation to 1:n
        do r2 = 1, n
            seq(r2) = r2
        end do

        !$omp parallel default(shared) num_threads(nthrs)
        !$omp do private(pos, r2, cmask, rmask, minor, c1_new, r2_new) schedule(dynamic)
        do c2 = 1, n

            ! Remove only columns not equal to c1 (in the case c1==c2 we have
            ! l1minor)
            if ( c2 /= c1 ) then
                ! Get the position in the new matrix with one column removed
                ! If c2 less then c1 then the index of c1 in the minor matrix is
                ! decreased by one, otherwise it remains
                if ( c2 < c1 ) then
                    c1_new = c1 - 1
                else
                    c1_new = c1
                end if

                ! Create a column mask without c1 and c2 column indices
                cmask = pack(seq, (seq/=c2 .and. seq/=c1))

                do r2 = 1, n
                    ! Remove only rows r2 not equal to r1
                    if ( r2 /= r1 ) then 

                        ! Create a row maks without r1 and r2 row indices
                        rmask = pack(seq, (seq/=r1 .and. seq/=r2))

                        ! Get a new position of row r2 in l2minor (after
                        ! removing row r1)
                        if ( r2 < r1 ) then 
                            r2_new = r2
                        else
                            r2_new = r2 - 1
                        end if

                        ! Compute the position in l2minors matrix to store
                        ! computed det of the minor
                        pos = (r2-1)*n + c2

                        ! Create a minor matrix using row and column masks
                        minor = ref(rmask, cmask)

                        ! Compute the determinant of the given minor and
                        ! multiply with the sign
                        l2minors(pos) = (-1)**(c1_new+r2_new) * mat_ge_det(minor)
                    end if
                end do
            end if
        end do
       !$omp end do
       !$omp end parallel

end subroutine getlvl2minors_blk
