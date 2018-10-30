module matrix_mod
    use global_defs, only : dp, errstop
    implicit none

    private
    public :: genmat_diag
! The following subroutines require LAPACK with Fortran95 interface.
    public :: mat_sy_ev
    public :: mat_sy_exp
    public :: mat_ge_det
    public :: mat_ge_mmm


contains


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: GenMat_Diag
    !> @brief Given a vector return diagonal matrix with values from vector on the diagonal.
    !----------------------------------------------------------------------------------------------
    function genmat_diag(vec) result(mat)
        real(dp), intent(in) :: vec(:)
        real(dp) :: mat(size(vec), size(vec))
        integer :: i

        mat = 0.0_dp
        do i = 1, size(vec)
            mat(i, i) = vec(i)
        end do
    end function genmat_diag



! The following subroutines require LAPACK with Fortran95 interface.
    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: mat_sy_ev
    !> @brief Eigensystem of a symmetric real matrix.
    !----------------------------------------------------------------------------------------------
    subroutine mat_sy_ev(a, evec, eval)
        use lapack95, only : syev
        real(dp), intent(in) :: a(:, :) !< Input matrix.
        real(dp), allocatable, intent(out) :: evec(:, :) !< Eigenvectors of the symmetric matrix.
        real(dp), allocatable, intent(out) :: eval(:) !< Eigenvalues of the symmetric matrix.
        integer :: n
        integer :: info

        n = size(a, 1)
        if (.not. allocated(evec)) allocate(evec, source=a)
        if (.not. allocated(eval)) allocate(eval(n))
 
        ! Calculate eigensystem.
        call syev(a=evec, w=eval, jobz='V', info=info)
        if (info /= 0) call errstop('mat_sy_ev', 'syev call failed.', info)
    end subroutine mat_sy_ev


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: mat_sy_exp
    !
    !> @brief Matrix exponential of a symmetric real matrix multiplied by a complex factor.
    !> @details
    !! The subroutine first calculates the eigenvectors/eigenvalues of the matrix. Then the matrix
    !! exponantial is calculated in the eigenvector basis and the matrix is converted back to the
    !! original basis.
    !----------------------------------------------------------------------------------------------
    function mat_sy_exp(a, mult) result (res)
        use lapack95, only : syev
        use blas95, only : gemm
        real(dp), intent(in) :: a(:, :) !< Input matrix.
        complex(dp), intent(in) :: mult !< Multiplier in exponent.
        complex(dp), allocatable :: res(:, :) !< Result matrix.
        real(dp), allocatable :: evec(:, :) !< Eigenvectors of the symmetric matrix.
        real(dp), allocatable :: eval(:) !< Eigenvalues of the symmetric matrix.
        complex(dp), dimension(:, :), allocatable :: tmp1, tmp2 !< Work matrices
        integer :: i, n
        integer :: info

        n = size(a, 1)
        allocate(evec, source=a)
        allocate(eval(n))
        allocate(tmp1(n, n))
        allocate(tmp2(n, n))
        if (.not. allocated(res)) allocate(res(n, n))
 
        ! Calculate eigensystem.
        call syev(a=evec, w=eval, jobz='V', info=info)
        if (info /= 0) call errstop('mat_sy_exp', 'syev call failed.', info)

        ! Matrix exponential in eigenvector basis.
        res = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, n
            res(i, i) = exp(eval(i) * mult)
        end do

        ! Conversion to original basis.
        tmp1 = cmplx(evec, 0.0_dp, dp)
        call gemm(tmp1, res, tmp2)
        call gemm(tmp2, tmp1, res, transb='C')
    end function mat_sy_exp


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: mat_ge_det
    !> @brief Compute determinant of general matrix.
    !> @details
    !! Uses LU factorization to compute the determinant. The input matrix is not changed.
    !----------------------------------------------------------------------------------------------
    pure function mat_ge_det(a) result(adet)
        use lapack95, only : getrf
        real(dp), intent(in) :: a(:, :) !< Matrix.
        real(dp) :: adet !< Determinant.
        real(dp), allocatable :: tmp(:, :) !< Work array.
        integer, allocatable :: ipiv(:)   ! pivot indices
        integer :: sgn !< Sign change based on pivots.
        integer :: info !< LAPACK exit status.
        integer :: i, n

        n = size(a, 1)
        allocate(tmp, source=a)
        allocate(ipiv(n))

        ! LU factorization using partial pivoting with row interchanges.
        call getrf(tmp, ipiv, info)
        if (info /= 0) then
            adet = 0.0_dp
            return
        end if

        sgn = 1
        adet = 1.0_dp
        do i = 1, n
            if (ipiv(i) /= i) sgn = -sgn
            adet = adet * tmp(i, i)
        end do
        adet = sgn * adet
    end function mat_ge_det


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: mat_ge_mmm
    !> @brief Multiply 3 matrices
    !----------------------------------------------------------------------------------------------
    subroutine mat_ge_mmm(a, b, c, d, transa, transc)
        use blas95, only : gemm
        real(dp) :: a(:, :) !< Input matrix 1.
        real(dp) :: b(:, :) !< Input matrix 2.
        real(dp) :: c(:, :) !< Input matrix 3.
        real(dp) :: d(:, :) !< Output matrix.
        character(len=1), intent(in), optional :: transa !< op(A)
        character(len=1), intent(in), optional :: transc !< op(C)

        real(dp), allocatable :: wrk(:, :)
        character(len=1) :: ta
        character(len=1) :: tc

        ta = 'N'
        tc = 'N'
        if (present(transa)) ta = transa
        if (present(transc)) tc = transc

        if (ta == 'N') then
            allocate(wrk(size(a, 1), size(b, 2)))
        else
            allocate(wrk(size(a, 2), size(b, 2)))
        end if

        call gemm(a, b, wrk, transa=ta)
        call gemm(wrk, c, d, transb=tc)
    end subroutine mat_ge_mmm



end module matrix_mod
