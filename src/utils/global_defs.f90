!--------------------------------------------------------------------------------------------------
! MODULE: Global_defs
!--------------------------------------------------------------------------------------------------
module global_defs
    use, intrinsic :: iso_fortran_env, only : &
        stderr => error_unit, &
        stdout => output_unit, &
        stdin => input_unit, &
        sp => real32, & 
        dp => real64, & 
        qp => real128 

! intel v15.0 doesn't print backtrace with call abort()
#if __INTEL_COMPILER
    interface 
        subroutine abort() bind(C, name="abort")
        end subroutine
    end interface
#endif


    !----------------------------------------------------------------------------------------------
    ! TYPE: IVec
    !> @brief Array of integer vectors of varying dimensions.
    !----------------------------------------------------------------------------------------------
    type ivec 
       integer, allocatable :: c(:) !< Coefficents of the vectors.
    end type ivec


    !----------------------------------------------------------------------------------------------
    ! TYPE: IMat
    !> @brief Array of integer matrices of varying dimensions.
    !----------------------------------------------------------------------------------------------
    type imat
        integer, allocatable :: c(:, :) !< Coefficients of the matrices.
    end type imat


    !----------------------------------------------------------------------------------------------
    ! TYPE: RVec
    !> @brief Array of real vectors of varying dimensions.
    !----------------------------------------------------------------------------------------------
    type rvec 
       real(dp), allocatable :: c(:) !< Coefficents of the vectors.
    end type rvec


    !----------------------------------------------------------------------------------------------
    ! TYPE: RMat
    !> @brief Array of real matrices of varying dimensions.
    !----------------------------------------------------------------------------------------------
    type rmat
        real(dp), allocatable :: c(:, :) !< Coefficients of the matrices.
    end type rmat


    !----------------------------------------------------------------------------------------------
    ! TYPE: LVec
    !> @brief Array of logical vectors of varying dimensions.
    !----------------------------------------------------------------------------------------------
    type lvec
       logical, allocatable :: l(:)
    end type lvec


    !----------------------------------------------------------------------------------------------
    ! TYPE: LMat
    !> @brief Array of logical matrices of varying dimensions.
    !----------------------------------------------------------------------------------------------
    type lmat
       logical, allocatable :: l(:, :)
    end type lmat


contains


    subroutine errstop(errsub, errmsg, errnum)
        character(len=*), intent(in) :: errsub !< Subroutine where the error occured.
        character(len=*), intent(in), optional :: errmsg !< Error message.
        integer, intent(in), optional :: errnum !< Error number.

        write(stderr, '(1x,a)') 'Error in '//errsub//'.'
        if (present(errmsg)) write(stderr, '(3x, a, a)')  'Error: ', errmsg
        if (present(errnum)) write(stderr, '(3x, a, i0)')  'Status: ', errnum
        call abort()
    end subroutine errstop


end module global_defs
