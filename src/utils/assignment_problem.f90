!--------------------------------------------------------------------------------------------------
! MODULE: assignment_problem_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date March, 2017
!
! DESCRIPTION: 
!> @brief Subroutine for solving the assignment problem.
!> @details
!! The (linear) assignment problem:
!!
!! The problem consists of n agents and n tasks. Any agent can be assigned to perform any task,
!! incurring a cost depending on the agent-task assignment. It is required to perform all tasks 
!! by assigning exactly one agent to each task and exactly one task to each agent in such a way 
!! that the total cost (sum of agent-task costs) of the assignment is minimized.
!!
!! A modified form of the Hungarian method (Kuhn–Munkres algorithm) is used.
!! Adapted from: http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html
!!
!! References:
!! J. Munkres, "Algorithms for the Assignment and Transportation Problems", J. SIAM, 5(1),
!!              32–38, 1957 March.
!--------------------------------------------------------------------------------------------------
module assignment_problem_mod
    use global_defs
    implicit none

    private
    public :: assignment_problem
    public :: getloops
    public :: getloop


    real(dp), parameter :: tinydp = 1.e-15_dp


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: assignment_problem
    !
    ! DESCRIPTION:
    !> @brief The Kuhn-Munkres algorithm for square matrices.
    !> @details
    !! Input is a real matrix with dimensions n x n. The output is a vector containing the optimal
    !! permutation of columns associated to each row of the matrix (e. g. element 1 of the vector
    !! is the column associated with row 1 of the matrix).
    !!
    !! The algorithm is split into steps, each contained in a separate helper subroutine.
    !----------------------------------------------------------------------------------------------
    subroutine assignment_problem(n, m, row2col, transform)
        integer, intent(in) :: n
        real(dp), intent(in) :: m(n, n)
        integer, intent(out) :: row2col(n)
        integer, intent(out), optional :: transform(n, n)

        integer :: step, i, j
        real(dp) :: cm(n, n)
        logical :: zero(n, n)
        logical :: rowc(n)
        logical :: colc(n)
        logical :: star(n, n)
        logical :: prim(n, n)
        integer :: p(2)

        cm = m ! The input matrix isn't modified by the suboutine.
        step = 1
        colc = .false.
        rowc = .false.
        prim = .false.
        star = .false.
        do 
            select case(step)
            case(1)
                call step1(step, n, cm, zero)
            case(2)
                call step2(step, n, zero, star)
            case(3)
                call step3(step, n, colc, star)
            case(4)
                call step4(step, n, zero, colc, rowc, star, prim, p)
            case(5)
                call step5(step, n, colc, rowc, star, prim, p)
            case(6)
                call step6(step, n, cm, zero, colc, rowc)
            case(7)
                call step7(n, star, row2col)
                exit
            end select
        end do

        if (present(transform)) then
            transform = 0
            do i=1,n
               do j=1,n
                  if (star(i,j)) transform(j,i) = 1
               end do
            end do
        end if
    end subroutine assignment_problem


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Step1
    !
    !> @details
    !! For each row of the matrix, find the smallest element and subtract it from every element in 
    !! its row. Repeat for columns. Create a logical matrix which is .true. where the value of the
    !! input matrix is zero. Go to Step 2. 
    !----------------------------------------------------------------------------------------------
    subroutine step1(step, n, cm, zero)
        integer, intent(out) :: step
        integer, intent(in) :: n
        real(dp), intent(inout) :: cm(n, n)
        logical, intent(out) :: zero(n, n)
        integer :: i

        do i = 1, n
            cm(i, :) = cm(i, :) - minval(cm(i, :))
        end do
        do i = 1, n
            cm(:, i) = cm(:, i) - minval(cm(:, i))
        end do
        zero = (cm < tinydp)
        step = 2
    end subroutine step1


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Step2
    !
    !> @details
    !! Cycle through the zeroes. If there is no starred zero in the row or column of a zero, star 
    !! it. Go to Step 3.
    !----------------------------------------------------------------------------------------------
    subroutine step2(step, n, zero, star)
        integer, intent(out) :: step
        integer, intent(in) :: n
        logical, intent(in) :: zero(n, n)
        logical, intent(out) :: star(n, n)
        integer :: i
        integer :: j

        star = .false.
        do i = 1, n
            if (any(star(i, :))) cycle
            do j = 1, n
                if (.not. zero(i,j)) cycle
                if (any(star(:, j))) cycle
                star(i, j) = .true.
                exit
            end do
        end do
        step = 3
    end subroutine step2


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Step3
    !
    !> @details
    !! Cover each column containing a starred zero. If n columns are covered, the starred zeros 
    !! describe a complete set of unique assignments (go to step 7). If not, go to step 4.  
    !----------------------------------------------------------------------------------------------
    subroutine step3(step, n, colc, star)
        integer, intent(out) :: step
        integer, intent(in) :: n
        logical, intent(inout) :: colc(n)
        logical, intent(inout) :: star(n, n)
        integer :: i
        integer :: k

        k = 0
        do i = 1, n
            if (.not. any(star(:, i))) cycle
            colc(i) = .true.
            k = k + 1
        end do

        if (k == n) then
            step = 7
        else
            step = 4
        end if
    end subroutine step3


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Step4
    !
    !> @details
    !! Find a noncovered zero and prime it.  If there is no starred zero in the row containing the 
    !! primed zero, mark the position of the primed zero and Go to Step 5.  Otherwise, cover this 
    !! row and uncover the column containing the starred zero. Repeat. If there are no uncovered 
    !! zeros left, go to Step 6.
    !----------------------------------------------------------------------------------------------
    subroutine step4(step, n, zero, colc, rowc, star, prim, p)
        integer, intent(out) :: step
        integer, intent(in) :: n
        logical, intent(in) :: zero(n, n)
        logical, intent(inout) :: colc(n)
        logical, intent(inout) :: rowc(n)
        logical, intent(in) :: star(n, n)
        logical, intent(inout) :: prim(n, n)
        integer, intent(out) :: p(2)
        integer :: i
        integer :: j
        logical :: starincol

        do
            p = [0, 0]
            starincol = .false.
outer:      do i = 1, n
                if (rowc(i)) cycle
inner:          do j = 1, n
                    if (colc(j)) cycle
                    if (.not. zero(i,j)) cycle
                    prim(i,j) = .true.
                    p = [i, j]
                    exit outer
                end do inner
            end do outer
            
            if (p(1) == 0) then
                step = 6
                return
            end if
            

            do i = 1, n
                if (.not. star(p(1), i)) cycle
                colc(i) = .false.
                rowc(p(1)) = .true.
                starincol = .true.
            end do

            if (.not. starincol) then
                step = 5
                return
            end if
        end do

    end subroutine step4


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Step5
    !
    !> @details
    !! Starting at the final primed zero of Step 4, move through the matrix in the following way: 
    !! Go to a starred zero in the same column as the primed zero (if possible). Go to a primed
    !! zero in the same row as the starred zero (always possible). Along the way, unstar all 
    !! starred zeros and star all primed zeros. When no starred zeros are present in the column of
    !! a primed zero, uncover all lines of the matrix, unprime all zeros and go to Step 3.
    !----------------------------------------------------------------------------------------------
    subroutine step5(step, n, colc, rowc, star, prim, p)
        integer, intent(out) :: step
        integer, intent(in) :: n
        logical, intent(out) :: colc(n)
        logical, intent(out) :: rowc(n)
        logical, intent(inout) :: star(n, n)
        logical, intent(inout) :: prim(n, n)
        integer, intent(inout) :: p(2)
        integer :: i
        logical :: starincol

        do
            starincol = .false.
            do i = 1, n
                if (.not. star(i, p(2))) cycle
                starincol = .true.
                star(p(1), p(2)) = .true.
                p(1) = i
                exit
            end do
            if (.not. starincol) then
                star(p(1), p(2)) = .true.
                prim = .false.
                rowc = .false.
                colc = .false.
                step = 3
                return
            end if
            do i = 1, n
                if (.not. prim(p(1), i)) cycle
                star(p(1), p(2)) = .false.
                p(2) = i
                exit
            end do
        end do
    end subroutine step5


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Step6
    !
    !> @details
    !! Find the smallest uncovered value. Add this value to every element of each covered row, and
    !! subtract it from every element of each uncovered column. Return to Step 4.
    !----------------------------------------------------------------------------------------------
    subroutine step6(step, n, cm, zero, colc, rowc)
        integer, intent(out) :: step
        integer, intent(in) :: n
        real(dp), intent(inout) :: cm(n, n)
        logical, intent(out) :: zero(n, n)
        logical, intent(in) :: colc(n)
        logical, intent(in) :: rowc(n)
        real(dp) :: mval
        logical :: mask(n, n)
        integer :: i

        mask = .true.
        do i = 1, n
            if (rowc(i)) mask(i, :) = .false.
            if (colc(i)) mask(:, i) = .false.
        end do
        mval = minval(cm, mask)

        do i = 1, n
            if (rowc(i)) cm(i, :) = cm(i, :) + mval
            if (.not. colc(i)) cm(:, i) = cm(:, i) - mval
        end do
        zero = (cm < tinydp)
        step = 4
    end subroutine step6


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: Step7
    !
    !> @details
    !! Construct the output vector from the matrix of starred zeros.
    !----------------------------------------------------------------------------------------------
    subroutine step7(n, star, row2col)
        integer, intent(in) :: n
        logical, intent(in) :: star(n, n)
        integer, intent(out) :: row2col(n)
        integer :: i
        integer :: j

        do i = 1, n
            do j = 1, n
                if (star(i, j)) row2col(i) = j
            end do
        end do
    end subroutine step7


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: GetLoops
    !> @brief Find all loops in given list of indexes.
    !----------------------------------------------------------------------------------------------
    subroutine getloops(list, loops)
        use ivec_mod
        integer, intent(in) :: list(:)
        type(ivec), allocatable, intent(out) :: loops(:)
        type(ivec) :: tmp(size(list))
        logical :: mask(size(list))
        integer :: i, n

        mask = .true.
        n = 0
        do i = 1, size(list)
            if (.not. mask(i)) cycle
            n = n + 1
            call getloop(list, i, tmp(n)%c, mask)
        end do

        if (allocated(loops)) deallocate(loops)
        allocate(loops(n))
        do i = 1, n
            loops(i) = tmp(i)
        end do
    end subroutine getloops


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: GetLoop
    !> @brief Follow list of indexes until it loops to initial index.
    !----------------------------------------------------------------------------------------------
    subroutine getloop(list, init, loop, mask)
        integer, intent(in) :: list(:) !< List of indexes to follow.
        integer, intent(in) :: init !< Initial position in list.
        integer, allocatable, intent(out) :: loop(:) !< List of indexes until init is found.
        logical :: mask(size(list))
        integer :: i, n

        n = 1
        i = init
        getn: do
            mask(i) = .false.
            if (list(i) == init) exit getn
            i = list(i)
            n = n + 1
            if (i > size(list)) then
                write(stderr, *) 'Error in getloop subroutine. Index out of bounds.'
                call abort()
            end if
            if (n > size(list)) then
                write(stderr, *) 'Error in getloop subroutine. Bad input list.'
                call abort()
            end if
        end do getn

        if (allocated(loop)) deallocate(loop)
        allocate(loop(n))
        i = init
        fill: do n = 1, size(loop)
            loop(n) = i
            if (list(i) == init) exit fill
            i = list(i)
        end do fill
    end subroutine getloop


end module assignment_problem_mod
