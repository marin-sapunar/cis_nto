module assignment_test
    use fruit
    use global_defs
    use assignment_problem_mod
    implicit none


contains


    subroutine test_assignment_square
        integer, parameter :: n = 4
        real(dp) :: test_mat(n, n)
        integer, allocatable :: res_row(:)
        integer, allocatable :: res_col(:)
        
        test_mat(1, :) = [1, 0, 0, 0]
        test_mat(2, :) = [0, 0, 0, 1]
        test_mat(3, :) = [0, 1, 0, 0]
        test_mat(4, :) = [0, 0, 1, 0]
        call assignment_problem(-test_mat, res_row, res_col)
        call assert_equals([1, 2, 3, 4], res_row, n, 'Assigned rows vector.')
        call assert_equals([1, 4, 2, 3], res_col, n, 'Assigned columns vector.')

        test_mat(1, :) = [1, 2, 3, 4]
        test_mat(2, :) = [2, 4, 6, 8]
        test_mat(3, :) = [3, 6, 9, 12]
        test_mat(4, :) = [4, 8, 12, 16]
        call assignment_problem(test_mat, res_row, res_col)
        call assert_equals([1, 2, 3, 4], res_row, n, 'Assigned rows vector.')
        call assert_equals([4, 3, 2, 1], res_col, n, 'Assigned columns vector.')
    end subroutine test_assignment_square


    subroutine test_assignment_rectangular_more_rows
        integer, parameter :: m = 5
        integer, parameter :: n = 3
        real(dp) :: test_mat(m, n)
        integer, allocatable :: res_row(:)
        integer, allocatable :: res_col(:)
        
        test_mat(1, :) = [1, 0, 0]
        test_mat(2, :) = [0, 0, 0]
        test_mat(3, :) = [0, 0, 1]
        test_mat(4, :) = [0, 1, 0]
        test_mat(5, :) = [0, 0, 0]
        call assignment_problem(-test_mat, res_row, res_col)
        call assert_equals([1, 3, 4], res_row, n, 'Assigned rows vector.')
        call assert_equals([1, 3, 2], res_col, n, 'Assigned columns vector.')
    end subroutine test_assignment_rectangular_more_rows


    subroutine test_assignment_rectangular_more_cols
        integer, parameter :: m = 3
        integer, parameter :: n = 5
        real(dp) :: test_mat(m, n)
        integer, allocatable :: res_row(:)
        integer, allocatable :: res_col(:)
        
        test_mat(1, :) = [1, 0, 0, 0, 0]
        test_mat(2, :) = [0, 0, 0, 1, 0]
        test_mat(3, :) = [0, 0, 1, 0, 0]
        call assignment_problem(-test_mat, res_row, res_col)
        call assert_equals([1, 2, 3], res_row, m, 'Assigned rows vector.')
        call assert_equals([1, 4, 3], res_col, m, 'Assigned columns vector.')
    end subroutine test_assignment_rectangular_more_cols


end module assignment_test
