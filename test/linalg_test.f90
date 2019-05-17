module linalg_test
    use fruit
    use global_defs
    implicit none

contains

    subroutine test_linalg_wrapper
        use linalg_wrapper_mod
        use matrix_mod
        real(dp) :: inp_mat(2, 2)
        real(dp) :: res_num
        real(dp) :: res_vec(2)
        real(dp) :: res_mat1(2, 2)
        real(dp) :: res_mat2(2, 2)
        real(dp), parameter :: tol = 1.E-12_dp
        real(dp), parameter :: rsq2_mat(2, 2) = reshape([rsq2, -rsq2, rsq2, rsq2], [2, 2])
        real(dp), parameter :: iden_mat(2, 2) = reshape([num1, num0, num0, num1], [2, 2])
        
        inp_mat = rsq2_mat
        
        res_num = dot(inp_mat(1, :), inp_mat(2, :))
        call assert_equals(num0, res_num, tol, 'Dot product')
        
        call gemm(inp_mat, inp_mat, res_mat1)
        res_mat2 = reshape([num0, -num1, num1, num0], [2, 2])
        call assert_equals(res_mat2, res_mat1, 2, 2, tol, 'Matrix multiplication')

        res_num = mat_ge_det(inp_mat)
        call assert_equals(num1, res_num, tol, 'Determinant')

        call gesvd(inp_mat, res_vec, res_mat1, res_mat2)
        call assert_equals([num1, num1], res_vec, 2, tol, 'Matrix SVD')
    end subroutine test_linalg_wrapper

end module linalg_test
