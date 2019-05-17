module ang_mom_defs_test
    use fruit
    use global_defs
    use ang_mom_defs
    implicit none

contains


    subroutine test_spherical_to_cartesian_gaussian
        real(dp), parameter :: tol = 1.e-12_dp
        real(dp), parameter :: s3o4 = sqrt(3.0_dp/4.0_dp)
        real(dp), parameter :: s3o8 = sqrt(3.0_dp/8.0_dp)
        real(dp), parameter :: s3o40 = sqrt(3.0_dp/40.0_dp)
        real(dp), parameter :: s5o8 = sqrt(5.0_dp/8.0_dp)
        real(dp), parameter :: s6o5 = sqrt(6.0_dp/5.0_dp)
        real(dp), parameter :: s9o8 = sqrt(9.0_dp/8.0_dp)
        real(dp), parameter :: s9o20 = sqrt(9.0_dp/20.0_dp)
        real(dp), parameter :: p_trans(3,3) = reshape( [ &
        &    num1, num0, num0, &
        &    num0, num1, num0, &
        &    num0, num0, num1  &
        &    ], shape(p_trans))
        real(dp), parameter :: d_trans(6,5) = reshape( [ &
        &   -f1o2,-f1o2, num1, num0, num0, num0, &
        &    num0, num0, num0, num0, num1, num0, &
        &    num0, num0, num0, num0, num0, num1, &
        &    s3o4,-s3o4, num0, num0, num0, num0, &
        &    num0, num0, num0, num1, num0, num0  &
        &    ], shape(d_trans))
        real(dp), parameter :: f_trans(10, 7) = reshape( [ &
        &     num0, num0, num1,  num0,-s9o20,  num0,-s9o20, num0, num0, num0, &
        &    -s3o8, num0, num0,  num0,  num0,-s3o40,  num0, s6o5, num0, num0, &
        &     num0,-s3o8, num0,-s3o40,  num0,  num0,  num0, num0, s6o5, num0, &
        &     num0, num0, num0,  num0,  s3o4,  num0, -s3o4, num0, num0, num0, &
        &     num0, num0, num0,  num0,  num0,  num0,  num0, num0, num0, num1, &
        &     s5o8, num0, num0,  num0,  num0, -s9o8,  num0, num0, num0, num0, &
        &     num0,-s5o8, num0,  s9o8,  num0,  num0,  num0, num0, num0, num0  &
        &     ], shape(f_trans))

        call amp%init(4)
        ! S
        call assert_equals(num1, amp%sphe_cart_trans(0)%c(1, 1), tol, 'S orbital transformation')
        ! P
        call assert_equals(p_trans, amp%sphe_cart_trans(1)%c, 3, 3, tol, 'P orbital transformation')
        ! D
        call assert_equals(d_trans, amp%sphe_cart_trans(2)%c, 5, 6, tol, 'D orbital transformation')
        ! F
        call assert_equals(f_trans, amp%sphe_cart_trans(3)%c, 10, 7, tol, 'F orbital transformation')
    end subroutine test_spherical_to_cartesian_gaussian


end module ang_mom_defs_test
