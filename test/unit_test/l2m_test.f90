module l2m_test
    use fruit
    use global_defs
    implicit none

contains

    subroutine test_l2m_getlvl2minors
        integer, parameter :: n = 4
        integer, parameter :: nv1 = 0
        integer, parameter :: nv2 = 0
        real(dp) :: csc(n+nv1, n+nv2)
        real(dp) :: l2minor(n**4)
        real(dp), parameter :: tol = 1.E-12_dp
        integer :: i
        external :: getlvl2minors_lu
        
        csc = reshape([ (i, i=1, (n+nv1)*(n+nv2)) ], shape(csc), order=[2, 1])
        
        call getlvl2minors_lu(n, nv2, csc, n+nv1, l2minor, n**3)
        call assert_equals(-4.0_dp, l2minor(18), tol, 'Level 2 minors')
        call assert_equals(-4.0_dp, l2minor(79), tol, 'Level 2 minors')
        call assert_equals(-4.0_dp, l2minor(226), tol, 'Level 2 minors')
        call assert_equals(-4.0_dp, l2minor(191), tol, 'Level 2 minors')
    end subroutine test_l2m_getlvl2minors

end module l2m_test
