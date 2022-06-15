!----------------------------------------------------------------------------------------------
! MODULE: molden_definitions_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date June, 2017
!
!> @brief Hard-coded definitions for molden subshell orderings.
!----------------------------------------------------------------------------------------------
module molden_definitions_mod
    use global_defs
    use file_mod, only : reader
    implicit none

    private

    public :: molden_sphe_order
    public :: molden_cart_order


contains


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: molden_cart_order
    !> @brief Return ordered ang. mom. numbers of cartesian bfs within subshell of given l.
    !----------------------------------------------------------------------------------------------
    function molden_cart_order(ll) result(lxyz)
        use ang_mom_defs, only : amp, cart_nums
        integer, intent(in) :: ll
        integer, allocatable :: lxyz(:, :)
        integer, parameter :: f_order(3,10) = reshape( [ &
        &      3, 0, 0, 1, 2, 2, 1, 0, 0, 1, &
        &      0, 3, 0, 2, 1, 0, 0, 1, 2, 1, &
        &      0, 0, 3, 0, 0, 1, 2, 2, 1, 1], &
        &     shape(f_order), order=[2,1])

        allocate(lxyz(3, amp%n_cart(ll)))
        select case(ll)
        case(3)
            lxyz = f_order
        case default
            lxyz = cart_nums(ll)
        end select
    end function molden_cart_order


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: molden_sphe_order
    !> @brief Return ordered ang. mom. numbers of spherical bfs within subshell of given l.
    !----------------------------------------------------------------------------------------------
    function molden_sphe_order(ll) result(lm)
        integer, intent(in) :: ll
        integer, allocatable :: lm(:, :)
        integer :: i

        allocate(lm(2, 2*ll+1))
        lm(1, :) = ll
        select case(ll)
        case(0)
            lm(2, :) = [0]
        case(1)
            lm(2, :) = [1, -1, 0]
        case default
            lm(2, 1) = 0
            do i = 1, ll
                lm(2, 2*i) = i
                lm(2, 2*i+1) = -i
            end do
        end select
    end function molden_sphe_order


end module molden_definitions_mod
