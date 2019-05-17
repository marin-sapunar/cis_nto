!--------------------------------------------------------------------------------------------------
! MODULE: math_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date March, 2019
!--------------------------------------------------------------------------------------------------
module math_mod
    use global_defs
    implicit none


    private
    public :: factorial
    public :: factorial2
    public :: binomial


contains


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: factorial
    !> @brief Product of all numbers up to n (n!).
    !----------------------------------------------------------------------------------------------
    elemental function factorial(n)
        integer, intent(in) :: n
        integer(j15) :: factorial
        integer(j15), parameter :: fact(0:14) = [                                                  &
        &                                        1_j15, 1_j15, 2_j15, 6_j15,                       &
        &                                        24_j15, 120_j15, 720_j15, 5040_j15,               &
        &                                        40320_j15, 362880_j15, 3628800_j15, 39916800_j15, &
        &                                        479001600_j15, 6227020800_j15, 87178291200_j15    &
        &                                        ]
        factorial = fact(n)
    end function factorial


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: factorial2
    !> @brief Product of all odd/even numbers up to n (n!!).
    !----------------------------------------------------------------------------------------------
    elemental function factorial2(n)
        integer, intent(in) :: n
        integer(j15) :: factorial2
        integer(j15), parameter :: dfact(-1:14) = [                                                &
        &                                          1_j15, 1_j15, 1_j15, 2_j15,                     &
        &                                          3_j15, 8_j15, 15_j15, 48_j15,                   &
        &                                          105_j15, 384_j15, 945_j15, 3840_j15,            &
        &                                          10395_j15, 46080_j15, 135135_j15, 645120_j15    &
        &                                          ]
        factorial2 = dfact(n)
    end function factorial2


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: binomial
    !> @brief Binomial coefficient for n over m.
    !----------------------------------------------------------------------------------------------
    function binomial(n, m)
        integer, intent(in) :: n
        integer, intent(in) :: m
        integer(j15) :: binomial
        if (m < 0) then
            binomial = 0
        else if (m > n) then
            binomial = 0
        else
            binomial = factorial(n) / factorial(m) / factorial(n-m)
        end if
    end function binomial


end module math_mod
