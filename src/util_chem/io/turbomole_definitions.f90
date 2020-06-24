!----------------------------------------------------------------------------------------------
! MODULE: turbomole_definitions_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @brief Hard-coded definitions for turbomole subshell orderings and available functionals.
!----------------------------------------------------------------------------------------------
module turbomole_definitions_mod
    use global_defs
    implicit none


    private
    public :: turbomole_sphe_order
    public :: turbomole_sphe_phase
    public :: turbomole_functional_type


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: turbomole_sphe_order
    !> @brief Return ordered ang. mom. numbers of spherical bfs within subshell of given l.
    !----------------------------------------------------------------------------------------------
    function turbomole_sphe_order(ll) result(lm)
        integer, intent(in) :: ll
        integer, allocatable :: lm(:, :)
        integer :: i

        allocate(lm(2, 2*ll+1))
        lm(1, :) = ll
        select case(ll)
        case(0)
            lm(2, :) = 0
        case(1)
            lm(2, :) = [1, -1, 0]
        case default
            lm(2, 1) = 0
            do i = 1, ll
                if (mod(i, 2) == 0) then
                    lm(2, 2*i) = -i
                    lm(2, 2*i+1) = i
                else
                    lm(2, 2*i) = i
                    lm(2, 2*i+1) = -i
                end if
            end do
        end select
    end function turbomole_sphe_order


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: turbomole_sphe_phase
    !> @brief Fix phase for spherical harmonincs where turbomole has non-standard sign.
    !----------------------------------------------------------------------------------------------
    function turbomole_sphe_phase(ll) result(phase)
        integer, intent(in) :: ll
        real(dp), allocatable :: phase(:)

        allocate(phase(2*ll+1))
        phase = num1
        select case(ll)
        case(3)
            phase(7) = -num1
        case(4)
            phase(5) = -num1
            phase(7) = -num1
        end select
    end function turbomole_sphe_phase


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: turbomole_functional_type
    !
    !> @brief Return type of DFT functional.
    !> @details
    !! The output is an integer based on the type of functional:
    !!  1 - LDA
    !!  2 - GGA
    !!  3 - MGGA
    !!  4 - Hybrid
    !!  5 - ODFT
    !!  6 - Double-hybrid
    !----------------------------------------------------------------------------------------------
    function turbomole_functional_type(func) result (functype)
        character(len=*), intent(in) :: func
        integer :: functype
        character(len=21), dimension(5), parameter :: ldafunc = ['slater-dirac-exchange', &
                                                                 's-vwn                ', &
                                                                 'vwn                  ', &
                                                                 's-vwn_Gaussian       ', &
                                                                 'pwlda                ']
        character(len=14), dimension(7), parameter :: ggafunc = ['becke-exchange', &
                                                                 'b-lyp         ', &
                                                                 'b-vwn         ', &
                                                                 'lyp           ', &
                                                                 'b-p           ', &
                                                                 'pbe           ', &
                                                                 'b97-d         ']
        character(len=4 ), dimension(1), parameter :: mggfunc = ['tpss']
        character(len=15), dimension(8), parameter :: hybfunc = ['bh-lyp         ', &
                                                                 'b3-lyp         ', &
                                                                 'b3-lyp_Gaussian', &
                                                                 'pbe0           ', &
                                                                 'tpshh          ', &
                                                                 'm06            ', &
                                                                 'm06-2x         ', &
                                                                 'pbeh-3c        ']
        character(len=3 ), dimension(2), parameter :: odffunc = ['lhf', &
                                                                 'oep']
        character(len=7 ), dimension(1), parameter :: dhyfunc = ['b2-plyp']
        functype = 0
        if (any(ldafunc == func)) functype = 1
        if (any(ggafunc == func)) functype = 2
        if (any(mggfunc == func)) functype = 3
        if (any(hybfunc == func)) functype = 4
        if (any(odffunc == func)) functype = 5
        if (any(dhyfunc == func)) functype = 6
        if (functype == 0) then
            write(stderr, *) 
            write(stderr, *) 'Error in turbomole_definitions module.'
            write(stderr, *) '  Unrecognized DFT functional.'
            stop
        end if
    end function turbomole_functional_type


end module turbomole_definitions_mod
