!--------------------------------------------------------------------------------------------------
! MODULE: one_el_op_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date November, 2016
!
!> @brief Contains the subroutine for computing a 1 el. operator over two sets of basis functions.
!--------------------------------------------------------------------------------------------------
module one_el_op_mod
    use global_defs
    implicit none

    private
    public :: one_el_op
    public :: cart_one_el_op


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: one_el_op
    !
    !> @brief Matrix of a simple one electron operator over two sets of basis functions.
    !> @details
    !! Calculates a matrix with elements:
    !!
    !!     S_12(i,j) = < B_1,i | op(x, y, z) | B_2,j >
    !!
    !! The operator to be calculated is given in the operator string either using a named operator
    !! or an analytic expression (eg: op = "x^1 * y^1 + z^2", see read_operator_string subroutine).
    !!
    !! The matrix is initially constructed in the basis of contracted Cartesian Gaussian functions
    !! ordered by 'internal' order (defined by functions given in ang_mom_defs). Afterwards, the
    !! basis is changed to the source format given in the basis_set variables.
    !----------------------------------------------------------------------------------------------
    subroutine one_el_op(geom1, geom2, bs1, bs2, operator_string, mat, ao_center)
        use linalg_wrapper_mod, only : gemm
        use basis_set_mod
        use basis_transform_mod
        real(dp), intent(in) :: geom1(:) !< Geometry 1.
        real(dp), intent(in) :: geom2(:) !< Geometry 2.
        type(basis_set), intent(inout) :: bs1 !< Basis set 1.
        type(basis_set), intent(inout) :: bs2 !< Basis set 2.
        character(len=*), intent(in) :: operator_string !< Operator to calculate.
        real(dp), allocatable, intent(out) :: mat(:, :) !< Operator matrix.
        integer, intent(in) :: ao_center !< Option for treating the atom centers:
                                         !!  -1 - Place all basis functions on geom1
                                         !!   0 - Calculate actual integrals (default)
                                         !!   1 - Place all basis functions on geom2
        real(dp), allocatable :: tmp(:, :)
        integer, allocatable :: lxyz(:, :)
        integer, allocatable :: add(:)
        integer :: clxyz(3), i
        type(basis_transform) :: trans1
        type(basis_transform) :: trans2
        character(len=:), allocatable :: tstring

        select case(operator_string)
        case('overlap')
            tstring = 'x^0'
        case('r^2')
            tstring = 'x^2 + y^2 + z^2'
        case default
            tstring = operator_string
        end select
        call read_operator_string(tstring, lxyz, add)

        ! Calculate matrix in basis of Cartesian GTOs in default order.
        do i = 1, size(lxyz, 1)
            clxyz = lxyz(i, :)
            select case(ao_center)
            case(-1)
                call cart_one_el_op(geom1, geom1, bs1, bs2, clxyz, tmp)
            case(1)
                call cart_one_el_op(geom2, geom2, bs1, bs2, clxyz, tmp)
            case default
                call cart_one_el_op(geom1, geom2, bs1, bs2, clxyz, tmp)
            end select
            if (.not. allocated(mat)) allocate(mat(1:size(tmp, 1), 1:size(tmp, 2)), source=0.0_dp)
            if (add(i) == 1) then
                mat = mat + tmp
            else
                mat = mat - tmp
            end if
        end do

        ! Convert to original basis format.
        call trans1%init(bs1, 'internal', bs1%source_format)
        call trans2%init(bs2, 'internal', bs2%source_format)
        call trans1%transform(mat, 1)
        call trans2%transform(mat, 2)
    end subroutine one_el_op


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_operator_string
    !
    !> @brief Read operator from analytical expression.
    !> @details
    !! Currently, the reader is limited to additions, subtractions and multiplications of terms
    !! of the form q^l, where l must be an integer. Multiplications have to be explicitly written.
    !! Examples of allowed operators:
    !!    x
    !!    x^2
    !!    x^2 + y^2
    !!    x^2 + y^2 * z^2
    !!    x^2 + y^2 * z^2 + ...
    !----------------------------------------------------------------------------------------------
    subroutine read_operator_string(operator_string, lxyz, add)
        use string_mod, only: string, string_parse
        character(len=*), intent(in) :: operator_string !< Operator to calculate.
        integer, allocatable :: lxyz(:, :) !< Exponents for the components of the operator.
        integer, allocatable, intent(out) :: add(:) !< Add (1) or subtract (2) each component.
        type(string), allocatable :: add_substrings(:)
        type(string), allocatable :: mult_substrings(:)
        type(string), allocatable :: exp_substrings(:)
        integer, allocatable :: found_signs(:)
        integer :: add_narg, mult_narg, exp_narg
        integer :: i, j, k, pow, chk
 
        call string_parse(operator_string, '+-', add_narg, add_substrings, found_signs)
        allocate(lxyz(add_narg, 3), source=0)
        allocate(add(add_narg), source=1)
        if (operator_string(1:1) == '-') add(1) = 2
        add(2:) = found_signs
        do i = 1, add_narg
            call string_parse(add_substrings(i)%s, '*', mult_narg, mult_substrings)
            do j = 1, mult_narg
                call string_parse(mult_substrings(j)%s, '^', exp_narg, exp_substrings)
                select case(exp_substrings(1)%s)
                case('x')
                    k = 1
                case('y')
                    k = 2
                case('z')
                    k = 3
                case default
                    if (exp_narg /= 2) then
                        write(stderr, *) 'Error in read_operator_string. Failed to read base.'
                        write(stderr, *) '    Substring: ', mult_substrings(j)%s
                        stop
                    end if
                end select
                select case(exp_narg)
                case(1)
                    pow = 1
                case(2)
                    read(exp_substrings(2)%s, *, iostat=chk) pow
                    if (chk /= 0) then
                        write(stderr, *) 'Error in read_operator_string. Failed to read power.'
                        write(stderr, *) '    Substring: ', mult_substrings(j)%s
                        stop
                    end if
                case default
                    write(stderr, *) 'Error in read_operator_string. Failed to parse.'
                    write(stderr, *) '    Substring: ', mult_substrings(j)%s
                    stop
                end select
                lxyz(i, k) = lxyz(i, k) + pow
            end do
        end do
    end subroutine read_operator_string


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cart_one_el_op
    !
    !> @brief Matrix of a simple one electron operator over two sets of cartesian basis functions.
    !> @details
    !! Calculates a matrix with elements:
    !!
    !!     S_12(i,j) = < B_1,i | x^lx*y^ly*z^lz | B_2,j >
    !!
    !! where < B_1,i | and | B_2,j > are contracnted Cartesian Gaussian functions. 
    !----------------------------------------------------------------------------------------------
    subroutine cart_one_el_op(geom1, geom2, bs1, bs2, lxyz, mat)
        use basis_set_mod
        use ang_mom_defs, only : amp
        use cgto_mod, only : cgto_one_el
        real(dp), intent(in) :: geom1(:) !< Geometry 1.
        real(dp), intent(in) :: geom2(:) !< Geometry 2.
        type(basis_set), intent(inout) :: bs1 !< Basis set 1.
        type(basis_set), intent(inout) :: bs2 !< Basis set 2.
        integer, intent(in) :: lxyz(3) !< List of exponents of the coordinates in the operator.
                                       !! (lx, ly, lz)
        real(dp), allocatable, intent(out) :: mat(:, :) !< Operator matrix.

        integer :: i, j, k, l, ck0, cl0, ck, cl, nbf1, nbf2, bi, bj
        real(dp) :: qi(3), qj(3), qc(3)

        if (allocated(mat)) deallocate(mat)
        allocate(mat(bs1%n_bf_cart, bs2%n_bf_cart), source=num0)

        ck = 0
        do i = 1, bs1%n_center
            bi = bs1%center_i_bs(i)
            ck0 = ck
            cl = 0
            do j = 1, bs2%n_center
                bj = bs2%center_i_bs(j)
                cl0 = cl
                qi = geom1(3*i-2:3*i)
                qj = geom2(3*j-2:3*j)
                qc = num0
                ck = ck0
                do k = 1, bs1%bs(bi)%n_subshell
                    nbf1 = amp%n_cart(bs1%bs(bi)%cg(k)%l)
                    cl = cl0
                    do l = 1, bs2%bs(bj)%n_subshell
                        nbf2 = amp%n_cart(bs2%bs(bj)%cg(l)%l)
                        mat(ck+1:ck+nbf1, cl+1:cl+nbf2) = cgto_one_el(bs1%bs(bi)%cg(k), bs2%bs(bj)%cg(l), qi, qc, qj, lxyz)
                        cl = cl + nbf2
                    end do
                    ck = ck + nbf1
                end do
            end do
        end do
    end subroutine cart_one_el_op


end module one_el_op_mod
