!--------------------------------------------------------------------------------------------------
! MODULE: basis_transform_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date April, 2019
!
!> @brief Contains the basis_transform type for transforming Gaussian basis sets.
!--------------------------------------------------------------------------------------------------
module basis_transform_mod
    use global_defs
    use cgto_mod
    use ang_mom_defs
    use basis_set_mod
    implicit none
 

    private
    public :: basis_transform


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: basis_transform
    !> @brief Transformation between different spherical/Cartesian Gaussian basis sets.
    !> @details
    !! Different formats expand the basis set using different ordering of the m (spherical) or lx,
    !! ly, lz (Cartesian) quantum numbers. This type is used to generate matrices for transforming
    !! vectors from one format to another.
    !!
    !! Calling the init procedure initializes the transformation matrix after which the transform
    !! procedure is used to transform vectors/arrays by multiplying them with the transformation
    !! matrix.
    !----------------------------------------------------------------------------------------------
    type basis_transform
        logical :: initialized = .false. !< Instance initialized?
        character(len=:), allocatable :: source_format !< Format from which to convert
        character(len=:), allocatable :: target_format !< Format to which to convert
        logical :: source_sphe !< Spherical/cartesian for source format
        logical :: target_sphe !< Spherical/cartesian for target format
        real(dp), allocatable :: trans(:, :) !< Transformation matrix
        real(dp), allocatable :: sphe_cart_trans(:, :) !< Sphe/cart basis transformation
        real(dp), allocatable :: source_trans(:, :) !< Transform source format to internal format
        real(dp), allocatable :: target_trans(:, :) !< Transform target format to internal format
        integer :: n1 = 0 !< Number of source basis functions
        integer :: n2 = 0 !< Number of target basis functions
    contains
        procedure :: init
        generic :: transform => transform_array, transform_arrays
        procedure, private :: transform_array, transform_arrays
    end type


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: init
    !> @brief Initialize basis set transformation matrix.
    !> @details
    !! The transformation matrix is generated from three separate transformations: from the input
    !! format to the internal format, (if needed) between internal Cartesian and spherical formats
    !! and finally from internal format to output format.
    !----------------------------------------------------------------------------------------------
    subroutine init(self, bset, source_format, target_format)
        use blas95, only : gemm
        use matrix_mod, only : mat_ge_mmm
        class(basis_transform), intent(inout) :: self
        type(basis_set), intent(in) :: bset !< Basis set for which the transformation is generated
        character(len=*), intent(in) :: source_format !< Format from which to convert
        character(len=*), intent(in) :: target_format !< Format to which to convert
        real(dp), allocatable :: wrk(:, :)

        if (self%initialized) then
            deallocate(self%sphe_cart_trans)
            deallocate(self%trans)
            self%initialized = .false.
        end if

        self%source_format = source_format
        self%target_format = target_format
        call transform_to_internal_format(bset, source_format, self%source_sphe, self%source_trans)
        call transform_to_internal_format(bset, target_format, self%target_sphe, self%target_trans)
        self%n1 = size(self%source_trans, 1)
        self%n2 = size(self%target_trans, 1)
        self%sphe_cart_trans = transform_sphe_cart(bset)

        allocate(self%trans(self%n1, self%n2))
        if (self%source_sphe .and. (.not. self%target_sphe)) then
           !allocate(wrk, source=transpose(self%sphe_cart_trans)) ! Major bug with gfortran v7.3
            allocate(wrk(1:self%n1, 1:self%n2))
            wrk = transpose(self%sphe_cart_trans)
            call mat_ge_mmm(self%source_trans, wrk, self%target_trans, self%trans, transc='T')
        else if ((.not. self%source_sphe) .and. self%target_sphe) then
            allocate(wrk, source=self%sphe_cart_trans)
            call mat_ge_mmm(self%source_trans, wrk, self%target_trans, self%trans, transc='T')
        else
            call gemm(self%source_trans, self%target_trans, self%trans, transb='T')
        end if

        self%initialized = .true.
    end subroutine init


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: transform_array
    !> @brief Transform 2D array where the basis coefficients are along dimension bf_dim.
    !----------------------------------------------------------------------------------------------
    subroutine transform_array(self, array, bf_dim)
        use blas95, only : gemm
        class(basis_transform), intent(inout) :: self
        real(dp), allocatable, intent(inout) :: array(:, :) !< Array to transform
        integer, intent(in) :: bf_dim !< Dimension along which to transform
        real(dp), allocatable :: wrk(:, :)

        if (bf_dim == 1) then
            allocate(wrk(self%n2, size(array, 2)))
            call gemm(self%trans, array, wrk, transa='T')
        else
            allocate(wrk(size(array, 1), self%n2))
            call gemm(array, self%trans, wrk)
        end if

        deallocate(array)
        allocate(array, source=wrk)
    end subroutine transform_array


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: transform_arrays
    !> @brief Transform 3D array where the basis coefficients are along dimension bf_dim.
    !----------------------------------------------------------------------------------------------
    subroutine transform_arrays(self, array, bf_dim)
        use blas95, only : gemm
        class(basis_transform), intent(inout) :: self
        real(dp), allocatable, intent(inout) :: array(:, :, :) !< Array to transform
        integer, intent(in) :: bf_dim !< Dimension along which to transform
        real(dp), allocatable :: wrk(:, :, :)
        integer :: i, lb(3), ub(3)

        lb = lbound(array)
        ub = ubound(array)
        if (bf_dim == 1) then
            allocate(wrk(self%n2, lb(2):ub(2), lb(3):ub(3)))
            do i = lbound(array, 3), ubound(array, 3)
                call gemm(self%trans, array(:, :, i), wrk(:, :, i), transa='T')
            end do
        else
            write(stderr, *) 'Error. basis_transform, transform_arrays'
            write(stderr, *) '    bf_dim /= 1 not implemented.'
            stop
        end if

        deallocate(array)
        allocate(array, source=wrk)
    end subroutine transform_arrays


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: transform_to_internal_format
    !> @brief Generate matrix for converting basis functions to internal format.
    !----------------------------------------------------------------------------------------------
    subroutine transform_to_internal_format(bset, basis_format, sphe, mat)
        use molden_definitions_mod, only : molden_sphe_order, molden_cart_order
        use turbomole_definitions_mod, only : turbomole_sphe_order, turbomole_sphe_phase
        use molpro_output_read_mod, only : molpro_sphe_order
        type(basis_set), intent(in) :: bset !<
        character(len=*), intent(in) :: basis_format
        logical, intent(out) :: sphe
        real(dp), allocatable, intent(out) :: mat(:, :)

        select case(basis_format)
        case('turbomole')
            sphe = .true.
            mat = transform_within_subshells(bset, sphe, turbomole_sphe_order, turbomole_sphe_phase)
        case('molden_sphe')
            sphe = .true.
            mat = transform_within_subshells(bset, sphe, molden_sphe_order)
        case('molden_cart')
            sphe = .false.
            mat = transform_within_subshells(bset, sphe, molden_cart_order)
        case('molpro')
            sphe = .true.
            mat = transform_within_subshells(bset, sphe, molpro_sphe_order)
        case('internal')
            sphe = .false.
            mat = transform_within_subshells(bset, sphe, amp%cart)
        case default
            write(stderr, *) 'Error in basis_transform_mod, transform_to_internal_format.'
            write(stderr, *) '    Format not implemented: '//trim(adjustl(basis_format))
            stop
        end select
    end subroutine transform_to_internal_format


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: transform_sphe_cart
    !> @brief Generate matrix for converting between cartesian/spherical basis functions.
    !----------------------------------------------------------------------------------------------
    function transform_sphe_cart(bset) result(wrk)
        type(basis_set), intent(in) :: bset !< Basis set for which the transformation is generated
        real(dp), allocatable :: wrk(:, :)
        integer :: i, j, bi, l, s, c

        allocate(wrk(bset%n_bf_cart, bset%n_bf_sphe), source=num0)
        s = 0
        c = 0
        do i = 1, bset%n_center
            bi = bset%center_i_bs(i)
            do j = 1, bset%bs(bi)%n_subshell
                l = bset%bs(bi)%cg(j)%l
                wrk(c+1:c+amp%n_cart(l), s+1:s+amp%n_sphe(l)) = amp%sphe_cart_trans(l)%c
                s = s + amp%n_sphe(l)
                c = c + amp%n_cart(l)
            end do
        end do
    end function transform_sphe_cart


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: transform_within_subshells
    !> @brief Generate matrix for reordering basis functions within subshells.
    !----------------------------------------------------------------------------------------------
    function transform_within_subshells(bset, sphe, order, phase) result(trans)
        type(basis_set), intent(in) :: bset !< Basis set for which the transformation is generated
        logical, intent(in) :: sphe !< Spherical/Cartesian basis set
        procedure(subshell_numbers) :: order !< Order of bfs within subshells for given l
        procedure(subshell_phases), optional :: phase !< Phase of bfs within subshells for given l
        real(dp), allocatable :: trans(:, :)
        integer :: i, bi, j, l, n(0:amp%max_l), c


        if (sphe) then
            allocate(trans(bset%n_bf_sphe, bset%n_bf_sphe), source=0.0_dp)
            n = amp%n_sphe
        else
            allocate(trans(bset%n_bf_cart, bset%n_bf_cart), source=0.0_dp)
            n = amp%n_cart
        end if


        c = 0
        do i = 1, bset%n_center
            bi = bset%center_i_bs(i)
            do j = 1, bset%bs(bi)%n_subshell
                l = bset%bs(bi)%cg(j)%l
                if (present(phase)) then
                    trans(c+1:c+n(l), c+1:c+n(l)) = amp%subshell_transform(l, order(l), phase(l))
                else
                    trans(c+1:c+n(l), c+1:c+n(l)) = amp%subshell_transform(l, order(l))
                end if
                c = c + n(l)
            end do
        end do
    end function transform_within_subshells


end module basis_transform_mod
