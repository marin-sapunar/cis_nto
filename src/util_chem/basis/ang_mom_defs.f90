!----------------------------------------------------------------------------------------------
! MODULE: ang_mom_defs
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date March, 2019
!
!> @brief Helper functions for gaussian basis set definitions and transformations.
!----------------------------------------------------------------------------------------------
module ang_mom_defs
    use global_defs
    implicit none


    private
    public :: ang_mom_parameters
    public :: amp
    public :: sphe_nums
    public :: cart_nums
    public :: matrix_sphe_to_cart_gaussian
    public :: subshell_numbers
    public :: subshell_phases

    character(len=1), parameter :: l_lett(0:11) = & ! Assumes max_l always < 11.
    &    ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n']

    interface
        function subshell_numbers(ll) result(nums)
            integer, intent(in) :: ll
            integer, allocatable :: nums(:, :)
        end function

        function subshell_phases(ll) result(phase)
            use global_defs, only : dp
            integer, intent(in) :: ll
            real(dp), allocatable :: phase(:)
        end function
    end interface



    !----------------------------------------------------------------------------------------------
    ! TYPE: ang_mom_parameters
    !> @brief Holds parameters for cartesian/spherical gaussian functions.
    !
    !> @details
    !! Parameters up to a given maximum angular momentum max_l are saved at the start by calling
    !! the %init(max_l) procedure.
    !----------------------------------------------------------------------------------------------
    type ang_mom_parameters
        logical :: initialized = .false. !< Instance initialized?
        integer :: max_l !< Chosen maximum l.
        integer, allocatable :: n_cart(:) !< Number of cartesian functions.
        integer, allocatable :: n_sphe(:) !< Number of spherical functions.
        character(len=1), allocatable :: l_lett(:) !< Letter for subshells with given l.
        type(rmat), allocatable :: sphe_cart_trans(:) !< Transformation matrices from sphe to cart.
        procedure(subshell_numbers), pointer, nopass :: sphe !< Order of l, m  combinations.
        procedure(subshell_numbers), pointer, nopass :: cart !< Order of lx, ly, lz combinations.
    contains
        procedure :: init => ang_mom_parameters_init
        procedure :: subshell_transform
        procedure, nopass :: l => subshell_l
    end type ang_mom_parameters


    type(ang_mom_parameters) :: amp !< Instance of ang_mom_parameters to hold default values.


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ang_mom_parameterers_init
    !> @brief Initialize ang_mom_parameterers instance.
    !----------------------------------------------------------------------------------------------
    subroutine ang_mom_parameters_init(self, max_l, sphe_nums_generator, cart_nums_generator)
        class(ang_mom_parameters), intent(inout) :: self
        procedure(subshell_numbers), optional :: sphe_nums_generator
        procedure(subshell_numbers), optional :: cart_nums_generator
        integer, intent(in) :: max_l
        integer :: l

        ! Deallocate if previously initialized.
        if (self%initialized) then
            deallocate(self%n_cart)
            deallocate(self%n_sphe)
            deallocate(self%l_lett)
            deallocate(self%sphe_cart_trans)
        end if

        ! Allocate arrays.
        self%max_l = max_l
        allocate(self%n_cart(0:max_l))
        allocate(self%n_sphe(0:max_l))
        allocate(self%l_lett(0:max_l))
        allocate(self%sphe_cart_trans(0:max_l))
        self%sphe => sphe_nums
        self%cart => cart_nums
        if (present(sphe_nums_generator)) self%sphe => sphe_nums_generator
        if (present(cart_nums_generator)) self%cart => cart_nums_generator

        ! Initialize arrays.
        self%n_cart = [( (l+1)*(l+2)/2, l = 0, max_l )]
        self%n_sphe = [( 2*l+1, l = 0, max_l )]
        self%l_lett = l_lett(0:max_l)
        do l = 0, max_l
            allocate(self%sphe_cart_trans(l)%c(self%n_sphe(l), self%n_cart(l)))
            self%sphe_cart_trans(l)%c = matrix_sphe_to_cart_gaussian(self%sphe(l), self%cart(l))
        end do
        self%initialized = .true.
    end subroutine ang_mom_parameters_init


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: subshell_l
    !> @brief Convert subshell letter (s/p/d...) to angular momentum number.
    !----------------------------------------------------------------------------------------------
    elemental function subshell_l(s_type) result(i)
        use string_mod, only : tolower
        character(len=1), intent(in) :: s_type
        integer :: i

        do i = 0, size(l_lett)
            if (l_lett(i) == tolower(s_type)) return
        end do
    end function subshell_l


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: subshell_transform
    !> @brief Return a matrix for reordering subshell functions
    !> @details
    !! The functions can also be multiplied by a scalar phase factor.
    !----------------------------------------------------------------------------------------------
    function subshell_transform(self, l, reorder, phase) result(trans_mat)
        class(ang_mom_parameters) :: self
        integer, intent(in) :: l
        integer, intent(in) :: reorder(:, :)
        real(dp), intent(in), optional :: phase(:)
        real(dp), allocatable :: trans_mat(:, :)
        integer, allocatable :: wrk(:, :)
        real(dp), allocatable :: wrk_phase(:)
        integer :: i, j, n

        if (size(reorder, 1) == 2) then
            n = self%n_sphe(l)
            allocate(wrk, source=self%sphe(l))
        else
            n = self%n_cart(l)
            allocate(wrk, source=self%cart(l))
        end if
        allocate(trans_mat(n, n), source=0.0_dp)
        allocate(wrk_phase(n), source=1.0_dp)
        if (present(phase)) wrk_phase = wrk_phase*phase

        do i = 1, n
            do j = 1, n
                if (any(wrk(:, i) /= reorder(:, j))) cycle
                trans_mat(j, i) = wrk_phase(j)
            end do
        end do
    end function subshell_transform


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: sphe_nums
    !> @brief Default generator function for order of spherical gaussian ang. mom. numbers.
    !----------------------------------------------------------------------------------------------
    function sphe_nums(ll) result(sphe)
        integer, intent(in) :: ll
        integer, allocatable :: sphe(:, :)
        integer :: i

        allocate(sphe(2, 2*ll+1))
        sphe(1, :) = ll
        select case(ll)
        case(0)
            sphe(2, :) = [0]
        case(1)
            sphe(2, :) = [1, -1, 0]
        case default
            sphe(2, 1) = 0
            do i = 1, ll
                sphe(2, 2*i) = i
                sphe(2, 2*i+1) = -i
            end do
        end select
    end function sphe_nums


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: cart_nums
    !> @brief Default generator function for order of cartesian gaussian ang. mom. numbers.
    !----------------------------------------------------------------------------------------------
    function cart_nums(ll) result(cart)
        integer, intent(in) :: ll
        integer, allocatable :: cart(:, :)
        integer :: i, l1, l2, l3

        allocate(cart(3, (ll+1)*(ll+2)/2))
        i = 0
        do l1 = ll, ceiling(ll/3.0), -1
            do l2 = min(ll-l1,l1), ceiling((ll-l1)/2.0), -1
                l3 = ll - l1 - l2
                i = i + 1
                cart(:, i) = [l1, l2, l3]
                if (l2 /= l3) then
                    i = i + 1
                    cart(:, i) = [l1, l3, l2]
                end if
                if (l1 /= l2) then
                    i = i + 1
                    cart(:, i) = [l2, l1, l3]
                end if
                if ((l1 /= l2) .and. (l1 /= l3) .and. (l2 /= l3)) then
                    i = i + 1
                    cart(:, i) = [l3, l1, l2]
                    i = i + 1
                    cart(:, i) = [l2, l3, l1]
                end if
                if (l1 /= l3) then
                    i = i + 1
                    cart(:, i) = [l3, l2, l1]
                end if
            end do
        end do
    end function cart_nums


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: matrix_sphe_to_cart_gaussian
    !> @brief Return a transformation matrix for converting from spherical to cartesian functions.
    !----------------------------------------------------------------------------------------------
    function matrix_sphe_to_cart_gaussian(sphe, cart) result(trans)
        integer, intent(in) :: sphe(:, :)
        integer, intent(in) :: cart(:, :)
        real(dp) :: trans(size(cart, 2), size(sphe, 2))
        integer :: i, j

        do i = 1, size(sphe, 2)
            do j = 1, size(cart, 2)
                trans(j, i) = sphe_to_cart_gaussian(sphe(1, i), sphe(2, i),                        &
                &                                   cart(1, j), cart(2, j), cart(3, j))
            end do
        end do
    end function matrix_sphe_to_cart_gaussian


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: sphe_to_cart_gaussian
    !> @brief Return a single spherical to cartesian gaussian transformation coefficient.
    !----------------------------------------------------------------------------------------------
    function sphe_to_cart_gaussian(ll, mm, lx, ly, lz) result(coeff)
        use math_mod
        integer, intent(in) :: ll
        integer, intent(in) :: mm
        integer, intent(in) :: lx
        integer, intent(in) :: ly
        integer, intent(in) :: lz
        real(dp) :: coeff
        real(dp) :: sq_lx, sq_ly, sq_lz, sq_ll, sq_lm, sq_tot
        real(dp) :: frac, sum1, sum2
        integer :: i, k, am, jj, ediff
        integer(j15) :: csum

        am = abs(mm)
        if (ll /= lx + ly + lz) then
            coeff = 0.0_dp
            return
        else if (abs(mod(lx+ly-am, 2)) == 1) then
            coeff = 0.0_dp
            return
        end if
        if (mm >= 0) then
            if (abs(mod(am-lx, 2)) == 1) then
                coeff = 0.0_dp
                return
            end if
            ediff = am - lx
        else if (mm < 0) then
            if (abs(mod(am-lx, 2)) == 0) then
                coeff = 0.0_dp
                return
            end if
            ediff = am - lx + 1
        end if
        jj = (lx + ly - am) / 2
        sq_lx = real(factorial(2*lx), kind=dp) / factorial(lx)
        sq_ly = real(factorial(2*ly), kind=dp) / factorial(ly)
        sq_lz = real(factorial(2*lz), kind=dp) / factorial(lz)
        sq_ll = real(factorial(ll), kind=dp) / factorial(2*ll)
        sq_lm = real(factorial(ll-am), kind=dp) / factorial(ll+am)
        sq_tot = sqrt(sq_lx*sq_ly*sq_lz*sq_ll*sq_lm)
        frac = 1.0_dp / 2**ll / factorial(ll)
        sum1 = 0.0_dp
        do i = 0, (ll-am)/2
            csum = binomial(ll, i) * binomial(i, jj)
            csum = csum * (-1)**i * factorial(2*ll - 2*i)
            csum = csum / factorial(ll - am - 2*i)
            sum1 = sum1 + csum
        end do
        sum2 = 0.0_dp
        do k = 0, jj
            csum = binomial(jj, k) * binomial(am, lx - 2*k)
            csum = csum * (-1)**(sign(1, mm) * (ediff + 2*k) / 2)
            sum2 = sum2 + csum
        end do
        if (mm > 0) then
            sum2 = sum2 * sqrt(2.0_dp)
        else if (mm < 0) then
            sum2 = - sum2 * sqrt(2.0_dp)
        end if
        coeff = sq_tot * frac * sum1 * sum2
    end function sphe_to_cart_gaussian


end module ang_mom_defs
