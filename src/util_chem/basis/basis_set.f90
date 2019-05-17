!--------------------------------------------------------------------------------------------------
! MODULE: basis_set_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date June, 2017
!
!> @brief Contains the basis_set type and accompanying subroutines.
!--------------------------------------------------------------------------------------------------
module basis_set_mod
    use global_defs
    use cgto_mod
    use ang_mom_defs, only : amp
    implicit none
 
    private
    public :: basis_set_single
    public :: basis_set


    !----------------------------------------------------------------------------------------------
    ! TYPE: basis_set_single
    !> @brief Contains all information about the basis set for a single atom.
    !----------------------------------------------------------------------------------------------
    type basis_set_single
        character(len=:), allocatable :: key !< Name of the basis set for the atom.
        integer :: n_subshell = 0 !< Number of basis functions for the atom.
        type(cgto_subshell), allocatable :: cg(:) !< Contracted Gaussian functions.
        integer :: n_cart !< Number of cartesian basis functions.
        integer :: n_sphe !< Number of spherical basis functions.
    contains
        procedure :: init => init_basis_set_single
    end type basis_set_single


    !----------------------------------------------------------------------------------------------
    ! TYPE: basis_set
    !> @brief Contains all information about the full basis set.
    !----------------------------------------------------------------------------------------------
    type basis_set
        integer :: n_bs = 0 !< Number of unique basis sets.
        type(basis_set_single), allocatable :: bs(:) !< Basis sets.
        integer :: n_center = 0 !< Number of centers(atoms) for basis functions.
        integer, allocatable :: center_i_bs(:) !< Index of basis set corresponding to each center.
        integer :: n_bf_cart = 0 !< Total number of cartesian basis functions.
        integer :: n_bf_sphe = 0 !< Total number of spherical basis functions.
        character(len=:), allocatable :: source_format !< Format from which the basis is read.
        logical :: sphe_mo
    contains
        procedure :: check_init => check_init_basis_set
    end type basis_set


contains

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: init_basis_set_single
    !> @brief Initialize instance of basis_set_single.
    !----------------------------------------------------------------------------------------------
    subroutine init_basis_set_single(self, nfunc, orbtype, nprim, zeta, beta, sort_l)
        class(basis_set_single), intent(inout) :: self
        integer, intent(in) :: nfunc
        character(len=1), intent(in) :: orbtype(:)
        integer, intent(in) :: nprim(:)
        real(dp), intent(in) :: zeta(:, :)
        real(dp), intent(in) :: beta(:, :)
        logical, intent(in), optional :: sort_l
        integer :: i, cfunc, l
        logical :: sort

        sort = .false.
        if (present(sort_l)) sort = sort_l

        cfunc = 0
        self%n_subshell = nfunc
        if (allocated(self%cg)) deallocate(self%cg)
        allocate(self%cg(nfunc))
        do l = 0, amp%max_l
            do i = 1, nfunc
                if ((amp%l(orbtype(i)) == l) .or. (.not. sort)) then
                    cfunc = cfunc + 1
                    allocate(self%cg(cfunc)%z(nprim(i)))
                    allocate(self%cg(cfunc)%b(nprim(i)))
                    self%cg(cfunc)%typ = orbtype(i)
                    self%cg(cfunc)%l = amp%l(orbtype(i))
                    self%cg(cfunc)%n_prim = nprim(i)
                    self%cg(cfunc)%z = zeta(i, 1:nprim(i))
                    self%cg(cfunc)%b = beta(i, 1:nprim(i))
                    call self%cg(cfunc)%norm_b()
                    self%n_sphe = self%n_sphe + amp%n_sphe(self%cg(cfunc)%l)
                    self%n_cart = self%n_cart + amp%n_cart(self%cg(cfunc)%l)
                end if
            end do
            if (cfunc == nfunc) exit
            if (l == amp%max_l) then
                write(stderr, *)
                write(stderr, *) 'Error in basis_set module, init_basis_set_single subroutine.'
                write(stderr, *) '    Basis functions with l>lmax present.'
               stop
            end if
        end do
    end subroutine init_basis_set_single


    subroutine check_init_basis_set(self)
        class(basis_set), intent(inout) :: self
        integer :: i, bi

        if (.not. allocated(self%bs)) then
            write(stderr, *) 'Error in check_init_basis_set. Basis set not defined.'
            stop
        end if
        if (self%n_bs == 0) self%n_bs = size(self%bs)
        if (.not. allocated(self%center_i_bs)) then
            allocate(self%center_i_bs(1:self%n_bs), source=[ (i, i = 1, self%n_bs) ])
        end if
        if (self%n_center == 0) self%n_center = size(self%center_i_bs)
        if (self%n_bf_cart == 0) then
            do i = 1, self%n_center
                bi = self%center_i_bs(i)
                self%n_bf_cart = self%n_bf_cart + self%bs(bi)%n_cart
            end do
        end if
        if (self%n_bf_sphe == 0) then
            do i = 1, self%n_center
                bi = self%center_i_bs(i)
                self%n_bf_sphe = self%n_bf_sphe + self%bs(bi)%n_sphe
            end do
        end if
    end subroutine check_init_basis_set


end module basis_set_mod
