!--------------------------------------------------------------------------------------------------
! MODULE: molecular_orbitals_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date May, 2019
!
!> @brief Contains the molecular_orbitals type and accompanying subroutines.
!--------------------------------------------------------------------------------------------------
module molecular_orbitals_mod
    use global_defs
    use cgto_mod
    use ang_mom_defs, only : amp
    implicit none
 
    private
    public :: molecular_orbitals


    !----------------------------------------------------------------------------------------------
    ! TYPE: basis_set
    !> @brief Contains all information about the full basis set.
    !----------------------------------------------------------------------------------------------
    type molecular_orbitals
        integer :: n_mo_a = 0 !< Number of alpha MOs
        integer :: n_mo_b = 0 !< Number of alpha MOs
        real(dp), allocatable :: ca(:, :) !< MO coefficients alpha
        real(dp), allocatable :: ea(:) !< MO energies alpha
        real(dp), allocatable :: oa(:) !< MO occupation alpha
        integer, allocatable :: na(:) !< MO numbers alpha
        real(dp), allocatable :: cb(:, :) !< MO coefficients beta
        real(dp), allocatable :: eb(:) !< MO energies beta
        real(dp), allocatable :: ob(:) !< MO occupation beta
        integer, allocatable :: nb(:) !< MO numbers beta
    contains
        procedure :: init => init_molecular_orbitals
        procedure :: to_unrestricted
    end type molecular_orbitals


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: init_molecular_orbitals
    !> @brief Initialize instance of molecular_orbitals
    !----------------------------------------------------------------------------------------------
    subroutine init_molecular_orbitals(self, ca, ea, oa, na, cb, eb, ob, nb)
        class(molecular_orbitals), intent(inout) :: self
        real(dp), allocatable, optional :: ca(:, :) !< MO coefficients alpha
        real(dp), allocatable, optional :: ea(:) !< MO energies alpha
        real(dp), allocatable, optional :: oa(:) !< MO occupation alpha
        integer, allocatable, optional :: na(:) !< MO numbers alpha
        real(dp), allocatable, optional :: cb(:, :) !< MO coefficients beta
        real(dp), allocatable, optional :: eb(:) !< MO energies beta
        real(dp), allocatable, optional :: ob(:) !< MO occupation beta
        integer, allocatable, optional :: nb(:) !< MO numbers beta
        integer :: i

        if (present(ca)) then
            if (allocated(ca)) then
                self%n_mo_a = size(ca, 2)
            end if
        end if
        if (present(cb)) then
            if (allocated(cb)) then
                self%n_mo_b = size(cb, 2)
            end if
        end if

        if (self%n_mo_a > 0) then
            ! Copy alpha arrays.
            if (present(ca)) self%ca = ca
            if (present(oa)) self%oa = oa
            if (present(ea)) self%ea = ea
            if (present(na)) self%na = na
            ! Initialize non-present arrays to default values.
            if (.not. allocated(self%oa)) allocate(self%oa(self%n_mo_a), source=0.0_dp)
            if (.not. allocated(self%ea)) allocate(self%ea(self%n_mo_a), source=0.0_dp)
            if (.not. allocated(self%na)) then
                allocate(self%na(self%n_mo_a))
                self%na = [ (i, i = 1, self%n_mo_a ) ]
            end if
        end if
        if (self%n_mo_b > 0) then
            ! Copy beta arrays.
            if (present(cb)) self%cb = cb
            if (present(ob)) self%ob = ob
            if (present(eb)) self%eb = eb
            if (present(nb)) self%nb = nb
            ! Initialize non-present arrays to default values.
            if (.not. allocated(self%ob)) allocate(self%ob(self%n_mo_b), source=0.0_dp)
            if (.not. allocated(self%eb)) allocate(self%eb(self%n_mo_b), source=0.0_dp)
            if (.not. allocated(self%nb)) then
                allocate(self%nb(self%n_mo_b))
                self%nb = [ (i, i = 1, self%n_mo_b ) ]
            end if
        end if
    end subroutine init_molecular_orbitals


    subroutine to_unrestricted(self)
        class(molecular_orbitals), intent(inout) :: self

        if (self%n_mo_b == self%n_mo_a) return
        self%n_mo_b = self%n_mo_a
        self%cb = self%ca
        self%ob = self%oa
        self%eb = self%ea
        self%nb = self%na
    end subroutine to_unrestricted


end module molecular_orbitals_mod
