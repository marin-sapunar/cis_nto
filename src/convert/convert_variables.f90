!--------------------------------------------------------------------------------------------------
! MODULE: convert_variables
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date April, 2019
!
!> @brief Hold global variables for MO convert program
!--------------------------------------------------------------------------------------------------
module convert_variables
    use global_defs
    use basis_set_mod, only : basis_set
    use basis_transform_mod, only : basis_transform
    use molecular_orbitals_mod, only : molecular_orbitals
   !use occupation_mod
    implicit none


    ! Options
    character(len=:), allocatable :: input_path
    character(len=:), allocatable :: input_format
    character(len=:), allocatable :: output_path
    character(len=:), allocatable :: output_format
    ! System
    real(dp), allocatable :: geom(:) !< Geometry
    integer, allocatable :: atnum(:) !< Atom numbers
    character(len=2), allocatable :: atsym(:) !< Atom symbols
    type(basis_set) :: bs !< Basis set
    type(molecular_orbitals) :: mos !< Molecular orbitals
    type(basis_transform) :: bs_trans


end module convert_variables
