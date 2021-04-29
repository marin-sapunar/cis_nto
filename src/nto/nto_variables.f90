!--------------------------------------------------------------------------------------------------
! MODULE: nto_variables
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date April, 2019
!
!> @brief Hold global variables for nto program.
!--------------------------------------------------------------------------------------------------
module nto_variables
    use global_defs
    use basis_set_mod
    use molecular_orbitals_mod
    use occupation_mod
    implicit none


    ! Options
    character(len=:), allocatable :: path
    character(len=:), allocatable :: input_format
    character(len=:), allocatable :: output_format
    character(len=:), allocatable :: outfile_nto
    real(dp) :: wf_threshold
    integer :: truncate_nex
    logical :: norm_states
    logical :: orth_states
    ! System
    real(dp), allocatable :: geom(:) !< Geometry
    integer, allocatable :: atnum(:) !< Atom numbers
    character(len=2), allocatable :: atsym(:) !< Atom symbols
    type(basis_set) :: bs !< Basis set
    type(molecular_orbitals) :: mos !< Molecular orbitals
    integer :: rhf !< Restricted (1) or unrestricted (2) calculation
    type(occupation_numbers) :: on !< Occupation numbers
    logical, allocatable :: occ(:, :) !< Occupied MO mask
    logical, allocatable :: act(:, :) !< Active MO mask from el. structure calculation
    real(dp), allocatable :: cisa(:, :, :) !< CIS matrix alpha
    real(dp), allocatable :: cisb(:, :, :) !< CIS matrix beta
    ! Results
    real(dp), allocatable :: nto_a(:, :, :) !< NTOs alpha.
    real(dp), allocatable :: nto_b(:, :, :) !< NTOs beta.
    real(dp), allocatable :: nto_c_a(:, :) !< NTO alpha coefficients.
    real(dp), allocatable :: nto_c_b(:, :) !< NTO beta coefficients.
    integer, allocatable :: na_a(:) !< Number of active alpha NTO orbitals.
    integer, allocatable :: na_b(:) !< Number of active beta NTO orbitals.
    type(molecular_orbitals) :: nto_mos !< NTOs for printing.
    ! Help
    integer, external :: omp_get_max_threads
    real(dp), external :: omp_get_wtime
    real(dp) :: time00, time0
    real(dp) :: time_in, time_ao, time_mo, time_wf, time_tot


end module nto_variables
