!--------------------------------------------------------------------------------------------------
! MODULE: overlap_variables
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @brief Hold global variables for overlap program.
!--------------------------------------------------------------------------------------------------
module overlap_variables
    use global_defs
    use basis_set_mod
    use molecular_orbitals_mod
    use occupation_mod
    implicit none


    ! Options
    character(len=:), allocatable :: path1
    character(len=:), allocatable :: path2
    character(len=:), allocatable :: input_format_1
    character(len=:), allocatable :: input_format_2
    logical :: ao_stop
    logical :: mo_stop
    logical :: dyson_c
    character(len=4) :: cis_algorithm
    real(dp) :: wf_threshold
    integer :: truncate_nex
    logical :: norm_states
    logical :: orth_states
    logical :: orth_overlap
    logical :: match_phase
    logical :: center_atoms
    logical :: center_pairs
    logical :: freeze_mo_norm
    real(dp) :: freeze_mo_norm_t
    character(len=:), allocatable :: outfile_ao
    character(len=:), allocatable :: outfile_mo
    character(len=:), allocatable :: outfile_wf
    character(len=:), allocatable :: prefix_dyson
    ! System
    real(dp), allocatable :: geom1(:) !< Geometry 1
    real(dp), allocatable :: geom2(:) !< Geometry 2
    character(len=2), allocatable :: atom_symbol(:) !< Atom symbols (should be same for both)
    integer, allocatable :: atom_number(:) !< Atom numbers (should be same for both)
    type(basis_set) :: bs1 !< Atomic orbitals 1
    type(basis_set) :: bs2 !< Atomic orbitals 2
    type(molecular_orbitals) :: mos1 !< Molecular orbitals 1
    type(molecular_orbitals) :: mos2 !< Molecular orbitals 2
    integer :: rhf !< Restricted (1) or unrestricted (2) calculation
    integer :: rhf1 !< Restricted (1) or unrestricted (2) calculation 1
    integer :: rhf2 !< Restricted (1) or unrestricted (2) calculation 2
    type(occupation_numbers) :: on1 !< Occupation numbers 1
    type(occupation_numbers) :: on2 !< Occupation numbers 2
    logical, allocatable :: occ1(:, :) !< Occupied MO mask 1
    logical, allocatable :: occ2(:, :) !< Occupied MO mask 2
    logical, allocatable :: act1(:, :) !< Active MO mask from el. structure calculation 1
    logical, allocatable :: act2(:, :) !< Active MO mask from el. structure calculation 2
    real(dp), allocatable :: cisa1(:, :, :) !< CIS matrix alpha 1
    real(dp), allocatable :: cisa2(:, :, :) !< CIS matrix alpha 2
    real(dp), allocatable :: cisb1(:, :, :) !< CIS matrix beta 1
    real(dp), allocatable :: cisb2(:, :, :) !< CIS matrix beta 2
    ! Results
    real(dp), allocatable :: s_ao(:, :) !< Atomic orbital overlaps.
                                        !! Dimensions: (n_ao1, n_ao2)
    real(dp), allocatable :: s_mo_a(:, :) !< Molecular orbital overlaps alpha
                                          !! Dimensions: (n_mo_a1, n_mo_a2)
    real(dp), allocatable :: s_mo_b(:, :) !< Molecular orbital overlaps beta
                                          !! Dimensions: (n_mo_b1, n_mo_b2)
    real(dp), allocatable :: s_wf(:, :) !< Wave function overlaps
                                        !! Dimensions: (0:n_ex_st1,  0:n_ex_st2)
    real(dp), allocatable :: dyson_ao(:, :, :) !< Dyson orbs. in terms of AO coefficients
                                               !! Dimensions: (n_ao1, 0:n_ex_st1,  0:n_ex_st2)
    real(dp), allocatable :: dyson_mo(:, :, :) !< Dyson orbs. in terms of MO coefficients
                                               !! Dimensions: (n_mo1, 0:n_ex_st1,  0:n_ex_st2)
    real(dp), allocatable :: dyson_norm(:, :) !< Norms of Dyson orbitals
                                                 !! Dimensions: (0:n_ex_st1,  0:n_ex_st2)
    ! Help
    integer, external :: omp_get_max_threads
    real(dp), external :: omp_get_wtime
    real(dp) :: time00, time0
    real(dp) :: time_in, time_ao, time_mo, time_wf, time_tot


end module overlap_variables
