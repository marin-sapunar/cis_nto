!----------------------------------------------------------------------------------------------
! MODULE: ao_overlap_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date March, 2019
!
!> @brief Subroutine for calculating overlap matrix between two sets of atomic orbitals.
!----------------------------------------------------------------------------------------------
module ao_overlap_mod
    use global_defs
    implicit none
    private


    public ao_overlap


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ao_overlap
    !> @brief Calculate overlap matrix between two sets of atomic orbitals.
    !----------------------------------------------------------------------------------------------
    subroutine ao_overlap(path1, format1, path2, format2, s_ao, center_atoms, center_pairs)
        use basis_set_mod
        use one_el_op_mod
        use read_all_mod
        character(len=*), intent(in) :: path1 !< Path for bra input
        character(len=*), intent(in) :: format1 !< Bra input format
        character(len=*), intent(in) :: path2 !< Path for ket input
        character(len=*), intent(in) :: format2 !< Ket input format
        real(dp), allocatable, intent(out) :: s_ao(:, :) !< Atomic orbital overlaps
        logical, intent(in), optional :: center_atoms !< Translate same atom to center.
        logical, intent(in), optional :: center_pairs !< Translate pairs of atoms to center.
        real(dp), allocatable :: geom1(:) !< Geometry 1
        real(dp), allocatable :: geom2(:) !< Geometry 2
        type(basis_set) :: bs1 !< Atomic orbitals 1
        type(basis_set) :: bs2 !< Atomic orbitals 2
        logical :: opt_atoms
        logical :: opt_pairs

        opt_atoms = .false.
        opt_pairs = .false.
        if (present(center_atoms)) opt_atoms = center_atoms
        if (present(center_pairs)) opt_pairs = center_pairs

        call read_geom(format1, path1, geom1)
        call read_geom(format2, path2, geom2)
        call read_basis(format1, path1, bs1)
        call read_basis(format2, path2, bs2)
        if (print_level >= 2) then
            write(stdout, *)
            write(stdout, '(1x,a)') 'Computing AO overlaps...'
        end if
        call one_el_op(geom1, geom2, bs1, bs2, [0,0,0], bs1%sphe_mo, bs2%sphe_mo, .true., .true.,  &
        &              s_ao, opt_atoms, opt_pairs)
    end subroutine ao_overlap


end module ao_overlap_mod
