!----------------------------------------------------------------------------------------------
! MODULE: mo_overlap_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date March, 2019
!
!> @brief Subroutine for calculating overlap matrix between two sets of molecular orbitals..
!----------------------------------------------------------------------------------------------
module mo_overlap_mod
    use global_defs
    implicit none
    private


    public mo_overlap


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: mo_overlap
    !> @brief Calculate overlap matrix between two sets of molecular orbitals.
    !----------------------------------------------------------------------------------------------
    subroutine mo_overlap(path1, format1, path2, format2, s_ao, s_mo_a, s_mo_b)
        use matrix_mod, only : mat_ge_mmm
        use read_all_mod
        character(len=*), intent(in) :: path1 !< Path for bra input
        character(len=*), intent(in) :: format1 !< Bra input format
        character(len=*), intent(in) :: path2 !< Path for ket input
        character(len=*), intent(in) :: format2 !< Ket input format
        real(dp), intent(in) :: s_ao(:, :) !< AO overlaps
        real(dp), allocatable, intent(out) :: s_mo_a(:, :) !< MO overlaps alpha
        real(dp), allocatable, intent(out) :: s_mo_b(:, :) !< MO overlaps beta
        real(dp), allocatable :: moa1(:, :) !< MO coefficients alpha 1
        real(dp), allocatable :: moa2(:, :) !< MO coefficients alpha 2
        real(dp), allocatable :: mob1(:, :) !< MO coefficients beta 1
        real(dp), allocatable :: mob2(:, :) !< MO coefficients beta 2
        logical :: uhf1, uhf2
        integer :: i

        if (allocated(s_mo_a)) deallocate(s_mo_a)
        if (allocated(s_mo_b)) deallocate(s_mo_b)

        ! Read MOs.
        call read_mo(format1, path1, moa_c=moa1, mob_c=mob1)
        call read_mo(format2, path2, moa_c=moa2, mob_c=mob2)

        ! Check unrestricted/restricted calculations.
        uhf1 = allocated(mob1)
        uhf2 = allocated(mob2)
        if (uhf1 .and. (.not. uhf2)) then
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(1x,a)') 'Unrestricted bra and restricted ket MOs...'
            end if
            allocate(mob2, source=moa2)
        else if (uhf2 .and. (.not. uhf1)) then
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(1x,a)') 'Restricted bra and unrestricted ket MOs...'
            end if
            allocate(mob1, source=moa1)
        end if

        ! Calculate overlaps.
        if (print_level >= 2) then
            write(stdout, *)
            write(stdout, '(1x,a)') 'Computing MO overlaps...'
        end if
        allocate(s_mo_a(size(moa1, 2), size(moa2, 2)))
        call mat_ge_mmm(moa1, s_ao, moa2, s_mo_a, transa='T')
        if (uhf1 .or. uhf2) then
            allocate(s_mo_b(size(mob1, 2), size(mob2, 2)))
            call mat_ge_mmm(mob1, s_ao, mob2, s_mo_b, transa='T')
        end if
    end subroutine mo_overlap


end module mo_overlap_mod
