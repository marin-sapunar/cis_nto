!----------------------------------------------------------------------------------------------
! MODULE: check_dimensions_overlap_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date January, 2019
!
! DESCRIPTION:
!> @brief Ensure dimensions for two sets of states are compatible before overlap calculations
!----------------------------------------------------------------------------------------------
module check_dimensions_overlap_mod
    ! General
    use global_defs
    ! Chem
    use occupation_mod
    implicit none


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: check_occ_orb
    !
    !> @brief Check match for number of (active) occupied orbitals for overlap calculation
    !> @details
    !! Number of occupied orbitals of each spin must be equal for both sets of states. If number
    !! of active occupied orbitals is not equal, frozen orbitals from both calculations are
    !! converted to active orbitals and CI arrays are expanded with excitations from the new active
    !! orbitals (with all coefficients equal to zero).
    !----------------------------------------------------------------------------------------------
    subroutine check_occ_orb(rhf1, rhf2, on1, on2, occ1, occ2, act1, act2, wfa1, wfa2, wfb1, wfb2)
        integer, intent(in) :: rhf1 !< Restricted (1) or unrestricted (2) calculation 1.
        integer, intent(in) :: rhf2 !< Restricted (1) or unrestricted (2) calculation 2.
        type(occupation_numbers), intent(inout) :: on1 !< Occupation numbers 1.
        type(occupation_numbers), intent(inout) :: on2 !< Occupation numbers 2.
        logical, intent(in) :: occ1(:, :) !< Occupied MO mask 1.
        logical, intent(in) :: occ2(:, :) !< Occupied MO mask 2.
        logical, intent(inout) :: act1(:, :) !< Active MO mask 1.
        logical, intent(inout) :: act2(:, :) !< Active MO mask 2.
        real(dp), allocatable, intent(inout) :: wfa1(:, :) !< CI vector alpha single 1.
        real(dp), allocatable, intent(inout) :: wfa2(:, :) !< CI vector alpha single 2.
        real(dp), allocatable, intent(inout) :: wfb1(:, :) !< CI vector beta single 1.
        real(dp), allocatable, intent(inout) :: wfb2(:, :) !< CI vector beta single 2.

        if ((on1%o(1) /= on2%o(1)) .or. (on1%o(rhf1) /= on2%o(rhf2))) then
            write(stderr, *)
            write(stderr, '(1x,a,a,a)') 'Error. Mismatch in number of occupied orbitals.'
            stop
        end if

        if ((on1%ao(1) /= on2%ao(1)) .or. (on1%ao(rhf1) /= on2%ao(rhf2))) then
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(1x,a,a,a)') 'Warning. Mismatch in number of active occupied orbitals.'
                write(stdout, '(3x,i0,a,i0,a)') on1%ao(1), ' alpha and ', on1%ao(rhf1), &
                &                               ' beta active occupied bra orbitals.'
                write(stdout, '(3x,i0,a,i0,a)') on2%ao(1), ' alpha and ', on2%ao(rhf2), &
                &                               ' beta active occupied ket orbitals.'
                write(stdout, '(3x,a)') 'Using previously frozen orbitals as active orbitals.'
            end if
            call wf_prepend_occ(on1%fo(1), on1%ao(1), on1%av(1), wfa1)
            call wf_prepend_occ(on2%fo(1), on2%ao(1), on2%av(1), wfa2)
            if (rhf1 == 2) call wf_prepend_occ(on1%fo(2), on1%ao(2), on1%av(2), wfb1)
            if (rhf2 == 2) call wf_prepend_occ(on2%fo(2), on2%ao(2), on2%av(2), wfb2)
            where(occ1) act1 = .true.
            where(occ2) act2 = .true.
            on1%ao(1) = on1%o(1)
            on2%ao(1) = on2%o(1)
            on1%ao(rhf1) = on1%o(rhf1)
            on2%ao(rhf2) = on2%o(rhf2)
        end if
    end subroutine check_occ_orb


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: check_rhf
    !
    !> @brief Check for overlap calculations between restricted and unrestricted sets of states.
    !> @details
    !! In case of mixed restricted/unrestricted calculation, the restricted MO and CI coefficients
    !! are used for both alpha and beta spin, and the CI coefficients are renormalized.
    !----------------------------------------------------------------------------------------------
    subroutine check_rhf(rhf1, rhf2, moa1, moa2, mob1, mob2, wfa1, wfa2, wfb1, wfb2)
        integer, intent(in) :: rhf1 !< Restricted (1) or unrestricted (2) calculation 1.
        integer, intent(in) :: rhf2 !< Restricted (1) or unrestricted (2) calculation 2.
        real(dp), allocatable :: moa1(:, :) !< Molecular orbital coefficients alpha 1.
        real(dp), allocatable :: moa2(:, :) !< Molecular orbital coefficients alpha 2.
        real(dp), allocatable :: mob1(:, :) !< Molecular orbital coefficients beta 1.
        real(dp), allocatable :: mob2(:, :) !< Molecular orbital coefficients beta 2.
        real(dp), allocatable, intent(inout) :: wfa1(:, :) !< CI vector alpha single 1.
        real(dp), allocatable, intent(inout) :: wfa2(:, :) !< CI vector alpha single 2.
        real(dp), allocatable, intent(inout) :: wfb1(:, :) !< CI vector beta single 1.
        real(dp), allocatable, intent(inout) :: wfb2(:, :) !< CI vector beta single 2.
        integer :: rhf

        rhf = max(rhf1, rhf2)
        if (rhf1 /= rhf) then
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(1x,a)') 'Unrestricted bra and restricted ket calculation...'
                write(stdout, '(5x,a)') 'Renormalizing bra CI coefficients...'
            end if
            wfa1 = wfa1 / sqrt(2.0_dp)
            allocate(mob1, source=moa1)
            allocate(wfb1, source=wfa1)
            wfb1 = -wfb1
        else if (rhf2 /= rhf) then
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(1x,a)') 'Restricted bra and unrestricted ket calculation...'
                write(stdout, '(5x,a)') 'Renormalizing ket CI coefficients...'
            end if
            wfa2 = wfa2 / sqrt(2.0_dp)
            allocate(mob2, source=moa2)
            allocate(wfb2, source=wfa2)
            wfb2 = -wfb2
        end if
    end subroutine check_rhf


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: wf_prepend_occ
    !
    !> @brief Add occupied orbitals to start of a wf array
    !----------------------------------------------------------------------------------------------
    subroutine wf_prepend_occ(nadd, no, nv, wf)
        integer, intent(in) :: nadd
        integer, intent(in) :: no
        integer, intent(in) :: nv
        real(dp), allocatable, intent(inout) :: wf(:, :)
        real(dp), allocatable :: tmp_wf(:, :)

        allocate(tmp_wf((nadd+no)*nv, size(wf,2)), source=0.0_dp)
        tmp_wf(nadd*nv+1:, :) = wf
        deallocate(wf)
        allocate(wf, source=tmp_wf)
    end subroutine wf_prepend_occ



end module check_dimensions_overlap_mod
