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
    ! SUBROUTINE: print_calculation_properties
    !
    !> @brief Print dimensions of an electronic structure calculation to stdout.
    !----------------------------------------------------------------------------------------------
    subroutine print_calculation_properties(msg, rhf, n_ex_state, n_sphe, n_cart, on)
        character(len=*), intent(in) :: msg
        integer, intent(in) :: rhf
        integer, intent(in) :: n_ex_state
        integer, intent(in) :: n_sphe
        integer, intent(in) :: n_cart
        type(occupation_numbers), intent(in) :: on

        if (print_level >= 1) then
            write(stdout, *)
            write(stdout, '(1x,a)') msg
            if (rhf == 1) write(stdout, '(5x, a)') 'Restricted calculation.'
            if (rhf == 2) write(stdout, '(5x, a)') 'Unrestricted calculation.'
            write(stdout, '(5x,a,2(1x,i0))') 'Number of excited states:           ', n_ex_state
            write(stdout, '(5x,a,2(1x,i0))') 'Number of basis functions (sphe):   ', n_sphe
            write(stdout, '(5x,a,2(1x,i0))') 'Number of basis functions (cart):   ', n_cart
            write(stdout, '(5x,a,2(1x,i0))') 'Number of orbitals:                 ', on%n
            write(stdout, '(5x,a,2(1x,i0))') 'Number of occupied orbitals:        ', on%o(1:rhf)
            write(stdout, '(5x,a,2(1x,i0))') 'Number of virtual orbitals:         ', on%v(1:rhf)
            write(stdout, '(5x,a,2(1x,i0))') 'Number of active occupied orbitals: ', on%ao(1:rhf)
            write(stdout, '(5x,a,2(1x,i0))') 'Number of active virtual orbitals:  ', on%av(1:rhf)
        end if
    end subroutine print_calculation_properties
        

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
    subroutine check_occ_orb(rhf1, rhf2, on1, on2, occ1, occ2, act1, act2, cisa1, cisa2, cisb1, cisb2)
        integer, intent(in) :: rhf1 !< Restricted (1) or unrestricted (2) calculation 1.
        integer, intent(in) :: rhf2 !< Restricted (1) or unrestricted (2) calculation 2.
        type(occupation_numbers), intent(inout) :: on1 !< Occupation numbers 1.
        type(occupation_numbers), intent(inout) :: on2 !< Occupation numbers 2.
        logical, intent(in) :: occ1(:, :) !< Occupied MO mask 1.
        logical, intent(in) :: occ2(:, :) !< Occupied MO mask 2.
        logical, intent(inout) :: act1(:, :) !< Active MO mask 1.
        logical, intent(inout) :: act2(:, :) !< Active MO mask 2.
        real(dp), allocatable, intent(inout) :: cisa1(:, :, :) !< CIS matrix alpha 1.
        real(dp), allocatable, intent(inout) :: cisa2(:, :, :) !< CIS matrix alpha 2.
        real(dp), allocatable, intent(inout) :: cisb1(:, :, :) !< CIS matrix beta 1.
        real(dp), allocatable, intent(inout) :: cisb2(:, :, :) !< CIS matrix beta 2.


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
            call cis_prepend_occ(on1%fo(1), on1%o(1), on1%av(1), cisa1)
            call cis_prepend_occ(on2%fo(1), on2%o(1), on2%av(1), cisa2)
            if (rhf1 == 2) call cis_prepend_occ(on1%fo(2), on1%o(2), on1%av(2), cisb1)
            if (rhf2 == 2) call cis_prepend_occ(on2%fo(2), on2%o(2), on2%av(2), cisb2)
            where(occ1) act1 = .true.
            where(occ2) act2 = .true.
            on1%ao(1) = on1%o(1)
            on2%ao(1) = on2%o(1)
            on1%ao(rhf1) = on1%o(rhf1)
            on2%ao(rhf2) = on2%o(rhf2)
        end if
    end subroutine check_occ_orb


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: remove_frozen_mo
    !
    !> @brief Remove MOs frozen during electronic structure calculation.
    !> @details
    !! These MOs are removed from the array of MO coefficients and are not used for the remainder
    !! of the calculation.
    !----------------------------------------------------------------------------------------------
    subroutine remove_frozen_mo(rhf1, rhf2, occ1, occ2, act1, act2, moa1, moa2, mob1, mob2)
        integer, intent(in) :: rhf1 !< Restricted (1) or unrestricted (2) calculation 1.
        integer, intent(in) :: rhf2 !< Restricted (1) or unrestricted (2) calculation 2.
        logical, intent(in) :: occ1(:, :) !< Occupied MO mask 1.
        logical, intent(in) :: occ2(:, :) !< Occupied MO mask 2.
        logical, intent(in) :: act1(:, :) !< Active MO mask 1.
        logical, intent(in) :: act2(:, :) !< Active MO mask 2.
        real(dp), allocatable, intent(inout) :: moa1(:, :) !< MO coefficients alpha 1.
        real(dp), allocatable, intent(inout) :: moa2(:, :) !< MO coefficients alpha 2.
        real(dp), allocatable, intent(inout) :: mob1(:, :) !< MO coefficients beta 1.
        real(dp), allocatable, intent(inout) :: mob2(:, :) !< MO coefficients beta 2.
        integer :: rhf

        rhf = max(rhf1, rhf2)
        if (print_level >= 1) then
            if ((.not. all(act1)) .or. (.not. all(act2))) then
                write(stdout, *)
                write(stdout, '(1x,a)') 'Removing MOs frozen in electronic structure calculation...'
            end if
        end if
        call sort_mo(occ1(:, 1), act1(:, 1), moa1, remove_inactive=.true.)
        call sort_mo(occ2(:, 1), act2(:, 1), moa2, remove_inactive=.true.)
        if (rhf == 2) call sort_mo(occ1(:, rhf1), act1(:, rhf1), mob1, remove_inactive=.true.)
        if (rhf == 2) call sort_mo(occ2(:, rhf2), act2(:, rhf2), mob2, remove_inactive=.true.)
    end subroutine remove_frozen_mo



    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: check_rhf
    !
    !> @brief Check for overlap calculations between restricted and unrestricted sets of states.
    !> @details
    !! In case of mixed restricted/unrestricted calculation, the restricted MO and CI coefficients
    !! are used for both alpha and beta spin, and the CI coefficients are renormalized.
    !----------------------------------------------------------------------------------------------
    subroutine check_rhf(rhf1, rhf2, moa1, moa2, mob1, mob2, cisa1, cisa2, cisb1, cisb2)
        integer, intent(in) :: rhf1 !< Restricted (1) or unrestricted (2) calculation 1.
        integer, intent(in) :: rhf2 !< Restricted (1) or unrestricted (2) calculation 2.
        real(dp), allocatable, intent(in) :: moa1(:, :) !< MO coefficients alpha 1.
        real(dp), allocatable, intent(in) :: moa2(:, :) !< MO coefficients alpha 2.
        real(dp), allocatable, intent(inout) :: mob1(:, :) !< MO coefficients beta 1.
        real(dp), allocatable, intent(inout) :: mob2(:, :) !< MO coefficients beta 2.
        real(dp), allocatable, intent(inout) :: cisa1(:, :, :) !< CIS matrix alpha 1.
        real(dp), allocatable, intent(inout) :: cisa2(:, :, :) !< CIS matrix alpha 2.
        real(dp), allocatable, intent(inout) :: cisb1(:, :, :) !< CIS matrix beta 1.
        real(dp), allocatable, intent(inout) :: cisb2(:, :, :) !< CIS matrix beta 2.
        integer :: rhf

        rhf = max(rhf1, rhf2)
        if (rhf1 /= rhf) then
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(1x,a)') 'Unrestricted bra and restricted ket calculation...'
                write(stdout, '(5x,a)') 'Renormalizing bra CI coefficients...'
            end if
            cisa1 = cisa1 / sqrt(2.0_dp)
            allocate(mob1, source=moa1)
            allocate(cisb1, source=cisa1)
            cisb1 = -cisb1
        else if (rhf2 /= rhf) then
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(1x,a)') 'Restricted bra and unrestricted ket calculation...'
                write(stdout, '(5x,a)') 'Renormalizing ket CI coefficients...'
            end if
            cisa2 = cisa2 / sqrt(2.0_dp)
            allocate(mob2, source=moa2)
            allocate(cisb2, source=cisa2)
            cisb2 = -cisb2
        end if
    end subroutine check_rhf


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: check_occ_mo_norms
    !
    !> @brief Check bra/ket MO norms in ket/bra basis and freeze MOs with low norms.
    !> @details
    !----------------------------------------------------------------------------------------------
    subroutine check_mo_norms(s_mo, mo1, mo2, cis1, cis2, freeze, freeze_threshold)
        real(dp), allocatable, intent(inout) :: s_mo(:, :)
        real(dp), allocatable, intent(inout) :: mo1(:, :) !< MO coefficients alpha 1.
        real(dp), allocatable, intent(inout) :: mo2(:, :) !< MO coefficients alpha 2.
        real(dp), allocatable, intent(inout) :: cis1(:, :, :) !< CIS matrix alpha 1.
        real(dp), allocatable, intent(inout) :: cis2(:, :, :) !< CIS matrix alpha 2.
        logical, intent(in) :: freeze
        real(dp), intent(in) :: freeze_threshold
        integer :: no, n_tot1, n_tot2, n_freeze, n_active
        integer, allocatable :: seq(:)
        integer, allocatable :: active_occ_bra(:)
        integer, allocatable :: active_occ_ket(:)
        integer, allocatable :: active_bra(:)
        integer, allocatable :: active_ket(:)
        logical, allocatable :: active_mask(:)
        real(dp), allocatable :: bra_norms(:)
        real(dp), allocatable :: ket_norms(:)
        real(dp), allocatable :: s2(:, :)
        real(dp), allocatable :: tmp_s_mo(:, :)
        real(dp), allocatable :: tmp_mo1(:, :) !< MO coefficients alpha 1.
        real(dp), allocatable :: tmp_mo2(:, :) !< MO coefficients alpha 2.
        real(dp), allocatable :: tmp_cis1(:, :, :) !< CIS matrix alpha 1.
        real(dp), allocatable :: tmp_cis2(:, :, :) !< CIS matrix alpha 2.
        integer :: i, bra_min(1)


        no = size(cis1, 2)
        allocate(seq(1:no), source=[ (i, i=1, no) ])
        allocate(s2(1:no, 1:no), source=s_mo(1:no, 1:no)**2)
        allocate(ket_norms(1:no), source=sum(s2, 1))
        allocate(bra_norms(1:no), source=sum(s2, 2))

        if (print_level >= 2) then
            write(stdout, *)
            write(stdout, '(3x,a)') 'Norms of ket occupied orbitals in bra basis:'
            write(stdout, '(3x,8es11.3)') ket_norms
            write(stdout, '(3x,a)') 'Norms of bra occupied orbitals in ket basis:'
            write(stdout, '(3x,8es11.3)') bra_norms
        end if

        if (freeze) then
            active_mask = (ket_norms > freeze_threshold)
            n_active = count(active_mask)
            n_freeze = no - n_active
            if (n_active == 0) then
                write(stderr, *)
                write(stderr, '(1x,a)') 'Error. Norms of all occupied orbitals below freeze threshold value.'
                stop
            end if
            active_occ_ket = pack(seq, active_mask)
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(3x,a,e10.4,a)') 'Freezing ket orbitals with norm in bra basis &
                                                &smaller than ', freeze_threshold,'.'
                write(stdout, '(5x,a)') 'Indices of frozen orbitals:'
                write(stdout, '(5x,20i4)') pack(seq, .not. active_mask)
                write(stdout, '(5x,a)') 'Norms of frozen orbitals:'
                write(stdout, '(5x,8es11.3)') ket_norms(pack(seq, .not. active_mask))
            end if
            active_mask = .true.
            do i = 1, n_freeze
                bra_min = minloc(bra_norms, mask=active_mask)
                active_mask(bra_min(1)) = .false.
            end do
            active_occ_bra = pack(seq, active_mask)
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(3x,a,i0,a)') 'Freezing same number (', n_freeze, ') of bra orbitals &
                                              &with smallest norms in ket basis.'
                write(stdout, '(5x,a)') 'Indices of frozen orbitals:'
                write(stdout, '(5x,20i4)') pack(seq, .not. active_mask)
                write(stdout, '(5x,a)') 'Norms of frozen orbitals:'
                write(stdout, '(5x,8es11.3)') bra_norms(pack(seq, .not. active_mask))
            end if
            n_tot1 = n_active+size(cis1, 1)
            n_tot2 = n_active+size(cis2, 1)
            allocate(active_bra(1:n_tot1))
            allocate(active_ket(1:n_tot2))
            active_bra(1:n_active) = active_occ_bra
            active_ket(1:n_active) = active_occ_ket
            active_bra(n_active+1:) = [ (i, i=no+1, no+size(cis1, 1)) ]
            active_ket(n_active+1:) = [ (i, i=no+1, no+size(cis2, 1)) ]
            allocate(tmp_s_mo(1:n_tot1, 1:n_tot2), source=s_mo(active_bra, active_ket))
            allocate(tmp_mo1(1:size(mo1, 1), 1:n_tot1), source=mo1(:, active_bra))
            allocate(tmp_mo2(1:size(mo2, 1), 1:n_tot2), source=mo2(:, active_ket))
            allocate(tmp_cis1(1:size(cis1, 1), 1:n_active, 1:size(cis1, 3)), &
            &                 source=cis1(:, active_occ_bra, :))
            allocate(tmp_cis2(1:size(cis2, 1), 1:n_active, 1:size(cis2, 3)), &
            &                 source=cis2(:, active_occ_ket, :))
            deallocate(s_mo)
            deallocate(mo1)
            deallocate(mo2)
            deallocate(cis1)
            deallocate(cis2)
            allocate(s_mo, source=tmp_s_mo)
            allocate(mo1, source=tmp_mo1)
            allocate(mo2, source=tmp_mo2)
            allocate(cis1, source=tmp_cis1)
            allocate(cis2, source=tmp_cis2)
        end if

    end subroutine check_mo_norms



    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_prepend_occ
    !
    !> @brief Add occupied orbitals to start of a wf array
    !----------------------------------------------------------------------------------------------
    subroutine cis_prepend_occ(nadd, no, nv, cis)
        integer, intent(in) :: nadd
        integer, intent(in) :: no
        integer, intent(in) :: nv
        real(dp), allocatable, intent(inout) :: cis(:, :, :)
        real(dp), allocatable :: tmp_cis(:, :, :)

        allocate(tmp_cis(nv, no, size(cis,3)), source=0.0_dp)
        tmp_cis(:, nadd+1:, :) = cis
        deallocate(cis)
        allocate(cis, source=tmp_cis)
    end subroutine cis_prepend_occ



end module check_dimensions_overlap_mod
