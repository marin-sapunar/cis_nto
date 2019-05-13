!----------------------------------------------------------------------------------------------
! MODULE: cis_util_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date March, 2019
!
!> @brief Helper subroutines for handling CIS type wave functions.
!----------------------------------------------------------------------------------------------
module cis_util_mod
    use global_defs
    use occupation_mod
    implicit none



contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: print_calculation_properties
    !> @brief Print dimensions of an electronic structure calculation to stdout.
    !----------------------------------------------------------------------------------------------
    subroutine print_calculation_properties(n_ex_state, on, set_name)
        integer, intent(in) :: n_ex_state
        type(occupation_numbers), intent(in) :: on
        character(len=*), intent(in) :: set_name

        if (print_level >= 1) then
            write(stdout, *)
            write(stdout, '(1x,a)') 'Properties for '//set_name//' states:'
            if (on%rhf == 1) write(stdout, '(5x, a)') 'Restricted calculation.'
            if (on%rhf == 2) write(stdout, '(5x, a)') 'Unrestricted calculation.'
            write(stdout, '(5x,a,2(1x,i0))') 'Number of excited states:          ', n_ex_state
            write(stdout, '(5x,a,2(1x,i0))') 'Number of orbitals:                ', on%n
            write(stdout, '(5x,a,2(1x,i0))') 'Number of occupied orbitals:       ', on%o(1:on%rhf)
            write(stdout, '(5x,a,2(1x,i0))') 'Number of virtual orbitals:        ', on%v(1:on%rhf)
            write(stdout, '(5x,a,2(1x,i0))') 'Number of active occupied orbitals:', on%ao(1:on%rhf)
            write(stdout, '(5x,a,2(1x,i0))') 'Number of active virtual orbitals: ', on%av(1:on%rhf)
        end if
    end subroutine print_calculation_properties


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: check_rhf
    !
    !> @brief Check for calculations between restricted and unrestricted sets of states.
    !> @details
    !! In case of mixed restricted/unrestricted calculation, the restricted MO and CI coefficients
    !! are used for both alpha and beta spin, and the CI coefficients are renormalized.
    !----------------------------------------------------------------------------------------------
    subroutine check_rhf(rhf, target_rhf, cisa, cisb, set_name)
        integer, intent(in) :: rhf
        integer, intent(in) :: target_rhf
        real(dp), allocatable :: cisa(:, :, :)
        real(dp), allocatable :: cisb(:, :, :)
        character(len=*), intent(in) :: set_name
        if ((target_rhf == 2) .and. (rhf == 1)) then
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(1x,a)') 'Restricted '//set_name//' MOs in unrestricted calculation'
                write(stdout, '(5x,a)') 'Renormalizing CI coefficients...'
            end if
            cisa = cisa / sqrt(num2)
            allocate(cisb, source=cisa)
            cisb = -cisb
        end if
    end subroutine check_rhf


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: sort_mo
    !> @brief Sort MO arrays based on occupation and active status.
    !----------------------------------------------------------------------------------------------
    subroutine sort_mo(occ, act, mo_dim, mo, s_mo)
        logical, intent(in) :: occ(:)
        logical, intent(in) :: act(:)
        integer, intent(in) :: mo_dim
        real(dp), allocatable, intent(inout) :: mo(:, :) !< MO coefficients.
        real(dp), allocatable, intent(inout) :: s_mo(:, :) !< MO overlap matrix.

        if (allocated(mo)) call mo_reshape(.true., .false., occ, act, 2, a2=mo)
        if (allocated(s_mo)) call mo_reshape(.true., .false., occ, act, mo_dim, a2=s_mo)
    end subroutine sort_mo


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: check_occ_orb
    !
    !> @brief Check match for number of (active) occupied orbitals for overlap/dyson calculation
    !> @details
    !! Number of occupied orbitals of each spin must be equal for overlap calculations. For Dyson
    !! orbital calculations, there should be one extra beta electron in the ket states. If number
    !! of active occupied orbitals is not equal, frozen orbitals from both calculations are
    !! converted to active orbitals and CI arrays are expanded with excitations from the new active
    !! orbitals (with all coefficients equal to zero).
    !----------------------------------------------------------------------------------------------
    subroutine check_occ_orb(n_el_d, s, on1, on2, cis1, cis2)
        integer, intent(in) :: n_el_d !< Required difference in number of electrons.
        integer, intent(in) :: s !< Current spin.
        type(occupation_numbers), intent(inout) :: on1 !< Bra occupation numbers.
        type(occupation_numbers), intent(inout) :: on2 !< Ket occupation numbers.
        real(dp), allocatable, intent(inout) :: cis1(:, :, :) !< Bra CIS matrices.
        real(dp), allocatable, intent(inout) :: cis2(:, :, :) !< Ket CIS matrices.
        character(len=:), allocatable :: ab

        ! Check number of occupied orbitals. This has to be compatible.
        if ((on1%o(s) + n_el_d /= on2%o(s))) then
            write(stderr, *)
            write(stderr, '(1x,a,a,a)') 'Error. Mismatch in number of occupied orbitals.'
            stop
        end if

        ! Check number of active occupied orbitals.
        if ((on1%ao(s) + n_el_d /= on2%ao(s))) then
            if (print_level >= 1) then
                if (s == 1) ab = 'alpha'
                if (s == 2) ab = 'beta'
                write(stdout, *)
                write(stdout, '(1x,a,a,a)') 'Warning. Mismatch in number of active orbitals.'
                write(stdout, '(3x,i0,a,i0,a)') on1%ao(s), ' bra and ', on2%ao(s), ' ket ', ab, &
                &                               ' active occupied orbitals.'
                write(stdout, '(3x,a)') 'Using previously frozen orbitals as active orbitals.'
            end if
            call prepend_zero_occ_cis(on1%fo(s), cis1)
            call prepend_zero_occ_cis(on2%fo(s), cis2)
            on1%ao(s) = on1%o(s)
            on2%ao(s) = on2%o(s)
            on1%fo(s) = 0
            on2%fo(s) = 0
        end if
    end subroutine check_occ_orb


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: remove_frozen_mo
    !> @brief Remove frozen MOs from arrays.
    !----------------------------------------------------------------------------------------------
    subroutine remove_frozen_mo(fo, fv, mo_dim, mo, s_mo)
        integer, intent(in) :: fo
        integer, intent(in) :: fv
        integer, intent(in) :: mo_dim
        real(dp), allocatable, intent(inout) :: mo(:, :) !< MO coefficients.
        real(dp), allocatable, intent(inout) :: s_mo(:, :) !< MO overlap matrix.
        logical, allocatable :: act(:)
        integer :: n

        if (.not. allocated(mo)) return
        n = size(mo, 2)
        allocate(act(n), source=.true.)
        if (fo > 0) act(1:fo) = .false.
        if (fv > 0) act(n-fv:n) = .false.
        if (allocated(mo)) call mo_reshape(.false., .true., act=act, mo_dim=2, a2=mo)
        if (allocated(s_mo)) call mo_reshape(.false., .true., act=act, mo_dim=mo_dim, a2=s_mo)
    end subroutine remove_frozen_mo


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: freeze_mo_by_norm
    !> @brief Select MOs with lowest norms and remove them from arrays.
    !----------------------------------------------------------------------------------------------
    subroutine freeze_mo_by_norm(threshold, spin, mo_dim, mo, s_mo, cis, n_freeze, first_n)
        real(dp), intent(in) :: threshold !< Freeze threshold
        integer, intent(in) :: spin !< Spin (just for print message).
        integer, intent(in) :: mo_dim !< Dimension of MOs in overlap matrix.
        real(dp), allocatable, intent(inout) :: mo(:, :) !< MO coefficients.
        real(dp), allocatable, intent(inout) :: s_mo(:, :) !< MO overlap matrix.
        real(dp), allocatable, intent(inout) :: cis(:, :, :) !< CIS matrix.
        integer, intent(inout) :: n_freeze !< Number of frozen orbitals.
        logical, intent(in) :: first_n !< Freeze first n_freeze orbitals instead of by threshold
        real(dp), allocatable :: norms(:)
        real(dp), allocatable :: s2(:, :)
        character(len=:), allocatable :: cset
        character(len=:), allocatable :: oset
        character(len=:), allocatable :: ab
        logical, allocatable :: act(:)
        integer, allocatable :: seq(:)
        integer :: no, nv, i, mpos(1)

        if (.not. allocated(s_mo)) return
        nv = size(cis, 1)
        no = size(cis, 2)
        allocate(act(1:no+nv), source=.true.)
        allocate(seq(1:no), source=[ (i, i = 1, no) ])
        allocate(s2, source=s_mo**2)
        if (mo_dim == 1) then
            cset = 'bra'
            oset = 'ket'
            allocate(norms(1:no), source=sum(s2(1:no, :), 2)) 
        else
            cset = 'ket'
            oset = 'bra'
            allocate(norms(1:no), source=sum(s2(:, 1:no), 1)) 
        end if

        if (print_level >= 2) then
            write(stdout, *)
            if (spin == 1) then
                ab = 'alpha'
            else
                ab = 'beta'
            end if
            write(stdout, '(3x,7a)') 'Norms of ', cset, ' spin ', ab, ' occupied orbitals in',     &
            &                        oset, ' basis:'
            write(stdout, '(3x,8es11.3)') norms
        end if

        if (threshold < 0.0_dp) return

        if (.not. first_n) then
            ! Find orbitals to freeze
            act(1:no) = (norms > threshold)
            n_freeze = count(.not. act)
            if (n_freeze == no) then
                write(stderr, *)
                write(stderr, '(1x,a)') 'Error. Norms of all orbitals below freeze threshold value.'
                stop
            end if
            ! Print frozen orbitals
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(5x,a,e10.4,a)') 'Freezing orbitals with norm below ', threshold,'.'
                if (all(act)) then
                    write(stdout, '(9x,a)') 'No orbitals frozen.'
                else
                    write(stdout, '(9x,a)') 'Indices of frozen orbitals:'
                    write(stdout, '(9x,20i4)') pack(seq, .not. act(1:no))
                    write(stdout, '(9x,a)') 'Norms of frozen orbitals:'
                    write(stdout, '(9x,8es11.3)') norms(pack(seq, .not. act(1:no)))
                end if
            end if
        else
            if (n_freeze == 0) return
            do i = 1, n_freeze
                mpos = minloc(norms, mask=act(1:no))
                act(mpos(1)) = .false.
            end do
            
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(5x,a,i0,a)') 'Freezing ', n_freeze, ' orbitals with lowest norms.'
                write(stdout, '(9x,a)') 'Indices of frozen orbitals:'
                write(stdout, '(9x,20i4)') pack(seq, .not. act)
                write(stdout, '(9x,a)') 'Norms of frozen orbitals:'
                write(stdout, '(9x,8es11.3)') norms(pack(seq, .not. act))
            end if

        end if

        if (allocated(mo)) call mo_reshape(.false., .true., act=act, mo_dim=2, a2=mo)
        if (allocated(s_mo)) call mo_reshape(.false., .true., act=act, mo_dim=mo_dim, a2=s_mo)
        if (allocated(cis)) call mo_reshape(.false., .true., act=act(1:no), mo_dim=1, a3=cis)
    end subroutine freeze_mo_by_norm


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: check_ci_norms
    !> @brief Check norms of CI vectors after freezing MOs. Warn if norm of any state becomes low.
    !----------------------------------------------------------------------------------------------
    subroutine check_ci_norms(rhf, cisa, cisb, set_name)
        integer, intent(in) :: rhf !< Restricted (1) or unrestricted (2) calculation.
        real(dp), allocatable, intent(in) :: cisa(:, :, :) !< CIS matrix alpha.
        real(dp), allocatable, intent(in) :: cisb(:, :, :) !< CIS matrix beta.
        character(len=*), intent(in) :: set_name
        real(dp), allocatable :: norms(:)
        real(dp), parameter :: thr = 0.5_dp

        allocate(norms(1:size(cisa, 3)))
        norms = sum(sum(cisa**2, 1), 1)
        if (rhf == 2) norms = norms + sum(sum(cisb**2, 1), 1)

        if (print_level >= 2) then
            write(stdout, *)
            write(stdout, '(3x,a,a,a)') 'Initial norms of ', set_name,' excited states:'
            write(stdout, '(5x,8es11.3)') norms
            if (any(norms < thr)) then
                write(stdout, '(5x,a)') 'Warning. States with very low initial norm present.'
                write(stdout, '(5x,a)') '         Check for excitations to/from frozen orbitals.'
            end if
        end if
    end subroutine check_ci_norms


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: prepend_zero_occ_cis
    !> @brief Add occupied orbitals (with all zero coefficients) to start of a wf array
    !----------------------------------------------------------------------------------------------
    subroutine prepend_zero_occ_cis(nadd, cis)
        integer, intent(in) :: nadd
        real(dp), allocatable, intent(inout) :: cis(:, :, :)
        real(dp), allocatable :: tmp_cis(:, :, :)

        allocate(tmp_cis(1:size(cis, 1), 1:size(cis, 2)+nadd, 1:size(cis, 3)), source=0.0_dp)
        tmp_cis(:, nadd+1:, :) = cis
        deallocate(cis)
        allocate(cis, source=tmp_cis)
    end subroutine prepend_zero_occ_cis


end module cis_util_mod
