!----------------------------------------------------------------------------------------------
! MODULE: overlap_cis_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date March, 2019
!
!> @brief Subroutine for calculating overlap matrix between two sets of molecular orbitals..
!----------------------------------------------------------------------------------------------
module overlap_cis_mod
    use global_defs
    implicit none


    private
    public overlap_cis


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: overlap_cis
    !> @brief Calculate overlap matrix between two sets of CIS wave functions.
    !----------------------------------------------------------------------------------------------
    subroutine overlap_cis()
        use overlap_variables
        use cis_overlap_nto_mod
        use cis_overlap_cis_mod
        use cis_overlap_l2m_mod
        use cis_dyson_nto_mod
        use cis_util_mod
        use read_all_mod
        integer :: n_freeze_a
        integer :: n_freeze_b
        integer :: diff_a
        integer :: diff_b

        diff_a = 0
        diff_b = 0
        if (dyson_c) diff_b = 1

        ! Read CIS WFs.
        time0 =  omp_get_wtime()

        call read_cis(input_format_1, path1, cisa=cisa1, cisb=cisb1, occ_mo=occ1, act_mo=act1,     &
        &             occ_num=on1, norm=norm_states, orthog=orth_states)
        call print_calculation_properties(size(cisa1, 3), on1, 'bra')

        call read_cis(input_format_2, path2, cisa=cisa2, cisb=cisb2, occ_mo=occ2, act_mo=act2,     &
        &             occ_num=on2, norm=norm_states, orthog=orth_states)
        call print_calculation_properties(size(cisa2, 3), on2, 'ket')

        ! Dimension checks
        rhf = max(on1%rhf, on2%rhf)
        ! -- Allocate beta WF coefficients for mixed restricted/unrestricted calculations.
        call check_rhf(on1%rhf, rhf, cisa1, cisb1, 'bra')
        call check_rhf(on2%rhf, rhf, cisa2, cisb2, 'ket')
        ! -- Sort MOs.
        call sort_mo(occ1(:, 1), act1(:, 1), 1, mos1%ca, s_mo_a)
        call sort_mo(occ2(:, 1), act2(:, 1), 2, mos2%ca, s_mo_a)
        call sort_mo(occ1(:, 2), act1(:, 2), 1, mos1%cb, s_mo_b)
        call sort_mo(occ2(:, 2), act2(:, 2), 2, mos2%cb, s_mo_b)
        deallocate(occ1)
        deallocate(occ2)
        deallocate(act1)
        deallocate(act2)
        ! -- Check number of occupied orbitals.
        call check_occ_orb(diff_a, 1, on1, on2, cisa1, cisa2)
        call check_occ_orb(diff_b, 2, on1, on2, cisa1, cisa2)
        ! -- Remove MOs frozen from el. struct. calculation.
        call remove_frozen_mo(on1%fo(1), on1%fv(1), 1, mos1%ca, s_mo_a)
        call remove_frozen_mo(on2%fo(1), on2%fv(1), 2, mos2%ca, s_mo_a)
        call remove_frozen_mo(on1%fo(2), on1%fv(2), 1, mos1%cb, s_mo_b)
        call remove_frozen_mo(on2%fo(2), on2%fv(2), 2, mos2%cb, s_mo_b)
        ! -- Freeze additional MOs based on norm threshold.
        call freeze_mo_by_norm(freeze_mo_norm_t, 1, 1, mos1%ca, s_mo_a, cisa1, n_freeze_a, .false.)
        call freeze_mo_by_norm(freeze_mo_norm_t, 1, 2, mos2%ca, s_mo_a, cisa2, n_freeze_a, .true.)
        call freeze_mo_by_norm(freeze_mo_norm_t, 2, 1, mos1%cb, s_mo_b, cisb1, n_freeze_b, .false.)
        call freeze_mo_by_norm(freeze_mo_norm_t, 2, 2, mos2%cb, s_mo_b, cisb2, n_freeze_b, .true.)
        ! -- Check CI vector norms.
        call check_ci_norms(rhf, cisa1, cisb1, 'bra')
        call check_ci_norms(rhf, cisa2, cisb2, 'ket')
        time_in =  time_in + omp_get_wtime() - time0

        ! Start calculation.
        time0 = omp_get_wtime()
        if (dyson_c) then
            if (print_level >= 2) then
                write(stdout, *)
                write(stdout, '(1x,a)') 'Computing Dyson orbitals overlaps...'
                write(stdout, *)
                write(stdout, '(5x,a,a,a)') 'Using ', trim(cis_algorithm), ' algorithm.'
            end if
            select case(cis_algorithm)
            case('NTO')
                call cis_dyson_nto(wf_threshold, truncate_nex, s_mo_a, s_mo_b, cisa1, cisa2, cisb1, cisb2, dyson_mo)
            case default
                write(stderr, *)
                write(stderr, *) 'Error. Selected algorithm not implemented.'
                stop
            end select
        else
            if (print_level >= 2) then
                write(stdout, *)
                write(stdout, '(1x,a)') 'Computing WF overlaps...'
                write(stdout, *)
                write(stdout, '(5x,a,a,a)') 'Using ', trim(cis_algorithm), ' algorithm.'
            end if
            select case(cis_algorithm)
            case('CIS')
                call cis_overlap_cis(rhf, wf_threshold, s_mo_a, s_mo_b, cisa1, cisa2, cisb1, cisb2, s_wf)
            case('L2M')
                call cis_overlap_l2m(rhf, s_mo_a, s_mo_b, cisa1, cisa2, cisb1, cisb2, s_wf)
            case('NTO')
                call cis_overlap_nto(rhf, wf_threshold, truncate_nex, s_mo_a, s_mo_b, cisa1, cisa2, cisb1, cisb2, s_wf)
            case default
                write(stderr, *)
                write(stderr, *) 'Error. Selected algorithm not implemented.'
                stop
            end select
        end if
        deallocate(cisa1)
        deallocate(cisa2)
        if (allocated(cisb1)) deallocate(cisb1)
        if (allocated(cisb2)) deallocate(cisb2)
        time_wf = omp_get_wtime() - time0
    end subroutine overlap_cis


end module overlap_cis_mod
