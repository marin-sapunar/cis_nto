!----------------------------------------------------------------------------------------------
! MODULE: nto_cis_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date March, 2019
!----------------------------------------------------------------------------------------------
module nto_cis_mod
    use global_defs
    implicit none


    private
    public nto_cis


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: nto_cis
    !> @brief Read input and generate NTOs from a set of CIS type wave functions.
    !----------------------------------------------------------------------------------------------
    subroutine nto_cis
        use nto_variables
        use read_all_mod
        use cis_util_mod
        use cis_nto_mod
        logical :: beta

        ! Read input
        time0 =  omp_get_wtime()
        call read_geom(input_format, path, geom, atsym, atnum)
        call read_basis(input_format, path, bs)
        call read_mo(input_format, path, mos, bs)
        call read_cis(input_format, path, cisa=cisa, cisb=cisb, occ_mo=occ, act_mo=act, occ_num=on)
        beta = .false.
        if (allocated(mos%cb)) then
            rhf = 2
            beta = .true.
        end if
        call print_calculation_properties(size(cisa, 3), on, '')
        ! -- Sort MOs
        call mo_reshape(.true., .true., occ(:, 1), act(:, 1), 2, a2=mos%ca)
        if (rhf == 2) call mo_reshape(.true., .true., occ(:, 2), act(:, 2), 2, a2=mos%ca)
        time_in =  time_in + omp_get_wtime() - time0

        ! Generate NTOs
        call cis_nto_and_convert_basis(mos%ca, cisa, nto_c_a, nto_a)
        if (rhf == 2) call cis_nto_and_convert_basis(mos%cb, cisb, nto_c_b, nto_b)
        call cis_nto_truncate(beta, wf_threshold, truncate_nex, nto_c_a, nto_c_b, na_a, na_b)
    end subroutine nto_cis


end module nto_cis_mod
