!----------------------------------------------------------------------------------------------
! MODULE: overlap_ao_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date March, 2019
!
!> @brief Subroutine for calculating overlap matrix between two sets of atomic orbitals.
!----------------------------------------------------------------------------------------------
module overlap_ao_mod
    use global_defs
    implicit none


    private
    public overlap_ao


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: overlap_ao
    !> @brief Calculate overlap matrix between two sets of atomic orbitals.
    !----------------------------------------------------------------------------------------------
    subroutine overlap_ao()
        use overlap_variables
        use basis_set_mod
        use one_el_op_mod
        use read_all_mod

        time0 = omp_get_wtime()
        call read_geom(input_format_1, path1, geom1)
        call read_geom(input_format_2, path2, geom2)
        call read_basis(input_format_1, path1, bs1)
        call read_basis(input_format_2, path2, bs2)
        time_in =  time_in + omp_get_wtime() - time0

        time0 = omp_get_wtime()
        if (print_level >= 2) then
            write(stdout, *)
            write(stdout, '(1x,a)') 'Computing AO overlaps...'
        end if
        call one_el_op(geom1, geom2, bs1, bs2, [0,0,0], bs1%sphe_mo, bs2%sphe_mo, .true., .true.,  &
        &              s_ao, center_atoms, center_pairs)
        time_ao = omp_get_wtime() - time0
    end subroutine overlap_ao


end module overlap_ao_mod
