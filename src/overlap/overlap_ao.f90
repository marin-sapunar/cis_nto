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
        integer :: i

        time0 = omp_get_wtime()
        call read_geom(input_format_1, path1, geom1, atom_symbol, atom_number)
        call read_geom(input_format_2, path2, geom2)
        call read_basis(input_format_1, path1, bs1)
        call read_basis(input_format_2, path2, bs2)
        time_in =  time_in + omp_get_wtime() - time0

        time0 = omp_get_wtime()
        if (print_level >= 2) then
            write(stdout, *)
            write(stdout, '(1x,a)') 'Computing AO overlaps...'
            write(stdout, '(5x,a,2(1x,i0))') 'Number of bra cartesian/spherical AOs:',             &
            &                                 bs1%n_bf_cart, bs1%n_bf_sphe
            write(stdout, '(5x,a,2(1x,i0))') 'Number of ket cartesian/spherical AOs:',             &
            &                                 bs2%n_bf_cart, bs2%n_bf_sphe
        end if
        if (print_level >= 4) then
            if (bs1%sphe_mo) then
                write(stdout, *) 'Expecting bra MOs written in terms of spherical AOs.'
            else
                write(stdout, *) 'Expecting bra MOs written in terms of cartesian AOs.'
            end if
            if (bs2%sphe_mo) then
                write(stdout, *) 'Expecting ket MOs written in terms of spherical AOs.'
            else
                write(stdout, *) 'Expecting ket MOs written in terms of cartesian AOs.'
            end if
        end if
        call one_el_op(geom1, geom2, bs1, bs2, operator_string, s_ao, center_atoms, center_pairs)
        if (print_level >= 3) then
            write(stdout, '(5x,a)') 'Diagonal of the AO overlap matrix: '
            write(stdout, '(6x,15f8.4)') [ ( s_ao(i, i), i=1, min(size(s_ao, 1), size(s_ao, 2)) ) ]
        end if
        time_ao = omp_get_wtime() - time0
    end subroutine overlap_ao


end module overlap_ao_mod
