!----------------------------------------------------------------------------------------------
! MODULE: overlap_mo_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date March, 2019
!
!> @brief Subroutine for calculating overlap matrix between two sets of molecular orbitals..
!----------------------------------------------------------------------------------------------
module overlap_mo_mod
    use global_defs
    implicit none


    private
    public overlap_mo


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: overlap_mo
    !> @brief Calculate overlap matrix between two sets of molecular orbitals.
    !----------------------------------------------------------------------------------------------
    subroutine overlap_mo()
        use overlap_variables
        use matrix_mod, only : mat_ge_mmm
        use read_all_mod
        logical :: uhf1, uhf2
        integer :: i

        if (allocated(s_mo_a)) deallocate(s_mo_a)
        if (allocated(s_mo_b)) deallocate(s_mo_b)

        ! Read MOs.
        time0 =  omp_get_wtime()
        call read_mo(input_format_1, path1, mos1, bs1)
        call read_mo(input_format_2, path2, mos2, bs2)

        ! Check unrestricted/restricted calculations.
        uhf1 = (mos1%n_mo_b > 0)
        uhf2 = (mos2%n_mo_b > 0)
        if (uhf1 .and. (.not. uhf2)) then
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(1x,a)') 'Unrestricted bra and restricted ket MOs...'
            end if
            call mos2%to_unrestricted()
        else if (uhf2 .and. (.not. uhf1)) then
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(1x,a)') 'Restricted bra and unrestricted ket MOs...'
            end if
            call mos1%to_unrestricted()
        end if
        time_in =  time_in + omp_get_wtime() - time0

        ! Calculate overlaps.
        time0 = omp_get_wtime()
        if (print_level >= 2) then
            write(stdout, *)
            write(stdout, '(1x,a)') 'Computing MO overlaps...'
            write(stdout, '(5x,a,2(1x,i0))') 'Number of bra of AOs/MOs:', shape(mos1%ca)
            write(stdout, '(5x,a,2(1x,i0))') 'Number of ket of AOs/MOs:', shape(mos2%ca)
        end if
        allocate(s_mo_a(mos1%n_mo_a, mos2%n_mo_a))
        call mat_ge_mmm(mos1%ca, s_ao, mos2%ca, s_mo_a, transa='T')
        if (uhf1 .or. uhf2) then
            allocate(s_mo_b(mos1%n_mo_b, mos2%n_mo_b))
            call mat_ge_mmm(mos1%cb, s_ao, mos2%cb, s_mo_b, transa='T')
        end if
        if (print_level >= 3) then
            write(stdout, '(5x,a)') 'Diagonal of the alpha MO overlap matrix: '
            write(stdout, '(6x,15f8.4)') [( s_mo_a(i, i), i=1, min(size(s_mo_a, 1), size(s_mo_a, 2)) )]
            if (allocated(s_mo_b)) then
                write(stdout, '(5x,a)') 'Diagonal of the beta MO overlap matrix: '
                write(stdout, '(10x,15(f10.5))') [( s_mo_b(i, i), i=1, min(size(s_mo_b, 1), size(s_mo_b, 2)) )]
            end if
        end if
        time_mo = omp_get_wtime() - time0
    end subroutine overlap_mo


end module overlap_mo_mod
