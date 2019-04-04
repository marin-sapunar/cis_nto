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
        call read_mo(input_format_1, path1, moa_c=moa1, mob_c=mob1)
        call read_mo(input_format_2, path2, moa_c=moa2, mob_c=mob2)

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
        time_in =  time_in + omp_get_wtime() - time0

        ! Calculate overlaps.
        time0 = omp_get_wtime()
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
        if (print_level >= 3) then
            write(stdout, '(5x,a)') 'Diagonal of the alpha MO overlap matrix: '
            write(stdout, '(6x,15f8.4)') [ ( s_mo_a(i, i), i=1, size(s_mo_a, 1) ) ]
            if (allocated(s_mo_b)) then
                write(stdout, '(5x,a)') 'Diagonal of the beta MO overlap matrix: '
                write(stdout, '(10x,15(f10.5))') [ ( s_mo_b(i, i), i=1, size(s_mo_b, 1) ) ]
            end if
        end if
        time_mo = omp_get_wtime() - time0
    end subroutine overlap_mo


end module overlap_mo_mod
