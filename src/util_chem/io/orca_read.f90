!----------------------------------------------------------------------------------------------
! MODULE: orca_read_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date 
!
!> @brief 
!----------------------------------------------------------------------------------------------
module orca_read_mod
    use global_defs
    use file_mod
    implicit none

    private
    public :: orca_read_cis

contains

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: orca_read_cis
    !> @brief Read CIS coefficients and necessary dimensions directly from binary ORCA .cis file.
    !----------------------------------------------------------------------------------------------
    subroutine orca_read_cis(path, wfa, occ_mo, act_mo)
        use occupation_mod
        character(len=*), intent(in) :: path
        real(dp), allocatable :: wfa(:, :)
        logical, allocatable :: occ_mo(:, :)
        logical, allocatable :: act_mo(:, :)
        character(len=4) :: statekey
        character(len=:), allocatable :: inputfile
        integer :: dummy, vdim, cdim
        integer :: i, j, inunit, nex, p_ind, firstvirt
        type(occupation_numbers) :: onum

        p_ind=index(path, '.molden')
        statekey = 'cis'
        inputfile=path(1:p_ind)//statekey
        if (print_level >= 2) then 
            write(stdout, *) '  '
            write(stdout, *) 'Reading CIS coefficients from:  ', inputfile
        endif
        
        !!call need_file(inputfile, 'Error in orca_read_cis subroutine.')
        open(newunit=inunit, file=inputfile, form='UNFORMATTED', access='STREAM', action='READ')
        read(unit=inunit) nex
        read(unit=inunit) onum%fo(1)        ! first non-frozen occupied orbital
        read(unit=inunit) onum%o(1)
        onum%o(1) = onum%o(1)+1             ! ORCA starts counting from 0        
        read(unit=inunit) firstvirt         ! first non-frozen virtual orbital
        read(unit=inunit) onum%n
        onum%av(1) = onum%n - firstvirt + 1
        onum%ao(1) = onum%o(1) - onum%fo(1)
        onum%n = onum%n + 1
        onum%v(1) = onum%n - onum%o(1)
        onum%fv(1) = onum%v(1) - onum%av(1)

      ! write(stdout, *) 'nex ', nex
      ! write(stdout, *) 'onum', onum
        do i = 1, 4
            read(unit=inunit) dummy
            !write(stdout, *) i, dummy
        end do
      ! write(stdout, *)  onum%ao(1), onum%av(1), onum%ao(1)*onum%av(1), nex
        vdim = onum%ao(1)*onum%av(1)
      ! write(stdout, *) vdim, nex, 'vdim, nex, check ints'
        if (allocated(wfa)) deallocate(wfa)
        allocate(wfa(vdim, nex))

        do i = 1, nex
            read(unit=inunit) cdim
            if (print_level >= 2) write(stdout, *) 'Reading state:  ', i
            if (cdim /= vdim) call errstop('orca_read_cis', 'Unexpected array dimensions!', 1)
            do j = 1, 9
                read(unit=inunit) dummy   ! if needed 4 = current state, 6 = energy (check)
              ! write(stdout, *) dummy
            end do
            read(unit=inunit) wfa(:, i)
        end do
        close(inunit)

        allocate(occ_mo(onum%n, 2), source=.false.)
        allocate(act_mo(onum%n, 2), source=.true.)
        do i = 1, onum%o(1)
            occ_mo(i, 1) = .true.
        end do
        do i = 1, onum%fo(1)
            act_mo(i, 1) = .false.
        end do
        if (onum%fv(1) .NE. 0) then
       !    write(stdout, *) 'There were ', onum%fv(1), ' frozen virtual orbitals!'
            do i = (onum%o(1)+1), (onum%o(1)+onum%fv(1))
                act_mo(j, 1) = .false.
            end do
        endif
    end subroutine orca_read_cis

end module orca_read_mod
