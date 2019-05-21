!----------------------------------------------------------------------------------------------
! MODULE: molden_write_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date June, 2017
!
!> @brief Contains subroutines for writing a Molden format file.
!> @details
!! The Molden file format is organized into sections marked by keywords (ex. [Atoms]). Each 
!! subroutine in this module writes a section of the file. 
!----------------------------------------------------------------------------------------------
module molden_write_mod
    use global_defs
    implicit none

    private
    public :: molden_write_atoms
    public :: molden_write_gto
    public :: molden_write_mo
    public :: molden_write_mo_single


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: molden_write_gto
    !> @brief Write the molden [Atoms] section.
    !----------------------------------------------------------------------------------------------
    subroutine molden_write_atoms(outunit, geom, at_symbol, at_number, angstrom)
        integer, intent(in) :: outunit
        real(dp), intent(in) :: geom(:)
        character(len=2), intent(in) :: at_symbol(:)
        integer, intent(in) :: at_number(:)
        logical, intent(in), optional :: angstrom
        character(len=*), parameter :: fatm = '(1x, a2, 2(1x,i4), 3(1x,e22.16))'
        logical :: tangst
        character(len=:), allocatable :: lenunit
        real(dp) :: lenconv
        integer :: i

        tangst = .false.
        if (present(angstrom)) tangst = angstrom
        if (tangst) then
            lenunit = 'Angs'
            lenconv = 1.88972612456 !< @todo: Add unit conversion constants module.
        else
            lenunit = 'AU'
            lenconv = 1.0_dp
        end if

        write(outunit, '(a)') '[Atoms] '//lenunit
        do i = 1, size(at_symbol)
            write(outunit, fatm) at_symbol(i), i, at_number(i), geom(3*i-2:3*i)*lenconv
        end do
    end subroutine molden_write_atoms


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: molden_write_gto
    !> @brief Write the molden [GTO] section.
    !----------------------------------------------------------------------------------------------
    subroutine molden_write_gto(outunit, bs)
        use basis_set_mod, only : basis_set
        integer, intent(in) :: outunit !< Unit in which the file is open.
        type(basis_set), intent(in) :: bs !< Basis sets.
        integer :: i, j, k, ii
        character(len=*), parameter :: fibas='(1x, i4, 1x, i1)'
        character(len=*), parameter :: fprim='(1x, a, 3x, i2, 1x, f4.2)'
        character(len=*), parameter :: fcoef='(2(1x, e20.14))'

        write(outunit, '(a)') '[GTO]'
        do i = 1, bs%n_center
            ii = bs%center_i_bs(i)
            write(outunit, fibas) i, 0
            do j = 1, bs%bs(ii)%n_subshell
                write(outunit, fprim) bs%bs(ii)%cg(j)%typ, bs%bs(ii)%cg(j)%n_prim, 1.00_dp
                do k = 1, bs%bs(ii)%cg(j)%n_prim
                    write(outunit, fcoef) bs%bs(ii)%cg(j)%z(k), bs%bs(ii)%cg(j)%b(k)
                end do
            end do
            write(outunit, '(a)') ' '
        end do
    end subroutine molden_write_gto

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: molden_write_mo
    !> @brief Write set of molecular orbitals to molden [MO] section.
    !----------------------------------------------------------------------------------------------
    subroutine molden_write_mo(outunit, mos)
        use molecular_orbitals_mod, only : molecular_orbitals
        integer, intent(in) :: outunit
        type(molecular_orbitals), intent(in) :: mos
        integer :: i

        write(outunit, '(a)') '[MO]'
        do i = 1, mos%n_mo_a
            call molden_write_mo_single(outunit, mos%na(i), 'a   ', 1, mos%ea(i), mos%oa(i),       &
            &                           mos%ca(:, i))
        end do
        do i = 1, mos%n_mo_b
            call molden_write_mo_single(outunit, mos%nb(i), 'a   ', 1, mos%eb(i), mos%ob(i),       &
            &                           mos%cb(:, i))
        end do
    end subroutine molden_write_mo


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: molden_write_mo_single
    !> @brief Write a single molecular orbital to the [MO] section.
    !----------------------------------------------------------------------------------------------
    subroutine molden_write_mo_single(outunit, n, sym, spin, ene, occup, coef)
        integer, intent(in) :: outunit !< Unit in which the file is open.
        integer, intent(in) :: n !< Number of the MO.
        character(len=*), intent(in) :: sym !< Symmetry group.
        integer, intent(in) :: spin !< Spin.
        real(dp), intent(in) :: ene !< Energy.
        real(dp), intent(in) :: occup !< Occupation number.
        real(dp), intent(in) :: coef(:) !< Coefficients.
        character(len=19), parameter :: fsym="(1x, 'Sym=', i6, a)"
        character(len=24), parameter :: fene="(1x, 'Ene=', 1x, e20.14)"
        character(len=20), parameter :: fspi="(1x, 'Spin=', 1x, a)"
        character(len=24), parameter :: focc="(1x, 'Occup=', 1x, f8.6)"
        character(len=16), parameter :: fcoe="(i6, 1x, e20.14)"
        integer :: i

        write(outunit, fsym) n, sym
        write(outunit, fene) ene
        if (spin == 1) write(outunit, fspi) 'Alpha'
        if (spin == 2) write(outunit, fspi) 'Beta'
        write(outunit, focc) occup
        do i = 1, size(coef)
            write(outunit, fcoe) i, coef(i)
        end do
    end subroutine molden_write_mo_single


end module molden_write_mod
