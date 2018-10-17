program cis_nto_prog
    use read_txt_mod
    use cis_nto_mod
    use write_molden_mod
    use blas95, only : gemm
    implicit none

    integer :: rhf = 0 !< Restricted (1) or unrestricted(2) calculation.
    integer :: nwf = 0 !< Number of wave functions 1.
    integer :: nao = 0 !< Number of atomic orbitals.
    integer :: nao_c = 0 !< Number of atomic orbitals (cartesian).
    integer :: nmo = 0 !< Number of molecular orbitals.
    integer :: noa = 0 !< Number of alpha occupied orbitals 1.
    integer :: nva = 0 !< Number of alpha virtual orbitals 1.
    integer :: nob = 0 !< Number of beta occupied orbitals 1.
    integer :: nvb = 0 !< Number of beta virtual orbitals 1.
    real(dp), allocatable :: mo(:, :, :) !< Molecular orbital coefficients.
                                         !! Dimensions: nao x nmo x rhf.
    real(dp), allocatable :: wfa(:, :, :) !< CI alpha single excitation coefficients.
                                           !! Dimensions: nva x noa x nwf.
    real(dp), allocatable :: wfb(:, :, :) !< CI beta single excitation coefficients.
                                           !! Dimensions: nvb x nob x nwf.
    real(dp), allocatable :: nto_a(:, :, :) !< NTOs alpha.
                                           !! Dimensions: nao x noa x nwf.
    real(dp), allocatable :: nto_b(:, :, :) !< NTOs beta.
                                           !! Dimensions: nao x nob x nwf.
    real(dp), allocatable :: nto_c_a(:, :) !< NTO alpha coefficients.
                                           !! Dimensions: noa x nwf.
    real(dp), allocatable :: nto_c_b(:, :) !< NTO beta coefficients.
                                           !! Dimensions: nob x nwf.
    integer, allocatable :: na_a(:) !< Number of active alpha NTO orbitals.
    integer, allocatable :: na_b(:) !< Number of active beta NTO orbitals.
    real(dp), allocatable :: c2s(:, :) !< Transformation between cartesian and spherical basis functions.
    real(dp), allocatable :: mo_c_a(:, :)
    real(dp), allocatable :: mo_c_b(:, :)
    logical :: beta = .false.
    integer :: i, j
    character(len=1000) :: temp
    real(dp) :: thr
    integer :: outunit
    character(len=100), parameter :: outfile='nto.molden'

    i = command_argument_count()
    thr = 1.0_dp
    if (i > 0) then
        call get_command_argument(1, temp)
        read(temp, *) thr
    end if

    call read_txt('mo', nao, nmo, rhf, mo)
    call read_txt('trans', nao, nao_c, c2s)
    if (rhf == 2)  beta = .true.
    call read_txt('wfa', nva, noa, nwf, wfa)
    allocate(mo_c_a(nao_c, nmo))
    call gemm(c2s, mo(:, :, 1), mo_c_a, transa='T')
    if (beta) then 
        call read_txt('wfb', nvb, nob, nwf, wfb)
        allocate(mo_c_b(nao_c, nmo))
        call gemm(c2s, mo(:, :, 2), mo_c_b, transa='T')
    end if

    call cis_nto_and_convert_basis(mo_c_a, wfa, nto_c_a, nto_a)
    if (beta) call cis_nto_and_convert_basis(mo_c_b, wfb, nto_c_b, nto_b)
    call cis_nto_truncate(beta, thr, nto_c_a, nto_c_b, na_a, na_b)
    open(newunit=outunit, file=outfile, action='write')
    write(outunit, '(a)') '[MO]'
    do i = 1, nwf
        do j = 1, na_a(i)
            call write_molden_mo_single(outunit, 1000*i+2*j-1, 'a   ', 1, 0.0_dp, nto_c_a(j, i),   &
            &                           nto_a(:, j, i))
            call write_molden_mo_single(outunit, 1000*i+2*j, 'a   ', 1, 0.0_dp, nto_c_a(j, i),     &
            &                           nto_a(:, noa+j, i))
        end do
        if (beta) then
            do j = 1, na_b(i)
                call write_molden_mo_single(outunit, 1000*i+2*j-1, 'a   ', 2, 0.0_dp, nto_c_b(j, i), &
                &                           nto_b(:, j, i))
                call write_molden_mo_single(outunit, 1000*i+2*j, 'a   ', 2, 0.0_dp, nto_c_b(j, i),   &
                &                           nto_b(:, nob+j, i))
            end do
        end if
    end do


end program cis_nto_prog
