program cis_dyson_prog
    use global_defs
    use cis_dyson_mod
    use read_txt_mod
    use orthog_mod
    use blas95, only : gemm, dot, gemv
    implicit none

    integer :: rhf = 0 !< Restricted (1) or unrestricted(2) calculation.
    integer :: rhf1 = 0 !< Restricted (1) or unrestricted (2) for 1 wave functions.
    integer :: rhf2 = 0 !< Restricted (1) or unrestricted (2) for 2 wave functions.
    integer :: nao1 = 0 !< Number of atomic orbitals 1.
    integer :: nao2 = 0 !< Number of atomic orbitals 2.
    integer :: nmo1 = 0 !< Number of molecular orbitals 1.
    integer :: nmo2 = 0 !< Number of molecular orbitals 2.
    integer :: nwf1 = 0 !< Number of wave functions 1.
    integer :: nwf2 = 0 !< Number of wave functions 2.
    integer :: noa1 = 0 !< Number of alpha occupied orbitals 1.
    integer :: noa2 = 0 !< Number of alpha occupied orbitals 2.
    integer :: nva1 = 0 !< Number of alpha virtual orbitals 1.
    integer :: nva2 = 0 !< Number of alpha virtual orbitals 1.
    integer :: nob1 = 0 !< Number of beta occupied orbitals 1.
    integer :: nob2 = 0 !< Number of beta occupied orbitals 2.
    integer :: nvb1 = 0 !< Number of beta virtual orbitals 1.
    integer :: nvb2 = 0 !< Number of beta virtual orbitals 1.
    real(dp), allocatable :: mo1(:, :, :) !< Molecular orbital coefficients 1.
                                          !! Dimensions: nao1 x nmo1 x rhf1.
    real(dp), allocatable :: mo2(:, :, :) !< Molecular orbital coefficients 2.
                                          !! Dimensions: nao2 x nmo2 x rhf2.
    real(dp), allocatable :: wfa1(:, :, :) !< CI alpha single excitation coefficients 1.
                                           !! Dimensions: nva1 x noa1 x nwf1.
    real(dp), allocatable :: wfa2(:, :, :) !< CI alpha single excitation coefficients 2.
                                           !! Dimensions: nva2 x noa2 x nwf2.
    real(dp), allocatable :: wfb1(:, :, :) !< CI beta single excitation coefficients 1.
                                           !! Dimensions: nvb1 x nob1 x nwf1.
    real(dp), allocatable :: wfb2(:, :, :) !< CI beta single excitation coefficients 2.
                                           !! Dimensions: nvb2 x nob2 x nwf2.
    real(dp), allocatable :: s_ao(:, :) !< Atomic orbital overlaps.
                                        !! Dimensions: nao1 x nao2.
    real(dp), allocatable :: s_mo(:, :, :) !< Molecular orbital overlaps.
                                           !! Dimensions: nmo1 x nmo2 x rhf
    real(dp), allocatable :: dys_mo(:, :, :) !< Dyson orbitals.
    real(dp), allocatable :: dys_ao(:, :, :) !< Dyson orbitals.
    real(dp), allocatable :: wrk(:, :)
    integer :: i, j, s
    character(len=1000) :: temp
    real(dp) :: thr
    integer :: outunit

    i = command_argument_count()
    thr = 1.0_dp
    if (i > 0) then
        call get_command_argument(1, temp)
        read(temp, *) thr
    end if

    call read_txt('s_ao', nao1, nao2, s_ao)
    call read_txt('mo1', nao1, nmo1, rhf1, mo1)
    call read_txt('mo2', nao2, nmo2, rhf2, mo2)

    call read_txt('wfa1', nva1, noa1, nwf1, wfa1)
    call read_txt('wfa2', nva2, noa2, nwf2, wfa2)
    if (rhf1 == 2) call read_txt('wfb1', nvb1, nob1, nwf1, wfb1)
    if (rhf2 == 2) call read_txt('wfb2', nvb2, nob2, nwf2, wfb2)
    rhf = max(rhf1, rhf2)
    if (rhf1 /= rhf) then
        allocate(wfb1(nva1, noa1, nwf1))
        wfa1 = wfa1 / sqrt(2.0_dp)
        wfb1 = -wfa1
    else if (rhf2 /= rhf) then
        allocate(wfb2(nva2, noa2, nwf2))
        wfa2 = wfa2 / sqrt(2.0_dp)
        wfb2 = -wfa2
    end if

    allocate(wrk(nmo1, nao2))
    allocate(s_mo(nmo1, nmo2, rhf))
    do s = 1, rhf
        call gemm(mo1(:, :, min(s, rhf1)), s_ao, wrk, transa='T')
        call gemm(wrk, mo2(:, :, min(s, rhf2)), s_mo(:, :, s))
    end do

    call cis_dyson(thr, s_mo, wfa1, wfa2, wfb1, wfb2, dys_mo)

    do i = 1, nwf1+1
        do j = 1, nwf2+1
            write(temp, '(a,i3.3,a,i3.3)') 'dys.', i-1, '.', j-1
            open(newunit=outunit, file=temp, action='write')
            write(outunit, *) dys_mo(:, i, j)
            close(outunit)
        end do
    end do


end program cis_dyson_prog
