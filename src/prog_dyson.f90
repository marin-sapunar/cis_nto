program cis_dyson_prog
    use global_defs
    use cis_dyson_mod
    use read_all_mod
    use write_molden_mod
    use write_txt_mod
    use orthog_mod
    use ccg_ao_mod
    use one_el_op_mod
    use blas95, only : gemm, dot, gemv
    implicit none

    integer :: rhf = 0 !< Restricted (1) or unrestricted(2) calculation.
    integer :: rhf1 = 0 !< Restricted (1) or unrestricted (2) for 1 wave functions.
    integer :: rhf2 = 0 !< Restricted (1) or unrestricted (2) for 2 wave functions.
    integer :: nwf1 !< Number of excited states 1.
    integer :: nwf2 !< Number of excited states 2.
    real(dp), allocatable :: trans1(:, :) !< Transformation matrix for basis functions 1.
    real(dp), allocatable :: trans2(:, :) !< Transformation matrix for basis functions 2.
    type(ccg), allocatable :: ccg1(:) !< Atomic orbitals 1.
    type(ccg), allocatable :: ccg2(:) !< Atomic orbitals 2.
    real(dp), allocatable :: moa1(:, :) !< Molecular orbital coefficients alpha 1.
    real(dp), allocatable :: moa2(:, :) !< Molecular orbital coefficients alpha 2.
    real(dp), allocatable :: mob1(:, :) !< Molecular orbital coefficients beta 1.
    real(dp), allocatable :: mob2(:, :) !< Molecular orbital coefficients beta 2.
    real(dp), allocatable :: wfa1(:, :, :) !< CI alpha single excitation coefficients 1.
    real(dp), allocatable :: wfa2(:, :, :) !< CI alpha single excitation coefficients 2.
    real(dp), allocatable :: wfb1(:, :, :) !< CI beta single excitation coefficients 1.
    real(dp), allocatable :: wfb2(:, :, :) !< CI beta single excitation coefficients 2.
    real(dp), allocatable :: s_ao(:, :) !< Atomic orbital overlaps.
    real(dp), allocatable :: s_mo(:, :, :) !< Molecular orbital overlaps.
    real(dp), allocatable :: dys_mo(:, :, :) !< Dyson orbitals in mo basis.
    real(dp), allocatable :: dys_ao(:, :, :) !< Dyson orbitals.
    real(dp), allocatable :: dys_norm(:, :) !< Norms of the dyson orbitals.
    real(dp), allocatable :: wrk(:, :)
    character(len=1000) :: temp
    character(len=:), allocatable :: input_format
    character(len=:), allocatable :: dir1
    character(len=:), allocatable :: dir2
    real(dp) :: thr
    integer :: outunit
    integer :: i, j

    input_format = 'turbomole'
    write(stdout, *) 'Path to N-1 electron calculation:'
    read(stdin, '(a)') temp
    dir1 = trim(adjustl(temp))
    write(stdout, *) 'Path to N electron calculation:'
    read(stdin, '(a)') temp
    dir2 = trim(adjustl(temp))
    write(stdout, *) 'Threshold for truncating wave functions:'
    read(stdin, *) thr

    call read_ccg_ao(input_format, dir1, ccg1, trans_ao=trans1)
    call read_ccg_ao(input_format, dir2, ccg2, trans_ao=trans2)
    call one_el_op(ccg1, ccg2, [0,0,0], trans1, trans2, s_ao)
    call write_txt('ao_ovl', s_ao)
    call write_txt('trans1', trans1)
    call write_txt('trans2', trans2)

    call read_mo(input_format, dir1, moa_c=moa1, mob_c=mob1)
    call read_mo(input_format, dir2, moa_c=moa2, mob_c=mob2)
    call write_txt('moa1', moa1)
    call write_txt('moa2', moa2)
    call write_txt('mob1', mob1)
    call write_txt('mob2', mob2)

    call read_cis(input_format, dir1, cisa=wfa1, cisb=wfb1, norm=.true.)
    call read_cis(input_format, dir2, cisa=wfa2, cisb=wfb2, norm=.true.)
    nwf1 = size(wfa1, 3)
    nwf2 = size(wfa2, 3)

    if (allocated(wfb1)) rhf1 = 2
    if (allocated(wfb2)) rhf2 = 2
    rhf = max(rhf1, rhf2)
    call write_txt('owfa1', wfa1)
    call write_txt('owfa2', wfa2)
    call write_txt('owfb1', wfb1)
    call write_txt('owfb2', wfb2)
    if (rhf1 /= rhf) then
        wfa1 = wfa1 / sqrt(2.0_dp)
        allocate(mob1, source=moa1)
        allocate(wfb1, source=wfa1)
        wfb1 = -wfb1
    else if (rhf2 /= rhf) then
        wfa2 = wfa2 / sqrt(2.0_dp)
        allocate(mob2, source=moa2)
        allocate(wfb2, source=wfa2)
        wfb2 = -wfb2
    end if
    call write_txt('nwfa1', wfa1)
    call write_txt('nwfa2', wfa2)
    call write_txt('nwfb1', wfb1)
    call write_txt('nwfb2', wfb2)

    allocate(wrk(size(moa1,2), size(moa2,1)))
    allocate(s_mo(size(moa1, 2), size(moa2, 2), rhf))
    call gemm(moa1, s_ao, wrk, transa='T')
    call gemm(wrk, moa2, s_mo(:, :, 1))
    if (rhf == 2) then
        call gemm(mob1, s_ao, wrk, transa='T')
        call gemm(wrk, mob2, s_mo(:, :, 2))
    end if
    call write_txt('s_mo', s_mo)

    call cis_dyson(thr, s_mo, wfa1, wfa2, wfb1, wfb2, dys_mo)
    allocate(dys_norm(nwf1+1, nwf2+1))
    dys_norm = sum(dys_mo(:, :, :)*dys_mo(:, :, :), dim=1)

    deallocate(wrk)
    allocate(wrk(size(trans2, 2), size(mob2, 2)))
    allocate(dys_ao(size(trans2, 2), nwf1, nwf2))
    call gemm(trans2, mob2(:, :), wrk, transa='T')
    do j = 1, nwf2
        call gemm(wrk, dys_mo(:, :, j), dys_ao(:, :, j))
    end do

    do i = 1, nwf1+1
        do j = 1, nwf2+1
            ! Write MO basis dyson orbital.
            write(temp, '(a,i3.3,a,i3.3,a)') 'dys.', i-1, '.', j-1, '.mo'
            open(newunit=outunit, file=temp, action='write')
            write(outunit, '(1e13.5)') dys_mo(:, i, j)
            close(outunit)
            ! Write AO basis dyson orbital.
            write(temp, '(a,i3.3,a,i3.3,a)') 'dys.', i-1, '.', j-1, '.ao'
            open(newunit=outunit, file=temp, action='write')
            call write_molden_mo_single(outunit, 1000*i+j, 'a   ', 1, 0.0_dp, dys_norm(i, j),      &
            &                           dys_ao(:, i, j))
            close(outunit)
        end do
    end do


end program cis_dyson_prog
