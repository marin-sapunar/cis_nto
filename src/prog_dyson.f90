program cis_dyson_prog
    use global_defs
    use cis_dyson_mod
    use read_all_mod
    use write_molden_mod
    use ccg_ao_mod
    use atom_basis_mod
    use one_el_op_mod
    use matrix_mod
    use occupation_mod
    use blas95, only : gemm
    implicit none

    integer :: rhf = 0 !< Restricted (1) or unrestricted(2) calculation.
    integer :: rhf1 = 1 !< Restricted (1) or unrestricted (2) for 1 wave functions.
    integer :: rhf2 = 1 !< Restricted (1) or unrestricted (2) for 2 wave functions.
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
    logical, allocatable :: occ1(:, :) !< Occupied MO mask 1.
    logical, allocatable :: occ2(:, :) !< Occupied MO mask 2.
    logical, allocatable :: act1(:, :) !< Active MO mask 1.
    logical, allocatable :: act2(:, :) !< Active MO mask 2.
    real(dp), allocatable :: s_ao(:, :) !< Atomic orbital overlaps.
    real(dp), allocatable :: s_mo(:, :, :) !< Molecular orbital overlaps.
    real(dp), allocatable :: dys_mo(:, :, :) !< Dyson orbitals in mo basis.
    real(dp), allocatable :: dys_ao(:, :, :) !< Dyson orbitals.
    real(dp), allocatable :: dys_norm(:, :) !< Norms of the dyson orbitals.
    type(atombasis), allocatable :: basis(:) !< Basis sets 2.
    integer, allocatable :: bindex(:) !< BS index for each atom 2.
    real(dp), allocatable :: geom(:) !< Geometry.
    character(len=2), allocatable :: atsym(:) !< Atom symbols.
    integer, allocatable :: atnum(:) !< Atom numbers.
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

    call read_geom(input_format, dir2, geom, atsym, atnum)
    call read_ccg_ao(input_format, dir1, ccg1, trans_ao=trans1)
    call read_ccg_ao(input_format, dir2, ccg2, trans_ao=trans2, basis=basis, basis_index=bindex)
    call one_el_op(ccg1, ccg2, [0,0,0], trans1, trans2, s_ao)

    call read_mo(input_format, dir1, moa_c=moa1, mob_c=mob1)
    call read_mo(input_format, dir2, moa_c=moa2, mob_c=mob2)

    call read_cis(input_format, dir1, cisa=wfa1, cisb=wfb1, occ_mo=occ1, act_mo=act1, norm=.true.)
    call read_cis(input_format, dir2, cisa=wfa2, cisb=wfb2, occ_mo=occ2, act_mo=act2, norm=.true.)
    nwf1 = size(wfa1, 3)
    nwf2 = size(wfa2, 3)

    if (allocated(wfb1)) rhf1 = 2
    if (allocated(wfb2)) rhf2 = 2
    rhf = max(rhf1, rhf2)
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
    call sort_mo(occ1(:, 1), act1(:, 1), moa1, remove_inactive = .true.)
    call sort_mo(occ2(:, 1), act2(:, 1), moa2, remove_inactive = .true.)
    if (rhf == 2) call sort_mo(occ1(:, rhf1), act1(:, rhf1), mob1, remove_inactive = .true.)
    if (rhf == 2) call sort_mo(occ2(:, rhf2), act2(:, rhf2), mob2, remove_inactive = .true.)

    allocate(s_mo(size(moa1, 2), size(moa2, 2), rhf))
    call mat_ge_mmm(moa1, s_ao, moa2, s_mo(:, :, 1), transa='T')
    if (rhf == 2) call mat_ge_mmm(mob1, s_ao, mob2, s_mo(:, :, 2), transa='T')

    call cis_dyson(thr, s_mo, wfa1, wfa2, wfb1, wfb2, dys_mo)
    allocate(dys_norm(nwf1+1, nwf2+1))
    dys_norm = sum(dys_mo(:, :, :)*dys_mo(:, :, :), dim=1)

    allocate(wrk(size(trans2, 2), size(mob2, 2)))
    allocate(dys_ao(size(trans2, 2), nwf1, nwf2))
    call gemm(trans2, mob2(:, :), wrk, transa='T')
    do j = 1, nwf2
        call gemm(wrk, dys_mo(:, :, j), dys_ao(:, :, j))
    end do

    open(newunit=outunit, file='dys.at', action='write')
    call write_molden_atoms(outunit, geom, atsym, atnum)
    close(outunit)
    open(newunit=outunit, file='dys.gto', action='write')
    call write_molden_gto(outunit, basis, bindex)
    close(outunit)
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
