program cis_olap_test
    use global_defs
    use cis_overlap_mod
    use read_all_mod
    use ccg_ao_mod
    use one_el_op_mod
    use matrix_mod
    use occupation_mod
    implicit none

    integer :: rhf = 0 !< Restricted (1) or unrestricted(2) calculation.
    integer :: rhf1 = 0 !< Restricted (1) or unrestricted (2) for 1 wave functions.
    integer :: rhf2 = 0 !< Restricted (1) or unrestricted (2) for 2 wave functions.
    integer :: nwf1 = 0 !< Number of wave functions 1.
    integer :: nwf2 = 0 !< Number of wave functions 2.
    type(ccg), allocatable :: ccg1(:) !< Atomic orbitals 1.
    type(ccg), allocatable :: ccg2(:) !< Atomic orbitals 2.
    real(dp), allocatable :: trans1(:, :) !< Transformation matrix for basis functions 1.
    real(dp), allocatable :: trans2(:, :) !< Transformation matrix for basis functions 2.
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
    real(dp), allocatable :: s_wf(:, :) !< Wave function overlaps.
    character(len=1000) :: temp
    character(len=:), allocatable :: input_format
    character(len=:), allocatable :: dir1
    character(len=:), allocatable :: dir2
    real(dp) :: thr
    integer :: i
    integer :: ounit

    input_format = 'turbomole'
    write(stdout, *) 'Path to input for bra states:'
    read(stdin, '(a)') temp
    dir1 = trim(adjustl(temp))
    write(stdout, *) 'Path to input for ket states:'
    read(stdin, '(a)') temp
    dir2 = trim(adjustl(temp))
    write(stdout, *) 'Threshold for truncating wave functions:'
    read(stdin, *) thr

    call read_ccg_ao(input_format, dir1, ccg1, trans_ao=trans1)
    call read_ccg_ao(input_format, dir2, ccg2, trans_ao=trans2)
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
    if (rhf == 2) call sort_mo(occ1(:, 2), act1(:, 2), mob1, remove_inactive = .true.)
    if (rhf == 2) call sort_mo(occ2(:, 2), act2(:, 2), mob2, remove_inactive = .true.)

    allocate(s_mo(size(moa1, 2), size(moa2, 2), rhf))
    call mat_ge_mmm(moa1, s_ao, moa2, s_mo(:, :, 1), transa='T')
    if (rhf == 2) call mat_ge_mmm(mob1, s_ao, mob2, s_mo(:, :, 2), transa='T')

    call cis_overlap(thr, s_mo, wfa1, wfa2, wfb1, wfb2, s_wf)

    open(newunit=ounit, file='omat', action='write')
    do i = 1, nwf2 + 1
        write(ounit,*) s_wf(:, i)
    end do
    close(ounit)

end program cis_olap_test
