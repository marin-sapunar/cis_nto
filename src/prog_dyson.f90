program cis_dyson_prog
    ! General
    use global_defs
    use file_mod
    use matrix_mod
    ! Chem
    use atom_basis_mod
    use occupation_mod
    use ccg_ao_mod
    use one_el_op_mod
    use cis_dyson_mod
    ! I/O
    use read_all_mod
    use write_molden_mod
    ! External
    use blas95, only : gemm
    implicit none

    ! System
    integer :: rhf = 0 !< Restricted (1) or unrestricted(2) calculation.
    integer :: rhf1 = 1 !< Restricted (1) or unrestricted (2) for 1 wave functions.
    integer :: rhf2 = 1 !< Restricted (1) or unrestricted (2) for 2 wave functions.
    integer :: nwf1 !< Number of excited states 1.
    integer :: nwf2 !< Number of excited states 2.
    real(dp), allocatable :: geom1(:) !< Geometry 1.
    real(dp), allocatable :: geom2(:) !< Geometry 2.
    character(len=2), allocatable :: atsym(:) !< Atom symbols.
    integer, allocatable :: atnum(:) !< Atom numbers.
    type(atombasis), allocatable :: basis(:) !< Basis sets.
    integer, allocatable :: bindex(:) !< Basis set index for each atom.
    type(occupation_numbers) :: on1 !< Occupation numbers 1.
    type(occupation_numbers) :: on2 !< Occupation numbers 2.
    type(ccg), allocatable :: ccg1(:) !< Atomic orbitals 1.
    type(ccg), allocatable :: ccg2(:) !< Atomic orbitals 2.
    real(dp), allocatable :: trans1(:, :) !< Transformation matrix for basis functions 1.
    real(dp), allocatable :: trans2(:, :) !< Transformation matrix for basis functions 2.
    real(dp), allocatable :: moa1(:, :) !< Molecular orbital coefficients alpha 1.
    real(dp), allocatable :: moa2(:, :) !< Molecular orbital coefficients alpha 2.
    real(dp), allocatable :: mob1(:, :) !< Molecular orbital coefficients beta 1.
    real(dp), allocatable :: mob2(:, :) !< Molecular orbital coefficients beta 2.
    real(dp), allocatable :: cisa1(:, :, :) !< CI alpha single excitation coefficients 1.
    real(dp), allocatable :: cisa2(:, :, :) !< CI alpha single excitation coefficients 2.
    real(dp), allocatable :: cisb1(:, :, :) !< CI beta single excitation coefficients 1.
    real(dp), allocatable :: cisb2(:, :, :) !< CI beta single excitation coefficients 2.
    logical, allocatable :: occ1(:, :) !< Occupied MO mask 1.
    logical, allocatable :: occ2(:, :) !< Occupied MO mask 2.
    logical, allocatable :: act1(:, :) !< Active MO mask 1.
    logical, allocatable :: act2(:, :) !< Active MO mask 2.
    ! Results
    real(dp), allocatable :: s_ao(:, :) !< Atomic orbital overlaps.
    real(dp), allocatable :: s_mo(:, :, :) !< Molecular orbital overlaps.
    real(dp), allocatable :: dys_mo(:, :, :) !< Dyson orbitals in mo basis.
    real(dp), allocatable :: dys_ao(:, :, :) !< Dyson orbitals.
    real(dp), allocatable :: dys_norm(:, :) !< Norms of the dyson orbitals.
    ! Options
    character(len=:), allocatable :: input_format
    character(len=:), allocatable :: dir1
    character(len=:), allocatable :: dir2
    character(len=:), allocatable :: prefix
    real(dp) :: thr
    logical :: norm
    logical :: orth
    ! Help
    integer :: i, j, narg, outunit
    character(len=1000) :: temp
    real(dp), allocatable :: wrk(:, :)
    integer, external :: omp_get_max_threads
    real(dp), external :: omp_get_wtime
    real(dp) :: time00, time0
    real(dp) :: time_io, time_ao, time_mo, time_wf, time_tot

    time00 = omp_get_wtime()

    ! Defaults.
    print_level = 2
    prefix = 'dys'
    input_format = 'turbomole'
    thr = 1.0_dp
    norm = .true.
    orth = .false.

    ! CLI
    narg = command_argument_count()
    i = 0
    do while (i < narg)
        i = i + 1
        call get_command_argument(i, temp)
        if (i == narg - 1) exit
        if (i == narg) temp = '--help'
        select case(temp)
        case('--help')
            write(stdout, *) 'usage: cis_overlap.exe [optional arguments] dir1 dir2'
            write(stdout, *)
            write(stdout, *) 'Calculate overlaps between two sets CIS type wave functions.'
            write(stdout, *)
            write(stdout, *) 'positional arguments:'
            write(stdout, *) '  dir1                  directory containing calculation for N-1 el. system    '
            write(stdout, *) '  dir2                  directory containing calculation for N el. system      '
            write(stdout, *)
            write(stdout, *) 'optional arguments:'
            write(stdout, *) '  --help                show this help message and exit                      '
            write(stdout, *) '  --(no-)norm-states    renormalize input states before calculation          '
            write(stdout, '(25x,a,l1)') 'default: ', norm
            write(stdout, *) '  --(no-)orth-states    reorthogonalize input states before calculation      '
            write(stdout, '(25x,a,l1)') 'default: ', orth
            write(stdout, *) '  --threshold t         truncate wave functions using given threshold        '
            write(stdout, *) '  --prefix pref         prefix for output files                              '
            write(stdout, *) '                        (geometry written to pref.at, basis set to pref.gto, '
            write(stdout, *) '                         and dyson orbitals to pref.#.#.ao and pref.#.#.mo)  '
            write(stdout, '(25x,a,a)') 'default: ', prefix
            write(stdout, *) '  --print-level p       control output level of program (0 = quiet)           '
            write(stdout, '(25x,a,i0)') 'default: ', print_level
            stop
        case('--norm-states')
            norm = .true.
        case('--no-norm-states')
            norm = .false.
        case('--orth-states')
            orth = .true.
        case('--no-orth-states')
            orth = .false.
        case('--threshold')
            i = i + 1
            call get_command_argument(i, temp)
            read(temp, *) thr
        case('--prefix')
            i = i + 1
            call get_command_argument(i, temp)
            prefix = trim(adjustl(temp))
        case('--print-level')
            i = i + 1
            call get_command_argument(i, temp)
            read(temp, *) print_level
        end select
    end do
    call get_command_argument(narg-1, temp)
    dir1 = trim(adjustl(temp))
    call get_command_argument(narg, temp)
    dir2 = trim(adjustl(temp))

    ! Begin run.
    if (print_level >= 1) then
        write(stdout, '(5x,a)') '-------------------------------------------------------------'
        write(stdout, '(5x,a)') '                      cis_dyson program                      '
        write(stdout, '(5x,a)') '                         version 1.0                         '
        write(stdout, '(5x,a)') '                                                             '
        write(stdout, '(5x,a)') ' Program compiled on '//__DATE__//' '//__TIME__//'.          '
        write(stdout, '(5x,a)') '-------------------------------------------------------------'
        write(stdout, *)
        write(stdout, '(1x,a,i0,a)') 'Using ', omp_get_max_threads(), ' threads.'
    end if

    ! Read input.
    time0 = omp_get_wtime()
    if (print_level >= 1) then
        write(stdout, *)
        write(stdout, '(1x,a)') 'Reading input from turbomole files...'
        write(stdout, '(5x,a,a)') ' Directory containing calculation for N-1 el. states: ', dir1
        write(stdout, '(5x,a,a)') ' Directory containing calculation for N el. states: ', dir2
    end if

    if (.not. is_dir(dir1)) then
        write(stderr, *)
        write(stderr, '(1x,a,a,a)') 'Error. Directory ', dir1, ' not found.'
        stop
    end if
    call read_ccg_ao(input_format, dir1, ccg1, geom=geom1, trans_ao=trans1)
    call read_mo(input_format, dir1, moa_c=moa1, mob_c=mob1)
    call read_cis(input_format, dir1, cisa=cisa1, cisb=cisb1, occ_mo=occ1, act_mo=act1, occ_num=on1, &
    &             norm=norm, orthog=orth)
    nwf1 = size(cisa1, 3)
    if (allocated(cisb1)) rhf1 = 2
    if (print_level >= 1) then
        write(stdout, *)
        write(stdout, '(1x,a)') 'Properties for N-1 el. states:'
        if (rhf1 == 1) write(stdout, '(5x, a)') 'Restricted calculation.'
        if (rhf1 == 2) write(stdout, '(5x, a)') 'Unrestricted calculation.'
        write(stdout, '(5x,a,2(1x,i0))') 'Number of excited states:           ', nwf1
        write(stdout, '(5x,a,2(1x,i0))') 'Number of basis functions (sphe):   ', size(trans1, 1)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of basis functions (cart):   ', size(trans1, 2)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of orbitals:                 ', on2%n
        write(stdout, '(5x,a,2(1x,i0))') 'Number of occupied orbitals:        ', on1%o(1:rhf1)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of virtual orbitals:         ', on1%v(1:rhf1)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of active occupied orbitals: ', on1%ao(1:rhf1)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of active virtual orbitals:  ', on1%av(1:rhf1)
    end if

    if (.not. is_dir(dir2)) then
        write(stderr, *)
        write(stderr, '(1x,a,a,a)') 'Error. Directory ', dir2, ' not found.'
        stop
    end if
    call read_geom(input_format, dir2, geom2, atsym, atnum)
    call read_ccg_ao(input_format, dir2, ccg2, geom=geom2, trans_ao=trans2, basis=basis, basis_index=bindex)
    call read_mo(input_format, dir2, moa_c=moa2, mob_c=mob2)
    call read_cis(input_format, dir2, cisa=cisa2, cisb=cisb2, occ_mo=occ2, act_mo=act2, occ_num=on2, &
    &             norm=norm, orthog=orth)
    nwf2 = size(cisa2, 3)
    if (allocated(cisb2)) rhf2 = 2
    if (print_level >= 1) then
        write(stdout, *)
        write(stdout, '(1x,a)') 'Properties for N el. states:'
        if (rhf2 == 1) write(stdout, '(5x, a)') 'Restricted calculation.'
        if (rhf2 == 2) write(stdout, '(5x, a)') 'Unrestricted calculation.'
        write(stdout, '(5x,a,2(1x,i0))') 'Number of excited states:           ', nwf2
        write(stdout, '(5x,a,2(1x,i0))') 'Number of basis functions (sphe):   ', size(trans2, 1)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of basis functions (cart):   ', size(trans2, 2)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of orbitals:                 ', on2%n
        write(stdout, '(5x,a,2(1x,i0))') 'Number of occupied orbitals:        ', on2%o(1:rhf2)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of virtual orbitals:         ', on2%v(1:rhf2)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of active occupied orbitals: ', on2%ao(1:rhf2)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of active virtual orbitals:  ', on2%av(1:rhf2)
    end if
    time_io = omp_get_wtime() - time0
    if ((on1%o(1) /= on2%o(1)) .or. (on1%o(rhf1) /= on2%o(rhf2)-1)) then
        write(stderr, *)
        write(stderr, '(1x,a,a,a)') 'Error. Mismatch in number of occupied orbitals.'
        stop
    end if
    if ((on1%ao(1) /= on2%ao(1)) .or. (on1%ao(rhf1) /= on2%ao(rhf2)-1)) then
        write(stderr, *)
        write(stderr, '(1x,a,a,a)') 'Error. Mismatch in number of active occupied orbitals.'
        stop
    end if

    ! Calculate AO overlaps.
    if (print_level >= 2) then
        write(stdout, *)
        write(stdout, '(1x,a)') 'Computing AO overlaps...'
    end if
    time0 = omp_get_wtime()
    call one_el_op(ccg1, ccg2, [0,0,0], geom1, geom2, trans1, trans2, s_ao)
    time_ao = omp_get_wtime() - time0

    ! Allocate arrays for mixed restricted/unrestricted calculation.
    rhf = max(rhf1, rhf2)
    if (rhf1 /= rhf) then
        if (print_level >= 1) then
            write(stdout, *)
            write(stdout, '(1x,a)') 'Unrestricted N-1 el. and restricted N el. calculation...'
            write(stdout, '(5x,a)') 'Renormalizing N-1 el. CI coefficients...'
        end if
        cisa1 = cisa1 / sqrt(2.0_dp)
        allocate(mob1, source=moa1)
        allocate(cisb1, source=cisa1)
        cisb1 = -cisb1
    else if (rhf2 /= rhf) then
        if (print_level >= 1) then
            write(stdout, *)
            write(stdout, '(1x,a)') 'Restricted N-1 el. and unrestricted N el. calculation...'
            write(stdout, '(5x,a)') 'Renormalizing N el. CI coefficients...'
        end if
        cisa2 = cisa2 / sqrt(2.0_dp)
        allocate(mob2, source=moa2)
        allocate(cisb2, source=cisa2)
        cisb2 = -cisb2
    end if

    ! Remove frozen mos and calculate MO overlap matrix.
    if (print_level >= 1) then
        if ((.not. all(act1)) .or. (.not. all(act2))) then
            write(stdout, *)
            write(stdout, '(1x,a)') 'Removing frozen MOs...'
        end if
    end if
    time0 = omp_get_wtime()
    call sort_mo(occ1(:, 1), act1(:, 1), moa1, remove_inactive = .true.)
    call sort_mo(occ2(:, 1), act2(:, 1), moa2, remove_inactive = .true.)
    if (rhf == 2) call sort_mo(occ1(:, rhf1), act1(:, rhf1), mob1, remove_inactive = .true.)
    if (rhf == 2) call sort_mo(occ2(:, rhf2), act2(:, rhf2), mob2, remove_inactive = .true.)
    if (print_level >= 2) then
        write(stdout, *)
        write(stdout, '(1x,a)') 'Computing MO overlaps..'
    end if
    allocate(s_mo(size(moa1, 2), size(moa2, 2), rhf))
    call mat_ge_mmm(moa1, s_ao, moa2, s_mo(:, :, 1), transa='T')
    if (rhf == 2) call mat_ge_mmm(mob1, s_ao, mob2, s_mo(:, :, 2), transa='T')
    time_mo = omp_get_wtime() - time0

    ! Calculate dyson orbitals
    time0 = omp_get_wtime()
    if (print_level >= 2) then
        write(stdout, *)
        write(stdout, '(1x,a)') 'Computing Dyson orbitals..'
    end if
    call cis_dyson(thr, s_mo, cisa1, cisa2, cisb1, cisb2, dys_mo)
    allocate(dys_norm(0:nwf1, 0:nwf2))
    dys_norm = sum(dys_mo(:, :, :)*dys_mo(:, :, :), dim=1)
    if (print_level >= 2) then
        write(stdout, '(5x,a)') 'Transforming Dyson orbitals to AO basis..'
    end if
    allocate(wrk(size(trans2, 2), size(mob2, 2)))
    allocate(dys_ao(size(trans2, 2), 0:nwf1, 0:nwf2))
    call gemm(trans2, mob2(:, :), wrk, transa='T')
    do j = 0, nwf2
        call gemm(wrk, dys_mo(:, :, j), dys_ao(:, :, j))
    end do
    time_wf = omp_get_wtime() - time0

    ! Output
    open(newunit=outunit, file=prefix//'.at', action='write')
    call write_molden_atoms(outunit, geom1, atsym, atnum)
    close(outunit)
    open(newunit=outunit, file=prefix//'.gto', action='write')
    call write_molden_gto(outunit, basis, bindex)
    close(outunit)
    do i = 0, nwf1
        do j = 0, nwf2
            ! Write MO basis dyson orbital.
            write(temp, '(a,a,i3.3,a,i3.3,a)') prefix, '.', i, '.', j, '.mo'
            open(newunit=outunit, file=temp, action='write')
            write(outunit, '(1e13.5)') dys_mo(:, i, j)
            close(outunit)
            ! Write AO basis dyson orbital.
            write(temp, '(a,a,i3.3,a,i3.3,a)') prefix, '.', i, '.', j, '.ao'
            open(newunit=outunit, file=temp, action='write')
            call write_molden_mo_single(outunit, 1000*i+j, 'a   ', 1, 0.0_dp, dys_norm(i, j),      &
            &                           dys_ao(:, i, j))
            close(outunit)
        end do
    end do

    time_tot = omp_get_wtime() - time00
    if (print_level >= 2) then
        write(stdout, *)
        write(stdout,'(1x,a)') 'Program time:'
        write(stdout, '(5x, a40, f14.4)') 'Input                     - time (sec):', time_io
        write(stdout, '(5x, a40, f14.4)') 'AO overlap                - time (sec):', time_ao
        write(stdout, '(5x, a40, f14.4)') 'MO overlap                - time (sec):', time_mo
        write(stdout, '(5x, a40, f14.4)') 'Dyson orbitals            - time (sec):', time_wf
        write(stdout, '(5x, 40x, a14)') '--------------'
        write(stdout, '(5x, a40, f14.4)') 'Total                     - time (sec):', time_tot
    end if

    if (print_level >= 1) then
        write(stdout, *)
        write(stdout, '(1x,a)') 'cis_dyson done                                               '
        write(stdout, '(1x,a)') '-------------------------------------------------------------'
    end if


end program cis_dyson_prog
