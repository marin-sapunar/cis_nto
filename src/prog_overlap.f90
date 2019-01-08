program cis_overlap_prog
    ! General
    use global_defs
    use file_mod
    use matrix_mod
    ! Chem
    use occupation_mod
    use ccg_ao_mod
    use one_el_op_mod
    use cis_overlap_nto_mod
    use cis_overlap_cis_mod
    use orthog_mod
    use phase_mod
    ! I/O
    use read_all_mod
    implicit none

    ! System
    integer :: rhf = 0 !< Restricted (1) or unrestricted(2) calculation.
    integer :: rhf1 = 1 !< Restricted (1) or unrestricted (2) for 1 wave functions.
    integer :: rhf2 = 1 !< Restricted (1) or unrestricted (2) for 2 wave functions.
    integer :: nwf1 !< Number of states 1.
    integer :: nwf2 !< Number of states 2.
    real(dp), allocatable :: geom1(:) !< Geometry 1.
    real(dp), allocatable :: geom2(:) !< Geometry 2.
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
    real(dp), allocatable :: wfa1(:, :) !< CI vector alpha single 1.
    real(dp), allocatable :: wfa2(:, :) !< CI vector alpha single 2.
    real(dp), allocatable :: wfb1(:, :) !< CI vector beta single 1.
    real(dp), allocatable :: wfb2(:, :) !< CI vector beta single 2.
    real(dp), allocatable :: cisa1(:, :, :) !< CIS matrix alpha 1.
    real(dp), allocatable :: cisa2(:, :, :) !< CIS matrix alpha 2.
    real(dp), allocatable :: cisb1(:, :, :) !< CIS matrix beta 1.
    real(dp), allocatable :: cisb2(:, :, :) !< CIS matrix beta 2.
    logical, allocatable :: occ1(:, :) !< Occupied MO mask 1.
    logical, allocatable :: occ2(:, :) !< Occupied MO mask 2.
    logical, allocatable :: act1(:, :) !< Active MO mask 1.
    logical, allocatable :: act2(:, :) !< Active MO mask 2.
    ! Results
    real(dp), allocatable :: s_ao(:, :) !< Atomic orbital overlaps.
    real(dp), allocatable :: s_mo(:, :, :) !< Molecular orbital overlaps.
    real(dp), allocatable :: s_wf(:, :) !< Wave function overlaps.
    real(dp), allocatable :: s_wf_raw(:, :) !< Non-orthogonalized wave function overlaps.
    real(dp) :: angle !< Angle between raw and orthogonalized overlap matrix.
    ! Options
    character(len=:), allocatable :: input_format
    character(len=:), allocatable :: dir1
    character(len=:), allocatable :: dir2
    character(len=:), allocatable :: outfile
    real(dp) :: thr
    logical :: norm
    logical :: orth
    logical :: orth_omat
    logical :: phase_omat
    logical :: center1
    logical :: center2
    integer :: alg
    ! Help
    integer :: i, narg, outunit
    character(len=1000) :: temp
    integer, external :: omp_get_max_threads
    real(dp), external :: omp_get_wtime
    real(dp) :: time00, time0
    real(dp) :: time_io, time_ao, time_mo, time_wf, time_tot

    time00 = omp_get_wtime()

    ! Defaults.
    print_level = 2
    outfile = 'omat'
    input_format = 'turbomole'
    alg = 2
    thr = 1.0_dp
    norm = .true.
    orth = .false.
    orth_omat = .false.
    phase_omat = .false.
    center1 = .false.
    center2 = .false.

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
            write(stdout, *) '  dir1                  directory containing calculation for bra states      '
            write(stdout, *) '  dir2                  directory containing calculation for ket states      '
            write(stdout, *)
            write(stdout, *) 'optional arguments:'
            write(stdout, *) '  --help                show this help message and exit                      '
            write(stdout, *) '  --cis                 use CIS algorithm for calculating overlaps           '
            write(stdout, *) '                        (very slow, not recommended)                         '
            write(stdout, *) '  --nto                 use NTO algorithm for calculating overlaps (default) '
            write(stdout, *) '  --(no-)norm-states    renormalize input states before calculation          '
            write(stdout, '(25x,a,l1)') 'default: ', norm
            write(stdout, *) '  --(no-)orth-states    reorthogonalize input states before calculation      '
            write(stdout, '(25x,a,l1)') 'default: ', orth
            write(stdout, *) '  --(no-)orth-overlap   orthogonalize overlap matrix                         '
            write(stdout, '(25x,a,l1)') 'default: ', orth_omat
            write(stdout, *) '  --(no-)phase-overlap  match phase between assigned bra/ket states          '
            write(stdout, '(25x,a,l1)') 'default: ', phase_omat
            write(stdout, *) '  --(no-)recenter-aos   attempt to remove effect of basis set translation by '
            write(stdout, *) '                        recentering pairs of AOs in AO overlap calculation   '
            write(stdout, *) '                        (untested, not recommended without further testing)  '
            write(stdout, '(25x,a,l1)') 'default: ', center2
            write(stdout, *) '  --threshold t         truncate wave functions using given threshold        '
            write(stdout, *) '  --outfile file        output final overlap matrix to file                  '
            write(stdout, '(25x,a,a)') 'default: ', outfile
            write(stdout, *) '  --print-level p       control output level of program (0 = quiet)           '
            write(stdout, '(25x,a,i0)') 'default: ', print_level
            stop
        case('--cis')
            alg = 1
        case('--nto')
            alg = 2
        case('--norm-states')
            norm = .true.
        case('--no-norm-states')
            norm = .false.
        case('--orth-states')
            orth = .true.
        case('--no-orth-states')
            orth = .false.
        case('--orth-overlap')
            orth_omat = .true.
        case('--no-orth-overlap')
            orth_omat = .false.
        case('--match-phase')
            phase_omat = .true.
        case('--no-match-phase')
            phase_omat = .false.
        case('--recenter-aos')
            center2 = .true.
        case('--no-recenter-aos')
            center2 = .false.
        case('--threshold')
            i = i + 1
            call get_command_argument(i, temp)
            read(temp, *) thr
        case('--outfile')
            i = i + 1
            call get_command_argument(i, temp)
            outfile = trim(adjustl(temp))
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
        write(stdout, '(5x,a)') '                     cis_overlap program                     '
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
        write(stdout, '(5x,a,a)') ' Directory containing calculation for bra states: ', dir1
        write(stdout, '(5x,a,a)') ' Directory containing calculation for ket states: ', dir2
    end if

    if (.not. is_dir(dir1)) then
        write(stderr, *) 
        write(stderr, '(1x,a,a,a)') 'Error. Directory ', dir1, ' not found.'
        stop
    end if
    call read_ccg_ao(input_format, dir1, ccg1, geom=geom1, trans_ao=trans1)
    call read_mo(input_format, dir1, moa_c=moa1, mob_c=mob1)
    call read_cis(input_format, dir1, wfa=wfa1, wfb=wfb1, occ_mo=occ1, act_mo=act1, occ_num=on1, &
    &             norm=norm, orthog=orth)
    nwf1 = size(wfa1, 2)
    if (allocated(wfb1)) rhf1 = 2
    if (print_level >= 1) then
        write(stdout, *) 
        write(stdout, '(1x,a)') 'Properties for bra states:'
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
    call read_ccg_ao(input_format, dir2, ccg2, geom=geom2, trans_ao=trans2)
    call read_mo(input_format, dir2, moa_c=moa2, mob_c=mob2)
    call read_cis(input_format, dir2, wfa=wfa2, wfb=wfb2, occ_mo=occ2, act_mo=act2, norm=norm, &
    &             orthog=orth, occ_num = on2)
    nwf2 = size(wfa2, 2)
    if (allocated(wfb2)) rhf2 = 2
    if (print_level >= 1) then
        write(stdout, *) 
        write(stdout, '(1x,a)') 'Properties for ket states:'
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
    if ((on1%o(1) /= on2%o(1)) .or. (on1%o(rhf1) /= on2%o(rhf2))) then
        write(stderr, *) 
        write(stderr, '(1x,a,a,a)') 'Error. Mismatch in number of occupied orbitals.'
        stop
    end if
    if ((on1%ao(1) /= on2%ao(1)) .or. (on1%ao(rhf1) /= on2%ao(rhf2))) then
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
    call one_el_op(ccg1, ccg2, [0,0,0], geom1, geom2, trans1, trans2, s_ao, &
    &              center_atom_pairs=center2, center_diagonal_block=center1)
    time_ao = omp_get_wtime() - time0

    ! Allocate arrays for mixed restricted/unrestricted calculation.
    rhf = max(rhf1, rhf2)
    if (rhf1 /= rhf) then
        if (print_level >= 1) then
            write(stdout, *) 
            write(stdout, '(1x,a)') 'Unrestricted bra and restricted ket calculation...'
            write(stdout, '(5x,a)') 'Renormalizing bra CI coefficients...'
        end if
        wfa1 = wfa1 / sqrt(2.0_dp)
        allocate(mob1, source=moa1)
        allocate(wfb1, source=wfa1)
        wfb1 = -wfb1
    else if (rhf2 /= rhf) then
        if (print_level >= 1) then
            write(stdout, *) 
            write(stdout, '(1x,a)') 'Restricted bra and unrestricted ket calculation...'
            write(stdout, '(5x,a)') 'Renormalizing ket CI coefficients...'
        end if
        wfa2 = wfa2 / sqrt(2.0_dp)
        allocate(mob2, source=moa2)
        allocate(wfb2, source=wfa2)
        wfb2 = -wfb2
    end if

    ! Remove frozen mos and calculate MO overlap matrix.
    if (print_level >= 1) then
        if ((.not. all(act1)) .or. (.not. all(act2))) then
            write(stdout, *) 
            write(stdout, '(1x,a)') 'Removing frozen MOs...'
        end if
    end if
    time0 = omp_get_wtime()
    call sort_mo(occ1(:, 1), act1(:, 1), moa1, remove_inactive=.true.)
    call sort_mo(occ2(:, 1), act2(:, 1), moa2, remove_inactive=.true.)
    if (rhf == 2) call sort_mo(occ1(:, rhf1), act1(:, rhf1), mob1, remove_inactive=.true.)
    if (rhf == 2) call sort_mo(occ2(:, rhf2), act2(:, rhf2), mob2, remove_inactive=.true.)
    if (print_level >= 2) then
        write(stdout, *) 
        write(stdout, '(1x,a)') 'Computing MO overlaps..'
    end if
    allocate(s_mo(size(moa1, 2), size(moa2, 2), rhf))
    call mat_ge_mmm(moa1, s_ao, moa2, s_mo(:, :, 1), transa='T')
    if (rhf == 2) call mat_ge_mmm(mob1, s_ao, mob2, s_mo(:, :, 2), transa='T')
    time_mo = omp_get_wtime() - time0

    ! Calculate WF overlap matrix.
    if (print_level >= 2) then
        write(stdout, *) 
        write(stdout, '(1x,a)') 'Computing WF overlaps..'
    end if
    time0 = omp_get_wtime()
    select case(alg)
    case(1)
        call cis_overlap_cis(thr, on1%ao(1), on1%ao(rhf1), s_mo, wfa1, wfa2, wfb1, wfb2, s_wf)
    case(2)
        cisa1 = reshape(wfa1, [on1%av(1), on1%ao(1), nwf1])
        deallocate(wfa1)
        cisa2 = reshape(wfa2, [on2%av(1), on2%ao(1), nwf2])
        deallocate(wfa2)
        if (rhf == 2) then
            cisb1 = reshape(wfb1, [on1%av(2), on1%ao(2), nwf1])
            deallocate(wfb1)
            cisb2 = reshape(wfb2, [on2%av(2), on2%ao(2), nwf2])
            deallocate(wfb2)
        end if
        call cis_overlap_nto(thr, s_mo, cisa1, cisa2, cisb1, cisb2, s_wf)
        allocate(s_wf_raw, source=s_wf)
    end select
    time_wf = omp_get_wtime() - time0

    ! Output.
    if (print_level >= 1) then
        write(stdout, *) 
        write(stdout, '(1x,a)') 'Raw overlap matrix:'
        do i = 0, nwf2
            write(stdout,'(5x,1000es24.16)') s_wf(i, :)
        end do
    end if

    if (orth_omat) then
        call orthog_lowdin(s_wf)
        if (print_level >= 1) then
            write(stdout, *) 
            write(stdout, '(1x,a)') 'Orthogonalized overlap matrix:'
            do i = 0, nwf2
                write(stdout,'(5x,1000es24.16)') s_wf(i, :)
            end do
            angle = acos(sum(s_wf*s_wf_raw) / mat_norm(s_wf) / mat_norm(s_wf_raw))
            write(stdout, *) 
            write(stdout, '(5x,a,f6.2)') 'Frobenius inner product angle: ', angle * 180 / 3.14159
            write(stdout, '(5x,a)') 'Norms of rows of raw overlap matrix: '
            write(stdout, '(9x,1000f10.4)') sum(s_wf_raw**2, 2)
            write(stdout, '(5x,a)') 'Norms of columns of raw overlap matrix: '
            write(stdout, '(9x,1000f10.4)') sum(s_wf_raw**2, 1)
        end if
    end if

    if (phase_omat) then
        call phasematch_assigned(s_wf)
        if (print_level >= 1) then
            write(stdout, *) 
            write(stdout, '(1x,a)') 'Overlap matrix with phase matching between assigned bra/ket states:'
            do i = 0, nwf2
                write(stdout, '(5x,1000es24.16)') s_wf(i, :)
            end do
        end if
    end if

    open(newunit=outunit, file=outfile, action='write')
    do i = 0, nwf2
        write(outunit, '(1000es24.16)') s_wf(i, :)
    end do
    close(outunit)

    time_tot = omp_get_wtime() - time00
    if (print_level >= 2) then
        write(stdout, *) 
        write(stdout,'(1x,a)') 'Program time:'
        write(stdout, '(5x, a40, f14.4)') 'Input                     - time (sec):', time_io
        write(stdout, '(5x, a40, f14.4)') 'AO overlap                - time (sec):', time_ao
        write(stdout, '(5x, a40, f14.4)') 'MO overlap                - time (sec):', time_mo
        write(stdout, '(5x, a40, f14.4)') 'WF overlap                - time (sec):', time_wf
        write(stdout, '(5x, 40x, a14)') '--------------'
        write(stdout, '(5x, a40, f14.4)') 'Total                     - time (sec):', time_tot
    end if

    if (print_level >= 1) then
        write(stdout, *) 
        write(stdout, '(1x,a)') 'cis_overlap done                                             '
        write(stdout, '(1x,a)') '-------------------------------------------------------------'
    end if

end program cis_overlap_prog
