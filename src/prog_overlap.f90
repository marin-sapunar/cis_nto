!--------------------------------------------------------------------------------------------------
! PROGRAM: cis_overlap_prog
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @brief Calculate overlap between two sets of CIS type wave functions.
!--------------------------------------------------------------------------------------------------
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
    use cis_overlap_l2m_mod
    use orthog_mod
    use phase_mod
    use cis_util_mod
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
    real(dp), allocatable :: cisa1(:, :, :) !< CIS matrix alpha 1.
    real(dp), allocatable :: cisa2(:, :, :) !< CIS matrix alpha 2.
    real(dp), allocatable :: cisb1(:, :, :) !< CIS matrix beta 1.
    real(dp), allocatable :: cisb2(:, :, :) !< CIS matrix beta 2.
    logical, allocatable :: occ1(:, :) !< Occupied MO mask 1.
    logical, allocatable :: occ2(:, :) !< Occupied MO mask 2.
    logical, allocatable :: act1(:, :) !< Active MO mask from el. structure calculation 1.
    logical, allocatable :: act2(:, :) !< Active MO mask from el. structure calculation 2.
    ! Results
    real(dp), allocatable :: s_ao(:, :) !< Atomic orbital overlaps.
    real(dp), allocatable :: s_mo_a(:, :) !< Molecular orbital overlaps alpha.
    real(dp), allocatable :: s_mo_b(:, :) !< Molecular orbital overlaps beta.
    real(dp), allocatable :: s_wf(:, :) !< Wave function overlaps.
    real(dp), allocatable :: s_wf_raw(:, :) !< Non-orthogonalized wave function overlaps.
    real(dp) :: angle !< Angle between raw and orthogonalized overlap matrix.
    ! Options
    character(len=:), allocatable :: out_fmt_s
    character(len=:), allocatable :: input_format
    character(len=:), allocatable :: dir1
    character(len=:), allocatable :: dir2
    character(len=:), allocatable :: outfile
    real(dp) :: thr
    real(dp) :: f_by_mo_norm_t
    logical :: norm
    logical :: orth
    logical :: orth_omat
    logical :: phase_omat
    logical :: center1
    logical :: center2
    logical :: f_by_mo_norm
    character(len=4) :: overlap_alg
    ! Help
    integer :: i, narg, outunit
    character(len=1000) :: temp
    integer, external :: omp_get_max_threads
    real(dp), external :: omp_get_wtime
    real(dp) :: time00, time0
    real(dp) :: time_io, time_ao, time_mo, time_wf, time_tot

    time00 = omp_get_wtime()

    ! Defaults.
    out_fmt_s = '(5x,1000f10.6)'
    print_level = 2
    outfile = 'omat'
    input_format = 'turbomole'
    overlap_alg = 'NTO'
    thr = 1.0_dp
    norm = .true.
    orth = .false.
    orth_omat = .false.
    phase_omat = .false.
    center1 = .false.
    center2 = .false.
    f_by_mo_norm = .false.
    f_by_mo_norm_t = -1.0_dp


    ! CLI
    narg = command_argument_count()
    i = 0
    if (narg == 0) call print_help()
    do while (i < narg)
        i = i + 1
        if (i == narg) call print_help()
        call get_command_argument(i, temp)
        if (i == narg - 1) exit
        select case(temp)
        case('--help', '-h')
            call print_help()
        case('--algorithm', '-alg')
            i = i + 1
            call get_command_argument(i, temp)
            overlap_alg = temp
        case('--norm-states', '-ns')
            norm = .true.
        case('--no-norm-states', '-nns')
            norm = .false.
        case('--orth-states', '-os')
            orth = .true.
        case('--no-orth-states', '-nos')
            orth = .false.
        case('--orth-overlap', '-oo')
            orth_omat = .true.
        case('--no-orth-overlap', '-noo')
            orth_omat = .false.
        case('--match-phase', '-mp')
            phase_omat = .true.
        case('--no-match-phase', '-nmp')
            phase_omat = .false.
        case('--recenter-aos', '-rao')
            center2 = .true.
        case('--no-recenter-aos', '-nrao')
            center2 = .false.
        case('--freeze-mo-norm', '-fmn')
            i = i + 1
            call get_command_argument(i, temp)
            read(temp, *) f_by_mo_norm_t
            if (f_by_mo_norm_t > 0.0_dp) f_by_mo_norm = .true.
        case('--threshold', '-t')
            i = i + 1
            call get_command_argument(i, temp)
            read(temp, *) thr
        case('--outfile', '-o')
            i = i + 1
            call get_command_argument(i, temp)
            outfile = trim(adjustl(temp))
        case('--print-level', '-p')
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

    call need_dir(dir1)
    call read_ccg_ao(input_format, dir1, ccg1, geom=geom1, trans_ao=trans1)
    call read_mo(input_format, dir1, moa_c=moa1, mob_c=mob1)
    call read_cis(input_format, dir1, cisa=cisa1, cisb=cisb1, occ_mo=occ1, act_mo=act1, occ_num=on1, &
    &             norm=norm, orthog=orth)
    nwf1 = size(cisa1, 3)
    if (allocated(cisb1)) rhf1 = 2
    call print_calculation_properties('Properties for bra states:', rhf1, nwf1, size(trans1, 1), &
    &                                 size(trans1, 2), on1)

    call need_dir(dir2)
    call read_ccg_ao(input_format, dir2, ccg2, geom=geom2, trans_ao=trans2)
    call read_mo(input_format, dir2, moa_c=moa2, mob_c=mob2)
    call read_cis(input_format, dir2, cisa=cisa2, cisb=cisb2, occ_mo=occ2, act_mo=act2, norm=norm, &
    &             orthog=orth, occ_num = on2)
    nwf2 = size(cisa2, 3)
    if (allocated(cisb2)) rhf2 = 2
    call print_calculation_properties('Properties for ket states:', rhf2, nwf2, size(trans2, 1), &
    &                                 size(trans2, 2), on2)

    ! Check input dimensions.
    rhf = max(rhf1, rhf2)
    call check_occ_orb(0, rhf1, rhf2, on1, on2, occ1, occ2, act1, act2, cisa1, cisa2, cisb1, cisb2)
    call check_rhf(rhf1, rhf2, moa1, moa2, mob1, mob2, cisa1, cisa2, cisb1, cisb2)
    call remove_pre_frozen_mo(rhf1, rhf2, occ1, occ2, act1, act2, moa1, moa2, mob1, mob2)
    time_io = omp_get_wtime() - time0

    ! Calculate AO overlaps.
    time0 = omp_get_wtime()
    if (print_level >= 2) then
        write(stdout, *) 
        write(stdout, '(1x,a)') 'Computing AO overlaps...'
    end if
    call one_el_op(ccg1, ccg2, [0,0,0], geom1, geom2, trans1, trans2, s_ao, &
    &              center_atom_pairs=center2, center_diagonal_block=center1)
    time_ao = omp_get_wtime() - time0

    ! Calculate MO overlap matrix.
    time0 = omp_get_wtime()
    if (print_level >= 2) then
        write(stdout, *) 
        write(stdout, '(1x,a)') 'Computing MO overlaps..'
    end if
    allocate(s_mo_a(size(moa1, 2), size(moa2, 2)))
    call mat_ge_mmm(moa1, s_ao, moa2, s_mo_a, transa='T')
    if (rhf == 2) then
        allocate(s_mo_b(size(mob1, 2), size(mob2, 2)))
        call mat_ge_mmm(mob1, s_ao, mob2, s_mo_b, transa='T')
    end if
    time_mo = omp_get_wtime() - time0

    ! Remove MO and CI coefficients based on freeze options.
    call check_mo_norms(s_mo_a, moa1, moa2, cisa1, cisa2, f_by_mo_norm, f_by_mo_norm_t)
    if (rhf == 2) then
        call check_mo_norms(s_mo_b, mob1, mob2, cisb1, cisb2, f_by_mo_norm, f_by_mo_norm_t)
    end if
    call check_ci_norms(rhf, cisa1, cisb1, 'bra')
    call check_ci_norms(rhf, cisa2, cisb2, 'ket')

    ! Calculate WF overlap matrix.
    if (print_level >= 2) then
        write(stdout, *) 
        write(stdout, '(1x,a)') 'Computing WF overlaps..'
        write(stdout, *) 
        write(stdout, '(5x,a,a,a)') 'Using ', trim(overlap_alg), ' algorithm.'
    end if
    time0 = omp_get_wtime()
    select case(overlap_alg)
    case('CIS')
        call cis_overlap_cis(rhf, thr, s_mo_a, s_mo_b, cisa1, cisa2, cisb1, cisb2, s_wf)
    case('L2M')
        call cis_overlap_l2m(rhf, s_mo_a, s_mo_b, cisa1, cisa2, cisb1, cisb2, s_wf)
    case('NTO')
        call cis_overlap_nto(rhf, thr, s_mo_a, s_mo_b, cisa1, cisa2, cisb1, cisb2, s_wf)
    end select
    allocate(s_wf_raw, source=s_wf)
    time_wf = omp_get_wtime() - time0

    ! Output.
    if (print_level >= 1) then
        write(stdout, *) 
        write(stdout, '(1x,a)') 'Raw overlap matrix:'
        do i = 0, nwf2
            write(stdout, out_fmt_s) s_wf(i, :)
        end do
    end if

    if (orth_omat) then
        call orthog_lowdin(s_wf)
        if (print_level >= 1) then
            write(stdout, *) 
            write(stdout, '(1x,a)') 'Orthogonalized overlap matrix:'
            do i = 0, nwf2
                write(stdout, out_fmt_s) s_wf(i, :)
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
                write(stdout, out_fmt_s) s_wf(i, :)
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


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: print_help
    !> @brief Print help message and exit program.
    !----------------------------------------------------------------------------------------------
    subroutine print_help()
        write(stdout, '(a)') 'usage: cis_overlap.exe [optional arguments] dir1 dir2'
        write(stdout, '(a)')
        write(stdout, '(a)') 'Calculate overlaps between two sets CIS type wave functions.'
        write(stdout, '(a)')
        write(stdout, '(a)') 'positional arguments:'
        write(stdout, '(a)') '  dir1                        directory containing calculation for bra states      '
        write(stdout, '(a)') '  dir2                        directory containing calculation for ket states      '
        write(stdout, '(a)')
        write(stdout, '(a)') 'optional arguments:'
        write(stdout, '(a)') '  -h, --help                  show this help message and exit                      '
        write(stdout, '(a)') '  -alg, --algorithm ALG       algorithm to use for the overlap calculation         '
        write(stdout, '(a)') '                              Available options:                                   '
        write(stdout, '(a)') '                                CIS (very slow, not recommended)                   '
        write(stdout, '(a)') '                                L2M                                                '
        write(stdout, '(a)') '                                NTO                                                '
        write(stdout, '(32x,a,a)') 'default: ', overlap_alg
        write(stdout, '(a)') '  -ns, --(no-)norm-states     renormalize input states before calculation          '
        write(stdout, '(32x,a,l1)') 'default: ', norm
        write(stdout, '(a)') '  -os, --(no-)orth-states     reorthogonalize input states before calculation      '
        write(stdout, '(32x,a,l1)') 'default: ', orth
        write(stdout, '(a)') '  -oo, --(no-)orth-overlap    orthogonalize overlap matrix                         '
        write(stdout, '(32x,a,l1)') 'default: ', orth_omat
        write(stdout, '(a)') '  -mp, --(no-)match-phase     match phase between assigned bra/ket states          '
        write(stdout, '(32x,a,l1)') 'default: ', phase_omat
        write(stdout, '(a)') '  -rao, --(no-)recenter-aos   attempt to remove effect of basis set translation by '
        write(stdout, '(a)') '                              recentering pairs of AOs in AO overlap calculation   '
        write(stdout, '(a)') '                              (untested, not recommended)                          '
        write(stdout, '(32x,a,l1)') 'default: ', center2
        write(stdout, '(a)') '  -fmn, --freeze-mo-norm t    freeze occupied ket MOs when their norm in bra basis '
        write(stdout, '(a)') '                              is smaller than given threshold. Same number of bra  '
        write(stdout, '(a)') '                              MOs with smallest norms in ket basis is also frozen. '
        write(stdout, '(a)') '                              Used when geometry of a small part of a system is    '
        write(stdout, '(a)') '                              significantly different between bra and ket states.  '
        write(stdout, '(a)') '                              (untested, not recommended)                          '
        write(stdout, '(32x,a,l1)') 'default: ', f_by_mo_norm
        write(stdout, '(a)') '  -t, --threshold t           truncate wave functions using given threshold        '
        write(stdout, '(a)') '  -o, --outfile file          output final overlap matrix to file                  '
        write(stdout, '(32x,a,a)') 'default: ',  outfile
        write(stdout, '(a)') '  -p, --print-level p         control output level of program (0 = quiet)          '
        write(stdout, '(32x,a,i0)') 'default: ', print_level
        stop
    end subroutine print_help


end program cis_overlap_prog
