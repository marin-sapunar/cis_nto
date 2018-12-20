program cis_olap_test
    use global_defs
    use read_all_mod
    use occupation_mod
    use ccg_ao_mod
    use one_el_op_mod
    use matrix_mod
    use orthog_mod
    use phase_mod
    use cis_overlap_nto_mod
    use cis_overlap_cis_mod
    implicit none

    integer :: rhf = 0 !< Restricted (1) or unrestricted(2) calculation.
    integer :: rhf1 = 1 !< Restricted (1) or unrestricted (2) for 1 wave functions.
    integer :: rhf2 = 1 !< Restricted (1) or unrestricted (2) for 2 wave functions.
    integer :: nwf1 !< Number of states 1.
    integer :: nwf2 !< Number of states 2.
    type(occupation_numbers) :: on1 !< Occupation numbers 1.
    type(occupation_numbers) :: on2 !< Occupation numbers 2.
    type(ccg), allocatable :: ccg1(:) !< Atomic orbitals 1.
    type(ccg), allocatable :: ccg2(:) !< Atomic orbitals 2.
    real(dp), allocatable :: geom1(:) !< Geometry 1.
    real(dp), allocatable :: geom2(:) !< Geometry 2.
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
    real(dp), allocatable :: s_ao(:, :) !< Atomic orbital overlaps.
    real(dp), allocatable :: s_mo(:, :, :) !< Molecular orbital overlaps.
    real(dp), allocatable :: s_wf(:, :) !< Wave function overlaps.
    character(len=1000) :: temp
    character(len=:), allocatable :: input_format
    character(len=:), allocatable :: dir1
    character(len=:), allocatable :: dir2
    real(dp) :: thr
    integer :: i, narg
    integer :: ounit
    logical :: norm
    logical :: orth
    logical :: orth_omat
    logical :: phase_omat
    logical :: center1
    logical :: center2
    integer :: alg
    real(dp), external :: omp_get_wtime
    real(dp) :: time00, time0
    real(dp) :: time_io, time_ao, time_mo, time_wf, time_tot

    print_level = 4
    time00 = omp_get_wtime()

    input_format = 'turbomole'
    narg = command_argument_count()
    if (narg /= 2) stop 'Call program with paths to bra and ket state calculations.'
    call get_command_argument(1, temp)
    dir1 = trim(adjustl(temp))
    call get_command_argument(2, temp)
    dir2 = trim(adjustl(temp))
    write(stdout, *) 'Algorithm for overlap calculation (1 - CIS, 2 - NTO):'
    read(stdin, *) alg
    write(stdout, *) 'Threshold for truncating wave functions:'
    read(stdin, *) thr
    write(stdout, *) 'Normalize wave functions before overlap calculation:'
    read(stdin, *) norm
    write(stdout, *) 'Orthogonalize wave functions before overlap calculation:'
    read(stdin, *) orth
    write(stdout, *) 'Orthogonalize overlap matrix after calculation:'
    read(stdin, *) orth_omat
    write(stdout, *) 'Match phase of states after calculation:'
    read(stdin, *) phase_omat
    write(stdout, *) 'Remove translation from diagonal blocks of atomic orbital overlap matrix:'
    read(stdin, *) center1
    write(stdout, *) 'Remove translation from atomic orbital overlap matrix:'
    read(stdin, *) center2

    ! Read input.
    time0 = omp_get_wtime()
    if (print_level >= 1) then
        write(stdout, *) 
        write(stdout, '(1x,a)') 'Reading input from turbomole files...'
        write(stdout, '(5x,a,a)') ' Directory containing calculation for bra states: ', dir1
        write(stdout, '(5x,a,a)') ' Directory containing calculation for ket states: ', dir2
    end if
    call read_ccg_ao(input_format, dir1, ccg1, geom=geom1, trans_ao=trans1)
    call read_mo(input_format, dir1, moa_c=moa1, mob_c=mob1)
    call read_cis(input_format, dir1, wfa=wfa1, wfb=wfb1, occ_mo=occ1, act_mo=act1, norm=norm, &
    &             orthog=orth, occ_num = on1)
    nwf1 = size(wfa1, 2)
    if (print_level >= 1) then
        write(stdout, *) 
        write(stdout, '(1x,a)') 'Properties for bra states:'
        if (rhf1 == 1) write(stdout, '(5x, a)') 'Restricted calculation.'
        if (rhf1 == 2) write(stdout, '(5x, a)') 'Unrestricted calculation.'
        write(stdout, '(5x,a,2(1x,i0))') 'Number of excited states:           ', nwf1
        write(stdout, '(5x,a,2(1x,i0))') 'Number of occupied orbitals:        ', on1%o(1:rhf1)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of occupied orbitals:        ', on1%o(1:rhf1)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of virtual orbitals:         ', on1%v(1:rhf1)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of active occupied orbitals: ', on1%ao(1:rhf1)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of active virtual orbitals:  ', on1%av(1:rhf1)
    end if

    call read_ccg_ao(input_format, dir2, ccg2, geom=geom2, trans_ao=trans2)
    call read_mo(input_format, dir2, moa_c=moa2, mob_c=mob2)
    call read_cis(input_format, dir2, wfa=wfa2, wfb=wfb2, occ_mo=occ2, act_mo=act2, norm=norm, &
    &             orthog=orth, occ_num = on2)
    nwf2 = size(wfa2, 2)
    if (print_level >= 1) then
        write(stdout, *) 
        write(stdout, '(1x,a)') 'Properties for ket states:'
        if (rhf2 == 1) write(stdout, '(5x, a)') 'Restricted calculation.'
        if (rhf2 == 2) write(stdout, '(5x, a)') 'Unrestricted calculation.'
        write(stdout, '(5x,a,2(1x,i0))') 'Number of excited states:           ', nwf2
        write(stdout, '(5x,a,2(1x,i0))') 'Number of occupied orbitals:        ', on2%o(1:rhf2)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of occupied orbitals:        ', on2%o(1:rhf2)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of virtual orbitals:         ', on2%v(1:rhf2)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of active occupied orbitals: ', on2%ao(1:rhf2)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of active virtual orbitals:  ', on2%av(1:rhf2)
    end if
    time_io = omp_get_wtime() - time0

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
    if (allocated(wfb1)) rhf1 = 2
    if (allocated(wfb2)) rhf2 = 2
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
            write(stdout, '(1x,a)') 'Unrestricted bra and restricted ket calculation...'
            write(stdout, '(5x,a)') 'Renormalizing ket CI coefficients...'
        end if
        wfa1 = wfa1 / sqrt(2.0_dp)
        wfa2 = wfa2 / sqrt(2.0_dp)
        allocate(mob2, source=moa2)
        allocate(wfb2, source=wfa2)
        wfb2 = -wfb2
    end if

    ! Remove frozen mos and calculate MO overlap matrix.
    if (print_level >= 1) then
        write(stdout, *) 
        write(stdout, '(1x,a)') 'Removing frozen MOs...'
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
    end select
    time_wf = omp_get_wtime() - time0

    ! Output.
    if (print_level >= 1) then
        write(stdout, *) 
        write(stdout, '(1x,a)') 'Raw overlap matrix:'
        do i = 0, nwf2
            write(ounit,'(5x,1000es24.16)') s_wf(i, :)
        end do
    end if

    if (orth_omat) then
        call orthog_lowdin(s_wf)
        if (print_level >= 1) then
            write(stdout, *) 
            write(stdout, '(1x,a)') 'Orthogonalized overlap matrix:'
            do i = 0, nwf2
                write(ounit,'(5x,1000es24.16)') s_wf(i, :)
            end do
        end if
    end if

    if (phase_omat) then
        call phasematch_assigned(s_wf)
        if (print_level >= 1) then
            write(stdout, *) 
            write(stdout, '(1x,a)') 'Overlap matrix with phase matching between assigned bra/ket states:'
            do i = 0, nwf2
                write(ounit,'(5x,1000es24.16)') s_wf(i, :)
            end do
        end if
    end if

    open(newunit=ounit, file='omat', action='write')
    do i = 0, nwf2
        write(ounit,*) s_wf(:, i)
    end do
    close(ounit)

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

end program cis_olap_test
