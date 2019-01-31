program cis_nto_prog
    ! General
    use global_defs
    use file_mod
    ! Chem
    use atom_basis_mod
    use occupation_mod
    use ccg_ao_mod
    use cis_nto_mod
    ! I/O
    use read_all_mod
    use write_molden_mod
    ! External
    use blas95, only : gemm
    implicit none

    ! System
    logical :: beta = .false.
    integer :: rhf = 1
    real(dp), allocatable :: geom(:) !< Geometry.
    character(len=2), allocatable :: atsym(:) !< Atom symbols.
    integer, allocatable :: atnum(:) !< Atom numbers.
    type(atombasis), allocatable :: basis(:) !< Basis sets.
    integer, allocatable :: bindex(:) !< Basis set index for each atom.
    type(occupation_numbers) :: on !< Occupation numbers.
    type(ccg), allocatable :: ccg1(:) !< Atomic orbitals.
    real(dp), allocatable :: trans(:, :) !< Transformation matrix for basis functions.
    real(dp), allocatable :: moa(:, :) !< Molecular orbital coefficients alpha.
    real(dp), allocatable :: mob(:, :) !< Molecular orbital coefficients beta.
    real(dp), allocatable :: moc_a(:, :) !< MO coefficients in cartesian basis alpha.
    real(dp), allocatable :: moc_b(:, :) !< MO coefficients in cartesian basis beta.
    real(dp), allocatable :: cisa(:, :, :) !< CIS matrix alpha.
    real(dp), allocatable :: cisb(:, :, :) !< CIS matrix beta
    logical, allocatable :: occ(:, :) !< Occupied MO mask.
    logical, allocatable :: act(:, :) !< Active MO mask.
    ! Results
    real(dp), allocatable :: nto_a(:, :, :) !< NTOs alpha.
    real(dp), allocatable :: nto_b(:, :, :) !< NTOs beta.
    real(dp), allocatable :: nto_c_a(:, :) !< NTO alpha coefficients.
    real(dp), allocatable :: nto_c_b(:, :) !< NTO beta coefficients.
    integer, allocatable :: na_a(:) !< Number of active alpha NTO orbitals.
    integer, allocatable :: na_b(:) !< Number of active beta NTO orbitals.
    ! Options
    character(len=:), allocatable :: input_format
    character(len=:), allocatable :: dir
    character(len=:), allocatable :: outfile
    real(dp) :: thr
    ! Help
    character(len=1000) :: temp
    integer :: i, j, narg, outunit

    ! Defaults
    print_level = 0
    outfile = 'nto.molden'
    input_format = 'turbomole'
    thr = 1.0_dp
    ! CLI
    narg = command_argument_count()
    i = 0
    do while (i < narg)
        i = i + 1
        call get_command_argument(i, temp)
        if (i > narg) temp = '--help'
        select case(temp)
        case('--help', '-h')
            write(stdout, *) 'usage: cis_nto.exe [optional arguments] dir'
            write(stdout, *)
            write(stdout, *) 'Return natural transition orbitals from a CIS/LR-TDDFT calculation.'
            write(stdout, *)
            write(stdout, *) 'positional arguments:'
            write(stdout, *) '  dir                   directory containing the calculation                 '
            write(stdout, *)
            write(stdout, *) 'optional arguments:'
            write(stdout, *) '  -h, --help                 show this help message and exit                      '
            write(stdout, *) '  -t, --threshold t          truncate wave functions using given threshold        '
            write(stdout, '(29x,a,f7.4)') 'default: ', thr
            write(stdout, *) '  -o, --outfile file         output final NTOs to file                            '
            write(stdout, '(29x,a,a)') 'default: ', outfile
            write(stdout, *) '  -p, --print-level p        control output level of program (0 = quiet)           '
            write(stdout, '(29x,a,i0)') 'default: ', print_level
            stop
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
    call get_command_argument(narg, temp)
    dir = trim(adjustl(temp))

    ! Begin run.
    if (print_level >= 1) then
        write(stdout, '(5x,a)') '-------------------------------------------------------------'
        write(stdout, '(5x,a)') '                       cis_nto program                       '
        write(stdout, '(5x,a)') '                         version 1.0                         '
        write(stdout, '(5x,a)') '                                                             '
        write(stdout, '(5x,a)') ' Program compiled on '//__DATE__//' '//__TIME__//'.          '
        write(stdout, '(5x,a)') '-------------------------------------------------------------'
        write(stdout, *)
    end if


    ! Read input.
    if (print_level >= 1) then
        write(stdout, *)
        write(stdout, '(1x,a)') 'Reading input from turbomole files...'
        write(stdout, '(5x,a,a)') ' Directory containing calculation: ', dir
    end if
    if (.not. is_dir(dir)) then
        write(stderr, *)
        write(stderr, '(1x,a,a,a)') 'Error. Directory ', dir, ' not found.'
        stop
    end if
    call read_geom(input_format, dir, geom, atsym, atnum)
    call read_ccg_ao(input_format, dir, ccg1, basis=basis, basis_index=bindex, trans_ao=trans)
    call read_mo(input_format, dir, moa_c=moa, mob_c=mob)
    call read_cis(input_format, dir, cisa=cisa, cisb=cisb, occ_mo=occ, act_mo=act, occ_num=on, &
    &             orthog=.false., norm=.true.)
    if (allocated(cisb)) beta = .true.
    if (allocated(cisb)) rhf = 2
    if (print_level >= 1) then
        write(stdout, *)
        write(stdout, '(1x,a)') 'Electronic structure calculation properties:'
        if (rhf == 1) write(stdout, '(5x, a)') 'Restricted calculation.'
        if (rhf == 2) write(stdout, '(5x, a)') 'Unrestricted calculation.'
        write(stdout, '(5x,a,2(1x,i0))') 'Number of excited states:           ', size(cisa, 3)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of basis functions (sphe):   ', size(trans, 1)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of basis functions (cart):   ', size(trans, 2)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of orbitals:                 ', on%n
        write(stdout, '(5x,a,2(1x,i0))') 'Number of occupied orbitals:        ', on%o(1:rhf)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of virtual orbitals:         ', on%v(1:rhf)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of active occupied orbitals: ', on%ao(1:rhf)
        write(stdout, '(5x,a,2(1x,i0))') 'Number of active virtual orbitals:  ', on%av(1:rhf)
    end if

    ! Remove frozen mos and transform mos to cartesian basis.
    if (print_level >= 1) then
        if (.not. all(act)) then
            write(stdout, *)
            write(stdout, '(1x,a)') 'Removing frozen MOs...'
        end if
    end if
    call sort_mo(occ(:, 1), act(:, 1), moa, remove_inactive=.true.)
    allocate(moc_a(size(trans, 2), size(moa, 2)))
    call gemm(trans, moa, moc_a, transa='T')
    if (rhf == 2) then 
        call sort_mo(occ(:, 2), act(:, 2), mob, remove_inactive=.true.)
        allocate(moc_b(size(trans, 2), size(mob, 2)))
        call gemm(trans, mob, moc_b, transa='T')
    end if

    ! Calculate NTOs.
    if (print_level >= 2) then
        write(stdout, *)
        write(stdout, '(1x,a)') 'Computing natural transition orbitals...'
    end if
    call cis_nto_and_convert_basis(moc_a, cisa, nto_c_a, nto_a)
    if (rhf == 2) call cis_nto_and_convert_basis(moc_b, cisb, nto_c_b, nto_b)
    call cis_nto_truncate(beta, thr, nto_c_a, nto_c_b, na_a, na_b)

    ! Output.
    if (print_level >= 1) then
        write(stdout, *)
        write(stdout, '(1x,a)') 'Printing NTOs to molden file '//outfile//'...'
    end if

    open(newunit=outunit, file=outfile, action='write')
    write(outunit, '(a)') '[Molden Format]'
    write(outunit, '(a)') '[Title]'
    write(outunit, '(a)') 'cis_nto output'
    call write_molden_atoms(outunit, geom, atsym, atnum)
    call write_molden_gto(outunit, basis, bindex)
    write(outunit, '(a)') '[MO]'
    do i = 1, size(nto_a, 3)
        do j = 1, na_a(i)
            call write_molden_mo_single(outunit, 1000*i+2*j-1, 'a   ', 1, 0.0_dp, nto_c_a(j, i),   &
            &                           nto_a(:, j, i))
            call write_molden_mo_single(outunit, 1000*i+2*j, 'a   ', 1, 0.0_dp, nto_c_a(j, i),     &
            &                           nto_a(:, on%ao(1)+j, i))
        end do
        if (rhf == 2) then
            do j = 1, na_b(i)
                call write_molden_mo_single(outunit, 1000*i+2*j-1, 'a   ', 2, 0.0_dp, nto_c_b(j, i), &
                &                           nto_b(:, j, i))
                call write_molden_mo_single(outunit, 1000*i+2*j, 'a   ', 2, 0.0_dp, nto_c_b(j, i),   &
                &                           nto_b(:, on%ao(2)+j, i))
            end do
        end if
    end do

    if (print_level >= 1) then
        write(stdout, *)
        write(stdout, '(1x,a)') 'cis_nto done                                             '
        write(stdout, '(1x,a)') '-------------------------------------------------------------'
    end if

end program cis_nto_prog
