program cis_nto_prog
    use atom_basis_mod
    use ccg_ao_mod
    use cis_nto_mod
    use occupation_mod
    use write_molden_mod
    use write_txt_mod
    use read_all_mod
    use blas95, only : gemm
    implicit none

    integer :: noa = 0 !< Number of occupied alpha orbitals.
    integer :: nob = 0 !< Number of occupied beta orbitals.
    real(dp), allocatable :: geom(:) !< Geometry.
    character(len=2), allocatable :: atsym(:) !< Atom symbols.
    integer, allocatable :: atnum(:) !< Atom numbers.
    type(atombasis), allocatable :: basis(:) !< Basis sets.
    integer, allocatable :: bindex(:) !< BS index for each atom.
    type(ccg), allocatable :: ccgao(:) !< Gaussian basis functions.
    real(dp), allocatable :: moa(:, :) !< Molecular orbital coefficients alpha.
    real(dp), allocatable :: mob(:, :) !< Molecular orbital coefficients beta.
    real(dp), allocatable :: wfa(:, :, :) !< CI alpha single excitation coefficients.
    real(dp), allocatable :: wfb(:, :, :) !< CI beta single excitation coefficients.
    real(dp), allocatable :: nto_a(:, :, :) !< NTOs alpha.
    real(dp), allocatable :: nto_b(:, :, :) !< NTOs beta.
    real(dp), allocatable :: nto_c_a(:, :) !< NTO alpha coefficients.
    real(dp), allocatable :: nto_c_b(:, :) !< NTO beta coefficients.
    integer, allocatable :: na_a(:) !< Number of active alpha NTO orbitals.
    integer, allocatable :: na_b(:) !< Number of active beta NTO orbitals.
    real(dp), allocatable :: c2s(:, :) !< Transformation between cartesian and spherical basis functions.
    real(dp), allocatable :: moc_a(:, :) !< MO coefficients in cartesian basis alpha.
    real(dp), allocatable :: moc_b(:, :) !< MO coefficients in cartesian basis beta.
    logical, allocatable :: occ(:, :) !< Occupied MO mask.
    logical, allocatable :: act(:, :) !< Active MO mask.

    logical :: beta = .false.
    integer :: i, j
    character(len=1000) :: temp
    real(dp) :: thr
    integer :: outunit
    character(len=:), allocatable :: input_format
    character(len=:), allocatable :: outfile
    character(len=:), allocatable :: dir

    input_format = 'turbomole'
    write(stdout, *) 'Path to electronic structure calculation:'
    read(stdin, '(a)') temp
    dir = trim(adjustl(temp))
    write(stdout, *) 'Output file:'
    read(stdin, '(a)') temp
    outfile = trim(adjustl(temp))
    write(stdout, *) 'Threshold for truncating wave functions:'
    read(stdin, *) thr

    call read_geom(input_format, dir, geom, atsym, atnum)
    call read_ccg_ao(input_format, dir, ccgao, trans_ao=c2s, basis=basis, basis_index=bindex)
    call read_mo(input_format, dir, moa_c=moa, mob_c=mob)
    call read_cis(input_format, dir, cisa=wfa, cisb=wfb, occ_mo=occ, act_mo=act, norm=.true.)
    call write_txt('wfa', wfa)
    if (allocated(wfb)) beta = .true.
    call sort_mo(occ(:, 1), act(:, 1), moa, remove_inactive = .true.)
    if (beta) call sort_mo(occ(:, 2), act(:, 2), mob, remove_inactive = .true.)

    allocate(moc_a(size(c2s, 2), size(moa, 2)))
    call gemm(c2s, moa, moc_a, transa='T')
    if (beta) then 
        allocate(moc_b(size(c2s, 2), size(mob, 2)))
        call gemm(c2s, mob, moc_b, transa='T')
    end if

    call cis_nto_and_convert_basis(moc_a, wfa, nto_c_a, nto_a)
    if (beta) call cis_nto_and_convert_basis(moc_b, wfb, nto_c_b, nto_b)
    call cis_nto_truncate(beta, thr, nto_c_a, nto_c_b, na_a, na_b)
    noa = size(nto_c_a, 1)
    nob = size(nto_c_b, 1)
    open(newunit=outunit, file=outfile, action='write')
    write(outunit, '(a)') '[Molden Format]'
    write(outunit, '(a)') '[Title]'
    write(outunit, '(a)') 'Dyson'
    call write_molden_atoms(outunit, geom, atsym, atnum)
    call write_molden_gto(outunit, basis, bindex)
    write(outunit, '(a)') '[MO]'
    do i = 1, size(nto_a, 3)
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
