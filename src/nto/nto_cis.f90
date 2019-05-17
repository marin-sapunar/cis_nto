!----------------------------------------------------------------------------------------------
! MODULE: nto_cis_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date March, 2019
!----------------------------------------------------------------------------------------------
module nto_cis_mod
    use global_defs
    implicit none


    private
    public nto_cis


contains

    subroutine nto_cis
        use nto_variables
        use read_all_mod
        use cis_util_mod
        use cis_nto_mod
        logical :: beta

        ! Read input
        time0 =  omp_get_wtime()
        call read_geom(input_format, path, geom, atsym, atnum)
        call read_basis(input_format, path, bs)
        call read_mo(input_format, path, mos)
        call read_cis(input_format, path, cisa=cisa, cisb=cisb, occ_mo=occ, act_mo=act, occ_num=on)
        beta = .false.
        if (allocated(mos%cb)) then
            rhf = 2
            beta = .true.
        end if
        call print_calculation_properties(size(cisa, 3), on, '')
        ! -- Sort MOs
        call mo_reshape(.true., .true., occ(:, 1), act(:, 1), 2, a2=mos%ca)
        if (rhf == 2) call mo_reshape(.true., .true., occ(:, 2), act(:, 2), 2, a2=mos%ca)
        time_in =  time_in + omp_get_wtime() - time0

        ! Generate NTOs
        call cis_nto_and_convert_basis(mos%ca, cisa, nto_c_a, nto_a)
        if (rhf == 2) call cis_nto_and_convert_basis(mos%cb, cisb, nto_c_b, nto_b)
        call cis_nto_truncate(beta, wf_threshold, truncate_nex, nto_c_a, nto_c_b, na_a, na_b)
    end subroutine nto_cis


!   subroutine nto_tmp()
!       use nto_variables
!       
!       use cis_nto_mod
!       use read_all_mod
!       use one_el_op_mod
!       use linalg_wrapper_mod, only : dot, gemm
!       use matrix_mod, only : mat_ge_mmm
!       logical :: beta
!       integer :: i, j
!       real(dp), allocatable :: x_wrk(:, :)
!       real(dp), allocatable :: y_wrk(:, :)
!       real(dp), allocatable :: z_wrk(:, :)
!       real(dp), allocatable :: r_wrk(:, :)
!       real(dp), allocatable :: wrk(:, :)
!       real(dp), allocatable :: x_ao(:, :)
!       real(dp), allocatable :: y_ao(:, :)
!       real(dp), allocatable :: z_ao(:, :)
!       real(dp) :: norm

!       ! AO overlap.
!       call one_el_op(geom, geom, bs, bs, [2, 0, 0], x_ao)
!       call one_el_op(geom, geom, bs, bs, [0, 2, 0], y_ao)
!       call one_el_op(geom, geom, bs, bs, [0, 0, 2], z_ao)

!       ! Calculate NTOs.
!     ! call cis_nto(cisa, nto_c_a, nto_a)
!       allocate(na_a(size(nto_a, 3)))
!       na_a = 1
!      !call cis_nto_truncate(beta, wf_threshold, nto_c_a, nto_c_b, na_a, na_b)

!       allocate(wrk(size(moa, 1), size(nto_a, 2)))
!       allocate(x_wrk(size(nto_a, 2), size(nto_a, 2)))
!       allocate(y_wrk(size(nto_a, 2), size(nto_a, 2)))
!       allocate(z_wrk(size(nto_a, 2), size(nto_a, 2)))
!       allocate(r_wrk(size(nto_a, 2), size(nto_a, 2)))
!       norm = 1.0_dp / sqrt(dot(geom, geom))
!       write(*, *) ' i,    coeff,   z_hole, z_particle'
!       do i = 1, size(nto_a, 3)
!       !   call gemm(moa, nto_a(:, :, i), wrk)
!           wrk = nto_a(:, :, i)
!           call mat_ge_mmm(wrk, x_ao, wrk, x_wrk, transa='T')
!           call mat_ge_mmm(wrk, y_ao, wrk, y_wrk, transa='T')
!           call mat_ge_mmm(wrk, z_ao, wrk, z_wrk, transa='T')
!           r_wrk = (x_wrk + y_wrk + z_wrk)!/ sqrt(dot(geom, geom))
!           x_wrk = x_wrk!/ sqrt(dot(geom(1::3), geom(1::3)))
!           y_wrk = y_wrk!/ sqrt(dot(geom(2::3), geom(2::3)))
!           z_wrk = z_wrk!/ sqrt(dot(geom(3::3), geom(3::3)))
!           r_wrk = (x_wrk + y_wrk + z_wrk)
!           write(*,'(i3,3f10.5)')  i, nto_c_a(1, i), z_wrk(1, 1), z_wrk(on%ao(1)+1, on%ao(1)+1)
!        !  do j = 1, na_a(i)
!        !      write(*,*) 'x', nto_c_a(j, i), x_wrk(j, j), x_wrk(on%ao(1)+j, on%ao(1)+j)
!        !      write(*,*) 'y', nto_c_a(j, i), y_wrk(j, j), y_wrk(on%ao(1)+j, on%ao(1)+j)
!        !      write(*,*) 'z', nto_c_a(j, i), z_wrk(j, j), z_wrk(on%ao(1)+j, on%ao(1)+j)
!        !      write(*,*) 'r', nto_c_a(j, i), r_wrk(j, j), r_wrk(on%ao(1)+j, on%ao(1)+j)
!        !    ! write(*,*) 'z/r', nto_c_a(j, i), z_wrk(j, j)/r_wrk(j, j), z_wrk(on%ao(1)+j, on%ao(1)+j)/r_wrk(on%ao(1)+j, on%ao(1)+j)
!        !    ! write(*,*) 'z*r', nto_c_a(j, i), z_wrk(j, j)*r_wrk(j, j), z_wrk(on%ao(1)+j, on%ao(1)+j)*r_wrk(on%ao(1)+j, on%ao(1)+j)
!        !  end do
!       end do

!   end subroutine nto_tmp


end module nto_cis_mod
