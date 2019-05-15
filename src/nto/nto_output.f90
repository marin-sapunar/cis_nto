!--------------------------------------------------------------------------------------------------
! PROGRAM: nto_output_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date March, 2019
!
!> @brief Subroutines for printing output of NTO program
!--------------------------------------------------------------------------------------------------
module nto_output_mod
    use global_defs
    use nto_variables
    use write_all_mod
    implicit none


    private
    public :: output_nto


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: output_nto
    !> @brief Print NTOs
    !----------------------------------------------------------------------------------------------
    subroutine output_nto()
        use basis_transform_mod
        real(dp), allocatable :: out_moa(:, :)
        real(dp), allocatable :: out_mob(:, :)
        real(dp), allocatable :: out_oa(:)
        real(dp), allocatable :: out_ob(:)
        integer, allocatable :: out_na(:)
        integer, allocatable :: out_nb(:)
        integer :: tot_n_a, tot_n_b, i, j, ca, cb, n_occ_a, n_occ_b
        type(basis_transform) :: trans


        call trans%init(bs, bs%source_format, output_format)

        tot_n_a = sum(na_a) * 2
        n_occ_a = size(nto_c_a, 1)
        call trans%transform(nto_a, 1)
        allocate(out_moa(size(nto_a, 1), tot_n_a))
        allocate(out_na(tot_n_a))
        allocate(out_oa(tot_n_a))
        if (rhf == 2) then
            tot_n_b = sum(na_b) * 2
            n_occ_b = size(nto_c_b, 1)
            call trans%transform(nto_b, 1)
            allocate(out_mob(size(nto_b, 1), tot_n_b))
            allocate(out_nb(tot_n_b))
            allocate(out_ob(tot_n_b))
        end if
        ca = 0
        cb = 0
        do i = 1, size(nto_a, 3)
            do j = 1, na_a(i)
                ca = ca + 1
                out_moa(:, ca) = nto_a(:, j, i)
                out_na(ca) = 1000*i+2*j-1
                out_oa(ca) = nto_c_a(j, i)
                ca = ca + 1
                out_moa(:, ca) = nto_a(:, j+n_occ_a, i)
                out_na(ca) = 1000*i+2*j
                out_oa(ca) = nto_c_a(j, i)
            end do
            if (rhf /= 2) cycle
            do j = 1, na_b(i)
                cb = cb + 1
                out_mob(:, cb) = nto_b(:, j, i)
                out_nb(cb) = 1000*i+2*j-1
                out_ob(cb) = nto_c_b(j, i)
                cb = cb + 1
                out_mob(:, cb) = nto_b(:, j+n_occ_b, i)
                out_nb(cb) = 1000*i+2*j
                out_ob(cb) = nto_c_b(j, i)
            end do
        end do

        call nto_mos%init(ca=out_moa, cb=out_mob, oa=out_oa, ob=out_ob, na=out_na, nb=out_nb)
        call write_all(output_format, outfile_nto, geom=geom, atom_symbol=atsym, atom_number=atnum,&
        &              basis=bs, mos=nto_mos)
    end subroutine output_nto


end module nto_output_mod
