!--------------------------------------------------------------------------------------------------
! MODULE: one_el_op_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date November, 2016
!
!> @brief Contains the subroutine for computing a 1 el. operator over two sets of basis functions.
!--------------------------------------------------------------------------------------------------
module one_el_op_mod
    use global_defs
    implicit none

    private
    public :: one_el_op
    public :: cart_one_el_op


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: one_el_op
    !
    !> @brief Matrix of a simple one electron operator over two sets of basis functions.
    !> @details
    !! Calculates a matrix with elements:
    !!
    !!     S_12(i,j) = < B_1,i | x^lx*y^ly*z^lz | B_2,j >
    !!
    !! The matrix is initially constructed in the basis of contracted Cartesian Gaussian functions
    !! ordered by 'internal' order (defined by functions given in ang_mom_defs). Afterwards, the
    !! basis is changed to the source format given in the basis_set variables.
    !! The center_atoms and center_pairs options can be used to partially remove the effect of
    !! translation from the overlaps. If center_atoms is .true., only overlaps between basis 
    !! functions on the same atom are calculated as if the atom is at the origin. If center_pairs
    !! is .true., overlaps between functions on two atoms are calculated as if the midpoint between
    !! the two atoms is at the origin.
    !----------------------------------------------------------------------------------------------
    subroutine one_el_op(geom1, geom2, bs1, bs2, lxyz, mat, center_atoms, center_pairs)
        use linalg_wrapper_mod, only : gemm
        use basis_set_mod
        use basis_transform_mod
        real(dp), intent(in) :: geom1(:) !< Geometry 1.
        real(dp), intent(in) :: geom2(:) !< Geometry 2.
        type(basis_set), intent(inout) :: bs1 !< Basis set 1.
        type(basis_set), intent(inout) :: bs2 !< Basis set 2.
        integer, intent(in) :: lxyz(3) !< List of exponents of the coordinates in the operator.
        real(dp), allocatable, intent(out) :: mat(:, :) !< Operator matrix.
        logical, intent(in), optional :: center_atoms !< Translate same atom to center.
        logical, intent(in), optional :: center_pairs !< Translate pairs of atoms to center.
        logical :: opt_atoms
        logical :: opt_pairs
        type(basis_transform) :: trans1
        type(basis_transform) :: trans2

        opt_atoms = .false.
        opt_pairs = .false.
        if (present(center_atoms)) opt_atoms = center_atoms
        if (present(center_pairs)) opt_pairs = center_pairs

        ! Calculate matrix in basis of Cartesian GTOs in default order.
        call cart_one_el_op(geom1, geom2, bs1, bs2, lxyz, mat, opt_atoms, opt_pairs)

        ! Convert to original basis format.
        call trans1%init(bs1, 'internal', bs1%source_format)
        call trans2%init(bs2, 'internal', bs2%source_format)
        call trans1%transform(mat, 1)
        call trans2%transform(mat, 2)

    end subroutine one_el_op


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cart_one_el_op
    !
    !> @brief Matrix of a simple one electron operator over two sets of cartesian basis functions.
    !> @details
    !! Calculates a matrix with elements:
    !!
    !!     S_12(i,j) = < B_1,i | x^lx*y^ly*z^lz | B_2,j >
    !!
    !! where < B_1,i | and | B_2,j > are contracnted Cartesian Gaussian functions. 
    !----------------------------------------------------------------------------------------------
    subroutine cart_one_el_op(geom1, geom2, bs1, bs2, lxyz, mat, center_atoms, center_pairs)
        use basis_set_mod
        use ang_mom_defs, only : amp
        use cgto_mod, only : cgto_one_el
        real(dp), intent(in) :: geom1(:) !< Geometry 1.
        real(dp), intent(in) :: geom2(:) !< Geometry 2.
        type(basis_set), intent(inout) :: bs1 !< Basis set 1.
        type(basis_set), intent(inout) :: bs2 !< Basis set 2.
        integer, intent(in) :: lxyz(3) !< List of exponents of the coordinates in the operator.
                                       !! (lx, ly, lz)
        real(dp), allocatable, intent(out) :: mat(:, :) !< Operator matrix.
        logical, intent(in), optional :: center_atoms !< Translate same atom to center.
        logical, intent(in), optional :: center_pairs !< Translate pairs of atoms to center.
        logical :: opt_atoms
        logical :: opt_pairs

        integer :: i, j, k, l, ck0, cl0, ck, cl, nbf1, nbf2, bi, bj
        real(dp) :: qi(3), qj(3), qc(3)

        if (allocated(mat)) deallocate(mat)
        allocate(mat(bs1%n_bf_cart, bs2%n_bf_cart), source=num0)

        opt_atoms = .false.
        opt_pairs = .false.
        if (present(center_atoms)) opt_atoms = center_atoms
        if (present(center_pairs)) opt_pairs = center_pairs

        ck = 0
        do i = 1, bs1%n_center
            bi = bs1%center_i_bs(i)
            ck0 = ck
            cl = 0
            do j = 1, bs2%n_center
                bj = bs1%center_i_bs(j)
                cl0 = cl
                if (opt_pairs) then
                    qi = (geom1(3*i-2:3*i) - geom1(3*j-2:3*j)) / 2
                    qj = (geom2(3*j-2:3*j) - geom2(3*i-2:3*i)) / 2
                else if (opt_atoms) then
                    if (i == j) then
                        qi = num0
                        qj = num0
                    end if
                else
                    qi = geom1(3*i-2:3*i)
                    qj = geom2(3*j-2:3*j)
                end if
                qc = num0
                ck = ck0
                do k = 1, bs1%bs(bi)%n_subshell
                    nbf1 = amp%n_cart(bs1%bs(bi)%cg(k)%l)
                    cl = cl0
                    do l = 1, bs2%bs(bj)%n_subshell
                        nbf2 = amp%n_cart(bs2%bs(bj)%cg(l)%l)
                        mat(ck+1:ck+nbf1, cl+1:cl+nbf2) = cgto_one_el(bs1%bs(bi)%cg(k), bs2%bs(bj)%cg(l), qi, qc, qj, lxyz)
                        cl = cl + nbf2
                    end do
                    ck = ck + nbf1
                end do
            end do
        end do
    end subroutine cart_one_el_op


end module one_el_op_mod
