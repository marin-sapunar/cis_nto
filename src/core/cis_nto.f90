!----------------------------------------------------------------------------------------------
! MODULE: cis_nto_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date October, 2018
!
! DESCRIPTION:
!> @brief Subroutines for working with natural transition orbitals from CIS wave functions.
!----------------------------------------------------------------------------------------------
module cis_nto_mod
    use global_defs
    implicit none


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_nto_1state
    !
    !> @brief Generate natural transition orbitals for a CIS wave function.
    !> @details
    !! Uses SVD decomposition of the CI coefficients matrix of a given spin to generate the natural
    !! transition orbitals for that spin. The NTOs are given as linear combinations of canonical 
    !! orbitals in the nat_orbs array. The array is filled with the no hole orbitals followed by 
    !! the first no particle orbitals.
    !
    !> @note The SVD decomposition is unique up to a sign for each pair of singular vectors (for 
    !! non-degenerate singular values).
    !
    !> @note When transforming the wave function to the NTO basis, the sign of can flip. This can
    !! be seen from the overlap between the reference determinant in the NTO and canonical basis. 
    !! If it is negative (-1), the sign of the wave function is fixed by changing the sign of the
    !! NTO pair corresponding to the smallest singular value.
    !----------------------------------------------------------------------------------------------
    subroutine cis_nto_1state(nv, no, wf, nto_c, nto_mo)
        use lapack95, only : gesvd
        use matrix_mod, only : mat_ge_det
        integer, intent(in) :: nv !< Number of virtual orbitals.
        integer, intent(in) :: no !< Number of occupied orbitals.
        real(dp), intent(in) :: wf(nv, no) !< CIS wave function coefficients.
        real(dp), intent(out) :: nto_c(no) !< NTO coefficients.
        real(dp), intent(out) :: nto_mo(no+nv, 2*no) !< NTOs in the MO basis.
        real(dp) :: wrk_c(nv, no) !< Array for CI coefficients.
        real(dp) :: wrk_o(no, no) !< Array for occupied NTOs.
        real(dp) :: wrk_v(nv, no) !< Array for virtual NTOs.
        wrk_c = wf(:, :)
        call gesvd(wrk_c, u=wrk_v, s=nto_c, vt=wrk_o)
        nto_mo(1:no, 1:no) = transpose(wrk_o)
        nto_mo(no+1:no+nv, no+1:2*no) = wrk_v
        ! Check/fix sign.
        if (mat_ge_det(wrk_o) < 0) then
            nto_mo(:, no) = -nto_mo(:, no)
            nto_mo(:, 2*no) = -nto_mo(:, 2*no)
        end if
    end subroutine cis_nto_1state


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_nto
    !> @brief Generate NTOs for a set of CIS wave functions.
    !----------------------------------------------------------------------------------------------
    subroutine cis_nto(wf, nto_c, nto_mo)
        real(dp), intent(in) :: wf(:, :, :) !< CIS wave function coefficients.
        real(dp), allocatable, intent(out) :: nto_c(:, :) !< NTO coefficients.
        real(dp), allocatable, intent(out) :: nto_mo(:, :, :) !< NTOs in the MO basis.
        integer :: nv !< Number of virtual orbitals.
        integer :: no !< Number of occupied orbitals.
        integer :: ns !< Number of states.
        integer :: i
        nv = size(wf, 1)
        no = size(wf, 2)
        ns = size(wf, 3)
        if (.not. allocated(nto_mo)) allocate(nto_mo(no+nv, 2*no, ns), source=0.0_dp)
        if (.not. allocated(nto_c)) allocate(nto_c(no, ns), source=0.0_dp)
        do i = 1, ns
            call cis_nto_1state(nv, no, wf(:, :, i), nto_c(:, i), nto_mo(:, :, i))
        end do
    end subroutine cis_nto


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_nto_and_convert_basis
    !> @brief Generate NTOs and convert to different basis.
    !----------------------------------------------------------------------------------------------
    subroutine cis_nto_and_convert_basis(mo, wf, nto_c, nto_bas)
        use blas95, only : gemm
        real(dp), intent(in) :: mo(:, :) !< MO coefficients in target basis.
        real(dp), intent(in) :: wf(:, :, :) !< CIS wave function coefficients.
        real(dp), allocatable, intent(out) :: nto_c(:, :) !< NTO coefficients.
        real(dp), allocatable, intent(out) :: nto_bas(:, :, :) !< NTOs in target basis.
        real(dp), allocatable :: nto_mo(:, :) !< NTOs in MO basis.
        integer :: nv !< Number of virtual orbitals.
        integer :: no !< Number of occupied orbitals.
        integer :: ns !< Number of states.
        integer :: i

        nv = size(wf, 1)
        no = size(wf, 2)
        ns = size(wf, 3)
        if (.not. allocated(nto_bas)) allocate(nto_bas(size(mo, 1), 2*no, ns), source=0.0_dp)
        if (.not. allocated(nto_c)) allocate(nto_c(no, ns), source=0.0_dp)
        allocate(nto_mo(no+nv, 2*no), source=0.0_dp)
        do i = 1, ns
            call cis_nto_1state(nv, no, wf(:, :, i), nto_c(:, i), nto_mo)
            call gemm(mo, nto_mo, nto_bas(:, :, i)) 
        end do
    end subroutine cis_nto_and_convert_basis


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_nto_truncate_1state_u
    !
    !> @brief Determine number of alpha/beta orbitals in truncated CIS wave function.
    !> @details
    !! Determine number of NTOs of each spin needed so the truncated CIS wave function has a square
    !! norm higher than trunc.
    !----------------------------------------------------------------------------------------------
    subroutine cis_nto_truncate_1state_u(trunc, c_a, c_b, n_a, n_b)
        real(dp), intent(in) :: trunc !< Threshold for truncating the wave functions.
        real(dp), intent(inout) :: c_a(:) !< NTO coefficients alpha.
        real(dp), intent(inout) :: c_b(:) !< NTO coefficients beta.
        integer, intent(out) :: n_a !< Number of alpha coefficients in truncated wf.
        integer, intent(out) :: n_b !< Number of beta coefficients in truncated wf.
        real(dp) :: norm !< Squared norm.
        real(dp) :: ia !< Current alpha excitation coefficient.
        real(dp) :: ib !< Current beta excitation coefficient.
        integer :: i

        if (trunc >= 1.0_dp) then
            n_a = size(c_a, 1)
            n_b = size(c_b, 1)
            return
        end if

        n_a = 0
        n_b = 0
        norm = 0.0_dp
        ia = c_a(n_a+1)
        ib = c_b(n_b+1)
        do i = 1, size(c_a, 1) + size(c_b, 1)
            if (ia > ib) then
                n_a = n_a + 1
                norm = norm + ia**2
                ia = c_a(n_a+1)
            else
                n_b = n_b + 1
                norm = norm + ib**2
                ib = c_b(n_b+1)
            end if
            if (norm > trunc) exit
        end do
        ! Renormalize
      ! c_a(1:n_a) = c_a(1:n_a) / sqrt(norm)
      ! c_b(1:n_b) = c_b(1:n_b) / sqrt(norm)
    end subroutine cis_nto_truncate_1state_u


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_nto_truncate_1state_r
    !
    !> @brief Determine number of orbitals in truncated CIS wave function.
    !> @details
    !! Determine number of NTOs needed so the truncated CIS wave function has a square norm higher 
    !! than trunc.
    !----------------------------------------------------------------------------------------------
    subroutine cis_nto_truncate_1state_r(trunc, c_a, n_a)
        real(dp), intent(in) :: trunc !< Threshold for truncating the wave functions.
        real(dp), intent(inout) :: c_a(:) !< NTO coefficients alpha.
        integer, intent(out) :: n_a !< Number of alpha coefficients in truncated wf.
        real(dp) :: norm !< Squared norm.
        integer :: i

        if (trunc >= 1.0_dp) then
            n_a = size(c_a, 1)
            return
        end if

        n_a = 0
        norm = 0.0_dp
        do i = 1, size(c_a, 1)
            n_a = n_a + 1
            norm = norm + c_a(n_a)**2
            if (norm > trunc) exit
        end do
    end subroutine cis_nto_truncate_1state_r


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_nto_truncate
    !
    !> @brief Determine number of alpha/beta orbitals in truncated CIS wave functions.
    !> @details
    !! Determine number of NTOs of each spin needed so the truncated CIS wave functions each have
    !! a square norm higher than trunc.
    !----------------------------------------------------------------------------------------------
    subroutine cis_nto_truncate(beta, trunc, c_a, c_b, n_a, n_b)
        logical, intent(in) :: beta !< Unrestricted calculation.
        real(dp), intent(in) :: trunc !< Threshold for truncating the wave functions.
        real(dp), intent(inout) :: c_a(:, :) !< NTO coefficients alpha.
        real(dp), intent(inout) :: c_b(:, :) !< NTO coefficients beta.
        integer, allocatable, intent(out) :: n_a(:) !< Number of alpha coefficients in truncated wf.
        integer, allocatable, intent(out) :: n_b(:) !< Number of beta coefficients in truncated wf.
        integer :: i
        if (.not. allocated(n_a)) allocate(n_a(size(c_a, 2)))
        if (beta .and. (.not. allocated(n_b))) allocate(n_b(size(c_b, 2)))
        do i = 1, size(c_a, 2)
            if (beta) then
                call cis_nto_truncate_1state_u(trunc, c_a(:, i), c_b(:, i), n_a(i), n_b(i))
            else
                call cis_nto_truncate_1state_r(trunc, c_a(:, i), n_a(i))
            end if
        end do
    end subroutine cis_nto_truncate


end module cis_nto_mod
