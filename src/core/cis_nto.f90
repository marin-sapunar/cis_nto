module cis_nto_mod
    use global_defs

    implicit none

contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_nto_single
    !
    !> @brief Generate natural transition orbitals for a CIS wave function.
    !> @details
    !! Uses SVD decomposition of the CI coefficients matrix of a given spin to generate the natural
    !! transition orbitals for that spin. The NTOs are given as linear combinations of canonical 
    !! orbitals in the nat_orbs array. The array is filled with the no hole orbitals followed by 
    !! the first no particle orbitals.
    !----------------------------------------------------------------------------------------------
    subroutine cis_nto_single(nv, no, wf, nat_coef, nat_orbs)
        use lapack95, only : gesvd
        integer, intent(in) :: nv !< Number of virtual orbitals.
        integer, intent(in) :: no !< Number of occupied orbitals.
        real(dp), intent(in) :: wf(nv, no) !< CIS wave function coefficients.
        real(dp), intent(out) :: nat_coef(no) !< NTO coefficients.
        real(dp), intent(out) :: nat_orbs(no+nv, 2*no) !< NTOs in the MO basis.
        real(dp) :: wrk_c(nv, no) !< Array for CI coefficients.
        real(dp) :: wrk_o(no, no) !< Array for occupied NTOs.
        real(dp) :: wrk_v(nv, no) !< Array for virtual NTOs.
        wrk_c = wf(:, :)
        call gesvd(wrk_c, u=wrk_v, s=nat_coef, vt=wrk_o)
        nat_orbs(1:no, 1:no) = transpose(wrk_o)
        nat_orbs(no+1:no+nv, no+1:2*no) = wrk_v
    end subroutine cis_nto_single


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_nto
    !
    !> @brief Generate natural transition orbitals for a set of CIS wave functions.
    !> @details
    !! Uses SVD decomposition of the CI coefficients matrix of a given spin to generate the natural
    !! transition orbitals for that spin. The NTOs are given as linear combinations of canonical 
    !! orbitals in the nat_orbs array. The array is filled with the no hole orbitals followed by 
    !! the first no particle orbitals.
    !----------------------------------------------------------------------------------------------
    subroutine cis_nto(nv, no, ns, wf, nat_coef, nat_orbs)
        use lapack95, only : gesvd
        integer, intent(in) :: nv !< Number of virtual orbitals.
        integer, intent(in) :: no !< Number of occupied orbitals.
        integer, intent(in) :: ns !< Number of states.
        real(dp), intent(in) :: wf(nv, no, ns) !< CIS wave function coefficients.
        real(dp), intent(out) :: nat_coef(no, ns) !< NTO coefficients.
        real(dp), intent(out) :: nat_orbs(no+nv, 2*no, ns) !< NTOs in the MO basis.
        real(dp) :: wrk_c(nv, no) !< Array for CI coefficients.
        real(dp) :: wrk_o(no, no) !< Array for occupied NTOs.
        real(dp) :: wrk_v(nv, no) !< Array for virtual NTOs.
        integer :: i
        do i = 1, ns
            call cis_nto_single(nv, no, wf(:, :, i), nat_coef(:, i), nat_orbs(:, :, i))
        end do
    end subroutine cis_nto


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cis_nto_truncate_single
    !
    !> @brief Determine number of alpha/beta orbitals in truncated CIS wave function.
    !> @details
    !! Determine number of NTOs of each spin needed so the truncated CIS wave function has a square
    !! norm higher than trunc.
    !----------------------------------------------------------------------------------------------
    subroutine cis_nto_truncate_single(beta, trunc, c_a, c_b, n_a, n_b)
        logical, intent(in) :: beta !< Unrestricted calculation.
        real(dp), intent(in) :: trunc !< Threshold for truncating the wave functions.
        real(dp), intent(inout) :: c_a(:) !< NTO coefficients alpha.
        real(dp), intent(inout) :: c_b(:) !< NTO coefficients beta.
        integer, intent(out) :: n_a !< Number of alpha coefficients in truncated wave function.
        integer, intent(out) :: n_b !< Number of beta coefficients in truncated wave function.
        real(dp) :: norm !< Squared norm.
        real(dp) :: ia !< Current alpha excitation coefficient.
        real(dp) :: ib !< Current beta excitation coefficient.
        integer :: i

        n_a = 0
        if (beta) n_b = 0
        if (trunc < 1.0_dp) then
            norm = 0.0_dp
            if (beta) then
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
            else
                do i = 1, size(c_a, 1)
                    n_a = n_a + 1
                    norm = norm + c_a(n_a)**2
                    if (norm > trunc) exit
                end do
            end if
            ! Renormalize
          ! c_a(1:n_a) = c_a(1:n_a) / sqrt(norm)
          ! if (beta) c_b(1:n_b) = c_b(1:n_b) / sqrt(norm)
        else
            n_a = size(c_a, 1)
            if (beta) n_b = size(c_b, 1)
        end if
    end subroutine cis_nto_truncate_single



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
        integer, intent(out) :: n_a(:) !< Number of alpha coefficients in truncated wave function.
        integer, intent(out) :: n_b(:) !< Number of beta coefficients in truncated wave function.
        integer :: st
        do st = 1, size(c_a, 2)
            call cis_nto_truncate_single(beta, trunc, c_a(:, st), c_b(:, st), n_a(st), n_b(st))
        end do
    end subroutine cis_nto_truncate


end module cis_nto_mod
