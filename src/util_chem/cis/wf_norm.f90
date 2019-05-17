!--------------------------------------------------------------------------------------------------
! MODULE: wf_norm_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date October, 2017
!
!> @brief Normalize wave function stored as separate alpha and beta coefficients.
!
!> @todo Cleanup/move somewhere else.
!--------------------------------------------------------------------------------------------------
module wf_norm_mod
    use global_defs
    implicit none
    
    private
    public :: wfab_normalize


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE wfab_normalize
    !> @brief Normalize wave function stored as separate alpha and beta coefficients.
    !----------------------------------------------------------------------------------------------
    subroutine wfab_normalize(beta, wfa, wfb)
        logical, intent(in) :: beta
        real(dp) :: wfa(:, :)
        real(dp) :: wfb(:, :)
        real(dp) :: norm
        integer :: i

        do i = 1, size(wfa, 2)
            norm = sum(wfa(:, i) * wfa(:, i))
            if (beta) then
                norm = norm + sum(wfb(:, i) * wfb(:, i))
                wfb(:, i) = wfb(:, i) / sqrt(norm)
            end if
            wfa(:, i) = wfa(:, i) / sqrt(norm)
        end do
    end subroutine wfab_normalize


end module wf_norm_mod
