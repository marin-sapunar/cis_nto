!----------------------------------------------------------------------------------------------
! MODULE: truncate_wf_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date November, 2018
!
!> @brief Subroutines for truncating a wave function.
!----------------------------------------------------------------------------------------------
module truncate_wf_mod
    use global_defs
    implicit none


    private
    public :: truncate_wf


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: truncate_wf
    !
    !> @brief Truncate a set of wave functions.
    !> @details
    !! Determine minimum set of coefficients for each state for the norm to surpass threshold.
    !! Positions of the coefficients are returned as logical arrays for masking the initial wave
    !! function coefficients. The amask_any and bmask_any arrays are returned as .true. for 
    !! coefficients which are selected in any state.
    !----------------------------------------------------------------------------------------------
    subroutine truncate_wf(thr, beta, wfa, wfb, amask, bmask)
        real(dp), intent(in) :: thr !< Threshold for truncating the wave function
        logical, intent(in) :: beta !< Restricted/unrestricted calculation.
        real(dp), intent(in) :: wfa(:, :) !< Wave function coefficients alpha.
        real(dp), intent(in) :: wfb(:, :) !< Wave function coefficients beta.
        logical, allocatable, intent(out) :: amask(:, :) !< Selected alpha coeffs. in each state.
        logical, allocatable, intent(out) :: bmask(:, :) !< Selected beta coeffs. in each state.
        real(dp), allocatable :: vec(:)
        logical, allocatable :: mask(:)
        integer :: i, dima, dimb, dimt, nwf

        nwf = size(wfa, 2)
        dima = size(wfa, 1)
        if (allocated(amask)) deallocate(amask)
        allocate(amask(dima, nwf))
        dimb = 0
        if (beta) then
            dimb = size(wfb, 1)
            if (allocated(bmask)) deallocate(bmask)
            allocate(bmask(dimb, nwf))
        end if
        dimt = dima + dimb
        allocate(vec(dimt))
        allocate(mask(dimt))

        do i = 1, nwf
            vec(1:dima) = wfa(:, i)
            if (beta) vec(dima+1:dimt) = wfb(:, i)
            call truncate_vector_norm(thr, vec, mask)
            amask(:, i) = mask(1:dima)
            if (beta) then
                bmask(:, i) = mask(dima+1:dimt)
            end if
        end do
    end subroutine truncate_wf


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: truncate_vector_norm
    !
    !> @brief Truncate a vector so it's norm is equal to or larger than a threshold.
    !> @details
    !! Determine minimum number of coefficients needed for norm of vector to surpass threshold.
    !! Positions of the coefficients are returned as a logical array for masking the initial 
    !! vector.
    !----------------------------------------------------------------------------------------------
    subroutine truncate_vector_norm(thr, vec, mask)
        real(dp), intent(in) :: thr !< Threshold for truncating the vector.
        real(dp), intent(in) :: vec(:) !< Vector coefficients.
        logical, intent(out) :: mask(:) !< Mask for selected coefficients.
        real(dp), allocatable :: v2(:) !< Squared coefficients.
        integer :: p(1) !< Position of current largest coefficient.
        integer :: i
        real(dp) :: norm

        allocate(v2(size(vec)), source=vec**2)
        mask = .true.
        norm = 0.0_dp

        do i = 1, size(vec)
            p = maxloc(v2, mask)
            norm = norm + v2(p(1)) 
            mask(p(1)) = .false.
            if (norm >= thr) exit
        end do
        mask = .not. mask
    end subroutine truncate_vector_norm


end module truncate_wf_mod
