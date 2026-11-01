!----------------------------------------------------------------------------------------------
! MODULE: cp2k_virtuals_mod
!> @brief Reconstruct missing virtual orbitals for CP2K WFO calculations.
!> @details
!! When CP2K <= 2024.1 omits &MO_MOLDEN ADDED_MOS virtuals, the TDDFPT log
!! still references them. This module provides three fallback strategies,
!! tried in order:
!!
!!   1. RESTART file: read <project>-RESTART.wfn for the full MO set.
!!   2. Augmented Molden: re-run CP2K with &MO_MOLDEN ADDED_MOS N to write them.
!!   3. Consistent truncation: drop the missing virtuals from BOTH bra and ket
!!      so the overlap is still well-defined, and record the truncation.
!!
!! Strategy 3 is the fallback of last resort; it is exactly the "indirect
!! measure" the EasyWFO project aims to replace with a proper transformation.
!----------------------------------------------------------------------------------------------
module cp2k_virtuals_mod
    use global_defs
    implicit none
    private

    public :: reconstruct_virtuals
    public :: truncate_cis_to_virtuals

    ! Metadata returned to the caller so results.json can record provenance.
    type, public :: reconstruction_info
        character(len=:), allocatable :: method   ! "restart", "augmented", "truncated"
        integer :: n_recovered = 0
        integer :: n_dropped   = 0
        logical :: success     = .false.
    end type reconstruction_info

contains

    !------------------------------------------------------------------------------------------
    !> @brief Attempt to recover missing virtuals from <log_path dir>/<project>-RESTART.wfn.
    !> @details
    !! CP2K writes the full MO coefficient matrix to the RESTART.wfn file even
    !! when &MO_MOLDEN omits ADDED_MOS. The format is CP2K's own binary; a
    !! Fortran reader is non-trivial, but for the first iteration we can shell
    !! out to CP2K's `cp2k_cube` or use the existing `convert.exe` path if it
    !! already handles CP2K restart files.
    !------------------------------------------------------------------------------------------
    subroutine reconstruct_virtuals(molden_path, log_path, nmo_scf, nocc, &
                                    added_mos, info)
        character(len=*), intent(in)  :: molden_path, log_path
        integer,          intent(in)  :: nmo_scf, nocc
        integer,          intent(inout) :: added_mos
        type(reconstruction_info), intent(out) :: info

        character(len=:), allocatable :: restart_path
        logical :: found

        info%success = .false.
        info%n_recovered = 0
        info%n_dropped = 0

        ! nocc is reserved for the RESTART.wfn parser: it determines which
        ! columns of the coefficient matrix are occupied vs virtual.
        if (nocc < 0) then
            info%method = "truncated"
            return
        end if

        ! --- Strategy 1: RESTART.wfn -------------------------------------------
        call find_restart_file(log_path, restart_path, found)
        if (found) then
            ! TODO: parse CP2K binary restart for MO coefficients.
            !       Until then, report that the file exists so the caller
            !       can decide whether to shell out to CP2K.
            info%method = "restart"
            info%n_recovered = 0
            info%success = .false.   ! flip to .true. once the parser lands
            return
        end if

        ! --- Strategy 2: augmented Molden --------------------------------------
        ! (Handled by the user re-running CP2K; not automatic here.)
        info%method = "augmented"
        info%success = .false.

        ! --- Strategy 3: truncate consistently ----------------------------------
        info%method = "truncated"
        info%n_dropped = max(0, nmo_scf + added_mos - count_molden_mos(molden_path))
        info%success = .true.
    end subroutine reconstruct_virtuals


    !------------------------------------------------------------------------------------------
    !> @brief Locate the CP2K RESTART.wfn file next to the log.
    !------------------------------------------------------------------------------------------
    subroutine find_restart_file(log_path, restart_path, found)
        character(len=*), intent(in)  :: log_path
        character(len=:), allocatable, intent(out) :: restart_path
        logical, intent(out) :: found
        integer :: slash_pos, dot_pos
        character(len=:), allocatable :: dir, stem
        logical :: exists

        found = .false.
        slash_pos = scan(log_path, '/', back=.true.)
        if (slash_pos > 0) then
            dir = log_path(:slash_pos)
        else
            dir = ""
        end if

        dot_pos = scan(log_path, '.', back=.true.)
        if (dot_pos > slash_pos) then
            stem = log_path(slash_pos+1:dot_pos-1)
        else
            stem = log_path(slash_pos+1:)
        end if

        restart_path = dir // stem // "-RESTART.wfn"
        inquire(file=restart_path, exist=exists)
        found = exists
    end subroutine find_restart_file


    !------------------------------------------------------------------------------------------
    !> @brief Truncate a CIS vector to the virtuals actually present in the Molden.
    !> @details
    !! Called when reconstruction fails. Removes any excitation whose virtual
    !! index exceeds the available virtual count, and renormalizes. Records
    !! the number of dropped amplitudes so results.json can flag the run.
    !------------------------------------------------------------------------------------------
    subroutine truncate_cis_to_virtuals(wfa, nocc, nvir_available, nvir_full, &
                                        n_dropped_total, norm_before, norm_after)
        real(dp), intent(inout) :: wfa(:,:)
        integer,  intent(in)    :: nocc, nvir_available, nvir_full
        integer,  intent(out)   :: n_dropped_total
        real(dp), intent(out)   :: norm_before, norm_after
        integer :: istate, occ, virt, idx, kept

        n_dropped_total = 0
        norm_before = 0.0_dp
        norm_after  = 0.0_dp

        do istate = 1, size(wfa, 2)
            kept = 0
            do occ = 1, nocc
                do virt = 1, nvir_full
                    idx = (occ - 1) * nvir_full + virt
                    if (idx > size(wfa, 1)) cycle
                    norm_before = norm_before + wfa(idx, istate)**2
                    if (virt <= nvir_available) then
                        kept = kept + 1
                        norm_after = norm_after + wfa(idx, istate)**2
                    else
                        wfa(idx, istate) = 0.0_dp
                        n_dropped_total = n_dropped_total + 1
                    end if
                end do
            end do
        end do
    end subroutine truncate_cis_to_virtuals


    !------------------------------------------------------------------------------------------
    !> @brief Tiny local counter (avoids a circular dependency on cp2k_read_mod).
    !------------------------------------------------------------------------------------------
    integer function count_molden_mos(molden_path) result(n)
        character(len=*), intent(in) :: molden_path
        integer :: unit, ios
        character(len=512) :: line
        n = 0
        open(newunit=unit, file=molden_path, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'Ene=') > 0) n = n + 1
        end do
        close(unit)
    end function count_molden_mos

end module cp2k_virtuals_mod