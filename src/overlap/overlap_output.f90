!--------------------------------------------------------------------------------------------------
! PROGRAM: overlap_output_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date March, 2019
!
!> @brief Subroutines for printing overlap program output.
!--------------------------------------------------------------------------------------------------
module overlap_output_mod
    use global_defs
    use write_txt_mod
    implicit none


    private
    public :: output_ao
    public :: output_mo
    public :: output_wf
    public :: output_orth
    public :: output_phase


    character(len=*), parameter :: out_fmt_s = '(5x,1000f10.6)'


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: output_ao
    !> @brief Print AO overlap matrix to file.
    !----------------------------------------------------------------------------------------------
    subroutine output_ao(outfile, s_ao)
        character(len=*), intent(in) :: outfile
        real(dp), intent(in) :: s_ao(:, :)
        if (outfile == 'None') return
        if (print_level >= 1) then
            write(stdout, '(1x,a)') 'Writing AO overlap matrix to '//outfile//'.'
        end if
        call write_txt(outfile, transpose(s_ao))
    end subroutine output_ao


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: output_mo
    !> @brief Print MO overlap matrix to file.
    !----------------------------------------------------------------------------------------------
    subroutine output_mo(outfile, s_mo_a, s_mo_b)
        character(len=*), intent(in) :: outfile
        real(dp), allocatable, intent(in) :: s_mo_a(:, :)
        real(dp), allocatable, intent(in) :: s_mo_b(:, :)
        if (outfile == 'None') return
        if (.not. allocated(s_mo_b)) then
            if (print_level >= 1) then
                write(stdout, '(1x,a)') 'Writing MO overlap matrix to '//outfile//'.'
            end if
            call write_txt(outfile, transpose(s_mo_a))
        else
            if (print_level >= 1) then
                write(stdout, '(1x,a)') 'Writing alpha MO overlap matrix to '//outfile//'_a.'
                write(stdout, '(1x,a)') 'Writing beta MO overlap matrix to '//outfile//'_b.'
            end if
            call write_txt(outfile//'_a', transpose(s_mo_a))
            call write_txt(outfile//'_b', transpose(s_mo_b))
        end if
    end subroutine output_mo


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: output_wf
    !> @brief Print WF overlap matrix to file.
    !----------------------------------------------------------------------------------------------
    subroutine output_wf(outfile, s_wf)
        character(len=*), intent(in) :: outfile
        real(dp), intent(in) :: s_wf(:, :)
        if (outfile == 'None') return
        if (print_level >= 1) then
            write(stdout, '(1x,a)') 'Writing WF overlap matrix to '//outfile//'.'
        end if
        call write_txt(outfile, transpose(s_wf))
    end subroutine output_wf


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: output_orth
    !> @brief Orthogonalize WF overlap matrix and print to stdout.
    !----------------------------------------------------------------------------------------------
    subroutine output_orth(orth_overlap, s_wf)
        use orthog_mod
        use matrix_mod
        logical, intent(in) :: orth_overlap
        real(dp), intent(inout) :: s_wf(:, :)
        real(dp), allocatable :: s_wf_orig(:, :)
        real(dp) :: angle
        integer :: i

        if (print_level >= 1) then
            write(stdout, *)
            write(stdout, '(1x,a)') 'Raw overlap matrix:'
            do i = lbound(s_wf, 1), ubound(s_wf, 1)
                write(stdout, out_fmt_s) s_wf(i, :)
            end do
        end if

        if (orth_overlap) then
            allocate(s_wf_orig, source=s_wf)
            call orthog_lowdin(s_wf)
            if (print_level >= 1) then
                write(stdout, *)
                write(stdout, '(1x,a)') 'Orthogonalized overlap matrix:'
                do i = lbound(s_wf, 1), ubound(s_wf, 1)
                    write(stdout, out_fmt_s) s_wf(i, :)
                end do
                angle = acos(sum(s_wf*s_wf_orig) / mat_norm(s_wf) / mat_norm(s_wf_orig))
                write(stdout, *)
                write(stdout, '(5x,a,f6.2)') 'Frobenius inner product angle: ', angle * 180 / 3.14159
                write(stdout, '(5x,a)') 'Norms of rows of raw overlap matrix: '
                write(stdout, '(9x,1000f10.4)') sum(s_wf_orig**2, 2)
                write(stdout, '(5x,a)') 'Norms of columns of raw overlap matrix: '
                write(stdout, '(9x,1000f10.4)') sum(s_wf_orig**2, 1)
            end if
        end if
    end subroutine output_orth


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: output_phase
    !> @brief Match phase of el. states in WF overlap matrix and print to stdout.
    !----------------------------------------------------------------------------------------------
    subroutine output_phase(match_phase, s_wf)
        use phase_mod
        use assignment_problem_mod
        logical, intent(in) :: match_phase
        real(dp), intent(inout) :: s_wf(:, :)
        integer, allocatable :: assigned_rows(:)
        integer, allocatable :: assigned_cols(:)
        integer :: i
        if (match_phase) then
            call phasematch_assigned_rotation(s_wf)
            if (print_level >= 1) then
                write(stdout, *) 
                write(stdout, '(1x,a)') 'Overlap matrix with phase matching between assigned bra/ket states:'
                do i = lbound(s_wf, 1), ubound(s_wf, 1)
                    write(stdout, out_fmt_s) s_wf(i, :)
                end do
            end if
            call assignment_problem(-s_wf**2, assigned_rows, assigned_cols)
            write(stdout, '(5x,a)') 'Assignment:'
            write(stdout, '(9x,a,1000(i0, 1x))') 'Rows: ', assigned_rows
            write(stdout, '(9x,a,1000(i0, 1x))') 'Cols: ', assigned_cols
        end if
    end subroutine output_phase


end module overlap_output_mod
