!----------------------------------------------------------------------------------------------
! MODULE: read_txt_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @brief Contains subroutines for reading arrays from simple ASCII files.
!----------------------------------------------------------------------------------------------
module read_txt_mod
    use global_defs
    implicit none


    private
    public :: read_txt


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_txt
    !
    !> @brief Read a real(dp) array from a text file.
    !> @details
    !! The first line of the file should contain the dimensions of the array.
    !! The subsequent lines loop over dimensions of the array from the outermost to the second 
    !! dimension. The first dimension is written as a vector in each line.
    !! If the array dimensions are given as input (have values different than zero on calling the 
    !! subroutine), an error is thrown if the dimensions given in the file don't match the input.
    !----------------------------------------------------------------------------------------------
    interface read_txt
        module procedure read_txt_dim2
        module procedure read_txt_dim3
    end interface read_txt


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_txt_dim2
    !> @brief Read a 2-dimensional array from a text file. See read_txt.
    !----------------------------------------------------------------------------------------------
    subroutine read_txt_dim2(fname, n1, n2, array)
        character(len=*), intent(in) :: fname
        integer, intent(inout) :: n1
        integer, intent(inout) :: n2
        real(dp), allocatable, intent(out) :: array(:, :)
        integer :: d1, d2, iunit, iotest

        call check_file_exists(fname)
        open(newunit=iunit, file=fname)
        read(iunit, *, iostat=iotest) d1, d2
        if (iotest /= 0) call read_error_dim(fname, 2, iotest)
        if ((n1 /= 0) .and. (n1 /= d1)) call dim_error(fname, 1, n1, d1)
        if ((n2 /= 0) .and. (n2 /= d2)) call dim_error(fname, 2, n2, d2)
        n1 = d1
        n2 = d2
        if (allocated(array)) deallocate(array)
        allocate(array(n1, n2))
        read(iunit, *, iostat=iotest) array
        if (iotest /= 0) call read_error_val(fname, iotest)
        close(iunit)
    end subroutine read_txt_dim2


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_txt_dim3
    !> @brief Read a 3-dimensional array from a text file. See read_txt.
    !----------------------------------------------------------------------------------------------
    subroutine read_txt_dim3(fname, n1, n2, n3, array)
        character(len=*), intent(in) :: fname
        integer, intent(inout) :: n1
        integer, intent(inout) :: n2
        integer, intent(inout) :: n3
        real(dp), allocatable, intent(out) :: array(:, :, :)
        integer :: d1, d2, d3, iunit, iotest

        call check_file_exists(fname)
        open(newunit=iunit, file=fname)
        read(iunit, *, iostat=iotest) d1, d2, d3
        if (iotest /= 0) call read_error_dim(fname, 3, iotest)
        if ((n1 /= 0) .and. (n1 /= d1)) call dim_error(fname, 1, n1, d1)
        if ((n2 /= 0) .and. (n2 /= d2)) call dim_error(fname, 2, n2, d2)
        if ((n3 /= 0) .and. (n3 /= d3)) call dim_error(fname, 3, n3, d3)
        n1 = d1
        n2 = d2
        n3 = d3
        if (allocated(array)) deallocate(array)
        allocate(array(n1, n2, n3))
        read(iunit, *, iostat=iotest) array
        if (iotest /= 0) call read_error_val(fname, iotest)
        close(iunit)
    end subroutine read_txt_dim3


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: check_file_exists
    !> @brief Print error message and abort program if file is not found.
    !----------------------------------------------------------------------------------------------
    subroutine check_file_exists(file_name, error_message)
        character(len=*), intent(in) :: file_name
        character(len=*), intent(in), optional :: error_message
        logical :: chk

        inquire(file=file_name, exist=chk)
        if (.not. chk) then
            write(stderr, *) ' Failed to find required file: '//file_name
            if (present(error_message)) write(stderr, *) error_message
            call abort()
        end if
    end subroutine check_file_exists


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: dim_error
    !> @brief Exit the program with an error if an array dimension doesn't match expected value.
    !----------------------------------------------------------------------------------------------
    subroutine dim_error(fname, cdim, n, d)
        character(len=*), intent(in) :: fname
        integer, intent(in) :: cdim
        integer, intent(in) :: n
        integer, intent(in) :: d

        write(stderr, '(3a)') 'Error reading array from file ', trim(fname), '.'
        write(stderr, '(1x, a, i0)') 'For array dimension: ', cdim
        write(stderr, '(3x, a, i0, a, i0, a)') 'Expected size ', n, ' instead of ', d, '.'
        call abort()
    end subroutine dim_error


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_error_dim
    !> @brief Exit the program with an error if array dimensions can't be read from the first line.
    !----------------------------------------------------------------------------------------------
    subroutine read_error_dim(fname, cdim, io)
        character(len=*), intent(in) :: fname
        integer, intent(in) :: cdim
        integer, intent(in) :: io

        write(stderr, '(3a)') 'Error reading array dimensions from file ', trim(fname), '.'
        write(stderr, '(1x, a, i0, a)') 'First line should contain ', cdim, &
        &                               ' integers with array dimensions.'
        write(stderr, '(1x,a,i0)') 'iostat: ', io
        call abort()
    end subroutine read_error_dim


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: read_error_val
    !> @brief Exit the program with an error if array dimensions can't be read from the first line.
    !----------------------------------------------------------------------------------------------
    subroutine read_error_val(fname, io)
        character(len=*), intent(in) :: fname
        integer, intent(in) :: io

        write(stderr, '(3a)') 'Error reading array values from file ', trim(fname), '.'
        write(stderr, '(1x,a,i0)') 'iostat: ', io
        call abort()
    end subroutine read_error_val


end module read_txt_mod
