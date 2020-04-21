!----------------------------------------------------------------------------------------------
! MODULE: write_txt_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @brief Contains subroutines for writing arrays to simple ASCII files.
!----------------------------------------------------------------------------------------------
module write_txt_mod
    use global_defs
    implicit none


    private
    public :: write_txt


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: write_txt
    !
    !> @brief Write a real(dp) array to a text file.
    !> @details
    !! The dimensions of the array are written to the first line of the file.
    !! The subsequent lines loop over dimensions of the array from the outermost to the second 
    !! dimension. The first dimension is written as a vector in each line.
    !----------------------------------------------------------------------------------------------
    interface write_txt
        module procedure write_txt_dim2
        module procedure write_txt_dim3
    end interface write_txt


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: write_txt_dim2
    !> @brief Write a 2-dimensional array to text file. See write_txt.
    !----------------------------------------------------------------------------------------------
    subroutine write_txt_dim2(fname, array)
        character(len=*), intent(in) :: fname
        real(dp), intent(in) :: array(:, :)
        integer :: i, ounit
        character(len=100) :: outfmt

        write(outfmt, '(a,i0,a)') '(', size(array, 1), 'e24.16)'

        open(newunit=ounit, file=fname, action='write')
        write(ounit, *) size(array, 1), size(array, 2)
        do i = 1, size(array, 2)
            write(ounit, outfmt) array(:, i)
        end do
        close(ounit)
    end subroutine write_txt_dim2


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: write_txt_dim3
    !> @brief Write a 3-dimensional array to text file. See write_txt.
    !----------------------------------------------------------------------------------------------
    subroutine write_txt_dim3(fname, array)
        character(len=*), intent(in) :: fname
        real(dp), intent(in) :: array(:, :, :)
        integer :: i, j, ounit
        character(len=100) :: outfmt

        write(outfmt, '(a,i0,a)') '(', size(array, 1), 'e24.16)'

        open(newunit=ounit, file=fname, action='write')
        write(ounit, *) size(array, 1), size(array, 2), size(array, 3)
        do i = 1, size(array, 3)
            do j = 1, size(array, 2)
                write(ounit, outfmt) array(:, j, i)
            end do
        end do
        close(ounit)
    end subroutine write_txt_dim3


end module write_txt_mod
