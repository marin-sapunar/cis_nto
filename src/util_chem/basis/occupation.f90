!--------------------------------------------------------------------------------------------------
! MODULE: occupation_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date August, 2018
!
!> @details Different ways of storing occupation numbers.
!> @todo Cleanup
!--------------------------------------------------------------------------------------------------
module occupation_mod
    use global_defs, only : dp
    implicit none
 
    private
    public :: occupation_numbers
    public :: occ_mask_to_nums
    public :: occ_nums_to_mask
    public :: occ_nums_to_vec
    public :: occ_mask_to_vec
    public :: mo_reshape


    character(len=1), parameter :: spin_char(0:3) = ['e', 'a', 'b', 'd'] !< MO occupation labels.


    !----------------------------------------------------------------------------------------------
    ! TYPE: occupation_numbers
    !> @brief Hold numbers of occupied/virtual active/frozen orbitals
    !----------------------------------------------------------------------------------------------
    type occupation_numbers
        integer :: rhf = 0 !< Restricted/unrestricted mos.
        integer :: n = 0 !< Number of mos (same for both spins).
        integer :: o(2) = 0 !< Number of occupied mos.
        integer :: v(2) = 0 !< Number of virtual mos.
        integer :: a(2) = 0 !< Number of active mos.
        integer :: f(2) = 0 !< Number of frozen mos.
        integer :: ao(2) = 0 !< Number of active occupied mos.
        integer :: av(2) = 0 !< Number of active virtual mos.
        integer :: fo(2) = 0 !< Number of frozen occupied mos.
        integer :: fv(2) = 0 !< Number of frozen virtual mos.
    end type occupation_numbers


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: occ_mask_to_nums
    !> @brief Count numbers of active/occupied orbitals in logical arrays.
    !----------------------------------------------------------------------------------------------
    subroutine occ_mask_to_nums(rhf, omask, amask, onum)
        integer, intent(in) :: rhf
        logical, intent(in) :: omask(:, :)
        logical, intent(in), optional :: amask(:, :)
        type(occupation_numbers), intent(out) :: onum

        onum%rhf = rhf
        onum%n = size(omask, 1)
        onum%o = count(omask, 1)
        onum%v = count(.not. omask, 1)
        if (present(amask)) then
            onum%a = count(amask, 1)
            onum%f = count(.not. amask, 1)
            onum%ao = count(omask .and. amask, 1)
            onum%av = count((.not. omask) .and. amask, 1)
            onum%fo = count(omask .and. (.not. amask), 1)
            onum%fv = count(((.not. omask) .and. (.not. amask)), 1)
        else
            onum%a = onum%n
            onum%ao = onum%o
            onum%av = onum%v
        end if
    end subroutine occ_mask_to_nums


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: occ_nums_to_mask
    !
    !> @brief Create logical arrays from occupation numbers.
    !> @details
    !! Order of the orbitals in the generated arrays: 1) frozen occupied, 2) active occupied, 
    !! 3) active virtual and 4) frozen virtual.
    !----------------------------------------------------------------------------------------------
    subroutine occ_nums_to_mask(onum, omask, amask)
        type(occupation_numbers), intent(in) :: onum
        logical, allocatable, intent(out) :: omask(:, :)
        logical, allocatable, intent(out), optional :: amask(:, :)
        integer :: i


        if (allocated(omask)) deallocate(omask)
        allocate(omask(onum%n, 2), source = .false.)
        do i = 1, 2
            omask(1:onum%o(i), i) = .true.
        end do
        if (present(amask)) then
            if (allocated(amask)) deallocate(amask)
            allocate(omask(onum%n, 2), source = .true.)
            do i = 1, 2
                if (onum%fo(i) > 0) amask(1:onum%fo(i), i) = .false.
                if (onum%fv(i) > 0) amask(onum%o(i)+onum%av(i)+1:, i) = .false.
            end do
        end if
    end subroutine occ_nums_to_mask


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: occ_nums_to_vec
    !> @brief Create integer vector of occupation numbers from occupation_numbers instance.
    !----------------------------------------------------------------------------------------------
    subroutine occ_nums_to_vec(onum, ovec)
        type(occupation_numbers), intent(in) :: onum
        integer, allocatable, intent(out) :: ovec(:)

        if (allocated(ovec)) deallocate(ovec)
        allocate(ovec(onum%n), source = 0)
        ovec(1:onum%o(1)) = ovec(1:onum%o(1)) + 1
        ovec(1:onum%o(2)) = ovec(1:onum%o(2)) + 2
    end subroutine occ_nums_to_vec


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: occ_mask_to_vec
    !> @brief Create integer vector of occupation numbers from logical array.
    !----------------------------------------------------------------------------------------------
    subroutine occ_mask_to_vec(omask, ovec)
        logical, intent(in) :: omask(:, :)
        integer, allocatable, intent(out) :: ovec(:)

        if (allocated(ovec)) deallocate(ovec)
        allocate(ovec(size(omask, 1)), source = 0)
        where (omask(:, 1)) ovec = ovec + 1
        where (omask(:, 2)) ovec = ovec + 2
    end subroutine occ_mask_to_vec


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: ovec_move_electron
    !> @brief Move single electron in occupation numbers vector.
    !----------------------------------------------------------------------------------------------
    function ovec_move_electron(ref, s, o, v) result(new)
        integer, intent(in) :: ref(:)
        integer, intent(in) :: s
        integer, intent(in) :: o
        integer, intent(in) :: v
        integer, allocatable :: new(:)
        allocate(new, source = ref)
        new(o) = new(o) - s
        new(v) = new(v) + s
    end function ovec_move_electron


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: sort_occ_act
    !
    !> @brief Get reordered vector of orbital indexes.
    !> @details
    !! Target order of the orbitals: 1) frozen occupied, 2) active occupied, 3) active virtual and 
    !! 4) frozen virtual.
    !----------------------------------------------------------------------------------------------
    subroutine sort_occ_act(occ, act, index_vector)
        logical, intent(in) :: occ(:)
        logical, intent(in) :: act(:)
        integer, allocatable, intent(inout) :: index_vector(:)
        integer :: seq(size(occ))
        logical :: mask(size(occ))
        integer :: i, n, c

        if (allocated(index_vector)) deallocate(index_vector)
        allocate(index_vector(size(occ)))
        seq = [ (i, i=1, size(occ)) ]
        c = 0
        ! Put inactive occupied orbitals first.
        mask = (occ .and. (.not. act))
        n = count(mask)
        index_vector(c+1:c+n) = pack(seq, mask)
        c = c + n
        ! Put active occupied orbitals second.
        mask = (occ .and. act)
        n = count(mask)
        index_vector(c+1:c+n) = pack(seq, mask)
        c = c + n
        ! Put active unoccupied orbitals third.
        mask = ((.not. occ) .and. act)
        n = count(mask)
        index_vector(c+1:c+n) = pack(seq, mask)
        c = c + n
        ! Put inactive unoccupied orbitals last.
        mask = ((.not. occ) .and. (.not. act))
        n = count(mask)
        index_vector(c+1:c+n) = pack(seq, mask)
        c = c + n
    end subroutine sort_occ_act


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: remove_inact
    !> @brief Remove inactive orbitals from orbital index vector.
    !----------------------------------------------------------------------------------------------
    subroutine remove_inact(act, index_vector)
        logical, intent(in) :: act(:)
        integer, allocatable, intent(inout) :: index_vector(:)
        integer, allocatable :: wrk(:)
        allocate(wrk, source=index_vector)
        deallocate(index_vector)
        allocate(index_vector, source=pack(wrk, act))
    end subroutine remove_inact


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: mo_reshape
    !
    !> @brief Reshape MO arrays.
    !> @details
    !! Reorder and/or remove inactive orbitals from MO arrays. For 2 or 3 dimensional arrays, the
    !! mo_dim variable tells the subroutine which dimension of the array corresponds to the MOs.
    !----------------------------------------------------------------------------------------------
    subroutine mo_reshape(sort, remove, occ, act, mo_dim, a1, a2, a3)
        logical :: sort
        logical :: remove
        logical, intent(in), optional :: occ(:)
        logical, intent(in) :: act(:)
        integer, intent(in) :: mo_dim
        real(dp), allocatable, optional :: a1(:)
        real(dp), allocatable, optional :: a2(:, :)
        real(dp), allocatable, optional :: a3(:, :, :)
        integer, allocatable :: sel(:)

        integer :: seq(size(act))
        integer, allocatable :: index_vector(:)
        real(dp), allocatable :: wrk_a1(:)
        real(dp), allocatable :: wrk_a2(:, :)
        real(dp), allocatable :: wrk_a3(:, :, :)
        integer :: i

        if (sort) then
            call sort_occ_act(occ, act, index_vector)
        else
            allocate(index_vector(size(act)))
            index_vector = [ (i, i=1, size(act)) ]
        end if
        if (remove) then
            call remove_inact(act, index_vector)
        end if

        if (present(a1)) then
            allocate(wrk_a1, source=a1(index_vector))
            deallocate(a1)
            allocate(a1(1:size(wrk_a1, 1)), source=wrk_a1)
        end if
        if (present(a2)) then
            if (mo_dim == 1) then
                allocate(wrk_a2, source=a2(index_vector, :))
            else if (mo_dim == 2) then
                allocate(wrk_a2, source=a2(:, index_vector))
            end if
            deallocate(a2)
            allocate(a2(1:size(wrk_a2, 1), 1:size(wrk_a2, 2)), source=wrk_a2)
        end if
        if (present(a3)) then
            if (mo_dim == 1) then
                allocate(wrk_a3, source=a3(index_vector, :, :))
            else if (mo_dim == 2) then
                allocate(wrk_a3, source=a3(:, index_vector, :))
            else if (mo_dim == 3) then
                allocate(wrk_a3, source=a3(:, :, index_vector))
            end if
            deallocate(a3)
            allocate(a3(1:size(wrk_a3, 1), 1:size(wrk_a3, 2), 1:size(wrk_a3, 3)), source=wrk_a3)
        end if
    end subroutine mo_reshape


end module occupation_mod
