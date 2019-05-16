!--------------------------------------------------------------------------------------------------
! MODULE: cgto_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date March, 2019
!
!> @brief Contains the cgto_subshell type and accompanying subroutines.
!--------------------------------------------------------------------------------------------------
module cgto_mod
    use global_defs
    use ang_mom_defs, only : amp
    implicit none
 
    private
    public :: cgto_subshell
    public :: cgto_one_el


    !----------------------------------------------------------------------------------------------
    ! TYPE: cgto_subshell
    !
    !> @brief Contracted Gaussian Type Orbital
    !----------------------------------------------------------------------------------------------
    type cgto_subshell
        integer :: l !< Angular momentum.
        character(len=1) :: typ !< Type of subshell.
        integer :: n_prim = 0 !< Number of primitive Gaussians.
        real(dp), allocatable :: z(:) !< Zeta coefficients. Dimensions: n_prim.
        real(dp), allocatable :: b(:) !< Beta coefficients. Dimensions: n_prim.
        real(dp), allocatable :: cb(:, :) !< Beta coefficients normalized for each cartesian function.
                                          !! Dimensions: n_prim x n_cart(l)
        logical :: init_cb = .false.
        logical :: inorm_b = .false.
    contains
        procedure :: norm_b => cgto_norm_b
        procedure :: ccgto_b
        procedure :: gen_all_ccgto_b
    end type cgto_subshell


contains


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: cgto_norm_b
    !
    !> @brief Normalize the cGTO.
    !> @details 
    !! Contraction coefficients are multiplied by the norm of the full contraction.
    !!   b'(i) = N * b(i)
    !! To fully normalize the basis set, the coefficients also need to be multiplied by the norms
    !! of the individual primitive (Cartesian) Gaussians.
    !----------------------------------------------------------------------------------------------
    subroutine cgto_norm_b(self)
        class(cgto_subshell) :: self
        integer:: i, j
        real(dp) :: bi, bj, zi, zj, pow
        real(dp):: nsum
        
        nsum = num0
        pow = f3o4 + self%l / num2
        do i = 1, self%n_prim
            bi = self%b(i)
            zi = self%z(i)
            do j = 1, self%n_prim
                bj = self%b(j)
                zj = self%z(j)
                nsum = nsum + bi * bj * zi**pow * zj**pow / (zi+zj)**(num2*pow)
            end do
        end do
        self%b = self%b / num2**pow / sqrt(nsum)
        self%inorm_b = .true.
    end subroutine cgto_norm_b


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: ccgto_b
    !> @brief Return cGTO coefficients multiplied by the ang. mom. dependent part of the norms.
    !----------------------------------------------------------------------------------------------
    function ccgto_b(self, lxyz) result(lxyz_b)
        use math_mod, only : factorial2
        class(cgto_subshell) :: self
        integer, intent(in) :: lxyz(3)
        real(dp) :: lxyz_b(size(self%b))
        integer :: lx, ly, lz, ll
        real(dp), parameter :: pi_3o2 = pi**(num3/num2)
        real(dp) :: pow, ggg, mul

        lx = lxyz(1)
        ly = lxyz(2)
        lz = lxyz(3)
        ll = sum(lxyz)
        pow = f3o4 + ll / num2
        ggg = factorial2(2*lx-1)*factorial2(2*ly-1)*factorial2(2*lz-1) / num2**ll * pi_3o2
        mul = 2**pow / sqrt(ggg)
        lxyz_b = self%b * mul * self%z**pow
    end function ccgto_b


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE: gen_all_ccgto_b
    !> @brief Return cGTO coefficients multiplied by the ang. mom. dependent part of the norms.
    !----------------------------------------------------------------------------------------------
    subroutine gen_all_ccgto_b(self)
        class(cgto_subshell) :: self
        integer :: i
        integer, allocatable :: wrk(:, :)

        if (.not. self%inorm_b) call self%norm_b
        if (allocated(self%cb)) deallocate(self%cb)
        allocate(self%cb(self%n_prim, amp%n_cart(self%l)))

        allocate(wrk(1:3,1:amp%n_cart(self%l)), source=amp%cart(self%l))
        do i = 1, amp%n_cart(self%l)
            self%cb(:, i) = self%ccgto_b(wrk(:, i))
        end do
        self%init_cb = .true.
    end subroutine gen_all_ccgto_b


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: cgto_one_el
    !> @brief Calculate one electron operator matrix between two cTGO subshells.
    !----------------------------------------------------------------------------------------------
    function cgto_one_el(cg1, cg2, q1, qc, q2, lc) result(op)
        type(cgto_subshell), intent(in) :: cg1
        type(cgto_subshell), intent(in) :: cg2
        real(dp), intent(in) :: q1(3)
        real(dp), intent(in) :: q2(3)
        integer, intent(in) :: lc(3)
        real(dp), allocatable :: op(:, :)
        real(dp) :: ssum, r2, z1, z2
        real(dp), parameter :: pi_3o2 = pi**(num3/num2)
        real(dp), parameter :: tinydp = 0.00000000000001_dp
        integer :: i, j, m, n
        real(dp) :: qc(3), os_int(3)
        integer, allocatable :: l1(:, :), l2(:, :)

        if (.not. cg1%init_cb) call cg1%gen_all_ccgto_b
        if (.not. cg2%init_cb) call cg2%gen_all_ccgto_b

        r2 = sum((q1-q2)**2)
        allocate(l1(1:3,1:amp%n_cart(cg1%l)), source=amp%cart(cg1%l))
        allocate(l2(1:3,1:amp%n_cart(cg2%l)), source=amp%cart(cg2%l))
        allocate(op(amp%n_cart(cg1%l), amp%n_cart(cg2%l)), source=num0)
        do i = 1, amp%n_cart(cg1%l)
            do j = 1, amp%n_cart(cg2%l)
                ssum = 0.0_dp
                do m = 1, cg1%n_prim
                    z1 = cg1%z(m)
                    do n = 1, cg2%n_prim
                        z2 = cg2%z(n)
                        os_int = osrec_multi(l1(:, i), lc, l2(:, j), q1, qc, q2, z1, z2)
                        ssum = ssum + cg1%cb(m, i) * cg2%cb(n, j) * product(os_int) * pi_3o2 * &
                        &      (z1+z2) ** (-num3/num2) * exp(-(z1*z2) * r2 / (z1+z2))
                    end do
                end do
                if (abs(ssum) > tinydp) op(i, j) = ssum
            end do
        end do
    end function cgto_one_el


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: osrec_multi
    !> @brief Calls osrec subroutine for each dimension
    !----------------------------------------------------------------------------------------------
    function osrec_multi(l1, lc, l2, q1, qc, q2, z1, z2) result(os_int)
        integer, intent(in) :: l1(:)
        integer, intent(in) :: lc(:)
        integer, intent(in) :: l2(:)
        real(dp), intent(in) :: q1(:)
        real(dp), intent(in) :: qc(:)
        real(dp), intent(in) :: q2(:)
        real(dp), intent(in) :: z1
        real(dp), intent(in) :: z2
        real(dp) :: os_int(size(l1))
        real(dp) :: qp, inv2p, inq(3)
        integer :: i, inl(3)

        inv2p = num1 / num2 / (z1 + z2)
        do i = 1, size(l1)
            inl = [l1(i), l2(i), lc(i)]
            inq = [q1(i), q2(i), qc(i)]
            qp = (z1*q1(i) + z2*q2(i)) / (z1+z2)
            os_int(i) = osrec(inl, inq, qp, inv2p)
        end do
    end function osrec_multi


    !----------------------------------------------------------------------------------------------
    ! FUNCTION: OSRec
    !
    !> @brief Obara-Saika recurrence relations for calculating integrals.
    !> @details
    !! The Obara-Saika recurrence relations for calculating the overlap and multipole integrals of
    !! Cartesian Gaussian functions of arbitrary quantum number:
    !! G(a,q,q0,lq) = (q-q0)^lq * exp(-a * (q-q0)^2)
    !! Integral: S(lq1,lqb,lq3) = < G(a,q,q0a,lqa) | (q-q03)^lq3 | G(b,q,q0b,lqb) >
    !! Recurrence relations:
    !!      S(i+1,j,k) = (qp-q0a) * S(i,j,k) + (i*S(i-1,j,k) + j*S(i,j-1,k) + k*S(i,j,k-1))/2(a+b)
    !!      S(i,j+1,k) = (qp-q0b) * S(i,j,k) + (i*S(i-1,j,k) + j*S(i,j-1,k) + k*S(i,j,k-1))/2(a+b)
    !!      S(i,j,k+1) = (qp-q03) * S(i,j,k) + (i*S(i-1,j,k) + j*S(i,j-1,k) + k*S(i,j,k-1))/2(a+b)
    !
    !> @param inl - The 3 quantum numbers: lqa, lqb and lq3.
    !> @param inq - Coordinates of the centres of the functions: q0a, q0b, q03
    !> @param qp - Position of the Gaussian overlap distribution: qp = (lqa*q0a + lqb*q0b)/(a+b)
    !> @param inv2p - 1/(2(a+b))
    !> @param sout - Result of applying the recurrence relations up to S(lqa,lqb,lq3).
    !!               For the value of the integral this should be multiplied by S(0,0,0).
    !----------------------------------------------------------------------------------------------
    function osrec(inl, inq, qp, inv2p) result(sout)
        integer,intent(in) :: inl(3)
        real(dp),intent(in) :: inq(3)
        real(dp),intent(in) :: qp
        real(dp),intent(in) :: inv2p
        real(dp) :: sout
        integer :: i, j, k, l1, l2, l3
        integer :: ord(3)
        real(dp) :: q1, q2, q3, qp1, qp2, qp3
        real(dp), allocatable :: s(:,:,:)

        if (inl(1)>=inl(2) .and. inl(1)>=inl(3)) then
            if (inl(2)>=inl(3)) then
                ord=[1,2,3]
            else
                ord=[1,3,2]
            end if
        else if (inl(2)>=inl(3)) then
            if (inl(1)>=inl(3)) then
                ord=[2,1,3]
            else
                ord=[2,3,1]
            end if
        else
            if (inl(1)>=inl(2)) then
                ord=[3,1,2]
            else
                ord=[3,2,1]
            end if
        end if
        l1=inl(ord(1))
        l2=inl(ord(2))
        l3=inl(ord(3))
        q1=inq(ord(1))
        q2=inq(ord(2))
        q3=inq(ord(3))

        allocate(s(0:l1,0:l2,0:l3))
        s = 0.0_dp
        qp1 = qp - q1
        qp2 = qp - q2
        qp3 = qp - q3
        s(0,0,0) = 1.0_dp
        if (l1>0) s(1,0,0) = qp1
        if (l2>0) then
            s(0,1,0) = qp2
            s(1,1,0) = qp1 * qp2 + inv2p
        end if
        do i=2,l1
            s(i,0,0) = qp1*s(i-1,0,0) + inv2p*(i-1)*s(i-2,0,0)
        end do
        if (l2>0) then
            do i=2,l1
                s(i,1,0) = qp2*s(i,0,0) + inv2p * i *s(i-1,0,0)
            end do
        end if
        do j=2,l2
            s(0,j,0) = qp1*s(0,j-1,0) + inv2p*(j-1)*s(0,j-2,0)
            s(1,j,0) = qp2*s(0,j,0)   + inv2p* j   *s(0,j-1,0)
        end do
        do j=2,l2
            do i=2,l1
                s(i,j,0) = qp2*s(i,j-1,0) + inv2p * (i*s(i-1,j-1,0) + (j-1)*s(i,j-2,0))
            end do
        end do
        do j=l2-l3+1,l2
            do i=l3-l2+j,l1
                s(i,j,1) = qp3*s(i,j,0) + inv2p * (i*s(i-1,j,0) + j*s(i,j-1,0))
            end do
        end do
        do k=2,l3
            do j=l2-l3+k,l2
                do i=l3-l2+j-k+1,l1
                    s(i,j,k)=qp3*s(i,j,k-1) + inv2p * (i*s(i-1,j,k-1) + j*s(i,j-1,k-1) + (k-1)*s(i,j,k-2))
                end do
            end do
        end do

        sout = s(l1,l2,l3)
        deallocate(s)
    end function osrec


end module cgto_mod
