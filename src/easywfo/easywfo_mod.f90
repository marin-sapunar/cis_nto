! Fortran wrapper for the EasyWFO C++ core (Shakiba & Akimov 2023).
! Original implementation.
module easywfo_mod
    use, intrinsic :: iso_c_binding
    implicit none

    ! Define the kind for complex double precision
    integer, parameter :: dp = c_double

    interface
        subroutine propagate_reprojection_c(n, C_old, H, S, dt, C_new) bind(C, name="propagate_reprojection_c")
            import :: c_int, c_double, c_double_complex
            integer(c_int) :: n
            complex(c_double_complex) :: C_old(n), C_new(n), H(n,n), S(n,n)
            real(c_double) :: dt
        end subroutine
    end interface

contains

    ! A clean Fortran wrapper to call from cis_nto
    subroutine easywfo_propagate(n, C_old, H, S, dt, C_new)
        integer, intent(in) :: n
        complex(dp), intent(in) :: C_old(n), H(n,n), S(n,n)
        real(dp), intent(in) :: dt
        complex(dp), intent(out) :: C_new(n)

        integer(c_int) :: n_c
        
        n_c = n
        call propagate_reprojection_c(n_c, C_old, H, S, dt, C_new)

    end subroutine easywfo_propagate

end module easywfo_mod
