program test_fortran_trajectory
    use iso_c_binding
    implicit none
    
    interface
        subroutine compute_cisd_matrix_c(n_states, n_dets, max_occ, &
                                         alpha_bra, beta_bra, coeffs_bra, &
                                         alpha_ket, beta_ket, coeffs_ket, &
                                         S_MO, n_mo, matrix_out) bind(C, name="compute_cisd_matrix_c")
            import :: c_int, c_double
            integer(c_int) :: n_states, n_dets, max_occ, n_mo
            integer(c_int) :: alpha_bra(*), beta_bra(*), alpha_ket(*), beta_ket(*)
            real(c_double) :: coeffs_bra(*), coeffs_ket(*), S_MO(*), matrix_out(*)
        end subroutine
    end interface

    integer(c_int), parameter :: n_states = 3, n_dets = 3, max_occ = 2, n_mo = 3
    integer(c_int) :: alpha_bra(n_dets * max_occ), beta_bra(n_dets * max_occ)
    integer(c_int) :: alpha_ket(n_dets * max_occ), beta_ket(n_dets * max_occ)
    real(c_double) :: coeffs_bra(n_states * n_dets), coeffs_ket(n_states * n_dets)
    real(c_double) :: S_MO(n_mo * n_mo), matrix_out(n_states * n_states)
    integer :: step, i, j
    
    character(len=50) :: filename

    ! Determinants (0-based indices in C++)
    alpha_bra = [1, 2, 1, 3, 1, 3]; beta_bra = [1, 2, 1, 2, 1, 3]
    alpha_ket = alpha_bra; beta_ket = beta_bra

    print *, "=== Generating NAMD Overlap Trajectory for Libra ==="
    
    ! Loop over 50 time steps
    do step = 0, 49
        ! 1. Update MO overlaps (1D array, column-major for Eigen)
        S_MO(1) = 1.0; S_MO(2) = 0.1 - (step * 0.002); S_MO(3) = 0.0
        S_MO(4) = S_MO(2); S_MO(5) = 1.0; S_MO(6) = 0.05 - (step * 0.001)
        S_MO(7) = 0.0; S_MO(8) = S_MO(6); S_MO(9) = 1.0
        
        ! 2. Update CI coefficients (1D array, row-major for C++)
        ! State 1: [0.95, 0.20, 0.15]
        coeffs_bra(1) = 0.95; coeffs_bra(2) = 0.20; coeffs_bra(3) = 0.15
        ! State 2: [0.0, 0.98, 0.1]
        coeffs_bra(4) = 0.0;  coeffs_bra(5) = 0.98; coeffs_bra(6) = 0.1
        ! State 3: [0.0, 0.1, 0.99]
        coeffs_bra(7) = 0.0;  coeffs_bra(8) = 0.1;  coeffs_bra(9) = 0.99
        
        ! Slightly perturb ket (t+dt)
        coeffs_ket(1) = 0.949; coeffs_ket(2) = 0.205; coeffs_ket(3) = 0.149
        coeffs_ket(4) = 0.01;  coeffs_ket(5) = 0.979; coeffs_ket(6) = 0.11
        coeffs_ket(7) = 0.01;  coeffs_ket(8) = 0.11;  coeffs_ket(9) = 0.989

        ! 3. Call C++ Matrix Engine
        call compute_cisd_matrix_c(n_states, n_dets, max_occ, &
                                   alpha_bra, beta_bra, coeffs_bra, &
                                   alpha_ket, beta_ket, coeffs_ket, &
                                   S_MO, n_mo, matrix_out)

        ! 4. Write to st_adi_X.txt for Libra
        write(filename, '(A,I0,A)') "st_adi_", step, ".txt"
        open(unit=10, file=trim(filename), status="replace", action="write")
        do i = 1, n_states
            write(10, '(3F12.8)') (matrix_out((i-1)*n_states + j), j=1,n_states)
        end do
        close(10)
    end do
    
    print *, "SUCCESS! Generated 50 Libra-compatible st_adi_*.txt files."
end program test_fortran_trajectory
