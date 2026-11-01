#include <ctime>
#include "easywfo_core.hpp"
#include <iostream>

// "extern C" tells the C++ compiler to use C-style naming conventions
// so Fortran can link to it.
extern "C" {
    
    // Wrapper for propagate_reprojection
    // Fortran passes pointers, so we use pointer arguments.
    void propagate_reprojection_c(
        int* n, 
        Complex* C_old_ptr, 
        Complex* H_ptr, 
        Complex* S_MCH_ptr, 
        double* dt, 
        Complex* C_new_ptr) 
    {
        // Map raw pointers to Eigen matrices/vectors without copying data
        Eigen::Map<VectorXcd> C_old(C_old_ptr, *n);
        Eigen::Map<MatrixXcd> H(H_ptr, *n, *n);
        Eigen::Map<MatrixXcd> S_MCH(S_MCH_ptr, *n, *n);
        
        // Call the C++ core
        VectorXcd C_new = EasyWFO::propagate_reprojection(C_old, H, S_MCH, *dt);
        
        // Map the result back to the Fortran pointer
        Eigen::Map<VectorXcd> out(C_new_ptr, *n);
        out = C_new;
    }
}
