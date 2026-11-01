// This file implements the Shakiba & Akimov (2023) reprojection method
// and SHARC-style wfthres truncation. Original implementation.
#ifndef EASYWFO_CORE_HPP
#define EASYWFO_CORE_HPP

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
#include <complex>
#include <iostream>
#include <algorithm>

using Complex = std::complex<double>;
using MatrixXcd = Eigen::MatrixXcd;
using VectorXcd = Eigen::VectorXcd;

class EasyWFO {
public:
    static MatrixXcd orthogonalized_T(const MatrixXcd& T, double tol = 1e-10, bool drop = false) {
        MatrixXcd S = T.adjoint() * T;
        Eigen::SelfAdjointEigenSolver<MatrixXcd> es(S);
        
        VectorXcd eigvals = es.eigenvalues();
        MatrixXcd eigvecs = es.eigenvectors();
        int n = T.rows();
        
        int rank = 0;
        for (int i = 0; i < n; ++i) {
            if (eigvals(i).real() > tol) rank++;
        }
        
        if (rank < n && !drop) {
            throw std::runtime_error("orthogonalized_T: T is rank-deficient. Set drop=true.");
        }
        
        MatrixXcd S_inv_sqrt = MatrixXcd::Zero(n, n);
        for (int i = 0; i < n; ++i) {
            if (!drop || eigvals(i).real() > tol) {
                double inv_sqrt = 1.0 / std::sqrt(eigvals(i).real());
                S_inv_sqrt += inv_sqrt * eigvecs.col(i) * eigvecs.col(i).adjoint();
            }
        }
        
        return T * S_inv_sqrt;
    }

    static VectorXcd truncate_wavefunction(VectorXcd coeffs, double wfthres = 0.99) {
        std::vector<std::pair<double, int>> mags;
        for (int i = 0; i < coeffs.size(); ++i) {
            mags.push_back({std::norm(coeffs(i)), i});
        }
        std::sort(mags.begin(), mags.end(), [](const auto& a, const auto& b) {
            return a.first > b.first;
        });
        
        double norm = 0.0;
        VectorXcd truncated = VectorXcd::Zero(coeffs.size());
        for (const auto& p : mags) {
            truncated(p.second) = coeffs(p.second);
            norm += p.first;
            if (norm >= wfthres) break;
        }
        return truncated;
    }

    static VectorXcd propagate_reprojection(
        VectorXcd C_old, 
        const MatrixXcd& H, 
        const MatrixXcd& S_MCH, 
        double dt) 
    {
        int n = C_old.size();
        MatrixXcd T_inv = S_MCH.inverse();
        MatrixXcd T = orthogonalized_T(T_inv, 1e-10, true);
        
        // Evaluate temporary matrix first to satisfy Eigen's lazy evaluation
        MatrixXcd M = -Complex(0, 1) * H * dt;
        MatrixXcd A = M.exp();
        
        VectorXcd C_new = T * A * C_old;
        return C_new;
    }
};
#endif
