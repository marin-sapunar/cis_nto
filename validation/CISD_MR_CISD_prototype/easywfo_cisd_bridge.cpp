#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <Eigen/Dense>
#include <algorithm>
#include <omp.h>

using namespace Eigen;

struct Determinant {
    std::vector<int> alpha_occ;
    std::vector<int> beta_occ;
    double coefficient;
    
    uint64_t get_alpha_key() const { 
        uint64_t key = 0; 
        for (int i : alpha_occ) key = key * 10000 + i; 
        return key; 
    }
    uint64_t get_beta_key() const { 
        uint64_t key = 0; 
        for (int i : beta_occ) key = key * 10000 + i; 
        return key; 
    }
};

extern "C" {
    void compute_cisd_matrix_c(int* n_states, int* n_dets, int* max_occ, 
                               int* alpha_idx_bra, int* beta_idx_bra, double* coeffs_bra,
                               int* alpha_idx_ket, int* beta_idx_ket, double* coeffs_ket,
                               double* S_MO_ptr, int* n_mo, double* matrix_out) {
        
        Map<MatrixXd> S_MO(S_MO_ptr, *n_mo, *n_mo);
        
        std::vector<std::vector<Determinant>> bra_states(*n_states);
        std::vector<std::vector<Determinant>> ket_states(*n_states);
        
        for (int s = 0; s < *n_states; ++s) {
            for (int i = 0; i < *n_dets; ++i) {
                Determinant d_bra, d_ket;
                d_bra.coefficient = coeffs_bra[s * (*n_dets) + i];
                d_ket.coefficient = coeffs_ket[s * (*n_dets) + i];
                
                for (int j = 0; j < *max_occ; ++j) {
                    int idx = i * (*max_occ) + j;
                    if (alpha_idx_bra[idx] > 0) d_bra.alpha_occ.push_back(alpha_idx_bra[idx] - 1);
                    if (beta_idx_bra[idx] > 0)  d_bra.beta_occ.push_back(beta_idx_bra[idx] - 1);
                    if (alpha_idx_ket[idx] > 0) d_ket.alpha_occ.push_back(alpha_idx_ket[idx] - 1);
                    if (beta_idx_ket[idx] > 0)  d_ket.beta_occ.push_back(beta_idx_ket[idx] - 1);
                }
                bra_states[s].push_back(d_bra);
                ket_states[s].push_back(d_ket);
            }
        }

        // PARALLELIZE OVER STATES
        #pragma omp parallel for schedule(dynamic)
        for (int I = 0; I < *n_states; ++I) {
            // Thread-local caches to avoid data races
            std::unordered_map<uint64_t, double> a_cache, b_cache;
            
            for (int J = 0; J < *n_states; ++J) {
                double total = 0.0;
                for (const auto& bra : bra_states[I]) {
                    for (const auto& ket : ket_states[J]) {
                        
                        uint64_t a_key = (bra.get_alpha_key() << 32) | ket.get_alpha_key();
                        double det_a = 1.0;
                        auto it_a = a_cache.find(a_key);
                        if (it_a == a_cache.end()) {
                            int n = bra.alpha_occ.size();
                            Matrix4d sub = Matrix4d::Identity(); // Static stack allocation
                            for(int r=0; r<n; ++r) for(int c=0; c<n; ++c) sub(r,c) = S_MO(bra.alpha_occ[r], ket.alpha_occ[c]);
                            det_a = sub.topLeftCorner(n, n).determinant();
                            a_cache[a_key] = det_a;
                        } else {
                            det_a = it_a->second;
                        }
                        
                        uint64_t b_key = (bra.get_beta_key() << 32) | ket.get_beta_key();
                        double det_b = 1.0;
                        auto it_b = b_cache.find(b_key);
                        if (it_b == b_cache.end()) {
                            int n = bra.beta_occ.size();
                            Matrix4d sub = Matrix4d::Identity();
                            for(int r=0; r<n; ++r) for(int c=0; c<n; ++c) sub(r,c) = S_MO(bra.beta_occ[r], ket.beta_occ[c]);
                            det_b = sub.topLeftCorner(n, n).determinant();
                            b_cache[b_key] = det_b;
                        } else {
                            det_b = it_b->second;
                        }
                        
                        total += bra.coefficient * ket.coefficient * (det_a * det_b);
                    }
                }
                matrix_out[I * (*n_states) + J] = total;
            }
        }
    }
}
