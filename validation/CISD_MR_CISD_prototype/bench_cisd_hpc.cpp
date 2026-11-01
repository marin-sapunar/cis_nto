#include <iostream>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <random>
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
        for (int i : alpha_occ) key = (key << 8) | i; 
        return key; 
    }
    uint64_t get_beta_key() const { 
        uint64_t key = 0; 
        for (int i : beta_occ) key = (key << 8) | i; 
        return key; 
    }
};

int main(int argc, char** argv) {
    int n_states = 10;
    int n_dets = 500;
    int max_occ = 4; 
    int n_mo = 20;

    std::cout << "=== EasyWFO HPC Benchmark (OpenMP + Static Matrices) ===" << std::endl;
    std::cout << "Threads Available: " << omp_get_max_threads() << std::endl;
    std::cout << "Total Pairs: " << n_states * n_states * n_dets * n_dets << std::endl;

    MatrixXd S_MO = MatrixXd::Identity(n_mo, n_mo);
    for(int i=0; i<n_mo; ++i) for(int j=i+1; j<n_mo; ++j) S_MO(i,j) = S_MO(j,i) = ((double)rand()/RAND_MAX) * 0.1;

    std::vector<std::vector<Determinant>> states(n_states);
    for (int s = 0; s < n_states; ++s) {
        for (int d = 0; d < n_dets; ++d) {
            Determinant det;
            det.coefficient = ((double)rand()/RAND_MAX) * 0.1;
            det.alpha_occ = {0, 1}; det.beta_occ = {0, 1};
            while(det.alpha_occ.size() < max_occ) {
                int v = 2 + (rand() % (n_mo - 2));
                if (std::find(det.alpha_occ.begin(), det.alpha_occ.end(), v) == det.alpha_occ.end()) det.alpha_occ.push_back(v);
            }
            while(det.beta_occ.size() < max_occ) {
                int v = 2 + (rand() % (n_mo - 2));
                if (std::find(det.beta_occ.begin(), det.beta_occ.end(), v) == det.beta_occ.end()) det.beta_occ.push_back(v);
            }
            std::sort(det.alpha_occ.begin(), det.alpha_occ.end());
            std::sort(det.beta_occ.begin(), det.beta_occ.end());
            states[s].push_back(det);
        }
    }

    auto start = std::chrono::high_resolution_clock::now();

    MatrixXd overlap_matrix = MatrixXd::Zero(n_states, n_states);

    // PARALLELIZE OVER STATES (I)
    #pragma omp parallel for schedule(dynamic)
    for (int I = 0; I < n_states; ++I) {
        // Thread-local caches avoid data races
        std::unordered_map<uint64_t, double> a_cache, b_cache;
        
        for (int J = 0; J < n_states; ++J) {
            double total = 0.0;
            for (const auto& bra : states[I]) {
                for (const auto& ket : states[J]) {
                    
                    uint64_t a_key = (bra.get_alpha_key() << 32) | ket.get_alpha_key();
                    double det_a = 1.0;
                    auto it_a = a_cache.find(a_key);
                    if (it_a == a_cache.end()) {
                        int n = bra.alpha_occ.size();
                        // STATIC MATRIX: No heap allocation! Max size 4x4
                        Matrix4d sub = Matrix4d::Identity();
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
            overlap_matrix(I, J) = total;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;

    std::cout << "Execution Time: " << elapsed.count() << " ms" << std::endl;
    
    return 0;
}
