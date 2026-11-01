#include <iostream>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>
#include <algorithm>

using namespace Eigen;

struct Determinant {
    std::vector<int> alpha_occ;
    std::vector<int> beta_occ;
    
    uint64_t get_alpha_key() const { 
        uint64_t key = 0; 
        for (int i : alpha_occ) key = key * 10000 + i; // FIXED
        return key; 
    }
};

int main() {
    int n_mo = 300; // Large molecule, > 255 orbitals!
    int n_dets = 1000;
    
    std::cout << "=== Stress Test: Unique Determinants with n_mo = " << n_mo << " ===" << std::endl;
    
    std::vector<Determinant> dets;
    for (int d = 0; d < n_dets; ++d) {
        Determinant det;
        // Generate exactly 1000 unique determinants with orbitals > 255
        det.alpha_occ = {0, 1, 256 + (d / 40), 256 + (d % 40)}; 
        dets.push_back(det);
    }

    std::unordered_map<uint64_t, int> test_cache;
    for (int i = 0; i < n_dets; ++i) {
        test_cache[dets[i].get_alpha_key()] = i;
    }

    std::cout << "Unique Determinants Inserted: " << n_dets << std::endl;
    std::cout << "Unique Cache Keys Generated  : " << test_cache.size() << std::endl;

    if (test_cache.size() == n_dets) {
        std::cout << "✅ PASS: No hash collisions for orbitals > 255!" << std::endl;
    } else {
        std::cout << "❌ FAIL: Hash collisions detected!" << std::endl;
    }
    return 0;
}
