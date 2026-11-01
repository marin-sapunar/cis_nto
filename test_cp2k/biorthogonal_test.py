#!/usr/bin/env python3
"""
Test script for the SVD-based biorthogonal transformation (EasyWFO).
Implements Löwdin orthogonalization (S^-1/2) to normalize raw MO overlaps,
then performs SVD to extract principal overlaps (λ).
"""
import numpy as np

def load_swf(filename):
    with open(filename, 'r') as f:
        n = int(f.readline().split()[0])
        M = np.zeros((n, n))
        for i in range(n):
            M[i] = [float(x) for x in f.readline().split()]
    return M

def main():
    S_MO = load_swf("s_wf")
    print(f"Loaded raw MO overlap matrix: {S_MO.shape[0]}x{S_MO.shape[1]}")
    
    # Löwdin Orthogonalization (Ref: Libra's get_Lowdin_general)
    vals, vecs = np.linalg.eigh(S_MO)
    vals_inv_sqrt = np.where(vals > 1e-12, 1.0 / np.sqrt(vals), 0.0)
    S_inv_sqrt = vecs @ np.diag(vals_inv_sqrt) @ vecs.T
    S_MO_ortho = S_inv_sqrt @ S_MO @ S_inv_sqrt
    
    # Biorthogonal SVD
    U, principal_overlaps, Vt = np.linalg.svd(S_MO_ortho)
    
    print("\nTop 5 Principal Overlaps (λ):")
    for i in range(min(5, len(principal_overlaps))):
        print(f"  Orbital {i+1}: {principal_overlaps[i]:.8f}")
        
    np.savetxt("principal_overlaps.txt", principal_overlaps, fmt="%.8e")
    print("\nSaved principal overlaps to principal_overlaps.txt")

if __name__ == "__main__":
    main()
