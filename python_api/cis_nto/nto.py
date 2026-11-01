import numpy as np

def compress_singles_nto(c_singles, threshold=0.999):
    """
    Compresses the singles block of a CI wavefunction into Natural Transition Orbitals (NTOs).
    
    Parameters
    ----------
    c_singles : ndarray, shape (n_occ, n_virt)
        The CIS singles coefficient matrix.
    threshold : float
        The cumulative squared norm threshold to retain (default 0.999).
        
    Returns
    -------
    U : ndarray, shape (n_occ, n_nto)
        The left singular vectors (occupied NTOs).
    Vt : ndarray, shape (n_nto, n_virt)
        The right singular vectors (virtual NTOs).
    S : ndarray, shape (n_nto,)
        The retained singular values.
    """
    # 1. Perform SVD: c_singles = U * diag(S) * Vt
    U, S, Vt = np.linalg.svd(c_singles, full_matrices=False)
    
    # 2. Determine how many NTO pairs are needed to reach the threshold
    norm_sq = S**2
    total_norm = np.sum(norm_sq)
    cumsum = np.cumsum(norm_sq) / total_norm
    
    # Find the index where cumsum exceeds the threshold
    n_keep = np.searchsorted(cumsum, threshold) + 1
    
    # 3. Truncate the matrices
    U_trunc = U[:, :n_keep]
    S_trunc = S[:n_keep]
    Vt_trunc = Vt[:n_keep, :]
    
    return U_trunc, Vt_trunc, S_trunc

def print_nto_summary(c_singles, threshold=0.999):
    """Utility function to print a summary of the NTO compression."""
    n_occ, n_virt = c_singles.shape
    U, Vt, S = compress_singles_nto(c_singles, threshold)
    
    print("=== NTO Compression Summary ===")
    print(f"Original Singles Determinants: {n_occ * n_virt}")
    print(f"NTO Pairs Kept: {len(S)}")
    print(f"Compression Factor: {(n_occ * n_virt) / len(S):.2f}x")
    print(f"Singular Values: {np.round(S, 4)}")
    
    norm_retained = np.sum(S**2) / np.sum(np.linalg.svd(c_singles, compute_uv=False)**2)
    print(f"Cumulative Norm Retained: {norm_retained * 100:.2f}%")
