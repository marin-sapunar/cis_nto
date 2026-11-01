import numpy as np
import pytest
from cis_nto.nto import compress_singles_nto, print_nto_summary

def test_nto_compression():
    # Create a mock singles block (10 occupied x 20 virtual)
    # dominated by 2 large excitations
    c_singles = np.zeros((10, 20))
    c_singles[0, 0] = 0.95
    c_singles[1, 1] = 0.30
    
    # Add tiny noise
    np.random.seed(42) # Fix seed for reproducibility
    c_singles += np.random.rand(10, 20) * 0.001

    U, Vt, S = compress_singles_nto(c_singles, threshold=0.999)

    # Check that the norm is conserved to within the threshold
    norm_sq = np.sum(S**2)
    total_norm = np.sum(c_singles**2)
    assert np.isclose(norm_sq, total_norm, rtol=1e-3)
    
    # Check that we compressed the matrix (rank should be <= 3)
    assert len(S) <= 3
    assert U.shape == (10, len(S))
    assert Vt.shape == (len(S), 20)

def test_nto_print(capsys):
    c_singles = np.random.rand(5, 5)
    print_nto_summary(c_singles)
    captured = capsys.readouterr()
    assert "NTO Compression Summary" in captured.out
