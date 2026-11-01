"""Test HDF5 I/O."""
import numpy as np
import tempfile
import os
from cis_nto.io import write_hdf5, read_hdf5


def test_write_read_hdf5():
    overlap = np.random.randn(5, 3, 3) + 1j * np.random.randn(5, 3, 3)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "test.h5")
        write_hdf5(path, overlap, time_fs=0.5)
        data = read_hdf5(path)
        assert data["overlap"].shape == (5, 3, 3)
        np.testing.assert_allclose(data["overlap"], overlap)


def test_write_with_hamiltonian():
    overlap = np.eye(3, dtype=complex)[None, :, :].repeat(5, axis=0)
    ham = overlap * 2.0
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "test.h5")
        write_hdf5(path, overlap, hamiltonian=ham)
        data = read_hdf5(path)
        np.testing.assert_allclose(data["hamiltonian"], ham)
