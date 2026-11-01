"""Test Libra schema adapter."""
import numpy as np
import tempfile
import os
from cis_nto.libra import write_libra_format, read_libra_format


def test_libra_roundtrip():
    overlap = np.random.randn(10, 4, 4) + 1j * np.random.randn(10, 4, 4)
    ham = overlap * 0.5
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "mem_data.hdf")
        write_libra_format(path, overlap, ham, time_fs=0.5, label="test")
        data = read_libra_format(path)
        np.testing.assert_allclose(data["overlap"], overlap)
        np.testing.assert_allclose(data["hamiltonian"], ham)
        assert data["attrs"]["source"] == "cis_nto"
        assert data["attrs"]["time_step_fs"] == 0.5
