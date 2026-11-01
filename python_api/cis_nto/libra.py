"""Adapter for Libra's mem_data.hdf schema."""
import numpy as np
from .io import write_hdf5, read_hdf5


def write_libra_format(path, overlap, hamiltonian=None, time_fs=0.5, label=""):
    """Write overlap data in Libra's mem_data.hdf schema."""
    metadata = {
        "source": "cis_nto",
        "schema": "libra_mem_data",
        "time_step_fs": time_fs,
        "label": label,
    }
    write_hdf5(path, overlap, hamiltonian, time_fs, metadata)
    print(f"Wrote Libra-compatible HDF5: {path}")
    print(f"  St_adi: {overlap.shape}")
    if hamiltonian is not None:
        print(f"  Ham_adi: {hamiltonian.shape}")


def read_libra_format(path):
    """Read cis_nto/Libra overlap data back from HDF5."""
    return read_hdf5(path)


def get_st_adi(path):
    return read_hdf5(path)["overlap"]


def get_ham_adi(path):
    return read_hdf5(path).get("hamiltonian")
