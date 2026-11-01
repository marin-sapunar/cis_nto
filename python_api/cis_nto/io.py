"""HDF5 I/O for cis_nto overlap data (Libra-compatible schema)."""
import h5py
import numpy as np
from .logging import get_logger

_log = get_logger("io")


def write_hdf5(path, overlap, hamiltonian=None, time_fs=0.5, metadata=None):
    """Write overlap (and optional Hamiltonian) to HDF5.

    Schema matches Libra's mem_data.hdf:
      /St_adi/data  (nsteps, nstates, nstates) complex128
      /Ham_adi/data (nsteps, nstates, nstates) complex128 (optional)
      /time/data    (nsteps,) float64 (atomic units)
    """
    overlap = np.atleast_3d(overlap).astype(np.complex128)
    nsteps = overlap.shape[0]

    with h5py.File(path, "w") as f:
        f.create_dataset("St_adi/data", data=overlap, dtype=np.complex128)
        if hamiltonian is not None:
            ham = np.atleast_3d(hamiltonian).astype(np.complex128)
            f.create_dataset("Ham_adi/data", data=ham, dtype=np.complex128)
        try:
            from libra_py import units
            dt_au = time_fs * units.fs2au
        except ImportError:
            dt_au = time_fs / 0.02418884326505
        f.create_dataset("time/data", data=np.arange(nsteps) * dt_au, dtype=np.float64)
        if metadata:
            for k, v in metadata.items():
                f.attrs[k] = v
    _log.debug(f"write_hdf5: {path} -> {overlap.shape}")


def read_hdf5(path):
    """Read overlap data from HDF5 (Libra or cis_nto format)."""
    result = {}
    with h5py.File(path, "r") as f:
        for key, h5path in [("overlap", "St_adi/data"),
                             ("hamiltonian", "Ham_adi/data"),
                             ("time_au", "time/data")]:
            if h5path in f:
                result[key] = np.array(f[h5path])
        result["attrs"] = dict(f.attrs)
    return result
