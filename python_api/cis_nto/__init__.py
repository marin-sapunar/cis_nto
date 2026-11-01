"""cis_nto Python API - overlap calculator for CI/MR-CISD wavefunctions."""
from .api import compute_overlap, compute_overlap_batch, compute_overlap_trajectory
from .io import read_hdf5, write_hdf5
from .libra import write_libra_format, read_libra_format
from .phase import correct_phase_matrix, verify_corrected
from .logging import get_logger, set_verbosity

__version__ = "1.1.0"
__all__ = [
    "compute_overlap", "compute_overlap_batch", "compute_overlap_trajectory",
    "read_hdf5", "write_hdf5", "write_libra_format", "read_libra_format",
    "correct_phase_matrix", "verify_corrected",
    "get_logger", "set_verbosity",
]
