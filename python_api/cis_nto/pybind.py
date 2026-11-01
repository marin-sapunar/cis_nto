"""Optional pybind11 acceleration."""
try:
    from ._core import orthogonalized_T, truncate_wavefunction, propagate_reprojection
    HAS_PYBIND = True
except ImportError:
    HAS_PYBIND = False
    orthogonalized_T = truncate_wavefunction = propagate_reprojection = None
