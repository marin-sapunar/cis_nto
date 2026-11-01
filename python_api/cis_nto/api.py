"""Main public API for cis_nto."""
import numpy as np
from pathlib import Path
from .wrapper import run_overlap_calculator
from .phase import correct_phase_matrix
from .io import write_hdf5
from .logging import get_logger

_log = get_logger("api")


def compute_overlap(wf_1, wf_2, format="cp2k", nstates=3,
                    method="reprojection", binary=None, correct_phase=True):
    """Compute overlap matrix between two CI wavefunctions.

    Parameters
    ----------
    wf_1, wf_2 : str
        Paths to wavefunction files
    format : str
        "cp2k" | "molden" | "molpro" | "columbus"
    nstates : int
        Number of excited states
    method : str
        "reprojection" | "finite_difference" | "raw"
    binary : str, optional
        Path to cis_overlap.exe
    correct_phase : bool
        Apply Hungarian permutation + sign flip (default True)

    Returns
    -------
    S : np.ndarray (nstates, nstates) complex
    """
    _log.info(f"compute_overlap: {wf_1} | {wf_2}")
    _log.debug(f"  format={format}, nstates={nstates}, method={method}, correct_phase={correct_phase}")

    stdout, _ = run_overlap_calculator(wf_1, wf_2, format, nstates, method, binary)
    S = _parse_overlap_output(stdout, nstates)

    if correct_phase:
        S, _, _ = correct_phase_matrix(S[None, :, :], return_swaps=True)
        S = S[0]
        _log.debug("  phase correction applied")

    return S


def compute_overlap_trajectory(wf_list, format="cp2k", nstates=3,
                                method="reprojection", correct_phase=True):
    """Process a list of wavefunctions sequentially, with batch phase correction.

    Returns
    -------
    overlaps : np.ndarray (nsteps, nstates, nstates) complex
    """
    if len(wf_list) < 2:
        raise ValueError(f"Need >= 2 wavefunctions, got {len(wf_list)}")

    _log.info(f"compute_overlap_trajectory: {len(wf_list)} wavefunctions, {nstates} states")

    overlaps = []
    for i in range(len(wf_list) - 1):
        S = compute_overlap(wf_list[i], wf_list[i + 1],
                            format, nstates, method, correct_phase=False)
        overlaps.append(S)
    overlaps = np.array(overlaps)

    if correct_phase:
        overlaps, swaps, flips = correct_phase_matrix(overlaps, return_swaps=True)
        _log.info(f"  phase correction: {sum(swaps)} swaps, {sum(flips)} flips total")

    return overlaps


def compute_overlap_batch(trajectory_dir, format="cp2k",
                          pattern="st_adi_*.txt", nstates=3,
                          output="overlaps.hdf5", correct_phase=True):
    """Process N geometries, write to HDF5 in one call."""
    files = sorted(Path(trajectory_dir).glob(pattern))
    if len(files) < 2:
        raise ValueError(f"Need >= 2 files matching '{pattern}' in {trajectory_dir}")

    _log.info(f"compute_overlap_batch: {len(files)} files in {trajectory_dir}")

    wf_list = [str(f) for f in files]
    overlaps = compute_overlap_trajectory(wf_list, format, nstates,
                                           correct_phase=correct_phase)

    metadata = {
        "source": "cis_nto",
        "trajectory_dir": str(trajectory_dir),
        "pattern": pattern,
        "nsteps": len(overlaps),
        "nstates": nstates,
        "phase_corrected": correct_phase,
    }
    write_hdf5(output, overlaps, metadata=metadata)
    _log.info(f"  wrote {output}")
    return output


def _parse_overlap_output(stdout, nstates):
    """Parse cis_overlap.exe text output into a complex numpy array."""
    lines = [l.strip() for l in stdout.split("\n") if l.strip()]
    S = np.zeros((nstates, nstates), dtype=complex)
    for i, line in enumerate(lines[-nstates:]):
        for j, val in enumerate(line.split()[:nstates]):
            try:
                S[i, j] = complex(val)
            except ValueError:
                _log.warning(f"  cannot parse token '{val}' at line {i}")
    return S
