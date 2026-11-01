# cis_nto Python API

Python-callable interface for the `cis_nto` overlap calculator.

## Quick Start

```python
import cis_nto
import numpy as np

# Write in Libra format
overlap = np.eye(3, dtype=complex)[None, :, :].repeat(5, axis=0)
cis_nto.write_libra_format("mem_data.hdf", overlap, time_fs=0.5)

# Read it back
data = cis_nto.read_libra_format("mem_data.h5")
print(data["overlap"].shape)  # (5, 3, 3)
```

## Installation

```bash
cd cis_nto/python_api
pip install -e .
```

## Usage

### Single overlap

```python
S = cis_nto.compute_overlap("geo1.molden:geo1.log", "geo2.molden:geo2.log",
                             format="cp2k", nstates=3)
```

### Batch processing

```python
cis_nto.compute_overlap_batch(trajectory_dir="results/", format="cp2k",
                           pattern="geo_*.molden", nstates=3,
                           output="overlaps.h5")
```

### Libra integration

```python
cis_nto.write_libra_format("mem_data.hdf", overlap, time_fs=0.5)
# Then in Libra:
#   import h5`y
#   St = h5py.File("mem_data.hdf")["St_adi/data"][:]
```

## API Reference

- `compute_overlap(wf_1, wf_2, format, nstates, method)` - single overlap
- `compute_overlap_batch(dir, format, pattern, nstates, output)` - batch -> HDF5
- `compute_overlap_trajectory(wf_list, format, nstates)` - list -> array
- `write_libra_format(path, overlap, hamiltonian, time_fs, label)` - write HDF5
- `read_libra_format(path)` - read HDF5

## Citation

- Sapunar et al. (2019) *JCTC* 15, 3461
- Shakiba & Akimov (2023) (reprojection matrix)
- Plasser et al. (2012) *JCP* 137, 22A514
