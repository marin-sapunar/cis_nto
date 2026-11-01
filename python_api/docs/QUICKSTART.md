# Quick Start (30 seconds)

```python
import cis_nto
import numpy as np

print(f"Version: {cis_nto.__version__}")

overlap = np.eye(3, dtype=complex)[None, :, :].repeat(5, axis=0)
cis_nto.write_libra_format("/tmp/test.hdf", overlap, time_fs=0.5)

data = cis_nto.read_libra_format("/tmp/test.hdf")
print(f"Read: {data['overlap'].shape}")
print(f"Attrs: {data['attrs']}")
```
