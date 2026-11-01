import h5py
import numpy as np
import matplotlib.pyplot as plt
from libra_py import units
import os

folders = { "LD": "FSSH_EasyWFO_LD", "NAC": "FSSH_Libra_NACs", "IC0": "FSSH_TUTORIAL_icond_0", "IC1": "FSSH_TUTORIAL_icond_1" }

def load_data(folder_name):
    filepath = f"{folder_name}/mem_data.hdf"
    if not os.path.exists(filepath): return None, None, None
    with h5py.File(filepath, 'r') as F:
        return np.array(F['time/data'][:]) * units.au2fs, np.array(F['se_pop_adi/data'][:]), np.array(F['sh_pop_adi/data'][:])

data = {k: load_data(v) for k, v in folders.items()}

fig, axs = plt.subplots(2, 2, figsize=(20, 16))
plt.rcParams.update({'font.size': 14})
colors = ['red', 'green', 'blue']

# TL: SE LD vs NAC
ax = axs[0, 0]
if data["LD"][0] is not None:
    for s in range(3): ax.plot(data["LD"][0], data["LD"][1][:, s], color=colors[s], linestyle='-', linewidth=3, label=f'LD State {s}')
if data["NAC"][0] is not None:
    for s in range(3): ax.plot(data["NAC"][0], data["NAC"][1][:, s], color=colors[s], linestyle='--', linewidth=3, label=f'NAC State {s}')
ax.set_title('SE Populations: EasyWFO (LD) vs Libra NACs', fontsize=20); ax.set_ylabel('Population', fontsize=18); ax.set_xlabel('Time (fs)', fontsize=18); ax.legend()

# TR: SH LD vs NAC
ax = axs[0, 1]
if data["LD"][0] is not None:
    for s in range(3): ax.plot(data["LD"][0], data["LD"][2][:, s], color=colors[s], linestyle='-', linewidth=3, label=f'LD State {s}')
if data["NAC"][0] is not None:
    for s in range(3): ax.plot(data["NAC"][0], data["NAC"][2][:, s], color=colors[s], linestyle='--', linewidth=3, label=f'NAC State {s}')
ax.set_title('SH Populations: EasyWFO (LD) vs Libra NACs', fontsize=20); ax.set_xlabel('Time (fs)', fontsize=18); ax.legend(ncol=2)

# BL: SE IC0 vs IC1
ax = axs[1, 0]
if data["IC0"][0] is not None:
    for s in range(3): ax.plot(data["IC0"][0], data["IC0"][1][:, s], color=colors[s], linestyle='-', linewidth=3, label=f'IC0 State {s}')
if data["IC1"][0] is not None:
    for s in range(3): ax.plot(data["IC1"][0], data["IC1"][1][:, s], color=colors[s], linestyle='--', linewidth=3, label=f'IC1 State {s}')
ax.set_title('SE Populations: IC 0 vs IC 1', fontsize=20); ax.set_ylabel('Population', fontsize=18); ax.set_xlabel('Time (fs)', fontsize=18); ax.legend(ncol=2)

# BR: SH IC0 vs IC1
ax = axs[1, 1]
if data["IC0"][0] is not None:
    for s in range(3): ax.plot(data["IC0"][0], data["IC0"][2][:, s], color=colors[s], linestyle='-', linewidth=3, label=f'IC0 State {s}')
if data["IC1"][0] is not None:
    for s in range(3): ax.plot(data["IC1"][0], data["IC1"][2][:, s], color=colors[s], linestyle='--', linewidth=3, label=f'IC1 State {s}')
ax.set_title('SH Populations: IC 0 vs IC 1', fontsize=20); ax.set_xlabel('Time (fs)', fontsize=18); ax.legend(ncol=2)

plt.tight_layout()
plt.savefig("all_populations_comparison.png", dpi=300)
print("Saved: all_populations_comparison.png")
