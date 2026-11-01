import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from libra_py import units
import glob
import os

def exp_funct(t, _x, beta):
    return np.exp(-np.power(t / _x, beta))

folders = sorted(glob.glob("FSSH_TUTORIAL_icond_*"))
if len(folders) == 0:
    print("ERROR: No FSSH_TUTORIAL_icond_* folders found.")
    exit()

ISTATE = 1
all_pops = []
time_fs = None

for folder in folders:
    filepath = f"{folder}/mem_data.hdf"
    if not os.path.exists(filepath): continue
        
    with h5py.File(filepath, 'r') as F:
        if time_fs is None:
            time_fs = np.array(F['time/data'][:]) * units.au2fs
        all_pops.append(np.array(F['sh_pop_adi/data'][:, ISTATE]))

all_pops = np.array(all_pops)
avg_pop = np.mean(all_pops, axis=0)

p0 = [1000.0, 0.1]
popt, _ = curve_fit(exp_funct, time_fs, avg_pop, p0=p0, bounds=([0.0, 0.0], [np.inf, 5.0]))
_tau, _beta = popt

residuals = avg_pop - exp_funct(time_fs, *popt)
r_squared = 1.0 - (np.sum(residuals**2) / np.sum((avg_pop - np.mean(avg_pop))**2))

fig, ax = plt.subplots(figsize=(10, 8))
for i in range(len(all_pops)):
    ax.plot(time_fs, all_pops[i], color='gray', alpha=0.15, linewidth=1)
ax.plot(time_fs, avg_pop, color='black', linewidth=4, label="Mg@C60 Ensemble Average")
ax.plot(time_fs, exp_funct(time_fs, _tau, _beta), color='red', linewidth=4, linestyle='--', label=f"Fit: Tau={_tau:.1f} fs, Beta={_beta:.2f}")

ax.set_xlabel("Time (fs)", fontsize=20)
ax.set_ylabel(f"Population (State {ISTATE})", fontsize=20)
ax.set_title("Mg@C60 Ensemble Population Decay", fontsize=24)
ax.legend(fontsize=16); ax.grid(True, alpha=0.3)

text_str = f"Tau = {_tau/1000:.3f} ps\nBeta = {_beta:.3f}\nR2 = {r_squared:.3f}"
ax.text(0.6, 0.8, text_str, transform=ax.transAxes, fontsize=16, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
plt.tight_layout()
plt.savefig("MgC60_ensemble_lifetime_fit.png", dpi=300)
print("Saved: MgC60_ensemble_lifetime_fit.png")
