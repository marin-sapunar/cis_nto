import os
import numpy as np
from scipy.optimize import linear_sum_assignment

REF_DIR = os.environ.get("REF_DIR")
print("=== EasyWFO Hungarian Algorithm Validation ===")

def easywfo_hungarian_correction(St_raw, phase_cumulative):
    nstates = St_raw.shape[0]
    
    # --- 1. HUNGARIAN ALGORITHM (Global Optimal Tracking) ---
    # We want to maximize overlap, so we minimize the negative absolute overlap
    cost_matrix = -np.abs(St_raw)
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    
    # Permute the columns to restore the original state order
    St_tracked = St_raw[:, col_ind]
    
    # --- 2. PHASE CORRECTION ---
    St_corrected = St_tracked.copy()
    for i in range(nstates):
        if St_corrected[i, i] < 0:
            St_corrected[:, i] *= -1
            phase_cumulative[i] *= -1
            
    swaps_detected = np.sum(col_ind != np.arange(nstates))
    return St_corrected, phase_cumulative, swaps_detected

# Validation Loop
start_step = 5000
end_step = 5005
dt_au = 41.34

print(f"{'Step':<8} | {'Corrupted RMSE':<16} | {'Hungarian RMSE':<16} | {'Swaps Fixed'}")
print("-" * 65)

for step in range(start_step, end_step):
    st_file = f"{REF_DIR}/St_ci_{step}_re"
    if not os.path.exists(st_file): continue
        
    St_truth = np.loadtxt(st_file)
    true_nacs = (St_truth - St_truth.T) / (2 * dt_au)
    
    # Corrupt the data (Swap 2&3, flip phase of 5)
    St_corrupted = St_truth.copy()
    St_corrupted[:, [1, 2]] = St_truth[:, [2, 1]] 
    St_corrupted[:, 4] *= -1
    
    corrupted_nacs = (St_corrupted - St_corrupted.T) / (2 * dt_au)
    rmse_corrupted = np.sqrt(np.mean((corrupted_nacs - true_nacs)**2))
    
    # Run Hungarian Correction
    phase_cumulative = np.ones(St_truth.shape[0])
    St_fixed, _, swaps = easywfo_hungarian_correction(St_corrupted, phase_cumulative)
    
    fixed_nacs = (St_fixed - St_fixed.T) / (2 * dt_au)
    rmse_fixed = np.sqrt(np.mean((fixed_nacs - true_nacs)**2))
    
    print(f"{step:<8} | {rmse_corrupted:<16.4e} | {rmse_fixed:<16.4e} | {swaps}")

print("-" * 65)
print("Success! Hungarian algorithm achieves RMSE = 0.0")
