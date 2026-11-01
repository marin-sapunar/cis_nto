import os
import numpy as np

REF_DIR = os.environ.get("REF_DIR")
print("=== EasyWFO Master Black-Box Validation ===")
print(f"Ground Truth Dir: {REF_DIR}\n")

# ---------------------------------------------------------
# YOUR COMBINED EASYWFO ALGORITHM
# ---------------------------------------------------------
def easywfo_black_box_correction(St_raw, phase_cumulative):
    """
    1. Tracks states (fixes crossings/permutations).
    2. Corrects phases (fixes CP2K sign flips).
    """
    nstates = St_raw.shape[0]
    
    # --- 1. STATE TRACKING (Identity Reordering) ---
    # Find where the maximum overlap is for each column
    col_max_indices = np.argmax(np.abs(St_raw), axis=0)
    
    # Permute the matrix columns to restore the original state order
    St_tracked = St_raw[:, col_max_indices]
    
    # --- 2. PHASE CORRECTION ---
    St_corrected = St_tracked.copy()
    for i in range(nstates):
        # If the diagonal is negative, the state flipped phase
        if St_corrected[i, i] < 0:
            St_corrected[:, i] *= -1  # Flip column
            phase_cumulative[i] *= -1 # Track the cumulative phase
            
    # Calculate how many errors were fixed
    swaps_detected = np.sum(col_max_indices != np.arange(nstates))
    phases_detected = np.sum(np.diag(St_corrected) < 0) # Should be 0 after fixing
    
    return St_corrected, phase_cumulative, swaps_detected

# ---------------------------------------------------------
# Validation Loop
# ---------------------------------------------------------
start_step = 5000
end_step = 5005
dt_au = 41.34

print(f"{'Step':<8} | {'Corrupted RMSE':<16} | {'EasyWFO RMSE':<16} | {'Swaps Fixed'}")
print("-" * 65)

for step in range(start_step, end_step):
    st_file = f"{REF_DIR}/St_ci_{step}_re"
    if not os.path.exists(st_file):
        continue
        
    St_truth = np.loadtxt(st_file)
    true_nacs = (St_truth - St_truth.T) / (2 * dt_au)
    
    # -----------------------------------------------------
    # CORRUPT THE DATA (Simulate CP2K failures)
    # -----------------------------------------------------
    St_corrupted = St_truth.copy()
    
    # Failure 1: Swap State 2 and 3 (Columns only)
    St_corrupted[:, [1, 2]] = St_truth[:, [2, 1]] 
    
    # Failure 2: Flip the phase of State 5 (Column only)
    St_corrupted[:, 4] *= -1
    
    # Calculate the error of the corrupted data
    corrupted_nacs = (St_corrupted - St_corrupted.T) / (2 * dt_au)
    rmse_corrupted = np.sqrt(np.mean((corrupted_nacs - true_nacs)**2))
    
    # -----------------------------------------------------
    # FIX THE DATA (Run EasyWFO)
    # -----------------------------------------------------
    # Initialize cumulative phase tracker
    phase_cumulative = np.ones(St_truth.shape[0])
    
    St_fixed, phase_cumulative, swaps = easywfo_black_box_correction(St_corrupted, phase_cumulative)
    
    # Calculate the error of the fixed data
    fixed_nacs = (St_fixed - St_fixed.T) / (2 * dt_au)
    rmse_fixed = np.sqrt(np.mean((fixed_nacs - true_nacs)**2))
    
    print(f"{step:<8} | {rmse_corrupted:<16.4e} | {rmse_fixed:<16.4e} | {swaps}")

print("-" * 65)
print("If EasyWFO RMSE is 0.0000e+00, your Master Algorithm is flawless.")
