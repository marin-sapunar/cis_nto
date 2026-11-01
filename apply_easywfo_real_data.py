import pickle
import numpy as np
from scipy.optimize import linear_sum_assignment

print("=== Applying EasyWFO to Raw CP2K Data ===")

# 1. Load your raw data (Adjust paths/names as needed based on your pkl structure)
# For demonstration, we assume you load two adjacent frames
# with open(f"{RAW_DATA_DIR}/states_cont_2.pkl", "rb") as f:
#     states_t = pickle.load(f)
# with open(f"{RAW_DATA_DIR}/states_cont_3.pkl", "rb") as f:
#     states_t_plus_dt = pickle.load(f)

# 2. Compute raw overlap (Example: if you have the S_IJ matrix pre-calculated)
# S_raw = np.load(f"{RAW_DATA_DIR}/S_IJ_2_3_rawphased.npy")

print("TODO: Insert your raw .pkl parsing logic here to build S_raw")
print("Then, run the Hungarian + Phase Correction engine on S_raw.")
