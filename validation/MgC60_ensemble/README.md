# Mg@C60 NAMD Validation (EasyWFO -> Libra FSSH)

This directory contains the validation scripts for running the Mg@C60 ensemble dynamics using the EasyWFO overlaps fed into Libra. 

## Reproducibility Steps

1. **Fetch the Data:** If you want to run the full 50-step dynamics, download `data3.tar.bz2` from the official Libra tutorial repo and extract it so `data3/res` is accessible. For this validation, a sample of the raw `.txt` files is already provided in the `results/` folder.

2. **Run the Ensemble:**
   python3 run_ensemble.py
   *Note: This script includes a `time.sleep(1.1)` delay between runs to ensure Libra's system-clock-based RNG generates unique seeds for all 20 initial conditions, preventing the "flat line" deterministic hopping artifact.*

3. **Generate the Lifetime Plot:**
   python3 plot_lifetime.py

4. **Generate the Comparison Plot:**
   python3 plot_comparison.py

## Expected Output
* `MgC60_ensemble_lifetime_fit.png` (Fit yields Tau ≈ 990.4 fs, Beta ≈ 0.16, R² ≈ 0.79).
* `all_populations_comparison.png` (Shows EasyWFO LD vs Libra NACs).
