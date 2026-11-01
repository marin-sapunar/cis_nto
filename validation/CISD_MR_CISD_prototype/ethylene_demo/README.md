# Ethylene Conical Intersection Demo (MR-CISD FSSH)

This demo reproduces the ensemble FSSH population dynamics through a conical intersection, a standard benchmark for nonadiabatic dynamics.

## Physics
We simulate 50 trajectories starting in the bright S2 state of ethylene. As the C=C double bond twists, the molecule passes through a conical intersection where S2 couples to the ground S1 state. The ensemble average shows the ultrafast decay from S2 to S1.

## Reproduce the Figure

    python3 generate_trajectory.py
    python3 run_fssh_ethylene.py

This generates 50 `st_adi_*.txt` matrices and the final `ethylene_ensemble_dynamics.png` plot.

## Expected Output
- 50 overlap matrices describing the CI passage.
- Ensemble FSSH dynamics showing population transfer from S2 (bright) to S1 (ground).
- The gray lines in the plot show individual trajectory hops; the bold lines are the ensemble average.
