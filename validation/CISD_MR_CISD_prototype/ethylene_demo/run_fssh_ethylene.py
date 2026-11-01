import numpy as np
import matplotlib.pyplot as plt
import glob

files = sorted(glob.glob("st_adi_*.txt"), key=lambda x: int(x.split('_')[-1].split('.')[0]))
print(f"Running FSSH Ensemble (50 trajectories) on Ethylene...")

n_states = 3
dt = 0.5
n_trajs = 50
n_steps = len(files)

# Store populations for all trajectories: shape (n_trajs, n_steps, n_states)
all_pops = np.zeros((n_trajs, n_steps, n_states))

for traj in range(n_trajs):
    np.random.seed(traj) # Unique seed for each trajectory
    current_state = 1 # Start in State 2 (bright state)
    all_pops[traj, 0, current_state] = 1.0
    
    for i in range(1, n_steps):
        S = np.loadtxt(files[i])
        NAC = (S - S.T) / (2.0 * dt)
        
        hopped = False
        for j in range(n_states):
            if j != current_state:
                prob = max(0.0, -2.0 * NAC[current_state, j] * dt)
                if np.random.rand() < prob:
                    all_pops[traj, i, current_state] = 0.0
                    current_state = j
                    all_pops[traj, i, current_state] = 1.0
                    hopped = True
                    break
        if not hopped:
            all_pops[traj, i, current_state] = 1.0

# Average over trajectories
avg_pops = np.mean(all_pops, axis=0)
time_axis = np.arange(n_steps) * dt

# Plotting
fig, ax = plt.subplots(figsize=(10, 6))
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

# Plot faint individual trajectories for State 1 (ground state)
for traj in range(n_trajs):
    ax.plot(time_axis, all_pops[traj, :, 0], color='gray', alpha=0.1, linewidth=1)

# Plot ensemble averages
for state in range(n_states):
    ax.plot(time_axis, avg_pops[:, state], linewidth=4, color=colors[state], label=f"State {state+1} (Avg)")

ax.set_xlabel("Time (fs)", fontsize=16)
ax.set_ylabel("Population", fontsize=16)
ax.set_title("Ethylene MR-CISD: Ensemble FSSH Dynamics", fontsize=18)
ax.legend(fontsize=14, loc='center right')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("ethylene_ensemble_dynamics.png", dpi=300)
print("✅ Saved final publication-quality plot: ethylene_ensemble_dynamics.png")
