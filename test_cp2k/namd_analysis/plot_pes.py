import matplotlib
matplotlib.use('Agg') # Headless mode
import matplotlib.pyplot as plt
import pickle

times = [i for i in range(20)]
energies = {s: [] for s in range(1, 16)}

for i in range(1, 21):
    with open(f'states_cont_{i}.pkl', 'rb') as f:
        states = pickle.load(f)
    for s in range(1, 16):
        if s in states:
            energies[s].append(states[s]['energy'])
        else:
            energies[s].append(float('nan'))

plt.figure(figsize=(10, 6))
for s in range(1, 16):
    # Highlight S1 and S2
    if s == 1: plt.plot(times, energies[s], 'r-o', label='S1', linewidth=2)
    elif s == 2: plt.plot(times, energies[s], 'b-o', label='S2', linewidth=2)
    else: plt.plot(times, energies[s], 'gray', alpha=0.3)

plt.xlabel('Time (fs)')
plt.ylabel('Energy (eV)')
plt.title('Potential Energy Surfaces (300K Production)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('PES_plot.png', dpi=300)
print("Saved PES_plot.png. Download this and look at it!")
