import matplotlib.pyplot as plt
import numpy as np

# Set a professional style
plt.rcParams.update({
    'font.size': 12, 
    'font.family': 'serif',
    'mathtext.fontset': 'dejavuserif' # Ensures math matches the serif font
})

# ==========================================
# Figure 1: Truncation Convergence
# ==========================================
n_dets = [5, 10, 20, 36]
overlaps = [0.9876055, 0.9876376, 0.9876979, 0.9876993]
exact_S = 0.9876993

fig1, ax1 = plt.subplots(figsize=(7, 4.5))
ax1.plot(n_dets, overlaps, 'o-', color='blue', linewidth=2, markersize=8, label='L4M-CISD Truncated')
ax1.axhline(y=exact_S, color='red', linestyle='--', linewidth=2, label='Full CISD (Exact)')

# Use raw strings (r'...') for LaTeX to prevent underscore escaping issues
ax1.set_xlabel(r'Number of Determinants ($N_{\mathrm{dets}}$)')
ax1.set_ylabel(r'Overlap $S(t_1, t_2)$')
ax1.set_title('Overlap Convergence Under Wavefunction Truncation', pad=15)
ax1.legend(loc='lower right')
ax1.grid(True, linestyle=':', alpha=0.7)

# Set y-axis limits to zoom in
ax1.set_ylim(0.98758, 0.98772)

# Ensure nothing gets cut off
fig1.tight_layout()
fig1.savefig('figure1_truncation.png', dpi=300, bbox_inches='tight')
print("Saved figure1_truncation.png")


# ==========================================
# Figure 2: Speedup Scaling
# ==========================================
configurations = ['nstates=2\n(n_occ=100)', 'nstates=2\n(n_occ=121)', 'nstates=16\n(n_occ=121)']
brute_times = [1800, 3510, 72021]
l4m_times = [22.3, 25.9, 7808]

x = np.arange(len(configurations))
width = 0.35

fig2, ax2 = plt.subplots(figsize=(8, 5))
rects1 = ax2.bar(x - width/2, brute_times, width, label='Brute-Force CISD', color='gray')
rects2 = ax2.bar(x + width/2, l4m_times, width, label='L4M-CISD', color='green')

# Add speedup annotations above the L4M bars
speedups = [81, 135, 9.2]
for i, rect in enumerate(rects2):
    height = rect.get_height()
    ax2.annotate(f'{speedups[i]}x',
                 xy=(rect.get_x() + rect.get_width() / 2, height),
                 xytext=(0, 3),
                 textcoords="offset points",
                 ha='center', va='bottom', fontweight='bold', color='darkgreen', fontsize=14)

ax2.set_ylabel('Time per Evaluation (ms)')
ax2.set_title('Performance: Brute-Force vs. L4M-CISD', pad=15)
ax2.set_xticks(x)
ax2.set_xticklabels(configurations)
ax2.set_yscale('log')
ax2.legend(loc='upper left')
ax2.grid(True, linestyle=':', alpha=0.7, axis='y')

fig2.tight_layout()
fig2.savefig('figure2_speedup.png', dpi=300, bbox_inches='tight')
print("Saved figure2_speedup.png")
