from read_molden_mo import read_molden_mo
import sys
import numpy as np

def read_ao_overlap(path):
    with open(path) as f:
        n = int(f.readline().split()[0])
        return np.loadtxt(f, max_rows=n)


def lowdin_inv_sqrt(S):
    vals, vecs = np.linalg.eigh(S)
    return vecs @ np.diag(1.0 / np.sqrt(np.maximum(vals, 1e-12))) @ vecs.T

fa, fb = int(sys.argv[1]), int(sys.argv[2])
s_ao_path = f's_ao_{fa}_{fb}'

S_AO = read_ao_overlap(s_ao_path)
C_A  = read_molden_mo(f'p4-di-PBE_{fa}-MOS-1_0.molden')
C_B  = read_molden_mo(f'p4-di-PBE_{fb}-MOS-1_0.molden')

print(f'=== Frame {fa} -> {fb} ===')
print(f'S_AO: {S_AO.shape}   C_A: {C_A.shape}   C_B: {C_B.shape}')

S_AA = C_A.T @ S_AO @ C_A
S_BB = C_B.T @ S_AO @ C_B
X_AB = C_A.T @ S_AO @ C_B

S_AA_m12 = lowdin_inv_sqrt(S_AA)
S_BB_m12 = lowdin_inv_sqrt(S_BB)
X_lowdin = S_AA_m12 @ X_AB @ S_BB_m12

print(f'S_AA raw err:      {np.abs(S_AA - np.eye(S_AA.shape[0])).max():.4e}')
print(f'S_BB raw err:      {np.abs(S_BB - np.eye(S_BB.shape[0])).max():.4e}')
print(f'S_AA Löwdin resid: {np.abs(S_AA_m12 @ S_AA @ S_AA_m12 - np.eye(S_AA.shape[0])).max():.2e}')
print(f'S_BB Löwdin resid: {np.abs(S_BB_m12 @ S_BB @ S_BB_m12 - np.eye(S_BB.shape[0])).max():.2e}')

sig = np.linalg.svd(X_lowdin, compute_uv=False)
print(f'sigma range:      [{sig.min():.6f}, {sig.max():.8f}]')

if sig.max() > 1 + 1e-6:
    print(f'⚠  sigma_max > 1 by {sig.max() - 1:.2e}')
else:
    print(f'✅ sigma_max = {sig.max():.8f}  (physical bound)')

np.savetxt(f'overlap_{fa}_{fb}.norm', X_lowdin, fmt='%.10e')
with open(f's_wf_{fa}_{fb}', 'w') as f:
    f.write(f'{X_lowdin.shape[0]}\n')
    np.savetxt(f, X_lowdin, fmt='%.10e')

# Also update the "current" s_wf so transform_ci.py works without editing
with open('s_wf', 'w') as f:
    f.write(f'{X_lowdin.shape[0]}\n')
    np.savetxt(f, X_lowdin, fmt='%.10e')

print(f'Saved overlap_{fa}_{fb}.norm and s_wf_{fa}_{fb}')
