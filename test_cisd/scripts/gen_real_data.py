import numpy as np
from pyscf import gto, scf, fci
from pyscf.fci import cistring

geom_a = '''
C  0.0000  0.0000  0.6670
C  0.0000  0.0000 -0.6670
H  0.9230  0.0000  1.2380
H -0.9230  0.0000  1.2380
H  0.9230  0.0000 -1.2380
H -0.9230  0.0000 -1.2380
'''
geom_b = geom_a.replace('C  0.0000  0.0000  0.6670', 'C  0.0000  0.0000  0.6750')

mol_a = gto.M(atom=geom_a, basis='sto-3g', verbose=0)
mol_b = gto.M(atom=geom_b, basis='sto-3g', verbose=0)
mf_a = scf.RHF(mol_a).run(verbose=0)
mf_b = scf.RHF(mol_b).run(verbose=0)

# 1. Real MO Overlap
S_ao = mol_a.intor_symmetric('int1e_ovlp')
S_mo = mf_a.mo_coeff.T @ S_ao @ mf_b.mo_coeff
for i in range(S_mo.shape[0]):
    if S_mo[i,i] < 0: S_mo[:,i] *= -1

# 2. Real Wavefunction (FCI in 4-orbital active space)
ncore = 6
nact = 4

# The trick: copy the molecule and manually set the number of electrons
mol_active = mol_a.copy()
mol_active.nelectron = 4  # 4 active electrons

cisolver = fci.FCI(mol_active, mf_a.mo_coeff[:, ncore:ncore+nact])
e, civec = cisolver.kernel()
civec = civec.flatten()

# Get top 10 determinants
idx = np.argsort(np.abs(civec))[::-1][:10]
coeffs = civec[idx]

strs_a = cistring.gen_strings4orblist(range(nact), 2)
strs_b = cistring.gen_strings4orblist(range(nact), 2)
nb_str = len(strs_b)

n_mo = S_mo.shape[0]  # 14 for ethylene STO-3G
n_occ = ncore + 2     # 6 core + 2 active = 8 occupied
n_dets = 10

with open("real_data.txt", "w") as f:
    f.write(f"{n_mo} {n_occ} {n_dets}\n")
    for i in range(n_mo):
        f.write(" ".join(f"{x:.15f}" for x in S_mo[i]) + "\n")

    for k in range(n_dets):
        i_a = idx[k] // nb_str
        i_b = idx[k] % nb_str
        sa = strs_a[i_a]
        sb = strs_b[i_b]

        # Core orbitals: 1 to ncore
        occ_a = [i+1 for i in range(ncore)]
        occ_b = [i+1 for i in range(ncore)]

        # Active orbitals: ncore+1 to ncore+nact
        for i in range(nact):
            if (sa & (1<<i)): occ_a.append(ncore + i + 1)
            if (sb & (1<<i)): occ_b.append(ncore + i + 1)

        f.write(" ".join(map(str, occ_a)) + "\n")
        f.write(" ".join(map(str, occ_b)) + "\n")
        f.write(f"{coeffs[k]:.15f}\n")

print(f"SUCCESS: Generated real_data.txt with n_mo={n_mo}, n_occ={n_occ}, n_dets={n_dets}")
