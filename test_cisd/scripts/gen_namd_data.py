import numpy as np
from pyscf import gto, scf, fci
from pyscf.fci import cistring

# Standard planar Benzene
geom_a = '''
C  1.3960  0.0000  0.0000
C  0.6980  1.2090  0.0000
C -0.6980  1.2090  0.0000
C -1.3960  0.0000  0.0000
C -0.6980 -1.2090  0.0000
C  0.6980 -1.2090  0.0000
H  2.4800  0.0000  0.0000
H  1.2400  2.1470  0.0000
H -1.2400  2.1470  0.0000
H -2.4800  0.0000  0.0000
H -1.2400 -2.1470  0.0000
H  1.2400 -2.1470  0.0000
'''
# Displace one carbon out of plane (simulating a B2u puckering vibration)
geom_b = geom_a.replace('C  1.3960  0.0000  0.0000', 'C  1.3960  0.0000  0.0500')

def run_cisd(geom, ncore, nact, n_act_elec):
    mol = gto.M(atom=geom, basis='sto-3g', verbose=0)
    mf = scf.RHF(mol).run(verbose=0)
    
    mol_act = mol.copy()
    mol_act.nelectron = n_act_elec
    
    cisolver = fci.FCI(mol_act, mf.mo_coeff[:, ncore:ncore+nact])
    e, civec = cisolver.kernel()
    return mf, civec.flatten()

# Active space: 6 electrons in 6 pi orbitals
ncore, nact, n_act_elec = 15, 6, 6
mf_a, civec_a = run_cisd(geom_a, ncore, nact, n_act_elec)
mf_b, civec_b = run_cisd(geom_b, ncore, nact, n_act_elec)

# Real MO Overlap between the two geometries
S_ao = mf_a.mol.intor_symmetric('int1e_ovlp')
S_mo = mf_a.mo_coeff.T @ S_ao @ mf_b.mo_coeff
for i in range(S_mo.shape[0]):
    if S_mo[i,i] < 0: S_mo[:,i] *= -1

# Map determinants (keep top 20 to simulate a truncated wavefunction)
idx = np.argsort(np.abs(civec_a))[::-1][:20]
coeffs_a = civec_a[idx]
coeffs_b = civec_b[idx] 

strs_a = cistring.gen_strings4orblist(range(nact), n_act_elec//2)
strs_b = cistring.gen_strings4orblist(range(nact), n_act_elec//2)
nb_str = len(strs_b)

n_mo = S_mo.shape[0]  # 36 for Benzene STO-3G
n_occ = ncore + 3     # 15 core + 3 active = 18 occupied
n_dets = 20

with open("real_data.txt", "w") as f:
    f.write(f"{n_mo} {n_occ} {n_dets}\n")
    for i in range(n_mo):
        f.write(" ".join(f"{x:.15f}" for x in S_mo[i]) + "\n")

    for k in range(n_dets):
        i_a = idx[k] // nb_str
        i_b = idx[k] % nb_str
        sa = strs_a[i_a]
        sb = strs_b[i_b]

        occ_a = [i+1 for i in range(ncore)] + [i+1+ncore for i in range(nact) if (sa & (1<<i))]
        occ_b = [i+1 for i in range(ncore)] + [i+1+ncore for i in range(nact) if (sb & (1<<i))]

        f.write(" ".join(map(str, occ_a)) + "\n")
        f.write(" ".join(map(str, occ_b)) + "\n")
        f.write(f"{coeffs_a[k]:.15f} {coeffs_b[k]:.15f}\n") # Write BOTH bra and ket coeffs

print(f"SUCCESS: True NAMD cross-overlap data for Benzene generated.")
print(f"n_mo={n_mo}, n_occ={n_occ}, n_dets={n_dets}")
