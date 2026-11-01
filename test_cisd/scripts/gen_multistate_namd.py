import numpy as np
from pyscf import gto, scf, fci
from pyscf.fci import cistring

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
geom_b = geom_a.replace('C  1.3960  0.0000  0.0000', 'C  1.3960  0.0000  0.0500')

def run_cisd(geom, ncore, nact, n_act_elec, nstates):
    mol = gto.M(atom=geom, basis='sto-3g', verbose=0)
    mf = scf.RHF(mol).run(verbose=0)
    mol_act = mol.copy()
    mol_act.nelectron = n_act_elec
    cisolver = fci.FCI(mol_act, mf.mo_coeff[:, ncore:ncore+nact])
    # Get the first nstates roots
    e, civec = cisolver.kernel(nroots=nstates)
    return mf, civec

ncore, nact, n_act_elec, nstates = 15, 6, 6, 3
mf_a, civecs_a = run_cisd(geom_a, ncore, nact, n_act_elec, nstates)
mf_b, civecs_b = run_cisd(geom_b, ncore, nact, n_act_elec, nstates)

# Real MO Overlap
S_ao = mf_a.mol.intor_symmetric('int1e_ovlp')
S_mo = mf_a.mo_coeff.T @ S_ao @ mf_b.mo_coeff
for i in range(S_mo.shape[0]):
    if S_mo[i,i] < 0: S_mo[:,i] *= -1

strs_a = cistring.gen_strings4orblist(range(nact), n_act_elec//2)
strs_b = cistring.gen_strings4orblist(range(nact), n_act_elec//2)
nb_str = len(strs_b)

n_mo = S_mo.shape[0]
n_occ = ncore + 3
n_dets = 20 # Keep top 20 dets per state to test truncation

with open("real_data_multi.txt", "w") as f:
    f.write(f"{n_mo} {n_occ} {n_dets} {nstates}\n")
    for i in range(n_mo):
        f.write(" ".join(f"{x:.15f}" for x in S_mo[i]) + "\n")
        
    for s in range(nstates):
        civec_a = civecs_a[s].flatten()
        civec_b = civecs_b[s].flatten()
        idx = np.argsort(np.abs(civec_a))[::-1][:n_dets]
        coeffs_a = civec_a[idx]
        coeffs_b = civec_b[idx]

        for k in range(n_dets):
            i_a = idx[k] // nb_str
            i_b = idx[k] % nb_str
            sa = strs_a[i_a]
            sb = strs_b[i_b]
            
            occ_a = [i+1 for i in range(ncore)] + [i+1+ncore for i in range(nact) if (sa & (1<<i))]
            occ_b = [i+1 for i in range(ncore)] + [i+1+ncore for i in range(nact) if (sb & (1<<i))]
            
            f.write(" ".join(map(str, occ_a)) + "\n")
            f.write(" ".join(map(str, occ_b)) + "\n")
            f.write(f"{coeffs_a[k]:.15f} {coeffs_b[k]:.15f}\n")

print("Generated real_data_multi.txt (3 states)")
