import numpy as np
from pyscf import gto, scf, fci
from pyscf.fci import cistring

# Benzene geometries
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

def run_cisd(geom, ncore, nact, n_act_elec):
    mol = gto.M(atom=geom, basis='sto-3g', verbose=0)
    mf = scf.RHF(mol).run(verbose=0)
    mol_act = mol.copy()
    mol_act.nelectron = n_act_elec
    cisolver = fci.FCI(mol_act, mf.mo_coeff[:, ncore:ncore+nact])
    e, civec = cisolver.kernel()
    return mf, civec.flatten()

ncore, nact, n_act_elec = 15, 6, 6
mf_a, civec_a = run_cisd(geom_a, ncore, nact, n_act_elec)
mf_b, civec_b = run_cisd(geom_b, ncore, nact, n_act_elec)

S_ao = mf_a.mol.intor_symmetric('int1e_ovlp')
S_mo = mf_a.mo_coeff.T @ S_ao @ mf_b.mo_coeff
for i in range(S_mo.shape[0]):
    if S_mo[i,i] < 0: S_mo[:,i] *= -1

strs_a = cistring.gen_strings4orblist(range(nact), n_act_elec//2)
strs_b = cistring.gen_strings4orblist(range(nact), n_act_elec//2)
nb_str = len(strs_b)

n_mo = S_mo.shape[0]
n_occ = ncore + 3

def write_file(filename, indices):
    with open(filename, "w") as f:
        f.write(f"{n_mo} {n_occ} {len(indices)} 1\n")
        for i in range(n_mo):
            f.write(" ".join(f"{x:.15f}" for x in S_mo[i]) + "\n")
            
        for k in indices:
            i_a = k // nb_str
            i_b = k % nb_str
            sa = strs_a[i_a]
            sb = strs_b[i_b]
            
            occ_a = [i+1 for i in range(ncore)] + [i+1+ncore for i in range(nact) if (sa & (1<<i))]
            occ_b = [i+1 for i in range(ncore)] + [i+1+ncore for i in range(nact) if (sb & (1<<i))]
            
            f.write(" ".join(map(str, occ_a)) + "\n")
            f.write(" ".join(map(str, occ_b)) + "\n")
            f.write(f"{civec_a[k]:.15f} {civec_b[k]:.15f}\n")

# Find the Hartree-Fock reference (all active orbitals filled up to n_act_elec/2)
ref_a = sum(1 << i for i in range(n_act_elec//2))
ref_b = ref_a
ref_idx_a = np.where(strs_a == ref_a)[0][0]
ref_idx_b = np.where(strs_b == ref_b)[0][0]
ref_idx = ref_idx_a * nb_str + ref_idx_b

singles_idx = []
doubles_idx = []

for k in range(len(civec_a)):
    if k == ref_idx: continue
    i_a = k // nb_str
    i_b = k % nb_str
    # Count excitations by comparing bits to reference
    exc_a = bin(strs_a[i_a] ^ ref_a).count('1') // 2
    exc_b = bin(strs_b[i_b] ^ ref_b).count('1') // 2
    total_exc = exc_a + exc_b
    
    if total_exc == 1:
        singles_idx.append(k)
    elif total_exc == 2:
        doubles_idx.append(k)

write_file("real_data_cisd.txt", np.argsort(np.abs(civec_a))[::-1][:36]) # Full
write_file("real_data_singles.txt", singles_idx) # Singles only
write_file("real_data_doubles.txt", doubles_idx) # Doubles only
print(f"Generated files: Full CISD (36 dets), Singles ({len(singles_idx)} dets), Doubles ({len(doubles_idx)} dets)")
