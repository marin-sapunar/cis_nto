import numpy as np
from pyscf import gto, scf, fci
from pyscf.fci import cistring

geom_a = """
H 0.0000000 0.0000000 0.0000000
N 0.0000000 0.0000000 1.0153460
C 1.2773601 0.0000000 1.7736663
C 1.2025642 -1.0803159 2.8827115
O 0.2174042 -1.1206722 3.6605502
H 1.5510343 0.9533542 2.2389295
N 2.2819439 -1.9130790 2.9853964
C 2.5527241 -2.7117443 4.1981942
C 3.8505248 -2.1878740 4.8692728
O 4.8487506 -1.9434037 4.1410359
N 3.8549146 -2.0297871 6.2222362
C 5.0559800 -1.4509892 6.8881543
C 5.2179200 -2.0915376 8.3009464
O 4.2949240 -2.0837731 9.1513927
H 5.0874371 -0.3508535 6.8495547
N 6.4347018 -2.6844644 8.5069308
H 6.6494282 -3.3332756 9.2744926
H 4.9525659 -0.5253686 -0.9845894
C 4.9256175 -1.2143869 -0.1096474
O 3.9266287 -1.2479954 0.6482254
N 6.0590835 -1.9509624 0.0999870
C 6.2919689 -2.7150761 1.3391207
C 7.6433932 -2.3165235 1.9862440
O 8.6614976 -2.1998167 1.2540089
N 7.6640198 -2.1538414 3.3392114
C 8.9253593 -1.8103360 4.0310941
C 9.0315228 -2.5605151 5.3812477
O 8.0108046 -2.9263712 6.0194565
H 9.1264944 -0.7274925 4.1194675
N 10.3162314 -2.8006441 5.8014575
C 10.6339476 -3.7870357 6.8616327
C 11.2097607 -3.1741987 8.1683830
O 12.3051961 -2.5970260 8.2127626
H 11.3766978 -4.5028005 6.4846698
O 10.4552836 -3.4531787 9.2345369
H -0.8799378 -0.1687248 1.4965178
H 2.0618704 -0.2519805 1.0442845
H 3.0562776 -1.7809528 2.3070137
H 1.6827465 -2.6557452 4.8733535
H 2.7255326 -3.7747180 3.9219594
H 2.9627443 -2.1196018 6.7429229
H 5.9109221 -1.7886832 6.2767265
H 7.0758518 -2.7133131 7.6979192
H 5.4458708 -2.5526737 2.0280533
H 6.7732631 -2.1791025 3.8686397
H 9.7352085 -2.1936400 3.3879564
H 11.0796122 -2.5424444 5.1505140
H 9.6946800 -4.3049793 7.1074725
H 6.8540612 -1.8218538 -0.5527054
H 6.3600408 -3.8014957 1.1082858
H 10.9398738 -3.1325161 10.1154717
"""

lines = geom_a.strip().split('\n')
lines[4] = "O 0.2174042 -1.1206722 3.7105502"
geom_b = '\n'.join(lines)

def run_fci(geom):
    mol = gto.M(atom=geom, basis='sto-3g', verbose=0)
    mf = scf.RHF(mol).run(verbose=0)
    
    nocc_mo = np.count_nonzero(mf.mo_occ > 0)
    ncore = nocc_mo - 2
    nact = 4 # 2 occupied + 2 virtual
    
    # Fix: copy mol and set nelectron manually
    mol_act = mol.copy()
    mol_act.nelectron = 4
    
    cisolver = fci.FCI(mol_act, mf.mo_coeff[:, ncore:ncore+nact])
    e, civec = cisolver.kernel()
    return mf, civec.flatten(), ncore, nocc_mo

print("Running RHF and FCI on the 77-atom 2B7 amyloid geometry...")
mf_a, civec_a, ncore, nocc_mo = run_fci(geom_a)
mf_b, civec_b, _, _ = run_fci(geom_b)

S_ao = mf_a.mol.intor_symmetric('int1e_ovlp')
S_mo = mf_a.mo_coeff.T @ S_ao @ mf_b.mo_coeff
for i in range(S_mo.shape[0]):
    if S_mo[i,i] < 0: S_mo[:,i] *= -1

nact_occ = 2
nact = 4

strs_a = cistring.gen_strings4orblist(range(nact), nact_occ)
strs_b = cistring.gen_strings4orblist(range(nact), nact_occ)
nb_str = len(strs_b)

n_mo = S_mo.shape[0]
n_occ = ncore + nact_occ
n_dets = len(civec_a)

with open("real_data_amyloid.txt", "w") as f:
    f.write(f"{n_mo} {n_occ} {n_dets} 1\n")
    for i in range(n_mo):
        f.write(" ".join(f"{x:.15f}" for x in S_mo[i]) + "\n")
        
    for k in range(n_dets):
        i_a = k // nb_str
        i_b = k % nb_str
        sa = strs_a[i_a]
        sb = strs_b[i_b]
        
        occ_a = [i+1 for i in range(ncore)] + [i+1+ncore for i in range(nact) if (sa & (1<<i))]
        occ_b = [i+1 for i in range(ncore)] + [i+1+ncore for i in range(nact) if (sb & (1<<i))]
        
        f.write(" ".join(map(str, occ_a)) + "\n")
        f.write(" ".join(map(str, occ_b)) + "\n")
        f.write(f"{civec_a[k]:.15f} {civec_b[k]:.15f}\n")

print(f"SUCCESS: Generated real_data_amyloid.txt (n_mo={n_mo}, n_occ={n_occ}, n_dets={n_dets})")
