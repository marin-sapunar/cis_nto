import numpy as np
import time
import os
from liblibra_core import *
from libra_py import units
import libra_py.dynamics.tsh.compute as tsh_dynamics

def read_files(data_path, istep, fstep, dt, nstates):
    NSTEPS = fstep - istep
    St_adi, Ham_adi, Hvib_adi, NAC_adi, Basis_transform = [], [], [], [], []
    for i in range(istep, fstep):
        st = CMATRIX(nstates, nstates); st.Load_Matrix_From_File(f"{data_path}/st_adi_{i}.txt")
        St_adi.append(st)
        ham = CMATRIX(nstates, nstates); ham.Load_Matrix_From_File(f"{data_path}/ham_adi_{i}.txt")
        Ham_adi.append(ham)
        
        nac = (st - st.H()) / (2.0 * dt * units.fs2au)
        NAC_adi.append(nac)
        Hvib_adi.append(ham - nac * (0.0 + 1.0j))
        bas_trans = CMATRIX(nstates, nstates); bas_trans.identity(); Basis_transform.append(bas_trans)
    return NSTEPS, nstates, St_adi, Ham_adi, Hvib_adi, NAC_adi, Basis_transform

def main():
    ISTATE = 1; NSTATES = 3; DT_FS = 0.5; istep = 0; fstep = 49; NSTEPS = fstep - istep
    data_path = "results" 
    
    if not os.path.exists(f"{data_path}/st_adi_0.txt"):
        print(f"ERROR: Cannot find {data_path}/st_adi_0.txt. Make sure the 'results' folder is here.")
        return

    NSTEPS, NSTATES, St_adi, Ham_adi, Hvib_adi, NAC_adi, Basis_transform = read_files(data_path, istep, fstep, DT_FS, NSTATES)

    class abstr_class: pass
    def compute_model(q, params, full_id):
        timestep = params["timestep"]
        obj = abstr_class()
        obj.ham_adi = Ham_adi[timestep]
        obj.nac_adi = NAC_adi[timestep]
        obj.basis_transform = Basis_transform[timestep]
        obj.hvib_adi = Hvib_adi[timestep]
        obj.time_overlap_adi = St_adi[timestep]
        return obj

    model_params = { "timestep":0, "icond":0, "model0":0, "nstates":NSTATES }
    
    dyn_general = { 
        "nsteps":NSTEPS, "ntraj":250, "nstates":NSTATES, "dt":DT_FS*units.fs2au, "nfiles": NSTEPS,
        "decoherence_rates": MATRIX(NSTATES, NSTATES), "ave_gaps": MATRIX(NSTATES, NSTATES),
        "progress_frequency":0.1, "which_adi_states":range(NSTATES), "which_dia_states":range(NSTATES),
        "mem_output_level":2, "properties_to_save":[ "timestep", "time", "se_pop_adi", "sh_pop_adi" ],
        "isNBRA":0, "ham_update_method":2, "ham_transform_method":0, "time_overlap_method":0, 
        "hvib_update_method":0, "force_method":0, "rep_force":1, "hop_acceptance_algo":31, "momenta_rescaling_algo":0,
        "decoherence_algo":-1, "decoherence_times_type":-1, "dephasing_informed":0, "do_ssy":0,
        "tsh_method":0, "rep_tdse":1, "electronic_integrator": 5, "nac_update_method":0,
        "state_tracking_algo": 21, "MK_scaling_function": 1, "MK_alpha": 1100.0, "do_phase_correction": 1, 
    }

    for icond in range(20):
        print(f"Running ICOND = {icond}...", flush=True)
        model_params.update({"icond": icond})
        dyn_general.update({"prefix": f"FSSH_TUTORIAL_icond_{icond}"})
        
        rnd = Random() 
        elec_params = {"ndia":NSTATES, "nadi":NSTATES, "verbosity":-1, "init_dm_type":0, "init_type":1, "rep":1, "istate":ISTATE}
        nucl_params = {"ndof":1, "init_type":3, "q":[-10.0], "p":[0.0], "mass":[2000.0], "force_constant":[0.01], "verbosity":-1}
        
        tsh_dynamics.generic_recipe(dyn_general, compute_model, model_params, elec_params, nucl_params, rnd)
        
        # CRITICAL FIX: Sleep for >1 second to guarantee a new system clock seed
        time.sleep(1.1) 

if __name__ == "__main__":
    main()
