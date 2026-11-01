import subprocess
import re

# Sweep over different numbers of determinants
n_dets_list = [5, 10, 20, 36, 50, 100]
overlaps = []

print("N_dets | Overlap S")
print("-" * 25)
for n in n_dets_list:
    # Generate file for this N
    # (Assuming you already have gen_truncation.py from earlier)
    # We will just run the Fortran code on the files we generated
    filename = f"real_data_n{n}.txt"
    try:
        result = subprocess.run(["./check_trunc", filename], capture_output=True, text=True)
        output = result.stdout
        match = re.search(r'S=\s+([-0-9.E+]+)', output)
        if match:
            s_val = float(match.group(1))
            overlaps.append((n, s_val))
            print(f"{n:6d} | {s_val:.9f}")
    except:
        pass
