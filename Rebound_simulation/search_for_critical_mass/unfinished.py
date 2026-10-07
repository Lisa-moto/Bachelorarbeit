import os, sys
#sys.path.insert(0, "/Pfad/zu/functions.py")   # Ordner anpassen
import functions

N_a = 200
n_masses = len(functions.build_mass_array(points_per_decade=100))
base = "/pfs/10/project/bw16g002/Lisa/Masterarbeit/critical_mass_search/progress_logs"

todo = []
for a in range(N_a):
    p = f"{base}/progress_a={a:03d}.log"
    entries = []
    if os.path.exists(p):
        for line in open(p).read().split("\n")[:-1]:
            parts = line.split(";")
            if len(parts) == 6:
                entries.append(parts)
    done = len(entries) >= n_masses or all(b is not None for b in functions.first_breaks(entries))
    if not done:
        todo.append(str(a))
print(",".join(todo))