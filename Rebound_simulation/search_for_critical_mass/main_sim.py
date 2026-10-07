from pathlib import Path
import sys
import numpy as np

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
	sys.path.insert(0, str(parent_dir))

# from sim_moon_auto import setupSimulation, simulation, safe_data

import sim_moon_auto
import functions

masses_to_test = functions.build_mass_array(points_per_decade=100)

start_a = 0.001
end_a = 0.5
N_a = 200
sma_values = np.linspace(start_a, end_a, N_a)

# Index des zu berechnenden a-Werts wird von aussen übergeben
# (z.B. vom SLURM Job Array via SLURM_ARRAY_TASK_ID)
a_index = int(sys.argv[1])
i = round(float(sma_values[a_index]), 6)

# "Masterarbeit" muss manuell existieren.
BASE_DIR = Path("/pfs/10/project/bw16g002/Lisa/Masterarbeit/critical_mass_search")
BASE_DIR.mkdir(exist_ok=True)  # legt nur "critical_mass_search" an, "Masterarbeit" muss schon da sein

readme_path = BASE_DIR / "README_columns.txt"
if not readme_path.exists():
    functions.atomic_write_text(
        readme_path,
        "Pro Simulation eine Datei:\n"
        "  simulation_data/a={a:.6f}/a={a:.6f}_m={m:.6e}.npz  (Datentyp: float32)\n"
        "  Beispiel: simulation_data/a=0.100000/a=0.100000_m=1.000000e-05.npz\n\n"
        "Laden mit:\n"
        "  d = np.load('.../data.npz')\n"
        "  ecc = d['ecc']\n\n"
        "Schluessel in data.npz:\n"
        "  ecc, sma (in AU), inc (in Grad), orbital_node (rad), omega (rad), l (rad):\n"
        "    Spalte 0: year\n"
        "    Spalte 1-7: b, c, d, e, f, g, moon\n"
        "  xyz_f, xyz_moon (in Meter):\n"
        "    Spalte 0: year\n"
        "    Spalte 1-3: x, y, z\n"
    )
    
    
# Eigene Unterordner für die Resonanz- und Instabilitaets-Dateien
resonance_dir = BASE_DIR / "resonance_breaking"
instability_dir = BASE_DIR / "instability_cases"
log_dir = BASE_DIR / "progress_logs"
for d in (resonance_dir, instability_dir, log_dir):
    d.mkdir(exist_ok=True)

resonance_path = resonance_dir / f"resonance_breaking_a={a_index:03d}.csv"
instability_path = instability_dir / f"instability_cases_a={a_index:03d}.csv"
log_path = log_dir / f"progress_a={a_index:03d}.log"

entries = functions.load_log(log_path, masses_to_test)
if entries:
    print(f"a_index={a_index} (a={i}): {len(entries)} Massen bereits erledigt, mache weiter.")

for idx in range(len(entries), len(masses_to_test)):
    if all(b is not None for b in functions.first_breaks(entries)):
        break  # alle drei Resonanzen gebrochen, fertig

    j = masses_to_test[idx]
    sim = sim_moon_auto.setupSimulation(a=i, m=j)

    try:
        ecc, sma, inc, omega, longitude, orbital_node, xyz_f, xyz_moon, break_year, end_year, early_stop = \
            sim_moon_auto.simulation(sim)
    except sim_moon_auto.SimulationInstabilityError as e:
        fields = [str(idx), repr(float(j)), "unstable",
                  str(e.reason).replace(";", ":").replace("\n", " "), f"{float(e.year):.6f}", "-"]
        functions.append_log(log_path, fields)
        entries.append(fields)
        print(f"a={i}, m={j}: instabil ({e.reason}) bei t={e.year:.2f} Jahren")
        continue

    sim_moon_auto.safe_data(
        ecc, sma, inc, omega, longitude, orbital_node, xyz_f, xyz_moon,
        i, j, break_year, end_year, early_stop, base_dir=BASE_DIR
    )

    fields = [str(idx), repr(float(j)), "ok"] + [
        "None" if break_year[k] is None else f"{float(break_year[k]):.6f}"
        for k in ("psi1", "psi2", "psi3")
    ]
    functions.append_log(log_path, fields)
    entries.append(fields)

functions.write_outputs_from_log(entries, i, resonance_path, instability_path)