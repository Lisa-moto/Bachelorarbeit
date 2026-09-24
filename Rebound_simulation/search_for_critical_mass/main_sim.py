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


# "Masterarbeit" muss manuell existieren.
BASE_DIR = Path("Masterarbeit") / "critical_mass_search"
BASE_DIR.mkdir(exist_ok=True)  # legt nur "critical_mass_search" an, "Masterarbeit" muss schon da sein

readme_path = BASE_DIR / "README_columns.txt"
if not readme_path.exists():
    readme_path.write_text(
        "Spaltenreihenfolge in den .npy-Dateien (Datentyp: float32):\n\n"
        "ecc.npy, sma.npy, inc.npy, orbital_node.npy, omega.npy, l.npy:\n"
        "  Spalte 0: year\n"
        "  Spalte 1-7: b, c, d, e, f, g, moon\n\n"
        "xyz_f.npy, xyz_moon.npy:\n"
        "  Spalte 0: year\n"
        "  Spalte 1-3: x, y, z\n"
    )

resonance_file = functions.init_resonance_file(BASE_DIR / "resonance_breaking.csv")
instability_file = functions.init_instability_file(BASE_DIR / "instability_cases.csv")

try:
    for i in sma_values:
        i = round(float(i), 6)

        break_mass_psi1, break_year_psi1 = None, None
        break_mass_psi2, break_year_psi2 = None, None
        break_mass_psi3, break_year_psi3 = None, None

        for j in masses_to_test:
            sim = sim_moon_auto.setupSimulation(a=i, m=j)

            try:
                ecc, sma, inc, omega, longitude, orbital_node, xyz_f, xyz_moon, break_year, end_year, early_stop = \
                    sim_moon_auto.simulation(sim)
            except sim_moon_auto.SimulationInstabilityError as e:
                functions.write_instability_row(instability_file, e.reason, e.year, i, j)
                print(f"a={i}, m={j}: instabil ({e.reason}) bei t={e.year:.2f} Jahren")
                continue

            sim_moon_auto.safe_data(
                ecc, sma, inc, omega, longitude, orbital_node, xyz_f, xyz_moon,
                i, j, break_year, end_year, early_stop, base_dir=BASE_DIR
            )

            if break_year["psi1"] is not None and break_mass_psi1 is None:
                break_mass_psi1 = j
                break_year_psi1 = break_year["psi1"]
            if break_year["psi2"] is not None and break_mass_psi2 is None:
                break_mass_psi2 = j
                break_year_psi2 = break_year["psi2"]
            if break_year["psi3"] is not None and break_mass_psi3 is None:
                break_mass_psi3 = j
                break_year_psi3 = break_year["psi3"]
            if break_mass_psi1 is not None and break_mass_psi2 is not None and break_mass_psi3 is not None:
                break

        functions.write_resonance_row(
            resonance_file, i,
            break_mass_psi1, break_year_psi1,
            break_mass_psi2, break_year_psi2,
            break_mass_psi3, break_year_psi3,
        )

finally:
    resonance_file.close()
    instability_file.close()