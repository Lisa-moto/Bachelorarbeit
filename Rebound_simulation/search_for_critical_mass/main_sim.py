from pathlib import Path
import sys

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
	sys.path.insert(0, str(parent_dir))

# from sim_moon_auto import setupSimulation, simulation, safe_data

import sim_moon_auto
import functions

#start_m = 0.00001
#end_m = 0.01
masses_to_test = functions.build_mass_array(points_per_decade=50)
start_a = 0.001
end_a = 0.5

resonance_file = functions.init_resonance_file("resonance_breaking.csv")
instability_file = functions.init_instability_file("instability_cases.csv")

try:
    i = start_a
    while i <= end_a:

        # Für jeden Winkel: erste Masse, bei der er bricht, und das Jahr,
        # in dem er bei dieser Masse gebrochen war
        break_mass_psi1, break_year_psi1 = None, None
        break_mass_psi2, break_year_psi2 = None, None
        break_mass_psi3, break_year_psi3 = None, None

        
        for j in masses_to_test:
            sim = sim_moon_auto.setupSimulation(a=i, m=j)

            try:
                ecc, sma, inc, omega, longitude, orbital_node, xyz_f, xyz_moon, break_year = \
                    sim_moon_auto.simulation(sim)
            except sim_moon_auto.SimulationInstabilityError as e:
                functions.write_instability_row(instability_file, e.reason, e.year, i, j)
                print(f"a={i}, m={j}: instabil ({e.reason}) bei t={e.year:.2f} Jahren")
                continue

            # nur beim ERSTEN Brechen Masse UND Jahr festhalten
            if break_year["psi1"] is not None and break_mass_psi1 is None:
                break_mass_psi1 = j
                break_year_psi1 = break_year["psi1"]
            if break_year["psi2"] is not None and break_mass_psi2 is None:
                break_mass_psi2 = j
                break_year_psi2 = break_year["psi2"]
            if break_year["psi3"] is not None and break_mass_psi3 is None:
                break_mass_psi3 = j
                break_year_psi3 = break_year["psi3"]


        functions.write_resonance_row(
            resonance_file, i,
            break_mass_psi1, break_year_psi1,
            break_mass_psi2, break_year_psi2,
            break_mass_psi3, break_year_psi3,
        )

        i += 0.001
        i = round(i, 3)

finally:
    resonance_file.close()
    instability_file.close()