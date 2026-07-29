from pathlib import Path
import sys
import functions

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
	sys.path.insert(0, str(parent_dir))

# from sim_moon_auto import setupSimulation, simulation, safe_data

import sim_moon_auto
import functions

start_m = 0.00001
end_m = 0.01
start_a = 0.001
end_a = 0.5

i = start_a
while i <= end_a:
    j = start_m
    while j<= end_m:
        sim = sim_moon_auto.setupSimulation(a=i, m=j)
        try:
            ecc,sma,inc,omega,longitude,orbital_node,xyz_f,xyz_moon = sim_moon_auto.simulation(sim)
        except sim_moon_auto.SimulationInstabilityError as e:
            print(f"a={i}, m={j}: instabil ({e.reason}) bei t={e.year:.2f} Jahren")
            j += 0.00001
            j = round(j, 5)
            continue

        psi1, psi2, psi3 = functions.laplace_angles(longitude)
        is_resonant_psi1 = functions.is_laplace_resonant(psi1)
        is_resonant_psi2 = functions.is_laplace_resonant(psi2)
        is_resonant_psi3 = functions.is_laplace_resonant(psi3)
        
        if not is_resonant_psi1 and not is_resonant_psi2 and not is_resonant_psi3:
            print(f"Critical mass found for a={i} and m={j}")
        j += 0.00001
        j = round(j, 5)
    i += 0.001
    i = round(i, 3)