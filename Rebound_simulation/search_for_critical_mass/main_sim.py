from pathlib import Path
import sys

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
	sys.path.insert(0, str(parent_dir))

# from sim_moon_auto import setupSimulation, simulation, safe_data

import sim_moon_auto
import functions
import pandas as pd

start_m = 0.00001
end_m = 0.01
start_a = 0.001
end_a = 0.5

resonance_file = pd.read_csv("resonance_breaking.csv")
instability_file = pd.read_csv("instability_cases.csv")

try:
    i = start_a
    while i <= end_a:
        
        # Merkt sich für dieses sma, bei welcher Masse die jeweilige
        # Resonanz zum ersten Mal gebrochen ist (None = noch nie gebrochen)
        break_mass_psi1 = None
        break_mass_psi2 = None
        break_mass_psi3 = None
        
        j = start_m
        while j<= end_m:
            sim = sim_moon_auto.setupSimulation(a=i, m=j)
            try:
                ecc,sma,inc,omega,longitude,orbital_node,xyz_f,xyz_moon = sim_moon_auto.simulation(sim)
            except sim_moon_auto.SimulationInstabilityError as e:
                functions.write_instability_row(instability_file, e.reason, e.year, i, j)
                print(f"a={i}, m={j}: instabil ({e.reason}) bei t={e.year:.2f} Jahren")
                j += 0.00001
                j = round(j, 5)
                continue

            psi1, psi2, psi3 = functions.laplace_angles(longitude)
            is_resonant_psi1 = functions.is_laplace_resonant(psi1)
            is_resonant_psi2 = functions.is_laplace_resonant(psi2)
            is_resonant_psi3 = functions.is_laplace_resonant(psi3)
            
            # nur beim ERSTEN Brechen die Masse festhalten
            if not is_resonant_psi1 and break_mass_psi1 is None:
                break_mass_psi1 = j
            if not is_resonant_psi2 and break_mass_psi2 is None:
                break_mass_psi2 = j
            if not is_resonant_psi3 and break_mass_psi3 is None:
                break_mass_psi3 = j
                
            j += 0.00001
            j = round(j, 5)
            
        functions.write_resonance_row(resonance_file, i, break_mass_psi1, break_mass_psi2, break_mass_psi3)
            
        i += 0.001
        i = round(i, 3)
        
finally:
    resonance_file.close()
    instability_file.close()