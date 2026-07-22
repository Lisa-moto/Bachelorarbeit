from pathlib import Path
import sys

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
        ecc,sma,inc,omega,longitude,orbital_node,xyz_f,xyz_moon = sim_moon_auto.simulation(sim)
        #sim_moon_auto.safe_data(ecc, sma, inc, omega, longitude, orbital_node, xyz_f, xyz_moon, a=i, m=j)
        #phi0, phi1, phi2, phi3, phi4, psi1, psi2, psi3 = resonance_plot_auto.resonant_angles(longitude, omega, orbital_node)
        #resonance_plot_auto.make_plot(phi0, phi1, phi2, phi3, phi4, psi1, psi2, psi3, a=i, m=j)
        is_resonant = functions.check_resonance(longitude, omega, orbital_node)
        if not is_resonant:
            print(f"Critical mass found for a={i} and m={j}")
            j += 0.00001
            j = round(j, 5)
    i += 0.001
    i = round(i, 3)