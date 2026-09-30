from pathlib import Path
import sys
import time
import resource
import numpy as np

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

import sim_moon_auto
import functions

masses_to_test = functions.build_mass_array(points_per_decade=100)

start_a = 0.001
end_a = 0.5
N_a = 200
sma_values = np.linspace(start_a, end_a, N_a)

# mittlere (mediane) Werte aus den jeweiligen Spektren
a_test = 0.3
m_test = 0.003

print(f"Teste a={a_test}, m={m_test}")

BASE_DIR = Path("/pfs/10/project/bw16g002/Lisa/Masterarbeit/critical_mass_search")
BASE_DIR.mkdir(exist_ok=True)

t0 = time.time()

sim = sim_moon_auto.setupSimulation(a=a_test, m=m_test)

try:
    ecc, sma, inc, omega, longitude, orbital_node, xyz_f, xyz_moon, break_year, end_year, early_stop = \
        sim_moon_auto.simulation(sim)

    sim_moon_auto.safe_data(
        ecc, sma, inc, omega, longitude, orbital_node, xyz_f, xyz_moon,
        a_test, m_test, break_year, end_year, early_stop, base_dir=BASE_DIR
    )

    print("break_year:", break_year)
    print("end_year:", end_year)
    print("early_stop:", early_stop)
    print(f"Save windows: {functions.compute_save_windows(break_year, end_year, early_stop, window_years=300)}")

except sim_moon_auto.SimulationInstabilityError as e:
    print(f"instabil ({e.reason}) bei t={e.year:.2f} Jahren")

t1 = time.time()
print(f"Laufzeit: {(t1 - t0)/60:.2f} Minuten")

# Spitzenspeicherverbrauch (Peak RSS) des Prozesses seit Start
peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print(f"Peak-Speicherverbrauch: {peak_kb / 1024:.1f} MB")

# Ergebnis:
""" 
The time is  4876 years 
The time is  4895 years 
The time is  4913 years 
The time is  4932 years 
The time is  4951 years 
The time is  4970 years 
The time is  4989 years 
break_year: {'psi1': np.float64(300.0154978397523), 'psi2': None, 'psi3': None}
end_year: 5000.0
early_stop: False
Save windows: [(np.float64(0.01549783975229957), np.float64(300.0154978397523)), (np.float64(4700.0), np.float64(5000.0))]
Laufzeit: 14.99 Minuten
Peak-Speicherverbrauch: 154.6 MB
"""