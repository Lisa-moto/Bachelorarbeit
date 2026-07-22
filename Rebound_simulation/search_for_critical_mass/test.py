import functions
import os
import numpy as np
from pathlib import Path
import sys

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
	sys.path.insert(0, str(parent_dir))

a_moon_short = 0.3
m_moon_short = 0.004

base_dir = Path(__file__).resolve().parents[1]
f_longitude = base_dir / 'data_aufgabe1_auto' / f'data_moon_a={a_moon_short}' / f'l_moon_a={a_moon_short}_m={m_moon_short}.txt'

longitude = np.loadtxt(f_longitude) # mean longitude (in rad)

psi1, psi2, psi3 = functions.laplace_angles(longitude)

is_resonant, diagnostics = functions.is_laplace_resonant(psi2, threshold_deg=177.0, return_diagnostics=True)

print(f"Is the system in Laplace resonance? {is_resonant}")
print(f"Amplitude of psi3: {diagnostics['amplitude_deg']:.2f} degrees")
print(f"Mean angle of psi3: {diagnostics['mean_angle_deg']:.2f} degrees")
print(f"Circular standard deviation of psi3: {diagnostics['circular_std_deg']:.2f} degrees")