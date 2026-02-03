import numpy as np
import matplotlib.pyplot as plt
import os

m_moon_short = 0.001  # mass of the moon relative to planet f

##########################################################################################
# load data
##########################################################################################

here = os.path.dirname(__file__)

# files without moon
f_sma = os.path.join(here, 'data_sim/sma_planets.txt')
f_longitude = os.path.join(here, 'data_sim/l_planets.txt')
f_omega = os.path.join(here, 'data_sim/omega_planets.txt')
f_node = os.path.join(here, 'data_sim/orbital_node_planets.txt')
f_ecc = os.path.join(here, 'data_sim/ecc_planets.txt')
f_inc = os.path.join(here, 'data_sim/inc_planets.txt')

"""
# filenames
f_sma = os.path.join(here, f'data_with_moon/sma_with_moon_m={m_moon_short}.txt')
f_longitude = os.path.join(here, f'data_with_moon/l_with_moon_m={m_moon_short}.txt')
f_omega = os.path.join(here, f'data_with_moon/omega_with_moon_m={m_moon_short}.txt')
f_node = os.path.join(here, f'data_with_moon/orbital_node_with_moon_m={m_moon_short}.txt')
f_ecc = os.path.join(here, f'data_with_moon/ecc_with_moon_m={m_moon_short}.txt')
f_inc = os.path.join(here, f'data_with_moon/inc_with_moon_m={m_moon_short}.txt')
"""

# load
# shape: (nsteps, nplanets)
sma = np.loadtxt(f_sma)        # semi major axis in AU
lamb = np.loadtxt(f_longitude) # mean longitude (in rad)
omega = np.loadtxt(f_omega)    # argument of pericenter (in rad)
Omega = np.loadtxt(f_node)     # longitude of ascending node (in rad)
ecc = np.loadtxt(f_ecc)       # eccentricity
inc = np.loadtxt(f_inc)       # inclination (in deg)
#res_angles = np.loadtxt(f_res_angles)  # resonant angles (n_angles, nsteps) in deg

##########################################################################################
# calculation of resonant angles
##########################################################################################

# Convert from radians to degrees in [0,360)
def rad_to_deg_0_360(arr):
  a = np.mod(arr, 2*np.pi)  # now in [0, 2pi)
  return a * 180.0 / np.pi

def resonant_angles(longitude):
  # return time series arrays (length n_steps) for each resonant angle

  # Leleu equations
  """
  phi0 = 3 * longitude[0, :] - 5 * longitude[1, :]
  phi1 = 1 * longitude[1, :] - 2 * longitude[2, :]
  phi2 = 2 * longitude[2, :] - 3 * longitude[3, :]
  phi3 = 2 * longitude[3, :] - 3 * longitude[4, :]
  phi4 = 3 * longitude[4, :] - 4 * longitude[5, :]

  psi1 = phi1 - phi2
  psi2 = phi2 - phi3
  psi3 = 0.5 * (phi3 - phi4)
  """

  # Murray and Dermott equations
  phi0 = 5 * longitude[:, 1] - 3 * longitude[:, 0] - 2 * (omega[:, 1] + Omega[:, 1])
  phi1 = 2 * longitude[:, 2] - longitude[:, 1] - (omega[:, 2] + Omega[:, 2])
  phi2 = 3 * longitude[:, 3] - 2 * longitude[:, 2] - (omega[:, 3] + Omega[:, 3])
  phi3 = 3 * longitude[:, 4] - 2 * longitude[:, 3] - (omega[:, 4] + Omega[:, 4])
  phi4 = 4 * longitude[:, 5] - 3 * longitude[:, 4] - (omega[:, 5] + Omega[:, 5])

  psi1 = longitude[:, 1] - 4 * longitude[:, 2] + 3 * longitude[:, 3]
  psi2 = 2 * longitude[:, 2] - 5 * longitude[:, 3] + 3 * longitude[:, 4]
  psi3 = longitude[:, 3] - 3 * longitude[:, 4] + 2 * longitude[:, 5]

  phi0 = rad_to_deg_0_360(phi0)
  phi1 = rad_to_deg_0_360(phi1)
  phi2 = rad_to_deg_0_360(phi2)
  phi3 = rad_to_deg_0_360(phi3)
  phi4 = rad_to_deg_0_360(phi4)

  psi1 = rad_to_deg_0_360(psi1)
  psi2 = rad_to_deg_0_360(psi2)
  psi3 = rad_to_deg_0_360(psi3)

  return phi0, phi1, phi2, phi3, phi4, psi1, psi2, psi3


nsteps, nplanets = sma.shape
# time axis from 0 to 500 years
time_years = np.linspace(0, 500, nsteps)

# compute resonant angles
phi0, phi1, phi2, phi3, phi4, psi1, psi2, psi3 = resonant_angles(lamb)
angles = np.array([phi0, phi1, phi2, phi3, phi4, psi1, psi2, psi3])

# number of angles and timesteps
n_angles = angles.shape[0]

# Create a 2x4 grid and plot first n_to_plot angles
fig, axes = plt.subplots(4, 2, figsize=(12, 9), sharex=True)
axes = axes.flatten()

for ax in axes:
  ax.set_ylim(0, 360)
  ax.grid(alpha=0.25)

for i in range(5):
  ax = axes[i]
  ax.plot(time_years, angles[i], lw=0.3)
  ax.set_title(f'$\phi_{i}$')
  ax.set_ylabel('resonant angle [deg]')

for i in range(3):
  ax = axes[i+5]
  ax.plot(time_years, angles[i+5], lw=0.8)
  ax.set_title(f'$\psi_{i+1}$')
  ax.set_ylabel('resonant angle [deg]')

# turn off any unused axes
for j in range(n_angles, len(axes)):
  axes[j].axis('off')

# label shared x-axis on bottom-left subplot (index 6 for our layout)
axes[6].set_xlabel('time [years]')

plt.tight_layout()
out = os.path.join(here, 'plots/all_resonances.png')
#out = os.path.join(here, f'plots_m={m_moon_short}/all_resonances_m={m_moon_short}.png')
plt.savefig(out, dpi=300, bbox_inches='tight')
print(f"Saved {out}")
plt.close()