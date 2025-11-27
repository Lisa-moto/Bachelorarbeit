import numpy as np
import matplotlib.pyplot as plt
import os


##########################################################################################
# load data
##########################################################################################

here = os.path.dirname(__file__)

# filenames
f_sma = os.path.join(here, 'sma_planets.txt')
f_longitude = os.path.join(here, 'l_planets.txt')
f_omega = os.path.join(here, 'omega_planets.txt')
f_node = os.path.join(here, 'orbital_node_planets.txt')
f_ecc = os.path.join(here, 'ecc_planets.txt')
f_inc = os.path.join(here, 'inc_planets.txt')
#f_res_angles = os.path.join(here, 'resonant_angles_planets.txt')

# load
sma = np.loadtxt(f_sma)        # shape (nplanets, nsteps) in AU
lamb = np.loadtxt(f_longitude) # mean longitude (rad)
omega = np.loadtxt(f_omega)    # argument of pericenter (rad)
Omega = np.loadtxt(f_node)     # longitude of ascending node (rad)
ecc = np.loadtxt(f_ecc)       # eccentricity
inc = np.loadtxt(f_inc)       # inclination (deg)
#res_angles = np.loadtxt(f_res_angles)  # resonant angles (n_angles, nsteps) in deg

##########################################################################################
# calculation of resonant angles
##########################################################################################

# Convert from radians to degrees in [0,360)
def rad_to_deg_0_360(arr):
  a = np.mod(arr, 2*np.pi)  # now in [0, 2pi)
  return a * 180.0 / np.pi

def resonant_angles(longitude):
  # longitude has shape (n_planets, n_timesteps)
  # build time series (use second axis = time)

  # Leleu equations
  #phi0 = 3 * longitude[0, :] - 5 * longitude[1, :]
  #phi1 = 1 * longitude[1, :] - 2 * longitude[2, :]
  #phi2 = 2 * longitude[2, :] - 3 * longitude[3, :]
  #phi3 = 2 * longitude[3, :] - 3 * longitude[4, :]
  #phi4 = 3 * longitude[4, :] - 4 * longitude[5, :]

  #psi1 = phi1 - phi2
  #psi2 = phi2 - phi3
  #psi3 = 0.5 * (phi3 - phi4)

  # Murray and0 Dermott equations
  phi0 = 5 * longitude[1, :] - 3 * longitude[0, :] - 2* (omega[1, :] + Omega[1, :])
  phi1 = 2 * longitude[2, :] - longitude[1, :] - (omega[2, :] + Omega[2, :])
  phi2 = 3 * longitude[3, :] - 2 * longitude[2, :] - (omega[3, :] + Omega[3, :])
  phi3 = 3 * longitude[4, :] - 2 * longitude[3, :] - (omega[4, :] + Omega[4, :])
  phi4 = 4 * longitude[5, :] - 3 * longitude[4, :] - (omega[5, :] + Omega[5, :])

  psi1 = longitude[1, :] - 4 * longitude[2, :] + 3 * longitude[3, :]
  psi2 = 2 * longitude[2, :] - 5 * longitude[3, :] + 3 * longitude[4, :]
  psi3 = longitude[3, :] - 3 * longitude[4, :] + 2 * longitude[5, :]


  phi0 = rad_to_deg_0_360(phi0)
  phi1 = rad_to_deg_0_360(phi1)
  phi2 = rad_to_deg_0_360(phi2)
  phi3 = rad_to_deg_0_360(phi3)
  phi4 = rad_to_deg_0_360(phi4)

  psi1 = rad_to_deg_0_360(psi1)
  psi2 = rad_to_deg_0_360(psi2)
  psi3 = rad_to_deg_0_360(psi3)

  return phi0, phi1, phi2, phi3, phi4, psi1, psi2, psi3


nplanets, nsteps = sma.shape
# time axis from 0 to 500 years
time_years = np.linspace(0, 500, nsteps)

# compute resonant angles
phi0, phi1, phi2, phi3, phi4, psi1,psi2,psi3 = resonant_angles(lamb)
angles = np.array([phi0, phi1, phi2, phi3, phi4, psi1,psi2,psi3])

# number of angles and timesteps
n_angles, n_steps = angles.shape

# Create a 2x4 grid and plot first n_to_plot angles
fig, axes = plt.subplots(4, 2, figsize=(12, 9), sharex=True)
axes = axes.flatten()

for ax in axes:
  ax.set_ylim(0, 360)
  ax.grid(alpha=0.25)

for i in range(n_angles):
  ax = axes[i]
  ax.plot(time_years, angles[i], lw=0.8)
  ax.set_title(f'angle {i}')
  ax.set_ylabel('deg')

# turn off any unused axes
for j in range(n_angles, len(axes)):
  axes[j].axis('off')

# label shared x-axis on bottom-left subplot (index 6 for our layout)
axes[6].set_xlabel('time [years]')

plt.tight_layout()
out = os.path.join(here, 'all_resonances.png')
plt.savefig(out, dpi=300, bbox_inches='tight')
print(f"Saved {out}")
plt.close()