import numpy as np
import matplotlib.pyplot as plt
import os

m_moon_short = 0.004  # mass of the moon relative to planet f

##########################################################################################
# load data
##########################################################################################

here = os.path.dirname(__file__)

# filenames
f_sma = os.path.join(here, f'data_with_moon/sma_with_moon_m={m_moon_short}.txt')
f_longitude = os.path.join(here, f'data_with_moon/l_with_moon_m={m_moon_short}.txt')
f_omega = os.path.join(here, f'data_with_moon/omega_with_moon_m={m_moon_short}.txt')
f_node = os.path.join(here, f'data_with_moon/orbital_node_with_moon_m={m_moon_short}.txt')
f_ecc = os.path.join(here, f'data_with_moon/ecc_with_moon_m={m_moon_short}.txt')
f_inc = os.path.join(here, f'data_with_moon/inc_with_moon_m={m_moon_short}.txt')

# load
# shape: (nsteps, nplanets)
sma = np.loadtxt(f_sma)        # semi major axis in AU
lamb = np.loadtxt(f_longitude) # mean longitude (in rad)
omega = np.loadtxt(f_omega)    # argument of pericenter (in rad)
Omega = np.loadtxt(f_node)     # longitude of ascending node (in rad)
ecc = np.loadtxt(f_ecc)       # eccentricity
inc = np.loadtxt(f_inc)       # inclination (in deg)


##########################################################################################
# plotting semi major axis and eccentricity e of the moon
##########################################################################################
nsteps, nplanets = sma.shape
# time axis from 0 to 500 years
time_years = np.linspace(0, 500, nsteps)

# Create a 1x2 grid
fig, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
axes = axes.flatten()

# Plot semi major axis
axes[0].plot(time_years, sma[:, 6], color='blue', marker='.', linestyle='none', markersize=2, label='sma of moon')
axes[0].set_ylabel('Semi Major Axis [AU]', fontsize=14)
axes[0].set_title('Semi Major Axis of the moon over time', fontsize=16)
axes[0].grid()
# y-Achse von 0 bis 0.2 AU
axes[0].set_ylim(0, 0.2)
# zum Vergleich plot der sma von Planet f und g
axes[0].plot(time_years, sma[:, 4], color='green', label='sma of Planet f')
axes[0].plot(time_years, sma[:, 5], color='orange', label='sma of Planet g')
axes[0].legend(fontsize=12)

# Plot eccentricity
axes[1].plot(time_years, ecc[:, 6], color='red', marker='.', linestyle='none', markersize=2, label='eccentricity of moon')
axes[1].set_ylabel('Eccentricity', fontsize=14)
axes[1].set_xlabel('Time [years]', fontsize=14)
axes[1].set_title('Eccentricity of the moon over time', fontsize=16)
axes[1].grid()
axes[1].legend(fontsize=12)

plt.tight_layout()
out = os.path.join(here, f'plots_m={m_moon_short}/moon_sma_ecc_m={m_moon_short}.png')
plt.savefig(out, dpi=300, bbox_inches='tight')
print(f"Saved {out}")
plt.close()