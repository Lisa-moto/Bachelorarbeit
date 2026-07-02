import numpy as np
import matplotlib.pyplot as plt
import os

m_moon_short = 0.003  # mass of the moon relative to planet f
a_moon_short = 0.2  # semi-major axis of the moon relative to planet f

moon_info = (f"(moon: $a={a_moon_short}\\,R_{{\\text{{Hill,f}}}}$, "
             f"$m={m_moon_short}\\,m_{{\\text{{f}}}}$)")

##########################################################################################
# load data
##########################################################################################

here = os.path.dirname(__file__)

# filenames
f_sma = os.path.join(here, f'data_aufgabe1_auto/data_moon_a={a_moon_short}/sma_moon_a={a_moon_short}_m={m_moon_short}.txt')
f_longitude = os.path.join(here, f'data_aufgabe1_auto/data_moon_a={a_moon_short}/l_moon_a={a_moon_short}_m={m_moon_short}.txt')
f_omega = os.path.join(here, f'data_aufgabe1_auto/data_moon_a={a_moon_short}/omega_moon_a={a_moon_short}_m={m_moon_short}.txt')
f_node = os.path.join(here, f'data_aufgabe1_auto/data_moon_a={a_moon_short}/orbital_node_moon_a={a_moon_short}_m={m_moon_short}.txt')
f_ecc = os.path.join(here, f'data_aufgabe1_auto/data_moon_a={a_moon_short}/ecc_moon_a={a_moon_short}_m={m_moon_short}.txt')
f_inc = os.path.join(here, f'data_aufgabe1_auto/data_moon_a={a_moon_short}/inc_moon_a={a_moon_short}_m={m_moon_short}.txt')


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

# Create a 1x3 grid
fig, axes = plt.subplots(3, 1, figsize=(12, 12), sharex=True)
axes = axes.flatten()

# Plot semi major axis
axes[0].plot(time_years, sma[:, 6], color='blue', zorder=10)
axes[0].set_ylabel('Semi Major Axis [AU]', fontsize=14)
axes[0].set_title(f'Semi Major Axes {moon_info}', fontsize=16)
axes[0].grid()
# y-Achse von 0 bis 0.2 AU
axes[0].set_ylim(0, 0.002371) # Hill Radius of planet f

#axes[0].legend(fontsize=12)

# Plot eccentricity
axes[1].plot(time_years, ecc[:, 6], color='red', marker='.', linestyle='none', markersize=2)
axes[1].set_ylabel('Eccentricity', fontsize=14)
axes[1].set_xlabel('Time [years]', fontsize=14)
axes[1].set_title('Eccentricity of the moon', fontsize=16)
axes[1].grid()
#axes[1].legend(fontsize=12)

# Plot inclination
axes[2].plot(time_years, inc[:, 6], color='green', marker='.', linestyle='none', markersize=2)
axes[2].set_ylabel('Inclination [deg]', fontsize=14)
axes[2].set_xlabel('Time [years]', fontsize=14)
axes[2].set_title('Inclination of the moon', fontsize=16)
axes[2].grid()
#axes[2].legend(fontsize=12)

plt.tight_layout()

##########################################################################################
# Save plot
##########################################################################################
output_dir = f'plots_a={a_moon_short}/plots_m={m_moon_short}'
os.makedirs(output_dir, exist_ok=True)

filename = f'moon_inc_sma_ecc_a={a_moon_short}_m={m_moon_short}.png'
plt.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches='tight')
print(f"Saved {filename}")
plt.close()

# print out the mean sma of the mooon in AU and km
mean_sma_moon = np.mean(sma[:, 6])
print(f"Mean semi major axis of the moon: {mean_sma_moon:.6f} AU")
au_to_km = 149597870.7
mean_sma_moon_km = mean_sma_moon * au_to_km
print(f"Mean semi major axis of the moon: {mean_sma_moon_km:.2f} km")