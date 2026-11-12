import os
import numpy as np
import plot_funtions as pf

m_short = 0.001  # mass of the moon relative to host planet
a_short = 0.3  # semi major axis of the moon relative to host planet

##########################################################################################
# load data
##########################################################################################

here = os.path.dirname(__file__)

# filenames
f_sma = os.path.join(here, f'data_moon_around_b/data_moon_a={a_short}/sma_moon_a={a_short}_m={m_short}.txt')
f_longitude = os.path.join(here, f'data_moon_around_b/data_moon_a={a_short}/l_moon_a={a_short}_m={m_short}.txt')
f_omega = os.path.join(here, f'data_moon_around_b/data_moon_a={a_short}/omega_moon_a={a_short}_m={m_short}.txt')
f_node = os.path.join(here, f'data_moon_around_b/data_moon_a={a_short}/orbital_node_moon_a={a_short}_m={m_short}.txt')
f_ecc = os.path.join(here, f'data_moon_around_b/data_moon_a={a_short}/ecc_moon_a={a_short}_m={m_short}.txt')
f_inc = os.path.join(here, f'data_moon_around_b/data_moon_a={a_short}/inc_moon_a={a_short}_m={m_short}.txt')
f_xyz_host_planet = os.path.join(here, f'data_moon_around_b/data_moon_a={a_short}/xyz_host_planet_a={a_short}_m={m_short}.txt')
f_xyz_moon = os.path.join(here, f'data_moon_around_b/data_moon_a={a_short}/xyz_moon_a={a_short}_m={m_short}.txt')

# load
# shape: (nsteps, nplanets)
sma = np.loadtxt(f_sma)        # semi major axis in AU
lamb = np.loadtxt(f_longitude) # mean longitude (in rad)
omega = np.loadtxt(f_omega)    # argument of pericenter (in rad)
Omega = np.loadtxt(f_node)     # longitude of ascending node (in rad)
ecc = np.loadtxt(f_ecc)       # eccentricity
inc = np.loadtxt(f_inc)       # inclination (in deg)
xyz_host_planet = np.loadtxt(f_xyz_host_planet)  # coordinates of host planet
xyz_moon = np.loadtxt(f_xyz_moon)  # coordinates of moon

phi0, phi1, phi2, phi3, phi4, psi1, psi2, psi3 = pf.resonant_angles(lamb, omega, Omega)
pf.make_resonant_angle_plot(psi1, psi2, psi3, a=a_short, m=m_short)

pf.make_sma_ecc_plot(sma, ecc, a=a_short, m=m_short)

pf.make_moon_xyz_plot(xyz_host_planet=xyz_host_planet, xyz_moon=xyz_moon, a=a_short, m=m_short)