import numpy as np
import matplotlib.pyplot as plt
import os

AU = 1.5e11

##########################################################################################
# Choose moon parameters!
##########################################################################################
m_moon_short = 0.003  # mass of the moon relative to mass of planet f
a_moon_short = 0.2  # semi major axis of the moon relative to planet f

##########################################################################################
# load data
##########################################################################################

here = os.path.dirname(__file__)

# filenames
f_sma = os.path.join(here, f'data/data_with_moon/data_moon_a={a_moon_short}/sma_moon_a={a_moon_short}_m={m_moon_short}.txt')
f_longitude = os.path.join(here, f'data/data_with_moon/data_moon_a={a_moon_short}/l_moon_a={a_moon_short}_m={m_moon_short}.txt')
f_omega = os.path.join(here, f'data/data_with_moon/data_moon_a={a_moon_short}/omega_moon_a={a_moon_short}_m={m_moon_short}.txt')
f_node = os.path.join(here, f'data/data_with_moon/data_moon_a={a_moon_short}/orbital_node_moon_a={a_moon_short}_m={m_moon_short}.txt')
f_ecc = os.path.join(here, f'data/data_with_moon/data_moon_a={a_moon_short}/ecc_moon_a={a_moon_short}_m={m_moon_short}.txt')
f_inc = os.path.join(here, f'data/data_with_moon/data_moon_a={a_moon_short}/inc_moon_a={a_moon_short}_m={m_moon_short}.txt')

f_sma_without_moon = os.path.join(here, f'data/data_without_moon/sma_planets.txt')
f_longitude_without_moon = os.path.join(here, f'data/data_without_moon/l_planets.txt')
f_omega_without_moon = os.path.join(here, 'data/data_without_moon/omega_planets.txt')
f_node_without_moon = os.path.join(here, 'data/data_without_moon/orbital_node_planets.txt')
f_ecc_without_moon = os.path.join(here, 'data/data_without_moon/ecc_planets.txt')
f_inc_without_moon = os.path.join(here, 'data/data_without_moon/inc_planets.txt')

# load
# shape: (nsteps, nplanets)
sma = np.loadtxt(f_sma)        # semi major axis in AU
longitude = np.loadtxt(f_longitude) # mean longitude (in rad)
omega = np.loadtxt(f_omega)    # argument of pericenter (in rad)
Omega = np.loadtxt(f_node)     # longitude of ascending node (in rad)
ecc = np.loadtxt(f_ecc)       # eccentricity
inc = np.loadtxt(f_inc)       # inclination (in deg)

sma_without_moon = np.loadtxt(f_sma_without_moon)  # semi major axis of planets without moon in AU
longitude_without_moon = np.loadtxt(f_longitude_without_moon)  # mean longitude without moon in rad
omega_without_moon = np.loadtxt(f_omega_without_moon)  # argument of pericenter without moon in rad
Omega_without_moon = np.loadtxt(f_node_without_moon)     # longitude of ascending node without moon in rad
ecc_without_moon = np.loadtxt(f_ecc_without_moon)  # eccentricity of planets without moon
inc_without_moon = np.loadtxt(f_inc_without_moon)  # inclination of planets without moon in deg


orbital_parameters = {
    'Semi Major Axis (AU)':              sma,
    'Eccentricity':                      ecc,
    'Mean Longitude (rad)':              longitude,
    'Argument of Pericenter (rad)':      omega,
    'Longitude of Ascending Node (rad)': Omega,
    'Inclination (deg)':                 inc,
}

orbital_parameters_without_moon = {
    'Semi Major Axis (AU)':              sma_without_moon,
    'Eccentricity':                      ecc_without_moon,
    'Mean Longitude (rad)':              longitude_without_moon,
    'Argument of Pericenter (rad)':      omega_without_moon,
    'Longitude of Ascending Node (rad)': Omega_without_moon,
    'Inclination (deg)':                 inc_without_moon,
}

short_names = {
    'sma':              sma,
    'ecc':              ecc,
    'longitude':        longitude,
    'omega':            omega,
    'Omega':            Omega,
    'inc':              inc
}

##########################################################################################
# Choose what to plot!
##########################################################################################
# What should be plotted? -> sma, ecc, longitude, omega, Omega oder inc
subject = longitude
# Which linestyle makes sense?
linestyle = ':'
# Which time period should be plotted? (in years, at most 500)
time_start = 0
time_end = 10
# Which planet should be plotted? (0-5 for planets b-g)
planet_index = 0


##########################################################################################
# Helper functions
##########################################################################################
def get_label(data_array, param_dict):
    """Gives the label corresponding to the data_array from the param_dict"""
    for name, arr in param_dict.items():
        if arr is data_array:
            return name
    return 'Unknown'


##########################################################################################
# Calculation/Preparation
##########################################################################################
subject_label = get_label(subject, orbital_parameters)
subject_short_name = get_label(subject, short_names)
subject_without_moon = orbital_parameters_without_moon[subject_label]

# Zeitachse: Zeilenindex = Jahr
n_steps = subject.shape[0]
time    = np.linspace(0, 500, n_steps)
mask    = (time >= time_start) & (time <= time_end)

planet_labels = ['Planet f']

moon_info = (f"(moon: $a={a_moon_short}\\,R_{{\\text{{Hill,f}}}}$, "
             f"$m={m_moon_short}\\,m_{{\\text{{f}}}}$)")


# Stuff for calculating the difference of the semi major axis of planet f with and without moon
a0_planet_f = 0.1039*AU
sim_diff_with_moon = a0_planet_f - sma[:, 4]*AU  # difference of semi major axis of planet f from initial value in m
sim_diff_without_moon = a0_planet_f - sma_without_moon[:, 4]*AU  # difference of semi major axis of planet f from initial value in m


##########################################################################################
# Plotting the chosen parameter (subject) for a planet with and without moon
##########################################################################################

def make_plot():
    # Create a 1x2 grid with shared y-axis
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True, sharey=True)
    axes = axes.flatten()
    
    # Plot subject with moon
    axes[0].plot(time[mask], subject[mask, planet_index], color='blue', linestyle=linestyle)
    axes[0].set_ylabel(subject_label, fontsize=14)
    axes[0].set_title(f'{subject_label} of {planet_labels[planet_index]} with Moon {moon_info}', fontsize=16)
    axes[0].grid()
    
    # Plot subject without moon
    axes[1].plot(time[mask], subject_without_moon[mask, planet_index], color='orange', linestyle=':')
    axes[1].set_ylabel(subject_label, fontsize=14)
    axes[1].set_xlabel('Time (years)', fontsize=14)
    axes[1].set_title(f'{subject_label} of {planet_labels[planet_index]} without Moon', fontsize=16)
    axes[1].grid()

    plt.tight_layout()
    
    return fig


fig = make_plot()

##########################################################################################
# Save plot
##########################################################################################
output_dir = f'plots/plots_a={a_moon_short}/plots_m={m_moon_short}'
os.makedirs(output_dir, exist_ok=True)

filename = f'{subject_short_name}_difference_a={a_moon_short}_m={m_moon_short}_{time_start}-{time_end}.png'
plt.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches='tight')
print(f"Saved {filename}")
plt.close()