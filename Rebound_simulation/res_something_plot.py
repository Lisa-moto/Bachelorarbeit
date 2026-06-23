import matplotlib.pyplot as plt
import numpy as np
import os

##########################################################################################
# Load data
##########################################################################################

here = os.path.dirname(__file__)

# filenames
f_sma = os.path.join(here, 'data_sim/sma_planets.txt')
f_ecc = os.path.join(here, 'data_sim/ecc_planets.txt')
f_longitude = os.path.join(here, 'data_sim/l_planets.txt')
f_omega = os.path.join(here, 'data_sim/omega_planets.txt')
f_node = os.path.join(here, 'data_sim/orbital_node_planets.txt')
f_inc = os.path.join(here, 'data_sim/inc_planets.txt')

# load
# shape: (nsteps, nplanets)
sma = np.loadtxt(f_sma)        # semi major axis in AU
ecc = np.loadtxt(f_ecc)       # eccentricity
lamb = np.loadtxt(f_longitude) # mean longitude (in rad)
omega = np.loadtxt(f_omega)    # argument of pericenter (in rad)
Omega = np.loadtxt(f_node)     # longitude of ascending node (in rad)
inc = np.loadtxt(f_inc)       # inclination (in deg)

orbital_parameters = {
    'Semi Major Axis (AU)':              sma,
    'Eccentricity':                      ecc,
    'Mean Longitude (rad)':              lamb,
    'Argument of Pericenter (rad)':      omega,
    'Longitude of Ascending Node (rad)': Omega,
    'Inclination (deg)':                 inc,
}

##########################################################################################
# Choose what to plot!
##########################################################################################
# What should be plotted? -> sma, ecc, lamb, omega, Omega oder inc
subject = ecc
# Which linestyle makes sense?
linestyle = ':'
# Which resonance is being considered? Psi 1, 2 or 3
psi_select = 2
# Which time period should be plotted? (in years, at most 500)
time_start = 300
time_end = 400

# 'single' -> ein Plot mit allen Planeten
# 'subplots' -> 6 kleine Plots, je ein Planet
plot_mode = 'subplots'

##########################################################################################
# Helper functions
##########################################################################################
def get_label(data_array, param_dict):
    """Gives the label corresponding to the data_array from the param_dict"""
    for name, arr in param_dict.items():
        if arr is data_array:
            return name
    return 'Unknown'

# Convert from radians to degrees in [0,360)
def rad_to_deg_0_360(arr):
  a = np.mod(arr, 2*np.pi)  # now in [0, 2pi)
  return a * 180.0 / np.pi

# Calculate resonant angles
def resonant_angles(longitude):
  # return time series arrays (length n_steps) for each resonant angle

  psi1 = longitude[:, 1] - 4 * longitude[:, 2] + 3 * longitude[:, 3]
  psi2 = 2 * longitude[:, 2] - 5 * longitude[:, 3] + 3 * longitude[:, 4]
  psi3 = longitude[:, 3] - 3 * longitude[:, 4] + 2 * longitude[:, 5]

  psi1 = rad_to_deg_0_360(psi1)
  psi2 = rad_to_deg_0_360(psi2)
  psi3 = rad_to_deg_0_360(psi3)

  return psi1, psi2, psi3

##########################################################################################
# Calculation/Preparation
##########################################################################################
psi1, psi2, psi3 = resonant_angles(lamb)
psi_arrays = {1: psi1, 2: psi2, 3: psi3}
psi_data   = psi_arrays[psi_select]

subject_label = get_label(subject, orbital_parameters)

# Zeitachse: Zeilenindex = Jahr
n_steps = subject.shape[0]
time    = np.linspace(0, 500, n_steps)
mask    = (time >= time_start) & (time <= time_end)

planet_labels = ['Planet b', 'Planet c', 'Planet d', 'Planet e', 'Planet f', 'Planet g']
colors        = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']


##########################################################################################
# Plot
##########################################################################################
def plot_single():
    fig, ax1 = plt.subplots(figsize=(12, 6))

    for i, (label, color) in enumerate(zip(planet_labels, colors)):
        ax1.plot(time[mask], subject[mask, i], label=label, color=color, linestyle=linestyle)

    ax1.set_xlabel('Time (years)')
    ax1.set_ylabel(subject_label)
    ax1.grid(alpha=0.3)
    ax1.legend(loc='upper left')

    ax2 = ax1.twinx()
    ax2.scatter(time[mask], psi_data[mask], s=1, color='black', alpha=0.4,
                label=rf'$\psi_{psi_select}$')
    ax2.set_ylabel(rf'$\psi_{psi_select}$ (deg)')
    ax2.set_ylim(0, 360)
    ax2.set_yticks([0, 90, 180, 270, 360])
    ax2.legend(loc='upper right')

    fig.suptitle(f'{subject_label}  &  $\\psi_{psi_select}$  over time\n')
    plt.tight_layout()
    return fig

def plot_subplots():
    fig, axes = plt.subplots(3, 2, figsize=(15, 8), sharex=True)
    axes = axes.flatten()
    
    #subject_max = subject[mask].max()
    #print(subject_max)
    #ylim_subject = (0, subject_max * 1.05)

    for i, (label, color) in enumerate(zip(planet_labels, colors)):
        ax1 = axes[i]
        ax1.plot(time[mask], subject[mask, i], color=color, linestyle=linestyle)
        ax1.set_title(label)
        ax1.set_ylabel(subject_label)
        ax1.grid(alpha=0.3)
        #ax1.set_ylim(ylim_subject)

        # x-Achsenbeschriftung nur in unterster Reihe
        if i >= 3:
            ax1.set_xlabel('Time (years)')

        ax2 = ax1.twinx()
        ax2.scatter(time[mask], psi_data[mask], s=1, color='black', alpha=0.4,
                    label=rf'$\psi_{psi_select}$')
        ax2.set_ylim(0, 360)
        ax2.set_yticks([0, 90, 180, 270, 360])

        # y-Achsenbeschriftung für psi nur rechts außen (Spalte 1)
        if i in (1, 3, 5):
            ax2.set_ylabel(rf'$\psi_{psi_select}$ (deg)')
        else:
            ax2.set_yticklabels([])

    fig.suptitle(f'{subject_label}  &  $\\psi_{psi_select}$  over time\n', y=1.01)
    plt.tight_layout()
    return fig

if plot_mode == 'single':
    fig = plot_single()
elif plot_mode == 'subplots':
    fig = plot_subplots()
else:
    raise ValueError(f"plot_mode must be 'single' or 'subplots', got '{plot_mode}'")

##########################################################################################
# Save plot
##########################################################################################
output_dir = f'resonance_break_plots'
os.makedirs(output_dir, exist_ok=True)

filename = f'{subject_label}_{plot_mode}_psi{psi_select}_{time_start}-{time_end}.png'
plt.savefig(os.path.join(output_dir, filename), dpi=300, bbox_inches='tight')
print(f"Saved {filename}")
plt.close()