import numpy as np
import matplotlib.pyplot as plt
import os

# Convert from radians to degrees in [0,360)
def rad_to_deg_0_360(arr):
  a = np.mod(arr, 2*np.pi)  # now in [0, 2pi)
  return a * 180.0 / np.pi

def resonant_angles(longitude, omega, Omega):
  # return time series arrays (length n_steps) for each resonant angle

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

def make_plot(phi0, phi1, phi2, phi3, phi4, psi1, psi2, psi3, a, m):
  n_steps = len(phi0)
  time_years = np.linspace(0, 500, n_steps)
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

  output_dir = f'plots_aufgabe1_auto/plots_a={a}'
  os.makedirs(output_dir, exist_ok=True)

  plt.tight_layout()
  plt.savefig(f'{output_dir}/all_resonances_a={a}_m={m}.png', dpi=300, bbox_inches='tight')
  print(f"Saved all_resonances_a={a}_m={m}")
  plt.close()