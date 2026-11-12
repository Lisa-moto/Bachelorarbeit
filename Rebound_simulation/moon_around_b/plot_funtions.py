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

def make_resonant_angle_plot(psi1, psi2, psi3, a, m):
    n_steps = len(psi1)
    time_years = np.linspace(0, 500, n_steps)
    angles = np.array([psi1, psi2, psi3])

    # number of angles and timesteps
    n_angles = angles.shape[0]

    # Create a 1x3 grid and plot first n_to_plot angles
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    axes = axes.flatten()

    for ax in axes:
        ax.set_ylim(0, 360)
        ax.grid(alpha=0.25)

    for i in range(3):
        ax = axes[i]
        ax.plot(time_years, angles[i], lw=0.8)
        ax.set_title(f'$\psi_{i+1}$', size=16)
        ax.set_ylabel('angle (deg)', size=15)

    # turn off any unused axes
    for j in range(n_angles, len(axes)):
        axes[j].axis('off')

    # label shared x-axis on bottom
    axes[2].set_xlabel('time (years)', size=15)
    
    # title
    fig.suptitle(f'Resonant Angles with Moon around Planet b; $a={a}\\,R_{{\\text{{Hill,f}}}}$, $m={m}\\,m_{{\\text{{f}}}}$', size=18)

    output_dir = f'plots_moon_around_b/a={a}'
    os.makedirs(output_dir, exist_ok=True)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/resonant_angles_a={a}_m={m}.png', dpi=300, bbox_inches='tight')
    print(f"Saved resonant_angles_a={a}_m={m}")
    plt.close()
  
def make_sma_ecc_plot(sma, ecc, a, m):
    n_steps, n_planets = sma.shape
    time_years = np.linspace(0, 500, n_steps)

    fig, axes = plt.subplots(2, 1, figsize=(12, 9))
    axes = axes.flatten()

    # Plot semi major axis
    axes[0].plot(time_years, sma[:, 6], color='blue', zorder=10, label='sma of moon')
    axes[0].set_ylabel('Semi Major Axis (AU)', fontsize=14)
    axes[0].set_xlabel('Time (years)', fontsize=14)
    axes[0].set_title(f'Semi Major Axes', fontsize=16)
    axes[0].grid()
    axes[0].set_ylim(0, 0.2)
    axes[0].plot(time_years, sma[:, 0], color='purple', label='sma of Planet b')
    axes[0].plot(time_years, sma[:, 1], color='green', label='sma of Planet c')
    axes[0].plot(time_years, sma[:, 2], color='orange', zorder=1, label='sma of Planet d')
    axes[0].legend(fontsize=12)

    # Plot eccentricity
    axes[1].plot(time_years, ecc[:, 6], color='red', marker='.', linestyle='none', markersize=2, label='eccentricity of moon')
    axes[1].set_ylabel('Eccentricity', fontsize=14)
    axes[1].set_xlabel('Time (years)', fontsize=14)
    axes[1].set_title('Eccentricity of the moon', fontsize=16)
    axes[1].grid()
    axes[1].legend(fontsize=12)
    
    plt.suptitle(f'Moon around Planet b; $a={a}\\,R_{{\\text{{Hill,f}}}}$, $m={m}\\,m_{{\\text{{f}}}}$', fontsize=18)

    output_dir = f'plots_moon_around_b/a={a}'
    os.makedirs(output_dir, exist_ok=True)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/sma_ecc_a={a}_m={m}.png', dpi=300, bbox_inches='tight')
    print(f"Saved sma_ecc_a={a}_m={m}")
    plt.close()
    
    
def make_moon_xyz_plot(xyz_host_planet, xyz_moon, a, m):
    # -----------------------------
    # compute moon orbit in frame of host planet
    # -----------------------------
    x_rel = xyz_moon[:, 0] - xyz_host_planet[:, 0]
    y_rel = xyz_moon[:, 1] - xyz_host_planet[:, 1]
    z_rel = xyz_moon[:, 2] - xyz_host_planet[:, 2]

    # only plot last N years
    Nsteps = xyz_host_planet.shape[0]
    steps_per_year = int(Nsteps/500) # as total simulation time is 500 years
    N_years = 5 # last N years
    x = x_rel[-steps_per_year*N_years:]
    y = y_rel[-steps_per_year*N_years:]
    z = z_rel[-steps_per_year*N_years:]

    # -----------------------------
    # Gemeinsame Achsengrenzen bestimmen
    # -----------------------------
    max_extent = np.max(np.abs([x, y, z]))
    lim = 1.05 * max_extent   # kleiner Rand

    # -----------------------------
    # Plot
    # -----------------------------
    fig, axes = plt.subplots(1, 3, figsize=(15,6))

    # xy plot
    axes[0].scatter(x, y, s=4, alpha=0.7)
    axes[0].scatter(0, 0, s=60, color="black")
    axes[0].set_xlabel("$x$ (m)", size=15)
    axes[0].set_ylabel("$y$ (m)", size=15)
    axes[0].set_title("$x$–$y$ projection", size=17)
    axes[0].set_xlim(-lim, lim)
    axes[0].set_ylim(-lim, lim)
    axes[0].set_aspect("equal")
    axes[0].grid(True)

    # xz plot
    axes[1].scatter(x, z, s=4, alpha=0.7)
    axes[1].scatter(0, 0, s=60, color="black")
    axes[1].set_xlabel("$x$ (m)", size=15)
    axes[1].set_ylabel("$z$ (m)", size=15)
    axes[1].set_title("$x$–$z$ projection", size=17)
    axes[1].set_xlim(-lim, lim)
    axes[1].set_ylim(-lim, lim)
    axes[1].set_aspect("equal")
    axes[1].grid(True)

    # yz plot
    axes[2].scatter(y, z, s=4, alpha=0.7)
    axes[2].scatter(0, 0, s=60, color="black")
    axes[2].set_xlabel("$y$ (m)", size=15)
    axes[2].set_ylabel("$z$ (m)", size=15)
    axes[2].set_title("$y$–$z$ projection", size=17)
    axes[2].set_xlim(-lim, lim)
    axes[2].set_ylim(-lim, lim)
    axes[2].set_aspect("equal")
    axes[2].grid(True)

    plt.suptitle(f"Moon Orbit around TOI-178b (last {N_years} simulation years)", fontsize=18)
    plt.tight_layout()
    plt.savefig(f"plots_moon_around_b/a={a}/moon_orbit_a={a}_m={m}.png", dpi=300, bbox_inches="tight")
    print(f"Saved moon_orbit_a={a}_m={m}")
    plt.close()