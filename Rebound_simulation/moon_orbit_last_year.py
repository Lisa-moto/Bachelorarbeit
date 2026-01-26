import numpy as np
import matplotlib.pyplot as plt

m_moon_short = 0.004  # mass of the moon relative to planet f

# -----------------------------
# load data
# -----------------------------
xyz_f = np.loadtxt(f"data_with_moon/xyz_f_with_moon_m={m_moon_short}.txt")
xyz_moon = np.loadtxt(f"data_with_moon/xyz_moon_m={m_moon_short}.txt")

# -----------------------------
# compute moon orbit in planet f frame
# -----------------------------
x_rel = xyz_moon[:, 0] - xyz_f[:, 0]
y_rel = xyz_moon[:, 1] - xyz_f[:, 1]
z_rel = xyz_moon[:, 2] - xyz_f[:, 2]

# only plot last year
Nsteps = xyz_f.shape[0]
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
axes[0].set_xlabel("x [m]")
axes[0].set_ylabel("y [m]")
axes[0].set_title("x–y projection")
axes[0].set_xlim(-lim, lim)
axes[0].set_ylim(-lim, lim)
axes[0].set_aspect("equal")
axes[0].grid(True)

# xz plot
axes[1].scatter(x, z, s=4, alpha=0.7)
axes[1].scatter(0, 0, s=60, color="black")
axes[1].set_xlabel("x [m]")
axes[1].set_ylabel("z [m]")
axes[1].set_title("x–z projection")
axes[1].set_xlim(-lim, lim)
axes[1].set_ylim(-lim, lim)
axes[1].set_aspect("equal")
axes[1].grid(True)

# yz plot
axes[2].scatter(y, z, s=4, alpha=0.7)
axes[2].scatter(0, 0, s=60, color="black")
axes[2].set_xlabel("y [m]")
axes[2].set_ylabel("z [m]")
axes[2].set_title("y–z projection")
axes[2].set_xlim(-lim, lim)
axes[2].set_ylim(-lim, lim)
axes[2].set_aspect("equal")
axes[2].grid(True)

plt.suptitle(f"Moon orbit around TOI-178 f (last {N_years} simulation years)", fontsize=14)
plt.tight_layout()
plt.savefig(f"plots_m={m_moon_short}/moon_orbit_3projections_last_{N_years}_years_m={m_moon_short}.png", dpi=300, bbox_inches="tight")
plt.close()