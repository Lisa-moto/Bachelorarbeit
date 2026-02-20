import numpy as np
import matplotlib.pyplot as plt
import os

m_moon_short = 0.001  # mass of the moon relative to planet f

# load transit times from file
here = os.path.dirname(__file__)

f_transit_times = os.path.join(here, "data_ttv/transit_times.txt")
f_transit_times_moon = os.path.join(here, f"data_ttv/transit_times_moon_m={m_moon_short}.txt")

transit_times = np.loadtxt(f_transit_times)
transit_times_moon = np.loadtxt(f_transit_times_moon)

# find the biggest ttv difference at a transit
ttv_diff = np.abs(transit_times - transit_times_moon)
max_diff = np.max(ttv_diff)
print(f"Max TTV difference between with and without moon: {max_diff/60:.2f} minutes")
# find the standard deviation of the ttv difference
std_diff = np.std(ttv_diff)
print(f"Standard deviation of TTV difference between with and without moon: {std_diff/60:.2f} minutes")

n = np.arange(len(transit_times))
# Plot TTVs
plt.figure(figsize=(8,5))
plt.plot(n, transit_times/60, lw=0.8, label="without moon")
plt.plot(n, transit_times_moon/60, lw=0.8, label="with moon")
plt.xlabel("transit number")
plt.ylabel("TTV [minutes]")
plt.title("TTV comparison with and without moon")
plt.legend()

# save plot
plt.savefig(f"TTV_plots/ttv_comparison_m={m_moon_short}.png", dpi=300)