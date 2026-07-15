import numpy as np
import matplotlib.pyplot as plt
import rebound

sa = rebound.Simulationarchive("sim_with_moon.bin")

sim = sa[0]


ob1 = rebound.OrbitPlot(sim, particles=[1,2,3,4,5,6])
ob2 = rebound.OrbitPlot(sim, particles=[7], primary=5, fig=ob1.fig, ax=ob1.ax, color='red')
plt.gca().set_aspect('equal', 'box')
plt.savefig('plots/test_orbitPlot.png', dpi=300, bbox_inches='tight')
plt.close()
