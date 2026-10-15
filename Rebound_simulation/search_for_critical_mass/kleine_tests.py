import numpy as np

shape = np.load("/pfs/10/project/bw16g002/Lisa/Masterarbeit/critical_mass_search_TEST/simulation_data/a=0.251754/a=0.251754_m=1.000000e-09.npz")["ecc"].shape
print(f"Shape of the eccentricity array: {shape}")