import time
import numpy as np
import matplotlib.pyplot as plt 
import rebound


# constants
M_E=5.972e24
M_S=1.989e30
R_sun = 696340000
AU = 1.5e11
Rstar = 0.651*R_sun

Ndays=500*365.25
orbit_time = 365.25
day_in_second = 60*60*24
Nsteps = 10000
times = np.linspace(0, Ndays*day_in_second, Nsteps)
timestep = (times[2]-times[1])
Nt = 6


# Data from the paper for each planet and of the star
# Masses with error
masses = np.zeros(7)
masses[0] = 0.65*M_S
masses[1] = 1.5*M_E
masses[2] = 4.77*M_E
masses[3] = 3.01*M_E
masses[4] = 3.86*M_E
masses[5] = 7.72*M_E
masses[6] = 3.94*M_E

mass_err1 = np.array([1-(0.44/1.50),1+(0.39/1.50)])
mass_err2 = np.array([1-(0.68/4.77),1+(0.55/4.77)])
mass_err3 = np.array([1-(1.03/3.01),1+(0.80/3.01)])
mass_err4 = np.array([1-(0.94/3.86),1+(1.25/3.86)])
mass_err5 = np.array([1-(1.52/7.72),1+(1.67/7.72)])
mass_err6 = np.array([1-(1.62/3.94),1+(1.31/3.94)])

# semi-major axes
sma = np.zeros(6)
sma[0] = 0.02607*AU
sma[1] = 0.0370*AU
sma[2] = 0.0592*AU
sma[3] = 0.0783*AU
sma[4] = 0.1039*AU
sma[5] = 0.1275*AU

# physical Radius
R = np.zeros(6)
R[0] = 0.01623*Rstar
R[1] = 0.0235*Rstar
R[2] = 0.03623*Rstar
R[3] = 0.0311*Rstar
R[4] = 0.0322*Rstar
R[5] = 0.0404*Rstar

# Hill- radii
Rh = np.zeros(6)
Rh[0] = (masses[1]/(3*masses[0]))**(1/3)*sma[0]
Rh[1] = (masses[2]/(3*masses[0]))**(1/3)*sma[1]
Rh[2] = (masses[3]/(3*masses[0]))**(1/3)*sma[2]
Rh[3] = (masses[4]/(3*masses[0]))**(1/3)*sma[3]
Rh[4] = (masses[5]/(3*masses[0]))**(1/3)*sma[4]
Rh[5] = (masses[6]/(3*masses[0]))**(1/3)*sma[5]

# Periods
P = np.zeros(6)
P[0] = 1.914558*day_in_second
P[1] = 3.23845*day_in_second
P[2] = 6.5577*day_in_second
P[3] = 9.961881*day_in_second
P[4] = 15.231915*day_in_second
P[5] = 20.70950*day_in_second

# Transit time in the observation
T0 = np.zeros(6)
T0[0] = 2458741.6365
T0[1] = 2458741.4783
T0[2] = 2458747.14623
T0[3] = 2458751.4658
T0[4] = 2458745.7178
T0[5] = 2458748.0302

# relative inclinations
incl = np.zeros(6)
incl[0] = np.deg2rad(0.4)
incl[1] = np.deg2rad(0)
incl[2] = np.deg2rad(0.18)
incl[3] = np.deg2rad(0.31)
incl[4] = np.deg2rad(0.323)
incl[5] = np.deg2rad(0.523)

"""
# inc (deg)
incl = np.zeros(6)
incl[0] = np.deg2rad(88.800)
incl[1] = np.deg2rad(88.400)
incl[2] = np.deg2rad(88.580)
incl[3] = np.deg2rad(88.710)
incl[4] = np.deg2rad(88.723)
incl[5] = np.deg2rad(88.823)

# lambda (deg)
lambd = np.zeros(6)
lambd[0] = 1.772610
lambd[1] = 0.632669
lambd[2] = -1.268245
lambd[3] = -0.934120
lambd[4] = 0.207905
lambd[5] = -1.737573
"""

### day zero date in JBD ###
date_ci = 2458354

### Calculation of initial position in orbit for each planet ###
lambd = np.zeros(6)
for i in range(6): 
  lambd[i] = -(2*np.pi/P[i])*((T0[i]-date_ci)*day_in_second)-np.pi/2


def setupSimulation():
# Setting up the Simulation
  sim = rebound.Simulation()
  sim.collision = 'direct'
  sim.collision_resolve = 'merge'
  sim.G = 6.6743e-11
  
  # placing the planets and the star

  TOI178 = rebound.Particle(m=masses[0],r=Rstar)
  sim.add(TOI178)
  TOI178_b = rebound.Particle(simulation=sim,primary=TOI178,m=masses[1], P=P[0], l=lambd[0],r=0.1*Rh[0],inc=incl[0])
  TOI178_c = rebound.Particle(simulation=sim,primary=TOI178,m=masses[2], P=P[1], l=lambd[1],r=0.1*Rh[1],inc=incl[1])
  TOI178_d = rebound.Particle(simulation=sim,primary=TOI178,m=masses[3], P=P[2], l=lambd[2],r=0.1*Rh[2],inc=incl[2])
  TOI178_e = rebound.Particle(simulation=sim,primary=TOI178,m=masses[4], P=P[3], l=lambd[3],r=0.1*Rh[3],inc=incl[3])
  TOI178_f = rebound.Particle(simulation=sim,primary=TOI178,m=masses[5], P=P[4], l=lambd[4],r=0.1*Rh[4],inc=incl[4])
  TOI178_g = rebound.Particle(simulation=sim,primary=TOI178,m=masses[6], P=P[5], l=lambd[5],r=0.1*Rh[5],inc=incl[5])
  
  sim.add(TOI178_b)
  sim.add(TOI178_c)
  sim.add(TOI178_d)
  sim.add(TOI178_e)
  sim.add(TOI178_f)
  sim.add(TOI178_g)
  

  sim.move_to_com()
  
  return sim


def simulation(sim):
### Definig simulation calculation and data entries ###
  N = sim.N
  # Store arrays with shape (nsteps, nplanets) so rows=time, columns=planet
  ecc = np.zeros((Nsteps, Nt))
  sma = np.zeros((Nsteps, Nt))
  inc = np.zeros((Nsteps, Nt))
  #mean_montion = np.zeros((Nsteps, Nt))
  omega = np.zeros((Nsteps, Nt))
  longitude = np.zeros((Nsteps, Nt))
  orbital_node = np.zeros((Nsteps, Nt))
  
  
  for i,t in enumerate(times):
  ### time step and data collection ###
    sim.integrate(t, exact_finish_time=0)
    N = sim.N
    
    ps = sim.particles
    
    for j in range(1,N):
      # store per-time-step in row i, planet index j-1 in column
      ecc[i, j-1] = ps[j].e
      sma[i, j-1] = ps[j].a / AU
      inc[i, j-1] = np.rad2deg(ps[j].inc)
      #mean_montion[i, j-1] = ps[j].n
      omega[i, j-1] = ps[j].pomega
      longitude[i, j-1] = ps[j].l
      orbital_node[i, j-1] = ps[j].Omega
    
   
    print("The time is %5d years "% (t/(60*60*24*365.25)))


  return ecc,sma,inc,omega,longitude,orbital_node


sim = setupSimulation()
ecc,sma,inc,omega,longitude,orbital_node = simulation(sim)

### saving data ###
# Save arrays directly: rows = timesteps, columns = planets
np.savetxt('data_sim/ecc_planets.txt', ecc)
np.savetxt('data_sim/sma_planets.txt', sma)
np.savetxt('data_sim/inc_planets.txt', inc)
np.savetxt('data_sim/orbital_node_planets.txt', orbital_node)
np.savetxt('data_sim/omega_planets.txt', omega)
np.savetxt('data_sim/l_planets.txt', longitude)


rebound.OrbitPlot(sim)
plt.gca().set_aspect('equal', 'box')
plt.savefig('plots/orbit_plot_rebound.png', dpi=300, bbox_inches='tight')
plt.close()