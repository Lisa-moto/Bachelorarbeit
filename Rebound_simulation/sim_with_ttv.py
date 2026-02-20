import numpy as np
import matplotlib.pyplot as plt 
import rebound


# constants
M_E=5.972e24
M_S=1.989e30
R_sun = 696340000
AU = 1.5e11
Rstar = 0.651*R_sun

Ndays=20*365.25
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


### day zero date in JBD ###
date_ci = 2458354

### Calculation of initial position in orbit for each planet ###
lambd = np.zeros(6)
for i in range(6): 
  lambd[i] = -(2*np.pi/P[i])*((T0[i]-date_ci)*day_in_second)-np.pi/2

# function for calculating the TTV
def compute_ttv(sim, n_transits, planet_index, check_step, post_step):
    N = n_transits
    transit_times = np.zeros(N)
    ps = sim.particles
    star = ps[0]
    planet = ps[planet_index]
    i = 0
    while i<N and sim.t < tmax:
        y_old = planet.y - star.y
        t_old = sim.t
        sim.integrate(sim.t+check_step) # integrate check step to check for transit
        t_new = sim.t
        if y_old*(planet.y-star.y)<0. and planet.x-star.x>0.: # sign changed (y_old*y<0), planet in front of star (x>0)
            while t_new - t_old > 1e-5: # bisect until prec of 1e-5 reached
                if y_old*(planet.y-star.y)<0.:
                    t_new = sim.t
                else:
                    t_old = sim.t
                sim.integrate((t_new+t_old)/2.)
            transit_times[i] = sim.t
            i += 1
            sim.integrate(sim.t+post_step)       # integrate post_step to be past the transit
            print(f"Found transit {i} at time %5.2f days" % (transit_times[i-1]/day_in_second))
    return transit_times

def setupSimulation():
# Setting up the Simulation
  sim = rebound.Simulation()
  sim.collision = 'direct'
  sim.collision_resolve = 'merge'
  sim.G = 6.6743e-11

  #sim.integrator = "whfast" # eigentlich nicht fuer close encounters und wenn Kollisionen auftreten
  #sim.dt = 0.001 * P[0] # fuer TTV-Genauigkeit kleineres timestep als fuer Stabilitaet
  
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


sim = setupSimulation()
tmax = Ndays * day_in_second
n_transits = int(tmax/P[4])
check_step = 0.002*P[4] # check for a transit in that interval
post_step = 0.005*P[4] # after finding a transit, integrate post_step to be past the transit to avoid finding the same transit again
transit_times = compute_ttv(sim, n_transits, 5, check_step, post_step)


# linear least square fit to remove the linear trend from the transit times
N = len(transit_times)
A = np.vstack([np.ones(N), range(N)]).T
c, m = np.linalg.lstsq(A, transit_times, rcond=-1)[0]

n = np.arange(len(transit_times))
ttv = (transit_times-m*np.array(range(N))-c) # in seconds

# save transit times to file
np.savetxt("data_ttv/transit_times.txt", ttv)

# fuer 10 Jahre 239 transits

# Plot TTVs
plt.figure(figsize=(8,5))
plt.plot(n, ttv/60, lw=0.5, marker='.')
plt.xlabel("Transit number")
plt.ylabel("TTV [minutes]")
plt.title("TTV of TOI-178 f")
plt.savefig("TTV_plots/TTV_TOI178_f.png", dpi=300)