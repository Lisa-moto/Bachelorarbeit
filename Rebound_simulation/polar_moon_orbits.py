import numpy as np
import matplotlib.pyplot as plt 
import rebound
import os


# constants
M_E=5.972e24
M_S=1.989e30
R_sun = 696340000
AU = 1.5e11
Rstar = 0.651*R_sun

Ndays=2*365.25
orbit_time = 365.25
day_in_second = 60*60*24
Nsteps = 10000
times = np.linspace(0, Ndays*day_in_second, Nsteps)
timestep = (times[2]-times[1])
Nt = 7 # 6 planets + 1 moon


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
sma[1] = 0.03700*AU
sma[2] = 0.05920*AU
sma[3] = 0.07830*AU
sma[4] = 0.10390*AU
sma[5] = 0.12750*AU

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
P[1] = 3.238450*day_in_second
P[2] = 6.557700*day_in_second
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

# moon parameters
a_moon = 0.2*Rh[4]    # semi-major axis of the moon around planet f
r_moon = 0.001*Rh[4]  # radius of the moon
moon_mass = 0.001*masses[5]
inc_moon = np.deg2rad(90) # inclination of the moon's orbit
a_moon_short = a_moon/Rh[4]
m_moon_short = round(moon_mass/masses[5], 3)

### day zero date in JBD ###
date_ci = 2458354

### Calculation of initial position in orbit for each planet ###
lambd = np.zeros(6)
for i in range(6): 
  lambd[i] = -(2*np.pi/P[i])*((T0[i]-date_ci)*day_in_second)-np.pi/2


def setupSimulation(a=a_moon_short, m=m_moon_short):
# Setting up the Simulation
  sim = rebound.Simulation()
  sim.collision = 'direct'
  sim.collision_resolve = 'merge'
  sim.G = 6.6743e-11

  a_moon = a*Rh[4]
  moon_mass = m*masses[5]
  
  # placing the planets and the star
  TOI178 = rebound.Particle(m=masses[0], r=Rstar, hash="star")
  sim.add(TOI178)
  TOI178_b = rebound.Particle(simulation=sim, primary=TOI178, m=masses[1], P=P[0], l=lambd[0], r=0.1*Rh[0], inc=incl[0], hash="b")
  TOI178_c = rebound.Particle(simulation=sim, primary=TOI178, m=masses[2], P=P[1], l=lambd[1], r=0.1*Rh[1], inc=incl[1], hash="c")
  TOI178_d = rebound.Particle(simulation=sim, primary=TOI178, m=masses[3], P=P[2], l=lambd[2], r=0.1*Rh[2], inc=incl[2], hash="d")
  TOI178_e = rebound.Particle(simulation=sim, primary=TOI178, m=masses[4], P=P[3], l=lambd[3], r=0.1*Rh[3], inc=incl[3], hash="e")
  TOI178_f = rebound.Particle(simulation=sim, primary=TOI178, m=masses[5], P=P[4], l=lambd[4], r=0.1*Rh[4], inc=incl[4], hash="f")
  TOI178_g = rebound.Particle(simulation=sim, primary=TOI178, m=masses[6], P=P[5], l=lambd[5], r=0.1*Rh[5], inc=incl[5], hash="g")
  
  sim.add(TOI178_b); sim.add(TOI178_c); sim.add(TOI178_d)
  sim.add(TOI178_e); sim.add(TOI178_f); sim.add(TOI178_g)
  
  # moon
  TOI178_f_moon = rebound.Particle(
      simulation=sim, primary=TOI178_f, m=moon_mass, r=r_moon, a=a_moon, inc=inc_moon, hash="moon"
  )
  sim.add(TOI178_f_moon)

  sim.move_to_com()
  return sim


def simulation(sim):
  N = sim.N
  ecc = np.zeros((Nsteps, Nt))
  sma = np.zeros((Nsteps, Nt))
  inc = np.zeros((Nsteps, Nt))
  omega = np.zeros((Nsteps, Nt))
  longitude = np.zeros((Nsteps, Nt))
  orbital_node = np.zeros((Nsteps, Nt))
  xyz_f = np.zeros((Nsteps,3))
  xyz_moon = np.zeros((Nsteps,3))

  # --- moon orbit counting: count true anomaly crossings ---
  prev_f = None  # true anomaly
  moon_periapsis_count = 0
  moon_ejection_time = None
  moon_period_initial_days = np.nan
  ejection_idx = -1

  for i, t in enumerate(times):
    sim.integrate(t, exact_finish_time=0)
    ps = sim.particles
    N = sim.N

    ecc[i, :] = np.nan
    sma[i, :] = np.nan
    inc[i, :] = np.nan
    omega[i, :] = np.nan
    longitude[i, :] = np.nan
    orbital_node[i, :] = np.nan
    xyz_f[i, :] = np.nan
    xyz_moon[i, :] = np.nan

    for j in range(1, N):
      if (j-1) < Nt:
        ecc[i, j-1] = ps[j].e
        sma[i, j-1] = ps[j].a / AU
        inc[i, j-1] = np.rad2deg(ps[j].inc)
        omega[i, j-1] = ps[j].pomega
        longitude[i, j-1] = ps[j].l
        orbital_node[i, j-1] = ps[j].Omega

    p_f = None
    p_moon = None
    try:
      p_f = sim.particles["f"]
      xyz_f[i, :] = [p_f.x, p_f.y, p_f.z]
    except Exception:
      pass

    try:
      p_moon = sim.particles["moon"]
      xyz_moon[i, :] = [p_moon.x, p_moon.y, p_moon.z]
    except Exception:
      if moon_ejection_time is None:
        moon_ejection_time = t
        ejection_idx = i

    if (p_f is not None) and (p_moon is not None):
      o = p_moon.orbit(primary=p_f)
      ecc[i, 6] = o.e
      sma[i, 6] = o.a / AU

      # first valid osculating period
      if np.isnan(moon_period_initial_days) and (o.a > 0) and (o.e < 1):
        moon_period_initial_days = o.P / day_in_second

      # count periapsis passages (true anomaly f crossing 0)
      if (o.a > 0) and (o.e < 1):
        f = o.f  # true anomaly
        
        if prev_f is not None:
          # detect crossing from negative to positive (periapsis at f=0)
          if (prev_f < 0) and (f >= 0):
            moon_periapsis_count += 1
          elif (prev_f > np.pi) and (f < np.pi):
            # crossing from apoapsis region back to periapsis
            if f < prev_f - np.pi/2:  # real crossing, not noise
              moon_periapsis_count += 1
        
        prev_f = f
      else:
        # unbound
        if moon_ejection_time is None:
          moon_ejection_time = t
          ejection_idx = i

    print("The time is %5d years " % (t/(60*60*24*365.25)))

  if moon_ejection_time is None:
    moon_ejection_time = times[-1]
    ejection_idx = Nsteps - 1

  # orbits = periapsis counts, refined from mean anomaly near ejection
  moon_orbits = moon_periapsis_count
  moon_time_days = moon_ejection_time / day_in_second
  moon_period_mean_days = moon_time_days / moon_orbits if moon_orbits > 0 else np.nan

  moon_stats = {
    "orbits_until_ejection": moon_orbits,
    "ejection_time_days": moon_time_days,
    "period_initial_days": moon_period_initial_days,
    "period_mean_until_ejection_days": moon_period_mean_days,
    "ejection_idx": ejection_idx,
  }

  return ecc, sma, inc, omega, longitude, orbital_node, xyz_f, xyz_moon, moon_stats

sim = setupSimulation(a=a_moon_short, m=m_moon_short)

ecc, sma, inc, omega, longitude, orbital_node, xyz_f, xyz_moon, moon_stats = simulation(sim)

print("Anzahl der particles nach Simulation: ", sim.N)
print("Moon stats:")
print(f"  Orbits until ejection: {moon_stats['orbits_until_ejection']:.4f}")
print(f"  Ejection time [days]:  {moon_stats['ejection_time_days']:.4f}")
print(f"  Initial period [days]: {moon_stats['period_initial_days']:.6f}")
print(f"  Mean period until ejection [days]: {moon_stats['period_mean_until_ejection_days']:.6f}")

# ...existing code...

# Nach dem simulation() call:
fig, axes = plt.subplots(2, 1, figsize=(10, 6))

# Moon eccentricity over time
time_years = times / (day_in_second * 365.25)
axes[0].plot(time_years, ecc[:, 6], 'r-', linewidth=1)
axes[0].axvline(moon_stats['ejection_time_days']/365.25, color='k', linestyle='--', label='Ejection')
axes[0].set_ylabel('Moon Eccentricity')
axes[0].set_xlabel('Time [years]')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Moon semi-major axis over time
axes[1].plot(time_years, sma[:, 6], 'b-', linewidth=1)
axes[1].axvline(moon_stats['ejection_time_days']/365.25, color='k', linestyle='--', label='Ejection')
axes[1].set_ylabel('Moon SMA [AU]')
axes[1].set_xlabel('Time [years]')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('moon_destabilization.png', dpi=150)