import sim_moon_auto
import resonance_plot_auto

start_m = 0.001
end_m = 0.013
start_a = 0.3
end_a = 0.4

i = start_a
while i <= end_a:
    j = start_m
    while j<= end_m:
        sim = sim_moon_auto.setupSimulation(a=i, m=j)
        ecc,sma,inc,omega,longitude,orbital_node,xyz_f,xyz_moon = sim_moon_auto.simulation(sim)
        sim_moon_auto.safe_data(ecc, sma, inc, omega, longitude, orbital_node, xyz_f, xyz_moon, a=i, m=j)
        phi0, phi1, phi2, phi3, phi4, psi1, psi2, psi3 = resonance_plot_auto.resonant_angles(longitude, omega, orbital_node)
        resonance_plot_auto.make_plot(phi0, phi1, phi2, phi3, phi4, psi1, psi2, psi3, a=i, m=j)
        j += 0.001
        j = round(j, 3)
    i += 0.1
    i = round(i, 1)