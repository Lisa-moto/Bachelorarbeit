#!/usr/bin/env python3
"""Generic resonance checker.

This script reads the simulation output files (sma, longitudes, nodes, omegas)
and can either scan planet pairs or test a single pair for candidate mean-motion
resonances. For each candidate m:n it computes common resonant angles and can plot them
for visual inspection.

Usage examples:
  # scan adjacent pairs, max denominator 6, produce plots
  python check_resonance.py --scan --max-den 6 --plot --outdir results

  # test pair (1,2) meaning planets with indices 1 and 2 (0-based); 0=b, 1=c, ...
  python check_resonance.py --pair 1 2 --max-den 8 --plot

Files expected (in same folder):
  sma_planets.txt, l_planets.txt, omega_planets.txt, orbital_node_planets.txt

Outputs: summary printed to stdout and optional PNG plots saved to --outdir.
"""
import os
import argparse
from fractions import Fraction
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load_2d(fn):
    a = np.loadtxt(fn)
    if a.ndim == 1:
        a = a[np.newaxis, :]
    return a


def best_rationals(x, max_den=8):
    """Return a list of (m,n,approx) candidate rationals close to x.

    We return the single best rational approximation using Fraction.limit_denominator.
    """
    f = Fraction(x).limit_denominator(max_den)
    return [(f.numerator, f.denominator, f.numerator / f.denominator)]


def compute_varpi(omega_arr, Omega_arr):
    return Omega_arr + omega_arr


def compute_phi_angles(lambda_inner, lambda_outer, varpi_inner, varpi_outer, m, n):
    """Compute two common resonant angle variants (rad):
    phi_outer = m*lambda_outer - n*lambda_inner - (m-n)*varpi_outer
    phi_inner = m*lambda_outer - n*lambda_inner - (m-n)*varpi_inner
    """
    phi_outer = m * lambda_outer - n * lambda_inner - (m - n) * varpi_outer
    phi_inner = m * lambda_outer - n * lambda_inner - (m - n) * varpi_inner
    return phi_outer, phi_inner


# The script computes and plots resonant angles; visual inspection
# of the plotted angles is used to assess resonance.


def pair_name(i, j):
    # planet labels: index 0 -> 'b', 1 -> 'c', ...
    def label(idx):
        return chr(ord('b') + idx)
    return f"{label(i)}-{label(j)}"


def run_pair_check(i, j, sma, lamb, omega, Omega, time_years, max_den=6, outdir=None, do_plot=True):
    nplanets = sma.shape[0]
    if i < 0 or j < 0 or i >= nplanets or j >= nplanets:
        raise ValueError('planet indices out of range')
    # inner = i, outer = j? We'll treat i as inner, j as outer if sma[i] < sma[j], else swap
    a_i = np.mean(sma[i])
    a_j = np.mean(sma[j])
    if a_i > a_j:
        inner, outer = j, i
    else:
        inner, outer = i, j

    a_inner = np.mean(sma[inner])
    a_outer = np.mean(sma[outer])
    period_ratio = (a_outer / a_inner) ** 1.5

    candidates = best_rationals(period_ratio, max_den=max_den)

    # prepare arrays (lambda in rad, varpi in rad)
    lambda_arr_inner = lamb[inner]
    lambda_arr_outer = lamb[outer]
    varpi_inner = compute_varpi(omega[inner], Omega[inner])
    varpi_outer = compute_varpi(omega[outer], Omega[outer])

    results = []
    for m, n, approx in candidates:
        phi_out, phi_in = compute_phi_angles(lambda_arr_inner, lambda_arr_outer, varpi_inner, varpi_outer, m, n)
        results.append({'m': m, 'n': n, 'approx': approx, 'phi_out': phi_out, 'phi_in': phi_in})

    # Print summary
    pair = pair_name(inner, outer)
    print(f"Pair {pair}: mean a_inner={a_inner:.6f} AU, a_outer={a_outer:.6f} AU, P_ratio~{period_ratio:.6f}")
    print("Candidates:")
    for r in results:
        print(f"  {r['m']}:{r['n']} (~{r['approx']:.4f})")

    # optional plotting: show sma and best phi (choose variant with smaller amp)
    if do_plot:
        fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
        ax[0].plot(time_years, sma[inner], label=f'planet inner ({pair.split("-")[0]})')
        ax[0].plot(time_years, sma[outer], label=f'planet outer ({pair.split("-")[1]})')
        ax[0].set_ylabel('a [AU]')
        ax[0].legend()
        ax[0].grid(alpha=0.3)

        # pick best candidate (first) and plot both angle variants (in degrees, 0-360)
        m, n, _ = candidates[0]
        phi_out, phi_in = compute_phi_angles(lambda_arr_inner, lambda_arr_outer, varpi_inner, varpi_outer, m, n)
        phi_out_deg = np.mod(np.degrees(phi_out), 360.0)
        phi_in_deg = np.mod(np.degrees(phi_in), 360.0)

        ax[1].plot(time_years, phi_out_deg, label=f'phi_out {m}:{n}')
        ax[1].set_ylabel('phi_out [deg]')
        ax[1].set_ylim(0, 360)
        ax[1].grid(alpha=0.3)
        ax[1].legend()

        ax[2].plot(time_years, phi_in_deg, label=f'phi_in {m}:{n}', color='C1')
        ax[2].set_ylabel('phi_in [deg]')
        ax[2].set_ylim(0, 360)
        ax[2].set_xlabel('time [years]')
        ax[2].grid(alpha=0.3)
        ax[2].legend()

        plt.tight_layout()
        if outdir is None:
            outdir = os.getcwd()
        os.makedirs(outdir, exist_ok=True)
        out = os.path.join(outdir, f'resonance_{pair}_{m}to{n}.png')
        plt.savefig(out, dpi=300)
        plt.close(fig)
        print(f"Saved plot to {out}")

    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pair', nargs=2, type=int, help='Two planet indices (0-based) to test, e.g. 1 2')
    parser.add_argument('--scan', action='store_true', help='Scan adjacent planet pairs')
    parser.add_argument('--max-den', type=int, default=6, help='Max denominator for rational approximation')
    parser.add_argument('--plot', action='store_true', help='Save plots for checked pairs')
    parser.add_argument('--outdir', default='resonance_results', help='Directory to save plots')
    args = parser.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    f_sma = os.path.join(here, 'sma_planets.txt')
    f_l = os.path.join(here, 'l_planets.txt')
    f_omega = os.path.join(here, 'omega_planets.txt')
    f_node = os.path.join(here, 'orbital_node_planets.txt')

    # load
    sma = load_2d(f_sma)
    lamb = load_2d(f_l)
    omega = load_2d(f_omega)
    Omega = load_2d(f_node)

    nplanets, nsteps = sma.shape
    time_years = np.linspace(0, 500, nsteps)

    checks = []
    if args.pair:
        i, j = args.pair
        checks.append((i, j))
    elif args.scan:
        # adjacent pairs by default
        for k in range(nplanets - 1):
            checks.append((k, k+1))
    else:
        parser.error('Either --pair or --scan must be specified')

    for (i, j) in checks:
        try:
            run_pair_check(i, j, sma, lamb, omega, Omega, time_years, max_den=args.max_den, outdir=args.outdir, do_plot=args.plot)
        except Exception as e:
            print(f"Error checking pair {i},{j}: {e}")


if __name__ == '__main__':
    main()
