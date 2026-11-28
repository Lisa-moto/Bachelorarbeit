import re
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import sys
#parent_dir = Path(__file__).resolve().parents[1]

SCRIPT_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPT_DIR) not in sys.path:
	sys.path.insert(0, str(SCRIPT_DIR))

import functions

BASE_DIR = Path("/pfs/10/project/bw16g002/Lisa/Masterarbeit/critical_mass_search")

START_A = 0.001
END_A = 0.5
N_A = 200
SMA_VALUES = np.linspace(START_A, END_A, N_A)   # exakt wie in der Simulation

# matcht a=0.100000_m=1.000000e-05.npz, aber NICHT ...tmp.npz
FILE_RE = re.compile(r"a=(?P<a>[\d.]+)_m=(?P<m>[\d.]+e[+-]\d+)\.npz$")


def resolve_a(a=None, a_index=None, base_dir=BASE_DIR):
    """
    Bestimmt das semi major axis a aus dem Simulationsgitter.
      a_index : 1-basierter Index (1 bis N_A)
      a       : beliebiger Wert -> nächstliegender Gitterwert
    Genau eines von beiden muss angegeben werden.
    """
    if (a is None) == (a_index is None):
        raise ValueError("Bitte genau eines von 'a' oder 'a_index' angeben.")

    if a_index is not None:
        if not 1 <= a_index <= N_A:
            raise ValueError(f"a_index muss zwischen 1 und {N_A} liegen.")
        idx = a_index - 1
    else:
        idx = int(np.argmin(np.abs(SMA_VALUES - a)))

    a_used = float(SMA_VALUES[idx])
    print(f"a_index={idx + 1}, a={a_used:.6f}"
          + (f" (nächstliegend zu angefragtem a={a})" if a is not None else ""))

    folder = Path(base_dir) / "simulation_data" / f"a={a_used:.6f}"
    if not folder.exists():
        raise FileNotFoundError(f"Ordner existiert (noch) nicht: {folder}")
    return a_used

# ----------------------------------------------------------------------------
# Dateien finden und Massen auswählen
# ----------------------------------------------------------------------------
def list_available_masses(a, base_dir=BASE_DIR):
    """Gibt sortierte Liste [(m, path), ...] für eine semi major axis a zurück."""
    folder = Path(base_dir) / "simulation_data" / f"a={a:.6f}"
    if not folder.exists():
        raise FileNotFoundError(f"Ordner nicht gefunden: {folder}")
    out = []
    for p in folder.glob("a=*_m=*.npz"):
        match = FILE_RE.match(p.name)
        if match:
            out.append((float(match.group("m")), p))
    out.sort(key=lambda x: x[0])
    return out


def select_masses(available, m_min=None, m_max=None, step=1, masses=None):
    """
    Wählt Massen aus.
      masses : explizite Liste -> jeweils nächstliegende vorhandene Masse (log)
      m_min, m_max : Bereich (in Einheiten von m_f), inklusive
      step : nimmt jede step-te Masse aus dem (gefilterten) Bereich
    """
    if masses is not None:
        m_avail = np.array([m for m, _ in available])
        chosen = []
        for m_req in masses:
            idx = int(np.argmin(np.abs(np.log10(m_avail) - np.log10(m_req))))
            if available[idx] not in chosen:
                chosen.append(available[idx])
        return sorted(chosen, key=lambda x: x[0])

    sel = [(m, p) for m, p in available
           if (m_min is None or m >= m_min) and (m_max is None or m <= m_max)]
    return sel[::step]


# ----------------------------------------------------------------------------
# Daten laden und Winkel berechnen
# ----------------------------------------------------------------------------
def rad_to_deg_0_360(arr):
    return np.mod(arr, 2 * np.pi) * 180.0 / np.pi


def compute_psi(longitude):
    """longitude: (n, 6) für b, c, d, e, f, g (ohne Jahr-Spalte, ohne Mond)."""
    psi1 = longitude[:, 1] - 4 * longitude[:, 2] + 3 * longitude[:, 3]
    psi2 = 2 * longitude[:, 2] - 5 * longitude[:, 3] + 3 * longitude[:, 4]
    psi3 = longitude[:, 3] - 3 * longitude[:, 4] + 2 * longitude[:, 5]
    return [rad_to_deg_0_360(p) for p in (psi1, psi2, psi3)]


def load_psi(path, window_years=300):
    """
    Lädt eine Simulationsdatei und gibt (years, [psi1, psi2, psi3]) zurück.
    window_years: nur die letzten N Jahre der gespeicherten Daten; None = alles.
    """
    d = np.load(path)
    l = d["l"].astype(np.float64)        # float64 für die Winkelkombination
    years = l[:, 0]
    lon = l[:, 1:7]                      # b..g, Mond (Spalte 7) ignorieren

    psis = compute_psi(lon)

    if window_years is not None:
        mask = years >= years.max() - window_years
        years = years[mask]
        psis = [p[mask] for p in psis]
    return years, psis


# ----------------------------------------------------------------------------
# Frame zeichnen
# ----------------------------------------------------------------------------
def render_frame(years, psis, a, m, t_end=5000.0, figsize=(12, 7), dpi=100,
                 threshold_deg=179.0):
    fig, axes = plt.subplots(3, 1, figsize=figsize, dpi=dpi, sharex=True)

    for i, ax in enumerate(axes):
        ax.plot(years, psis[i], ls="none", marker=".", ms=2)
        ax.set_ylim(0, 360)
        ax.set_yticks([0, 90, 180, 270, 360])
        ax.grid(alpha=0.25)
        ax.set_ylabel("angle (deg)", size=12)

        resonant = functions.is_laplace_resonant(psis[i], threshold_deg=threshold_deg)
        if resonant:
            ax.set_title(rf"$\psi_{i+1}$, resonant", size=14, color="black")
        else:
            ax.set_title(rf"$\psi_{i+1}$, not resonant", size=14, color="red")
    axes[-1].set_xlabel("time (years)", size=12)

    early = years.max() < t_end - 1
    note = f"  (early stop, last data at {years.max():.0f} yr)" if early else ""
    fig.suptitle(
        rf"Resonant angles, moon around planet f; "
        rf"$a={a:.4f}\,R_{{\mathrm{{Hill,f}}}}$, $m={m:.3e}\,m_{{\mathrm{{f}}}}${note}",
        size=15,
    )
    fig.tight_layout()

    fig.canvas.draw()
    frame = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
    plt.close(fig)
    return frame


# ----------------------------------------------------------------------------
# Video
# ----------------------------------------------------------------------------
def make_video(a=None, a_index=None, out_dir=None, m_min=None, m_max=None,
               step=1, masses=None, window_years=300, seconds_per_frame=1.0,
               fps=10, t_end=5000.0, base_dir=BASE_DIR):
    a = resolve_a(a, a_index, base_dir)          # <-- neu
    available = list_available_masses(a, base_dir)
    chosen = select_masses(available, m_min, m_max, step, masses)
    if not chosen:
        raise ValueError("Keine Massen ausgewählt / keine Dateien gefunden.")
    print(f"{len(chosen)} von {len(available)} Massen werden geplottet.")

    out_dir = Path(out_dir) if out_dir is not None else SCRIPT_DIR / "analysis_of_results/resonance_videos"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"video_psi_a={a:.6f}.mp4"

    n_repeat = max(1, int(round(seconds_per_frame * fps)))
    writer = imageio.get_writer(str(out_path), fps=fps, codec="libx264",
                                macro_block_size=1)
    try:
        for k, (m, path) in enumerate(chosen, 1):
            years, psis = load_psi(path, window_years)
            frame = render_frame(years, psis, a, m, t_end=t_end)
            for _ in range(n_repeat):
                writer.append_data(frame)
            print(f"[{k}/{len(chosen)}] m={m:.3e}")
    finally:
        writer.close()
    print(f"Gespeichert: {out_path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", type=float, default=None)          # nächstliegendes a
    ap.add_argument("--a_index", type=int, default=None)      # 1 bis 200
    ap.add_argument("--m_min", type=float, default=None)
    ap.add_argument("--m_max", type=float, default=None)
    ap.add_argument("--step", type=int, default=1)
    ap.add_argument("--masses", type=float, nargs="+", default=None)
    ap.add_argument("--window_years", type=float, default=300)  # 0 = alles
    ap.add_argument("--seconds_per_frame", type=float, default=1.0)
    ap.add_argument("--out_dir", type=str, default=None)  # Standard: Skript-Ordner
    args = ap.parse_args()

    make_video(
        a=args.a, a_index=args.a_index, m_min=args.m_min, m_max=args.m_max, step=args.step,
        masses=args.masses,
        window_years=None if args.window_years == 0 else args.window_years,
        seconds_per_frame=args.seconds_per_frame, out_dir=args.out_dir,
    )