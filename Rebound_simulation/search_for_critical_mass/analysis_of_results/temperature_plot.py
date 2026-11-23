import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, LogFormatterSciNotation
from pathlib import Path


# -----------------------------------------------------------------------------
# Pfade und Konstanten
# -----------------------------------------------------------------------------

BASE_DIR = Path("/pfs/10/project/bw16g002/Lisa/Masterarbeit/critical_mass_search")

LOG_DIR = BASE_DIR / "progress_logs"                 # Ordner mit den Logs
LOG_PATTERN = "progress_a={a:03d}.log"                # Dateiname pro a_index
LOG_DELIM = ";"                             # Trennzeichen im Log

N_a = 200
sma_values = np.linspace(0.001, 0.5, N_a)

ANGLES = ("psi1", "psi2", "psi3")


def load_all_logs():
    """Alle (a, m)-Kombinationen aus den Logs: a_index, sma, mass, status, psi1-3."""
    rows = []
    for a in range(N_a):
        files = list(LOG_DIR.glob(LOG_PATTERN.format(a=a)))
        if not files:
            continue
        with open(files[0], newline="") as f:
            for r in csv.reader(f, delimiter=LOG_DELIM):
                if len(r) < 6 or not r[0].strip().isdigit():
                    continue  # Header/Leerzeilen überspringen
                def num(x):
                    try:
                        return float(x)
                    except ValueError:   # "None" oder "-"
                        return np.nan
                rows.append(dict(
                    a_index=a, sma=sma_values[a], mass=float(r[1]), status=r[2],
                    psi1=num(r[3]), psi2=num(r[4]), psi3=num(r[5]),
                ))
    return pd.DataFrame(rows)


def plot_break_years(df, save_path=None, log_color=False,
                     cmap="viridis", c_nobreak="orange", c_unstable="grey"):
    """
    Drei Panels (psi1, psi2, psi3): sma (x) gegen Masse (y, log).
    Farbe = Jahr des Bruchs; kein Bruch und instabil haben eigene Farben.
    Nicht simulierte Kombinationen (Abbruch nach Bruch aller drei) bleiben leer.
    """
    ok = df["status"] == "ok"
    unstable = df["status"] == "unstable"

    # gemeinsame Farbskala für alle drei Panels
    years = np.concatenate([df.loc[ok, a].dropna().values for a in ANGLES])
    vmin, vmax = years.min(), years.max()
    norm = LogNorm(vmin=max(vmin, 1e-3), vmax=vmax) if log_color \
        else Normalize(vmin=vmin, vmax=vmax)

    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True, sharey=True,
                             constrained_layout=True)

    for ax, angle, lab in zip(axes, ANGLES, (r"$\Psi_1$", r"$\Psi_2$", r"$\Psi_3$")):
        broke = ok & df[angle].notna()
        nobreak = ok & df[angle].isna()

        ax.scatter(df.loc[nobreak, "sma"], df.loc[nobreak, "mass"],
                   s=6, marker="s", color=c_nobreak, linewidths=0, rasterized=True)
        ax.scatter(df.loc[unstable, "sma"], df.loc[unstable, "mass"],
                   s=8, marker="x", color=c_unstable, linewidths=0.5, rasterized=True)
        sc = ax.scatter(df.loc[broke, "sma"], df.loc[broke, "mass"],
                c=df.loc[broke, angle], cmap=cmap, norm=norm,
                s=6, marker="s", linewidths=0, rasterized=True)

        ax.set_yscale("log")
        ax.yaxis.set_major_locator(LogLocator(base=10, numticks=20))
        ax.yaxis.set_major_formatter(LogFormatterSciNotation(base=10))
        ax.set_ylabel(r"mass in $m_{\text{f}}$")
        ax.set_title(f"Resonance breaking {lab}")
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel(r"semi major axis in $R_{\text{Hill,f}}$")

    cbar = fig.colorbar(sc, ax=axes, pad=0.02)
    cbar.set_label("year of resonance break")

    handles = [
        Line2D([], [], marker="s", ls="", color=c_nobreak, label="no break"),
        Line2D([], [], marker="x", ls="", color=c_unstable, label="unstable"),
    ]
    axes[0].legend(handles=handles, loc="upper left", fontsize=8, framealpha=0.9)

    if save_path is not None:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
    return fig, axes


df_all = load_all_logs()
fig, axes = plot_break_years(
    df_all, save_path="Plots_critical_mass/break_year_temperature_plot.png"
)