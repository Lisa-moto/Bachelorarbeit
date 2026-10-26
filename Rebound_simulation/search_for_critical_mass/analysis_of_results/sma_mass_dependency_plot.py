import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import sys

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
	sys.path.insert(0, str(parent_dir))
 
import functions


# -----------------------------------------------------------------------------
# Pfade und Konstanten
# -----------------------------------------------------------------------------

BASE_DIR = Path("/pfs/10/project/bw16g002/Lisa/Masterarbeit/critical_mass_search")

masses_to_test = functions.build_mass_array(points_per_decade=100)
# mass grid nach m={m:.6e}
MASS_GRID = [round(float(m), 6) for m in masses_to_test]
print(f"Mass grid: {MASS_GRID[:3]}")

N_a = 200
sma_values = np.linspace(0.001, 0.5, N_a)
SMA_GRID = [round(float(sma), 6) for sma in sma_values]


# Namen der drei Resonanzwinkel und erwartete CSV-Spalten

ANGLES = ("psi1", "psi2", "psi3")
EXPECTED_COLUMNS = [
    "sma",
    "psi1_break_mass", "psi1_break_year",
    "psi2_break_mass", "psi2_break_year",
    "psi3_break_mass", "psi3_break_year",
]

# Schlüssel in den .npz-Dateien
NPZ_KEYS = ("ecc", "sma", "inc", "orbital_node", "omega", "l",
                  "xyz_f", "xyz_moon")


# -----------------------------------------------------------------------------
# Einlesen
# -----------------------------------------------------------------------------

df = pd.read_csv(BASE_DIR / "resonance_breaking_ALL.csv", na_values="None")

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------

def plot_resonance_breaking(df, save_path=None, log_y=True):
    """
    Plottet die kritische Masse, bei der die Resonanz bricht, gegen die
    große Halbachse (sma) für psi1, psi2, psi3.

    Punkte ohne Resonanzbruch (None/NaN) werden in einem schmalen
    Zusatzplot unterhalb des Hauptplots als farbige Striche markiert.
    """
    sma = df["sma"].values

    colors  = {"psi1": "tab:blue", "psi2": "tab:orange", "psi3": "tab:green"}
    markers = {"psi1": "o", "psi2": "s", "psi3": "^"}
    labels  = {"psi1": r"$\Psi_1$", "psi2": r"$\Psi_2$", "psi3": r"$\Psi_3$"}

    fig, (ax, ax_none) = plt.subplots(
        2, 1, sharex=True, figsize=(10, 6.5),
        gridspec_kw={"height_ratios": [6, 1], "hspace": 0.08},
    )

    for i, angle in enumerate(ANGLES):
        mass = df[f"{angle}_break_mass"].values.astype(float)
        has_break = np.isfinite(mass)

        # Hauptplot: nur Punkte, keine Verbindungslinien
        ax.scatter(
            sma[has_break], mass[has_break],
            s=18, marker=markers[angle], color=colors[angle],
            alpha=0.8, label=labels[angle], linestyle="none",
        )

        # Zusatzplot: Striche dort, wo kein Bruch auftritt
        ax_none.scatter(
            sma[~has_break], np.full((~has_break).sum(), i),
            s=40, marker="|", color=colors[angle], linewidths=1.5,
        )

    # Hauptplot formatieren
    if log_y:
        ax.set_yscale("log")
    ax.set_ylabel("critical mass in $m_{\text{f}}$")
    ax.set_title("Resonance Breaking: Dependency on Semi Major Axis and Mass of the Moon")
    ax.grid(True, alpha=0.3)
    ax.legend(title="Resonanzwinkel", loc="best")
    ax.tick_params(labelbottom=False)

    # Zusatzplot formatieren
    ax_none.set_yticks(range(len(ANGLES)))
    ax_none.set_yticklabels([labels[a] for a in ANGLES])
    ax_none.set_ylim(-0.7, len(ANGLES) - 0.3)
    ax_none.set_xlabel("semi major axis in $R_{\text{Hill,f}}$")
    ax_none.set_ylabel("no break", fontsize=9)
    ax_none.grid(True, axis="x", alpha=0.3)

    if save_path is not None:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
    return fig, (ax, ax_none)


fig, axes = plot_resonance_breaking(df, save_path="Plots_critical_mass/resonance_breaking_plot.png")