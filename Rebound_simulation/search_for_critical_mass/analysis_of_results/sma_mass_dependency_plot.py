import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import sys
from matplotlib.ticker import LogLocator, LogFormatterSciNotation, NullFormatter

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
df_ins = pd.read_csv(BASE_DIR / "instability_cases_ALL.csv")

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
            alpha=0.8, label=labels[angle]
        )

        # Zusatzplot: Striche dort, wo kein Bruch auftritt
        ax_none.scatter(
            sma[~has_break], np.full((~has_break).sum(), i),
            s=40, marker="|", color=colors[angle], linewidths=1.5,
        )

    # Hauptplot formatieren
    if log_y:
        ax.set_yscale("log")

        # Achsenbereich aus den Daten, mit etwas Rand
        all_mass = np.concatenate([
            df[f"{a}_break_mass"].values.astype(float) for a in ANGLES
        ])
        all_mass = all_mass[np.isfinite(all_mass) & (all_mass > 0)]
        if all_mass.size:
            ax.set_ylim(all_mass.min() / 1.5, all_mass.max() * 1.5)

        # Major-Ticks: jede Dekade beschriftet
        ax.yaxis.set_major_locator(LogLocator(base=10, numticks=20))
        ax.yaxis.set_major_formatter(LogFormatterSciNotation(base=10))

        # Minor-Ticks: 2,3,...,9 mal 10^k, mit Label
        ax.yaxis.set_minor_locator(
            LogLocator(base=10, subs=(2, 3, 4, 5, 6, 7), numticks=50)
        )
        ax.yaxis.set_minor_formatter(
            LogFormatterSciNotation(base=10, labelOnlyBase=False,
                                    minor_thresholds=(10, 10))
        )
        ax.tick_params(axis="y", which="minor", labelsize=7)
    ax.set_ylabel(r"critical mass in $m_{\text{f}}$")
    ax.set_title("Resonance Breaking: Dependency on Semi Major Axis and Mass of the Moon")
    ax.grid(True, alpha=0.3)
    ax.legend(title="Resonance Angles", loc="best")
    ax.tick_params(labelbottom=False)

    # Zusatzplot formatieren
    ax_none.set_yticks(range(len(ANGLES)))
    ax_none.set_yticklabels([labels[a] for a in ANGLES])
    ax_none.set_ylim(-0.7, len(ANGLES) - 0.3)
    ax_none.set_xlabel(r"semi major axis in $R_{\text{Hill,f}}$")
    ax_none.set_ylabel("no break", fontsize=9)
    ax_none.grid(True, axis="x", alpha=0.3)

    if save_path is not None:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
    return fig, (ax, ax_none)


fig, axes = plot_resonance_breaking(df, save_path="Plots_critical_mass/resonance_breaking_plot.png")


def plot_instability_cases(df, save_path=None, log_y=True):
    """
    Scatterplot der instabilen Läufe: große Halbachse (sma) gegen Masse,
    mit unterschiedlichen Markern je nach Instabilitätsgrund.
    """
    reason = df["reason"].astype(str).str.lower()

    # Gruppen anhand des Grund-Strings (enthält "collision" bzw. "ejection")
    groups = {
        "collision": reason.str.contains("collision"),
        "ejection":  reason.str.contains("ejection"),
    }
    # Alles, was keiner der beiden Gruppen zugeordnet werden kann
    groups["other"] = ~(groups["collision"] | groups["ejection"])

    style = {
        "collision": dict(marker="o", color="tab:red",    label="Collision"),
        "ejection":  dict(marker="^", color="tab:purple", label="Ejection"),
        "other":     dict(marker="x", color="gray",       label="Other"),
    }

    fig, ax = plt.subplots(figsize=(10, 6))

    for name, mask in groups.items():
        if not mask.any():
            continue
        ax.scatter(
            df.loc[mask, "sma"], df.loc[mask, "mass"],
            s=18, alpha=0.7, **style[name],
        )

    if log_y:
        ax.set_yscale("log")
        mass = df["mass"].values.astype(float)
        mass = mass[np.isfinite(mass) & (mass > 0)]
        if mass.size:
            ax.set_ylim(mass.min() / 1.5, mass.max() * 1.5)
        ax.yaxis.set_major_locator(LogLocator(base=10, numticks=20))
        ax.yaxis.set_major_formatter(LogFormatterSciNotation(base=10))

    ax.set_xlabel(r"semi major axis in $R_{\text{Hill,f}}$")
    ax.set_ylabel(r"mass in $m_{\text{f}}$")
    ax.set_title("Unstable Cases: Dependency on Semi Major Axis and Mass of the Moon")
    ax.grid(True, alpha=0.3)
    ax.legend(title="Reason", loc="best")

    if save_path is not None:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
    return fig, ax


fig, ax = plot_instability_cases(
    df_ins, save_path="Plots_critical_mass/instability_cases_plot.png"
)