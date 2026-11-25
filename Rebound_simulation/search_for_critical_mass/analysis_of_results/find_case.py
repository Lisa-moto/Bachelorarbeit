import numpy as np
import pandas as pd
from pathlib import Path
import sys

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))
import functions

BASE_DIR = Path("/pfs/10/project/bw16g002/Lisa/Masterarbeit/critical_mass_search")
N_a = 200
sma_values = np.linspace(0.001, 0.5, N_a)
N_MASSES = len(functions.build_mass_array(points_per_decade=100))
ANGLES = ("psi1", "psi2", "psi3")

# --- Daten einlesen ----------------------------------------------------------
df = pd.read_csv(BASE_DIR / "resonance_breaking_ALL.csv", na_values="None")
df_ins = pd.read_csv(BASE_DIR / "instability_cases_ALL.csv")

# a_index aus sma rekonstruieren (nächster Gitterpunkt)
df["a_index"] = [int(np.argmin(np.abs(sma_values - s))) for s in df["sma"]]
# Kontrolle: passt das Gitter wirklich?
max_dev = np.abs(sma_values[df["a_index"]] - df["sma"]).max()
print(f"Max. Abweichung sma zu Gitter: {max_dev:.2e} (sollte ~1e-6 oder kleiner sein)")

# Anzahl instabiler Läufe pro a
df_ins["a_index"] = [int(np.argmin(np.abs(sma_values - s))) for s in df_ins["sma"]]
n_unstable = df_ins.groupby("a_index").size()
df["n_unstable"] = df["a_index"].map(n_unstable).fillna(0).astype(int)
df["n_stable_run"] = N_MASSES - df["n_unstable"]

# --- Kein Bruch bei ALLEN drei Winkeln --------------------------------------
no_break_all = df[ANGLES[0] + "_break_mass"].isna()
for a in ANGLES[1:]:
    no_break_all &= df[f"{a}_break_mass"].isna()

sel = df[no_break_all].sort_values("a_index")
print(f"\nAnzahl a-Werte ohne Bruch bei allen drei Winkeln: {len(sel)}")
print(sel[["a_index", "sma", "n_unstable", "n_stable_run"]].to_string(index=False))

# --- Sind es Nachbarn? -------------------------------------------------------
if len(sel) > 1:
    idx = sel["a_index"].values
    gaps = np.diff(idx)
    print(f"\nAbstände der a_index: {gaps}")
    # zusammenhängende Gruppen (Differenz == 1)
    groups = np.split(idx, np.where(gaps != 1)[0] + 1)
    print(f"Zusammenhängende Gruppen: {[list(g) for g in groups]}")
    print(f"Gitterabstand in sma: {sma_values[1] - sma_values[0]:.6f}")
elif len(sel) == 1:
    print("\nGenau ein a-Wert, kein Nachbar betroffen.")
else:
    print("\nKein a-Wert ohne Bruch bei allen drei Winkeln.")

# --- Teilweise ohne Bruch: pro Winkel einzeln --------------------------------
print("\nAnzahl a-Werte ohne Bruch pro Winkel:")
for a in ANGLES:
    print(f"  {a}: {df[f'{a}_break_mass'].isna().sum()}")

# --- Bereits komplett instabile a's (falls noch in der Datei) ----------------
fully = df[df["n_unstable"] == N_MASSES]["a_index"].tolist()
print(f"\nKomplett instabile a_index in der Resonanz-CSV: {fully}")