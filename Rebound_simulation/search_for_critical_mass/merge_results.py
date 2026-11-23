from pathlib import Path
import sys
import functions

BASE = Path("/pfs/10/project/bw16g002/Lisa/Masterarbeit/critical_mass_search")
N_a = 200
N_MASSES = len(functions.build_mass_array(points_per_decade=100))


def read_csv_lines(path):
    """Gibt (header, datenzeilen) zurück."""
    lines = path.read_text().splitlines()
    return lines[0], lines[1:]


def merge(subdir, prefix, out_name, skip_a=frozenset()):
    """Fügt die Einzeldateien zusammen. a-Indizes in skip_a werden ausgelassen."""
    header, rows, missing = None, [], []
    for a in range(N_a):
        if a in skip_a:
            continue
        p = BASE / subdir / f"{prefix}_a={a:03d}.csv"
        if not p.exists():
            missing.append(a)
            continue
        file_header, data = read_csv_lines(p)
        if header is None:
            header = file_header
        elif file_header != header:
            raise RuntimeError(f"{p}: Header weicht ab")
        rows += data
    if header is not None:
        (BASE / out_name).write_text(
            header + "\n" + "\n".join(rows) + ("\n" if rows else "")
        )
    return len(rows), missing


def find_fully_unstable():
    """a-Indizes, bei denen ALLE Massen instabil waren."""
    fully = set()
    for a in range(N_a):
        p = BASE / "instability_cases" / f"instability_cases_a={a:03d}.csv"
        if not p.exists():
            continue
        _, data = read_csv_lines(p)
        if len(data) == N_MASSES:
            fully.add(a)
    return fully


fully_unstable = find_fully_unstable()

# Resonanz-Datei ohne komplett instabile a's, Instabilitäts-Datei vollständig
n_res, miss_res = merge("resonance_breaking", "resonance_breaking",
                        "resonance_breaking_ALL.csv", skip_a=fully_unstable)
n_ins, miss_ins = merge("instability_cases", "instability_cases",
                        "instability_cases_ALL.csv")

print(f"Resonanz:     {n_res} Zeilen, fehlende a_index: {miss_res}")
print(f"Instabilität: {n_ins} Zeilen, fehlende a_index: {miss_ins}")
print(f"Komplett instabile a_index (nicht in Resonanz-CSV): {sorted(fully_unstable)}")