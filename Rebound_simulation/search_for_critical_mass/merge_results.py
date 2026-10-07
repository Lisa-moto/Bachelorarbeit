from pathlib import Path

BASE = Path("/pfs/10/project/bw16g002/Lisa/Masterarbeit/critical_mass_search")
N_a = 200


def merge(subdir, prefix, out_name):
    header, rows, missing = None, [], []
    for a in range(N_a):
        p = BASE / subdir / f"{prefix}_a={a:03d}.csv"
        if not p.exists():
            missing.append(a)
            continue
        lines = p.read_text().splitlines()
        if header is None:
            header = lines[0]
        elif lines[0] != header:
            raise RuntimeError(f"{p}: Header weicht ab")
        rows += lines[1:]
    if header is not None:
        (BASE / out_name).write_text(header + "\n" + "\n".join(rows) + ("\n" if rows else ""))
    return len(rows), missing


n_res, miss_res = merge("resonance_breaking", "resonance_breaking", "resonance_breaking_ALL.csv")
n_ins, miss_ins = merge("instability_cases", "instability_cases", "instability_cases_ALL.csv")

print(f"Resonanz:     {n_res} Zeilen, fehlende a_index: {miss_res}")
print(f"Instabilität: {n_ins} Zeilen, fehlende a_index: {miss_ins}")