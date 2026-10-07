import numpy as np
import csv, os

# Convert from radians to degrees in [0,360)
def rad_to_deg_0_360(arr):
  a = np.mod(arr, 2*np.pi)  # now in [0, 2pi)
  return a * 180.0 / np.pi

def laplace_angles(longitude):
  # return time series arrays (length n_steps) for each resonant angle

  psi1 = longitude[:, 1] - 4 * longitude[:, 2] + 3 * longitude[:, 3]
  psi2 = 2 * longitude[:, 2] - 5 * longitude[:, 3] + 3 * longitude[:, 4]
  psi3 = longitude[:, 3] - 3 * longitude[:, 4] + 2 * longitude[:, 5]

  psi1 = rad_to_deg_0_360(psi1)
  psi2 = rad_to_deg_0_360(psi2)
  psi3 = rad_to_deg_0_360(psi3)

  return psi1, psi2, psi3


def is_laplace_resonant(phi_deg, threshold_deg=179.0, return_diagnostics=False):
    """
    Prüft, ob ein Laplace-Winkel (Zeitreihe, in Grad) noch resoniert,
    d.h. librativ statt zirkulierend ist.

    Parameters
    ----------
    phi_deg : array_like
        Zeitreihe des Laplace-Winkels in Grad (beliebig gewrappt,
        z.B. [0,360) oder (-180,180]).
    threshold_deg : float
        Amplituden-Grenze in Grad, unterhalb derer Libration angenommen wird.
        Standard 179.0° statt genau 180°, um numerisches Rauschen an der
        Separatrix nicht als "gerade noch resonant" fehlzuinterpretieren.
    return_diagnostics : bool
        Wenn True, zusätzlich ein dict mit Amplitude, zirkulärem Mittelwert
        und Std zurückgeben.

    Returns
    -------
    bool  (oder (bool, dict) falls return_diagnostics=True)
    """
    phi = np.deg2rad(np.asarray(phi_deg, dtype=float))

    # zirkulärer Mittelwert (robust gegen 0/360-Wrap)
    mean_angle = np.arctan2(np.mean(np.sin(phi)), np.mean(np.cos(phi)))

    # Abweichung vom Mittelwert, sauber in (-pi, pi] gewrappt
    delta = (phi - mean_angle + np.pi) % (2 * np.pi) - np.pi

    amplitude_deg = np.degrees(np.max(np.abs(delta)))

    # zirkuläre Std als Zusatzinfo (nicht als Entscheidungskriterium)
    R = np.hypot(np.mean(np.sin(phi)), np.mean(np.cos(phi)))
    circ_std_deg = np.degrees(np.sqrt(-2 * np.log(R))) if R > 0 else np.inf

    is_resonant = amplitude_deg < threshold_deg

    if return_diagnostics:
        return is_resonant, {
            "amplitude_deg": amplitude_deg,
            "mean_angle_deg": np.degrees(mean_angle) % 360,
            "circular_std_deg": circ_std_deg,
        }
    return is_resonant


### für die CSV-Ausgabe der Simulationsergebnisse (irrelevant) ###

# def init_resonance_file(path):
#     """Öffnet die Resonanz-Ausgabedatei neu und schreibt den Header."""
#     f = open(path, 'w')
#     f.write("sma,psi1_break_mass,psi1_break_year,psi2_break_mass,psi2_break_year,psi3_break_mass,psi3_break_year\n")
#     f.flush()
#     return f


# def init_instability_file(path):
#     """Öffnet die Instabilitäts-Ausgabedatei neu und schreibt den Header."""
#     f = open(path, 'w')
#     f.write("reason,year,sma,mass\n")
#     f.flush()
#     return f
    
# def write_resonance_row(f, sma, psi1_mass, psi1_year, psi2_mass, psi2_year, psi3_mass, psi3_year):
#     f.write(
#         f"{_fmt(sma, 3)},"
#         f"{_fmt(psi1_mass)},{_fmt(psi1_year, 2)},"
#         f"{_fmt(psi2_mass)},{_fmt(psi2_year, 2)},"
#         f"{_fmt(psi3_mass)},{_fmt(psi3_year, 2)}\n"
#     )
#     f.flush()


# def write_instability_row(f, reason, year, sma, mass):
#     f.write(f"{reason},{year:.4f},{_fmt(sma, 3)},{_fmt(mass)}\n")
#     f.flush()


def _fmt(x, decimals=9):
    """Formatiert einen Float als String mit fester Nachkommastellenzahl,
    oder gibt 'None' zurück, falls x None ist."""
    return "None" if x is None else f"{x:.{decimals}f}"

def _fmt_mass(m):
    return "None" if m is None else f"{m:.6e}"


### für die Massen-Array-Berechnung im Resonanz-Suchlauf ###

def build_mass_array(discrete_low=(1e-9, 1e-8, 1e-7),
                      log_start=1e-6, log_end=0.01,
                      points_per_decade=100):
    """
    Baut das Massen-Array: feste Einzelwerte im untersten Bereich,
    danach logarithmisch gleichverteilt.
    """
    n_decades = np.log10(log_end / log_start)
    n_log_points = int(round(n_decades * points_per_decade)) + 1
    log_part = np.logspace(np.log10(log_start), np.log10(log_end), n_log_points)
    return np.concatenate([np.array(discrete_low), log_part])

### für die Bestimmung der zu speichernden Zeitfenster ###

def compute_save_windows(break_year, end_year, early_stop, window_years=300):
    """
    Bestimmt die zu speichernden Zeitfenster (in Jahren), gemergt und ohne
    Überlappung.

    - Für jeden Winkel, der VOR Simulationsende gebrochen ist (break_year < end_year):
      Fenster [break_year - window_years, break_year].
    - Zusätzlich das Fenster [end_year - window_years, end_year], AUSSER die
      Simulation wurde früh abgebrochen, weil alle drei Winkel bereits
      gebrochen waren (early_stop=True).
    """
    windows = []
    for key in ("psi1", "psi2", "psi3"):
        y = break_year.get(key)
        if y is not None and y <= end_year:
            start = max(0.0, y - window_years)
            windows.append((start, y))

    if not early_stop:
        start = max(0.0, end_year - window_years)
        windows.append((start, end_year))

    if not windows:
        return []

    windows.sort(key=lambda w: w[0])
    merged = [windows[0]]
    for start, stop in windows[1:]:
        last_start, last_stop = merged[-1]
        if start <= last_stop:
            merged[-1] = (last_start, max(last_stop, stop))
        else:
            merged.append((start, stop))

    return merged

# neue Funktionen für Fortschritt-Log zur Ausfallsicherheit

def atomic_write_text(path, text):
    path = str(path)
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w") as f:
        f.write(text)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def append_log(path, fields):
    """Hängt genau eine Zeile an und erzwingt das Schreiben auf die Platte."""
    with open(path, "a") as f:
        f.write(";".join(fields) + "\n")
        f.flush()
        os.fsync(f.fileno())


def load_log(path, masses):
    """Liest das Log. Unvollständige/kaputte Zeilen am Ende (Abbruch beim
    Schreiben) werden verworfen. Gibt eine Liste von Feldlisten zurück,
    Eintrag k gehört zu masses[k]."""
    entries = []
    if os.path.exists(path):
        with open(path) as f:
            raw = f.read()
        # letztes Element von split ist entweder "" oder eine halbe Zeile -> weg
        for k, line in enumerate(raw.split("\n")[:-1]):
            p = line.split(";")
            if len(p) != 6 or p[0] != str(k):
                break
            if k >= len(masses) or not np.isclose(float(p[1]), masses[k], rtol=1e-9, atol=0):
                raise RuntimeError(f"{path}: Zeile {k} passt nicht zum aktuellen "
                                   f"Massen-Array. Bitte prüfen/löschen.")
            entries.append(p)
        # Datei bereinigt neu schreiben, damit das nächste Anhängen sauber ist
        atomic_write_text(path, "".join(";".join(p) + "\n" for p in entries))
    return entries


def first_breaks(entries):
    """Pro psi die erste Masse mit gebrochener Resonanz: [(mass, year) oder None]*3."""
    result = [None, None, None]
    for p in entries:
        if p[2] != "ok":
            continue
        for n in range(3):
            if result[n] is None and p[3 + n] != "None":
                result[n] = (float(p[1]), float(p[3 + n]))
    return result


def write_outputs_from_log(entries, sma, resonance_path, instability_path):
    """Erzeugt beide CSVs komplett aus dem Log (keine Duplikate möglich)."""
    b = first_breaks(entries)
    cols = []
    for x in b:
        cols += [_fmt_mass(None if x is None else x[0]), _fmt(None if x is None else x[1], 2)]
    atomic_write_text(resonance_path,
        "sma,psi1_break_mass,psi1_break_year,psi2_break_mass,psi2_break_year,"
        "psi3_break_mass,psi3_break_year\n" + f"{_fmt(sma, 6)}," + ",".join(cols) + "\n")

    lines = ["reason,year,sma,mass\n"]
    for p in entries:
        if p[2] == "unstable":
            reason = p[3].replace(",", ";")
            lines.append(f"{reason},{float(p[4]):.4f},{_fmt(sma, 6)},{_fmt_mass(float(p[1]))}\n")
    atomic_write_text(instability_path, "".join(lines))