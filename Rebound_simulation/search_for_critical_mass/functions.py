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
        Standard 177° statt genau 180°, um numerisches Rauschen an der
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






def init_resonance_file(path):
    """Öffnet die Resonanz-Ausgabedatei neu und schreibt den Header."""
    f = open(path, 'w')
    f.write("sma,psi1,psi1_break_year,psi2,psi2_break_year,psi3,psi3_break_year\n")
    f.flush()
    return f


def init_instability_file(path):
    """Öffnet die Instabilitäts-Ausgabedatei neu und schreibt den Header."""
    f = open(path, 'w')
    f.write("reason,year,sma,mass\n")
    f.flush()
    return f


def _fmt(x, decimals=5):
    """Formatiert einen Float als String mit fester Nachkommastellenzahl,
    oder gibt 'None' zurück, falls x None ist."""
    return "None" if x is None else f"{x:.{decimals}f}"


def write_resonance_row(f, sma, psi1_mass, psi1_year, psi2_mass, psi2_year, psi3_mass, psi3_year):
    f.write(
        f"{_fmt(sma, 3)},"
        f"{_fmt(psi1_mass)},{_fmt(psi1_year, 2)},"
        f"{_fmt(psi2_mass)},{_fmt(psi2_year, 2)},"
        f"{_fmt(psi3_mass)},{_fmt(psi3_year, 2)}\n"
    )
    f.flush()


def write_instability_row(f, reason, year, sma, mass):
    f.write(f"{reason},{year:.4f},{_fmt(sma, 3)},{_fmt(mass)}\n")
    f.flush()
