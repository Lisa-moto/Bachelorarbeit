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


def is_laplace_resonant(phi_deg, threshold_deg=177.0, return_diagnostics=False):
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



def libration_score(psi_deg):
    """
    Gibt zurück: (is_librating: bool, circ_std_deg: float)
    psi_deg: 1D array, Winkel in Grad über die Zeit
    """
    psi_rad = np.deg2rad(psi_deg)
    R = np.sqrt(np.mean(np.sin(psi_rad))**2 + np.mean(np.cos(psi_rad))**2)
    circ_std = np.inf if R <= 1e-12 else np.rad2deg(np.sqrt(-2*np.log(R)))
    is_librating = circ_std < 90  # Schwellwert, siehe unten
    return is_librating, circ_std


def init_results_file(path='results_moon_scan.csv'):
    file_exists = os.path.isfile(path)
    f = open(path, 'a', newline='')
    writer = csv.writer(f)
    if not file_exists:
        writer.writerow(['a', 'm', 'psi1_librating', 'psi1_std',
                          'psi2_librating', 'psi2_std',
                          'psi3_librating', 'psi3_std'])
    return f, writer


def append_result(writer, f, a, m, results):
    # results: dict {'psi1': (bool, std), 'psi2': (...), 'psi3': (...)}
    writer.writerow([
        a, m,
        results['psi1'][0], results['psi1'][1],
        results['psi2'][0], results['psi2'][1],
        results['psi3'][0], results['psi3'][1],
    ])
    f.flush()
    os.fsync(f.fileno())
