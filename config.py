"""Costanti, percorsi e proprieta' dei solventi."""
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW      = os.path.join(BASE, "00_raw")
MEDIAN   = os.path.join(BASE, "01_median")
DETREND  = os.path.join(BASE, "02_detrended")
WINDOWS  = os.path.join(BASE, "03_windows")
FITS     = os.path.join(BASE, "04_fits")
TABLES   = os.path.join(BASE, "05_tables")
PLOTS    = os.path.join(BASE, "06_plots")

# --- acquisizione -----------------------------------------------------
NPTS    = 4000          # campioni per traccia
DT_MS   = 0.025         # passo di campionamento (40 kHz)
G0      = 7.748091729e-5
T1_MS   = 8.0           # ancoraggi del detrending lineare
T2_MS   = 92.0
WIN_MS  = 10.0          # semiperiodo a 50 Hz

# --- convenzione sulla tensione ---------------------------------------
# "step":    Rs = |dV al fronte| / |A|  (2 V nei semiperiodi interni)
# "plateau": Rs = |V_plateau| / |A|     (~1 V)
# Rm usa sempre |V_plateau|.
VOLTAGE_MODE = "step"

# Finestre escluse: 1 e 10 sono i segmenti a bias DC.
SKIP_WINDOWS = {1, 10}
# Anisolo: le finestre 2 e 3 precedono la stabilizzazione dell'onda quadra.
EXTRA_SKIP = {"anisole": {2, 3}}

# --- solventi ---------------------------------------------------------
ORDER = ["mesitylene", "decane", "anisole", "tcb", "dcb", "gvl", "gbl", "pc"]
GROUP = {"mesitylene": "A", "decane": "A",
         "anisole": "B", "tcb": "B", "dcb": "B",
         "gvl": "C", "gbl": "C", "pc": "C"}

# Proprieta' a 25 C. Fonti in 05_tables/proprieta.bib.
PROPS = {
 "mesitylene": dict(eps_r=2.279,  eta_mPa_s=0.661, mu_D=0.00,   rho_g_mL=0.8615),
 "decane":     dict(eps_r=1.9853, eta_mPa_s=0.838, mu_D=0.00,   rho_g_mL=0.7266),
 "anisole":    dict(eps_r=4.30,   eta_mPa_s=1.056, mu_D=1.2623, rho_g_mL=0.9940),
 "tcb":        dict(eps_r=2.24,   eta_mPa_s=1.926, mu_D=1.26,   rho_g_mL=1.4590),
 "dcb":        dict(eps_r=10.12,  eta_mPa_s=1.324, mu_D=2.50,   rho_g_mL=1.3059),
 "gvl":        dict(eps_r=36.9,   eta_mPa_s=1.820, mu_D=4.71,   rho_g_mL=1.0551),
 "gbl":        dict(eps_r=39.0,   eta_mPa_s=1.720, mu_D=4.27,   rho_g_mL=1.1296),
 "pc":         dict(eps_r=66.14,  eta_mPa_s=2.500, mu_D=4.90,   rho_g_mL=1.2047),
}

# --- elettronegativita' media (Pauling) --------------------------------
PAULING = {"H": 2.20, "C": 2.55, "O": 3.44, "Cl": 3.16}
FORMULA = {
 "mesitylene": dict(C=9,  H=12),
 "decane":     dict(C=10, H=22),
 "anisole":    dict(C=7,  H=8,  O=1),
 "tcb":        dict(C=6,  H=3,  Cl=3),
 "dcb":        dict(C=6,  H=4,  Cl=2),
 "gvl":        dict(C=5,  H=8,  O=2),
 "gbl":        dict(C=4,  H=6,  O=2),
 "pc":         dict(C=4,  H=6,  O=3),
}

def chi_mean(solvent):
    """Media aritmetica dell'elettronegativita' di Pauling su tutti gli atomi."""
    f = FORMULA[solvent]
    return sum(PAULING[e] * c for e, c in f.items()) / sum(f.values())
