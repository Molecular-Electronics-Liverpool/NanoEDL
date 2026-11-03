"""02 - detrending lineare della corrente mediana.

Sottrae la retta passante per (T1_MS, I) e (T2_MS, I), poi riancora il segnale
al valore originale in T2 e converte in nA. Identico a 02_detrend.py originale.
Scrive 02_detrended/detrended_<solv>.txt con t(ms), I(nA), V(V).
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *

os.makedirs(DETREND, exist_ok=True)
for s in ORDER:
    f = os.path.join(MEDIAN, f"median_{s}.txt")
    if not os.path.exists(f): continue
    d = np.loadtxt(f, skiprows=1)
    t, V, I = d[:, 0], d[:, 1], d[:, 4]
    i1 = int(np.argmin(np.abs(t - T1_MS))); i2 = int(np.argmin(np.abs(t - T2_MS)))
    base = I[i1] + (I[i2] - I[i1]) * (t - t[i1]) / (t[i2] - t[i1])
    I_nA = (I - base + I[i2]) * 1e9
    np.savetxt(os.path.join(DETREND, f"detrended_{s}.txt"),
               np.column_stack((t, I_nA, V)),
               header="t(ms)\tI_detrended(nA)\tV(V)", fmt="%.8e", delimiter="\t")
    print(f"[ok] {s}", flush=True)
print("02 fatto.")
