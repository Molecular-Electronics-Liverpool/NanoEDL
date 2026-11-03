"""03 - taglio in semiperiodi e misura della tensione di ogni finestra.

Ogni finestra e' un semiperiodo da WIN_MS. I nomi sono a due cifre
(window_01 ... window_10) cosi' l'ordinamento alfabetico coincide con quello
temporale: negli script originali sorted(glob(...)) metteva window_10 fra
window_1 e window_2, scambiando le etichette.

Per ogni finestra misura sui dati due tensioni, che servono al passo 04:
  V_plateau  mediana di V nella seconda meta' della finestra
  dV_edge    salto di tensione al fronte, V_plateau(w) - V_plateau(w-1)
E' questo che distingue il primo semiperiodo, dove l'escursione parte dal bias
DC e vale circa 1.2 V, dai semiperiodi interni, dove vale 2 V.
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *

os.makedirs(WINDOWS, exist_ok=True)
meta = []
for s in ORDER:
    f = os.path.join(DETREND, f"detrended_{s}.txt")
    if not os.path.exists(f): continue
    d = np.loadtxt(f, skiprows=1)
    t, I, V = d[:, 0], d[:, 1], d[:, 2]
    outdir = os.path.join(WINDOWS, s); os.makedirs(outdir, exist_ok=True)
    npw = int(round(WIN_MS / (t[1] - t[0])))
    nw  = len(I) // npw
    prev = 0.0
    for w in range(1, nw + 1):
        sl = slice((w-1)*npw, w*npw)
        np.savetxt(os.path.join(outdir, f"window_{w:02d}.txt"),
                   np.column_stack((t[sl], I[sl])),
                   header="t(ms)\tI(nA)", fmt="%.8e", delimiter="\t")
        vw   = V[sl]
        plat = float(np.median(vw[len(vw)//2:]))
        skip = SKIP_WINDOWS | EXTRA_SKIP.get(s, set())
        meta.append(dict(solvent=s, window=w, V_plateau=plat,
                         dV_edge=plat - prev, used=(w not in skip)))
        prev = plat
    print(f"[ok] {s}: {nw} finestre", flush=True)

pd.DataFrame(meta).to_csv(os.path.join(WINDOWS, "windows_meta.csv"), index=False)
print("03 fatto.")
