"""01 - traccia mediana per ogni solvente.

Legge 00_raw/newfile_<solv>.txt (colonne: indice, V, log10(G/G0)), taglia in
tracce da NPTS campioni, prende la mediana punto per punto su tutte le tracce,
e ricostruisce la corrente I = G0 * 10^log(G/G0) * V.
Scrive 01_median/median_<solv>.txt con colonne t(ms), V(V), log(G/G0), G(S), I(A).
"""
import os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *

os.makedirs(MEDIAN, exist_ok=True)
counts = {}
for s in ORDER:
    p = os.path.join(RAW, f"newfile_{s}.txt")
    if not os.path.exists(p):
        print(f"[skip] {s}: manca {p}", flush=True); continue
    t0 = time.time()
    Vs, Ls = [], []
    for ch in pd.read_csv(p, sep="\t", header=None, usecols=[1, 2],
                          dtype=np.float32, chunksize=4_000_000, engine="c"):
        Vs.append(ch[1].to_numpy()); Ls.append(ch[2].to_numpy())
    V = np.concatenate(Vs); L = np.concatenate(Ls); del Vs, Ls
    n = len(L) // NPTS
    V = V[:n*NPTS].reshape(n, NPTS); L = L[:n*NPTS].reshape(n, NPTS)
    Vm = np.median(V, axis=0).astype(np.float64)
    Lm = np.median(L, axis=0).astype(np.float64)
    G  = G0 * 10.0**Lm
    I  = G * Vm
    t  = np.arange(NPTS) * DT_MS
    np.savetxt(os.path.join(MEDIAN, f"median_{s}.txt"),
               np.column_stack((t, Vm, Lm, G, I)),
               header="t(ms)\tV(V)\tlog(G/G0)\tG(S)\tI(A)", fmt="%.8e", delimiter="\t")
    counts[s] = n
    print(f"[ok] {s}: {n} tracce, {time.time()-t0:.1f}s", flush=True)

pd.DataFrame([{"solvent": k, "n_traces": v} for k, v in counts.items()]) \
  .to_csv(os.path.join(MEDIAN, "n_traces.csv"), index=False)
print("01 fatto.")
