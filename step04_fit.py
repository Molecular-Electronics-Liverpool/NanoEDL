"""04 - fit di ogni finestra e parametri del circuito equivalente.

Modello, con le lettere del SI:
    i(t) = A * exp(-t/tau) + B        A = transitorio EDL, B = plateau molecolare

tau e B vengono dalla return map I(k+1) = n*I(k) + c stimata con Theil-Sen:
    tau = -dt/ln(n)        B = c/(1-n)
A e' poi il minimo quadrato lineare con tau e B fissi.

Conversione in resistenze, con le tensioni MISURATE finestra per finestra:
    Rs  = V_ref      / |A|      V_ref secondo config.VOLTAGE_MODE
    Rm  = |V_plateau|/ |B|
    Ceq = tau / Rs
Con VOLTAGE_MODE = "plateau" entrambe usano la tensione applicata al giunto nel
semiperiodo, circa 1 V, misurata finestra per finestra. Con "step" Rs usa invece
il gradino al fronte, 2 V nei semiperiodi interni e circa 1.2 V nel primo.
Le due scelte differiscono di un fattore due su Rs e, in senso opposto, su Ceq.
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *

def theil_sen(x, y):
    sl = []
    for i in range(len(x) - 1):
        dx = x[i+1:] - x[i]; dy = y[i+1:] - y[i]
        m = dx != 0
        if np.any(m): sl.append(dy[m] / dx[m])
    sl = np.concatenate(sl)
    a = float(np.median(sl))
    return a, float(np.median(y - a * x))

os.makedirs(FITS, exist_ok=True)
meta = pd.read_csv(os.path.join(WINDOWS, "windows_meta.csv"))
rows = []
for s in ORDER:
    dirw = os.path.join(WINDOWS, s)
    if not os.path.isdir(dirw): continue
    for fn in sorted(os.listdir(dirw)):
        if not fn.startswith("window_"): continue
        w = int(fn.split("_")[1].split(".")[0])
        d = np.loadtxt(os.path.join(dirw, fn), skiprows=1)
        t, I = d[:, 0] - d[0, 0], d[:, 1]
        n, c = theil_sen(I[:-1], I[1:])
        rec = dict(solvent=s, group=GROUP[s], window=w, n_slope=n, c_intercept=c)
        m = meta[(meta.solvent == s) & (meta.window == w)].iloc[0]
        rec.update(V_plateau=float(m.V_plateau), dV_edge=float(m.dV_edge),
                   used=bool(m.used))
        if not (0.0 < n < 1.0):
            rec.update(tau_ms=np.nan, A_nA=np.nan, B_nA=np.nan, edge="",
                       Rs_MOhm=np.nan, Rm_MOhm=np.nan, Ceq_pF=np.nan, note="n fuori (0,1)")
            rows.append(rec); continue
        tau = -(t[1] - t[0]) / np.log(n)
        B   = c / (1.0 - n)
        e   = np.exp(-t / tau)
        A   = float(np.sum((I - B) * e) / np.sum(e * e))
        V_ref = abs(rec["V_plateau"]) if VOLTAGE_MODE == "plateau" else abs(rec["dV_edge"])
        Rs  = V_ref                 / abs(A * 1e-9) / 1e6
        Rm  = abs(rec["V_plateau"]) / abs(B * 1e-9) / 1e6
        rec.update(V_ref=V_ref, tau_ms=tau, A_nA=A, B_nA=B,
                   edge=("rising" if A > 0 else "falling"),
                   Rs_MOhm=Rs, Rm_MOhm=Rm,
                   Ceq_pF=(tau * 1e-3) / (Rs * 1e6) * 1e12, note="")
        rows.append(rec)
    print(f"[ok] {s}", flush=True)

df = pd.DataFrame(rows)
cols = ["solvent","group","window","used","edge","V_plateau","dV_edge","V_ref",
        "n_slope","c_intercept","tau_ms","A_nA","B_nA",
        "Rs_MOhm","Rm_MOhm","Ceq_pF","note"]
df[cols].to_csv(os.path.join(FITS, "fit_per_window.csv"), index=False)
print("04 fatto.")
