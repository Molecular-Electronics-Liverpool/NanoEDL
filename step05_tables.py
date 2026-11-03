"""05 - aggregazione per solvente e tabelle finali."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *

os.makedirs(TABLES, exist_ok=True)
df = pd.read_csv(os.path.join(FITS, "fit_per_window.csv"))
df = df[df.used & df.tau_ms.notna()]
ntr = pd.read_csv(os.path.join(MEDIAN, "n_traces.csv")).set_index("solvent").n_traces.to_dict()

def build(sel, tag):
    rows = []
    for s in ORDER:
        sub = sel(df[df.solvent == s])
        if not len(sub): continue
        rows.append(dict(solvent=s, group=GROUP[s], filepath=f"00_raw/newfile_{s}.txt",
                         Rm_MOhm=sub.Rm_MOhm.mean(), Rs_MOhm=sub.Rs_MOhm.mean(),
                         Ceq_nF=sub.Ceq_pF.mean()/1000.0, tau_ms=sub.tau_ms.mean(),
                         Rm_sem=sub.Rm_MOhm.sem(), Rs_sem=sub.Rs_MOhm.sem(),
                         Ceq_sem_nF=sub.Ceq_pF.sem()/1000.0, tau_sem=sub.tau_ms.sem(),
                         n_windows_used=len(sub), n_traces=ntr.get(s),
                         **PROPS[s], chi_mean=chi_mean(s)))
    cols = ["solvent","filepath","Rm_MOhm","Rs_MOhm","Ceq_nF","tau_ms",
            "n_windows_used","eps_r","eta_mPa_s","mu_D","chi_mean"]      # come gli originali
    out = pd.DataFrame(rows)
    out[cols].to_csv(os.path.join(TABLES, f"solvent_table_{tag}.csv"), index=False)
    return out

full = build(lambda d: d, "total")
build(lambda d: d[d.edge == "rising"],  "rising")
build(lambda d: d[d.edge == "falling"], "falling")
full.to_csv(os.path.join(TABLES, "solvent_table_total_full.csv"), index=False)

with open(os.path.join(TABLES, "SUMMARY.txt"), "w") as f:
    f.write("TABELLA FINALE - media +/- SEM sulle finestre ritenute valide\n")
    f.write("(1 e 10 sono i segmenti a bias DC; esclusioni per solvente in config.EXTRA_SKIP)\n")
    f.write(f"VOLTAGE_MODE = {VOLTAGE_MODE}\n")
    f.write("Rs = V_ref/|A|, Rm = |V plateau|/|B|, Ceq = tau/Rs, con le tensioni\n")
    if VOLTAGE_MODE == "step":
        f.write("misurate finestra per finestra. V_ref = gradino al fronte: 2.00 V nei\n")
        f.write("semiperiodi interni, circa 1.2 V nel primo, che parte dal bias DC.\n\n")
    else:
        f.write("misurate finestra per finestra. V_ref = plateau del semiperiodo, ~1 V.\n\n")
    f.write(f"{'solvente':11s} {'gr':>2} {'n_tracce':>9} {'n_fin':>5} "
            f"{'Rs(MOhm)':>16} {'tau(ms)':>15} {'Ceq(pF)':>15} {'Rm(MOhm)':>15}\n")
    for _, r in full.iterrows():
        f.write(f"{r.solvent:11s} {r.group:>2} {int(r.n_traces):>9} {int(r.n_windows_used):>5} "
                f"{r.Rs_MOhm:>9.1f}+-{r.Rs_sem:<5.1f} {r.tau_ms:>9.3f}+-{r.tau_sem:<4.3f} "
                f"{r.Ceq_nF*1000:>9.1f}+-{r.Ceq_sem_nF*1000:<4.1f} "
                f"{r.Rm_MOhm:>9.2f}+-{r.Rm_sem:<4.2f}\n")
    f.write("\nRipartizione dei fronti (una risposta lineare ne da' meta' e meta'):\n")
    d2 = pd.read_csv(os.path.join(FITS, "fit_per_window.csv"))
    d2 = d2[d2.used & d2.tau_ms.notna()]
    for s in ORDER:
        sub = d2[d2.solvent == s]
        if not len(sub): continue
        r = int((sub.edge == "rising").sum()); fa = int((sub.edge == "falling").sum())
        flag = "  <-- tutte dello stesso segno" if min(r, fa) == 0 else ""
        f.write(f"  {s:11s} gruppo {GROUP[s]}   salita {r}   discesa {fa}{flag}\n")
print(open(os.path.join(TABLES, "SUMMARY.txt")).read())
print("05 fatto.")
