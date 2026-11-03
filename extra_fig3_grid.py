"""EXTRA - Figura 3, pannelli c-f in un'unica griglia 2x2 larga 8 cm.

   (c) alto-sinistra : R_m  vs densita'
   (d) alto-destra   : R_s  vs viscosita'
   (e) basso-sinistra: C_eq (asse sx) e tau (asse dx) vs momento di dipolo
   (f) basso-destra  : C_eq (asse sx) e tau (asse dx) vs elettronegativita'
                       media di Pauling

   In basso C_eq e tau sono etichettati una sola volta: C_eq sul bordo
   sinistro di (e), tau sul bordo destro di (f). Cambia solo l'ascissa.

   Stile, colori e fit Theil-Sen identici a step06_plots.py; assi lineari.
   Output in 06_plots/: fig3_grid_cdef_8cm_clean|with_text .png e .svg
"""
import os, sys, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import TABLES, PLOTS
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, MultipleLocator

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica"],
    "axes.titlesize": 8, "axes.labelsize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})
DPI_OUT = 800
cm = 1 / 2.54

# ------------------------------------------------------------------ geometria
FIG_W_CM   = 8.0          # larghezza fissata dalle linee guida
LEFT_CM    = 1.35         # ticks + nome dell'asse y di sinistra
RIGHT_CM   = 1.20         # ticks + nome dell'asse y di destra (tau in f)
GAP_X_CM   = 1.15         # spazio fra le due colonne (serve a (d) per il suo asse y)
TOP_CM     = 0.45         # spazio per le lettere (c) e (d)
GAP_Y_CM   = 1.40         # ticks + ascissa di c,d, piu' le lettere (e) e (f)
BOTTOM_CM  = 1.00
AX_W_CM    = (FIG_W_CM - LEFT_CM - RIGHT_CM - GAP_X_CM) / 2.0
AX_H_CM    = 2.30
FIG_H_CM   = TOP_CM + AX_H_CM + GAP_Y_CM + AX_H_CM + BOTTOM_CM

# ------------------------------------------------------------------- stile
def _c(rgb):
    return tuple(np.array(rgb) / 255.0)

RED, BLUE = _c((220, 50, 47)), _c((49, 54, 149))
RED_SOLVENTS  = {"mesitylene", "decane"}
BLUE_SOLVENTS = {"anisole", "tcb", "dcb"}
KEEP = ["mesitylene", "decane", "anisole", "tcb", "dcb"]
CODE = {"mesitylene": "1A", "decane": "2A", "anisole": "1B", "tcb": "2B", "dcb": "3B"}

X_LABELS = {
    "rho_g_mL":  r"$\rho$ (g mL$^{-1}$)",
    "eta_mPa_s": r"$\eta$ (mPa s)",
    "mu_D":      r"$\mu$ (D)",
    "chi_mean":  r"$\overline{\chi}$ (Pauling)",
}
Y_LABELS = {
    "Rm_MOhm": r"$R_m$ (M$\Omega$)",
    "Rs_MOhm": r"$R_s$ (M$\Omega$)",
    "Ceq_nF":  r"$C_{eq}$ (nF)",
    "tau_ms":  r"$\tau$ (ms)",
}
YLIM = {"Rm_MOhm": (0.0, 20.0), "Rs_MOhm": (0.0, 550.0),
        "Ceq_nF": (0.0, 0.06),  "tau_ms": (0.0, 1.1)}
# quantita' per cui il fit Theil-Sen si fa su log10(y), come in step06_plots.py
LOG_FIT = {"Rm_MOhm", "Rs_MOhm", "tau_ms"}

FIT_POINTS, FIT_LW = 200, 0.9

# etichette dei solventi: offset in punti, pannello per pannello
LABEL_OFFSETS = {
    "rho_g_mL":  {"mesitylene": (4.5, 1.5), "decane": (4.5, 1.5),
                  "anisole": (4.5, 1.5), "tcb": (-13.5, 3.5), "dcb": (4.5, 1.5)},
    "eta_mPa_s": {"mesitylene": (4.5, 2.5), "decane": (4.5, -8.5),
                  "anisole": (4.5, 3.5), "tcb": (-13.5, 3.5), "dcb": (4.5, 2.0)},
    "mu_D":      {"mesitylene": (6.0, 0.0), "decane": (6.0, 0.0),
                  "anisole": (6.0, -1.0), "tcb": (6.0, -1.0), "dcb": (-14.0, 5.0)},
    "chi_mean":  {"mesitylene": (2.5, 12.0), "decane": (-2.0, 12.0),
                  "anisole": (5.0, -1.0), "tcb": (5.0, -1.0), "dcb": (-15.0, -9.0)},
}
LABEL_OVERRIDE = {"mu_D": {"mesitylene": "1A, 2A", "decane": None}}

# ------------------------------------------------- Theil-Sen (come step06)
def _theil_sen(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 2 or len(np.unique(x)) < 2:
        return None
    sl = []
    for i in range(len(x) - 1):
        dx = x[i+1:] - x[i]; dy = y[i+1:] - y[i]
        v = dx != 0
        if np.any(v):
            sl.extend((dy[v] / dx[v]).tolist())
    if not sl:
        return None
    a = float(np.median(sl)); b = float(np.median(y - a * x))
    return a, b, float(x.min()), float(x.max())

def _overlay(ax, x, y, colour, ykey, linestyle):
    logfit = ykey in LOG_FIT
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if logfit:
        m &= (y > 0)
    x, y = x[m], y[m]
    fit = _theil_sen(x, np.log10(y) if logfit else y)
    if fit is None:
        return None
    a, b, x0, x1 = fit
    if x1 <= x0:
        return None
    xx = np.linspace(x0, x1, FIT_POINTS)
    ax.plot(xx, 10 ** (a * xx + b) if logfit else a * xx + b,
            linestyle=linestyle, linewidth=FIT_LW, color=colour, zorder=1)
    yp = 10 ** (a * x + b) if logfit else a * x + b
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - float(np.sum((y - yp) ** 2)) / ss_tot if ss_tot > 0 else float("nan")
    return a, b, r2

def _frame(ax):
    ax.tick_params(axis="both", labelsize=6, width=1.0, length=4)
    for sp in ax.spines.values():
        sp.set_linewidth(1.0)

def _annotate(ax, x, y, solv, xkey):
    OFF = LABEL_OFFSETS[xkey]
    OVR = LABEL_OVERRIDE.get(xkey, {})
    for xi, yi, k in zip(x, y, solv):
        lab = OVR.get(k, CODE[k])
        if lab is None:
            continue
        ax.annotate(lab, (xi, yi), textcoords="offset points",
                    xytext=OFF.get(k, (4.5, 2.0)), fontsize=6, color="k")

def _xlim(ax, x, right_factor=2.4):
    pad = 0.10 * (x.max() - x.min()) if x.max() > x.min() else 0.5
    ax.set_xlim(x.min() - pad, x.max() + right_factor * pad)

# ---------------------------------------------------------------- pannelli
def panel_single(ax, rows, xkey, ykey, with_text):
    """Pannelli (c) e (d): una sola quantita' sull'asse y."""
    solv = [r["solvent"].strip().lower() for r in rows]
    x = np.array([float(r[xkey]) for r in rows])
    y = np.array([float(r[ykey]) for r in rows])
    col = [RED if k in RED_SOLVENTS else BLUE for k in solv]

    ax.scatter(x, y, s=22, c=col, edgecolors="none", zorder=3)
    rep = {}
    for name, group, colour in (("A", RED_SOLVENTS, RED), ("B", BLUE_SOLVENTS, BLUE)):
        m = np.array([k in group for k in solv])
        rep[name] = _overlay(ax, x[m], y[m], colour, ykey, "--")

    ax.set_ylim(*YLIM[ykey])
    _xlim(ax, x)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    _frame(ax)
    if with_text:
        ax.set_xlabel(X_LABELS[xkey])
        ax.set_ylabel(Y_LABELS[ykey])
        _annotate(ax, x, y, solv, xkey)
    return rep

def panel_dual(axL, axR, rows, xkey, with_text):
    """Pannelli (e) e (f): C_eq a sinistra, tau a destra."""
    solv = [r["solvent"].strip().lower() for r in rows]
    x   = np.array([float(r[xkey]) for r in rows])
    ceq = np.array([float(r["Ceq_nF"]) for r in rows])
    tau = np.array([float(r["tau_ms"]) for r in rows])
    col = [RED if k in RED_SOLVENTS else BLUE for k in solv]

    axL.scatter(x, ceq, s=22, c=col, marker="o", edgecolors="none", zorder=3)
    axR.scatter(x, tau, s=22, facecolors="none", edgecolors=col, marker="s",
                linewidths=0.9, zorder=3)
    rep = {}
    for name, group, colour in (("A", RED_SOLVENTS, RED), ("B", BLUE_SOLVENTS, BLUE)):
        m = np.array([k in group for k in solv])
        rep[name] = (_overlay(axL, x[m], ceq[m], colour, "Ceq_nF", "--"),
                     _overlay(axR, x[m], tau[m], colour, "tau_ms", ":"))

    axL.set_ylim(*YLIM["Ceq_nF"])
    axR.set_ylim(*YLIM["tau_ms"])
    _xlim(axL, x, right_factor=2.6 if xkey == "chi_mean" else 2.4)
    axL.yaxis.set_major_locator(MaxNLocator(nbins=4))
    axR.yaxis.set_major_locator(MaxNLocator(nbins=4))
    if xkey == "chi_mean":
        axL.xaxis.set_major_locator(MultipleLocator(0.1))
    else:
        axL.xaxis.set_major_locator(MaxNLocator(nbins=4))
    _frame(axL); _frame(axR)
    if with_text:
        axL.set_xlabel(X_LABELS[xkey])
        _annotate(axL, x, ceq, solv, xkey)
    return rep

# --------------------------------------------------------------------- IO
def read_rows():
    with open(os.path.join(TABLES, "solvent_table_total_full.csv"), newline="") as f:
        rows = [r for r in csv.DictReader(f)]
    rows = [r for r in rows if r["solvent"].strip().lower() in KEEP]
    rows.sort(key=lambda r: KEEP.index(r["solvent"].strip().lower()))
    return rows

def _axes_rect(col, row):
    """row 0 = riga superiore."""
    left = LEFT_CM + col * (AX_W_CM + GAP_X_CM)
    top  = TOP_CM + row * (AX_H_CM + GAP_Y_CM)
    bottom = FIG_H_CM - top - AX_H_CM
    return [left / FIG_W_CM, bottom / FIG_H_CM,
            AX_W_CM / FIG_W_CM, AX_H_CM / FIG_H_CM]

def make_grid(rows, with_text, outpath):
    fig = plt.figure(figsize=(FIG_W_CM * cm, FIG_H_CM * cm))

    ax_c = fig.add_axes(_axes_rect(0, 0))
    ax_d = fig.add_axes(_axes_rect(1, 0))
    ax_e = fig.add_axes(_axes_rect(0, 1)); ax_e_r = ax_e.twinx()
    ax_f = fig.add_axes(_axes_rect(1, 1)); ax_f_r = ax_f.twinx()

    rep = {}
    rep["c"] = panel_single(ax_c, rows, "rho_g_mL",  "Rm_MOhm", with_text)
    rep["d"] = panel_single(ax_d, rows, "eta_mPa_s", "Rs_MOhm", with_text)
    rep["e"] = panel_dual(ax_e, ax_e_r, rows, "mu_D",     with_text)
    rep["f"] = panel_dual(ax_f, ax_f_r, rows, "chi_mean", with_text)

    # C_eq e tau etichettati una sola volta
    ax_f.tick_params(labelleft=False)     # C_eq: numeri solo in (e)
    ax_e_r.tick_params(labelright=False)  # tau:  numeri solo in (f)
    if with_text:
        ax_e.set_ylabel(Y_LABELS["Ceq_nF"])
        ax_f_r.set_ylabel(Y_LABELS["tau_ms"])
        for ax, letter in ((ax_c, "c"), (ax_d, "d"), (ax_e, "e"), (ax_f, "f")):
            ax.text(-0.30, 1.04, f"({letter})", transform=ax.transAxes,
                    fontsize=7, fontweight="bold", ha="left", va="bottom")
        handles = [Line2D([], [], marker="o", ls="none", color="0.35", ms=3.5,
                          label=r"$C_{eq}$ (left)"),
                   Line2D([], [], marker="s", ls="none", mfc="none", mec="0.35",
                          ms=3.5, label=r"$\tau$ (right)")]
        ax_e.legend(handles=handles, loc="center left", frameon=False,
                    fontsize=5.5, handletextpad=0.3, borderpad=0.1,
                    labelspacing=0.25)

    fig.savefig(outpath, dpi=DPI_OUT, transparent=True)
    fig.savefig(os.path.splitext(outpath)[0] + ".svg", transparent=True)
    plt.close(fig)
    print("written", os.path.basename(outpath), f"({FIG_W_CM} x {FIG_H_CM:.2f} cm)")
    return rep

def main():
    rows = read_rows()
    os.makedirs(PLOTS, exist_ok=True)
    for wt, tag in ((False, "clean"), (True, "with_text")):
        rep = make_grid(rows, wt, os.path.join(PLOTS, f"fig3_grid_cdef_8cm_{tag}.png"))
    print("\n--- Theil-Sen (a, b, R2) ---")
    for panel in ("c", "d", "e", "f"):
        print(f"  ({panel}) {rep[panel]}")

if __name__ == "__main__":
    main()
