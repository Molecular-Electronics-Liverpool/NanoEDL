# NanoEDL

Post-processing code for

> **Probing Solvent Dynamics at the Single-Molecule Limit**, *Nano Letters* (2026).
> DOI: [10.1021/acs.nanolett.6c03272](https://doi.org/10.1021/acs.nanolett.6c03272)

The scripts take the raw AC-STMBJ traces, extract the parameters of the equivalent circuit (R<sub>s</sub>, R<sub>m</sub>, C<sub>eq</sub>, τ) for each solvent, and generate Figure 3c–f.

## Raw data

All raw data acquired in Liverpool are available free of charge from the University of Liverpool Data Catalogue (deposit no. 3083, DOI: [10.17638/datacat.liverpool.ac.uk/3083](https://doi.org/10.17638/datacat.liverpool.ac.uk/3083)).

Download the traces and place them in a folder named `00_raw/` at the root of this repository, one file per solvent:

```
00_raw/newfile_<solvent>.txt
```

with `<solvent>` one of `mesitylene`, `decane`, `anisole`, `tcb`, `dcb`, `gvl`, `gbl`, `pc`. Each file is tab-separated with three columns: sample index, bias V (V), and log<sub>10</sub>(G/G<sub>0</sub>). Consecutive blocks of 4000 samples (40 kHz, 100 ms) form one trace.

## Requirements

Python ≥ 3.10 with

```
pip install numpy pandas matplotlib
```

## Running the analysis

Run the scripts in order from the repository root:

```
python3 step01_median.py
python3 step02_detrend.py
python3 step03_windows.py
python3 step04_fit.py
python3 step05_tables.py
python3 extra_fig3_grid.py
```

Each step reads the output of the previous one and writes to its own folder, created on first run:

| Step | Script | What it does | Output |
|---|---|---|---|
| 1 | `step01_median.py` | Splits the raw record into traces and takes the point-by-point median over all traces; reconstructs the current I = G<sub>0</sub>·10<sup>log(G/G<sub>0</sub>)</sup>·V | `01_median/` |
| 2 | `step02_detrend.py` | Removes a linear baseline from the median current | `02_detrended/` |
| 3 | `step03_windows.py` | Cuts the trace into 10 ms half-periods of the 50 Hz square wave and measures the plateau voltage and the voltage step at each edge | `03_windows/` |
| 4 | `step04_fit.py` | Fits i(t) = A·exp(−t/τ) + B in each window and converts the fit into circuit parameters | `04_fits/fit_per_window.csv` |
| 5 | `step05_tables.py` | Averages the valid windows for each solvent and adds the solvent properties | `05_tables/` |
| 6 | `extra_fig3_grid.py` | Plots Figure 3c–f | `06_plots/` |

Step 1 reads about 3.3 GB of raw data and is by far the slowest; steps 4–6 run in seconds.

## Method and conventions

**Fit.** τ and the plateau B come from the return map I(k+1) = n·I(k) + c, fitted with the Theil–Sen estimator: τ = −Δt / ln n and B = c / (1 − n). The transient amplitude A is then a linear least-squares fit with τ and B fixed.

**Circuit parameters.** All voltages are measured window by window from the bias trace:

```
Rs  = V_ref       / |A|
Rm  = |V_plateau| / |B|
Ceq = tau / Rs
```

With `VOLTAGE_MODE = "step"` (the setting used in the paper), V<sub>ref</sub> is the voltage step at the edge that triggers the transient — 2 V for the inner half-periods. At the edge, the voltage across R<sub>s</sub> is the difference between the applied bias and the voltage already on the capacitor, which has fully charged to −V<sub>0</sub> during the previous half-period. Setting `VOLTAGE_MODE = "plateau"` uses the applied bias (≈ 1 V) instead, which doubles C<sub>eq</sub> and halves R<sub>s</sub>; R<sub>m</sub> and τ do not change.

**Windows.** Windows 1 and 10 are the DC-bias segments before and after the square wave and are excluded. For anisole, windows 2 and 3 are also excluded, since the response has not yet settled. Exclusions are set in `config.py` (`SKIP_WINDOWS`, `EXTRA_SKIP`).

**Solvent properties.** Dielectric constant, viscosity, dipole moment and density at 25 °C are collected in `config.py` (`PROPS`). The mean Pauling electronegativity used in Figure 3e–f is the arithmetic mean over all atoms of the molecule (`chi_mean`).

All physical constants, acquisition settings and solvent properties live in `config.py`; no other file needs editing to change them.

## Citation

If you use this code, please cite the paper above.

## Licence

See `LICENSE`.
