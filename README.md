# Wickert2026-hydroRaVENS-model-design

Model-selection calibration experiments for the hydroRaVENS paper.
Seven Dakota optimisations test whether individual process additions
improve fit over a 3-reservoir + snowpack baseline, using AIC to
penalise added parameters.

## Experiment design

Each experiment adds one process to the M01 baseline (3-reservoir + snowpack).
All models are calibrated against the 1991–2011 Cannon River record using
the `logKGE_logFDC_BFI` composite metric.

| ID  | Description                              | Added parameter(s)              | k |
|-----|------------------------------------------|---------------------------------|---|
| M00 | 3-reservoir only (no snowpack)           | —                               | 5 |
| M01 | 3-reservoir + snowpack *(baseline)*      | PDD melt factor                 | 6 |
| M02 | M01 + Nash-cascade routing               | K_route                         | 7 |
| M03 | M01 + saturation-excess threshold        | H_max (shallow reservoir)       | 7 |
| M04 | M01 + frozen-ground index                | FDD threshold                   | 7 |
| M05 | M01 + frozen-ground index + insulation   | FDD threshold + k_insulation    | 8 |
| M06 | M01 + regional baseflow import           | Q_base                          | 7 |

The five parameters shared by every experiment are the three e-folding
reservoir times (τ_shallow, τ_soil, τ_karst) and the two exfiltration
fractions (f_shallow, f_soil).

## Repository layout

```
cannon_river/
├── Cannon1991-2011Input.csv          input time series (1991-10-01 – 2011-09-30)
├── cannon_cfg_1991_2011_et_tc.yml    hydroRaVENS config template (Thornthwaite-Chang ET)
├── driver.py                         Dakota driver (reads params.yml)
├── run_driver.sh                     wrapper: calls driver.py in the dakota-env
├── generate_dakota_in.py             generates dakota.in from params.yml
├── plot_best.py                      best-fit diagnostic plot
├── archive_run.sh                    archives a completed run to runs/<name>/
├── run.sh                            single-experiment entry point
├── run_all.sh                        runs M00–M06 in series
├── compare_aic.py                    collects AIC + metrics table across experiments
└── experiments/
    ├── M00/   params.yml  dakota.in  runs/  [symlinks to shared files]
    ├── M01/   ...
    ├── M02/   ...
    ├── M03/   ...
    ├── M04/   ...
    ├── M05/   ...
    └── M06/   ...
```

Shared infrastructure lives at `cannon_river/`; each experiment directory
contains only its own `params.yml` and `dakota.in`, with symlinks to
everything else.

## Running the experiments

### All experiments in series (recommended)

```bash
cd cannon_river
bash run_all.sh v1
```

Each experiment runs to completion before the next starts.
Results are archived to `experiments/<MXX>/runs/<timestamp>_v1/`.

### Single experiment

```bash
cd cannon_river/experiments/M01
bash run.sh v1
```

`run.sh` regenerates `dakota.in` from `params.yml` automatically before
calling Dakota, so editing `params.yml` is all that is needed to change
the run configuration.

### Regenerate `dakota.in` without running

```bash
cd cannon_river/experiments/M01
python generate_dakota_in.py
```

## Collecting results

After all experiments have at least one archived run:

```bash
cd cannon_river
python compare_aic.py              # prints table to stdout
python compare_aic.py --csv results_table.csv
python compare_aic.py --experiments M01 M03 M05   # subset
```

`compare_aic.py` finds the most-recently archived run in each experiment,
re-runs hydroRaVENS with the best-fit parameters, and reports logKGE, NSE,
KGE, KGE_logFDC, AIC, BFI_obs, and BFI_mod.

## Dependencies

- [hydroRaVENS](https://github.com/MNiMORPH/hydroRaVENS) — installed in the
  `dakota-env` conda environment
- [Dakota](https://dakota.sandia.gov/) ≥ 6.18 — likewise in `dakota-env`
