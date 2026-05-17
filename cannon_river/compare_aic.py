#!/usr/bin/env python3
"""
compare_aic.py — Collect best-fit metrics from all model-selection experiments
and print a summary table suitable for the paper.

For each experiment with the given prefix (e.g. M00–M06 or W00–W02):
  1. Locate the most-recently archived run.
  2. Read the best-fit row from evaluations.dat (lowest neg_kge).
  3. Re-run hydroRaVENS with those parameters (same logic as plot_best.py).
  4. Extract logKGE, NSE, KGE, KGE_logFDC, AIC, BFI_obs, BFI_mod.

Reservoir structure is driven by driver.reservoir_order in each params.yml.
Defaults to ['shallow', 'soil', 'karst'][:n_reservoirs] for backward
compatibility with experiments that predate that key.

Usage (from cannon_river/):
    python compare_aic.py W
    python compare_aic.py W --experiments W01 W02   # subset
    python compare_aic.py W --csv results_table.csv
"""

import argparse
import sys
import warnings
import yaml
import pandas as pd
import numpy as np
from pathlib import Path
from hydroravens import HydrographSeparation, run_and_score
from hydroravens.calibration import _nse, _kge, _log_kge

# Intentional: et_scale carries explicit responsibility for the water balance
# in Y-series experiments that use enforce_water_balance='none'.
warnings.filterwarnings('ignore', message=r"enforce_water_balance='none'",
                        category=UserWarning)

ROUTING_N = 2

_MODULE_PARAMS_MAP = {
    'snowpack':      ['PDD_melt_factor'],
    'frozen_ground': ['log__fdd_threshold', 'snow_insulation_k'],
    'direct_runoff': ['f_direct_runoff'],
    'rain_on_snow':  [],
}

# Short labels used in annotations; extend as new reservoir types are added.
_RES_LABEL = {
    'shallow': 'sh',
    'soil':    'soil',
    'karst':   'karst',
    'deep':    'deep',
}


def _load_params_yml(path):
    """
    Return (metric, modules, params, cfg_template,
            reservoir_order, fixed_rec, enforce_wb, spin_up).

    reservoir_order is a list of reservoir labels in cascade order, e.g.
    ['shallow', 'soil', 'karst'] or ['soil', 'karst', 'deep'].  It drives
    all parameter-name → reservoir-index mappings.

    fixed_rec, if not None, is a list of pre-fixed recession exponents
    (one per reservoir) read from driver.recession_exponents.
    """
    with open(path) as f:
        pcfg = yaml.safe_load(f)
    modules = pcfg.get('modules', {})
    params  = pcfg['parameters']
    for mod, names in _MODULE_PARAMS_MAP.items():
        if not modules.get(mod, True):
            for name in names:
                if name in params:
                    params[name]['active'] = False
    drv = pcfg['driver']
    cfg_template = drv['config_template']
    metric       = drv['metric']
    n_res        = drv.get('n_reservoirs', 3)
    enforce_wb   = drv.get('enforce_water_balance', 'water-year')
    spin_up      = drv.get('spin_up_cycles', 0)
    fixed_rec    = drv.get('recession_exponents', None)
    # Default preserves the historical shallow/soil/karst ordering.
    reservoir_order = drv.get('reservoir_order',
                               ['shallow', 'soil', 'karst'][:n_res])
    return (metric, modules, params, cfg_template,
            reservoir_order, fixed_rec, enforce_wb, spin_up)


def _is_active(params, name):
    return params.get(name, {}).get('active', False)


def _get(row_or_dict, params, name):
    """Return active param value from dat row, or fixed fallback from params."""
    if _is_active(params, name):
        return float(row_or_dict[name])
    return float(params.get(name, {}).get('fixed', 0.0))


def _best_row(evaluations_dat):
    df = pd.read_csv(evaluations_dat, sep=r'\s+')
    df = df.rename(columns={'%eval_id': 'eval_id'})
    for col in df.columns:
        if col != 'interface':
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df.loc[df['neg_kge'].idxmin()]


def _latest_run_dir(exp_dir):
    runs_path = Path(exp_dir, 'runs')
    if not runs_path.is_dir():
        return None
    runs = sorted(runs_path.iterdir())
    runs = [r for r in runs if r.is_dir()]
    return runs[-1] if runs else None


# ---------------------------------------------------------------------------
# Reservoir-order–aware parameter helpers
# ---------------------------------------------------------------------------

def _pdm_list(row, params, reservoir_order):
    """PDM H0 list in reservoir_order; None entry where PDM is inactive."""
    names = [f'log__H0_pdm_{label}' for label in reservoir_order]
    vals  = [10 ** _get(row, params, n) if _is_active(params, n) else None
             for n in names]
    return vals if any(v is not None for v in vals) else None


def _tile_list(row, params, reservoir_order):
    """
    Tile-drain fraction list in reservoir_order order, or None if no tile
    parameters are present.

    Legacy compatibility: when 'shallow' is in reservoir_order but
    f_tile_shallow is absent while f_tile_soil is present, shallow and soil
    share the same tile network (W/V/N/M series convention).
    """
    names = [f'f_tile_{label}' for label in reservoir_order]
    if not any(n in params for n in names):
        return None
    vals = [_get(row, params, n) if n in params else 0.0 for n in names]
    if ('shallow' in reservoir_order
            and 'f_tile_shallow' not in params
            and 'f_tile_soil' in params):
        idx_sh = reservoir_order.index('shallow')
        idx_so = reservoir_order.index('soil')
        vals[idx_sh] = vals[idx_so]
    return vals


def _tau_tile(row, params):
    if 'log__tau_tile' not in params:
        return None
    return 10 ** _get(row, params, 'log__tau_tile')


def _hmax(row, params, reservoir_order):
    """Hmax list for the first reservoir; None if not in params."""
    key = f'log__Hmax_{reservoir_order[0]}'
    if key not in params:
        return None
    return [10 ** _get(row, params, key)]


def _et_scale(row, params):
    return _get(row, params, 'et_scale') if 'et_scale' in params else None


def _et_alpha(row, params):
    return _get(row, params, 'et_alpha') if 'et_alpha' in params else None


def _et_scale(row, params):
    return _get(row, params, 'et_scale') if 'et_scale' in params else None


def _wp_soil(row, params):
    return _get(row, params, 'wp_soil') if 'wp_soil' in params else None


def _wp_soil_sigma(row, params):
    return _get(row, params, 'wp_soil_sigma') if 'wp_soil_sigma' in params else None


def _recession_exponents(row, params, reservoir_order, fixed_rec):
    """
    Build (exponents, n_calibrated) for run_and_score.

    Priority:
      1. fixed_rec from driver.recession_exponents (e.g. W00)
      2. shared recession_b parameter (applies to all non-shallow reservoirs)
      3. per-reservoir recession_b_{label} parameters
      4. default 1.0 (linear) for unlisted or 'shallow' reservoirs
    Returns (None, 0) if all exponents are 1.0.
    """
    if fixed_rec is not None:
        return list(fixed_rec), 0

    has_shared = 'recession_b' in params
    exponents  = []
    n_cal      = 0
    for label in reservoir_order:
        label_key = f'recession_b_{label}'
        if label == 'shallow':
            exponents.append(1.0)
        elif has_shared:
            exponents.append(_get(row, params, 'recession_b'))
        elif label_key in params:
            exponents.append(_get(row, params, label_key))
            if params[label_key].get('active', False):
                n_cal += 1
        else:
            exponents.append(1.0)
    if has_shared:
        n_cal = 1 if any(l != 'shallow' for l in reservoir_order) else 0
    if all(e == 1.0 for e in exponents):
        return None, 0
    return exponents, n_cal


def _initial_states(cfg_template, reservoir_order, exp_dir):
    """Compute data-driven initial reservoir states via HydrographSeparation."""
    cfg_path = Path(exp_dir) / cfg_template
    cfg_dir  = cfg_path.parent
    with open(cfg_path) as f:
        model_cfg = yaml.safe_load(f)
    area_km2   = model_cfg['catchment']['drainage_basin_area__km2']
    datafile   = model_cfg['timeseries']['datafile']
    data_path  = (Path(datafile) if Path(datafile).is_absolute()
                  else cfg_dir / datafile)
    df         = pd.read_csv(data_path, parse_dates=['Date'])
    Q_specific = df['Discharge [m^3/s]'].values * 86400.0 / (area_km2 * 1e3)
    precip     = df['Precipitation [mm/day]'].values
    hs         = HydrographSeparation(Q_specific,
                                       n_reservoirs=len(reservoir_order),
                                       precip=precip)
    hs.fit()
    return {'reservoirs': hs.get_initial_conditions()['H0']}


# ---------------------------------------------------------------------------
# Model runner
# ---------------------------------------------------------------------------

def _run_model(row, params, modules, metric, cfg_template, exp_dir,
               reservoir_order, fixed_rec, enforce_wb='water-year', spin_up=0):
    rec_exp, rec_k = _recession_exponents(row, params, reservoir_order, fixed_rec)
    cfg_path = str(Path(exp_dir) / cfg_template)
    init     = _initial_states(cfg_template, reservoir_order, exp_dir)
    return run_and_score(
        cfg_path,
        t_efold                = [10 ** _get(row, params, f'log__t_efold_{l}')
                                   for l in reservoir_order],
        f_to_discharge         = [_get(row, params, f'f_exfiltration_{l}')
                                   for l in reservoir_order[:-1]],
        melt_factor            =  _get(row, params, 'PDD_melt_factor'),
        fdd_threshold          =  10 ** _get(row, params, 'log__fdd_threshold'),
        snow_insulation_k      =  _get(row, params, 'snow_insulation_k'),
        Hmax                   =  _hmax(row, params, reservoir_order),
        pdm_H0                 =  _pdm_list(row, params, reservoir_order),
        f_tile                 =  _tile_list(row, params, reservoir_order),
        tau_tile               =  _tau_tile(row, params),
        direct_runoff_fraction =  _get(row, params, 'f_direct_runoff'),
        baseflow_Q             =  _get(row, params, 'baseflow_Q'),
        et_scale               =  _et_scale(row, params),
        et_alpha               =  _et_alpha(row, params),
        wp_soil                =  _wp_soil(row, params),
        wp_soil_sigma          =  _wp_soil_sigma(row, params),
        recession_exponents            = rec_exp,
        recession_exponents_calibrated = rec_k,
        routing_K              =  10 ** _get(row, params, 'log__routing_K'),
        routing_N              =  ROUTING_N,
        modules                =  modules,
        metric                 =  metric,
        enforce_water_balance  =  enforce_wb,
        spin_up_cycles         =  spin_up,
        initial_states         =  init,
    )


def _n_active(params):
    return sum(1 for p in params.values() if p.get('active', False))


# ---------------------------------------------------------------------------
# Per-experiment processing
# ---------------------------------------------------------------------------

def process_experiment(exp_dir):
    exp_name    = Path(exp_dir).name
    params_path = Path(exp_dir, 'params.yml')
    if not params_path.exists():
        return None

    (metric, modules, params, cfg_template,
     reservoir_order, fixed_rec,
     enforce_wb, spin_up) = _load_params_yml(params_path)
    k = _n_active(params)

    run_dir = _latest_run_dir(exp_dir)
    if run_dir is None:
        print(f'  {exp_name}: no archived runs — skipping', file=sys.stderr)
        return None
    eval_file = run_dir / 'evaluations.dat'
    if not eval_file.exists():
        print(f'  {exp_name}: evaluations.dat missing in {run_dir.name} — skipping',
              file=sys.stderr)
        return None

    print(f'  {exp_name}: reading {run_dir.name} ...', file=sys.stderr)
    best = _best_row(eval_file)

    result = _run_model(best, params, modules, metric, cfg_template, exp_dir,
                        reservoir_order, fixed_rec, enforce_wb, spin_up)
    b    = result.buckets
    mask = (b.hydrodata['Specific Discharge (modeled) [mm/day]'].notna()
            & b.hydrodata['Specific Discharge [mm/day]'].notna())
    m    = np.asarray(b.hydrodata.loc[mask, 'Specific Discharge (modeled) [mm/day]'])
    o    = np.asarray(b.hydrodata.loc[mask, 'Specific Discharge [mm/day]'])

    return {
        'Model':      exp_name,
        'k':          k,
        'logKGE':     round(_log_kge(m, o), 3),
        'NSE':        round(_nse(m, o),     3),
        'KGE':        round(_kge(m, o),     3),
        'KGE_logFDC': round(result.kge_logfdc, 3),
        'AIC':        round(result.aic,         1),
        'BFI_obs':    round(result.bfi_obs,     3),
        'BFI_mod':    round(result.bfi_mod,     3),
        'run':        run_dir.name,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('prefix',
                        help='Experiment prefix letter (e.g. M or W)')
    parser.add_argument('--experiments', nargs='+', default=None,
                        help='Subset of experiment names (default: all with prefix)')
    parser.add_argument('--csv', default=None, help='Save table to CSV')
    args = parser.parse_args()

    experiments_root = Path('experiments')
    if not experiments_root.is_dir():
        sys.exit('Run from cannon_river/ — experiments/ directory not found.')

    if args.experiments:
        exp_dirs = [experiments_root / e for e in args.experiments]
    else:
        exp_dirs = sorted(experiments_root.iterdir())
        exp_dirs = [d for d in exp_dirs
                    if d.is_dir() and d.name.startswith(args.prefix)]

    print('Collecting results...\n', file=sys.stderr)
    rows = []
    for d in exp_dirs:
        r = process_experiment(d)
        if r:
            rows.append(r)

    if not rows:
        sys.exit('No results found.')

    df   = pd.DataFrame(rows).set_index('Model')
    cols = ['k', 'logKGE', 'NSE', 'KGE', 'KGE_logFDC',
            'AIC', 'BFI_obs', 'BFI_mod', 'run']
    df   = df[cols]
    print('\n' + df.to_string())

    if args.csv:
        df.to_csv(args.csv)
        print(f'\nSaved to {args.csv}')


if __name__ == '__main__':
    main()
