#!/usr/bin/env python3
"""
compare_aic.py — Collect best-fit metrics from all model-selection experiments
and print a summary table suitable for the paper.

For each experiment with the given prefix (e.g. M00–M06 or N00–N04):
  1. Locate the most-recently archived run.
  2. Read the best-fit row from evaluations.dat (lowest neg_kge).
  3. Re-run hydroRaVENS with those parameters (same logic as plot_best.py).
  4. Extract logKGE, NSE, KGE, KGE_logFDC, AIC, BFI_obs, BFI_mod.

Usage (from cannon_river/):
    python compare_aic.py M
    python compare_aic.py N --experiments N01 N03   # subset
    python compare_aic.py M --csv results_table.csv
"""

import argparse
import os
import sys
import math
import yaml
import pandas as pd
import numpy as np
from pathlib import Path
from hydroravens import run_and_score
from hydroravens.calibration import _nse, _kge, _log_kge

ROUTING_N = 2

_MODULE_PARAMS_MAP = {
    'snowpack':      ['PDD_melt_factor'],
    'frozen_ground': ['log__fdd_threshold', 'snow_insulation_k'],
    'direct_runoff': ['f_direct_runoff'],
    'rain_on_snow':  [],
}


def _load_params_yml(path):
    with open(path) as f:
        pcfg = yaml.safe_load(f)
    modules = pcfg.get('modules', {})
    params  = pcfg['parameters']
    for mod, names in _MODULE_PARAMS_MAP.items():
        if not modules.get(mod, True):
            for name in names:
                if name in params:
                    params[name]['active'] = False
    cfg_template = pcfg['driver']['config_template']
    metric       = pcfg['driver']['metric']
    n_reservoirs = pcfg['driver'].get('n_reservoirs', 3)
    return metric, modules, params, cfg_template, n_reservoirs


def _is_active(params, name):
    return params.get(name, {}).get('active', False)


def _get(row_or_dict, params, name):
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
    runs = sorted(Path(exp_dir, 'runs').iterdir())
    runs = [r for r in runs if r.is_dir()]
    if not runs:
        return None
    return runs[-1]


def _run_model(row, params, modules, metric, cfg_template, exp_dir, n_reservoirs=3):
    g = lambda name: _get(row, params, name)
    t_efold_all      = [10 ** g('log__t_efold_shallow'),
                        10 ** g('log__t_efold_soil'),
                        10 ** g('log__t_efold_karst')]
    f_discharge_all  = [g('f_exfiltration_shallow'),
                        g('f_exfiltration_soil')]
    cfg_path = str(Path(exp_dir, cfg_template))
    return run_and_score(
        cfg_path,
        t_efold               = t_efold_all[:n_reservoirs],
        f_to_discharge        = f_discharge_all[:n_reservoirs - 1],
        melt_factor           =  g('PDD_melt_factor'),
        fdd_threshold         =  10 ** g('log__fdd_threshold'),
        snow_insulation_k     =  g('snow_insulation_k'),
        Hmax                  = [10 ** g('log__Hmax_shallow')],
        direct_runoff_fraction=  g('f_direct_runoff'),
        baseflow_Q            =  g('baseflow_Q'),
        routing_K             =  10 ** g('log__routing_K'),
        routing_N             =  ROUTING_N,
        modules               =  modules,
        metric                =  metric,
    )


def _n_active(params):
    return sum(1 for p in params.values() if p.get('active', False))


def process_experiment(exp_dir):
    exp_name = Path(exp_dir).name
    params_path = Path(exp_dir, 'params.yml')
    if not params_path.exists():
        return None

    metric, modules, params, cfg_template, n_reservoirs = _load_params_yml(params_path)
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

    result = _run_model(best, params, modules, metric, cfg_template, exp_dir, n_reservoirs)
    b    = result.buckets
    mask = (b.hydrodata['Specific Discharge (modeled) [mm/day]'].notna()
            & b.hydrodata['Specific Discharge [mm/day]'].notna())
    m    = np.asarray(b.hydrodata.loc[mask, 'Specific Discharge (modeled) [mm/day]'])
    o    = np.asarray(b.hydrodata.loc[mask, 'Specific Discharge [mm/day]'])

    return {
        'Model':       exp_name,
        'k':           k,
        'logKGE':      round(_log_kge(m, o), 3),
        'NSE':         round(_nse(m, o),     3),
        'KGE':         round(_kge(m, o),     3),
        'KGE_logFDC':  round(result.kge_logfdc, 3),
        'AIC':         round(result.aic,         1),
        'BFI_obs':     round(result.bfi_obs,     3),
        'BFI_mod':     round(result.bfi_mod,     3),
        'run':         run_dir.name,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('prefix',
                        help='Experiment prefix letter (e.g. M or N)')
    parser.add_argument('--experiments', nargs='+',
                        default=None, help='Subset of experiment names (default: all with prefix)')
    parser.add_argument('--csv', default=None, help='Save table to CSV')
    args = parser.parse_args()

    experiments_root = Path('experiments')
    if not experiments_root.is_dir():
        sys.exit('Run from cannon_river/ — experiments/ directory not found.')

    if args.experiments:
        exp_dirs = [experiments_root / e for e in args.experiments]
    else:
        exp_dirs = sorted(experiments_root.iterdir())
        exp_dirs = [d for d in exp_dirs if d.is_dir() and d.name.startswith(args.prefix)]

    print('Collecting results...\n', file=sys.stderr)
    rows = []
    for d in exp_dirs:
        r = process_experiment(d)
        if r:
            rows.append(r)

    if not rows:
        sys.exit('No results found.')

    df = pd.DataFrame(rows).set_index('Model')
    cols = ['k', 'logKGE', 'NSE', 'KGE', 'KGE_logFDC', 'AIC', 'BFI_obs', 'BFI_mod', 'run']
    df = df[cols]

    print('\n' + df.to_string())

    if args.csv:
        df.to_csv(args.csv)
        print(f'\nSaved to {args.csv}')


if __name__ == '__main__':
    main()
