#!/usr/bin/env python3
"""
Find the best-fit parameters from a completed Dakota calibration run,
re-run hydroRaVENS with those parameters, and produce a diagnostic plot.

Figure layout
-------------
Left column  : precipitation (top, inverted) + observed/modelled discharge
Right column : flow duration curve (log scale) with observed BFI annotated

Reservoir structure is driven by driver.reservoir_order in params.yml.
Defaults to ['shallow', 'soil', 'karst'][:n_reservoirs] for backward
compatibility with experiments that predate that key.

Usage (from the experiment directory):
    python plot_best.py                      # uses dakota.dat, saves best_fit.png
    python plot_best.py --dat dakota_test.dat --save test_fit.png
"""

import argparse
import sys
from pathlib import Path
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from hydroravens import HydrographSeparation, run_and_score
from hydroravens.calibration import _nse, _kge, _log_kge, _kge_logfdc

try:
    with open('params.yml') as _f:
        CFG_TEMPLATE = yaml.safe_load(_f)['driver']['config_template']
except (FileNotFoundError, KeyError):
    CFG_TEMPLATE = 'cannon_cfg_template.yml'

OBJECTIVE_COL = 'neg_kge'
ROUTING_N     = 2      # Nash-cascade shape; must match driver.py

_MODULE_PARAMS_MAP = {
    'snowpack':      ['PDD_melt_factor'],
    'frozen_ground': ['log__fdd_threshold', 'snow_insulation_k'],
    'direct_runoff': ['f_direct_runoff'],
    'rain_on_snow':  [],
}

_RES_LABEL = {
    'shallow': 'sh',
    'soil':    'soil',
    'karst':   'karst',
    'deep':    'deep',
}


def _load_params_yml(path):
    """
    Return (metric, modules, params, reservoir_order, fixed_rec,
            enforce_wb, spin_up).
    """
    try:
        with open(path) as f:
            pcfg = yaml.safe_load(f)
        modules = pcfg.get('modules', {})
        params  = pcfg['parameters']
        for mod, names in _MODULE_PARAMS_MAP.items():
            if not modules.get(mod, True):
                for name in names:
                    if name in params:
                        params[name]['active'] = False
        drv       = pcfg['driver']
        n_res     = drv.get('n_reservoirs', 3)
        res_order = drv.get('reservoir_order',
                             ['shallow', 'soil', 'karst'][:n_res])
        fixed_rec  = drv.get('recession_exponents', None)
        enforce_wb = drv.get('enforce_water_balance', 'water-year')
        spin_up    = drv.get('spin_up_cycles', 0)
        return (drv['metric'], modules, params,
                res_order, fixed_rec, enforce_wb, spin_up)
    except FileNotFoundError:
        return ('KGE_logKGE', {}, {},
                ['shallow', 'soil', 'karst'], None, 'water-year', 0)


# ---------------------------------------------------------------------------
# Module-level globals — populated in __main__ after arg parsing
# ---------------------------------------------------------------------------
METRIC           = None
MODULES          = None
_PARAMS          = None
RESERVOIR_ORDER  = ['shallow', 'soil', 'karst']
_FIXED_RECESSION = None
ENFORCE_WB       = 'water-year'
SPIN_UP_CYCLES   = 0
INITIAL_STATES   = None


# ---------------------------------------------------------------------------
# Parameter helpers (use module-level globals)
# ---------------------------------------------------------------------------

def _is_active(name):
    return (_PARAMS or {}).get(name, {}).get('active', False)


def _get(row, name):
    """Return active param value from dat row, or fixed fallback."""
    if _is_active(name):
        return float(row[name])
    return float((_PARAMS or {}).get(name, {}).get('fixed', 0.0))


def read_best_params(dat_file):
    try:
        df = pd.read_csv(dat_file, sep=r'\s+')
    except FileNotFoundError:
        sys.exit(f'Error: {dat_file} not found. Run Dakota first.')
    df = df.rename(columns={'%eval_id': 'eval_id'})
    for col in df.columns:
        if col != 'interface':
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df.loc[df[OBJECTIVE_COL].idxmin()]


def _pdm_list(row):
    names = [f'log__H0_pdm_{l}' for l in RESERVOIR_ORDER]
    vals  = [10 ** _get(row, n) if _is_active(n) else None for n in names]
    return vals if any(v is not None for v in vals) else None


def _tile_list(row):
    """
    Tile-drain fractions in reservoir_order; None if no tile params present.
    Preserves the W/V/N/M series shared-tile convention when 'shallow' is in
    the order but f_tile_shallow is absent.
    """
    params = _PARAMS or {}
    names  = [f'f_tile_{l}' for l in RESERVOIR_ORDER]
    if not any(n in params for n in names):
        return None
    vals = [_get(row, n) if n in params else 0.0 for n in names]
    if ('shallow' in RESERVOIR_ORDER
            and 'f_tile_shallow' not in params
            and 'f_tile_soil' in params):
        idx_sh = RESERVOIR_ORDER.index('shallow')
        idx_so = RESERVOIR_ORDER.index('soil')
        vals[idx_sh] = vals[idx_so]
    return vals


def _tau_tile(row):
    if 'log__tau_tile' not in (_PARAMS or {}):
        return None
    return 10 ** _get(row, 'log__tau_tile')


def _hmax(row):
    key = f'log__Hmax_{RESERVOIR_ORDER[0]}'
    if key not in (_PARAMS or {}):
        return None
    return [10 ** _get(row, key)]


def _et_alpha(row):
    if 'et_alpha' not in (_PARAMS or {}):
        return None
    return _get(row, 'et_alpha')


def _wp_soil(row):
    if 'wp_soil' not in (_PARAMS or {}):
        return None
    return _get(row, 'wp_soil')


def _wp_soil_sigma(row):
    if 'wp_soil_sigma' not in (_PARAMS or {}):
        return None
    return _get(row, 'wp_soil_sigma')


def _recession_exponents(row):
    """Return (exponents, n_calibrated) using module-level globals."""
    params = _PARAMS or {}
    if _FIXED_RECESSION is not None:
        return list(_FIXED_RECESSION), 0
    has_shared = 'recession_b' in params
    exponents  = []
    n_cal      = 0
    for label in RESERVOIR_ORDER:
        label_key = f'recession_b_{label}'
        if label == 'shallow':
            exponents.append(1.0)
        elif has_shared:
            exponents.append(_get(row, 'recession_b'))
        elif label_key in params:
            exponents.append(_get(row, label_key))
            n_cal += 1
        else:
            exponents.append(1.0)
    if has_shared:
        n_cal = 1 if any(l != 'shallow' for l in RESERVOIR_ORDER) else 0
    if all(e == 1.0 for e in exponents):
        return None, 0
    return exponents, n_cal


def run_model(row):
    rec_exp, rec_k = _recession_exponents(row)
    return run_and_score(
        CFG_TEMPLATE,
        t_efold                = [10 ** _get(row, f'log__t_efold_{l}')
                                   for l in RESERVOIR_ORDER],
        f_to_discharge         = [_get(row, f'f_exfiltration_{l}')
                                   for l in RESERVOIR_ORDER[:-1]],
        melt_factor            =  _get(row, 'PDD_melt_factor'),
        fdd_threshold          =  10 ** _get(row, 'log__fdd_threshold'),
        snow_insulation_k      =  _get(row, 'snow_insulation_k'),
        Hmax                   =  _hmax(row),
        pdm_H0                 =  _pdm_list(row),
        f_tile                 =  _tile_list(row),
        tau_tile               =  _tau_tile(row),
        direct_runoff_fraction =  _get(row, 'f_direct_runoff'),
        baseflow_Q             =  _get(row, 'baseflow_Q'),
        et_alpha               =  _et_alpha(row),
        wp_soil                =  _wp_soil(row),
        wp_soil_sigma          =  _wp_soil_sigma(row),
        recession_exponents            = rec_exp,
        recession_exponents_calibrated = rec_k,
        routing_K              =  10 ** _get(row, 'log__routing_K'),
        routing_N              =  ROUTING_N,
        modules                =  MODULES,
        metric                 =  METRIC,
        enforce_water_balance  =  ENFORCE_WB,
        spin_up_cycles         =  SPIN_UP_CYCLES,
        initial_states         =  INITIAL_STATES,
    )


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def make_plot(result, params_row, save_path):
    b     = result.buckets
    aic   = result.aic

    mask  = (b.hydrodata['Specific Discharge (modeled) [mm/day]'].notna()
             & b.hydrodata['Specific Discharge [mm/day]'].notna())
    m_all = np.asarray(b.hydrodata.loc[mask, 'Specific Discharge (modeled) [mm/day]'])
    o_all = np.asarray(b.hydrodata.loc[mask, 'Specific Discharge [mm/day]'])
    nse        = _nse(m_all, o_all)
    kge        = _kge(m_all, o_all)
    log_kge    = _log_kge(m_all, o_all)
    kge_logfdc = result.kge_logfdc
    dates      = b.hydrodata['Date']

    fig = plt.figure(figsize=(14, 7))
    gs  = fig.add_gridspec(2, 2, width_ratios=[3, 1], height_ratios=[1, 2.5],
                           hspace=0.05, wspace=0.25)
    ax_p   = fig.add_subplot(gs[0, 0])
    ax_q   = fig.add_subplot(gs[1, 0], sharex=ax_p)
    ax_fdc = fig.add_subplot(gs[:, 1])

    ax_p.bar(dates, b.hydrodata['Precipitation [mm/day]'],
             width=1, color='steelblue', alpha=0.7)
    ax_p.set_ylabel('Precip.\n[mm/day]')
    ax_p.invert_yaxis()
    ax_p.yaxis.set_label_position('right')
    ax_p.yaxis.tick_right()
    plt.setp(ax_p.get_xticklabels(), visible=False)

    ax_q.plot(dates, b.hydrodata['Specific Discharge [mm/day]'],
              color='royalblue', lw=1.5, label='Observed')
    ax_q.plot(dates, b.hydrodata['Specific Discharge (modeled) [mm/day]'],
              color='k', lw=1.5, label='Modelled')
    ax_q.set_ylabel('Specific discharge [mm/day]')
    ax_q.set_xlabel('Date')
    ax_q.set_ylim(bottom=0)
    ax_q.legend(loc='upper right', fontsize=9)
    ax_q.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    plt.setp(ax_q.get_xticklabels(), rotation=30, ha='right')

    # ---- annotation text ----
    score_str = (f'logKGE = {log_kge:.3f}   NSE = {nse:.3f}   KGE = {kge:.3f}'
                 f'   KGE$_{{logFDC}}$ = {kge_logfdc:.3f}   AIC = {aic:.1f}')

    tau_parts = []
    for label in RESERVOIR_ORDER:
        val = 10 ** _get(params_row, f'log__t_efold_{label}')
        lbl = _RES_LABEL.get(label, label)
        fmt = '.1f' if val < 10 else '.0f'
        tau_parts.append(f'$\\tau_{{{lbl}}}$ = {val:{fmt}} d')
    tau_str = ',  '.join(tau_parts)

    f_parts = []
    for label in RESERVOIR_ORDER[:-1]:
        val = _get(params_row, f'f_exfiltration_{label}')
        lbl = _RES_LABEL.get(label, label)
        f_parts.append(f'$f_{{{lbl}}}$ = {val:.3f}')
    f_str = ',  '.join(f_parts)

    param_lines = (f'BFI: obs = {result.bfi_obs:.3f},  mod = {result.bfi_mod:.3f}\n'
                   + tau_str + '\n' + f_str)

    if _is_active('PDD_melt_factor'):
        param_lines += (f',  PDD = {_get(params_row, "PDD_melt_factor"):.2f}'
                        f' mm °C$^{{-1}}$ d$^{{-1}}$')

    hmax_key = f'log__Hmax_{RESERVOIR_ORDER[0]}'
    if _is_active(hmax_key):
        param_lines += f',  $H_{{max}}$ = {10**_get(params_row, hmax_key):.0f} mm'

    for label in RESERVOIR_ORDER:
        name = f'log__H0_pdm_{label}'
        if _is_active(name):
            lbl = _RES_LABEL.get(label, label)
            param_lines += (f',  $H_{{0,\\mathrm{{{lbl}}}}}$'
                            f' = {10**_get(params_row, name):.0f} mm')

    if _is_active('log__fdd_threshold'):
        param_lines += (f',  FDD$_{{thresh}}$'
                        f' = {10**_get(params_row, "log__fdd_threshold"):.0f} °C·d')
    if _is_active('snow_insulation_k'):
        param_lines += (f',  $k_{{ins}}$'
                        f' = {_get(params_row, "snow_insulation_k"):.4f}'
                        f' mm$^{{-1}}$ SWE')

    # Tile: shared (W/V/N/M convention) or per-reservoir
    params = _PARAMS or {}
    _shared_tile = ('shallow' in RESERVOIR_ORDER
                    and 'f_tile_shallow' not in params
                    and 'f_tile_soil' in params
                    and _is_active('f_tile_soil'))
    if _shared_tile:
        param_lines += f',  $f_{{tile}}$ = {_get(params_row, "f_tile_soil"):.3f}'
    else:
        for label in RESERVOIR_ORDER:
            name = f'f_tile_{label}'
            if _is_active(name):
                lbl = _RES_LABEL.get(label, label)
                param_lines += f',  $f_{{tile,{lbl}}}$ = {_get(params_row, name):.3f}'
    if _is_active('log__tau_tile'):
        param_lines += (f',  $\\tau_{{tile}}$'
                        f' = {10**_get(params_row, "log__tau_tile"):.1f} d')

    if _is_active('et_alpha'):
        param_lines += f',  $\\alpha_{{ET}}$ = {_get(params_row, "et_alpha"):.3f}'
    if _is_active('wp_soil'):
        param_lines += f',  $H_{{wp}}$ = {_get(params_row, "wp_soil"):.1f} mm'
    if _is_active('wp_soil_sigma'):
        param_lines += f',  $\\sigma_{{wp}}$ = {_get(params_row, "wp_soil_sigma"):.1f} mm'

    rec_exp, _ = _recession_exponents(params_row)
    if rec_exp is not None:
        labeled   = [(l, bv) for l, bv in zip(RESERVOIR_ORDER, rec_exp) if bv != 1.0]
        unique_b  = list(dict.fromkeys(f'{bv:.6f}' for _, bv in labeled))
        if len(unique_b) == 1:
            param_lines += f',  $b$ = {labeled[0][1]:.3f}'
        else:
            for label, bv in labeled:
                lbl = _RES_LABEL.get(label, label)
                param_lines += f',  $b_{{{lbl}}}$ = {bv:.3f}'

    routing_K = 10 ** _get(params_row, 'log__routing_K')
    param_lines += f'\n$K_{{route}}$ = {routing_K:.2f} d  (N={ROUTING_N})'
    if _is_active('f_direct_runoff'):
        param_lines += (f',  $\\gamma_{{direct}}$'
                        f' = {_get(params_row, "f_direct_runoff"):.3f}')
    if _is_active('baseflow_Q'):
        param_lines += f',  $Q_{{base}}$ = {_get(params_row, "baseflow_Q"):.4f} mm/d'

    ann = score_str + '\n' + param_lines
    ax_q.text(0.02, 0.97, ann, transform=ax_q.transAxes,
              va='top', fontsize=8.5,
              bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))

    ax_fdc.semilogy(result.fdc_obs.index, result.fdc_obs.values,
                    color='royalblue', lw=1.5, label='Observed')
    ax_fdc.semilogy(result.fdc_mod.index, result.fdc_mod.values,
                    color='k', lw=1.5, label='Modelled')
    ax_fdc.set_xlabel('Exceedance probability [%]')
    ax_fdc.set_ylabel('Specific discharge [mm/day]')
    ax_fdc.set_xlim(0, 100)
    ax_fdc.legend(fontsize=9)
    ax_fdc.set_title('Flow duration curve', fontsize=10)
    ax_fdc.grid(True, which='both', alpha=0.3)

    fig.suptitle('hydroRaVENS – Cannon River best-fit calibration', fontsize=13)
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f'Figure saved to {save_path}')


# ---------------------------------------------------------------------------
# __main__
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dat',     default='dakota.dat',   help='Dakota tabular data file')
    parser.add_argument('--save',    default='best_fit.png', help='Output figure path')
    parser.add_argument('--params',  default='params.yml',   help='params.yml config file')
    parser.add_argument('--no-show', action='store_true',    help='Save only; skip plt.show()')
    args = parser.parse_args()

    (METRIC, MODULES, _PARAMS,
     RESERVOIR_ORDER, _FIXED_RECESSION,
     ENFORCE_WB, SPIN_UP_CYCLES) = _load_params_yml(args.params)

    # Data-driven initial reservoir states — mirrors driver.py exactly
    cfg_path   = Path(CFG_TEMPLATE)
    cfg_dir    = cfg_path.parent
    with open(cfg_path) as _f:
        _mcfg = yaml.safe_load(_f)
    _datafile   = _mcfg['timeseries']['datafile']
    _data_path  = (Path(_datafile) if Path(_datafile).is_absolute()
                   else cfg_dir / _datafile)
    _area_km2   = _mcfg['catchment']['drainage_basin_area__km2']
    _df_raw     = pd.read_csv(_data_path, parse_dates=['Date'])
    _Q_specific = _df_raw['Discharge [m^3/s]'].values * 86400.0 / (_area_km2 * 1e3)
    _precip     = _df_raw['Precipitation [mm/day]'].values
    _hs         = HydrographSeparation(_Q_specific,
                                        n_reservoirs=len(RESERVOIR_ORDER),
                                        precip=_precip)
    _hs.fit()
    INITIAL_STATES = {'reservoirs': _hs.get_initial_conditions()['H0']}

    best = read_best_params(args.dat)

    # Summary to stdout
    print(f'\nBest evaluation: {int(best["eval_id"])}')
    print(f'  metric           = {METRIC}')
    print(f'  reservoir_order  = {RESERVOIR_ORDER}')
    for label in RESERVOIR_ORDER:
        val = 10 ** _get(best, f'log__t_efold_{label}')
        lbl = _RES_LABEL.get(label, label)
        print(f'  t_efold_{lbl:<8} = {val:.1f} days')
    for label in RESERVOIR_ORDER[:-1]:
        val = _get(best, f'f_exfiltration_{label}')
        lbl = _RES_LABEL.get(label, label)
        print(f'  f_exfilt_{lbl:<7} = {val:.4f}')
    if _is_active('PDD_melt_factor'):
        print(f'  PDD_melt_factor  = {_get(best, "PDD_melt_factor"):.4f} mm/°C/day')
    hmax_key = f'log__Hmax_{RESERVOIR_ORDER[0]}'
    if _is_active(hmax_key):
        print(f'  Hmax_{RESERVOIR_ORDER[0]:<10} = {10**_get(best, hmax_key):.1f} mm')
    for label in RESERVOIR_ORDER:
        name = f'log__H0_pdm_{label}'
        if _is_active(name):
            print(f'  {name:<22} = {10**_get(best, name):.1f} mm')
    if _is_active('log__fdd_threshold'):
        print(f'  fdd_threshold    = {10**_get(best, "log__fdd_threshold"):.1f} °C·day')
    if _is_active('snow_insulation_k'):
        print(f'  snow_insulation_k= {_get(best, "snow_insulation_k"):.4f} mm⁻¹ SWE')
    params = _PARAMS or {}
    _shared_tile = ('shallow' in RESERVOIR_ORDER
                    and 'f_tile_shallow' not in params
                    and 'f_tile_soil' in params
                    and _is_active('f_tile_soil'))
    if _shared_tile:
        print(f'  f_tile (shared)  = {_get(best, "f_tile_soil"):.4f}')
    else:
        for label in RESERVOIR_ORDER:
            name = f'f_tile_{label}'
            if _is_active(name):
                print(f'  {name:<22} = {_get(best, name):.4f}')
    if _is_active('log__tau_tile'):
        print(f'  tau_tile         = {10**_get(best, "log__tau_tile"):.2f} days')
    if _is_active('et_alpha'):
        print(f'  et_alpha         = {_get(best, "et_alpha"):.4f}')
    rec_exp, _ = _recession_exponents(best)
    if rec_exp is not None:
        labeled  = [(l, bv) for l, bv in zip(RESERVOIR_ORDER, rec_exp) if bv != 1.0]
        unique_b = list(dict.fromkeys(f'{bv:.6f}' for _, bv in labeled))
        if len(unique_b) == 1:
            print(f'  recession_b      = {labeled[0][1]:.4f}')
        else:
            for label, bv in labeled:
                print(f'  recession_b_{label:<8} = {bv:.4f}')
    if _is_active('wp_soil'):
        print(f'  wp_soil          = {_get(best, "wp_soil"):.2f} mm')
    if _is_active('wp_soil_sigma'):
        print(f'  wp_soil_sigma    = {_get(best, "wp_soil_sigma"):.2f} mm')
    if _is_active('baseflow_Q'):
        print(f'  baseflow_Q       = {_get(best, "baseflow_Q"):.4f} mm/day')
    rK = 10 ** _get(best, 'log__routing_K')
    print(f'  routing_K        = {rK:.3f} days  (N={ROUTING_N},'
          f' mean travel time = {ROUTING_N * rK:.2f} days)')
    if _is_active('f_direct_runoff'):
        print(f'  f_direct_runoff  = {_get(best, "f_direct_runoff"):.4f}')

    result = run_model(best)
    b_     = result.buckets
    mask   = (b_.hydrodata['Specific Discharge (modeled) [mm/day]'].notna()
              & b_.hydrodata['Specific Discharge [mm/day]'].notna())
    m_all  = np.asarray(b_.hydrodata.loc[mask, 'Specific Discharge (modeled) [mm/day]'])
    o_all  = np.asarray(b_.hydrodata.loc[mask, 'Specific Discharge [mm/day]'])
    print(f'  logKGE           = {_log_kge(m_all, o_all):.4f}')
    print(f'  NSE              = {_nse(m_all, o_all):.4f}')
    print(f'  KGE              = {_kge(m_all, o_all):.4f}')
    print(f'  KGE_logFDC       = {result.kge_logfdc:.4f}')
    print(f'  AIC              = {result.aic:.2f}')
    print(f'  BFI obs          = {result.bfi_obs:.4f}')
    print(f'  BFI mod          = {result.bfi_mod:.4f}')

    make_plot(result, best, save_path=args.save)
    if not args.no_show:
        plt.show()
