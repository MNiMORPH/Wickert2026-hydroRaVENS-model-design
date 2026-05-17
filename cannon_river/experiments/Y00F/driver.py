#!/usr/bin/env python3
"""
Generic Dakota driver for hydroRaVENS calibration.
Reads all reservoir structure and parameters from params.yml — no hard-coded names.

Compatible with any reservoir_order (2-res, 3-res, ...) and any combination of
active/fixed parameters. Designed for the X/Y-series experiments.
"""

import yaml
import pandas as pd
import dakota.interfacing as di
import numpy as np
from hydroravens import HydrographSeparation, run_and_score

with open('params.yml') as f:
    _cfg = yaml.safe_load(f)

_driver    = _cfg['driver']
_param_cfg = _cfg['parameters']

CONFIG_TEMPLATE = _driver['config_template']
METRIC          = _driver['metric']
SPIN_UP_CYCLES  = _driver['spin_up_cycles']
ROUTING_N       = _driver['routing_N']
DECADE_START    = _driver.get('decade_start', None)
DECADE_END      = _driver.get('decade_end', None)
N_RESERVOIRS    = _driver.get('n_reservoirs', 3)
ENFORCE_WB      = _driver.get('enforce_water_balance', 'water-year')
MODULES         = _cfg.get('modules', {})
RESERVOIR_ORDER = _driver.get('reservoir_order',
                               ['shallow', 'soil', 'karst'][:N_RESERVOIRS])
_FIXED_REC      = _driver.get('recession_exponents', None)

_MODULE_PARAMS = {
    'snowpack':      ['PDD_melt_factor'],
    'frozen_ground': ['log__fdd_threshold', 'snow_insulation_k'],
    'direct_runoff': ['f_direct_runoff'],
    'rain_on_snow':  [],
}
for _mod, _names in _MODULE_PARAMS.items():
    if not MODULES.get(_mod, True):
        for _name in _names:
            if _name in _param_cfg:
                _param_cfg[_name]['active'] = False

with open(CONFIG_TEMPLATE) as _f:
    _model_cfg = yaml.safe_load(_f)
_area_km2 = _model_cfg['catchment']['drainage_basin_area__km2']
_datafile  = _model_cfg['timeseries']['datafile']
_df        = pd.read_csv(_datafile, parse_dates=['Date'])
_Q_spec    = _df['Discharge [m^3/s]'].values * 86400.0 / (_area_km2 * 1e3)
_precip    = _df['Precipitation [mm/day]'].values
_hs        = HydrographSeparation(_Q_spec, n_reservoirs=len(RESERVOIR_ORDER),
                                   precip=_precip)
_hs.fit()
INITIAL_STATES = {'reservoirs': _hs.get_initial_conditions()['H0']}

PENALTY = 2.0

params, results = di.read_parameters_file()


def get(name):
    p = _param_cfg[name]
    return params[name] if p['active'] else p['fixed']


def _is_active(name):
    return _param_cfg.get(name, {}).get('active', False)


def _pdm_list():
    names = [f'log__H0_pdm_{l}' for l in RESERVOIR_ORDER]
    vals  = [10 ** get(n) if n in _param_cfg and _param_cfg[n]['active'] else None
             for n in names]
    return vals if any(v is not None for v in vals) else None


def _tile_list():
    names = [f'f_tile_{l}' for l in RESERVOIR_ORDER]
    if not any(n in _param_cfg for n in names):
        return None
    vals = [get(n) if n in _param_cfg else 0.0 for n in names]
    # Legacy shared-tile for W/V/N/M series: soil tile fraction duplicated to shallow
    if ('shallow' in RESERVOIR_ORDER
            and 'f_tile_shallow' not in _param_cfg
            and 'f_tile_soil' in _param_cfg):
        idx_sh = RESERVOIR_ORDER.index('shallow')
        idx_so = RESERVOIR_ORDER.index('soil')
        vals[idx_sh] = vals[idx_so]
    return vals


def _tau_tile():
    if 'log__tau_tile' not in _param_cfg:
        return None
    return 10 ** get('log__tau_tile')


def _hmax():
    key = f'log__Hmax_{RESERVOIR_ORDER[0]}'
    if key not in _param_cfg:
        return None
    return [10 ** get(key)]


def _et_alpha():
    if 'et_alpha' not in _param_cfg:
        return None
    p = _param_cfg['et_alpha']
    return params['et_alpha'] if p['active'] else float(p['fixed'])


def _et_scale():
    if 'et_scale' not in _param_cfg:
        return None
    p = _param_cfg['et_scale']
    return params['et_scale'] if p['active'] else float(p['fixed'])


def _recession_exponents():
    if _FIXED_REC is not None:
        return list(_FIXED_REC), 0
    has_shared = 'recession_b' in _param_cfg
    exponents  = []
    n_cal      = 0
    for label in RESERVOIR_ORDER:
        label_key = f'recession_b_{label}'
        if label == 'shallow':
            exponents.append(1.0)
        elif has_shared:
            exponents.append(get('recession_b'))
        elif label_key in _param_cfg:
            exponents.append(get(label_key))
            n_cal += 1
        else:
            exponents.append(1.0)
    if has_shared:
        n_cal = 1 if any(l != 'shallow' for l in RESERVOIR_ORDER) else 0
    if all(e == 1.0 for e in exponents):
        return None, 0
    return exponents, n_cal


try:
    rec_exp, rec_k = _recession_exponents()
    result = run_and_score(
        CONFIG_TEMPLATE,
        t_efold                = [10 ** get(f'log__t_efold_{l}') for l in RESERVOIR_ORDER],
        f_to_discharge         = [get(f'f_exfiltration_{l}') for l in RESERVOIR_ORDER[:-1]],
        melt_factor            =  get('PDD_melt_factor'),
        fdd_threshold          =  10 ** get('log__fdd_threshold'),
        snow_insulation_k      =  get('snow_insulation_k'),
        Hmax                   =  _hmax(),
        pdm_H0                 =  _pdm_list(),
        f_tile                 =  _tile_list(),
        tau_tile               =  _tau_tile(),
        direct_runoff_fraction =  get('f_direct_runoff'),
        baseflow_Q             =  get('baseflow_Q'),
        et_scale               =  _et_scale(),
        et_alpha               =  _et_alpha(),
        recession_exponents            = rec_exp,
        recession_exponents_calibrated = rec_k,
        modules                =  MODULES,
        routing_K              =  10 ** get('log__routing_K'),
        routing_N              =  ROUTING_N,
        enforce_water_balance  =  ENFORCE_WB,
        initial_states         =  INITIAL_STATES,
        start                  =  DECADE_START,
        end                    =  DECADE_END,
        spin_up_cycles         =  SPIN_UP_CYCLES,
        metric                 =  METRIC,
    )
    neg_score = 1.0 - result.score if np.isfinite(result.score) else PENALTY

except Exception:
    neg_score = PENALTY

results['neg_kge'].function = neg_score
results.write()
