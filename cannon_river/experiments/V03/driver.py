#!/usr/bin/env python3
"""
Dakota driver for decade-by-decade hydroRaVENS calibration.

Run settings and active parameters are read from params.yml.
Returns (1 - score) so Dakota minimisation is equivalent to metric maximisation.
"""

import yaml
import pandas as pd
import dakota.interfacing as di
import numpy as np
from hydroravens import HydrographSeparation, run_and_score

with open('params.yml') as f:
    _cfg = yaml.safe_load(f)

_driver     = _cfg['driver']
_param_cfg  = _cfg['parameters']

CONFIG_TEMPLATE = _driver['config_template']
METRIC         = _driver['metric']
SPIN_UP_CYCLES = _driver['spin_up_cycles']
ROUTING_N      = _driver['routing_N']
DECADE_START   = _driver['decade_start']
DECADE_END     = _driver['decade_end']
N_RESERVOIRS   = _driver.get('n_reservoirs', 3)
ENFORCE_WB     = _driver.get('enforce_water_balance', 'water-year')
MODULES        = _cfg.get('modules', {})

# Data-driven initial conditions from HydrographSeparation.
# Fitted once per process on the observed record; H0 is fixed across all
# Dakota evaluations so the optimizer sees a consistent starting point.
with open(CONFIG_TEMPLATE) as _f:
    _model_cfg = yaml.safe_load(_f)
_area_km2   = _model_cfg['catchment']['drainage_basin_area__km2']
_datafile   = _model_cfg['timeseries']['datafile']
_df         = pd.read_csv(_datafile, parse_dates=['Date'])
_Q_specific = _df['Discharge [m^3/s]'].values * 86400.0 / (_area_km2 * 1e3)
_precip     = _df['Precipitation [mm/day]'].values
_hs         = HydrographSeparation(_Q_specific, n_reservoirs=N_RESERVOIRS,
                                    precip=_precip)
_hs.fit()
INITIAL_STATES = {'reservoirs': _hs.get_initial_conditions()['H0']}

# Mirror generate_dakota_in.py's module auto-fix so active flags match dakota.in.
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

PENALTY = 2.0   # returned on model failure; safely above any real 1 - score

params, results = di.read_parameters_file()


def get(name):
    """Return the Dakota parameter value if active, else the fixed fallback."""
    p = _param_cfg[name]
    return params[name] if p['active'] else p['fixed']


_T_NAMES    = ['log__t_efold_shallow', 'log__t_efold_soil', 'log__t_efold_karst']
_F_NAMES    = ['f_exfiltration_shallow', 'f_exfiltration_soil']
_PDM_NAMES  = ['log__H0_pdm_shallow', 'log__H0_pdm_soil', 'log__H0_pdm_karst']
_TILE_NAMES = ['f_tile_shallow', 'f_tile_soil', 'f_tile_karst']

def _pdm_list():
    """Return pdm_H0 list if any reservoir has an active PDM parameter, else None."""
    vals = [10 ** get(n) if n in _param_cfg and _param_cfg[n]['active'] else None
            for n in _PDM_NAMES[:N_RESERVOIRS]]
    return vals if any(v is not None for v in vals) else None

def _tile_list():
    """Return f_tile list if any tile parameter exists in params.yml, else None."""
    if not any(n in _param_cfg for n in _TILE_NAMES[:N_RESERVOIRS]):
        return None
    vals = [get(n) if n in _param_cfg else 0.0
            for n in _TILE_NAMES[:N_RESERVOIRS]]
    # Shallow shares f_tile with soil: same tile network, one calibrated parameter.
    if 'f_tile_shallow' not in _param_cfg and 'f_tile_soil' in _param_cfg:
        vals[0] = vals[1]
    return vals

def _tau_tile():
    """Return tile residence time [days] if log__tau_tile exists, else None."""
    if 'log__tau_tile' not in _param_cfg:
        return None
    return 10 ** get('log__tau_tile')

def _et_alpha():
    """Return et_alpha if active (calibrated), else None."""
    if 'et_alpha' not in _param_cfg or not _param_cfg['et_alpha']['active']:
        return None
    return get('et_alpha')

def _wp_soil():
    """Return wp_soil if active, else None (model default 0.0 = no threshold)."""
    if 'wp_soil' not in _param_cfg or not _param_cfg['wp_soil']['active']:
        return None
    return get('wp_soil')

try:
    result = run_and_score(
        CONFIG_TEMPLATE,
        t_efold               = [10 ** get(n) for n in _T_NAMES[:N_RESERVOIRS]],
        f_to_discharge        = [get(n) for n in _F_NAMES[:N_RESERVOIRS - 1]],
        melt_factor           =  get('PDD_melt_factor'),
        fdd_threshold         =  10 ** get('log__fdd_threshold'),
        snow_insulation_k     =  get('snow_insulation_k'),
        Hmax                  = [10 ** get('log__Hmax_shallow')],
        pdm_H0                =  _pdm_list(),
        f_tile                =  _tile_list(),
        tau_tile              =  _tau_tile(),
        direct_runoff_fraction=  get('f_direct_runoff'),
        baseflow_Q            =  get('baseflow_Q'),
        et_alpha              =  _et_alpha(),
        wp_soil               =  _wp_soil(),
        modules               =  MODULES,
        routing_K             =  10 ** get('log__routing_K'),
        routing_N             =  ROUTING_N,
        enforce_water_balance =  ENFORCE_WB,
        initial_states        =  INITIAL_STATES,
        start                 =  DECADE_START,
        end                   =  DECADE_END,
        spin_up_cycles        =  SPIN_UP_CYCLES,
        metric                =  METRIC,
    )
    neg_score = 1.0 - result.score if np.isfinite(result.score) else PENALTY

except Exception:
    neg_score = PENALTY

results['neg_kge'].function = neg_score
results.write()
