#!/usr/bin/env python3
"""
Dakota driver for decade-by-decade hydroRaVENS calibration.

Run settings and active parameters are read from params.yml.
Returns (1 - score) so Dakota minimisation is equivalent to metric maximisation.
"""

import yaml
import dakota.interfacing as di
import numpy as np
from hydroravens import run_and_score

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
MODULES        = _cfg.get('modules', {})
INITIAL_STATES = None   # set to a CalibResult.final_states dict for chained decades

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


_T_NAMES = ['log__t_efold_shallow', 'log__t_efold_soil', 'log__t_efold_karst']
_F_NAMES = ['f_exfiltration_shallow', 'f_exfiltration_soil']

try:
    result = run_and_score(
        CONFIG_TEMPLATE,
        t_efold               = [10 ** get(n) for n in _T_NAMES[:N_RESERVOIRS]],
        f_to_discharge        = [get(n) for n in _F_NAMES[:N_RESERVOIRS - 1]],
        melt_factor           =  get('PDD_melt_factor'),
        fdd_threshold         =  10 ** get('log__fdd_threshold'),
        snow_insulation_k     =  get('snow_insulation_k'),
        Hmax                  = [10 ** get('log__Hmax_shallow')],
        direct_runoff_fraction=  get('f_direct_runoff'),
        baseflow_Q            =  get('baseflow_Q'),
        modules               =  MODULES,
        routing_K             =  10 ** get('log__routing_K'),
        routing_N             =  ROUTING_N,
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
