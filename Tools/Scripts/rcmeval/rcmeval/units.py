#
#    This file is part of ICTP RegCM.
#
#    Use of this source code is governed by an MIT-style license that can
#    be found in the LICENSE file or at
#
#         https://opensource.org/licenses/MIT.
#
"""Units of the evaluation: every field, model or reference, is brought to
the target unit of its variable before any statistic is computed."""

import numpy as np

# Target units by variable (CMIP names)
TARGET = {
    'pr': 'mm/day', 'prc': 'mm/day', 'evspsbl': 'mm/day',
    'tas': 'degC', 'tasmax': 'degC', 'tasmin': 'degC', 'ts': 'degC',
    'clt': '%', 'cll': '%', 'clm': '%', 'clh': '%', 'hurs': '%',
    'hur': '%', 'cl': '%',
    'clwvi': 'g/m2', 'clivi': 'g/m2', 'prw': 'kg/m2',
    'hfls': 'W/m2', 'hfss': 'W/m2', 'rsut': 'W/m2', 'rlut': 'W/m2',
    'rsdt': 'W/m2', 'rsds': 'W/m2', 'rlds': 'W/m2', 'rsus': 'W/m2',
    'rlus': 'W/m2', 'rsns': 'W/m2', 'rsnl': 'W/m2', 'rsutcs': 'W/m2',
    'rlutcs': 'W/m2',
    'sfcWind': 'm/s', 'uas': 'm/s', 'vas': 'm/s', 'ua': 'm/s', 'va': 'm/s',
    'psl': 'hPa', 'ps': 'hPa',
    'ta': 'K', 'hus': 'g/kg', 'wap': 'Pa/s', 'zg': 'm',
}

# Canonical spelling of the units found in the files
_ALIAS = {
    'kg m-2 s-1': 'kg/m2/s', 'kg m**-2 s**-1': 'kg/m2/s',
    'kg/m^2/s': 'kg/m2/s', 'kg m^-2 s^-1': 'kg/m2/s',
    'mm/day': 'mm/day', 'mm d-1': 'mm/day', 'mm day-1': 'mm/day',
    'mm/d': 'mm/day',
    'k': 'K', 'kelvin': 'K', 'degc': 'degC', 'celsius': 'degC',
    'degrees_celsius': 'degC', 'deg_c': 'degC', 'c': 'degC',
    '%': '%', 'percent': '%', '1': '1', '(0 - 1)': '1', 'fraction': '1',
    'kg m-2': 'kg/m2', 'kg m**-2': 'kg/m2', 'kg/m2': 'kg/m2',
    'g/m2': 'g/m2', 'g m-2': 'g/m2',
    'w m-2': 'W/m2', 'w m**-2': 'W/m2', 'w/m2': 'W/m2', 'w/m^2': 'W/m2',
    'm s-1': 'm/s', 'm/s': 'm/s', 'm s**-1': 'm/s',
    'pa': 'Pa', 'hpa': 'hPa', 'mb': 'hPa',
    'kg kg-1': 'kg/kg', 'kg kg**-1': 'kg/kg', 'kg/kg': 'kg/kg',
    'g/kg': 'g/kg', 'pa s-1': 'Pa/s', 'pa s**-1': 'Pa/s', 'pa/s': 'Pa/s',
    'm': 'm',
}

# (from, to): (factor, offset), value_to = value_from*factor + offset
_CONV = {
    ('kg/m2/s', 'mm/day'): (86400.0, 0.0),
    ('K', 'degC'): (1.0, -273.15),
    ('degC', 'K'): (1.0, 273.15),
    ('1', '%'): (100.0, 0.0),
    ('kg/m2', 'g/m2'): (1000.0, 0.0),
    ('Pa', 'hPa'): (0.01, 0.0),
    ('kg/kg', 'g/kg'): (1000.0, 0.0),
}


def canonical(units):
    u = (units or '').strip()
    return _ALIAS.get(u.lower(), _ALIAS.get(u, u))


def conversion(units, var):
    """Factor and offset from the file units to the target unit of var."""
    tgt = TARGET.get(var)
    src = canonical(units)
    if tgt is None or src == tgt:
        return 1.0, 0.0
    if (src, tgt) in _CONV:
        return _CONV[(src, tgt)]
    raise ValueError(f'No conversion of {var} from "{units}" to "{tgt}"')


def to_target(da, var, units=None):
    """DataArray in the target unit of var (units from attrs if not given)."""
    f, o = conversion(units if units is not None
                      else da.attrs.get('units', ''), var)
    if f != 1.0 or o != 0.0:
        attrs = dict(da.attrs)
        da = da * f + o
        da.attrs = attrs
    da.attrs['units'] = TARGET.get(var, da.attrs.get('units', ''))
    return da


def qsat(t, p):
    """Saturation specific humidity over liquid water (kg/kg), t in K,
    p in Pa."""
    es = 611.2 * np.exp(17.67 * (t - 273.15) / (t - 29.65))
    return 0.622 * es / (p - 0.378 * es)
