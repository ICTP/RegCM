#
#    This file is part of ICTP RegCM.
#
#    Use of this source code is governed by an MIT-style license that can
#    be found in the LICENSE file or at
#
#         https://opensource.org/licenses/MIT.
#
"""Periods, seasons, regions and the statistics of the comparison."""

import calendar
import numpy as np

_MON = 'JFMAMJJASOND'
_NAMES = [m.upper() for m in calendar.month_abbr[1:]]


def season_months(season):
    """Months (1-12) of a season: ANN, a month name (SEP), or consecutive
    month initials (DJF, JJAS, ONDJF)."""
    s = season.upper()
    if s == 'ANN':
        return list(range(1, 13))
    if s in _NAMES:
        return [_NAMES.index(s) + 1]
    two = _MON + _MON
    i = two.find(s)
    if i < 0 or len(s) > 12:
        raise ValueError(f'Unknown season {season}')
    return [(i + k) % 12 + 1 for k in range(len(s))]


def parse_years(spec):
    """'1990-1999', [1990, 1991], 1990 or '1990:1999:2'."""
    if isinstance(spec, int):
        return [spec]
    if isinstance(spec, (list, tuple)):
        return [int(y) for y in spec]
    s = str(spec)
    for sep in ('-', ':'):
        if sep in s:
            p = [int(x) for x in s.split(sep)]
            return list(range(p[0], p[1] + 1, p[2] if len(p) > 2 else 1))
    return [int(s)]


def season_mean(monthly, months):
    """Mean over the given calendar months of a monthly series (time, ...),
    weighted by the length of the months; None without data."""
    t = monthly['time']
    sel = monthly.sel(time=t['time.month'].isin(months))
    if sel.sizes['time'] == 0:
        return None
    w = sel['time'].dt.days_in_month.astype(float)
    return (sel * w).sum('time') / w.sum()


def months_present(monthly):
    return sorted(set(int(m) for m in monthly['time.month'].values))


def region_mask(lat, lon, land, valid, region):
    """Cells of a region: box [lat0, lat1, lon0, lon1] and surface
    (all, land, ocean) from the land fraction."""
    la0, la1, lo0, lo1 = region.get('box', [-90, 90, -180, 180])
    m = valid & (lat >= la0) & (lat <= la1) & (lon >= lo0) & (lon <= lo1)
    surf = region.get('surface', 'all')
    if surf != 'all' and land is not None:
        m &= (land >= 0.5) if surf == 'land' else (land < 0.5)
    return m


# Fewest cells of a region for a pattern correlation: with fewer the
# correlation is noise and is not given.
MINCORR = 10


def compare(model, ref, area, mask):
    """Area weighted means, bias, RMSE and pattern correlation of two fields
    on the same grid over the cells of mask."""
    ok = mask & np.isfinite(model) & np.isfinite(ref)
    if ok.sum() < 1:
        return None
    w = area[ok] / area[ok].sum()
    a = model[ok]
    b = ref[ok]
    ma = float((w * a).sum())
    mb = float((w * b).sum())
    rmse = float(np.sqrt((w * (a - b) ** 2).sum()))
    if ok.sum() >= MINCORR:
        cov = (w * (a - ma) * (b - mb)).sum()
        sa = np.sqrt((w * (a - ma) ** 2).sum())
        sb = np.sqrt((w * (b - mb) ** 2).sum())
        corr = float(cov / (sa * sb)) if sa > 0 and sb > 0 else np.nan
    else:
        corr = np.nan
    return {'model': ma, 'refmean': mb, 'bias': ma - mb, 'rmse': rmse,
            'corr': corr, 'cells': int(ok.sum())}


def region_mean(field, area, mask):
    ok = mask & np.isfinite(field)
    if ok.sum() < 1:
        return np.nan
    return float((field[ok] * area[ok]).sum() / area[ok].sum())
