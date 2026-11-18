#
#    This file is part of ICTP RegCM.
#
#    Use of this source code is governed by an MIT-style license that can
#    be found in the LICENSE file or at
#
#         https://opensource.org/licenses/MIT.
#
"""Figures: bias maps, region profiles and annual cycles."""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt          # noqa: E402
import matplotlib.colors as mcolors     # noqa: E402

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAVE_CARTOPY = True
except ImportError:
    HAVE_CARTOPY = False

# Default bias levels and colour maps (dry/cold to wet/warm)
BIAS_LEVELS = {
    'pr': [-6, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 6],
    'evspsbl': [-3, -2, -1, -0.5, -0.25, 0.25, 0.5, 1, 2, 3],
    'tas': [-8, -6, -4, -2, -1, 1, 2, 4, 6, 8],
    'clt': [-50, -40, -30, -20, -10, 10, 20, 30, 40, 50],
}
BIAS_CMAP = {'pr': 'BrBG', 'evspsbl': 'BrBG', 'clt': 'PuOr_r',
             'cll': 'PuOr_r', 'clm': 'PuOr_r', 'clh': 'PuOr_r',
             'prw': 'BrBG'}


def _levels(var, data, nlev=10):
    if var in BIAS_LEVELS:
        return BIAS_LEVELS[var]
    a = np.nanpercentile(np.abs(np.concatenate(
        [np.ravel(d) for d in data])), 95)
    if not np.isfinite(a) or a == 0:
        a = 1.0
    e = 10 ** np.floor(np.log10(a))
    top = np.ceil(a / e) * e
    return list(np.linspace(-top, top, nlev + 1))


def bias_maps(var, run, panels, path, title=''):
    """panels: rows of (row label, [(col label, lat, lon, bias, mean)])."""
    nr = len(panels)
    nc = max(len(r[1]) for r in panels)
    levels = _levels(var, [p[3] for r in panels for p in r[1]])
    cmap = plt.get_cmap(BIAS_CMAP.get(var, 'RdBu_r'))
    norm = mcolors.BoundaryNorm(levels, cmap.N, extend='both')
    lat0 = panels[0][1][0][1]
    lon0 = panels[0][1][0][2]
    asp = (np.nanmax(lat0) - np.nanmin(lat0)) / max(
        np.nanmax(lon0) - np.nanmin(lon0), 1.0)
    kw = {'projection': ccrs.PlateCarree()} if HAVE_CARTOPY else {}
    fig, axs = plt.subplots(nr, nc, figsize=(3.2 * nc + 1.0,
                                             3.2 * nr * max(asp, 0.4) + 0.6),
                            subplot_kw=kw, squeeze=False)
    im = None
    for i, (rlab, cols) in enumerate(panels):
        for j in range(nc):
            ax = axs[i, j]
            if j >= len(cols):
                ax.set_visible(False)
                continue
            clab, lat, lon, bias, mean = cols[j]
            tk = {'transform': ccrs.PlateCarree()} if HAVE_CARTOPY else {}
            im = ax.pcolormesh(lon, lat, np.ma.masked_invalid(bias),
                               cmap=cmap, norm=norm, shading='auto', **tk)
            if HAVE_CARTOPY:
                ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
                ax.add_feature(cfeature.BORDERS, linewidth=0.3,
                               edgecolor='grey')
                ax.set_extent([np.nanmin(lon), np.nanmax(lon),
                               np.nanmin(lat), np.nanmax(lat)],
                              crs=ccrs.PlateCarree())
            ax.set_title(f'{rlab} {run} - {clab}\nmean bias {mean:+.2f}',
                         fontsize=8)
    cb = fig.colorbar(im, ax=axs, shrink=0.8, extend='both')
    cb.set_label(title, fontsize=9)
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


def profiles(region, var, units, plev, lines, path):
    """lines: [(label, values, style)] with values on plev (hPa)."""
    fig, ax = plt.subplots(figsize=(3.6, 4.6))
    for lab, v, st in lines:
        ax.plot(v, plev, st, label=lab, lw=1.6)
    ax.set_ylim(max(plev), min(plev))
    ax.set_ylabel('hPa')
    ax.set_xlabel(f'{var} [{units}]')
    ax.set_title(region, fontsize=9)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


def annual_cycle(region, var, units, lines, path):
    """lines: [(label, months(1-12), values, style)]."""
    fig, ax = plt.subplots(figsize=(5.0, 3.4))
    for lab, m, v, st in lines:
        ax.plot(m, v, st, label=lab, lw=1.6, marker='o', ms=3)
    ax.set_xticks(range(1, 13))
    ax.set_xticklabels([c for c in 'JFMAMJJASOND'])
    ax.set_ylabel(f'{var} [{units}]')
    ax.set_title(region, fontsize=9)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


def outname(outdir, *parts):
    os.makedirs(outdir, exist_ok=True)
    return os.path.join(outdir, '_'.join(str(p).replace(' ', '-')
                                         for p in parts) + '.png')
