#
#    This file is part of ICTP RegCM.
#
#    Use of this source code is governed by an MIT-style license that can
#    be found in the LICENSE file or at
#
#         https://opensource.org/licenses/MIT.
#
"""Driver: reads the configuration, computes the statistics of every
variable, region and season for every run against its references, and
writes the figures, a Markdown report and a CSV table."""

import argparse
import csv
import os
import sys
import time
import numpy as np
import xarray as xr
import yaml

from . import __version__
from .cache import Cache
from .catalog import load_catalog, Reference
from .grid import Remapper
from .model import open_run
from .stats import (season_months, parse_years, season_mean, months_present,
                    region_mask, compare, region_mean)
from .units import TARGET
from . import plots

HERE = os.path.dirname(os.path.abspath(__file__))


def log(msg):
    print(time.strftime('%H:%M:%S'), msg, flush=True)


class Evaluation:

    def __init__(self, config, cfgdir):
        self.cfg = config
        cat = config.get('catalog', 'ictp')
        if not os.path.exists(cat):
            for c in (os.path.join(cfgdir, cat),
                      os.path.join(HERE, '..', 'catalog', cat + '.yaml')):
                if os.path.exists(c):
                    cat = c
                    break
        self.catalog = load_catalog(cat)
        self.cache = Cache(config.get('cache', './rcmeval_cache'))
        self.out = config.get('output', './rcmeval_output')
        os.makedirs(self.out, exist_ok=True)
        self.years = parse_years(config['years'])
        self.seasons = config.get('seasons', ['DJF', 'MAM', 'JJA', 'SON'])
        self.regions = config.get('regions', {})
        self.runs = {n: open_run(n, s) for n, s in config['runs'].items()}
        self.rows = []
        self._monthly = {}

    # ---------------------------------------------------------------- data
    def model_monthly(self, run, var):
        k = (run.name, var)
        if k not in self._monthly:
            key = self.cache.key(['model-monthly', var, self.years,
                                  run.sources(var, self.years)])
            got = self.cache.get(key)
            if got is None:
                log(f'  {run.name}: monthly {var}')
                da = run.monthly(var, self.years)
                self.cache.put(key, da.rename(var), f'{run.name} {var}')
            else:
                da = got[var]
            self._monthly[k] = da
        return self._monthly[k]

    def ref_monthly(self, ref, years, bbox, months):
        key = self.cache.key(['ref-monthly', ref.name, ref.var, years,
                              [round(b, 3) for b in bbox], months,
                              ref.vs, ref.sources(years, months)])
        got = self.cache.get(key)
        if got is not None:
            return got[ref.var]
        log(f'  {ref.name}: monthly {ref.var}')
        da = ref.monthly(years, bbox, months)
        self.cache.put(key, da, f'{ref.name} {ref.var}')
        return da

    def _ref(self, entry, var):
        if isinstance(entry, str):
            entry = {'name': entry}
        name = entry['name']
        if name not in self.catalog:
            raise KeyError(f'Dataset {name} not in the catalog')
        ref = Reference(name, self.catalog[name], entry.get('var', var))
        ref.var = var
        yrs = parse_years(entry['years']) if 'years' in entry else self.years
        label = name if 'years' not in entry else \
            f'{name} {yrs[0]}-{yrs[-1]}'
        return ref, yrs, label

    def _regions(self, lat, lon, land, valid):
        out = {'domain': valid}
        for rn, r in self.regions.items():
            out[rn] = region_mask(lat, lon, land, valid, r)
        return out

    # ---------------------------------------------------------- 2D fields
    def surface_variable(self, var, spec):
        refs = spec.get('refs', [])
        for run in self.runs.values():
            if not run.has(var):
                log(f'{run.name} has no {var}')
                continue
            mm = self.model_monthly(run, var)
            mmonths = months_present(mm)
            grid = run.grid()
            bbox = grid.bbox(margin=1.0)
            panels = {s: [] for s in self.seasons}
            for entry in refs:
                try:
                    ref, ryears, label = self._ref(entry, var)
                    rmon = mmonths if ryears == self.years else None
                    rm = self.ref_monthly(ref, ryears, bbox, rmon)
                except (KeyError, FileNotFoundError) as e:
                    log(f'  skip {entry}: {e}')
                    continue
                rem = Remapper(grid, ref.grid(rm), self.cache)
                com = rem.common
                masks = self._regions(com.lat, com.lon, rem.land, rem.valid)
                for s in self.seasons:
                    months = [m for m in season_months(s) if m in mmonths]
                    if not months:
                        continue
                    a = season_mean(mm, months)
                    b = season_mean(rm, months)
                    if a is None or b is None:
                        continue
                    ac = rem.to_common(a.values, 'model')
                    bc = rem.to_common(b.values, 'ref')
                    if ref.surface != 'all':
                        surf = rem.land >= 0.5 if ref.surface == 'land' \
                            else rem.land < 0.5
                        ac = np.where(surf, ac, np.nan)
                        bc = np.where(surf, bc, np.nan)
                    for rn, msk in masks.items():
                        st = compare(ac, bc, com.area, msk)
                        if st:
                            self.rows.append(dict(
                                var=var, stat='mean', season=s, region=rn,
                                run=run.name, ref=label, **st))
                    dom = compare(ac, bc, com.area, masks['domain'])
                    panels[s].append((label, com.lat, com.lon, ac - bc,
                                      dom['bias'] if dom else np.nan))
            if spec.get('maps', True):
                rows = [(s, panels[s]) for s in self.seasons if panels[s]]
                if rows:
                    plots.bias_maps(var, run.name, rows,
                                    plots.outname(self.out, 'bias', var,
                                                  run.name),
                                    f'{var} bias [{TARGET.get(var, "")}]')
            if spec.get('annual_cycle', True) and len(mmonths) >= 3:
                self.cycle(var, run, mm, refs)

    def cycle(self, var, run, mm, refs):
        grid = run.grid()
        masks = self._regions(grid.lat, grid.lon, grid.land, grid.valid)
        for rn, msk in masks.items():
            lines = [(run.name, months_present(mm),
                      [region_mean(season_mean(mm, [m]).values, grid.area,
                                   msk) for m in months_present(mm)], '-')]
            plots.annual_cycle(rn, var, TARGET.get(var, ''), lines,
                               plots.outname(self.out, 'cycle', var,
                                             run.name, rn))

    # --------------------------------------------------------- extremes
    def p99(self, var, refs):
        for run in self.runs.values():
            key = self.cache.key(['model-daily', var, self.years,
                                  run.sources(var, self.years, 'daily')])
            got = self.cache.get(key)
            if got is None:
                log(f'  {run.name}: daily {var}')
                md = run.daily(var, self.years)
                self.cache.put(key, md.rename(var), f'{run.name} {var} day')
            else:
                md = got[var]
            grid = run.grid()
            bbox = grid.bbox(margin=1.0)
            panels = []
            for entry in refs:
                ref, ryears, label = self._ref(entry, var)
                key = self.cache.key(['ref-daily', ref.name, var, ryears,
                                      [round(b, 3) for b in bbox],
                                      ref.sources(ryears)])
                got = self.cache.get(key)
                if got is None:
                    log(f'  {ref.name}: daily {var}')
                    rd = ref.daily(ryears, bbox)
                    self.cache.put(key, rd, f'{ref.name} {var} day')
                else:
                    rd = got[var]
                rem = Remapper(grid, ref.grid(rd), self.cache)
                com = rem.common
                # quantiles of the daily fields on the common grid
                qkey = self.cache.key(['p99', var, self.years, ryears,
                                       run.sources(var, self.years, 'daily'),
                                       ref.sources(ryears), rem.fine_is,
                                       com.fingerprint()])
                got = self.cache.get(qkey)
                if got is None:
                    ac = np.nanpercentile(rem.to_common(md.values, 'model'),
                                          99, axis=0)
                    bc = np.nanpercentile(rem.to_common(rd.values, 'ref'),
                                          99, axis=0)
                    self.cache.put(qkey, xr.Dataset(
                        {'model': (('y', 'x'), ac), 'ref': (('y', 'x'), bc)}),
                        f'{run.name} {label} p99 {var}')
                else:
                    ac = got['model'].values
                    bc = got['ref'].values
                masks = self._regions(com.lat, com.lon, rem.land, rem.valid)
                for rn, msk in masks.items():
                    st = compare(ac, bc, com.area, msk)
                    if st:
                        self.rows.append(dict(var=var, stat='p99',
                                              season='ANN', region=rn,
                                              run=run.name, ref=label, **st))
                dom = compare(ac, bc, com.area, masks['domain'])
                panels.append((label, com.lat, com.lon, ac - bc,
                               dom['bias'] if dom else np.nan))
            if panels:
                plots.bias_maps(var, run.name, [('p99', panels)],
                                plots.outname(self.out, 'p99', var,
                                              run.name),
                                f'{var} p99 bias [{TARGET.get(var, "")}]')

    # --------------------------------------------------------- profiles
    def profiles(self, spec):
        plev = [float(p) for p in spec.get('levels', [
            925, 850, 700, 600, 500, 400, 300, 250, 200])]
        regions = spec.get('regions', list(self.regions))
        for var in spec.get('variables', ['ta', 'hus', 'ua', 'wap']):
            lines = {(rn, s): [] for rn in regions for s in self.seasons}
            for run in self.runs.values():
                key = self.cache.key(['model-profile', var, self.years, plev,
                                      run.sources(var, self.years)])
                got = self.cache.get(key)
                if got is None:
                    log(f'  {run.name}: profiles {var}')
                    mp = run.profile_monthly(var, self.years, plev)
                    self.cache.put(key, mp.rename(var), f'{run.name} {var}')
                else:
                    mp = got[var]
                mmonths = months_present(mp)
                grid = run.grid()
                bbox = grid.bbox(margin=1.0)
                for entry in spec.get('refs', ['ERA5']):
                    ref, ryears, label = self._ref(entry, var)
                    rmon = mmonths if ryears == self.years else None
                    rm = self.ref_monthly(ref, ryears, bbox, rmon)
                    rm = rm.sel(plev=plev, method='nearest')
                    rem = Remapper(grid, ref.grid(rm), self.cache)
                    com = rem.common
                    masks = self._regions(com.lat, com.lon, rem.land,
                                          rem.valid)
                    for s in self.seasons:
                        months = [m for m in season_months(s)
                                  if m in mmonths]
                        if not months:
                            continue
                        a = rem.to_common(season_mean(mp, months).values,
                                          'model')
                        b = rem.to_common(season_mean(rm, months).values,
                                          'ref')
                        for rn in regions:
                            msk = masks[rn]
                            av = [region_mean(a[k], com.area, msk)
                                  for k in range(len(plev))]
                            bv = [region_mean(b[k], com.area, msk)
                                  for k in range(len(plev))]
                            lst = lines[(rn, s)]
                            if not any(lab == label for lab, _, _ in lst):
                                lst.append((label, bv, 'k--'))
                            lst.append((run.name, av, '-'))
                            for k, p in enumerate(plev):
                                self.rows.append(dict(
                                    var=f'{var}{int(p)}', stat='profile',
                                    season=s, region=rn, run=run.name,
                                    ref=label, model=av[k], refmean=bv[k],
                                    bias=av[k] - bv[k], rmse=np.nan,
                                    corr=np.nan, cells=0))
            for (rn, s), ln in lines.items():
                if ln:
                    plots.profiles(f'{rn} {s}', var, TARGET.get(var, ''),
                                   plev, ln, plots.outname(
                                       self.out, 'profile', var, rn, s))

    # ----------------------------------------------------------- report
    def report(self):
        keys = ['var', 'stat', 'season', 'region', 'run', 'ref', 'model',
                'refmean', 'bias', 'rmse', 'corr', 'cells']
        path = os.path.join(self.out, 'stats.csv')
        with open(path, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in self.rows:
                w.writerow({k: r.get(k, '') for k in keys})
        md = [f'# RegCM evaluation (rcmeval {__version__})', '',
              f'Runs: {", ".join(self.runs)}; years '
              f'{self.years[0]}-{self.years[-1]}.', '']
        groups = {}
        for r in self.rows:
            groups.setdefault((r['var'], r['stat'], r['season']),
                              []).append(r)

        def fmt(x):
            return '' if x is None or not np.isfinite(x) else f'{x:.2f}'
        for (var, stat, s), rows in groups.items():
            refs = list(dict.fromkeys(r['ref'] for r in rows))
            runs = list(dict.fromkeys(r['run'] for r in rows))
            regs = list(dict.fromkeys(r['region'] for r in rows))
            base = var.rstrip('0123456789') if stat == 'profile' else var
            md += [f'## {var} {stat} {s} [{TARGET.get(base, "")}]', '',
                   '| region | ' + ' | '.join(refs + runs) + ' |',
                   '|---' * (1 + len(refs) + len(runs)) + '|']
            for rg in regs:
                cells = []
                for rf in refs:
                    v = [r for r in rows if r['region'] == rg and
                         r['ref'] == rf]
                    cells.append(fmt(v[0]['refmean']) if v else '')
                for rn in runs:
                    v = [r for r in rows if r['region'] == rg and
                         r['run'] == rn]
                    cells.append(fmt(v[0]['model']) if v else '')
                md.append(f'| {rg} | ' + ' | '.join(cells) + ' |')
            md.append('')
            if stat != 'profile':
                md += ['| region | run | ref | bias | RMSE | corr | cells |',
                       '|---|---|---|---|---|---|---|']
                for r in rows:
                    md.append(f'| {r["region"]} | {r["run"]} | {r["ref"]} '
                              f'| {r["bias"]:+.2f} | {r["rmse"]:.2f} | '
                              f'{fmt(r["corr"])} | {r["cells"]} |')
                md.append('')
        with open(os.path.join(self.out, 'report.md'), 'w') as f:
            f.write('\n'.join(md) + '\n')
        log(f'Report: {os.path.join(self.out, "report.md")}, {path}')

    def run(self):
        for var, spec in self.cfg.get('variables', {}).items():
            log(f'{var}')
            self.surface_variable(var, spec or {})
        for var, refs in self.cfg.get('p99', {}).items():
            log(f'{var} p99')
            self.p99(var, refs)
        if 'profiles' in self.cfg:
            log('profiles')
            self.profiles(self.cfg['profiles'])
        self.report()


def main(argv=None):
    import warnings
    warnings.filterwarnings('ignore', message='All-NaN slice')
    warnings.filterwarnings('ignore', message='Mean of empty slice')
    ap = argparse.ArgumentParser(
        prog='rcmeval', description='RegCM evaluation against observations')
    ap.add_argument('-c', '--config', required=True, help='YAML configuration')
    ap.add_argument('--clean-cache', action='store_true',
                    help='remove all the cached results first')
    a = ap.parse_args(argv)
    with open(a.config) as f:
        cfg = yaml.safe_load(f)
    ev = Evaluation(cfg, os.path.dirname(os.path.abspath(a.config)))
    if a.clean_cache:
        log(f'removed {ev.cache.clean()} cache entries')
    ev.run()
    return 0


if __name__ == '__main__':
    sys.exit(main())
