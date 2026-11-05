#!/usr/bin/env python3
"""
Tests of the assembly of archive-size chunks from monthly segments
(utilities/chunks.py), on small synthetic CORDEX files.

Run from the pycordexer directory: python3 tests/test_chunks.py
"""

import os
import random
import shutil
import sys
import tempfile
import time
import uuid

import numpy as np
from netCDF4 import Dataset, date2num, num2date
from cftime import datetime as cfdatetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(
    __file__))))
from utilities import chunks

UNITS = 'days since 1950-01-01 00:00:00'
PREFIX = 'tas_DOM_GCM_historical_r1i1p1f1_ICTP_RegCM_v1-r1_{}'
NFAIL = 0


def check(cond, msg):
    global NFAIL
    print(('PASS: ' if cond else 'FAIL: ') + msg)
    if not cond:
        NFAIL += 1


def value(t):
    """The field value of a time (days): a function of time only"""
    return np.float32(np.sin(t / 10.0) * 10.0 + 280.0)


def write_segment(var_dir, freq, dates, calendar='standard', shift=0.0):
    """A CORDEX-like segment with the records at dates"""
    times = np.asarray(date2num(dates, UNITS, calendar), dtype=np.float64)
    if freq == 'mon':
        fmt = '{0.year:04d}{0.month:02d}'
    elif freq == 'day':
        fmt = '{0.year:04d}{0.month:02d}{0.day:02d}'
    else:
        fmt = '{0.year:04d}{0.month:02d}{0.day:02d}{0.hour:02d}{0.minute:02d}'
    name = '{}_{}-{}.nc'.format(PREFIX.format(freq), fmt.format(dates[0]),
                                fmt.format(dates[-1]))
    seg_dir = os.path.join(var_dir, chunks.SEGMENT_DIR)
    os.makedirs(seg_dir, exist_ok=True)
    path = os.path.join(seg_dir, name)
    with Dataset(path, 'w', format='NETCDF4_CLASSIC') as ncf:
        ncf.setncattr('tracking_id', 'hdl:21.14103/' + str(uuid.uuid1()))
        ncf.setncattr('creation_date', '2000-01-01T00:00:00Z')
        ncf.setncattr('frequency', freq)
        ncf.createDimension('time', None)
        ncf.createDimension('bnds', 2)
        ncf.createDimension('y', 3)
        ncf.createDimension('x', 4)
        tvar = ncf.createVariable('time', 'f8', ('time',))
        tvar.units = UNITS
        tvar.calendar = calendar
        tvar.bounds = 'time_bnds'
        tvar[:] = times
        bvar = ncf.createVariable('time_bnds', 'f8', ('time', 'bnds'))
        bvar[:] = np.stack((times - 0.1, times + 0.1), axis=-1)
        lat = ncf.createVariable('lat', 'f8', ('y', 'x'))
        lat[:] = np.arange(12.0).reshape(3, 4)
        tas = ncf.createVariable('tas', 'f4', ('time', 'y', 'x'), zlib=True,
                                 complevel=1, fill_value=np.float32(1.e20))
        tas.units = 'K'
        tas[:] = np.array([np.full((3, 4), value(t) + shift) for t in times])
    return path


def month_dates(year, month, freq, calendar='standard', overlap=False):
    """The record dates of a month for a frequency"""
    t0 = cfdatetime(year, month, 1, calendar=calendar)
    ny, nm = (year + 1, 1) if month == 12 else (year, month + 1)
    t1 = cfdatetime(ny, nm, 1, calendar=calendar)
    if freq == 'mon':
        return [cfdatetime(year, month, 15, calendar=calendar)]
    step = 24.0 if freq == 'day' else float(freq.replace('hr', ''))
    half = 12.0 if freq == 'day' else 0.0
    d0 = date2num(t0, UNITS, calendar) + half / 24.0
    d1 = date2num(t1, UNITS, calendar)
    days = np.arange(d0, d1 - 1.e-9, step / 24.0)
    if overlap:
        days = np.append(days, d1 + half / 24.0)
    return list(num2date(days, UNITS, calendar))


def chunk_files(var_dir):
    return sorted(f for f in os.listdir(var_dir) if f.endswith('.nc'))


def segments(var_dir):
    seg_dir = os.path.join(var_dir, chunks.SEGMENT_DIR)
    if not os.path.isdir(seg_dir):
        return []
    return sorted(f for f in os.listdir(seg_dir) if f.endswith('.nc'))


def read(path):
    with Dataset(path) as ncf:
        return (np.asarray(ncf['time'][:]), np.asarray(ncf['tas'][:]),
                ncf.tracking_id, ncf['time'].calendar)


def test_rules():
    check(chunks.chunk_years(1983, 'day') == (1981, 1985) and
          chunks.chunk_years(1986, 'day') == (1986, 1990) and
          chunks.chunk_years(1990, 'mon') == (1981, 1990) and
          chunks.chunk_years(1991, 'mon') == (1991, 2000) and
          chunks.chunk_years(1984, '6hr') == (1984, 1984) and
          chunks.chunk_years(1983, 'day', first_year=1982) == (1982, 1985),
          'chunk years of the CORDEX-CMIP6 rules')
    check(not chunks.chunked_frequency('fx') and
          chunks.chunked_frequency('1hr') and
          chunks.chunked_frequency('day'), 'chunked frequencies')


def test_year_in_random_order(root):
    var_dir = os.path.join(root, 'a')
    os.makedirs(var_dir)
    months = list(range(1, 13))
    random.Random(7).shuffle(months)
    done = []
    for i, m in enumerate(months):
        seg = write_segment(var_dir, '6hr', month_dates(1981, m, '6hr'))
        done += chunks.assemble_after_write(seg, '6hr')
        if i < 11:
            check(not done, 'no chunk before month {} of 12'.format(i + 1)) \
                if i in (0, 10) else None
    files = chunk_files(var_dir)
    check(files == [PREFIX.format('6hr') +
                    '_198101010000-198112311800.nc'],
          'one year of 6-hourly segments in random order gives the chunk')
    check(segments(var_dir) == [], 'the segments are removed')
    times, data, tid, _ = read(os.path.join(var_dir, files[0]))
    check(times.size == 1460 and np.all(np.diff(times) > 0),
          'the chunk has the 1460 records in order')
    check(np.allclose(data[:, 0, 0], [value(t) for t in times]),
          'the chunk has the values of the segments')
    check(tid.startswith('hdl:21.14103/'), 'the chunk has a tracking_id')
    return var_dir


def test_missing_month(root):
    var_dir = os.path.join(root, 'b')
    os.makedirs(var_dir)
    for m in range(1, 13):
        if m == 7:
            continue
        seg = write_segment(var_dir, 'day', month_dates(1982, m, 'day'))
        chunks.assemble_after_write(seg, 'day', last_year=1982)
    check(chunk_files(var_dir) == [] and len(segments(var_dir)) == 11,
          'a missing month leaves the segments')
    written = chunks.finalize(var_dir, last_year=1982)
    check(written == [] and len(segments(var_dir)) == 11,
          'finalize does not assemble a chunk with a gap')


def test_overlap_and_rerun(root):
    var_dir = os.path.join(root, 'c')
    os.makedirs(var_dir)
    for m in range(1, 13):
        seg = write_segment(var_dir, '6hr',
                            month_dates(1983, m, '6hr', overlap=True))
        chunks.assemble_after_write(seg, '6hr')
    files = chunk_files(var_dir)
    times, _, tid, _ = read(os.path.join(var_dir, files[0]))
    check(len(files) == 1 and times.size == 1460,
          'months overlapping by one record: no repeated record')
    left = segments(var_dir)
    check(len(left) == 1 and '198312' in left[0],
          'the segment with the first record of 1984 is kept')
    os.remove(os.path.join(var_dir, chunks.SEGMENT_DIR, left[0]))
    # A month processed again after the assembly replaces its records
    seg = write_segment(var_dir, '6hr', month_dates(1983, 6, '6hr'),
                        shift=100.0)
    chunks.assemble_after_write(seg, '6hr')
    files2 = chunk_files(var_dir)
    times2, data2, tid2, _ = read(os.path.join(var_dir, files2[0]))
    june = np.array([d.month == 6 for d in num2date(times2, UNITS)])
    check(files2 == files and times2.size == 1460 and
          np.allclose(data2[june, 0, 0],
                      [value(t) + 100.0 for t in times2[june]]) and
          np.allclose(data2[~june, 0, 0],
                      [value(t) for t in times2[~june]]),
          'a month processed again replaces its records in the chunk')
    check(tid2 != tid and segments(var_dir) == [],
          'the rebuilt chunk has a new tracking_id, no segment left')


def test_daily_first_year(root):
    var_dir = os.path.join(root, 'd')
    os.makedirs(var_dir)
    cal = 'noleap'
    for y in (1983, 1984, 1985):
        for m in range(1, 13):
            seg = write_segment(var_dir, 'day',
                                month_dates(y, m, 'day', cal), cal)
            chunks.assemble_after_write(seg, 'day', first_year=1983)
    files = chunk_files(var_dir)
    check(files == [PREFIX.format('day') + '_19830101-19851231.nc'],
          'daily chunk clipped at the first year (noleap calendar)')
    times, _, _, calendar = read(os.path.join(var_dir, files[0]))
    check(times.size == 3 * 365 and calendar == 'noleap',
          'the daily chunk has 3 x 365 records and keeps the calendar')


def test_monthly(root):
    var_dir = os.path.join(root, 'e')
    os.makedirs(var_dir)
    for y in range(1991, 2001):
        for m in range(1, 13):
            seg = write_segment(var_dir, 'mon', month_dates(y, m, 'mon'))
            chunks.assemble_after_write(seg, 'mon')
    files = chunk_files(var_dir)
    times = read(os.path.join(var_dir, files[0]))[0] if files else []
    check(files == [PREFIX.format('mon') + '_199101-200012.nc'] and
          len(times) == 120, 'ten years of monthly means give the chunk')


def test_finalize_partial(root):
    var_dir = os.path.join(root, 'f')
    os.makedirs(var_dir)
    for m in (3, 4, 5):
        seg = write_segment(var_dir, '6hr', month_dates(1990, m, '6hr'))
        chunks.assemble_after_write(seg, '6hr')
    check(chunk_files(var_dir) == [], 'a partial year is not assembled')
    dry = chunks.finalize(root, dry_run=True)
    check(any(d[0] == var_dir for d in dry) and
          len(segments(var_dir)) == 3, 'finalize --dry-run writes nothing')
    chunks.finalize(var_dir)
    check(chunk_files(var_dir) == [PREFIX.format('6hr') +
                                   '_199003010000-199005311800.nc'] and
          segments(var_dir) == [],
          'finalize assembles a partial chunk at the experiment edge')


def point_dates(year, month):
    """6-hourly point values of a RegCM month: from the 1st at 06:00 to the
    1st of the next month at 00:00"""
    t0 = date2num(cfdatetime(year, month, 1), UNITS) + 0.25
    ny, nm = (year + 1, 1) if month == 12 else (year, month + 1)
    t1 = date2num(cfdatetime(ny, nm, 1), UNITS)
    return list(num2date(np.arange(t0, t1 + 1.e-9, 0.25), UNITS))


def test_point_values_edges(root):
    var_dir = os.path.join(root, 'h')
    os.makedirs(var_dir)
    months = [(y, m) for y in (1980, 1981) for m in range(1, 13)]
    order = [(1980, 12)] + [ym for ym in months if ym != (1980, 12)]
    random.Random(3).shuffle(order[1:])
    for i, (y, m) in enumerate(order):
        seg = write_segment(var_dir, '6hr', point_dates(y, m))
        chunks.assemble_after_write(seg, '6hr', first_year=1980,
                                    last_year=1981)
        if i == 0:
            check(chunk_files(var_dir) == [],
                  'one month out of order does not complete the first year')
    files = chunk_files(var_dir)
    check(files == [PREFIX.format('6hr') + '_198001010600-198012311800.nc',
                    PREFIX.format('6hr') + '_198101010000-198112311800.nc'],
          'point values: the first year starts at 06:00, the next at 00:00')
    n = [read(os.path.join(var_dir, f))[0].size for f in files]
    check(n == [366 * 4 - 1, 365 * 4] and segments(var_dir) == [],
          'all records once, the record after the last year dropped')


def test_stale_lock(root):
    var_dir = os.path.join(root, 'g')
    os.makedirs(var_dir)
    segs = [write_segment(var_dir, 'mon', month_dates(2001, m, 'mon'))
            for m in range(1, 13)]
    lock = os.path.join(var_dir, chunks.SEGMENT_DIR,
                        '.lock_{}_2001-2001'.format(PREFIX.format('mon')))
    open(lock, 'w').close()
    check(chunks.assemble_after_write(segs[-1], 'mon', last_year=2001) == [],
          'a recent lock stops the assembly')
    old = time.time() - 2 * chunks.LOCK_MAX_AGE
    os.utime(lock, (old, old))
    done = chunks.assemble_after_write(segs[-1], 'mon', last_year=2001)
    check(len(done) == 1 and not os.path.exists(lock),
          'a stale lock is removed and the chunk assembled')


def main():
    root = tempfile.mkdtemp(prefix='test_chunks_')
    try:
        test_rules()
        test_year_in_random_order(root)
        test_missing_month(root)
        test_overlap_and_rerun(root)
        test_daily_first_year(root)
        test_monthly(root)
        test_finalize_partial(root)
        test_stale_lock(root)
        test_point_values_edges(root)
    finally:
        shutil.rmtree(root)
    print('test_chunks: {}'.format('all passed' if NFAIL == 0 else
                                   '{} failed'.format(NFAIL)))
    return 1 if NFAIL else 0


if __name__ == '__main__':
    sys.exit(main())
