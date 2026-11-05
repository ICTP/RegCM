"""
Archive-size CORDEX-CMIP6 files assembled from monthly segments.

Each RegCM output file (one month) produces, for every variable and
frequency, a segment: a complete CORDEX file with the records of that
month, saved in the hidden directory ".segments" of the variable directory.
When the segments cover a whole chunk of the archive, they are merged into
the chunk file, which gets its final name only when complete; the segments
are then removed. The chunks follow the CORDEX-CMIP6 archiving
specifications:

    sub-daily (1hr, 3hr, 6hr): one year
    day:  five years, starting with the years ending in 1 or 6
    mon:  ten years, starting with the years ending in 1

clipped to the first and last year of the experiment when they are given.

No file is ever modified in place: a chunk is written under a temporary
name and renamed, so an interrupted run leaves the segments and the previous
files intact and a rerun assembles them again. Records found in more than
one source are taken once: segments come before an existing chunk file, so
a month processed again replaces its old records. A chunk with missing
records is not assembled.
"""

import glob
import logging
import os
import re
import time
import uuid

import numpy as np
from netCDF4 import Dataset, num2date
from cftime import datetime as cfdatetime

LOGGER = logging.getLogger(__name__)

SEGMENT_DIR = '.segments'
# A lock older than this (seconds) was left by an interrupted run
LOCK_MAX_AGE = 7200

_SUBDAILY = re.compile(r'^(\d+)hr(Pt)?$')
_RANGE = re.compile(r'^(.*)_(\d{6,12})-(\d{6,12})\.nc$')


def chunked_frequency(freq):
    """True if the files of this frequency are assembled in chunks"""
    return freq in ('day', 'mon') or _SUBDAILY.match(freq) is not None


def time_step(freq):
    """The time step of a frequency: hours, or 'month'"""
    if freq == 'mon':
        return 'month'
    if freq == 'day':
        return 24.0
    match = _SUBDAILY.match(freq)
    if match is None:
        raise ValueError('No time step for frequency {}'.format(freq))
    return float(match.group(1))


def chunk_years(year, freq, first_year=None, last_year=None):
    """The first and last year of the chunk that contains year"""
    if freq == 'mon':
        start = year - (year - 1) % 10
        end = start + 9
    elif freq == 'day':
        start = year - (year - 1) % 5
        end = start + 4
    else:
        start = end = year
    if first_year is not None:
        start = max(start, first_year)
    if last_year is not None:
        end = min(end, last_year)
    return start, end


def in_experiment(years, first_year=None, last_year=None):
    """Which of the years are in the experiment"""
    years = np.asarray(years)
    keep = np.ones(years.shape, dtype=bool)
    if first_year is not None:
        keep &= years >= first_year
    if last_year is not None:
        keep &= years <= last_year
    return keep


def segment_path(final_path):
    """The segment path of a file that would be saved as final_path"""
    seg_dir = os.path.join(os.path.dirname(final_path), SEGMENT_DIR)
    os.makedirs(seg_dir, exist_ok=True)
    return os.path.join(seg_dir, os.path.basename(final_path))


def temporary_path(path):
    """A hidden temporary name, in the same directory, for path"""
    return os.path.join(os.path.dirname(path),
                        '.' + os.path.basename(path) + '.tmp')


def _file_prefix(basename):
    """The file name without its time range, or None"""
    match = _RANGE.match(basename)
    return None if match is None else match.group(1)


def _read_times(path):
    """Time values (days), units and calendar of a CORDEX file"""
    with Dataset(path) as ncf:
        tvar = ncf.variables['time']
        times = np.asarray(tvar[:], dtype=np.float64)
        return times, tvar.units, getattr(tvar, 'calendar', 'standard')


def _date_string(date, freq):
    if freq == 'mon':
        return '{0.year:04d}{0.month:02d}'.format(date)
    if freq == 'day':
        return '{0.year:04d}{0.month:02d}{0.day:02d}'.format(date)
    return '{0.year:04d}{0.month:02d}{0.day:02d}{0.hour:02d}{0.minute:02d}' \
        .format(date)


def check_records(dates, freq, start, end, complete=True, late_start=0.0):
    """
    Check that the dates form a continuous series of the frequency from the
    start of year start to the end of year end. With complete=False the
    series may begin after the start and stop before the end (the edges of
    an experiment). late_start (days) lets a sub-daily or daily series begin
    up to that much after the start of the chunk: the first chunk of an
    experiment has no record at its very first time when the point values
    start one step after the initial time. Returns None if the check passes,
    else the reason.
    """
    if len(dates) == 0:
        return 'no records'
    calendar = dates[0].calendar or 'standard'
    step = time_step(freq)
    if step == 'month':
        index = np.array([d.year * 12 + d.month - 1 for d in dates])
        if np.any(np.diff(index) != 1):
            return 'missing or repeated months'
        first_ok = index[0] == start * 12
        last_ok = index[-1] == end * 12 + 11
    else:
        step_days = step / 24.0
        t0 = cfdatetime(start, 1, 1, calendar=calendar)
        t1 = cfdatetime(end + 1, 1, 1, calendar=calendar)
        days = np.array([(d - t0).total_seconds() for d in dates]) / 86400.0
        if np.any(np.abs(np.diff(days) - step_days) > 1.0e-3):
            return 'missing or repeated time steps'
        first_ok = 0.0 <= days[0] < max(step_days, late_start + 1.0e-6)
        span = (t1 - t0).total_seconds() / 86400.0
        last_ok = days[-1] < span <= days[-1] + step_days + 1.0e-6
    if complete and not first_ok:
        return 'the chunk does not start at its first time step'
    if complete and not last_ok:
        return 'the chunk does not end at its last time step'
    return None


class _Source(object):
    """A file with records for a chunk: a segment or a chunk file"""

    def __init__(self, path):
        self.path = path
        self.times, self.units, self.calendar = _read_times(path)
        self.dates = num2date(self.times, self.units, self.calendar)
        self.years = np.array([d.year for d in self.dates])


def _sources(var_dir, prefix):
    """Readable segments of prefix, in time order, then chunk files"""
    segments = []
    pattern = os.path.join(var_dir, SEGMENT_DIR, prefix + '_*.nc')
    for path in sorted(glob.glob(pattern)):
        try:
            segments.append(_Source(path))
        except (OSError, KeyError, ValueError) as err:
            LOGGER.warning('Skipping unreadable segment %s: %s', path, err)
    segments.sort(key=lambda s: s.times[0] if s.times.size else np.inf)
    finals = []
    for path in sorted(glob.glob(os.path.join(var_dir, prefix + '_*.nc'))):
        try:
            finals.append(_Source(path))
        except (OSError, KeyError, ValueError) as err:
            LOGGER.warning('Skipping unreadable file %s: %s', path, err)
    return segments, finals


def _lock(lock_path):
    """Take the lock of a chunk; False if another process holds it"""
    for _ in range(2):
        try:
            fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.write(fd, '{} {}\n'.format(os.uname()[1], os.getpid())
                     .encode())
            os.close(fd)
            return True
        except FileExistsError:
            try:
                age = time.time() - os.path.getmtime(lock_path)
            except OSError:
                continue
            if age < LOCK_MAX_AGE:
                return False
            LOGGER.warning('Removing the stale lock %s', lock_path)
            try:
                os.remove(lock_path)
            except OSError:
                pass
    return False


def _write_chunk(records, out_path):
    """
    Write the records, a list of (source, indices) in time order, in a new
    file at out_path. The variables without the time dimension and all the
    attributes come from the first source; the file gets a new tracking_id
    and creation_date.
    """
    template = records[0][0].path
    tmp_path = temporary_path(out_path)
    with Dataset(template) as src, \
            Dataset(tmp_path, 'w', format=src.data_model) as dst:
        attrs = {a: src.getncattr(a) for a in src.ncattrs()}
        if 'tracking_id' in attrs:
            prefix = attrs['tracking_id'].rsplit('/', 1)[0]
            attrs['tracking_id'] = prefix + '/' + str(uuid.uuid1())
        if 'creation_date' in attrs:
            attrs['creation_date'] = time.strftime('%Y-%m-%dT%H:%M:%SZ',
                                                   time.gmtime())
        dst.setncatts(attrs)
        for name, dim in src.dimensions.items():
            unlimited = dim.isunlimited() or name == 'time'
            dst.createDimension(name, None if unlimited else len(dim))
        timed = []
        for name, var in src.variables.items():
            filters = var.filters() or {}
            chunking = var.chunking()
            fill = var.getncattr('_FillValue') \
                if '_FillValue' in var.ncattrs() else False
            out = dst.createVariable(
                name, var.datatype, var.dimensions,
                zlib=filters.get('zlib', False),
                complevel=filters.get('complevel', 4),
                shuffle=filters.get('shuffle', False),
                fletcher32=filters.get('fletcher32', False),
                chunksizes=chunking if isinstance(chunking, list) else None,
                fill_value=fill,
            )
            # All the attributes but the fill value (set at the creation):
            # the quantization one records the already quantized data
            out.setncatts({a: var.getncattr(a) for a in var.ncattrs()
                           if a != '_FillValue'})
            if 'time' in var.dimensions:
                if var.dimensions[0] != 'time':
                    raise ValueError('Variable {} of {} has time not as its '
                                     'first dimension'.format(name, template))
                timed.append(name)
            else:
                out[...] = var[...]
        start = 0
        for source, indices in records:
            with Dataset(source.path) as ncf:
                for name in timed:
                    data = ncf.variables[name][indices, ...]
                    dst.variables[name][start:start + len(indices), ...] = \
                        data
            start += len(indices)
    os.replace(tmp_path, out_path)


def assemble_chunk(var_dir, prefix, freq, start, end, complete=True,
                   first_year=None, last_year=None):
    """
    Assemble the chunk of years start-end of the files named prefix_* in
    var_dir from its segments (and an existing chunk file). With
    complete=True the chunk is written only if it is complete; otherwise
    (at the end of the experiment) a continuous partial chunk is accepted.
    The records outside the years first_year-last_year of the experiment
    belong to no chunk. Returns the path of the chunk file, or None if it
    was not written.
    """
    seg_dir = os.path.join(var_dir, SEGMENT_DIR)
    lock_path = os.path.join(seg_dir, '.lock_{}_{}-{}'.format(prefix, start,
                                                             end))
    if not _lock(lock_path):
        LOGGER.info('Chunk %s %d-%d is being assembled by another process',
                    prefix, start, end)
        return None
    try:
        segments, finals = _sources(var_dir, prefix)
        in_chunk = []
        for source in segments + finals:
            sel = np.nonzero((source.years >= start) &
                             (source.years <= end))[0]
            if sel.size:
                in_chunk.append((source, sel))
        if not any(s in segments for s, _ in in_chunk):
            return None
        # First source wins: segments in time order, then chunk files
        chosen = {}
        for source, sel in in_chunk:
            for i in sel:
                key = round(float(source.times[i]), 6)
                if key not in chosen:
                    chosen[key] = (source, i)
        keys = sorted(chosen)
        dates = [chosen[k][0].dates[chosen[k][1]] for k in keys]
        late = 1.0 if first_year is not None and start == first_year \
            else 0.0
        problem = check_records(dates, freq, start, end, complete, late)
        if problem is not None:
            if complete:
                LOGGER.debug('Chunk %s %d-%d not assembled: %s', prefix,
                             start, end, problem)
            else:
                LOGGER.error('Chunk %s %d-%d not assembled: %s', prefix,
                             start, end, problem)
            return None
        # Consecutive records of the same source are written together
        records = []
        for k in keys:
            source, i = chosen[k]
            if records and records[-1][0] is source:
                records[-1][1].append(i)
            else:
                records.append((source, [i]))
        out_path = os.path.join(var_dir, '{}_{}-{}.nc'.format(
            prefix, _date_string(dates[0], freq),
            _date_string(dates[-1], freq)))
        LOGGER.info('Assembling %s from %d files', out_path, len(records))
        _write_chunk(records, out_path)
        # The previous chunk files of these years are replaced
        for source in finals:
            if source.path != out_path and source in [s for s, _ in in_chunk]:
                if np.all((source.years >= start) & (source.years <= end)):
                    os.remove(source.path)
        # A segment used here is removed when all its records of the
        # experiment are in chunk files (this one or others already done)
        assembled = set()
        for path in glob.glob(os.path.join(var_dir, prefix + '_*.nc')):
            assembled.update(np.round(_read_times(path)[0], 6).tolist())
        used = [s for s, _ in in_chunk]
        for source in segments:
            if source not in used:
                continue
            mine = in_experiment(source.years, first_year, last_year)
            if all(t in assembled
                   for t in np.round(source.times[mine], 6).tolist()):
                os.remove(source.path)
        return out_path
    finally:
        try:
            os.remove(lock_path)
        except OSError:
            pass


def assemble_after_write(seg_path, freq, first_year=None, last_year=None):
    """
    After a segment has been written, assemble the chunks it touches if
    they are complete.
    """
    var_dir = os.path.dirname(os.path.dirname(seg_path))
    prefix = _file_prefix(os.path.basename(seg_path))
    if prefix is None:
        return []
    years = np.unique(_Source(seg_path).years)
    years = years[in_experiment(years, first_year, last_year)]
    done = []
    for year in years:
        start, end = chunk_years(int(year), freq, first_year, last_year)
        path = assemble_chunk(var_dir, prefix, freq, start, end,
                              complete=True, first_year=first_year,
                              last_year=last_year)
        if path is not None:
            done.append(path)
    return done


def finalize(root_dir, first_year=None, last_year=None, dry_run=False):
    """
    Assemble all the chunks left in the segment directories under root_dir,
    also the partial ones at the edges of the experiment. Returns the list
    of chunk files written (or that would be written with dry_run).
    """
    written = []
    for seg_dir, _, files in os.walk(root_dir):
        if os.path.basename(seg_dir) != SEGMENT_DIR:
            continue
        var_dir = os.path.dirname(seg_dir)
        chunks = set()
        for name in files:
            prefix = _file_prefix(name)
            if prefix is None or name.startswith('.'):
                continue
            freq = prefix.rsplit('_', 1)[-1]
            if not chunked_frequency(freq):
                continue
            try:
                years = np.unique(_Source(os.path.join(seg_dir, name)).years)
            except (OSError, KeyError, ValueError) as err:
                LOGGER.warning('Skipping unreadable segment %s: %s', name,
                               err)
                continue
            for year in years[in_experiment(years, first_year, last_year)]:
                chunks.add((prefix, freq) +
                           chunk_years(int(year), freq, first_year,
                                       last_year))
        for prefix, freq, start, end in sorted(chunks):
            if dry_run:
                LOGGER.info('Would assemble %s %d-%d in %s', prefix, start,
                            end, var_dir)
                written.append((var_dir, prefix, start, end))
                continue
            path = assemble_chunk(var_dir, prefix, freq, start, end,
                                  complete=False, first_year=first_year,
                                  last_year=last_year)
            if path is not None:
                written.append(path)
    return written



def remaining_segments(root_dir):
    """The segment files left under root_dir"""
    left = []
    for seg_dir, _, files in os.walk(root_dir):
        if os.path.basename(seg_dir) != SEGMENT_DIR:
            continue
        left += [os.path.join(seg_dir, f) for f in sorted(files)
                 if f.endswith('.nc') and not f.startswith('.')]
    return left
