#
#    This file is part of ICTP RegCM.
#
#    Use of this source code is governed by an MIT-style license that can
#    be found in the LICENSE file or at
#
#         https://opensource.org/licenses/MIT.
#
"""Content addressed cache of intermediate results.

The key of an entry is a hash of everything the result depends on: the
source files with their size and modification time, the variable, the
period, the operation, the grids and the version of this package. A
changed input gives a new key, so stale entries are never used; old
entries can be removed at any time (rcmeval --clean-cache)."""

import hashlib
import json
import os
import numpy as np
import scipy.sparse as sp
import xarray as xr

from . import __version__


def _code_hash():
    """Hash of the source of this package: a code change invalidates the
    cached results."""
    h = hashlib.sha1(__version__.encode())
    here = os.path.dirname(os.path.abspath(__file__))
    for f in sorted(os.listdir(here)):
        if f.endswith('.py'):
            with open(os.path.join(here, f), 'rb') as fh:
                h.update(fh.read())
    return h.hexdigest()[:12]


CODE = _code_hash()


def files_signature(files):
    """Paths with size and modification time."""
    out = []
    for f in sorted(files):
        st = os.stat(f)
        out.append([os.path.abspath(f), st.st_size, st.st_mtime_ns])
    return out


class Cache:

    def __init__(self, path):
        self.path = os.path.abspath(os.path.expanduser(path))
        os.makedirs(self.path, exist_ok=True)

    def key(self, parts):
        txt = json.dumps([CODE] + list(parts), sort_keys=True,
                         default=str)
        return hashlib.sha1(txt.encode()).hexdigest()[:24]

    def _file(self, key, ext):
        return os.path.join(self.path, key + ext)

    def get(self, key):
        f = self._file(key, '.nc')
        if not os.path.exists(f):
            return None
        with xr.open_dataset(f) as ds:
            return ds.load()

    def put(self, key, ds, label=''):
        f = self._file(key, '.nc')
        tmp = f + '.tmp'
        if isinstance(ds, xr.DataArray):
            ds = ds.to_dataset(name=ds.name or 'data')
        ds = ds.assign_attrs(rcmeval_label=label)
        ds.to_netcdf(tmp)
        os.replace(tmp, f)

    def get_npz(self, key):
        f = self._file(key, '.npz')
        if not os.path.exists(f):
            return None
        z = np.load(f)
        w = sp.csr_matrix((z['data'], z['indices'], z['indptr']),
                          shape=tuple(z['shape']))
        return {'w': w, 'cover': z['cover']}

    def put_npz(self, key, w, cover):
        f = self._file(key, '.npz')
        tmp = f + '.tmp.npz'
        np.savez(tmp, data=w.data, indices=w.indices, indptr=w.indptr,
                 shape=np.array(w.shape), cover=cover)
        os.replace(tmp, f)

    def clean(self):
        n = 0
        for f in os.listdir(self.path):
            if f.endswith(('.nc', '.npz')):
                os.unlink(os.path.join(self.path, f))
                n += 1
        return n
