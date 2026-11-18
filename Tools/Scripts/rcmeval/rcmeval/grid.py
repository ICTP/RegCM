#
#    This file is part of ICTP RegCM.
#
#    Use of this source code is governed by an MIT-style license that can
#    be found in the LICENSE file or at
#
#         https://opensource.org/licenses/MIT.
#
"""Grids and the remapping of a field from the finer to the coarser of two
grids.

The remapping aggregates: every valid cell of the fine grid is assigned to
the coarse cell whose centre is nearest (KD-tree on unit vectors), and a
coarse value is the area weighted mean of its fine cells. A coarse cell is
valid when its fine cells cover at least half of its area; coarse cells
that get no fine cell (similar resolutions) take the value of the nearest
fine cell within reach. This approximates the conservative remapping of
fine fields onto a coarser grid, for any curvilinear grid, with scipy only.
The sparse weights are computed once for a pair of grids and cached.
"""

import hashlib
import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree

EARTH_RADIUS = 6371000.0


def _xyz(lat, lon):
    la = np.radians(lat)
    lo = np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo),
                     np.sin(la)], axis=-1)


class Grid:
    """Cell centres (2D lat, lon), cell areas (m2), valid cells, and the
    land fraction if known."""

    def __init__(self, lat, lon, area, valid=None, land=None, name=''):
        self.lat = np.asarray(lat, float)
        self.lon = ((np.asarray(lon, float) + 180.0) % 360.0) - 180.0
        if self.lat.ndim == 1:
            self.lon, self.lat = np.meshgrid(self.lon, self.lat)
        self.area = np.asarray(area, float)
        self.valid = (np.ones(self.lat.shape, bool) if valid is None
                      else np.asarray(valid, bool))
        self.land = None if land is None else np.asarray(land, float)
        self.name = name

    @property
    def shape(self):
        return self.lat.shape

    @property
    def resolution(self):
        """Typical cell size (km)."""
        return float(np.sqrt(np.nanmedian(self.area[self.valid]))) / 1000.0

    def bbox(self, margin=0.0):
        la = self.lat[self.valid]
        lo = self.lon[self.valid]
        return (la.min() - margin, la.max() + margin,
                lo.min() - margin, lo.max() + margin)

    def fingerprint(self):
        h = hashlib.sha1()
        for a in (self.lat, self.lon, self.valid):
            h.update(np.ascontiguousarray(a).tobytes())
        return h.hexdigest()[:16]


def regular(lat1d, lon1d, name=''):
    """Grid of a regular latitude-longitude dataset."""
    lat1d = np.asarray(lat1d, float)
    lon1d = np.asarray(lon1d, float)

    def edges(c):
        e = np.empty(len(c) + 1)
        e[1:-1] = 0.5 * (c[1:] + c[:-1])
        e[0] = c[0] - 0.5 * (c[1] - c[0])
        e[-1] = c[-1] + 0.5 * (c[-1] - c[-2])
        return e
    ely = np.clip(edges(lat1d), -90.0, 90.0)
    elx = edges(lon1d)
    dsin = np.abs(np.diff(np.sin(np.radians(ely))))
    dlon = np.abs(np.radians(np.diff(elx)))
    area = EARTH_RADIUS**2 * np.outer(dsin, dlon)
    return Grid(lat1d, lon1d, area, name=name)


def remap_weights(src, dst):
    """Sparse matrix W (dst cells x src cells) of source cell areas, and the
    fraction of each destination cell covered by valid source cells."""
    nsrc = src.lat.size
    ndst = dst.lat.size
    dvalid = np.flatnonzero(dst.valid.ravel())
    svalid = np.flatnonzero(src.valid.ravel())
    dxyz = _xyz(dst.lat.ravel()[dvalid], dst.lon.ravel()[dvalid])
    sxyz = _xyz(src.lat.ravel()[svalid], src.lon.ravel()[svalid])
    dtree = cKDTree(dxyz)
    dist, idx = dtree.query(sxyz)
    # reach: about half the diagonal of the destination cell
    dsize = np.sqrt(dst.area.ravel()[dvalid]) / EARTH_RADIUS
    ok = dist < 0.75 * dsize[idx]
    rows = dvalid[idx[ok]]
    cols = svalid[ok]
    sarea = src.area.ravel()[cols]
    w = sp.csr_matrix((sarea, (rows, cols)), shape=(ndst, nsrc))
    cover = np.asarray(w.sum(axis=1)).ravel() / dst.area.ravel()
    # destination cells without source cells: nearest source cell within
    # reach (similar resolutions)
    empty = dvalid[np.asarray(w.sum(axis=1)).ravel()[dvalid] == 0.0]
    if empty.size:
        stree = cKDTree(sxyz)
        d2, i2 = stree.query(_xyz(dst.lat.ravel()[empty],
                                  dst.lon.ravel()[empty]))
        ok2 = d2 < 0.75 * np.sqrt(dst.area.ravel()[empty]) / EARTH_RADIUS
        if ok2.any():
            extra = sp.csr_matrix((np.ones(ok2.sum()),
                                   (empty[ok2], svalid[i2[ok2]])),
                                  shape=(ndst, nsrc))
            w = w + extra
            cover[empty[ok2]] = 1.0
    return w.tocsr(), cover.reshape(dst.shape)


def apply(w, field, cover=None, mincover=0.5):
    """Remap field (..., ny, nx) of the source grid with the weights w;
    missing values are left out of the means."""
    a = np.asarray(field, float)
    lead = a.shape[:-2]
    flat = a.reshape((-1, a.shape[-2] * a.shape[-1])).T
    good = np.isfinite(flat)
    num = w @ np.where(good, flat, 0.0)
    den = w @ good.astype(float)
    with np.errstate(invalid='ignore', divide='ignore'):
        out = np.where(den > 0.0, num / den, np.nan)
    nshape = w.shape[0]
    out = out.T.reshape(lead + (-1,))
    if cover is not None:
        out[..., cover.ravel() < mincover] = np.nan
    return out.reshape(lead + cover.shape if cover is not None
                       else lead + (nshape,))


class Remapper:
    """Remapping between a model grid and a reference grid onto the coarser
    of the two, with the weights cached on disk."""

    def __init__(self, model, ref, cache=None, mincover=0.5):
        self.model = model
        self.ref = ref
        self.mincover = mincover
        if model.resolution >= ref.resolution:
            self.common, fine, self.fine_is = model, ref, 'ref'
        else:
            self.common, fine, self.fine_is = ref, model, 'model'
        key = None
        if cache is not None:
            key = cache.key(['remap', fine.fingerprint(),
                             self.common.fingerprint()])
            got = cache.get_npz(key)
        else:
            got = None
        if got is None:
            self.w, self.cover = remap_weights(fine, self.common)
            if cache is not None:
                cache.put_npz(key, w=self.w, cover=self.cover)
        else:
            self.w, self.cover = got['w'], got['cover']
        valid = self.common.valid & (self.cover >= mincover)
        if self.fine_is == 'model' and self.common.land is None:
            # land fraction of the reference cells from the model mask
            land = apply(self.w, model.land if model.land is not None
                         else np.zeros(model.shape), self.cover, mincover)
        else:
            land = self.common.land
        self.valid = valid
        self.land = land

    def to_common(self, field, which):
        """Field of the model or of the reference on the common grid."""
        if which == self.fine_is:
            return apply(self.w, field, self.cover, self.mincover)
        out = np.array(field, float)
        out[..., ~self.valid] = np.nan
        return out
