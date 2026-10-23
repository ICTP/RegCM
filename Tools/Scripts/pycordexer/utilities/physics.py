"""
Vertical interpolation, geopotential height and CAPE/CIN computations.

This module replaces the former f2py Fortran modules vertint.F90, hgt.F90
and capecin.F90. The functions keep the same argument order of the Fortran
routines, but are vectorized over the horizontal grid. All 3D arrays are
expected with shape (kz, iy, jx), with the vertical index going from the
model top to the bottom, while 2D arrays have shape (iy, jx).
"""

import numpy as np


__copyright__ = 'Copyright (C) 2026 ICTP'
__author__ = 'Graziano Giuliani <ggiulian@ictp.it>'
__credits__ = ["Graziano Giuliani"]


# Standard Gravity (m/sec**2) 3rd CGPM
EGRAV = 9.80665
# Gas constant for dry air in Joules/kg/K
RGAS = 287.05823
# Standard atmosphere ICAO 1993
LRATE = 0.00649
RGLROG = RGAS * LRATE / EGRAV
# Standard surface pressure
P00 = 100000.0
ROVG = 29.2716599
CPOVR = 3.5


def _scalar(x):
    return float(np.ravel(x)[0])


def _take(f, k):
    return np.take_along_axis(f, k[np.newaxis], axis=0)[0]


def _first_above(levels, value):
    """
    Index of the first level strictly greater than value, limited to the
    interval [1, km-1] so that both k and k-1 are valid indexes.
    levels has shape (kz,) or (kz, iy, jx) and must be monotonically
    increasing along the first axis.
    """
    if levels.ndim == 1:
        levels = levels[:, np.newaxis, np.newaxis]
    kx = np.count_nonzero(levels <= value, axis=0)
    return np.clip(kx, 1, levels.shape[0] - 1)


# Vertical interpolation (former mod_vertint)

def _interp_lin(f, lev, target):
    f = np.asarray(f)
    kx = _first_above(lev, target)
    if lev.ndim == 1:
        lk, ln = lev[kx], lev[kx - 1]
    else:
        lk, ln = _take(lev, kx), _take(lev, kx - 1)
    wp = (lk - target) / (lk - ln)
    return (1.0 - wp) * _take(f, kx) + wp * _take(f, kx - 1)


def _interp_log(f, lev, target):
    f = np.asarray(f)
    kx = _first_above(lev, target)
    if lev.ndim == 1:
        lk, ln = lev[kx], lev[kx - 1]
    else:
        lk, ln = _take(lev, kx), _take(lev, kx - 1)
    wp = np.log(lk / target) / np.log(lk / ln)
    return (1.0 - wp) * _take(f, kx) + wp * _take(f, kx - 1)


def intlin_hy(f, pstar, sig, ptop, p):
    """
    Linear in pressure interpolation of f at pressure p (Pa) for the
    hydrostatic core. Above the top and below the bottom level the field is
    considered to have 0 vertical derivative.
    """
    ptp = _scalar(ptop) * 100.0
    sig = np.asarray(sig, dtype=np.float32)
    sigp = (p - ptp) / (pstar - ptp)
    with np.errstate(divide='ignore', invalid='ignore'):
        fp = _interp_lin(f, sig, sigp)
    fp = np.where(sigp <= sig[0], f[0], fp)
    fp = np.where(sigp >= sig[-1], f[-1], fp)
    return fp.astype(np.float32)


def intlog_hy(f, pstar, sig, ptop, p):
    """
    Linear in log(pressure) interpolation of f at pressure p (Pa) for the
    hydrostatic core. Below the bottom level a standard lapse rate is used.
    """
    ptp = _scalar(ptop) * 100.0
    sig = np.asarray(sig, dtype=np.float32)
    sigp = (p - ptp) / (pstar - ptp)
    with np.errstate(divide='ignore', invalid='ignore'):
        fp = _interp_log(f, sig, sigp)
        fx = 0.5 * (f[-1] + f[-2]) * np.exp(RGLROG * np.log(sigp / sig[-1]))
    fp = np.where(sigp <= sig[0], f[0], fp)
    fp = np.where(sigp > sig[-1], fx, fp)
    return fp.astype(np.float32)


def _normalized_pressure(mp, p):
    # Pressure is scaled to [0,1] between the top and the bottom level
    tp = mp[0]
    bp = mp[-1] - tp
    return (mp - tp) / bp, (p - tp) / bp


def intlin_nonhy(f, mp, p):
    """
    Linear in pressure interpolation of f at pressure p (Pa) given the full
    3D pressure field mp (Pa).
    """
    spp, sp = _normalized_pressure(mp, p)
    with np.errstate(divide='ignore', invalid='ignore'):
        fp = _interp_lin(f, spp, sp)
    fp = np.where(sp <= spp[0], f[0], fp)
    fp = np.where(sp >= spp[-1], f[-1], fp)
    return fp.astype(np.float32)


def intlog_nonhy(f, mp, p):
    """
    Linear in log(pressure) interpolation of f at pressure p (Pa) given the
    full 3D pressure field mp (Pa).
    """
    spp, sp = _normalized_pressure(mp, p)
    with np.errstate(divide='ignore', invalid='ignore'):
        fp = _interp_log(f, spp, sp)
        fx = 0.5 * (f[-1] + f[-2]) * np.exp(RGLROG * np.log(sp / spp[-1]))
    fp = np.where(sp <= spp[0], f[0], fp)
    fp = np.where(sp > spp[-1], fx, fp)
    return fp.astype(np.float32)


# Geopotential height (former mod_hgt)

def _height_on_pressure(t, psig, htsig, psfc, zsfc, p):
    """
    Height of the pressure level p given temperature t, pressure psig and
    height htsig on model levels.
    For upward extrapolation, t is considered to have 0 vertical derivative.
    For downward extrapolation, t has a lapse rate of LRATE.
    """
    t = np.asarray(t)
    kt = np.clip(np.count_nonzero(psig < p, axis=0) - 1, 0, psig.shape[0] - 2)
    kb = kt + 1
    pt = _take(psig, kt)
    pb = _take(psig, kb)
    tb = _take(t, kb)
    with np.errstate(divide='ignore', invalid='ignore'):
        wt = np.log(pb / p) / np.log(pb / pt)
        wb = np.log(p / pt) / np.log(pb / pt)
        temp = 0.5 * ((wt * _take(t, kt) + wb * tb) + tb)
        h_mid = _take(htsig, kb) + ROVG * temp * np.log(pb / p)
    h_top = htsig[0] + ROVG * t[0] * np.log(psig[0] / p)
    h_low = zsfc + ROVG * t[-1] * np.log(psfc / p)
    tm = 0.5 * (t[-1] + t[-2])
    h_blw = zsfc + (tm / LRATE) * (1.0 - np.exp(ROVG * LRATE * np.log(p / psfc)))
    hp = np.where(p <= psig[0], h_top,
           np.where(p < psig[-1], h_mid,
             np.where(p <= psfc, h_low, h_blw)))
    return hp.astype(np.float32)


def _hydrostatic_height(t, psig, psfc, topo):
    htsig = np.empty(psig.shape, dtype=np.float32)
    htsig[-1] = topo + ROVG * t[-1] * np.log(psfc / psig[-1])
    dh = ROVG * 0.5 * (t[:-1] + t[1:]) * np.log(psig[1:] / psig[:-1])
    htsig[:-1] = htsig[-1] + np.cumsum(dh[::-1], axis=0)[::-1]
    return htsig


def height_hydro(t, pstar, topo, sig, ptop, p):
    """Height (m) of the pressure level p (Pa) for the hydrostatic core."""
    ptp = _scalar(ptop) * 100.0
    sig = np.asarray(sig, dtype=np.float32)[:, np.newaxis, np.newaxis]
    psig = sig * (pstar - ptp) + ptp
    htsig = _hydrostatic_height(t, psig, pstar, topo)
    return _height_on_pressure(t, psig, htsig, pstar, topo, p)


def height_nonhydro(t, pstar, ps, topo, sig, ptop, pp, p):
    """Height (m) of the pressure level p (Pa) for the non-hydrostatic core."""
    ptp = _scalar(ptop) * 100.0
    sig = np.asarray(sig, dtype=np.float32)[:, np.newaxis, np.newaxis]
    psig = sig * (pstar - ptp) + ptp + pp
    htsig = _hydrostatic_height(t, psig, ps, topo)
    return _height_on_pressure(t, psig, htsig, ps, topo, p)


def height_moloch(t, pai, ps, topo, a, b, p):
    """Height (m) of the pressure level p (Pa) for the MOLOCH core."""
    a = np.asarray(a, dtype=np.float32)[:, np.newaxis, np.newaxis]
    b = np.asarray(b, dtype=np.float32)[:, np.newaxis, np.newaxis]
    htsig = a + b * topo
    psig = P00 * pai ** CPOVR
    return _height_on_pressure(t, psig, htsig, ps, topo, p)


# CAPE and CIN (former mod_capecin)
#
# The parcel computation is the getcape routine by George H. Bryan
# (NCAR, version 1.02), aligned with Share/mod_capecin.F90 of RegCM.
# It is compiled with numba on first use.

_CAPE_KERNEL = None


def _cape_kernel():
    global _CAPE_KERNEL
    if _CAPE_KERNEL is not None:
        return _CAPE_KERNEL

    import numba

    pinc = 100.0      # Pressure increment (Pa)
    adiabat = 3       # 1 = pseudoadiabatic, liquid only
                      # 2 = reversible, liquid only
                      # 3 = pseudoadiabatic, with ice
                      # 4 = reversible, with ice
    converge = 0.002

    wlhv = 2.50080e6
    wlhs = 2.83435e6
    cpv = 1846.0932676
    cpw = 4186.95
    cpi = 2117.27
    rgas = 287.0569248
    egrav = 9.80665
    tzero = 273.15
    cpd = 1004.6992368
    rwat = 461.5233169
    p00 = 1.0e5
    ep2 = 0.6219770795
    lv1 = wlhv + (cpw - cpv) * tzero
    lv2 = cpw - cpv
    ls1 = wlhs + (cpi - cpv) * tzero
    ls2 = cpi - cpv
    reps = rwat / rgas
    rddcp = rgas / cpd
    cpdg = cpd / egrav
    ice = adiabat == 3 or adiabat == 4
    pseudo = adiabat == 1 or adiabat == 3

    @numba.njit(inline='always')
    def getqvs(p, t):
        es = 611.2 * np.exp(17.67 * (t - 273.15) / (t - 29.65))
        return ep2 * es / (p - es)

    @numba.njit(inline='always')
    def getqvi(p, t):
        es = 611.2 * np.exp(21.8745584 * (t - 273.15) / (t - 7.66))
        return ep2 * es / (p - es)

    @numba.njit(inline='always')
    def getdewp(tc, rh):
        gm = np.log(rh * np.exp((18.678 - tc / 234.5) * (tc / (257.14 + tc))))
        return tzero + 257.14 * gm / (18.678 - gm)

    @numba.njit(inline='always')
    def getthe(p, t, td, q):
        if (td - t) >= -0.1:
            tlcl = t
        else:
            # Bolton (1980) Eq. 22
            e_hpa = p * q / ((ep2 + q) * 100.0)
            tlcl = 2840.0 / (3.5 * np.log(t) - np.log(e_hpa) - 4.805) + 55.0
        return (t * ((100000.0 / p) ** (0.2854 * (1.0 - 0.28 * q))) *
                np.exp(((3376.0 / tlcl) - 2.54) * q * (1.0 + 0.81 * q)))

    @numba.njit
    def getcape(p, t, rh):
        # p, t, rh are ordered from the bottom to the top
        nk = p.size
        pi = np.empty(nk)
        td = np.empty(nk)
        q = np.empty(nk)
        th = np.empty(nk)
        thv = np.empty(nk)
        for k in range(nk):
            pi[k] = (p[k] / p00) ** rddcp
            td[k] = getdewp(t[k] - tzero, rh[k])
            q[k] = getqvs(p[k], td[k])
            th[k] = t[k] / pi[k]
            thv[k] = th[k] * (1.0 + reps * q[k]) / (1.0 + q[k])

        # Most unstable parcel (max theta-e) below 500 hPa
        kmax = 0
        if p[0] >= 50000.0:
            maxthe = 0.0
            for k in range(nk):
                if p[k] >= 50000.0:
                    the = getthe(p[k], t[k], td[k], q[k])
                    if the > maxthe:
                        maxthe = the
                        kmax = k
        k = kmax
        th2 = th[kmax]
        p2 = p[kmax]
        t2 = t[kmax]
        qv2 = q[kmax]
        b2 = 0.0

        narea = 0.0
        ql2 = 0.0
        qi2 = 0.0
        qt = qv2
        cape = 0.0
        cin = 0.0
        doit = True

        while doit and k < nk - 1:
            k += 1
            b1 = b2
            dp = p[k-1] - p[k]
            if dp < pinc:
                nloop = 1
            else:
                nloop = 1 + int(dp / pinc)
                dp = dp / nloop
            for n in range(nloop):
                p1 = p2
                t1 = t2
                th1 = th2
                qv1 = qv2
                ql1 = ql2
                qi1 = qi2
                p2 = p2 - dp
                pi2 = (p2 / p00) ** rddcp
                thlast = th1
                it = 0
                not_converged = True
                while not_converged and it < 100:
                    it += 1
                    t2 = thlast * pi2
                    if ice:
                        fliq = max(min((t2 - 233.15) / (tzero - 233.15), 1.0), 0.0)
                        fice = 1.0 - fliq
                    else:
                        fliq = 1.0
                        fice = 0.0
                    qv2 = min(qt, fliq * getqvs(p2, t2) + fice * getqvi(p2, t2))
                    qi2 = max(fice * (qt - qv2), 0.0)
                    ql2 = max(qt - qv2 - qi2, 0.0)
                    tbar = 0.5 * (t1 + t2)
                    qvbar = 0.5 * (qv1 + qv2)
                    qlbar = 0.5 * (ql1 + ql2)
                    qibar = 0.5 * (qi1 + qi2)
                    lhv = lv1 - lv2 * tbar
                    lhs = ls1 - ls2 * tbar
                    rm = rgas + rwat * qvbar
                    cpm = cpd + cpv * qvbar + cpw * qlbar + cpi * qibar
                    th2 = th1 * np.exp(lhv * (ql2 - ql1) / (cpm * tbar) +
                                       lhs * (qi2 - qi1) / (cpm * tbar) +
                                       (rm / cpm - rgas / cpd) * np.log(p2 / p1))
                    if abs(th2 - thlast) < converge:
                        not_converged = False
                    else:
                        thlast = thlast + 0.3 * (th2 - thlast)
                if pseudo:
                    qt = qv2
                    ql2 = 0.0
                    qi2 = 0.0

            thv2 = th2 * (1.0 + reps * qv2) / (1.0 + qv2 + ql2 + qi2)
            b2 = egrav * (thv2 - thv[k]) / thv[k]
            dz = -cpdg * 0.5 * (thv[k] + thv[k-1]) * (pi[k] - pi[k-1])

            if b2 >= 0.0 and b1 <= 0.0:
                # first trip into positive area
                frac = b2 / (b2 - b1) if b2 != b1 else 0.0
                parea = 0.5 * b2 * dz * frac
                narea = narea - 0.5 * b1 * dz * (1.0 - frac)
                cin += narea
                narea = 0.0
            elif b2 < 0.0 and b1 > 0.0:
                # first trip into negative area
                frac = b1 / (b1 - b2)
                parea = 0.5 * b1 * dz * frac
                narea = -0.5 * b2 * dz * (1.0 - frac)
            elif b2 < 0.0:
                # still collecting negative buoyancy
                parea = 0.0
                narea = narea - 0.5 * dz * (b1 + b2)
            else:
                # still collecting positive buoyancy
                parea = 0.5 * dz * (b1 + b2)
                narea = 0.0

            cape += max(0.0, parea)

            # stop if b < 0 and p < 100 hPa
            if p[k] <= 10000.0 and b2 < 0.0:
                doit = False
        return cape, cin

    @numba.njit(parallel=True, cache=True)
    def kernel(p, t, rh):
        kz, ny, nx = t.shape
        cape = np.empty((ny, nx), dtype=np.float32)
        cin = np.empty((ny, nx), dtype=np.float32)
        for ij in numba.prange(ny * nx):
            j = ij // nx
            i = ij % nx
            pa = np.empty(kz)
            ta = np.empty(kz)
            ra = np.empty(kz)
            for k in range(kz):
                kk = kz - 1 - k
                pa[kk] = p[k, j, i]
                ta[kk] = t[k, j, i]
                ra[kk] = rh[k, j, i] * 0.01
            cape[j, i], cin[j, i] = getcape(pa, ta, ra)
        return cape, cin

    _CAPE_KERNEL = kernel
    return kernel


def _getcape(p, t, rh):
    kernel = _cape_kernel()
    return kernel(np.ascontiguousarray(p, dtype=np.float64),
                  np.ascontiguousarray(t, dtype=np.float64),
                  np.ascontiguousarray(rh, dtype=np.float64))


def getcape_hy(ps, t, rh, sigma, ptop):
    """CAPE and CIN (J/kg) for the hydrostatic core. rh is in %."""
    ptp = _scalar(ptop) * 100.0
    sigma = np.asarray(sigma, dtype=np.float32)[:, np.newaxis, np.newaxis]
    return _getcape(sigma * (ps - ptp) + ptp, t, rh)


def getcape_nhy(ps, t, p0, rh, sigma, ptop, pp):
    """CAPE and CIN (J/kg) for the non-hydrostatic core. rh is in %."""
    ptp = _scalar(ptop) * 100.0
    sigma = np.asarray(sigma, dtype=np.float32)[:, np.newaxis, np.newaxis]
    return _getcape(sigma * (p0 - ptp) + ptp + pp, t, rh)


def getcape_moloch(ps, t, rh, pai):
    """CAPE and CIN (J/kg) for the MOLOCH core. rh is in %."""
    return _getcape(P00 * pai ** CPOVR, t, rh)
