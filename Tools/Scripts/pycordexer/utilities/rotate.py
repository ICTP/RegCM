import numpy as np


__copyright__ = 'Copyright (C) 2018 ICTP'
__author__ = 'Graziano Giuliani <ggiulian@ictp.it>'
__credits__ = ["Stefano Piani", "Graziano Giuliani"]


def get_x(lon, clon, cone):
    if clon >= 0.0 and lon >= 0.0 or clon < 0.0 and lon < 0.0:
        return np.radians(clon - lon) * cone

    elif clon >= 0.0:
        if abs(clon - lon + 360.0) < abs(clon - lon):
            return np.radians(clon - lon + 360) * cone
        else:
            return np.radians(clon - lon) * cone

    elif abs(clon - lon - 360.0) < abs(clon - lon):
        return np.radians(clon - lon - 360) * cone

    else:
        return np.radians(clon - lon) * cone


def grid_to_earth_uvrotate(proj, lon, lat, clon, clat, cone=None,
                           plon=None, plat=None):
    if proj == 'NORMER':
        return 1, 0

    elif proj == 'ROTLLR':
        if np.abs(plat - 90.0) < 0.001:
            return 1, 0
        plam = np.radians(plon)
        pphi = np.radians(plat)
        phi = np.radians(lat)
        lam = np.where(abs(lat)>89.99999, 0.0, np.radians(lon))
        dlam = plam - lam
        f1 = np.cos(pphi)*np.sin(dlam)
        f2 = np.cos(phi)*np.sin(pphi) - np.sin(phi)*np.cos(pphi)*np.cos(dlam)
        delta = np.arctan(f1/f2)
        return np.cos(delta), np.sin(delta)

    elif proj == 'ROTMER':
        zphi = np.radians(lat)
        zrla = np.radians(lon)
        zrla = np.where(abs(lat) > 89.99999, 0.0, zrla)
        if plat > 0.0:
            pollam = plon + 180.0
            polphi = 90.0 - plat
        else:
            pollam = plon
            polphi = 90.0 + plat
        if pollam > 180.0:
            pollam = pollam - 360.0

        polcphi = np.cos(np.radians(polphi))
        polsphi = np.sin(np.radians(polphi))
        zrlap = np.radians(pollam) - zrla
        zarg1 = polcphi * np.sin(zrlap)
        zarg2 = polsphi*np.cos(zphi) - polcphi*np.sin(zphi)*np.cos(zrlap)
        znorm = 1.0/np.sqrt(zarg1**2+zarg2**2)
        sindel = zarg1*znorm
        cosdel = zarg2*znorm

        return cosdel, sindel

    else:
        if np.isscalar(lon):
            x = get_x(lon, clon, cone)
        else:
            c = np.vectorize(get_x, excluded=['clon', 'cone'])
            x = c(lon, clon, cone)
        xc = np.cos(x)
        xs = np.sin(x)

        if clat >= 0:
            xs *= -1

        return xc, xs
