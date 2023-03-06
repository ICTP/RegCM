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
                           plon=None, plat=None, ds=None):
    if proj == 'NORMER':
        return 1, 0

    elif proj == 'ROTLLR':
        drot = ds/6371229.0
        plam = 0.0
        pphi = 0.0
        if np.radians(plon) < 0.0:
            plam = np.radians(plon) + np.pi
        elif np.radians(plon) > 0.0:
            plam = np.radians(plon) - np.pi
        if np.radians(plat) > 0.0:
            pphi = np.pi*0.5 - np.radians(plat)
        elif np.radians(plat) < 0.0:
            pphi = - (np.pi*0.5 + np.radians(plat))
        if clat > 90.0:
            cphi = np.radians(90.0-clat)
        elif clat < -90.0:
            cphi = np.radians(clat+90.0)
        else:
            cphi = np.radians(clat)
        if clon > 180.0:
            clam = np.radians(clon-360.0)
        elif clon < -180.0:
            clam = np.radians(clon+360.0)
        else:
            clam = np.radians(clon)
        rphi0 = np.arcsin(-np.cos(cphi)*np.sin(pphi)*np.cos(clam-plam) +
                           np.sin(cphi)*np.cos(pphi))
        if np.abs(np.abs(rphi0)-np.pi*0.5) > 1.0e-7 and np.abs(pphi) > 1.0e-7:
            rlam0 = ((np.sin(cphi)-np.cos(pphi)*np.sin(rphi0)) /
                     (np.sin(pphi)*np.cos(rphi0)))
            if rlam0 < -1.0 and rlam0 > -1.00001:
                rlam0 = -1.0
            if rlam0 > 1.0 and rlam0 < 1.00001:
                rlam0 = 1.0
            rlam0 = np.arccos(rlam0)
            if clam < plam:
                rlam0 = -rlam0
        else:
            rlam0 = clam
        phi = np.radians(lat)
        lam = np.where(abs(lat)>89.99999, 0.0, np.radians(lon))
        ny, nx = np.shape(lat)
        yy = np.arange(-ny/2+0.5,ny/2)
        xx = np.tile(np.array([ x 
                        for x in np.arange(-nx/2+0.5,nx/2)] ),(ny,1))
        yy = np.tile(np.array([[y 
                        for y in np.arange(-ny/2+0.5,ny/2)]]).transpose(),nx)
        rotlam = xx * drot + rlam0
        rotphi = yy * drot + rphi0
        dlam = lam - plam
        f5 = np.where(np.abs(np.cos(phi))>1.0e-7,
              -(np.sin(pphi)*np.sin(rotlam))/np.cos(phi), np.nan)
        f6 = np.where(np.abs(np.cos(phi))>1.0e-7,
              (np.cos(pphi)*np.cos(rotphi) -
              np.sin(pphi)*np.sin(rotphi)*np.cos(rotlam))/np.cos(phi),np.nan)
        f7 = np.where(np.abs(np.sin(dlam))>0.075e-1,
              np.cos(rotphi)/(np.sin(dlam)*np.sin(pphi)),np.nan)
        f8 = np.where(np.abs(np.cos(rotphi))>1.0e-7,
              -(np.cos(dlam)*np.sin(pphi)*np.sin(phi) +
               np.cos(pphi)*np.cos(phi))/np.cos(rotphi),np.nan)
        return f5, f6, f7, f8

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
