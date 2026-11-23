#!/usr/bin/env python3
#
#    This file is part of ICTP RegCM.
#
#    Use of this source code is governed by an MIT-style license that can
#    be found in the LICENSE file or at
#
#         https://opensource.org/licenses/MIT.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# Download the CAMS global reanalysis (EAC4) on pressure levels for the
# chemistry and aerosol boundary conditions of chem_icbc (chemtyp = 'CAMSR'),
# in the layout read by PreProc/ICBC/mod_ch_cams.F90:
#
#     <inpglob>/CAMS/<year>/<NAME>_<year>_<month>.nc
#
# one file per variable and month, 6 hourly (00, 06, 12, 18 UTC). All the
# variables of a month are retrieved with one request (one wait in the
# queue of the data store) and split per variable (needs netCDF4).
# Needs the cdsapi package and an Atmosphere Data Store key: the
# ~/.cdsapirc key of the Copernicus services works, the url is set here.
# The CAMS licence must be accepted once on the ADS web page of the dataset.
#
# Examples:
#   cams_eac4_download.py --dest $INPGLOB/CAMS --start 2023-09 --end 2023-10
#   cams_eac4_download.py --dest CAMS --start 2023-09 --end 2023-09 \
#       --area 25 -120 -65 0 --only AE
#
# The area (north west south east, degrees) should cover the domain with a
# margin of a few degrees. Without --area the whole globe is downloaded.

import argparse
import calendar
import os
import sys
import zipfile

ADS_URL = 'https://ads.atmosphere.copernicus.eu/api'
DATASET = 'cams-global-reanalysis-eac4'

# file name -> ADS variable (aerosols, AE)
AEROSOLS = {
    'SSLT1': 'sea_salt_aerosol_0.03-0.5um_mixing_ratio',
    'SSLT2': 'sea_salt_aerosol_0.5-5um_mixing_ratio',
    'SSLT3': 'sea_salt_aerosol_5-20um_mixing_ratio',
    'DUST1': 'dust_aerosol_0.03-0.55um_mixing_ratio',
    'DUST2': 'dust_aerosol_0.55-0.9um_mixing_ratio',
    'DUST3': 'dust_aerosol_0.9-20um_mixing_ratio',
    'OCHL': 'hydrophilic_organic_matter_aerosol_mixing_ratio',
    'OCHB': 'hydrophobic_organic_matter_aerosol_mixing_ratio',
    'BCHL': 'hydrophilic_black_carbon_aerosol_mixing_ratio',
    'BCHB': 'hydrophobic_black_carbon_aerosol_mixing_ratio',
    'SO4': 'sulphate_aerosol_mixing_ratio',
    'SO2': 'sulphur_dioxide',
}

# file name -> ADS variable (gas phase chemistry, CH)
GASES = {
    'O3': 'ozone',
    'CO': 'carbon_monoxide',
    'ETH': 'ethane',
    'CH2O': 'formaldehyde',
    'C3H8': 'propane',
    'NO': 'nitrogen_monoxide',
    'NO2': 'nitrogen_dioxide',
    'H2O2': 'hydrogen_peroxide',
    'ISOP': 'isoprene',
    'SO2': 'sulphur_dioxide',
    'HNO3': 'nitric_acid',
    'CH4': 'methane_chemistry',
    'PAN': 'peroxyacetyl_nitrate',
}

# file name -> variable name in the netCDF files of the data store
SHORT = {'SSLT1': 'aermr01', 'SSLT2': 'aermr02', 'SSLT3': 'aermr03',
         'DUST1': 'aermr04', 'DUST2': 'aermr05', 'DUST3': 'aermr06',
         'OCHL': 'aermr07', 'OCHB': 'aermr08', 'BCHL': 'aermr09',
         'BCHB': 'aermr10', 'SO4': 'aermr11', 'SO2': 'so2', 'O3': 'go3',
         'CO': 'co', 'ETH': 'c2h6', 'CH2O': 'hcho', 'C3H8': 'c3h8',
         'NO': 'no', 'NO2': 'no2', 'H2O2': 'h2o2', 'ISOP': 'c5h8',
         'HNO3': 'hno3', 'CH4': 'ch4', 'PAN': 'pan'}

LEVELS = ['1', '2', '3', '5', '7', '10', '20', '30', '50', '70', '100',
          '150', '200', '250', '300', '400', '500', '600', '700', '800',
          '850', '900', '925', '950', '1000']


def months(start, end):
    y, m = (int(x) for x in start.split('-'))
    ye, me = (int(x) for x in end.split('-'))
    while (y, m) <= (ye, me):
        yield y, m
        m += 1
        if m > 12:
            y, m = y + 1, 1


def unpack(path, odir):
    """the retrieved file: one netCDF, or a zip of netCDF files"""
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as z:
            names = [os.path.join(odir, n) for n in z.namelist()]
            z.extractall(odir)
        os.remove(path)
        return names
    return [path]


def split(files, targets):
    """write each wanted variable to its own file (name -> target)"""
    import netCDF4 as nc
    for path in files:
        with nc.Dataset(path) as src:
            for name, target in targets.items():
                v = SHORT[name]
                if v not in src.variables and name == 'CH4':
                    v = 'ch4_c' if 'ch4_c' in src.variables else v
                if v not in src.variables:
                    continue
                with nc.Dataset(target + '.part', 'w') as dst:
                    dst.setncatts(src.__dict__)
                    for d, dim in src.dimensions.items():
                        dst.createDimension(
                            d, None if dim.isunlimited() else len(dim))
                    for c in src.variables.values():
                        if c.name != v and c.name not in src.dimensions:
                            continue
                        fill = c.__dict__.get('_FillValue')
                        o = dst.createVariable(c.name, c.dtype, c.dimensions,
                                               fill_value=fill, zlib=True)
                        o.setncatts({k: x for k, x in c.__dict__.items()
                                     if k != '_FillValue'})
                        o.set_auto_maskandscale(False)
                        c.set_auto_maskandscale(False)
                        o[...] = c[...]
                os.replace(target + '.part', target)
    missing = [t for t in targets.values() if not os.path.isfile(t)]
    if missing:
        sys.exit('variables not found in the retrieved data: ' +
                 ' '.join(missing))
    for path in files:
        os.remove(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.
                                 RawDescriptionHelpFormatter)
    ap.add_argument('--dest', required=True, help='CAMS directory')
    ap.add_argument('--start', required=True, help='first month YYYY-MM')
    ap.add_argument('--end', required=True, help='last month YYYY-MM')
    ap.add_argument('--area', nargs=4, type=float,
                    metavar=('N', 'W', 'S', 'E'))
    ap.add_argument('--only', choices=['AE', 'CH'],
                    help='only aerosols or gases')
    ap.add_argument('--days', help='day range of each month, e.g. 1-2 '
                    '(tests)')
    a = ap.parse_args()
    try:
        import cdsapi
    except ImportError:
        sys.exit('The cdsapi package is needed (pip install cdsapi)')
    client = cdsapi.Client(url=ADS_URL)
    want = {}
    if a.only != 'CH':
        want.update(AEROSOLS)
    if a.only != 'AE':
        want.update(GASES)
    for y, m in months(a.start, a.end):
        d0, d1 = 1, calendar.monthrange(y, m)[1]
        if a.days:
            d0, d1 = (int(x) for x in a.days.split('-'))
        odir = os.path.join(a.dest, '%04d' % y)
        os.makedirs(odir, exist_ok=True)
        targets = {}
        for name in want:
            target = os.path.join(odir, '%s_%04d_%02d.nc' % (name, y, m))
            if not os.path.isfile(target):
                targets[name] = target
        if not targets:
            continue
        req = {
            'variable': sorted(set(want[n] for n in targets)),
            'pressure_level': LEVELS,
            'date': ['%04d-%02d-%02d/%04d-%02d-%02d' %
                     (y, m, d0, y, m, d1)],
            'time': ['00:00', '06:00', '12:00', '18:00'],
            'data_format': 'netcdf',
        }
        if a.area:
            req['area'] = a.area
        tmp = os.path.join(odir, 'cams_%04d_%02d.download' % (y, m))
        print('Retrieving %d variables for %04d-%02d' % (len(targets), y, m),
              flush=True)
        client.retrieve(DATASET, req, tmp)
        split(unpack(tmp, odir), targets)


if __name__ == '__main__':
    main()
