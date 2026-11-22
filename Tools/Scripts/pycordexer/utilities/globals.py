import glob
import json
import os
import re

import numpy as np

# The current version of RegCM. This information will be used only if there is
# no way to get it from the RegCM files
ICTP_Institution_fallback = 'ICTP'
ICTP_Model_fallback = 'RegCM5-0'
ICTP_Model_Version_fallback = 'v1-r1'

# Get the path of the directory that contains the software
MAIN_DIR = os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + '/..')

# The directory where the CORDEX file will be saved
OUTPUTDIR = 'CORDEX-CMIP6'

# NetCDF output options
DATATYPE_AUXILIARIES = 'f8'
DATATYPE_MAIN = 'f'
COMPRESSION = 'zlib'
SHUFFLE = True
FLETCHER32 = True
COMPRESSION_LEVEL = 1
# Remove from the CORDEX files the lateral boundary zone (nspgx-2 points
# per side, the boundary_points attribute of the RegCM files)
CROP_BOUNDARY = True
# RegCM horizontal dimensions cropped with the boundary zone
HORIZONTAL_DIMS = ('jx', 'iy')

# Default number of significant decimal digits kept by the netCDF quantization
# (BitGroom) of single precision variables
SIGNIFICANT_DIGITS = 4

# Significant digits per variable. The values follow the GRIB2 packing
# precision of the NCEP Unified Post Processor (UPP, parm/post_avblflds.xml,
# the <scale> of each field, also used for the WRF GRIB2 output). GRIB2 simple
# packing keeps "scale" significant digits of the field range, an absolute
# precision; the netCDF quantization keeps significant digits of every value,
# a relative precision. The numbers below give, at the typical magnitude of
# each variable, a precision at least as fine as the UPP GRIB2 one:
#   temperatures 4 (UPP TMP 2 m/surface 4, isobaric 3: 0.1 K, here ~0.02 K)
#   pressures 6 (PRES, PRMSL 6: 0.1 Pa)
#   geopotential height and orography 6 (HGT 6: 0.01 m)
#   boundary layer height 5 (HPBL 6 over a 3 km range: 0.01 m)
#   specific humidity 4 (SPFH 4), relative humidity 3 (RH 2 m 3)
#   horizontal wind 4 (UGRD/VGRD 4), gusts 3 (GUST 3)
#   vertical wind 3 (DZDT 3)
#   precipitation and evaporation fluxes 4 (CPRAT, EVP 4)
#   runoff 4 (SSRUN, BGRUN 4), soil moisture 3 (SOILW, SOILL 3)
#   snow 3 (WEASD, SNOD, SNOWC 3), snow melt 4 (SNOM 4)
#   turbulent and radiative fluxes 4 (SHTFL, LHTFL, DSWRF... 4)
#   momentum flux 3 (UFLX 3)
#   cloud cover 3 (TCDC 3), cloud water and ice path 3 (TCOLW 3)
#   precipitable water 3 (PWAT 3)
#   CAPE and CIN 4 (CAPE, CIN 4), lifted index 3 (LFTX 3)
#   sea ice 3 (ICEC, ICETK 3), roughness length 3 (SFCR 2.7)
#   total ozone 4 (TOZNE 4)
#   land and sea fractions 3, cell areas 6 (no UPP field: area weighted
#   sums need the extra digits)
# The variables on pressure or height levels (ta850, ua100m, ...) take the
# value of their family (ta, ua, ...). Not listed variables (aerosols,
# chemistry, ocean) use SIGNIFICANT_DIGITS.
SIGNIFICANT_DIGITS_VAR = {
    'tas': 4, 'tasmax': 4, 'tasmin': 4, 'ts': 4, 'tsl': 4, 'ta': 4,
    'tos': 4, 'tob': 4, 'sitemptop': 4,
    'ps': 6, 'psl': 6,
    'zg': 6, 'orog': 6,
    'zmla': 5,
    'huss': 4, 'hus': 4, 'hurs': 3,
    'uas': 4, 'vas': 4, 'ua': 4, 'va': 4, 'sfcWind': 4, 'sfcWindmax': 4,
    'wsgsmax': 3,
    'wa': 3,
    'pr': 4, 'prc': 4, 'prsn': 4, 'prhmax': 4,
    'evspsbl': 4, 'evspsblpot': 4,
    'mrro': 4, 'mrros': 4,
    'mrso': 3, 'mrsos': 3, 'mrsol': 3, 'mrsfl': 3, 'mrfso': 3, 'mrfsos': 3,
    'mrsofc': 3,
    'snw': 3, 'snd': 3, 'snc': 3, 'snm': 4,
    'hfls': 4, 'hfss': 4, 'hfds': 4,
    'rsds': 4, 'rsus': 4, 'rlds': 4, 'rlus': 4, 'rsdt': 4, 'rsut': 4,
    'rlut': 4, 'rsdscs': 4, 'rsuscs': 4, 'rldscs': 4, 'rluscs': 4,
    'rsutcs': 4, 'rlutcs': 4, 'rsdsdir': 4, 'rsntds': 4,
    'rsdsaf': 4, 'rsusaf': 4, 'rldsaf': 4, 'rlusaf': 4, 'rsutaf': 4,
    'rlutaf': 4, 'rsdscsaf': 4, 'rsuscsaf': 4, 'rldscsaf': 4,
    'rluscsaf': 4, 'rsutcsaf': 4, 'rlutcsaf': 4,
    'tauu': 3, 'tauv': 3,
    'clt': 3, 'cll': 3, 'clm': 3, 'clh': 3,
    'clwvi': 3, 'clivi': 3, 'prw': 3,
    'cape': 4, 'capemax': 4, 'cin': 4, 'cinmax': 4, 'li': 3, 'limax': 3,
    'siconc': 3, 'siconca': 3, 'sithick': 3,
    'z0': 3,
    'toz': 4,
    'sftlf': 3, 'sftgif': 3, 'sftlaf': 3, 'sftof': 3, 'sfturf': 3,
    'areacella': 6, 'areacellr': 6, 'areacello': 6,
}

# Families of the variables on pressure or height levels
LEVEL_VARIABLE_FAMILY = re.compile(r'^(ta|ua|va|wa|zg|hus)\d+m?$')


def significant_digits(var_name):
    """
    Number of significant decimal digits kept by the quantization of a
    single precision variable (see SIGNIFICANT_DIGITS_VAR)
    """
    if var_name in SIGNIFICANT_DIGITS_VAR:
        return SIGNIFICANT_DIGITS_VAR[var_name]
    family = LEVEL_VARIABLE_FAMILY.match(var_name)
    if family:
        return SIGNIFICANT_DIGITS_VAR[family.group(1)]
    return SIGNIFICANT_DIGITS

# Some dimensions in the RegCM files must have a different name in the CORDEX
# files. The following dictionary stores the needed replacements
REPLACE_DIMS = {
    'jx': 'x',
    'iy': 'y',
    'ntimes': 'time_bnds',
}
REPLACE_DIMS_ROT = {
    'jx': 'rlon',
    'iy': 'rlat',
    'ntimes': 'time_bnds',
}

# The same for the variables
REPLACE_VARS = {
    'jx': 'x',
    'iy': 'y',
    'xlat': 'lat',
    'xlon': 'lon',
    'xlat_bnds': 'lat_bnds',
    'xlon_bnds': 'lon_bnds',
}
REPLACE_VARS_ROT = {
    'jx': 'rlon',
    'iy': 'rlat',
    'xlat': 'lat',
    'xlon': 'lon',
    'xlat_bnds': 'lat_bnds',
    'xlon_bnds': 'lon_bnds',
}

# These attributes will not be copied from the RegCM files into the CORDEX ones
EXCLUDED_ATTRIBUTES = [
    '_CoordinateAxisType',
    '_CoordinateAxisTypes',
    '_CoordinateTransformType'
]

# Read the JSON files with the definitions of the cordex variables
CORDEX_VARS = {}
CORDEX_VARS_DIR = os.path.join(MAIN_DIR, 'variables/cordex_vars')
for json_file in os.listdir(CORDEX_VARS_DIR):
    if not json_file.lower().endswith('json'):
        continue

    with open(os.path.join(CORDEX_VARS_DIR, json_file), 'r') as f:
        try:
            for var_name, var_description in json.load(f).items():
                CORDEX_VARS[var_name.lower()] = var_description
        except Exception:
            print('Error reading '+json_file)
            raise(SyntaxError)


# Read the JSON file with the association among files and vars
ASSOCIATION_FILE = os.path.join(MAIN_DIR, 'variables/file_associations.json')
with open(ASSOCIATION_FILE, 'r') as f:
    try:
        ASSOCIATIONS = json.load(f)
    except Exception:
        print('Error reading '+ASSOCIATION_FILE)
        raise(SyntaxError)


# The values in the mask field are regular expressions. They will be replaced
# with their compiled version
for data_file in ASSOCIATIONS:
    data_file['mask'] = re.compile(data_file['mask'])


# Read the JSON file with the association among RegCM variables and their names
# inside the RegCM NetCDF files
REGCM_VAR_FILE = os.path.join(MAIN_DIR, 'variables/regcm_vars.json')
with open(REGCM_VAR_FILE, 'r') as f:
    try:
        REGCM_VARS = json.load(f)
    except Exception:
        print('Error reading '+REGCM_VAR_FILE)
        raise(SyntaxError)

CORDEX_CMIP6_DEFINITIONS = { }
jsonfiles = glob.glob(os.path.join(MAIN_DIR,'variables','CORDEX-CMIP6*json'))
for jf in jsonfiles:
    with open(jf, 'r') as f:
        try:
            elements = json.load(f)
        except Exception:
            print('Error reading ' + jf)
            raise(SyntaxError)
        for k in elements:
            CORDEX_CMIP6_DEFINITIONS[k] = elements[k]

# The variables of the CORDEX-CMIP6 tables of each frequency. They give the
# cell_methods of the variables written at that frequency
CORDEX_TABLE_VARS = {}
for frequency in ('1hr', '6hr', 'day', 'mon', 'fx'):
    table_file = f'CORDEX-CMIP6_{frequency}.json'
    with open(os.path.join(MAIN_DIR, 'variables', table_file)) as f:
        CORDEX_TABLE_VARS[frequency] = json.load(f)['variable_entry']
del frequency, table_file

# The variables of the monthly table, used to select and describe the monthly
# means computed from the daily variables
CORDEX_MON_VARS = CORDEX_TABLE_VARS['mon']

# Default fill values for NETCDF. These are the same values that you can find
# from the netCDF4._default_fillvals dictionary
NETCDF_DEFAULT_FILL_VALUES = {
    'i8': -9223372036854775806,
    'f4': 9.969209968386869e+36,
    'f': 9.969209968386869e+36,
    'u8': 18446744073709551614,
    'i1': -127,
    'u4': 4294967295,
    'S1': '\x00',
    'i2': -32767,
    'u1': 255,
    'i4': -2147483647,
    'u2': 65535,
    'f8': 9.969209968386869e+36,

    np.int64: -9223372036854775806,
    np.float32: 9.969209968386869e+36,
    np.uint64: 18446744073709551614,
    np.int8: -127,
    np.uint32: 4294967295,
    np.char: '\x00',
    np.int16: -32767,
    np.uint8: 255,
    np.int32: -2147483647,
    np.uint16: 65535,
    np.float64: 9.969209968386869e+36,

    np.int64(0).dtype: -9223372036854775806,
    np.float32(0).dtype: 9.969209968386869e+36,
    np.uint64(0).dtype: 18446744073709551614,
    np.int8(0).dtype: -127,
    np.uint32(0).dtype: 4294967295,
    np.int16(0).dtype: -32767,
    np.uint8(0).dtype: 255,
    np.int32(0).dtype: -2147483647,
    np.uint16(0).dtype: 65535,
    np.float64(0).dtype: 9.969209968386869e+36
}
