import os
import json
import re
import numpy as np

# The current version of RegCM. This information will be used only if there is
# no way to get it from the RegCM files
ICTP_Model_fallback = 'ICTP-RegCM5-0'
ICTP_Model_Version_fallback = 'v0'

# Get the path of the directory that contains the software
MAIN_DIR = os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + '/..')

# The directory where the CORDEX file will be saved
OUTPUTDIR = 'CORDEX'

# NetCDF output options
DATATYPE_AUXILIARIES = 'f8'
DATATYPE_MAIN = 'f'
COMPRESSION = True
SHUFFLE = True
FLETCHER32 = True
COMPRESSION_LEVEL = 1

# Some dimensions in the RegCM files must have a different name in the CORDEX
# files. The following dictionary stores the needed replacements
REPLACE_DIMS = {
    'jx': 'x',
    'iy': 'y',
    'ntimes': 'time_bnds',
}

# The same for the variables
REPLACE_VARS = {
    'jx': 'x',
    'iy': 'y',
    'xlat': 'lat',
    'xlon': 'lon',
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
        except:
            print('Error reading '+json_file)
            raise(SyntaxError)


# Read the JSON file with the association among files and vars
ASSOCIATION_FILE = os.path.join(MAIN_DIR, 'variables/file_associations.json')
with open(ASSOCIATION_FILE, 'r') as f:
    try:
        ASSOCIATIONS = json.load(f)
    except:
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
    except:
        print('Error reading '+REGCM_VAR_FILE)
        raise(SyntaxError)

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
