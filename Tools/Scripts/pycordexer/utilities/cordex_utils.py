import ctypes
import time
import uuid
from logging import getLogger
from multiprocessing import Array
from os import path
from numpy.ma import is_masked
from netCDF4 import Dataset, num2date
from copy import copy
from traceback import format_exc

from utilities.globals import *
from utilities.variables import Variable, NetCDFData, get_var_with_name


__copyright__ = 'Copyright (C) 2018 ICTP'
__author__ = 'Graziano Giuliani <ggiulian@ictp.it>'
__credits__ = ["Stefano Piani", "Graziano Giuliani"]


LOGGER = getLogger(__name__)


class RegcmSimulation(object):
    """
    A wrapper that contains the metadata that identify a particular run of RegCM
    and that are saved inside the CORDEX files
    """

    def __init__(self, mail='esp@ictp.it', domain='NONE', global_model='NONE',
                 experiment='none', ensemble='NN', notes='none'):
        self.mail = mail
        self.domain = domain
        self.product = 'output'
        self.global_model = global_model
        self.experiment = experiment
        self.ensemble = ensemble
        self.notes = notes


def dtype_2_ctype(dtype):
    """
    This function is a hack to get to find a c_type that is big enough to store
    a single entry of a numpy array with a particular dtype.

    It returns a ctype that is exactly as big as the numpy datatype passed as
    input.

    There is no guarantee of coherence with the type (integers can become
    floats, for example). Only the size is preserved.

    :param dtype: A numpy dtype object
    :return: a c_type object
    """

    # This function is a really ugly hack! It would be better to find something
    # better

    # This create an array of one element using the dtype and then it returns
    # the size of its entry
    dtype_size = np.array([0, ], dtype=dtype).itemsize

    ctypes_by_size = [
        ctypes.c_byte,
        ctypes.c_short,
        ctypes.c_float,
        ctypes.c_double,
        ctypes.c_longdouble,
    ]

    ctypes_by_size.sort(key=ctypes.sizeof)
    for c in ctypes_by_size:
        if ctypes.sizeof(c) == dtype_size:
            return c
        elif ctypes.sizeof(c) > dtype_size:
            break

    raise ValueError(
        'No ctype found with size {} (suitable for dtype {})'
        .format(dtype_size, dtype)
    )


def read_rcm_map(netcdf_file_pointer):
    LOGGER.debug('Trying to read rcm_map variable from rcm_map variable')
    try:
        rcm_map = netcdf_file_pointer.variables['rcm_map']
        LOGGER.debug('Variable rcm_map read!')
        orcmap = 'rcm_map'
        return rcm_map, orcmap
    except KeyError:
        LOGGER.debug('No rcm_map variable found')

    LOGGER.debug('Trying to read rcm_map variable from crs variable')
    try:
        rcm_map = netcdf_file_pointer.variables['crs']
        LOGGER.debug('Variable crs read!')
        orcmap = 'crs'
        return rcm_map, orcmap
    except KeyError:
        LOGGER.debug('No crs variable found')

    LOGGER.debug('No rcm_map found. rcm_map value will be set to "None"')
    rcm_map = None
    orcmap = None
    return rcm_map, orcmap


def get_first_and_last_date_str(dates, frequency):
    if frequency == 'month':
        dd1 = '{:04d}{:02d}{:02d}'.format(
            dates[0].year,
            dates[0].month,
            dates[0].day
        )
        dd2 = '{:04d}{:02d}{:02d}'.format(
            dates[-2].year,
            dates[-2].month,
            dates[-2].day
        )
    elif frequency == 'day':
        # Notice the trailing '12'.
        dd1 = '{:04d}{:02d}{:02d}'.format(
            dates[0].year,
            dates[0].month,
            dates[0].day
        )
        try:
            dd2 = '{:04d}{:02d}{:02d}'.format(
                dates[-2].year,
                dates[-2].month,
                dates[-2].day
            )
        except:
            dd2 = dd1
    else:
        dd1 = '{:04d}{:02d}{:02d}00'.format(
            dates[0].year,
            dates[0].month,
            dates[0].day
        )
        dd2 = '{:04d}{:02d}{:02d}00'.format(
            dates[-1].year,
            dates[-1].month,
            dates[-1].day
        )
    return dd1, dd2


class RegcmOutputFile(object):
    """
    This object is a wrapper around a NetCDF file written by RegCM as output
    It contains all the general attributes of the files that are useful
    regardless the specific CORDEX variable that must be saved

    :param ncf: an open Dataset object from the library netCDF4
    :param datafile (optional): the path of the NetCDF file (used just for logs)
    """

    def __init__(self, ncf, datafile=None, regcm_version=None,
                 regcm_version_id=None, regcm_nest_tag=None):

        # in_file and _of_file are variables used just inside the log lines to
        # report the name of the file where the operations have been performed
        if datafile is None:
            in_file = ''
            of_file = ''
            self._name = 'RegcmFile'
        else:
            in_file = ' in file {}'.format(path.basename(datafile))
            of_file = ' of file {}'.format(path.basename(datafile))
            self._name = 'RegcmFile({})'.format(path.basename(datafile))

        # Save a pointer to the netcdf file inside this object
        # When init ends, this pointer will be destroyed
        self.__ncf = ncf

        times = ncf.variables['time']
        if len(times) < 1:
            raise ValueError('No timesteps found{}'.format(in_file))

        LOGGER.debug('Reading attributes%s', of_file)
        self.attributes = {}
        for attribute in ncf.ncattrs():
            LOGGER.debug('Reading attribute "%s"%s', attribute, in_file)
            attribute_value = ncf.getncattr(attribute)
            LOGGER.debug(
                'Saving value "%s" as value for attribute "%s"%s',
                attribute_value,
                attribute,
                in_file
            )
            self.attributes[attribute] = attribute_value

        LOGGER.debug('Copying values for times%s in memory', of_file)
        self._save_in_memory(times, 'times', np.float32)

        LOGGER.debug('Cleaning attribute of variable times')
        if 'calendar' in self.times_attributes:
            current_val = self.times_attributes['calendar']
            if 'gregorian' in current_val.lower():
                if current_val != 'proleptic_gregorian':
                    LOGGER.debug(
                        'Changing "calendar" attribute from "%s" to'
                        '"proleptic_gregorian"',
                        current_val
                    )
                self.times_attributes['calendar'] = 'proleptic_gregorian'
        if 'bounds' in self.times_attributes:
            LOGGER.debug('Removing "bounds" attribute from times')
            del self.times_attributes['bounds']
        if 'units' in self.times_attributes:
            LOGGER.debug('Removing "UTC" string from attribute "units"')
            new_value = self.times_attributes['units'].replace(' UTC', '')
            self.times_attributes['units'] = new_value
        LOGGER.debug('Attribute cleaned')

        xlat = get_var_with_name(['xlat', 'lat'], ncf)
        xlon = get_var_with_name(['xlon', 'lon'], ncf)

        LOGGER.debug('Copying values for latitudes%s in memory', of_file)
        self._save_in_memory(xlat, 'xlat', np.float32)

        LOGGER.debug('Copying values for longitudes%s in memory', of_file)
        self._save_in_memory(xlon, 'xlon', np.float32)

        iy = get_var_with_name(['iy', 'y', 'rlat'], ncf)
        jx = get_var_with_name(['jx', 'x', 'rlon'], ncf)

        LOGGER.debug('Copying values of variable iy%s in memory', of_file)
        self._save_in_memory(iy, 'iy', np.float32)

        LOGGER.debug('Copying values of variable jx%s in memory', of_file)
        self._save_in_memory(jx, 'jx', np.float32)

        LOGGER.debug(
            'Try to recognize the file type%s from its variables',
            of_file
        )
        self.ftype = self._get_file_type()

        rcm_map, orcmap = read_rcm_map(ncf)
        if orcmap is None:
            self._contains_map = False
            self._map_type = None
            LOGGER.debug(
                'No map found. The script will ignore this variable from now on'
            )
        else:
            if orcmap == 'rcm_map':
                self._contains_map = True
                self._map_type = 'rcm_map'
            elif orcmap == 'crs':
                self._contains_map = True
                self._map_type = 'crs'
            else:
                raise ValueError(
                    'Unexpected type of map type ({}). Accepted only "rcm_map" '
                    'and "crs".'.format(orcmap)
                )
            LOGGER.debug('Saving map in memory')
            self._save_in_memory(rcm_map, 'map', rcm_map.datatype)

        self._has_timebounds = True
        try:
            LOGGER.debug('Checking for time bounds%s', in_file)
            timebnds = ncf.variables['time_bnds']
            LOGGER.debug('Time bounds found!')
        except KeyError:
            LOGGER.debug(
                'No time bounds found (no variable time_bnds found%s)',
                in_file
            )
            self._has_timebounds = False
            timebnds = None

        if self._has_timebounds:
            LOGGER.debug('Copying values for time bounds in memory')
            self._save_in_memory(timebnds, 'timebounds')

        LOGGER.debug('Reading the domain attribute%s', of_file)
        self._domain = getattr(ncf, 'experiment')
        self._product = 'output'
        LOGGER.debug('The domain is %s', self._domain)
        LOGGER.debug('The product is %s', self._product)

        LOGGER.debug('Reading the RegCM revision attribute%s', of_file)
        try:
            rev_temp = getattr(ncf, 'model_revision')
            LOGGER.debug('The model revision is %s', rev_temp)
        except:
            LOGGER.warning(
                'No attribute "model_revision" found%s. "0000" will be used as '
                'placeholder!', in_file
            )
            rev_temp = 'rev0000'

        if regcm_version is not None:
            self._revision = regcm_version
            if regcm_version_id is not None:
                self._rev_version = 'v' + str(regcm_version_id)
                self._nest_tag = None
            elif regcm_nest_tag is not None:
                self._rev_version = regcm_nest_tag
                self._nest_tag = regcm_nest_tag
        else:
            try:
                if rev_temp.lower().startswith('tag'):
                    LOGGER.debug('Found a tagged version of RegCM')
                    if rev_temp.lower().startswith('tag-'):
                        cleaned_rev_temp = rev_temp[4:]
                    else:
                        cleaned_rev_temp = rev_temp[3:]
                    LOGGER.debug(
                        'Removed the initial "tag" from the name, now the '
                        'string is %s', cleaned_rev_temp
                    )
                    LOGGER.debug(
                        'Splitting string by points and read the result as a '
                        'list of integers'
                    )
                    rev_numbers = [int(i) for i in cleaned_rev_temp.split('.')]
                    rev_numbers_str = [str(i) for i in rev_numbers]
                    self._revision = '-'.join(rev_numbers_str)
                    self._rev_version = 'v{}'.format(rev_numbers[-1])
                else:
                    LOGGER.debug('Found an untagged version of RegCM')
                    if rev_temp.lower().startswith('rev'):
                        rev_temp = rev_temp[3:]
                    self._revision = '4-git' + rev_temp
                    self._rev_version = 'v0'
                LOGGER.debug(
                    'The model revision is %s; the version is %s',
                    self._revision,
                    self._rev_version
                )
                self._nest_tag = None
            except Exception:
                LOGGER.warning(
                    'Unable to understand the revision "%s"%s for the following '
                    'reason:\n{}'.format(format_exc()),
                    rev_temp,
                    of_file
                )
                self._revision = None
                self._rev_version = None

        self.__ncf = None

    def _get_file_type(self):
        def check_variable(var_name):
            LOGGER.debug('Looking for variable "%s"', var_name)
            for v_name in get_regcm_variable_names(var_name):
                if v_name in self.__ncf.variables:
                    LOGGER.debug(
                        'Variable "%s" found with name "%s"',
                        v_name,
                        var_name
                    )
                    return True
            else:
                LOGGER.debug('Variable "%s" not found', var_name)

        if check_variable('tsmax'):
            LOGGER.debug('Variable "tsmax" found: the type is "STS"')
            return 'STS'
        if check_variable('tsmin'):
            LOGGER.debug('Variable "tsmin" found: the type is "STS"')
            return 'STS'

        if check_variable('lakets'):
            LOGGER.debug('Variable "lakets" found: the type is "LAK"')
            return 'LAK'

        if check_variable('ua'):
            LOGGER.debug('Variable "ua" found: the type is "ATM"')
            return 'ATM'
        if check_variable('u'):
            LOGGER.debug('Variable "u" found: the type is "ATM"')
            return 'ATM'

        if check_variable('qrs'):
            LOGGER.debug('Variable "qrs" found: the type is "RAD"')
            return 'RAD'
        if check_variable('firtp'):
            LOGGER.debug('Variable "firtp" found: the type is "RAD"')
            return 'RAD'
        if check_variable('rts'):
            LOGGER.debug('Variable "rts" found: the type is "RAD"')
            return 'RAD'
        if check_variable('rsut'):
            LOGGER.debug('Variable "rts" found: the type is "RAD"')
            return 'RAD'

        if check_variable('prhmax'):
            LOGGER.debug('Variable "prhmax" found: the type is "SHF"')
            return 'SHF'

        LOGGER.debug('No known variable found: assuming SRF')
        return 'SRF'

    def _save_in_memory(self, ncf_var, var_name, dtype=None):
        """
        Save the content of a netcdf variable in the shared memory (taking
        care of converting km in meters if needed)
        """
        # If no dtype is specified, use the one of the NetCDF variable
        if dtype is None:
            dtype = ncf_var.datatype

        # Save the shape
        setattr(self, '_' + var_name + '_shape', ncf_var.shape)

        # Save the data type
        setattr(self, '_' + var_name + '_dtype', dtype)

        # Save the data in a multiprocessing array
        LOGGER.debug('Finding a suitable ctype to save a %s type', dtype)
        c_datatype = dtype_2_ctype(dtype)
        LOGGER.debug('Using %s to save a %s type', c_datatype, dtype)

        LOGGER.debug('Allocating space in memory')
        setattr(
            self,
            '_' + var_name + '_data',
            Array(c_datatype, int(ncf_var.size), lock=False)
        )

        change_in_meters = False

        # Save the attributes
        var_attributes = {}
        for attrib in ncf_var.ncattrs():
            var_attributes[attrib] = ncf_var.getncattr(attrib)
            if attrib == 'units' and var_attributes[attrib] == 'km':
                var_attributes[attrib] = 'm'
                change_in_meters = True
        setattr(self, var_name + '_attributes', var_attributes)

        raw_data = np.array(ncf_var[:], dtype=dtype).flatten()
        if change_in_meters:
            raw_data *= 1000
            LOGGER.debug(
                'Variable "%s" has been converted from km to m and '
                'therefore has been multiplied by 1000',
                var_name
            )

        # This is a terrible hack, which I hope is also portable
        # The idea is to save inside the multiprocessing Array the data
        # saved inside the raw_data as they are, simply copying the bytes
        # and bypassing completely the numpy interface.
        # The idea is that the data of this class must be shared among all
        # the process that this script could spawn. Unfortunately, using
        # multiprocessing you can only share "C arrays". Numpy array are a
        # more compless objects with more features. Therefore, we use a C
        # array to dump all the data of the numpy arrays and then we rebuild
        # them in each process starting from the C array.
        # So, we use the multiprocessing Array just as a place where the
        # data can be dumped. It might happen that integers numbers where
        # saved in a float multiprocessing Array. Therefore the values that
        # you can find inside these vector are completely meaningless.
        # Nevertheless, when we rebuild the numpy array, the values should
        # be reconstructed in the right way
        mem_pointer = raw_data.ctypes.data_as(ctypes.POINTER(c_datatype))
        data_size = int(int(ncf_var.size))
        getattr(self, '_' + var_name + '_data')[:] = mem_pointer[:data_size]

        # Save the dimensions
        var_dimensions = []
        for dim_name in ncf_var.dimensions:
            dim_value = len(self.__ncf.dimensions[dim_name])
            new_dim_data = (dim_name, dim_value)
            var_dimensions.append(new_dim_data)
        setattr(self, var_name + '_dimensions', var_dimensions)

    @property
    def has_timebounds(self):
        return self._has_timebounds

    @property
    def times(self):
        raw_data = np.frombuffer(self._times_data, dtype=self._times_dtype)
        return raw_data.reshape(self._times_shape)

    @property
    def dates(self):
        """
        Dates is an array with the same content of times but in another format
        (netCDF4 dates instead of floats)
        """
        return num2date(
            self.times,
            units=self.times_attributes['units'],
            calendar=self.times_attributes['calendar']
        )

    @property
    def xlat(self):
        raw_data = np.frombuffer(self._xlat_data, dtype=self._xlat_dtype)
        return raw_data.reshape(self._xlat_shape)

    @property
    def xlon(self):
        raw_data = np.frombuffer(self._xlon_data, dtype=self._xlon_dtype)
        return raw_data.reshape(self._xlon_shape)

    @property
    def iy(self):
        raw_data = np.frombuffer(self._iy_data, dtype=self._iy_dtype)
        return raw_data.reshape(self._iy_shape)

    @property
    def jx(self):
        raw_data = np.frombuffer(self._jx_data, dtype=self._jx_dtype)
        return raw_data.reshape(self._jx_shape)

    @property
    def timebounds(self):
        if not self.has_timebounds:
            raise AttributeError(
                'Trying to access to the timebounds of a RegCM output file '
                'that does not have the timebounds'
            )
        else:
            raw_data = np.frombuffer(
                self._timebounds_data,
                dtype=self._timebounds_dtype,
            )
            return raw_data.reshape(self._timebounds_shape)

    @property
    def domain(self):
        return self._domain

    @property
    def product(self):
        return self._product

    @property
    def revision(self):
        return self._revision

    @property
    def revision_version(self):
        return self._rev_version

    @property
    def nesting_tag(self):
        return self._nest_tag

    @property
    def contains_map(self):
        return self._contains_map

    @property
    def map_type(self):
        if self.contains_map:
            return self._map_type
        else:
            raise AttributeError(
                'Trying to access to the type of the map of a RegCM output '
                'file that does not have a map'
            )

    @property
    def map(self):
        if self.contains_map:
            raw_data = np.frombuffer(self._map_data, dtype=self._map_dtype)
            return raw_data.reshape(self._map_shape)
        else:
            raise AttributeError(
                'Trying to access to the type of the map of a RegCM output '
                'file that does not have a map'
            )

    @property
    def grid_size(self):
        return self.attributes['grid_size_in_meters']

    @property
    def map_projection(self):
        return self.attributes['projection']

    @property
    def map_projection_latitude_origin(self):
        return self.attributes['latitude_of_projection_origin']

    @property
    def map_projection_longitude_origin(self):
        return self.attributes['longitude_of_projection_origin']

    @property
    def map_projection_longitude_origin(self):
        return self.attributes['longitude_of_projection_origin']

    @property
    def map_projection_standard_parallel(self):
        if 'standard_parallel' in self.attributes:
            return self.attributes['standard_parallel']
        else:
            return None

    @property
    def lock(self):
        return None

    def __str__(self):
        return self._name


class CordexDataset(Dataset):
    """
    This object represents one of the output file of this program. It inherits
    from the Dataset class of the NetCDF4 module, but, as soon as it is
    created, it also contains the default variables that are common in all the
    cortex files

    This class is intended only for writing files (the mode is forced to "w").

    :param output_file_path: a str that contains the path of the file that will
    be created
    :param regcm_file: an object from the RegcmOutputFile class
    :param simulation: an object from the RegCMSimulation class
    """

    def __init__(self, output_file_path, regcm_file, simulation):

        LOGGER.debug('Creating file %s', output_file_path)
        super(CordexDataset, self).__init__(
            output_file_path,
            mode='w',
            format='NETCDF4_CLASSIC'
        )

        if regcm_file.revision is not None:
            ICTP_Model = 'ICTP-RegCM{}'.format(regcm_file.revision)
            ICTP_Model_Version = regcm_file.revision_version
        else:
            LOGGER.warning('Using fallback values for RegCM version...')
            ICTP_Model = ICTP_Model_fallback
            ICTP_Model_Version = ICTP_Model_Version_fallback

        scenario = simulation.experiment.upper().replace('.', '')
        newattr = {
            'project_id': 'CORDEX',
            'ipcc_scenario_code': scenario,
            'institute_id': 'ICTP',
            'note': 'The domain is larger than ' + simulation.domain,
            'comment': 'RegCM CORDEX {} run'.format(regcm_file.domain),
            'experiment': simulation.domain,
            'experiment_id': simulation.experiment.replace('.', ''),
            'driving_experiment': '{}, {}, {}'.format(
                simulation.global_model,
                simulation.experiment.replace('.', ''),
                simulation.ensemble
            ),
            'driving_model_id': simulation.global_model,
            'driving_model_ensemble_member': simulation.ensemble,
            'driving_experiment_name': simulation.experiment.replace('.', ''),
            'institution': 'International Centre for Theoretical Physics',
            'model_id': ICTP_Model,
            'creation_date': time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                           time.localtime(time.time())),
            'CORDEX_domain': simulation.domain,
            'rcm_version_id': ICTP_Model_Version,
            'ICTP_version_note': simulation.notes,
            'contact': simulation.mail,
            'product': 'output',
            'Conventions': 'CF-1.7',
            'tracking_id': str(uuid.uuid1()),
        }

        if regcm_file.nesting_tag is not None:
            newattr['nesting_tag'] = regcm_file.nesting_tag

        file_name = path.basename(output_file_path)

        # Write Global attributes
        for attribute, value in newattr.items():
            LOGGER.debug(
                'Adding attribute "%s" with value "%s" to file %s',
                attribute,
                value,
                file_name
            )
            self.setncattr(attribute, value)

        for attribute, value in regcm_file.attributes.items():
            if attribute not in newattr.keys():
                LOGGER.debug(
                    'Adding attribute "%s" with value "%s" to file %s',
                    attribute,
                    value,
                    file_name
                )
                self.setncattr(attribute, value)
            else:
                LOGGER.debug(
                    'Attribute "%s" found in the RegCM file is not copied in '
                    'file %s because it is already set',
                    attribute,
                    file_name
                )

        xlat = self.create_var_from_data(
            REPLACE_VARS['xlat'],
            regcm_file.xlat,
            regcm_file.xlat_dimensions,
            attributes=regcm_file.xlat_attributes
        )
        #LOGGER.debug(
        #    'Adding attribute named "grid_mapping" with value "crs" to %s',
        #    REPLACE_VARS['xlat']
        #)
        #xlat.setncattr('grid_mapping', 'crs')

        xlon = self.create_var_from_data(
            REPLACE_VARS['xlon'],
            regcm_file.xlon,
            regcm_file.xlon_dimensions,
            attributes=regcm_file.xlon_attributes
        )
        #LOGGER.debug(
        #    'Adding attribute named "grid_mapping" with value "crs" to %s',
        #    REPLACE_VARS['xlon']
        #)
        #xlon.setncattr('grid_mapping', 'crs')

        x = self.create_var_from_data(
            REPLACE_VARS['jx'],
            regcm_file.jx,
            regcm_file.jx_dimensions,
            attributes=regcm_file.jx_attributes
        )
        LOGGER.debug(
            'To identify variable x as an axis, an attribute named "axis" will '
            'be added with value "X"'
        )
        x.setncattr('axis', 'X')

        y = self.create_var_from_data(
            REPLACE_VARS['iy'],
            regcm_file.iy,
            regcm_file.iy_dimensions,
            attributes=regcm_file.iy_attributes
        )
        LOGGER.debug(
            'To identify variable y as an axis, an attribute named "axis" will '
            'be added with value "Y"'
        )
        y.setncattr('axis', 'Y')

        if regcm_file.contains_map:
            self._copy_map_from_regcm_file(regcm_file)
        else:
            LOGGER.debug(
                'Map will not be generated because the regcm file does not '
                'contain it'
            )

    def create_var_from_data(self, var_name, data, dims,
                             datatype=DATATYPE_AUXILIARIES, attributes=None,
                             replace_dims=REPLACE_DIMS):
        """
        Create a new variable in the netcdf file starting from a numpy array
        and some metadata (like the name of the dimensions and the attributes).

        The dimensions that do not already exists (in the root group) will be
        created.

        If a dimension is unlimited, it will be saved as a fixed length
        dimension and the unlimited attributed will be lost.

        :param var_name: A string that will be used as the name of the new var
        :param data: A numpy array with the data that must be saved into the
        variable
        :param dims: A list of couples (tuples with two element). The first
        element of the tuple is the name of the dimension and the second one
        is an integer with its length. The order of the list must be such that
        [i[1] for i in dimensions] is the shape of data
        :param datatype: A string that represent the type of data of the
        variable. For example, "f8" means double precision
        :param attributes: A dictionary-like item
        :param replace_dims: A dictionary. If the name of a dimension is inside
        this dictionary, it will be replaced with its corresponding one. This is
        useful because some dimensions in RegCM are saved with another name in
        the CORDEX files (for example, jx becomes x).
        """
        LOGGER.debug('Saving variable %s', var_name)
        LOGGER.debug('dimensions %s', dims)

        if replace_dims is None:
            replace_dims = {}
        if attributes is None:
            attributes = {}

        # Copy dims
        dim_name_list = []
        for dim_name, dim_len in dims:
            LOGGER.debug(
                'A dimension named "%s" is required to create the variable %s',
                dim_name,
                var_name
            )
            if dim_name in replace_dims:
                LOGGER.debug(
                    'The name "%s" for a dimension is in the replace dict. It '
                    'will be called "%s" instead',
                    dim_name,
                    replace_dims[dim_name],
                )
                dim_name = replace_dims[dim_name]

            dim_name_list.append(dim_name)

            if dim_name in self.dimensions:
                LOGGER.debug(
                    'Dimension "%s" will not be created because it has already '
                    'been created on the file',
                    dim_name
                )
                dim_current_len = len(self.dimensions[dim_name])
                if dim_current_len != dim_len:
                    raise ValueError(
                        'The length of dimension {} is already set to {}. '
                        'To save variable "{}", it should be {}!'
                        .format(dim_name, dim_len, var_name, dim_current_len)
                    )
            else:
                LOGGER.debug(
                    'Creating dimension %s of length %s',
                    dim_name,
                    dim_len
                )
                if dim_name == 'time':
                    self.createDimension(dim_name, None)
                else:
                    self.createDimension(dim_name, dim_len)

        # Finding the fill value
        if '_FillValue' in attributes:
            fill_value = attributes['_FillValue']
            LOGGER.debug('Using %s as fill value (as requested)', fill_value)
        else:
            if is_masked(data):
                if datatype in NETCDF_DEFAULT_FILL_VALUES:
                    fill_value = NETCDF_DEFAULT_FILL_VALUES[datatype]
                    LOGGER.debug(
                        'Using %s as fill value (default value)',
                        fill_value
                    )
                else:
                    raise ValueError(
                        'Data is masked but not appropriated fill_value has '
                        'been found for datatype %s', datatype
                    )
            else:
                LOGGER.debug(
                    'No _FillValue specified and data is unmasked. Setting '
                    'fill_value flag as False'
                )
                fill_value = False

        # Create the variable
        LOGGER.debug(
            'Creating variable %s, with datatype "%s" and dimensions %s',
            var_name,
            datatype,
            tuple(dim_name_list)
        )
        ncdf_variable = self.createVariable(
            var_name,
            datatype,
            tuple(dim_name_list),
            fill_value=fill_value,
            zlib=COMPRESSION,
            complevel=COMPRESSION_LEVEL,
            shuffle=SHUFFLE,
            fletcher32=FLETCHER32
        )

        # Copy attributes
        for attr, attr_val in attributes.items():
            if attr == '_FillValue' or attr == 'missing_value':
                ncdf_variable.setncattr('missing_value', attr_val)
                continue

            if attr in EXCLUDED_ATTRIBUTES:
                LOGGER.debug(
                    'Avoiding to copy the attribute %s because it is in the '
                    'EXCLUDED_ATTRIBUTES list (file globals.py)', attr
                )
                continue

            LOGGER.debug(
                'Adding attribute "%s" with value "%s" for variable "%s"',
                attr,
                attr_val,
                var_name,
            )
            ncdf_variable.setncattr(attr, attr_val)

        if is_masked(data):
            LOGGER.debug('Data is masked')
        else:
            LOGGER.debug('Data is not masked')

        # Copy the values into the variable
        LOGGER.debug('Copying data inside the variable')
        ncdf_variable[:] = data

        return ncdf_variable

    def _copy_map_from_regcm_file(self, regcm_file):
        """
        Copy the crs (or regcm_map) variable into this cordex file.
        This function is a wrapper around the create_var_from_data function
        that also adds some metadata
        :param regcm_file: A RegCM file that contains a map
        """
        LOGGER.debug('Copying crs values into the CORDEX file')

        if not regcm_file.contains_map:
            raise ValueError('The RegCM file does not contain any map')

        crs_var = self.create_var_from_data(
            'crs',
            regcm_file.map,
            regcm_file.map_dimensions,
            datatype=regcm_file.map.dtype,
            attributes=regcm_file.map_attributes
        )

        if 'grid_mapping_name' not in regcm_file.map_attributes:
            LOGGER.debug('Setting attributes about the projection')

            crs_var.setncattr('semi_major_axis', 6371229.0)
            crs_var.setncattr('inverse_flattening', 0.0)
            crs_var.setncattr('false_easting', -regcm_file.grid_size / 2.0)
            crs_var.setncattr('false_northing', -regcm_file.grid_size / 2.0)

            if '_CoordinateTransformType' not in EXCLUDED_ATTRIBUTES:
                crs_var.setncattr('_CoordinateTransformType', 'Projection')
            if '_CoordinateAxisTypes' not in EXCLUDED_ATTRIBUTES:
                crs_var.setncattr('_CoordinateAxisTypes', 'GeoX GeoY')

            if regcm_file.map_projection == 'LAMCON':
                LOGGER.debug('Setting attributes for LAMCON proj')
                crs_var.setncattr(
                    'grid_mapping_name',
                    'lambert_conformal_conic'
                )
                crs_var.setncattr(
                    'standard_parallel',
                    regcm_file.map_projection_standard_parallel
                )
                crs_var.setncattr(
                    'latitude_of_projection_origin',
                    regcm_file.map_projection_latitude_origin,
                )
                crs_var.setncattr(
                    'false_easting',
                    -regcm_file.grid_size / 2.0,
                )
                crs_var.setncattr(
                    'false_northing',
                    -regcm_file.grid_size / 2.0,
                )
            elif regcm_file.map_projection == 'POLSTR':
                LOGGER.debug('Setting attributes for POLSTR proj')
                crs_var.setncattr(
                    'grid_mapping_name',
                    'stereographic'
                )
                crs_var.setncattr(
                    'latitude_of_projection_origin',
                    regcm_file.map_projection_latitude_origin,
                )
                crs_var.setncattr(
                    'longitude_of_projection_origin',
                    regcm_file.map_projection_longitude_origin,
                )
                crs_var.setncattr(
                    'scale_factor_at_projection_origin',
                    1.0
                )
                crs_var.setncattr(
                    'false_easting',
                    -regcm_file.grid_size / 2.0,
                )
                crs_var.setncattr(
                    'false_northing',
                    -regcm_file.grid_size / 2.0,
                )
            elif regcm_file.map_projection == 'NORMER':
                LOGGER.debug('Setting attributes for NORMER proj')
                crs_var.setncattr(
                    'grid_mapping_name',
                    'mercator'
                )
                crs_var.setncattr(
                    'standard_parallel',
                    regcm_file.map_projection_latitude_origin,
                )
                crs_var.setncattr(
                    'latitude_of_projection_origin',
                    regcm_file.map_projection_latitude_origin,
                )
                crs_var.setncattr(
                    'longitude_of_projection_origin',
                    regcm_file.map_projection_longitude_origin,
                )
                crs_var.setncattr(
                    'false_easting',
                    -regcm_file.grid_size / 2.0,
                )
                crs_var.setncattr(
                    'false_northing',
                    -regcm_file.grid_size / 2.0,
                )
            elif regcm_file.map_projection == 'ROTMER':
                LOGGER.debug('Setting attributes for ROTMER proj')
                crs_var.setncattr(
                    'grid_mapping_name',
                    'oblique_mercator'
                )
                crs_var.setncattr(
                    'azimuth_of_central_line',
                    89.999999,
                )
                crs_var.setncattr(
                    'latitude_of_projection_origin',
                    regcm_file.map_projection_latitude_origin
                )
                crs_var.setncattr(
                    'longitude_of_projection_origin',
                    regcm_file.map_projection_longitude_origin
                )
                crs_var.setncattr(
                    'scale_factor_at_projection_origin',
                    1.0
                )
                crs_var.setncattr(
                    'false_easting',
                    -regcm_file.grid_size / 2.0,
                )
                crs_var.setncattr(
                    'false_northing',
                    -regcm_file.grid_size / 2.0,
                )
            else:
                raise ValueError(
                    'Unknown projection type: {}'
                    .format(regcm_file.map_projection)
                )
        else:
            # It is expected that copying the file attributes all the
            # information about the projection will be also copied (if it was
            # saved in the original RegCM file)
            LOGGER.debug(
                'semi_major_axis attribute is already present in the RegCM '
                'file. Skipping the procedure that set the information about '
                'the projection'
            )

    def save_cordex_variable(self, cordex_var):
        attributes = copy(cordex_var.attributes)

        LOGGER.debug('Setting parameter "frequency" for the netCDF file')
        self.setncattr('frequency', cordex_var.frequency)

        if 'coordinates' in attributes:
            # The "coordinates" attribute refers to variables that could have
            # a different name now. Therefore, it must be updated
            LOGGER.debug('Updating "coordinates" attribute (still in memory)')
            for dim_old_name, dim_new_name in REPLACE_VARS.items():
                coord_val = attributes['coordinates']
                LOGGER.debug(
                    'Replacing any reference to "%s" with "%s" for "%s"',
                    dim_old_name,
                    dim_new_name,
                    coord_val
                )
                attributes['coordinates'] = re.sub(
                    r'\b{}\b'.format(dim_old_name),
                    dim_new_name,
                    coord_val
                )

        if cordex_var.depends_on_time:
            LOGGER.debug('Saving variable "time"')
            time_attributes = cordex_var.times_attributes
            if cordex_var.needs_time_bounds:
                time_attributes['bounds'] = 'time_bnds'

            times = cordex_var.times
            time_step = cordex_var.time_step_size
            if cordex_var.time_step_size >= 24:
                LOGGER.debug('Switching time units to day instead of hours')
                times = times / 24
                time_step /= 24
                if 'units' in time_attributes:
                    new_val = time_attributes['units'].replace('hours', 'days')
                    LOGGER.debug(
                        'Replaced attr "units" in "time" from "%s" to "%s"',
                        time_attributes['units'],
                        new_val
                    )
                    time_attributes['units'] = new_val

            self.create_var_from_data(
                'time',
                times,
                (('time', cordex_var.times.size),),
                attributes=time_attributes,
            )
        else:
            LOGGER.debug(
                'Variable "time" not created (the cordex var does not depend '
                'on time)'
            )

        if cordex_var.needs_time_bounds:
            LOGGER.debug('Creating time bounds variable')
            dims = (('time', times.size), ('bnds', 2))
            bounds = np.empty([d[1] for d in dims], dtype=times.dtype)
            bounds[:, 0] = times - (time_step / 2)
            bounds[:, 1] = times + (time_step / 2)

            bounds_attributes = {
                'units': time_attributes['units'],
                'calendar': time_attributes['calendar']
            }
            self.create_var_from_data(
                'time_bnds',
                bounds,
                dims,
                attributes=bounds_attributes
            )
        else:
            LOGGER.debug('time_bnds creation skipped because it was not needed')

        for auxiliary_var in cordex_var.auxiliary_variables:
            LOGGER.debug('Saving variable %s', auxiliary_var.name)
            with auxiliary_var.data:
                auxiliary_data = auxiliary_var.data()

            self.create_var_from_data(
                auxiliary_var.name,
                auxiliary_data,
                auxiliary_var.dimensions,
                attributes=auxiliary_var.attributes,
                datatype=auxiliary_data.dtype,
            )

        with cordex_var.data:
            data = cordex_var.data()

        self.create_var_from_data(
            cordex_var.name,
            data,
            cordex_var.dimensions,
            attributes=attributes,
            datatype=data.dtype,
        )


def get_regcm_variable_names(var_name):
    LOGGER.debug(
        'Checking if "%s" is a known variable of the regcm files',
        var_name
    )

    if var_name in REGCM_VARS:
        regcm_variables = REGCM_VARS[var_name]
        LOGGER.debug(
            'Variable "%s" is known! Its names are: %s',
            var_name,
            ', '.join(regcm_variables)
        )
    else:
        regcm_variables = [var_name, ]
        LOGGER.debug(
            '"%s" is not in the regcm_vars file. Trying to open a variable '
            'with the same name',
            var_name
        )

    if len(regcm_variables) == 0:
        raise KeyError(
            'Invalid regcm_vars file! No available names for variable "{}"'
            .format(var_name)
        )
    return regcm_variables


def read_variable_from_regcmfile(var_name, regcm_file, regcm_file_path,
                                 need_time_bounds=False):

    regcm_variables = get_regcm_variable_names(var_name)

    LOGGER.debug('Opening file %s', regcm_file_path)
    with Dataset(regcm_file_path, 'r') as ncdf_file:
        var_in_file = None
        for var_alternative in regcm_variables:
            LOGGER.debug(
                'Trying to read variable "%s" using name "%s"',
                var_name,
                var_alternative
            )
            if var_alternative in ncdf_file.variables:
                LOGGER.debug(
                    'Var found! Using "%s" as var name',
                    var_alternative
                )
                var_in_file = var_alternative
                break
            else:
                LOGGER.debug('Variable "%s" not found', var_name)
        else:
            raise KeyError(
                'Variable "{}" not found in file {} (tried with names: {})'
                .format(
                    var_name,
                    regcm_file_path,
                    ', '.join(regcm_variables)
                )
            )

        regcm_var = ncdf_file.variables[var_in_file]

        LOGGER.debug('Reading dimensions')
        dimensions = []
        for dim_name in regcm_var.dimensions:
            dim_value = len(ncdf_file.dimensions[dim_name])
            LOGGER.debug(
                'Found dimension "%s" of length %s',
                dim_name,
                dim_value
            )
            dimensions.append((dim_name, dim_value))

        LOGGER.debug('Read attributes')
        attributes = {}
        for attr_name in regcm_var.ncattrs():
            attr_value = regcm_var.getncattr(attr_name)
            attributes[attr_name] = attr_value
            LOGGER.debug(
                'Read attribute "%s" with values "%s"',
                attr_name,
                attr_value
            )

    # Create a callable object that can be retrieve the data of this
    # variable
    data = NetCDFData(
        var_in_file,
        regcm_file_path,
        DATATYPE_MAIN,
        lock=regcm_file.lock
    )

    if 'time' in [i[0] for i in dimensions]:
        times = regcm_file.times
        LOGGER.debug('Copying time from the RegCM file for the variable')
    else:
        LOGGER.debug(
            'Time not saved inside the variable (no dimension "time" found'
        )
        times = None

    read_var = Variable(
        name=var_name,
        data=data,
        dimensions=dimensions,
        attributes=attributes,
        times=times,
        needs_time_bounds=need_time_bounds,
        times_attributes=regcm_file.times_attributes.copy(),
    )

    LOGGER.debug('Variable read')
    return read_var


def prepare_cordex_file_dir(var_name, var_dates, var_freq, simul, regcm_file,
                            cordex_dir):
    """
    Create the directory tree needed to save the new NetCDf file with the CORDEX
    variable. Moreover, return the path of the CORDEX file.

    :param var_name: The name of the CORDEX variable as a string
    :param var_dates: an array of netCDF4 dates for the variable
    :param: var_freq: a string that represent the time step of the variable
    (like "3h" or "day")
    :param simul: An object of the Simulation class
    :param cordex_dir: The main dir where the CORDEX files will be saved
    :return: the path where the new cordex file will be created
    """

    if regcm_file.revision is not None:
        ICTP_Model = 'ICTP-RegCM{}'.format(regcm_file.revision)
        ICTP_Model_Version = regcm_file.revision_version
    else:
        LOGGER.warning('Using fallback values for RegCM version...')
        ICTP_Model = ICTP_Model_fallback
        ICTP_Model_Version = ICTP_Model_Version_fallback

    # This is the directory where the NetCDf file will be saved
    cordex_last_dir = os.path.join(
        cordex_dir,
        simul.product,
        simul.domain,
        'ICTP',
        simul.global_model,
        simul.experiment.translate({None: '.'}),
        simul.ensemble,
        ICTP_Model,
        ICTP_Model_Version,
        var_freq,
        var_name
    )
    LOGGER.debug('Checking if directory {} exists'.format(cordex_last_dir))
    if os.path.exists(cordex_last_dir):
        LOGGER.debug('Path already exists. Checking if it is a directory')
        if os.path.isdir(cordex_last_dir):
            LOGGER.debug('Path "%s" exists and is a directory', cordex_last_dir)
        else:
            raise IOError(
                'Path "{}" is not a directory!'.format(cordex_last_dir)
            )
    else:
        LOGGER.debug('Path does not exist. Creating it...')
        try:
            os.makedirs(cordex_last_dir)
        except FileExistsError:
            pass
        LOGGER.debug('Path created!')

    LOGGER.debug('Preparing strings for the first and last date of the var')
    dd1, dd2 = get_first_and_last_date_str(var_dates, var_freq)
    LOGGER.debug('First date: %s    Last date: %s', dd1, dd2)

    LOGGER.debug('Creating netCDF file basename')
    cordex_basename = '_'.join([
        var_name,
        simul.domain,
        simul.global_model,
        simul.experiment.translate({None: '.'}),
        simul.ensemble,
        ICTP_Model,
        ICTP_Model_Version,
        var_freq,
        dd1 + '-' + dd2 + '.nc',
    ])
    LOGGER.debug('The new file will be called "%s"', cordex_basename)

    return os.path.join(cordex_last_dir, cordex_basename)


def mrtosph(x, oname):
    if oname == 'humidity_mixing_ratio':
        sph = 1.0 / (1.0 - x)
    else:
        sph = x
    return sph


def correct_evaporation(x):
    corrected = np.where(x < 0.0, 0.0, x)
    return corrected


def unit_factor_correction(old, new):
    mm_per_day_units = (
        'mm/day', 'mm day-1', 'Kg m-2 day-1', 'kg m-2 day-1', 'kg m-2 d-1'
    )

    # Convert rain units
    if old in mm_per_day_units and new == 'kg m-2 s-1':
        return 1 / 86400.0

    # convert pressure units
    elif old == 'mb' and new == 'Pa':
        return 100.0

    elif old == 'kg m-2' and new == 'm':
        return 1 / 1000.0

    elif old == 'kg m-2 s-1' and new == 'W m-2':
        return 2.50080e6

    return 1.0
