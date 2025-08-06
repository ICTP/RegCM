from netCDF4 import Dataset, num2date
from logging import getLogger
from traceback import format_exc
import numpy as np
from numpy.ma import is_masked, masked_array
import numbers
import sys

from utilities.globals import REGCM_VARS


__copyright__ = 'Copyright (C) 2017-2018 ICTP'
__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__credits__ = ["Stefano Piani"]


LOGGER = getLogger(__name__)


def step_postponer(slice_data):
    """
    Remove the steps from the slices, and put them in another slices.
    In other words, return a couple sl1 and sl2 such that arr[slice_data] is
    equal to arr[sl1][sl2] and the slice of sl1 do not have steps while sl2
    has only steps and not start and stop (therefore is like [::step])

    Check https://github.com/Unidata/netcdf4-python/issues/680 for the
    motivation

    :param slice_data: a slice or an iterable of slices
    :return: a tuple of slices
    """
    slice1 = []
    slice2 = []
    used_step = None

    try:
        slice_list = [s for s in slice_data]
    except TypeError:
        slice_list = [slice_data]

    for s in slice_list:
        if isinstance(s, int):
            slice1.append(s)
            slice2.append(np.index_exp[:][0])
            continue

        if s.step is not None:
            used_step = True

        slice1.append(np.index_exp[s.start:s.stop][0])
        slice2.append(np.index_exp[::s.step][0])

    if used_step:
        return slice1, slice2
    else:
        return slice1, None


class SliceStr(object):
    """
    A wrapper around an object like a tuple of slices, a list of slices or a
    slice.
    It allows to represent the object in a human readable format.
    """
    def __init__(self, slice_data):
        self._slice = slice_data

    @staticmethod
    def _slice_str(slice_obj):
        if isinstance(slice_obj, int):
            return str(slice_obj)
        output = ''
        start = slice_obj.start
        stop = slice_obj.stop
        step = slice_obj.step
        if start is None and stop is None:
            output += 'All elements'
        else:
            if start is None:
                output += 'From beginning'
            else:
                output += 'From element {}'.format(start)
            if stop is None:
                output += ' to end'
            else:
                output += ' to element {}'.format(stop)
        if step is not None:
            output += ' using a step of {}'.format(step)
        return output

    def __str__(self):
        try:
            slice_str_list = [SliceStr._slice_str(s) for s in self._slice]
        except TypeError:
            return SliceStr._slice_str(self._slice)

        return '[' + ', '.join(slice_str_list) + ']'


class Variable(object):
    """
    A Variable is an object that can save all the information that a NetCDF
    variable inside a file can store

    """
    def __init__(self, name, data, dimensions, attributes, times=None,
                 times_attributes=None, auxiliary_vars=None,
                 needs_time_bounds=False):
        self.name = name
        self._data = data
        self._dimensions = dimensions
        self._attributes = attributes
        self._times = times
        if times_attributes is None:
            self._times_attributes = {}
        else:
            self._times_attributes = times_attributes
        if auxiliary_vars is None:
            self._auxiliary_variables = []
        else:
            self._auxiliary_variables = auxiliary_vars

        # In some situation override the flag needs_time_bounds
        if self.time_step_size == 24 or self.time_step_size > 670:
            self._needs_time_bounds = True
        elif not self.depends_on_time:
            self._needs_time_bounds = False
        else:
            self._needs_time_bounds = needs_time_bounds

    @property
    def data(self):
        return self._data

    @property
    def dimensions(self):
        return self._dimensions

    @property
    def attributes(self):
        return self._attributes

    @property
    def auxiliary_variables(self):
        return self._auxiliary_variables

    @property
    def depends_on_time(self):
        if self._times is None:
            return False
        else:
            return True

    @property
    def needs_time_bounds(self):
        return self._needs_time_bounds

    @property
    def times(self):
        return self._times

    @property
    def times_attributes(self):
        return self._times_attributes

    @property
    def time_step_size(self):
        if not self.depends_on_time:
            return -1
        if len(self.times) >= 2:
            return self.times[1] - self.times[0]
        return 24.0

    @property
    def frequency(self):
        if not self.depends_on_time:
            return 'fx'

        if self.time_step_size == 24:
            return 'day'

        if self.time_step_size > 670:
            return 'month'

        return '{}hr'.format(int(self.time_step_size))

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

    def __str__(self):
        return 'variable {}'.format(self.name)


class Data(object):
    def __init__(self):
        raise NotImplementedError

    def __enter__(self):
        return self

    def __exit__(self, dtype, value, tb):
        pass


class NetCDFData(Data):
    """
    This class represent data that must be read from a single variable into a
    NetCDf file
    """
    def __init__(self, var_name, netcdf_file_path, datatype, lock=None):
        self.var_name = var_name
        self._path = netcdf_file_path
        self.__f_pointer = None
        self.__datatype = datatype
        self.__lock = lock

    def __enter__(self):
        if self.__f_pointer is not None:
            raise IOError('Trying to open again file {}'.format(self._path))

        if self.__lock is not None:
            LOGGER.debug('Acquiring lock on file')
            self.__lock.acquire()
            LOGGER.debug('Lock acquired')

        self.__f_pointer = Dataset(self._path, 'r')
        return self

    def __call__(self, slice_data=None):
        if self.__f_pointer is None:
            raise IOError('File {} is still closed'.format(self._path))

        if slice_data is None:
            slice_data = np.index_exp[:]

        LOGGER.debug(
            'Reading data for variable %s from file %s (slice: %s)',
            self.var_name,
            self._path,
            SliceStr(slice_data)
        )
        data_var = self.__f_pointer.variables[self.var_name]

        # The data_raw has the same data type of the netCDF file and could be
        # a masked array
        slice_start_stop, slice_steps = step_postponer(slice_data)
        if slice_steps is not None:
            data_raw = data_var[tuple(slice_start_stop)][tuple(slice_steps)]
        else:
            data_raw = data_var[tuple(slice_start_stop)]
        LOGGER.debug('Data read')

        # Now we cast the data
        data = np.array(data_raw, dtype=self.__datatype)

        if not is_masked(data_raw):
            return data
        else:
            LOGGER.debug(
                'Some values are masked (because they are outside the valid '
                'range or because they coincide with the fill_value). These '
                'values will set to nan (if datatype is float)'
            )
            data_mask = data_raw.mask

            # If data is float, put nan on the masked values
            data_flat = data.reshape((data.size,))
            data_element = data_flat[0]
            if isinstance(data_element, numbers.Real):
                data[data_mask] = np.nan

            return masked_array(data, mask=data_mask)

    def __exit__(self, dtype, value, tb):
        if self.__f_pointer is not None:
            try:
                self.__f_pointer.close()
            except:
                LOGGER.warning(
                    'Unable to close file {}:\n{}'.format(
                        self._path,
                        format_exc()
                    )
                )
        if self.__lock is not None:
            LOGGER.debug('Releasing file lock')
            self.__lock.release()

        self.__f_pointer = None

    @property
    def source(self):
        return 'disk'


class MemoryData(Data):
    """
    This class represent data that are saved in memory as a numpy array (or a
    numpy masked array)
    """
    def __init__(self, data_array):
        self.__data_array = data_array

    def __call__(self, slice_data=None):
        if slice_data is None:
            slice_data = np.index_exp[:]
        return self.__data_array[tuple(slice_data)]

    @property
    def source(self):
        return 'memory'


class FilterData(Data):
    """
    This class represent data that are elaborated starting from another Data
    object by calling a function on them
    """
    def __init__(self, prev_data, f, f_preserves_data_shape=False):
        self.__prev_data = prev_data
        self.__f = f
        self.preserve_data_shape = f_preserves_data_shape

    def __enter__(self):
        self.__prev_data.__enter__()
        return self

    def __call__(self, slice_data=None):
        if self.preserve_data_shape:
            return self.__f(self.__prev_data(slice_data))
        else:
            if slice_data is None:
                slice_data = np.index_exp[:]
            return self.__f(self.__prev_data())[tuple(slice_data)]

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__prev_data.__exit__(exc_type, exc_val, exc_tb)

    @property
    def source(self):
        return self.__prev_data.source


class SubselectedData(Data):
    """
    This class represent data that are elaborated starting from another Data
    object by executing a subselection (with a slice) on them
    """
    def __init__(self, prev_data, data_slice):
        self.__prev_data = prev_data
        self.__slice = data_slice

    def __enter__(self):
        self.__prev_data.__enter__()
        return self

    def __call__(self, slice_data=None):
        if slice_data is None:
            slice_data = np.index_exp[:]
        return self.__prev_data(self.__slice)[tuple(slice_data)]

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__prev_data.__exit__(exc_type, exc_val, exc_tb)

    @property
    def source(self):
        return self.__prev_data.source



class SliceData(Data):
    """
    This class represent data that are elaborated starting from another Data
    object by reading a slice where a specific dimension is fixed
    """
    def __init__(self, prev_data, dimension, fixed_index):
        self.__prev_data = prev_data
        self.dimension = dimension
        self.index = fixed_index

    def __enter__(self):
        self.__prev_data.__enter__()
        return self

    def __call__(self, slice_data=None):
        if slice_data is None:
            slice_data = np.index_exp[:]

        try:
            len(slice_data)
        except TypeError:
            slice_data = (slice_data,)

        slice_data = list(slice_data)
        while len(slice_data) < self.dimension - 1:
            slice_data.append(slice(None))

        new_slice_data = slice_data[:self.dimension] + [self.index] + \
                         slice_data[self.dimension + 1:]

        LOGGER.debug(
            'Transforming slice from %s to %s',
            SliceStr(slice_data),
            SliceStr(new_slice_data)
        )

        return self.__prev_data(new_slice_data)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__prev_data.__exit__(exc_type, exc_val, exc_tb)

    @property
    def source(self):
        return self.__prev_data.source


class NetCDFMultiData(Data):
    """
    This class represent data that must be read from several different variables
    from one single NetCDf file and must be combined using a function
    """
    def __init__(self, var_names, combiner, netcdf_file_path, datatype,
                 lock=None, slice_as_first_argument=False):
        self.var_names = var_names
        self.combiner = combiner
        self._path = netcdf_file_path
        self.__f_pointer = None
        self.__datatype = datatype
        self.__lock = lock
        self.__slice_as_first_argument=slice_as_first_argument

    def __enter__(self):
        if self.__f_pointer is not None:
            raise IOError('Trying to open again file {}'.format(self._path))

        if self.__lock is not None:
            LOGGER.debug('Acquiring lock on file')
            self.__lock.acquire()
            LOGGER.debug('Lock acquired')

        self.__f_pointer = Dataset(self._path, 'r')
        return self

    def __call__(self, slice_data=None):
        if self.__f_pointer is None:
            raise IOError('File {} is still closed'.format(self._path))

        if slice_data is None:
            slice_data = np.index_exp[:]

        regcm_vars = []
        for var_name in self.var_names:
            LOGGER.debug(
                'Reading variable %s from file %s',
                var_name,
                self._path
            )
            if var_name in REGCM_VARS:
                var_name_list = REGCM_VARS[var_name]
            else:
                var_name_list = [var_name,]
            regcm_vars.append(
                get_var_with_name(var_name_list, self.__f_pointer)
            )

        combiner = self.combiner
        if self.__slice_as_first_argument:
            output_val = np.array(
                combiner(*[slice_data] + [v[slice_data] for v in regcm_vars]),
                dtype=self.__datatype
            )
        else:
            output_val = np.array(
                combiner(*[v[slice_data] for v in regcm_vars]),
                dtype=self.__datatype
            )

        return output_val

    def __exit__(self, dtype, value, tb):
        if self.__f_pointer is not None:
            try:
                self.__f_pointer.close()
            except:
                LOGGER.warning(
                    'Unable to close file {}:\n{}'.format(
                        self._path,
                        format_exc()
                    )
                )

        if self.__lock is not None:
            LOGGER.debug('Releasing file lock')
            self.__lock.release()

        self.__f_pointer = None


def get_var_with_name(variable_names, netcdf_file_pointer):
    """
    Sometimes, the same variable can be saved with different names among several
    netcdf files.
    This function takes a netcdf file and a list of strings as input and
    returns a pointer to a variable that is called with the first string in the
    iterable for which exists a variable with that name.

    :param variable_names: a iterable of strings
    :param netcdf_file_pointer: an open Dataset object from the library netCDF4
    :return: A netcdf variable
    """
    for var_name in variable_names:
        try:
            LOGGER.debug('Trying to read variable %s', var_name)
            netcdf_var = netcdf_file_pointer.variables[var_name]
            LOGGER.debug('Read success!')
            return netcdf_var
        except KeyError:
            LOGGER.debug('No variable %s. Retry again', var_name)

    LOGGER.error('No other names for the variable available. Read failed.')
    raise KeyError(
        'No NetCDF variable found with one of these names: ' +
        ', '.join(['"' + v + '"' for v in variable_names])
    )
