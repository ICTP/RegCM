import re
import numpy as np
from numpy.ma import is_masked
import logging
from importlib import import_module
from inspect import isclass, getmembers
from copy import copy
from netCDF4 import Dataset

from utilities.globals import DATATYPE_MAIN, DATATYPE_AUXILIARIES
from utilities.rotate import grid_to_earth_uvrotate
from utilities.repeater import Repeater
from utilities.variables import Variable, MemoryData, FilterData,\
    SubselectedData, SliceData, NetCDFMultiData, SliceStr
from utilities.cordex_utils import CordexDataset, prepare_cordex_file_dir, \
    read_variable_from_regcmfile, get_regcm_variable_names, get_var_with_name, \
    unit_factor_correction

# f2py modules
from utilities.vertint import mod_vertint
from utilities.hgt import mod_hgt
from utilities.capecin import mod_capecin


__copyright__ = 'Copyright (C) 2018 ICTP'
__author__ = 'Stefano Piani <stefano.piani@exact-lab.it>'
__credits__ = ["Stefano Piani", "Graziano Giuliani"]


LOGGER = logging.getLogger(__name__)

class fake_var(np.ndarray):

    def __new__(cls, input_array, dimensions=[None,], name='noname'):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.dimensions = dimensions
        obj.name = name
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self.dimensions = getattr(obj, 'dimensions', [None,])
        self.name = getattr(obj, 'name', 'noname')

class Method(object):
    """A Method is a callable object that performs actions on a variable.

    It is invoked in the JSON file that describes the actions to perform on a
    variable.

    This class is an abstract class, that make sure that at least the important
    methods have been set in the descendant classes.
    """

    def __init__(self, *args, **kwargs):
        raise NotImplementedError

    def __call__(self, *args, **kwargs):
        raise NotImplementedError


class Action(object):
    """An Action is a chain of Methods that retrieves a variable.

    Each Method (beside the first one) receives the output of its ancestor.
    """

    def __init__(self, methods_list, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):
        self.simulation = simulation
        self.regcm_file = regcm_file
        self.regcm_file_path = regcm_file_path
        self._methods = methods_list
        self.corrflag = corrflag
        self.cordex_dir = cordex_dir
        if not isinstance(self._methods[0], ActionStarter):
            raise ValueError(
                'The first Method of the list must be an ActionStarter'
            )

    def execute(self):
        LOGGER.debug('Executing method %s', self._methods[0].__class__.__name__)
        current_output = self._methods[0](
            regcm_file=self.regcm_file,
            regcm_file_path=self.regcm_file_path,
            simulation=self.simulation,
            corrflag=self.corrflag,
            cordex_dir=self.cordex_dir,
        )
        LOGGER.debug('Method %s executed', self._methods[0].__class__.__name__)
        for method in self._methods[1:]:
            LOGGER.debug('Executing method %s', method.__class__.__name__)
            current_output = method(
                current_output,
                regcm_file=self.regcm_file,
                regcm_file_path=self.regcm_file_path,
                simulation=self.simulation,
                corrflag=self.corrflag,
                cordex_dir=self.cordex_dir
            )
            LOGGER.debug('Method %s executed', method.__class__.__name__)


class ActionStarter(Method):
    """An ActionStarter does not require any input from another method."""

    def __init__(self, *args, **kwargs):
        raise NotImplementedError

    def __call__(self, regcm_file, regcm_file_path, simulation, corrflag,
                 cordex_dir):
        raise NotImplementedError


class Filter(Method):
    """A Filter requires as input the output of another method to execute."""

    def __init__(self, *args, **kwargs):
        raise NotImplementedError

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):
        raise NotImplementedError


def get_methods():
    """Retrieve the classes inheriting from Method, stored in the
    ``pycordexer_methods`` module.

    Returns:

        dict: a dictionary that associates each name (as a string) with its
        corresponding class.
    """
    methods_avail = {}

    methods_module = import_module('utilities.pycordexer_methods')

    for obj_name, obj in getmembers(methods_module):
        if obj in (Method, Filter, Action, ActionStarter):
            continue
        if isclass(obj) and issubclass(obj, Method):
            methods_avail[obj_name] = obj

    return methods_avail


# From now on there are the real methods that can be used in the cordex_vars
# JSON files


class ReadVariableFromFile(ActionStarter):
    """Read a variable from a file.

    Parameters:

        var_name (string): the name of the variable to extract from the file.

        need_time_bounds (boolean, optional): whether to use an auxiliary
            variable to store the interval range. Defaults to ``False``.
    """

    def __init__(self, var_name, need_time_bounds=False):
        self.var_name = var_name
        self.need_time_bounds = need_time_bounds

    def __call__(self, regcm_file, regcm_file_path, simulation, corrflag,
                 cordex_dir):
        return read_variable_from_regcmfile(
            self.var_name,
            regcm_file,
            regcm_file_path,
            self.need_time_bounds
        )


class SaveVariableToDisk(Filter):
    """Save a variable to disk.

    Parameters:

        var_name (string, optional): if provided, it overrides the name of the
            variable stored in the previous Method of the chain.

        fill_value (float, optional): if provided, fill the variable with this
            value.

        new_attributes (dict ``string -> object``, optional): if provided,
            contains attributes that can override those from the previous
            Method object. netCDF takes care of applying the object value
            correctly.
    """

    def __init__(self, var_name=None, fill_value=None, new_attributes=None):
        self.var_name = var_name
        if fill_value is None:
            self.fill_value = None
        else:
            # Save the fill value in the same format that is used to save the
            # main variable. The conversion is made creating an array with just
            # one value (the fill_value) and then retrieving it back from the
            # array
            self.fill_value = np.array([fill_value, ], dtype=DATATYPE_MAIN)[0]
        if new_attributes is None:
            self.attributes = {}
        else:
            self.attributes = new_attributes

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):

        var_name = prev_result.name
        if self.var_name is not None:
            var_name = self.var_name

        cordex_path = prepare_cordex_file_dir(
            var_name,
            regcm_file.dates,
            prev_result.frequency,
            simulation,
            regcm_file,
            cordex_dir
        )

        attributes = copy(prev_result.attributes)

        # Change the units. If needed, multiply the data for the right
        # conversion number
        new_data = prev_result.data
        if 'units' in attributes and 'units' in self.attributes:
            old_units = attributes['units']
            new_units = self.attributes['units']
            if old_units != new_units:
                LOGGER.debug(
                    'Changing units from "%s" to "%s"',
                    old_units,
                    new_units
                )
                units_factor = unit_factor_correction(old_units, new_units)
                LOGGER.debug('The correction factor is %s', units_factor)
                new_data = FilterData(
                    prev_result.data,
                    lambda x: x * units_factor,
                    f_preserves_data_shape=True
                )
            else:
                LOGGER.debug('No conversion applied for the unit change')

        for attr_name, attr_value in self.attributes.items():
            LOGGER.debug(
                'Preparing to add attribute "%s" with value "%s" to var "%s"',
                attr_name,
                attr_value,
                var_name
            )
            attributes[attr_name] = attr_value

        if self.fill_value is not None:
            attributes['_FillValue'] = self.fill_value

        with new_data:
            cordex_var = Variable(
                name=var_name,
                data=new_data,
                dimensions=prev_result.dimensions,
                attributes=attributes,
                times=prev_result.times,
                times_attributes=prev_result.times_attributes,
                auxiliary_vars=prev_result.auxiliary_variables,
                needs_time_bounds=prev_result.needs_time_bounds,
            )

        LOGGER.info('Writing on file %s', cordex_path)
        with CordexDataset(cordex_path, regcm_file, simulation) as cdf:
            cdf.save_cordex_variable(cordex_var)
        LOGGER.info('Finished writing on file %s', cordex_path)

        return prev_result


class SaveMultipleVariablesToDisk(Filter):
    """Save multiple variables to disk.

    Parameters:

        fill_value (float, optional): if provided, fill the variable with this
            value.

        new_attributes (dict ``string -> object``, optional): if provided,
            contains attributes that can override those from the previous
            Method object. netCDF takes care of applying the object value
            correctly.
    """

    def __init__(self, fill_value=None, new_attributes=None):

        self.__writer = SaveVariableToDisk(
            fill_value=fill_value,
            new_attributes=new_attributes
        )

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):

        for cordex_var in prev_result:
            LOGGER.debug('Saving variable %s', cordex_var.name)

            self.__writer(
                cordex_var,
                regcm_file,
                regcm_file_path,
                simulation,
                corrflag,
                cordex_dir
            )

        return prev_result


class CorrectTime(Filter):
    """Apply an offset to correct the time of the variable.

    Parameters:

        offset (dict): a dict containing an ``ftype`` field to apply as offset.
    """

    def __init__(self, offset):
        self._offset = offset

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):
        LOGGER.debug('Executing time correction')
        if corrflag is False:
            LOGGER.debug('Time correction skipped because of the corrflag')
            return prev_result

        ftype = regcm_file.ftype
        LOGGER.debug('Looking for %s in the offset dictionary', ftype)
        if ftype not in self._offset:
            raise KeyError(
                'No offset value found that is suitable for {} files'
                .format(ftype)
            )

        raw_offset = str(self._offset[ftype])
        sub_offset = raw_offset.replace(
            'HALFFREQUENCY',
            str(prev_result.time_step_size / 2.)
        )
        offset = float(sub_offset)
        LOGGER.debug('Offset value is %s', offset)

        LOGGER.debug('Applying offset')
        correct_var = Variable(
            name=prev_result.name,
            data=prev_result.data,
            dimensions=prev_result.dimensions,
            attributes=prev_result.attributes,
            times=prev_result.times + offset,
            needs_time_bounds=prev_result.needs_time_bounds,
            times_attributes=prev_result.times_attributes,
            auxiliary_vars=prev_result.auxiliary_variables,
        )

        return correct_var


class AddDephtVariable(Filter):
    """Add the ground depth, discarding empty data."""

    def __init__(self):
        pass

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):

        with Dataset(regcm_file_path, 'r') as f:
            soil_names = get_regcm_variable_names('soil_layer')
            bound_names = get_regcm_variable_names('soil_bounds')
            LOGGER.debug('Reading soil_layer variable from regcm file')
            try:
                soil_layer = get_var_with_name(soil_names, f)
                depth_var = Variable(
                    name='soil_layer',
                    data=MemoryData(copy(f.variables['soil_layer'][:])),
                    dimensions=list(zip(f.variables['soil_layer'].dimensions,
                                   np.shape(f.variables['soil_layer']))),
                    attributes={
                        'standard_name': 'depth',
                        'long_name': 'Soil layer depth',
                        'positive': 'down',
                        'bounds': 'soil_bounds',
                        'units': 'm',
                        'axis': 'Z'
                    }
                )
                LOGGER.debug('soil_layer variable found in regcm file')
            except:
                soil_layer = np.array([0.00710063541719354, 0.0279250004153169,
                 0.062258573936546, 0.118865066900143, 0.212193395908963,
                 0.366065797104704, 0.619758497929827, 1.03802705000157,
                 1.7276353086672, 2.86460711317969])
                depth_var = Variable(
                    name='soil_layer',
                    data=MemoryData(soil_layer),
                    dimensions=(('soil_layer',10),),
                    attributes={
                        'standard_name': 'depth',
                        'long_name': 'Soil layer depth',
                        'positive': 'down',
                        'bounds': 'soil_bounds',
                        'units': 'm',
                        'axis': 'Z'
                    }
                )
            LOGGER.debug('Reading soil_bounds variable from regcm file')
            try:
                soil_bounds = get_var_with_name(bound_names, f)
                bounds_var = Variable(
                    name='soil_bounds',
                    data=MemoryData(copy(f.variables['soil_bounds'][:])),
                    dimensions=list(zip(f.variables['soil_bounds'].dimensions,
                                    np.shape(f.variables['soil_bounds']))),
                    attributes={
                        'standard_name': 'depth',
                        'long_name': 'Soil layer depth',
                        'units': 'm',
                    }
                )
                LOGGER.debug('soil_bounds variable found in regcm file')
            except:
                soil_bounds = np.array([[0, 0.0175128179162552],
                      [0.0175128179162552, 0.0450917871759315],
                      [0.0450917871759315, 0.0905618204183447],
                      [0.0905618204183447, 0.165529231404553],
                      [0.165529231404553, 0.289129596506834],
                      [0.289129596506834, 0.492912147517265],
                      [0.492912147517265, 0.828892773965698],
                      [0.828892773965698, 1.38283117933438],
                      [1.38283117933438, 2.29612121092344],
                      [2.29612121092344, 3.80188191232272]])
                bounds_var = Variable(
                    name='soil_bounds',
                    data=MemoryData(soil_bounds),
                    dimensions=(('soil_layer',10),('bnds', 2)),
                    attributes={
                        'standard_name': 'depth',
                        'long_name': 'Soil layer depth',
                        'units': 'm',
                    }
                )

        attributes = copy(prev_result.attributes)
        attributes['coordinates'] = 'soil_layer lat lon'

        depth_added_var = Variable(
            name=prev_result.name,
            data=prev_result.data,
            dimensions=prev_result.dimensions,
            attributes=attributes,
            times=prev_result.times,
            needs_time_bounds=prev_result.needs_time_bounds,
            times_attributes=prev_result.times_attributes,
            auxiliary_vars = (prev_result.auxiliary_variables +
                   [depth_var, bounds_var, ]),
        )

        return depth_added_var


class ExtractGroundHeight(Filter):
    """Extract the height from the ground, discarding empty data."""

    def __init__(self):
        pass

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):

        dimensions = prev_result.dimensions
        LOGGER.debug('Looking into the dimensions for a variable called mDDD')
        position = None
        height_val = None
        for i, dim in enumerate(dimensions):
            dim_name, dim_value = dim
            dim_match = re.match(r'^m(?P<height>\d+)', dim_name)
            if dim_match:
                position = i
                height_val = int(dim_match.groups('height')[0])
                LOGGER.debug(
                    'Found dimension %s in position %s (height value %s)',
                    dim_name,
                    position,
                    height_val,
                )
        if position is None:
            raise KeyError('No dimension found to extract height')

        LOGGER.debug('Removing dimension %s from the variable', position)
        reduced_dims = dimensions[0:position] + dimensions[position + 1:]

        attributes = copy(prev_result.attributes)
        if 'coordinates' in attributes:
            LOGGER.debug('Updating attribute "coordinates" adding "height"')
            attributes['coordinates'] = 'height ' + attributes['coordinates']

        data = SliceData(prev_result.data, position, 0)

        height_var = Variable(
            name='height',
            data=MemoryData(np.array([height_val], dtype=DATATYPE_AUXILIARIES)),
            dimensions=(),
            attributes={
                'standard_name': 'height',
                'long_name': 'height',
                'positive': 'up',
                'units': 'm',
                'axis': 'Z'
            }
        )

        remove_height_var = Variable(
            name=prev_result.name,
            data=data,
            dimensions=reduced_dims,
            attributes=attributes,
            times=prev_result.times,
            needs_time_bounds=prev_result.needs_time_bounds,
            times_attributes=prev_result.times_attributes,
            auxiliary_vars=prev_result.auxiliary_variables + [height_var, ],
        )

        return remove_height_var


class InterpolateHeight(Filter):
    """
    Apply a vertical interpolation with the height expressed as a pressure.

    Parameters:

        pressure_level (float): the height, expressed as pressure, to
            interpolate at.

        method (string in ``["linear", "logarithmic"]``): the method to use to
            interpolate.
    """

    def __init__(self, pressure_level, method, core):
        self.pressure_level = float(pressure_level)
        self.method = method.lower()
        if method not in ('linear', 'logarithmic'):
            raise ValueError(
                'Unknown method "{}" for vertical interpolation'.format(method)
            )
        self.core = core

        if method == 'linear':
            if core == 'hydrostatic':
                self.__interpolator = mod_vertint.intlin_hy
            elif core == 'moloch':
                self.__interpolator = mod_vertint.intlin_nonhy
            else:
                self.__interpolator = mod_vertint.intlin_nonhy
        else:
            if core == 'hydrostatic':
                self.__interpolator = mod_vertint.intlog_hy
            elif core == 'moloch':
                self.__interpolator = mod_vertint.intlog_nonhy
            else:
                self.__interpolator = mod_vertint.intlog_nonhy

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):

        dim_names = [i[0] for i in prev_result.dimensions]

        if 'kz' not in dim_names:
            raise ValueError(
                'No height dimension (kz) found in variable "{}"'
                .format(prev_result.name)
            )

        kz_index = dim_names.index('kz')
        LOGGER.debug('kz dimension is in position %s', kz_index)

        new_dimensions = [i for i in prev_result.dimensions if i[0] != 'kz']
        new_data_shape = [i[1] for i in prev_result.dimensions if i[0] != 'kz']
        LOGGER.debug('The output will be an array of shape %s', new_data_shape)

        LOGGER.debug('Allocating memory for the interpolation')
        interpolation_result = np.empty(new_data_shape, dtype=DATATYPE_MAIN)

        ps_names = get_regcm_variable_names('ps')
        ptop_names = get_regcm_variable_names('ptop')
        sigma_names = get_regcm_variable_names('sigma')

        p0_names = get_regcm_variable_names('p0')
        pp_names = get_regcm_variable_names('pp')
        pai_names = get_regcm_variable_names('pai')

        with Dataset(regcm_file_path, 'r') as f:
            if self.core == 'hydrostatic':
                LOGGER.debug('Reading ps variable from regcm file')
                ps = get_var_with_name(ps_names, f)
                LOGGER.debug(
                    'Reading ptop variable from regcm file and saving its content '
                    'in memory'
                )
                try:
                    ptop_array = np.array(
                        get_var_with_name(ptop_names, f)[:],
                        dtype=DATATYPE_AUXILIARIES,
                    )
                except:
                    ptop_array = np.array([100.0,])
            elif self.core == 'moloch':
                LOGGER.debug('Reading pai variable from regcm file')
                pai = get_var_with_name(pai_names, f)
            else:
                LOGGER.debug('Reading ps variable from regcm file')
                ps = get_var_with_name(ps_names, f)

                LOGGER.debug(
                    'Reading p0 variable from regcm file and saving its content '
                    'in memory'
                )
                p0 = np.array(
                    get_var_with_name(p0_names, f)[:],
                    dtype=DATATYPE_AUXILIARIES,
                )
                LOGGER.debug('Reading pp variable from regcm file')
                pp = get_var_with_name(pp_names, f)
                LOGGER.debug(
                    'Reading ptop variable from regcm file and saving its content '
                    'in memory'
                )
                try:
                    ptop_array = np.array(
                        get_var_with_name(ptop_names, f)[:],
                        dtype=DATATYPE_AUXILIARIES,
                    )
                except:
                    ptop_array = np.array([100.0,])
            LOGGER.debug(
                'Reading sigma variable from regcm file and saving its content '
                'in memory'
            )
            sigma_array = np.array(
                get_var_with_name(sigma_names, f)[:],
                dtype=DATATYPE_AUXILIARIES,
            )

            time_steps = prev_result.times.size

            with prev_result.data:
                for t in range(time_steps):

                    data_slice = [slice(None) for _ in prev_result.dimensions]
                    data_slice[0] = t

                    if self.core == 'hydrostatic':
                        LOGGER.debug('Reading time step %s of variable ps', t)
                        # PS, PTOP, PLEV must have same units. Use Pa
                        ps_t = np.array(ps[t, :], 
                                        dtype=DATATYPE_AUXILIARIES)
                        LOGGER.debug(
                            'Interpolating time step %s of %s',
                            t,
                            time_steps
                        )
                        interpolation_result[t, :] = self.__interpolator(
                            prev_result.data(data_slice),
                            ps_t,
                            sigma_array,
                            ptop_array,
                            self.pressure_level
                        )
                    elif self.core == 'moloch':
                        LOGGER.debug('Reading time step %s of variable pai', t)
                        pai_t = np.array(pai[t, :],
                                        dtype=DATATYPE_AUXILIARIES)
                        p_t = np.empty_like(prev_result.data(data_slice))
                        nk = np.shape(p_t)[0]
                        LOGGER.debug('Computing pressure')
                        p_t = np.power(pai_t,3.5)*100000.0
                        LOGGER.debug(
                            'Interpolating time step %s of %s',
                            t,
                            time_steps
                        )
                        interpolation_result[t, :] = self.__interpolator(
                            prev_result.data(data_slice),
                            p_t,
                            self.pressure_level
                        )
                    else:
                        # PTOP, PLEV must have same units. Use Pa
                        LOGGER.debug('Reading time step %s of variable pp', t)
                        # PS0, PP, PLEV must have same units. Use Pa
                        pp_t = np.array(pp[t, :], 
                                        dtype=DATATYPE_AUXILIARIES)
                        p_t = np.empty_like(prev_result.data(data_slice))
                        nk = np.shape(p_t)[0]
                        LOGGER.debug('Computing pressure')
                        ptpa = ptop_array * 100.0
                        for k in range (0,nk):
                            p_t[k,:] = ((p0[:]-ptpa)*sigma_array[k] + 
                                    ptpa) + pp_t[k,:]
                        LOGGER.debug(
                            'Interpolating time step %s of %s',
                            t,
                            time_steps
                        )
                        interpolation_result[t, :] = self.__interpolator(
                            prev_result.data(data_slice),
                            p_t,
                            self.pressure_level
                        )

            new_data = MemoryData(interpolation_result)

            attributes = copy(prev_result.attributes)
            if 'coordinates' in attributes:
                LOGGER.debug('Updating attribute "coordinates" adding "plev"')
                attributes['coordinates'] = 'plev ' + attributes['coordinates']

            pval = np.array([1], dtype=DATATYPE_AUXILIARIES)
            pval[0] = self.pressure_level
            height_var = Variable(
                name='plev',
                data=MemoryData(copy(pval)),
                dimensions=(),
                attributes={
                    'standard_name': 'air_pressure',
                    'long_name': 'pressure level',
                    'positive': 'down',
                    'units': 'Pa',
                    'axis': 'Z'
                }
            )

            return Variable(
                name=prev_result.name,
                data=new_data,
                dimensions=new_dimensions,
                attributes=attributes,
                times=prev_result.times,
                times_attributes=prev_result.times_attributes,
                auxiliary_vars=prev_result.auxiliary_variables + [height_var, ],
                needs_time_bounds=prev_result.needs_time_bounds,
            )


class InterpolateOnMultipleHeights(Filter):
    """
    Apply a vertical interpolation at many heights, expressed as pressures.

    Parameters:

        pressure_levels (iterable of floats): an iterable of heights, expressed
            as pressure, to interpolate at.

        method (string in ``["linear", "logarithmic"]``): the method to use to
            interpolate.

    Returns:

        list of ``Variable``: list of interpolated ``Variable`` objects.
    """

    def __init__(self, pressure_levels, method):
        self.pressure_levels = [float(i*100) for i in pressure_levels]
        self.method = method.lower()
        if method not in ('linear', 'logarithmic'):
            raise ValueError(
                'Unknown method "{}" for vertical interpolation'.format(method)
            )

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):

        interp_vars = []
        for pressure in self.pressure_levels:
            LOGGER.debug('Interpolating on pressure %s', pressure)

            LOGGER.debug(
                'Creating an InterpolateHeight filter for pressure %s',
                pressure
            )
            mcore = 'hydrostatic'
            try:
                core = Dataset(regcm_file_path).dynamical_core
                if core == 1:
                    mcore = 'hydrostatic'
                elif core == 2:
                    mcore = 'non-hydrostatic'
                elif core == 3:
                    mcore = 'moloch'
                else:
                    raise ValueError(
                            'Unknown dynamical_core : "{}"'.format(core))
            except:
                pass
            interp_filter = InterpolateHeight(pressure, self.method, mcore)

            LOGGER.debug('Preforming interpolation')
            interp_var = interp_filter(
                prev_result,
                regcm_file,
                regcm_file_path,
                simulation,
                corrflag,
                cordex_dir
            )

            interp_vars.append(interp_var)

        LOGGER.debug('Appending pressure level to the names of the variables')
        for interp_var, pressure in zip(interp_vars, self.pressure_levels):
            # The format .0f removes the ".0" at the end of the float
            interp_var.name += '{:.0f}'.format(pressure*0.01)

        return interp_vars


class ExtractLevel(Filter):
    """
    Given a level index and the name of a dimension, move the returned
    Variable object on that specific level.

    Parameters:

        level (integer): the index of the level to extract.

        dimension (string): the name of the dimension to extract.
    """

    def __init__(self, level, dimension):
        self.level = int(level)
        self.dimension = str(dimension)

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):

        dim_names = [i[0] for i in prev_result.dimensions]

        if self.dimension not in dim_names:
            raise ValueError('No dimension "{}" found in variable "{}"'
                             .format(self.dimension, prev_result.name))

        dim_index = dim_names.index(self.dimension)
        LOGGER.debug('%s is in position %s', self.dimension, dim_index)

        new_dims = [i for i in prev_result.dimensions if i[0] != self.dimension]
        LOGGER.debug(
            'Changing dimensions of variable %s from %s to %s',
            prev_result.name,
            prev_result.dimensions,
            new_dims,
        )

        LOGGER.debug('Preparing new data')
        new_data = SliceData(
            prev_result.data,
            dim_index,
            self.level
        )

        LOGGER.debug('Returning new variable with extracted level')
        return Variable(
            name=prev_result.name,
            data=new_data,
            dimensions=new_dims,
            attributes=prev_result.attributes,
            times=prev_result.times,
            times_attributes=prev_result.times_attributes,
            auxiliary_vars=prev_result.auxiliary_variables,
            needs_time_bounds=prev_result.needs_time_bounds,
        )


class NormFromCoords(ActionStarter):
    """
    Compute norms of a list of variables and merge them into a single one.

    Parameters:

        var_name (string): the name of the variable to create.

        compute_from (list of strings): the list of variables to calculate the
            norms on, and to collect as a single variable to be returned.

        copy_attributes_from (string, optional): if provided, copy the variable
            attributes from that variable, otherwise use the previous variable
            attributes.

        new_attributes (dict ``string -> object``, optional): if provided,
            contains attributes that can override those from the previous
            Method object. netCDF takes care of applying the object value
            correctly.

        need_time_bounds (boolean, optional): whether to use an auxiliary
            variable to store the interval range. Defaults to ``False``.
    """

    def __init__(self, var_name, compute_from, copy_attributes_from=None,
                 new_attributes=None, need_time_bounds=False):
        self.var_name = var_name
        self.compute_from = compute_from

        if (copy_attributes_from is not None and
                copy_attributes_from not in compute_from):
            raise ValueError(
                '{} is not in the list of variables to be read'
                .format(copy_attributes_from)
            )
        self.attributes_from = copy_attributes_from

        self.new_attributes = new_attributes
        self.need_time_bounds = need_time_bounds

    def __call__(self, regcm_file, regcm_file_path, simulation, corrflag,
                 cordex_dir):

        if self.attributes_from is None:
            k = 0
        else:
            k = self.compute_from.index(self.attributes_from)

        LOGGER.debug('Reading variable %s as reference', self.compute_from[k])
        reference_var = read_variable_from_regcmfile(
            self.compute_from[k],
            regcm_file,
            regcm_file_path,
            self.need_time_bounds
        )

        if self.attributes_from is None:
            attributes = {}
        else:
            LOGGER.debug('Copying attributes from %s', reference_var.name)
            attributes = reference_var.attributes

        if self.new_attributes is not None:
            for attr_name, attr_value in self.new_attributes.items():
                LOGGER.debug(
                    'Adding attribute "%s" with value "%s"',
                    attr_name,
                    attr_value
                )
                attributes[attr_name] = attr_value

        LOGGER.debug('Preparing data for the new variable')

        # def norm_f(*args):
        #     return np.sqrt(sum([x * x for x in args]))

        var_data = NetCDFMultiData(
            self.compute_from,
            np.linalg.norm,
            regcm_file_path,
            DATATYPE_MAIN,
            lock=regcm_file.lock,
        )

        return Variable(
            name=self.var_name,
            data=var_data,
            dimensions=reference_var.dimensions,
            attributes=attributes,
            times=reference_var.times,
            times_attributes=reference_var.times_attributes,
            auxiliary_vars=[],
            needs_time_bounds=self.need_time_bounds,
        )


class Thin(Filter):
    """
    Change the time-step of a variable.

    Note that if the timestep is the same as the actual one (compared with an
    epsilon to avoid floating-point comparison), the same input ``Variable`` is
    returned.

    Parameters:

        new_time_step (float): the new time-step to apply to the variable.
    """

    def __init__(self, new_time_step):
        self.new_time_step = float(new_time_step)

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):
        if not prev_result.depends_on_time:
            raise ValueError('Variable {} does not depend on time'.format(
                prev_result.name
            ))

        current_step = prev_result.time_step_size
        if abs(current_step - self.new_time_step) < 0.01:
            LOGGER.debug(
                'No thinning needed (the timestep is already %s hours)',
                self.new_time_step
            )
            return prev_result

        LOGGER.debug(
            'Trying to extract from a variable that has a value every %s hours'
            ' a new variable with a value every %s',
            current_step,
            self.new_time_step
        )

        if current_step > self.new_time_step:
            LOGGER.debug(
                'Cannot operate: (the timestep is %s hours)',
                current_step
            )
            return prev_result

        thin_step = self.new_time_step / current_step
        if abs(thin_step - int(thin_step)) > 0.1:
            raise ValueError(
                'The old time step ({}h) and the new one ({}h) are not '
                'multiples'.format(current_step, self.new_time_step)
            )
        thin_step = int(thin_step)
        LOGGER.debug('Taking a value every %s', thin_step)

        LOGGER.debug('Computing the starting point for the thinning')
        starting_point = 0
        for i, t in enumerate(prev_result.times[:thin_step]):
            if int(t) % thin_step == 0:
                LOGGER.debug('Starting from %s', i)
                starting_point = i
                break
        else:
            LOGGER.warning('Unable to find a correct starting point for '
                           'thinning. Starting from 0')

        dims_names = [d[0] for d in prev_result.dimensions]
        time_index = dims_names.index('time')
        new_data_slice = [slice(None) for _ in prev_result.dimensions]
        new_data_slice[time_index] = slice(starting_point, None, thin_step)

        var_data = SubselectedData(prev_result.data, new_data_slice)
        new_times = prev_result.times[starting_point::thin_step]

        new_dimensions = copy(prev_result.dimensions)
        new_dimensions[time_index] = ('time', new_times.size)
        LOGGER.debug(
            'Changing variable dimension from %s to %s',
            prev_result.dimensions,
            new_dimensions
        )

        return Variable(
            name=prev_result.name,
            data=var_data,
            dimensions=new_dimensions,
            attributes=prev_result.attributes,
            times=new_times,
            times_attributes=prev_result.times_attributes,
            auxiliary_vars=prev_result.auxiliary_variables,
            needs_time_bounds=prev_result.needs_time_bounds,
        )

class ComputeMaximum(Filter):
    """
    Compute the maximum of the variable using a given time-step as interval.

    Note that if the timestep is the same as the actual one (compared with an
    epsilon to avoid floating-point comparison), the same input ``Variable`` is
    returned.

    Parameters:

        new_time_step (float): the new time-step to apply to the variable.
    """

    def __init__(self, new_time_step):
        self.new_time_step = float(new_time_step)

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):
        if not prev_result.depends_on_time:
            raise ValueError('Variable {} does not depend on time'.format(
                prev_result.name
            ))

        current_step = prev_result.time_step_size
        if abs(current_step - self.new_time_step) < 0.01:
            LOGGER.debug(
                'No maximum needed (the timestep is already %s hours)',
                self.new_time_step
            )
            return prev_result

        LOGGER.debug(
            'Trying to create maximums of length %s from a variable that has a '
            'value every %s hours',
            self.new_time_step,
            current_step,
        )

        if current_step > self.new_time_step:
            LOGGER.debug(
                'Cannot operate: (the timestep is %s hours)',
                current_step
            )
            return prev_result

        maximum_step = self.new_time_step / current_step
        if abs(maximum_step - int(maximum_step)) > 0.1:
            raise ValueError(
                'The old time step ({}h) and the new one ({}h) are not '
                'multiples'.format(current_step, self.new_time_step)
            )
        maximum_step = int(maximum_step)
        LOGGER.debug(
            'A new value will be the maximum of %s old values',
            maximum_step
        )

        LOGGER.debug('Computing the starting point for the maximum')
        # This is a tricky point. Let say that we have timesteps every 6 hours
        # and that we want to perform a daily means. Ideally, we would like to
        # have times at the hours 3, 9, 15 and 21. In this case, we simply have
        # to perform the mean among all four timesteps. If the simulation
        # starts between two days, otherwise, we could have something like
        # 15, 21, 3, 9, 15, 21, ...
        # In that case, we have to discard the first two values. Therefore, the
        # variable starting_point will assume the value 2.
        # Finally, there is another problem to address.
        # If the values do not refer to some temporal means but just to some
        # temporal instant (and, therefore, they do not have temporal bounds),
        # the values could be like 6, 12, 18, 0, ...
        # In this case, we will subtract 3 to all the time steps. This allows
        # to have a correct mean considering 00.00 as the last point of one day
        # (and not as the first point of the following day)
        starting_point = 0

        if prev_result.needs_time_bounds:
            corr_factor = 0
        else:
            # NEED TO INVESTIGATE. WHY THIS?
            corr_factor = -prev_result.time_step_size / 2.
        LOGGER.debug(
            'Using a correction factor for times of value %s',
            corr_factor
        )

        for i in range(maximum_step):
            LOGGER.debug('Trying starting from %s', i)
            j = i + maximum_step
            times_to_maximum = prev_result.times[i:j] + corr_factor
            time_mean = np.mean(times_to_maximum)

            start_interval = time_mean - (self.new_time_step / 2)
            if start_interval % maximum_step == 0:
                LOGGER.debug('Starting from %s', i)
                starting_point = i
                break
            else:
                LOGGER.debug('Trying again!')
        else:
            LOGGER.warning('Unable to find a correct starting point for '
                           'maximums. Starting from 0')

        maximums_num = (prev_result.times.size - starting_point) // maximum_step
        LOGGER.debug(
            'There are %s time step with a frequency of %s. The means start '
            'the timestep number %s. Therefore, %s means will be computed',
            prev_result.times.size,
            prev_result.frequency,
            starting_point,
            maximums_num
        )

        LOGGER.debug('Looking for a dimension called "time"')
        dims_names = [d[0] for d in prev_result.dimensions]
        time_index = dims_names.index('time')
        LOGGER.debug('Time is the dimension number %s', time_index)

        LOGGER.debug('Computing the new dimensions of the variables')
        new_dimensions = copy(prev_result.dimensions)
        new_dimensions[time_index] = ('time', maximums_num)
        LOGGER.debug('The new dimensions will be %s', new_dimensions)

        LOGGER.debug('Allocating space for the new time vector')
        new_times = np.empty(maximums_num, dtype=prev_result.times.dtype)

        LOGGER.debug('Allocating space for the new data of the variable')
        new_data_array = np.empty(
            [d[1] for d in new_dimensions],
            dtype=DATATYPE_MAIN,
        )

        LOGGER.debug('Preparing a mask to hide NaN values')
        new_data_mask = np.zeros(
            [d[1] for d in new_dimensions],
            dtype=bool
        )

        shape_template = [slice(None) for _ in prev_result.dimensions]
        old_data_slice = shape_template[:]
        new_data_slice = shape_template[:]
        LOGGER.debug('Computing maximums')
        with prev_result.data:
            for i in range(maximums_num):
                start = starting_point + maximum_step * i
                end = start + maximum_step
                LOGGER.debug(
                    'Performing maximum of timesteps from %s to %s',
                    start,
                    end
                )
                new_times[i] = np.mean(
                    prev_result.times[start:end] + corr_factor
                )

                old_data_slice[time_index] = slice(start, end)
                new_data_slice[time_index] = slice(i,i)

                prev_data = prev_result.data(old_data_slice)
                new_data_array[i] = np.nanmax(
                    prev_data,
		    axis=time_index,
		    keepdims=True,
		    )

                if is_masked(prev_data):
                    new_data_mask[tuple(new_data_slice)] = np.prod(
                        prev_data.mask,
                        axis=time_index
                    )

        if np.any(new_data_mask):
            LOGGER.debug(
                'Converting array in a masked array to remove invalid '
                'values'
            )
            new_data_array = np.ma.masked_array(
                new_data_array,
                new_data_mask
            )

        var_data = MemoryData(new_data_array)

        LOGGER.debug('Going to save the variable with name %s', prev_result.name)

        attributes = copy(prev_result.attributes)
        attributes['cell_methods'] = 'time: maximum'

        return Variable(
            name=prev_result.name,
            data=var_data,
            dimensions=new_dimensions,
            attributes=attributes,
            times=new_times,
            times_attributes=prev_result.times_attributes,
            auxiliary_vars=prev_result.auxiliary_variables,
            needs_time_bounds=True,
        )


class IfNeededMaximumAndSave(Filter):
    """
    Compute the maximum  of the variable in a new interval, and save it to
    disk if needed.

    Parameters:

        new_time_step (float): the new time-step to apply to the variable.

        var_name (string): the name of the variable to extract from the file.

        fill_value (float, optional): if provided, fill the variable with this
            value.

        new_attributes (dict ``string -> object``, optional): if provided,
            contains attributes that can override those from the previous
            Method object. netCDF takes care of applying the object value
            correctly.
    """

    def __init__(self, new_time_step, var_name=None, fill_value=None,
                 new_attributes=None):
        self.__maximum = ComputeMaximum(new_time_step)
        self.__savedisk = SaveVariableToDisk(var_name, fill_value,
                                             new_attributes)

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):
        LOGGER.debug('Performing maximum')
        averaged_var = self.__maximum(
            prev_result,
            regcm_file,
            regcm_file_path,
            simulation,
            corrflag,
            cordex_dir
        )

        # Since Variable doesn't provide an __eq__ method, Python defaults to
        # use the identity equality. Moreover, the average is a function that
        # may return the same identical object in case no average operation
        # has to be performed (if the time-step is the same as the one
        # requested).
        if averaged_var is not prev_result:
            LOGGER.debug('Saving on disk')
            self.__savedisk(
                averaged_var,
                regcm_file,
                regcm_file_path,
                simulation,
                corrflag,
                cordex_dir
            )

        return prev_result

class IfNeededThinAndSave(Filter):
    """
    Change the time-step of a variable and save it to disk if needed.

    Parameters:

        new_time_step (float): the new time-step to apply to the variable.

        var_name (string): the name of the variable to extract from the file.

        fill_value (float, optional): if provided, fill the variable with this
            value.

        new_attributes (dict ``string -> object``, optional): if provided,
            contains attributes that can override those from the previous
            Method object. netCDF takes care of applying the object value
            correctly.
    """

    def __init__(self, new_time_step, var_name=None, fill_value=None,
                 new_attributes=None):
        self.__thin = Thin(new_time_step)
        self.__savedisk = SaveVariableToDisk(var_name, fill_value,
                                             new_attributes)

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):
        LOGGER.debug('Performing thinning')
        thinned_var = self.__thin(
            prev_result,
            regcm_file,
            regcm_file_path,
            simulation,
            corrflag,
            cordex_dir
        )

        # Since Variable doesn't provide an __eq__ method, Python defaults to
        # use the identity equality. Moreover, the thinning is a function that
        # may return the same identical object in case no thinning operation
        # has to be performed (if the time-step is the same as the one
        # requested).
        if thinned_var is not prev_result:
            LOGGER.debug('Saving on disk')
            self.__savedisk(
                thinned_var,
                regcm_file,
                regcm_file_path,
                simulation,
                corrflag,
                cordex_dir
            )

        return prev_result

class IfNeededThinAndSaveForEach(Filter):
    """
    Change the time-step of multiple variable and save them to disk if needed.

    Parameters:

        new_time_step (float): the new time-step to apply to the variable.

        var_name (string): the name of the variable to extract from the file.

        fill_value (float, optional): if provided, fill the variable with this
            value.

        new_attributes (dict ``string -> object``, optional): if provided,
            contains attributes that can override those from the previous
            Method object. netCDF takes care of applying the object value
            correctly.
    """

    def __init__(self, new_time_step, var_name=None, fill_value=None,
                 new_attributes=None):
        self.__thin = Thin(new_time_step)
        self.__savedisk = SaveVariableToDisk(var_name, fill_value,
                                             new_attributes)

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):
        thin_variables = []
        for cordex_var in prev_result:
            LOGGER.debug('Performing thinning of variable %s', cordex_var.name)
            new_thin = self.__thin(
                cordex_var,
                regcm_file,
                regcm_file_path,
                simulation,
                corrflag,
                cordex_dir
            )
            # Since Variable doesn't provide an __eq__ method, Python defaults to
            # use the identity equality. Moreover, the thinning is a function that
            # may return the same identical object in case no thinning operation
            # has to be performed (if the time-step is the same as the one
            # requested).
            if new_thin is not prev_result:
                LOGGER.debug('Saving on disk')
                self.__savedisk(
                    new_thin,
                    regcm_file,
                    regcm_file_path,
                    simulation,
                    corrflag,
                    cordex_dir
                )

                thin_variables.append(new_thin)

        return thin_variables


class ComputeAverage(Filter):
    """
    Compute the average of the variable using a given time-step as interval.

    Note that if the timestep is the same as the actual one (compared with an
    epsilon to avoid floating-point comparison), the same input ``Variable`` is
    returned.

    Parameters:

        new_time_step (float): the new time-step to apply to the variable.
    """

    def __init__(self, new_time_step):
        self.new_time_step = float(new_time_step)

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):
        if not prev_result.depends_on_time:
            raise ValueError('Variable {} does not depend on time'.format(
                prev_result.name
            ))

        current_step = prev_result.time_step_size
        if abs(current_step - self.new_time_step) < 0.01:
            LOGGER.debug(
                'No average needed (the timestep is already %s hours)',
                self.new_time_step
            )
            return prev_result

        LOGGER.debug(
            'Trying to create averages of length %s from a variable that has a '
            'value every %s hours',
            self.new_time_step,
            current_step,
        )

        if current_step > self.new_time_step:
            LOGGER.debug(
                'Cannot operate: (the timestep is %s hours)',
                current_step
            )
            return prev_result

        average_step = self.new_time_step / current_step
        if abs(average_step - int(average_step)) > 0.1:
            raise ValueError(
                'The old time step ({}h) and the new one ({}h) are not '
                'multiples'.format(current_step, self.new_time_step)
            )
        average_step = int(average_step)
        LOGGER.debug(
            'A new value will be the average of %s old values',
            average_step
        )

        LOGGER.debug('Computing the starting point for the average')
        # This is a tricky point. Let say that we have timesteps every 6 hours
        # and that we want to perform a daily means. Ideally, we would like to
        # have times at the hours 3, 9, 15 and 21. In this case, we simply have
        # to perform the mean among all four timesteps. If the simulation
        # starts between two days, otherwise, we could have something like
        # 15, 21, 3, 9, 15, 21, ...
        # In that case, we have to discard the first two values. Therefore, the
        # variable starting_point will assume the value 2.
        # Finally, there is another problem to address.
        # If the values do not refer to some temporal means but just to some
        # temporal instant (and, therefore, they do not have temporal bounds),
        # the values could be like 6, 12, 18, 0, ...
        # In this case, we will subtract 3 to all the time steps. This allows
        # to have a correct mean considering 00.00 as the last point of one day
        # (and not as the first point of the following day)
        starting_point = 0

        if prev_result.needs_time_bounds:
            corr_factor = 0
        else:
            # NEED TO INVESTIGATE. WHY THIS?
            corr_factor = -prev_result.time_step_size / 2.
        LOGGER.debug(
            'Using a correction factor for times of value %s',
            corr_factor
        )

        for i in range(average_step):
            LOGGER.debug('Trying starting from %s', i)
            j = i + average_step
            times_to_average = prev_result.times[i:j] + corr_factor
            time_mean = np.mean(times_to_average)

            start_interval = time_mean - (self.new_time_step / 2)
            if start_interval % average_step == 0:
                LOGGER.debug('Starting from %s', i)
                starting_point = i
                break
            else:
                LOGGER.debug('Trying again!')
        else:
            LOGGER.warning('Unable to find a correct starting point for '
                           'averages. Starting from 0')

        averages_num = (prev_result.times.size - starting_point) // average_step
        LOGGER.debug(
            'There are %s time step with a frequency of %s. The means start '
            'the timestep number %s. Therefore, %s means will be computed',
            prev_result.times.size,
            prev_result.frequency,
            starting_point,
            averages_num
        )

        LOGGER.debug('Looking for a dimension called "time"')
        dims_names = [d[0] for d in prev_result.dimensions]
        time_index = dims_names.index('time')
        LOGGER.debug('Time is the dimension number %s', time_index)

        LOGGER.debug('Computing the new dimensions of the variables')
        new_dimensions = copy(prev_result.dimensions)
        new_dimensions[time_index] = ('time', averages_num)
        LOGGER.debug('The new dimensions will be %s', new_dimensions)

        LOGGER.debug('Allocating space for the new time vector')
        new_times = np.empty(averages_num, dtype=prev_result.times.dtype)

        LOGGER.debug('Allocating space for the new data of the variable')
        new_data_array = np.empty(
            [d[1] for d in new_dimensions],
            dtype=DATATYPE_MAIN,
        )

        LOGGER.debug('Preparing a mask to hide NaN values')
        new_data_mask = np.zeros(
            [d[1] for d in new_dimensions],
            dtype=bool
        )

        shape_template = [ slice(None) for _ in prev_result.dimensions]
        old_data_slice = shape_template[:]
        new_data_slice = shape_template[:]
        LOGGER.debug('Computing means')
        with prev_result.data:
            for i in range(averages_num):
                start = starting_point + average_step * i
                end = start + average_step
                LOGGER.debug(
                    'Performing mean of timesteps from %s to %s',
                    start,
                    end
                )
                new_times[i] = np.mean(
                    prev_result.times[start:end] + corr_factor
                )

                old_data_slice[time_index] = slice(start, end)
                new_data_slice[time_index] = i

                prev_data = prev_result.data(old_data_slice)
                new_data_array[tuple(new_data_slice)] =  np.ma.mean(
                    prev_data,
                    axis=time_index
                )

                if is_masked(prev_data):
                    new_data_mask[tuple(new_data_slice)] = np.prod(
                        prev_data.mask,
                        axis=time_index
                    )

        if np.any(new_data_mask):
            LOGGER.debug(
                'Converting array in a masked array to remove invalid '
                'values'
            )
            new_data_array = np.ma.masked_array(
                new_data_array,
                new_data_mask
            )

        var_data = MemoryData(new_data_array)

        LOGGER.debug('Going to save the variable with name %s', prev_result.name)

        attributes = copy(prev_result.attributes)
        attributes['cell_methods'] = 'time: mean'

        return Variable(
            name=prev_result.name,
            data=var_data,
            dimensions=new_dimensions,
            attributes=attributes,
            times=new_times,
            times_attributes=prev_result.times_attributes,
            auxiliary_vars=prev_result.auxiliary_variables,
            needs_time_bounds=True,
        )


class IfNeededAverageAndSave(Filter):
    """
    Compute the average of the variable in a new interval, and save it to
    disk if needed.

    Parameters:

        new_time_step (float): the new time-step to apply to the variable.

        var_name (string): the name of the variable to extract from the file.

        fill_value (float, optional): if provided, fill the variable with this
            value.

        new_attributes (dict ``string -> object``, optional): if provided,
            contains attributes that can override those from the previous
            Method object. netCDF takes care of applying the object value
            correctly.
    """

    def __init__(self, new_time_step, var_name=None, fill_value=None,
                 new_attributes=None):
        self.__averg = ComputeAverage(new_time_step)
        self.__savedisk = SaveVariableToDisk(var_name, fill_value,
                                             new_attributes)

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):
        LOGGER.debug('Performing average')
        averaged_var = self.__averg(
            prev_result,
            regcm_file,
            regcm_file_path,
            simulation,
            corrflag,
            cordex_dir
        )

        # Since Variable doesn't provide an __eq__ method, Python defaults to
        # use the identity equality. Moreover, the average is a function that
        # may return the same identical object in case no average operation
        # has to be performed (if the time-step is the same as the one
        # requested).
        if averaged_var is not prev_result:
            LOGGER.debug('Saving on disk')
            self.__savedisk(
                averaged_var,
                regcm_file,
                regcm_file_path,
                simulation,
                corrflag,
                cordex_dir
            )

        return prev_result


class LinearCombinationOfVars(ActionStarter):
    """
    Combine the Variables using a linear combination.

    Parameters:

        var_name (string): the name of the variable to extract from the file.

        combination (iterable of iterables ``(string, number)``): each
            pair has to be an iterable of ``(variable_name, coefficient)``.

        copy_attributes_from (string, optional): if provided, copy the variable
            attributes from that variable, otherwise use the previous variable
            attributes.

        new_attributes (dict ``string -> object``, optional): if provided,
            contains attributes that can override those from the previous
            Method object. netCDF takes care of applying the object value
            correctly.

        need_time_bounds (boolean, optional): whether to use an auxiliary
            variable to store the interval range. Defaults to ``False``.
    """

    def __init__(self, var_name, combination, copy_attributes_from=None,
                 new_attributes=None, need_time_bounds=False):
        self.var_name = var_name
        self.var_list = combination

        if copy_attributes_from is not None:
            if copy_attributes_from not in [v[0] for v in combination]:
                raise ValueError(
                    '{} is not in the list of variables to be read'
                    .format(copy_attributes_from)
                )
        self.attributes_from = copy_attributes_from

        self.new_attributes = new_attributes
        self.need_time_bounds = need_time_bounds

    def __call__(self, regcm_file, regcm_file_path, simulation, corrflag,
                 cordex_dir):

        var_names = [v[0] for v in self.var_list]
        var_coefficients = [v[1] for v in self.var_list]

        if self.attributes_from is None:
            k = 0
        else:
            k = var_names.index(self.attributes_from)

        LOGGER.debug('Reading variable %s as reference', var_names[k])
        reference_var = read_variable_from_regcmfile(
            var_names[k],
            regcm_file,
            regcm_file_path,
            self.need_time_bounds
        )

        if self.attributes_from is None:
            attributes = {}
        else:
            LOGGER.debug('Copying attributes from %s', reference_var.name)
            attributes = reference_var.attributes

        if self.new_attributes is not None:
            for attr_name, attr_value in self.new_attributes.items():
                LOGGER.debug(
                    'Adding attribute "%s" with value "%s"',
                    attr_name,
                    attr_value
                )
                attributes[attr_name] = attr_value

        LOGGER.debug('Preparing data for the new variable')

        def linear_comb(*args):
            return sum([c * x for c, x in zip(var_coefficients, args)])

        var_data = NetCDFMultiData(
            var_names,
            linear_comb,
            regcm_file_path,
            DATATYPE_MAIN,
            lock=regcm_file.lock
        )

        return Variable(
            name=self.var_name,
            data=var_data,
            dimensions=reference_var.dimensions,
            attributes=attributes,
            times=reference_var.times,
            times_attributes=reference_var.times_attributes,
            auxiliary_vars=[],
            needs_time_bounds=self.need_time_bounds,
        )


class SumOnDimension(Filter):
    """Sum variables on a given dimension.

    Parameters:

        dimension (string): the name of the dimension to extract.
    """

    def __init__(self, dimension):
        self.dimension = dimension

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):
        dimensions = prev_result.dimensions
        LOGGER.debug(
            'Looking into the dimensions for one called %s',
            self.dimension
        )

        position = None
        for i, dim in enumerate(dimensions):
            dim_name, dim_value = dim
            if dim_name == self.dimension:
                position = i
                LOGGER.debug(
                    'Found dimension %s in position %s',
                    dim_name,
                    position,
                )
        if position is None:
            raise KeyError('No dimension found named {}'.format(self.dimension))

        LOGGER.debug('Removing dimension %s from the variable', position)
        reduced_dims = dimensions[0:position] + dimensions[position + 1:]
        LOGGER.debug('New dimensions are %s', reduced_dims)

        LOGGER.debug('Copying attributes from the previous variable')
        attributes = copy(prev_result.attributes)

        LOGGER.debug('Allocating space in memory to save the sum')
        new_data_array = np.empty(
            [d[1] for d in reduced_dims],
            dtype=DATATYPE_MAIN
        )

        LOGGER.debug('Preparing a mask to hide NaN values')
        new_data_mask = np.zeros(
            [d[1] for d in reduced_dims],
            dtype=bool
        )

        if prev_result.depends_on_time:
            LOGGER.debug('Checking position of time dimension')
            time_position = [d[0] for d in prev_result.dimensions].index('time')
            LOGGER.debug('Time is the dimension number %s', time_position)

            slice_template_new = [slice(None) for _ in reduced_dims]
            slice_template_old = [slice(None) for _ in prev_result.dimensions]

            if time_position < position:
                position_no_time = position - 1
            else:
                position_no_time = position
            LOGGER.debug(
                'When cycling on the time dimension, the position of variable '
                '%s will be %s',
                self.dimension,
                position_no_time
            )

            LOGGER.debug('Reading data from the previous variable')
            with prev_result.data:
                for i in range(prev_result.times.size):
                    LOGGER.debug(
                        'Summing timestep %s of %s',
                        i + 1,
                        prev_result.times.size
                    )

                    current_slice_new = slice_template_new[:]
                    current_slice_new[time_position] = i

                    current_slice_old = slice_template_old[:]
                    current_slice_old[time_position] = i

                    prev_data = prev_result.data(current_slice_old)
                    new_data_array[tuple(current_slice_new)] = np.ma.sum(
                        prev_data,
                        axis=position_no_time
                    )

                    if is_masked(prev_data):
                        new_data_mask[tuple(current_slice_new)] = np.prod(
                            prev_data.mask,
                            axis=position_no_time
                        )

        else:
            LOGGER.debug('Reading data from the previous variable')
            with prev_result.data:
                prev_data = prev_result.data()
                new_data_array[:] = np.ma.sum(prev_data, axis=position)
                if is_masked(prev_data):
                    new_data_mask[:] = np.prod(prev_data.mask, axis=position)

        if np.any(new_data_mask):
            LOGGER.debug(
                'Converting array in a masked array to remove invalid '
                'values'
            )
            new_data_array = np.ma.masked_array(
                new_data_array,
                new_data_mask
            )

        new_var = Variable(
            name=prev_result.name,
            data=MemoryData(new_data_array),
            dimensions=reduced_dims,
            attributes=attributes,
            times=prev_result.times,
            needs_time_bounds=prev_result.needs_time_bounds,
            times_attributes=prev_result.times_attributes,
            auxiliary_vars=prev_result.auxiliary_variables,
        )

        return new_var


class ComputeAverageOfEachVariable(Filter):
    """Compute the average of each variable, given a certain time-step.

    Note that if the timestep is the same as the actual one (compared with an
    epsilon to avoid floating-point comparison), the same input ``Variable`` is
    returned.

    Parameters:

        new_time_step (float): the new time-step to apply to the variable.

    Returns:

        list of ``Variable``: list of ``Variable`` objects.
    """

    def __init__(self, new_time_step):
        self.__compute_average = ComputeAverage(new_time_step=new_time_step)

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):

        mean_variables = []
        for cordex_var in prev_result:
            LOGGER.debug('Computing average of variable %s', cordex_var.name)
            new_mean = self.__compute_average(
                cordex_var,
                regcm_file,
                regcm_file_path,
                simulation,
                corrflag,
                cordex_dir
            )

            mean_variables.append(new_mean)

        return mean_variables


class ComputeMaximumOfEachVariable(Filter):
    """Compute the maximum of each variable, given a certain time-step.

    Note that if the timestep is the same as the actual one (compared with an
    epsilon to avoid floating-point comparison), the same input ``Variable`` is
    returned.

    Parameters:

        new_time_step (float): the new time-step to apply to the variable.

    Returns:

        list of ``Variable``: list of ``Variable`` objects.
    """

    def __init__(self, new_time_step):
        self.__compute_maximum = ComputeMaximum(new_time_step=new_time_step)

    def __call__(self, prev_result, regcm_file, regcm_file_path, simulation,
                 corrflag, cordex_dir):

        max_variables = []
        for cordex_var in prev_result:
            LOGGER.debug('Computing maximum of variable %s', cordex_var.name)
            new_max = self.__compute_maximum(
                cordex_var,
                regcm_file,
                regcm_file_path,
                simulation,
                corrflag,
                cordex_dir
            )

            max_variables.append(new_max)

        return max_variables


class ComputeGeopotentialHeight(ActionStarter):
    """
    Compute the geopotential at a certain given height expressed as pressure.

    Parameters:

        height (float): the height at which to compute the geopotential,
            expressed as pressure.

    Returns:

        ``Variable``: a ``Variable`` object, containing the height expressed as
        meters.
    """

    def __init__(self, height):
        self.height = float(height * 100.0) #Pa

    def __call__(self, regcm_file, regcm_file_path, simulation, corrflag,
                 cordex_dir):

        mcore = 'hydrostatic'
        try:
            core = Dataset(regcm_file_path).dynamical_core
            if core == 1:
                mcore = 'hydrostatic'
            elif core == 2:
                mcore = 'non-hydrostatic'
            elif core == 3:
                mcore = 'moloch'
            else:
                raise ValueError(
                        'Unknown dynamical_core : "{}"'.format(core))
        except:
            pass

        # Prepare a function to read exactly one timestep of a netcdf var
        def read_timestep(var_pointer, t_step):
            if 'time' not in var_pointer.dimensions:
                return no_time_vars_data[var_pointer.name]

            t_dim = var_pointer.dimensions.index('time')
            slicer = [slice(None) for _ in var_pointer.shape]
            slicer[t_dim] = t_step
            LOGGER.debug(
                'Reading timestep %s of variable %s',
                t_step,
                var_pointer.name
            )
            return np.array(var_pointer[slicer], dtype=DATATYPE_MAIN)

        if mcore == 'hydrostatic':
            needed_var_names = ['ta', 'ps', 'topo', 'sigma', 'ptop']
            funcname = mod_hgt.height_hydro
        elif mcore == 'moloch':
            needed_var_names = ['ta', 'pai', 'ps', 'topo', 'a', 'b']
            funcname = mod_hgt.height_moloch
        else:
            needed_var_names = ['ta', 'p0', 'ps', 'topo', 'sigma', 'ptop', 'pp']
            funcname = mod_hgt.height_nonhydro

        with Dataset(regcm_file_path, 'r') as f:
            needed_vars = []
            for v in needed_var_names:
                v_names = get_regcm_variable_names(v)
                if v_names == ['ptop']:
                    try:
                        v_pointer = get_var_with_name(v_names, f)
                    except:
                        v_pointer = fake_var((100.0,))
                else:
                    v_pointer = get_var_with_name(v_names, f)
                needed_vars.append(v_pointer)

            # This is the name of the variable that will be used to copy the
            # time steps, the time attributes and so on...
            for v_name, v in zip(needed_var_names, needed_vars):
                if 'time' in v.dimensions:
                    reference_var_name = v_name
                    LOGGER.debug(
                        'Variable %s will be used as reference',
                        reference_var_name
                    )
                    break
            else:
                raise ValueError('No variable found to use as a reference')

            LOGGER.debug('Check the total number of steps')
            time_dim = needed_vars[0].dimensions.index('time')
            time_steps = needed_vars[0].shape[time_dim]
            LOGGER.debug('%s timesteps found', time_steps)

            kz_dim = needed_vars[0].dimensions.index('kz')
            old_shape = needed_vars[0].shape
            new_shape = old_shape[:kz_dim] + old_shape[kz_dim + 1:]
            old_dims = needed_vars[0].dimensions
            new_dims = old_dims[:kz_dim] + old_dims[kz_dim + 1:]
            dimensions = list(zip(new_dims, new_shape))
            LOGGER.debug(
                'The output data will have the following shape: %s',
                new_shape
            )
            LOGGER.debug(
                'The output data will have the following dimensions: %s',
                dimensions
            )
            LOGGER.debug('Allocating memory')
            data = np.empty(new_shape, dtype=DATATYPE_MAIN)

            # Read the data of the variables that do not depend on time
            no_time_vars_data = {}
            for v in needed_vars:
                if 'time' not in v.dimensions:
                    LOGGER.debug('Reading data of variable %s', v.name)
                    if v.name == 'ptop':
                        v_data = np.array(v[:], dtype=DATATYPE_AUXILIARIES)
                    else:
                        v_data = np.array(v[:], dtype=DATATYPE_MAIN)
                    no_time_vars_data[v.name] = v_data

            for t in range(time_steps):
                LOGGER.debug(
                    'Computing geopotential height for %s Pa at timestep %s',
                    self.height,
                    t
                )
                hgt_args = [read_timestep(data_v, t) for data_v in needed_vars]
                hgt_args.append(self.height)

                LOGGER.debug('Calling Fortran function mod_hgt.height')
                data[t, :] = funcname(*hgt_args)

            LOGGER.debug('Reading reference var')
            reference_var = read_variable_from_regcmfile(
                reference_var_name,
                regcm_file,
                regcm_file_path,
                False,
            )

            attributes = copy(reference_var.attributes)
            if 'coordinates' in attributes:
                LOGGER.debug('Updating attribute "coordinates" adding "plev"')
                attributes['coordinates'] = 'plev ' + attributes['coordinates']

            pval = np.array([1], dtype=DATATYPE_AUXILIARIES)
            pval[0] = self.height
            height_var = Variable(
                name='plev',
                data=MemoryData(copy(pval)),
                dimensions=(),
                attributes={
                    'standard_name': 'air_pressure',
                    'long_name': 'pressure level',
                    'positive': 'down',
                    'units': 'Pa',
                    'axis': 'Z'
                }
            )

            return Variable(
                name="zg{:.0f}".format(self.height*0.01),
                data=MemoryData(data),
                dimensions=dimensions,
                attributes=attributes,
                times=reference_var.times,
                times_attributes=reference_var.times_attributes,
                auxiliary_vars=(height_var,),
                needs_time_bounds=False,
            )


class ComputeGeopotentialHeightOnSeveralPressures(ActionStarter):
    """
    Compute the geopotential at several pressure levels.

    Parameters:

        pressures (list of float): compute the geopotential at each pressure
            level provided.

    Returns:

        list of ``Variable``: a list of ``Variable`` objects, containing the
        height expressed as meters.
    """

    def __init__(self, pressures):
        pressures = [float(l) for l in pressures]
        self.__height_gens = [ComputeGeopotentialHeight(l) for l in pressures]

    def __call__(self, regcm_file, regcm_file_path, simulation, corrflag,
                 cordex_dir):

        new_vars = []
        for height_generator in self.__height_gens:
            LOGGER.debug(
                'Computing geopotential height on %s Pa',
                height_generator.height
            )
            new_var = height_generator(
                regcm_file,
                regcm_file_path,
                simulation,
                corrflag,
                cordex_dir
            )
            new_vars.append(new_var)

        return new_vars


class ComputeCapeCinFromVerticalProfile(ActionStarter):
    """
    Compute the integrated CAPE and CIN from vertical model profile

    Parameters:

    None

    Returns:

    list of ``Variable``: a list of ``Variable`` objects, containing the
    CAPE and CIN for each model point
    """

    def __init__(self):
        self.__generator = ComputeCapeCin( )

    def __call__(self, regcm_file, regcm_file_path, simulation, corrflag,
         cordex_dir):

        LOGGER.debug(
            'Computing CAPE and CIN'
        )
        new_vars = self.__generator(
            regcm_file,
            regcm_file_path,
            simulation,
            corrflag,
            cordex_dir
        )

        return new_vars

class ComputeCapeCin(ActionStarter):
    """
    Compute the integrated CAPE and CIN form a vertical profile

    Parameters:

    None

    Returns:

    ``Variable``: a ``Variable`` object list, containing CAPE and CIN
    meters.
    """

    def __init__(self):
        pass

    def __call__(self, regcm_file, regcm_file_path, simulation, corrflag,
                 cordex_dir):

        mcore = 'hydrostatic'
        try:
            core = Dataset(regcm_file_path).dynamical_core
            if core == 1:
                mcore = 'hydrostatic'
            elif core == 2:
                mcore = 'non-hydrostatic'
            elif core == 3:
                mcore = 'moloch'
            else:
                raise ValueError(
                   'Unknown dynamical_core : "{}"'.format(core))
        except:
            pass

        # Prepare a function to read exactly one timestep of a netcdf var
        def read_timestep(var_pointer, t_step):
            if 'time' not in var_pointer.dimensions:
                return no_time_vars_data[var_pointer.name]

            t_dim = var_pointer.dimensions.index('time')
            slicer = [slice(None) for _ in var_pointer.shape]
            slicer[t_dim] = t_step
            LOGGER.debug(
                'Reading timestep %s of variable %s',
                t_step,
                var_pointer.name
            )
            return np.array(var_pointer[slicer], dtype=DATATYPE_MAIN)

        if mcore == 'hydrostatic':
            needed_var_names = ['ps', 'ta', 'rh', 'sigma', 'ptop']
            funcname = mod_capecin.getcape_hy
        elif mcore == 'moloch':
            needed_var_names = ['ps', 'ta', 'rh', 'pai']
            funcname = mod_capecin.getcape_moloch
        else:
            needed_var_names = ['ps', 'ta', 'p0', 'rh', 'sigma', 'ptop', 'pp']
            funcname = mod_capecin.getcape_nhy

        with Dataset(regcm_file_path, 'r') as f:
            needed_vars = []
            for v in needed_var_names:
                v_names = get_regcm_variable_names(v)
                try:
                    v_pointer = get_var_with_name(v_names, f)
                except:
                    v_pointer = fake_var((100.0,))
                needed_vars.append(v_pointer)

            # This is the name of the variable that will be used to copy the
            # time steps, the time attributes and so on...
            reference_var_name = 'ps'

            LOGGER.debug('Check the total number of steps')
            time_dim = needed_vars[0].dimensions.index('time')
            time_steps = needed_vars[0].shape[time_dim]
            LOGGER.debug('%s timesteps found', time_steps)

            new_shape = needed_vars[0].shape
            new_dims = needed_vars[0].dimensions
            dimensions = list(zip(new_dims, new_shape))
            LOGGER.debug(
                'The 2 output data will have the following shape: %s',
                new_shape
            )
            LOGGER.debug(
                'The 2 output data will have the following dimensions: %s',
                dimensions
            )
            LOGGER.debug('Allocating memory')
            cape = np.empty(new_shape, dtype=DATATYPE_MAIN)
            cin = np.empty(new_shape, dtype=DATATYPE_MAIN)

            # Read the data of the variables that do not depend on time
            no_time_vars_data = {}
            for v in needed_vars:
                if 'time' not in v.dimensions:
                    LOGGER.debug('Reading data of variable %s', v.name)
                    if v.name == 'ptop':
                        v_data = np.array(v[:], dtype=DATATYPE_AUXILIARIES)
                    else:
                        v_data = np.array(v[:], dtype=DATATYPE_MAIN)
                    no_time_vars_data[v.name] = v_data

            for t in range(time_steps):
                LOGGER.debug(
                    'Computing CAPE and CIN at timestep %s',
                    t
                )
                getc_args = [read_timestep(data_v, t) for data_v in needed_vars]

                LOGGER.debug('Calling Fortran function mod_capecin.getcape')
                cape[t, :],cin[t, :] = funcname(*getc_args)

        LOGGER.debug('Reading reference var')
        reference_var = read_variable_from_regcmfile(
            reference_var_name,
            regcm_file,
            regcm_file_path,
            False,
        )

        attributes = copy(reference_var.attributes)
        attributes['standard_name'] = 'atmosphere_convective_available_potential_energy_wrt_surface'
        attributes['long_name'] = 'Convective Available Potential Energy'
        attributes['units'] = 'J kg-1'

        v1 = Variable(
            name="CAPE",
            data=MemoryData(cape),
            dimensions=dimensions,
            attributes=attributes,
            times=reference_var.times,
            times_attributes=reference_var.times_attributes,
            needs_time_bounds=False,
        )
        attributes = copy(reference_var.attributes)
        attributes['standard_name'] = 'atmosphere_convective_inhibition_wrt_surface'
        attributes['long_name'] = 'Convective INhibition'
        attributes['units'] = 'J kg-1'
        v2 = Variable(
            name="CIN",
            data=MemoryData(cin),
            dimensions=dimensions,
            attributes=attributes,
            times=reference_var.times,
            times_attributes=reference_var.times_attributes,
            needs_time_bounds=False,
        )

        return [v1, v2]


class ComputeGeoCoordinateFromGridCoordinate(ActionStarter):
    def __init__(self, var_name, grid_eastward, grid_northward, direction,
                 new_attributes=None, need_time_bounds=False):

        if direction.lower() not in ('eastward', 'northward'):
            raise ValueError(
                'The only allowed directions are: eastward, northward'
            )

        self.var_name = var_name
        self.direction = direction.lower()

        self.__x = grid_eastward
        self.__y = grid_northward

        self.new_attributes = new_attributes
        self.need_time_bounds = need_time_bounds

    def __call__(self, regcm_file, regcm_file_path, simulation, corrflag,
                 cordex_dir):

        if self.direction == 'eastward':
            v = self.__x
        else:
            v = self.__y

        LOGGER.debug('Reading variable %s as reference', v)
        reference_var = read_variable_from_regcmfile(
            v,
            regcm_file,
            regcm_file_path,
            self.need_time_bounds
        )

        if "grid_" not in reference_var.attributes['standard_name']:
            return reference_var

        attributes = copy(reference_var.attributes)

        if self.new_attributes is not None:
            for attr_name, attr_value in self.new_attributes.items():
                LOGGER.debug(
                    'Adding attribute "%s" with value "%s"',
                    attr_name,
                    attr_value
                )
                attributes[attr_name] = attr_value

        dimensions = reference_var.dimensions

        LOGGER.debug(
            'Searching the index that individuates the latitude in the variable'
        )
        lat_index = -1
        for i, dim in enumerate(dimensions):
            if dim[0] in ('iy', 'y'):
                lat_index = i
                LOGGER.debug(
                    'Latitude is dimension %s (position %s)',
                    dim[0],
                    lat_index
                )
                break
        else:
            raise ValueError(
                'Latitude dimension ("y") not found in variable {}'.format(v)
            )

        LOGGER.debug(
            'Searching the index that individuates the longitude in the '
            'variable'
        )
        lon_index = -1
        for i, dim in enumerate(dimensions):
            if dim[0] in ('jx', 'x'):
                lon_index = i
                LOGGER.debug(
                    'Longituide is dimension %s (position %s)',
                    dim[0],
                    lon_index
                )
                break
        else:
            raise ValueError(
                'Longitude dimension ("x") not found in variable {}'.format(v)
            )

        LOGGER.debug('Reading the cone from the attribute "grid_factor"')
        if 'grid_factor' in regcm_file.attributes:
            cone = float(regcm_file.attributes['grid_factor'])
            LOGGER.debug('Cone is %s', cone)
        else:
            cone = None
            LOGGER.debug('Cone not found: missing attribute')

        LOGGER.debug(
            'Reading the plon attribute from the attribute '
            '"grid_north_pole_longitude"'
        )
        if 'grid_north_pole_longitude' in regcm_file.attributes:
            plon = float(regcm_file.attributes['grid_north_pole_longitude'])
            LOGGER.debug('Plon is %s', plon)
        else:
            plon = None
            LOGGER.debug('Plon not found: missing attribute')

        if 'grid_north_pole_latitude' in regcm_file.attributes:
            plat = float(regcm_file.attributes['grid_north_pole_latitude'])
            LOGGER.debug('Plat is %s', plat)
        else:
            plat = None
            LOGGER.debug('Plat not found: missing attribute')

        # This is the function that computes the coefficients for the rotation
        def compute_rotation_coefficients(lat_lon_slice):
            LOGGER.debug('Computing transformation coefficients')
            ff = grid_to_earth_uvrotate(
                regcm_file.map_projection,
                regcm_file.xlon[lat_lon_slice],
                regcm_file.xlat[lat_lon_slice],
                regcm_file.map_projection_longitude_origin,
                regcm_file.map_projection_latitude_origin,
                cone,
                plon,
                plat
            )
            if self.direction == 'eastward':
                return ff[0], ff[1]
            elif self.direction == 'northward':
                return -ff[1], ff[0]
            else:
                raise ValueError('Invalid direction: {}'.format(self.direction))

        # Avoid repetition with the same slice
        compute_rotation_coefficients = Repeater(compute_rotation_coefficients)

        # Create now the new function to rotate the coordinates
        def transform_coords(slice_data, x, y):
            try:
                len(slice_data)
            except TypeError:
                slice_data = (slice_data,)

            slice_data = list(slice_data)

            while len(slice_data) < max(lat_index, lon_index) + 1:
                slice_data.append(slice(None))
            LOGGER.debug('Slicing data using %s', SliceStr(slice_data))

            lat_lon_slice = (slice_data[lat_index], slice_data[lon_index])
            LOGGER.debug(
                'Slicing lat and lon using %s',
                SliceStr(lat_lon_slice)
            )

            f1, f2 = compute_rotation_coefficients(lat_lon_slice)
            return f1 * x + f2 * y

        # Return the variable that uses the previous function to compute
        # its data
        data = NetCDFMultiData(
            [self.__x, self.__y],
            transform_coords,
            regcm_file_path,
            DATATYPE_MAIN,
            lock=regcm_file.lock,
            slice_as_first_argument=True
        )
        return Variable(
            name=self.var_name,
            data=data,
            dimensions=reference_var.dimensions,
            attributes=attributes,
            times=reference_var.times,
            times_attributes=reference_var.times_attributes,
            auxiliary_vars=[],
            needs_time_bounds=self.need_time_bounds,
        )
