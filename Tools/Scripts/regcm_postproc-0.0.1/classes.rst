Modules within regcm_postproc
=============================

.. py:module:: read

.. py:class:: Reader(pattern)

    A generic class that implements input from netCDF file, requiring the name of the variable in netCDF and (optionally) imposed limits on the data. Can be subclassed, representing netCDF files from different sources.

    :param str pattern: a glob pattern that can be expanded to one or several netCDF files

    .. py:attribute:: crd_names

    The dictionary of variables representing coordinates in the netCDF file. `crd_names` keys should be 'lat' and 'lon'

    .. py:method:: get_value(self, var_name, [imposed_limits=None, [latlon_limits=None]])

    Gets the data from netCDF files, optionally within given time and coordinate limits. The limits must be the dictionaries of tuples of 2 elements, representing minimum and maximum limits.

    :param str var_name: The name of a variable, as given in the netCDF file
    :param dict imposed_limits: A dictionary of various constraints imposed on data (not including coordinates)
    :param dict latlon_limits: A dictionary of coordinate constraints (keys must be 'lat' and 'lon')
    :return: the Value object
    :rtype: Value

.. py:class:: RegCMReader(pattern)

    A class that can read from RegCM files

.. py:class:: CRUReader(pattern)

    A class that can read from CRU files

.. py:module:: value

.. py:class:: Value([data=None[, title = None[, units=None[, dims=None[, dim_names=None[, latlon=None[, limits=None[, latlon_limits=None]]]]]]]])

    A class representing the data read from netCDF files. If all parameters passed to class constructor are equal to `None`, this will create an empty value.

    :param ndarray data: actual data read from netCDF files
    :param str title: the name of the data variable (usually is stored in `long_name` netCDF variable attribute)
    :param str units: measurement units of the data variable
    :param tuple dims: dimensions of data variable, a tuple of arrays
    :param tuple dim_names: dimension names of data variable, a tuple of strings
    :param list latlon: a list of 2 2-d arrays representing latitude-longitude mesh
    :param dict limits: a dictionary of limits (lists of min and max values) of the data dimensions
    :param dict latlon_limits: a dictionary of coordinate limits (keys are 'lat' and 'lon')

    .. py:method:: update(self, value)

    Updates itself with the data from another instance of a Value class

    :param Value value: data for update
    :return: the updated Value

    .. py:method:: regrid(self, latlon)

    Regrids itself to provided latlon grid using RectBivariateSpline from scipy

    :param list latlon: a list of 2 2-d arrays representing new latitude-longitude mesh

    .. py:method:: get_limits(self, name):

    Returns data limits by name (can't be coordinate limits)

    :param str name: the dimension name for which the limits are requested
    :return: the limits for the dimension


    .. py:method:: get_latlonlimits(self)

    Returns coordinate limits for the value

    .. py:method:: mean(self):

    Returns mean value over time axis

    :rtype: Value

    .. py:method:: __sub__(self, other)

    Returns a value with value.data subtracted from self.data

    :param Value other: a Value with data to subtract
    :rtype: Value

    .. py:method:: __abs__(self)

    Returns a value with absolute values of data

    :rtype: Value

    .. py:method:: to_C(self)

    Converts its data to degrees Celsius

    .. py:method:: to_K(self)

    Converts its data to Kelvins

    .. py:method:: m_average(self, names, data, times):

    Returns monthly average





