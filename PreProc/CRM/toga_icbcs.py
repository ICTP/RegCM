import netCDF4 as nc
import pickle
from numpy import *


with nc.Dataset("input/crm_test_nh_DOMAIN000.nc","r+") as fin:
    # set 0 coriolis force
    fin.variables["coriol"][:] = 0.0
    # set cartesian map
    fin.variables["dmap"][:] = 1.0
    # set cartesian map
    fin.variables["xmap"][:] = 1.0
    # set longitude/latitude
    fin.variables["xlon"][:] = 0.0
    fin.variables["dlon"][:] = 0.0
    fin.variables["xlat"][:] = 0.0
    fin.variables["dlat"][:] = 0.0


with nc.Dataset("input/crm_test_nh_ICBC.1990060100.nc","r+") as fin:
    with open('../regcm-crm/CRM/toga_sounding_interpolators.pk','rb') as fpk:
        toga_interpolators = pickle.load(fpk)

    # load the vertical coordinate
    sigma = fin.variables['sigma']

    # surface temperature
    fin.variables["ts"][:] = 305 #toga_interpolators['t'](1)

    # surface pressure
    #fin.variables["ps"][:] = (toga_interpolators['pp'](1) + 101325)/100
    fin.variables["ps"][:] = 101325

    # temperature
    fin.variables["t"][:] = ones(shape(fin.variables['t'][:]))*toga_interpolators['t'](sigma)[newaxis,:,newaxis,newaxis]

    # humidity
    fin.variables["qv"][:] = ones(shape(fin.variables['qv'][:]))*toga_interpolators['q'](sigma)[newaxis,:,newaxis,newaxis]

    # zonal wind
    fin.variables["u"][:] = ones(shape(fin.variables['u'][:]))*toga_interpolators['u'](sigma)[newaxis,:,newaxis,newaxis]

    # meridional wind
    fin.variables["v"][:] = ones(shape(fin.variables['v'][:]))*toga_interpolators['v'](sigma)[newaxis,:,newaxis,newaxis]

    # pressure perturbation
    #fin.variables["pp"][:] = ones(shape(fin.variables['pp'][:]))*toga_interpolators['pp'](sigma)[newaxis,:,newaxis,newaxis]


    # set longitude/latitude
    fin.variables["xlon"][:] = 0.0
    fin.variables["dlon"][:] = 0.0
    fin.variables["xlat"][:] = 0.0
    fin.variables["dlat"][:] = 0.0

    vars_to_zero = ["w","pp"]
    for var in vars_to_zero:
        fin.variables[var][:] = 1e-7

