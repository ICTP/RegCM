""" Creates bare-bones ICBC files for RegCM """
from __future__ import division
import os
import datetime
import netCDF4 as nc
from numpy import *

def initialize_icbc(    icbc_directory, \
                        domain_name, \
                        iy, \
                        jx, \
                        kz, \
                        dx, \
                        ptop, \
                        current_file_start, \
                        bc_start_date, \
                        nstep, \
                        bc_frequency, \
                        do_clobber = False, \
                        be_verbose = True, \
                        do_hydrostatic = False):
    """ Initializes a RegCM CRM domain file.

        input:
        ------

            icbc_directory      : the directory to place ICBC files (should already exist)

            domain_name         : the name of the given domain

            iy, jx, kz          : the number of domain grid points in the zonal, meridional and
                                  vertical directions

            dx                  : the grid resolution [m]

            ptop                : the model top pressure [cb]

            current_file_start  : the date of the start of the given file (a datetime.date object)

            bc_start_date       : the date of the start of the ICBCs (a datetime.date object)

            nstep               : the number of timesteps in the given file

            bc_frequency        : the timespan between boundary conditions (hours)

            do_clobber          : flags whether to clobber an existing file

            be_verbose          : flags whether to print status statements

            do_hydrostatic      : flags whether to generate hydrostatic boundary conditions

        returns:
        --------

            icbc_file_path      : the name and path of the created icbc file

    """
    def vprint(msg):
        """ Prints a message only if in verbose mode """
        if be_verbose:
            print(msg)
            
    # generate the ICBC file name and path
    icbc_file_name = "{}_ICBC.{:04}{:02}{:02}00.nc".format(domain_name, \
                                                           current_file_start.year, \
                                                           current_file_start.month, \
                                                           1)
    icbc_file_path = os.path.abspath("{}/{}".format(icbc_directory,icbc_file_name))


    # Open the icbc file for writing
    vprint("Opening {} for writing".format(icbc_file_path))
    with nc.Dataset(icbc_file_path,"w",clobber = do_clobber) as fout:
        # Create dimensions
        fout.createDimension('jx',jx)
        fout.createDimension('iy',iy)
        fout.createDimension('kz',kz)
        fout.createDimension('time',None)
        fout.createDimension('soil_layer',3)
        fout.createDimension('ntex',17)

        #******************
        # Define variables
        #******************
        vjx = fout.createVariable("jx",'f4',('jx'))
        vjx.long_name = "x-coordinate in Cartesian system"
        vjx.standard_name = "projection_x_coordinate"
        vjx.units = "m"
        vjx.axis = "X"
        vjx._CoordinateAxisType = "GeoX"

        viy = fout.createVariable("iy",'f4',('iy'))
        viy.long_name = "y-coordinate in Cartesian system"
        viy.standard_name = "projection_y_coordinate"
        viy.units = "m"
        viy.axis = "Y"
        viy._CoordinateAxisType = "GeoY"

        vsigma = fout.createVariable("sigma",'f4',("kz",))
        vsigma.long_name = "Sigma at full model layers"
        vsigma.standard_name = "atmosphere_sigma_coordinate"
        vsigma.units = "1"
        vsigma.axis = "Z"
        vsigma.positive = "down"
        vsigma.formula = "p(n,k,j,i) = ptop+sigma(k)*(p0(j,i)-ptop)+ppa(n,k,j,i)"
        vsigma._CoordinateAxisType = "GeoZ"

        vptop  = fout.createVariable("ptop",'f4')
        vptop.long_name = "Pressure at model top"
        vptop.standard_name = "air_pressure"
        vptop.units = "hPa"

        vxlon = fout.createVariable("xlon",'f4',("iy", "jx"))
        vxlon.long_name = "Longitude on Cross Points"
        vxlon.standard_name = "longitude"
        vxlon.units = "degrees_east"
        vxlon.grid_mapping = "crs"

        vxlat = fout.createVariable("xlat",'f4',("iy", "jx"))
        vxlat.long_name = "Latitude on Cross Points"
        vxlat.standard_name = "latitude"
        vxlat.units = "degrees_north"
        vxlat.grid_mapping = "crs"

        vdlon = fout.createVariable("dlon",'f4',("iy", "jx"))
        vdlon.long_name = "Longitude on Dot Points"
        vdlon.standard_name = "longitude"
        vdlon.units = "degrees_east"
        vdlon.grid_mapping = "crs"

        vdlat = fout.createVariable("dlat",'f4',("iy", "jx"))
        vdlat.long_name = "latitude"
        vdlat.standard_name = "Latitude on Dot Points"
        vdlat.units = "degrees_north"
        vdlat.grid_mapping = "crs"

        vmask = fout.createVariable("mask",'f4',("iy", "jx"))
        vmask.long_name = "Land Mask"
        vmask.standard_name = "land_binary_mask"
        vmask.units = "1"
        vmask.coordinates = "xlat xlon"
        vmask.grid_mapping = "crs"

        vtopo = fout.createVariable("topo",'f4',("iy", "jx"))
        vtopo.long_name = "Surface Model Elevation"
        vtopo.standard_name = "surface_altitude"
        vtopo.units = "m"
        vtopo.coordinates = "xlat xlon"
        vtopo.grid_mapping = "crs"

        vps = fout.createVariable("ps",'f4',("time", "iy", "jx"))
        vps.long_name = "Surface pressure"
        vps.standard_name = "surface_air_pressure"
        vps.units = "hPa"
        vps.coordinates = "xlat xlon"
        vps.grid_mapping = "crs"
        vps.cell_methods = "time: point"

        vts = fout.createVariable("ts",'f4',("time", "iy", "jx"))
        vts.long_name = "Surface Temperature"
        vts.standard_name = "surface_temperature"
        vts.units = "K"
        vts.coordinates = "xlat xlon"
        vts.grid_mapping = "crs"
        vts.cell_methods = "time: point"
        
        vt = fout.createVariable("t",'f4',("time", "kz", "iy", "jx"))
        vt.long_name = "Temperature"
        vt.standard_name = "air_temperature"
        vt.units = "K"
        vt.coordinates = "xlat xlon"
        vt.grid_mapping = "crs"
        vt.cell_methods = "time: point"

        vqv = fout.createVariable("qv",'f4',("time", "kz", "iy", "jx"))
        vqv.long_name = "Water vapor mixing ratio"
        vqv.standard_name = "humidity_mixing_ratio"
        vqv.units = "kg kg-1"
        vqv.coordinates = "xlat xlon"
        vqv.grid_mapping = "crs"
        vqv.cell_methods = "time: point"

        vu = fout.createVariable("u",'f4',("time", "kz", "iy", "jx"))
        vu.long_name = "Zonal component (westerly) of wind"
        vu.standard_name = "eastward_wind"
        vu.units = "m s-1"
        vu.coordinates = "dlat dlon"
        vu.grid_mapping = "crs"
        vu.cell_methods = "time: point"

        vv = fout.createVariable("v",'f4',("time", "kz", "iy", "jx"))
        vv.long_name = "Meridional component (southerly) of wind"
        vv.standard_name = "northward_wind"
        vv.units = "m s-1"
        vv.coordinates = "dlat dlon"
        vv.grid_mapping = "crs"
        vv.cell_methods = "time: point"

        if not do_hydrostatic:
            vw = fout.createVariable("w",'f4',("time", "kz", "iy", "jx"))
            vw.long_name = "Vertical wind"
            vw.standard_name = "upward_air_velocity"
            vw.units = "m s-1"
            vw.coordinates = "xlat xlon"
            vw.grid_mapping = "crs"
            vw.cell_methods = "time: point"

            vwtop = fout.createVariable("wtop",'f4',("time", "iy", "jx"))
            vwtop.long_name = "Model top vertical velocity"
            vwtop.standard_name = "upward_air_velocity"
            vwtop.units = "m s-1"
            vwtop.coordinates = "xlat xlon"
            vwtop.grid_mapping = "crs"
            vwtop.cell_methods = "time: point"

            vpp = fout.createVariable("pp",'f4',("time", "kz", "iy", "jx"))
            vpp.long_name = "Pressure perturbation"
            vpp.standard_name = "difference_of_air_pressure_from_model_reference"
            vpp.units = "Pa"
            vpp.coordinates = "xlat xlon"
            vpp.grid_mapping = "crs"
            vpp.cell_methods = "time: point"

        vtime = fout.createVariable("time",'f4',("time"))
        vtime.long_name = "time"
        vtime.standard_name = "time"
        vtime.units = "hours since 1949-12-01 00:00:00 UTC"
        vtime.calendar = "gregorian"
        
        vcrs  = fout.createVariable("crs",'c')
        vcrs.proj4_params = "+proj=merc +lat_ts=0.00 +lon_0=0.0 +x_0=-{:0.0f} +y_0=-{:0.0f} +ellps=sphere +a=6371229. +b=6371229. +units=m +no_defs".format(dx/2,dx/2)
        vcrs.grid_mapping_name = "mercator"
        vcrs.standard_parallel = 0.
        vcrs.latitude_of_projection_origin = 0.
        vcrs.longitude_of_projection_origin = 0.
        vcrs.semi_major_axis = 6371229.
        vcrs.inverse_flattening = 0.
        vcrs.false_easting = -dx/2
        vcrs.false_northing = -dx/2
        vcrs._CoordinateTransformType = "Projection"
        vcrs._CoordinateAxisTypes = "GeoX GeoY"

        fout.title = "ICTP Regional Climatic model V5"
        fout.institution = "ICTP"
        fout.source = "RegCM Model output file"
        fout.Conventions = "CF-1.4"
        fout.references = "http://gforge.ictp.it/gf/project/regcm"
        fout.model_revision = "tag-4.6.1"
        fout.history = "{} : Created by RegCM crm_icbc program".format(datetime.datetime.now())
        fout.experiment = domain_name
        fout.projection = "NORMER"
        fout.grid_size_in_meters = dx
        fout.latitude_of_projection_origin = 0.
        fout.longitude_of_projection_origin = 180.
        fout.grid_factor = 0.
        fout.atm_source = "TOGA-COARE"
        
        #************************
        # Write simple variables
        #************************
        vptop[:] = ptop*10
        vxlon[:] = 0.0
        vxlat[:] = 0.0
        vdlon[:] = 0.0
        vdlat[:] = 0.0
        vmask[:] = 0.0
        vtopo[:] = 0.0
        vu[:] = 1.e-7
        vv[:] = 1.e-7
        if not do_hydrostatic:
            vw[:] = 1.e-7
            vpp[:] = 1.e-7
            vwtop[:] = 1.e-7

        # create the timestamps for each timestep in the file
        dates = [current_file_start + datetime.timedelta(hours=n*bc_frequency) for n in range(nstep)]
        # convert these to a julian type date
        times = nc.date2num(dates,vtime.units)
        # write the times
        vtime[:] = times




    return icbc_file_path
