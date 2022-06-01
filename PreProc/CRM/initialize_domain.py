""" Initialize a RegCM CRM domain file """
from __future__ import division
import os
import datetime
import netCDF4 as nc
from numpy import *

# import definitions from mod_constants.F90
import constants as c
# set constants as defined on mod_nhinterp.F90
p0 = c.stdp
ts0 = c.stdt
tlp = 50.0

def initialize_domain(  icbc_directory,
                        domain_name,
                        iy,
                        jx,
                        kz,
                        dx,
                        ptop,
                        ts0, \
                        do_clobber = False,
                        be_verbose = True, \
                        do_hydrostatic = False):
    """ Initializes a RegCM CRM domain file.

        input:
        ------

            icbc_directory  : the directory to place ICBC files (should already exist)

            domain_name     : the name of the given domain

            iy, jx, kz      : the number of domain grid points in the zonal, meridional and
                              vertical directions

            dx              : the grid resolution [m]

            ptop            : the model top pressure [cb]

            ts0             : the base state surface temperature [K]

            do_clobber      : flags whether to clobber an existing file

            be_verbose      : flags whether to print status statements

            do_hydrostatic      : flags whether to generate hydrostatic boundary conditions

        returns:
        --------

            domain_file_path    : the name and path of the created domain file

    """
    def vprint(msg):
        """ Prints a message only if in verbose mode """
        if be_verbose:
            print(msg)

    # Set the domain file's name and path
    domain_file_name = "{}_DOMAIN000.nc".format(domain_name)
    domain_file_path = os.path.abspath("{}/{}".format(icbc_directory,domain_file_name))

    
    # Open the domain file for writing
    vprint("Opening {} for writing".format(domain_file_path))
    with nc.Dataset(domain_file_path,"w",clobber = do_clobber) as fout:
        # Create dimensions
        fout.createDimension('jx',jx)
        fout.createDimension('iy',iy)
        fout.createDimension('kz',kz+1)
        fout.createDimension('soil_layer',3)
        fout.createDimension('ntex',17)

        #******************
        # Define variables
        #******************
        vjx = fout.createVariable("jx",'f8',('jx'))
        vjx.long_name = "x-coordinate in Cartesian system"
        vjx.standard_name = "projection_x_coordinate"
        vjx.units = "m"
        vjx.axis = "X"
        vjx._CoordinateAxisType = "GeoX"

        viy = fout.createVariable("iy",'f8',('iy'))
        viy.long_name = "y-coordinate in Cartesian system"
        viy.standard_name = "projection_y_coordinate"
        viy.units = "m"
        viy.axis = "Y"
        viy._CoordinateAxisType = "GeoY"

        vsigma = fout.createVariable("sigma",'f8',("kz",))
        vsigma.long_name = "Sigma at full model layers"
        vsigma.standard_name = "atmosphere_sigma_coordinate"
        vsigma.units = "1"
        vsigma.axis = "Z"
        vsigma.positive = "down"
        vsigma.formula = "p(n,k,j,i) = ptop+sigma(k)*(p0(j,i)-ptop)+ppa(n,k,j,i)"
        vsigma._CoordinateAxisType = "GeoZ"

        vptop  = fout.createVariable("ptop",'f8')
        vptop.long_name = "Pressure at model top"
        vptop.standard_name = "air_pressure"
        vptop.units = "hPa"

        vxlon = fout.createVariable("xlon",'f8',("iy", "jx"))
        vxlon.long_name = "Longitude on Cross Points"
        vxlon.standard_name = "longitude"
        vxlon.units = "degrees_east"
        vxlon.grid_mapping = "crs"

        vxlat = fout.createVariable("xlat",'f8',("iy", "jx"))
        vxlat.long_name = "Latitude on Cross Points"
        vxlat.standard_name = "latitude"
        vxlat.units = "degrees_north"
        vxlat.grid_mapping = "crs"

        vdlon = fout.createVariable("dlon",'f8',("iy", "jx"))
        vdlon.long_name = "Longitude on Dot Points"
        vdlon.standard_name = "longitude"
        vdlon.units = "degrees_east"
        vdlon.grid_mapping = "crs"

        vdlat = fout.createVariable("dlat",'f8',("iy", "jx"))
        vdlat.long_name = "latitude"
        vdlat.standard_name = "Latitude on Dot Points"
        vdlat.units = "degrees_north"
        vdlat.grid_mapping = "crs"

        vxmap = fout.createVariable("xmap",'f8',("iy", "jx"))
        vxmap.long_name = "Map Factor on Cross Points"
        vxmap.standard_name = "map_factor"
        vxmap.units = "1"
        vxmap.coordinates = "xlat xlon"
        vxmap.grid_mapping = "crs"

        vdmap = fout.createVariable("dmap",'f8',("iy", "jx"))
        vdmap.long_name = "Map Factor on Dot Points"
        vdmap.standard_name = "map_factor"
        vdmap.units = "1"
        vdmap.coordinates = "xlat xlon"
        vdmap.grid_mapping = "crs"

        vcoriol = fout.createVariable("coriol",'f8',("iy", "jx"))
        vcoriol.long_name = "Coriolis Parameter"
        vcoriol.standard_name = "coriolis_parameter"
        vcoriol.units = "s-1"
        vcoriol.coordinates = "xlat xlon"
        vcoriol.grid_mapping = "crs"

        vmask = fout.createVariable("mask",'f8',("iy", "jx"))
        vmask.long_name = "Land Mask"
        vmask.standard_name = "land_binary_mask"
        vmask.units = "1"
        vmask.coordinates = "xlat xlon"
        vmask.grid_mapping = "crs"

        vtopo = fout.createVariable("topo",'f8',("iy", "jx"))
        vtopo.long_name = "Surface Model Elevation"
        vtopo.standard_name = "surface_altitude"
        vtopo.units = "m"
        vtopo.coordinates = "xlat xlon"
        vtopo.grid_mapping = "crs"

        vlanduse = fout.createVariable("landuse",'f8',("iy", "jx"))
        vlanduse.long_name = "Landuse category as defined in BATS1E"
        vlanduse.standard_name = "land_type"
        vlanduse.units = "1"
        vlanduse.coordinates = "xlat xlon"
        vlanduse.grid_mapping = "crs"
        vlanduse.legend = [ "1  => Crop/mixed farming\n", \
                            "2  => Short grass\n", \
                            "3  => Evergreen needleleaf tree\n", \
                            "4  => Deciduous needleleaf tree\n", \
                            "5  => Deciduous broadleaf tree\n", \
                            "6  => Evergreen broadleaf tree\n", \
                            "7  => Tall grass\n", \
                            "8  => Desert\n", \
                            "9  => Tundra\n", \
                            "10 => Irrigated Crop\n", \
                            "11 => Semi-desert\n", \
                            "12 => Ice cap/glacier\n", \
                            "13 => Bog or marsh\n", \
                            "14 => Inland water\n", \
                            "15 => Ocean\n", \
                            "16 => Evergreen shrub\n", \
                            "17 => Deciduous shrub\n", \
                            "18 => Mixed Woodland\n", \
                            "19 => Forest/Field mosaic\n", \
                            "20 => Water and Land mixture\n", \
                            "21 => Urban\n", \
                            "22 => Sub-Urban" ]

        vsnowam = fout.createVariable("snowam",'f8',("iy", "jx"), fill_value = 1.e20)
        vsnowam.long_name = "Snow initial LWE in mm"
        vsnowam.standard_name = "snowfall_amount"
        vsnowam.units = "mm"
        vsnowam.coordinates = "xlat xlon"
        vsnowam.grid_mapping = "crs"

        vsmoist = fout.createVariable("smoist",'f8',("iy", "jx"), fill_value = 1.e20)
        vsmoist.long_name = "Soil Moisture"
        vsmoist.standard_name = "volume_fraction_of_water_in_soil"
        vsmoist.units = "1"
        vsmoist.coordinates = "xlat xlon"
        vsmoist.grid_mapping = "crs"

        vtexture = fout.createVariable("texture",'f8',("iy", "jx"))
        vtexture.long_name = "Texture dominant category"
        vtexture.standard_name = "soil_type"
        vtexture.units = "1"
        vtexture.coordinates = "xlat xlon"
        vtexture.grid_mapping = "crs"
        vtexture.legend = ( "1  => Sand\n", \
                            "2  => Loamy Sand\n", \
                            "3  => Sandy Loam\n", \
                            "4  => Silt Loam\n", \
                            "5  => Silt\n", \
                            "6  => Loam\n", \
                            "7  => Sandy Clay Loam\n", \
                            "8  => Silty Clay Loam\n", \
                            "9  => Clay Loam\n", \
                            "10 => Sandy Clay\n", \
                            "11 => Silty Clay\n", \
                            "12 => Clay\n", \
                            "13 => OM\n", \
                            "14 => Water\n", \
                            "15 => Bedrock\n", \
                            "16 => Other\n", \
                            "17 => No data" )

        vrmoist = fout.createVariable("rmoist",'f8',("soil_layer", "iy", "jx"), fill_value = 1.e20)
        vrmoist.long_name = "Soil Moisture"
        vrmoist.standard_name = "volume_fraction_of_water_in_soil"
        vrmoist.units = "kg m-2"
        vrmoist.coordinates = "xlat xlon"
        vrmoist.grid_mapping = "crs"

        vtexture_fraction = fout.createVariable("texture_fraction",'f8',("ntex", "iy", "jx"), fill_value = 1.e20)
        vtexture_fraction.long_name = "Texture category fraction"
        vtexture_fraction.standard_name = "soil_type_fraction"
        vtexture_fraction.units = "1"
        vtexture_fraction.coordinates = "xlat xlon"
        vtexture_fraction.grid_mapping = "crs"

        if not do_hydrostatic:
            vps0 = fout.createVariable("ps0",'f8',("iy", "jx"))
            vps0.long_name = "Reference State Surface Pressure"
            vps0.standard_name = "air_pressure"
            vps0.units = "Pa"
            vps0.coordinates = "xlat xlon"
            vps0.grid_mapping = "crs"

            vpr0 = fout.createVariable("pr0",'f8',("kz", "iy", "jx"))
            vpr0.long_name = "Reference State Pressure"
            vpr0.standard_name = "air_pressure"
            vpr0.units = "Pa"
            vpr0.coordinates = "xlat xlon"
            vpr0.grid_mapping = "crs"

            vt0 = fout.createVariable("t0",'f8',("kz", "iy", "jx"))
            vt0.long_name = "Reference State Temperature"
            vt0.standard_name = "air_temperature"
            vt0.units = "K"
            vt0.coordinates = "xlat xlon"
            vt0.grid_mapping = "crs"

            vrho0 = fout.createVariable("rho0",'f8',("kz", "iy", "jx"))
            vrho0.long_name = "Reference State Density"
            vrho0.standard_name = "air_density"
            vrho0.units = "kg m-3"
            vrho0.coordinates = "xlat xlon"
            vrho0.grid_mapping = "crs"

            vz0 = fout.createVariable("z0",'f8',("kz", "iy", "jx"))
            vz0.long_name = "Reference State Elevation"
            vz0.standard_name = "height"
            vz0.units = "m"
            vz0.coordinates = "xlat xlon"
            vz0.grid_mapping = "crs"

        vcrs  = fout.createVariable("crs",'c')
        vcrs.proj4_params = "+proj=merc +lat_ts=0.00 +lon_0=0.0 +x_0=-{:0.0f} +y_0=-{:0.0f} +ellps=sphere +a=6371229. +b=6371229. +units=m +no_defs".format(dx/2,dx/2)
        vcrs.grid_mapping_name = "mercator"
        vcrs.standard_parallel = 0.
        vcrs.latitude_of_projection_origin = 0.
        vcrs.longitude_of_projection_origin = 0.
        vcrs.semi_major_axis = 6371229.
        vcrs.inverse_flattening = 0.
        vcrs.false_easting = -dx/2.
        vcrs.false_northing = -dx/2.
        vcrs._CoordinateTransformType = "Projection"
        vcrs._CoordinateAxisTypes = "GeoX GeoY"

        vsoil_layer = fout.createVariable("soil_layer",'f8',("soil_layer",))
        vsoil_layer.long_name = "Soil layer levels"
        vsoil_layer.standard_name = "root_depth"
        vsoil_layer.units = "m"

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
        fout.boundary_smoothing = "No"
        fout.minimum_h2o_pct_for_water = 50.
        fout.smoothing_level = 1
        fout.h2o_hgt_over_water = "Yes"
        fout.intermediate_resolution = 0
        fout.landuse_fudging = "No"
        fout.texture_fudging = "No"
        fout.initialized_soil_moisture = "No"
        if not do_hydrostatic:
            fout.base_state_surface_temperature = ts0
        
        #************************
        # Write simple variables
        #************************
        vptop[:] = ptop*10
        vxlon[:] = 0.0
        vxlat[:] = 0.0
        vdlon[:] = 0.0
        vdlat[:] = 0.0
        vxmap[:] = 1.0
        vdmap[:] = 1.0
        vmask[:] = 0.0
        vtopo[:] = 0.0
        vlanduse[:] = 15
        vsnowam[:] = 0
        vsmoist[:] = 1e20
        vtexture[:] = 14
        vrmoist[:] = -1
        vtexture_fraction[:] = 0.0
        vsoil_layer[:] = 0.0

        return domain_file_path

def temppres(pressure):
    """ Given a pressure, returns the corresponding temperature for the nonhydrostatic base state.

        input:
        ------
            pressure    : atmospheric pressure [Pa]

        returns:
        -------
            temperature : nonhydrostatic reference temperature [K]

        Note: this function is modeled after the temppres function in Share/mod_nhinterp.F90

    """

    pr = pressure

    if ( pr > 22632.0): 
        tp = max(ts0 + tlp*log(pr/p0), c.tiso)
    else: 
        if ( pr > 5474.9): 
            tp = c.tiso + min(log(5474.9/pr),228.65)
        else:
            tp = 228.65+ 2.0*min(log(5474.9/pr),271.15)

    temperature = tp

    return temperature

def reference_profiles(sigma,ptoppa,ps = p0, zs = 0):
    """ Return reference profiles of temperature, pressure, and density; to be used as
        reference states for non-hydrostatic RegCM

        input:
        ------
            sigma   : sigma levels (interface levels; including 0 and 1)

            ptoppa  : model top pressure [Pa]

            ps      : the surface pressure [Pa]; defaults to 101325 Pa

            zs      : height of the lowest model level

        returns:
        --------
            pr0, t0, rho0, z0 : reference profiles of pressure, temperature, density, and height

        This code is modeled after code in nhbase() in mod_nhinterp.F90

    """

    # set some constants
    tlp = 47.7 # logp lapse rate K/ln(Pa)
    st0 = 288.15 # K
    rovg = c.rgas/c.egrav

    # set the sigma levels
    sigmah = sigma
    kz = len(sigmah)

    # set the relative surface pressure)
    psp = ps - ptoppa

    # intialize the reference vectors
    pr0 = zeros(kz)
    t0 = zeros(kz)
    rho0 = zeros(kz)
    z0 = zeros(kz)

    for k in range(kz):
        pr0[k] = psp * sigmah[k] + ptoppa
        t0[k] = temppres(pr0[k])
        rho0[k] = pr0[k]/c.rgas/t0[k]
        alnp = log(pr0[k]/(psp+ptoppa))
        z0[k] = max([-(0.5*rovg*tlp*alnp*alnp + rovg*st0*alnp),\
                       0.0])



    return pr0,t0,rho0,z0

def nonhydrostatic_pressure(sigma,ptoppa,t,q,pr0,t0,ps,ps0=p0):
    """ Return the non-hydrostatic pressure perturbation required to have an atmosphere in hydrostatic balance.

        input:
        ------
            sigma   : sigma levels (interface levels; including 0 and 1)

            ptoppa  : model top pressure [Pa]

            ps0     : the surface pressure [Pa]; defaults to 101325 Pa

            t       : temperature [K]

            q       : specific humidity [kg/kg]

            pr0     : reference pressure [Pa]

            t0      : reference temperature [K]

            ps      : surface pressure [Pa]

        returns:
        --------
            pp : pressure perturbation [Pa]

        This code is modeled after code in nhpp() in mod_nhinterp.F90.
    """
    d_one = 1.
    kxs = len(t)-1

    # calculate virtual temperature
    tv = t*(1. + c.ep1*q)

    pp = zeros(len(t))

    #
    #  Correct pressure perturbation field (mass) to be consistent
    #  with psa assuming density perturbation remains constant in
    #  lowest half-layer.   Start with pp at surface.
    #
    p0surf = ps0 + ptoppa
    psp = ps - ps0
    delp0 = p0surf - pr0[kxs]
    tvk = tv[kxs]
    tk = t[kxs]
    tvpot = (tvk - t0[kxs]) / tk
    pp[kxs] = (tvpot*delp0 + psp) / (d_one + delp0/pr0[kxs])
    for k in range(kxs-1,-1,-1):
        tvkp1 = tvk
        tvk = tv[k]
        tkp1 = tk
        tk = t[k]
        wtl = (sigma[k+1] - sigma[k]) / (sigma[k+2] - sigma[k])
        wtu = d_one - wtl
        aa = c.egrav / (pr0[k+1] - pr0[k])
        bb = c.egrav * wtl / pr0[k+1] * t0[k+1] / tkp1
        cc = c.egrav * wtu / pr0[k] * t0[k] / tk
        tvpot = wtl * ((tvkp1 - t0[k+1]) / tkp1) + wtu * ((tvk   - t0[k]) / tk  )
        pp[k] = (c.egrav * tvpot + pp[k+1] * (aa - bb)) / (aa + cc)

    print(pp)

    return pp

def nonhydrostatic_pressure2(sigma,ptoppa,t,q,pr0,t0,ps,ps0=p0):
    """ Return the non-hydrostatic pressure perturbation required to have an atmosphere in hydrostatic balance.

        input:
        ------
            sigma   : sigma levels (interface levels; including 0 and 1)

            ptoppa  : model top pressure [Pa]

            ps0     : the surface pressure [Pa]; defaults to 101325 Pa

            t       : temperature [K]

            q       : specific humidity [kg/kg]

            pr0     : reference pressure [Pa]

            t0      : reference temperature [K]

            ps      : surface pressure [Pa]

        returns:
        --------
            pp : pressure perturbation [Pa]

    """

    # initialize the pressure perturbation
    pp = zeros(len(pr0))

    # set virtual temperature
    tv = t*(1. + c.ep1*q)

    # define the linear system matrix and RHS vector
    A = zeros([len(pr0),len(pr0)])
    R = zeros(len(pr0))

    kz = len(pr0)-1

    # set the top boundary condition (P' = 0)
    A[0,0] = 1
    R[0] = 0

    # set the bottom boundary condition (P' = ps-ps0)
    A[kz,kz] = 1
    R[kz] = ps-ps0
    

    # calculate the alpha terms
    alpha = zeros(len(pr0))
    beta = zeros(len(pr0))
    pr0mid = 0.5*(pr0[1:] + pr0[:-1])
    t0mid = 0.5*(t0[1:] + t0[:-1])
    alpha[1:] = pr0mid*tv/t0mid
    alpha[0] = pr0[0]*tv[0]/t0[0]
    beta[1:] = -(pr0[:-1] - pr0[1:])
    beta[0] = beta[1]

    # set the solution
    for k in range(1,len(pr0)-1):
        A[k,k] = - 0.5/alpha[k] - 1./beta[k]
        A[k,k+1] = 1./beta[k] - 0.5/alpha[k]

        R[k] = (pr0mid[k] - alpha[k])/alpha[k]
        
    # solve the system
    pp = linalg.solve(A,R)#[::-1]

    ppmid = 0.5*(pp[1:] + pp[:-1])

    print(ppmid)

    return ppmid

