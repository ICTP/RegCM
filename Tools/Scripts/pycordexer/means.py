#!/usr/bin/env python3

"""
This small program is computing time means of netcdf file variables
of the RegCM output file and writing a netcdf file with the means.
"""

def compute_mean(datafile,window='mon'):
    """
Make the mean of an input regcm data output file.
Output file will be created in the current directory, and its name 
will follow the following convention:

Its format will be NETCDF4_CLASSIC, and all 2D+ variables will be
compressed in disk.
    """
    from netCDF4 import Dataset , num2date
    import numpy as np
    import time
    from netcdftime import datetime , utime
    import os

    names = { 'day' : 'day',
              'mon' : 'month' }

    ncf = Dataset(datafile)

    times = ncf.variables['time']
    if len(times) < 1:
        print('No timesteps in file !')
        sys.exit(0)
    if times.units.find('hours since') == 0:
        xt = times[:]-0.01
        dates = num2date(xt, units=times.units, calendar=times.calendar)
    else:
        dates = num2date(times[:], units=times.units, calendar=times.calendar)
    d1 = datetime(dates[0].year,dates[0].month,dates[0].day)
    d2 = datetime(dates[-1].year,dates[-1].month,dates[-1].day)
    f1 = (repr(dates[0].year).zfill(4)+repr(dates[0].month).zfill(2)+
          repr(dates[0].day).zfill(2))
    f2 = (repr(dates[-1].year).zfill(4)+repr(dates[-1].month).zfill(2)+
          repr(dates[-1].day).zfill(2))

    creation_date = time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                  time.localtime(time.time())),

    pieces = os.path.basename(os.path.splitext(datafile)[0]).split('_')
 
    if window == 'day':
        if ncf.frequency == 'day' or ncf.frequency == 'mon':
            print('How to make daily mean on day or monthly dataset?')
            sys.exit(-1)
        try:
            nco = Dataset(str.join('_',pieces[0:8])+
                          '_'+window+'_'+f1+'12-'+f2+'12.nc',
                          'w', format='NETCDF4_CLASSIC')
        except:
            print(str.join('_',pieces[0:8])+'_'+window+'_'+f1+'12-'+f2+'12.nc')
            raise RuntimeError('Cannot open output file')
        tunit = 'days since 1950-01-01 00:00:00'
    elif window == 'mon':
        if ncf.frequency == 'mon':
            print('How to make monthly mean on monthly dataset?')
            sys.exit(-1)
        try:
            nco = Dataset(str.join('_',pieces[0:8])+'_'+
                          window+'_'+f1+'-'+f2+'.nc',
                          'w', format='NETCDF4_CLASSIC')
        except:
            raise RuntimeError('Cannot open output file')
        tunit = 'days since 1950-01-01 00:00:00'
    else:
        raise RuntimeError(
                'Unsupported time window. Only day and mon implemented')

    cdftime = utime(tunit,calendar=times.calendar)        

    for attr in ncf.ncattrs():
        if attr == 'frequency':
            nco.setncattr(attr,window)
        elif attr == 'creation_date':
            nco.setncattr(attr,creation_date)
        else:
            nco.setncattr(attr,getattr(ncf,attr))

    for dim in ncf.dimensions:
        if ( ncf.dimensions[dim].isunlimited() or dim == 'time' ):
            nco.createDimension(dim)
        else:
            nco.createDimension(dim,len(ncf.dimensions[dim]))

    if 'bnds' not in ncf.dimensions:
        nco.createDimension('bnds',2)

    tbnds = nco.createVariable('time_bnds','f8',['time','bnds'])
    tbnds.setncattr('units',tunit)
    tbnds.setncattr('calendar',times.calendar)

    for var in ncf.variables:
        nctype = ncf.variables[var].datatype
        if ('x' in ncf.variables[var].dimensions and
            'y' in ncf.variables[var].dimensions):
            nco.createVariable(var,nctype,ncf.variables[var].dimensions,
                               shuffle=True,significant_digits=4,
                               fletcher32=True,zlib=True,complevel=9,
                               fill_value=1.0e+20)
        else:
            if var == 'time_bnds':
                pass
            else:
                nco.createVariable(var,nctype,ncf.variables[var].dimensions)
        if var == 'time':
            hasbnds = False
            for attr in ncf.variables[var].ncattrs():
              if attr == 'units':
                  nco.variables[var].setncattr('units',tunit)
              elif attr == 'bounds':
                  nco.variables[var].setncattr('bounds','time_bnds')
                  hasbnds = True
              else:
                  nco.variables[var].setncattr(attr,
                                    getattr(ncf.variables[var],attr))
            if not hasbnds:
                nco.variables[var].setncattr('bounds','time_bnds')
        elif var == 'time_bnds':
            pass
        else:
            if 'time' in ncf.variables[var].dimensions:
                if 'max' in var:
                    if window == 'day':
                        nco.variables[var].setncattr('cell_methods',
                                               'time: maximum')
                    else:
                        nco.variables[var].setncattr('cell_methods',
                                               'time: maximum within days '+
                                               'time: mean over days')
                elif 'min' in var:
                    if window == 'day':
                        nco.variables[var].setncattr('cell_methods',
                                               'time: minimum')
                    else:
                        nco.variables[var].setncattr('cell_methods',
                                               'time: minimum within days '+
                                               'time: mean over days')
                else:
                    nco.variables[var].setncattr('cell_methods',  
                                                 'time: mean')
            for attr in ncf.variables[var].ncattrs():
                if attr != 'cell_methods':
                    if attr != '_FillValue':
                        nco.variables[var].setncattr(attr,
                              getattr(ncf.variables[var],attr))

    if window == 'day':
        rc = cdftime.date2num(dates)
        ic = rc.astype(int)
        dic = np.unique(ic)
        nco.variables['time'][:] = dic+0.5
        tbnds[:,0] = dic+0.0
        tbnds[:,1] = dic+1.0
    else:
        ic = np.zeros(len(times),dtype='i4')
        for it in range(0,len(times)):
            ic[it] = int(dates[it].year*100+dates[it].month)
        dic = np.unique(ic)
        for it in range(0,len(dic)):
            indx = (ic == dic[it]) 
            if times.units.find('hours since') == 0:
                hh = (times[1]-times[0])/48.0
                nco.variables['time'][it] = np.median(times[indx])/24.0 - hh
            else:
                nco.variables['time'][it] = np.median(times[indx])
        if times.units.find('hours since') == 0:
            diff = (times[1]-times[0])/2.0
        else:
            diff = 0.5
        for it in range(0,len(dic)):
            indx = (ic == dic[it]) 
            if times.units.find('hours since') == 0:
                hh = (times[1]-times[0])/48.0
                tbnds[it,0] = (np.min(times[indx])-diff)/24.0 - hh
                tbnds[it,1] = (np.max(times[indx])+diff)/24.0 - hh
            else:
                tbnds[it,0] = np.min(times[indx])-diff
                tbnds[it,1] = np.max(times[indx])+diff

    for var in ncf.variables:
        if var == 'crs':
            continue
        if 'time' not in ncf.variables[var].dimensions:
            nco.variables[var][:] = ncf.variables[var][:]
        else:
            if var == 'time' or var == 'time_bnds':
              pass
            else:
                for it in range(0,len(dic)):
                    indx = (ic == dic[it]) 
                    nco.variables[var][it,Ellipsis] = (
                        np.nanmean(ncf.variables[var][indx,Ellipsis],axis=0,
                                 keepdims=True) )
    ncf.close()
    nco.close()

if ( __name__ == '__main__' ):
    import sys
    if len(sys.argv) < 2:
        print('Need file name')
        sys.exit(-1)
    infile = sys.argv[1]
    if ( len(sys.argv) < 3 ):
        timewin = 'mon'
    else:
        timewin = sys.argv[2]
    compute_mean(infile, timewin)
