#!/usr/bin/env python

"""
This small program is computing time means of netcdf file variables
of the RegCM output file and writing a netcdf file with the means.
"""

def compute_mean(datafile,window='day'):
  """
Make the mean of an input regcm data output file.
Output file will be created in the current directory, and its name 
will follow the following convention:

Its format will be NETCDF4_CLASSIC, and all 2D+ variables will be
compressed in disk.
  """
  from netCDF4 import Dataset
  import numpy as np
  import time
  from netcdftime import datetime , num2date , utime
  import os
  from string import join

  ncf = Dataset(datafile)

  times = ncf.variables['time']
  if len(times) < 1:
    print('No timesteps in file !')
    sys.exit(0)
  dates = num2date(times[:]-0.01, units=times.units, calendar=times.calendar)
  d1 = datetime(dates[0].year,dates[0].month,dates[0].day)
  d2 = datetime(dates[-1].year,dates[-1].month,dates[-1].day)
  f1 = (repr(dates[0].year).zfill(4)+repr(dates[0].month).zfill(2)+
        repr(dates[0].day).zfill(2))
  f2 = (repr(dates[-1].year).zfill(4)+repr(dates[-1].month).zfill(2)+
        repr(dates[-1].day).zfill(2))

  pieces = os.path.basename(os.path.splitext(datafile)[0]).split('_')
 
  if window == 'day':
    if ncf.frequency == 'day' or ncf.frequency == 'mon':
      print('How to make daily mean on day or monthly dataset?')
      sys.exit(-1)
    try:
      nco = Dataset(join(pieces[0:7],'_')+'_'+window+'_'+f1+'12-'+f2+'12.nc',
                    'w', format='NETCDF4_CLASSIC')
    except:
      raise RuntimeError('Cannot open output file')
    tunit = 'days since 1949-12-01 00:00:00 UTC'
  elif window == 'mon':
    if ncf.frequency == 'mon':
      print('How to make monthly mean on monthly dataset?')
      sys.exit(-1)
    try:
      nco = Dataset(join(pieces[0:7],'_')+'_'+window+'_'+f1+'-'+f2+'.nc',
                    'w', format='NETCDF4_CLASSIC')
    except:
      raise RuntimeError('Cannot open output file')
    tunit = 'days since 1949-12-01 00:00:00 UTC'
  else:
    raise RuntimeError(
                'Unsupported time window. Only day and mon implemented')

  cdftime = utime(tunit,calendar=times.calendar)        

  for attr in ncf.ncattrs():
    if attr == 'frequency':
      nco.setncattr(attr,window)
    else:
      nco.setncattr(attr,getattr(ncf,attr))

  for dim in ncf.dimensions:
    if ( ncf.dimensions[dim].isunlimited() ):
      nco.createDimension(dim)
    else:
      nco.createDimension(dim,len(ncf.dimensions[dim]))

  if 'time_bnds' not in ncf.dimensions:
    nco.createDimension('time_bnds',2)

  tbnds = nco.createVariable('time_bnds','f8',['time','time_bnds'])
  tbnds.setncattr('units',tunit)
  tbnds.setncattr('calendar',times.calendar)

  for var in ncf.variables:
    nctype = ncf.variables[var].datatype
    if 'jx' in ncf.variables[var].dimensions:
      nco.createVariable(var,nctype,ncf.variables[var].dimensions,
                         zlib=True,complevel=9)
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
          nco.variables[var].setncattr(attr,getattr(ncf.variables[var],attr))
      if not hasbnds:
        nco.variables[var].setncattr('bounds','time_bnds')
    elif var == 'time_bnds':
      pass
    else:
      if 'time' in ncf.variables[var].dimensions:
        if 'cell_methods' in ncf.variables[var].ncattrs():
          attvalue = (getattr(ncf.variables[var],'cell_methods') +
                ' within '+ncf.frequency+' time: mean over '+window)
          nco.variables[var].setncattr('cell_methods',attvalue)
        else:
          nco.variables[var].setncattr('cell_methods',
                                 'time: mean over'+window)
      for attr in ncf.variables[var].ncattrs():
        if attr != 'cell_methods':
          nco.variables[var].setncattr(attr,getattr(ncf.variables[var],attr))

  if window == 'day':
    rc = cdftime.date2num(dates)
    ic = rc.astype(int)
    dic = np.unique(ic)
    nco.variables['time'][:] = dic+0.5
    tbnds[:,0] = dic+0.0
    tbnds[:,1] = dic+1.0
  else:
    # Assume monthly data up to now!
    rc = cdftime.date2num(num2date(np.median(times[:]),
                           units=times.units, calendar=times.calendar))
    ic = np.arange(0,len(times),dtype=np.int)
    dic = np.array([0,])
    nco.variables['time'][0] = rc
    if times.units.find('hours') > 0:
      diff = 12.0
    else:
      diff = 0.5
    x1 = num2date(times[0]-diff, units=times.units, calendar=times.calendar)
    x2 = num2date(times[-1]+diff, units=times.units, calendar=times.calendar)
    tbnds[0,0] = round(cdftime.date2num(x1))
    tbnds[0,1] = round(cdftime.date2num(x2))

  for var in ncf.variables:
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
  compute_mean(sys.argv[1], sys.argv[2])
