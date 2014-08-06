#!/usr/bin/env python2

from netCDF4 import num2date, date2num , Dataset
from vertint import mod_vertint
from humid import mod_humid
from hgt import mod_hgt
import numpy as np
import sys
import time
import os

ICTP_Model_id = 'RegCM4'
ICTP_Model = 'ICTP-'+ICTP_Model_id
ICTP_Model_Version = 'v4'

def wspeed(u,v):
  spd = np.zeros(np.shape(u))
  spd = u*u+v*v
  return spd

def integrate(x,r):
  inr = np.zeros(np.shape(x))
  inr = np.sum(r,axis=0)
  return inr

def searchvar(vnames,nc):
  v = None
  indx = 0
  for vname in vnames:
    try:
      v = nc.variables[vname]
      break
    except:
      indx = indx + 1
      pass
  return (v , indx)

def wrapper(func,args,it):
  valist = [ ]
  for arg in args:
    if np.size(arg.shape) == 4:
      valist.append(arg[it,:,:,:])
    elif np.size(arg.shape) == 3:
      valist.append(arg[it,:,:])
    elif np.size(arg.shape) <= 2 and np.size(arg.shape) >= 1:
      valist.append(arg[:])
    else:
      valist.append(arg)
  return wrapper2(func,valist)

def wrapper2(func,args):
  return func(*args)

def copyvar(nc,name,ovar,nnvar=None,bnds=True,stsf=False,ds=0.0):
  ssdim = list(ovar.dimensions)
  try:
    ssdim[ssdim.index('iy')] = 'y'
  except:
    pass
  try:
    ssdim[ssdim.index('jx')] = 'x'
  except:
    pass
  if nnvar == None:
    nvar = name
  else:
    nvar = nnvar
  if 'ntimes' in ssdim:
    ssdim[ssdim.index('ntimes')] = 'time_bounds'
  if name == 'rcm_map' or name == 'crs':
    xvar = nc.createVariable(nvar,ovar.datatype,tuple(ssdim))
  else:
    xvar = nc.createVariable(nvar,'f8',tuple(ssdim))
  attrlist = [ ]
  for attr in ovar.ncattrs():
    attrlist.append(attr)
    if attr == 'bounds':
      if bnds:
        xvar.setncattr(attr,getattr(ovar,attr))
      else:
        pass
    else:
      if name.find('time') >= 0 and getattr(ovar,attr) == 'gregorian':
	xvar.setncattr(attr,'proleptic_gregorian')
      elif name.find('time') >= 0 and attr == 'units':
	if stsf:
          attval = 'days since 1949-12-01 00:00:00'
	else:
          attval = getattr(ovar,attr).replace(' UTC','')
	xvar.setncattr(attr,attval)
      else:
        xvar.setncattr(attr,getattr(ovar,attr))
  if nnvar == 'crs':
    if 'semi_major_axis' not in attrlist:
      xvar.setncattr('semi_major_axis',6371229.0)
    if 'inverse_flattening' not in attrlist:
      xvar.setncattr('inverse_flattening',0.0)
    if '_CoordinateTransformType' not in attrlist:
      xvar.setncattr('_CoordinateTransformType','Projection')
    if '_CoordinateAxisTypes' not in attrlist:
      xvar.setncattr('_CoordinateAxisTypes','GeoX GeoY')
    if 'false_easting' not in attrlist:
      xvar.setncattr('false_easting',-ds/2.0)
    if 'false_northing' not in attrlist:
      xvar.setncattr('false_northing',-ds/2.0)
    if 'scale_factor_at_projection_origin' not in attrlist:
      xvar.setncattr('scale_factor_at_projection_origin',1.0)
  return xvar

def unitcorrect(old,new):
  if ( old == 'mm/day' or old == 'mm day-1' or 
       old == 'Kg m-2 day-1' or old == 'kg m-2 day-1' or
       old == 'kg m-2 d-1'):
    if new == 'kg m-2 s-1':
      return 1.0/86400.0
  elif old == 'mb':
    if new == 'Pa':
      return 100.0
  elif old == 'kg m-2':
    if new == 'm':
      return 1.0/1000.0
  return 1.0

if len(sys.argv) < 9:
  print('Missing argments:')
  print('Usage:')
  print(sys.argv[0]+' datafile variable mail domain global_model'+
                    ' experiment ensemble notes [corrflag]')
  print('')
  sys.exit(-1)

datafile = sys.argv[1]
variable = sys.argv[2]
mail = sys.argv[3]
domain = sys.argv[4]
global_model = sys.argv[5]
experiment = sys.argv[6]
ensemble = sys.argv[7]
notes = sys.argv[8]
try:
  corrflag = sys.argv[9]
except:
  corrflag = 1

lookup = { 'tas' : { 'name' : ['t2m','t2avg'],
                     'long_name' : 'Near-Surface Air Temperature',
                     'units' : 'K',
                     'needbound' : False,
                     'timecorr' : { 'SRF' : 0.0,
                                    'STS' : -12.0,
                                  }
                   },
           'pr' :  { 'name' : ['tpr','pcpavg'],
                     'units' : 'kg m-2 s-1',
                     'long_name' : 'Surface Precipitation',
                     'needbound' : True,
                     'timecorr' : { 'SRF' : -1.5,
                                    'STS' : -12.0,
                                  },
                   },
           'prc' :  { 'name' : ['prc','prcv'],
                     'units' : 'kg m-2 s-1',
                     'standard_name' : 'convective_precipitation_flux',
                     'long_name' : 'Convective Precipitation',
                     'needbound' : True,
                     'timecorr' : { 'SRF' : -1.5, },
                   },
           'huss' :  { 'name' : ['qas','q2m'],
                     'units' : '1',
                     'long_name' : 'Near Surface Specific Humidity',
                     'needbound' : True,
                     'timecorr' : { 'SRF' : -1.5,
                                    'STS' : -12.0,
                                  },
                   },
           'evspsbl' :  { 'name' : ['evp'],
                     'units' : 'kg m-2 s-1',
                     'long_name' : 'Evaporation',
                     'needbound' : True,
                     'timecorr' : { 'SRF' : -1.5 },
                   },
           'mrros' :  { 'name' : ['runoff'],
                     'units' : 'kg m-2 s-1',
                     'long_name' : 'Surface Runoff',
                     'needbound' : True,
                     'timecorr' : { 'SRF' : -1.5 },
                   },
           'ps' :  { 'name' : ['ps'],
                     'units' : 'Pa',
                     'long_name' : 'Surface Air Pressure',
                     'needbound' : True,
                     'timecorr' : { 'SRF' : -1.5 ,
                                    'STS' : -12.0,
                                  },
                   },
           'tasmax' :  { 'name' : ['t2max'],
                     'units' : 'K',
                     'long_name' : 'Daily-Maximum Near-Surface Air Temperature',
                     'needbound' : True,
                     'timecorr' : { 'STS' : -12.0, },
                   },
           'tasmin' :  { 'name' : ['t2min'],
                     'units' : 'K',
                     'long_name' : 'Daily-Minimum Near-Surface Air Temperature',
                     'needbound' : True,
                     'timecorr' : { 'STS' : -12.0, },
                   },
           'sfcWindmax' :  { 'name' : ['w10max'],
                     'units' : 'm s-1',
                     'long_name' : 'Daily-Maximum Near-Surface Wind Speed',
                     'needbound' : True,
                     'timecorr' : { 'STS' : -12.0, },
                   },
           'mrro' :  { 'name' : ['mrro'],
                     'units' : 'kg m-2 s-1',
                     'long_name' : 'Total Runoff',
                     'standard_name' : 'runoff_flux',
                     'needbound' : True,
                     'timecorr' : { 'SRF' : -1.5, },
                     'formula' : [['ps',],['mrro',]],
                     'tocall' : { 'method' : integrate ,
                                  'dimension' : 2, } ,
                   },
           'sfcWind' :  { 'name' : ['sfcWind'],
                     'units' : 'm s-1',
                     'long_name' : 'Near-Surface Wind Speed',
                     'standard_name' : 'wind_speed',
                     'needbound' : True,
                     'timecorr' : { 'SRF' : -1.5, },
                     'formula' : [['uas','u10m',],['vas','v10m',]],
                     'tocall' : { 'method' : wspeed ,
                                  'dimension' : 3, } ,
                   },
           'ua850' :  { 'name' : ['ua','u'],
                     'units' : 'm s-1',
                     'long_name' : 'Eastward Wind',
                     'needbound' : False,
                     'timecorr' : { 'ATM' : 0.0, },
                     'vertint' : 850.0,
                     'vertmod' : 'linear',
                   },
           'va850' :  { 'name' : ['va','v'],
                     'units' : 'm s-1',
                     'long_name' : 'Northward Wind',
                     'needbound' : False,
                     'timecorr' : { 'ATM' : 0.0, },
                     'vertint' : 850.0,
                     'vertmod' : 'linear',
                   },
           'ta850' :  { 'name' : ['ta','t'],
                     'units' : 'K',
                     'long_name' : 'Temperature',
                     'needbound' : False,
                     'timecorr' : { 'ATM' : 0.0, },
                     'vertint' : 850.0,
                     'vertmod' : 'log',
                   },
           'hus850' :  { 'name' : ['qas','qv'],
                     'units' : '1',
                     'long_name' : 'Specific Humidity',
                     'needbound' : False,
                     'timecorr' : { 'ATM' : 0.0, },
                     'vertint' : 850.0,
                     'vertmod' : 'log',
                   },
           'hurs' :  { 'name' : ['hurs'],
                     'units' : '%',
                     'long_name' : 'Near-Surface Relative Humidity',
                     'standard_name' : 'relative_humidity',
                     'needbound' : False,
                     'timecorr' : { 'SRF' : -1.5, },
                     'formula' : [['tas','t2m',],['qas','q2m'],['ps',]],
                     'tocall' : { 'method' : mod_humid.humid2,
                                  'dimension' : 3, } ,
                   },
           'psl' :  { 'name' : ['mslp'],
                     'units' : 'Pa',
                     'long_name' : 'Sea Level Pressure',
                     'standard_name' : 'air_pressure_at_sea_level',
                     'needbound' : False,
                     'timecorr' : { 'ATM' : 0.0, },
                     'formula' : [['ps',],['t','ta'],['topo',]],
                     'tocall' : { 'method' : mod_hgt.mslp,
                                  'dimension' : 2, } ,
                   },
           'ua500' :  { 'name' : ['ua','u'],
                     'units' : 'm s-1',
                     'long_name' : 'Eastward Wind',
                     'needbound' : False,
                     'timecorr' : { 'ATM' : 0.0, },
                     'vertint' : 500.0,
                     'vertmod' : 'linear',
                   },
           'va500' :  { 'name' : ['va','v'],
                     'units' : 'm s-1',
                     'long_name' : 'Northward Wind',
                     'needbound' : False,
                     'timecorr' : { 'ATM' : 0.0, },
                     'vertint' : 500.0,
                     'vertmod' : 'linear',
                   },
           'ta500' :  { 'name' : ['ta','t'],
                     'units' : 'K',
                     'long_name' : 'Temperature',
                     'needbound' : False,
                     'timecorr' : { 'ATM' : 0.0, },
                     'vertint' : 500.0,
                     'vertmod' : 'log',
                   },
           'zg500' :  { 'name' : ['hgt'],
                     'units' : 'm',
                     'long_name' : 'Geopotential Elevation',
                     'standard_name' : 'height',
                     'needbound' : False,
                     'timecorr' : { 'ATM' : 0.0, },
                     'formula' : [['t','ta'],['ps',],['topo',],
                                  ['sigma',],['ptop',]],
                     'tocall' : { 'method' : mod_hgt.height,
                                  'dimension' : 3, } ,
                     'vertint' : 500.0,
                     'vertmod' : 'compute',
                   },
           'ua200' :  { 'name' : ['ua','u'],
                     'units' : 'm s-1',
                     'long_name' : 'Eastward Wind',
                     'needbound' : False,
                     'timecorr' : { 'ATM' : 0.0, },
                     'vertint' : 200.0,
                     'vertmod' : 'linear',
                   },
           'va200' :  { 'name' : ['va','v'],
                     'units' : 'm s-1',
                     'long_name' : 'Northward Wind',
                     'needbound' : False,
                     'timecorr' : { 'ATM' : 0.0, },
                     'vertint' : 200.0,
                     'vertmod' : 'linear',
                   },
           'ta200' :  { 'name' : ['ta','t'],
                     'units' : 'K',
                     'long_name' : 'Temperature',
                     'needbound' : False,
                     'timecorr' : { 'ATM' : 0.0, },
                     'vertint' : 200.0,
                     'vertmod' : 'log',
                   },
           'zg200' :  { 'name' : ['hgt'],
                     'units' : 'm',
                     'long_name' : 'Geopotential Elevation',
                     'standard_name' : 'height',
                     'needbound' : False,
                     'timecorr' : { 'ATM' : 0.0, },
                     'formula' : [['t','ta'],['ps',],['topo',],
                                  ['sigma',],['ptop',]],
                     'tocall' : { 'method' : mod_hgt.height,
                                  'dimension' : 3, } ,
                     'vertint' : 200.0,
                     'vertmod' : 'compute',
                   },
           'rh850' :  { 'name' : ['rh'],
                     'units' : '%',
                     'long_name' : 'Relative Humidity',
                     'standard_name' : 'relative_humidity',
                     'needbound' : False,
                     'timecorr' : { 'ATM' : 0.0, },
                     'formula' : [['ta','t'],['qv','qas'],['ps',],
                                  ['sigma',],['ptop',]],
                     'tocall' : { 'method' : mod_humid.humid1,
                                  'dimension' : 3, } ,
                     'vertint' : 850.0,
                     'vertmod' : 'compute',
                   },
         }

if variable not in lookup.keys():
  print('Nothing to do for variable '+variable)
  sys.exit(0)

lookupdim = {
          'jx' : 'x',
          'iy' : 'y',
          'm2' : 'height',
          'm10' : 'height',
          'sigma' : 'plev',
          'kz' : 'plev',
          'plev' : 'plev',
          'ntimes' : 'time_bnds',
}

lookupvar = {
          'm2' :  { 'name' : 'height',
                    'attributes' : { 'standard_name' : 'height',
                                     'long_name' : 'height',
                                     'positive' : 'up',
                                     'units' : 'm',
                                     'axis' : 'Z',
                        },
                    'default_value' : 2.0,
                 },
          'm10' : { 'name' : 'height',
                    'attributes' : { 'standard_name' : 'height',
                                     'long_name' : 'height',
                                     'positive' : 'up',
                                     'units' : 'm',
                                     'axis' : 'Z',
                        },
                    'default_value' : 10.0,
                 },
          'plev' : { 'name' : 'plev',
                    'attributes' : { 'standard_name' : 'air_pressure',
                                     'long_name' : 'pressure',
                                     'positive' : 'down',
                                     'units' : 'Pa',
                                     'axis' : 'Z',
                        },
                    'default_value' : 'plev',
                 },
          'sigma' : { 'name' : 'plev',
                    'attributes' : { 'standard_name' : 'air_pressure',
                                     'long_name' : 'pressure',
                                     'positive' : 'down',
                                     'units' : 'Pa',
                                     'axis' : 'Z',
                        },
                    'default_value' : 'plev',
                 },
          'kz' : { 'name' : 'plev',
                    'attributes' : { 'standard_name' : 'air_pressure',
                                     'long_name' : 'pressure',
                                     'positive' : 'down',
                                     'units' : 'Pa',
                                     'axis' : 'Z',
                        },
                    'default_value' : 'plev',
                 },
}

# Open input dataset
ncf = Dataset(datafile)

times = ncf.variables['time']
if len(times) < 1:
  print('No timesteps in file !')
  sys.exit(0)

xlat = ncf.variables['xlat']
xlon = ncf.variables['xlon']
iy = ncf.variables['iy']
jx = ncf.variables['jx']

stsf = False
ftype = 'SRF'
if ncf.variables.has_key('tsmax') or ncf.variables.has_key('t2max'):
  ftype = 'STS'
  stsf = True
if ncf.variables.has_key('ua') or ncf.variables.has_key('u'):
  ftype = 'ATM'
if ncf.variables.has_key('qrs') or ncf.variables.has_key('firtp'):
  ftype = 'RAD'

try:
  rcm_map = ncf.variables['rcm_map']
  orcmap = 'rcm_map'
except:
  try:
    rcm_map = ncf.variables['crs']
    orcmap = 'crs'
  except:
    rcm_map = None
    orcmap = None
try:
  timebnds = ncf.variables['time_bnds']
except:
  timebnds = None

use_formula = False

[var,idtc] = searchvar([variable,]+lookup[variable]['name'],ncf)
if var is None:
  if lookup[variable].has_key('formula'):
    compvars = [ ]
    for name in lookup[variable]['formula']:
      vv = searchvar(name,ncf)
      if vv[0] is None:
        print ('Variable not in file and missing one in '+repr(name))
        sys.exit(-1)
      if 'time' in vv[0].dimensions:
        compvars.append(vv[0])
      else:
        compvars.append(vv[0][:])
    var = compvars[0]
    use_formula = True
    if lookup[variable].has_key('vertint'):
      compvars.append(np.array([lookup[variable]['vertint'],]))
  else:
    print('Variable not in file !')
    sys.exit(-1)

correct_time = times[:]
if corrflag == 1:
  correct_time = correct_time + lookup[variable]['timecorr'][ftype]
dates = num2date(correct_time, units=times.units, calendar=times.calendar)

ff = correct_time[1]-correct_time[0]
if ff == 24:
  frequency = 'day'
  needbound = True
elif ff > 670:
  frequency = 'month'
  needbound = True
else:
  frequency = repr(ff.astype(int))+'hr'
  needbound = lookup[variable]['needbound']

olddom = getattr(ncf,'experiment')

newattr = {
  'project_id' : 'CORDEX',
  'ipcc_scenario_code' : experiment.upper().replace('.',''),
  'institute_id' : 'ICTP',
  'note' : 'The domain is larger than '+domain,
  'comment' : 'RegCM CORDEX '+olddom+' run',
  'experiment' : domain,
  'experiment_id' : experiment.replace('.',''),
  'driving_experiment' : (global_model + ', ' +
                          experiment.replace('.','') + ', ' + ensemble),
  'frequency' : frequency,
  'driving_model_id' : global_model,
  'driving_model_ensemble_member' : ensemble,
  'driving_experiment_name' : experiment.replace('.',''),
  'institution' : 'International Centre for Theoretical Physics',
  'model_id' : ICTP_Model,
  'creation_date' : time.asctime(time.localtime(time.time()-4000000)),
  'CORDEX_domain' : domain,
  'rcm_version_id' : '4.3-rc15',
  'ICTP_version_note' : notes,
  'contact' : mail,
  'product' : 'output',
}

if frequency == 'month':
  dd1 = (repr(dates[0].year).zfill(4)+repr(dates[0].month).zfill(2)+
         repr(dates[0].day).zfill(2))
  dd2 = (repr(dates[-1].year).zfill(4)+repr(dates[-1].month).zfill(2)+
         repr(dates[-1].day).zfill(2))
elif frequency == 'day':
  dd1 = (repr(dates[0].year).zfill(4)+repr(dates[0].month).zfill(2)+
         repr(dates[0].day).zfill(2)+'12')
  dd2 = (repr(dates[-1].year).zfill(4)+repr(dates[-1].month).zfill(2)+
         repr(dates[-1].day).zfill(2)+'12')
else:
  dd1 = (repr(dates[0].year).zfill(4)+repr(dates[0].month).zfill(2)+
         repr(dates[0].day).zfill(2)+repr(dates[0].hour).zfill(2)+
         repr(dates[0].minute).zfill(2))
  dd2 = (repr(dates[-1].year).zfill(4)+repr(dates[-1].month).zfill(2)+
         repr(dates[-1].day).zfill(2)+repr(dates[-1].hour).zfill(2)+
         repr(dates[-1].minute).zfill(2))

cordexname = (variable+'_'+domain+'_'+global_model+
              '_'+experiment.translate(None,'.')+'_'+
              ensemble+'_'+ICTP_Model+'_'+ICTP_Model_Version+'_'+
              frequency+'_'+dd1+'-'+dd2+'.nc')

# Open output filename
nco = Dataset(cordexname, 'w' , format='NETCDF4_CLASSIC')

# Write Global attributes
nco.setncatts(newattr)
for attr in ncf.ncattrs():
  if attr in newattr.keys():
    pass
  else:
    nco.setncattr(attr,getattr(ncf,attr))

# Define needed dimensions

newvardims = []
reducedims = None
for dim in var.dimensions:
  if ( ncf.dimensions[dim].isunlimited() ):
    nco.createDimension(dim)
    newvardims.append(dim)
  else:
    if dim == 'sigma' or dim == 'plev' or dim == 'kz':
      nco.createDimension(lookupdim[dim], 1)
      reducedims = 'plev'
      #newvardims.append(lookupdim[dim])
    elif dim == 'time_bnds' or dim == 'ntimes':
      pass
    else:
      try:
        nco.createDimension(lookupdim[dim], len(ncf.dimensions[dim]))
	if ( lookupdim[dim] != 'height' ):
          newvardims.append(lookupdim[dim])
	else:
          reducedims = 'height'
      except:
        nco.createDimension(dim,len(ncf.dimensions[dim]))
        newvardims.append(dim)

if timebnds is not None and needbound:
  try:
    nco.createDimension('time_bounds',len(ncf.dimensions['time_bounds']))
  except:
    nco.createDimension('time_bounds',len(ncf.dimensions['ntimes']))
else:
  if needbound:
    nco.createDimension('time_bnds',2)

# Define Variables

newxlat = copyvar(nco,'xlat',xlat,nnvar='lat')
newxlon = copyvar(nco,'xlon',xlon,nnvar='lon')
newjx = copyvar(nco,'jx',jx,nnvar='x')
newjx.setncattr('axis','X')
newiy = copyvar(nco,'iy',iy,nnvar='y')
newiy.setncattr('axis','Y')

newrcm_map = None
if rcm_map is not None:
  newrcm_map = copyvar(nco,orcmap,rcm_map,nnvar='crs', 
		       ds=getattr(ncf,'grid_size_in_meters'))

newtime = copyvar(nco,'time',times,bnds=needbound,stsf=stsf)

if timebnds is not None and needbound:
  newtimebnds = copyvar(nco,'time_bnds',timebnds)
else:
  if needbound:
    newtimebnds = nco.createVariable('time_bnds','f8',('time','time_bnds'))

newvar = nco.createVariable(variable,'f',tuple(newvardims),
                            shuffle=True,fletcher32='true',
			    zlib=True,complevel=9)
oldunits = ''
for attr in var.ncattrs():
  if attr in lookup[variable].keys():
    newvar.setncattr(attr,lookup[variable][attr])
    if attr == 'units':
      oldunits = getattr(var,attr)
  else:
    if attr == 'coordinates':
      newattr = getattr(var,attr)
      newattr = newattr.replace('xlat','lat')
      newattr = newattr.replace('xlon','lon')
      if reducedims is not None:
        newvar.setncattr(attr,reducedims+" "+newattr)
      else:
        newvar.setncattr(attr,newattr)
    elif attr == 'grid_mapping':
      newvar.setncattr(attr,'crs')
    else:
      newvar.setncattr(attr,getattr(var,attr))

# Search for dimension variables to be added

if lookup[variable].has_key('vertint'):
  if 'plev' in var.dimensions:
    plevs = ncf.variables['plev'][:]
    mask = np.where(plevs == lookup[variable]['vertint'])
  else:
    ps = ncf.variables['ps']
    ptop = ncf.variables['ptop'][:]
    sigma = ncf.variables['sigma'][:]

avars = []
addvar = []
for dim in var.dimensions:
  if dim in lookupvar.keys():
    addvar.append(dim)
for avar in addvar:
  bbvar = nco.createVariable(lookupvar[avar]['name'],'f8',lookupdim[avar])
  for attr in lookupvar[avar]['attributes'].keys():
    bbvar.setncattr(attr,lookupvar[avar]['attributes'][attr])
  avars.append(bbvar)

# Write variables in out file

newxlat[Ellipsis] = xlat[Ellipsis]
newxlon[Ellipsis] = xlon[Ellipsis]
newiy[Ellipsis] = iy[Ellipsis]
newjx[Ellipsis] = jx[Ellipsis]

if newrcm_map is not None:
  newrcm_map[Ellipsis] = rcm_map[Ellipsis]

for avar , avname in zip(avars,addvar):
  try:
    avar[Ellipsis] = ncf.variables[avname][Ellipsis]
  except:
    if lookupvar[avname]['default_value'] == 'plev':
      avar[Ellipsis] = np.array( lookup[variable]['vertint'], )*100.0
    else:
      avar[Ellipsis] = np.array( lookupvar[avname]['default_value'], )

if stsf:
  newtime[:] = correct_time[:]/24.0
else:
  newtime[:] = correct_time[:]
if timebnds is not None and needbound:
  if stsf:
    newtimebnds[:] = timebnds[:]/24.0
  else:
    newtimebnds[:] = timebnds[:]
else:
  if needbound:
    if stsf:
      newtimebnds[:,0] = (times[:] - ff/2)/24.0
      newtimebnds[:,1] = (times[:] + ff/2)/24.0
    else:
      newtimebnds[:,0] = times[:] - ff/2
      newtimebnds[:,1] = times[:] + ff/2

xfac = unitcorrect(oldunits,lookup[variable]['units'])

for it in range(0,np.size(correct_time)):
  if lookup[variable].has_key('vertint'):
    if 'plev' in var.dimensions:
      newvar[it,Ellipsis] = var[it,mask[0],Ellipsis] * xfac
    else:
      if lookup[variable]['vertmod'] == 'linear':
        newvar[it,Ellipsis] = mod_vertint.intlin(var[it,:,:,:],
                                 ps[it,:,:],sigma,ptop,
                                 lookup[variable]['vertint']) * xfac
      elif lookup[variable]['vertmod'] == 'log':
        newvar[it,Ellipsis] = mod_vertint.intlog(var[it,:,:,:],
                                 ps[it,:,:],sigma,ptop,
                                 lookup[variable]['vertint']) * xfac
      elif lookup[variable]['vertmod'] == 'compute':
        if not lookup[variable].has_key('formula'):
          raise RuntimeError('Cannot compute without formula!')
        func = lookup[variable]['tocall']['method']
        if lookup[variable]['tocall']['dimension'] == 2:
          newvar[it,Ellipsis] = wrapper(func,compvars,it) * xfac
        else:
          newvar[it,Ellipsis] = wrapper(func,compvars,it) * xfac
      else:
        raise RuntimeError('Unknow vertical interpolation method!')
  elif lookup[variable].has_key('formula') and use_formula:
    func = lookup[variable]['tocall']['method']
    if lookup[variable]['tocall']['dimension'] == 2:
      if variable == 'mslp':
        tmp = wrapper(func,compvars,it) * xfac
        mod_hgt.gs_filter(tmp,compvars[0][it,Ellipsis])
        newvar[it,Ellipsis] = tmp
      else:
        newvar[it,Ellipsis] = wrapper(func,compvars,it) * xfac
    else:
      newvar[it,Ellipsis] = wrapper(func,compvars,it) * xfac
  else:
    newvar[it,Ellipsis] = var[it,Ellipsis] * xfac

nco.close()
ncf.close()
