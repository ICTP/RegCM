#!/usr/bin/env python3

import sys
import os
import glob
import yaml
import xarray as xr
import xesmf as xe
from pathlib import Path

def season_string_to_monthlist(str):
    months = "DJFMAMJJASONDJ"
    idx = months.find(str)
    res = list(range(idx,idx+len(str)))
    if res[0] == 0:
        res[0] = 12
    if res[0] == -1:
        res[0] = 11
        res[1] = 12
    if res[-1] == 13:
        res[-1] = 1
    return res

def load_model_data(files):
    nds = xr.open_mfdataset(files, combine='nested',
                concat_dim="time", chunks={"time": 100}).unify_chunks()
    try:
        nds = nds.rename({"xlat": "lat", "xlon": "lon"})
    except:
        pass
    return nds

def load_obs_data(files):
    if ( len(files) > 1 ):
        return xr.open_mfdataset(files, combine='by_coords',
                                 chunks={"time": 100})
    else:
        return xr.open_dataset(files[0], engine='netcdf4',
                                 chunks={"time": 100})

class data_processor:

    ds = None

    def __init__(self,ds):
        if not (isinstance(ds,xr.Dataset) or isinstance(ds,xr.DataArray)):
            raise TypeError("data_processor needs Xarray Dataset or DataArray")
        self.ds = ds

    def quantile(self,q=99):
        pq_val = self.ds.quantile(q/100, dim = "time", keep_attrs = True)
        return pq_val

    def seasonal_mean(self, season):
        xtimes = self.ds['time.month'].isin(season_string_to_monthlist(season))
        ds_season = self.ds.sel(time = xtimes)
        nds = ds_season.mean(dim = "time", keep_attrs = True)
        nds = nds.assign_attrs(seasonal_average = season)
        return nds

class observation_reader:

    config = None
    cache = None
    obs = None
    var = None
    dom = None
    years = [ ]
    files = [ ]

    def __init__(self, cache, years, obs, var,
                 dom=None, config= "obsdata.yaml"):
        cpath = os.path.join(os.path.dirname(
            os.path.abspath(__file__)), config)
        with open(cpath,"r") as f:
            self.config = yaml.safe_load(f)['obsdata']
        try:
            fvar = self.config[obs][var]
        except:
            raise FileNotFoundError(
                       f"No files found for {var} in {obs}")
        path = self.config[obs]["path"]
        fvname = self.config[obs][var]["varfile"]
        cyears = (str(y) for y in years)
        if obs == "EOBS":
            names = [os.path.join(path,
                fvname+"_ens_mean_0.1deg_reg_v30.0e.nc"),]
        elif obs == "CRU":
            names = [os.path.join(path,fvname+".dat.nc"),]
        elif obs == "CPC":
            names = (os.path.join(path,
                fvname,fvname+"."+y+".nc") for y in cyears)
        elif obs == "MSWEP":
            names = (os.path.join(path,
                "monthly",y+"[0-1][0-9].nc") for y in cyears)
        elif obs == "ERA5":
            names = (os.path.join(path,"monthly",
                y,fvname+"_"+y+"_[0-1][0-9].nc") for y in cyears)
        elif obs == "EURO4M":
            names = (os.path.join(path,'EURO4M-APGD-1971-2008.nc'),)
        else:
            names = [ ]
        self.files = sorted(glob.glob(n) for n in names)[0]
        if len(self.files) == 0:
            raise FileNotFoundError(
                       f"No files found for {var} in {obs}")
        self.cache = cache
        self.years = years
        self.var = var
        self.obs = obs
        if dom is None:
            self.dom = 'TEST'
        else:
            self.dom = dom

    def listfiles(self):
        return self.files

    def load_data(self,newgrid=None):
        factor = self.config[self.obs][self.var]["factor"]
        offset = self.config[self.obs][self.var]["offset"]
        oname = self.config[self.obs][self.var]["var_ds"]
        ds = load_obs_data(self.files).rename({oname: self.var})
        if self.obs == 'CPC' or self.obs == 'ERA5' or self.obs == 'EOBS':
            try:
                ds = ds.rename({"latitude": "lat", "longitude": "lon"})
            except:
                pass
            ds = ds.assign_coords(
                        lon=((ds.lon + 180) % 360) - 180).sortby("lon")
        if newgrid is not None:
            method = "bilinear"
            method = "nearest_s2d"
            regridder = xe.Regridder(ds, newgrid, method)
            dr = ds[self.var].sel(time=ds.time.dt.year.isin(self.years))
            nds = regridder(dr, keep_attrs=True)
        else:
            nds = ds[self.var].sel(time=ds.time.dt.year.isin(self.years))
        nds *= factor
        nds += offset
        return nds

    def seasonal_data(self,seasons,newgrid=None):
        fname = (self.dom+"_"+self.obs + "_" + self.var +
            '.' + "-".join((str(x) for x in (self.years[0],self.years[-1]))))
        fseas = [ ]
        for seas in seasons:
            fseas.append(fname+"_"+seas+".nc")
        ds = self.cache.retrieve(fseas)
        if ds is None:
            proc = data_processor(self.load_data(newgrid))
            for s in seasons:
                xtmp = proc.seasonal_mean(s)
                sname = fname+"_"+s+".nc"
                self.cache.store(xtmp,sname,self.var)
                print('Observation '+self.obs + ' for ' + self.var +
                      ' season '+s+' mean created.')
                del(xtmp)
            ds = self.cache.retrieve(fseas)
        return ds

    def quantile_data(self,q=99,newgrid=None):
        fname = (self.dom+'_'+self.obs + "_q" + repr(q) + "_" + self.var +
         '.' + "-".join((str(x) for x in (self.years[0],self.years[-1])))+".nc")
        ds = self.cache.retrieve(fname)
        if ds is None:
            proc = data_processor(self.load_data(newgrid))
            xtmp = proc.quantile(q)
            self.cache.store(xtmp,fname,self.var)
            print('Observation ' + self.obs + ' for ' + self.var +
                  ' quantile '+ repr(q)+' created.')
            del(xtmp)
            ds = self.cache.retrieve(fname)
        return ds

class regcm_format:

    sid = None
    fid = None
    model_dir = None

    def __init__(self,sid,fid,model_dir):
        self.sid = sid
        self.fid = fid
        self.model_dir = model_dir

class cordex_format:

    cordex_dir = None

    def __init__(self,cordex_dir):
        self.model_dir = cordex_dir

class model_reader:

    cache = None
    model_dir = None
    var = None
    sid = None
    fid = None
    years = [ ]
    files = [ ]

    def __init__(self,cache,years,var,file_format=None):
        self.var = var
        self.cache = cache
        self.years = years
        if isinstance(file_format,regcm_format):
            self.sid = file_format.sid
            self.fid = file_format.fid[var]
            self.model_dir = file_format.model_dir
            search = list(
                f"{self.sid}_{self.fid}.{y}[0-9][0-9][0-9][0-9][0-9][0-9].nc"
                for y in years)
            for x in search:
                self.files = (self.files +
                        sorted(glob.glob(os.path.join(self.model_dir,x))))
        elif isinstance(file_format,cordex_format):
            # <activity>/<product>/<Domain>/<Institution>/<GCMModelName>/
            #    <CMIP5ExperimentName>/<CMIP5EnsembleMember>/<RCMModelName>/
            #    <BiasAdjustment>/<Frequency>/<VariableName>
            pp = Path(file_format.model_dir)
            base = os.path.join('day',var,'**',var+'_*day_')
            search = list(f"{base}{y}*.nc" for y in years)
            for x in search:
                self.files = (self.files +
                        sorted(str(p) for p in pp.rglob(x)))
            self.model_dir = file_format.model_dir
            if len(self.files) > 0:
                from netCDF4 import Dataset
                ds = Dataset(self.files[0])
                self.sid = ds.domain_id
                ds.close( )
            else:
                print(var,years)
                print('NO SID!!!')
                sys.exit(-1)
        else:
            raise NotImplementedError("Not yet implemented")

    def listfiles(self):
        return self.files

    def load_data(self):
        tmp = load_model_data(self.files)[self.var]
        # Apply unit conversion
        if self.var == "pr":
            tmp *= 86400  # Convert from kg/m²/s to mm/day
            tmp = tmp.assign_attrs(units='mm/day')
        elif self.var == "tas":
            tmp -= 273.15  # Convert from K to °C
            tmp = tmp.assign_attrs(units='Celsius')
        return tmp

    def seasonal_data(self,seasons):
        fname = (self.sid + "_" + self.var +
           '.' + "-".join((str(x) for x in (self.years[0],self.years[-1]))))
        fseas = [ ]
        for seas in seasons:
            fseas.append(fname+"_"+seas+".nc")
        ds = self.cache.retrieve(fseas)
        if ds is None:
            proc = data_processor(self.load_data())
            for s in seasons:
                xtmp = proc.seasonal_mean(s)
                sname = fname+"_"+s+".nc"
                self.cache.store(xtmp,sname,self.var)
                print('Model ' + self.sid + ' for ' + self.var +
                      ' season ' + s + ' mean created.')
            ds = self.cache.retrieve(fseas)
        return ds

    def quantile_data(self,q=99):
        fname = (self.sid + "_q" + repr(q) + "_" + self.var +
         '.' + "-".join((str(x) for x in (self.years[0],self.years[-1])))+".nc")
        ds = self.cache.retrieve(fname)
        if ds is None:
            proc = data_processor(self.load_data())
            xtmp = proc.quantile(q)
            self.cache.store(xtmp,fname,self.var)
            print('Model ' + self.sid + ' for ' + self.var +
                  ' quantile ' + repr(q) + ' created.')
            del(xtmp)
            ds = self.cache.retrieve(fname)
        return ds

if __name__ == "__main__":
    from storage_class import CacheDirectory
    for s in ["DJF","MAM","JJA","SON", "JJAS", "DJFM"]:
        print(s,season_string_to_monthlist(s))

    if False:
        cache = CacheDirectory("./tmp")
        file_format = regcm_format('MED-12',{'pr' : 'STS'},
          '/leonardo_scratch/large/userexternal/ggiulian/run/output')
        mdl_accessor = model_reader(cache,[1980],'tas',file_format)
        mds = mdl_accessor.seasonal_data(['JJA',])
        mds1 = mdl_accessor.quantile_data(99)
        obs_accessor = observation_reader(cache,[1980],'CRU','tas')
        ods = obs_accessor.seasonal_data(['JJA',], mds)
        ods1 = obs_accessor.quantile_data(99, mds)
        print(cache.content( ))
        cache.cleanup( )
