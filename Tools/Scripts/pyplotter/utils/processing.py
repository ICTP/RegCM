#!/usr/bin/env python3

import os
import re
import numpy as np
import xarray as xr

def var_cache_dataset(ds,var,path):
    ds.to_netcdf(path=path, format = 'NETCDF4', engine = 'netcdf4',
                 encoding = {var : {"zlib": True,
                                    "complevel": 9}},
                 unlimited_dims = ['time',])

def obsconf( ):
    import os
    import yaml
    cpath = os.path.join(os.path.dirname(
        os.path.abspath(__file__)), 'obsdata.yaml')
    with open(cpath,"r") as f:
        config = yaml.safe_load(f)
    return config

def season_string_to_monthlist(str):
    months = "DJFMAMJJASONDJ"
    idx = months.find(str)
    res = list(range(idx,idx+len(str)))
    if res[0] == 0:
        res[0] = 12
    if res[-1] == 13:
        res[-1] = 1
    return res

def getmodelfnames(outdir,sid,fid,years):
    import os
    import glob
    search = list(f"{sid}_{fid}.{y}[0-9][0-9][0-9][0-9][0-9][0-9].nc"
            for y in years)
    files = [ ]
    for x in search:
       files = files+sorted(glob.glob(os.path.join(outdir,x)))
    return files

def getobsfnames(config,obs,var,years):
    import os
    import glob
    try:
        fvar = config[obs][var]
    except:
        return None
    path = config[obs]["path"]
    fvname = config[obs][var]["varfile"]
    cyears = (str(y) for y in years)
    if obs == "EOBS":
        names = [os.path.join(path,fvname+"_ens_mean_0.1deg_reg_v30.0e.nc"),]
    elif obs == "CRU":
        names = [os.path.join(path,fvname+".dat.nc"),]
    elif obs == "CPC":
        names = (os.path.join(path,fvname,fvname+"."+y+".nc") for y in cyears)
    elif obs == "MSWEP":
        names = (os.path.join(path,"monthly",y+"[0-1][0-9].nc") for y in cyears)
    elif obs == "ERA5":
        names = (os.path.join(path,"monthly",
            y,fvname+"_"+y+"_[0-1][0-9].nc") for y in cyears)
    else:
        return None
    return sorted(glob.glob(n) for n in names)[0]

def load_model_data(files, var):
    ds = xr.open_mfdataset(files, combine='nested',
            concat_dim="time", chunks={"time": 100}).unify_chunks()
    return ds[var]

def load_obs_data(files, var):
    return xr.open_mfdataset(files, combine='by_coords', chunks={"time": 100})

def compute_p99(ds):
    p99_val = ds.quantile(0.99, dim="time", keep_attrs=True)
    p99_val = p99_val.assign_attrs(quantile = 99)
    return p99_val

def compute_seasonal_mean(ds, season):
    xtimes = ds['time.month'].isin(season_string_to_monthlist(season))
    ds_season = ds.sel(time = xtimes)
    nds = ds_season.mean(dim="time",keep_attrs=True)
    nds = nds.assign_attrs(seasonal_average = season)
    return nds

def load_model_seasonal_data(cachedir,mout,sid,fid,years,var,seasons):
    fname = sid + "_" + var + '.' + "-".join((str(x) for x in years))
    fseas = [ ]
    for seas in seasons:
       sname = fname+"_"+seas+".nc"
       fseas.append(os.path.join(cachedir,sname))
    if all(os.path.exists(target) for target in fseas):
        ds = xr.open_mfdataset(fseas, combine='nested', concat_dim="time")
        print('Model seasonal mean loaded.')
    else:
        files = getmodelfnames(mout,sid,fid,years)
        if not files:
            raise FileNotFoundError(
                       f"No model files found for {var} in {mout}")
        tmp = load_model_data(files, var)
        # Apply unit conversion
        if var == "pr":
            tmp *= 86400  # Convert from kg/m²/s to mm/day
            tmp = tmp.assign_attrs(units='mm/day')
        elif var == "tas":
            tmp -= 273.15  # Convert from K to °C
            tmp = tmp.assign_attrs(units='Celsius')
        for s in seasons:
            xtmp = compute_seasonal_mean(tmp, s)
            sname = fname+"_"+s+".nc"
            var_cache_dataset(xtmp,var,os.path.join(cachedir,sname))
            print('Model season '+s+' mean created.')
        tmp.close( )
        ds = xr.open_mfdataset(fseas, combine='nested', concat_dim="time")
    return ds

def load_model_full_data(mout,sid,fid,years,var):
    files = getmodelfnames(mout,sid,fid,years)
    if not files:
       raise FileNotFoundError(
                       f"No model files found for {var} in {mout}")
    tmp = load_model_data(files, var)
    # Apply unit conversion
    if var == "pr":
        tmp *= 86400  # Convert from kg/m²/s to mm/day
        tmp = tmp.assign_attrs(units='mm/day')
    elif var == "tas":
        tmp -= 273.15  # Convert from K to °C
        tmp = tmp.assign_attrs(units='Celsius')
    return tmp

def load_model_p99_data(cachedir,mout,sid,fid,years,var):
    fname = sid + "_p99_" + var + '.' + "-".join((str(x) for x in years))
    target = os.path.join(cachedir,fname+'.nc')
    if os.path.exists(target):
        return xr.open_dataset(target)
    else:
        ds = load_model_full_data(mout,sid,fid,years,var)
        tmp = compute_p99(ds)
        var_cache_dataset(tmp,var,target)
        return tmp.to_dataset( )

def load_obs_seasonal_data(cachedir, obs, var, years, seasons, model_grid):
    fname = obs + "_" + var + '.' + "-".join((str(x) for x in years))
    fseas = [ ]
    for seas in seasons:
       sname = fname+"_"+seas+".nc"
       fseas.append(os.path.join(cachedir,sname))
    if all(os.path.exists(target) for target in fseas):
        ds = xr.open_mfdataset(fseas, combine='nested', concat_dim="time")
        print('Observation '+obs+' seasonal mean loaded.')
    else:
        config = obsconf( )['obsdata']
        files = getobsfnames(config,obs,var,years)
        if files:
            oname = config[obs][var]["var_ds"]
            ds = load_obs_data(files, oname)
            ds = ds.rename({oname: var})
            factor = config[obs][var]["factor"]
            offset = config[obs][var]["offset"]
            nds = ds[var].sel(time=ds.time.dt.year.isin(years))
            nds *= factor
            nds += offset
            if obs == 'CPC' or obs == 'ERA5' or obs == "EOBS":
                try:
                    nds = nds.rename({"latitude": "lat", "longitude": "lon"})
                except:
                    pass
                nds = nds.assign_coords(
                        lon=((nds.lon + 180) % 360) - 180).sortby("lon")
            if var == "pr":
                nds = nds.interp(lat = model_grid["xlat"],
                                 lon = model_grid["xlon"], method="nearest")
            else:
                nds = nds.interp(lat = model_grid["xlat"],
                                 lon = model_grid["xlon"], method="linear")
            for s in seasons:
                xtmp = compute_seasonal_mean(nds, s)
                sname = fname+"_"+s+".nc"
                var_cache_dataset(xtmp,var,os.path.join(cachedir,sname))
                xtmp.close( )
                print('Observation '+obs+' season '+s+' mean created.')
            nds.close( )
            ds.close( )
            ds = xr.open_mfdataset(fseas, combine='nested', concat_dim="time")
        else:
            print(obs)
            print(var)
            print(years)
            raise ValueError('No dataset found')
    return ds

def load_obs_full_data(cachedir, obs, var, years, model_grid):
    config = obsconf( )['obsdata']
    files = getobsfnames(config,obs,var,years)
    if files:
        oname = config[obs][var]["var_ds"]
        ds = load_obs_data(files, oname)
        ds = ds.rename({oname: var})
        factor = config[obs][var]["factor"]
        offset = config[obs][var]["offset"]
        nds = ds[var].sel(time=ds.time.dt.year.isin(years))
        nds *= factor
        nds += offset
        if obs == 'CPC' or obs == 'ERA5' or obs == "EOBS":
            try:
                nds = nds.rename({"latitude": "lat", "longitude": "lon"})
            except:
                pass
            nds = nds.assign_coords(
                        lon=((nds.lon + 180) % 360) - 180).sortby("lon")
        if var == "pr":
            nds = nds.interp(lat = model_grid["xlat"],
                             lon = model_grid["xlon"], method="nearest")
        else:
            nds = nds.interp(lat = model_grid["xlat"],
                                 lon = model_grid["xlon"], method="linear")
    else:
        print(obs)
        print(var)
        print(years)
        raise ValueError('No dataset found')
    return nds

def load_obs_p99_data(cachedir,obs,var,years,model_grid):
    fname = obs + "_p99_" + var + '.' + "-".join((str(x) for x in years))
    target = os.path.join(cachedir,fname+'.nc')
    if os.path.exists(target):
        return xr.open_dataset(target)
    else:
        ds = load_obs_full_data(cachedir,obs,var,years,model_grid)
        tmp = compute_p99(ds)
        var_cache_dataset(tmp,var,target)
        return tmp.to_dataset( )

if __name__ == "__main__":
    print(obsconf( ))
