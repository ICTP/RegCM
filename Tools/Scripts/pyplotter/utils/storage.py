#!/usr/bin/env python3

import os
import xarray as xr

class CacheDirectory:

    location = None
    create = True
    verbose = False
    valid = True

    def __init__(self, location = "./tmp", create = True,
                       valid = True, verbose = False):
        self.verbose = verbose
        if self.verbose:
            print('Opening Cache directory...')
        if self.create:
            if not os.path.exists(location):
                try:
                    os.makedirs(location,exist_ok=True)
                    if self.verbose:
                        print('Created cache directory '+location)
                except:
                    print('Cannot crate cache directory '+location)
                    print('Please check path!')
        if not os.path.isdir(location):
            print('Not a directory : '+location)
            print('Please check path!')
        self.location = location
        self.valid = valid
        if not self.valid:
            self.cleanup( )
        if self.verbose:
            print(self)

    def __str__(self):
        output = "Cache location : "+self.location+os.linesep
        isize = self.size( )
        if (self.valid):
            output = output + "Status is valid."+os.linesep
            if isize > 0:
              output = output + "Current size is "+self.ssize( )+os.linesep
              output = output + "Content : "+os.linesep
              for f in os.listdir(self.location):
                  if os.path.isfile(f):
                      output = output+f+os.linesep
              output = output+"####################################"
            else:
                output = output + 'Cache is empty.'
        else:
            output = output + "Status is invalid."
        return output

    def content(self):
        isize = self.size( )
        res = [ ]
        if (self.valid):
            if isize > 0:
              for f in os.listdir(self.location):
                  if os.path.isfile(os.path.join(self.location,f)):
                      res.append(f)
        return res

    def size(self):
        if self.valid:
            return sum(os.path.getsize(os.path.join(self.location,f)) for f in
                       os.listdir(self.location)
                   if os.path.isfile(os.path.join(self.location,f)))
        else:
            return 0

    def ssize(self):
        return str(ByteSize(self.size( )))

    def cleanup(self):
        if self.verbose:
            print('Requested to invalidate directory content...')
        try:
            with os.scandir(self.location) as entries:
                for entry in entries:
                    if entry.is_file():
                        os.unlink(entry.path)
            if self.verbose:
                print("All files deleted successfully.")
            self.valid = True
        except OSError:
            print("Error occurred while deleting files.")

    def __del__(self):
        if self.verbose:
            print("Closed cache directory "+self.location)
            print("Total cache size is now: "+self.ssize( ))
        self.location = None
        self.valid = False

    def store(self,ds,path,var=None):
        if isinstance(ds,list):
            for x in ds:
                self.store_one(x,path,var)
        else:
            self.store_one(ds,path,var)

    def store_one(self,ds,path,var=None):
        if ds is not None:
            cachepath = os.path.join(self.location,path)
            xds = ds
            if isinstance(ds,xr.DataArray):
                xds = ds.to_dataset( )
            encoding = None
            if var is not None:
                encoding = {var : {"zlib": True,
                                   "complevel": 4}}
            if self.verbose:
                print('Caching dataset :')
                print(cachepath,ds.info)
            ds.to_netcdf(path = cachepath,
                         format = 'NETCDF4',
                         engine = 'netcdf4',
                         encoding = encoding,
                         unlimited_dims = ['time',])

    def retrieve(self,path):
        try:
            if isinstance(path,list):
                cachepath = (os.path.join(self.location,x) for x in path)
                ds = xr.open_mfdataset(cachepath,
                                       combine='nested',
                                       concat_dim="time")
            else:
                cachepath = os.path.join(self.location,path)
                ds = xr.open_dataset(cachepath)
            return ds
        except:
            return None

class ByteSize(int):

    _KB = 1024
    _suffixes = 'B', 'KB', 'MB', 'GB', 'PB'

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        self.bytes = self.B = int(self)
        self.kilobytes = self.KB = self / self._KB**1
        self.megabytes = self.MB = self / self._KB**2
        self.gigabytes = self.GB = self / self._KB**3
        self.petabytes = self.PB = self / self._KB**4
        *suffixes, last = self._suffixes
        suffix = next((
            suffix
            for suffix in suffixes
            if 1 < getattr(self, suffix) < self._KB
        ), last)
        self.readable = suffix, getattr(self, suffix)

        super().__init__()

    def __str__(self):
        return self.__format__('.2f')

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, super().__repr__())

    def __format__(self, format_spec):
        suffix, val = self.readable
        return '{val:{fmt}} {suf}'.format(val=val, fmt=format_spec, suf=suffix)

    def __sub__(self, other):
        return self.__class__(super().__sub__(other))

    def __add__(self, other):
        return self.__class__(super().__add__(other))

    def __mul__(self, other):
        return self.__class__(super().__mul__(other))

    def __rsub__(self, other):
        return self.__class__(super().__sub__(other))

    def __radd__(self, other):
        return self.__class__(super().__add__(other))

    def __rmul__(self, other):
        return self.__class__(super().__rmul__(other))

if __name__ == "__main__":
    import numpy as np
    import pandas as pd
    cachedir = CacheDirectory("./tmp",create=True,valid=True,verbose=True)
    np.random.seed(0)
    temperature = 15 + 8 * np.random.randn(2, 3, 4)
    precipitation = 10 * np.random.rand(2, 3, 4)
    lon = [-99.83, -99.32]
    lat = [42.25, 42.21]
    instruments = ["manufac1", "manufac2", "manufac3"]
    time = pd.date_range("2014-09-06", periods=4)
    reference_time = pd.Timestamp("2014-09-05")
    ds = xr.Dataset(
            data_vars=dict(
              temperature=(["loc", "instrument", "time"], temperature),
              precipitation=(["loc", "instrument", "time"], precipitation),
            ),
            coords=dict(
              lon=("loc", lon),
              lat=("loc", lat),
              instrument=instruments,
              time=time,
              reference_time=reference_time,
            ),
            attrs=dict(description="Weather related data."),
         )
    cachedir.store([ds,],'pippo.nc')
    cachedir.store(None,'peppe.nc')
    print('Cache content: ',cachedir.content( ))
    del(cachedir)
    print('Deleted cache object, creating new object.')
    cachedir = CacheDirectory("./tmp",create=False,valid=True,verbose=False)
    print('Retrieving dataset from existing cache.')
    ds = cachedir.retrieve('pippo.nc')
    print(ds.info)
    if cachedir.retrieve('pappo.nc') is None:
        print('Good, cannot retrieve non exixting resource.')
    cachedir.cleanup( )
    del(cachedir)
    os.rmdir("./tmp")
