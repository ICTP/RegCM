#!/usr/bin/env python3

import os
import sys
import yaml
import glob
import numpy as np
import xarray as xr
from utils.storage import getsize, ByteSize, delete_files_in_directory
from utils.processing import load_model_seasonal_data, load_obs_seasonal_data
from utils.processing import load_model_p99_data, load_obs_p99_data
from utils.biasplot import biasplot
from utils.p99plot import p99plot
from os.path import expanduser
import cartopy
import argparse

plot_types = [ "biasplot", "p99plot" ]

parser = argparse.ArgumentParser(description = "Basic RegCM plotting",
        epilog = "Edit YAML file to configure.", allow_abbrev = True)
parser.add_argument("-c","--config", metavar = "config.yaml",
        default = "config.yaml",help = "Configuration file")
parser.add_argument("-i", "--invalidate_cache",
        help = "Clean up the cache directory",
        action = "store_true")
args = parser.parse_args()

try:
    with open(args.config,"r") as f:
        config = yaml.safe_load(f)
except:
    print("Usage:")
    print(sys.argv[0]+" [-c yaml_configuration_file] [--invalidate_cache]")
    sys.exit(-1)

modelout = config["datapath"]
rcmformat = config["format"]
what = config["plot"]
cachedir = config["cache"]
plotdir = config["plotpath"]

if any(x not in plot_types for x in what.keys( )):
    print("Unrecognized plotting requests: "+what)
    sys.exit(-1)

if not os.path.exists(cachedir):
    try:
        os.mkdir(cachedir)
        print('Created cache directory ',cachedir)
    except:
        print('Cannot crate cache directory: please check path in '+cfile)
        sys.exit(-1)
else:
    if not os.path.isdir(cachedir):
        print('Not a directory: '+cachedir)
        print('Please check cache path in '+cfile)
        sys.exit(-1)
    if args.invalidate_cache:
        delete_files_in_directory(cachedir)
    else:
        size = getsize(cachedir)
        print('Storage used so far: ',ByteSize(size))

for plot in what.keys( ):
    if plot == "biasplot":
        dpi = what[plot]["dpi"]
        pname = what[plot]["name"]
        if rcmformat == "original":
            yys = what[plot]["years"]
            if isinstance(yys[0],str) and ':' in yys[0]:
                xys = list(int(x) for x in yys[0].split(':'))
                if ( len(xys) == 2 ):
                    years = range(xys[0],xys[1]+1)
                else:
                    years = range(xys[0],xys[1]+1,xys[2])
            else:
                years = list(int(x) for x in yys)
            for var in what[plot]["variables"].keys( ):
                sid = config["original"]["simulation"]
                fid = config["original"]["variables"][var]
                sea = what[plot]["seasons"]
                msd = load_model_seasonal_data(cachedir,modelout,sid,fid,
                        years,var,sea)
                osd = [ ]
                oblist = what[plot]["variables"][var]["obs"]
                for obs in oblist:
                    osd.append(load_obs_seasonal_data(cachedir,
                        obs,var,years,sea,msd))
                biasplot(var,years,sea,oblist,msd,osd,dpi,pname,plotdir)
        else:
            print("To be implemented.")
            sys.exit(0)
    elif plot == "p99plot":
        dpi = what[plot]["dpi"]
        pname = what[plot]["name"]
        if rcmformat == "original":
            yys = what[plot]["years"]
            if isinstance(yys[0],str) and ':' in yys[0]:
                xys = list(int(x) for x in yys[0].split(':'))
                if ( len(xys) == 2 ):
                    years = range(xys[0],xys[1]+1)
                else:
                    years = range(xys[0],xys[1]+1,xys[2])
            else:
                years = list(int(x) for x in yys)
            for var in what[plot]["variables"].keys( ):
                sid = config["original"]["simulation"]
                fid = config["original"]["variables"][var]
                msd = load_model_p99_data(cachedir,modelout,sid,fid,years,var)
                osd = [ ]
                oblist = what[plot]["variables"][var]["obs"]
                for obs in oblist:
                    osd.append(load_obs_p99_data(cachedir,obs,var,years,msd))
                p99plot(var,years,oblist,msd,osd,dpi,pname,plotdir)
        else:
            print("To be implemented.")
            sys.exit(0)
