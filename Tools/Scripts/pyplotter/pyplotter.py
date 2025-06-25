#!/usr/bin/env python3

import os
import sys
import argparse
import yaml

from utils.storage import CacheDirectory
from utils.processing import regcm_format, cordex_format
from utils.processing import model_reader, observation_reader
from utils.plot import biasplot, p99plot

def parse_years(yys):
    if isinstance(yys[0],str) and ':' in yys[0]:
        xys = list(int(x) for x in yys[0].split(':'))
        if ( len(xys) == 2 ):
            years = range(xys[0],xys[1]+1)
        else:
            years = range(xys[0],xys[1]+1,xys[2])
    else:
        years = list(int(x) for x in yys)
    return years

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

rcmformat = config["format"]
datadir = config["datapath"]

valid_cache = not args.invalidate_cache

file_format = None
if rcmformat == "original":
    file_format = regcm_format(config["original"]["simulation"],
                               config["original"]["variables"],
                               datadir)
elif rcmformat == "cordex":
    file_format = cordex_format(datadir)
else:
    print("File reader for ",rcmformat," yet to be implemented.")
    sys.exit(0)

try:
    cache = CacheDirectory(config["cache"],valid=valid_cache)
except:
    print('Cannot access data cache.')
    sys.exit(-1)

plotdir = config["plotpath"]
what = config["plot"]
for plot in what.keys( ):
    if plot == "biasplot":
        dpi = what[plot]["dpi"]
        pname = what[plot]["name"]
        seasons = what[plot]["seasons"]
        if "bounds" in what[plot]:
            bounds = what[plot]["bounds"]
        else:
            bounds = None
        years = parse_years(what[plot]["years"])
        for var in what[plot]["variables"].keys( ):
            observations = what[plot]["variables"][var]["obs"]
            mdl_accessor = model_reader(cache,years,var,file_format)
            mds = mdl_accessor.seasonal_data(seasons)
            ods = [ ]
            for obs in observations:
                obs_accessor = observation_reader(cache,years,obs,var,
                                                  mdl_accessor.sid)
                ods.append(obs_accessor.seasonal_data(seasons, mds))
            bias_plot = biasplot(dpi=dpi,dest=plotdir,bounds=bounds)
            bias_plot.biasplot(var, years, seasons, observations,
                               mds, ods, pname)
    elif plot == "p99plot":
        dpi = what[plot]["dpi"]
        pname = what[plot]["name"]
        years = parse_years(what[plot]["years"])
        if "bounds" in what[plot]:
            bounds = what[plot]["bounds"]
        else:
            bounds = None
        for var in what[plot]["variables"].keys( ):
            observations = what[plot]["variables"][var]["obs"]
            mdl_accessor = model_reader(cache,years,var,file_format)
            mds = mdl_accessor.quantile_data(99)
            ods = [ ]
            for obs in observations:
                obs_accessor = observation_reader(cache,years,obs,var,
                                                  mdl_accessor.sid)
                ods.append(obs_accessor.quantile_data(99, mds))
            p99_plot = p99plot(dpi=dpi,dest=plotdir,bounds=bounds)
            p99_plot = p99_plot.p99plot(var,years,observations,mds,ods,pname)
    else:
        print("Plot ", plot, 'not available!')
        print("To be implemented.")
        continue

print('All Done!')
sys.exit(0)
