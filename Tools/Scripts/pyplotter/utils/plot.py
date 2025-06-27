#!/usr/bin/env pyton3

import os
import yaml
import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors

class p99plot:

    config = None
    dpi = 100
    dest = './plots'

    def __init__(self,config='plotspec.yaml',
                      dpi=100,
                      bounds=None,
                      dest='./plots',
                      cartopy_cache='~/.cartopy-data'):
        cartopy.config['data_dir'] = os.path.expanduser(cartopy_cache)
        cpath = os.path.join(os.path.dirname(
            os.path.abspath(__file__)), config)
        with open(cpath,"r") as f:
            self.config = yaml.safe_load(f)
        self.dpi = dpi
        self.dest = dest
        self.bounds = bounds
        os.makedirs(dest,exist_ok=True)

    def p99plot(self,var,years,obs,msd,osd,pname):
        nrows = 1
        ncols = len(obs)
        print(f"P99 : Processing variable: {var}")

        colors = self.config["p99_hex_colors"][var]
        levels = self.config["p99_levels"].get(var, np.linspace(-8, 8, 11))

        label = ('P99 Bias '+msd[var].long_name + ' [' +
                 msd[var].units + ']')

        # Define colormap
        if isinstance(colors, list):
            custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom",
                          colors)
        else:
            # Predefined Matplotlib colormap
            custom_cmap = plt.get_cmap(colors)

        # Plotting biases
        if self.bounds:
            bounds = self.bounds
            if bounds[0] > 90.0 and bounds[1] < 0.0:
                bounds[1] = 360.0 + bounds[1]
            lonrange = (bounds[1]-bounds[0])
            latrange = (bounds[3]-bounds[2])
            xratio = abs(lonrange/latrange)
            yratio = 1.0/xratio
        else:
            bounds = [msd.lon.min(), msd.lon.max(),
                      msd.lat.min(), msd.lat.max()]
            if bounds[0] > 90.0 and bounds[1] < 0.0:
                bounds[1] = 360.0 + bounds[1]
            lonrange = (bounds[1]-bounds[0])
            latrange = (bounds[3]-bounds[2])
            xratio = lonrange/28
            yratio = latrange/28

        xpos = bounds[0] + 0.20*lonrange
        ypos = bounds[2] + 0.20*latrange
        ydim = yratio*nrows*3
        xdim = xratio*ncols*3
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(xdim,ydim),
                           subplot_kw={"projection": ccrs.PlateCarree()})

        for col, obs_name in enumerate(obs):
            bias = msd[var] - osd[col][var]
            if ncols == 1:
                ax = axes
            else:
                ax = axes[col]
            ax.set_extent(bounds, crs=ccrs.PlateCarree())
            cs = ax.contourf(msd.lon, msd.lat, bias.squeeze( ),
                    levels=levels, cmap=custom_cmap, extend="both")
            ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor="grey")
            ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
            if col == 0:
                ax.gridlines(draw_labels={"left":"y","bottom":"x"},
                    dms=True, x_inline=False, y_inline=False)
            else:
                ax.gridlines(draw_labels={"bottom":"x"},
                    dms=True, x_inline=False, y_inline=False)
            ax.set_title(f"P99 Bias {pname} - {obs_name}", fontsize=8*xratio)

        plt.subplots_adjust(wspace=0.02, hspace=0.02)
        cbar = fig.colorbar(cs, ax=axes,
               orientation="vertical", shrink=0.8,
               fraction = 0.02*yratio, pad = 0.03)
        cbar.set_label(label,size=10*yratio)

        cyears = "-".join(str(y) for y in (years[0],years[-1]))
        pfile = pname.replace(" ","_")
        out_file = os.path.join(self.dest, f"{pfile}_p99_{var}_{cyears}.png")
        plt.savefig(out_file, bbox_inches="tight", dpi=int(self.dpi))
        plt.close(fig)

class biasplot:

    config = None
    dpi = 100
    dest = './plots'

    def __init__(self,config='plotspec.yaml',
                      dpi=100,
                      bounds = None,
                      dest='./plots',
                      cartopy_cache='~/.cartopy-data'):
        cartopy.config['data_dir'] = os.path.expanduser(cartopy_cache)
        cpath = os.path.join(os.path.dirname(
            os.path.abspath(__file__)), config)
        with open(cpath,"r") as f:
            self.config = yaml.safe_load(f)
        self.dpi = dpi
        self.dest = dest
        self.bounds = bounds
        os.makedirs(dest,exist_ok=True)

    def biasplot(self,var,years,seas,obs,msd,osd,pname):
        nrows = len(seas)
        ncols = len(obs)
        print(f"BIAS : Processing variable: {var}")

        colors = self.config["bias_hex_colors"][var]
        levels = self.config["bias_levels"].get(var, np.linspace(-8, 8, 11))
        label = (msd[var].long_name + ' Bias [' +
                 msd[var].units + ']')

        # Define colormap
        if isinstance(colors, list):
            custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom",
                          colors)
        else:
            # Predefined Matplotlib colormap
            custom_cmap = plt.get_cmap(colors)

        # Plotting biases
        if self.bounds:
            bounds = self.bounds
            if bounds[0] > 90.0 and bounds[1] < 0.0:
                bounds[1] = 360.0 + bounds[1]
            lonrange = (bounds[1]-bounds[0])
            latrange = (bounds[3]-bounds[2])
            xratio = abs(lonrange/latrange)
            yratio = 1.0/xratio
        else:
            bounds = [msd.lon.min(), msd.lon.max(),
                      msd.lat.min(), msd.lat.max()]
            if (bounds[0] > 90.0 and bounds[1] < 0.0):
                bounds[1] = 360.0 + bounds[1]
            lonrange = (bounds[1]-bounds[0])
            latrange = (bounds[3]-bounds[2])
            xratio = lonrange/28
            yratio = latrange/28

        xpos = bounds[0] + 0.20*lonrange
        ypos = bounds[2] + 0.20*latrange
        ydim = yratio*nrows*3
        xdim = xratio*ncols*3
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(xdim,ydim),
                           subplot_kw={"projection": ccrs.PlateCarree()})

        for row, season in enumerate(seas):
            rcm_seasonal = msd[var][row]
            for col, obs_name in enumerate(obs):
                obs_seasonal = osd[col][var][row]
                bias = rcm_seasonal - obs_seasonal
                if nrows == 1 and ncols == 1:
                    ax = axes
                elif ncols == 1:
                    ax = axes[row]
                else:
                    ax = axes[row, col]
                ax.set_extent(bounds, crs=ccrs.PlateCarree())
                ax.add_feature(cfeature.BORDERS,
                        linewidth=0.5, edgecolor='grey')
                ax.add_feature(cfeature.COASTLINE, linewidth=0.5) 
                cs = ax.contourf(rcm_seasonal.lon, rcm_seasonal.lat,
                                 bias.squeeze(), levels=levels,
                                 cmap=custom_cmap, extend="both")
                if row == 0:
                    ax.set_title(f"{pname} - {obs_name}", fontsize=8*xratio)
                # Add rectangle with season label
                if col == 0:
                    ax.text(xpos, ypos, f"{season}",
                        fontsize=10, va="center", ha="center", color='grey',
                        bbox=dict(facecolor="white", alpha=0.5,
                        edgecolor="grey", boxstyle="round,pad=0.3"))
                    if row == (len(seas)-1):
                        ax.gridlines(draw_labels={"left":"y","bottom":"x"},
                            dms=True, x_inline=False, y_inline=False)
                    else:
                        ax.gridlines(draw_labels={"left":"y"},
                            dms=True, x_inline=False, y_inline=False)
                else:
                    if row == (len(seas)-1):
                        ax.gridlines(draw_labels={"bottom":"x"},
                            dms=True, x_inline=False, y_inline=False)
                    else:
                        ax.gridlines(draw_labels=False,
                            dms=True, x_inline=False, y_inline=False)

        plt.subplots_adjust(wspace=0.02, hspace=0.02)
        if nrows == 1 and ncols == 1:
            cbar = fig.colorbar(cs, ax = axes, shrink = 0.8,
                    orientation = "vertical",
                    fraction = 0.05*yratio, pad = 0.03)
            cbar.set_label(label,size=7*yratio)
        else:
            if ncols == 1:
                cbar = fig.colorbar(cs, ax = axes, shrink = 0.8,
                    orientation = "vertical",
                    fraction = 0.04*yratio, pad = 0.03)
            else:
                cbar = fig.colorbar(cs, ax = axes, shrink = 0.8,
                    orientation = "vertical",
                    fraction = 0.02*yratio, pad = 0.03)
            cbar.set_label(label,size=10*yratio)

        cyears = "-".join(str(y) for y in (years[0],years[-1]))
        pfile = pname.replace(" ","_")
        out_file = os.path.join(self.dest,
                f"{pfile}_bias_{var}_{cyears}_4s.png")
        plt.savefig(out_file, bbox_inches="tight", dpi=int(self.dpi))
        plt.close(fig)

if __name__ == "__main__":
    from processing_class import regcm_format
    file_format = regcm_format('MED-12',{"pr": 'STS'},
          '/leonardo_scratch/large/userexternal/ggiulian/run/output')
    if False:
        from storage_class import CacheDirectory
        from processing_class import model_reader, observation_reader
        cache = CacheDirectory("./tmp")
        years = [1980,]
        seasons = ['DJF','MAM','JJA','SON']
        observations = ['ERA5','EOBS','CRU','MSWEP','CPC']
        mdl_accessor = model_reader(cache,years,'pr',file_format)
        mds = mdl_accessor.seasonal_data(seasons)
        ods = [ ]
        obs_accessor = { }
        for obs in observations:
            obs_accessor[obs] = observation_reader(cache,years,obs,'pr')
            ods.append(obs_accessor[obs].seasonal_data(seasons, mds))
        bias_plot = biasplot( )
        bias_plot.biasplot('pr', years, seasons, observations,
                mds, ods, 'RegCM5')
        mds = mdl_accessor.quantile_data(99)
        ods = [ ]
        for obs in observations:
            ods.append(obs_accessor[obs].quantile_data(99, mds))
        p99_plot = p99plot( )
        p99_plot = p99_plot.p99plot('pr',years,observations,mds,ods,'RegCM5')
        #cache.cleanup( )

