#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors

levels = [-50, -40, -30, -20, -10, -5, 5, 10, 20, 30, 40, 50]
#hex_colors = "BrBG", #"CBR_drywet
hex_colors = ["#543005", "#8c510a", "#bf812d", "#d8b365", "#ebd9b2",
              "#ffffff", "#c7eae5", "#80cdc1", "#35978f", "#01665e",
              "#014636"]

def p99plot(var,years,obs,msd,osd,dpi,pname,dest):
    os.makedirs(dest,exist_ok=True)

    ncols = len(obs)
    print(f"P99 : Processing variable: {var}")

    fig, axes = plt.subplots(nrows=1, ncols=ncols,
            figsize=(ncols * 3, 11.5),
            subplot_kw={"projection": ccrs.PlateCarree()})
    custom_cmap = mcolors.LinearSegmentedColormap.from_list(
            "custom", hex_colors)
    # Predefined Matplotlib colormap
    #custom_cma = plt.get_cmap(hex_colors[var])

    for col, obs_name in enumerate(obs):
        bias = msd[var] - osd[col][var]
        ax = axes[col]
        cs = ax.contourf(msd.xlon.values, msd.xlat.values, bias.squeeze(),
                levels=levels, cmap=custom_cmap, extend="both")
        ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor="grey")
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax.set_title(f"P99 Bias ({pname} - {obs_name})", fontsize=10)

    cbar = fig.colorbar(cs, ax=ax, orientation="vertical", pad=0.05, shrink=0.8)
    cbar.set_label("P99 Precipitation Bias (mm/day)")

    cyears = "-".join(str(y) for y in years)
    pfile = pname.replace(" ","_")
    out_file = os.path.join(dest, f"{pfile}_p99_{var}_{cyears}.png")
    plt.savefig(out_file, bbox_inches="tight", dpi=int(dpi))
    plt.close(fig)

if __name__ == "__main__":
    print("Nothing to test.")
