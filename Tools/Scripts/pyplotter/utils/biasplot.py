#!/usr/bin/env pyton3

import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors

bar_labels = {"pr": "Precipitation Bias (mm/day)",
              "tas": "Temperature Bias (°C)",
              "clt": "Cloud Cover Bias (%)",
              "rsnl": "Surface Net Upward Longwave Flux Bias (W/m²)"}

    # Colormap settings 
hex_colors = {
    "pr": ["#543005", "#8c510a", "#bf812d", "#d8b365", "#ebd9b2",
           "#ffffff", "#c7eae5", "#80cdc1", "#35978f", "#01665e",
           "#014636"],
    "tas": ["#08306b", "#287bb9", "#5aaecf", "#9cd2e9", "#c1e4f2",
            "#ffffff", "#fddbc7", "#fc9272", "#fc6a54", "#d73027",
            "#7a1420"],
    "clt": "BrBG", #"CBR_drywet",
    "rsnl": ["#08306b", "#287bb9", "#5aaecf", "#9cd2e9", "#c1e4f2",
             "#ffffff", "#fddbc7", "#fc9272", "#fc6a54", "#d73027",
             "#7a1420"]
}
# Plotting levels
var_levels = {
    "pr":   [-8, -6, -4, -2, -0.5, 0.5, 2, 4, 6, 8],
    "tas":  [-8, -6, -4, -2, -1, 1, 2, 4, 6, 8],
    "clt":  [-50, -40, -30, -20, -10, 10, 20, 30, 40, 50],
    "rsnl": [-50, -25, -15, -10, -5, 5, 10, 15, 25, 50]
}

def biasplot(var,years,seas,obs,msd,osd,dpi,pname,dest):

    os.makedirs(dest,exist_ok=True)

    nrows = len(seas)
    ncols = len(obs)

    print(f"Processing variable: {var}")
    
    # Define colormap
    if isinstance(hex_colors[var], list):
        custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom",
                hex_colors[var])
    else:
        # Predefined Matplotlib colormap
        custom_cmap = plt.get_cmap(hex_colors[var])

    # Plotting biases
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
            figsize=(ncols * 3, 11.5),
            subplot_kw={"projection": ccrs.PlateCarree()})
    
    for row, season in enumerate(seas):
        rcm_seasonal = msd[var][row]
        for col, obs_name in enumerate(obs):
            obs_seasonal = osd[col][var][row]
            bias = rcm_seasonal - obs_seasonal
            ax = axes[row] if ncols == 1 else axes[row, col]
            ax.set_extent([rcm_seasonal.xlon.min(), rcm_seasonal.xlon.max(),
                           rcm_seasonal.xlat.min(), rcm_seasonal.xlat.max()],
                           crs=ccrs.PlateCarree())
            ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='grey')
            ax.add_feature(cfeature.COASTLINE, linewidth=0.5) 
 
            levels = var_levels.get(var, np.linspace(-8, 8, 11))  
            cs = ax.contourf(rcm_seasonal.xlon,
                             rcm_seasonal.xlat,
                             bias.squeeze(),
                             levels=levels,
                             cmap=custom_cmap,
                             extend="both")
            
            if row == 0:
                ax.set_title(f"{pname} - {obs_name}", fontsize=10)

            # Add rectangle with season label
            if col == 0:
                ax.text(-14.5, -40.5, f"{season}",
                        fontsize=10, va="center", ha="center", color='grey',
                        bbox=dict(facecolor="white", alpha=0.5,
                            edgecolor="grey", boxstyle="round,pad=0.3"))
    
    plt.subplots_adjust(wspace=0.02, hspace=0.02)
    if var == 'rsnl': #(single column)
        cbar = fig.colorbar(cs, ax=axes,
                orientation="vertical", fraction=0.08, pad=0.04)
    else:
        cbar = fig.colorbar(cs, ax=axes,
                orientation="vertical", fraction=0.04, pad=0.04)
    cbar.set_label(bar_labels[var])

    cyears = "-".join(str(y) for y in years)
    pfile = pname.replace(" ","_")
    out_file = os.path.join(dest, f"{pfile}_bias_{var}_{cyears}_4s.png")
    plt.savefig(out_file, bbox_inches="tight", dpi=int(dpi))
    plt.close(fig)

if __name__ == "__main__":
    print("Nothing to test.")
