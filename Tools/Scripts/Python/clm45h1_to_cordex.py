#!/usr/bin/env python3
"""
clm45_to_cordex.py
==================
Convert CLM4.5 RegCM-coupled output (1-D gridcell dimension) to a set of
CORDEX-CMIP6 compliant netCDF4 files on the full 2-D grid, one output file
per source variable.

Algorithm
---------
1. Read the source CLM4.5 file:
   - Load `regcm_mask(iy, jx)`: the 2-D land mask whose True cells, scanned
     south-to-north and west-to-east (i.e. row-major starting at [0,0] = SW),
     define the order of the 1-D `gridcell` dimension.
   - Build the scatter index arrays (iy_idx, jx_idx) from the mask.
   - Convert the time axis from "hours since 1950-01-01" to
     "days since 1950-01-01" and derive time bounds.
   - Identify every variable with dimensions (time, gridcell) — these are
     the fields to export.

2. Read the destination template:
   - Copy all coordinate/ancillary variables verbatim
   - Use global attributes as the baseline for every output file,
     overriding only the file-specific ones.

3. For each source variable write one output file:
   - Scatter the 1-D data into the 2-D grid; ocean/outside-mask cells
     keep the fill value (1.e+20f).
   - Variable name is the source name lowercased.
   - Filename follows the CORDEX-CMIP6 convention (see build_filename()).
   - Variable attributes and global attributes are adjusted to reflect
     the CLM4.5 land-component origin and the urban/rural surface class.

Usage
-----
    python clm45_to_cordex.py <source_clm45.nc> <destination_template.nc> \\
                              [--outdir <output_directory>]

Arguments
---------
    source_clm45.nc          CLM4.5 RegCM-coupled daily output file.
    destination_template.nc  CORDEX-CMIP6 template (e.g. a tas file).
    --outdir                 Directory for output files (default: current dir).

Dependencies
------------
    numpy, netCDF4  (both standard in any scientific Python stack)
"""

import sys
import os
import uuid
import argparse
import datetime
import numpy as np
import netCDF4 as nc

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

FILL_VALUE = 1.e+20          # canonical CORDEX fill value
EPOCH_DAYS = "days since 1950-01-01 00:00:00"
CLM_NOTE   = (
    "This variable is a diagnostic from the CLM4.5 land-surface component "
    "of the RegCM coupled model. It is NOT part of the official CORDEX-CMIP6 "
    "data request. Suffix _R denotes the rural surface class; "
    "suffix _U denotes the urban surface class."
)

# Metadata overrides per variable name (lowercase).
# Keys: long_name, standard_name (if CF-compliant one exists), units,
#       cell_methods, comment.
# Anything not listed here falls back to what is in the source file.
VAR_META = {
    # latent heat
    "eflxlhtotr": dict(standard_name="surface_upward_latent_heat_flux",
                          cell_methods="area: time: mean"),
    "eflxlhtotu": dict(standard_name="surface_upward_latent_heat_flux",
                          cell_methods="area: time: mean"),
    # ground heat flux
    "fgrr": dict(standard_name="downward_heat_flux_at_ground_level_in_soil",
                  cell_methods="area: time: mean"),
    "fgru": dict(standard_name="downward_heat_flux_at_ground_level_in_soil",
                  cell_methods="area: time: mean"),
    # net LW radiation
    "firar": dict(standard_name="surface_net_downward_longwave_flux",
                   cell_methods="area: time: mean"),
    "firau": dict(standard_name="surface_net_downward_longwave_flux",
                   cell_methods="area: time: mean"),
    # absorbed SW radiation
    "fsar": dict(standard_name="surface_net_downward_shortwave_flux",
                  cell_methods="area: time: mean"),
    "fsau": dict(standard_name="surface_net_downward_shortwave_flux",
                  cell_methods="area: time: mean"),
    # sensible heat
    "fshr": dict(standard_name="surface_upward_sensible_heat_flux",
                  cell_methods="area: time: mean"),
    "fshu": dict(standard_name="surface_upward_sensible_heat_flux",
                  cell_methods="area: time: mean"),
    # snow melt heat flux — no single CF name; keep source long_name
    "fsmr": dict(cell_methods="area: time: mean"),
    "fsmu": dict(cell_methods="area: time: mean"),
    # runoff
    "qrunoffr": dict(standard_name="surface_runoff_flux",
                      cell_methods="area: time: mean"),
    "qrunoffu": dict(standard_name="surface_runoff_flux",
                      cell_methods="area: time: mean"),
    # relative humidity
    "rh2mr": dict(standard_name="relative_humidity",
                   cell_methods="area: time: mean"),
    "rh2mu": dict(standard_name="relative_humidity",
                   cell_methods="area: time: mean"),
    # soil alpha — no CF name
    "soilalphau": dict(cell_methods="area: time: mean"),
    # ground temperature
    "tgr": dict(standard_name="soil_temperature",
                 cell_methods="area: time: mean"),
    "tgu": dict(standard_name="soil_temperature",
                 cell_methods="area: time: mean"),
    # daily min T2m
    "trefmnavr": dict(standard_name="air_temperature",
                       cell_methods="area: time: minimum within days time: mean over days"),
    "trefmnavu": dict(standard_name="air_temperature",
                       cell_methods="area: time: minimum within days time: mean over days"),
    # daily max T2m
    "trefmxavr": dict(standard_name="air_temperature",
                       cell_methods="area: time: maximum within days time: mean over days"),
    "trefmxavu": dict(standard_name="air_temperature",
                       cell_methods="area: time: maximum within days time: mean over days"),
    # mean T2m
    "tsar": dict(standard_name="air_temperature",
                  cell_methods="area: time: mean"),
    "tsau": dict(standard_name="air_temperature",
                  cell_methods="area: time: mean"),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("source",      help="CLM4.5 RegCM-coupled source netCDF file")
    p.add_argument("destination", help="CORDEX-CMIP6 template netCDF file")
    p.add_argument("--outdir", default=".", help="Output directory (default: .)")
    return p.parse_args()


def hours_to_days(hours_values, time_units_hours):
    """
    Convert a time array in 'hours since <epoch>' to
    'days since 1950-01-01 00:00:00'.
    Returns the day-valued array and the matching bounds array.
    """
    import netCDF4 as nc
    # Decode to datetime objects
    dts = nc.num2date(hours_values, time_units_hours,
                      calendar="proleptic_gregorian")
    # Re-encode as days since 1950-01-01
    days = nc.date2num(dts, EPOCH_DAYS, calendar="proleptic_gregorian")
    return np.array(days, dtype=np.float64)


def build_mask_index(regcm_mask):
    """
    Given regcm_mask(iy, jx) with 1=land, build two 1-D arrays
    iy_idx and jx_idx of length == number of True cells, in the
    south-to-north, west-to-east scan order (row-major from [0,0]).

    Returns
    -------
    iy_idx, jx_idx : np.ndarray of int, shape (n_land,)
    """
    mask_bool = np.array(regcm_mask[:], dtype=bool)   # (iy, jx)
    # np.argwhere returns row-major order -> south-to-north, west-to-east
    indices = np.argwhere(mask_bool)                   # (n_land, 2)
    iy_idx  = indices[:, 0]
    jx_idx  = indices[:, 1]
    return iy_idx, jx_idx


def scatter_to_2d(data_1d, iy_idx, jx_idx, ny, nx):
    """
    Scatter a 2-D array (time, gridcell) to (time, ny, nx).
    Non-land cells are set to FILL_VALUE.
    """
    nt = data_1d.shape[0]
    out = np.full((nt, ny, nx), FILL_VALUE, dtype=np.float32)
    out[:, iy_idx, jx_idx] = data_1d
    return out


def build_filename(template_path, var_name_lower, date_start, date_end, outdir):
    """
    Build a CORDEX-CMIP6 compliant filename from the template filename,
    substituting the variable name and date range.

    Template example:
      tas_MED-12_ERA5_evaluation_r1i1p1f1_ICTP_RegCM-ES1-1_v1-r1_day_19800101-19800131.nc

    Pattern (underscore-separated tokens):
      <var>_<domain>_<driver>_<experiment>_<ensemble>_<institution>_<model>_<version>_<freq>_<daterange>.nc
    """
    basename = os.path.basename(template_path)
    # Strip .nc extension
    if basename.endswith(".nc"):
        basename = basename[:-3]
    tokens = basename.split("_")
    # Replace first token (variable name)
    tokens[0] = var_name_lower
    # Replace last token (date range)
    d0 = date_start.strftime("%Y%m%d")
    d1 = date_end.strftime("%Y%m%d")
    tokens[-1] = f"{d0}-{d1}"
    fname = "_".join(tokens) + ".nc"
    opathd = os.path.join(outdir,var_name_lower)
    os.makedirs(opathd,exist_ok=True)
    return os.path.join(opathd, fname)


def copy_coord_var(src_var, dst_ds):
    """Copy a coordinate/ancillary variable verbatim."""
    kwargs = {}
    if hasattr(src_var, '_FillValue'):
        kwargs['fill_value'] = src_var._FillValue
    dst_var = dst_ds.createVariable(
        src_var.name, src_var.dtype, src_var.dimensions, **kwargs)
    skip = {'_FillValue'}
    dst_var.setncatts(
        {k: src_var.getncattr(k) for k in src_var.ncattrs() if k not in skip})
    dst_var[:] = src_var[:]
    return dst_var


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Open source and destination simultaneously (read-only)
    # ------------------------------------------------------------------
    print(f"Opening source      : {args.source}")
    print(f"Opening destination : {args.destination}")

    with nc.Dataset(args.source, 'r') as src, \
         nc.Dataset(args.destination, 'r') as dst_tmpl:

        # ---------- grid dimensions from destination ----------
        try:
            rlat_size = dst_tmpl.dimensions['rlat'].size
            rlon_size = dst_tmpl.dimensions['rlon'].size
            rotllr = True
        except:
            rlat_size = dst_tmpl.dimensions['iy'].size
            rlon_size = dst_tmpl.dimensions['jx'].size
            rotllr = False
        try:
            nv_size   = dst_tmpl.dimensions['nv'].size
        except:
            nv_size = 0

        # ---------- mask & scatter indices ----------
        if 'regcm_mask' not in src.variables:
            sys.exit("ERROR: 'regcm_mask' not found in source file.")
        regcm_mask = src.variables['regcm_mask']
        iy_idx, jx_idx = build_mask_index(regcm_mask)
        n_land = len(iy_idx)
        print(f"  Land cells found  : {n_land}  (gridcell dimension = "
              f"{src.dimensions['gridcell'].size})")
        if n_land != src.dimensions['gridcell'].size:
            sys.exit(
                f"ERROR: mask land cells ({n_land}) != "
                f"gridcell dimension ({src.dimensions['gridcell'].size})."
            )

        # ---------- time conversion ----------
        src_time_var  = src.variables['time']
        src_time_vals = src_time_var[:]
        src_time_unit = src_time_var.units            # hours since 1950-01-01
        # Registered time is the end of the average interval. Need to shift.
        dst_time_vals = hours_to_days(src_time_vals, src_time_unit) - 0.5

        # Time bounds: shift by ±0.5 day around each daily value
        dst_time_bnds = np.column_stack(
            [dst_time_vals - 0.5, dst_time_vals + 0.5])

        # Date objects for filename
        dts = nc.num2date(src_time_vals, src_time_unit,
                          calendar="proleptic_gregorian")
        date_start = datetime.date(dts[0].year,  dts[0].month,  dts[0].day)
        date_end   = datetime.date(dts[-1].year, dts[-1].month, dts[-1].day)

        # ---------- identify target variables ----------
        target_vars = {}
        for vname, var in src.variables.items():
            if var.dimensions == ('time', 'gridcell'):
                target_vars[vname] = var
        print(f"  Variables to export: {sorted(target_vars.keys())}")

        # ---------- baseline global attributes from destination ----------
        base_global = {k: dst_tmpl.getncattr(k) for k in dst_tmpl.ncattrs()}

        # ---------- destination template coordinate variables ----------
        coord_vars_to_copy = [
            'lat', 'lon', 'lat_bnds', 'lon_bnds', 'rlat', 'rlon',
            'iy', 'jx', 'crs'
        ]

        # ------------------------------------------------------------------
        # 2. Write one file per variable
        # ------------------------------------------------------------------
        now_str = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")

        for src_vname, src_var in sorted(target_vars.items()):
            var_lower = src_vname.lower().replace('_','')
            outpath   = build_filename(
                args.destination, var_lower, date_start, date_end, args.outdir)

            if os.path.exists(outpath):
                print(f"  SKIP (exists): {outpath}")
                continue

            print(f"  Writing : {os.path.basename(outpath)}")

            # Load 1-D data and scatter to 2-D
            raw_1d = src_var[:]                            # masked array (time, gridcell)
            # Replace any masked/fill values with FILL_VALUE before scatter
            data_1d = np.where(
                raw_1d.mask if np.ma.is_masked(raw_1d) else False,
                FILL_VALUE, np.array(raw_1d, dtype=np.float32))
            data_2d = scatter_to_2d(data_1d, iy_idx, jx_idx,
                                    rlat_size, rlon_size)   # (time, lat, lon)

            with nc.Dataset(outpath, 'w', format='NETCDF4') as out:

                # ---- dimensions ----
                if rotllr:
                    out.createDimension('rlat', rlat_size)
                    out.createDimension('rlon', rlon_size)
                else:
                    out.createDimension('iy', rlat_size)
                    out.createDimension('jx', rlon_size)
                if nv_size > 0:
                    out.createDimension('nv',   nv_size)
                out.createDimension('time', None)   # unlimited
                out.createDimension('bnds', 2)

                # ---- coordinate variables from destination template ----
                for cv in coord_vars_to_copy:
                    if cv in dst_tmpl.variables:
                        copy_coord_var(dst_tmpl.variables[cv], out)

                # ---- time ----
                tv = out.createVariable('time', 'f8', ('time',))
                tv.long_name    = "time"
                tv.standard_name = "time"
                tv.units        = EPOCH_DAYS
                tv.calendar     = "proleptic_gregorian"
                tv.leap_seconds = "none"
                tv.bounds       = "time_bnds"
                tv[:]           = dst_time_vals

                tb = out.createVariable('time_bnds', 'f8', ('time', 'bnds'))
                tb[:]           = dst_time_bnds

                # ---- data variable ----
                meta = VAR_META.get(var_lower, {})

                if rotllr:
                    out_var = out.createVariable(
                        var_lower, 'f4', ('time', 'rlat', 'rlon'),
                        fill_value=FILL_VALUE,
                        compression='zlib',
                        complevel=6,
                        significant_digits=4)
                else:
                    out_var = out.createVariable(
                        var_lower, 'f4', ('time', 'iy', 'jx'),
                        fill_value=FILL_VALUE,
                        compression='zlib',
                        complevel=6,
                        significant_digits=4)

                # Start from source attributes, then apply overrides
                src_attrs = {k: src_var.getncattr(k)
                             for k in src_var.ncattrs()
                             if k != '_FillValue'}
                # Rename cell_method -> cell_methods (CF >= 1.3)
                if 'cell_method' in src_attrs:
                    src_attrs['cell_methods'] = src_attrs.pop('cell_method')

                out_var.long_name    = src_attrs.get('long_name', src_vname)
                out_var.units        = src_attrs.get('units', '1')
                out_var.cell_methods = src_attrs.get('cell_methods',
                                                     'area: time: mean')
                out_var.coordinates  = "lat lon"
                out_var.grid_mapping = "crs"
                out_var.missing_value = np.float32(FILL_VALUE)
                out_var.cell_measures = "area: areacella"

                # Apply VAR_META overrides
                if 'standard_name' in meta:
                    out_var.standard_name = meta['standard_name']
                if 'cell_methods' in meta:
                    out_var.cell_methods = meta['cell_methods']

                # Provenance and class annotation
                suffix = src_vname[-2:]   # '_R' or '_U'
                if suffix == '_R':
                    surface_class = "rural"
                elif suffix == '_U':
                    surface_class = "urban"
                else:
                    surface_class = "unspecified"

                out_var.surface_class = surface_class
                out_var.source_variable = src_vname
                out_var.source_component = "CLM4.5 land-surface model"
                out_var.comment = (
                    f"Diagnostic from the CLM4.5 land-surface component of "
                    f"the RegCM coupled model ({surface_class} surface class). "
                    f"Source variable: {src_vname}. "
                    f"NOT part of the official CORDEX-CMIP6 data request."
                )

                out_var[:] = data_2d

                # ---- global attributes ----
                gattrs = dict(base_global)  # copy baseline from template

                # Update / add fields specific to this file
                gattrs['history'] = (
                    f"{now_str} : Converted from CLM4.5 source by "
                    f"clm45_to_cordex.py (source variable: {src_vname})\n"
                    + gattrs.get('history', '')
                )
                gattrs['variable_id']  = var_lower
                gattrs['table_id']     = "day"
                gattrs['frequency']    = "day"
                gattrs['clm45_note']   = CLM_NOTE
                gattrs['source']       = (
                    gattrs.get('source', '') +
                    " CLM4.5 land-surface component output."
                )
                gattrs['Conventions']  = "CF-1.11"
                # Explicitly flag as non-standard CORDEX variable
                gattrs['non_standard_variable'] = (
                    f"{var_lower}: CLM4.5 land-component diagnostic, "
                    f"{surface_class} surface class, not in CORDEX-CMIP6 data request"
                )
                gattrs['tracking_id'] = 'hdl:21.14103/'+str(uuid.uuid4( ))

                out.setncatts(gattrs)

        print(f"\nDone. Output files written to: {os.path.abspath(args.outdir)}")


if __name__ == '__main__':
    main()
