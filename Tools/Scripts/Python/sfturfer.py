#!/usr/bin/env python3
"""
sfturfer.py
==============
Create a CORDEX-CMIP6 compliant netCDF4 file containing the variable `sfturf`
(Surface Urban Area Fraction, in %) by reading the urb_2d data from CLM4.5
RegCM input file and writing in a new netCDF file copying all the required
bits from the pycordexed orog fx file.

Algorithm
---------
1. Read the CLM4.5 source file (read-only).
   - Load `urb_2d(numurbl, iy, jx)`.
   - Sum over the `numurbl` dimension → shape (iy, jx).
   - Trim the border: keep [1:-2, 1:-2] → shape (iy-3, jx-3).
     This must match (rlat, rlon) of the target template.

2. Open the CORDEX-CMIP6 target template file (read-only).
   - Copy every dimension, every coordinate/ancillary variable, and
     every global attribute.
   - Replace the `orog` variable with `sfturf` (same shape, metadata as in
     the CORDEX-CMIP6 fx table).
   - Set the global attributes that belong to the new file: variable_id,
     tracking_id, creation_date, history, and external_variables.

3. Write the output file as NETCDF4_CLASSIC (required by CORDEX-CMIP6),
   compressed with deflate level 1 and shuffle.

Usage
-----
    python sfturfer.py <source_clm45.nc> <target_template.nc> <output.nc>

All three paths are positional arguments; the output file must not exist
beforehand (the script will refuse to overwrite it for safety).

Changes (October 2026, Chen Lu), so the output passes the compliance checker
(wcrp_cordex_cmip6:1.0) with no Required flags:
  - NETCDF4_CLASSIC instead of NETCDF4                               [FILE002]
  - variable_id = "sfturf" instead of the copied "orog"     [ATTR005, PATH001]
  - a new tracking_id (the copied one belonged to the orog file) and a new
    creation_date
  - sfturf attributes as in the CORDEX-CMIP6 fx table: standard_name
    "area_fraction", long_name, cell_methods "area: mean", and cell_measures
    "area: areacella" (with the global attribute external_variables =
    "areacella"). The table's scalar coordinate type = "urban" is left out,
    because the checker requires all coordinate variables to be float64
    (a character coordinate gives a VAR011 Required flag).
  - missing_value stored as float, like the data; deflate level 1 + shuffle
"""

import sys
import os
import uuid
import datetime
import numpy as np
import netCDF4 as nc


# Compression recommended by the CORDEX-CMIP6 archive specifications
COMPRESSION = dict(zlib=True, complevel=1, shuffle=True)


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def copy_variable(src_var, dst_ds, extra_attrs=None):
    """
    Copy a variable from an open source Dataset/Group into dst_ds,
    preserving datatype, dimensions, fill value, and all attributes.
    `extra_attrs` (dict) overrides specific attribute values.
    """
    # Collect kwargs for createVariable
    kwargs = dict(COMPRESSION) if src_var.ndim > 0 else {}
    if hasattr(src_var, '_FillValue'):
        kwargs['fill_value'] = src_var._FillValue

    dst_var = dst_ds.createVariable(
        src_var.name,
        src_var.dtype,
        src_var.dimensions,
        **kwargs,
    )

    # Copy attributes (skip _FillValue – already set via fill_value kwarg)
    skip = {'_FillValue'}
    attrs = {k: src_var.getncattr(k)
               for k in src_var.ncattrs() if k not in skip}
    if extra_attrs:
        attrs.update(extra_attrs)
    dst_var.setncatts(attrs)

    # Copy data
    src_var.set_auto_maskandscale(False)
    dst_var.set_auto_maskandscale(False)
    dst_var[...] = src_var[...]
    return dst_var


def new_tracking_id(old):
    """Same prefix as the template (e.g. 'hdl:21.14103/'), new UUID."""
    old = str(old or '')
    prefix = old.rsplit('/', 1)[0] + '/' if '/' in old else 'hdl:21.14103/'
    return prefix + str(uuid.uuid4())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(
            "Usage: python "+sys.argv[0]+" <source_clm45.nc> "
            "<target_template.nc> <output.nc>"
        )

    src_path    = sys.argv[1]
    target_path = sys.argv[2]
    out_path    = sys.argv[3]

    # Refuse to overwrite
    if os.path.exists(out_path):
        sys.exit(f"ERROR: output file already exists: {out_path}\n"
                 "Remove it first or choose a different name.")

    # Read and process the CLM4.5 source
    print(f"Opening source file : {src_path}")
    with nc.Dataset(src_path, 'r') as src:
        if 'urb_2d' not in src.variables:
            sys.exit("ERROR: variable 'urb_2d' not found in source file.")

        urb_2d = src.variables['urb_2d']   # (numurbl, iy, jx)
        raw    = urb_2d[:]                  # masked array

        print(f"  urb_2d shape      : {raw.shape}  (numurbl, iy, jx)")

        # Replace fill values with 0 before summing so they don't
        # contaminate valid urban fractions.
        if hasattr(urb_2d, '_FillValue'):
            fill_val = urb_2d._FillValue
        else:
            fill_val = -9999.0
        data = np.where(np.ma.getdata(raw) == fill_val, 0.0, np.ma.getdata(raw))

        # Sum over numurbl axis (axis=0)
        urb_sum = np.sum(data, axis=0)          # (iy, jx)
        print(f"  after sum         : {urb_sum.shape}  (iy, jx)")

        # Trim border: [1:-2, 1:-2]
        sfturf_data = urb_sum[1:-2, 1:-2]       # (iy-3, jx-3)
        print(f"  after trim [1:-2,1:-2]: {sfturf_data.shape}  (rlat, rlon)")

    # Open target template and build output
    print(f"Opening target file : {target_path}")
    with nc.Dataset(target_path, 'r') as tgt, \
         nc.Dataset(out_path, 'w', format='NETCDF4_CLASSIC') as dst:

        # ---- verify shape agreement ----
        rlat = tgt.dimensions['rlat'].size
        rlon = tgt.dimensions['rlon'].size
        if sfturf_data.shape != (rlat, rlon):
            sys.exit(
                f"ERROR: processed data shape {sfturf_data.shape} does not "
                f"match target grid ({rlat}, {rlon}).\n"
                "Check the trim slice or source/target file pair."
            )

        # ---- dimensions ----
        for dim_name, dim in tgt.dimensions.items():
            size = None if dim.isunlimited() else len(dim)
            dst.createDimension(dim_name, size)
            print(f"  dimension: {dim_name} = {len(dim)}")

        # ---- coordinate and ancillary variables (all except the main one) ----
        skip_vars = {'orog', str(getattr(tgt, 'variable_id', 'orog'))}
        for var_name, var in tgt.variables.items():
            if var_name in skip_vars:
                continue
            copy_variable(var, dst)
            print(f"  copied variable  : {var_name}")

        # ---- sfturf variable (replaces orog) ----
        FILL_VALUE = np.float32(1.e+20)
        sfturf_var = dst.createVariable(
            'sfturf',
            'f4',
            ('rlat', 'rlon'),
            fill_value=FILL_VALUE,
            **COMPRESSION,
        )
        sfturf_var.standard_name = "area_fraction"
        sfturf_var.long_name     = "Percentage of the Grid Cell Occupied by Urban Area"
        sfturf_var.units         = "%"
        sfturf_var.coordinates   = "lat lon"
        sfturf_var.grid_mapping  = "crs"
        sfturf_var.missing_value = FILL_VALUE
        sfturf_var.cell_methods  = "area: mean"
        sfturf_var.cell_measures = "area: areacella"
        sfturf_var.comment       = (
            "Sum of CLM4.5 urb_2d urban class fractions "
            "(TBD=Tall Building District, HD=High Density, MD=Medium Density)"
        )
        sfturf_var[...] = sfturf_data.astype('f4')
        print("  created variable : sfturf")

        # ---- global attributes (copy all from target, then set the new ones) ----
        global_attrs = {k: tgt.getncattr(k) for k in tgt.ncattrs()}

        now = datetime.datetime.now(datetime.timezone.utc)
        global_attrs['variable_id']        = 'sfturf'
        global_attrs['tracking_id']        = new_tracking_id(global_attrs.get('tracking_id'))
        global_attrs['creation_date']      = now.strftime("%Y-%m-%dT%H:%M:%SZ")
        global_attrs['external_variables'] = 'areacella'

        # Update history with today's action
        old_history = global_attrs.get('history', '')
        global_attrs['history'] = (
            now.strftime("%Y-%m-%d %H:%M:%S")
            + " : sfturf created from CLM4.5 urb_2d by "
            + os.path.basename(sys.argv[0]) + "\n" + old_history
        )

        dst.setncatts(global_attrs)
        print("  global attributes : copied; variable_id, tracking_id, "
              "creation_date, external_variables, and history set")

    print(f"\nDone. Output written to: {out_path}")


if __name__ == '__main__':
    main()
