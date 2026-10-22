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
     every global attribute exactly as-is.
   - Replace the `orog` variable with `sfturf` (same shape, new
     metadata appropriate for urban fraction).

3. Write the output file.

Usage
-----
    python sfturfer.py <source_clm45.nc> <target_template.nc> <output.nc>

All three paths are positional arguments; the output file must not exist
beforehand (the script will refuse to overwrite it for safety).
"""

import sys
import os
import datetime
import numpy as np
import netCDF4 as nc


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def copy_variable(src_var, dst_ds, extra_attrs=None):
    """
    Copy a variable from an open source Dataset/Group into dst_ds,
    preserving datatype, dimensions, fill value, chunk hints, and all
    attributes.  `extra_attrs` (dict) overrides specific attribute values.
    """
    # Collect kwargs for createVariable
    kwargs = {}
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
    dst_var[:] = src_var[:]
    return dst_var


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
        else
            fill_val = -9999.0
        data = np.where(raw.data == fill_val, 0.0, raw.data)

        # Sum over numurbl axis (axis=0)
        urb_sum = np.sum(data, axis=0)          # (iy, jx)
        print(f"  after sum         : {urb_sum.shape}  (iy, jx)")

        # Trim border: [1:-2, 1:-2]
        sfturf_data = urb_sum[1:-2, 1:-2]       # (iy-3, jx-3)
        print(f"  after trim [1:-2,1:-2]: {sfturf_data.shape}  (rlat, rlon)")

    # Open target template and build output
    print(f"Opening target file : {target_path}")
    with nc.Dataset(target_path, 'r') as tgt, \
         nc.Dataset(out_path, 'w', format='NETCDF4') as dst:

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

        # ---- coordinate and ancillary variables (all except 'orog') ----
        skip_vars = {'orog'}
        for var_name, var in tgt.variables.items():
            if var_name in skip_vars:
                continue
            copy_variable(var, dst)
            print(f"  copied variable  : {var_name}")

        # ---- sfturf variable (replaces orog) ----
        FILL_VALUE = 1.e+20
        sfturf_var = dst.createVariable(
            'sfturf',
            'f4',
            ('rlat', 'rlon'),
            fill_value=FILL_VALUE,
        )
        sfturf_var.long_name    = "Surface Urban Area Fraction"
        sfturf_var.standard_name = "surface_urban_area_fraction"
        sfturf_var.units        = "%"
        sfturf_var.coordinates  = "lat lon"
        sfturf_var.grid_mapping = "crs"
        sfturf_var.missing_value = FILL_VALUE
        sfturf_var.cell_methods  = "area: sum"
        sfturf_var.comment       = (
            "Sum of CLM4.5 urb_2d urban class fractions "
            "(TBD=Tall Building District, HD=High Density, MD=Medium Density)"
        )
        sfturf_var[:] = sfturf_data.astype('f4')
        print("  created variable : sfturf")

        # ---- global attributes (copy all from target, update history) ----
        global_attrs = {k: tgt.getncattr(k) for k in tgt.ncattrs()}

        # Update history with today's action
        now = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
        old_history = global_attrs.get('history', '')
        global_attrs['history'] = (
            f"{now} : sfturf created from CLM4.5 urb_2d by "+sys.argv[0]+"\n"
            + old_history
        )

        dst.setncatts(global_attrs)
        print("  global attributes : copied and history updated")

    print(f"\nDone. Output written to: {out_path}")


if __name__ == '__main__':
    main()
