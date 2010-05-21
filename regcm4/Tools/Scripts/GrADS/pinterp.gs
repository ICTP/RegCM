function pinterp(field,pgrid,plev)
*----------------------------------------------------------------------
* Bob Hart (hart@ems.psu.edu) /  PSU Meteorology
* 2/26/1999
*
* 2/26/1999 - Fixed a bug that caused the script to crash on 
*             certain machines.
*
* GrADS function to interpolate within a 3-D grid to a specified
* pressure level.  Can also be used on non-pressure level data, such
* as sigma or eta-coordinate output where pressure is a function
* of time and grid level.
* 
* Advantages:  Easy to use, no UDFs.  Disadvantages:  Can take a few secs.
*
* Arguments:
*    field = name of 3-D grid to interpolate
*
*    pgrid = name of 3-D grid holding pressure values at each gridpoint
*            If you are using regular pressure-level data, this should be
*            set to the builtin GrADS variable 'lev'.
*
*    plev  = pressure level at which to interpolate
*
* Function returns:  defined grid interp holding interpolated values
* 
* NOTE:  YOU NEED TO INCLUDE A COPY OF THIS FUNCTION IN YOUR SCRIPT
*
* Note:  Areas having plev below bottom level or above upper level 
*        will be undefined in output field. Extrapolation is NOT
*        performed!!
*
*---------------------------------------------------------------------
*
* EXAMPLE FUNCTION CALLS:
*
* Sample variables: u = u-wind in m/s
*                   v = v-wind in m/s
*                   t = temperature in K
*                  PP = pressure data in mb
*
* 1) Display a temperature field interpolated to 925mb
*
*      "d "pinterp(t,PP,925)
*
* 2) Display a streamline field interpolated to 550mb
*
*      "define u550="pinterp(u,PP,550)
*      "define v550="pinterp(v,PP,550)
*      "set gxout stream"
*      "d u550;v550"
*
*
* PROBLEMS:  Send email to Bob Hart (hart@ems.psu.edu)
*
*-----------------------------------------------------------------------
*-------------------- BEGINNING OF FUNCTION ----------------------------
*-----------------------------------------------------------------------

* Get initial dimensions of dataset so that exit dimensions will be
* same

"q dims"
rec=sublin(result,4)
ztype=subwrd(rec,3)
if (ztype = "fixed") 
   zmin=subwrd(rec,9)
   zmax=zmin
else
   zmin=subwrd(rec,11)
   zmax=subwrd(rec,13)
endif

* epsilon here will take care of rare instance where user requests 
* an interpolated level that already exists in the file.   Necessary
* to make maskout functions below work; otherwise, the 
* routine finds no set of levels that is valid for interpolation
* Essentially, adding epsilon below 'tips' the
* level selection to one side when two choices are available without
* altering the output of the interpolated results.

epsilon=0.0001

* Get full vertical dimensions of dataset

"q file"
rec=sublin(result,5)
zsize=subwrd(rec,9)
"set z 1 "zsize
"q dims"
rec=sublin(result,4)
pmax=subwrd(rec,6)
pmin=subwrd(rec,8)

* Shift epsilon tip to other side in case user requests lowest level

If (plev = pmax) 
    epsilon=-epsilon
Endif

* Determine spatially varying bounding pressure levels for p-surface
* pabove = pressure-value at level above ; pbelow = pressure value at level
* below for each gridpoint

"set z 1 "zsize-1
"define pabove=0.5*maskout("pgrid","plev+epsilon"-"pgrid")+0.5*maskout("pgrid","pgrid"(z-1)-"plev+epsilon")"

"set z 2 "zsize
"define pbelow=0.5*maskout("pgrid","plev+epsilon"-"pgrid"(z+1))+0.5*maskout("pgrid","pgrid"-"plev+epsilon")"

* Isolate field values at bounding pressure levels
* fabove = requested field value above pressure surface
* fbelow = requested field value below pressure surface

"define fabove=pabove*0+"field
"define fbelow=pbelow*0+"field

* Turn this 3-D grid of values (mostly undefined) into a 2-D press layer.
* mean is used here only for simplicity.  

"set z 1"
"define pabove=mean(pabove,z=1,z="zsize")"
"define fabove=mean(fabove,z=1,z="zsize")"
"define pbelow=mean(pbelow,z=1,z="zsize")"
"define fbelow=mean(fbelow,z=1,z="zsize")"

* Finally, interpolate linearly in log-pressure and create surface.

"set z "zmin " " zmax

"define interp=fbelow+log(pbelow/"plev")*(fabove-fbelow)/log(pbelow/pabove)"

say "Done.  Newly defined variable interp has "plev"mb "field"-field."

return(interp)
