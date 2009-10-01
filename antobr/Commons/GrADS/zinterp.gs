function zinterp(field,zgrid,zlev)
*----------------------------------------------------------------------
* 
* Bob Hart (hart@ems.psu.edu) /  PSU Meteorology
* 3/4/1999
*
* GrADS function to interpolate within a 3-D grid to a specified
* height level.  Can also be used on non-pressure level data, such
* as sigma or eta-coordinate output where height is a function
* of time and grid level.
* 
* Advantages:  Easy to use, no UDFs.  Disadvantages:  Can take 3-10 secs.
*
* Arguments:
*    field = name of 3-D grid to interpolate
*
*    zgrid = name of 3-D grid holding height values at each gridpoint
*            Units of measurement are arbitrary.
*
*    zlev  = height level at which to interpolate (having same units as zgrid)
*
* Function returns:  defined grid interp holding interpolated values
*
* NOTE:  Areas having zlev below bottom level or above upper level 
*        in output will be undefined.
*
* NOTE:  No distinction in the function is made between height above
*        sea level, height above ground surface, or geopotential. The
*        function will give output regardless of which is sent.  
*        It is up to the user to be aware which height variable is
*        being passed to the function and treat the output accordingly.
*
* Example function calls:
*
*      "d "zinterp(temp,height,5000)
*
* Would display a temperature field interpolated to 5000.
*      
*      "define t1000="zinterp(temp,height,1000)
*
* Would define a new variable, t1000, as a temp field at 1000.
*
*      "d p1000="zinterp(lev,height,1000)
*
* Would display a field of the pressure at a height of 1000.
*
* PROBLEMS:  Send email to Bob Hart (hart@ems.psu.edu)
* 
*-----------------------------------------------------------------------
*-------------------- BEGINNING OF FUNCTION ----------------------------
*-----------------------------------------------------------------------
*
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

* Get full z-dimensions of dataset.

"q file"
rec=sublin(result,5)
zsize=subwrd(rec,9)

* Determine spatially varying bounding height levels for height surface
* zabove = height-value at level above ; zbelow = height value at level
* below for each gridpoint

"set z 2 "zsize
"define zabove=0.5*maskout("zgrid","zgrid"-"zlev")+0.5*maskout("zgrid","zlev"-"zgrid"(z-1))"
"set z 1 "zsize-1
"define zbelow=0.5*maskout("zgrid","zgrid"(z+1)-"zlev")+0.5*maskout("zgrid","zlev"-"zgrid")"

* Isolate field values at bounding height levels
* fabove = requested field value above height surface
* fbelow = requested field value below height surface

"set z 2 "zsize
"define fabove=zabove*0+"field
"set z 1 "zsize-1
"define fbelow=zbelow*0+"field

* Turn this 3-D grid of values (mostly undefined) into a 2-D height layer
* mean is used here below for simplicity, since mean ignores undefined
* values.

"set z 1"
"define zabove=mean(zabove,z=2,z="zsize")"
"define fabove=mean(fabove,z=2,z="zsize")"
"define zbelow=mean(zbelow,z=1,z="zsize-1")"
"define fbelow=mean(fbelow,z=1,z="zsize-1")"

* Finally, interpolate linearly in height and create surface.

"set z "zmin " " zmax

"define slope=(fabove-fbelow)/(zabove-zbelow)"
"define b=fbelow-slope*zbelow"
"define interp=slope*"zlev"+b"

* variable interp now holds height field and its named it returned
* for use by the user.

say "Done.  Newly defined variable interp has "zlev" "field"-field."

Return(interp)
