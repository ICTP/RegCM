function isen(field,tgrid,pgrid,tlev)
*----------------------------------------------------------------------
* Bob Hart (hart@ems.psu.edu) /  PSU Meteorology
* 2/26/1999
*
* 2/26/99 - Fixed a bug that caused the script to crash on
*           certain machines.  
*
* GrADS function to interpolate within a 3-D grid to a specified
* isentropic level.  Can also be used on non-pressure level data, such
* as sigma or eta-coordinate output where pressure is a function
* of time and grid level.  Can be used to create isentropic PV surfaces
* (examples are given at end of documentation just prior to
* function.)
* 
* Advantages:  Easy to use, no UDFs.  Disadvantages:  Can take 5-20 secs.
*
* Arguments:
*    field = name of 3-D grid to interpolate
*
*    tgrid = name of 3-D grid holding temperature values (deg K) at each
*            gridpoint.
*
*    pgrid = name of 3-D grid holding pressure values (mb) at each gridpoint
*            If you are using regular pressure-level data, this should be
*            set to the builtin GrADS variable 'lev'.
*
*    tlev  = theta-level (deg K) at which to interpolate
*
* Function returns:  defined grid interp holding interpolated values
*
* NOTE:  YOU NEED TO INCLUDE A COPY OF THIS FUNCTION IN YOUR SCRIPT
*
* NOTE:  Areas having tlev below bottom level or above upper level 
*        will be undefined in output field. Extrapolation is NOT
*        performed!!
*
*------------------------------------------------------------------------
*
* EXAMPLE FUNCTION CALLS:
*
* Sample variables: u = u-wind in m/s
*                   v = v-wind in m/s
*                   w = vertical velocity
*                   t = temperature in K
*                  PP = pressure data in mb
*
* 1) Display vertical velocity field on 320K surface:
* 
*    "d "isen(w,t,PP,320)
*
* 2) Create & Display colorized streamlines on 320K surface:
*
*    "define u320="isen(u,t,PP,320)
*    "define v320="isen(v,t,PP,320)
*    "set z 1"
*    "set gxout stream"
*    "d u320;v320;mag(u320,v320)"
*
* 3) Create & display a 320K isentropic PV surface:
*
*    "set lev 1050 150"
*    "define coriol=2*7.29e-5*sin(lat*3.1415/180)"
*    "define dudy=cdiff(u,y)/(111177*cdiff(lat,y))"
*    "define dvdx=cdiff(v,x)/(111177*cdiff(lon,x)*cos(lat*3.1415/180))"
*    "define dt=t(z-1)*pow(1000/PP(z-1),0.286)-t(z+1)*pow(1000/PP(z+1),0.286)"
*    "define dp=100*(PP(z-1)-PP(z+1))"
*    "define dtdp=dt/dp"
*    "define part1="isen(dvdx,t,PP,320)
*    "define part2="isen(dudy,t,PP,320)
*    "define part3="isen(dtdp,t,PP,320)
*    "define pv320=-9.8*(coriol+part1-part2)*part3"
*    "set z 1"
*    "d pv320"
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

* Get full z-dimensions of dataset.

"q file"
rec=sublin(result,5)
zsize=subwrd(rec,9)

* Determine spatially varying bounding pressure levels for isen surface
* tabove = theta-value at level above ; tbelow = theta value at level
* below for each gridpoint

"set z 1 "zsize
"define theta="tgrid"*pow(1000/"pgrid",0.286)"
"set z 2 "zsize
"define thetam="tgrid"(z-1)*pow(1000/"pgrid"(z-1),0.286)"
"set z 1 "zsize-1
"define thetap="tgrid"(z+1)*pow(1000/"pgrid"(z+1),0.286)"

"define tabove=0.5*maskout(theta,theta-"tlev")+0.5*maskout(theta,"tlev"-thetam)"
"define tbelow=0.5*maskout(theta,thetap-"tlev")+0.5*maskout(theta,"tlev"-theta)"

* Isolate field values at bounding pressure levels
* fabove = requested field value above isen surface
* fbelow = requested field value below isen surface

"define fabove=tabove*0+"field
"define fbelow=tbelow*0+"field

"set z 1"

* Turn this 3-D grid of values (mostly undefined) into a 2-D isen layer

* If more than one layer is valid (rare), take the mean of all the
* valid levels. Not the best way to deal with the multi-layer issue,
* but works well, rarely if ever impacts output, and is quick.
* Ideally, only the upper most level would be used.  However, this
* is not easily done using current GrADS intrinsic functions.

"define fabove=mean(fabove,z=1,z="zsize")"
"define fbelow=mean(fbelow,z=1,z="zsize")"
"define tabove=mean(tabove,z=1,z="zsize")"
"define tbelow=mean(tbelow,z=1,z="zsize")"

* Finally, interpolate linearly in theta and create isen surface.
* Linear interpolation in theta works b/c it scales as height,
* or log-P, from Poisson equation for pot temp.

"set z "zmin " " zmax

"define slope=(fabove-fbelow)/(tabove-tbelow)"
"define b=fbelow-slope*tbelow"
"define interp=slope*"tlev"+b"

* variable interp now holds isentropic field and its named it returned
* for use by the user.

say "Done.  Newly defined variable interp has "tlev"K "field"-field."

return(interp)
