* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* meteogram_eta.gs
*
* Customized for ETA meteograms on COLA's web site --JMA May 2003
* Data taken from GDS at http://cola8.iges.org:9191
*
* This script draws a meteogram based on NCEP forecast data.
* The data is available through the GrADS-DODS Server at COLA.
* You MUST be using a DODS-enabled version of GrADS.
*
* Usage:   meteogram_eta <Name> <yyyymmddhh> <lon> <lat> <e>
* Example: meteogram_eta Boston  2003022100   -71    42   e 
*
* The 'e' argument is for British units. Default is metric.
*
* Note: This script must be run in a directory in which 
* you have write permission because intermediate files 
* are written out to disk in order to speed up the display
* and minimize the number of hits to the data server. 
*  
* Originally written by Paul Dirmeyer
* Modification History:
* J.M. Adams   Oct 2001
* Jim Kinter   Oct 2001
* J.M. Adams   Dec 2001
* Joe Wielgosz Jan 2002
* J.M. Adams   Jul 2002
* J.M. Adams   May 2003
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
function main(args)

* Make sure GrADS is in portrait mode
'q gxinfo'
pagesize = sublin(result,2)
xsize = subwrd(pagesize,4)
ysize = subwrd(pagesize,6)
if (xsize != 8.5 & ysize != 11)
  say 'You must be in PORTRAIT MODE to draw a meteogram'
  return
endif

* Parse the arguments: date, longitude, latitude, units
if (args = '')
  prompt 'Enter forecast date (yyyymmddhh) --> ' 
  pull date
  prompt 'Enter longitude --> '
  pull hilon
  prompt 'Enter latitude --> '
  pull hilat
  prompt 'Metric units? [y/n] --> '
  pull metric
  if (metric='n' | metric='N') ; _units='e' ; endif
  prompt 'Enter location Name --> '
  pull name
else 
  name  = subwrd(args,1)
  date  = subwrd(args,2)
  hilon = subwrd(args,3)
  hilat = subwrd(args,4)
 _units = subwrd(args,5)
endif


* Open the data file
'reinit'
_baseurl = 'http://cola8.iges.org:9191/dods/eta/'
_dataset = 'eta'date
'sdfopen '_baseurl%_dataset
check1 = sublin(result,2)
check2 = subwrd(check1,1)
if (check2 != 'Found') ; return ; endif

* Get info from the descriptor file
'q ctlinfo'
_ctl = result
_undef = getctl(undef)
_tdef = getctl(tdef)
_zdef = getctl(zdef)

* Get the Time axis info
tsize = subwrd(_tdef,2)
_t1 = 1
_t2 = tsize
'set t '_t1' '_t2
'q dims'
times  = sublin(result,5)
_time1 = subwrd(times,6)  
_time2 = subwrd(times,8)
_tdim = _time1' '_time2
tincr = subwrd(_tdef,5)
_tdef = 'tdef 'tsize' linear '_time1' 'tincr

* Get Vertical grid info 
zsize = subwrd(_zdef,2)
z = 1
zlevs = ''
rhzlevs = ''
while (z <= zsize)
  'set z 'z
  lev = subwrd(result,4)
  if lev = 500 ; z500 = z ; endif
  zlevs = zlevs%lev%' '
* ETA relative humidity grid values are every 50 mb
  rem = math_fmod(lev,50)
  if (rem = 0)
   rhzlevs = rhzlevs%lev%' '  
  endif
  z = z + 1
endwhile

* Find the grid point closest to requsted location
'set lon 'hilon
hilon = subwrd(result,4)
'set lat 'hilat
hilat = subwrd(result,4)
_xdim = hilon' 'hilon 
_ydim = hilat' 'hilat

* Determine pressure range for hovmoellers
getseries(ps,pshov,1000)
'set lon 'hilon
'set lat 'hilat
'd ave(pshov,t='_t1',t='_t2')*0.01-15.0'

data = sublin(result,2)
mmm = subwrd(data,4)
meanps = math_nint(mmm)
cnt = 1
while (cnt<zsize)
  el1 = subwrd(zlevs,cnt)
  el2 = subwrd(zlevs,cnt+1)
  if (meanps > el2)
    elb = el1
    elt = subwrd(zlevs,z500+cnt-1)
    break
  endif
  cnt=cnt+1
endwhile
if (elt < 500) ; elt = 500 ; endif
_zbot = elb
_ztop = elt
_zgrd = _zbot' '_ztop

* Get pressure range for ETA relative humidity grid
rem = math_fmod(_zbot,50)
if (rem != 0)
  _rhzbot = _zbot - 25
else 
  _rhzbot = _zbot
endif
rem = math_fmod(_ztop,50)
if (rem != 0)
  _rhztop = _ztop + 25
else 
  _rhztop = _ztop
endif
_rhzgrd = _rhzbot' '_rhztop

* Get number of pressure levels between _zbot and _ztop
n = 1
p = subwrd(zlevs,n)
while (p != _zbot)
  n = n + 1
  p = subwrd(zlevs,n)
endwhile
levs = p
nlevs = 1
while (p > _ztop)
  n = n + 1
  p = subwrd(zlevs,n)
  levs = levs%' 'p
  nlevs = nlevs + 1
endwhile
_newzsize = nlevs
_zdef = 'zdef '_newzsize' levels 'levs

* Get TOTAL number of ETA pressure levels between _rhzbot and _rhztop
* including those with no rh value
n = 1
p = subwrd(zlevs,n)
while (p != _rhzbot)
  n = n + 1
  p = subwrd(zlevs,n)
endwhile
tlevs = p
nlevs = 1
while (p > _rhztop)
  n = n + 1
  p = subwrd(zlevs,n)
  tlevs = tlevs%' 'p
  nlevs = nlevs + 1
endwhile
_trhzsize = nlevs
_trhzdef = 'zdef '_trhzsize' levels 'tlevs

* Get number of pressure levels for ETA relative humidity grid
* that contain rh data
n = 1
p = subwrd(rhzlevs,n)
while (p != _rhzbot)
  n = n + 1
  p = subwrd(rhzlevs,n)
endwhile
_rhlevs = p
nlevs = 1
while (p > _rhztop)
  n = n + 1
  p = subwrd(rhzlevs,n)
  _rhlevs = _rhlevs%' 'p
  nlevs = nlevs + 1
endwhile
_newrhzsize = nlevs
_rhzdef = 'zdef '_newrhzsize' levels '_rhlevs

* Get the Z,T grids for the upper air variables
getgrid(t,t)
getgrid(u,u)
getgrid(v,v) 
getgrid('z',hgt) 
getetarh(rh,rhum)

* Set up a few preliminary characteristics
setcols(1)
'set display color white'
'c'

* Determine the plot areas for each panel
* Panels: Xsect, thickness, stability, slp, u10, t2m, rh2m, precip
npanels = 8  
x1 =  1.20
x2 =  8.15
y1 =  7.50
y2 = 10.10
panel.npanels = x1' 'x2' 'y1' 'y2   ;* hovmoeller panel
ytop = 7.5  ;* y boundaries for rest of panels except precip
ybot = 1.7
int = (ytop-ybot)/(npanels-2)     ;* get height of middle panels
int = 0.001 * (math_nint(1000*int))
n=npanels-1
y2 = ytop
while (n >= 2)
  y2 = ytop - (npanels-n-1)*int
  y1 = ytop - (npanels-n)*int
  panel.n = x1' 'x2' 'y1' 'y2        ;* coords of middle panels
  n = n - 1
endwhile
xincr = (8.15 - 1.2)/tsize           ;* size of one time step
xincr = 0.01 * math_nint(100*xincr)
panel.1 = x1+xincr' 'x2' 0.55 'y1     ;* coords of precip panel

* Set the Plot Area for the Upper Air Panel
p = npanels
'set parea 'panel.p
'set vpage off'
'set grads off'
'set grid on'

* Draw the Relative Humidity Shading
'set gxout shaded'
'set csmooth on'
'set clevs  30 50 70 90 100'
'set ccols 0 20 21 23 25 26'
'set xlopts 1 4 0.16'
'set xlpos 0 t'
*'set ylab `1%g'
'set ylab %g'
'set ylint 100'
if (_units = 'e')
  'define temp = (t-273.16)*1.8+32'
else
  'define temp = (t-273.16)'
endif
'set t '_t1-0.5' '_t2+0.5
'set lev '_zbot+50' '_ztop-50
'd rhum'
'set gxout contour'
'set grid off'
'set ccolor 15'
'set clab off'
'set clevs 10 30 50 70 90'
'd rhum'
'set ccolor 0'
'set clab on'
'set cstyle 5'
'set clopts 15'
'set clevs 10 30 50 70 90'
'd rhum'

* Draw the Temperature Contours
'set clopts -1'
'set cstyle 1'
'set ccolor rainbow'
'set rbcols 9 14 4 11 5 13 12 8 2 6'
if (_units = 'e')
  'set cint 10'
  'set cthick 6'
  'd temp'
  'set clevs 32'
  'set cthick 12'
  'set ccolor 1'
  'set clab off'
  'd temp'
  'set background 1'
  'set ccolor 20'
  'set clevs 32'
  'set cthick 4'
  'set clab on'
  'set clab `4FR'
else
  'set cint 5'
  'set cthick 6'
  'd temp'
  'set clevs 0'
  'set cthick 12'
  'set ccolor 1'
  'set clab off'
  'd temp'
  'set background 1'
  'set ccolor 20'
  'set clevs 0'
  'set cthick 4'
  'set clab on'
endif
'd temp'

* Draw the Wind Barbs
'set background 0'
'set gxout barb'
'set ccolor 1'
'set xlab off'
'set ylab off'
*'d uwnd;vwnd'
'd urh;vrh'

* Draw a rectangle over the year to clear the area for a title
'set line 0'
'draw recf 0.5 10.6 2.1 11.0'

* Define Thickness
'set lev 1000'
'set t '_t1' '_t2
*getseries('z',hgt,1000)
getseries('z',z5,500)
getseries('z',z10,1000)
'define thk1 = (z5-z10)/10'

* Next Panel: 1000-500 thickness
p = p - 1
'set parea 'panel.p
'set gxout line'
'set vpage off'
'set grads off'
'set grid on'
'set xlab on'
'set ylab on'
vrng(thk1, thk1)
'set ccolor 5'
'set cmark 4'
'set t '_t1-0.5' '_t2+0.5
'd thk1'

* Draw a rectangle over the x-axis labels
xlo = subwrd(panel.p,1)
xhi = subwrd(panel.p,2)
ylo = subwrd(panel.p,3)
yhi = subwrd(panel.p,4)
'set line 0'
'draw recf 'xlo-0.4' 'ylo-0.8' 'xhi+0.4' 'ylo-0.02

* Next Panel: Stability Indices
p = p - 1
'set parea 'panel.p
ylo = subwrd(panel.p,3)
yhi = subwrd(panel.p,4)
ymid = ylo + (yhi-ylo)/2

'set t '_t1' '_t2
'set lev 1000'
'define rh8 = rhum(lev=850)'
'define t8 = t(lev=850)'
'define t5 = t(lev=500)'
'set vpage off'
'set grads off'
'set grid on'
'set xlab on'
'set gxout bar'
'set barbase 40'
'set bargap 50'
'define toto = 1/(1/t8-log(rh8*0.01)*461/2.5e6)-t5*2+t8'
'set axlim 11 69'
'set yaxis 11 69 10'
'set ccolor 8'
'set t '_t1-0.5' '_t2+0.5
'set grid on'
'd (toto-40+abs(toto-40))*0.5+40'
'set grid off'
'set ccolor 7'
'd (toto-40-abs(toto-40))*0.5+40'

* draw a rectangle over 'toto' yaxis labels
'set line 0'
'draw recf 0.2 'ylo' 1.175 'yhi-0.07

* Lifted Index
getseries(li,li,1000)
'set gxout line'
'set grid off'
'set vrange 5.9 -5.9'
'set yaxis 5.9 -5.9 2'
'set ccolor 2'
'set cstyle 3'
'set cmark 7'
'set cmax 0'
'set datawarn off'
'd li'

* draw a zero line
'set ccolor 15'
'set cmark 0'
'set cstyle 3'
'd const(li,0)'

* Draw a rectangle over the x-axis labels
xlo = subwrd(panel.p,1)
xhi = subwrd(panel.p,2)
ylo = subwrd(panel.p,3)
yhi = subwrd(panel.p,4)
'set line 0'
'draw recf 'xlo-0.4' 'ylo-0.8' 'xhi+0.4' 'ylo-0.02

* Next Panel: SLP
getseries(slpe,slp,1000)
'define slp = slp*0.01'
p = p - 1
'set parea 'panel.p
'set vpage off'
'set lon 'hilon
'set lat 'hilat
'set grid on'
'set gxout contour'
vrng(slp,slp)
'set ccolor 11'
'set cmark 0'
'set t '_t1-0.5' '_t2+0.5
'd slp'

* Draw a rectangle over the x-axis labels
xlo = subwrd(panel.p,1)
xhi = subwrd(panel.p,2)
ylo = subwrd(panel.p,3)
yhi = subwrd(panel.p,4)
'set line 0'
'draw recf 'xlo-0.4' 'ylo-0.8' 'xhi+0.4' 'ylo-0.02

* Next Panel: Surface Wind Speed
p = p - 1
getseries(u10m,ubot,1000)
getseries(v10m,vbot,1000)
'set parea 'panel.p
'set vpage off'
'set grads off'
if (_units = 'e')
  'define ubot = 2.2374*ubot'
  'define vbot = 2.2374*vbot'
endif
'define wind = mag(ubot,vbot)'
vrng(wind,wind)
'set ccolor 26'
'set cmark 7'
'set grid on'
'set t '_t1-0.5' '_t2+0.5
'set gxout contour'
'd wind'

* Draw a rectangle over the x-axis labels
xlo = subwrd(panel.p,1)
xhi = subwrd(panel.p,2)
ylo = subwrd(panel.p,3)
yhi = subwrd(panel.p,4)
'set line 0'
'draw recf 'xlo-0.4' 'ylo-0.8' 'xhi+0.4' 'ylo-0.02

* Next Panel: 2m Temperatures
getseries(t2m,t2m,1000)
getseries(rh2m,rh2m,1000)
p = p - 1
'set parea 'panel.p
'set vpage off'
'set frame on'
'set grads off'
'set ylab on'
'set gxout line'
'set grid off'
if (_units = 'e')
  'define t2mc = const((t2m-273.16),0,-u)'
  'define t2m  = const((t2m-273.16)*9/5+32,0,-u)'
  'define dewpt = t2mc-((14.55+0.114*t2mc)*(1-0.01*rh2m)+pow((2.5+0.007*t2mc)*(1-0.01*rh2m),3)+(15.9+0.117*t2mc)*pow((1-0.01*rh2m),14))'
  'define dewpt = dewpt*9/5+32'
else
  'define t2mf = const((t2m-273.16)*1.8+32,0,-u)'
  'define t2m  = const((t2m-273.16),0,-u)'
  'define dewpt = t2m-((14.55+0.114*t2m)*(1-0.01*rh2m)+pow((2.5+0.007*t2m)*(1-0.01*rh2m),3)+(15.9+0.117*t2m)*pow((1-0.01*rh2m),14))'

endif
vrng(t2m,dewpt)
'set t '_t1-0.5' '_t2+0.5
if (_units = 'e')
  'set ylint 10'
  'set gxout linefill'
  expr = 't2m;const(t2m'
  'set lfcols  9 0' ; 'd 'expr',-60,-a)'
  'set lfcols  9 0' ; 'd 'expr',-60,-a)'
  'set lfcols 14 0' ; 'd 'expr',-10,-a)'
  'set lfcols  4 0' ; 'd 'expr',0,-a)'
  'set lfcols 11 0' ; 'd 'expr',10,-a)'
  'set lfcols  5 0' ; 'd 'expr',20,-a)'
  'set lfcols 13 0' ; 'd 'expr',30,-a)'
  'set lfcols  3 0' ; 'd 'expr',40,-a)'
  'set lfcols 10 0' ; 'd 'expr',50,-a)'
  'set lfcols  7 0' ; 'd 'expr',60,-a)'
  'set lfcols 12 0' ; 'd 'expr',70,-a)'
  'set lfcols  8 0' ; 'd 'expr',80,-a)'
  'set lfcols  2 0' ; 'd 'expr',90,-a)'
  'set lfcols  6 0' ; 'd 'expr',100,-a)'
  'set gxout line'
  'set ccolor 15'
  'set cstyle 3'
  'set cmark 0'
  'd t2m'
else
  'set ylint 5'
  'set gxout linefill'
  expr = 't2m;const(t2m'
  'set lfcols  9 0' ; 'd 'expr',-60,-a)'
  'set lfcols 14 0' ; 'd 'expr',-25,-a)'
  'set lfcols  4 0' ; 'd 'expr',-20,-a)'
  'set lfcols 11 0' ; 'd 'expr',-15,-a)'
  'set lfcols  5 0' ; 'd 'expr',-10,-a)'
  'set lfcols 13 0' ; 'd 'expr',-5,-a)'
  'set lfcols  3 0' ; 'd 'expr',0,-a)'
  'set lfcols 10 0' ; 'd 'expr',5,-a)'
  'set lfcols  7 0' ; 'd 'expr',10,-a)'
  'set lfcols 12 0' ; 'd 'expr',15,-a)'
  'set lfcols  8 0' ; 'd 'expr',20,-a)'
  'set lfcols  2 0' ; 'd 'expr',25,-a)'
  'set lfcols  6 0' ; 'd 'expr',30,-a)'
  'set gxout line'
  'set ccolor 15'
  'set cstyle 3'
  'set cmark 0'
  'd t2m'
endif
'set grid on'
'set cmark 8'
'set ccolor 2'
'd t2m'
'set ccolor 97'
'set cmark 9'
'd dewpt'

* Draw a rectangle over the x-axis labels
xlo = subwrd(panel.p,1)
xhi = subwrd(panel.p,2)
ylo = subwrd(panel.p,3)
yhi = subwrd(panel.p,4)
'set line 0'
'draw recf 'xlo-0.4' 'ylo-0.8' 'xhi+0.4' 'ylo-0.02

* Back up to Previous Panel: 10m Wind Barbs
p = p + 1
'set parea 'panel.p
'set ccolor 1'
lap1 = hilat + 0.1
lam1 = hilat - 0.1
'set lon 'hilon ;* ??
'set lat 'lam1' 'lap1
'set frame off'
'set grid off'
'set gxout barb'
'set xyrev on'
'set xlab off'
'set ylab off'
if (_units = 'e')
  'd 2.2374*u10m.1;2.2374*v10m.1'
else
  'd u10m.1;v10m.1'
endif

* Reset dimension and graphics parameters
'set lat 'hilat
'set lon 'hilon
'set vpage off'
'set frame on'
'set grads off'
'set ylab on'
'set xlab on'
'set gxout line'
'set grid off'

* Skip to Next Panel: 2m Relative Humidity
p = p - 2
'set parea 'panel.p
*'set vpage off'
*'set grads off'
rh2vrng(rh2m)
'set gxout linefill'
'set lfcols 20 0' ; 'd rh2m;const(rh2m,00.01,-a)'
'set lfcols 21 0' ; 'd rh2m;const(rh2m,20.01,-a)'
'set lfcols 22 0' ; 'd rh2m;const(rh2m,30.01,-a)'
'set lfcols 23 0' ; 'd rh2m;const(rh2m,40.01,-a)'
'set lfcols 24 0' ; 'd rh2m;const(rh2m,50.01,-a)'
'set lfcols 25 0' ; 'd rh2m;const(rh2m,60.01,-a)'
'set lfcols 26 0' ; 'd rh2m;const(rh2m,70.01,-a)'
'set lfcols 27 0' ; 'd rh2m;const(rh2m,80.01,-a)'
'set lfcols 28 0' ; 'd rh2m;const(rh2m,90.01,-a)'
'set ccolor 28'
'set gxout line'
'set grid on'
'set cmark 2'
'd rh2m'

* Draw a rectangle over the x-axis labels
xlo = subwrd(panel.p,1)
xhi = subwrd(panel.p,2)
ylo = subwrd(panel.p,3)
yhi = subwrd(panel.p,4)
'set line 0'
'draw recf 'xlo-0.4' 'ylo-0.8' 'xhi+0.4' 'ylo-0.02

* Final Panel: Precipitation
getseries('p',pt,1000)
getseries(pc,pc,1000)
p = p - 1
'set parea 'panel.p
'set vpage off'
'set grid on'
'set grads off'
'define ptot  = 0.5*(pt+abs(pt))'
'define pconv = 0.5*(pc+abs(pc))'
if (_units = 'e')
  'define ptot  = const(ptot,0,-u)/25.4'
  'define pconv = const(pconv,0,-u)/25.4'
else
  'define ptot  = const(ptot,0,-u)'
  'define pconv = const(pconv,0,-u)'
endif

* Adjust the accumulation periods to be every three hours
fixprecip()

* Get Total Precipitation Range
'set gxout stat'
'd total'
data = sublin(result,8)
pmx = subwrd(data,5)
if (_units = 'e')
  if (pmx < 0.05) 
    pmx = 0.0499
  else
    pmx = pmx + (0.05*pmx)
  endif
else
  if (pmx < 1.0) 
    pmx = 0.99
  else
    pmx = pmx + (0.05*pmx)
  endif
endif
'set vrange 0 'pmx
incr = 0.01 * (math_nint(100*pmx/5))
'set ylint 'incr
'set t '_t1+0.5' '_t2+0.5

* Rain (Total Precipitation)
'set gxout bar'
'set barbase 0'
'set bargap 20'
'set ccolor 42'
'd total'

* Convective Precipitation
'set bargap 80'
'set ccolor 2'
'd conv'

* Draw all the Y-axis labels

* First panel
'set line 21' ; 'draw recf 0.4 7.65 0.62  8.18'
'set line 22' ; 'draw recf 0.4 7.65 0.58  8.18'
'set line 23' ; 'draw recf 0.4 7.65 0.535 8.18'
'set line 25' ; 'draw recf 0.4 7.65 0.49  8.18'
'set line 26' ; 'draw recf 0.4 7.65 0.445 8.18'
'set string 0 c 4 90' ; 'draw string 0.5 7.93 RH (%)'
'set string 2 l 4 90' ; 'draw string 0.5 8.36 T'
'set string 8 l 4 90' ; 'draw string 0.5 8.43 e'
'set string 5 l 4 90' ; 'draw string 0.5 8.50 m'
'set string 4 l 4 90' ; 'draw string 0.5 8.62 p'
'set string 9 l 4 90' ; 'draw string 0.5 8.69 .'
if (_units = 'e')
  'set string 2 l 4 90' ; 'draw string 0.5 8.79 (F)'
  'set string 1 c 4 90' ; 'draw string 0.5 9.53 Wind (mph)'
else
  'set string 2 l 4 90' ; 'draw string 0.5 8.79 (C)'
  'set string 1 c 4 90' ; 'draw string 0.5 9.53 Wind (m/s)'
endif
'draw string 0.75 8.63 `1m i l l i b a r s'

* Next Panel
'set strsiz 0.08 0.12'
p = npanels - 1 
ylo = subwrd(panel.p,3)
yhi = subwrd(panel.p,4)
ymid = ylo + (yhi-ylo)/2
'set string  5 c 4 90'  
'draw string 0.5 'ymid' Thickness'
'draw string 0.3 'ymid' 1000-500mb'
'set string  1 c 4 90'  
'draw string 0.74 'ymid' (dm)'

* Next Panel
p = p - 1 
ylo = subwrd(panel.p,3)
yhi = subwrd(panel.p,4)
ymid = ylo + (yhi-ylo)/2
'set string 8 c 4 90' ; 'draw string 0.15 'ymid' Total-totals'

* Total-Totals Y-axis Legend
'set strsiz 0.08 0.11'
'set string 15 r 4 0' ; 'draw string 0.45 'ymid' 40'
'set string  7 r 4 0' ; 'draw string 0.45 'ymid-0.133' 30'
'set string  7 r 4 0' ; 'draw string 0.45 'ymid-0.266' 20'
'set string  8 r 4 0' ; 'draw string 0.45 'ymid+0.133' 50'
'set string  8 r 4 0' ; 'draw string 0.45 'ymid+0.266' 60'
'set strsiz 0.08 0.12'
'set string 2 c 4 90' ; 'draw string 0.69 'ymid' Lifted Index'

* Next Panel
p = p - 1 
ylo = subwrd(panel.p,3)
yhi = subwrd(panel.p,4)
ymid = ylo + (yhi-ylo)/2
'set string 11 c 4 90' ; 'draw string 0.3 'ymid' SLP'
'set string  1 c 4 90' ; 'draw string 0.6 'ymid' (mb)'

* Next Panel
p = p - 1 
ylo = subwrd(panel.p,3)
yhi = subwrd(panel.p,4)
ymid = ylo + (yhi-ylo)/2
'set string 26 c 4 90' ; 'draw string 0.15 'ymid' 10m Wind'
'set string 26 c 4 90' ; 'draw string 0.35 'ymid' Speed'
'set string  1 c 4 90' ; 'draw string 0.55 'ymid' & Barbs'
if (_units = 'e')
  'draw string 0.75 'ymid' (mph)'
else
  'draw string 0.75 'ymid' (m/s)'
endif

* Next Panel
p = p - 1 
ylo = subwrd(panel.p,3)
yhi = subwrd(panel.p,4)
ymid = ylo + (yhi-ylo)/2
'set string  2 c 4 90' ; 'draw string 0.15 'ymid' 2m Temp '
'set string 97 c 4 90' ; 'draw string 0.35 'ymid' 2m DewPt '
if (_units = 'e')
  'set string 1 c 4 90'
  'draw string 0.75 'ymid' (F)'
else
  'set string 1 c 4 90'
  'draw string 0.75 'ymid' (C)'
endif

* Next Panel
p = p - 1 
ylo = subwrd(panel.p,3)
yhi = subwrd(panel.p,4)
ymid = ylo + (yhi-ylo)/2
'set string 26 c 4 90' ; 'draw string 0.35 'ymid' 2m RH'
'set string  1 c 4 90' ; 'draw string 0.75 'ymid' (%)'

* Bottom Panel
p = p - 1 
ylo = subwrd(panel.p,3)
yhi = subwrd(panel.p,4)
ymid = ylo + (yhi-ylo)/2
if (_units = 'e')
  'set string 1 c 4 90' ; 'draw string .85 'ymid' (inches)'
else
  'set string 1 c 4 90' ; 'draw string .85 'ymid' (mm)'
endif

'set string  1 c 4 90' ; 'draw string 0.15 'ymid' 3hr Precip'
'set string 42 c 4 90' ; 'draw string 0.35 'ymid' Total'
'set string  2 c 4 90' ; 'draw string 0.55 'ymid' Convective'

* Draw Labels at the top of the page
'set string 1 r 1 0'
'set strsiz 0.14 .17'
label = 'ETA 0-84hr Forecast Meteogram for ('
if (hilon < 0)  ; label = label%hilon*(-1.0)'W, ' ; endif
if (hilon >= 0) ; label = label%hilon'E, ' ; endif
if (hilat < 0)  ; label = label%hilat*(-1.0)'S)'; endif
if (hilat >= 0) ; label = label%hilat'N)' ; endif
'draw string 8.15 10.75 'label

* Draw the station label
'set strsiz 0.18 0.22'
'set string 21 l 12 0' ; 'draw string 0.12 10.79 `1'name
'set string  1 l  8 0' ; 'draw string 0.10 10.81 `1'name

* Remove the dummy files
'!rm -f dummy.ctl'
'!rm -f dummy.dat'

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* END OF MAIN SCRIPT
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

function setcols(args)
'set rgb 20 234 245 234'
'set rgb 21 200 215 200'
'set rgb 22 160 205 160'
'set rgb 23 120 215 120'
'set rgb 24  80 235  80'
'set rgb 25   0 255   0'
'set rgb 26   0 195   0'
'set rgb 27   0 160   0'
'set rgb 28   0 125   0'

'set rgb 30 255 160 120'
'set rgb 31 160 120 255'
'set rgb 32 160 180 205'

'set rgb 42  32 208  32'
'set rgb 43 208  32 208'
'set rgb 44  64  64 255'
'set rgb 45 255 120  32'
'set rgb 46  32 208 208'
'set rgb 47 240 240   0'

'set rgb 96 139 115  85'
'set rgb 97 100 100 100'
'set rgb 98  64  64  96'
'set rgb 99 254 254 254'
return

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
function vrng(f1,f2)
'set gxout stat'
'd 'f1
data = sublin(result,8)
ymx = subwrd(data,5)
ymn = subwrd(data,4)
'd 'f2
data = sublin(result,8)
zmx = subwrd(data,5)
zmn = subwrd(data,4)
if (zmx > ymx) ; ymx = zmx ; endif
if (zmn < ymn) ; ymn = zmn ; endif
dy = ymx-ymn
ymx = ymx + 0.08 * dy
ymn = ymn - 0.08 * dy
if ((ymx-ymn)/2.2 < 1)
  incr = (ymx-ymn)/4
  incr = 0.01 * (math_nint(100*incr))
else
  incr = math_nint((ymx-ymn)/4)
endif
'set vrange 'ymn' 'ymx
'set ylint 'incr
*say 'vrng: 'ymn' 'ymx' 'incr
if (ymn=0 & ymx=0 & incr=0)
*  say 'vrng: resetting zeros to -.9 .9 1'
  'set vrange -.9 .9'
  'set ylint 1'
endif
'set gxout line'
return

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
function rh2vrng(f1)
'set gxout stat'
'd 'f1
data = sublin(result,8)
ymn = subwrd(data,4)
ymx = subwrd(data,5)
if (ymn < 20) 
  miny = 0 
  'set ylevs 20 40 60 80'
endif
if (ymn >= 20 & ymn < 30) 
  miny = 20 
  'set ylevs 30 50 70 90'
endif
if (ymn >= 30 & ymn < 40) 
  miny = 30 
  'set ylevs 40 50 60 70 80 90'
endif
if (ymn >= 40 & ymn < 50) 
  miny = 40 
  'set ylevs 50 60 70 80 90'
endif
if (ymn >= 50 & ymn < 60) 
  miny = 50
  'set ylevs 60 70 80 90'
endif
if (ymn >= 60) 
  miny = 60
  'set ylevs 70 80 90'
endif
'set vrange 'miny' 'ymx+3
return

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
function getctl(handle)
line = 1
found = 0
while (!found)
  info = sublin(_ctl,line)
  if (subwrd(info,1)=handle)
    _handle = info
    found = 1
  endif
  line = line + 1
endwhile
return _handle

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
function getgrid(dodsvar,myvar)

'set lon '_xdim
'set lat '_ydim
'set lev '_zgrd
'set time '_tdim

* Write the variable to a file
'set gxout fwrite'
'set fwrite dummy.dat'
'd 'dodsvar
'disable fwrite'

* Write a descriptor file 
rc = write(dummy.ctl,'dset ^dummy.dat')
rc = write(dummy.ctl,_undef,append)
rc = write(dummy.ctl,'xdef 1 linear 1 1',append)
rc = write(dummy.ctl,'ydef 1 linear 1 1',append)
rc = write(dummy.ctl,_zdef,append)
rc = write(dummy.ctl,_tdef,append)
rc = write(dummy.ctl,'vars 1',append)
rc = write(dummy.ctl,'dummy '_newzsize' -999 dummy',append)
rc = write(dummy.ctl,'endvars',append)
rc = close (dummy.ctl)

* Open the dummy file, define variable, close dummy file
'open dummy.ctl'
line = sublin(result,2)
dummyfile = subwrd(line,8)
'set dfile 'dummyfile
'set lon 1'
'set lat 1'
'set lev '_zbot' '_ztop
'set time '_time1' '_time2
'define 'myvar' = dummy.'dummyfile
'close 'dummyfile
'set dfile 1'
return

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
function getetarh(dodsvar,myvar)

* swap out original pressure vars
tmpzgrd = _zgrd
tmpzdef = _zdef
tmpzbot = _zbot
tmpztop = _ztop
tmpzsize = _newzsize

* retrieve rh data over the rh pressure range
_zgrd = _rhzgrd
_zdef = _trhzdef
_ztop = _rhztop
_zbot = _rhzbot
_newzsize = _trhzsize
getgrid(dodsvar,tmprh)

* swap in original pressure vars
_zgrd = tmpzgrd
_zdef = tmpzdef
_zbot = tmpzbot
_ztop = tmpztop
_newzsize = tmpzsize

'set lon '_xdim
'set lat '_ydim
'set lev '_rhzgrd
'set time '_tdim

* Write the variable to a file
'set gxout fwrite'
'set fwrite dummy.dat'
t = _t1
while (t <= _t2)
  'set t 't
  z = 1
  while (z <= _newrhzsize)
    level = subwrd(_rhlevs,z) 
    'set lev 'level
    'd tmprh'
    z = z + 1
  endwhile
  z = 1
  while (z <= _newrhzsize)
    level = subwrd(_rhlevs,z) 
    'set lev 'level
    'd u'
    z = z + 1
  endwhile
  z = 1
  while (z <= _newrhzsize)
    level = subwrd(_rhlevs,z) 
    'set lev 'level
    'd v'
    z = z + 1
  endwhile
  t = t + 1
endwhile
'disable fwrite'

* Write a descriptor file 
rc = write(dummy.ctl,'dset ^dummy.dat')
rc = write(dummy.ctl,_undef,append)
rc = write(dummy.ctl,'xdef 1 linear 1 1',append)
rc = write(dummy.ctl,'ydef 1 linear 1 1',append)
rc = write(dummy.ctl,_rhzdef,append)
rc = write(dummy.ctl,_tdef,append)
rc = write(dummy.ctl,'vars 3',append)
rc = write(dummy.ctl,'dummy '_newrhzsize' -999 dummy',append)
rc = write(dummy.ctl,'dummyu '_newrhzsize' -999 dummy',append)
rc = write(dummy.ctl,'dummyv '_newrhzsize' -999 dummy',append)
rc = write(dummy.ctl,'endvars',append)
rc = close (dummy.ctl)

* Open the dummy file, define variable, close dummy file
'open dummy.ctl'
line = sublin(result,2)
dummyfile = subwrd(line,8)
'set dfile 'dummyfile
'set lon 1'
'set lat 1'
'set lev '_rhzbot' '_rhztop
'set time '_time1' '_time2
'q dims'
'define 'myvar' = dummy.'dummyfile
if (_units = 'e')
  'define urh = dummyu.'dummyfile'*2.2374'
  'define vrh = dummyv.'dummyfile'*2.2374'
else
  'define urh = dummyu.'dummyfile
  'define vrh = dummyv.'dummyfile
endif
'close 'dummyfile
'set dfile 1'

return

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
function getseries(dodsvar,myvar,level)

'set lon '_xdim
'set lat '_ydim
'set lev 'level' 'level
'set time '_tdim

* Write the variable to a file
'set fwrite dummy.dat'
'set gxout fwrite'
'd 'dodsvar
'disable fwrite'

* Write a descriptor file 
rc = write(dummy.ctl,'dset ^dummy.dat')
rc = write(dummy.ctl,_undef,append)
rc = write(dummy.ctl,'xdef 1 linear 1 1',append)
rc = write(dummy.ctl,'ydef 1 linear 1 1',append)
rc = write(dummy.ctl,'zdef 1 linear 1 1',append)
rc = write(dummy.ctl,_tdef,append)
rc = write(dummy.ctl,'vars 1',append)
rc = write(dummy.ctl,'dummy 0 -999 dummy',append)
rc = write(dummy.ctl,'endvars',append)
rc = close(dummy.ctl)

* Open the dummy file, define variable, close dummy file
'open dummy.ctl'
line = sublin(result,2)
dummyfile = subwrd(line,8)
'set dfile 'dummyfile
'set lon 1'
'set lat 1'
'set lev 'level
'set time '_time1' '_time2
'define 'myvar' = dummy.'dummyfile
'close 'dummyfile
'set dfile 1'
'set gxout contour'

return

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
function fixprecip

'set lon '_xdim
'set lat '_ydim
'set z 1'

* Write the precip variables to a file
'set fwrite dummy.dat'
'set gxout fwrite'
'set t '_t1
'd ptot'; 'd pconv'
index = 2
while (index <= 26)  ;* 26 is 3 timesteps (or 9 hours) before the final time 
  'set t 'index
  'd ptot'; 'd pconv'
  'set t 'index+1
  'define tmp = ptot-ptot(t-1)  '; 'd tmp'
  'define tmp = pconv-pconv(t-1)'; 'd tmp'
  'set t 'index+2
  'define tmp = ptot-ptot(t-1)  '; 'd tmp'
  'define tmp = pconv-pconv(t-1)'; 'd tmp'
  'set t 'index+3
  'define tmp = ptot-ptot(t-1)  '; 'd tmp'
  'define tmp = pconv-pconv(t-1)'; 'd tmp'
  index = index + 4
endwhile
'disable fwrite'

* Write a descriptor file 
rc = write(dummy.ctl,'dset ^dummy.dat')
rc = write(dummy.ctl,_undef,append)
rc = write(dummy.ctl,'xdef 1 linear 1 1',append)
rc = write(dummy.ctl,'ydef 1 linear 1 1',append)
rc = write(dummy.ctl,'zdef 1 linear 1 1',append)
rc = write(dummy.ctl,_tdef,append)
rc = write(dummy.ctl,'vars 2',append)
rc = write(dummy.ctl,'dummy1 0 -999 dummy',append)
rc = write(dummy.ctl,'dummy2 0 -999 dummy',append)
rc = write(dummy.ctl,'endvars',append)
rc = close(dummy.ctl)

* Open the dummy file, define variable, close dummy file
'open dummy.ctl'
line = sublin(result,2)
dummyfile = subwrd(line,8)
'set dfile 'dummyfile
'set lon 1'
'set lat 1'
'set z 1'
'set time '_time1' '_time2
'define 'total' = dummy1.'dummyfile
'define 'conv'  = dummy2.'dummyfile
'close 'dummyfile
'set dfile 1'
'set gxout contour'

return

