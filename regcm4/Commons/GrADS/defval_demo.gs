* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* defval_demo.gs
*
* Illustrates the use of the 'q defval' and 'set defval' commands
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

'open model.ctl'

* restrict the dimension environment to see things more clearly
xmin =  0.5
xmax = 15.5
ymin =  0.5
ymax = 10.5
'set x 'xmin' 'xmax
'set y 'ymin' 'ymax

* the "defval" commands will only work with 2-D defined variables 
'define var = ps'

* Display the variable with colored grid cells and their values
'c'
'set grid off'
'set mproj off'
'set xaxis 'xmin' 'xmax
'set yaxis 'ymin' 'ymax
'set gxout grfill'
'd var'
'set gxout grid'
'set digsiz .14'
'd var'

* Use the mouse to click on a grid point to change
say 'Click on any grid point'
'q pos'
xscreen = subwrd(result,3)
yscreen = subwrd(result,4)

* Convert screen positions to grid coordinates
'q xy2gr 'xscreen' 'yscreen
xgrid = subwrd(result,3)
ygrid = subwrd(result,6)

* Round the grid values to the nearest integer
gx = math_nint(xgrid)
gy = math_nint(ygrid)

* Get the value of the defined variable 
'q defval var 'gx' 'gy
val = subwrd(result,3)
say 'The value at grid point ('gx','gy') is --> 'val

* Ask for a new replacement value and assign it
prompt 'Enter a new value --> ' 
pull newval
'set defval var 'gx' 'gy' 'newval

* Display the newly updated variable
'c'
'set xaxis 'xmin' 'xmax
'set yaxis 'ymin' 'ymax
'set gxout grfill'
'd var'
'set gxout grid'
'd var'

