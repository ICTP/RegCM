* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* basemap.gs
*
* This script overlays a land or ocean mask that exactly fits 
* the low-resolution coastal outlines. The masks are composed 
* of hundreds of polygons that are specified in the accompanying 
* ascii files lpoly.asc and opoly.asc. It will work with any 
* scaled or latlon map projection. 
*
* If you are using Grads version 1.8 or higher, this script
* will also work properly with the robinson projection and 
* polar stereographic projections from 0-90, 15-90, and 20-90 
* (North and South). 
*
* Other projections will work but are not guaranteed because 
* GrADS may not clip the basemap properly. A solution to this 
* problem is to use "set mpvals" to override the dimension 
* environment limits. For example:
*    set mproj nps
*    set lon -180 180
*    set lat 0 90
*    set mpvals -180 180 60 90 
*    display <something>
*    basemap L
* The resulting plot will be a properly clipped square.
*
* A new high-resolution polygon data set is available covering the
* region between 15N-53N, 130W-60W. To use this option, you must 
* first download the required ascii files containing the 
* hires polygons (opoly_hires.asc and lpoly_hires.asc) and then put 
* an additional argument ("hi" or "lo") and the end of the complete
* argument string. For example:
* 
*    set mpdset hires
*    display <something>
*    basemap L 15 0 hi
* 
* An additional option is to mask out the Mexican and
* Canadian land regions surrounding the US, so that only the
* conterminous states are seen. To use this feature, change
* your land polygon file from lpoly.asc to lpoly_US.asc and
* then run basemap twice:
*    basemap o 0 0  (<- that's O zero zero)
*    basemap L 0 0
* This will only work properly if your domain is within the
* boundaries 20N-50N, 130W-60W. Low-res maps only. 
* 
* Download all polygon files from ftp://grads.iges.org/grads/scripts/
*      
* Written 12/2000 by Jennifer M. Adams, jma@cola.iges.org
* Updated 01/2001, 05/2001, 09/2001
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
function main(args)

* There are defaults for the colors and resolution
* User must specify which mask
  if (args='') 
    say 'Usage: basemap L(and)/O(cean) <fill_color> <outline_color> <hi/lo>'
    return
  else 
    type    = subwrd(args,1)
    color   = subwrd(args,2)
    outline = subwrd(args,3)
    res     = subwrd(args,4)
    if (color = '')   ; color = 15  ; endif
    if (outline = '') ; outline = 0 ; endif
    if (res = 'HI' | res = 'hi') ; hires = 1 ; lores = 0 ; endif
    if (res = 'LO' | res = 'lo') ; hires = 0 ; lores = 1 ; endif
    if (res = '') ; hires = 0 ; lores = 1  ; endif
  endif 

  if (type = 'L' | type = 'l') 
    if (lores) ; file = 'lpoly.asc' ; endif
    if (hires) ; file = 'lpoly_hires.asc' ; endif
  endif
  if (type = 'O' | type = 'o') 
    if (lores) ; file = 'opoly.asc' ; endif
    if (hires) ; file = 'opoly_hires.asc' ; endif
  endif

* Make sure there's a plot already drawn
  'q gxinfo'
  line5 = sublin(result,5)
  line6 = sublin(result,6)
  xaxis = subwrd(line5,3)
  yaxis = subwrd(line5,6)
  proj  = subwrd(line6,3)
  if (xaxis = 'None' | yaxis = 'None') 
    say 'Error: Please display a variable before using basemap'
    return
  endif

* See what version of Grads is running 
  'q config'
  line = sublin(result,1)
  word = subwrd(line,2)
  version = substr(word,2,3)
  if (version >= 1.8) 
    newgrads = 1 
  else 
    newgrads = 0
  endif
   
* See if map projection will be supported
  if (newgrads = 0) 
    if (proj != 1 & proj != 2)
      say 'Error: Only scaled or latlon projections are supported with GrADS v'version
      return
    endif
  endif

* Clip image accordingly
  'q gxinfo'
  line3 = sublin(result,3)
  line4 = sublin(result,4)
  x1 = subwrd(line3,4)
  x2 = subwrd(line3,6)
  y1 = subwrd(line4,4)
  y2 = subwrd(line4,6)
  'set clip 'x1' 'x2' 'y1' 'y2

* Read the first record from the polygon file
  result = read(file)
  rc = sublin(result,1)
  rc = subwrd(rc,1)
  if (rc!=0)
    say 'Error reading 'file
    return
  endif
  nwcmd = sublin(result,2)

* Read subsequent records, allowing for read input buffer overflow
  flag = 1
  while (flag)
    ignore = 0
    wcmd = nwcmd
    while(1)
      result = read(file)
      rc = sublin(result,1)
      rc = subwrd(rc,1)
      if (rc!=0)
        flag = 0
        break
      else 
        nwcmd = sublin(result,2)
        if (subwrd(nwcmd,5) != 'draw') 
          wcmd = wcmd % nwcmd
        else
          break
        endif
      endif
    endwhile

*   Get the lat/lon range of the current dimension environment
    'q dims'
    line1 = sublin(result,2)
    line2 = sublin(result,3)
    minlon = subwrd(line1,6)
    maxlon = subwrd(line1,8)
    minlat = subwrd(line2,6)
    maxlat = subwrd(line2,8)

*   The range of the polygon is specified in the first four words of the record
    minwx = subwrd(wcmd,1)
    maxwx = subwrd(wcmd,2)
    minwy = subwrd(wcmd,3)
    maxwy = subwrd(wcmd,4)

*   If the polygon is outside the current dimension, ignore it
    if (minwx >= maxlon) ; ignore = 1 ; endif 
    if (maxwx <= minlon) ; ignore = 1 ; endif 
    if (minwy >= maxlat) ; ignore = 1 ; endif 
    if (maxwy <= minlat) ; ignore = 1 ; endif 

    if (!ignore)    
      count = 7
      nvert = 1
      if (newgrads)  
        cmd = 'draw mappoly ' 
      else 
        cmd = 'draw polyf '   
      endif 
      while (1)
        countx = count
        county = count + 1
        wx = subwrd(wcmd,countx)
        wy = subwrd(wcmd,county)
        if ((wx = '') | (wy = ''))
          break 
        endif

*       Convert world coordinates to screen coordinates if necessary
        if (newgrads)  
          sx = wx
          sy = wy
        else 
          'q w2xy 'wx' 'wy
          sx = subwrd(result,3)
          sy = subwrd(result,6)
        endif

*       Append the coordinates to the draw command
        cmd = cmd%sx' 'sy' '
        count = count + 2
      endwhile   

*     Draw the polygon
      'set line 'color
      cmd
    endif
  endwhile

* Draw the continental outline
  if (lores) ; 'set mpdset lowres' ; endif 
  if (hires) ; 'set mpdset hires'  ; endif 
  'set mpt * 'outline
  'draw map'

* Close the polygon file 
  rc = close(file)
  return

* THE END *
