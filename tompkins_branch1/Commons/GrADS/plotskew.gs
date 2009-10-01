function plotskew(sndtemp,snddewp,sndspd,snddir)
*************************************************************************
*
* GrADS Script to Plot a SkewT/LogP Diagram      
*
* Bob Hart 
* Penn State University / Dept of Meteorology
* Last Update:  January 23, 2001
*
* Recent Changes:
*
* 01/23/01 - Fixed a small bug in the theta-e calculation.  
*            Errors averaged 0.5-3K.  Thank you George Bryan.
*
* 11/10/99 - Change in calculation method for CAPE/CIN.  Trapezoid 
*            integration method is now used.  Speeds up execution
*            by 25%, and increases accuracy by 5-10%.
*
* 10/18/99 - Minor glitch fixed that occasionally caused crash.
*
*  8/26/99 - Datasets with missing data can now be used.
*
* Features:
*   - All features of standard skewt/logp plot
*   - RH sounding 
*   - LCL location
*   - Parcel trajectory for both sfc based convection and elevated from
*     most unstable level (highest theta-e level reported)
*   - Stability indices and precipitable water calculations
*   - CAPE & CIN Calculations
*   - Wind Profile
*   - Hodograph / Hodograph scaling
*   - Helicity and SR Helicity Calculations and Display
*   - Color aspects of output
*   - Line Thickness, style aspects of output
*   - Can be run in either PORTRAIT or LANDSCAPE mode.
*
* There are numerous tunable parameters below to change the structure
* and output for the diagram.
*
* Function Arguments:
*    sndtemp - temperature data (Celsius) as a function of pressure
*    snddewp - dewpoint data (Celsius) as a function of pressure      
*    sndspd  - wind speed data (knots) as a function of pressure
*    snddir  - wind direction data as a function of pressure
*
* Use '-1' for any of the above 4 arguments to indicate that you
* are not passing that variable.  The appropriate options will 
* be ignored based on your specifying '-1' for that variable.
*
* NOTE:  Make sure to set the vertical range of the plot before running.
*        I.e., "SET LEV 1050 150", for example.   This does not have to 
*        be limited to the pressure range of your data.
*
* Labelling:  Pressure/Height is labelled along left side.  Temperature is 
*             labelled along bottom.  Mixing ratio is labelled along right
*             side/top.
*
* 
* PROBLEMS:  First check out the web page for the script (which also
*            has a link to a FAQ with answers to many common questions
*            about using the script):
*            http://www.ems.psu.edu/~hart/skew.html
*
* Please send any further problems, comments, or suggestions to 
* <hart@ems.psu.edu>
*
* ACKNOWLEDGMENTS:  Thanks go to the innumerable users who have helped
* fine tune the script from the horrible mess from which it began.
* In particular, thanks go out to Steve Lord (NCEP), Mike Fiorino (ECMWF),
* George Bryan (PSU), Davide Sacchetti (CMIRL), and Enrico Minguzzi (CMIRL).
*
**************************************************************************
*           !!!!!   BEGINNING OF USER-SPECIFIED OPTIONS  !!!!!!
**************************************************************************
*
* --------------------- Initialization options  ----------------------
*
* ClrScrn = Whether to clear the screen before drawing diagram
*           [1 = yes, 0 = no]

ClrScrn = 1

*
* ------------------- Define Skew-T Diagram Shape/Slope-----------------
*
* (P1,T1) = Pres, Temp of some point on left-most side
* (P2,T2) = Pres, Temp of some point on right-most side
* (P3,T3) = Pres, Temp of some point in diagram which is mid-point 
*           in the horizontal between 1 and 2.
*
* P1, P2, P3 are in mb ; T1, T2, T3 are in Celsius
*
* These define the SLOPE and WIDTH of the diagram as you see it but DO NOT
* DEFINE THE HEIGHT of the diagram as you see it.  In other words,
* 1 and 2 do NOT necessarily need to be at the bottom of the diagram and
* 3 does NOT necessarily need to be at the top.  THE VERTICAL PRESSURE 
* RANGE OF THE SKEWT AS YOU SEE IT IS DETERMINED BY YOUR 'SET Z ...'  
* COMMAND OR THE 'SET LEV ...' COMMAND BEFORE RUNNING THIS SCRIPT.
*
*    _______________________
*   |                       |
*   |                       |
*   |           3           |
*   |                       |
*   |                       |
*   |                       |
*   |                       |
*   |                       |
*   |                       |
*   |                       |
*   |                       |
*   |1                     2|
*   |                       |
*   |_______________________|
*   
*
* A good set of defining points are given below.   Feel free
* to experiment with variations.


P1 = 1000
T1 = -40

P2 = 1000
T2 = 40

P3 = 200
T3 = -50

* Another good set of defining points suggested by George Bryan (PSU)
* are:
*
* P1 = 1000
* T1 = -30
*
* P2 = 1000
* T2 = 40
*
* P3 = 500
* T3 = -18

* ------------------- Contour Intervals / Levels --------------------------
*
* All variables below are contour intervals/levels for diagram
*
* Thetaint = interval for potential temperature lines
* Thetwint = interval for moist pseudo adiabats
* tempint  = interval for temperature lines
* wsclevs  = contour LEVELS for mixing ratio lines
*
*
thetaint= 10
thetwint= 5
tempint = 10
wsclevs = "1 2 3 4 6 8 10 15 20 25 30 35 40"
*
*
* ------------------------ Output Options --------------------------------
*
* All variables below are logical .. 1=yes, 0=no, unless otherwise
* specified.
*
* DrawBarb = Draw wind barbs along right side of plot
* DrawThet = Draw dry adiabats
* DrawThtw = Draw moist pseudo-adiabats
* DrawTemp = Draw temperature lines 
* DrawMix  = Draw mixing ratio lines
* DrawTSnd = Draw temperature sounding
* DrawDSnd = Draw dewpoint sounding
* DrawRH   = Draw relative humidity sounding
* DrawPrcl = Draw parcel path from surface upward
* DrawPMax = Draw parcel path from most unstable level upward
* DrawIndx = Display stability indices & CAPE
* DrawHeli = Calculate and display absolute and storm-relative helicity
* DrawHodo = Draw hodograph
* DrawPLev = Draw Pressure Levels 
* DrawZLev = Draw height levels and lines 
*            0 = no lines
*            1 = above ground level (AGL)
*            2 = above sea level (ASL)
* DrawZSTD = Draw Height levels using standard atm lapse rate           
* LblAxes  = Label the x,y axes (temperature, pressure,mixing ratio)
*
* ThtwStop = Pressure level at which to stop drawing Theta-w lines
* MixStop  = Pressure level at which to stop drawing Mixratio lines

DrawBarb= 1
DrawThet= 1
DrawThtw= 0
DrawTemp= 1
DrawMix = 1
DrawTSnd= 1
DrawDSnd= 1
DrawRH  = 0
DrawPrcl= 0
DrawPMax= 0
DrawIndx= 1
DrawHeli= 1
DrawHodo= 1
DrawPLev= 1
DrawZLev= 0
DrawZSTD= 0
LblAxes = 1

ThtwStop = 200
MixStop  = 600


* 
* -----------------  Sounding Geography options ------------------------
*
* SfcElev = Elevation above sea-level (meters) of lowest level reported
*           in sounding.  Used only if DrawZLev = 2

SfcElev = 0


*
* ------------------ Thermodynamic Index Options --------------------
* 
* All variables here are in inches.  Use -1 for the default values.
*
*  Text1XC = X-location of midpoint of K,TT,PW output box
*  Text1YC = Y-location of midpoint of K,TT,PW output box
*  Text2XC = X-Location of midpoint of surface indices output box
*  Text2YC = Y-location of midpoint of surface indices output box
*  Text3XC = X-Location of midpoint of most unstable level-based indices
*            output box
*  Text3YC = Y-location of midpoint of most unstable level-based indices
*            output box

Text1XC = -1
Text1YC = -1
Text2XC = -1
Text2YC = -1
Text3XC = -1
Text3YC = -1

*
* ----------------- Wind Barb Profile Options ----------------------------
*
* All variables here are in units of inches, unless otherwise specified
*
*  barbint = Interval for plotting barbs (in units of levels)
*  poleloc = X-Location of profile.  Choose -1 for the default.
*  polelen = Length of wind-barb pole
*  Len05   = Length of each 5-knot barb 
*  Len10   = Length of each 10-knot barb
*  Len50   = Length of each 50-knot flag
*  Wid50   = Width of base of 50-knot flag 
*  Spac50  = Spacing between 50-knot flag and next barb/flag 
*  Spac10  = Spacing between 10-knot flag and next flag
*  Spac05  = Spacing between 5-knot flag and next flag
*  Flagbase= Draw flagbase (filled circle) for each windbarb [1=yes, 0 =no] 
*  Fill50  = Solid-fill 50-knot flag [1=yes, 0=no]
*  barbline= Draw a vertical line connecting all the wind barbs [1=yes, 0=no]
*
barbint = 1
poleloc = -1
polelen = 0.35
len05   = 0.07
len10   = 0.15
len50   = 0.15
wid50   = 0.06
spac50  = 0.07
spac10  = 0.05
spac05  = 0.05
Fill50  = 1
flagbase= 1
barbline= 1

*
*
*---------------- Hodograph Options -------------------------------------
*
* All variables here are in units of inches, unless otherwise specified
*
* HodXcent= x-location of hodograph center.  Use -1 for default location.
* HodYcent= y-location of hodograph center.  Use -1 for default location.
* HodSize = Size of hodograph in inches 
* NumRing = Number of rings to place in hodograph (must be at least 1)
* HodRing = Wind speed increment of each hodograph ring
* HodoDep = Depth (above lowest level in mb) of end of hodograph trace
* TickInt = Interval (in kts) at which tick marks are drawn along the axes
*           Use 0 for no tick marks.
* TickSize= Size of tick mark in inches
* Text4XC = X-location of midpoint of hodograph text output. Use -1 for default.
* Text4YC = Y-location of midpoint of hodograph text output. Use -1 for default.

HodXcent= -1
HodYcent= -1
HodSize = 2
NumRing = 3
HodRing = 15
HodoDep = 300
TickInt = 5
TickSize= 0.05
Text4XC = -1
Text4YC = -1

*--------------- Helicity Options ---------------------------------------
*
* MeanVTop = Top pressure level (mb) of mean-wind calculation
* MeanVBot = Bottom pressure level (mb) of mean-wind calculation
* HelicDep = Depth in mb (above ground) of helicity integration
* StormMot = Type of storm motion estimation scheme.  Use following:  
*            0 = No departure from mean wind.
*            1 = Davies-Jones (1990) approach
* FillArrw = Whether to fill the arrowhead of the storm motion vector
*            [1 = yes, 0 = no]

MeanVTop= 300 
MeanVBot= 850
HelicDep= 300
StormMot= 1
FillArrw= 1

*
*---------------- Color Options ------------------------------------------
*
* ThetCol = Color of dry adiabats
* TempCol = Color of temperature lines
* MixCol  = Color of mixing ratio lines
* ThtwCol = Color of moist adiabats
* TSndCol = Color of Temperature Sounding 
* DSndCol = Color of Dewpoint Sounding
* RHCol   = Color of RH Sounding
* PrclCol = Color of parcel trace
* BarbCol = Color of wind barbs (choose -1 for color according to speed)
* HodoCol = Color of hodograph trace

ThetCol = 2
TempCol = 4
MixCol  = 7
ThtwCol = 3
TSndCol = 1 
DSndCol = 1
RHCol   = 3
PrclCol = 5
BarbCol = -1
HodoCol = 1

*
*-------------------- Line Style Options ------------------------------------
*
* GrADS Styles: 1=solid;2=long dash;3=short dash;4=long,short dashed;
*               5=dotted;6=dot dash;7=dot dot dash
*
* ThetLine = Line Style of dry adiabats
* TempLine = Line Style of temperature lines
* MixLine  = Line Style of mixing ratio lines
* ThtwLine = Line Style of moist adiabats
* TSndLine = Line Style of Temperature Sounding
* DSndLine = Line Style of Dewpoint Sounding
* RHLine   = Line Style of RH sounding
* PrclLine = Line Style of parcel trace
* HodoLine = Line Style of hodograph trace
*

ThetLine = 1
TempLine = 1
MixLine  = 5
ThtwLine = 3
TSndLine = 1
DSndLine = 1
RHLine   = 1
PrclLine = 1
HodoLine = 1

*
*------------------- Line Thickness Options---------------------------------
* GrADS Line Thickness: increases with increasing number. Influences 
*                       hardcopy output more strongly than screen output.
*
*
* ThetThk = Line Thickness of dry adiabats
* TempThk = Line Thickness of temperature lines
* MixThk  = Line Thickness of mixing ratio lines
* ThtwThk = Line Thickness of moist adiabats
* TSndThk = Line Thickness of temperature sounding
* DSndThk = Line thickness of dewpoint sounding
* RHThk   = Line thickness of RH sounding
* PrclThk = Line thickness of parcel trace
* HodoThk = Line thickness of hodograph trace
* BarbThk = Line thickness of wind barbs

ThetThk = 3
TempThk = 3
MixThk  = 3
ThtwThk = 3
TSndThk = 8 
DSndThk = 8
RHThk   = 8
PrclThk = 6 
HodoThk = 3 
BarbThk = 2

*
*------------------- Data Point Marker Options -----------------------------
* GrADS Marker Types: 0 = none ; 1 = cross ; 2 = open circle ; 
*                     3 = closed circle ; 4 = open square ; 5 = closed square
*                     6 = X ; 7 = diamond ; 8 = triangle ; 9 = none
*                    10 = open circle with vertical line ; 11 = open oval
*
* TSndMrk = Mark type of data point marker for temperature sounding
* DSndMrk = Mark type of data point marker for dewpoint sounding
* RHMrk   = Mark type of data point marker for relative humidity sounding
* MrkSize = Mark size (inches) of each data marker

TSndMrk = 0
DSndMrk = 0
RHMrk   = 0
MrkSize = 0.1


* !!!!! YOU SHOULD NOT NEED TO CHANGE ANYTHING BELOW HERE !!!!!
****************************************************************************

*-------------------------------------------
* grab user-specified environment dimensions
*-------------------------------------------

"q dims"
rec=sublin(result,2)
_xtype=subwrd(rec,3)
_xval=subwrd(rec,9)
rec=sublin(result,3)
_yval=subwrd(rec,9)
_ytype=subwrd(rec,3)
rec=sublin(result,4)
_ptype=subwrd(rec,3)
_pmax=subwrd(rec,6)
_pmin=subwrd(rec,8)
_zmin=subwrd(rec,11)
_zmax=subwrd(rec,13)
rec=sublin(result,5)
_ttype=subwrd(rec,3)
_tval=subwrd(rec,9)

"q file"
rec=sublin(result,5)
_zmaxfile=subwrd(rec,9)

*-------------------------------------------------------------
* Check to ensure that dimensions are valid.  Warn & exit if not.
*--------------------------------------------------------------

dimrc=0
If (_xtype != "fixed")
  say "X-Dims Error:  Not fixed.  Use 'set lon' or 'set x' to specify a value."
  dimrc=-1
Endif

If (_ytype != "fixed")
  say "Y-Dims Error:  Not fixed.  Use 'set lat' or 'set y' to specify a value"
  dimrc=-1
Endif

If (_ptype != "varying")
   say "Z-Dims Error:  Not varying.  Use 'set lev' or 'set z' to specify a range."
   dimrc=-1
Endif

If (_ttype != "fixed")
  say "Time Error:     Not fixed.  Use 'set time' or 'set t' to specify a value"
  dimrc=-1
Endif


If (dimrc < 0)
  Return(-1)
Endif


*
* A few global variables used in units conversion
*

_pi=3.14159265
_dtr=_pi/180
_rtd=1/_dtr 
_ktm=0.514444
_mtk=1/_ktm

* A few global constants used in thermo calcs

_C0=0.99999683 
_C1=-0.90826951/100
_C2= 0.78736169/10000
_C3=-0.61117958/1000000
_C4= 0.43884187/pow(10,8) 
_C5=-0.29883885/pow(10,10)
_C6= 0.21874425/pow(10,12)
_C7=-0.17892321/pow(10,14)
_C8= 0.11112018/pow(10,16)          
_C9=-0.30994571/pow(10,19)

* A pressure array of power calculations which should be performed
* only once to reduce execution time.

zz=1100
while (zz > 10)
    subscr=zz/10
    _powpres.subscr=pow(zz,0.286)
    zz=zz-10
endwhile

*
* Turn off options not available due to user data limitations
*

If (ClrScrn = 1) 
  "clear"
Endif

If (sndspd = -1 | snddir = -1) 
  DrawBarb = 0
  DrawHodo = 0
  DrawHeli = 0
Endif

If (snddewp = -1) 
  DrawDSnd = 0
  DrawRH   = 0
  DrawPrcl = 0
  DrawPMax = 0
  DrawIndx = 0
Endif

If (sndtemp = -1)
  DrawTSnd = 0
  DrawRH   = 0
  DrawPrcl = 0
  DrawPMax = 0
  DrawIndx = 0
  DrawZLev = 0
Endif

If (NumRing < 1) 
  DrawHodo = 0
Endif 
  
"q gxinfo"
rec=sublin(result,2)
xsize=subwrd(rec,4)

If (xsize = 11) 
   PageType = "Landscape"
Else
   PageType = "Portrait"
Endif
  
*------------------------------------------------------
* calculate constants determining slope/shape of diagram
* based on temp/pressure values given by user
*-------------------------------------------------------

"set x 1"
"set y 1"
"set z 1"
"set t 1"
_m1=(T1+T2-2*T3)/(2*log10(P2/P3))
_m2=(T2-T3-_m1*log10(P2/P3))/50
_m3=(T1-_m1*log10(P1))

"set z "_zmin" "_zmax           
"set zlog on"
"set xlab off"

*-------------------------------------------------
* perform coordinate transformation to Skew-T/LogP
*-------------------------------------------------

"set gxout stat"
"set x "_xval
"set y "_yval
"set t "_tval
"define tempx=("sndtemp"-"_m1"*log10(lev)-"_m3")/"_m2
"define dewpx=("snddewp"-"_m1"*log10(lev)-"_m3")/"_m2

If (PageType = "Portrait") 
   "set parea 0.7 7 0.75 10.5"
Else
   "set parea 0.7 6.5 0.5 8"
Endif
  
"set axlim 0 100"
"set lon 0 100"
"set grid on 1 1"

"set z "_zmin " " _zmax
"set lon 0 100"
"set clevs -900"
"set gxout contour"

*-------------------------------------
* Draw pressure lines 
*-------------------------------------

If (DrawPLev = 0) 
   "set ylab off"
Else
   "set ylab on"
   "set ylopts 1 3 0.10"
   "set xlopts 1 3 0.125"
Endif

"d lon"

*--------------------------------------
* Determine corners of skewt/logp frame
*--------------------------------------

"q w2xy 100 "_pmin
rxloc=subwrd(result,3)
tyloc=subwrd(result,6)
"q w2xy 0 "_pmax
lxloc=subwrd(result,3)
byloc=subwrd(result,6)

If (DrawPLev = 1 & LblAxes = 1)
   "set strsiz 0.10"
   "set string 1 c 3 0"
   If (PageType = "Portrait") 
      "draw string 0.5 10.85 mb"
   Else
      "draw string 0.5 8.35 mb"
   Endif
Endif

*---------------------------------------------------
* Calculate & draw actual height lines using temp data
*---------------------------------------------------

If (DrawZLev > 0)
   say "Calculating observed height levels from temp/pressure data."
   zz=1
   "set gxout stat"
   "set x "_xval
   "set y "_yval
   "set t "_tval
   count=0
   while (zz < _zmax)
      "set z "zz
      pp.zz=subwrd(result,4)
      lpp.zz=log(pp.zz)
      "d "sndtemp
      rec=sublin(result,8)
      tt=subwrd(rec,4)
      if (tt > -900) 
         tk=tt+273.15
         count=count+1
         zzm=zz-1
         If (count = 1) 
            If (DrawZLev = 2)
               htlb="ASL"
               height.zz=SfcElev
            Else
               htlb="AGL"
               height.zz=0
            Endif
            sfcz=height.zz
         Else
            DZ=29.2857*(lpp.zzm-lpp.zz)*(lpp.zz*tk+lpp.zzm*tkold)/(lpp.zz+lpp.zzm) 
            height.zz=height.zzm+DZ
            highz=height.zz
         Endif
      else
         height.zz = -9999
      endif
      tkold=tk
      zz=zz+1
   endwhile

   maxht=int(highz/1000)
   if (int(sfcz/1000) = sfcz/1000)
      minht=int(sfcz/1000)
   else
      minht=1+int(sfcz/1000)
   endif

   ht=minht
   "set line 1 3 1"
   "set strsiz 0.10"
   "set string 1 l 3 0"
   while (ht <= maxht)
       zz=1
       while (height.zz/1000 <= ht)
          zz=zz+1
       endwhile
       zzm=zz-1
       PBelow=pp.zzm
       PAbove=pp.zz
       HBelow=height.zzm
       HAbove=height.zz
       DZ=HAbove-HBelow
       DP=PAbove-PBelow
       Del=ht*1000-HBelow
       Est=PBelow+Del*DP/DZ
       If (Est >= _pmin & Est <= _pmax)
          "q w2xy 1 " Est
          yloc=subwrd(result,6)
          "draw line " lxloc " " yloc " " rxloc " " yloc
          "draw string 0.22 "yloc-0.05" "ht
       Endif
       ht=ht+1
   endwhile
   "set strsiz 0.10"
   "set string 1"
   If (LblAxes = 1)
      If (PageType = "Portrait") 
         "draw string 0.25 10.85 km"
         "draw string 0.25 10.75 "htlb
         "draw string 0.25 10.65 OBS"
      Else
         "draw string 0.25 8.35 km"
         "draw string 0.25 8.25 "htlb
         "draw string 0.25 8.15 OBS"
      Endif
   Endif
Endif


*---------------------------------------------------
* Draw height levels (height above MSL using Std Atm)
*---------------------------------------------------

If (DrawZSTD = 1)
   "set strsiz 0.10"
   minht=30.735*(1-pow(_pmax/1013.26,0.287))
   minht=int(minht+0.5)
   maxht=30.735*(1-pow(_pmin/1013.26,0.287))
   maxht=int(maxht)
   "set gxout stat"
   zcount=minht        
   while (zcount <= maxht) 
      plev=1013.26*pow((1-zcount/30.735),3.4843)
      "q w2xy 0 "plev 
      yloc=subwrd(result,6)
      "draw string 0 "yloc-0.05" "zcount
      zcount=zcount+1
   endwhile
   "set strsiz 0.10"
   If (LblAxes = 1)
      If (PageType = "Portrait") 
         "draw string 0 10.85 km"
         "draw string 0 10.75 ASL"
         "draw string 0 10.65 STD"
      Else
         "draw string 0 8.35 km"
         "draw string 0 8.25 ASL"
         "draw string 0 8.15 STD"
      Endif
  Endif
Endif


*-----------------------
* Plot temperature lines 
*-----------------------

If (DrawTemp = 1)
   "set strsiz 0.1"
   "set z "_zmin " " _zmax
   "set line "TempCol " " TempLine " "TempThk
   "set string 1 c 3 0"
   "set gxout stat"
   maxtline=GetTemp(100,_pmax)
   mintline=GetTemp(0,_pmin)

   maxtline=tempint*int(maxtline/tempint)
   mintline=tempint*int(mintline/tempint)

   tloop=mintline
   While (tloop <= maxtline) 
       Botxtemp=GetXLoc(tloop,_pmax)
       "q w2xy "Botxtemp " " _pmax
       Botxloc=subwrd(result,3)
       Botyloc=byloc           
       Topxtemp=GetXLoc(tloop,_pmin)
        "q w2xy "Topxtemp " " _pmin
       Topxloc=subwrd(result,3)
       Topyloc=tyloc            
       If (Botxtemp <= 100 | Topxtemp <= 100) 
          If (Topxtemp > 100)
             Slope=(Topyloc-Botyloc)/(Topxtemp-Botxtemp)
             b=Topyloc-Slope*Topxtemp
             Topyloc=Slope*100+b
             Topxloc=rxloc         
          Endif
          If (Botxtemp < 0)
             Slope=(Topyloc-Botyloc)/(Topxtemp-Botxtemp)
             b=Botyloc-Slope*Botxtemp
             Botyloc=b
             Botxloc=lxloc 
          Else
             "draw string " Botxloc-0.05 " " Botyloc-0.15 " " tloop
          Endif
          "draw line "Botxloc " " Botyloc " " Topxloc " " Topyloc
       Endif
       tloop=tloop+tempint
   EndWhile
   If (LblAxes = 1) 
      "set strsiz 0.15"
      "set string 1 c"
      If (PageType = "Portrait")
         "draw string 4.0 0.35 Temperature (`3.`0C)"
      Else
         "draw string 3.5 0.15 Temperature (`3.`0C)"
      Endif
   Endif
Endif


*------------------
* Plot dry adiabats
*------------------

If (DrawThet = 1)
   temp=GetTemp(100,_pmin)
   maxtheta=GetThet2(temp,-100,_pmin)
   maxtheta=thetaint*int(maxtheta/thetaint)
   temp=GetTemp(0,_pmax)
   mintheta=GetThet2(temp,-100,_pmax)
   mintheta=thetaint*int(mintheta/thetaint)
   
   "set lon 0 100"
   "set y 1"
   "set z 1"
   tloop=mintheta
   "set line "ThetCol" "ThetLine " "ThetThk
   While (tloop <= maxtheta)
     PTemp=LiftDry(tloop,1000,_pmin,1,_pmin,_pmax)     
     tloop=tloop+thetaint
   Endwhile
Endif

*------------------------
* Plot mixing ratio lines
*------------------------

If (DrawMix = 1)
   If (MixStop < _pmin) 
      MixStop = _pmin
   Endif
   "set string 1 l"
   "set z "_zmin " " _zmax
   "set cint 1"
   "set line "MixCol" " MixLine " "MixThk
   cont = 1
   mloop=subwrd(wsclevs,1)
   count = 1
   While (cont = 1) 
       BotCoef=log(mloop*_pmax/3801.66)
       BotTval=-245.5*BotCoef/(BotCoef-17.67)
       Botxtemp=GetXLoc(BotTval,_pmax)
       "q w2xy "Botxtemp " " _pmax
       Botxloc=subwrd(result,3)
       Botyloc=byloc            
       TopCoef=log(mloop*MixStop/3801.66)
       TopTval=-245.5*TopCoef/(TopCoef-17.67)
       Topxtemp=GetXLoc(TopTval,MixStop)
       "q w2xy "Topxtemp " " MixStop
       Topxloc=subwrd(result,3)
       Topyloc=subwrd(result,6) 
       "set string "MixCol" l 3"
       "set strsiz 0.09"
       If (Botxtemp <= 100 | Topxtemp <= 100) 
          If (Topxtemp > 100)
             Slope=(Topyloc-Botyloc)/(Topxtemp-Botxtemp)
             b=Topyloc-Slope*Topxtemp
             Topyloc=Slope*100+b
             Topxloc=rxloc 
             "draw string " Topxloc+0.05 " " Topyloc  " " mloop
          Else
             "draw string " Topxloc " " Topyloc+0.1 " " mloop
          Endif
          If (Botxtemp < 0)
             Slope=(Topyloc-Botyloc)/(Topxtemp-Botxtemp)
             b=Botyloc-Slope*Botxtemp
             Botyloc=b
             Botxloc=lxloc        
          Endif
          "draw line "Botxloc " " Botyloc " " Topxloc " " Topyloc
       Endif
       count=count+1
       mloop=subwrd(wsclevs,count)
       If (mloop = "" | count > 50) 
          cont = 0
       Endif
   EndWhile
   If (LblAxes = 1)
      "set strsiz 0.15"
      "set string 1 c 3 90"
      If (PageType = "Portrait")
         "draw string 7.40 4.75 Mixing Ratio (g/kg)"
      Else
         "draw string 6.90 4.25 Mixing Ratio (g/kg)"
      Endif
      "set string 1 c 3 0"
   Endif
Endif

*-----------------------------
* Plot moist (pseudo) adiabats
*-----------------------------

If (DrawThtw = 1)
   "set lon 0 100"
   "set y 1"
   "set z 1"
   "set gxout stat"
   tloop=80
   "set line "ThtwCol" "ThtwLine " "ThtwThk
   While (tloop > -80)
     PTemp=LiftWet(tloop,1000,ThtwStop,1,_pmin,_pmax)     
     tloop=tloop-thetwint
   Endwhile
Endif


*-----------------------------------------------------
* Plot transformed user-specified temperature sounding
*-----------------------------------------------------

If (DrawTSnd = 1)
   say "Drawing temperature sounding."
   "set gxout line"
   "set x "_xval
   "set y "_yval
   "set z "_zmin" "_zmax     
   "set ccolor "TSndCol
   "set cstyle "TSndLine
   "set cmark "TSndMrk
   "set digsiz "MrkSize
   "set cthick "TSndThk 
   "set missconn on"
   "d tempx"
Endif

*---------------------------------------------------
* Plot transformed user-specified dewpoint sounding
*---------------------------------------------------

If (DrawDSnd = 1)
   say "Drawing dewpoint sounding."
   "set gxout line"
   "set x "_xval
   "set y "_yval
   "set z "_zmin" "_zmax
   "set cmark "DSndMrk 
   "set digsiz "MrkSize
   "set ccolor "DSndCol
   "set cstyle "DSndLine
   "set cthick "DSndThk
   "set missconn on"
   "d dewpx"
Endif

*----------------------------------------
* Determine lowest level of reported  data
*----------------------------------------

If (DrawTSnd = 1) 
   field=sndtemp
Else
   field=sndspd
Endif

"set gxout stat"
"set x "_xval
"set y "_yval
"set t "_tval
"set lev " _pmax " " _pmin
"d maskout(lev,"field"+300)"
rec=sublin(result,1)
check=substr(rec,1,6)
If (check = "Notice") 
    rec=sublin(result,9)
Else
    rec=sublin(result,8)
Endif
SfcPlev=subwrd(rec,5)

If (DrawTSnd = 1 & DrawDSnd = 1)
   "set lev "SfcPlev
   "d "sndtemp
   rec=sublin(result,8)
   Sfctemp=subwrd(rec,4)
   "d "snddewp
   rec=sublin(result,8)
   Sfcdewp=subwrd(rec,4)
   SfcThee=Thetae(Sfctemp,Sfcdewp,SfcPlev)

*------------------------------------------
* Calculate temperature and pressure of LCL
*------------------------------------------

   TLcl=Templcl(Sfctemp,Sfcdewp)
   PLcl=Preslcl(Sfctemp,Sfcdewp,SfcPlev)
Endif

*----------------------------------------------------------
* Plot parcel path from surface to LCL and up moist adiabat
*----------------------------------------------------------

If (DrawPrcl = 1)
   say "Drawing parcel path from surface upward."
   If (PageType = "Portrait")
      xloc=7.15
   Else
      xloc=6.65
   Endif
   "q w2xy 1 "PLcl
   rec=sublin(result,1)
   yloc=subwrd(rec,6)
   "set strsiz 0.1"
   If (PLcl < _pmax) 
      "set string 1 l"
      "draw string "xloc" "yloc" LCL"
      "set line 1 1 1"
      "draw line "xloc-0.15" "yloc" "xloc-0.05" "yloc
   Endif
   "set lon 0 100"
   "set gxout stat"
   "set line "PrclCol" "PrclLine " " PrclThk
   PTemp=LiftDry(Sfctemp,SfcPlev,PLcl,1,_pmin,_pmax)
   Ptemp=LiftWet(TLcl,PLcl,_pmin,1,_pmin,_pmax)
Endif

*-------------------------------------------------------
* Determine level within lowest 250mb having highest
* theta-e value
*-------------------------------------------------------

If (DrawTSnd = 1 & DrawDSnd = 1)
  "set x "_xval
  "set y "_yval
  "set t "_tval
   zz=1
   MaxThee=-999
   "set gxout stat"
   while (zz <= _zmax & pp > _pmax-250)
       "set z "zz
       pp=subwrd(result,4)
       "d "sndtemp
       rec=sublin(result,8)
       tt=subwrd(rec,4)
       "d "snddewp
       rec=sublin(result,8)
       dd=subwrd(rec,4)
       If (abs(tt) < 130 & abs(dd) < 130) 
          Thee=Thetae(tt,dd,pp)
          If (Thee > MaxThee) 
             MaxThee=Thee
             TMaxThee=tt
             DMaxThee=dd 
             PMaxThee=pp
          Endif
       endif
       zz=zz+1
   Endwhile
   If (PMaxThee = SfcPlev-250) 
      PMaxThee = SfcPlev
   Endif
*------------------------------------------------------
* Calculate temperature and pressure of LCL from highest   
* theta-e level
*------------------------------------------------------
   If (SfcPlev != PMaxThee) 
      TLclMax=Templcl(TMaxThee,DMaxThee)
      PLclMax=Preslcl(TMaxThee,DMaxThee,PMaxThee)
   Endif
Endif

*----------------------------------------------------------
* Plot parcel path from highest theta-e level to LCL and up
* moist adiabat
*----------------------------------------------------------

If (DrawPMax = 1 & SfcPlev != PMaxThee)
   say "Drawing parcel path from most unstable level upward."
   If (PageType = "Portrait")
      xloc=7.15
   Else
      xloc=6.65
   Endif
   "q w2xy 1 "PLclMax
   rec=sublin(result,1)
   yloc=subwrd(rec,6)
   "set strsiz 0.1"
   If (PLclMax < _pmax) 
      "set string 1 l"
      "draw string "xloc" "yloc" LCL"
      "set line 1 1 1"
      "draw line "xloc-0.15" "yloc" "xloc-0.05" "yloc
   Endif
   "set lon 0 100"
   "set gxout stat"
   "set line "PrclCol" "PrclLine " " PrclThk
   PTemp=LiftDry(TMaxThee,PMaxThee,PLclMax,1,_pmin,_pmax)
   Ptemp=LiftWet(TLclMax,PLclMax,_pmin,1,_pmin,_pmax)
Endif

*--------------------------------
* Draw thermodynamic indices
*--------------------------------

If (DrawIndx = 1) 
   "set string 1 l"
   "set strsiz 0.10"
   "set x "_xval
   "set y "_yval
   "set t "_tval
   say "Calculating precipitable water."
   pw=precipw(sndtemp,snddewp,_pmax,_pmin)
   say "Calculating thermodynamic indices."
   Temp850=interp(sndtemp,850)
   Temp700=interp(sndtemp,700)
   Temp500=interp(sndtemp,500)
   Dewp850=interp(snddewp,850)
   Dewp700=interp(snddewp,700) 
   Dewp500=interp(snddewp,500)
   If (Temp850>-900 & Dewp850>-900 & Dewp700>-900 & Temp700>-900 & Temp500>-900) 
      K=Temp850+Dewp850+Dewp700-Temp700-Temp500
   Else 
      K=-999
   Endif
   If (Temp850 > -900 & Dewp850 > -900 & Temp500 > -900)
      tt=Temp850+Dewp850-2*Temp500
   Else
      tt=-999
   Endif
   Temp500V=virtual2(Temp500+273.15,Dewp500+273.15,500)-273.15
   PclTemp=LiftWet(TLcl,PLcl,500,0)
   PclTempV=virtual2(PclTemp+273.15,PclTemp+273.15,500)-273.15
   SLI=Temp500V-PclTempV
   rec=CAPE(TLcl,PLcl,100,sndtemp,snddewp)
   Pos=subwrd(rec,1)
   CIN=subwrd(rec,2)

   If (SfcPlev != PMaxThee) 
      PclTemp=LiftWet(TLclMax,PLclMax,500,0)
      PclTempV=virtual2(PclTemp+273.15,PclTemp+273.15,500)-273.15
      LIMax=Temp500V-PclTempV
      rec=CAPE(TLclMax,PLclMax,100,sndtemp,snddewp)
      PosMax=subwrd(rec,1)
      CINMax=subwrd(rec,2)
   Else
      LIMax=SLI
      PosMax=Pos
      CINMax=CIN
      MaxThee=SfcThee
   Endif
 
   If (PageType = "Portrait") 
      If (Text1XC = -1)
         Text1XC=rxloc-0.75
      Endif
      If (Text1YC = -1)
         Text1YC=tyloc-2.25
      Endif
      If (Text2XC = -1)
         Text2XC=rxloc-0.75
      Endif
      If (Text2YC = -1)
         Text2YC=tyloc-3.25
      Endif
      If (Text3XC = -1)
          Text3XC=rxloc-0.75
      Endif
      If (Text3YC = -1)
         Text3YC=tyloc-4.40
      Endif
   Else
      If (Text1XC = -1)
         Text1XC=rxloc+2.50
      Endif
      If (Text1YC = -1)
         Text1YC=tyloc-3.00
      Endif
      If (Text2XC = -1)
         Text2XC=rxloc+2.50
      Endif
      If (Text2YC = -1)
         Text2YC=tyloc-4.00
      Endif
      If (Text3XC = -1)
         Text3XC=rxloc+2.50
      Endif
      If (Text3YC = -1)
         Text3YC=tyloc-5.10
      Endif
   Endif
   "set string 1 l 3"
   "set line 0 1 3"
   "draw recf  "Text1XC-0.75 " " Text1YC-0.40 " " Text1XC+0.75 " " Text1YC+0.25
   "set line 1 1 3"
   "draw rec  "Text1XC-0.75 " " Text1YC-0.40 " " Text1XC+0.75 " " Text1YC+0.25
   "draw string "Text1XC-0.70 " " Text1YC+0.10"  K" 
   "draw string "Text1XC+0.25 " " Text1YC+0.10" " int(K)      
   "draw string "Text1XC-0.70 " " Text1YC-0.10 "  TT" 
   "draw string "Text1XC+0.25 " " Text1YC-0.10 " " int(tt)
   "draw string "Text1XC-0.70 " " Text1YC-0.25 "  PW(cm)" 
   "draw string "Text1XC+0.25 " " Text1YC-0.25 " " int(pw*100)/100
   "set line 0 1 3"
   "draw recf  "Text2XC-0.75 " " Text2YC-0.60 " " Text2XC+0.75 " " Text2YC+0.60
   "set line 1 1 3"
   "draw rec  "Text2XC-0.75 " " Text2YC-0.60 " " Text2XC+0.75 " " Text2YC+0.60
   "draw string "Text2XC-0.35 " " Text2YC+0.50 " Surface"
   "draw string "Text2XC-0.70 " " Text2YC+0.30 "  Temp(`3.`0C)" 
   "draw string "Text2XC+0.25 " " Text2YC+0.30 " " int(Sfctemp*10)/10
   "draw string "Text2XC-0.70 " " Text2YC+0.15 "  Dewp(`3.`0C)" 
   "draw string "Text2XC+0.25 " " Text2YC+0.15 " " int(Sfcdewp*10)/10
   "draw string "Text2XC-0.70 " " Text2YC "   `3z`0`bE`n(K)"
   "draw string "Text2XC+0.25 " " Text2YC " " int(SfcThee) 
   "draw string "Text2XC-0.70 " " Text2YC-0.15 "  LI"
   "draw string "Text2XC+0.25 " " Text2YC-0.15 " " round(SLI)
   "draw string "Text2XC-0.70 " " Text2YC-0.30 "  CAPE(J)"
   "draw string "Text2XC+0.25 " " Text2YC-0.30 " " int(Pos)   
   "draw string "Text2XC-0.70 " " Text2YC-0.45 "  CIN(J)"
   "draw string "Text2XC+0.25 " " Text2YC-0.45 " " int(CIN)      
   "set line 0 1 3"
   "draw recf  "Text3XC-0.75 " " Text3YC-0.55 " "  Text3XC+0.75 " " Text3YC+0.55
   "set line 1 1 3"
   "draw rec  "Text3XC-0.75 " " Text3YC-0.55 " "  Text3XC+0.75 " " Text3YC+0.55
   "draw string "Text3XC-0.60 " " Text3YC+0.45 " Most Unstable"
   "draw string "Text3XC-0.70 " " Text3YC+0.20 "  Press(mb)" 
   "draw string "Text3XC+0.25 " " Text3YC+0.20 " " int(PMaxThee)
   "draw string "Text3XC-0.70 " " Text3YC+0.05 " `3z`0`bE`n(K)"
   "draw string "Text3XC+0.25 " " Text3YC+0.05 " " int(MaxThee)
   "draw string "Text3XC-0.70 " " Text3YC-0.10 " LI" 
   "draw string "Text3XC+0.25 " " Text3YC-0.10 " "round(LIMax)
   "draw string "Text3XC-0.70 " " Text3YC-0.25 " CAPE(J)" 
   "draw string "Text3XC+0.25 " " Text3YC-0.25 " "int(PosMax) 
   "draw string "Text3XC-0.70 " " Text3YC-0.40 " CIN(J)"
   "draw string "Text3XC+0.25 " " Text3YC-0.40 " " int(CINMax) 
Endif

*-----------------------------
* Draw wind profile along side
*-----------------------------

If (DrawBarb = 1) 
   say "Drawing Wind Profile."
   If (poleloc = -1) 
      If (PageType = "Portrait")
         poleloc = 8.0
      Else
         poleloc = 7.5 
      Endif
   Endif
   If (barbline = 1)
      "set line 1 1 3"
      "draw line "poleloc " " byloc " " poleloc " " tyloc
   Endif
   If (BarbCol = -1) 
      'set rgb 41 255 0 132'
      'set rgb 42 255 0 168'
      'set rgb 43 255 0 204'
      'set rgb 44 255 0 240'
      'set rgb 45 255 0 255'
      'set rgb 46 204 0 255'
      'set rgb 47 174 0 255'
      'set rgb 48 138 0 255'
      'set rgb 49 108 0 255'
      'set rgb 50 84 0 255'
      'set rgb 51 40 0 255'
      'set rgb 52 0 0  255'
      'set rgb 53 0 42 255'
      'set rgb 54 0 84 255'
      'set rgb 55 0 120 255'
      'set rgb 56 0 150 255'
      'set rgb 57 0 192 255'
      'set rgb 58 0 240 255'
      'set rgb 59 0 255 210'
      'set rgb 60 0 255 160' 
      'set rgb 61 0 255 126'
      'set rgb 62 0 255 78'
      'set rgb 63 84 255 0'
      'set rgb 64 138 255 0'
      'set rgb 65 188 255 0'
      'set rgb 66 236 255 0'
      'set rgb 67 255 255 0'
      'set rgb 68 255 222 0'
      'set rgb 69 255 192 0'
      'set rgb 70 255 162 0'
      'set rgb 71 255 138 0'
      'set rgb 72 255 108 0'
      'set rgb 73 255 84 0'
      'set rgb 74 255 54 0'
      'set rgb 75 255 12 0'
      'set rgb 76 255 0 34'
      'set rgb 77 255 0 70'
      'set rgb 78 255 0 105'
      'set rgb 79 255 0 140'
      'set rgb 80 255 0 175'
      'set rgb 81 255 0 215'
      'set rgb 82 255 0 255'
      'set rgb 83 255 255 255'

      col1='83 83 83 83 83 83 83 83 83 83 82 81 80 79 78'
      col2='77 76 75 74 73 72 71 70 69 68 67 66 65 64 63'
      col3='62 61 60 59 58 57 56 55 54 53 52 51 50 49 48'
      'set rbcols 'col1' 'col2' 'col3
   Endif
   "set z "_zmin" "_zmax
   "set gxout stat"
   zz=1
   wspd=-999
   cont=1
   While (cont = 1 & zz < _zmax) 
      "set z "zz
      pres=subwrd(result,4)
      "d "sndspd
      rec=sublin(result,8)
      wspd=subwrd(rec,4)
      if (wspd < 0 | pres > _pmax) 
          zz=zz+1
      else
          cont=0
      Endif
   Endwhile
   While (zz <= _zmax)
      "d "sndspd"(z="zz")"
      rec=sublin(result,8)
      wspd=subwrd(rec,4)
      If (BarbCol >= 0)
         "set line "BarbCol " 1 "BarbThk
      Else
         tempbcol=55+wspd/5     
         If (tempbcol > 83) 
            tempbcol = 83
         Endif
         "set line "tempbcol " 1 "BarbThk
      Endif
      "d "snddir"(z="zz")"
      rec=sublin(result,8)
      wdir=subwrd(rec,4)
      xwind=GetUWnd(wspd,wdir)
      ywind=GetVWnd(wspd,wdir)
      "query gr2xy 5 "zz
      y1=subwrd(result,6) 
      if (wspd > 0) 
         cc=polelen/wspd
         xendpole=poleloc-xwind*cc
         yendpole=y1-ywind*cc
      endif
      if (xendpole>0 & wspd >= 0.5)
        if (flagbase = 1) 
           "draw mark 3 "poleloc " " y1 " 0.05"
        endif
        "draw line " poleloc " " y1 "  " xendpole " " yendpole
        flagloop=wspd/10
        windcount=wspd
        flagcount=0
        xflagstart=xendpole
        yflagstart=yendpole
        dx=cos((180-wdir)*_dtr)
        dy=sin((180-wdir)*_dtr)
        while (windcount > 47.5)
           flagcount=flagcount+1
           dxflag=-len50*dx
           dyflag=-len50*dy
           xflagend=xflagstart+dxflag
           yflagend=yflagstart+dyflag
           windcount=windcount-50
           x1=xflagstart+0.5*wid50*xwind/wspd
           y1=yflagstart+0.5*wid50*ywind/wspd
           x2=xflagstart-0.5*wid50*xwind/wspd
           y2=yflagstart-0.5*wid50*ywind/wspd
           If (Fill50 = 1) 
              "draw polyf "x1" "y1" "x2" "y2" "xflagend" "yflagend" "x1" "y1
           Else
              "draw line "x1 " "y1 " " xflagend " " yflagend " "  
              "draw line "x2 " "y2 " " xflagend " " yflagend
              "draw line "x1 " "y1 " " x2 " " y2
           Endif
           xflagstart=xflagstart+spac50*xwind/wspd
           yflagstart=yflagstart+spac50*ywind/wspd
        endwhile
        while (windcount > 7.5 ) 
           flagcount=flagcount+1
           dxflag=-len10*dx
           dyflag=-len10*dy
           xflagend=xflagstart+dxflag
           yflagend=yflagstart+dyflag
           windcount=windcount-10
           "draw line " xflagstart " " yflagstart " " xflagend " " yflagend
           xflagstart=xflagstart+spac10*xwind/wspd
           yflagstart=yflagstart+spac10*ywind/wspd
        endwhile
        if (windcount > 2.5) 
           flagcount=flagcount+1
           if (flagcount = 1) 
              xflagstart=xflagstart+spac05*xwind/wspd
              yflagstart=yflagstart+spac05*ywind/wspd
           endif
           dxflag=-len05*dx
           dyflag=-len05*dy
           xflagend=xflagstart+dxflag
           yflagend=yflagstart+dyflag
           windcount=windcount-5
           "draw line " xflagstart " " yflagstart " " xflagend " " yflagend
        endif
      else
        if (wspd < 0.5 & wspd >= 0) 
           "draw mark 2 " poleloc " " y1 " 0.08"
        endif
      endif
      zz=zz+barbint
   endwhile
Endif

*----------------
* Draw Hodograph
*----------------

If (DrawHodo = 1)
   say "Drawing Hodograph."

   If (HodXcent = -1 | HodYcent = -1) 
      If (PageType = "Portrait") 
         HodXcent=6
         HodYcent=9.5
      Else
         HodXcent=9
         HodYcent=7.0
      Endif
   Endif
   HodL=HodXcent-HodSize/2.0
   HodR=HodXcent+HodSize/2.0
   HodT=HodYcent+HodSize/2.0
   HodB=HodYcent-HodSize/2.0
   RingSpac=HodSize/(NumRing*2)
   "set line 0"
   "draw recf "HodL" "HodB" "HodR" "HodT
   "set line "HodoCol" 1 6"
   "draw rec "HodL" "HodB" "HodR" "HodT
   "set line 1 1 3"
   "set string 1 c"
   "draw mark 1 "HodXcent " "HodYcent " " HodSize
   i=1
   While (i <= NumRing)
     "set strsiz 0.10"
     "set string 1 c 3 45"
     uwnd=-i*HodRing*cos(45*_dtr)
     xloc=HodXcent+uwnd*RingSpac/HodRing
     yloc=HodYcent+uwnd*RingSpac/HodRing
  
     "draw mark 2 "HodXcent " " HodYcent " " i*HodSize/NumRing
     "draw string "xloc " " yloc " " HodRing*i
     i=i+1
   Endwhile
   "set string 1 l 3 0"
   If (TickInt > 0) 
      i=0
      while (i < HodRing*NumRing) 
         dist=i*RingSpac/HodRing
         hrxloc=HodXcent+dist                
         hlxloc=HodXcent-dist                      
         htyloc=HodYcent+dist
         hbyloc=HodYcent-dist
         "set line 1 1 3"
         "draw line "hrxloc " " HodYcent-TickSize/2 " " hrxloc " " HodYcent+TickSize/2
         "draw line "hlxloc " " HodYcent-TickSize/2 " " hlxloc " " HodYcent+TickSize/2
         "draw line "HodXcent+TickSize/2 " " htyloc " " HodXcent-TickSize/2 " " htyloc
         "draw line "HodXcent+TickSize/2 " " hbyloc " " HodXcent-TickSize/2 " " hbyloc
         i=i+TickInt
      endwhile
   Endif
   "set line "HodoCol " " HodoLine " "HodoThk
   "draw string "HodL+0.05 " " HodT-0.1 " knots"
   zloop=_zmin
   xold=-999
   yold=-999
   count=0
   Depth=0
   While (zloop < _zmax & Depth < HodoDep)
      "set z "zloop
      pres=subwrd(result,4)
      "d "sndspd
      rec=sublin(result,8)
      wspd=subwrd(rec,4)
      "d "snddir
      rec=sublin(result,8)
      wdir=subwrd(rec,4)
      uwnd=GetUWnd(wspd,wdir)
      vwnd=GetVWnd(wspd,wdir)
      If (wspd >= 0) 
         xloc=HodXcent+uwnd*RingSpac/HodRing
         yloc=HodYcent+vwnd*RingSpac/HodRing
         If (xloc > 0 & yloc > 0 & xold > 0 & yold > 0) 
            Depth=Depth+pold-pres
            count=count+1
            If (count = 1) 
               "draw mark 3 "xold " " yold " 0.05"
            Endif
            "draw line "xold" "yold" "xloc" "yloc
         Endif
         xold=xloc
         yold=yloc
      Endif
      zloop=zloop+1
      pold=pres
   EndWhile

   If (count > 0) 
      "draw mark 3 "xold " " yold " 0.05"
   Endif
Endif

*----------------------------------------------
* Calculate and Display Absolute & S-R Helicity
*----------------------------------------------

If (DrawHeli = 1) 
   say "Calculating Helicity & SR Helicity."
   delp=10
   UTotal=0
   VTotal=0

* First, calculate mass-weighted mean wind 
* Since delp is a constant, and mass is proportional to 
* delp, this is a simple sum.

   "set lev "_pmax " " _pmin
   "define uwndarr="sndspd"*cos((270-"snddir")*"_dtr")"
   "define vwndarr="sndspd"*sin((270-"snddir")*"_dtr")"

   pres=MeanVBot
   While (pres >= MeanVTop)
      uwnd=interp(uwndarr,pres)*_ktm
      vwnd=interp(vwndarr,pres)*_ktm
      If (uwnd > -900 & vwnd > -900) 
         UTotal=UTotal+uwnd
         VTotal=VTotal+vwnd
      Endif
      pres=pres-delp
   EndWhile
   vcount=1+(MeanVBot-MeanVTop)/delp 
   Umean=UTotal/vcount
   Vmean=VTotal/vcount
   Spdmean=GetWSpd(Umean,Vmean)
   MeanDir=GetWDir(Umean,Vmean)

* Now, rotate and reduce mean wind to get storm motion

   If (StormMot = 1) 
      If (Spdmean < 15) 
         Reduct=0.25
         Rotate=30
      Else
         Reduct=0.20
         Rotate=20
      Endif
   Else
      Reduct=0.0
      Rotate=0.0
   Endif

   UReduce=(1-Reduct)*Umean
   VReduce=(1-Reduct)*Vmean
   StormSpd=GetWSpd(UReduce,VReduce)

   StormDir=GetWDir(UReduce,VReduce)+Rotate
   If (StormDir >= 360) 
      StormDir=StormDir-360
   Endif

   StormU=GetUWnd(StormSpd,StormDir)
   StormV=GetVWnd(StormSpd,StormDir)

* Draw Storm Motion Vector

   xloc=HodXcent+_mtk*StormU*RingSpac/HodRing
   yloc=HodYcent+_mtk*StormV*RingSpac/HodRing

   "set line 1 1 4"
   "draw line "HodXcent " " HodYcent " " xloc " " yloc
   Arr1U=GetUWnd(HodRing/10,StormDir+30)
   Arr1V=GetVWnd(HodRing/10,StormDir+30)
   Arr2U=GetUWnd(HodRing/10,StormDir-30)
   Arr2V=GetVWnd(HodRing/10,StormDir-30)

   xloc2=xloc-Arr1U/HodRing
   xloc3=xloc-Arr2U/HodRing
   yloc2=yloc-Arr1V/HodRing
   yloc3=yloc-Arr2V/HodRing

   "set line 1 1 3"

   If (FillArrw = 0) 
      "draw line "xloc" "yloc" "xloc2" "yloc2
      "draw line "xloc" "yloc" "xloc3" "yloc3
   Else
      "draw polyf "xloc" "yloc" "xloc2" "yloc2" "xloc3" "yloc3" "xloc" "yloc
   Endif

 
* Now, calculate SR and Environmental Helicity
 
   helic=0
   SRhelic=0
   MinP=SfcPlev-HelicDep
   pres=SfcPlev
   uwndold=-999
   vwndold=-999
   While (pres >= MinP)
      uwnd=interp(uwndarr,pres)*_ktm
      vwnd=interp(vwndarr,pres)*_ktm
      If (uwnd > -900 & uwndold > -900)
          du=uwnd-uwndold
          dv=vwnd-vwndold
          ubar=0.5*(uwnd+uwndold)
          vbar=0.5*(vwnd+vwndold)
          uhelic=-dv*ubar                   
          vhelic=du*vbar                   
          SRuhelic=-dv*(ubar-StormU)
          SRvhelic=du*(vbar-StormV)
          SRhelic=SRhelic+SRuhelic+SRvhelic
          helic=helic+uhelic+vhelic
      Endif
      uwndold=uwnd
      vwndold=vwnd
      pres=pres-delp
   EndWhile

   "set strsiz 0.1"
   "set string 1 l 3"
   If (PageType = "Portrait") 
      If (Text4XC = -1)
         Text4XC=rxloc-0.75
      Endif
      If (Text4YC = -1) 
         Text4YC=tyloc-5.45
      Endif
   Else
      If (Text4XC = -1)
         Text4XC=rxloc+2.50
      Endif
      If (Text4YC = -1)
         Text4YC=tyloc-6.10
      Endif
   Endif
   "set line 0 1 3"
   "draw recf  "Text4XC-0.75 " "Text4YC-0.5 " " Text4XC+0.75 " " Text4YC+0.5
   "set line 1 1 3"
   "draw rec  "Text4XC-0.75 " "Text4YC-0.5 " " Text4XC+0.75 " " Text4YC+0.5
   "draw string "Text4XC-0.45 " " Text4YC+0.40 " Hodograph"
   "draw string "Text4XC-0.70 " " Text4YC+0.20 " EH"
   "draw string "Text4XC+0.25 " " Text4YC+0.20 " "int(helic)
   "draw string "Text4XC-0.70 " " Text4YC+0.05 " SREH"
   "draw string "Text4XC+0.25 " " Text4YC+0.05 " " int(SRhelic)
   "draw string "Text4XC-0.70 " " Text4YC-0.20 " StmDir"
   "draw string "Text4XC+0.25 " " Text4YC-0.20 " " int(StormDir)"`3.`0"
   "draw string "Text4XC-0.70 " " Text4YC-0.35 " StmSpd(kt)"
   "draw string "Text4XC+0.25 " " Text4YC-0.35 " " int(_mtk*StormSpd)
Endif

*---------------------------------------------------
* Plot RH profile.
*---------------------------------------------------

If (DrawRH = 1)
  "set z "_zmin" "_zmax
  "set x "_xval
  "set y "_yval
  "set t "_tval
  "set zlog on"
  "set gxout line"
  "set ccolor "RHCol
  "set cstyle "RHLine
  "set cmark "RHMrk
  "set digsiz "MrkSize
  "set missconn on"
  "set xlab on"
  "set frame off"
  "set vrange 0 350"
  "set xlpos 0 t"
  "set xlevs 25 50 75 100"
  "set grid vertical 5"
  "define rh=100*exp((17.2694*"snddewp")/("snddewp"+237.3)-(17.2694*"sndtemp")/("sndtemp"+237.3))"
  "d rh"
   If (LblAxes = 1)
     "set string 1 c 3 0"
     "set strsiz 0.125"
     If (PageType = "Portrait")
       "draw string 1.5 10.85 RH (%)"
     Else
       "draw string 1.75 8.35 RH (%)"
     Endif
   Endif
Endif

*------------------------------------------
* Reset environment to original dimensions
*------------------------------------------

"set t "_tval
"set x "_xval 
"set y "_yval 
"set z "_zmin " "_zmax

say "Done."

Return(0)

*************************************************************************
function Templcl(temp,dewp)

*------------------------------------------------------
* Calculate the temp at the LCL given temp & dewp in C
*------------------------------------------------------

tempk=temp+273.15
dewpk=dewp+273.15
Parta=1/(dewpk-56)
Partb=log(tempk/dewpk)/800
Tlcl=1/(Parta+Partb)+56
return(Tlcl-273.15)

**************************************************************************

function Preslcl(temp,dewp,pres)

*-------------------------------------------------------
* Calculate press of LCL given temp & dewp in C and pressure
*-------------------------------------------------------

Tlclk=Templcl(temp,dewp)+273.15
tempk=temp+273.15
theta=tempk*pow(1000/pres,0.286)
plcl=1000*pow(Tlclk/theta,3.48)
return(plcl)

**************************************************************************
function LiftWet(startt,startp,endp,display,Pmin,Pmax)

*--------------------------------------------------------------------
* Lift a parcel moist adiabatically from startp to endp.
* Init temp is startt in C.  If you wish to see the parcel's
* path plotted, display should be 1.  Returns temp of parcel at endp.
*--------------------------------------------------------------------

temp=startt
pres=startp
cont = 1
delp=10
While (pres >= endp & cont = 1) 
    If (display = 1) 
       xtemp=GetXLoc(temp,pres)
       "q w2xy "xtemp" "pres
       xloc=subwrd(result,3)
       yloc=subwrd(result,6)
       If (xtemp < 0 | xtemp > 100)
          cont=0
       Else
          If (pres >= Pmin & pres < Pmax & pres < startp)  
             "draw line "xold" "yold" "xloc" "yloc 
          Endif
       Endif
    Endif
    xold=xloc
    yold=yloc
    temp=temp-100*delp*gammaw(temp,pres-delp/2,100)
    pres=pres-delp
EndWhile
return(temp)


**************************************************************************
function LiftDry(startt,startp,endp,display,Pmin,Pmax)

*--------------------------------------------------------------------
* Lift a parcel dry adiabatically from startp to endp.
* Init temp is startt in C.  If you wish to see the parcel's
* path plotted, display should be 1.  Returns temp of parcel at endp.
*--------------------------------------------------------------------

starttk=startt+273.15
cont = 1
delp=10
round=int(startp/10)*10
subscr=0.1*round
powstart=pow(startp,-0.286)
temp=starttk*_powpres.subscr*powstart-273.15
pres=round-10
While (pres >= endp & cont = 1) 
    subscr=0.1*pres
    temp=starttk*_powpres.subscr*powstart-273.15
    If (display = 1) 
       xtemp=GetXLoc(temp,pres)
       "q w2xy "xtemp" "pres
       xloc=subwrd(result,3)
       yloc=subwrd(result,6)
       If (xtempold > 0 & xtempold < 100 & xtemp > 0 & xtemp < 100) 
          If (pres >= Pmin & pres < Pmax & pres < startp)  
             "draw line "xold" "yold" "xloc" "yloc 
          Endif
       Endif
    Endif
    xold=xloc
    xtempold=xtemp
    yold=yloc
    pres=pres-delp
EndWhile
return(temp)

**************************************************************************
function CAPE(startt,startp,endp,sndtemp,snddewp)

*---------------------------------------------------------------------
* Returns all postive area and convective inhibition above LCL.
* Parcel is lifted from LCL at startt,startp and is halted
* at endp.   Integration method used is trapezoid method.
*---------------------------------------------------------------------

pres=startp
PclTemp=startt
PclTempV=virtual2(PclTemp+273.15,PclTemp+273.15,pres)-273.15
delp=10
Pos=0
Neg=0
Neg2=0

count=0
While (pres >= endp)
   EnvTemp=interp(sndtemp,pres)
   EnvDewp=interp(snddewp,pres)
   EnvTempV=virtual2(EnvTemp+273.15,EnvDewp+273.15,pres)-273.15
   DelT=PclTempV-EnvTempV
   If (abs(EnvTempV) < 130 & abs(PclTempV) < 130)
     count=count+1
     If (count > 1) 
       Val=DelT/pres+Prev 
       If (Val > 0)
          Pos=Pos+Val
          Neg2=0
       Else
          Neg=Neg+abs(Val)
          Neg2=Neg2+abs(Val)
       Endif
     Endif
     Prev=DelT/pres
   Endif
   pres=pres-delp
   PclTemp=PclTemp-100*delp*gammaw(PclTemp,pres,100)
   PclTempV=virtual2(PclTemp+273.15,PclTemp+273.15,pres)-273.15
Endwhile

Pos=0.5*Pos*287*delp
CIN=0.5*(Neg-Neg2)*287*delp

return(Pos" "CIN)

***************************************************************************
function gammaw(tempc,pres,rh)

*-----------------------------------------------------------------------
* Function to calculate the moist adiabatic lapse rate (deg C/Pa) based
* on the temperature, pressure, and rh of the environment.
*----------------------------------------------------------------------

tempk=tempc+273.15
es=satvap2(tempc)
ws=mixratio(es,pres)
w=rh*ws/100
tempv=virtual(tempk,w)
latent=latentc(tempc)

A=1.0+latent*ws/(287*tempk)
B=1.0+0.622*latent*latent*ws/(1005*287*tempk*tempk)
Density=100*pres/(287*tempv)
lapse=(A/B)/(1005*Density)
return(lapse)

*************************************************************************
function latentc(tempc)

*-----------------------------------------------------------------------
* Function to return the latent heat of condensation in J/kg given
* temperature in degrees Celsius.
*-----------------------------------------------------------------------

val=2502.2-2.43089*tempc

return(val*1000)

*************************************************************************
function precipw(sndtemp,snddewp,startp,endp)

*-----------------------------------------------------------------------
* Function to calculate the precipitable water (cm) in a sounding
* starting at pressure level startp and ending at pressure level endp.
*-----------------------------------------------------------------------

ppold=-999
ttold=-999
ddold=-999
delp=10
Int=0
mix=0 
pres=startp
logpp=log(pres)
logppm=log(pres-delp)
while (pres >= endp)
   tt=interp(sndtemp,pres)
   dd=interp(snddewp,pres)
   if (tt>-900 & ttold>-900 & dd>-900 & ddold>-900) 
      e=satvap2(dd)
      mix=mixratio(e,pres)
      mixavg=(logpp*mix+logppm*mixold)/(logpp+logppm)
      Int=Int+1.020408*mixavg*delp
   endif
   ttold=tt
   ddold=dd
   ppold=pp
   mixold=mix
   pres=pres-delp
   logpp=logppm
   logppm=log(pres-delp)
endwhile

return(Int)

*************************************************************************

function virtual(temp,mix)

*------------------------------------------------------------
* Function to return virtual temperature given temperature in 
* kelvin and mixing ratio in g/g.
*-------------------------------------------------------------

tempv=temp*(1.0+0.6*mix)

return (tempv)

************************************************************************

function virtual2(temp,dewp,pres)
  
*------------------------------------------------------------
* Function to return virtual temperature in kelvin given temperature in
* kelvin and dewpoint in kelvin and pressure in mb
*-------------------------------------------------------------

if (temp > 0 & dewp > 0) 
  vap=satvap2(dewp-273.15)
  mix=mixratio(vap,pres)
  tempv=virtual(temp,mix)
else
  tempv=-9999
endif

return (tempv)
  
************************************************************************

function satvapor(temp)

*---------------------------------------------------------------
* Given temp in Celsius, returns saturation vapor pressure in mb
*---------------------------------------------------------------

pol=_C0+temp*(_C1+temp*(_C2+temp*(_C3+temp*(_C4+temp*(_C5+temp*(_C6+temp*(_C7+temp*(_C8+temp*(_C9)))))))))

return(6.1078/pow(pol,8))

************************************************************************

function satvap2(temp)

*---------------------------------------------------------------
* Given temp in Celsius, returns saturation vapor pressure in mb
*---------------------------------------------------------------

es=6.112*exp(17.67*temp/(temp+243.5))

return(es)

*************************************************************************

function mixratio(e,p)

*------------------------------------------------------
* Given vapor pressure and pressure, return mixing ratio
*-------------------------------------------------------

mix=0.622*e/(p-e)

return(mix)

************************************************************************

function getrh(temp,dewp,pres)

tempk=temp+273.15
dewpk=dewp+273.15

es=satvap2(temp)

If (temp > 0) 
   A=2.53*pow(10,9)
   B=5420
Else
   A=3.41*pow(10,10)
   B=6130
Endif

w=A*0.622*exp(-B/dewpk)/pres
ws=mixratio(es,pres)

return(100*w/ws)

************************************************************************
function interp(array,pres)

*------------------------------------------------------------------------
* Interpolate inside array for pressure level pres.
* Returns estimated value of array at pressure pres.
*------------------------------------------------------------------------

"set gxout stat"
"set lev "pres
altpres=subwrd(result,4)
"q dims"
rec=sublin(result,4)
zlev=subwrd(rec,9)

If (zlev < 2 | zlev > _zmaxfile)
  Vest = -9999.0
Else
  If (altpres > pres)
    zlev=zlev+1
  Endif
  "set z "zlev
  PAbove=subwrd(result,4)
  "d "array"(lev="PAbove")"
  rec=sublin(result,8)
  VAbove=subwrd(rec,4)
  "set z "zlev-1
  PBelow=subwrd(result,4)
  "d "array"(lev="PBelow")"
  rec=sublin(result,8)
  VBelow=subwrd(rec,4)

* Now if we are in a region of missing data, find next good level.

  If (abs(VAbove) > 130 & zlev > 1 & zlev < _zmaxfile)
     zz=zlev
     While (abs(VAbove) > 130 & zz < _zmaxfile)
       zz=zz+1
       "set z "zz
       PAbove=subwrd(result,4)
       "d "array"(lev="PAbove")"
       rec=sublin(result,8)
       VAbove=subwrd(rec,4)
     EndWhile
  Endif

  If (abs(VBelow) > 130 & zlev > 1 & zlev < _zmaxfile)
      zz=zlev-1
      While (abs(VBelow) > 130 & zz > 1) 
        zz=zz-1
        "set z "zz
        PBelow=subwrd(result,4)
        "d "array"(lev="PBelow")"
        rec=sublin(result,8)
        VBelow=subwrd(rec,4)
      EndWhile
  Endif

  If (abs(VAbove) < 130 & abs(VBelow) < 130)
     Vest=VBelow+log(PBelow/pres)*(VAbove-VBelow)/log(PBelow/PAbove)
  Else
     Vest=-9999.0
  Endif

Endif

Return(Vest)

****************************************************************************

function GetUWnd(wspd,wdir)

*------------------------
* Get x-component of wind. 
*------------------------


If (wspd >= 0) 
   xwind=wspd*cos((270-wdir)*_dtr)
Else
   xwind = -9999.0
Endif
return(xwind)

**************************************************************************

function GetVWnd(wspd,wdir)

*-----------------------
* Get y-component of wind
*------------------------

If (wspd >= 0) 
   ywind=wspd*sin((270-wdir)*_dtr)
Else
   ywind = -9999.0
Endif
return(ywind)


*************************************************************************

function GetWSpd(xwind,ywind)


"set gxout stat"
"d mag("xwind","ywind")"
rec=sublin(result,8)
val=subwrd(rec,4)

return (val)

*************************************************************************

function GetWDir(xwind,ywind)

* Return wind direction given x and y components

"set gxout stat"
"define theta=270-"_rtd"*atan2("ywind","xwind")"
"d theta"
rec=sublin(result,8)
Dir=subwrd(rec,4)

If (Dir > 360)
   Dir=Dir-360
Endif

If (Dir < 0)
   Dir=360+Dir
Endif

return(Dir)

*************************************************************************

function GetXLoc(temp,pres)

*-------------------------------------------------
* Get x-location on skew-t based on temp, pressure
*-------------------------------------------------

xloc=(temp-_m1*log10(pres)-_m3)/_m2
return(xloc)

*************************************************************************
 
function GetTemp(xloc,pres) 

*------------------------------------------------- 
* Return temperature at location given by xloc,pres
*-------------------------------------------------

tempval=_m1*log10(pres)+_m2*xloc+_m3
return(tempval)

**************************************************************************

function GetTheta(temp,pres)         

*---------------------------------------------------
* Calculate theta for a given temperature and pressure
*---------------------------------------------------

theta=(temp+273.15)*pow(1000/pres,0.286)-273.15
return(theta)


*************************************************************************

function GetThet2(temp,dewp,pres)         

*---------------------------------------------------
* Calculate theta for a given temperature,dewp, and pressure
*---------------------------------------------------

tempk=273.15+temp
dewpk=273.15+dewp

es=satvap2(temp)
ws=mixratio(es,pres)

mix=10*getrh(temp,dewp,pres)*ws

exponent=0.2854*(1.0-0.00028*mix)
theta=(temp+273.15)*pow(1000/pres,exponent)-273.15
return(theta)

*************************************************************************

function Thetae(temp,dewp,pres)

*--------------------------------------------------------------
* Return equiv. pot. temp in Kelvin given temp, dewp in celsius
* From Bolton (1980) Mon Wea Rev
*--------------------------------------------------------------

es=satvap2(temp)
ws=mixratio(es,pres)
mix=10*getrh(temp,dewp,pres)*ws
theta=GetThet2(temp,dewp,pres)+273.15
TLcl=Templcl(temp,dewp)+273.15
thetae=theta*exp((3.376/TLcl-0.00254)*mix*(1.0+0.00081*mix))

return(thetae)

**************************************************************************


function int(i0)

*--------------------------
* Return integer of i0
*--------------------------
  i=0
  while(i<12)
    i=i+1
    if(substr(i0,i,1)='.')
      i0=substr(i0,1,i-1)
      break
    endif
  endwhile
return(i0)

*************************************************************************

function abs(i)

*----------------------------
* return absolute value of i
*----------------------------

  if (i < 0) 
     absval=-i
  else 
     absval=i
  endif

return(absval)

*************************************************************************

function log(i)

*---------------------------
* return natural log of i
*---------------------------

"set gxout stat"
"d log("i")"
rec=sublin(result,8)
val=subwrd(rec,4)
return(val)

*************************************************************************

function log10(i)

*--------------------------
* return log base 10 of i
*--------------------------

"set gxout stat"
"d log10("i")"
rec=sublin(result,8)
val=subwrd(rec,4)
return(val)

*************************************************************************

function pow(i,j)

*-------------------------------
* return power of i raised to j
*-------------------------------

"set gxout stat"
"d pow("i","j")"
rec=sublin(result,8)
val=subwrd(rec,4)
return(val)

************************************************************************

function cos(i)

*-----------------------------------------
* return cosine of i, where i is in radians
*------------------------------------------

"set gxout stat"
"d cos("i")"
rec=sublin(result,8)
val=subwrd(rec,4)
return(val)

************************************************************************

function sin(i)

*------------------------------------------
* return sine of i, where i is in radians
*------------------------------------------

"set gxout stat"
"d sin("i")"
rec=sublin(result,8)
val=subwrd(rec,4)
return(val)

************************************************************************

function exp(i)

*------------------------------------------
* return exponential of i
*------------------------------------------

"set gxout stat"
"d exp("i")"
rec=sublin(result,8)
val=subwrd(rec,4)
return(val)

***********************************************************************
function round(i)

rr=abs(1.0*i)
rr=int(rr+0.5)
if (i < 0)
   rr=-1*rr      
endif
return(rr)
