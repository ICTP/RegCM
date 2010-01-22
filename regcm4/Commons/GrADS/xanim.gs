* xanim.gs
* A GrADS script for Xwindows animation
*
* Features:
*	 Uses double buffering for smooth animations 
*	 Contour/no-contour option
*	 User-specified delay or mouse-controlled animation
*        The default contour levels and colors from the first display 
*        of the animation will be used for the remaining displays.
*
* Setup: To override the default contour levels and/or colors, 
*        execute 'set clevs' and 'set ccols' and 'set dbuff on' 
*        before running xanim.gs with the -levs option. 
*
*        Set the variable "colorbar" to point to your version of cbar.
*        i.e., change the line: colorbar = 'cbarn.gs'
*
* Usage: xanim [options] variable title
*
* Available options are:
*   -contour      Displays a shaded plot with contours overlaid (default)
*   -nocontour    Displays only the shaded plot
*   -notitle      No title will be drawn
*   -levs         Use pre-specified levs and ccols, not the defaults.
*                 You must 'set dbuff on' before calling xanim.gs
*   -pause        Allows mouse button control of animation:
*                 left=forward, right=backwards, middle=quit
*   -sec n        Pauses n seconds after each frame
*   -repeat n     Loops through images n+1 times 
*   -skip n       Displays every n-th frame (1, n+1, 2n+1, ...)
*   -script name  Executes named script for each display
*
* Notes:  'variable' can be any GrADS expression as long as
*         it does not contain blanks. e.g., (tmax-273)*1.8+32
*
*         If 'title' is omitted, the 'variable' expression is used.
*
*         If the '-script name' option is used, then the 'variable' 
*         and 'title' arguments are ignored. 
*         
* Written by Wesley Ebisuzaki June 2001
* Comments and cosmetic adjustments added by J.M.Adams 
* v1.0b
*
function main(args)

* Default settings
contour  = 1
notitle  = 0
pause    = 0
sec      = 0
repeat   = 0
skip     = 1
dbuff    = 0
script   = ''
colorbar = 'cbarn.gs'
colorbar = 'cbarb.gs'

* Parse all the options
i = 1
wrd = subwrd(args,i)
while (substr(wrd,1,1) = '-')
   if (wrd = '-contour')
      contour = 1
   endif
   if (wrd = '-nocontour')
      contour = 0
   endif
   if (wrd = '-notitle')
      notitle = 1
   endif
   if (wrd = '-levs')
      dbuff = 1
   endif
   if (wrd = '-pause')
      pause = 1
      sec = 0
   endif
   if (wrd = '-sec')
      i = i + 1
      sec = subwrd(args,i)
   endif
   if (wrd = '-repeat')
      i = i + 1
      repeat = subwrd(args,i)
   endif
   if (wrd = '-skip')
      i = i + 1
      skip = subwrd(args,i)
   endif
   if (wrd = '-script')
      i = i + 1
      script = subwrd(args,i)
      notitle = 1
      contour = 0
   endif
   i = i + 1
   wrd = subwrd(args,i)
endwhile

* If scripts are not used, get the variable expression
if (script = '')
  var = subwrd(args,i)
  i = i + 1
endif

* Get the title which may contain more than one word
* If no title is specified, it will be the variable expression
title = ''
while (subwrd(args,i) != '')
   title = title%' '%subwrd(args,i)
   i = i+1
endwhile
if (title = '' & notitle = 0)
   title = var
endif

* Get the dimension environment
'q dim'
diminfo = result
line5 = sublin(diminfo,5)
time1 = subwrd(line5,11)
time2 = subwrd(line5,13)
if (time2 = '')
  say 'The time dimension must be varying'
  exit 8
endif

* Start plotting
if (dbuff = 0)
  'set dbuff on'
endif
qshade = 0
while (repeat >= 0)
   it = time1
   while (it <= time2)
      'set t 'it

*     Add display commands here
      'set gxout shaded'
      'set grads off'

      if (qshade = 1)
         'set clevs 'lev
         'set ccols 'color
      endif
      if (script = '')
         'd 'var
         r1 = subwrd(result,1)
         if (r1 != 'Cannot')
            'run 'colorbar
         endif
      else
         'run 'script
      endif

*     Get the shaded contour levels and colors if not pre-specified
      if (qshade = 0)
         'query shades'
         shdinfo = result
         r1 = subwrd(shdinfo,1)
         if (r1 != 'None')
            nlevs = subwrd(shdinfo,5)
            rec = sublin(shdinfo,2)
            color = subwrd(rec,1)
            lev = ''
            n = 2
            while (n <= nlevs)
               rec = sublin(shdinfo,n+1)
               color = color%' '%subwrd(rec,1)
               lev = lev%' '%subwrd(rec,2)
               n = n + 1
            endwhile 
            qshade = 1
         endif
      endif

*     Settings for contour overlay
      if (contour = 1 & qshade = 1)
        'set gxout contour'
        'set clevs 'lev
        'set ccolor 15'
        'set ccols 15'
        'd 'var
      endif

*     Time stamp is added to the title 
      'q dim'
      rec = sublin(result,5)
      time = subwrd(rec,6)
      if (notitle = 0)
        'draw title 'title'   'time
      endif

*     Swap buffers and then pause or sleep
      'swap'
      if (pause != 0)
         'q pos'
         mousekey = subwrd(result,5)
         if (mousekey = 2) 
            it = time2
            repeat = 0
         endif
         if (mousekey = 3)
            if (it < time1)
               it = time2+skip
            endif
            it = it-skip-skip
         endif
      else
         if (sec != 0)
           '!sleep ' sec
         endif
      endif

*     Move to next time step
      it = it+skip
   endwhile

   repeat = repeat - 1
   if (repeat >= 0 & sec != 0)
      '!sleep ' sec
      '!sleep ' sec
   endif

endwhile

* Clean up
'set t 'time1' 'time2
if (dbuff = 0) 
  'set dbuff off'
endif


