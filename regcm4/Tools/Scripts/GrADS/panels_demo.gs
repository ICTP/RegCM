* panels_demo.gs 
* 
* This script demonstrates the use of panels.gsf, a dynamically 
* loaded script function that sets up the virtual page commands 
* for a multi-panel plot. 
* 
* Written by JMA March 2001
*
function main(args)
  rc = gsfallow("on")
  if (args='') 
    say 'Two arguments are required: the # of rows and # of columns'
    return 
  else 
    nrows = subwrd(args,1)
    ncols = subwrd(args,2)
  endif

  'use model.ctl'
  'reset'
  panels(args)
  p = 1
  ptot = nrows * ncols
  'set mproj scaled'

* Loop through each panel and draw a plot
  while (p <= ptot)
    _vpg.p
    'set t 'p
    'set grads off'
    'd t'
    p = p + 1
  endwhile
