function main(args)

  input  = subwrd(args, 1)
  varp  = subwrd(args, 2)

  outname = input'_'varp
  caption = 'Plot of last status of var 'varp' from 'input
  titl    = 'Himet - Cetemps'

  'open 'input
  'set t last'
  'q time'
  dow   = subwrd(result, 6)
  res   = subwrd(result, 3)
  hour  = substr(res, 1, 2)
  day   = substr(res, 4, 2)
  month = substr(res, 6, 3)
  year  = substr(res, 9, 5)
  leadstr = dow', 'day' 'month' 'year'  'hour' UTC'

  'set csmooth on'
  'set mpdset hires'
  'set datawarn off'
  'set grads off'
  'set gxout shaded'

  'd 'varp
  'draw map'

  'set font 0'
  'set string 1 tl 8'
  'set strsiz 0.18'
  'draw string 0.44 8.49 'titl
  'set strsiz 0.10'
  'draw string 0.5 0.50 'caption
  'set strsiz 0.18'
  'set string 4 tl 8'
  'draw string 4.80 8.49 'leadstr
  'printim 'outname'.png x751 y580 white'
  'quit'

return
