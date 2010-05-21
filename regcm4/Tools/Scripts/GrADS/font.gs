*
*	usage (in GrADS)
*
* run font # 
*	where # is the font number
*
*
function font(arg)
if (arg='')
  arg = 1
endif

fch = '!"#$%&''()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOP'
fch = fch % 'QRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

'set font 1'
'set strsiz 0.3'
'set string 1 c 3'
'draw string 5.5 8.0 Font Set 'arg
x = 1.0
y = 7.0
i = 1
'set font 0'
'set strsiz 0.2'
'set string 1 c 1'
while (i<95)
  'draw string 'x' 'y' 'substr(fch,i,1)' `'%arg%substr(fch,i,1)
  x = x + 1.5
  if (x>10.5)
    x = 1.0
    y = y - 0.4
  endif
  i = i + 1
endwhile
