*  Draws all wx symbols
wx = 1
x = 1
y = 7.2
'set string 1 c 6'
'set strsiz 0.25'
while (wx<42) 
  'draw wxsym 'wx' 'x' 'y' 0.7 -1 6'
  'draw string 'x' '%(y-0.55)%' 'wx
  wx = wx + 1
  x = x + 1
  if (x > 10) 
    x = 1
    y = y - 1.5
  endif
endwhile 
