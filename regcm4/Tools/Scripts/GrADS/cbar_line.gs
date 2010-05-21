*
*  Script to plot a legend for line graphs
*
*  Written MC, LSCE, Mars 2000
*
function cbarline (args)
*
****** run script to create a legend
******    cbar_line -x X -y Y -c color -m mark -l linestyle -t text -p
******
******  cbar_line parameters:
******    x        x position on graph
******    y        y position on graph
******    color    colors of the line and symbols
******    mark     symbols type
******    text     text to be written for each line
******             must be double quoted
******    p        if set, user can click on graphic to place legend.
******
******  Use cbar_line without parameters to get a short discription.
******  Function call and switches can be uppercase and lowercase.
******
******  Defaults are no symbol, no lines, color white and no text.
******  If one keyword misses a value, 
******  cbar_line uses the default for the remaining values.
******  If you wish to include zeros (e.g. no symbol for one line),
******  zeros has to be double quoted (s. 2. example below).
******
****** Example calls of cbar_line:
******
******    cbar_line -x 1 -y 7 -c 2 3 -l 1 1 -m 1 2 -t "1.Text" "2.Text"
******
******    cbar_line -x 1 -y 7 -c 2 3 -l 1 1 -m 1 "0" -t "Symbol and line" "Only Line"
******  to place it manually 
******    cbar_line -x 1 -y 7 -c 2 3 -l 1 1 -m 1 "0" -t "Symbol and line" "Only Line" -p
******
****** Example calls of cbar_line from a script:
******
******    'cbar_line -x 1 -y 7 -c 2 3 -l 1 1 -m 1 "0" -t "Symbol and line" "Only Line" '
******
******   Use defined variables
******   'cbar_line -x 3 -y 7.8 -c 'color1' 'color2' 'color3' -m 'mark1' 'mark2' 'mark3' -l 'line1' 'line2' 'line3' -t "Monat'month1'" "Monat'month2'" "Monat'month3'"'
******    
*
***********************   STANDARDS   **************************************************
i=1
f=subwrd(args,i)
xerr=0
yerr=0
cerr=0
merr=0
lerr=0
terr=0
xx=0
yy=0
k=0
l=0
m=0
t=0
placeit=0
***********************   USAGE   **************************************************
usage1="usage: cbar_line -x X -y Y -c color -m mark -l linestyle -t text -p"
usage2="-x x        position on graph (Default 1)"
usage3="-y y        position on graph (Default 1)"
usage4="-c color    colors of the line and symbols (Default 1=white)"
usage5="-l line     line type"
usage6="-m mark     symbols type (Default 0=no symbol)"
usage7="-t text     text for each line/symbol, must be double quoted"
usage8="-p          if set, user can click on graphic to set legend."
usage9=" "
usage10="If one keyword misses a value, uses default for the remaining"
usage11="To include zeros, zeros has to be double quoted"
usage12='Ex.: cbar_line -c 1 2 3 -m 1 2 3 -l 1 2 3 -t "line1" "line2" "line3" -p'

if (f='')
 say usage1
 say usage2
 say usage3
 say usage4
 say usage5
 say usage6
 say usage7
 say usage8
 say usage9
 say usage10
 say usage11
 say usage12
 return
endif

while (f!='')
  opt = substr(f,2,1)
  if (opt='x' | opt='y' | opt='X' | opt='Y')
   i=i+1
   arg = subwrd(args,i)
   ttt=substr(arg,1,1)
   if (ttt='-' | ttt='')
    i=i-1
    if (opt='x' | opt='X')
     xerr=1
     say '-x error: X position missing!'
     say '          x=1 will be used.'
    endif
    if (opt='y' | opt='Y')
     yerr=1
     say '-y error: Y position missing!'
     say '          y=1 will be used.'
    endif
   endif
   if ((opt='x' | opt='X') & xerr=0); xx=stripstr(arg); endif;
   if ((opt='y' | opt='Y') & yerr=0); yy=stripstr(arg); endif;
  endif
  if (opt='p'); placeit=1; endif;
  if (opt='c' | opt='C')
   i=i+1
   arg = subwrd(args,i)
   ttt=substr(arg,1,1)
   if (ttt='-' | ttt='')
    cerr=1
    say '-c error: colors missing!'
    say '          take 1=white for all'
   endif
   if (cerr!=1)
    cerr=2
    while (ttt!='' & ttt!='-')
     k=k+1
     i=i+1
     cols.k=stripstr(arg)
     arg = subwrd(args,i)
     ttt=substr(arg,1,1)
    endwhile
    i=i-1
   endif
  endif
  if (opt='m' | opt='M')
   i=i+1
   arg = subwrd(args,i)
   ttt=substr(arg,1,1)
   if (ttt='-' | ttt='')
    merr=1
    say '-m error: marks missing!'
    say '          take 0= no marks'
   endif
   if (merr!=1)
    merr=2
    while (ttt!='' & ttt!='-')
     m=m+1
     i=i+1
     marks.m=stripstr(arg)
     arg = subwrd(args,i)
     ttt=substr(arg,1,1)
    endwhile
    i=i-1
   endif
  endif
  if (opt='l' | opt='L')
   i=i+1
   arg = subwrd(args,i)
   ttt=substr(arg,1,1)
   if (ttt='-' | ttt='')
    lerr=1
    say '-l error: lines missing!'
    say '          take 0= no lines'
   endif
   if (lerr!=1)
    lerr=2
    while (ttt!='' & ttt!='-')
     l=l+1
     i=i+1
     lines.l=stripstr(arg)
     arg = subwrd(args,i)
     ttt=substr(arg,1,1)
    endwhile
    i=i-1
   endif
  endif
  if (opt='t' | opt='T')
   i=i+1
   arg = subwrd(args,i)
   ttt=substr(arg,1,1)
   ta=substr(arg,1,1)
   if (ttt != '')
     te=substr(arg,wrdlen(arg),1)
   endif
   if (ttt='-' | ttt='')
    terr=1
    say '-t error: text missing!'
    say 'ABORT ! ! !'
    return
   endif
   if (terr!=1)
    terr=2
    if (ta!='"')
      say '-t error: text must be double quoted!'
      say 'ABORT ! ! !'
      return
    endif
    while (ttt!='' & ttt!='-')
     t=t+1
     i=i+1
     texts.t=stripstr(arg)
     while (te!='"')
      arg = subwrd(args,i)
      ttt=substr(arg,1,1)
      ta=substr(arg,1,1)
      if (ttt != '')
       te=substr(arg,wrdlen(arg),1)
      endif
      if (ttt='-' | ttt='')
        say '-t error: text must be quoted!'
        say 'ABORT ! ! !'
        return
      endif
      if (te='"') 
        arg=stripstr(arg)
      endif
      texts.t=texts.t%' '%arg
      i=i+1
     endwhile
     arg = subwrd(args,i)
     ttt=substr(arg,1,1)
     ta=substr(arg,1,1)
     if (ttt != '')
      te=substr(arg,wrdlen(arg),1)
     endif
    endwhile
    i=i-1
   endif
  endif
  i=i+1
  f=subwrd(args,i)
endwhile

if (k!=m | k!=l | k!=t | m!=l | m!=t | l!=t)
   say '# of colors, marks, lines and/or text not equal.'
   say '# colors 'k
   say '# marks  'm
   say '# lines  'l
   say '# texts  't
endif

tmp=max(k,l)
tmp1=max(tmp,m)
maximum=max(tmp1,t)

p=0 
while (p<maximum)
  p=p+1
  if (p>k); cols.p=1; endif
  if (p>l); lines.p=0; endif
  if (p>m); marks.p=0; endif
  if (p>t); texts.p=''; endif
endwhile


'query gxinfo'
****** example: rec2 =>   Page Size = 11 by 8.5
rec2 = sublin(result,2)
******          rec3 =>   X Limits = 1.3 to 10.2
rec3 = sublin(result,3)
******          rec4 =>   Y Limits = 1.58 to 6.92
rec4 = sublin(result,4)

xsiz = subwrd(rec2,4)
ysiz = subwrd(rec2,6)
yhi = subwrd(rec4,6)
ylo = subwrd(rec4,4)
xhi = subwrd(rec3,6)
xlo = subwrd(rec3,4)
xd = xsiz - xhi
yd = ysiz - ylo

if (xx=0); xx=xlo; endif;
if (yy=0); yy=yhi; endif;
if (placeit=1)
  say ''
  say 'Click where you want the left upper corner of the legend'
  'query bpos'
  x = subwrd(result,3)
  y = subwrd(result,4)
  say 'Print legend at X Y: 'x' 'y
  xx=x
  yy=y
endif

xl=xx
xwid = xsiz/20
xr=xl+xwid
y=yy
ywid=yd/20
y=y+ywid

'set strsiz 0.12 0.13'

o=0
while (o<maximum)
  o=o+1
  y=y-ywid
  'set line 'cols.o' 'lines.o
  'draw mark 'marks.o' 'xl' 'y' 0.1'
  'draw line 'xl' 'y' 'xr' 'y
  'draw mark 'marks.o' 'xr' 'y' 0.1'
  'set string 1 l 5'
  'draw string '%(xr+0.1)%' 'y' 'texts.o
endwhile


exit







*************************   MAXIMUM NUMBER   *************************************************
function max(wert1, wert2)
* liefert die groessere von zwei Zahlen
if (wert1 <= wert2)
 return (wert2)
else
 return (wert1)
endif
return

*************************   LENGTH OF WORD   *************************************************
function wrdlen (arg)
i=1
s=subwrd(arg,1)
if (s='') 
 return 0
else
 t=substr(s,i,i)
 while (t!='')
  i=i+1
  t=substr(s,i,i)
 endwhile 
endif

return (i-1)


*************************   STRIP STRING OF DOUBLE QUOTES  *************************************************
function stripstr (arg)

i=1
len=0
s=subwrd(arg,i)
if (s='') 
 return 0
else
 a=substr(s,1,1)
 if (a='"')
  s=substr(s,2,wrdlen(s))
 endif
 e=substr(s,wrdlen(s),1)
 if (e='"')
  s=substr(s,1,wrdlen(s)-1)
 endif
endif

return s
