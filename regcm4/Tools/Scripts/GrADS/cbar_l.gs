*
*  Script to plot a legend for line graphs
*
*  Written MC, LSCE, Mars 2000
*  Modified CBAR_LINE script
*
function cbarline (args)
*
****** run script to create a legend when grads choosed line style color and symbols
******    cbar_l -x X -y Y -n number -t text -p
******
******  cbar_l parameters: 
******     x        x position on graph
******     y        y position on graph
******     number   # of line graphs, optional (max. 10)
******     text     text to be written for each line
******     must be double quoted (max 10 texts)
******     p        if set, user can click on graphic to place legend.
******                      
******  Use cbar_l without parameters to get a short discription.
******  Function call and switches can be uppercase and lowercase.
******  Defaults are 0 and no text.
******  If one keyword misses a value, cbar_l uses the default for the remaining values.
******
****** Example calls of cbar_l:
******
******    cbar_l -x 1 -y 7 -n 2 -t "1.Text" "2.Text"
******  same as
******    cbar_l -x 1 -y 7 -t "1.Text" "2.Text"
******  to place it manually 
******    cbar_l -t "1.Text" "2.Text" -p
******
****** Example calls of cbar_l from a script:
******
******   'cbar_l -x 1 -y 7 -c 2 -t "1.Text" "2.Text" '
******
****** Use defined variables
******   'cbar_l -x 3 -y 7.8 -t "Monat'month1'" "Monat'month2'" "Monat'month3'"'
******    
*
**********************************************************************************
**********************************************************************************
***********************   STANDARDS   **************************************************
i=1
f=subwrd(args,i)
xerr=0
yerr=0
nerr=0
terr=0
xx=0
yy=0
nn=0
k=0
l=0
m=0
t=0
placeit=0
***********************   USAGE   **************************************************
usage1="usage: cbar_l -x X -y Y -n # -t text -p"
usage2="-x x        position on graph (Default 1)"
usage3="-y y        position on graph (Default 1)"
usage4="-n number   # number of line (Default # of text)"
usage5="-t text     text for each line/symbol, must be double quoted"
usage6="-p          if set, user can click on graphic to set legend."
usage7=" "
usage8="If one keyword misses a value, uses default for the remaining"
usage9="To include zeros, zeros has to be double quoted"
usage10='Ex.: cbar_l -t "line1" "line2" "line3" -p'

***********************   HELP   **************************************************
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
 return
endif

***********************   OPTIIONS   **************************************************
while (f!='')
  opt = substr(f,2,1)
  if (opt='x' | opt='y' | opt='X' | opt='Y' | opt='n' | opt='N')
   i=i+1
   arg = subwrd(args,i)
   ttt=substr(arg,1,1)
   if (ttt='-' | ttt='')
    i=i-1
    if (opt='x' | opt='X')
     xerr=1
     say '-x error: X position missing!'
     say '          left side will be used.'
    endif
    if (opt='y' | opt='Y')
     yerr=1
     say '-y error: Y position missing!'
     say '          upper side will be used.'
    endif
    if (opt='n' | opt='N')
     nerr=1
     say '-n error: number not specified!'
     say '          take # of texts.'
    endif
   endif
   if ((opt='x' | opt='X') & xerr=0); xx=stripstr(arg); endif;
   if ((opt='y' | opt='Y') & yerr=0); yy=stripstr(arg); endif;
   if ((opt='n' | opt='N') & nerr=0); nn=stripstr(arg); endif;
  endif
  if (opt='p'); placeit=1; endif;
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
    say '          Draw lines and symbols but no text.'
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

***********************   CONSISTENCY   **************************************************
if (nn=0 & t=0)
   say 'nothing specified!'
   return
endif
maximum=max(nn,t)
if (nn!=0 & nn!=t)
   if (nn>t)
     say '# of lines greater # of texts'
     say 'Rest will be lines with empty strings'
   else
     say '# of lines lower # of texts'
     say 'discard # of lines; take # of texts'
   endif
endif
if (maximum>10)
  say '# limit to 10'
  say 'discard the rest'
  maximum=10
endif

***********************   FILL TEXT   **************************************************
p=0 
while (p<maximum)
  p=p+1
  if (p>t); texts.p=''; endif
endwhile

***********************   GRADS STANDARDS   **************************************************
cols.1=1
cols.2=3
cols.3=7
cols.4=2
cols.5=6
cols.6=9
cols.7=10
cols.8=11
cols.9=12
cols.10=15
lines.1=1;lines.2=1;lines.3=1;lines.4=1;lines.5=1;lines.6=1;lines.7=1;lines.8=1;lines.9=1;lines.10=1;
marks.1=2;marks.2=3;marks.3=4;marks.4=5;marks.5=1;marks.6=2;marks.7=3;marks.8=4;marks.9=5;marks.10=1;

***********************   GRAPHICS   **************************************************
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

***********************   DRAW LEGEND   **************************************************
o=0
while (o<maximum)
  o=o+1
  y=y-ywid
  'set line 'cols.o' 'lines.o
  'draw mark 'marks.o' 'xl' 'y' 0.1'
  'draw line 'xl' 'y' 'xr' 'y
  'draw mark 'marks.o' 'xr' 'y' 0.1'
  'set string 1 l 5'
  if (texts.o!='')
    'draw string '%(xr+0.1)%' 'y' 'texts.o
  endif
endwhile


exit
***********************   END   **************************************************
**********************************************************************************
**********************************************************************************
**********************************************************************************







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
