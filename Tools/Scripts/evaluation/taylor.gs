***************************************************************************************
*       $Id: taylor.gs,v 1.34 2012/07/17 23:55:28 bguan Exp $
*       Copyright (C) 2012 Bin Guan.
*       Distributed under GNU/GPL.
***************************************************************************************
function taylor(arg)
*
* Draw a Taylor diagram.
*
rc=gsfallow('on')

tmpdir='/tmp'
whoamifile='./.whoami.bGASL'
'!whoami>'whoamifile
whoami=sublin(read(whoamifile),2)
rc=close(whoamifile)
'!unlink 'whoamifile
mytmpdir=tmpdir'/bGASL-'whoami
'!mkdir -p 'mytmpdir

pi=3.14159
small_spacing=0.011
tick_length=0.05
'query pp2xy 0 0'
tmpxa=subwrd(result,3)
'query pp2xy 1 1'
tmpxb=subwrd(result,3)
rvratio=tmpxb-tmpxa
small_spacing=small_spacing*rvratio
tick_length=tick_length*rvratio

*
* Parse -s and -r options (standard deviation and correlation).
*
num_var=parseopt(arg,'-','s','std')
num_var2=parseopt(arg,'-','r','corr')

*
* Parse -i option (input file; will overwrite -s and -r).
*
_.file.1=''
_.columnTEXT.1=''
rc=parseopt(arg,'-','i','file')
if(_.file.1!='')
  flag=1
  cnt=0
  while(flag)
    result=read(_.file.1)
    status=sublin(result,1)
    if(status!=0)
      flag=0
    else
      cnt=cnt+1
      line.cnt=sublin(result,2)
    endif
  endwhile
  num_line=cnt
  _.rowstart.1=1
  _.rowend.1=num_line
  _.columnstd.1=1
  _.columncorr.1=2
  rc=parseopt(arg,'-','rows','rowstart')
  rc=parseopt(arg,'-','rowe','rowend')
  rc=parseopt(arg,'-','cols','columnstd')
  rc=parseopt(arg,'-','colr','columncorr')
  rc=parseopt(arg,'-','colT','columnTEXT')
  num_var=_.rowend.1-_.rowstart.1+1
  cnt=1
  while(cnt<=num_var)
    tmpcnt=cnt+_.rowstart.1-1
    _.std.cnt=subwrd(line.tmpcnt,_.columnstd.1)
    _.corr.cnt=subwrd(line.tmpcnt,_.columncorr.1)
    if(_.columnTEXT.1!='')
      _.TEXT.cnt=subwrd(line.tmpcnt,_.columnTEXT.1)
    endif
    cnt=cnt+1
  endwhile
endif

*
* -s/-r/-i options were parsed first to know # of input variables.
*
*if(num_var=0 | (_.file.1='' & num_var2=0))
*  usage()
*  return
*endif

*
* Initialize other options.
*
cnt=1
while(cnt<=num_var)
  _.mark.cnt=cnt+1
  _.marksize.cnt=0.11
  _.color.cnt=cnt
  if(_.columnTEXT.1='')
    _.TEXT.cnt='Variable 'cnt
  endif
  cnt=cnt+1
endwhile
_.xlim.1=2.5
_.append.1=0
_.text.1='A'
_.text.2='B'
_.text.3='C'
_.text.4='D'
_.text.5='E'
_.text.6='F'
_.text.7='G'
_.text.8='H'
_.text.9='I'
_.text.10='J'
_.text.11='K'
_.text.12='L'
_.text.13='M'
_.text.14='N'
_.text.15='O'
_.text.16='P'
_.text.17='Q'
_.text.18='R'
_.text.19='S'
_.text.20='T'
_.text.21='U'
_.text.22='V'
_.text.23='W'
_.text.24='X'
_.text.25='Y'
_.text.26='Z'
num_levs=1
_.levs.1=1
num_levr=0
_.levc.1=1

*
* Parse other options.
*
rc=parseopt(arg,'-','l','xlim')
if(_.xlim.1<0)
  halfcircle=1
  _.xlim.1=math_abs(_.xlim.1)
else
  halfcircle=0
endif
rc=parseopt(arg,'-','m','mark')
if(rc=1)
  cnt=2
  while(cnt<=num_var)
    _.mark.cnt=_.mark.1
    cnt=cnt+1
  endwhile
endif
rc=parseopt(arg,'-','c','color')
if(rc=1)
  cnt=2
  while(cnt<=num_var)
    _.color.cnt=_.color.1
    cnt=cnt+1
  endwhile
endif
rc=parseopt(arg,'-','z','marksize')
if(rc=1)
  cnt=2
  while(cnt<=num_var)
    _.marksize.cnt=_.marksize.1
    cnt=cnt+1
  endwhile
endif
rc=parseopt(arg,'-','t','text')
if(rc=1)
  cnt=2
  while(cnt<=num_var)
    _.text.cnt=_.text.1
    cnt=cnt+1
  endwhile
endif
rc=parseopt(arg,'-','T','TEXT')
if(rc=1)
  cnt=2
  while(cnt<=num_var)
    _.TEXT.cnt=_.TEXT.1
    cnt=cnt+1
  endwhile
endif
num_levs=parseopt(arg,'-','levs','levs')
if(num_levs=0)
  num_levs=1
endif
num_levr=parseopt(arg,'-','levr','levr')
rc=parseopt(arg,'-','levc','levc')
rc=parseopt(arg,'-','append','append')

*
* Make a .ctl file with no data.
*
ctllines=10
ctlline.1='DSET ^%y4.dat'
ctlline.2='UNDEF -9999'
ctlline.3='options template'
ctlline.4='xdef 2 levels 0 216'
ctlline.5='ydef 2 levels -90 90'
ctlline.6='zdef 1 levels 1000'
ctlline.7='tdef 1 linear 01jan0001 1dy'
ctlline.8='VARS 1'
ctlline.9='var 0 99 var'
ctlline.10='ENDVARS'
cnt=1
while(cnt<=ctllines)
  status=write(mytmpdir'/taylor.ctl~',ctlline.cnt)
  cnt=cnt+1
endwhile
status=close(mytmpdir'/taylor.ctl~')
'open 'mytmpdir'/taylor.ctl~'
* At this point, at least one file is opened.

qdims()

*
* My trick to get ACTUAL plot area (gxinfo does not give true plot area before a plot is actually made).
*
'set lon 0 216'
'set lat -90 90'
'set z 1'
'set t 1'
'set mproj scaled'
* Above: 'set mproj latlon' may not work properly since the requested plot area (e.g., by 'set parea' elsewhere) may not be optimal.
'set grid off'
'set grads off'
'set frame off'
'set xlab off'
'set ylab off'
'set mpdraw off'
'set gxout contour'
'set clevs -1e9'
'display lat'
'q gxinfo'
line3=sublin(result,3)
line4=sublin(result,4)
xa=subwrd(line3,4)
xb=subwrd(line3,6)
ya=subwrd(line4,4)
yb=subwrd(line4,6)
xlength=xb-xa
ylength=yb-ya

*
* Draw frame.
*
if(!halfcircle)
  if(xlength/ylength<1)
    radius=xlength
  else
    radius=ylength
  endif
  xa=xa+(xlength-radius)/2
  xb=xb-(xlength-radius)/2
  ya=ya+(ylength-radius)/2
  yb=yb-(ylength-radius)/2
  xorig=xa
  'draw line 'xa' 'ya' 'xb' 'ya
  'draw line 'xa' 'ya' 'xa' 'yb
  deg=0
  while(deg<=89.9)
    x1=xorig+radius*math_cos(deg/180*pi)
    y1=ya+radius*math_sin(deg/180*pi)
    x2=xorig+radius*math_cos((deg+0.1)/180*pi)
    y2=ya+radius*math_sin((deg+0.1)/180*pi)
    'draw line 'x1' 'y1' 'x2' 'y2
    deg=deg+0.1
  endwhile
endif

*
* Draw frame.
*
if(halfcircle)
  if(xlength/ylength<2)
    radius=xlength/2
  else
    radius=ylength
  endif
  xa=xa+(xlength-2*radius)/2
  xb=xb-(xlength-2*radius)/2
  ya=ya+(ylength-radius)/2
  yb=yb-(ylength-radius)/2
  xorig=xa+radius
  'draw line 'xa' 'ya' 'xb' 'ya
  deg=0
  while(deg<=179.9)
    x1=xorig+radius*math_cos(deg/180*pi)
    y1=ya+radius*math_sin(deg/180*pi)
    x2=xorig+radius*math_cos((deg+0.1)/180*pi)
    y2=ya+radius*math_sin((deg+0.1)/180*pi)
    'draw line 'x1' 'y1' 'x2' 'y2
    deg=deg+0.1
  endwhile
endif

unit=radius/_.xlim.1

*
* Draw CORR tick marks.
*
cnt=0
while(cnt<=0.9)
  x1=xorig+_.xlim.1*unit*math_cos(math_acos(cnt))
  y1=ya+_.xlim.1*unit*math_sin(math_acos(cnt))
  x2=xorig+(_.xlim.1*unit-tick_length)*math_cos(math_acos(cnt))
  y2=ya+(_.xlim.1*unit-tick_length)*math_sin(math_acos(cnt))
  x3=xorig+(_.xlim.1*unit+tick_length)*math_cos(math_acos(cnt))
  y3=ya+(_.xlim.1*unit+tick_length)*math_sin(math_acos(cnt))
  'draw line 'x1' 'y1' 'x2' 'y2
  'set string 1 l 4 'math_acos(cnt)/pi*180
  'draw string 'x3' 'y3' 'cnt
  cnt=cnt+0.1
endwhile
cnt=0.91
while(cnt<=1.0)
  x1=xorig+_.xlim.1*unit*math_cos(math_acos(cnt))
  y1=ya+_.xlim.1*unit*math_sin(math_acos(cnt))
  x2=xorig+(_.xlim.1*unit-tick_length)*math_cos(math_acos(cnt))
  y2=ya+(_.xlim.1*unit-tick_length)*math_sin(math_acos(cnt))
  x3=xorig+(_.xlim.1*unit+tick_length)*math_cos(math_acos(cnt))
  y3=ya+(_.xlim.1*unit+tick_length)*math_sin(math_acos(cnt))
  'draw line 'x1' 'y1' 'x2' 'y2
  if(cnt=0.95 | cnt=0.99 | cnt=1.00)
    'set string 1 l 4 'math_acos(cnt)/pi*180
    'draw string 'x3' 'y3' 'cnt
  endif
  cnt=cnt+0.01
endwhile
if(halfcircle)
  cnt=0
  while(cnt>=-0.9)
    x1=xorig+_.xlim.1*unit*math_cos(math_acos(cnt))
    y1=ya+_.xlim.1*unit*math_sin(math_acos(cnt))
    x2=xorig+(_.xlim.1*unit-tick_length)*math_cos(math_acos(cnt))
    y2=ya+(_.xlim.1*unit-tick_length)*math_sin(math_acos(cnt))
    x3=xorig+(_.xlim.1*unit+tick_length)*math_cos(math_acos(cnt))
    y3=ya+(_.xlim.1*unit+tick_length)*math_sin(math_acos(cnt))
    'draw line 'x1' 'y1' 'x2' 'y2
    'set string 1 r 4 'math_acos(cnt)/pi*180-180
    'draw string 'x3' 'y3' 'cnt
    cnt=cnt-0.1
  endwhile
  cnt=-0.91
  while(cnt>=-1.0)
    x1=xorig+_.xlim.1*unit*math_cos(math_acos(cnt))
    y1=ya+_.xlim.1*unit*math_sin(math_acos(cnt))
    x2=xorig+(_.xlim.1*unit-tick_length)*math_cos(math_acos(cnt))
    y2=ya+(_.xlim.1*unit-tick_length)*math_sin(math_acos(cnt))
    x3=xorig+(_.xlim.1*unit+tick_length)*math_cos(math_acos(cnt))
    y3=ya+(_.xlim.1*unit+tick_length)*math_sin(math_acos(cnt))
    'draw line 'x1' 'y1' 'x2' 'y2
    if(cnt=-0.95 | cnt=-0.99 | cnt=-1.00)
      'set string 1 r 4 'math_acos(cnt)/pi*180-180
      'draw string 'x3' 'y3' 'cnt
    endif
    cnt=cnt-0.01
  endwhile
endif

*
* Draw x tick marks.
*
tick_step=math_nint(_.xlim.1*10/5)/10
if(tick_step=0.3 | tick_step=0.4)
  tick_step=0.25
endif
if(tick_step>=0.6 & tick_step<=0.9)
  tick_step=0.5
endif
if(tick_step>1 & tick_step<2)
  tick_step=1
endif
if(tick_step>2 & tick_step<5)
  tick_step=2
endif
cnt=0
while(cnt<=_.xlim.1)
  x1=xorig+cnt*unit
  y1=ya
  x2=xorig+cnt*unit
  y2=ya+tick_length
  x3=xorig+cnt*unit
  y3=ya-tick_length
  'draw line 'x1' 'y1' 'x2' 'y2
  'set string 1 tc'
  'draw string 'x3' 'y3' 'cnt
  cnt=cnt+tick_step
endwhile
if(halfcircle)
  cnt=0
  while(cnt<=_.xlim.1)
    x1=xorig-cnt*unit
    y1=ya
    x2=xorig-cnt*unit
    y2=ya+tick_length
    x3=xorig-cnt*unit
    y3=ya-tick_length
    'draw line 'x1' 'y1' 'x2' 'y2
    'set string 1 tc'
    'draw string 'x3' 'y3' 'cnt
    cnt=cnt+tick_step
  endwhile
endif

*
* Draw y tick marks.
*
if(!halfcircle)
  cnt=0
  while(cnt<=_.xlim.1)
    x1=xorig
    y1=ya+cnt*unit
    x2=xorig+tick_length
    y2=ya+cnt*unit
    x3=xorig-tick_length
    y3=ya+cnt*unit
    'draw line 'x1' 'y1' 'x2' 'y2
    'set string 1 r'
    'draw string 'x3' 'y3' 'cnt
    cnt=cnt+tick_step
  endwhile
endif

*
* Draw STD levels.
*
if(num_levs!=2 | num_levr!=2)
  'set line '_.levc.1' 1'
  cnt=1
  while(cnt<=num_levs)
    dash_step=(5/_.levs.cnt)*(_.xlim.1/2.5)*rvratio
    if(halfcircle);dash_step=dash_step*2;endif
    if(!halfcircle)
      deg=0
      while(deg<=90-dash_step/2)
        x1=xorig+unit*_.levs.cnt*math_cos(deg/180*pi)
        y1=ya+unit*_.levs.cnt*math_sin(deg/180*pi)
        x2=xorig+unit*_.levs.cnt*math_cos((deg+dash_step/2)/180*pi)
        y2=ya+unit*_.levs.cnt*math_sin((deg+dash_step/2)/180*pi)
        'draw line 'x1' 'y1' 'x2' 'y2
        deg=deg+dash_step
      endwhile
*     begin draw last partial dash if needed
      if(deg<=90)
        x1=xorig+unit*_.levs.cnt*math_cos(deg/180*pi)
        y1=ya+unit*_.levs.cnt*math_sin(deg/180*pi)
        x2=xorig+unit*_.levs.cnt*math_cos((90)/180*pi)
        y2=ya+unit*_.levs.cnt*math_sin((90)/180*pi)
        'draw line 'x1' 'y1' 'x2' 'y2
      endif
*     end draw last partial dash if needed
    else
      deg=0
      while(deg<=180-dash_step/2)
        x1=xorig+unit*_.levs.cnt*math_cos(deg/180*pi)
        y1=ya+unit*_.levs.cnt*math_sin(deg/180*pi)
        x2=xorig+unit*_.levs.cnt*math_cos((deg+dash_step/2)/180*pi)
        y2=ya+unit*_.levs.cnt*math_sin((deg+dash_step/2)/180*pi)
        'draw line 'x1' 'y1' 'x2' 'y2
        deg=deg+dash_step
      endwhile
*     begin draw last partial dash if needed
      if(deg<=180)
        x1=xorig+unit*_.levs.cnt*math_cos(deg/180*pi)
        y1=ya+unit*_.levs.cnt*math_sin(deg/180*pi)
        x2=xorig+unit*_.levs.cnt*math_cos((180)/180*pi)
        y2=ya+unit*_.levs.cnt*math_sin((180)/180*pi)
        'draw line 'x1' 'y1' 'x2' 'y2
      endif
*     end draw last partial dash if needed
    endif
    cnt=cnt+1
  endwhile
endif

*
* Draw CORR levels.
*
if(num_levs!=2 | num_levr!=2)
  'set line '_.levc.1' 2'
  cnt=1
  while(cnt<=num_levr)
    x1=xorig
    y1=ya
    x2=xorig+(_.xlim.1*unit)*math_cos(math_acos(_.levr.cnt))
    y2=ya+(_.xlim.1*unit)*math_sin(math_acos(_.levr.cnt))
    'draw line 'x1' 'y1' 'x2' 'y2
    cnt=cnt+1
  endwhile
  'set line 1 1'
endif

*
* Draw STD and CORR "box": STD boundaries.
*
if(num_levs=2 & num_levr=2)
  'set line '_.levc.1' 1'
  cnt=1
  while(cnt<=num_levs)
    dash_step=(5/_.levs.cnt)*(_.xlim.1/2.5)*rvratio
    if(halfcircle);dash_step=dash_step*2;endif
    deg_s=math_acos(_.levr.1)/pi*180
    deg_e=math_acos(_.levr.2)/pi*180
    if(deg_s>deg_e)
      deg_tmp=deg_s
      deg_s=deg_e
      deg_e=deg_tmp
    endif
    deg=deg_s
    while(deg<=deg_e-dash_step)
      x1=xorig+unit*_.levs.cnt*math_cos(deg/180*pi)
      y1=ya+unit*_.levs.cnt*math_sin(deg/180*pi)
      x2=xorig+unit*_.levs.cnt*math_cos((deg+dash_step)/180*pi)
      y2=ya+unit*_.levs.cnt*math_sin((deg+dash_step)/180*pi)
      'draw line 'x1' 'y1' 'x2' 'y2
      deg=deg+dash_step
    endwhile
*   begin draw last partial dash if needed
    if(deg<=deg_e)
      x1=xorig+unit*_.levs.cnt*math_cos(deg/180*pi)
      y1=ya+unit*_.levs.cnt*math_sin(deg/180*pi)
      x2=xorig+unit*_.levs.cnt*math_cos((deg_e)/180*pi)
      y2=ya+unit*_.levs.cnt*math_sin((deg_e)/180*pi)
      'draw line 'x1' 'y1' 'x2' 'y2
    endif
*   end draw last partial dash if needed
    cnt=cnt+1
  endwhile
endif

*
* Draw STD and CORR "box": CORR boundaries.
*
if(num_levs=2 & num_levr=2)
  'set line '_.levc.1' 1'
  cnt=1
  while(cnt<=num_levr)
    x1=xorig+(_.levs.1*unit)*math_cos(math_acos(_.levr.cnt))
    y1=ya+(_.levs.1*unit)*math_sin(math_acos(_.levr.cnt))
    x2=xorig+(_.levs.2*unit)*math_cos(math_acos(_.levr.cnt))
    y2=ya+(_.levs.2*unit)*math_sin(math_acos(_.levr.cnt))
    'draw line 'x1' 'y1' 'x2' 'y2
    cnt=cnt+1
  endwhile
  'set line 1 1'
endif

*
* Plot data.
*
cnt=1
while(cnt<=num_var)
  x=xorig+_.std.cnt*unit*math_cos(math_acos(_.corr.cnt))
  y=ya+_.std.cnt*unit*math_sin(math_acos(_.corr.cnt))
  'set line '_.color.cnt
  'draw mark '_.mark.cnt' 'x' 'y' '_.marksize.cnt*rvratio
  'set line 1'
  cnt=cnt+1
endwhile
cnt=1
while(cnt<=num_var)
  x=xorig+_.std.cnt*unit*math_cos(math_acos(_.corr.cnt))
  y=ya+_.std.cnt*unit*math_sin(math_acos(_.corr.cnt))
  if(_.text.cnt!='')
    'set string 1 bc'
    'draw string 'x' 'y+_.marksize.cnt*rvratio/2+small_spacing' '_.text.cnt
  endif
  cnt=cnt+1
endwhile

*
* Write legend file.
*
if(_.append.1!=1)
  rc=write(mytmpdir'/legend.txt~','taylor')
endif
cnt=1
while(cnt<=num_var)
  line=_.mark.cnt' '_.marksize.cnt' '_.color.cnt' | '_.text.cnt' | '_.TEXT.cnt
  if(_.append.1!=1)
    rc=write(mytmpdir'/legend.txt~',line)
  else
    rc=write(mytmpdir'/legend.txt~',line,append)
  endif
  cnt=cnt+1
endwhile
rc=close(mytmpdir'/legend.txt~')

*
* Restore original dimension environment.
*
'set mproj latlon'
'set grid on'
'set grads on'
'set frame on'
'set xlab on'
'set ylab on'
'set mpdraw on'
'set gxout contour'
'set x '_xs' '_xe
'set y '_ys' '_ye
'set z '_zs' '_ze
'set t '_ts' '_te
'close 'file_number()

return
***************************************************************************************
function file_number()
*
* Get the number of files opened.
*
'q files'
line1=sublin(result,1)
if(line1='No files open')
  return 0
endif

lines=1
while(sublin(result,lines+1)!='')
  lines=lines+1
endwhile

return lines/3
***************************************************************************************
function parseopt(instr,optprefix,optname,outname)
*
* Parse an option, store argument(s) in a global variable array.
*
rc=gsfallow('on')
cnt=1
cnt2=0
while(subwrd(instr,cnt)!='')
  if(subwrd(instr,cnt)=optprefix''optname)
    cnt=cnt+1
    word=subwrd(instr,cnt)
    while(word!='' & (valnum(word)!=0 | substr(word,1,1)''999!=optprefix''999))
      cnt2=cnt2+1
      _.outname.cnt2=parsestr(instr,cnt)
      if(_end_wrd_idx=-9999);return cnt2;endif
      cnt=_end_wrd_idx+1
      word=subwrd(instr,cnt)
    endwhile
  endif
  cnt=cnt+1
endwhile
return cnt2
***************************************************************************************
function usage()
*
* Print usage information.
*
say '  Draw a Taylor diagram.'
say ''
say '  USAGE: taylor -s <STD1> [<STD2>...] -r <CORR1> [<CORR2>...] [-i <file>] [-rows <row_start>] [-rowe <row_end>]'
say '         [-cols <column_STD>] [-colr <column_CORR>] [-colT <column_TEXT>] [-l <limit>] [-m <mark1> [<mark2>...]]'
say '         [-z <size1> [<size2>...]] [-c <color1> [<color2>...]] [-t <text1> [<text2>...]] [-T <TEXT1> [<TEXT2>...]]'
say '         [-levs <level1> [<level2>...]] [-levr <level1> [<level2>...]] [-levc <color>] [-append 1]'
say '    <STD>: standard deviation.'
say '    <CORR>: correlation.'
say '    <file>: input file containing standard deviations (in one column) and correlations (in another column).'
say '    <row_start>: first row to read. Default=1.'
say '    <row_end>: last row to read. Default=number of rows in <file>.'
say '    <column_STD>: which column contains standard deviations. Default=1.'
say '    <column_CORR>: which column contains correlations. Default=2.'
say '    <column_TEXT>: which column contains text to be shown in the legend (use "legend.gs").'
say '    <limit>: limit of x/y axis. Default=2.5. A quarter circle is drawn if positive; a half circle is drawn if negative.'
say '    <mark>: defaults to "2 3 4...", i.e., open circle, closed circle, open square, closed square, and so on.'
say '    <size>: mark size. Default=0.11,'
say '    <color>: defaults to "1 2 3...", i.e., foreground color, red, green, dark blue, and so on.'
say '    <text>: text to be shown above each mark. Text beginning with a minus sign or containing spaces must be double quoted.'
say '    <TEXT>: text to be shown in the legend (use "legend.gs"). Text beginning with a minus sign or containing spaces must be double quoted.'
say '    <level>: contour levels to be drawn for <STD> and/or <CORR>. Contour level 1 is drawn by default for standard deviation.'
say '    <color>: line color for <STD> and/or <CORR> contours. Default=1.'
say '    -append 1: use if appending to an existing plot. (Run "legend.gs" only once after all data are plotted.)'
say ''
say '  EXAMPLE 1:'
say '    subplot 1 1 1 -xy 1'
say '    taylor -s 0.8 1.25 -r 0.8 0.9 -t A B -T "Model A" "Model B" -c 2 -m 2'
say '    drawstr -p 6 9 corr -t "Standard Deviation" "Standard Deviation" Correlation'
say '    legend'
say ''
say '  EXAMPLE 2:'
say '    subplot 1 1 1 -xy 2'
say '    taylor -s 0.8 1.25 -r -0.8 -0.9 -l -2.5 -t 1 2 -T "Model A" "Model B" -c 4 -m 3'
say '    drawstr -p 6 corr -t "Standard Deviation" Correlation'
say '    legend'
say ''
say '  EXAMPLE 3 (input given by a text file "myinput.txt", where 1st column gives standard deviations, and 2nd column correlations):'
say '    subplot 1 1 1 -xy 1'
say '    taylor -i myinput.txt -T "Model A" "Model B" -c 2 -m 2'
say '    drawstr -p 6 9 corr -t "Standard Deviation" "Standard Deviation" Correlation'
say '    legend'
say ''
say '  Dependencies: parsestr.gsf, qdims.gsf'
say ''
say '  See also: legend.gs'
say ''
say '  Copyright (C) 2012 Bin Guan.'
say '  Distributed under GNU/GPL.'
return
