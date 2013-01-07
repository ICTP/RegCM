* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* mapPre.gs
* 
* This script creates Precipitation's maps for each seasons 
*
* Written by L.Mariotti, May 2011
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


month.1=jan
month.2=apr
month.3=jul
month.4=oct

title.1='DJF'
title.2='MAM'
title.3='JJA'
title.4='SON'


'open RegCM4monthly_SRF.1998_2003.nc.ctl'


*seasonal loop
k=1
kmax=4
while(k<=kmax)
  'c'
  'set gxout shaded'
  'set grads off'
  'set grid off'
  'set mpdset hires'
  'set map 1 1 8'
  'set xlopts 1 4 0.15'
  'set ylopts 1 4 0.15'
  'set clopts -1 -1 0.16'

  'define mean1=ave(ave(pr.1,t-1,t+1),time='month.k'1998,time='month.k'2003,1yr)'

  'run colors.gs'
  'set clevs  1  2  4  8 12 16 24'
  'set ccols 24 26 30 31 32 33 34 35 '
  
  'd mean1'

  'draw title PRE RegCM4 'title.k
  'run cbarn.gs'

  'enable print mapPRE_'title.k'.gmf'
  'print'
  'disable print'
  '!gxeps -c -i mapPRE_'title.k'.gmf -o mapPRE_'title.k'.eps'
 
k=k+1
endwhile


