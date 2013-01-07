* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* mapPre.gs
* 
* This script creates Precipitation's maps of the differences 
* of RegCM4 minus CRU in mm/day and in % for each seasons 
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
  'set clopts -1 -1 0.13'

  'open RegCM4monthly_SRF.1998_2003.nc.ctl'
  'define rcmmean=ave(ave(pr.1,t-1,t+1),time='month.k'1998,time='month.k'2002,1yr)'
  'close 1'

  'sdfopen /home/netapp-clima/shared/OBS/CRU/TS2.1/CRUPRE.CDF'
  'set lon 14.65 120.25'
  'set lat  -20.05 50.13'
  'define crumean=ave(ave(pre.1,t-1,t+1),time='month.k'1998,time='month.k'2002,1yr)'


* RegCM-CRU 
  'define obscru=lterp(crumean,rcmmean)'
  'define diffcru=(rcmmean-obscru)'
  'define percentcru=(rcmmean-obscru)*100/(obscru+0.01)'
  'define percentcru=maskout(percentcru,abs(obscru)-1.0)'

  'run subpg l2 1'
  'run colors.gs'
  'set clevs -10  -5 -2  -1 -0.5 0.5  1  2  5 10'
  'set ccols  21  22  23  24  26  0  30 31 32 33 34 35 '
  'd diffcru'

  'draw title PRE RegCM4 minus CRU mm/day 'title.k
  'run cbarn.gs'


  
  'run subpg l2 2'
  'run colors.gs'
  'set clevs -100 -75 -50 -25 -10 10 25 50 75 100'
  'set ccols   21  22  23  24  26  0 30 31 32 33 34 35 '
  'd percentcru'

  'draw title PRE RegCM4 minus CRU % 'title.k
  'run cbarn.gs'

  'enable print biasPRE_'title.k'.gmf'
  'print'
  'disable print'
  '!gxeps -c -i biasPRE_'title.k'.gmf -o biasPRE_'title.k'.eps'
  'close 1'
k=k+1
endwhile


