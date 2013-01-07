* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* mapPre.gs
* 
* This script creates Temperature's maps of the differences 
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
  'define rcmmean=ave(ave(s01tas.1,t-1,t+1),time='month.k'1998,time='month.k'2002,1yr)-273.15'
  'close 1'

  'sdfopen /home/netapp-clima/shared/OBS/CRU/TS2.1/CRUTMP.CDF'
  'set lon 14.65 120.25'
  'set lat  -20.05 50.13'
  'define crumean=ave(ave(tmp.1,t-1,t+1),time='month.k'1998,time='month.k'2002,1yr)'


* RegCM-CRU 
  'define obscru=lterp(crumean,rcmmean)'
  'define diffcru=(rcmmean-obscru)'
  'define percentcru=(rcmmean-obscru)*100/(obscru+0.01)'
  'define percentcru=maskout(percentcru,abs(obscru)-1.0)'

  'run colors.gs'
  'set clevs -8 -5 -2 -1 -0.5 0.5  1  2  4   8'
  'set ccols 42 43 44 45  46   0   51 52 53 54 55'
  'd diffcru'

  'draw title TMP RegCM4 minus CRU degC 'title.k
  'run cbarn.gs'


  'enable print biasTMP_'title.k'.gmf'
  'print'
  'disable print'
  '!gxeps -c -i biasTMP_'title.k'.gmf -o biasTMP_'title.k'.eps'
  'close 1'
k=k+1
endwhile


