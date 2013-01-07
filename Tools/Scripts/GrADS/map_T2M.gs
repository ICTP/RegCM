* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* mapPre.gs
* 
* This script creates 2meters Temperature's maps for each seasons 
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

  'define mean1=ave(ave(s01tas.1,t-1,t+1),time='month.k'1998,time='month.k'2003,1yr)-273.15'

  'run colors.gs'
  'set clevs  1  4  8 12 16 20 24 28 32 '
  'set ccols 41 42 43 44 45 46 51 52 53 54 55'

  'd mean1'

  'draw title TMP RegCM4 'title.k
  'run cbarn.gs'

  'enable print mapTMP_'title.k'.gmf'
  'print'
  'disable print'
  '!gxeps -c -i mapTMP_'title.k'.gmf -o mapTMP_'title.k'.eps'
 
k=k+1
endwhile


