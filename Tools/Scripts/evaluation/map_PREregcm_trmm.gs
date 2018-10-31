stringa='comeTaleena'
stringa2='gulland=0.80; cevaplnd=0.1e-4'
stringa3='cevapoce=0.1e-4; qck1land=0.1e-4; qck1oce=0.5e-3'
stringa4='rprc_lnd=6e-2'

year1='1998'
year2='2002'

month.1=jan
month.2=apr
month.3=jul
month.4=oct
month.5=feb

title.1='DJF'
title.2='MAM'
title.3='JJA'
title.4='SON'
title.5=year1'-'year2

k=1
kmax=4
while (k<=kmax)
'reinit'
  'set string 2 l 5 90'
  'set strsiz 0.1'
  'draw string 0.1 5.8  ' stringa
  'draw string 0.3 5.6  ' stringa2
  'draw string 0.5 5.4  ' stringa3
  'draw string 0.7 5.2  ' stringa4

'open /home/netapp-clima/users/mariotti/25km_regcm4/africa_fra/monthly_comeTaleena_subex8mod-daMarconi/Africa_SRF.1998-2002.nc.ctl'

'set display color white'
'set grads off'
'set grid off'
'set map 1 1 8'
'set xlopts 1 5 0.15'
'set ylopts 1 5 0.15'
'set clopts -1 -1 0.15'
'set lon -40 60'
'set lat -40 40'
*'set mpdset hires'

'define mean1=ave(ave(pr.1,t-1,t+1),time='month.k''year1',time='month.k''year2',1yr)*86400'

'set gxout shaded'
'run subpg p3 2'
'/home/netapp-clima-users1/fraffael/mycolors.gs'
*'precip_mmdiff_light.gs'
'set clevs 1 2 4 8 12 16 24'
'set ccols  0 49 50 51 52 53 54 47'
'd mean1'
'draw title RegCM CLM45 PRE 'title.5' ' title.k
'run /home/netapp-clima/users/mariotti/RegCM4/cbarn.gs'

*'set sdfwrite TiKfclm45_pbl1-'title.k'.nc'
*'sdfwrite mean1'
'close 1'

* TRMM data
'sdfopen /home/esp-shared-a/Observations/TRMM/MonthlyTRMM_3B43/3B43.1998_2008.nc'

'set lon -40 60'
'set lat -40 40'
'define mean3=ave(ave(precip.1*24,t-1,t+1),time='month.k''year1',time='month.k''year2',1yr)'

'run subpg p3 1'
*'precip_mmdiff_light.gs'
'/home/netapp-clima-users1/fraffael/mycolors.gs'
'set clevs 1 2 4 8 12 16 24'
'set ccols  0 49 50 51 52 53 54 47'
'd mean3'
'draw title TRMM PRE 'title.5' 'title.k
'run /home/netapp-clima/users/mariotti/RegCM4/cbarn.gs'


* RegCM-TRMM
  'define obstrmm=lterp(mean3,mean1)'
  'define difftrmm=((mean1-obstrmm))'
  'define percenttrmm=(mean1-obstrmm)*100/(obstrmm+0.01)'
  'define percenttrmm=maskout(percenttrmm,abs(difftrmm)-0.5)'

'run subpg p3 3'
'/home/netapp-clima-users1/fraffael/precip_mmdiff_light.gs'
*'/home/netapp-clima-users1/fraffael/mycolors.gs'
'set clevs  -16 -8 -4 -2 -1 1 2 4 8 16 '
'set ccols  41 42 43 44 45 20 21 22 23 25 27 29 30'
'd difftrmm'
'draw title BIAS_TRMM 'title.k' (mm/day)'
'run /home/netapp-clima/users/mariotti/RegCM4/cbarn.gs'
 

'gxprint africaband-'title.k'.eps white'

'close 1'
k=k+1
endwhile
'!gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=africaband-comeTaleena-sub8mod-1998-2002.pdf africaband-DJF.eps africaband-MAM.eps africaband-JJA.eps africaband-SON.eps'

'quit'
