stringa='25km.TiKf_comeTaleena'
*stringa2='comeCordex'
*stringa3='exe47mod'
*stringa4='NogherottoMoisture scheme - WET'

year1='1979'
year2='1981'

month.1=jan
month.2=apr
month.3=jul
month.4=oct

title.1='DJF'
title.2='MAM'
title.3='JJA'
title.4='SON'
title.5=year1'-'year2

k=1
kmax=4
while (k<=kmax)
'reinit'
  'set string 2 tl 5 0'
  'set strsiz 0.1'
  'draw string 0.1 10.8  ' stringa
*  'draw string 0.1 10.6  ' stringa2
* 'draw string 0.1 10.4  ' stringa3
*  'draw string 0.1 10.2  ' stringa4

*'open /home/netapp-clima/users/mariotti/50km_regcm4/africa_clm45_sep2017/monthly_EEclm45_pbl1_comeCordex_newmicro/Africa_SRF.nc.ctl'
*'open /home/netapp-clima/users/mariotti/50km_regcm4/africa_clm45_sep2017/monthly_EEclm45_pbl1_comeCordex/Africa_SRF.nc.ctl'
*'open /home/netapp-clima/users/mariotti/50km_regcm4/africa_clm45/monthly_EEclm45_pbl1_comeCordex/Africa_SRF.nc.ctl'
'open /home/clima-archive4/CORDEX2/monthly_original/AFRICA/eraint/Africa_SRF.1979-1981.nc.ctl'

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

'define mean1=ave(ave(s01tas.1,t-1,t+1),time='month.k''year1',time='month.k''year2',1yr)'
'define mean1=mean1-273.15'

'set gxout shaded'
'run subpg p3 2'
'/home/netapp-clima-users1/fraffael/temp_diff_light.gs'
'set clevs -5 0 5 10 12 15 20 25 30'
'set ccols  21 22 23 24 25 26 27 28 29 30 31'
'd mean1'
'draw title RegCM CLM45 TEMP 'title.5' ' title.k
'run /home/netapp-clima/users/mariotti/RegCM4/cbarn.gs'

'close 1'

* TRMM data
'sdfopen /home/esp-shared-a/Observations/CRU/CRU_TS4.00/cru_ts4.00.1901.2015.tmp.dat.nc'
'set lon -40 60'
'set lat -40 40'
'define mean3=ave(ave(tmp.1,t-1,t+1),time='month.k''year1',time='month.k''year2',1yr)'

'run subpg p3 1'
'/home/netapp-clima-users1/fraffael/temp_diff_light.gs'
'set clevs -5 0 5 10 12 15 20 25 30'
'set ccols  21 22 23 24 25 26 27 28 29 30 31'
'd mean3'
'draw title CRU TEMP 'title.5' 'title.k
'run /home/netapp-clima/users/mariotti/RegCM4/cbarn.gs'


* RegCM-TRMM
  'define obstrmm=lterp(mean3,mean1)'
  'define difftrmm=((mean1-obstrmm))'
*  'define percenttrmm=(mean1-obstrmm)*100/(obstrmm+0.01)'
*  'define percenttrmm=maskout(percenttrmm,abs(difftrmm)-0.5)'

'run subpg p3 3'
'/home/netapp-clima-users1/fraffael/temp_diff_light.gs'
'set clevs -6 -4 -2 -1 1 2 4 6'
'set ccols  21 22 23  24  0  27 28 36 37'
'd difftrmm'
'draw title BIAS_CRU 'title.k' (degrees)'
'run /home/netapp-clima/users/mariotti/RegCM4/cbarn.gs'
 

'gxprint africaband-'title.k'.eps white'

'close 1'
k=k+1
endwhile
'!gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=africaband-ERAINT1979-1981.pdf africaband-DJF.eps africaband-MAM.eps africaband-JJA.eps africaband-SON.eps'

'quit'
