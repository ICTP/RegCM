'reinit'

'open /home/clima-archive4/CORDEX2/monthly_original/AFRICA/eraint/Africa_SRF.1979-1981avg.nc.ctl'
*'open /home/netapp-clima/users/mariotti/25km_regcm4/africa_fra/monthly_comeTaleena_subex8mod-daMarconi/Africa_SRF.1998-2002avg.nc.ctl'

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

'define mean1=s01tas-273.15'

'set gxout shaded'
'run subpg p2 1'
'temp_diff_light.gs'
'set clevs -5 0 5 10 12 15 20 25 30'
'set ccols  21 22 23 25 24 26 27 28 29 30 31'
'd mean1'
'draw title RegCM CLM45 TEMP 1979-1981'
'run cbarn.gs'

'close 1'

* TRMM data
'sdfopen CRU-1979-1981avg.nc'
'set lon -40 60'
'set lat -40 40'
'define mean3=tmp'


* RegCM-TRMM
  'define obstrmm=lterp(mean3,mean1)'
  'define difftrmm=((mean1-obstrmm))'

'run subpg p2 2'
'temp_diff_light.gs'
'set clevs -6 -4 -2 -1 1 2 4 6'
'set ccols  21 22 23  24  0  27 28 36 37'
'd difftrmm'
'draw title BIAS_CRU (degrees)'
'cbarn.gs'
 

'gxprint africaband-ERAINTannualT.eps white'
'quit'
