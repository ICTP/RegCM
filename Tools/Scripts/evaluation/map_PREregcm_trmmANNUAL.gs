stampa.1='/home/netapp-clima-users1/fraffael/subpg p2 '1
stampa.2='/home/netapp-clima-users1/fraffael/subpg p2 '2

'open /home/clima-archive4/CORDEX2/monthly_original/AFRICA/eraint/Africa_SRF.1979-1981avg.nc.ctl'

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

'define meanP=pr*86400'


'set gxout shaded'
'run subpg p2 1'
'/home/netapp-clima-users1/fraffael/mycolors.gs'
*'precip_mmdiff_light.gs'
'set clevs 1 2 4 8 12 16 24'
'set ccols  0 49 50 51 52 53 54 47'
'd meanP'
'draw title RegCM CLM45 PRE 1979-1981'
'run /home/netapp-clima/users/mariotti/RegCM4/cbarn.gs'
'close 1'

* TRMM data
'sdfopen TRMM-1998_2008avg.nc'

'set lon -40 60'
'set lat -40 40'
'define meanT=precip*24'


* RegCM-TRMM
  'define obstrmm=lterp(meanT,meanP)'
  'define difftrmm=((meanP-obstrmm))'
  'define percenttrmm=(meanP-obstrmm)*100/(obstrmm+0.01)'
  'define percenttrmm=maskout(percenttrmm,abs(difftrmm)-0.5)'

'run subpg p2 2'
'/home/netapp-clima-users1/fraffael/precip_mmdiff_light.gs'
*'/home/netapp-clima-users1/fraffael/mycolors.gs'
'set clevs  -16 -8 -4 -2 -1 1 2 4 8 16 '
'set ccols  41 42 43 44 45 20 21 22 23 25 27 29 30'
'd difftrmm'
'draw title BIAS_TRMM (mm/day)'
'run /home/netapp-clima/users/mariotti/RegCM4/cbarn.gs'
 

'gxprint africaband-ERAINTannual.eps white'
'quit'
