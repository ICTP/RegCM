stampa.1='/home/netapp-clima-users1/fraffael/subpg p6 '1
stampa.2='/home/netapp-clima-users1/fraffael/subpg p6 '2
stampa.3='/home/netapp-clima-users1/fraffael/subpg p6 '3
stampa.4='/home/netapp-clima-users1/fraffael/subpg p6 '4
stampa.5='/home/netapp-clima-users1/fraffael/subpg p6 '5
stampa.6='/home/netapp-clima-users1/fraffael/subpg p6 '6
stampa.7='/home/netapp-clima-users1/fraffael/subpg l8 '7
stampa.8='/home/netapp-clima-users1/fraffael/subpg l8 '8

  'reinit'
*TRMM-monthly1998-2000ymon.nc = cdo ymonmean TRMM-monthly1998-2000.nc TRMM-monthly1998-2000ymon.nc
 'sdfopen TRMM-monthly1998-2000ymon.nc'
 'set t 1 12'
 'define trmm=precip*24'
 'define trmmNWSA=tloop(aave(trmm,lon=-18,lon=8,lat=10,lat=20))'
 'define trmmESA=tloop(aave(trmm,lon=30,lon=50,lat=5,lat=20))'
 'define trmmEQA=tloop(aave(trmm,lon=8,lon=42,lat=-15,lat=5))'
 'define trmmSOA=tloop(aave(trmm,lon=10,lon=40,lat=-35,lat=-15))'
 'define trmmSWSA=tloop(aave(trmm,lon=-18,lon=8,lat=5,lat=10))'
 'define trmmCSA=tloop(aave(trmm,lon=8,lon=30,lat=5,lat=20))'
'close 1'


  'open Africa_STS-prOK_1998-2000-comeTaleena-ymon.nc.ctl'
 'set t 1 12'
 'define modNWSA=tloop(aave(pr,lon=-18,lon=8,lat=10,lat=20))'
 'define modESA=tloop(aave(pr,lon=30,lon=50,lat=5,lat=20))'
 'define modEQA=tloop(aave(pr,lon=8,lon=42,lat=-15,lat=5))'
 'define modSOA=tloop(aave(pr,lon=10,lon=40,lat=-35,lat=-15))'
 'define modSWSA=tloop(aave(pr,lon=-18,lon=8,lat=5,lat=10))'
 'define modCSA=tloop(aave(pr,lon=8,lon=30,lat=5,lat=20))'

* 'close 1'


  'run 'stampa.1
  'set x 1'
  'set y 1'
  'set t 1 12'
  'set vrange 0 10'
  'd modNWSA'
  'd trmmNWSA'
*  'd udelNWSA'
  'draw title Precip Annual Cycle (NWSA)'
  'set string 1 bl 5 0'
  'set strsiz 0.2'
  'draw string 3 6 RegCM'
  'set string 3 bl 5 0'
  'draw string 3 5.6 TRMM'
*  'set string 7 bl 5 0'
*  'draw string 3 5.2 UDEL'




  'run 'stampa.2
  'set x 1'
  'set y 1'
  'set t 1 12'
  'set vrange 0 10'
  'd modSOA'
  'd trmmSOA'
*  'd udelSOA'
  'draw title Precip Annual Cycle (SOA)'
  'set string 1 bl 5 0'
  'set strsiz 0.2'
  'draw string 3 6 RegCM'
  'set string 3 bl 5 0'
  'draw string 3 5.6 TRMM'
*  'set string 7 bl 5 0'
*  'draw string 3 5.2 UDEL'


  'run 'stampa.3
  'set x 1'
  'set y 1'
  'set t 1 12'
  'set vrange 0 10'
  'd modESA'
  'd trmmESA'
*  'd udelESA'
  'draw title Precip Annual Cycle (ESA)'
  'set string 1 bl 5 0'
  'set strsiz 0.2'
  'draw string 3 6 RegCM'
  'set string 3 bl 5 0'
  'draw string 3 5.6 TRMM'
*  'set string 7 bl 5 0'
*  'draw string 3 5.2 UDEL'



  'run 'stampa.4
  'set x 1'
  'set y 1'
  'set t 1 12'
  'set vrange 0 10'
  'd modSWSA'
  'd trmmSWSA'
*  'd udelSWSA'
  'draw title Precip Annual Cycle (SWSA)'
  'set string 1 bl 5 0'
  'set strsiz 0.2'
  'draw string 3 6 RegCM'
  'set string 3 bl 5 0'
  'draw string 3 5.6 TRMM'
*  'set string 7 bl 5 0'
*  'draw string 3 5.2 UDEL'



  'run 'stampa.5
  'set x 1'
  'set y 1'
  'set t 1 12'
  'set vrange 0 10'
  'd modEQA'
  'd trmmEQA'
*  'd udelEQA'
  'draw title Precip Annual Cycle (EQA)'
  'set string 1 bl 5 0'
  'set strsiz 0.2'
  'draw string 3 6 RegCM'
  'set string 3 bl 5 0'
  'draw string 3 5.6 TRMM'
*  'set string 7 bl 5 0'
*  'draw string 3 5.2 UDEL'



  'run 'stampa.6
  'set x 1'
  'set y 1'
  'set t 1 12'
  'set vrange 0 10'
  'd modCSA'
  'd trmmCSA'
*  'd udelCSA'
  'draw title Precip Annual Cycle (CSA)'
  'set string 1 bl 5 0'
  'set strsiz 0.2'
  'draw string 3 6 RegCM'
  'set string 3 bl 5 0'
  'draw string 3 5.6 TRMM'
*  'set string 7 bl 5 0'
*  'draw string 3 5.2 UDEL'


 'gxprint annualcyWestAfricaProva.eps white'
  '!ps2pdf -DEPSCrop annualcyWestAfricaProva.eps annualcyWestAfricaProva.pdf'
  '!convert -density 600 -rotate 90 annualcyWestAfricaProva.eps annualcyWestAfricaProva.jpg'



