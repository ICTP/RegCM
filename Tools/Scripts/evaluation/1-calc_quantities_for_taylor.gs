******************************************************************************
* To plot Taylor diagram                                                     *
* 29/07/2015. Csaba Torma, Trieste                                           *
* /home/netapp-clima/users/ctorma/Scripts_scratch/Taylor/CORDEX_paper/Italy  *
******************************************************************************
*Resolution
region='italy'
obspath='/home/netapp-clima-users1/fraffael/ERAI-VALIDATION/OBS/ITALY/REMAP11/int-durFREQ.nc'

* Region: 'region'
lon1='10.0'
lon2='19.0'
lat1='37.0'
lat2='43.4820'



***********************************
model1='CCLM'
model2='EuroRegCM4'
model3='HIRHAM5'
model4='WRF331F'
model5='RACMO22E'
model6='RCA4'
model7='ALADIN'
model8='RegCM4'
model9='GUF-CCLM'
model10='PROMES'

*EURO-CORDEX
filename1='sdii_m_CLMcom-CCLM4-8-17'
filename2='sdii_m_Euro-RegCM4'
filename3='sdii_m_DMI-HIRHAM5'
filename4='sdii_m_IPSL-INERIS-WRF331F'
filename5='sdii_m_KNMI-RACMO22E'
filename6='sdii_m_SMHI-RCA4'

*Med-CORDEX
filename7='sdii_m_CNRM-ALADIN'  
filename8='sdii_m_ICTP-RegCM4' 
filename9='sdii_m_GUF-CCLM'
filename10='sdii_m_UCLM-PROMES'


*EURO-CORDEX       
filename11='sdii_h_CLMcom-CCLM4-8-17'
filename12='sdii_h_Euro-RegCM4'
filename13='sdii_h_DMI-HIRHAM5'      
filename14='sdii_h_IPSL-INERIS-WRF331F'          
filename15='sdii_h_KNMI-RACMO22E'
filename16='sdii_h_SMHI-RCA4'

*Med-CORDEX
filename17='sdii_h_CNRM-ALADIN' 
filename18='sdii_h_ICTP-RegCM4'
filename19='sdii_h_GUF-CCLM'
filename20='sdii_h_UCLM-PROMES'
filename21='sdii_m_RCM' 
filename22='sdii_h_RCM' 
filename23='sdii_m_medRCM' 
filename24='sdii_h_medRCM'
filename25='sdii_ERAINT'


  'reinit'

***************************************************************
*Reading data from OBS, ERA-Interim and RCMs on 0.44 and 0.11 *
***************************************************************


*1  Italy
  'sdfopen 'obspath''
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
    'define meanobs=pre'
  'close 1'

* 0.44 CLMcom-CCLM4-8-17
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/CCLM_044_sdii.nc'
*  'define dum=lterp(pr,meanobs)'
*  'define dur=maskout(dum,meanobs)'
*  'define mean1annual=dur'

   'define a=lterp(pr,meanobs)'
   'define b=maskout(a,meanobs)'
   'define meanobs=maskout(meanobs,b)'
   'define mean1annual=b'

  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.44 0.44 Euro-RegCM4
* 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/EuroRegCM4_044_sdii.nc'
*  'define dum=lterp(pr,meanobs)'
*  'define dur=maskout(dum,meanobs)'
*  'define mean2annual=dur'
*  'set LON 'lon1' 'lon2''
*  'set LAT 'lat1' 'lat2''
*  'close 1'

* 0.44 DMI-HIRHAM5
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/HIRHAM5_044_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean3annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.44 IPSL-INERIS-WRF331F
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/WRF331F_044_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean4annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.44 KNMI-RACMO22E
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/RACMO22E_044_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean5annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.44 SMHI-RCA4
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/RCA4_044_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean6annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.44 CNRM-ALADIN
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/ALADIN_044_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean7annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.44 ICTP-RegCM4
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/RegCM4_044_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean8annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.44 GUF-CCLM
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/GUF-CCLM_044_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean9annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.44 UCLM-PROMES
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/PROMES_044_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean10annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.11 CLMcom-CCLM4-8-17
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/CCLM_011_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean11annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.11 Euro-RegCM4
* 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/EuroRegCM4_011_sdii.nc'
*  'define dum=lterp(pr,meanobs)'
*  'define dur=maskout(dum,meanobs)'
*  'define mean12annual=dur'
*  'set LON 'lon1' 'lon2''
*  'set LAT 'lat1' 'lat2''
*  'close 1'

* 0.11 DMI-HIRHAM5
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/HIRHAM5_011_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean13annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.11 IPSL-INERIS-WRF331F
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/WRF331F_011_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean14annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.11 KNMI-RACMO22E
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/RACMO22E_011_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean15annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.11 SMHI-RCA4
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/RCA4_011_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean16annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.11 CNRM-ALADIN
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/ALADIN_011_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean17annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.11 ICTP-RegCM4
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/RegCM4_011_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean18annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.11 GUF-CCLM
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/GUF-CCLM_011_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean19annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.11 UCLM-PROMES
 'sdfopen /home/netapp-clima/users/ctorma/Scripts/SDII/'region'/PROMES_011_sdii.nc'
  'define dum=lterp(pr,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean20annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''
  'close 1'

* 0.11 ERA-Interim
 'sdfopen /home/netapp-clima-users1/fraffael/ERAI-VALIDATION/ERAINT/REMAP11/int-durFREQ_'region'.nc'
  'define dum=lterp(tp,meanobs)'
  'define dur=maskout(dum,meanobs)'
  'define mean25annual=dur'
  'set LON 'lon1' 'lon2''
  'set LAT 'lat1' 'lat2''


*50 km RCMs mean
  'define mean21annual=(mean1annual+mean3annual+mean4annual+mean5annual+mean6annual+mean7annual+mean8annual+mean9annual+mean10annual)/9'


*12 km RCMs mean
  'define mean22annual=(mean11annual+mean13annual+mean14annual+mean15annual+mean16annual+mean17annual+mean18annual+mean19annual+mean20annual)/9'


*50 km RCMs mean - Only Med-CORDEX
  'define mean23annual=(mean7annual+mean8annual+mean9annual+mean10annual)/4'


*12 km RCMs mean - Only Med-CORDEX
  'define mean24annual=(mean17annual+mean18annual+mean19annual+mean20annual)/4'



********************************************************
* Calculating spatial correlation and write into files *
********************************************************
* 1 0.44 CLMcom-CCLM4-8-17 -> OBS

  'set gxout stat'

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean1annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean1annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename1, value,append)
  res=close (filename1)

* 2 0.44 DMI-HIRHAM5 -> OBS

*  'set gxout stat'

*  'd meanobs'
*  lineacru = sublin(result,13)
*  sigmacru = subwrd(lineacru,2)
*  'd mean2annual'
*  linea = sublin(result,13)
*  sigma = subwrd(linea,2)
*  sigmaannual=sigma/sigmacru
*  'd scorr(mean2annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
*  linea2 = sublin(result,8)
*  value2 = subwrd(linea2,4)
*  value = sigmaannual ' ' value2
*  res = write(filename2, value,append)
*  res=close (filename2)


****************************
* 3 0.44 IPSL-INERIS-WRF331F -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean3annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean3annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename3, value,append)
  res=close (filename3)

****************************
* 4 0.44 KNMI-RACMO22E -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean4annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean4annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename4, value,append)
  res=close (filename4)

****************************
* 5 0.44 SMHI-RCA4 -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean5annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean5annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename5, value,append)
  res=close (filename5)


****************************
* 6 0.44 CNRM-ALADIN -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean6annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean6annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename6, value,append)
  res=close (filename6)

****************************
* 7 0.44 ICTP-RegCM4 -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean7annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean7annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename7, value,append)
  res=close (filename7)

****************************
* 8 0.44 GUF-CCLM -> OBS


  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean8annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean8annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename8, value,append)
  res=close (filename8)

****************************
* 9 0.44 UCLM-PROMES -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean9annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean9annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename9, value,append)
  res=close (filename9)

****************************
* 10 0.11 CLMcom-CCLM4-8-17 -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean10annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean10annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename10, value,append)
  res=close (filename10)

****************************
* 11 0.11 DMI-HIRHAM5 -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean11annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean11annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename11, value,append)
  res=close (filename11)

****************************
* 12 0.11 IPSL-INERIS-WRF331F -> OBS

*  'd meanobs'
*  lineacru = sublin(result,13)
*  sigmacru = subwrd(lineacru,2)
*  'd mean12annual'
*  linea = sublin(result,13)
*  sigma = subwrd(linea,2)
*  sigmaannual=sigma/sigmacru
*  'd scorr(mean12annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
*  linea2 = sublin(result,8)
*  value2 = subwrd(linea2,4)
*  value = sigmaannual ' ' value2
*  res = write(filename12, value,append)
*  res=close (filename12)

****************************
* 13 0.11 KNMI-RACMO22E -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean13annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean13annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename13, value,append)
  res=close (filename13)

****************************
* 14 0.11 SMHI-RCA4 -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean14annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean14annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename14, value,append)
  res=close (filename14)

****************************
* 15 0.11 CNRM-ALADIN -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean15annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean15annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename15, value,append)
  res=close (filename15)

****************************
* 16 0.11 ICTP-RegCM4 -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean16annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean16annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename16, value,append)
  res=close (filename16)

****************************
* 17 0.11 GUF-CCLM -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean17annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean17annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename17, value,append)
  res=close (filename17)

****************************
* 18 0.11 UCLM-PROMES -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean18annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean18annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename18, value,append)
  res=close (filename18)

****************************
* 19 0.44 RCMs -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean19annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean19annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename19, value,append)
  res=close (filename19)

****************************
* 20 0.11 RCMs -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean20annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean20annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename20, value,append)
  res=close (filename20)

****************************
* 21 ERA-Int -> OBS

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean21annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean21annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename21, value,append)
  res=close (filename21)

****************************

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean22annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean22annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename22, value,append)
  res=close (filename22)

****************************

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean23annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean23annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename23, value,append)
  res=close (filename23)

****************************

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean24annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean24annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename24, value,append)
  res=close (filename24)

****************************

  'd meanobs'
  lineacru = sublin(result,13)
  sigmacru = subwrd(lineacru,2)
  'd mean25annual'
  linea = sublin(result,13)
  sigma = subwrd(linea,2)
  sigmaannual=sigma/sigmacru
  'd scorr(mean25annual,meanobs,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
  linea2 = sublin(result,8)
  value2 = subwrd(linea2,4)
  value = sigmaannual ' ' value2
  res = write(filename25, value,append)
  res=close (filename25)


