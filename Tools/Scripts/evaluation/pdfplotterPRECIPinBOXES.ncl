;function createplot(wks,imdmask,regcmmask,aphmask,namereg)

  function createplot(wks,aphmask,namereg,region)
  local dir50km,dirimd,dir20km,dirtrmm,v1,v2,v3,v4,km50,imd,km20,trmm,m2,mas2,kmres50,imdres,kmres20,trmmres,a50km,aimd,a20km,atrmm,nvar,nbin 
  begin

  fobs     = addfile ("TRMM.day.mean.1998_2009_BOXafrica.nc", "r")
  femn     = addfile ("Africa_STS-prOK_1998-2002-remap.nc", "r")
  femnc     = addfile ("Africa_SRF-prOK_1998-2002-remap.nc", "r")
  m1      = addfile("mask25km-remap.nc", "r")

  if (region.eq.1) then
    pobs1  = fobs->hrf(:1095,{10:20},{-18:8})
    pemn1  = femn->pr(:1095,{10:20},{-18:8})
    pemnc1  = femnc->prc(:1095,{10:20},{-18:8})
    mas1   = m1->mask({10:20},{-18:8})
  else if (region.eq.2) then
    pobs1  = fobs->hrf(:1095,{5:20},{30:50})
    pemn1  = femn->pr(:1095,{5:20},{30:50})
    pemnc1  = femnc->prc(:1095,{5:20},{30:50})
    mas1   = m1->mask({5:20},{30:50})
  else if (region.eq.3) then
    pobs1  = fobs->hrf(:1095,{-15:5},{8:42})
    pemn1  = femn->pr(:1095,{-15:5},{8:42})
    pemnc1  = femnc->prc(:1095,{-15:5},{8:42})
    mas1   = m1->mask({-15:5},{8:42})
  else if (region.eq.4) then
    pobs1  = fobs->hrf(:1095,{-35:-15},{10:40})
    pemn1  = femn->pr(:1095,{-35:-15},{10:40})
    pemnc1  = femnc->prc(:1095,{-35:-15},{10:40})
    mas1   = m1->mask({-35:-15},{10:40})
  else if (region.eq.5) then
    pobs1  = fobs->hrf(:1095,{5:10},{-18:8})
    pemn1  = femn->pr(:1095,{5:10},{-18:8})
    pemnc1  = femnc->prc(:1095,{5:10},{-18:8})
    mas1   = m1->mask({5:10},{-18:8})
  else
    pobs1  = fobs->hrf(:1095,{5:20},{8:30})
    pemn1  = femn->pr(:1095,{5:20},{8:30})
    pemnc1  = femnc->prc(:1095,{5:20},{8:30})
    mas1   = m1->mask({5:20},{8:30})
  end if
  end if
  end if
  end if
  end if

;  imd1@_FillValue = default_fillvalue("float")
;  trmm1@_FillValue = default_fillvalue("float")
;
  pobs1@_FillValue = default_fillvalue("float")
  pobs1 = where(pobs1.lt.1, default_fillvalue("float"), pobs1)
  pemn1 = where(pemn1.lt.1, default_fillvalue("float"), pemn1)
  pemnc1 = where(pemnc1.lt.1, default_fillvalue("float"), pemnc1)
;  km50 = where(km501.lt.1, default_fillvalue("float"), km501)
;  imd = where(imd1.lt.1, default_fillvalue("float"), imd1)
;  km20 = where(km201.lt.1, default_fillvalue("float"), km201)
;  trmm = where(trmm1.lt.1, default_fillvalue("float"), trmm1)

   
 ;m2     = addfile("/home/netapp-clima-users1/sdas/old_20km_50km/daily_20km/new_20km_lmsk/new20km_mask/"+aphmask, "r")

  
 ;  if (namereg .eq. "All Indian regions") then
 ;  kmres50     = mask(km50, mas2(:,:), 1)
 ;  imdres     = mask(imd, mas2(:,:), 1)
 ;  kmres20    = mask(km20, mas2(:,:), 1)
 ;  trmmres    = mask(trmm, mas2(:,:),1)
 ;  else
 
   
  pemn_m    = mask(pemn1, mas1(:,:),2)
  pemnc_m    = mask(pemnc1, mas1(:,:),2)
  pobs_m    = mask(pobs1, mas1(:,:),2)
 ;  end if


  opt = True
  opt@bin_min = 0.
  opt@bin_max = 800.

  aemn= pdfx(pemn_m, 200, opt)
  aemnc = pdfx(pemnc_m, 200, opt)
  aobs = pdfx(pobs_m, 200, opt)
  nvar = 3

  nbin = aemn@nbins
;  print(a50km@bin_bound_min)
;  print(a50km@bin_center)
;  print(a50km@bin_center)
;  print(a50km)
  
  xx = new((/nvar, nbin/), typeof(aemn))
  xx(0,:) = aemn@bin_center
  xx(1,:) = aemnc@bin_center
  xx(2,:) = aobs@bin_center
;  xx@_FillValue = default_fillvalue(typeof(xx))
;  xx = where(xx.le.0.000001,xx@_FillValue,xx)
  yy = new((/nvar, nbin/), typeof(aemn))
  yy(0,:) = (/aemn/)/100 
  yy(1,:) = (/aemnc/)/100 
  yy(2,:) = (/aobs/)/100 
  yy@_FillValue = default_fillvalue(typeof(yy))
  yy = where(yy.le.0.00000001,yy@_FillValue,yy)
  oarry = exp(fspan(-15,1,30))
  arry = 10^(fspan(-8,0,200))
;  arry = 10^(fspan(-7,0,100))
;  print(oarry)
;  print(arry)
  
  res  = True
  res@gsnDraw   = False                             ; don't draw yet
  res@gsnFrame  = False
  res@trXMinF = 0
  res@trXMaxF = 900
  res@trYMinF = 0
  res@trYMaxF = 1
  res@xyMarkLineModes          = (/"Markers","Markers","Markers","Markers","Markers","Markers"/)
  res@xyMarkers                = (/16,16,16,16,16,16/)
  res@xyMarkerColors            = (/"blue","red","green","black","blueviolet","orange"/)
  res@trYAxisType              = "IrregularAxis"
  res@xyYIrregularPoints       = arry
 
  res@tmYLMode                 = "Explicit"
;  res@tmYLValues               = (/0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,1/)
;  res@tmYLLabels               = (/"0.0000001","0.000001","0.00001","0.0001","0.001","0.01","0.1","1"/)
  res@tmYLValues               = (/0.000000001,0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,1/)
  res@tmYLLabels               = (/"10~S~-9~N~","10~S~-8~N~","10~S~-7~N~","10~S~-6~N~","10~S~-5~N~","10~S~-4~N~","10~S~-3~N~","10~S~-2~N~","10~S~-1~N~","10~S~0~N~"/)

  res@pmLegendDisplayMode    = "Always"            ; turn on legend
  res@pmLegendSide           = "Top"               ; Change location of 
  res@pmLegendParallelPosF   = .75                 ; move units right
  res@pmLegendOrthogonalPosF = -0.4                ; move units down
  res@pmLegendWidthF         = 0.125               ; Change width and
  res@pmLegendHeightF        = 0.15                ; height of legend.
  res@lgPerimOn              = True                ; turn off/on box around
  res@lgLabelFontHeightF     = .015                ; label font height
  res@xyExplicitLegendLabels = (/"RegCM PR","RegCM PRC","TRMM 1998-2002","R20km","TRMM","TDKF","EMNF"/)  ; create explicit labels
  res@tiMainString           = namereg + " 200 Bins" ; changed from 100 bins
  res@gsnMaximize = True
  plot = gsn_csm_xy (wks, xx, yy, res)
  return(plot)
  end
  begin

  wks  = gsn_open_wks ("eps","TiKf_25km-comeTaleena-sub8mod-1998-2002")
  

 ;  listaimd = (/"CNEmsk.nc.66-100-6-39_imd.nc","NEmsk.nc.66-100-6-39_imd.nc","NWmsk.nc.66-100-6-39_imd.nc","PENINmsk.nc.66-100-6-39_imd.nc","WCmsk.nc.66-100-6-39_imd.nc","IMD55kmsk.nc.66-100-6-39_imd.nc"/)
 
 ;  listareg = (/"CNE50kmRegCM.nc","NE50kmRegCM.nc","NW50kmRegCM.nc","PC50kmRegCM.nc","WC50kmRegCM.nc","IMD50kmRegCM.nc"/)

  listaaph = (/"CNEmsk.nc_box_remapnn50","NEmsk.nc_box_remapnn50","NWmsk.nc_box_remapnn50","PENINmsk.nc_box_remapnn50","WCmsk.nc_box_remapnn50","IMD55kmsk.nc_box_remapnn50"/)

  regionname = (/"NWSA","ESA","EQA","SOA","SWSA","CSA"/) 

 ;  plot1 = createplot(wks,listaimd(0),listareg(0),listaaph(0),regionname(0))
 ;  plot2 = createplot(wks,listaimd(1),listareg(1),listaaph(1),regionname(1))
 ;  plot3 = createplot(wks,listaimd(2),listareg(2),listaaph(2),regionname(2))
 ;  plot4 = createplot(wks,listaimd(3),listareg(3),listaaph(3),regionname(3))
 ;  plot5 = createplot(wks,listaimd(4),listareg(4),listaaph(4),regionname(4))
 ;  plot6 = createplot(wks,listaimd(5),listareg(5),listaaph(5),regionname(5))

  plot1 = createplot(wks,listaaph(0),regionname(0),1)
  plot2 = createplot(wks,listaaph(0),regionname(1),2)
  plot3 = createplot(wks,listaaph(0),regionname(2),3)
  plot4 = createplot(wks,listaaph(0),regionname(3),4)
  plot5 = createplot(wks,listaaph(0),regionname(4),5)
  plot6 = createplot(wks,listaaph(0),regionname(5),6)

 ;---Panel all four plots
 pres                    = True
 pres@gsnMaximize        = True
 pres@gsnPanelLabelBar   = True
 pres@pmLabelBarWidthF   = 0.8
 pres@lbLabelFontHeightF = 0.008
 pres@txString = "PDF TiKf comeTaleena sub8mod"
 ; pres@txString           = "chi200_ud_smooth.nc"
 pres@gsnPanelLabelBar = False 
gsn_panel(wks,(/plot1,plot2,plot3,plot4,plot5,plot6/),(/2,3/),pres)
; gsn_panel(wks,(/plot1/),(/1,1/),pres)
 end
