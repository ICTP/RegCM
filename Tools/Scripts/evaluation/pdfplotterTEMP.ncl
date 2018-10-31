;function createplot(wks,imdmask,regcmmask,aphmask,namereg)

  function createplot(wks,aphmask,namereg,region)
  local dir50km,dirimd,dir20km,dirtrmm,v1,v2,v3,v4,km50,imd,km20,trmm,m2,mas2,kmres50,imdres,kmres20,trmmres,a50km,aimd,a20km,atrmm,nvar,nbin 
  begin

  fref     = addfile ("temp-REFmaskLANDbox.nc", "r")


    pref  = fref->temp(0:11320,:,:)


  pref@_FillValue = default_fillvalue("float")

   
   
  opt = True
  opt@bin_min = 0.
  opt@bin_max = 45.

;20 e' il numero di classi in cui suddivido i valori di precip e x le quali mi calcolo la frequenza

  aobs = pdfx(pref, 200, opt)
  nvar = 1

;  printVarSummary(aobs)
;  printVarSummary(aemnc)
;  exit()

  nbin = aobs@nbins
  
  xx = new((/nvar, nbin/), typeof(aobs))
  xx(0,:) = aobs@bin_center
  yy = new((/nvar, nbin/), typeof(aobs))
  yy(0,:) = (/aobs/)/100 
  yy@_FillValue = default_fillvalue(typeof(yy))
;  yy = where(yy.le.0.00000001,yy@_FillValue,yy)
  oarry = exp(fspan(-15,1,30))
  arry = 10^(fspan(-8,0,200))
 
  print(aobs@bin_center+" "+aobs) 
  res  = True
  res@gsnDraw   = False                             ; don't draw yet
  res@gsnFrame  = False
  res@trXMinF = 0
  res@trXMaxF = 50
;  res@trYMinF = 0
;  res@trYMaxF = 1
  res@xyMarkLineModes          = (/"Markers","Markers","Markers","Markers","Markers","Markers"/)
  res@xyMarkers                = (/16,16,16,16,16,16/)
  res@xyMonoMarkerSize             = True
  res@xyMarkerSizeF            = 0.007    ; default = 0.01
  res@xyMarkerColors            = (/"blue","green3","red","black","blueviolet","orange"/)
;  res@trYAxisType              = "IrregularAxis"
;  res@xyYIrregularPoints       = arry
 
;  res@tmYLMode                 = "Explicit"
;  res@tmYLValues               = (/0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,1/)
;  res@tmYLLabels               = (/"0.0000001","0.000001","0.00001","0.0001","0.001","0.01","0.1","1"/)
;  res@tmYLValues               = (/0.000000001,0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,1/)
;  res@tmYLLabels               = (/"10~S~-9~N~","10~S~-8~N~","10~S~-7~N~","10~S~-6~N~","10~S~-5~N~","10~S~-4~N~","10~S~-3~N~","10~S~-2~N~","10~S~-1~N~","10~S~0~N~"/)

;  res@pmLegendDisplayMode    = "Always"            ; turn on legend
  res@pmLegendSide           = "Top"               ; Change location of 
  res@pmLegendParallelPosF   = .75                 ; move units right
  res@pmLegendOrthogonalPosF = -0.4                ; move units down
  res@pmLegendWidthF         = 0.1               ; Change width and
  res@pmLegendHeightF        = 0.1                ; height of legend.
  res@lgPerimOn              = True                ; turn off/on box around
  res@lgLabelFontHeightF     = .015                ; label font height
;  res@xyExplicitLegendLabels = (/"1980-2010","2011-2040","2041-2070","2071-2100","TRMM","TDKF","EMNF"/)  ; create explicit labels
  res@xyExplicitLegendLabels = (/"1970-2000"/)  ; create explicit labels
;  res@tiMainString           = namereg + " 200 Bins" ; changed from 100 bins
  res@tiMainString           = "Temperature PDF - MPI - Tropics (Land Only)" ; changed from 100 bins
  res@gsnMaximize = True
  plot = gsn_csm_xy (wks, xx, yy, res)
  return(plot)
  end
  begin

  wks  = gsn_open_wks ("eps","TropicsLAND-TEMP")
  

 ;  listaimd = (/"CNEmsk.nc.66-100-6-39_imd.nc","NEmsk.nc.66-100-6-39_imd.nc","NWmsk.nc.66-100-6-39_imd.nc","PENINmsk.nc.66-100-6-39_imd.nc","WCmsk.nc.66-100-6-39_imd.nc","IMD55kmsk.nc.66-100-6-39_imd.nc"/)
 
 ;  listareg = (/"CNE50kmRegCM.nc","NE50kmRegCM.nc","NW50kmRegCM.nc","PC50kmRegCM.nc","WC50kmRegCM.nc","IMD50kmRegCM.nc"/)

  listaaph = (/"CNEmsk.nc_box_remapnn50","NEmsk.nc_box_remapnn50","NWmsk.nc_box_remapnn50","PENINmsk.nc_box_remapnn50","WCmsk.nc_box_remapnn50","IMD55kmsk.nc_box_remapnn50"/)

;  regionname = (/"NWSA","ESA","EQA","SOA","SWSA","CSA"/) 
  regionname = (/"GLOBAL"/) 

  plot1 = createplot(wks,listaaph(0),regionname(0),1)
;  plot2 = createplot(wks,listaaph(0),regionname(1),2)
;  plot3 = createplot(wks,listaaph(0),regionname(2),3)
;  plot4 = createplot(wks,listaaph(0),regionname(3),4)
;  plot5 = createplot(wks,listaaph(0),regionname(4),5)
;  plot6 = createplot(wks,listaaph(0),regionname(5),6)

 ;---Panel all four plots
 pres                    = True
 pres@gsnMaximize        = True
 pres@gsnPanelLabelBar   = True
 pres@pmLabelBarWidthF   = 0.8
 pres@lbLabelFontHeightF = 0.008
; pres@txString = "PDF Global (MPI) rcp8.5"
 pres@gsnPanelLabelBar = False 
; gsn_panel(wks,(/plot1,plot2,plot3,plot4,plot5,plot6/),(/2,3/),pres)
 gsn_panel(wks,(/plot1/),(/1,1/),pres)
 end
