;********************************************************************
; Plot terrain NetCDF output file
;********************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"    
begin
;*************************
; read in precip data (units cm)
;*************************
  a = addfile("EUROPE_DOMAIN000.nc","r")
  topo = a->topo

  xlat = a->xlat
  xlon = a->xlon
  dims = dimsizes(xlat)
  nlat = dims(1)                            
  nlon = dims(0)

;***********************
; plot
;***********************
  wks = gsn_open_wks("x11","Topography")    ; open a workstation
  gsn_define_colormap(wks,"WhViBlGrYeOrReWh")
  
  res            = True                    ; plot mods desired
; !!!!! any plot of data that is on a native grid, must use the "corners"
; method of zooming in on map.
  
  res@mpLimitMode       = "Corners"            ; choose range of map
  res@mpLeftCornerLatF  = xlat(0,0)
  res@mpLeftCornerLonF  = xlon(0,0)
  res@mpRightCornerLatF = xlat(nlon-1,nlat-1)
  res@mpRightCornerLonF = xlon(nlon-1,nlat-1)

; The following 4 pieces of information are REQUIRED to properly display
; data on a native grid.  This data should be specified somewhere in the
; model itself.
  
  prj = a@projection
  trlats = a@standard_parallel
  clon = a@longitude_of_projection_origin

  if (prj .eq. "LAMCON") then
    res@mpProjection        = "LambertConformal"
    res@mpLambertParallel1F = trlats(0)
    res@mpLambertParallel2F = trlats(1)
    res@mpLambertMeridianF  = clon
  end if
  
; usually, when data is placed onto a map, it is TRANSFORMED to the specified
; projection. Since this model is already on a native lambert conformal grid,
; we want to turn OFF the tranformation.
  
  res@tfDoNDCOverlay         = True        ; do not transform
  res@cnFillOn               = True        ; color plot desired
  res@cnLinesOn              = False	   ; no contour lines
  
  res@mpGeophysicalLineColor = "red"       ; color of continental outlines
  res@mpPerimOn              = True        ; draw box around map
  res@mpGridLineDashPattern  = 2           ; lat/lon lines as dashed
  res@cnLevelSelectionMode        = "ManualLevels"       ; manual levels
  res@cnMinLevelValF              = 0.
  res@cnMaxLevelValF              = 3000.
  res@cnLevelSpacingF             = 300.
  
  res@tiMainString           = a@title     ;

  res@gsnAddCyclic           = False       ; regional data don't add
  
  plot = gsn_csm_contour_map(wks,topo,res) 

end
