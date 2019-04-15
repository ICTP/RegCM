;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
;
;    This file is part of ICTP RegCM.
;
;    ICTP RegCM is free software: you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation, either version 3 of the License, or
;    (at your option) any later version.
;
;    ICTP RegCM is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
;
;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"    
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/ut_string.ncl"

begin
;*************************
; read in data
;*************************
  a = addfile("CAS_ATM.1990010100.nc","r")
  t = a->ta

  xlat = a->xlat
  xlon = a->xlon
  dims = dimsizes(xlat)
  nlat = dims(1)                            
  nlon = dims(0)

  tdims = dimsizes(t)
  ntimes = tdims(0)
  nsigma = tdims(1)

  times = a->time;

;***********************
; plot
;***********************
  wks = gsn_open_wks("ps","RegCM_Temperature_at_lowest_sigma")
  gsn_define_colormap(wks,"gui_default")
  
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
  clon = a@longitude_of_projection_origin

  if (prj .eq. "LAMCON") then
    trlats = a@standard_parallel
    res@mpProjection        = "LambertConformal"
    res@mpLambertParallel1F = trlats(0)
    res@mpLambertParallel2F = trlats(1)
    res@mpLambertMeridianF  = clon
  end if
  if (prj .eq. "NORMER") then
    res@mpProjection        = "Mercator"
  end if
  if (prj .eq. "POLSTR") then
    clat = a@latitude_of_projection_origin
    res@mpProjection        = "Stereographic"
    res@mpRelativeCenterLon = True
    res@mpCenterLonF = clon
    res@mpRelativeCenterLat = True
    res@mpCenterLatF = clat
  end if
  if (prj .eq. "ROTMER") then
    clat = a@latitude_of_projection_origin
    res@mpProjection        = "Mercator"
    res@mpCenterLonF = clon
    res@mpCenterLatF = clat
  end if
  
; usually, when data is placed onto a map, it is TRANSFORMED to the specified
; projection. Since this model is already on a native lambert conformal grid,
; we want to turn OFF the tranformation.
  
  res@tfDoNDCOverlay         = True      ; do not transform
  res@cnFillOn               = True      ; color plot desired
  res@cnLinesOn              = False	   ; no contour lines
  
  res@mpGeophysicalLineColor = "red"     ; color of continental outlines
  res@mpPerimOn              = True      ; draw box around map
  res@mpGridLineDashPattern  = 2         ; lat/lon lines as dashed
  
  res@gsnAddCyclic           = False     ; regional data don't add

  res@pmTickMarkDisplayMode  = "Always"  ;
  res@lbAutoManage           = True      ;

  res@tiMainFuncCode = "~" ;
  i = 0
  kz = nsigma - 1
  do while (i .lt. ntimes)
    res@tiMainString = a@title + "~C~" + ut_string(times(i), "");
    plot = gsn_csm_contour_map(wks,t(i,kz,:,:),res)
    i = i + 1
  end do

end
