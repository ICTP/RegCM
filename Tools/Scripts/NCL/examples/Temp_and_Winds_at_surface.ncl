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
  a = addfile("CAS_SRF.1990010100.nc","r")
  t2m  = a->tas(:,0,:,:)
  u10m = a->uas(:,0,:,:)
  v10m = a->vas(:,0,:,:)

  xlat = a->xlat
  xlon = a->xlon
  dims = dimsizes(xlat)
  nlat = dims(1)                            
  nlon = dims(0)

  tdims = dimsizes(t2m)
  ntimes = tdims(0)

  times = a->time;

;***********************
; plot
;***********************
  wks = gsn_open_wks("ps","RegCM_Temperature_and_Winds")
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
  
  res@cnFillOn               = True
  res@cnLinesOn              = False
  res@cnLineLabelsOn         = False
  res@gsnScalarContour       = True
  res@gsnSpreadColors        = True

  res@mpGeophysicalLineColor = "red"     ; color of continental outlines
  res@mpPerimOn              = True      ; draw box around map
  res@mpFillOn               = False      ; turn off map fill
  res@mpGridLineDashPattern  = 2         ; lat/lon lines as dashed
  
  res@gsnAddCyclic           = False     ; regional data don't add

  res@pmTickMarkDisplayMode  = "Always"  ;
  res@lbAutoManage           = True      ;

  res@vcRefMagnitudeF         = 10.0             ; define vector ref mag
  res@vcRefLengthF            = 0.045            ; define length of vec ref
  res@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
  res@vcMinDistanceF          = 0.017
  res@vcVectorDrawOrder       = "PostDraw"       ; draw vectors last

  res@tiMainFuncCode = "~" ;
  i = 1
  do while (i .lt. ntimes)
    res@tiMainString = "Temperature and winds" +  \
                       "~C~" + ut_string(times(i), "");
    plot = gsn_csm_vector_scalar_map(wks,u10m(i,:,:), \
                                     v10m(i,:,:),t2m(i,:,:),res)
    i = i + 1
  end do

end
