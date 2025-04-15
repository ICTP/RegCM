;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
;
;    This file is part of ICTP RegCM.
;    
;    Use of this source code is governed by an MIT-style license that can
;    be found in the LICENSE file or at
;
;         https://opensource.org/licenses/MIT.
;
;    ICTP RegCM is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
;
;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

; XY line plot for RegCM4 output

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

; Set path to input file
  diri = "/home/graziano/test/output/";
  fili = "CAS_ATM.1990010100.nc"

; Pointers to input files
  f  = addfile(diri+fili,"r")

; Read in variables
  var    = f->ta(:,17,:,:)      ; (time,lev,lat,lon)
  var1   = f->time

; Calculate area average temperature over the entire domain
  varave = wgt_areaave(var,1.0,1.0,0)

  wks  = gsn_open_wks ("ps", "xy")            ; open workstation
 
  res                        = True               ; plot mods desired
 
  res@xyLineColors           = (/"blue"/)
  res@xyLineThicknesses      = (/"4.0"/)

  res@tiMainString           = "RegCM4 xy_plot"

  plot = gsn_csm_xy(wks,var1,varave,res)

end
