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
