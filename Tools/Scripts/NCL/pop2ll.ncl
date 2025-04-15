;
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
;
;**************************************************************************
;
; Remap POP gx1v3 grid to regular latlon 1dx1d grid
; POP gx1v3 grid is used in CCSM3.
; To build the map_gx1v3_to_1x1d_aave_da_050401.nc use SCRIP software.
;
; Reference at http://www-pcmdi.llnl.gov/ipcc/model_documentation/CCSM3.htm
;
;**************************************************************************

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

begin

;*************************
; read in data
;*************************
  a = addfile("gx1v3_mn.sst.nc","r")
  
  sst = a->SST

;************************************************
; convert pop to a 1x1 degree grid 
;************************************************

  sst1d = PopLatLon(sst,"gx1v3","1x1d","aave","da","050401") 

;************************************************
; write output
;************************************************ 

  diro = "./" ; output dir
  filo = "ccsm_mn.sst.nc" ; output file name
  system ("/bin/rm -f "+diro+filo) 
  fout = addfile (diro+filo , "c")
  fout->SST = sst1d

end
