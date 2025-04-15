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

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

begin

  rad_files = systemfunc("ls output/*RAD*.nc")
  nfiles = dimsizes(rad_files)
  do n = 0 , nfiles-1
    rfile = rad_files(n)
    f = addfile(rfile, "rw")
    setfileoption(f,"DefineMode",True)
    filevardef(f, "cldlo", typeof(f->ps), getvardims(f->ps))
    filevarattdef(f, "cldlo", f->cld)
    filevardef(f, "cldmi", typeof(f->ps), getvardims(f->ps))
    filevarattdef(f, "cldmi", f->cld)
    filevardef(f, "cldhi", typeof(f->ps), getvardims(f->ps))
    filevarattdef(f, "cldhi", f->cld)
    filevardef(f, "cldtt", typeof(f->ps), getvardims(f->ps))
    filevarattdef(f, "cldtt", f->cld)
    setfileoption(f,"DefineMode",False)
    f->cldhi(:,:,:) = dim_product_n((1.0-f->cld(:,0:5,:,:)), 1)
    f->cldmi(:,:,:) = dim_product_n((1.0-f->cld(:,5:10,:,:)), 1)
    f->cldlo(:,:,:) = dim_product_n((1.0-f->cld(:,10:17,:,:)), 1)
    f->cldtt(:,:,:) = dim_product_n((1.0-f->cld(:,:,:,:)), 1)
  end do

end
