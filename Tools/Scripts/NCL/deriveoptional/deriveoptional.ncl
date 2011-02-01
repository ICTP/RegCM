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
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$REGCMSRCDIR/Tools/Scripts/NCL/deriveoptional/calculateFog.ncl"

begin

  calcvar = True
  calcvar@Pressure = True
  calcvar@GeoHeight = True
  calcvar@Density = True
  calcvar@Fog = True
  calcvar@Theta = True

  bWriteOut = True

; Constant definitions

  ep1 = 0.608
  ep2 = 0.62197
  rgas = 287.0058
  gti = 9.080665
  rwat = 461.90
  cpd = 1005.46
  cpv = 1869.46
  cpw = 4186.95
  wlhf = 0.3336e6
  wlhv = 2.51040e6

  p0exner = 1000.0
  rovcp = rgas/cpd



  fhatm = addfile(atmfilename,"r")
  fhsrf = addfile(srffilename,"r")
  fhrad = addfile(radfilename,"r")
  fhopt = addfile(optfilename,"w")


  ;Run through dependencies
  ;TODO: Put this in a function
;Functions that depend on a level 2 variable at highest
;Functions that depend on a level 1 variable at highest
  if(calcvar@Fog)then
    calcvar@Geoheight = True
  end if
;Functions that only depend on a level 0 variable
  if(calcvar@Theta)then
    calcvar@Pressure = True
  end if
  if(calcvar@Density)then
    calcvar@Pressure = True
  end if
  if(calcvar@GeoHeight)then
    calcvar@Pressure = True
  end if


  if(calcvar@Pressure)then
    if(bWriteOut)then
      print("Calculating Pressure...")
    end if
    ps    = fhatm->ps
    sigma = fhatm->sigma
    ptop  = fhatm->ptop
    time = fhatm->time

    t = fhatm->t

    ;Create the pressure variable
    dim4d = dimsizes(t)
    p = t

    ;Create a temporary 4d surface pressure variable
    ;and a 4d sigma variable
    ps4d = conform_dims(dim4d,ps,(/0,2,3/))
    sig4d = conform_dims(dim4d,sigma,(/1/))

    p = (sig4d*(ps4d-ptop) + ptop)

    delete(ps4d)
    delete(sig4d)


    p@long_name = "Pressure"
    p@units = "hPa"

    fhopt->pres=p
    if(bWriteOut)then
      print("...Finished Writing Pressure")
    end if
  end if

  if(calcvar@Density)then
    if(bWriteOut)then
      print("Calculating Density...")
    end if
    t = fhatm->t
    qv = fhatm->qv
    qc = fhatm->qc

    tv = t*(1 + ep1*qv - qc)

    rho = t
    rho = 100.0*p/(rgas*tv)
    rho@long_name = "Density"
    rho@standard_name = "density"
    rho@units = "kg m-3"

    fhopt->rho = rho

    if(bWriteOut)then
      print("...Finished Writing Density")
    end if
    delete(tv)
  end if

  if(calcvar@GeoHeight)then
    if(bWriteOut)then
      print("Calculating Geopotential Height...")
    end if
    ps    = fhatm->ps
    topo = fhatm->topo
    t = fhatm->t
    qv = fhatm->qv
    qc = fhatm->qc
    ts = fhsrf->t2m(:,0,:,:)
    qvs = fhsrf->q2m(:,0,:,:)
    qcs = qvs
    qcs = 0.0

    dim4d = dimsizes(t)
    nlev = dim4d(1)

    ;Make temporary variables that merge the surface data in
    tfull = merge_levels_sfc(t,ts,1)
    qvfull = merge_levels_sfc(qv,qvs,1)
    qcfull = merge_levels_sfc(qc,qcs,1)
    pfull = merge_levels_sfc(p,ps,1)

    tvfull = tfull
    tvfull = tfull*(1 + ep1*qvfull - qcfull)

    topo3d = conform_dims(dimsizes(ps),topo,(/1,2/))

    zgeotmp = hydro( pfull(time|:,iy|:,jx|:,KZ|::-1),  \
                  tvfull(time|:,iy|:,jx|:,KZ|::-1),  \
                  topo3d)

    copy_VarCoords(tvfull(time|:,iy|:,jx|:,KZ|:),zgeotmp)
    copy_VarAtts(tvfull(time|:,iy|:,jx|:,KZ|:),zgeotmp)
    zgeo = zgeotmp(time|:,KZ|1:nlev,iy|:,jx|:)
    zgeo!1 = "kz"

    delete(tfull)
    delete(qvfull)
    delete(qcfull)
    delete(pfull)
    delete(tvfull)
    delete(topo3d)
    delete(zgeotmp)


    zgeo@long_name = "Geopotential Height"
    zgeo@standard_name = "geopotential_height"
    zgeo@units = "m"

    fhopt->hgt = zgeo(:,::-1,:,:)
    if(bWriteOut)then
      print("...Finished Writing Height")
    end if
  end if

  if(calcvar@Theta)then
    if(bWriteOut)then
      print("Calculating Theta...")
    end if
      t = fhatm->t
      theta = t
      theta = t*(p0exner/p)^rovcp

      theta@long_name = "Potential Temperature" 
      theta@long_name = "potential_temperature"
      fhopt->theta = theta

    if(bWriteOut)then
      print("...Finished Writing Theta")
    end if
  end if

  if(calcvar@Fog)then
    if(bWriteOut)then
      print("Calculating Fog...")
    end if
;    cld = fhrad->cld
    qc = fhatm->qc
    relgeo = zgeo
    relgeo = zgeo - conform_dims(dimsizes(zgeo),topo,(/2,3/))
    ;Use the calculateFog() function to calculate fog
    resfog = True
    resfog@UseLWC = True
    resfog@LWCthreshold = 5e-5
    ;fog = calculateFog(cld(time|:,iy|:,jx|:,kz|:),relgeo(time|:,iy|:,jx|:,kz|::-1),resfog) 
    fog = calculateFog(qc(time|:,iy|:,jx|:,kz|:),relgeo(time|:,iy|:,jx|:,kz|::-1),resfog) 
    fhopt->fog = fog
    if(bWriteOut)then
      print("...Finished Writing Fog")
    end if
  end if

  
end
