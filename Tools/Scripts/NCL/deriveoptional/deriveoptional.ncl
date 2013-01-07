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
  calcvar@InCldLWP = False
  calcvar@Subsidence = True
  calcvar@ZDens = True
  calcvar@rho2m = True
  calcvar@occzda2000 = True
  calcvar@CurlTau = True
  calcvar@SnowFreeAlbedoSW = True

;;  ;Temporary turn-off
;  calcvar@Pressure = False
;  calcvar@GeoHeight = False
;  calcvar@Density = False
;  calcvar@Fog = False
;  calcvar@Theta = False
;  calcvar@InCldLWP = False
;  calcvar@ZDens = False
;  calcvar@occzda2000 = False


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

  gammastd = 6.5e-3
  t0std = 288.15
  p0std = 101325
  rho0std = p0std/(rgas*t0std);
  

  p0exner = 1000.0
  rovcp = rgas/cpd

  areathresh = 0.5 ;The threshold for aeral coverage of cloud (in km^2)

  fhatm = addfile(atmfilename,"r")
  fhsrf = addfile(srffilename,"r")
  fhrad = addfile(radfilename,"r")
  fhopt = addfile(optfilename,"w")


  ;Run through dependencies
  ;TODO: Put this in a function
;Level 3: Functions that depend on a level 2 variable at highest
;Level 2: Functions that depend on a level 1 variable at highest
  if(calcvar@Fog)then
    calcvar@Geoheight = True
  end if
  if(calcvar@Subsidence)then
    calcvar@Density = True
  end if
  if(calcvar@occzda2000)then
    calcvar@ZDens = True
  end if
;Level 1: Functions that only depend on a level 0 variable
  if(calcvar@Theta)then
    calcvar@Pressure = True
  end if
  if(calcvar@Density)then
    calcvar@Pressure = True
  end if
  if(calcvar@GeoHeight)then
    calcvar@Pressure = True
  end if
  if(calcvar@ZDens)then
    calcvar@rho2m = True
  end if
  if(calcvar@CurlTau)then
    calcvar@rho2m = True
  end if
;Level 0: Variables that depend on nothing:
  ; Pressure
  ; Density Altitude


  if(calcvar@Pressure)then
    if(bWriteOut)then
      print("Calculating Pressure...")
    end if
    ps    = fhatm->ps
    sigma = fhatm->sigma
    ptop  = fhatm->ptop
    time = fhatm->time

    t = fhatm->ta

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
    t = fhatm->ta
    qv = fhatm->qas
    qc = fhatm->clw

    qv = 1/(1+qv)
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
    t = fhatm->ta
    qv = fhatm->qas
    qc = fhatm->clw
    ts = fhsrf->tas(:,0,:,:)
    qvs = fhsrf->qas(:,0,:,:)
    qv = 1 / (1+qv)
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
      t = fhatm->ta
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
;    cld = fhrad->cl
    qc = fhatm->clw
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

  if(calcvar@InCldLWP)then
    if(bWriteOut)then
      print("Calculating In-Cloud Liquid Water Path...")
    end if
    ;Read in the grid-cell cloud fraction
    cld = fhrad->cl
    ;Read in the grid-cell average liquid water path
    lwp = fhrad->clw

    totcld = 1 - dim_product(1 - cld(time|:,iy|:,jx|:,kz|:))
    totlwp = dim_sum_Wrap(lwp(time|:,iy|:,jx|:,kz|:))

    ds = fhrad@grid_size_in_meters/1000 ; Get the grid resolution (in meters)
    cloudarea = totcld*ds^2 ;Calculate the approximate areal coverage of cloud
    

    incldlwp = totlwp
    incldlwp@_FillValue = 1e-30
    ;Get the in-cloud liquid water path for clouds whose area is greater
    ; than the threshold (approximates having a lower limit to a satellite's 
    ; resolution)
    incldlwp = totlwp/where(cloudarea.lt.areathresh,incldlwp@_FillValue,totcld)
    delete(cloudarea)
    delete(ds)

    incldlwp@long_name = "Integrated In-Cloud Liquid Water Path"
    incldlwp@standard_name = ""
    
    fhopt->incldlwp = incldlwp
    if(bWriteOut)then
      print("...Finished Writing In-Cloud Liquid Water Path")
    end if
  end if

  if(calcvar@Subsidence)then
    if(bWriteOut)then
      print("Calculating Subsidence Rate...")
    end if
      omega = fhatm->omega
      
      w = omega
      w = -omega/(rho*gti)
      w@long_name = "Vertical Velocity"
      w@standard_name = "upward_air_velocity"
      w@units = "m s-1"

      fhopt->w = w

    if(bWriteOut)then
      print("...Finished writing Subsidence Rate")
    end if
  end if

  if(calcvar@rho2m)then
    if(bWriteOut)then
      print("Calculating 2m Density...")
    end if

      alpha = gti/(rgas*gammastd)
      t2m = fhsrf->tas(:,0,:,:)
      q2m = fhsrf->qas(:,0,:,:)
      ps = fhsrf->ps

      tv2m = t2m*(1 + ep1*q2m)

      rho2m = t2m
      rho2m = 100.0*ps/(rgas*tv2m)
      rho2m@long_name = "2m Air Density"
      rho2m@standard_name = "air_density_at_2m"
      rho2m@units = "kg m-3"

      fhopt->rho2m = rho2m

    if(bWriteOut)then
      print("...Finished writing 2m Density") 
    end if
  end if


  if(calcvar@ZDens)then
    if(bWriteOut)then
      print("Calculating Density Altitude...")
    end if

      zdens = rho2m
      zdens = t0std/gammastd*(1 - (rho2m/rho0std)^(1/(alpha-1)));

      zdens@long_name = "Density Altitude"
      zdens@standard_name = "density_equivalent_altitude_in_standard_atmosphere"
      zdens@units = "m"

      fhopt->zdens = zdens

    if(bWriteOut)then
      print("...Finished writing Density Altitude") 
    end if
  end if

  if(calcvar@occzda2000)then
    if(bWriteOut)then
      print("Calculating Density Altitude Anomaly > 2000 m...")
    end if

      topo = fhatm->topo
      topo3d = conform_dims(dimsizes(zdens),topo,(/1,2/))

      occzda2000 = zdens
      occzda2000 = where((zdens-topo3d).ge.2000,1.0,0.0)
      occzda2000@long_name = "Dens. Alt. Anomaly > 2000m Occurrence"
      occzda2000@standard_name = "density_altitude_anomaly_occurrence"
      occzda2000@units = "fraction"

      fhopt->occzda2000 = occzda2000 

    if(bWriteOut)then
      print("...Finished writing Density Altitude Anomaly > 2000 m...")
    end if
  end if

  if(calcvar@CurlTau)then
    if(bWriteOut)then
      print("Calculating Wind Stress Curl...")
    end if

      uvdrag = fhsrf->drag
      u10m = fhsrf->uas(:,0,:,:)
      v10m = fhsrf->vas(:,0,:,:)
      ds = fhsrf@grid_size_in_meters; Get the grid resolution (in meters)

      taux = u10m
      tauy = v10m
      taux = rho2m*uvdrag*u10m*abs(u10m)
      tauy = rho2m*uvdrag*v10m*abs(v10m)

      ;TODO: Add in map scale factor
      dtauxdy = center_finite_diff(taux(time|:,jx|:,iy|:),ds,False,0)
      dtauydx = center_finite_diff(tauy(time|:,iy|:,jx|:),ds,False,0)

      dtauxdy!0 = "time"
      dtauxdy!1 = "jx"
      dtauxdy!2 = "iy"
      dtauydx!0 = "time"
      dtauydx!2 = "jx"
      dtauydx!1 = "iy"

      dtauxdy&time = u10m&time
      dtauxdy&jx = u10m&jx
      dtauxdy&iy = u10m&iy
      dtauydx&time = v10m&time
      dtauydx&jx = v10m&jx
      dtauydx&iy = v10m&iy

      curltau = dtauydx
      curltau = dtauydx - dtauxdy(time|:,iy|:,jx|:)
      curltau@long_name = "Wind Stress Curl"
      curltau@units = "Pa m-1"

      fhopt->curltau = curltau 

    if(bWriteOut)then
      print("...Finished writing Wind Stress Curl") 
    end if
  end if

  if(calcvar@SnowFreeAlbedoSW)then
    if(bWriteOut)then
      print("Calculating Snow Free Land Albedo (SW)...")
    end if

      aldirs = fhsrf->aldirs
      scv = fhsrf->snv

      aldirs@_FillValue = -1e-20
      sfalbedosw = aldirs

      sfalbedosw = where(scv.eq.0.0,aldirs,aldirs@_FillValue)
      sfalbedosw@long_name = "Surface Albedo (Snow Free)"
      sfalbedosw@standard_name = ""
      sfalbedosw@units = "fraction"

      fhopt->sfaldirs = sfalbedosw 

    if(bWriteOut)then
      print("...Finished writing Snow Free Land Albedo (SW)...")
    end if
  end if




end
