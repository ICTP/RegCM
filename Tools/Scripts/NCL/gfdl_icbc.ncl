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
;This script uses NCL to generate ICBC's for RegCM4.1 from GFDL model output
;Created 04/2011 by Travis A. O'Brien for RegCM v4.1

;load assistive scripts
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$REGCMSRCDIR/Tools/Scripts/NCL/rcmdatefunctions.ncl"

;*******************************************************************************
;*******************************************************************************
;************************** Start of GFDL Script *******************************
;*******************************************************************************
;*******************************************************************************

begin

  pi = 3.1415926
  gammastd = 6.5e-3
  t0std = 288.15
  g = 9.81
  rgas = 287.
  cpd = 1004.
  alpha = g/(rgas*gammastd)

;*******************************************************************************
;*******************************************************************************
;************************** Initialize ICBC program ****************************
;*******************************************************************************
;*******************************************************************************
;These variables are set by the wrapper script gfdl_ncl
;DOMNAME = "GFDL.hist"
;istartdate = 1990010100
;ienddate = 1990010200
;bcdir = "/storage/Data/GFDL/modern"
;inpdir = "./input"

domfile = inpdir + "/" + DOMNAME + "_DOMAIN000.nc"

dtbc = 6 ;The GCM time step
sdate = parse_date(istartdate)
edate = parse_date(ienddate)

edatep1 = increment_date(edate,6,False)

bMethod1 = True

;*******************************************************************************
;*******************************************************************************
;************************** Read RCM Domain File *******************************
;*******************************************************************************
;*******************************************************************************

  frcmdom = addfile(domfile,"r")
  topo = frcmdom->topo
  xmap = frcmdom->xmap
  dmap = frcmdom->dmap
  xlat = frcmdom->xlat
  xlon = frcmdom->xlon
  dlat = frcmdom->dlat
  dlon = frcmdom->dlon
  landuse = frcmdom->landuse
  sigma = frcmdom->sigma
  ptop = frcmdom->ptop
  clon = dble2flt(frcmdom@longitude_of_projection_origin)
  xvec = frcmdom->jx
  yvec = frcmdom->iy

;*******************************************************************************
;*******************************************************************************
;*************************** Pre-set File Parameters****************************
;*******************************************************************************
;*******************************************************************************

  kz = dimsizes(sigma)-1
  dsize = dimsizes(xlat)
  iy = dsize(0)
  jx = dsize(1)
  delete(dsize)

  sigmahalf = sigma(1:)
  sigmahalf = (sigma(0:kz-1) + sigma(1:kz))/2

  ;Determine the sinx and cosx parameters for rotating u and v
  if(clon.ge.0)then
    clonshift = clon
  else
    clonshift = clon +360.
  end if
  dlonshift = where(dlon.ge.0.0,dlon,dlon+360)
  x = (clonshift - dlonshift)*(pi/180.)*dmap
  sinx = sin(x)
  cosx = cos(x)


  fAtt = True
  fAtt@coordinates = "xlat xlon"

  dimNames = (/"time","iy","jx","kz"/)
  dimSizes = (/-1,iy,jx,kz/)
  dimUnlim = (/True,False,False,False/)

  varNames4D = (/"t","u","v","qv"/)
  varTypes4D = (/"float","float","float","float"/)
  dims4D = (/"time","kz","iy","jx"/)

  varNames3D = (/"ts","ps"/)
  varTypes3D = (/"float","float"/)
  dims3D = (/"time","iy","jx"/)

  varNames2D = (/"xlat","xlon","topo"/)
  varTypes2D = (/"float","float","float"/)
  dims2D = (/"iy","jx"/)

  ptop = ptop*100 ;Convert to Pa

  timeunits = "hours since " + date2str(sdate)

;*******************************************************************************
;*******************************************************************************
;*************************** Main loop over dates ******************************
;*******************************************************************************
;*******************************************************************************

bProcessingDates = True
curdate = sdate
itimegcm = 0
itimeicbc = 0
do while(bProcessingDates)

  ;*****************************************************************************
  ;************************** Initialize I/O ***********************************
  ;*****************************************************************************

  ;Check if a new ICBC file needs to be opened
  startmonth = curdate
  startmonth@day = 1
  startmonth@hour = 0 
  if(equiv_date(curdate,startmonth))then
    itimeicbc = 0
  ;Open a new output file
    ficbcname = inpdir+"/"+DOMNAME+"_ICBC."+idatestr(startmonth)+".nc"
    print("filename: " + ficbcname)
    system("rm " + ficbcname + " >& /dev/null")
    ficbc = addfile(ficbcname,"c")
    ;put the file in define mode
    setfileoption(ficbc,"DefineMode",True)
    ;Define the dimensions of the file
    filedimdef(ficbc,dimNames,dimSizes,dimUnlim)
    ;Define the variables of the file
    filevardef(ficbc,varNames4D,varTypes4D,dims4D)
    filevardef(ficbc,varNames3D,varTypes3D,dims3D)
    filevardef(ficbc,varNames2D,varTypes2D,dims2D)
    filevardef(ficbc,"sigma","float","kz")
    filevardef(ficbc,"iy","float","iy")
    filevardef(ficbc,"jx","float","jx")
    filevardef(ficbc,"time","double","time")

    dumAtt = fAtt
    dumAtt@standard_name = "surface_air_pressure" ;
    dumAtt@long_name = "Surface pressure" ;
    dumAtt@units = "hPa" ;
    filevarattdef(ficbc,"ps",dumAtt)
    dumAtt = fAtt
    dumAtt@standard_name = "surface_temperature" ;
    dumAtt@long_name = "Surface Temperature" ;
    dumAtt@units = "K" ;
    filevarattdef(ficbc,"ts",dumAtt)
    dumAtt = fAtt
    dumAtt@standard_name = "eastward_wind" ;
    dumAtt@long_name = "U component (westerly) of wind" ;
    dumAtt@units = "m s-1" ;
    filevarattdef(ficbc,"u",dumAtt)
    dumAtt = fAtt
    dumAtt@standard_name = "northward_wind" ;
    dumAtt@long_name = "V component (southerly) of wind" ;
    dumAtt@units = "m s-1" ;
    filevarattdef(ficbc,"v",dumAtt)
    dumAtt = fAtt
    dumAtt@standard_name = "air_temperature" ;
    dumAtt@long_name = "Temperature" ;
    dumAtt@units = "K" ;
    filevarattdef(ficbc,"t",dumAtt)
    dumAtt = fAtt
    dumAtt@standard_name = "humidity_mixing_ratio" ;
    dumAtt@long_name = "Water vapor mixing ratio" ;
    dumAtt@units = "kg kg-1" ;
    filevarattdef(ficbc,"qv",dumAtt)

    dumAtt = True
    dumAtt@units = timeunits
    filevarattdef(ficbc,"time",dumAtt)

    ;Grab the attributes from the domain file
    copy_VarAtts(frcmdom,ficbc)
    ;Add a custom comment
    ficbc@method = "Created from gfdl_icbc.ncl (written by TA O'Brien)"

    ;Take the file out of define mode
    setfileoption(ficbc,"DefineMode",False)
    ;Write out static variables
    ficbc->sigma = sigmahalf
    ficbc->jx = xvec
    ficbc->iy = yvec
    ficbc->xlat = xlat
    ficbc->xlon = xlon
    ficbc->topo = topo
  end if

  startyear = startmonth
  startyear@month = 1
  if(equiv_date(curdate,startmonth))then
  ;Open the new input files
    itimegcm = 0
    endyear = curdate
    endyear@month = 12
    endyear@day = 31
    endyear@hour = 23
    fmidstr = idatestr(startyear) + "-" + idatestr(endyear)

    ;Set all the variable file names
    uname = bcdir+"/ua_A6."+fmidstr+".nc"
    vname = bcdir+"/va_A6."+fmidstr+".nc"
    tname = bcdir+"/ta_A6."+fmidstr+".nc"
    qname = bcdir+"/hus_A6."+fmidstr+".nc"
    psname = bcdir+"/ps_A6."+fmidstr+".nc"
    sstname = bcdir+"/sst_O6."+fmidstr+"_interp.nc"
    oroname = bcdir + "/../orog_A1_subset.static.nc"
    
    ;Open all the variable files
    fu = addfile(uname,"r")
    fv = addfile(vname,"r")
    ft = addfile(tname,"r")
    fq = addfile(qname,"r")
    fps = addfile(psname,"r")
    fsst = addfile(sstname,"r")
    foro = addfile(oroname,"r")

    ;Create 3d sinx and cosx arrays for rotating winds
    ;also pre-create temporary variables
    kzgcm = dimsizes(fu->ua(0,:,0,0))
    iygcm = dimsizes(fu->ua(0,0,:,0))
    jxgcm = dimsizes(fu->ua(0,0,0,:))
    ds3d = (/kzgcm,iy,jx/)
    sinx3d = conform_dims(ds3d,sin(x),(/1,2/))
    cosx3d = conform_dims(ds3d,cos(x),(/1,2/))
    utmp = sinx3d
    vtmp = sinx3d

    ;Create a super-variable for interpolating
    supervar = new((/2*kzgcm+1,iygcm,jxgcm/),float)

  end if

  ;Print user output
  print("Processing: " + date2str(curdate))

  ;*****************************************************************************
  ;************************* Read GCM T,U,V,Q,PS *******************************
  ;*****************************************************************************
  t_gcm = ft->ta(itimegcm,:,:,:)
  qv_gcm = fq->hus(itimegcm,:,:,:)
  u_gcm = fu->ua(itimegcm,:,:,:)
  v_gcm = fv->va(itimegcm,:,:,:)
  ps_gcm = fps->ps(itimegcm,:,:)
  sst_gcm = fsst->sst(itimegcm,:,:)
  oro_gcm = foro->orog(:,:)

  ;Shift the longitudes (rgrid2rcm won't work with 0:360 longitudes; it needs 
  ; -180:180 longitudes)
  t_gcm&lon = t_gcm&lon - 360
  qv_gcm&lon = qv_gcm&lon - 360
  u_gcm&lon = u_gcm&lon - 360
  v_gcm&lon = v_gcm&lon - 360
  ps_gcm&lon = ps_gcm&lon - 360
  sst_gcm&lon = sst_gcm&lon - 360
  oro_gcm&lon = oro_gcm&lon - 360


  ;*****************************************************************************
  ;************************ Convert surface pressure to SLP ********************
  ;*****************************************************************************

  slp_gcm = ps_gcm
  slp_gcm = ps_gcm/((t0std - gammastd*oro_gcm)/t0std)^alpha

  ;Create a set of GCM pressures
  p_gcm = conform_dims(dimsizes(t_gcm),dble2flt(t_gcm&plev),(/0/))

  theta_gcm = t_gcm
  theta_gcm = t_gcm*(100000/p_gcm)^(rgas/cpd)

  ;*****************************************************************************
  ;************************ Vertically interpolate to P-Levs *******************
  ;*****************************************************************************
  ;TAO: GFDL data are already on pressure levels

  ;*****************************************************************************
  ;********************** Horizontally interpolate to RCM grid *****************
  ;*****************************************************************************
  if(bMethod1)then
    theta_hint = rgrid2rcm(theta_gcm&lat,theta_gcm&lon,theta_gcm,xlat,xlon,1)
    qv_hint = rgrid2rcm(qv_gcm&lat,qv_gcm&lon,qv_gcm,xlat,xlon,1)
    ;Interpolate the SSTs
    sst_rcm = rgrid2rcm(sst_gcm&lat,sst_gcm&lon,sst_gcm,xlat,xlon,1)
    ;Interpolate slp to the cross and dot grids
    slpx_hint = rgrid2rcm(slp_gcm&lat,slp_gcm&lon,slp_gcm,xlat,xlon,1)
    slpd_hint = rgrid2rcm(slp_gcm&lat,slp_gcm&lon,slp_gcm,dlat,dlon,1)
    u_hint = rgrid2rcm(u_gcm&lat,u_gcm&lon,u_gcm,dlat,dlon,1)
    v_hint = rgrid2rcm(v_gcm&lat,v_gcm&lon,v_gcm,dlat,dlon,1)
  else
    ;Interpolate cross-grid vars
    supervar(0:kzgcm-1,:,:) = theta_gcm
    supervar(kzgcm:(2*kzgcm-1),:,:) = qv_gcm
    supervar(2*kzgcm,:,:) = slp_gcm
    sghint = rgrid2rcm(theta_gcm&lat,theta_gcm&lon,supervar,xlat,xlon,1)
    theta_hint = sghint(0:kzgcm-1,:,:)
    qv_hint = sghint(kzgcm:(2*kzgcm-1),:,:)
    slpx_hint = sghint(2*kzgcm,:,:)
    ;Interpolate dot-grid vars
    supervar(0:kzgcm-1,:,:) = u_gcm
    supervar(kzgcm:(2*kzgcm-1),:,:) = v_gcm
    sghint = rgrid2rcm(theta_gcm&lat,theta_gcm&lon,supervar,dlat,dlon,1)
    u_hint = sghint(0:kzgcm-1,:,:)
    v_hint = sghint(kzgcm:(2*kzgcm-1),:,:)
    slpd_hint = sghint(2*kzgcm,:,:)
    sst_rcm = rgrid2rcm(sst_gcm&lat,sst_gcm&lon,sst_gcm,xlat,xlon,1)
  end if

  ;*****************************************************************************
  ;********************** Calculate surface P, T *******************************
  ;*****************************************************************************
  psx_rcm = slpx_hint*((t0std - gammastd*topo)/t0std)^alpha
  psd_rcm = slpx_hint*((t0std - gammastd*topo)/t0std)^alpha

  ;*****************************************************************************
  ;********************** Rotate wind vectors **********************************
  ;*****************************************************************************
  utmp = u_hint*cosx3d + v_hint*sinx3d
  vtmp = -u_hint*sinx3d + v_hint*cosx3d
  u_hint = utmp
  v_hint = vtmp

  ;*****************************************************************************
  ;******************* Vertically interpolate to RCM grid **********************
  ;*****************************************************************************
  ;Calculate the pressure levels of the rcm grid
  dsrcm = dimsizes(theta_hint)
  dsrcm(0) = dimsizes(sigmahalf)
  psx3d = conform_dims(dsrcm,psx_rcm,(/1,2/))
  psd3d = conform_dims(dsrcm,psd_rcm,(/1,2/))
  sig3d = conform_dims(dsrcm,sigmahalf,(/0/))

  px_rcm = sig3d*(psx3d-ptop) + ptop
  pd_rcm = sig3d*(psd3d-ptop) + ptop
  
  theta_rcm = int2p_n(dble2flt(theta_gcm&plev),theta_hint,px_rcm,-1,0)
  qv_rcm = int2p_n(dble2flt(qv_gcm&plev),qv_hint,px_rcm,-2,0) ;Do log interpolation
  u_rcm = int2p_n(dble2flt(u_gcm&plev),u_hint,pd_rcm,-1,0)
  v_rcm = int2p_n(dble2flt(v_gcm&plev),v_hint,pd_rcm,-1,0)

  t_rcm = theta_rcm
  t_rcm = theta_rcm*(px_rcm/100000)^(rgas/cpd)
  ts_rcm = theta_rcm(kz-1,:,:)
  ts_rcm = theta_rcm(kz-1,:,:)*(psx_rcm/100000)^(rgas/cpd)

  ;*****************************************************************************
  ;*********************** Set SSTs ********************************************
  ;*****************************************************************************
;  ts_rcm = where(landuse.eq.14.or.landuse.eq.15,sst_rcm,ts_rcm)
  ts_rcm = where(landuse.eq.15,sst_rcm,ts_rcm)

  ;*****************************************************************************
  ;*********************** Pre-output conversion *******************************
  ;*****************************************************************************
  psx_rcm = psx_rcm/100
  iyear = round(curdate@year,3) 
  imonth = round(curdate@month,3) 
  iday = round(curdate@day,3) 
  ihour = round(curdate@hour,3) 
  curtime = ut_inv_calendar(iyear,imonth,iday,ihour,0,0,timeunits,0)

  ;*****************************************************************************
  ;*************************** Write output ************************************
  ;*****************************************************************************
  ficbc->time(itimeicbc) = (/curtime/)
  ficbc->t(itimeicbc,:,:,:) = (/t_rcm/)
  ficbc->u(itimeicbc,:,:,:) = (/u_rcm/)
  ficbc->v(itimeicbc,:,:,:) = (/v_rcm/)
  ficbc->qv(itimeicbc,:,:,:) = (/qv_rcm/)
  ficbc->ps(itimeicbc,:,:) = (/psx_rcm/)
  ficbc->ts(itimeicbc,:,:) = (/ts_rcm/)

  ;Increment the date
  prevdate = curdate
  curdate = increment_date(prevdate,6,False)
  itimegcm = itimegcm + 1 ;Increment the GCM data time counter
  itimeicbc = itimeicbc + 1 ;Increment the ICBC data time counter
  ;Check if we reached the end
  if(equiv_date(curdate,edatep1))then
    bProcessingDates = False
  end if
end do


end





