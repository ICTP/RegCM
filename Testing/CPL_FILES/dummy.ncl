;-----------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;-----------------------------------------------------------
;***********************************************************
;*** This script is used to create dummy forcing file ******
;*** for the coupled RegCM-ROMS model **********************
;***********************************************************
;***********************************************************
;*** Ufuk Turuncoglu ***************************************
;*** Send bug reports to turuncu@be.itu.edu.tr *************
;***********************************************************
begin
  ;--- parameters ---
  dstr = "1995-01-01"
  clen = 2 ;365
  cval = 1.0d20
  ;time = fspan(0,364.75,365*4)
  time = fspan(0,clen-0.25,tointeger(clen)*4)

  ;--- open grid file --- 
  grid = addfile("initial-v1/CAS_grd_v1.nc", "r")

  ;--- read variables ---
  x_rho = grid->x_rho
  y_rho = grid->y_rho
  x_u = grid->x_u
  y_u = grid->y_u
  x_v = grid->x_v
  y_v = grid->y_v

  ;--- get dimension sizes ---
  dims = dimsizes(x_rho)
  eta_rho = dims(0)
  xi_rho = dims(1)
  dims = dimsizes(x_u)
  eta_u = dims(0)
  xi_u = dims(1)
  dims = dimsizes(x_v)
  eta_v = dims(0)
  xi_v = dims(1)
  ocean_time = dimsizes(time)

  ;--- get variables ---
  x_rho = grid->x_rho
  y_rho = grid->y_rho
  x_u = grid->x_u
  y_u = grid->y_u
  x_v = grid->x_v
  y_v = grid->y_v

  ;--- create forcing file ---
  fname = "CAS_forcing_dummy.nc"
  system("/bin/rm -f "+fname)
  fout = addfile(fname, "c")

  ;--- define mode on ---
  setfileoption(fout, "DefineMode", True)

  ;--- add global attributes ---
  attr = True
  attr@title = "RegCM-ROMS dummy forcing file"
  attr@source_file = "roms_dummy.ncl"
  attr@creation_date = systemfunc ("date")
  fileattdef(fout, attr)
  delete(attr)

  ;--- define dimensions ---
  dimNames = (/ "xi_rho"   , "xi_u"     , "xi_v"     , "eta_rho"  , \
                "eta_v"    , "eta_u"    , "sms_time" , "shf_time" , \
                "srf_time" , "swf_time" , "wind_time", "pair_time", \
                "tair_time", "qair_time", "rain_time", "lrf_time" /)
  dimSizes = (/ xi_rho, xi_u, xi_v, eta_rho, eta_v, eta_u, ocean_time, ocean_time, \
                ocean_time, ocean_time, ocean_time, ocean_time, ocean_time, \
                ocean_time, ocean_time, ocean_time /)
  dimUnlim = (/ False, False, False, False, False, False, False, False, \
                False, False, False, False, False, False, False, False /)
  filedimdef(fout, dimNames, dimSizes, dimUnlim)

  filevardef(fout, "sustr" , typeof(x_rho), (/ "sms_time", "eta_u"  , "xi_u"   /))
  filevardef(fout, "svstr" , typeof(y_rho), (/ "sms_time", "eta_v"  , "xi_v"   /))
  filevardef(fout, "shflux" , typeof(y_rho), (/ "shf_time", "eta_rho", "xi_rho" /))
  filevardef(fout, "swrad" , typeof(y_rho), (/ "srf_time", "eta_rho", "xi_rho" /))
  filevardef(fout, "swflux" , typeof(y_rho), (/ "swf_time", "eta_rho", "xi_rho" /))
  filevardef(fout, "Pair" , typeof(y_rho), (/ "pair_time", "eta_rho", "xi_rho" /))
  filevardef(fout, "Tair" , typeof(y_rho), (/ "tair_time", "eta_rho", "xi_rho" /))
  filevardef(fout, "Qair" , typeof(y_rho), (/ "qair_time", "eta_rho", "xi_rho" /))
  filevardef(fout, "lwrad" , typeof(y_rho), (/ "lrf_time", "eta_rho", "xi_rho" /))
  filevardef(fout, "lwrad_down" , typeof(y_rho), (/ "lrf_time", "eta_rho", "xi_rho" /))
  filevardef(fout, "Uwind" , typeof(y_rho), (/ "wind_time", "eta_rho", "xi_rho" /))
  filevardef(fout, "Vwind" , typeof(y_rho), (/ "wind_time", "eta_rho", "xi_rho" /))
  filevardef(fout, "rain" , typeof(y_rho), (/ "rain_time", "eta_rho", "xi_rho" /))
  filevardef(fout, "sms_time" , typeof(x_rho), "sms_time" )
  filevardef(fout, "shf_time" , typeof(x_rho), "shf_time" )
  filevardef(fout, "srf_time" , typeof(x_rho), "srf_time" )
  filevardef(fout, "swf_time" , typeof(x_rho), "swf_time" )
  filevardef(fout, "wind_time", typeof(x_rho), "wind_time")
  filevardef(fout, "pair_time", typeof(x_rho), "pair_time")
  filevardef(fout, "tair_time", typeof(x_rho), "tair_time")
  filevardef(fout, "qair_time", typeof(x_rho), "qair_time")
  filevardef(fout, "rain_time", typeof(x_rho), "rain_time")
  filevardef(fout, "lrf_time" , typeof(x_rho), "lrf_time" )

  ;--- add attributes ---
  ;--- time dimensions ---
  attr = True
  attr@long_name = "surface momentum stress time"
  attr@units = "days since "+dstr+" 00:00:00" 
  attr@calendar = "standard"
  if (clen .gt. 0) then
  attr@cycle_length = clen
  end if
  attr@field = "sms_time, scalar, series" ;
  filevarattdef(fout, "sms_time", attr)
  delete(attr)

  attr = True
  attr@long_name = "surface heat flux time"
  attr@units = "days since "+dstr+" 00:00:00"
  attr@calendar = "standard"
  if (clen .gt. 0) then
  attr@cycle_length = clen
  end if
  attr@field = "shf_time, scalar, series" ;
  filevarattdef(fout, "shf_time", attr)
  delete(attr)

  attr = True
  attr@long_name = "solar shortwave radiation time"
  attr@units = "days since "+dstr+" 00:00:00"
  attr@calendar = "standard"
  if (clen .gt. 0) then
  attr@cycle_length = clen
  end if
  attr@field = "srf_time, scalar, series" ;
  filevarattdef(fout, "srf_time", attr)
  delete(attr)

  attr = True
  attr@long_name = "surface freshwater flux time"
  attr@units = "days since "+dstr+" 00:00:00"
  attr@calendar = "standard"
  if (clen .gt. 0) then
  attr@cycle_length = clen
  end if
  attr@field = "swf_time, scalar, series" ;
  filevarattdef(fout, "swf_time", attr)
  delete(attr)

  attr = True
  attr@long_name = "surface wind time"
  attr@units = "days since "+dstr+" 00:00:00"
  attr@calendar = "standard"
  if (clen .gt. 0) then
  attr@cycle_length = clen
  end if
  attr@field = "wind_time, scalar, series" ;
  filevarattdef(fout, "wind_time", attr)
  delete(attr)

  attr = True
  attr@long_name = "surface pressure time"
  attr@units = "days since "+dstr+" 00:00:00"
  attr@calendar = "standard"
  if (clen .gt. 0) then
  attr@cycle_length = clen
  end if
  attr@field = "pair_time, scalar, series" ;
  filevarattdef(fout, "pair_time", attr)
  delete(attr)

  attr = True
  attr@long_name = "surface temperature time"
  attr@units = "days since "+dstr+" 00:00:00"
  attr@calendar = "standard"
  if (clen .gt. 0) then
  attr@cycle_length = clen
  end if
  attr@field = "tair_time, scalar, series" ;
  filevarattdef(fout, "tair_time", attr)
  delete(attr)

  attr = True
  attr@long_name = "surface humidity time"
  attr@units = "days since "+dstr+" 00:00:00"
  attr@calendar = "standard"
  if (clen .gt. 0) then
  attr@cycle_length = clen
  end if
  attr@field = "qair_time, scalar, series" ;
  filevarattdef(fout, "qair_time", attr)
  delete(attr)

  attr = True
  attr@long_name = "precipitation time"
  attr@units = "days since "+dstr+" 00:00:00"
  attr@calendar = "standard"
  if (clen .gt. 0) then
  attr@cycle_length = clen
  end if
  attr@field = "rain_time, scalar, series" ;
  filevarattdef(fout, "rain_time", attr)
  delete(attr)

  attr = True
  attr@long_name = "longwave radiation time"
  attr@units = "days since "+dstr+" 00:00:00"
  attr@calendar = "standard"
  if (clen .gt. 0) then
  attr@cycle_length = clen
  end if
  attr@field = "lrf_time, scalar, series" ;
  filevarattdef(fout, "lrf_time", attr)
  delete(attr)

  ;--- variables ---
  attr = True
  attr@long_name = "surface u-momentum stress"
  attr@units = "Newton meter-2"
  attr@time = "sms_time"
  filevarattdef(fout, "sustr", attr)
  delete(attr)

  attr = True
  attr@long_name = "surface v-momentum stress"
  attr@units = "Newton meter-2" 
  attr@time = "sms_time"
  filevarattdef(fout, "svstr", attr)
  delete(attr)

  attr = True
  attr@long_name = "surface net heat flux"
  attr@units = "Watts meter-2"
  attr@time = "shf_time"
  filevarattdef(fout, "shflux", attr)
  delete(attr)

  attr = True
  attr@long_name = "solar shortwave radiation"
  attr@units = "Watts meter-2"
  attr@time = "srf_time"
  filevarattdef(fout, "swrad", attr)
  delete(attr)

  attr = True
  attr@long_name = "surface freshwater flux (E-P)"
  attr@units = "centimeter day-1"
  attr@time = "swf_time"
  filevarattdef(fout, "swflux", attr)
  delete(attr)

  attr = True
  attr@long_name = "surface air pressure"
  attr@units = "pascal"
  attr@time = "pair_time"
  filevarattdef(fout, "Pair", attr)
  delete(attr)

  attr = True
  attr@long_name = "surface air temperature"
  attr@units = "Celsius"
  attr@time = "tair_time"
  filevarattdef(fout, "Tair", attr)
  delete(attr)

  attr = True
  attr@long_name = "surface air specific humidity"
  attr@units = "kg/kg"
  attr@time = "qair_time"
  filevarattdef(fout, "Qair", attr)
  delete(attr)

  attr = True
  attr@long_name = "surface u-wind component"
  attr@units = "meter second-1"
  attr@time = "wind_time"
  filevarattdef(fout, "Uwind", attr)
  delete(attr)

  attr = True
  attr@long_name = "surface v-wind component"
  attr@units = "meter second-1"
  attr@time = "wind_time"
  filevarattdef(fout, "Vwind", attr)
  delete(attr)

  attr = True
  attr@long_name = "rain fall rate"
  attr@units = "kilogram meter-2 second-1"
  attr@time = "rain_time"
  filevarattdef(fout, "rain", attr)
  delete(attr)

  attr = True
  attr@long_name = "net longwave radiation flux"
  attr@units = "Watt meter-2"
  attr@time = "lrf_time"
  filevarattdef(fout, "lwrad", attr)
  delete(attr)

  attr = True
  attr@long_name = "downward longwave radiation flux"
  attr@units = "Watt meter-2"
  attr@time = "lrf_time"
  filevarattdef(fout, "lwrad_down", attr)
  delete(attr)

  ;--- exit file definition mode ---
  setfileoption(fout, "DefineMode", False)

  ;--- fill data ---
  fout->sms_time  = (/ time /)
  fout->shf_time  = (/ time /)
  fout->srf_time  = (/ time /)
  fout->swf_time  = (/ time /)
  fout->wind_time = (/ time /)
  fout->pair_time = (/ time /)
  fout->tair_time = (/ time /)
  fout->qair_time = (/ time /)
  fout->lrf_time  = (/ time /)
  fout->rain_time  = (/ time /)

  zeros_u = new((/ ocean_time, eta_u, xi_u /), typeof(x_u))
  zeros_u = 0.0+cval
  fout->sustr = (/ zeros_u /)

  zeros_v = new((/ ocean_time, eta_v, xi_v /), typeof(x_v))
  zeros_v = 0.0+cval
  fout->svstr = (/ zeros_v /)

  zeros_r = new((/ ocean_time, eta_rho, xi_rho /), typeof(x_rho))
  zeros_r = 0.0+cval
  fout->shflux = (/ zeros_r /)
  fout->swrad  = (/ zeros_r /)
  fout->swflux = (/ zeros_r /)
  fout->Pair = (/ zeros_r /)
  fout->Tair = (/ zeros_r /)
  fout->Qair = (/ zeros_r /)
  fout->lwrad = (/ zeros_r /)
  fout->lwrad_down = (/ zeros_r /)
  fout->Uwind = (/ zeros_r /)
  fout->Vwind = (/ zeros_r /)
  fout->rain = (/ zeros_r /)
end
