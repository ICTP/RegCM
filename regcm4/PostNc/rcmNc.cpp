/***************************************************************************
 *   Copyright (C) 2008-2009 by Graziano Giuliani                          *
 *   graziano.giuliani at aquila.infn.it                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details. (see COPYING)            *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 *                                                                         *
 *   LIC: GPL                                                              *
 *                                                                         *
 ***************************************************************************/

#include <rcmNc.h>
#include <ctime>

using namespace rcm;

rcmNc::rcmNc(char *fname, char *experiment, header_data &outhead)
{
  f = new NcFile(fname, NcFile::Replace);

  f->add_att("title", "ICTP Regional Climatic model V4 output");
  f->add_att("institution", "ICTP Trieste");
  f->add_att("source", "RegCM Model simulation");
  f->add_att("Conventions", "CF-1.4");
  char buffer[256];
  time_t xtime = time(&xtime);
  snprintf(buffer, 256, "%s", ctime(&xtime));
  char *p = strchr(buffer, '\n');
  *p = 0;
  strncat(buffer, ": Created from model output", 256);
  f->add_att("history", buffer);
  f->add_att("references", "http://users.ictp.it/RegCNET");
  snprintf(buffer, 256, "Experiment Name is : %s", experiment);
  f->add_att("comment", buffer);
  f->add_att("projection", outhead.proj);
  f->add_att("domain_spacing_km", outhead.ds);
  f->add_att("center_latitude_degrees_north", outhead.clat);
  f->add_att("center_longitude_degrees_east", outhead.clon);
  f->add_att("pole_latitude_degrees_north", outhead.xplat);
  f->add_att("pole_longitude_degrees_east", outhead.xplat);
  f->add_att("pole_longitude_degrees_east", outhead.xplon);
  f->add_att("model_boundary_conditions", outhead.lbcs());
  f->add_att("model_cumulous_convection_scheme", outhead.cums());
  f->add_att("model_boundary_layer_scheme", outhead.pbls());
  f->add_att("model_moist_physics_scheme", outhead.moists());

  iy = f->add_dim("iy", outhead.nx);
  jx = f->add_dim("jx", outhead.ny);
  kz = f->add_dim("kz", outhead.kz);
  tt = f->add_dim("time");

  // Add time variable
  timevar = f->add_var("time", ncDouble, tt);
  timevar->add_att("long_name", "time");
  timevar->add_att("standard_name", "time");
  timevar->add_att("calendar", "julian");
  timevar->add_att("units", "hours since 1970-01-01 00:00:00 UTC");

  // Add time independent variables and time itself
  NcVar *sigfvar = f->add_var("level", ncFloat, kz);
  sigfvar->add_att("standard_name", "atmosphere_sigma_coordinate");
  sigfvar->add_att("long_name", "Sigma at layer midpoints");
  sigfvar->add_att("positive", "down");
  sigfvar->add_att("units", "1");
  sigfvar->add_att("axis", "Z");
  sigfvar->add_att("formula_terms", "sigma: level ps: psa ptop: ptop");
  NcVar *xlatvar = f->add_var("xlat", ncFloat, iy, jx);
  xlatvar->add_att("standard_name", "latitude");
  xlatvar->add_att("long_name", "Latitude");
  xlatvar->add_att("units", "degrees_north");
  NcVar *xlonvar = f->add_var("xlon", ncFloat, iy, jx);
  xlonvar->add_att("standard_name", "longitude");
  xlonvar->add_att("long_name", "Longitude");
  xlonvar->add_att("units", "degrees_east");
  NcVar *htvar = f->add_var("ht", ncFloat, iy, jx);
  htvar->add_att("standard_name", "surface_altitude");
  htvar->add_att("long_name", "Domain surface elevation");
  htvar->add_att("coordinates", "xlon xlat");
  htvar->add_att("units", "m");
  NcVar *htsdvar = f->add_var("htsd", ncFloat, iy, jx);
  htsdvar->add_att("standard_name", "surface_altitude");
  htsdvar->add_att("long_name", "Domain elevation stantard deviation");
  htsdvar->add_att("coordinates", "xlon xlat");
  htsdvar->add_att("units", "m");
  htsdvar->add_att("cell_method", "area: standard_deviation");
  NcVar *veg2dvar = f->add_var("veg2d", ncFloat, iy, jx);
  veg2dvar->add_att("long_name", "Vegetation category as defined in BATS");
  veg2dvar->add_att("units", "1");
  veg2dvar->add_att("coordinates", "xlon xlat");
  NcVar *landusevar = f->add_var("landuse", ncFloat, iy, jx);
  landusevar->add_att("long_name", "Landuse category as defined in BATS");
  landusevar->add_att("standard_name", "soil_type");
  landusevar->add_att("units", "1");
  landusevar->add_att("coordinates", "xlon xlat");
  NcVar *xmapvar = f->add_var("xmap", ncFloat, iy, jx);
  xmapvar->add_att("long_name", "Map factor in domain cross points");
  xmapvar->add_att("units", "1");
  xmapvar->add_att("coordinates", "xlon xlat");
  NcVar *dmapvar = f->add_var("dmap", ncFloat, iy, jx);
  dmapvar->add_att("long_name", "Map factor in domain dot points");
  dmapvar->add_att("units", "1");
  dmapvar->add_att("coordinates", "xlon xlat");
  NcVar *coriolvar = f->add_var("coriol", ncFloat, iy, jx);
  coriolvar->add_att("long_name", "Coriolis parameter");
  coriolvar->add_att("standard_name", "coriolis_parameter");
  coriolvar->add_att("units", "s-1");
  coriolvar->add_att("coordinates", "xlon xlat");
  NcVar *maskvar = f->add_var("mask", ncFloat, iy, jx);
  maskvar->add_att("long_name", "Land Sea mask");
  maskvar->add_att("standard_name", "land_binary_mask");
  maskvar->add_att("units", "1");
  maskvar->add_att("coordinates", "xlon xlat");
  NcVar *ptopvar = f->add_var("ptop", ncFloat);
  ptopvar->add_att("long_name", "Pressure at model top");
  ptopvar->add_att("standard_name", "air_pressure");
  ptopvar->add_att("units", "hPa");

  sigfvar->put(outhead.hsigf, outhead.kz);
  xlatvar->put(outhead.xlat, outhead.nx, outhead.ny);
  xlonvar->put(outhead.xlon, outhead.nx, outhead.ny);
  htvar->put(outhead.ht, outhead.nx, outhead.ny);
  htsdvar->put(outhead.htsd, outhead.nx, outhead.ny);
  veg2dvar->put(outhead.veg2d, outhead.nx, outhead.ny);
  landusevar->put(outhead.landuse, outhead.nx, outhead.ny);
  xmapvar->put(outhead.xmap, outhead.nx, outhead.ny);
  dmapvar->put(outhead.dmap, outhead.nx, outhead.ny);
  coriolvar->put(outhead.coriol, outhead.nx, outhead.ny);
  maskvar->put(outhead.mask, outhead.nx, outhead.ny);
  ptopvar->put(&outhead.ptop, 1);
  f->sync();

  // Manage time setup
  struct tm ref;
  memset(&ref, 0, sizeof(struct tm));
  ref.tm_year = 70;
  time_t tzero = mktime(&ref);
  memset(&ref, 0, sizeof(struct tm));
  unsigned int base = outhead.mdate0;
  unsigned int basey = outhead.mdate0/1000000;
  base = base-basey*1000000;
  ref.tm_year = basey - 1900;
  unsigned int basem = base/10000;
  base = base-basem*10000;
  ref.tm_mon = basem-1;
  unsigned int based = base/100;
  base = base-based*100;
  ref.tm_mday = based-1;
  ref.tm_hour = base;
  time_t tref = mktime(&ref);
  reference_time = (unsigned long) (difftime(tref,tzero)/3600.0);
}

rcmNc::~rcmNc( )
{
  f->close( );
}

rcmNcAtmo::rcmNcAtmo(char *fname, char *experiment, header_data &h)
  : rcmNc(fname, experiment, h)
{
  float fillv = -1e+34;
  psvar = f->add_var("psa", ncFloat, tt, iy, jx);
  psvar->add_att("standard_name", "surface_air_pressure");
  psvar->add_att("long_name", "Surface pressure");
  psvar->add_att("coordinates", "xlon xlat");
  psvar->add_att("units", "hPa");
  tprvar = f->add_var("tpr", ncFloat, tt, iy, jx);
  tprvar->add_att("standard_name", "precipitation_flux");
  tprvar->add_att("long_name", "Total daily precipitation rate");
  tprvar->add_att("coordinates", "xlon xlat");
  tprvar->add_att("units", "kg m-2 day-1");
  tgbvar = f->add_var("tgb", ncFloat, tt, iy, jx);
  tgbvar->add_att("standard_name", "soil_temperature");
  tgbvar->add_att("long_name", "Lower groud temperature in BATS");
  tgbvar->add_att("coordinates", "xlon xlat");
  tgbvar->add_att("units", "K");
  swtvar = f->add_var("swt", ncFloat, tt, iy, jx);
  swtvar->add_att("standard_name", "moisture_content_of_soil_layer");
  swtvar->add_att("_FillValue", fillv);
  swtvar->add_att("long_name", "Total soil water");
  swtvar->add_att("coordinates", "xlon xlat");
  swtvar->add_att("units", "kg m-2");
  rnovar = f->add_var("rno", ncFloat, tt, iy, jx);
  rnovar->add_att("standard_name", "runoff_amount");
  rnovar->add_att("_FillValue", fillv);
  rnovar->add_att("long_name", "Runoff accumulated infiltration");
  rnovar->add_att("coordinates", "xlon xlat");
  rnovar->add_att("units", "kg m-2");
  uvar = f->add_var("u", ncFloat, tt, kz, iy, jx);
  uvar->add_att("standard_name", "eastward_wind");
  uvar->add_att("long_name", "U component (westerly) of wind");
  uvar->add_att("coordinates", "xlon xlat");
  uvar->add_att("units", "m/s");
  vvar = f->add_var("v", ncFloat, tt, kz, iy, jx);
  vvar->add_att("standard_name", "northward_wind");
  vvar->add_att("long_name", "V component (southerly) of wind");
  vvar->add_att("coordinates", "xlon xlat");
  vvar->add_att("units", "m/s");
  ovar = f->add_var("omega", ncFloat, tt, kz, iy, jx);
  ovar->add_att("standard_name", "lagrangian_tendency_of_air_pressure");
  ovar->add_att("long_name", "Pressure velocity");
  ovar->add_att("coordinates", "xlon xlat");
  ovar->add_att("units", "hPa s-1");
  tvar = f->add_var("t", ncFloat, tt, kz, iy, jx);
  tvar->add_att("standard_name", "air_temperature");
  tvar->add_att("long_name", "Temperature");
  tvar->add_att("coordinates", "xlon xlat");
  tvar->add_att("units", "K");
  qvvar = f->add_var("qv", ncFloat, tt, kz, iy, jx);
  qvvar->add_att("standard_name", "humidity_mixing_ratio");
  qvvar->add_att("long_name", "Water vapor mixing ratio");
  qvvar->add_att("coordinates", "xlon xlat");
  qvvar->add_att("units", "kg kg-1");
  qcvar = f->add_var("qc", ncFloat, tt, kz, iy, jx);
  qcvar->add_att("standard_name", "cloud_liquid_water_mixing_ratio");
  qcvar->add_att("long_name", "Cloud water mixing ratio");
  qcvar->add_att("coordinates", "xlon xlat");
  qcvar->add_att("units", "kg kg-1");
  count = 0;
}

void rcmNcAtmo::put_rec(atmodata &a)
{
  double xtime = reference_time + count*a.dt;
  timevar->put_rec(&xtime, count);
  psvar->put_rec(a.psa, count);
  tprvar->put_rec(a.tpr, count);
  tgbvar->put_rec(a.tgb, count);
  swtvar->put_rec(a.swt, count);
  rnovar->put_rec(a.rno, count);
  uvar->put_rec(a.u, count);
  vvar->put_rec(a.v, count);
  ovar->put_rec(a.omega, count);
  tvar->put_rec(a.t, count);
  qvvar->put_rec(a.qv, count);
  qcvar->put_rec(a.qc, count);
  count ++;
  return;
}

rcmNcSrf::rcmNcSrf(char *fname, char *experiment, header_data &h)
  : rcmNc(fname, experiment, h)
{
  float fillv = -1e+34;
  // Add two more vertical dimensions for 10m and 2m hgts
  NcDim *m10 = f->add_dim("m10", 1);
  NcDim *m2 = f->add_dim("m2", 1);
  NcDim *soil = f->add_dim("soil_layer", 2);

  NcVar *m10var = f->add_var("m10", ncFloat, m10);
  m10var->add_att("standard_name", "altitude");
  m10var->add_att("long_name", "Convenience 10 m elevation level");
  m10var->add_att("positive", "up");
  m10var->add_att("units", "m");
  m10var->add_att("axis", "Z");
  NcVar *m2var = f->add_var("m2", ncFloat, m2);
  m2var->add_att("standard_name", "altitude");
  m2var->add_att("long_name", "Convenience 2 m elevation level");
  m2var->add_att("positive", "up");
  m2var->add_att("units", "m");
  m2var->add_att("axis", "Z");
  NcVar *soilvar = f->add_var("layer", ncFloat, soil);
  soilvar->add_att("standard_name", "model_level_number");
  soilvar->add_att("long_name", "Surface and root zone");
  soilvar->add_att("positive", "down");
  soilvar->add_att("units", "1");
  soilvar->add_att("axis", "Z");
  float val[2];
  val[0] = 10.0;
  m10var->put(val, 1);
  val[0] = 2.0;
  m2var->put(val, 1);
  val[0] = 0.0;
  val[1] = 1.0;
  soilvar->put(val, 2);

  // Setup variables
  u10mvar = f->add_var("u10m", ncFloat, tt, m10, iy, jx);
  u10mvar->add_att("standard_name", "eastward_wind");
  u10mvar->add_att("long_name", "10 meters U component (westerly) of wind");
  u10mvar->add_att("coordinates", "xlon xlat");
  u10mvar->add_att("units", "m/s");
  v10mvar = f->add_var("v10m", ncFloat, tt, m10, iy, jx);
  v10mvar->add_att("standard_name", "northward_wind");
  v10mvar->add_att("long_name", "10 meters V component (southerly) of wind");
  v10mvar->add_att("coordinates", "xlon xlat");
  v10mvar->add_att("units", "m/s");
  uvdragvar = f->add_var("uvdrag", ncFloat, tt, iy, jx);
  uvdragvar->add_att("standard_name", "surface_drag_coefficient_in_air");
  uvdragvar->add_att("long_name", "Surface drag stress");
  uvdragvar->add_att("coordinates", "xlon xlat");
  uvdragvar->add_att("units", "1");
  tgvar = f->add_var("tg", ncFloat, tt, iy, jx);
  tgvar->add_att("standard_name", "surface_temperature");
  tgvar->add_att("long_name", "Ground temperature");
  tgvar->add_att("coordinates", "xlon xlat");
  tgvar->add_att("units", "K");
  tlefvar = f->add_var("tlef", ncFloat, tt, iy, jx);
  tlefvar->add_att("standard_name", "canopy_temperature");
  tlefvar->add_att("long_name", "Foliage temperature");
  tlefvar->add_att("coordinates", "xlon xlat");
  tlefvar->add_att("_FillValue", fillv);
  tlefvar->add_att("units", "K");
  t2mvar = f->add_var("t2m", ncFloat, tt, m2, iy, jx);
  t2mvar->add_att("standard_name", "air_temperature");
  t2mvar->add_att("long_name", "2 meters temperature");
  t2mvar->add_att("coordinates", "xlon xlat");
  t2mvar->add_att("units", "K");
  q2mvar = f->add_var("q2m", ncFloat, tt, m2, iy, jx);
  q2mvar->add_att("standard_name", "humidity_mixing_ratio");
  q2mvar->add_att("long_name", "2 meters vapour mixing ratio");
  q2mvar->add_att("coordinates", "xlon xlat");
  q2mvar->add_att("units", "kg kg-1");
  smwvar = f->add_var("smw", ncFloat, tt, soil, iy, jx);
  smwvar->add_att("standard_name", "soil_moisture_content");
  smwvar->add_att("long_name", "Moisture content");
  smwvar->add_att("coordinates", "xlon xlat");
  smwvar->add_att("_FillValue", fillv);
  smwvar->add_att("units", "kg m-2");
  tprvar = f->add_var("tpr", ncFloat, tt, iy, jx);
  tprvar->add_att("standard_name", "precipitation_amount");
  tprvar->add_att("long_name", "Total precipitation");
  tprvar->add_att("coordinates", "xlon xlat");
  tprvar->add_att("units", "kg m-2");
  evpvar = f->add_var("evp", ncFloat, tt, iy, jx);
  evpvar->add_att("standard_name", "water_evaporation_amount");
  evpvar->add_att("long_name", "Total evapotranspiration");
  evpvar->add_att("coordinates", "xlon xlat");
  evpvar->add_att("units", "kg m-2");
  runoffvar = f->add_var("runoff", ncFloat, tt, iy, jx);
  runoffvar->add_att("standard_name", "surface_runoff_flux");
  runoffvar->add_att("long_name", "surface runoff");
  runoffvar->add_att("coordinates", "xlon xlat");
  runoffvar->add_att("_FillValue", fillv);
  runoffvar->add_att("units", "kg m-2 day-1");
  scvvar = f->add_var("scv", ncFloat, tt, iy, jx);
  scvvar->add_att("standard_name", "snowfall_amount");
  scvvar->add_att("long_name", "Snow precipitation");
  scvvar->add_att("coordinates", "xlon xlat");
  scvvar->add_att("_FillValue", fillv);
  scvvar->add_att("units", "kg m-2");
  senavar = f->add_var("sena", ncFloat, tt, iy, jx);
  senavar->add_att("standard_name", "surface_downward_sensible_heat_flux");
  senavar->add_att("long_name", "Sensible heat flux");
  senavar->add_att("coordinates", "xlon xlat");
  senavar->add_att("units", "W m-2");
  flwvar = f->add_var("flw", ncFloat, tt, iy, jx);
  flwvar->add_att("standard_name", "net_upward_longwave_flux_in_air");
  flwvar->add_att("long_name", "Net infrared energy flux");
  flwvar->add_att("coordinates", "xlon xlat");
  flwvar->add_att("units", "W m-2");
  fswvar = f->add_var("fsw", ncFloat, tt, iy, jx);
  fswvar->add_att("standard_name", "surface_downwelling_shortwave_flux_in_air");
  fswvar->add_att("long_name", "Solar absorbed energy flux");
  fswvar->add_att("coordinates", "xlon xlat");
  fswvar->add_att("units", "W m-2");
  sinavar = f->add_var("sina", ncFloat, tt, iy, jx);
  sinavar->add_att("standard_name",
                   "net_downward_radiative_flux_at_top_of_atmosphere_model");
  sinavar->add_att("long_name", "Incident solar energy flux");
  sinavar->add_att("coordinates", "xlon xlat");
  sinavar->add_att("units", "W m-2");
  prcvvar = f->add_var("prcv", ncFloat, tt, iy, jx);
  prcvvar->add_att("standard_name", "convective_rainfall_flux");
  prcvvar->add_att("long_name", "Convective precipitation");
  prcvvar->add_att("coordinates", "xlon xlat");
  prcvvar->add_att("units", "kg m-2 day-1");
  psbvar = f->add_var("psa", ncFloat, tt, iy, jx);
  psbvar->add_att("standard_name", "air_pressure");
  psbvar->add_att("long_name", "Surface pressure");
  psbvar->add_att("coordinates", "xlon xlat");
  psbvar->add_att("units", "hPa");
  zpblvar = f->add_var("zpbl", ncFloat, tt, iy, jx);
  zpblvar->add_att("standard_name", "atmosphere_boundary_layer_thickness");
  zpblvar->add_att("long_name", "PBL layer thickness");
  zpblvar->add_att("coordinates", "xlon xlat");
  zpblvar->add_att("units", "m");
  tgmaxvar = f->add_var("tgmax", ncFloat, tt, iy, jx);
  tgmaxvar->add_att("standard_name", "surface_temperature");
  tgmaxvar->add_att("long_name", "Maximum surface temperature");
  tgmaxvar->add_att("coordinates", "xlon xlat");
  tgmaxvar->add_att("cell_methods", "time: maximum within days");
  tgmaxvar->add_att("units", "K");
  tgminvar = f->add_var("tgmin", ncFloat, tt, iy, jx);
  tgminvar->add_att("standard_name", "surface_temperature");
  tgminvar->add_att("long_name", "Maximum surface temperature");
  tgminvar->add_att("coordinates", "xlon xlat");
  tgminvar->add_att("cell_methods", "time: minimum within days");
  tgminvar->add_att("units", "K");
  t2maxvar = f->add_var("t2max", ncFloat, tt, m2, iy, jx);
  t2maxvar->add_att("standard_name", "air_temperature");
  t2maxvar->add_att("long_name", "Maximum 2 meters temperature");
  t2maxvar->add_att("coordinates", "xlon xlat");
  t2maxvar->add_att("cell_methods", "time: maximum within days");
  t2maxvar->add_att("units", "K");
  t2minvar = f->add_var("t2min", ncFloat, tt, m2, iy, jx);
  t2minvar->add_att("standard_name", "air_temperature");
  t2minvar->add_att("long_name", "Minimum 2 meters temperature");
  t2minvar->add_att("coordinates", "xlon xlat");
  t2minvar->add_att("cell_methods", "time: minimum within days");
  t2minvar->add_att("units", "K");
  w10maxvar = f->add_var("w10max", ncFloat, tt, m10, iy, jx);
  w10maxvar->add_att("standard_name", "wind_speed");
  w10maxvar->add_att("long_name", "Maximum speed of 10m wind");
  w10maxvar->add_att("coordinates", "xlon xlat");
  w10maxvar->add_att("cell_methods", "time: maximum within days");
  w10maxvar->add_att("units", "m/s");
  ps_minvar = f->add_var("ps_min", ncFloat, tt, iy, jx);
  ps_minvar->add_att("standard_name", "air_pressure");
  ps_minvar->add_att("long_name", "Surface pressure");
  ps_minvar->add_att("coordinates", "xlon xlat");
  ps_minvar->add_att("cell_methods", "time: minimum within days");
  ps_minvar->add_att("units", "hPa");
  count = 0;
}

void rcmNcSrf::put_rec(srfdata &s)
{
  double xtime = reference_time + count*s.dt;
  timevar->put_rec(&xtime, count);
  u10mvar->put_rec(s.u10m, count);
  v10mvar->put_rec(s.v10m, count);
  uvdragvar->put_rec(s.uvdrag, count);
  tgvar->put_rec(s.tg, count);
  tlefvar->put_rec(s.tlef, count);
  t2mvar->put_rec(s.t2m, count);
  q2mvar->put_rec(s.q2m, count);
  smwvar->put_rec(s.smw, count);
  tprvar->put_rec(s.tpr, count);
  evpvar->put_rec(s.evp, count);
  runoffvar->put_rec(s.runoff, count);
  scvvar->put_rec(s.scv, count);
  senavar->put_rec(s.sena, count);
  flwvar->put_rec(s.flw, count);
  fswvar->put_rec(s.fsw, count);
  sinavar->put_rec(s.sina, count);
  prcvvar->put_rec(s.prcv, count);
  psbvar->put_rec(s.psb, count);
  zpblvar->put_rec(s.zpbl, count);
  tgmaxvar->put_rec(s.tgmax, count);
  tgminvar->put_rec(s.tgmin, count);
  t2maxvar->put_rec(s.t2max, count);
  t2minvar->put_rec(s.t2min, count);
  w10maxvar->put_rec(s.w10max, count);
  ps_minvar->put_rec(s.ps_min, count);
  count ++;
  return;
}

