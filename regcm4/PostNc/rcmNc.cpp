/***************************************************************************
 *   Copyright (C) 2010 Graziano Giuliani                                  *
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
#include <map>
#include <algorithm>
#include <string>
#include <locale.h>

using namespace rcm;

void regcmvar::addvar(std::string name)
{
  std::string::size_type lp = name.find_first_not_of(",",0);
  std::string::size_type pos = name.find_first_of(",",lp);
  while (std::string::npos != pos || std::string::npos != lp)
  {
    rcmname.push_back(name.substr(lp,pos-lp));
    lp = name.find_first_not_of(",",pos);
    pos = name.find_first_of(",", lp);
  }
  return;
}

bool regcmvar::isthere(std::string name)
{
  std::list<std::string>::iterator result = rcmname.begin( );
  result = std::find(result, rcmname.end(), name);
  if (result != rcmname.end()) return true;
  return false;
}

rcmNc::rcmNc(regcmout &fnc, header_data &outhead, bool full)
{
  f = new NcFile(fnc.fname.c_str(), NcFile::Replace);

  f->add_att("title", "ICTP Regional Climatic model V4 output");
  f->add_att("institution", "ICTP");
  f->add_att("source", "RegCM V4 Model simulation");
  f->add_att("Conventions", "CF-1.4");
  char buffer[256];
  time_t xtime = time(&xtime);
  snprintf(buffer, 256, "%s", ctime(&xtime));
  char *p = strchr(buffer, '\n');
  *p = 0;
  strncat(buffer, ": Created from model output", 256);
  f->add_att("history", buffer);
  f->add_att("references", "http://users.ictp.it/RegCNET");
  snprintf(buffer, 256, "Experiment Name is : %s", fnc.experiment.c_str());
  f->add_att("comment", buffer);
  f->add_att("projection", outhead.proj);
  f->add_att("grid_size_in_meters", outhead.ds);
  f->add_att("latitude_of_projection_origin", outhead.clat);
  f->add_att("longitude_of_projection_origin", outhead.clon);
  if (strcmp(outhead.proj, "POLSTR") == 0)
  {
    f->add_att("latitude_of_projection_pole", outhead.xplat);
    f->add_att("longitude_of_projection_pole", outhead.xplon);
  }
  if (strcmp(outhead.proj, "LAMCON") == 0)
  {
    double stp[2];
    stp[0] = outhead.trlat1;
    stp[1] = outhead.trlat2;
    f->add_att("standard_parallel", 2, stp);
  }
  f->add_att("model_boundary_conditions", outhead.lbcs());
  f->add_att("model_cumulous_convection_scheme", outhead.cums());
  if (outhead.icup == 2)
    f->add_att("model_convective_closure_assumption", outhead.cumsclos());
  f->add_att("model_boundary_layer_scheme", outhead.pbls());
  f->add_att("model_moist_physics_scheme", outhead.moists());
  f->add_att("model_ocean_flux_scheme", outhead.ocnflxs());
  f->add_att("model_pressure_gradient_force_scheme", outhead.pgfs());
  f->add_att("model_use_emission_factor", outhead.emisss());
  f->add_att("model_use_lake_model", outhead.lakes());
  f->add_att("model_chemistry_dust_aerosol", outhead.chems());
  f->add_att("model_simulation_initial_start", outhead.mdate0);
  f->add_att("model_simulation_start", outhead.idate1);
  f->add_att("model_simulation_expected_end", outhead.idate2);
  f->add_att("model_simulation_is_a_restart", outhead.ifrest);
  f->add_att("model_timestep_in_seconds", outhead.dt);
  f->add_att("model_timestep_in_minutes_solar_rad_calc", outhead.radfrq);
  f->add_att("model_timestep_in_seconds_bats_calc", outhead.abatm);
  f->add_att("model_timestep_in_hours_radiation_calc", outhead.abemh);
  f->add_att("model_timestep_in_hours_boundary_input", outhead.ibdyfrq);

  if (full)
  {
    iy = f->add_dim("iy", outhead.nx);
    jx = f->add_dim("jx", outhead.ny);
  }
  kz = f->add_dim("kz", outhead.kz);
  tt = f->add_dim("time");

  // Add time variable
  timevar = f->add_var("time", ncDouble, tt);
  timevar->add_att("long_name", "time");
  timevar->add_att("standard_name", "time");
  timevar->add_att("calendar", "standard");
  timevar->add_att("units", "hours since 1970-01-01 00:00:00 UTC");

  // Add time independent variables and time itself
  NcVar *iyvar = f->add_var("iy", ncFloat, iy);
  iyvar->add_att("long_name", "y-coordinate in Cartesian system");
  iyvar->add_att("standard_name", "projection_y_coordinate");
  iyvar->add_att("axis", "Y");
  iyvar->add_att("units", "m");
  NcVar *jxvar = f->add_var("jx", ncFloat, jx);
  jxvar->add_att("long_name", "x-coordinate in Cartesian system");
  jxvar->add_att("standard_name", "projection_x_coordinate");
  jxvar->add_att("axis", "X");
  jxvar->add_att("units", "m");
  NcVar *sigfvar = f->add_var("level", ncFloat, kz);
  sigfvar->add_att("standard_name", "atmosphere_sigma_coordinate");
  sigfvar->add_att("long_name", "Sigma at model layer midpoints");
  sigfvar->add_att("positive", "down");
  sigfvar->add_att("units", "1");
  sigfvar->add_att("axis", "Z");
  sigfvar->add_att("formula_terms", "sigma: level ps: psa ptop: ptop");
  NcVar *ptopvar = f->add_var("ptop", ncFloat);
  ptopvar->add_att("long_name", "Pressure at model top");
  ptopvar->add_att("standard_name", "air_pressure");
  ptopvar->add_att("units", "hPa");
  ctl = &(fnc.ctl);
  if (ctl->doit)
  {
    ctl->head("ICTP Regional Climatic model V4 output", -1e+34);
    ctl->set_grid(outhead);
  }

  if (full)
  {
    NcVar *xlatvar = f->add_var("xlat", ncFloat, iy, jx);
    xlatvar->add_att("standard_name", "latitude");
    xlatvar->add_att("long_name", "Latitude");
    xlatvar->add_att("units", "degrees_north");
    NcVar *xlonvar = f->add_var("xlon", ncFloat, iy, jx);
    xlonvar->add_att("standard_name", "longitude");
    xlonvar->add_att("long_name", "Longitude");
    xlonvar->add_att("units", "degrees_east");
    NcVar *htvar = f->add_var("topo", ncFloat, iy, jx);
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
    landusevar->add_att("legend",
           "1  => Crop/mixed farming\n2  => Short grass\n"
           "3  => Evergreen needleleaf tree\n4  => Deciduous needleleaf tree\n"
           "5  => Deciduous broadleaf tree\n6  => Evergreen broadleaf tree\n"
           "7  => Tall grass\n8  => Desert\n9  => Tundra\n"
           "10 => Irrigated Crop\n11 => Semi-desert\n"
           "12 => Ice cap/glacier\n13 => Bog or marsh\n"
           "14 => Inland water\n15 => Ocean\n16 => Evergreen shrub\n"
           "17 => Deciduous shrub\n18 => Mixed Woodland\n"
           "19 => Forest/Field mosaic\n20 => Water and Land mixture");
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
    if (ctl->doit)
    {
      gradsvar gv;
      gv.set("xlat","xlat","Latitude (degrees_north)",0,false);
      ctl->add_var(gv);
      gv.set("xlon","xlon","Longitude (degrees_east)",0,false);
      ctl->add_var(gv);
      gv.set("topo","ht","Surface elevation (m)",0,false);
      ctl->add_var(gv);
      gv.set("htsd","htsd","Surface elevation std. dev. (m)",0,false);
      ctl->add_var(gv);
      gv.set("veg2d","veg2d","Vegetation BATS category. (1)",0,false);
      ctl->add_var(gv);
      gv.set("landuse","landuse","Landuse BATS category. (1)",0,false);
      ctl->add_var(gv);
      gv.set("xmap","xmap","Map factor cross points. (1)",0,false);
      ctl->add_var(gv);
      gv.set("dmap","dmap","Map factor dot points. (1)",0,false);
      ctl->add_var(gv);
      gv.set("coriol","coriol","Coriolis parameter. (s-1)",0,false);
      ctl->add_var(gv);
      gv.set("mask","mask","Land/Sea mask. (1)",0,false);
      ctl->add_var(gv);
    }
  }
  sigfvar->put(outhead.hsigm, outhead.kz);
  float *tmp;
  if (ctl->doit)
  {
    tmp = new float[outhead.kz];
    for (unsigned int i = 0; i < outhead.kz; i ++)
      tmp[i] = outhead.hsigm[i]*1000.0f;
    ctl->set_levs(tmp, outhead.kz);
    delete [] tmp;
  }
  ptopvar->put(&outhead.ptop, 1);
  tmp = new float[outhead.nx];
  tmp[0] = -(((float) outhead.nx-1)/2.0f) * outhead.ds;
  for (unsigned int i = 1; i < outhead.nx; i ++)
    tmp[i] = tmp[i-1] + outhead.ds;
  iyvar->put(tmp, outhead.nx);
  delete [] tmp;
  tmp = new float[outhead.ny];
  tmp[0] = -(((float) outhead.ny-1)/2.0f) * outhead.ds;
  for (unsigned int i = 1; i < outhead.ny; i ++)
    tmp[i] = tmp[i-1] + outhead.ds;
  jxvar->put(tmp, outhead.ny);
  delete [] tmp;

  f->sync();

  rcount = 0;
  tcount = 0;
}

rcmNc::~rcmNc( )
{
  f->close( );
  delete f;
}

rcmNcAtmo::rcmNcAtmo(regcmout &fnc, header_data &h)
  : rcmNc(fnc, h, true)
{
  // Check if var is on wanted list
  // Special all key
  for (int i = 0; i < 10; i ++)
    varmask[i] = false;
  if (fnc.vl.isthere("all"))
  {
    for (int i = 0; i < 10; i ++)
      varmask[i] = true;
  }
  else
  {
    if (fnc.vl.isthere("u")) varmask[0] = true;
    if (fnc.vl.isthere("v")) varmask[1] = true;
    if (fnc.vl.isthere("o")) varmask[2] = true;
    if (fnc.vl.isthere("dv")) varmask[3] = true;
    if (fnc.vl.isthere("vr")) varmask[4] = true;
    if (fnc.vl.isthere("p")) varmask[5] = true;
    if (fnc.vl.isthere("hg")) varmask[6] = true;
    if (fnc.vl.isthere("t")) varmask[7] = true;
    if (fnc.vl.isthere("td")) varmask[8] = true;
    if (fnc.vl.isthere("tp")) varmask[9] = true;
    if (fnc.vl.isthere("rh")) varmask[10] = true;
    if (fnc.vl.isthere("qv")) varmask[11] = true;
    if (fnc.vl.isthere("qc")) varmask[12] = true;
    if (fnc.vl.isthere("tpr")) varmask[13] = true;
    if (fnc.vl.isthere("tgb")) varmask[14] = true;
    if (fnc.vl.isthere("swt")) varmask[15] = true;
    if (fnc.vl.isthere("rno")) varmask[16] = true;
  }

  // Manage time setup
  rcmdate d(h.idate1);
  reference_time = d.unixtime( );
  if (ctl->doit)
    ctl->set_time(ctl->gradstime(d.basey,d.basem,d.based,d.baseh,(int) h.dto));

  float fillv = -1e+34;
  psvar = f->add_var("psa", ncFloat, tt, iy, jx);
  psvar->add_att("standard_name", "surface_air_pressure");
  psvar->add_att("long_name", "Surface pressure");
  psvar->add_att("coordinates", "xlon xlat");
  psvar->add_att("units", "hPa");

  if (varmask[0] || varmask[1])
  {
    uvar = f->add_var("u", ncFloat, tt, kz, iy, jx);
    uvar->add_att("standard_name", "eastward_wind");
    uvar->add_att("long_name", "U component (westerly) of wind");
    uvar->add_att("coordinates", "xlon xlat");
    uvar->add_att("units", "m s-1");
    vvar = f->add_var("v", ncFloat, tt, kz, iy, jx);
    vvar->add_att("standard_name", "northward_wind");
    vvar->add_att("long_name", "V component (southerly) of wind");
    vvar->add_att("coordinates", "xlon xlat");
    vvar->add_att("units", "m s-1");
  }
  if (varmask[2])
  {
    ovar = f->add_var("omega", ncFloat, tt, kz, iy, jx);
    ovar->add_att("standard_name", "lagrangian_tendency_of_air_pressure");
    ovar->add_att("long_name", "Pressure velocity");
    ovar->add_att("coordinates", "xlon xlat");
    ovar->add_att("units", "hPa s-1");
  }
  if (varmask[3])
  {
    dvvar = f->add_var("dv", ncFloat, tt, kz, iy, jx);
    dvvar->add_att("standard_name", "divergence_of_wind");
    dvvar->add_att("long_name", "Wind divergence");
    dvvar->add_att("coordinates", "xlon xlat");
    dvvar->add_att("units", "s-1");
  }
  if (varmask[4])
  {
    vrvar = f->add_var("vr", ncFloat, tt, kz, iy, jx);
    vrvar->add_att("standard_name", "atmosphere_relative_vorticity");
    vrvar->add_att("long_name", "Vorticity");
    vrvar->add_att("coordinates", "xlon xlat");
    vrvar->add_att("units", "s-1");
  }
  if (varmask[5])
  {
    pvar = f->add_var("p", ncFloat, tt, kz, iy, jx);
    pvar->add_att("standard_name", "air_pressure");
    pvar->add_att("long_name", "Pressure");
    pvar->add_att("coordinates", "xlon xlat");
    pvar->add_att("units", "hPa");
  }
  if (varmask[6])
  {
    hgtvar = f->add_var("hgt", ncFloat, tt, kz, iy, jx);
    hgtvar->add_att("standard_name", "geopotential_height");
    hgtvar->add_att("long_name", "Geopotential Height");
    hgtvar->add_att("coordinates", "xlon xlat");
    hgtvar->add_att("units", "m");
  }
  if (varmask[7])
  {
    tvar = f->add_var("t", ncFloat, tt, kz, iy, jx);
    tvar->add_att("standard_name", "air_temperature");
    tvar->add_att("long_name", "Temperature");
    tvar->add_att("coordinates", "xlon xlat");
    tvar->add_att("units", "K");
  }
  if (varmask[8])
  {
    tdvar = f->add_var("td", ncFloat, tt, kz, iy, jx);
    tdvar->add_att("standard_name", "dew_point_temperature");
    tdvar->add_att("long_name", "Dewpoint Temperature");
    tdvar->add_att("coordinates", "xlon xlat");
    tdvar->add_att("units", "K");
  }
  if (varmask[9])
  {
    tpvar = f->add_var("tp", ncFloat, tt, kz, iy, jx);
    tpvar->add_att("standard_name", "air_potential_temperature");
    tpvar->add_att("long_name", "Potential Temperature");
    tpvar->add_att("coordinates", "xlon xlat");
    tpvar->add_att("units", "K");
  }
  if (varmask[10])
  {
    rhvar = f->add_var("rh", ncFloat, tt, kz, iy, jx);
    rhvar->add_att("standard_name", "relative_humidity");
    rhvar->add_att("long_name", "Relative Humidity");
    rhvar->add_att("coordinates", "xlon xlat");
    rhvar->add_att("units", "1");
  }
  if (varmask[11])
  {
    qvvar = f->add_var("qv", ncFloat, tt, kz, iy, jx);
    qvvar->add_att("standard_name", "humidity_mixing_ratio");
    qvvar->add_att("long_name", "Water vapor mixing ratio");
    qvvar->add_att("coordinates", "xlon xlat");
    qvvar->add_att("units", "kg kg-1");
  }
  if (varmask[12])
  {
    qcvar = f->add_var("qc", ncFloat, tt, kz, iy, jx);
    qcvar->add_att("standard_name", "cloud_liquid_water_mixing_ratio");
    qcvar->add_att("long_name", "Cloud water mixing ratio");
    qcvar->add_att("coordinates", "xlon xlat");
    qcvar->add_att("units", "kg kg-1");
  }
  if (varmask[13])
  {
    tprvar = f->add_var("tpr", ncFloat, tt, iy, jx);
    tprvar->add_att("standard_name", "precipitation_flux");
    tprvar->add_att("long_name", "Total daily precipitation rate");
    tprvar->add_att("coordinates", "xlon xlat");
    tprvar->add_att("units", "kg m-2 day-1");
  }
  if (varmask[14])
  {
    tgbvar = f->add_var("tgb", ncFloat, tt, iy, jx);
    tgbvar->add_att("standard_name", "soil_temperature");
    tgbvar->add_att("long_name", "Lower groud temperature in BATS");
    tgbvar->add_att("coordinates", "xlon xlat");
    tgbvar->add_att("units", "K");
  }
  if (varmask[15])
  {
    swtvar = f->add_var("swt", ncFloat, tt, iy, jx);
    swtvar->add_att("standard_name", "moisture_content_of_soil_layer");
    swtvar->add_att("_FillValue", fillv);
    swtvar->add_att("long_name", "Total soil water");
    swtvar->add_att("coordinates", "xlon xlat");
    swtvar->add_att("units", "kg m-2");
  }
  if (varmask[16])
  {
    rnovar = f->add_var("rno", ncFloat, tt, iy, jx);
    rnovar->add_att("standard_name", "runoff_amount");
    rnovar->add_att("_FillValue", fillv);
    rnovar->add_att("long_name", "Runoff accumulated infiltration");
    rnovar->add_att("coordinates", "xlon xlat");
    rnovar->add_att("units", "kg m-2");
  }

  if (ctl->doit)
  {
    gradsvar gv;
    gv.set("psa","psa","Surface pressure (hPa)",0,true);
    ctl->add_var(gv);
    if (varmask[0] || varmask[1])
    {
      ctl->addentry("vectorpairs u,v");
      gv.set("u","u","U component (westerly) of wind (m s-1)",h.nz,true);
      ctl->add_var(gv);
      gv.set("v","v","V component (southerly) of wind (m s-1)",h.nz,true);
      ctl->add_var(gv);
    }
    if (varmask[2])
    {
      gv.set("omega","omega","Pressure velocity (hPa s-1)",h.nz,true);
      ctl->add_var(gv);
    }
    if (varmask[3])
    {
      gv.set("dv","dv","Wind divergence (s-1)",h.nz,true);
      ctl->add_var(gv);
    }
    if (varmask[4])
    {
      gv.set("vr","vr","Wind vorticity (s-1)",h.nz,true);
      ctl->add_var(gv);
    }
    if (varmask[5])
    {
      gv.set("p","p","Pressure (hPa)",h.nz,true);
      ctl->add_var(gv);
    }
    if (varmask[6])
    {
      gv.set("hgt","hgt","Geopotential height (m)",h.nz,true);
      ctl->add_var(gv);
    }
    if (varmask[7])
    {
      gv.set("t","t","Temperature (K)",h.nz,true);
      ctl->add_var(gv);
    }
    if (varmask[8])
    {
      gv.set("td","td","Dew point temperature (K)",h.nz,true);
      ctl->add_var(gv);
    }
    if (varmask[9])
    {
      gv.set("tp","tp","Potential temperature (K)",h.nz,true);
      ctl->add_var(gv);
    }
    if (varmask[10])
    {
      gv.set("rh","rh","Relative humidity (1)",h.nz,true);
      ctl->add_var(gv);
    }
    if (varmask[11])
    {
      gv.set("qv","qv","Water vapor mixing ratio (kg kg-1)",h.nz,true);
      ctl->add_var(gv);
    }
    if (varmask[12])
    {
      gv.set("qc","qc","Cloud water mixing ratio (kg kg-1)",h.nz,true);
      ctl->add_var(gv);
    }
    if (varmask[13])
    {
      gv.set("tpr","tpr","Total precipitation (kg m-2 day-1)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[14])
    {
      gv.set("tgb","tgb","Soil temperature (K)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[15])
    {
      gv.set("swt","swt","Soil moisture (kg m-2)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[16])
    {
      gv.set("rno","rno","Runoff amount (kg m-2)",0,true);
      ctl->add_var(gv);
    }
  }
}

void rcmNcAtmo::put_rec(atmodata &a, t_atm_deriv &d)
{
  double xtime = reference_time + tcount*a.dt;
  timevar->put_rec(&xtime, rcount);
  psvar->put_rec(a.psa, rcount);
  if (varmask[0] || varmask[1])
  {
    uvar->put_rec(a.u, rcount);
    vvar->put_rec(a.v, rcount);
  }
  if (varmask[2]) ovar->put_rec(a.omega, rcount);
  if (varmask[3]) dvvar->put_rec(d.dv, rcount);
  if (varmask[4]) vrvar->put_rec(d.vr, rcount);
  if (varmask[5]) pvar->put_rec(d.p, rcount);
  if (varmask[6]) hgtvar->put_rec(d.hg, rcount);
  if (varmask[7]) tvar->put_rec(a.t, rcount);
  if (varmask[8]) tdvar->put_rec(d.td, rcount);
  if (varmask[9]) tpvar->put_rec(d.tp, rcount);
  if (varmask[10]) rhvar->put_rec(d.rh, rcount);
  if (varmask[11]) qvvar->put_rec(a.qv, rcount);
  if (varmask[12]) qcvar->put_rec(a.qc, rcount);
  if (varmask[13]) tprvar->put_rec(a.tpr, rcount);
  if (varmask[14]) tgbvar->put_rec(a.tgb, rcount);
  if (varmask[15]) swtvar->put_rec(a.swt, rcount);
  if (varmask[16]) rnovar->put_rec(a.rno, rcount);
  rcount ++;
  tcount ++;
  if (ctl->doit)
    ctl->add_time( );
  return;
}

rcmNcSrf::rcmNcSrf(regcmout &fnc, header_data &h)
  : rcmNc(fnc, h, true)
{
  // Check if var is on wanted list
  // Special all key
  for (int i = 0; i < 10; i ++)
    varmask[i] = false;
  if (fnc.vl.isthere("all"))
  {
    for (int i = 0; i < 10; i ++)
      varmask[i] = true;
  }
  else
  {
    if (fnc.vl.isthere("u10")) varmask[0] = true;
    if (fnc.vl.isthere("v10")) varmask[1] = true;
    if (fnc.vl.isthere("udg")) varmask[2] = true;
    if (fnc.vl.isthere("tg")) varmask[3] = true;
    if (fnc.vl.isthere("tfl")) varmask[4] = true;
    if (fnc.vl.isthere("t2")) varmask[5] = true;
    if (fnc.vl.isthere("r2")) varmask[6] = true;
    if (fnc.vl.isthere("q2")) varmask[7] = true;
    if (fnc.vl.isthere("sm")) varmask[8] = true;
    if (fnc.vl.isthere("tpr")) varmask[9] = true;
    if (fnc.vl.isthere("evp")) varmask[10] = true;
    if (fnc.vl.isthere("rno")) varmask[11] = true;
    if (fnc.vl.isthere("scv")) varmask[12] = true;
    if (fnc.vl.isthere("sen")) varmask[13] = true;
    if (fnc.vl.isthere("flw")) varmask[14] = true;
    if (fnc.vl.isthere("fsw")) varmask[15] = true;
    if (fnc.vl.isthere("fld")) varmask[16] = true;
    if (fnc.vl.isthere("sin")) varmask[17] = true;
    if (fnc.vl.isthere("prc")) varmask[18] = true;
    if (fnc.vl.isthere("zpb")) varmask[19] = true;
    if (fnc.vl.isthere("tga")) varmask[20] = true;
    if (fnc.vl.isthere("tgi")) varmask[21] = true;
    if (fnc.vl.isthere("t2a")) varmask[22] = true;
    if (fnc.vl.isthere("t2i")) varmask[23] = true;
    if (fnc.vl.isthere("wma")) varmask[24] = true;
    if (fnc.vl.isthere("psi")) varmask[25] = true;
  }

  // Manage time setup
  rcmdate d(h.idate1);
  reference_time = d.unixtime( );
  if (ctl->doit)
    ctl->set_time(ctl->gradstime(d.basey,d.basem,d.based,d.baseh,(int) h.dtb));

  float fillv = -1e+34;
  // Add three more vertical dimensions for 10m, 2m hgts and soil layers
  NcDim *m10 = f->add_dim("m10", 1);
  NcDim *m2 = f->add_dim("m2", 1);
  NcDim *soil = f->add_dim("soil_layer", 2);
  NcDim *nv = f->add_dim("nv", 2);

  timevar->add_att("bounds", "time_bnds");

  tbnd = f->add_var("time_bnds", ncFloat, tt, nv);
  tbnd->add_att("calendar", "standard");
  tbnd->add_att("units", "hours since 1970-01-01 00:00:00 UTC");

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
  char cell_method[64];
  psbvar = f->add_var("psa", ncFloat, tt, iy, jx);
  psbvar->add_att("standard_name", "air_pressure");
  psbvar->add_att("long_name", "Surface pressure");
  psbvar->add_att("coordinates", "xlon xlat");
  psbvar->add_att("units", "hPa");

  if (varmask[0] || varmask[1])
  {
    u10mvar = f->add_var("u10m", ncFloat, tt, m10, iy, jx);
    u10mvar->add_att("standard_name", "eastward_wind");
    u10mvar->add_att("long_name", "10 meters U component (westerly) of wind");
    u10mvar->add_att("coordinates", "xlon xlat");
    u10mvar->add_att("units", "m s-1");
    v10mvar = f->add_var("v10m", ncFloat, tt, m10, iy, jx);
    v10mvar->add_att("standard_name", "northward_wind");
    v10mvar->add_att("long_name", "10 meters V component (southerly) of wind");
    v10mvar->add_att("coordinates", "xlon xlat");
    v10mvar->add_att("units", "m s-1");
  }
  if (varmask[2])
  {
    uvdragvar = f->add_var("uvdrag", ncFloat, tt, iy, jx);
    uvdragvar->add_att("standard_name", "surface_drag_coefficient_in_air");
    uvdragvar->add_att("long_name", "Surface drag stress");
    uvdragvar->add_att("coordinates", "xlon xlat");
    uvdragvar->add_att("units", "1");
  }
  if (varmask[3])
  {
    tgvar = f->add_var("tg", ncFloat, tt, iy, jx);
    tgvar->add_att("standard_name", "surface_temperature");
    tgvar->add_att("long_name", "Ground temperature");
    tgvar->add_att("coordinates", "xlon xlat");
    tgvar->add_att("units", "K");
  }
  if (varmask[4])
  {
    tlefvar = f->add_var("tlef", ncFloat, tt, iy, jx);
    tlefvar->add_att("standard_name", "canopy_temperature");
    tlefvar->add_att("long_name", "Foliage temperature");
    tlefvar->add_att("coordinates", "xlon xlat");
    tlefvar->add_att("_FillValue", fillv);
    tlefvar->add_att("units", "K");
  }
  if (varmask[5])
  {
    t2mvar = f->add_var("t2m", ncFloat, tt, m2, iy, jx);
    t2mvar->add_att("standard_name", "air_temperature");
    t2mvar->add_att("long_name", "2 meters temperature");
    t2mvar->add_att("coordinates", "xlon xlat");
    t2mvar->add_att("units", "K");
  }
  if (varmask[6])
  {
    r2mvar = f->add_var("r2m", ncFloat, tt, m2, iy, jx);
    r2mvar->add_att("standard_name", "relative_humidity");
    r2mvar->add_att("long_name", "2 meters relative humidity");
    r2mvar->add_att("coordinates", "xlon xlat");
    r2mvar->add_att("units", "1");
  }
  if (varmask[7])
  {
    q2mvar = f->add_var("q2m", ncFloat, tt, m2, iy, jx);
    q2mvar->add_att("standard_name", "humidity_mixing_ratio");
    q2mvar->add_att("long_name", "2 meters vapour mixing ratio");
    q2mvar->add_att("coordinates", "xlon xlat");
    q2mvar->add_att("units", "kg kg-1");
  }
  if (varmask[8])
  {
    smwvar = f->add_var("smw", ncFloat, tt, soil, iy, jx);
    smwvar->add_att("standard_name", "soil_moisture_content");
    smwvar->add_att("long_name", "Moisture content");
    smwvar->add_att("coordinates", "xlon xlat");
    smwvar->add_att("_FillValue", fillv);
    smwvar->add_att("units", "kg m-2");
  }
  if (varmask[9])
  {
    tprvar = f->add_var("tpr", ncFloat, tt, iy, jx);
    tprvar->add_att("standard_name", "precipitation_amount");
    tprvar->add_att("long_name", "Total precipitation");
    tprvar->add_att("coordinates", "xlon xlat");
    tprvar->add_att("units", "kg m-2");
  }
  if (varmask[10])
  {
    evpvar = f->add_var("evp", ncFloat, tt, iy, jx);
    evpvar->add_att("standard_name", "water_evaporation_amount");
    evpvar->add_att("long_name", "Total evapotranspiration");
    evpvar->add_att("coordinates", "xlon xlat");
    evpvar->add_att("units", "kg m-2");
  }
  if (varmask[11])
  {
    runoffvar = f->add_var("runoff", ncFloat, tt, iy, jx);
    runoffvar->add_att("standard_name", "surface_runoff_flux");
    runoffvar->add_att("long_name", "surface runoff");
    runoffvar->add_att("coordinates", "xlon xlat");
    runoffvar->add_att("_FillValue", fillv);
    runoffvar->add_att("units", "kg m-2 day-1");
  }
  if (varmask[12])
  {
    scvvar = f->add_var("scv", ncFloat, tt, iy, jx);
    scvvar->add_att("standard_name", "snowfall_amount");
    scvvar->add_att("long_name", "Snow precipitation");
    scvvar->add_att("coordinates", "xlon xlat");
    scvvar->add_att("_FillValue", fillv);
    scvvar->add_att("units", "kg m-2");
  }
  if (varmask[13])
  {
    senavar = f->add_var("sena", ncFloat, tt, iy, jx);
    senavar->add_att("standard_name", "surface_downward_sensible_heat_flux");
    senavar->add_att("long_name", "Sensible heat flux");
    senavar->add_att("coordinates", "xlon xlat");
    senavar->add_att("units", "W m-2");
  }
  if (varmask[14])
  {
    flwvar = f->add_var("flw", ncFloat, tt, iy, jx);
    flwvar->add_att("standard_name", "net_upward_longwave_flux_in_air");
    flwvar->add_att("long_name", "Net infrared energy flux");
    flwvar->add_att("coordinates", "xlon xlat");
    flwvar->add_att("units", "W m-2");
  }
  if (varmask[15])
  {
    fswvar = f->add_var("fsw", ncFloat, tt, iy, jx);
    fswvar->add_att("standard_name",
                    "surface_downwelling_shortwave_flux_in_air");
    fswvar->add_att("long_name", "Solar absorbed energy flux");
    fswvar->add_att("coordinates", "xlon xlat");
    fswvar->add_att("units", "W m-2");
  }
  if (varmask[16])
  {
    flwdvar = f->add_var("fld", ncFloat, tt, iy, jx);
    flwdvar->add_att("standard_name",
                    "surface_downwelling_longwave_flux_in_air");
    flwdvar->add_att("long_name", "Downward LW flux");
    flwdvar->add_att("coordinates", "xlon xlat");
    flwdvar->add_att("units", "W m-2");
  }
  if (varmask[17])
  {
    sinavar = f->add_var("sina", ncFloat, tt, iy, jx);
    sinavar->add_att("standard_name",
                     "net_downward_radiative_flux_at_top_of_atmosphere_model");
    sinavar->add_att("long_name", "Incident solar energy flux");
    sinavar->add_att("coordinates", "xlon xlat");
    sinavar->add_att("units", "W m-2");
  }
  if (varmask[18])
  {
    prcvvar = f->add_var("prcv", ncFloat, tt, iy, jx);
    prcvvar->add_att("standard_name", "convective_rainfall_flux");
    prcvvar->add_att("long_name", "Convective precipitation");
    prcvvar->add_att("coordinates", "xlon xlat");
    prcvvar->add_att("units", "kg m-2 day-1");
  }
  if (varmask[19])
  {
    zpblvar = f->add_var("zpbl", ncFloat, tt, iy, jx);
    zpblvar->add_att("standard_name", "atmosphere_boundary_layer_thickness");
    zpblvar->add_att("long_name", "PBL layer thickness");
    zpblvar->add_att("coordinates", "xlon xlat");
    zpblvar->add_att("units", "m");
  }
  if (varmask[20])
  {
    tgmaxvar = f->add_var("tgmax", ncFloat, tt, iy, jx);
    tgmaxvar->add_att("standard_name", "surface_temperature");
    tgmaxvar->add_att("long_name", "Maximum surface temperature");
    tgmaxvar->add_att("coordinates", "xlon xlat");
    sprintf(cell_method, "time: maximum (interval: %d hour)", (int) h.dtb);
    tgmaxvar->add_att("cell_methods", cell_method);
    tgmaxvar->add_att("units", "K");
  }
  if (varmask[21])
  {
    tgminvar = f->add_var("tgmin", ncFloat, tt, iy, jx);
    tgminvar->add_att("standard_name", "surface_temperature");
    tgminvar->add_att("long_name", "Maximum surface temperature");
    tgminvar->add_att("coordinates", "xlon xlat");
    sprintf(cell_method, "time: minimum (interval: %d hour)", (int) h.dtb);
    tgminvar->add_att("cell_methods", cell_method);
    tgminvar->add_att("units", "K");
  }
  if (varmask[22])
  {
    t2maxvar = f->add_var("t2max", ncFloat, tt, m2, iy, jx);
    t2maxvar->add_att("standard_name", "air_temperature");
    t2maxvar->add_att("long_name", "Maximum 2 meters temperature");
    t2maxvar->add_att("coordinates", "xlon xlat");
    sprintf(cell_method, "time: maximum (interval: %d hour)", (int) h.dtb);
    t2maxvar->add_att("cell_methods", cell_method);
    t2maxvar->add_att("units", "K");
  }
  if (varmask[23])
  {
    t2minvar = f->add_var("t2min", ncFloat, tt, m2, iy, jx);
    t2minvar->add_att("standard_name", "air_temperature");
    t2minvar->add_att("long_name", "Minimum 2 meters temperature");
    t2minvar->add_att("coordinates", "xlon xlat");
    sprintf(cell_method, "time: minimum (interval: %d hour)", (int) h.dtb);
    t2minvar->add_att("cell_methods", cell_method);
    t2minvar->add_att("units", "K");
  }
  if (varmask[24])
  {
    w10maxvar = f->add_var("w10max", ncFloat, tt, m10, iy, jx);
    w10maxvar->add_att("standard_name", "wind_speed");
    w10maxvar->add_att("long_name", "Maximum speed of 10m wind");
    w10maxvar->add_att("coordinates", "xlon xlat");
    sprintf(cell_method, "time: maximum (interval: %d hour)", (int) h.dtb);
    w10maxvar->add_att("cell_methods", cell_method);
    w10maxvar->add_att("units", "m s-1");
  }
  if (varmask[25])
  {
    ps_minvar = f->add_var("ps_min", ncFloat, tt, iy, jx);
    ps_minvar->add_att("standard_name", "air_pressure");
    ps_minvar->add_att("long_name", "Surface pressure");
    ps_minvar->add_att("coordinates", "xlon xlat");
    sprintf(cell_method, "time: minimum (interval: %d hour)", (int) h.dtb);
    ps_minvar->add_att("cell_methods", cell_method);
    ps_minvar->add_att("units", "hPa");
  }
  last_time = reference_time;
  if (ctl->doit)
  {
    gradsvar gv;
    gv.set("psa","psa","Surface pressure (hPa)",0,true);
    ctl->add_var(gv);
    if (varmask[1] || varmask[2])
    {
      ctl->addentry("vectorpairs u10,v10");
      gv.set("u10","u10","10m U component (westerly) of wind (m s-1)",0,true);
      ctl->add_var(gv);
      gv.set("v10","v10","10m V component (southerly) of wind (m s-1)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[2])
    {
      gv.set("uvdrag","uvdrag","Surface drag stress (1)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[3])
    {
      gv.set("tg","tg","Ground temperature (K)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[4])
    {
      gv.set("tlef","tlef","Foliage temperature (K)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[5])
    {
      gv.set("t2m","t2m","2m air temperature (K)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[6])
    {
      gv.set("r2m","r2m","2m relative humidity (1)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[7])
    {
      gv.set("q2m","q2m","2m vapor mixing ratio (kg kg-1)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[8])
    {
      std::cout << "Warning: Variable smw is not plottable by GrADS."
                << std::endl
                << "Is not added to CTL file, albeit present in NetCDF"
                << std::endl;
    }
    if (varmask[9])
    {
      gv.set("tpr","tpr","Total precipitation amount (kg m-2)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[10])
    {
      gv.set("evp","evp","Total evapotranspiration (kg m-2)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[11])
    {
      gv.set("runoff","runoff","Surface runoff flux (kg m-2 day-1)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[12])
    {
      gv.set("scv","scv","Snowfall amount (kg m-2)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[13])
    {
      gv.set("sena","sena","Sensible heat flux (W m-2)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[14])
    {
      gv.set("flw","flw","Net infrared energy flux (W m-2)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[15])
    {
      gv.set("fsw","fsw","Net solar absorbed energy flux (W m-2)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[16])
    {
      gv.set("fld","fld","Net downward LW flux (W m-2)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[17])
    {
      gv.set("sina","sina","Incident solar energy flux (W m-2)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[18])
    {
      gv.set("prcv","prcv","Convective precipitation (kg m-2 day-1)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[19])
    {
      gv.set("zpbl","zpbl","PBL layer thickness (m)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[20])
    {
      gv.set("tgmax","tgmax","Maximum surface temperature (K)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[21])
    {
      gv.set("tgmin","tgmin","Minimum surface temperature (K)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[22])
    {
      gv.set("t2max","t2max","Maximum 2m temperature (K)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[23])
    {
      gv.set("t2min","t2min","Minimum 2m temperature (K)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[24])
    {
      gv.set("w10max","w10max","Maximum wind speed at 10m (m s-1)",0,true);
      ctl->add_var(gv);
    }
    if (varmask[25])
    {
      gv.set("ps_min","ps_min","Minimum surface pressure (hPa)",0,true);
      ctl->add_var(gv);
    }
  }
}

void rcmNcSrf::put_rec(srfdata &s, t_srf_deriv &d)
{
  double xtime[2];
  xtime[0] = reference_time + tcount*s.dt;
  xtime[1] = (double) last_time;
  timevar->put_rec(xtime, rcount);
  tbnd->put_rec(xtime, rcount);
  last_time = (unsigned long) xtime[0];
  psbvar->put_rec(s.psb, rcount);
  if (varmask[0] || varmask[1])
  {
    u10mvar->put_rec(s.u10m, rcount);
    v10mvar->put_rec(s.v10m, rcount);
  }
  if (varmask[2]) uvdragvar->put_rec(s.uvdrag, rcount);
  if (varmask[3]) tgvar->put_rec(s.tg, rcount);
  if (varmask[4]) tlefvar->put_rec(s.tlef, rcount);
  if (varmask[5]) t2mvar->put_rec(s.t2m, rcount);
  if (varmask[6]) r2mvar->put_rec(d.r2, rcount);
  if (varmask[7]) q2mvar->put_rec(s.q2m, rcount);
  if (varmask[8]) smwvar->put_rec(s.smw, rcount);
  if (varmask[9]) tprvar->put_rec(s.tpr, rcount);
  if (varmask[10]) evpvar->put_rec(s.evp, rcount);
  if (varmask[11]) runoffvar->put_rec(s.runoff, rcount);
  if (varmask[12]) scvvar->put_rec(s.scv, rcount);
  if (varmask[13]) senavar->put_rec(s.sena, rcount);
  if (varmask[14]) flwvar->put_rec(s.flw, rcount);
  if (varmask[15]) fswvar->put_rec(s.fsw, rcount);
  if (varmask[16]) flwdvar->put_rec(s.flwd, rcount);
  if (varmask[17]) sinavar->put_rec(s.sina, rcount);
  if (varmask[18]) prcvvar->put_rec(s.prcv, rcount);
  if (varmask[19]) zpblvar->put_rec(s.zpbl, rcount);
  if (varmask[20]) tgmaxvar->put_rec(s.tgmax, rcount);
  if (varmask[21]) tgminvar->put_rec(s.tgmin, rcount);
  if (varmask[22]) t2maxvar->put_rec(s.t2max, rcount);
  if (varmask[23]) t2minvar->put_rec(s.t2min, rcount);
  if (varmask[24]) w10maxvar->put_rec(s.w10max, rcount);
  if (varmask[25]) ps_minvar->put_rec(s.ps_min, rcount);
  rcount ++;
  tcount ++;
  if (ctl->doit)
    ctl->add_time( );
  return;
}

rcmNcRad::rcmNcRad(regcmout &fnc, header_data &h)
  : rcmNc(fnc, h, true)
{
  // Manage time setup
  rcmdate d(h.idate1);
  reference_time = d.unixtime( );
  if (ctl->doit)
    ctl->set_time(ctl->gradstime(d.basey,d.basem,d.based,d.baseh,(int) h.dtr));

  // Setup variables
  psvar = f->add_var("psa", ncFloat, tt, iy, jx);
  psvar->add_att("standard_name", "surface_air_pressure");
  psvar->add_att("long_name", "Surface pressure");
  psvar->add_att("coordinates", "xlon xlat");
  psvar->add_att("units", "hPa");
  cldvar = f->add_var("cld", ncFloat, tt, kz, iy, jx);
  cldvar->add_att("standard_name", "cloud_area_fraction_in_atmosphere_layer");
  cldvar->add_att("long_name", "Cloud fractional cover");
  cldvar->add_att("coordinates", "xlon xlat");
  cldvar->add_att("units", "1");
  clwpvar = f->add_var("clwp", ncFloat, tt, kz, iy, jx);
  clwpvar->add_att("standard_name",
                   "atmosphere_optical_thickness_due_to_cloud");
  clwpvar->add_att("long_name", "Cloud liquid water path");
  clwpvar->add_att("coordinates", "xlon xlat");
  clwpvar->add_att("units", "1");
  qrsvar = f->add_var("qrs", ncFloat, tt, kz, iy, jx);
  qrsvar->add_att("standard_name",
                   "tendency_of_air_temperature_due_to_shortwave_heating");
  qrsvar->add_att("long_name", "Solar heating rate");
  qrsvar->add_att("coordinates", "xlon xlat");
  qrsvar->add_att("units", "K s-1");
  qrlvar = f->add_var("qrl", ncFloat, tt, kz, iy, jx);
  qrlvar->add_att("standard_name",
                   "tendency_of_air_temperature_due_to_longwave_heating");
  qrlvar->add_att("long_name", "Longwave cooling rate");
  qrlvar->add_att("coordinates", "xlon xlat");
  qrlvar->add_att("units", "K s-1");
  frsavar = f->add_var("frsa", ncFloat, tt, iy, jx);
  frsavar->add_att("standard_name",
                   "surface_downwelling_shortwave_flux_in_air");
  frsavar->add_att("long_name", "Surface absorbed solar flux");
  frsavar->add_att("coordinates", "xlon xlat");
  frsavar->add_att("units", "W m-2");
  frlavar = f->add_var("frla", ncFloat, tt, iy, jx);
  frlavar->add_att("standard_name", "downwelling_longwave_flux_in_air");
  frlavar->add_att("long_name", "Longwave cooling of surface flux");
  frlavar->add_att("coordinates", "xlon xlat");
  frlavar->add_att("units", "W m-2");
  clrstvar = f->add_var("clrst", ncFloat, tt, iy, jx);
  clrstvar->add_att("standard_name",
                   "downwelling_shortwave_flux_in_air_assuming_clear_sky");
  clrstvar->add_att("long_name", "clearsky total column absorbed solar flux");
  clrstvar->add_att("coordinates", "xlon xlat");
  clrstvar->add_att("units", "W m-2");
  clrssvar = f->add_var("clrss", ncFloat, tt, iy, jx);
  clrssvar->add_att("standard_name",
                   "net_downward_shortwave_flux_in_air_assuming_clear_sky");
  clrssvar->add_att("long_name", "clearsky surface absorbed solar flux");
  clrssvar->add_att("coordinates", "xlon xlat");
  clrssvar->add_att("units", "W m-2");
  clrltvar = f->add_var("clrlt", ncFloat, tt, iy, jx);
  clrltvar->add_att("standard_name",
                   "toa_net_upward_longwave_flux_assuming_clear_sky");
  clrltvar->add_att("long_name", "clearsky net upward LW flux at TOA");
  clrltvar->add_att("coordinates", "xlon xlat");
  clrltvar->add_att("units", "W m-2");
  clrlsvar = f->add_var("clrls", ncFloat, tt, iy, jx);
  clrlsvar->add_att("standard_name",
                   "net_upward_longwave_flux_in_air_assuming_clear_sky");
  clrlsvar->add_att("long_name", "clearsky LW cooling at surface");
  clrlsvar->add_att("coordinates", "xlon xlat");
  clrlsvar->add_att("units", "W m-2");
  solinvar = f->add_var("solin", ncFloat, tt, iy, jx);
  solinvar->add_att("standard_name", "toa_instantaneous_shortwave_forcing");
  solinvar->add_att("long_name", "Instantaneous incident solar");
  solinvar->add_att("coordinates", "xlon xlat");
  solinvar->add_att("units", "W m-2");
  sabtpvar = f->add_var("sabtp", ncFloat, tt, iy, jx);
  sabtpvar->add_att("standard_name",
                  "atmosphere_net_rate_of_absorption_of_shortwave_energy");
  sabtpvar->add_att("long_name", "Total column absorbed solar flux");
  sabtpvar->add_att("coordinates", "xlon xlat");
  sabtpvar->add_att("units", "W m-2");
  firtpvar = f->add_var("firtp", ncFloat, tt, iy, jx);
  firtpvar->add_att("standard_name",
                  "atmosphere_net_rate_of_absorption_of_longwave_energy");
  firtpvar->add_att("long_name", "net upward LW flux at TOA");
  firtpvar->add_att("coordinates", "xlon xlat");
  firtpvar->add_att("units", "W m-2");
  if (ctl->doit)
  {
    gradsvar gv;
    gv.set("psa","psa","Surface pressure (hPa)",0,true);
    ctl->add_var(gv);
    gv.set("cld","cld","Cloud fractional cover (1)",h.nz,true);
    ctl->add_var(gv);
    gv.set("clwp","clwp","Cloud liquid water path (1)",h.nz,true);
    ctl->add_var(gv);
    gv.set("qrs","qrs","Solar heating rate (K s-1)",h.nz,true);
    ctl->add_var(gv);
    gv.set("qrl","qrl","Longwave cooling rate (K s-1)",h.nz,true);
    ctl->add_var(gv);
    gv.set("frsa","frsa","Surface absorbed solar flux (W m-2)",0,true);
    ctl->add_var(gv);
    gv.set("frla","frla","Longwave cooling flux (W m-2)",0,true);
    ctl->add_var(gv);
    gv.set("clrst","clrst","Clearsky total abs. solar flux (W m-2)",0,true);
    ctl->add_var(gv);
    gv.set("clrss","clrss","Clearsky surface abs. solar flux (W m-2)",0,true);
    ctl->add_var(gv);
    gv.set("clrlt","clrlt","Clearsky net upward LW flux TOA (W m-2)",0,true);
    ctl->add_var(gv);
    gv.set("clrls","clrls","Clearsky LW surface cooling flux (W m-2)",0,true);
    ctl->add_var(gv);
    gv.set("solin","solin","Instantaneous incident solar flux (W m-2)",0,true);
    ctl->add_var(gv);
    gv.set("sabtp","sabtp","Total column absorbed solar flux (W m-2)",0,true);
    ctl->add_var(gv);
    gv.set("firtp","firtp","Net upward LW flux TOA (W m-2)",0,true);
    ctl->add_var(gv);
  }
}

void rcmNcRad::put_rec(raddata &r)
{
  double xtime = reference_time + tcount*r.dt;
  timevar->put_rec(&xtime, rcount);
  psvar->put_rec(r.psa, rcount);
  cldvar->put_rec(r.cld, rcount);
  clwpvar->put_rec(r.clwp, rcount);
  qrsvar->put_rec(r.qrs, rcount);
  qrlvar->put_rec(r.qrl, rcount);
  frsavar->put_rec(r.frsa, rcount);
  frlavar->put_rec(r.frla, rcount);
  clrstvar->put_rec(r.clrst, rcount);
  clrssvar->put_rec(r.clrss, rcount);
  clrltvar->put_rec(r.clrlt, rcount);
  clrlsvar->put_rec(r.clrls, rcount);
  solinvar->put_rec(r.solin, rcount);
  sabtpvar->put_rec(r.sabtp, rcount);
  firtpvar->put_rec(r.firtp, rcount);
  rcount ++;
  tcount ++;
  if (ctl->doit)
    ctl->add_time( );
  return;
}

rcmNcChe::rcmNcChe(regcmout &fnc, header_data &h)
  : rcmNc(fnc, h, true)
{
  // Manage time setup
  rcmdate d(h.idate1);
  reference_time = d.unixtime( );
  if (ctl->doit)
    ctl->set_time(ctl->gradstime(d.basey,d.basem,d.based,d.baseh,(int) h.dtc));

  nx = h.nx;
  ny = h.ny;

  // This guy does not respect any coventions
  f->add_att("Conventions", "None");
  f->add_att("tracer_model_names", h.trnames.c_str());

  // Add tracer numbers dimension
  trc = f->add_dim("tracer", h.ntr);

  // Setup variables
  psvar = f->add_var("psa", ncFloat, tt, iy, jx);
  psvar->add_att("standard_name", "surface_air_pressure");
  psvar->add_att("long_name", "Surface pressure");
  psvar->add_att("coordinates", "xlon xlat");
  psvar->add_att("units", "hPa");

  trac3Dvar = f->add_var("trac", ncFloat, tt, trc, kz, iy, jx);
  trac3Dvar->add_att("standard_name", "atmosphere_mixing_ratio_of_tracer");
  trac3Dvar->add_att("long_name", "Tracers mixing ratios");
  trac3Dvar->add_att("coordinates", "xlon xlat");
  trac3Dvar->add_att("units", "kg kg-1");
  aext8var = f->add_var("aext8", ncFloat, tt, kz, iy, jx);
  aext8var->add_att("standard_name", "aerosol_extincion_coefficient");
  aext8var->add_att("long_name", "aer mix. ext. coef");
  aext8var->add_att("coordinates", "xlon xlat");
  aext8var->add_att("units", "1");
  assa8var = f->add_var("assa8", ncFloat, tt, kz, iy, jx);
  assa8var->add_att("standard_name", "aerosol_single_scattering_albedo");
  assa8var->add_att("long_name", "aer mix. sin. scat. alb");
  assa8var->add_att("coordinates", "xlon xlat");
  assa8var->add_att("units", "1");
  agfu8var = f->add_var("agfu8", ncFloat, tt, kz, iy, jx);
  agfu8var->add_att("standard_name", "aerosol_asymmetry_parameter");
  agfu8var->add_att("long_name", "aer mix. ass. par");
  agfu8var->add_att("coordinates", "xlon xlat");
  agfu8var->add_att("units", "1");
  colbvar = f->add_var("colb", ncFloat, tt, trc, iy, jx);
  colbvar->add_att("standard_name", "instantaneous_deposition_of_tracer");
  colbvar->add_att("long_name", "columnburden inst");
  colbvar->add_att("coordinates", "xlon xlat");
  colbvar->add_att("units", "mg m-2");
  wdlscvar = f->add_var("wdlsc", ncFloat, tt, trc, iy, jx);
  wdlscvar->add_att("standard_name",
      "tendency_of_wet_deposition_of_tracer_due_to_large_scale_precipitation");
  wdlscvar->add_att("long_name", "wet dep lgscale");
  wdlscvar->add_att("coordinates", "xlon xlat");
  wdlscvar->add_att("units", "mg m-2 day-1");
  wdcvcvar = f->add_var("wdcvc", ncFloat, tt, trc, iy, jx);
  wdcvcvar->add_att("standard_name",
      "tendency_of_wet_deposition_of_tracer_due_to_convective_precipitation");
  wdcvcvar->add_att("long_name", "wet dep convect");
  wdcvcvar->add_att("coordinates", "xlon xlat");
  wdcvcvar->add_att("units", "mg m-2 day-1");
  sdrdpvar = f->add_var("sdrdp", ncFloat, tt, trc, iy, jx);
  sdrdpvar->add_att("standard_name", "tendency_of_dry_deposition_of_tracer");
  sdrdpvar->add_att("long_name", "surf dry depos");
  sdrdpvar->add_att("coordinates", "xlon xlat");
  sdrdpvar->add_att("units", "mg m-2 day-1");
  xgascvar = f->add_var("xgasc", ncFloat, tt, trc, iy, jx);
  xgascvar->add_att("standard_name", "tendency_of_gas_conversion_of_tracer");
  xgascvar->add_att("long_name", "chem gas conv");
  xgascvar->add_att("coordinates", "xlon xlat");
  xgascvar->add_att("units", "mg m-2 day-1");
  xaqucvar = f->add_var("xaquc", ncFloat, tt, trc, iy, jx);
  xaqucvar->add_att("standard_name",
                    "tendency_of_aqueous_conversion_of_tracer");
  xaqucvar->add_att("long_name", "chem aqu conv");
  xaqucvar->add_att("coordinates", "xlon xlat");
  xaqucvar->add_att("units", "mg m-2 day-1");
  emissvar = f->add_var("emiss", ncFloat, tt, trc, iy, jx);
  emissvar->add_att("standard_name",
                    "tendency_of_surface_emission_of_tracer");
  emissvar->add_att("long_name", "surf emission");
  emissvar->add_att("coordinates", "xlon xlat");
  emissvar->add_att("units", "mg m-2 day-1");
  acstoarfvar = f->add_var("acstoarf", ncFloat, tt, iy, jx);
  acstoarfvar->add_att("standard_name",
                    "toa_instantaneous_radiative_forcing");
  acstoarfvar->add_att("long_name", "TOArad forcing av.");
  acstoarfvar->add_att("coordinates", "xlon xlat");
  acstoarfvar->add_att("units", "W m-2");
  acstsrrfvar = f->add_var("acstsrrf", ncFloat, tt, iy, jx);
  acstsrrfvar->add_att("standard_name",
                    "surface_instantaneous_radiative_forcing");
  acstsrrfvar->add_att("long_name", "SRFrad forcing av.");
  acstsrrfvar->add_att("coordinates", "xlon xlat");
  acstsrrfvar->add_att("units", "W m-2");
  if (ctl->doit)
  {
    gradsvar gv;
    gv.set("psa","psa","Surface pressure (hPa)",0,true);
    ctl->add_var(gv);
    gv.set("acstoarf","acstoarf","TOA rad forcing av. (W m-2)",0,true);
    ctl->add_var(gv);
    gv.set("acstsrrf","acstsrrf","SRF rad forcing av. (W m-2)",0,true);
    ctl->add_var(gv);
  }
}

void rcmNcChe::put_rec(chedata &c)
{
  double xtime = reference_time + tcount*c.dt;
  timevar->put_rec(&xtime, rcount);
  psvar->put_rec(c.psa, rcount);
  trac3Dvar->put_rec(c.trac3D, rcount);
  aext8var->put_rec(c.aext8, rcount);
  assa8var->put_rec(c.assa8, rcount);
  agfu8var->put_rec(c.agfu8, rcount);
  float *colb;
  float *wdlsc;
  float *wdcvc;
  float *sdrdp;
  float *xgasc;
  float *xaquc;
  float *emiss;
  float *base;
  for (int i = 0; i < c.ntr; i++)
  {
    base = c.trac2D+(i*(7*c.size2D));

    colb = base;
    colbvar->set_cur(rcount , i, 0, 0);
    colbvar->put(colb, 1, 1, nx, ny);
    wdlsc = colb+c.size2D;
    wdlscvar->set_cur(rcount , i, 0, 0);
    wdlscvar->put(wdlsc, 1, 1, nx, ny);
    wdcvc = wdlsc+c.size2D;
    wdcvcvar->set_cur(rcount , i, 0, 0);
    wdcvcvar->put(wdcvc, 1, 1, nx, ny);
    sdrdp = wdcvc+c.size2D;
    sdrdpvar->set_cur(rcount , i, 0, 0);
    sdrdpvar->put(sdrdp, 1, 1, nx, ny);
    xgasc = sdrdp+c.size2D;
    xgascvar->set_cur(rcount , i, 0, 0);
    xgascvar->put(xgasc, 1, 1, nx, ny);
    xaquc = xgasc+c.size2D;
    xaqucvar->set_cur(rcount , i, 0, 0);
    xaqucvar->put(xaquc, 1, 1, nx, ny);
    emiss = xaquc+c.size2D;
    emissvar->set_cur(rcount , i, 0, 0);
    emissvar->put(emiss, 1, 1, nx, ny);
  }
  acstoarfvar->put_rec(c.acstoarf, rcount);
  acstsrrfvar->put_rec(c.acstsrrf, rcount);
  rcount ++;
  tcount ++;
  if (ctl->doit)
    ctl->add_time( );
  return;
}

rcmNcSub::rcmNcSub(regcmout &fnc, header_data &h, subdom_data &s)
  : rcmNc(fnc, h, false)
{
  // Manage time setup
  rcmdate d(h.idate1);
  reference_time = d.unixtime( );
  if (ctl->doit)
    ctl->set_time(ctl->gradstime(d.basey,d.basem,d.based,d.baseh,(int) h.dtb));

  f->add_att("domain_name", s.name);
  f->add_att("input_dataset_resolution_in_minutes", s.ntypec_s);
  if (s.fudge_lnd_s)
    f->add_att("landuse_fudging", "Yes");
  else
    f->add_att("landuse_fudging", "No");
  if (s.fudge_tex_s)
    f->add_att("texture_fudging", "Yes");
  else
    f->add_att("texture_fudging", "No");
  float fillv = -1e+34;

  // Add subgrid dimensions
  NcDim *iys = f->add_dim("iys", s.nx);
  NcDim *jxs = f->add_dim("jxs", s.ny);

  // Add three more vertical dimensions for 10m, 2m hgts and soil layers
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
  NcVar *iyvar = f->add_var("iy", ncFloat, iy);
  iyvar->add_att("long_name", "y-coordinate in Cartesian system");
  iyvar->add_att("standard_name", "projection_y_coordinate");
  iyvar->add_att("axis", "Y");
  iyvar->add_att("units", "m");
  NcVar *jxvar = f->add_var("jx", ncFloat, jx);
  jxvar->add_att("long_name", "x-coordinate in Cartesian system");
  jxvar->add_att("standard_name", "projection_x_coordinate");
  jxvar->add_att("axis", "X");
  jxvar->add_att("units", "m");
  NcVar *xlatvar = f->add_var("xlatsub", ncFloat, iys, jxs);
  xlatvar->add_att("standard_name", "latitude");
  xlatvar->add_att("long_name", "Latitude");
  xlatvar->add_att("units", "degrees_north");
  NcVar *xlonvar = f->add_var("xlonsub", ncFloat, iys, jxs);
  xlonvar->add_att("standard_name", "longitude");
  xlonvar->add_att("long_name", "Longitude");
  xlonvar->add_att("units", "degrees_east");
  NcVar *htvar = f->add_var("topo", ncFloat, iys, jxs);
  htvar->add_att("standard_name", "surface_altitude");
  htvar->add_att("long_name", "Domain surface elevation");
  htvar->add_att("coordinates", "xlonsub xlatsub");
  htvar->add_att("units", "m");
  NcVar *htsdvar = f->add_var("htsdsub", ncFloat, iys, jxs);
  htsdvar->add_att("standard_name", "surface_altitude");
  htsdvar->add_att("long_name", "Domain elevation stantard deviation");
  htsdvar->add_att("coordinates", "xlonsub xlatsub");
  htsdvar->add_att("units", "m");
  htsdvar->add_att("cell_method", "area: standard_deviation");
  NcVar *landusevar = f->add_var("landusesub", ncFloat, iys, jxs);
  landusevar->add_att("long_name", "Landuse category as defined in BATS");
  landusevar->add_att("standard_name", "soil_type");
  landusevar->add_att("units", "1");
  landusevar->add_att("coordinates", "xlonsub xlatsub");
  landusevar->add_att("legend",
         "1  => Crop/mixed farming\n2  => Short grass\n"
         "3  => Evergreen needleleaf tree\n4  => Deciduous needleleaf tree\n"
         "5  => Deciduous broadleaf tree\n6  => Evergreen broadleaf tree\n"
         "7  => Tall grass\n8  => Desert\n9  => Tundra\n"
         "10 => Irrigated Crop\n11 => Semi-desert\n"
         "12 => Ice cap/glacier\n13 => Bog or marsh\n"
         "14 => Inland water\n15 => Ocean\n16 => Evergreen shrub\n"
         "17 => Deciduous shrub\n18 => Mixed Woodland\n"
         "19 => Forest/Field mosaic\n20 => Water and Land mixture");
  NcVar *xmapvar = f->add_var("xmapsub", ncFloat, iys, jxs);
  xmapvar->add_att("long_name", "Map factor in domain cross points");
  xmapvar->add_att("units", "1");
  xmapvar->add_att("coordinates", "xlonsub xlatsub");
  NcVar *coriolvar = f->add_var("coriolsub", ncFloat, iys, jxs);
  coriolvar->add_att("long_name", "Coriolis parameter");
  coriolvar->add_att("standard_name", "coriolis_parameter");
  coriolvar->add_att("units", "s-1");
  coriolvar->add_att("coordinates", "xlonsub xlatsub");
  NcVar *maskvar = f->add_var("masksub", ncFloat, iys, jxs);
  maskvar->add_att("long_name", "Land Sea mask");
  maskvar->add_att("standard_name", "land_binary_mask");
  maskvar->add_att("units", "1");
  maskvar->add_att("coordinates", "xlonsub xlatsub");

  xlatvar->put(s.xlat, s.nx, s.ny);
  xlonvar->put(s.xlon, s.nx, s.ny);
  htvar->put(s.ht, s.nx, s.ny);
  htsdvar->put(s.htsd, s.nx, s.ny);
  landusevar->put(s.landuse, s.nx, s.ny);
  xmapvar->put(s.xmap, s.nx, s.ny);
  coriolvar->put(s.coriol, s.nx, s.ny);
  maskvar->put(s.mask, s.nx, s.ny);
  float *tmp = new float[s.nx];
  float incr = s.ds/((float) s.nsg);
  tmp[0] = -(((float) s.nx-1)/2.0f) * incr;
  for (unsigned int i = 1; i < s.nx; i ++)
    tmp[i] = tmp[i-1] + incr;
  iyvar->put(tmp, s.nx);
  delete [] tmp;
  tmp = new float[s.ny];
  tmp[0] = -(((float) s.ny-1)/2.0f) * incr;
  for (unsigned int i = 1; i < s.ny; i ++)
    tmp[i] = tmp[i-1] + incr;
  jxvar->put(tmp, s.ny);
  delete [] tmp;

  f->sync();

  u10mvar = f->add_var("u10m", ncFloat, tt, m10, iys, jxs);
  u10mvar->add_att("standard_name", "eastward_wind");
  u10mvar->add_att("long_name", "10 meters U component (westerly) of wind");
  u10mvar->add_att("coordinates", "xlonsub xlatsub");
  u10mvar->add_att("units", "m s-1");
  v10mvar = f->add_var("v10m", ncFloat, tt, m10, iys, jxs);
  v10mvar->add_att("standard_name", "northward_wind");
  v10mvar->add_att("long_name", "10 meters V component (southerly) of wind");
  v10mvar->add_att("coordinates", "xlonsub xlatsub");
  v10mvar->add_att("units", "m s-1");
  uvdragvar = f->add_var("uvdrag", ncFloat, tt, iys, jxs);
  uvdragvar->add_att("standard_name", "surface_drag_coefficient_in_air");
  uvdragvar->add_att("long_name", "Surface drag stress");
  uvdragvar->add_att("coordinates", "xlonsub xlatsub");
  uvdragvar->add_att("units", "1");
  tgvar = f->add_var("tg", ncFloat, tt, iys, jxs);
  tgvar->add_att("standard_name", "surface_temperature");
  tgvar->add_att("long_name", "Ground temperature");
  tgvar->add_att("coordinates", "xlonsub xlatsub");
  tgvar->add_att("units", "K");
  tlefvar = f->add_var("tlef", ncFloat, tt, iys, jxs);
  tlefvar->add_att("standard_name", "canopy_temperature");
  tlefvar->add_att("long_name", "Foliage temperature");
  tlefvar->add_att("coordinates", "xlonsub xlatsub");
  tlefvar->add_att("_FillValue", fillv);
  tlefvar->add_att("units", "K");
  t2mvar = f->add_var("t2m", ncFloat, tt, m2, iys, jxs);
  t2mvar->add_att("standard_name", "air_temperature");
  t2mvar->add_att("long_name", "2 meters temperature");
  t2mvar->add_att("coordinates", "xlonsub xlatsub");
  t2mvar->add_att("units", "K");
  q2mvar = f->add_var("q2m", ncFloat, tt, m2, iys, jxs);
  q2mvar->add_att("standard_name", "humidity_mixing_ratio");
  q2mvar->add_att("long_name", "2 meters vapour mixing ratio");
  q2mvar->add_att("coordinates", "xlonsub xlatsub");
  q2mvar->add_att("units", "kg kg-1");
  r2mvar = f->add_var("r2m", ncFloat, tt, m2, iys, jxs);
  r2mvar->add_att("standard_name", "relative_humidity");
  r2mvar->add_att("long_name", "2 meters relative humidity");
  r2mvar->add_att("coordinates", "xlonsub xlatsub");
  r2mvar->add_att("units", "1");
  smwvar = f->add_var("smw", ncFloat, tt, soil, iys, jxs);
  smwvar->add_att("standard_name", "soil_moisture_content");
  smwvar->add_att("long_name", "Moisture content");
  smwvar->add_att("coordinates", "xlonsub xlatsub");
  smwvar->add_att("_FillValue", fillv);
  smwvar->add_att("units", "kg m-2");
  tprvar = f->add_var("tpr", ncFloat, tt, iys, jxs);
  tprvar->add_att("standard_name", "precipitation_amount");
  tprvar->add_att("long_name", "Total precipitation");
  tprvar->add_att("coordinates", "xlonsub xlatsub");
  tprvar->add_att("units", "kg m-2");
  evpvar = f->add_var("evp", ncFloat, tt, iys, jxs);
  evpvar->add_att("standard_name", "water_evaporation_amount");
  evpvar->add_att("long_name", "Total evapotranspiration");
  evpvar->add_att("coordinates", "xlonsub xlatsub");
  evpvar->add_att("units", "kg m-2");
  runoffvar = f->add_var("runoff", ncFloat, tt, iys, jxs);
  runoffvar->add_att("standard_name", "surface_runoff_flux");
  runoffvar->add_att("long_name", "surface runoff");
  runoffvar->add_att("coordinates", "xlonsub xlatsub");
  runoffvar->add_att("_FillValue", fillv);
  runoffvar->add_att("units", "kg m-2 day-1");
  scvvar = f->add_var("scv", ncFloat, tt, iys, jxs);
  scvvar->add_att("standard_name", "snowfall_amount");
  scvvar->add_att("long_name", "Snow precipitation");
  scvvar->add_att("coordinates", "xlonsub xlatsub");
  scvvar->add_att("_FillValue", fillv);
  scvvar->add_att("units", "kg m-2");
  senavar = f->add_var("sena", ncFloat, tt, iys, jxs);
  senavar->add_att("standard_name", "surface_downward_sensible_heat_flux");
  senavar->add_att("long_name", "Sensible heat flux");
  senavar->add_att("coordinates", "xlonsub xlatsub");
  senavar->add_att("units", "W m-2");
  prcvvar = f->add_var("prcv", ncFloat, tt, iys, jxs);
  prcvvar->add_att("standard_name", "convective_rainfall_flux");
  prcvvar->add_att("long_name", "Convective precipitation");
  prcvvar->add_att("coordinates", "xlonsub xlatsub");
  prcvvar->add_att("units", "kg m-2 day-1");
  psbvar = f->add_var("psa", ncFloat, tt, iys, jxs);
  psbvar->add_att("standard_name", "air_pressure");
  psbvar->add_att("long_name", "Surface pressure");
  psbvar->add_att("coordinates", "xlonsub xlatsub");
  psbvar->add_att("units", "hPa");
  if (ctl->doit)
  {
    gradsvar gv;
    ctl->addentry("vectorpairs u10m,v10m");
  }
}

void rcmNcSub::put_rec(subdata &s, t_srf_deriv &d)
{
  double xtime = reference_time + tcount*s.dt;
  timevar->put_rec(&xtime, rcount);
  u10mvar->put_rec(s.u10m, rcount);
  v10mvar->put_rec(s.v10m, rcount);
  uvdragvar->put_rec(s.uvdrag, rcount);
  tgvar->put_rec(s.tg, rcount);
  tlefvar->put_rec(s.tlef, rcount);
  t2mvar->put_rec(s.t2m, rcount);
  q2mvar->put_rec(s.q2m, rcount);
  r2mvar->put_rec(d.r2, rcount);
  smwvar->put_rec(s.smw, rcount);
  tprvar->put_rec(s.tpr, rcount);
  evpvar->put_rec(s.evp, rcount);
  runoffvar->put_rec(s.runoff, rcount);
  scvvar->put_rec(s.scv, rcount);
  senavar->put_rec(s.sena, rcount);
  prcvvar->put_rec(s.prcv, rcount);
  psbvar->put_rec(s.psb, rcount);
  rcount ++;
  tcount ++;
  if (ctl->doit)
    ctl->add_time( );
  return;
}

domNc::domNc(regcmout &fnc)
{
  f = new NcFile(fnc.fname.c_str(), NcFile::Replace);
  f->add_att("title", "ICTP Regional Climatic model V4 domain");
  f->add_att("institution", "ICTP");
  f->add_att("source", "RegCM Model simulation DOMAIN Input");
  f->add_att("Conventions", "CF-1.4");
  char buffer[256];
  time_t xtime = time(&xtime);
  snprintf(buffer, 256, "%s", ctime(&xtime));
  char *p = strchr(buffer, '\n');
  *p = 0;
  strncat(buffer, ": Created from DOMAIN input", 256);
  f->add_att("history", buffer);
  f->add_att("references", "http://users.ictp.it/RegCNET");
  snprintf(buffer, 256, "Experiment Name is : %s", fnc.experiment.c_str());
  f->add_att("comment", buffer);
  ctl = &(fnc.ctl);
  if (ctl->doit)
    ctl->head("ICTP Regional Climatic model V4 domain", -1e+34);
}

domNc::~domNc()
{
  f->close( );
  delete f;
}

void domNc::write(domain_data &d)
{
  f->add_att("projection", d.proj);
  f->add_att("grid_size_in_meters", d.ds);
  f->add_att("latitude_of_projection_origin", d.clat);
  f->add_att("longitude_of_projection_origin", d.clon);
  if (strcmp(d.proj, "POLSTR") == 0)
  {
    f->add_att("latitude_of_projection_pole", d.xplat);
    f->add_att("longitude_of_projection_pole", d.xplon);
  }
  if (strcmp(d.proj, "LAMCON") == 0)
  {
    double stp[2];
    stp[0] = d.trlat1;
    stp[1] = d.trlat2;
    f->add_att("standard_parallel", 2, stp);
  }
  f->add_att("domain_name", d.name);
  f->add_att("input_dataset_resolution_in_minutes", d.ntypec);
  if (d.anal)
    f->add_att("data_interpolation", "Cressman type objective analysis");
  else
    f->add_att("data_interpolation", "Overlapping parabolic 16 points");
  if (d.smthbdy) 
    f->add_att("boundary_smoothing", "No");
  else
    f->add_att("boundary_smoothing", "Yes");
  if (d.lakadj)
    f->add_att("great_lakes_adjustment", "Yes");
  else
    f->add_att("great_lakes_adjustment", "No");
  if (d.fudge_lnd)
    f->add_att("landuse_fudging", "Yes");
  else
    f->add_att("landuse_fudging", "No");
  if (d.fudge_tex)
    f->add_att("texture_fudging", "Yes");
  else
    f->add_att("texture_fudging", "No");
  f->add_att("number_of_textures", d.ntex);
  f->add_att("minimum_h2o_pct_for_water", d.h2opct);

  NcDim *iy = f->add_dim("iy", d.nx);
  NcDim *jx = f->add_dim("jx", d.ny);

  // Add time independent variables
  NcVar *iyvar = f->add_var("iy", ncFloat, iy);
  iyvar->add_att("long_name", "y-coordinate in Cartesian system");
  iyvar->add_att("standard_name", "projection_y_coordinate");
  iyvar->add_att("axis", "Y");
  iyvar->add_att("units", "m");
  NcVar *jxvar = f->add_var("jx", ncFloat, jx);
  jxvar->add_att("long_name", "x-coordinate in Cartesian system");
  jxvar->add_att("standard_name", "projection_x_coordinate");
  jxvar->add_att("axis", "X");
  jxvar->add_att("units", "m");
  NcVar *xlatvar = f->add_var("xlat", ncFloat, iy, jx);
  xlatvar->add_att("standard_name", "latitude");
  xlatvar->add_att("long_name", "Latitude");
  xlatvar->add_att("units", "degrees_north");
  NcVar *xlonvar = f->add_var("xlon", ncFloat, iy, jx);
  xlonvar->add_att("standard_name", "longitude");
  xlonvar->add_att("long_name", "Longitude");
  xlonvar->add_att("units", "degrees_east");
  NcVar *dlatvar = f->add_var("dlat", ncFloat, iy, jx);
  dlatvar->add_att("standard_name", "latitude");
  dlatvar->add_att("long_name", "Dot Point Latitude");
  dlatvar->add_att("units", "degrees_north");
  NcVar *dlonvar = f->add_var("dlon", ncFloat, iy, jx);
  dlonvar->add_att("standard_name", "longitude");
  dlonvar->add_att("long_name", "Dot Point Longitude");
  dlonvar->add_att("units", "degrees_east");
  NcVar *htvar = f->add_var("topo", ncFloat, iy, jx);
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
  NcVar *landusevar = f->add_var("landuse", ncFloat, iy, jx);
  landusevar->add_att("long_name", "Landuse category as defined in BATS");
  landusevar->add_att("standard_name", "soil_type");
  landusevar->add_att("units", "1");
  landusevar->add_att("coordinates", "xlon xlat");
  landusevar->add_att("legend",
         "1  => Crop/mixed farming\n2  => Short grass\n"
         "3  => Evergreen needleleaf tree\n4  => Deciduous needleleaf tree\n"
         "5  => Deciduous broadleaf tree\n6  => Evergreen broadleaf tree\n"
         "7  => Tall grass\n8  => Desert\n9  => Tundra\n"
         "10 => Irrigated Crop\n11 => Semi-desert\n"
         "12 => Ice cap/glacier\n13 => Bog or marsh\n"
         "14 => Inland water\n15 => Ocean\n16 => Evergreen shrub\n"
         "17 => Deciduous shrub\n18 => Mixed Woodland\n"
         "19 => Forest/Field mosaic\n20 => Water and Land mixture");
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
  NcVar *snowvar = f->add_var("snowam", ncFloat, iy, jx);
  snowvar->add_att("long_name", "Snow amount");
  snowvar->add_att("standard_name", "snowfall_amount");
  snowvar->add_att("units", "kg m-2");
  snowvar->add_att("coordinates", "xlon xlat");
  NcVar *maskvar = f->add_var("mask", ncFloat, iy, jx);
  maskvar->add_att("long_name", "Land Sea mask");
  maskvar->add_att("standard_name", "land_binary_mask");
  maskvar->add_att("units", "1");
  maskvar->add_att("coordinates", "xlon xlat");
  if (ctl->doit)
  {
    ctl->set_grid(d);
    gradsvar gv;
    gv.set("xlat","xlat","Latitude (degrees_north)",0,false);
    ctl->add_var(gv);
    gv.set("xlon","xlon","Longitude (degrees_east)",0,false);
    ctl->add_var(gv);
    gv.set("dlat","dlat","Dot point Latitude (degrees_north)",0,false);
    ctl->add_var(gv);
    gv.set("dlon","dlon","Dot point Longitude (degrees_east)",0,false);
    ctl->add_var(gv);
    gv.set("topo","ht","Surface elevation (m)",0,false);
    ctl->add_var(gv);
    gv.set("htsd","htsd","Surface elevation standard dev. (m)",0,false);
    ctl->add_var(gv);
    gv.set("landuse","landuse","Landuse BATS category. (1)",0,false);
    ctl->add_var(gv);
    gv.set("xmap","xmap","Map factor cross points. (1)",0,false);
    ctl->add_var(gv);
    gv.set("dmap","dmap","Map factor dot points. (1)",0,false);
    ctl->add_var(gv);
    gv.set("coriol","coriol","Coriolis parameter. (s-1)",0,false);
    ctl->add_var(gv);
    gv.set("snowam","snowam","Snow amount. (kg m-2)",0,false);
    ctl->add_var(gv);
    gv.set("mask","mask","Land/Sea mask. (1)",0,false);
    ctl->add_var(gv);
  }

  float *tmp = new float[d.nx];
  tmp[0] = -(((float) d.nx-1)/2.0f) * d.ds;
  for (unsigned int i = 1; i < d.nx; i ++)
    tmp[i] = tmp[i-1] + d.ds;
  iyvar->put(tmp, d.nx);
  delete [] tmp;
  tmp = new float[d.ny];
  tmp[0] = -(((float) d.ny-1)/2.0f) * d.ds;
  for (unsigned int i = 1; i < d.ny; i ++)
    tmp[i] = tmp[i-1] + d.ds;
  jxvar->put(tmp, d.ny);
  delete [] tmp;

  xlatvar->put(d.xlat, d.nx, d.ny);
  xlonvar->put(d.xlon, d.nx, d.ny);
  dlatvar->put(d.dlat, d.nx, d.ny);
  dlonvar->put(d.dlon, d.nx, d.ny);
  htvar->put(d.ht, d.nx, d.ny);
  htsdvar->put(d.htsd, d.nx, d.ny);
  landusevar->put(d.landuse, d.nx, d.ny);
  xmapvar->put(d.xmap, d.nx, d.ny);
  dmapvar->put(d.dmap, d.nx, d.ny);
  coriolvar->put(d.coriol, d.nx, d.ny);
  snowvar->put(d.snowam, d.nx, d.ny);
  maskvar->put(d.mask, d.nx, d.ny);

  if (ctl->doit)
    ctl->add_time( );
  return;
}

bcNc::bcNc(regcmout &fnc, domain_data &d)
{
  // Check if var is on wanted list
  // Special all key
  for (int i = 0; i < 10; i ++)
    varmask[i] = false;
  if (fnc.vl.isthere("all"))
  {
    for (int i = 0; i < 10; i ++)
      varmask[i] = true;
  }
  else
  {
    if (fnc.vl.isthere("u")) varmask[0] = true;
    if (fnc.vl.isthere("v")) varmask[1] = true;
    if (fnc.vl.isthere("t")) varmask[2] = true;
    if (fnc.vl.isthere("qv")) varmask[3] = true;
    if (fnc.vl.isthere("ts")) varmask[4] = true;
    if (fnc.vl.isthere("so4")) varmask[5] = true;
    if (fnc.vl.isthere("sm")) varmask[6] = true;
    if (fnc.vl.isthere("it")) varmask[7] = true;
    if (fnc.vl.isthere("st")) varmask[8] = true;
    if (fnc.vl.isthere("sn")) varmask[9] = true;
  }

  f = new NcFile(fnc.fname.c_str(), NcFile::Replace);
  f->add_att("title", "ICTP Regional Climatic model V4 BC input");
  f->add_att("institution", "ICTP");
  f->add_att("source", "RegCM Model simulation DOMAIN Boundary Conditions");
  f->add_att("Conventions", "CF-1.4");
  char buffer[256];
  time_t xtime = time(&xtime);
  snprintf(buffer, 256, "%s", ctime(&xtime));
  char *p = strchr(buffer, '\n');
  *p = 0;
  strncat(buffer, ": Created from model boundary conditions", 256);
  f->add_att("history", buffer);
  f->add_att("references", "http://users.ictp.it/RegCNET");
  snprintf(buffer, 256, "Experiment Name is : %s", fnc.experiment.c_str());
  f->add_att("comment", buffer);
  f->add_att("projection", d.proj);
  f->add_att("grid_size_in_meters", d.ds);
  f->add_att("latitude_of_projection_origin", d.clat);
  f->add_att("longitude_of_projection_origin", d.clon);
  if (strcmp(d.proj, "POLSTR") == 0)
  {
    f->add_att("latitude_of_projection_pole", d.xplat);
    f->add_att("longitude_of_projection_pole", d.xplon);
  }
  if (strcmp(d.proj, "LAMCON") == 0)
  {
    double stp[2];
    stp[0] = d.trlat1;
    stp[1] = d.trlat2;
    f->add_att("standard_parallel", 2, stp);
  }

  iy = f->add_dim("iy", d.nx);
  jx = f->add_dim("jx", d.ny);
  kz = f->add_dim("kz", d.nz);
  tt = f->add_dim("time");

  ctl = &(fnc.ctl);
  if (ctl->doit)
  {
    ctl->head("ICTP Regional Climatic model V4 domain", -1e+34);
    ctl->set_grid(d);
  }
  // Add time variable
  timevar = f->add_var("time", ncDouble, tt);
  timevar->add_att("long_name", "time");
  timevar->add_att("standard_name", "time");
  timevar->add_att("calendar", "standard");
  timevar->add_att("units", "hours since 1970-01-01 00:00:00 UTC");

  // Add time independent variables
  NcVar *iyvar = f->add_var("iy", ncFloat, iy);
  iyvar->add_att("long_name", "y-coordinate in Cartesian system");
  iyvar->add_att("standard_name", "projection_y_coordinate");
  iyvar->add_att("axis", "Y");
  iyvar->add_att("units", "m");
  NcVar *jxvar = f->add_var("jx", ncFloat, jx);
  jxvar->add_att("long_name", "x-coordinate in Cartesian system");
  jxvar->add_att("standard_name", "projection_x_coordinate");
  jxvar->add_att("axis", "X");
  jxvar->add_att("units", "m");
  NcVar *sigfvar = f->add_var("level", ncFloat, kz);
  sigfvar->add_att("standard_name", "atmosphere_sigma_coordinate");
  sigfvar->add_att("long_name", "Sigma at model layer midpoints");
  sigfvar->add_att("positive", "down");
  sigfvar->add_att("units", "1");
  sigfvar->add_att("axis", "Z");
  sigfvar->add_att("formula_terms", "sigma: level ps: psa ptop: ptop");
  NcVar *ptopvar = f->add_var("ptop", ncFloat);
  ptopvar->add_att("long_name", "Pressure at model top");
  ptopvar->add_att("standard_name", "air_pressure");
  ptopvar->add_att("units", "hPa");
  NcVar *xlatvar = f->add_var("xlat", ncFloat, iy, jx);
  xlatvar->add_att("standard_name", "latitude");
  xlatvar->add_att("long_name", "Latitude");
  xlatvar->add_att("units", "degrees_north");
  NcVar *xlonvar = f->add_var("xlon", ncFloat, iy, jx);
  xlonvar->add_att("standard_name", "longitude");
  xlonvar->add_att("long_name", "Longitude");
  xlonvar->add_att("units", "degrees_east");
  NcVar *dlatvar = f->add_var("dlat", ncFloat, iy, jx);
  dlatvar->add_att("standard_name", "latitude");
  dlatvar->add_att("long_name", "Dot Point Latitude");
  dlatvar->add_att("units", "degrees_north");
  NcVar *dlonvar = f->add_var("dlon", ncFloat, iy, jx);
  dlonvar->add_att("standard_name", "longitude");
  dlonvar->add_att("long_name", "Dot Point Longitude");
  dlonvar->add_att("units", "degrees_east");
  NcVar *htvar = f->add_var("topo", ncFloat, iy, jx);
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
  NcVar *landusevar = f->add_var("landuse", ncFloat, iy, jx);
  landusevar->add_att("long_name", "Landuse category as defined in BATS");
  landusevar->add_att("standard_name", "soil_type");
  landusevar->add_att("units", "1");
  landusevar->add_att("coordinates", "xlon xlat");
  landusevar->add_att("legend",
         "1  => Crop/mixed farming\n2  => Short grass\n"
         "3  => Evergreen needleleaf tree\n4  => Deciduous needleleaf tree\n"
         "5  => Deciduous broadleaf tree\n6  => Evergreen broadleaf tree\n"
         "7  => Tall grass\n8  => Desert\n9  => Tundra\n"
         "10 => Irrigated Crop\n11 => Semi-desert\n"
         "12 => Ice cap/glacier\n13 => Bog or marsh\n"
         "14 => Inland water\n15 => Ocean\n16 => Evergreen shrub\n"
         "17 => Deciduous shrub\n18 => Mixed Woodland\n"
         "19 => Forest/Field mosaic\n20 => Water and Land mixture");
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
  NcVar *snowvar = f->add_var("snowam", ncFloat, iy, jx);
  snowvar->add_att("long_name", "Snow amount");
  snowvar->add_att("standard_name", "snowfall_amount");
  snowvar->add_att("units", "kg m-2");
  snowvar->add_att("coordinates", "xlon xlat");
  NcVar *maskvar = f->add_var("mask", ncFloat, iy, jx);
  maskvar->add_att("long_name", "Land Sea mask");
  maskvar->add_att("standard_name", "land_binary_mask");
  maskvar->add_att("units", "1");
  maskvar->add_att("coordinates", "xlon xlat");
  if (ctl->doit)
  {
    gradsvar gv;
    gv.set("xlat","xlat","Latitude (degrees_north)",0,false);
    ctl->add_var(gv);
    gv.set("xlon","xlon","Longitude (degrees_east)",0,false);
    ctl->add_var(gv);
    gv.set("dlat","dlat","Latitude dot points (degrees_north)",0,false);
    ctl->add_var(gv);
    gv.set("dlon","dlon","Longitude dot points (degrees_east)",0,false);
    ctl->add_var(gv);
    gv.set("topo","ht","Surface elevation (m)",0,false);
    ctl->add_var(gv);
    gv.set("htsd","htsd","Surface elevation std. dev. (m)",0,false);
    ctl->add_var(gv);
    gv.set("landuse","landuse","Landuse BATS categories. (1)",0,false);
    ctl->add_var(gv);
    gv.set("xmap","xmap","Map factor. (1)",0,false);
    ctl->add_var(gv);
    gv.set("dmap","dmap","Map factor dot points. (1)",0,false);
    ctl->add_var(gv);
    gv.set("coriol","coriol","Coriolis parameter. (1)",0,false);
    ctl->add_var(gv);
    gv.set("snowam","snowam","Snow amount. (kg m-2)",0,false);
    ctl->add_var(gv);
    gv.set("mask","mask","Land/Sea mask. (1)",0,false);
    ctl->add_var(gv);
  }

  float *tmp = new float[d.nz];
  for (unsigned int i = 0; i < d.nz; i ++)
    tmp[i] = d.hsigm[d.nz-1-i];
  sigfvar->put(tmp, d.nz);
  if (ctl->doit)
  {
    for (unsigned int i = 0; i < d.nz; i ++)
      tmp[i] *= 1000.0;
    ctl->set_levs(tmp, d.nz);
  }
  delete [] tmp;
  tmp = new float[d.nx];
  tmp[0] = -(((float) d.nx-1)/2.0f) * d.ds;
  for (unsigned int i = 1; i < d.nx; i ++)
    tmp[i] = tmp[i-1] + d.ds;
  iyvar->put(tmp, d.nx);
  delete [] tmp;
  tmp = new float[d.ny];
  tmp[0] = -(((float) d.ny-1)/2.0f) * d.ds;
  for (unsigned int i = 1; i < d.ny; i ++)
    tmp[i] = tmp[i-1] + d.ds;
  jxvar->put(tmp, d.ny);
  delete [] tmp;
  ptopvar->put(&d.ptop, 1);
  xlatvar->put(d.xlat, d.nx, d.ny);
  xlonvar->put(d.xlon, d.nx, d.ny);
  dlatvar->put(d.dlat, d.nx, d.ny);
  dlonvar->put(d.dlon, d.nx, d.ny);
  htvar->put(d.ht, d.nx, d.ny);
  htsdvar->put(d.htsd, d.nx, d.ny);
  landusevar->put(d.landuse, d.nx, d.ny);
  xmapvar->put(d.xmap, d.nx, d.ny);
  dmapvar->put(d.dmap, d.nx, d.ny);
  coriolvar->put(d.coriol, d.nx, d.ny);
  snowvar->put(d.snowam, d.nx, d.ny);
  maskvar->put(d.mask, d.nx, d.ny);

  f->sync();

  notinit = true;
}

bcNc::~bcNc()
{
  f->close();
  delete f;
}

void bcNc::put_rec(bcdata &b)
{
  double xtime;
  if (notinit)
  {
    f->add_att("input_sst_type", b.ssttyp);
    f->add_att("input_analysis_type", b.dattyp);
    // Manage time setup
    rcmdate d(b.rdate);
    reference_time = d.unixtime( );
    rcount = 0;
    tcount = 0;

    if (ctl->doit)
      ctl->set_time(ctl->gradstime(d.basey,d.basem,d.based,d.baseh,(int) b.dt));

    // This will allways be here for vertical level calculation
    psvar = f->add_var("psa", ncFloat, tt, iy, jx);
    psvar->add_att("standard_name", "surface_air_pressure");
    psvar->add_att("long_name", "Surface pressure");
    psvar->add_att("coordinates", "xlon xlat");
    psvar->add_att("units", "hPa");

    if (varmask[0] || varmask[1])
    {
      uvar = f->add_var("u", ncFloat, tt, kz, iy, jx);
      uvar->add_att("standard_name", "eastward_wind");
      uvar->add_att("long_name", "U component (westerly) of wind");
      uvar->add_att("coordinates", "xlon xlat");
      uvar->add_att("units", "m s-1");
      vvar = f->add_var("v", ncFloat, tt, kz, iy, jx);
      vvar->add_att("standard_name", "northward_wind");
      vvar->add_att("long_name", "V component (southerly) of wind");
      vvar->add_att("coordinates", "xlon xlat");
      vvar->add_att("units", "m s-1");
    }
    if (varmask[2])
    {
      tvar = f->add_var("t", ncFloat, tt, kz, iy, jx);
      tvar->add_att("standard_name", "air_temperature");
      tvar->add_att("long_name", "Temperature");
      tvar->add_att("coordinates", "xlon xlat");
      tvar->add_att("units", "K");
    }
    if (varmask[3])
    {
      qvvar = f->add_var("qv", ncFloat, tt, kz, iy, jx);
      qvvar->add_att("standard_name", "humidity_mixing_ratio");
      qvvar->add_att("long_name", "Water vapor mixing ratio");
      qvvar->add_att("coordinates", "xlon xlat");
      qvvar->add_att("units", "kg kg-1");
    }
    if (varmask[4])
    {
      tsvar = f->add_var("ts", ncFloat, tt, iy, jx);
      tsvar->add_att("standard_name", "soil_temperature");
      tsvar->add_att("long_name", "Temperature");
      tsvar->add_att("coordinates", "xlon xlat");
      tsvar->add_att("units", "K");
    }
    if (ctl->doit)
    {
      gradsvar gv;
      gv.set("psa","psa","Surface pressure (hPa)",0,true);
      ctl->add_var(gv);
      if (varmask[0] || varmask[1])
      {
        ctl->addentry("vectorpais u,v");
        gv.set("u","u","U component (westerly) of wind (m s-1)",b.nz,true);
        ctl->add_var(gv);
        gv.set("v","v","V component (southerly) of wind (m s-1)",b.nz,true);
        ctl->add_var(gv);
      }
      if (varmask[2])
      {
        gv.set("t","t","Temperature (K)",b.nz,true);
        ctl->add_var(gv);
      }
      if (varmask[3])
      {
        gv.set("qv","qv","Humidity mixing ratio (kg kg-1)",b.nz,true);
        ctl->add_var(gv);
      }
      if (varmask[4])
      {
        gv.set("ts","ts","Soil temperature (K)",0,true);
        ctl->add_var(gv);
      }
    }

    if (b.ehso4 && varmask[5])
    {
      so4var = f->add_var("so4", ncFloat, tt, iy, jx);
      so4var->add_att("standard_name", "atmosphere_sulfate_content");
      so4var->add_att("long_name", "Sulfate");
      so4var->add_att("coordinates", "xlon xlat");
      so4var->add_att("units", "kg m-2");
      if (ctl->doit)
      {
        gradsvar gv;
        gv.set("so4","so4","Sulfate content (kg m-1)",0,true);
        ctl->add_var(gv);
      }
    }

    if (b.usgs)
    {
      soil = f->add_dim("soil", 4);
      if (varmask[6])
      {
        smvar = f->add_var("sm", ncFloat, tt, soil, iy, jx);
        smvar->add_att("standard_name", "moisture_content_of_soil_layer");
        smvar->add_att("long_name", "Soil Moisture");
        smvar->add_att("coordinates", "xlon xlat");
        smvar->add_att("units", "kg m-2");
      }
      if (varmask[7])
      {
        itvar = f->add_var("it", ncFloat, tt, soil, iy, jx);
        itvar->add_att("standard_name", "ice_temperature_of_soil_layer");
        itvar->add_att("long_name", "Soil ice temperature");
        itvar->add_att("coordinates", "xlon xlat");
        itvar->add_att("units", "K");
      }
      if (varmask[8])
      {
        stvar = f->add_var("ts", ncFloat, tt, soil, iy, jx);
        stvar->add_att("standard_name", "soil_temperaturei_of_soil_layer");
        stvar->add_att("long_name", "Soil Temperature");
        stvar->add_att("coordinates", "xlon xlat");
        stvar->add_att("units", "K");
      }
      if (varmask[9])
      {
        snam = f->add_var("sa", ncFloat, tt, iy, jx);
        snam->add_att("standard_name", "thickness_of_snowfall_amount");
        snam->add_att("long_name", "Snow amount");
        snam->add_att("coordinates", "xlon xlat");
        snam->add_att("units", "m");
      }
      if (ctl->doit)
      {
        std::cout <<
          "Warning: Soil levels variables will not be included in CTL file." <<
          std::endl << "GrADS is unable to plot different verical levels."  <<
          std::endl;
      }
    }

    xtime = (double) reference_time;
    notinit = false;
  }
  else
    xtime = (double) reference_time + b.dt*tcount;
  timevar->put_rec(&xtime, rcount);
  for (unsigned int i = 0; i < b.size2D; i ++) b.px[i] *= 10.0;
  psvar->put_rec(b.px, rcount);
  if (varmask[0] || varmask[1])
  {
    uvar->put_rec(b.u, rcount);
    vvar->put_rec(b.v, rcount);
  }
  if (varmask[2]) tvar->put_rec(b.t, rcount);
  if (varmask[3]) qvvar->put_rec(b.q, rcount);
  if (varmask[4]) tsvar->put_rec(b.ts, rcount);
  if (b.ehso4 && varmask[5])
    so4var->put_rec(b.so4, rcount);
  if (b.usgs)
  {
    if (varmask[6]) smvar->put_rec(b.sm, rcount);
    if (varmask[7]) itvar->put_rec(b.icet, rcount);
    if (varmask[8]) stvar->put_rec(b.soilt, rcount);
    if (varmask[9]) snam->put_rec(b.snowd, rcount);
  }
  rcount ++;
  tcount ++;
  if (ctl->doit)
    ctl->add_time( );
  return;
}
