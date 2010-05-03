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

#include <iostream>

#include <netcdf.hh>

#include <cstring>
#include <cstdio>
#include <ctime>
#include <libgen.h>
#include <unistd.h>
#include <getopt.h>
#include <rcmio.h>

using namespace rcm;

const static char version[] = SVN_REV;

void help(char *pname)
{
  std::cerr << std::endl
      << "                 RegCM V4 ICTP NetCDF Postprocessor." << std::endl
      << std::endl
      << "This simple programs converts binary domain files from RegCM V4 "
      << "into NetCDF" << std::endl << "CF-1.4 convention compliant data files."
      << std::endl << std::endl << "I need two mandatory arguments:"
      << std::endl << std::endl
      << "    DOMAIN.INFO - path to DOMAIN.INFO of RegCM model v4" << std::endl
      << "    expname     - a (meaningful) name for this expertiment"
      << std::endl << std::endl << "Example:" << std::endl
      << std::endl << "     " << pname
      << " [options] /home/regcm/Input/DOMAIN.INFO ACWA_reference"
      << std::endl << std::endl
      << "where options can be in:" << std::endl << std::endl
  << "   --sequential              : Set I/O non direct (direct access default)"
      << std::endl
  << "   --little_endian           : Set I/O endianess to LITTLE (BIG default)"
      << std::endl
  << "   --help/-h                 : Print this help"
      << std::endl
  << "   --version/-V              : Print versioning information"
      << std::endl << std::endl;
   return;
}

int main(int argc, char *argv[])
{
  bool ldirect, lbigend;
  int iseq, ilittle;
  ldirect = true;
  lbigend = true;
  iseq = 0;
  ilittle = 0;

  while (1)
  {
    static struct option long_options[] = {
      { "sequential", no_argument, &iseq, 1},
      { "little_endian", no_argument, &ilittle, 1},
      { "help", no_argument, 0, 'h'},
      { "version", no_argument, 0, 'V'},
      { 0, 0, 0, 0 }
    };
    int optind, c = 0;
    c = getopt_long (argc, argv, "hV",
                     long_options, &optind);
    if (c == -1) break;
    switch (c)
    {
      case 0:
        if (long_options[optind].flag != 0) break;
      case 'h':
        help(argv[0]);
        return 0;
        break;
      case 'V':
        std::cerr << "This is " << basename(argv[0]) << " version " << version
                  << std::endl;
        return 0;
      case '?':
        break;
      default:
        std::cerr << "Unknown switch '" << (char) c << "' discarded."
                  << std::endl;
        break;
    }
  }

  if (argc - optind != 2)
  {
    std::cerr << std::endl << "Howdy there, not enough arguments." << std::endl;
    help(argv[0]);
    return -1;
  }

  if (iseq == 1) ldirect = false;
  if (ilittle == 1) lbigend = false;

  try
  {
    char *dominfo = strdup(argv[optind++]);
    char *experiment = strdup(argv[optind++]);
    
    domain_data d;
    rcmio rcmout(0, lbigend, ldirect);
    rcmout.read_domain(dominfo, d);

    char fname[PATH_MAX];
    sprintf(fname, "DOMAIN_%s.nc", experiment);

    NcFile *f = new NcFile(fname, NcFile::Replace);
    f->add_att("title", "ICTP Regional Climatic model V4 output");
    f->add_att("institution", "ICTP Trieste");
    f->add_att("source", "RegCM Model simulation DOMAIN Input");
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

    NcDim *iy = f->add_dim("iy", d.nx);
    NcDim *jx = f->add_dim("jx", d.ny);
    NcDim *kz = f->add_dim("kz", d.nz);

    // Add time independent variables and time itself
    NcVar *sigfvar = f->add_var("level", ncFloat, kz);
    sigfvar->add_att("standard_name", "atmosphere_sigma_coordinate");
    sigfvar->add_att("long_name", "Sigma at model layer midpoints");
    sigfvar->add_att("positive", "down");
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

    sigfvar->put(d.hsigm, d.nz);
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

    f->close( );
  }
  catch (const char *e)
  {
    std::cerr << "Error : " << e << std::endl;
    return -1;
  }
  catch (...)
  {
    return -1;
  }

  std::cout << "Successfully completed processing." << std::endl;
  return 0;
}
