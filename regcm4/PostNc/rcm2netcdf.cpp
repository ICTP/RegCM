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

#include <netcdf.hh>
#include <cstring>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <rcmio.h>
#include <rcmNc.h>

using namespace rcm;

int main(int argc, char *argv[])
{
  if (argc != 3)
  {
    std::cerr << std::endl
        << "Howdy there, wrong number of arguments." << std::endl
        << std::endl << "I need two arguments:" << std::endl
        << "    outdir  - output directory of RegCM model v4" << std::endl
        << "    expname - a (meaningful) name for this expertiment" << std::endl
        << std::endl << "Example:" << std::endl << "     " << argv[0]
        << " /home/regcm/Run/output ACWA_reference" << std::endl << std::endl;
   return -1;
  }

  try
  {
    rcmio rcmout(argv[1], true, true);
    header_data outhead;
    rcmout.read_header(outhead);

    if (rcmout.has_atm)
    {
      atmodata a(outhead.nx, outhead.ny, outhead.nz, 
                 outhead.mdate0, outhead.dto);
      char fname[PATH_MAX];
      sprintf(fname, "ATM_%s_%d.nc", argv[2], outhead.mdate0);
      rcmNcAtmo atmnc(fname, argv[2], outhead);
      // Add Atmospheric variables
      while ((rcmout.atmo_read_tstep(a)) == 0)
        atmnc.put_rec(a);
    }

    if (rcmout.has_srf)
    {
      srfdata s(outhead.nx, outhead.ny, outhead.mdate0, outhead.dtb);
      char fname[PATH_MAX];
      sprintf(fname, "SRF_%s_%d.nc", argv[2], outhead.mdate0);
      rcmNcSrf srfnc(fname, argv[2], outhead);
      // Add Surface variables
      while ((rcmout.srf_read_tstep(s)) == 0)
        srfnc.put_rec(s);
    }

    if (rcmout.has_rad)
    {
      raddata r(outhead.nx, outhead.ny, outhead.nz,
                outhead.mdate0, outhead.dtr);
      char fname[PATH_MAX];
      sprintf(fname, "RAD_%s_%d.nc", argv[2], outhead.mdate0);
      rcmNcRad radnc(fname, argv[2], outhead);
      // Add Radiation variables
      while ((rcmout.rad_read_tstep(r)) == 0)
        radnc.put_rec(r);
    }

    if (rcmout.has_che)
    {
      chedata c(outhead.nx, outhead.ny, outhead.nz,
                outhead.mdate0, outhead.dtc);
      char fname[PATH_MAX];
      sprintf(fname, "CHE_%s_%d.nc", argv[2], outhead.mdate0);
      rcmNcChe chenc(fname, argv[2], outhead);
      // Add Chemical tracers variables
      while ((rcmout.che_read_tstep(c)) == 0)
        chenc.put_rec(c);
    }

    outhead.free_space( );

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

  std::cout << "Done" << std::endl;
  return 0;
}
