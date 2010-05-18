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

#include <iostream>

#include <netcdf.hh>

#include <cstring>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <libgen.h>
#include <unistd.h>
#include <getopt.h>

#include <rcmio.h>
#include <rcmNc.h>

using namespace rcm;

const static char version[] = SVN_REV;

void help(char *pname)
{
  std::cerr << std::endl
      << "                 RegCM V4 ICTP NetCDF Postprocessor." << std::endl
      << std::endl
      << "This simple program converts binary ICBC files from RegCM V4 "
      << "into NetCDF" << std::endl << "CF-1.4 convention compliant data files."
      << std::endl << std::endl << "I need four mandatory arguments:"
      << std::endl << std::endl
      << "    regcm.in       - path to regcm.in of RegCM model v4" << std::endl
      << "    DOMAIN.INFO    - path to DOMAIN.INFO of RegCM model v4"
      << std::endl
      << "    ICBCYYYYMMDDHH - path to ICBC of RegCM model v4" << std::endl
      << "    expname        - a (meaningful) name for this expertiment"
      << std::endl << std::endl << "Example:" << std::endl
      << std::endl << "     " << pname
      << " [options] DOMAIN.INFO ICBC1990060100 ACWA_reference"
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

  if (argc - optind != 4)
  {
    std::cerr << std::endl << "Howdy there, not enough arguments." << std::endl;
    help(argv[0]);
    return -1;
  }

  if (iseq == 1) ldirect = false;
  if (ilittle == 1) lbigend = false;

  try
  {
    char *regcmin = strdup(argv[optind++]);
    char *dominfo = strdup(argv[optind++]);
    char *arg2 = strdup(argv[optind++]);
    char *brg2 = strdup(arg2);
    char *icbcf = basename(arg2);
    char *icbcd = dirname(brg2);
    char *experiment = strdup(argv[optind++]);
 
    rcminp inpf(regcmin);
    domain_data d(inpf);
    rcmio rcmout(icbcd, lbigend, ldirect);
    rcmout.read_domain(dominfo, d);

    bcdata b(d, inpf);

    char fname[PATH_MAX];
    sprintf(fname, "ICBC_%s.nc", experiment);
    bcNc bcnc(fname, experiment, d);

    std::cout << "Processing ICBC";
    while ((rcmout.bc_read_tstep(b)) == 0)
    {
      std::cout << ".";
      std::cout.flush();
      bcnc.put_rec(b);
    }
    std::cout << " Done." << std::endl;

    free(dominfo);
    free(arg2);
    free(brg2);
    free(experiment);
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
