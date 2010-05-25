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
#include <string>

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
#include <gradsctl.h>

using namespace rcm;

static void help(char *pname);
static const char version[] = SVN_REV;

void help(char *pname)
{
  std::cerr << std::endl
      << "                 RegCM V4 ICTP NetCDF Postprocessor." << std::endl
      << std::endl
      << "This simple program converts binary ICBC files from RegCM V4 "
      << "into NetCDF" << std::endl << "CF-1.4 convention compliant data files."
      << std::endl << std::endl << "I need ONE mandatory argument:"
      << std::endl << std::endl
      << "    regcm.in       - path to regcm.in of RegCM model v4" << std::endl
      << std::endl << "Example:" << std::endl << std::endl << "     " << pname
      << " [options] regcm.in"
      << std::endl << std::endl
      << "where options can be in:" << std::endl << std::endl
  << "   --sequential              : Set I/O sequential (direct access default)"
      << std::endl
  << "   --little_endian           : Set I/O endianess to LITTLE (BIG default)"
      << std::endl
  << "   --grads                   : Produce a CTL file for GrADS"
      << std::endl
  << "   --var/-v[all|name[,name]] : Include only some vars (all default)"
      << std::endl
  << "   --list/-l                 : Output list of names for -v option"
      << std::endl
  << "   --help/-h                 : Print this help"
      << std::endl
  << "   --version/-V              : Print versioning information"
      << std::endl << std::endl;
   return;
}

int main(int argc, char *argv[])
{
  bool ldirect, lbigend, lgrads;
  int iseq, ilittle, igrads;
  std::string vnames = "all";
  ldirect = true;
  lbigend = true;
  lgrads = false;
  iseq = 0;
  ilittle = 0;
  igrads = 0;

  char *pname = basename(argv[0]);
  while (1)
  {
    static struct option long_options[] = {
      { "sequential", no_argument, &iseq, 1},
      { "little_endian", no_argument, &ilittle, 1},
      { "grads", no_argument, &igrads, 1},
      { "help", no_argument, 0, 'h'},
      { "list", no_argument, 0, 'l'},
      { "var", required_argument, 0, 'v'},
      { "version", no_argument, 0, 'V'},
      { 0, 0, 0, 0 }
    };
    int optind, c = 0;
    c = getopt_long (argc, argv, "hVlv:",
                     long_options, &optind);
    if (c == -1) break;
    switch (c)
    {
      case 0:
        if (long_options[optind].flag != 0) break;
      case 'l':
        std::cout << std::endl
                  << "List of names for " << pname << " is:" << std::endl
                  << std::endl
                  << "  u  :  U component of wind (westerly)" << std::endl
                  << "  v  :  V component of wind (southerly)" << std::endl
                  << "  t  :  Temperature" << std::endl
                  << "  qv :  Water Vapor mixing ratio" << std::endl
                  << "  ts :  Soil temperature" << std::endl
                  << "  so4:  Sulfate (if avail)" << std::endl
                  << "  sm :  Soil moisture (if avail)" << std::endl
                  << "  it :  Ice temperature (if avail)" << std::endl
                  << "  st :  Soil temperature (if avail)" << std::endl
                  << "  sn :  Snow thickness (if avail)" << std::endl
                  << std::endl;
        return 0;
      case 'v':
        vnames = optarg;
        break;
      case 'h':
        help(pname);
        return 0;
        break;
      case 'V':
        std::cerr << "This is " << pname << " version " << version
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

  if (argc - optind != 1)
  {
    std::cerr << std::endl << "Howdy there, wrong number of arguments."
              << std::endl;
    help(pname);
    return -1;
  }

  if (iseq == 1) ldirect = false;
  if (ilittle == 1) lbigend = false;
  if (igrads == 1) lgrads = true;

  try
  {
    char *regcmin = strdup(argv[optind++]);
 
    rcminp inpf(regcmin);
    domain_data d(inpf);
    char *datadir = strdup(inpf.valuec("dirglob"));
    rcmio rcmout(datadir, lbigend, ldirect);

    char *experiment = strdup(inpf.valuec("domname"));
    char dominfo[PATH_MAX];

    sprintf(dominfo, "%s%s%s.INFO", inpf.valuec("dirter"),
            separator, experiment);
    rcmout.read_domain(dominfo, d);

    bcdata b(d, inpf);

    regcmout outnc;
    outnc.experiment = experiment;
    outnc.fname = "ICBC_"+outnc.experiment+".nc";
    if (lgrads)
    {
      char ctlname[PATH_MAX];
      sprintf(ctlname, "ICBC_%s.ctl", experiment);
      outnc.ctl.open(ctlname, (char *) outnc.fname.c_str());
    }
    outnc.vl.addvar(vnames);

    bcNc bcnc(outnc, d);

    std::cout << "Processing ICBC";
    while ((rcmout.bc_read_tstep(b)) == 0)
    {
      std::cout << ".";
      std::cout.flush();
      bcnc.put_rec(b);
    }
    std::cout << " Done." << std::endl;

    if (lgrads)
      outnc.ctl.finalize( );
    free(datadir);
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
