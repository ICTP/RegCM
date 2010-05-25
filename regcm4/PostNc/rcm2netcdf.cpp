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
#include <rcminp.h>
#include <rcmNc.h>
#include <calc.h>

using namespace rcm;

static void help(char *pname);
static const char version[] = SVN_REV;

void help(char *pname)
{
  std::cerr << std::endl
      << "                 RegCM V4 ICTP NetCDF Postprocessor." << std::endl
      << std::endl
      << "This simple program converts binary output files from RegCM V4 "
      << "into NetCDF" << std::endl << "CF-1.4 convention compliant data files."
      << std::endl << std::endl << "I need ONE mandatory argument:"
      << std::endl << std::endl
      << "    regcm.in  - path to regcm.in of RegCM model v4" << std::endl
      << std::endl << "Example:" << std::endl << std::endl << "     " << pname
      << " [options] /home/regcm/Run/regcm.in"
      << std::endl << std::endl
      << "where options can be in:" << std::endl << std::endl
  << "   --sequential              : Set I/O non direct (direct access default)"
      << std::endl
  << "   --little_endian           : Set I/O endianess to LITTLE (BIG default)"
      << std::endl
  << "   --onlyatm/-a              : Process ATM file (default do all)"
      << std::endl
  << "   --onlysrf/-s              : Process SRF file (default do all)"
      << std::endl
  << "   --onlysub/-u              : Process SUB file (default do all)"
      << std::endl
  << "   --onlyrad/-r              : Process RAD file (default do all)"
      << std::endl
  << "   --onlyche/-c              : Process CHE file (default do all)"
      << std::endl
  << "   --startstep/-t [number]   : Start at timestep number 'number'"
      << std::endl
  << "   --nsteps/-n [number]      : Process just 'number' timesteps"
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
  int iseq, ilittle, istart, nsteps;
  bool onlyatm, onlysrf, onlysub, onlyrad, onlyche;
  ldirect = true;
  lbigend = true;
  onlyatm = false;
  onlysrf = false;
  onlysub = false;
  onlyrad = false;
  onlyche = false;
  iseq = 0;
  ilittle = 0;
  istart = 0;
  nsteps = -1;

  char *pname = basename(argv[0]);
  while (1)
  {
    static struct option long_options[] = {
      { "sequential", no_argument, &iseq, 1},
      { "little_endian", no_argument, &ilittle, 1},
      { "onlyatm", no_argument, 0, 'a'},
      { "onlysrf", no_argument, 0, 's'},
      { "onlysub", no_argument, 0, 'u'},
      { "onlyrad", no_argument, 0, 'r'},
      { "onlyche", no_argument, 0, 'c'},
      { "startstep", required_argument, 0, 't'},
      { "nsteps", required_argument, 0, 'n'},
      { "help", no_argument, 0, 'h'},
      { "version", no_argument, 0, 'V'},
      { 0, 0, 0, 0 }
    };
    int optind, c = 0;
    c = getopt_long (argc, argv, "asruct:n:hV",
                     long_options, &optind);
    if (c == -1) break;
    switch (c)
    {
      case 0:
        if (long_options[optind].flag != 0) break;
      case 'h':
        help(pname);
        return 0;
        break;
      case 'V':
        std::cerr << "This is " << pname << " version " << version
                  << std::endl;
        return 0;
      case 'a':
        onlyatm = true;
        break;
      case 's':
        onlysrf = true;
        break;
      case 'r':
        onlyrad = true;
        break;
      case 'u':
        onlysub = true;
        break;
      case 'c':
        onlyche = true;
        break;
      case 't':
        sscanf(optarg, "%d", &istart);
        break;
      case 'n':
        sscanf(optarg, "%d", &nsteps);
        break;
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

  if (onlyatm == false && onlysrf == false && onlyrad == false &&
      onlysub == false && onlyche == false)
  {
    onlyatm = onlysrf = onlyrad = onlysub = onlyche = true;
  }
  if (iseq == 1) ldirect = false;
  if (ilittle == 1) lbigend = false;

  try
  {
    char *regcmin = strdup(argv[optind++]);
    rcminp inpf(regcmin);

    char *experiment = strdup(inpf.valuec("domname"));
    char *outdir = strdup(inpf.valuec("dirout"));

    rcmio rcmout(outdir, lbigend, ldirect);
    header_data outhead(inpf);
    rcmout.read_header(outhead);

    char fname[PATH_MAX];
    char ctlname[PATH_MAX];

    if (rcmout.has_atm && onlyatm)
    {
      int recnum = 1;
      int astart = istart - 1;
      std::cout << "Found Atmospheric data ATM and processing";
      atmodata a(outhead);
      sprintf(fname, "ATM_%s_%d.nc", experiment, outhead.idate1);
      sprintf(ctlname, "ATM_%s_%d.ctl", experiment, outhead.idate1);
      gradsctl ctl(ctlname, fname);
      rcmNcAtmo atmnc(fname, experiment, outhead, ctl);
      atmcalc c(outhead);
      t_atm_deriv d;
      // Add Atmospheric variables
      while ((rcmout.atmo_read_tstep(a)) == 0)
      {
        if (nsteps > 0)
          if (recnum > nsteps) break;
        if (astart > 0)
        {
          astart --;
          atmnc.increment_time();
          continue;
        }
        std::cout << ".";
        std::cout.flush();
        c.do_calc(a, d);
        atmnc.put_rec(a, d, ctl);
        recnum ++;
      }
      ctl.finalize( );
      std::cout << " Done." << std::endl;
    }

    if (rcmout.has_srf && onlysrf)
    {
      int recnum = 1;
      int sstart = istart - 1;
      std::cout << "Found Surface data SRF and processing";
      srfdata s(outhead);
      sprintf(fname, "SRF_%s_%d.nc", experiment, outhead.idate1);
      sprintf(ctlname, "SRF_%s_%d.ctl", experiment, outhead.idate1);
      gradsctl ctl(ctlname, fname);
      rcmNcSrf srfnc(fname, experiment, outhead, ctl);
      // Add Surface variables
      srfcalc c(outhead);
      t_srf_deriv d;
      while ((rcmout.srf_read_tstep(s)) == 0)
      {
        if (nsteps > 0)
          if (recnum > nsteps) break;
        if (sstart > 0)
        {
          sstart --;
          srfnc.increment_time();
          continue;
        }
        std::cout << ".";
        std::cout.flush();
        c.do_calc(s, d);
        srfnc.put_rec(s, d, ctl);
        recnum ++;
      }
      ctl.finalize( );
      std::cout << " Done." << std::endl;
    }

    if (rcmout.has_rad && onlyrad)
    {
      int recnum = 1;
      int rstart = istart - 1;
      std::cout << "Found Radiation data RAD and processing";
      raddata r(outhead);
      sprintf(fname, "RAD_%s_%d.nc", experiment, outhead.idate1);
      sprintf(ctlname, "RAD_%s_%d.ctl", experiment, outhead.idate1);
      gradsctl ctl(ctlname, fname);
      rcmNcRad radnc(fname, experiment, outhead, ctl);
      // Add Radiation variables
      while ((rcmout.rad_read_tstep(r)) == 0)
      {
        if (nsteps > 0)
          if (recnum > nsteps) break;
        if (rstart > 0)
        {
          rstart --;
          radnc.increment_time();
          continue;
        }
        std::cout << ".";
        std::cout.flush();
        radnc.put_rec(r, ctl);
        recnum ++;
      }
      ctl.finalize( );
      std::cout << " Done." << std::endl;
    }

    if (rcmout.has_che && onlyche)
    {
      int recnum = 1;
      int cstart = istart - 1;
      std::cout << "Found Chemical data CHE and processing";
      chedata c(outhead);
      sprintf(fname, "CHE_%s_%d.nc", experiment, outhead.idate1);
      sprintf(ctlname, "CHE_%s_%d.ctl", experiment, outhead.idate1);
      gradsctl ctl(ctlname, fname);
      rcmNcChe chenc(fname, experiment, outhead, ctl);
      // Add Chemical tracers variables
      while ((rcmout.che_read_tstep(c)) == 0)
      {
        if (nsteps > 0)
          if (recnum > nsteps) break;
        if (cstart > 0)
        {
          cstart --;
          chenc.increment_time();
          continue;
        }
        std::cout << ".";
        std::cout.flush();
        chenc.put_rec(c, ctl);
        recnum ++;
      }
      ctl.finalize( );
      std::cout << " Done." << std::endl;
    }

    if (rcmout.has_sub && onlysub)
    {
      int recnum = 1;
      int ustart = istart - 1;
      std::cout << "Found Surface Subgrid data SUB and processing";
      subdom_data subdom(inpf);
      sprintf(fname, "%s%s%s%03d.INFO", inpf.valuec("dirter"),
              separator, inpf.valuec("domname"), inpf.valuei("nsg"));
      rcmout.read_subdom(outhead, subdom, fname);
      subdata s(outhead, subdom);
      sprintf(fname, "SUB_%s_%d.nc", experiment, outhead.idate1);
      sprintf(ctlname, "SUB_%s_%d.ctl", experiment, outhead.idate1);
      gradsctl ctl(ctlname, fname);
      rcmNcSub subnc(fname, experiment, outhead, subdom, ctl);
      subcalc c(subdom);
      t_srf_deriv d;
      // Add Subgrid variables
      while ((rcmout.sub_read_tstep(s)) == 0)
      {
        if (nsteps > 0)
          if (recnum > nsteps) break;
        if (ustart > 0)
        {
          ustart --;
          subnc.increment_time();
          continue;
        }
        std::cout << ".";
        std::cout.flush();
        c.do_calc(s, d);
        subnc.put_rec(s, d, ctl);
        recnum ++;
      }
      ctl.finalize( );
      std::cout << " Done." << std::endl;
    }

    outhead.free_space( );
    free(outdir);
    free(experiment);
    free(regcmin);
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
