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
  << "   --sequential                 : Set sequential (direct access default)"
      << std::endl
  << "   --little_endian              : Set endianess to LITTLE (BIG default)"
      << std::endl
  << "   --grads                      : Produce a CTL file for GrADS"
      << std::endl
  << "   --list/-l                    : Output list of names for -v option"
      << std::endl
  << "   --var/-v[all|name[,name]]    : Include only some vars (all default)"
      << std::endl
  << "   --tstart/-t YYYY[MM[DD[HH]]] : Start processing at this date"
      << std::endl
  << "   --tend/-e YYYY[MM[DD[HH]]]   : Stop processing at this date"
      << std::endl
  << "   --onlyatm/-a                 : Process ATM file (default do all)"
      << std::endl
  << "   --onlysrf/-s                 : Process SRF file (default do all)"
      << std::endl
  << "   --onlysub/-u                 : Process SUB file (default do all)"
      << std::endl
  << "   --onlyrad/-r                 : Process RAD file (default do all)"
      << std::endl
  << "   --onlyche/-c                 : Process CHE file (default do all)"
      << std::endl
  << "   --help/-h                    : Print this help"
      << std::endl
  << "   --version/-V                 : Print versioning information"
      << std::endl << std::endl;
   return;
}

int main(int argc, char *argv[])
{
  bool ldirect, lbigend, lgra;
  int iseq, ilittle, igra;
  bool onlyatm, onlysrf, onlysub, onlyrad, onlyche;
  int date1 = -1, date2= -1;
  std::string vnames = "all";
  ldirect = true;
  lbigend = true;
  lgra = false;
  onlyatm = false;
  onlysrf = false;
  onlysub = false;
  onlyrad = false;
  onlyche = false;
  iseq = 0;
  ilittle = 0;
  igra = 0;

  char *pname = basename(argv[0]);
  while (1)
  {
    static struct option long_options[] = {
      { "sequential", no_argument, &iseq, 1},
      { "little_endian", no_argument, &ilittle, 1},
      { "grads", no_argument, &igra, 1},
      { "onlyatm", no_argument, 0, 'a'},
      { "onlysrf", no_argument, 0, 's'},
      { "onlysub", no_argument, 0, 'u'},
      { "onlyrad", no_argument, 0, 'r'},
      { "onlyche", no_argument, 0, 'c'},
      { "list", no_argument, 0, 'l'},
      { "var", required_argument, 0, 'v'},
      { "tstart", required_argument, 0, 't'},
      { "tend", required_argument, 0, 'e'},
      { "help", no_argument, 0, 'h'},
      { "version", no_argument, 0, 'V'},
      { 0, 0, 0, 0 }
    };
    int optind, c = 0;
    c = getopt_long (argc, argv, "asruct:e:hVlv:",
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
      case 'l':
        std::cout << std::endl
                  << "List of names for " << pname << " is:" << std::endl
                  << std::endl
                  << "###############" << std::endl
                  << "ATM file source" << std::endl
                  << "###############" << std::endl
                  << "  u  :  U component of wind (westerly)" << std::endl
                  << "  v  :  V component of wind (southerly)" << std::endl
                  << "  o  :  Vertical pressure velocity" << std::endl
                  << "  dv :  Wind divergence" << std::endl
                  << "  vr :  Wind vorticity" << std::endl
                  << "  p  :  Pressure" << std::endl
                  << "  hg :  Geopotential Height" << std::endl
                  << "  t  :  Temperature" << std::endl
                  << "  td :  Dew point temperature" << std::endl
                  << "  tp :  Potential temperature" << std::endl
                  << "  rh :  Relative humidity" << std::endl
                  << "  qv :  Water Vapor mixing ratio" << std::endl
                  << "  qc :  Cloud water mixing ratio" << std::endl
                  << "  tpr:  Total precipitation" << std::endl
                  << "  tgb:  Temperature of ground" << std::endl
                  << "  swt:  Soil water content" << std::endl
                  << "  rno:  Surface runoff" << std::endl
                  << "###############" << std::endl
                  << "SRF file source" << std::endl
                  << "###############" << std::endl
                  << "  u10:  10m wind U component (westerly)" << std::endl
                  << "  v10:  10m wind V component (southerly)" << std::endl
                  << "  udg:  Wind drag" << std::endl
                  << "  tg :  Ground temperature" << std::endl
                  << "  tfl:  Foliage temperature" << std::endl
                  << "  t2 :  2m air temperature" << std::endl
                  << "  r2 :  2m air relative humidity" << std::endl
                  << "  q2 :  2m air water vapor mixing ratio" << std::endl
                  << "  sm :  Soil moisture" << std::endl
                  << "  tpr:  Total precipitation" << std::endl
                  << "  evp:  Total evapotranspiration" << std::endl
                  << "  rno:  Surface runoff" << std::endl
                  << "  scv:  Snow precipitation" << std::endl
                  << "  sen:  Sensible heat" << std::endl
                  << "  flw:  Upward LW" << std::endl
                  << "  fsw:  Downward SW" << std::endl
                  << "  fld:  Downward LW" << std::endl
                  << "  sin:  Solar input energy" << std::endl
                  << "  prc:  Convective precipitation" << std::endl
                  << "  zpb:  PBL height" << std::endl
                  << "  tga:  Maximum ground temperature" << std::endl
                  << "  tgi:  Minimum ground temperature" << std::endl
                  << "  t2a:  Maximum 2m air temperature" << std::endl
                  << "  t2i:  Minimum 2m air temperature" << std::endl
                  << "  wma:  Maximum speed 10m wind" << std::endl
                  << "  psi:  Minimum pressure" << std::endl
                  << "###############" << std::endl
                  << "RAD file source" << std::endl
                  << "###############" << std::endl
                  << "  cld:  Cloud fractional cover" << std::endl
                  << "  clw:  Cloud liquid water path" << std::endl
                  << "  qrs:  Solar heating rate" << std::endl
                  << "  qrl:  LW cooling rate" << std::endl
                  << "  frs:  Surface absorbed solar flux" << std::endl
                  << "  frl:  Longwave cooling of surface flux" << std::endl
                  << "  crs:  Clearsky SW solar column abs. flux" << std::endl
                  << "  css:  Clearsky SW solar surface flux" << std::endl
                  << "  crl:  Clearsky LW column flux at TOA" << std::endl
                  << "  csl:  Clearsky LW flux at surface" << std::endl
                  << "  soi:  Incident solar flux" << std::endl
                  << "  sab:  SW absorption rate" << std::endl
                  << "  fir:  LW absorption rate" << std::endl
                  << "###############" << std::endl
                  << "CHE file source" << std::endl
                  << "###############" << std::endl
                  << "  trc:  Tracers mixing ratio" << std::endl
                  << "  axt:  Aerosol extinction coefficient" << std::endl
                  << "  asa:  Aerosol single scattering albedo" << std::endl
                  << "  agf:  Aerosol asymmetry parameter" << std::endl
                  << "  cbd:  Tracer deposition" << std::endl
                  << "  wdl:  Wet deposition large scale" << std::endl
                  << "  wdc:  Wet deposition convective" << std::endl
                  << "  drd:  Dry deposition" << std::endl
                  << "  xgs:  Gas conversion" << std::endl
                  << "  xaq:  Aqueous conversion" << std::endl
                  << "  tem:  Emission of tracer from surface" << std::endl
                  << "  fst:  TOA SW radiative forcing" << std::endl
                  << "  fss:  Surface SW radiative forcing" << std::endl
                  << "  flt:  TOA LW radiative forcing" << std::endl
                  << "  fls:  Surface LW radiative forcing" << std::endl
                  << "###############" << std::endl
                  << "SUB file source" << std::endl
                  << "###############" << std::endl
                  << "  u10:  10m wind U component (westerly)" << std::endl
                  << "  v10:  10m wind V component (southerly)" << std::endl
                  << "  udg:  Wind drag" << std::endl
                  << "  tg :  Ground temperature" << std::endl
                  << "  tfl:  Foliage temperature" << std::endl
                  << "  t2 :  2m air temperature" << std::endl
                  << "  r2 :  2m air relative humidity" << std::endl
                  << "  q2 :  2m air water vapor mixing ratio" << std::endl
                  << "  sm :  Soil moisture" << std::endl
                  << "  tpr:  Total precipitation" << std::endl
                  << "  evp:  Total evapotranspiration" << std::endl
                  << "  rno:  Surface runoff" << std::endl
                  << "  scv:  Snow precipitation" << std::endl
                  << "  sen:  Sensible heat" << std::endl
                  << "  prc:  Convective precipitation" << std::endl
                  << std::endl;
        return 0;
      case 'v':
        vnames = optarg;
        break;
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
        sscanf(optarg, "%d", &date1);
        break;
      case 'e':
        sscanf(optarg, "%d", &date2);
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
  if (igra == 1) lgra = true;

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
    if (date1 < 1) date1 = inpf.valuei("idate1");
    if (date2 < 1) date2 = inpf.valuei("idate2");
    // Year only
    if (date1 < 10000) date1 = date1*1000000+10100;
    if (date2 < 10000) date2 = date2*1000000+10100;
    // Year+month only
    if (date1 < 1000000) date1 = date1*10000+100;
    if (date2 < 1000000) date2 = date2*10000+100;
    // Year+month+day only
    if (date1 < 100000000) date1 = date1*100;
    if (date2 < 100000000) date2 = date2*100;

    t_time_interval t;
    t.idate0 = date1;
    t.idate1 = date2;

    if (rcmout.has_atm && onlyatm)
    {
      std::cout << "Found Atmospheric data ATM and processing";
      atmodata a(outhead, t);
      regcmout outnc;
      outnc.experiment = experiment;
      sprintf(fname, "ATM_%s_%d-%d.nc", experiment, a.date0, a.date1);
      outnc.fname = fname;
      if (lgra)
      {
        char ctlname[PATH_MAX];
        sprintf(ctlname, "ATM_%s_%d-%d.ctl", experiment, a.date0, a.date1);
        outnc.ctl.open(ctlname, (char *) outnc.fname.c_str());
	outnc.ctl.doit = true;
      }
      outnc.vl.addvar(vnames);

      rcmNcAtmo atmnc(outnc, outhead, t);
      atmcalc c(outhead);
      t_atm_deriv d;
      // Add Atmospheric variables
      while ((rcmout.atmo_read_tstep(a)) == 0)
      {
        std::cout << ".";
        std::cout.flush();
        c.do_calc(a, d);
        atmnc.put_rec(a, d);
      }
      if (lgra)
        outnc.ctl.finalize( );
      std::cout << " Done." << std::endl;
    }

    if (rcmout.has_srf && onlysrf)
    {
      std::cout << "Found Surface data SRF and processing";
      srfdata s(outhead, t);
      regcmout outnc;
      outnc.experiment = experiment;
      sprintf(fname, "SRF_%s_%d-%d.nc", experiment, s.date0, s.date1);
      outnc.fname = fname;
      if (lgra)
      {
        char ctlname[PATH_MAX];
        sprintf(ctlname, "SRF_%s_%d-%d.ctl", experiment, s.date0, s.date1);
        outnc.ctl.open(ctlname, (char *) outnc.fname.c_str());
	outnc.ctl.doit = true;
      }
      outnc.vl.addvar(vnames);

      rcmNcSrf srfnc(outnc, outhead, t);
      // Add Surface variables
      srfcalc c(outhead);
      t_srf_deriv d;
      while ((rcmout.srf_read_tstep(s)) == 0)
      {
        std::cout << ".";
        std::cout.flush();
        c.do_calc(s, d);
        srfnc.put_rec(s, d);
      }
      if (lgra)
        outnc.ctl.finalize( );
      std::cout << " Done." << std::endl;
    }

    if (rcmout.has_rad && onlyrad)
    {
      std::cout << "Found Radiation data RAD and processing";
      raddata r(outhead, t);
      regcmout outnc;
      outnc.experiment = experiment;
      sprintf(fname, "RAD_%s_%d-%d.nc", experiment, r.date0, r.date1);
      outnc.fname = fname;
      if (lgra)
      {
        char ctlname[PATH_MAX];
        sprintf(ctlname, "RAD_%s_%d-%d.ctl", experiment, r.date0, r.date1);
        outnc.ctl.open(ctlname, (char *) outnc.fname.c_str());
	outnc.ctl.doit = true;
      }
      outnc.vl.addvar(vnames);

      rcmNcRad radnc(outnc, outhead, t);
      // Add Radiation variables
      while ((rcmout.rad_read_tstep(r)) == 0)
      {
        std::cout << ".";
        std::cout.flush();
        radnc.put_rec(r);
      }
      if (lgra)
        outnc.ctl.finalize( );
      std::cout << " Done." << std::endl;
    }

    if (rcmout.has_che && onlyche)
    {
      std::cout << "Found Chemical data CHE and processing";
      chedata c(outhead, t);
      regcmout outnc;
      outnc.experiment = experiment;
      sprintf(fname, "CHE_%s_%d-%d.nc", experiment, c.date0, c.date1);
      outnc.fname = fname;
      if (lgra)
      {
        char ctlname[PATH_MAX];
        sprintf(ctlname, "CHE_%s_%d-%d.ctl", experiment, c.date0, c.date1);
        outnc.ctl.open(ctlname, (char *) outnc.fname.c_str());
	outnc.ctl.doit = true;
      }
      outnc.vl.addvar(vnames);

      rcmNcChe chenc(outnc, outhead, t);
      // Add Chemical tracers variables
      while ((rcmout.che_read_tstep(c)) == 0)
      {
        std::cout << ".";
        std::cout.flush();
        chenc.put_rec(c);
      }
      if (lgra)
        outnc.ctl.finalize( );
      std::cout << " Done." << std::endl;
    }

    if (rcmout.has_sub && onlysub)
    {
      std::cout << "Found Surface Subgrid data SUB and processing";
      subdom_data subdom(inpf);
      sprintf(fname, "%s%s%s%03d.INFO", inpf.valuec("dirter"),
              separator, inpf.valuec("domname"), inpf.valuei("nsg"));
      rcmout.read_subdom(outhead, subdom, fname);
      subdata s(outhead, subdom, t);
      regcmout outnc;
      outnc.experiment = experiment;
      sprintf(fname, "SUB_%s_%d-%d.nc", experiment, s.date0, s.date1);
      outnc.fname = fname;
      if (lgra)
      {
        char ctlname[PATH_MAX];
        sprintf(ctlname, "SUB_%s_%d-%d.ctl", experiment, s.date0, s.date1);
        outnc.ctl.open(ctlname, (char *) outnc.fname.c_str());
	outnc.ctl.doit = true;
      }
      outnc.vl.addvar(vnames);

      rcmNcSub subnc(outnc, outhead, subdom, t);
      subcalc c(subdom);
      t_srf_deriv d;
      // Add Subgrid variables
      while ((rcmout.sub_read_tstep(s)) == 0)
      {
        std::cout << ".";
        std::cout.flush();
        c.do_calc(s, d);
        subnc.put_rec(s, d);
      }
      if (lgra)
        outnc.ctl.finalize( );
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
