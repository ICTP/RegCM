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

#include <rcmio.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/stat.h>
#include <cmath>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

using namespace rcm;

static const int mlen[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
static bool fexist(char *fname);
static void reduce(float *a, float *b, int ni, int nj, int nsg);
static int calcnfiles(int date0, int date1);
static int calcnsteps(int date0, int date1, int dt);
static int nextmonth(int rdate);
static bool isleap(int year);
static int mdays(int mon, int year);
static int idatei(int y, int m, int d, int h);

bool fexist(char *fname)
{
  struct stat finfo;
  int retval = stat(fname, &finfo);
  if (retval == 0) return true;
  return false;
}

void reduce(float *a, float *b, int ni, int nj, int nsg)
{
  int count = 0;
  for (int i = nsg; i < ni-nsg; i ++)
  {
    for (int j = nsg; j < nj-nsg; j ++)
    {
       b[count] = a[i*nj+j];
       count ++;
    }
  }
  return;
}

int calcnfiles(int date0, int date1)
{
  // Calculate number of months (i.e. files) to read
  unsigned int base = date0;
  unsigned int basey = base/1000000;
  int year1 = (int) basey;
  base = base-basey*1000000;
  unsigned int basem = base/10000;
  int month1 = (int) basem;
  base = date1;
  basey = base/1000000;
  int year2 = (int) basey;
  base = base-basey*1000000;
  basem = base/10000;
  int month2 = (int) basem;
  return (year2-year1)*12+(month2-month1);
}

int nextmonth(int rdate)
{
  rcmdate d(rdate);
  d.basem += 1;
  if (d.basem > 12)
  {
    d.basey += 1;
    d.basem = 1;
  }
  return d.idate( );
}

int calcnsteps(int date0, int date1, int dt)
{
  rcmdate d(date0);
  unsigned long hdiff = d.datediffh(date1);
  return (int) (hdiff/(unsigned long) dt);
}

typedef struct {
  size_t len;
  char *data;
} t_fortran_record;

header_data::header_data(rcminp &in)
{
  hsigf = hsigm = ht = htsd = veg2d = landuse = 0;
  xlat = xlon = xmap = dmap = coriol = mask = 0;
  idate0 = in.valuei("idate0");
  idate1 = in.valuei("idate1");
  idate2 = in.valuei("idate2");
  ifrest = in.valueb("ifrest");
  radfrq = in.valuef("radfrq");
  abatm = in.valuef("abatm");
  abemh = in.valuef("abemh");
  dt = in.valuef("dt");
  ibdyfrq = in.valuef("ibdyfrq");
  igcc = (unsigned short) in.valuei("igcc");
  iocnflx = (unsigned short) in.valuei("iocnflx");
  ipgf = (unsigned short) in.valuei("ipgf");
  iemiss = (unsigned short) in.valuei("iemiss");
  lakemod = (unsigned short) in.valuei("lakemod");
  ichem = (unsigned short) in.valuei("ichem");
  if (ichem == 1)
  {
    ntr = in.valuei("ntr");
    nbin = in.valuei("nbin");
    trnames = in.valuec("inpchtrname");
  }
}

header_data::~header_data( )
{
  free_space( );
}

void header_data::free_space( )
{
  if (hsigf) delete [ ] hsigf;
  if (hsigm) delete [ ] hsigm;
  if (ht) delete [ ] ht;
  if (htsd) delete [ ] htsd;
  if (veg2d) delete [ ] veg2d;
  if (landuse) delete [ ] landuse;
  if (xlat) delete [ ] xlat;
  if (xlon) delete [ ] xlon;
  if (xmap) delete [ ] xmap;
  if (dmap) delete [ ] dmap;
  if (coriol) delete [ ] coriol;
  if (mask) delete [ ] mask;
  hsigf = hsigm = ht = htsd = veg2d = landuse = 0;
  xlat = xlon = xmap = dmap = coriol = mask = 0;
}

const char *header_data::lbcs( )
{
  switch (iboudy)
  {
    case 0:
     return "Fixed";
     break;
    case 1:
     return "Relaxation, linear technique.";
     break;
    case 2:
      return "Time-dependent";
      break;
    case 3:
      return "Time and inflow/outflow dependent.";
      break;
    case 4:
      return "Sponge (Perkey & Kreitzberg, MWR 1976).";
      break;
    case 5:
      return "Relaxation, exponential technique.";
      break;
    default:
      return "Unknown or not specified";
  }
}

const char *header_data::cums( )
{
  switch(icup)
  {
    case 1:
      return "Kuo";
      break;
    case 2:
      return "Grell";
      break;
    case 3:
      return "Betts-Miller (1986)";
      break;
    case 4:
      return "Emanuel (1991)";
      break;
    default:
      return "Unknown or not specified";
  }
}

const char *header_data::cumsclos( )
{
  switch(igcc)
  {
    case 1:
      return "Arakawa & Schubert (1974)";
      break;
    case 2:
      return "Fritsch & Chappell (1980)";
      break;
    default:
      return "Unknown or not specified";
  }
}

const char *header_data::pbls( )
{
  switch (ibltyp)
  {
    case 0:
      return "Frictionless";
      break;
    case 1:
      return "Holtslag PBL (Holtslag, 1990)";
      break;
    default:
      return "Unknown or not specified";
  }
}

const char *header_data::moists( )
{
  switch(imoist)
  {
    case 1:
     return "Explicit moisture (SUBEX; Pal et al 2000)";
     break;
    default:
      return "Unknown or not specified";
  }
}

const char *header_data::ocnflxs( )
{
  switch(iocnflx)
  {
    case 1:
      return "Use BATS1e Monin-Obukhov";
      break;
    case 2:
      return "Zeng et al (1998)";
      break;
    default:
      return "Unknown or not specified";
  }
}

const char *header_data::pgfs( )
{
  switch(ipgf)
  {
    case 0:
      return "Use full fields";
      break;
    case 1:
      return "Hydrostatic deduction with perturbation temperature";
      break;
    default:
      return "Unknown or not specified";
  }
}

const char *header_data::emisss( )
{
  switch(iemiss)
  {
    case 0:
      return "No";
      break;
    case 1:
      return "Yes";
      break;
    default:
      return "Unknown or not specified";
  }
}

const char *header_data::lakes( )
{
  switch(lakemod)
  {
    case 0:
      return "No";
      break;
    case 1:
      return "Yes";
      break;
    default:
      return "Unknown or not specified";
  }
}

const char *header_data::chems( )
{
  switch(ichem)
  {
    case 0:
      return "Not active";
      break;
    case 1:
      return "Active";
      break;
    default:
      return "Unknown or not specified";
  }
}

bcdata::bcdata(domain_data &d, rcminp &in, t_time_interval &t)
{
  try
  {
    dt = in.valuei("ibdyfrq");
  }
  catch (...)
  {
    std::cout << std::endl << "Boundary frequency not specified." << std::endl
              << "Defaulting to 6 hours" << std::endl;
    dt = 6;
  }
  date0 = in.valuei("globidate1");
  date1 = in.valuei("globidate2");
  rcmdate d1(date0);
  rdate = date0;
  if (t.idate0 < date0)
  {
    std::cerr << std::endl << "Warning: requested from " << t.idate0
              << ", will start at " << date0 << std::endl;
  }
  else if (t.idate0 > date0)
  {
    rcmdate d2(t.idate0);
    // Same month
    if (d1.idatendh( ) != d2.idatendh( )) rdate = d2.idatendh( )+100;
    date0 = t.idate0;
  }
  if (t.idate1 < date1) date1 = t.idate1;
  else if (t.idate1 > date1)
  {
    std::cerr << std::endl << "Warning: requested up to " << t.idate1
              << ", will stop at " << date1 << std::endl;
  }
  nfiles = calcnfiles(date0, date1);
  nsteps = calcnsteps(date0, date1, (int) dt);
  size2D = d.nx*d.ny;
  size3D = size2D*d.nz;

  nz = d.nz;
  ehso4 = in.valueb("ehso4");
  if (strcmp(in.valuec("lsmtyp"), "BATS") == 0)
    usgs = false;
  else if (strcmp(in.valuec("lsmtyp"), "USGS") == 0)
    usgs = true;
  buffer = 0;
  strncpy(name, in.valuec("domname"), 256);
  strncpy(ssttyp, in.valuec("ssttyp"), 16);
  strncpy(dattyp, in.valuec("dattyp"), 16);
}

bcdata::~bcdata()
{
  if (buffer) delete [] buffer;
}

atmodata::atmodata(header_data &h, t_time_interval &t)
{
  dt = h.dto;
  date0 = h.idate1;
  date1 = h.idate2;
  rcmdate d1(date0);
  rdate = date0;
  if (t.idate0 < date0)
  {
    std::cerr << std::endl << "Warning: requested from " << t.idate0
              << ", will start at " << date0 << std::endl;
  }
  else if (t.idate0 > date0)
  {
    rcmdate d2(t.idate0);
    // Same month
    if (d1.idatendh( ) != d2.idatendh( )) rdate = d2.idatendh( )+100;
    date0 = t.idate0;
  }
  if (t.idate1 < date1) date1 = t.idate1;
  else if (t.idate1 > date1)
  {
    std::cerr << std::endl << "Warning: requested up to " << t.idate1
              << ", will stop at " << date1 << std::endl;
  }
  nfiles = calcnfiles(date0, date1);
  nsteps = calcnsteps(date0, date1, (int) dt);

  n3D = 6;
  n2D = 5;
  size2D = h.nx*h.ny;
  size3D = size2D*h.nz;
  nvals = n3D*size3D + n2D*size2D;
  datasize = nvals*sizeof(float);
  buffer = new char[datasize];
  u = (float *) buffer;
  v = u + size3D;
  omega = v + size3D;
  this->t = omega + size3D;
  qv = this->t + size3D;
  qc = qv + size3D;
  ps = qc + size3D;
  tpr = ps + size2D;
  tgb = tpr + size2D;
  swt = tgb + size2D;
  rno = swt + size2D;
}

atmodata::~atmodata( )
{
  delete [] buffer;
}

raddata::raddata(header_data &h, t_time_interval &t)
{
  dt = h.dtr;
  date0 = h.idate1;
  date1 = h.idate2;
  rcmdate d1(date0);
  rdate = date0;
  if (t.idate0 < date0)
  {
    std::cerr << std::endl << "Warning: requested from " << t.idate0
              << ", will start at " << date0 << std::endl;
  }
  else if (t.idate0 > date0)
  {
    rcmdate d2(t.idate0);
    // Same month
    if (d1.idatendh( ) != d2.idatendh( )) rdate = d2.idatendh( )+100;
    date0 = t.idate0;
  }
  if (t.idate1 < date1) date1 = t.idate1;
  else if (t.idate1 > date1)
  {
    std::cerr << std::endl << "Warning: requested up to " << t.idate1
              << ", will stop at " << date1 << std::endl;
  }
  nfiles = calcnfiles(date0, date1);
  nsteps = calcnsteps(date0, date1, (int) dt);

  n3D = 4;
  n2D = 10;
  size2D = h.nx*h.ny;
  size3D = size2D*h.nz;
  nvals = n3D*size3D + n2D*size2D;
  datasize = nvals*sizeof(float);
  buffer = new char[datasize];
  cld = (float *) buffer;
  clwp = cld + size3D;
  qrs = clwp + size3D;
  qrl = qrs + size3D;
  frsa = qrl + size3D;
  frla = frsa + size2D;
  clrst = frla + size2D;
  clrss = clrst + size2D;
  clrlt = clrss + size2D;
  clrls = clrlt + size2D;
  solin = clrls + size2D;
  sabtp = solin + size2D;
  firtp = sabtp + size2D;
  ps = firtp + size2D;
}

raddata::~raddata( )
{
  delete [] buffer;
}

chedata::chedata(header_data &h, t_time_interval &t)
{
  dt = h.dtc;
  date0 = h.idate1;
  date1 = h.idate2;
  rcmdate d1(date0);
  rdate = date0;
  if (t.idate0 < date0)
  {
    std::cerr << std::endl << "Warning: requested from " << t.idate0
              << ", will start at " << date0 << std::endl;
  }
  else if (t.idate0 > date0)
  {
    rcmdate d2(t.idate0);
    // Same month
    if (d1.idatendh( ) != d2.idatendh( )) rdate = d2.idatendh( )+100;
    date0 = t.idate0;
  }
  if (t.idate1 < date1) date1 = t.idate1;
  else if (t.idate1 > date1)
  {
    std::cerr << std::endl << "Warning: requested up to " << t.idate1
              << ", will stop at " << date1 << std::endl;
  }
  nfiles = calcnfiles(date0, date1);
  nsteps = calcnsteps(date0, date1, (int) dt);

  ntr = h.ntr;
  n3D = 13;
  n2D = 75;
  size2D = h.nx*h.ny;
  size3D = size2D*h.nz;
  nvals = n3D*size3D + n2D*size2D;
  datasize = nvals*sizeof(float);
  buffer = new char[datasize];
  trac3D = (float *) buffer;
  aext8 = trac3D + 10*size3D;
  assa8 = aext8 + size3D;
  agfu8 = assa8 + size3D;
  trac2D = agfu8 + size3D;
  acstoarf = trac2D + 70*size2D;
  acstsrrf = acstoarf + size2D;
  acstalrf = acstsrrf + size2D;
  acssrlrf = acstalrf + size2D;
  ps = acssrlrf + size2D;
}

chedata::~chedata( )
{
  delete [] buffer;
}

srfdata::srfdata(header_data &h, t_time_interval &t)
{
  dt = h.dtb;
  date0 = h.idate1;
  date1 = h.idate2;
  rcmdate d1(date0);
  rdate = date0;
  if (t.idate0 < date0)
  {
    std::cerr << std::endl << "Warning: requested from " << t.idate0
              << ", will start at " << date0 << std::endl;
  }
  else if (t.idate0 > date0)
  {
    rcmdate d2(t.idate0);
    // Same month
    if (d1.idatendh( ) != d2.idatendh( )) rdate = d2.idatendh( )+100;
    date0 = t.idate0;
  }
  if (t.idate1 < date1) date1 = t.idate1;
  else if (t.idate1 > date1)
  {
    std::cerr << std::endl << "Warning: requested up to " << t.idate1
              << ", will stop at " << date1 << std::endl;
  }
  nfiles = calcnfiles(date0, date1);
  nsteps = calcnsteps(date0, date1, (int) dt);

  n2D = 27;
  size2D = h.nx*h.ny;
  nvals = n2D*size2D;
  datasize = nvals*sizeof(float);
  buffer = new char[datasize];
  u10m = (float *) buffer;
  v10m = u10m + size2D;
  uvdrag = v10m + size2D;
  tg = uvdrag + size2D;
  tlef = tg + size2D;
  t2m = tlef + size2D;
  q2m = t2m + size2D;
  smw = q2m + size2D;
  tpr = smw + 2*size2D;
  evp = tpr + size2D;
  runoff = evp + size2D;
  scv = runoff + size2D;
  sena = scv + size2D;
  flw = sena + size2D;
  fsw = flw + size2D;
  flwd = fsw + size2D;
  sina = flwd + size2D;
  prcv = sina + size2D;
  psb = prcv + size2D;
  zpbl = psb + size2D;
  tgmax = zpbl + size2D;
  tgmin = tgmax + size2D;
  t2max = tgmin + size2D;
  t2min = t2max + size2D;
  w10max = t2min + size2D;
  ps_min = w10max + size2D;
}

srfdata::~srfdata( )
{
  delete [] buffer;
}

domain_data::domain_data(rcminp &in)
{
  n2D = 14;
  nx = ny = nz = 0;
  buffer = 0;
  hsigf = hsigm = 0;
  strncpy(name, in.valuec("domname"), 256);
  ntypec = in.valuei("ntypec");
  anal = in.valueb("ifanal");
  smthbdy = in.valueb("smthbdy");
  lakadj = in.valueb("lakadj");
  fudge_lnd = in.valueb("fudge_lnd");
  fudge_tex = in.valueb("fudge_tex");
  ntex = in.valuei("ntex");
  h2opct = in.valuef("h2opct");
}

domain_data::~domain_data( )
{
  if (buffer) delete [ ] buffer;
  if (hsigf) delete [ ] hsigf;
  if (hsigm) delete [ ] hsigm;
}

subdata::subdata(header_data &h, subdom_data &s, t_time_interval &t)
{
  dt = h.dtb;
  date0 = h.idate1;
  date1 = h.idate2;
  rcmdate d1(date0);
  rdate = date0;
  if (t.idate0 < date0)
  {
    std::cerr << std::endl << "Warning: requested from " << t.idate0
              << ", will start at " << date0 << std::endl;
  }
  else if (t.idate0 > date0)
  {
    rcmdate d2(t.idate0);
    // Same month
    if (d1.idatendh( ) != d2.idatendh( )) rdate = d2.idatendh( )+100;
    date0 = t.idate0;
  }
  if (t.idate1 < date1) date1 = t.idate1;
  else if (t.idate1 > date1)
  {
    std::cerr << std::endl << "Warning: requested up to " << t.idate1
              << ", will stop at " << date1 << std::endl;
  }
  nfiles = calcnfiles(date0, date1);
  nsteps = calcnsteps(date0, date1, (int) dt);

  n2D = 16;
  size2D = s.nx*s.ny;
  nvals = n2D*size2D;
  datasize = nvals*sizeof(float);
  buffer = new char[datasize];
  u10m = (float *) buffer;
  v10m = u10m + size2D;
  uvdrag = v10m + size2D;
  tg = uvdrag + size2D;
  tlef = tg + size2D;
  t2m = tlef + size2D;
  q2m = t2m + size2D;
  smw = q2m + size2D;
  tpr = smw + 2*size2D;
  evp = tpr + size2D;
  runoff = evp + size2D;
  scv = runoff + size2D;
  sena = scv + size2D;
  prcv = sena + size2D;
  psb = prcv + size2D;
}

subdata::~subdata( )
{
  delete [] buffer;
}

rcmio::rcmio(char *directory, bool lbig, bool ldirect)
{
  // big or little endian swapping ?
  if ((littlearch() && lbig) || (bigarch() && ! lbig)) doswap = true;
  else doswap = false;

  // direct access or sequential file ?
  doseq = ! ldirect;
  if (directory != NULL)
    strncpy(outdir, directory, PATH_MAX);

  // All to false
  initatm = false;
  has_atm = false;
  initche = false;
  has_che = false;
  initsrf = false;
  has_srf = false;
  initrad = false;
  has_rad = false;
  initsub = false;
  has_sub = false;
  initbc  = false;
  storage = 0;
}

subdom_data::subdom_data(rcminp &in)
{
  nsg = 1;
  ht = htsd = landuse = xlat = xlon = xmap = coriol = mask = 0;
  strncpy(name, in.valuec("domname"), 256);
  ntypec_s = in.valuei("ntypec_s");
  fudge_lnd_s = in.valueb("fudge_lnd_s");
  fudge_tex_s = in.valueb("fudge_tex_s");
}

subdom_data::~subdom_data( )
{
  if (ht) delete [] ht;
  if (htsd) delete [] htsd;
  if (landuse) delete [] landuse;
  if (xlat) delete [] xlat;
  if (xlon) delete [] xlon;
  if (xmap) delete [] xmap;
  if (coriol) delete [] coriol;
  if (mask) delete [] mask;
}

void rcmio::read_domain(char *name, domain_data &d)
{
  std::ifstream rcmf;
  rcmf.open(name, std::ios::binary);
  if (! rcmf.good()) throw "Invalid Input. Cannot open.";

  // Read entire file
  size_t len;
  rcmf.seekg (0, std::ios::end);
  len = rcmf.tellg();
  rcmf.seekg (0, std::ios::beg);

  // Check for minimal size for reading in basic informations
  size_t minsize = 3*sizeof(int);
  if (doseq) minsize += sizeof(int);

  if (len <= minsize)
  {
    throw "DOMAIN.INFO file size is too short.";
  }

  char *buffer = new char[len];
  rcmf.read(buffer, len);
  rcmf.close();

  char *buf = buffer;
  if (doseq) buf = buf+sizeof(int);

  // Read dimensions
  d.nx = intvalfrombuf(buf); buf = buf + sizeof(int);
  d.ny = intvalfrombuf(buf); buf = buf + sizeof(int);
  d.nz = intvalfrombuf(buf); buf = buf + sizeof(int);
  if (d.nx <= 0 || d.ny <= 0 || d.nz <= 0)
  {
    throw "DOMAIN.INFO dims inconsistent";
  }
  // Now I have record len of fortran I/O
  int nvals = d.nx*d.ny;
  d.ds = rvalfrombuf(buf); buf = buf + sizeof(float);
  d.dsmeters = d.ds;
  d.clat = rvalfrombuf(buf); buf = buf + sizeof(float);
  d.clon = rvalfrombuf(buf); buf = buf + sizeof(float);
  d.xplat = rvalfrombuf(buf); buf = buf + sizeof(float);
  d.xplon = rvalfrombuf(buf); buf = buf + sizeof(float);
  d.grdfac = rvalfrombuf(buf); buf = buf + sizeof(float);
  d.proj[6] = 0;
  memcpy(d.proj, buf, 6); buf = buf + 6;
  // Calculate and store half sigma levels
  if (d.hsigf) delete [] d.hsigf;
  if (d.hsigm) delete [] d.hsigm;
  d.hsigf = new float[d.nz+1];
  d.hsigm = new float[d.nz];
  for (unsigned int i = 0; i < d.nz+1; i ++)
  {
    d.hsigf[i] = rvalfrombuf(buf); buf = buf + sizeof(float);
  }
  for (unsigned int i = 0; i < d.nz; i ++)
    d.hsigm[i] = d.hsigf[i]+0.5f*(d.hsigf[i+1]-d.hsigf[i]);
  d.ptop = rvalfrombuf(buf); buf = buf + sizeof(float);
  d.ptop *= 10.0f; // Put in hPa from cbar
  // skip I/O flags
  buf = buf + sizeof(int);
  buf = buf + sizeof(int);
  d.trlat1 = rvalfrombuf(buf); buf = buf + sizeof(float);
  d.trlat2 = rvalfrombuf(buf); buf = buf + sizeof(float);
  
  // First record is mostly empty except above infos
  buf = buffer + nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  d.buffer = new char[nvals*d.n2D*sizeof(float)];
  char *p = (char *) d.buffer;

  d.ht = (float *) p;
  vectorfrombuf(buf, d.ht, nvals);
  buf += nvals*sizeof(float);
  p += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  d.htsd = (float *) p;
  vectorfrombuf(buf, d.htsd, nvals);
  buf += nvals*sizeof(float);
  p += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  d.landuse = (float *) p;
  vectorfrombuf(buf, d.landuse, nvals);
  buf += nvals*sizeof(float);
  p += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  d.xlat = (float *) p;
  vectorfrombuf(buf, d.xlat, nvals);
  buf += nvals*sizeof(float);
  p += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  d.xlon = (float *) p;
  vectorfrombuf(buf, d.xlon, nvals);
  buf += nvals*sizeof(float);
  p += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  d.dlat = (float *) p;
  vectorfrombuf(buf, d.dlat, nvals);
  buf += nvals*sizeof(float);
  p += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  d.dlon = (float *) p;
  vectorfrombuf(buf, d.dlon, nvals);
  buf += nvals*sizeof(float);
  p += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  d.xmap = (float *) p;
  vectorfrombuf(buf, d.xmap, nvals);
  buf += nvals*sizeof(float);
  p += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  d.dmap = (float *) p;
  vectorfrombuf(buf, d.dmap, nvals);
  buf += nvals*sizeof(float);
  p += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  d.coriol = (float *) p;
  vectorfrombuf(buf, d.coriol, nvals);
  buf += nvals*sizeof(float);
  p += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  d.snowam = (float *) p;
  vectorfrombuf(buf, d.snowam, nvals);
  buf += nvals*sizeof(float);
  p += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  d.mask = (float *) p;
  vectorfrombuf(buf, d.mask, nvals);
  buf += nvals*sizeof(float);
  p += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  delete [ ] buffer;
  return;
}

void rcmio::read_subdom(header_data &h, subdom_data &s, char *fname)
{
  std::ifstream rcmf;
  rcmf.open(fname, std::ios::binary);
  if (! rcmf.good()) throw "Invalid Input. Cannot open.";

  // Read entire file
  t_fortran_record header;
  rcmf.seekg (0, std::ios::end);
  header.len = rcmf.tellg();
  rcmf.seekg (0, std::ios::beg);

  // Check for minimal size for reading in basic informations
  size_t minsize = 3*sizeof(int);
  if (doseq) minsize += sizeof(int);

  if (header.len <= minsize)
  {
    throw "OUT_HEAD file size is too short.";
  }

  header.data = new char[header.len];
  rcmf.read(header.data, header.len);
  rcmf.close();

  char *buf = header.data;
  if (doseq) buf = buf+sizeof(int);
  int iy = intvalfrombuf(buf); buf = buf + sizeof(int);
  int nsg = iy/h.iy;
  s.nx = (h.iy-2)*nsg;
  int jx = intvalfrombuf(buf); buf = buf + sizeof(int);
  s.ny = (h.jx-2)*nsg;
  // We do not need any vertical information
  int kz = intvalfrombuf(buf); buf = buf + sizeof(int);

  // Now I have record len of fortran I/O
  int nvals = iy*jx;
  int rdc = s.nx*s.ny;
  size_t size2D = nvals*sizeof(float);

  // Double check file size if correct before memcopy
  // We need to read only up to mask
  size_t chk = size2D*14;
  if (doseq) chk += 28*sizeof(int);
  if (chk < header.len)
  {
    delete [] header.data;
    throw "ERROR IN READING DOMAIN INFO DATA";
  }
  // We will take projection informations from here
  s.ds = rvalfrombuf(buf); buf = buf + sizeof(float);
  s.dsmeters = s.ds;
  s.clat = rvalfrombuf(buf); buf = buf + sizeof(float);
  s.clon = rvalfrombuf(buf); buf = buf + sizeof(float);
  s.xplat = rvalfrombuf(buf); buf = buf + sizeof(float);
  s.xplon = rvalfrombuf(buf); buf = buf + sizeof(float);
  // Skip grdfac
  buf = buf + sizeof(float);
  s.proj[6] = 0;
  memcpy(s.proj, buf, 6); buf = buf + 6;
  // Skip siglevs
  buf = buf + (kz+1)*sizeof(float);
  // Skip ptop
  buf = buf + sizeof(float);
  // Skip igrads, ibigend...
  buf = buf + 2*sizeof(int);
  // ... but read truelats.
  s.trlat1 = rvalfrombuf(buf); buf = buf + sizeof(float);
  s.trlat2 = rvalfrombuf(buf); buf = buf + sizeof(float);

  // First record is mostly empty except above infos
  buf = header.data + nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  float *space = new float[nvals];

  // Read terrain elevation
  vectorfrombuf(buf, space, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);
  if (s.ht) delete [] s.ht;
  s.ht = new float[rdc];
  reduce(space, s.ht, iy, jx, nsg);

  // Read terrain elevation std dev
  vectorfrombuf(buf, space, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);
  if (s.htsd) delete [] s.htsd;
  s.htsd = new float[rdc];
  reduce(space, s.htsd, iy, jx, nsg);

  // Read surface landuse
  vectorfrombuf(buf, space, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);
  if (s.landuse) delete [] s.landuse;
  s.landuse = new float[rdc];
  reduce(space, s.landuse, iy, jx, nsg);

  // Read latitude cross points
  vectorfrombuf(buf, space, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);
  if (s.xlat) delete [] s.xlat;
  s.xlat = new float[rdc];
  reduce(space, s.xlat, iy, jx, nsg);

  // Read longitude cross points
  vectorfrombuf(buf, space, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);
  if (s.xlon) delete [] s.xlon;
  s.xlon = new float[rdc];
  reduce(space, s.xlon, iy, jx, nsg);

  // Skip dlat, dlon
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  // Read xmap
  vectorfrombuf(buf, space, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);
  if (s.xmap) delete [] s.xmap;
  s.xmap = new float[rdc];
  reduce(space, s.xmap, iy, jx, nsg);

  // Skip dmap
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  // Read Coriolis force
  vectorfrombuf(buf, space, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);
  if (s.coriol) delete [] s.coriol;
  s.coriol = new float[rdc];
  reduce(space, s.coriol, iy, jx, nsg);

  // Skip snowam
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  // Read Land/Sea mask
  vectorfrombuf(buf, space, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);
  if (s.mask) delete [] s.mask;
  s.mask = new float[rdc];
  reduce(space, s.mask, iy, jx, nsg);

  delete [] header.data;
  delete [] space;
  return;
}

void rcmio::read_header(header_data &h)
{
  char fname[PATH_MAX];
  sprintf(fname, "%s%s%s", outdir, separator, "OUT_HEAD");

  std::ifstream rcmf;
  rcmf.open(fname, std::ios::binary);
  if (! rcmf.good()) throw "Invalid Input. Cannot open.";

  // Read entire file
  t_fortran_record header;
  rcmf.seekg (0, std::ios::end);
  header.len = rcmf.tellg();
  rcmf.seekg (0, std::ios::beg);

  // Check for minimal size for reading in basic informations
  size_t minsize = 8*sizeof(int);
  if (doseq) minsize += sizeof(int);

  if (header.len <= minsize)
  {
    throw "OUT_HEAD file size is too short.";
  }

  header.data = new char[header.len];
  rcmf.read(header.data, header.len);
  rcmf.close();

  char *buf = header.data;
  if (doseq) buf = buf+sizeof(int);
  h.mdate0 = intvalfrombuf(buf); buf = buf + sizeof(int);
  if (h.mdate0 != h.idate0)
  {
    std::cerr << "idate0 in input regcm is : " << h.idate0 << std::endl;
    std::cerr << "mdate0 in " << fname << " is : " << h.mdate0 << std::endl;
    throw "Regcm input file is not relative to OUT_HEAD file in output dir";
  }
  h.ibltyp = (short int) intvalfrombuf(buf); buf = buf + sizeof(int);
  h.icup = (short int) intvalfrombuf(buf); buf = buf + sizeof(int);
  h.imoist = (short int) intvalfrombuf(buf); buf = buf + sizeof(int);
  h.iboudy = (short int) intvalfrombuf(buf); buf = buf + sizeof(int);
  h.iy = intvalfrombuf(buf); buf = buf + sizeof(int);
  h.nx = h.iy-2;
  h.jx = intvalfrombuf(buf); buf = buf + sizeof(int);
  h.ny = h.jx-2;
  h.kz = intvalfrombuf(buf); buf = buf + sizeof(int);
  h.nz = h.kz;

  // Now I have record len of fortran I/O
  int nvals = (h.nx)*(h.ny);

  // Double check file size if correct before memcopy
  size_t chk = nvals*11*sizeof(float);
  if (doseq) chk += 22*sizeof(int);
  if (chk > header.len)
  {
    std::cerr << "CHK : " << chk << std::endl;
    std::cerr << "EXP : " << header.len << std::endl;
    std::cerr << "JX  : " << h.jx << std::endl;
    std::cerr << "IY  : " << h.iy << std::endl;
    delete [] header.data;
    throw "ERROR IN READING HEADER DATA";
  }

  // Calculate and store half sigma levels
  if (h.hsigf) delete [] h.hsigf;
  if (h.hsigm) delete [] h.hsigm;
  h.hsigf = new float[h.kz+1];
  h.hsigm = new float[h.kz];
  for (unsigned int i = 0; i < h.kz+1; i ++)
  {
    h.hsigf[i] = rvalfrombuf(buf); buf = buf + sizeof(float);
  }
  for (unsigned int i = 0; i < h.kz; i ++)
    h.hsigm[i] = h.hsigf[i]+0.5f*(h.hsigf[i+1]-h.hsigf[i]);
  h.ds = rvalfrombuf(buf); buf = buf + sizeof(float);
  h.dsmeters = h.ds*1000.0;
  h.ptop = rvalfrombuf(buf); buf = buf + sizeof(float);
  h.ptop *= 10.0; // Put in hPa from cbar
  h.clat = rvalfrombuf(buf); buf = buf + sizeof(float);
  h.clon = rvalfrombuf(buf); buf = buf + sizeof(float);
  h.xplat = rvalfrombuf(buf); buf = buf + sizeof(float);
  h.xplon = rvalfrombuf(buf); buf = buf + sizeof(float);
  h.proj[6] = 0;
  memcpy(h.proj, buf, 6); buf = buf + 6;
  h.dto = rvalfrombuf(buf); buf = buf + sizeof(float);
  h.dtb = rvalfrombuf(buf); buf = buf + sizeof(float);
  h.dtr = rvalfrombuf(buf); buf = buf + sizeof(float);
  h.dtc = rvalfrombuf(buf); buf = buf + sizeof(float);
  h.idirect = (unsigned short) intvalfrombuf(buf); buf = buf + sizeof(int);
  h.trlat1 = rvalfrombuf(buf); buf = buf + sizeof(float);
  h.trlat2 = rvalfrombuf(buf); buf = buf + sizeof(float);

  // First record is mostly empty except above infos
  buf = header.data + nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  // Read terrain elevation
  if (h.ht) delete [] h.ht;
  h.ht = new float[nvals];
  vectorfrombuf(buf, h.ht, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  // Read terrain elevation std dev
  if (h.htsd) delete [] h.htsd;
  h.htsd = new float[nvals];
  vectorfrombuf(buf, h.htsd, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  // Read vegetation type
  if (h.veg2d) delete [] h.veg2d;
  h.veg2d = new float[nvals];
  vectorfrombuf(buf, h.veg2d, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  // Read surface landuse
  if (h.landuse) delete [] h.landuse;
  h.landuse = new float[nvals];
  vectorfrombuf(buf, h.landuse, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  // Read latitude cross points
  if (h.xlat) delete [] h.xlat;
  h.xlat = new float[nvals];
  vectorfrombuf(buf, h.xlat, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  // Read longitude cross points
  if (h.xlon) delete [] h.xlon;
  h.xlon = new float[nvals];
  vectorfrombuf(buf, h.xlon, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  // Read map factor cross points
  if (h.xmap) delete [] h.xmap;
  h.xmap = new float[nvals];
  vectorfrombuf(buf, h.xmap, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  // Read map factor dot points
  if (h.dmap) delete [] h.dmap;
  h.dmap = new float[nvals];
  vectorfrombuf(buf, h.dmap, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  // Read Coriolis force
  if (h.coriol) delete [] h.coriol;
  h.coriol = new float[nvals];
  vectorfrombuf(buf, h.coriol, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  // Read Land/Sea mask
  if (h.mask) delete [] h.mask;
  h.mask = new float[nvals];
  vectorfrombuf(buf, h.mask, nvals);
  buf += nvals*sizeof(float);
  if (doseq) buf = buf + 2*sizeof(int);

  delete [] header.data;

  // Check what output is present for that header
  sprintf(fname, "%s%sATM.%d", outdir, separator, h.idate1);
  has_atm = fexist(fname);
  sprintf(fname, "%s%sSRF.%d", outdir, separator, h.idate1);
  has_srf = fexist(fname);
  sprintf(fname, "%s%sCHE.%d", outdir, separator, h.idate1);
  has_che = fexist(fname);
  sprintf(fname, "%s%sRAD.%d", outdir, separator, h.idate1);
  has_rad = fexist(fname);
  sprintf(fname, "%s%sSUB.%d", outdir, separator, h.idate1);
  has_sub = fexist(fname);
}

int rcmio::atmo_read_tstep(atmodata &a)
{
  if (! has_atm) return 1;
  if (! initatm)
  {
    if (doseq)
      readsize = a.datasize+(a.n3D+a.n2D)*2*sizeof(int);
    else
      readsize = a.datasize;
    char fname[PATH_MAX];
    sprintf(fname, "%s%sATM.%d", outdir, separator, a.rdate);
    atmf.open(fname, std::ios::binary);
    if (! atmf.good()) return -1;
    storage = new char[readsize];
    atmf.seekg (0, std::ios::end);
    atmsize = atmf.tellg();
    rcmdate df(a.rdate);
    int ss = df.datediffh(a.date0)/a.dt;
    atmf.seekg (ss*readsize, std::ios::beg);
    a.rdate = a.date0;
    initatm = true;
  }
  size_t pos = atmf.tellg( );
  if (pos+readsize > atmsize)
  {
    if (a.nfiles > 0)
    {
      atmf.close( );
      // Add one month
      a.rdate = nextmonth(a.rdate);
      // Open BC file
      char fname[PATH_MAX];
      sprintf(fname, "%s%sATM.%d", outdir, separator, a.rdate);
      atmf.open(fname, std::ios::binary);
      if (! atmf.good()) return -1;

      atmf.seekg (0, std::ios::end);
      atmsize = atmf.tellg();
      atmf.seekg (readsize, std::ios::beg);

      if (atmsize < readsize)
      {
        a.nfiles = 0;
        delete [] storage;
        storage = 0;
        atmf.close();
        return 1;
      }
      a.nfiles--;
    }
    else
    {
      delete [] storage;
      storage = 0;
      atmf.close();
      return 1;
    }
  }
  if (a.nsteps < 0)
  {
    delete [] storage;
    storage = 0;
    atmf.close();
    return 1;
  }
  atmf.read(storage, readsize);
  a.nsteps --;
  if (doseq)
  {
    char *p1 = storage;
    char *p2 = storage;
    for (int i = 0; i < a.n3D; i ++)
    {
      p1 += sizeof(int);
      memcpy(p2, p1, a.size3D*sizeof(float));
      p2 += a.size3D*sizeof(float);
      p1 += a.size3D*sizeof(float)+sizeof(int);
    }
    for (int i = 0; i < a.n2D; i ++)
    {
      p1 += sizeof(int);
      memcpy(p2, p1, a.size2D*sizeof(float));
      p2 += a.size2D*sizeof(float);
      p1 += a.size2D*sizeof(float)+sizeof(int);
    }
  }
  vectorfrombuf(storage, (float *) a.buffer, a.nvals);
  return 0;
}

int rcmio::srf_read_tstep(srfdata &s)
{
  if (! has_srf) return 1;
  if (! initsrf)
  {
    if (doseq)
      readsize = s.datasize+s.n2D*2*sizeof(int);
    else
      readsize = s.datasize;
    char fname[PATH_MAX];
    sprintf(fname, "%s%sSRF.%d", outdir, separator, s.rdate);
    srff.open(fname, std::ios::binary);
    if (! srff.good()) return -1;
    storage = new char[readsize];
    srff.seekg (0, std::ios::end);
    srfsize = srff.tellg();
    rcmdate df(s.rdate);
    int ss = df.datediffh(s.date0)/s.dt;
    srff.seekg (ss*readsize, std::ios::beg);
    s.rdate = s.date0;
    initsrf = true;
  }
  size_t pos = srff.tellg( );
  if (pos+readsize > srfsize)
  {
    if (s.nfiles > 0)
    {
      srff.close( );
      // Add one month
      s.rdate = nextmonth(s.rdate);
      // Open BC file
      char fname[PATH_MAX];
      sprintf(fname, "%s%sSRF.%d", outdir, separator, s.rdate);
      srff.open(fname, std::ios::binary);
      if (! srff.good()) return -1;

      srff.seekg (0, std::ios::end);
      srfsize = srff.tellg();
      srff.seekg (readsize, std::ios::beg);

      if (srfsize < readsize)
      {
        s.nfiles = 0;
        delete [] storage;
        storage = 0;
        srff.close();
        return 1;
      }
      s.nfiles--;
    }
    else
    {
      delete [] storage;
      storage = 0;
      srff.close();
      return 1;
    }
  }
  if (s.nsteps < 0)
  {
    delete [] storage;
    storage = 0;
    srff.close();
    return 1;
  }
  srff.read(storage, readsize);
  s.nsteps --;
  if (doseq)
  {
    char *p1 = storage;
    char *p2 = storage;
    for (int i = 0; i < s.n2D; i ++)
    {
      p1 += sizeof(int);
      memcpy(p2, p1, s.size2D*sizeof(float));
      p2 += s.size2D*sizeof(float);
      p1 += s.size2D*sizeof(float)+sizeof(int);
    }
  }
  vectorfrombuf(storage, (float *) s.buffer, s.nvals);
  return 0;
}

int rcmio::sub_read_tstep(subdata &u)
{
  if (! has_sub) return 1;
  if (! initsub)
  {
    if (doseq)
      readsize = u.datasize+u.n2D*2*sizeof(int)+3*sizeof(int);
    else
      readsize = u.datasize;
    char fname[PATH_MAX];
    sprintf(fname, "%s%sSUB.%d", outdir, separator, u.rdate);
    subf.open(fname, std::ios::binary);
    if (! subf.good()) return -1;
    storage = new char[readsize];
    subf.seekg (0, std::ios::end);
    subsize = subf.tellg();
    rcmdate df(u.rdate);
    int ss = df.datediffh(u.date0)/u.dt;
    subf.seekg (ss*readsize, std::ios::beg);
    u.rdate = u.date0;
    initsub = true;
  }
  size_t pos = subf.tellg( );
  if (pos+readsize > subsize)
  {
    if (u.nfiles > 0)
    {
      subf.close( );
      // Add one month
      u.rdate = nextmonth(u.rdate);
      // Open BC file
      char fname[PATH_MAX];
      sprintf(fname, "%s%sSUB.%d", outdir, separator, u.rdate);
      subf.open(fname, std::ios::binary);
      if (! subf.good()) return -1;

      subf.seekg (0, std::ios::end);
      subsize = subf.tellg();
      subf.seekg (readsize, std::ios::beg);

      if (subsize < readsize)
      {
        u.nfiles = 0;
        delete [] storage;
        storage = 0;
        subf.close();
        return 1;
      }
      u.nfiles--;
    }
    else
    {
      delete [] storage;
      storage = 0;
      subf.close();
      return 1;
    }
  }
  if (u.nsteps < 0)
  {
    delete [] storage;
    storage = 0;
    subf.close();
    return 1;
  }
  subf.read(storage, readsize);
  u.nsteps --;
  if (doseq)
  {
    char *p1 = storage+3*sizeof(int);
    char *p2 = storage+3*sizeof(int);
    for (int i = 0; i < u.n2D; i ++)
    {
      p1 += sizeof(int);
      memcpy(p2, p1, u.size2D*sizeof(float));
      p2 += u.size2D*sizeof(float);
      p1 += u.size2D*sizeof(float)+sizeof(int);
    }
  }
  vectorfrombuf(storage, (float *) u.buffer, u.nvals);
  return 0;
}

int rcmio::rad_read_tstep(raddata &r)
{
  if (! has_rad) return 1;
  if (! initrad)
  {
    if (doseq)
      readsize = r.datasize+(r.n3D+r.n2D)*2*sizeof(int);
    else
      readsize = r.datasize;

    char fname[PATH_MAX];
    sprintf(fname, "%s%sRAD.%d", outdir, separator, r.rdate);
    radf.open(fname, std::ios::binary);
    if (! radf.good()) return -1;
    storage = new char[readsize];
    radf.seekg (0, std::ios::end);
    radsize = radf.tellg();
    rcmdate df(r.rdate);
    int ss = df.datediffh(r.date0)/r.dt;
    radf.seekg (ss*readsize, std::ios::beg);
    r.rdate = r.date0;
    initrad = true;
  }
  size_t pos = radf.tellg( );
  if (pos+readsize > radsize)
  {
    if (r.nfiles > 0)
    {
      radf.close( );
      // Add one month
      r.rdate = nextmonth(r.rdate);
      // Open BC file
      char fname[PATH_MAX];
      sprintf(fname, "%s%sRAD.%d", outdir, separator, r.rdate);
      radf.open(fname, std::ios::binary);
      if (! radf.good()) return -1;

      radf.seekg (0, std::ios::end);
      radsize = radf.tellg();
      radf.seekg (readsize, std::ios::beg);

      if (radsize < readsize)
      {
        r.nfiles = 0;
        delete [] storage;
        storage = 0;
        radf.close();
        return 1;
      }
      r.nfiles--;
    }
    else
    {
      delete [] storage;
      storage = 0;
      radf.close();
      return 1;
    }
  }
  if (r.nsteps < 0)
  {
    delete [] storage;
    storage = 0;
    radf.close();
    return 1;
  }
  radf.read(storage, readsize);
  r.nsteps --;
  if (doseq)
  {
    char *p1 = storage;
    char *p2 = storage;
    for (int i = 0; i < r.n3D; i ++)
    {
      p1 += sizeof(int);
      memcpy(p2, p1, r.size3D*sizeof(float));
      p2 += r.size3D*sizeof(float);
      p1 += r.size3D*sizeof(float)+sizeof(int);
    }
    for (int i = 0; i < r.n2D; i ++)
    {
      p1 += sizeof(int);
      memcpy(p2, p1, r.size2D*sizeof(float));
      p2 += r.size2D*sizeof(float);
      p1 += r.size2D*sizeof(float)+sizeof(int);
    }
  }
  vectorfrombuf(storage, (float *) r.buffer, r.nvals);
  return 0;
}

int rcmio::che_read_tstep(chedata &c)
{
  if (! has_che) return 1;
  if (! initche)
  {
    if (doseq)
      readsize = c.datasize+(c.n3D+c.n2D)*2*sizeof(int);
    else
      readsize = c.datasize;
    char fname[PATH_MAX];
    sprintf(fname, "%s%sCHE.%d", outdir, separator, c.rdate);
    chef.open(fname, std::ios::binary);
    if (! chef.good()) return -1;
    storage = new char[readsize];
    chef.seekg (0, std::ios::end);
    chesize = chef.tellg();
    rcmdate df(c.rdate);
    int ss = df.datediffh(c.date0)/c.dt;
    chef.seekg (ss*readsize, std::ios::beg);
    c.rdate = c.date0;
    initche = true;
  }
  size_t pos = chef.tellg( );
  if (pos+readsize > chesize)
  {
    if (c.nfiles > 0)
    {
      chef.close( );
      // Add one month
      c.rdate = nextmonth(c.rdate);
      // Open BC file
      char fname[PATH_MAX];
      sprintf(fname, "%s%sCHE.%d", outdir, separator, c.rdate);
      chef.open(fname, std::ios::binary);
      if (! chef.good()) return -1;

      chef.seekg (0, std::ios::end);
      chesize = chef.tellg();
      chef.seekg (readsize, std::ios::beg);

      if (chesize < readsize)
      {
        c.nfiles = 0; 
        delete [] storage;
        storage = 0;
        chef.close();
        return 1;
      }
      c.nfiles--;
    }
    else
    {
      delete [] storage;
      storage = 0;
      chef.close();
      return 1;
    }
  }
  if (c.nsteps < 0)
  {
    delete [] storage;
    storage = 0;
    chef.close();
    return 1;
  }
  chef.read(storage, readsize);
  c.nsteps --;
  if (doseq)
  {
    char *p1 = storage;
    char *p2 = storage;
    for (int i = 0; i < c.n3D; i ++)
    {
      p1 += sizeof(int);
      memcpy(p2, p1, c.size3D*sizeof(float));
      p2 += c.size3D*sizeof(float);
      p1 += c.size3D*sizeof(float)+sizeof(int);
    }
    for (int i = 0; i < c.n2D; i ++)
    {
      p1 += sizeof(int);
      memcpy(p2, p1, c.size2D*sizeof(float));
      p2 += c.size2D*sizeof(float);
      p1 += c.size2D*sizeof(float)+sizeof(int);
    }
  }
  vectorfrombuf(storage, (float *) c.buffer, c.nvals);
  return 0;
}

int rcmio::find_nextstep(std::ifstream &f, size_t recsize, size_t end)
{
  char check[3*sizeof(int)];
  char tbc[3*sizeof(int)];

  // Goto begin, skip time integer
  f.seekg(0, std::ios::beg);
  if (doseq) f.seekg(sizeof(int), std::ios::cur);
  f.seekg(sizeof(int), std::ios::cur);

  // Read three integers (nx, ny, nz)
  f.read(check, 3*sizeof(int));
  // Goto next record
  f.seekg(recsize, std::ios::beg);
  if (doseq) f.seekg(2*sizeof(int), std::ios::cur);

  int irec = 1;
  size_t pos;
  // Start search
  while ((pos = f.tellg( )) < end)
  {
    f.seekg (sizeof(int), std::ios::cur);
    f.read(tbc, 3*sizeof(int));
    if (memcmp((const void *) tbc, (const void *) check, 3*sizeof(int)) == 0)
    {
      // Found again nx, ny, nz
      f.seekg (0, std::ios::beg);
      return irec;
    }
    // This is not a new time. Goto next
    f.seekg(pos, std::ios::beg);
    f.seekg(recsize, std::ios::cur);
    if (doseq) f.seekg (2*sizeof(int), std::ios::cur);
    irec ++;
  }
  f.seekg (0, std::ios::beg);
  return irec;
}

int rcmio::bc_read_tstep(bcdata &b)
{
  if (! initbc)
  {
    // Open BC file
    char fname[PATH_MAX];
    sprintf(fname, "%s%s%s_ICBC%d", outdir, separator, b.name, b.rdate);
    bcf.open(fname, std::ios::binary);
    if (! bcf.good()) return -1;

    bcf.seekg (0, std::ios::end);
    bcsize = bcf.tellg();
    bcf.seekg (0, std::ios::beg);

    b.n3D = 4;
    int baserec = b.n3D*b.nz+3;
    int baseso4 = baserec+1;
    int baseusg = baserec+13;
    int allrecs = baserec+14;

    int nrecs = find_nextstep(bcf, b.size2D*sizeof(float), bcsize);
    if (nrecs == baserec)
    {
      b.n2D = 3;
      b.ehso4 = false;
      b.usgs = false;
    }
    else if (nrecs == baseso4)
    {
      b.n2D = 4;
      b.ehso4 = true;
      b.usgs = false;
    }
    else if (nrecs == baseusg)
    {
      b.n2D = 16;
      b.ehso4 = false;
      b.usgs = true;
    }
    else if (nrecs == allrecs)
    {
      b.n2D = 17;
      b.ehso4 = true;
      b.usgs = true;
    } 
    else
    {
      std::cerr << "Found NREC = " << nrecs << std::endl;
      throw "Cannot find lsmtyp and dattyp";
    }
    b.nvals = b.n3D*b.size3D + (b.n2D-1)*b.size2D;
    b.datasize = b.nvals*sizeof(float);
    b.buffer = new char[b.datasize];
    readsize = b.datasize+b.size2D*sizeof(float);
    if (doseq)
      readsize = readsize+(b.n3D+b.n2D)*2*sizeof(int);
    storage = new char[readsize];

    b.u = (float *) b.buffer;
    b.v = b.u + b.size3D;
    b.t = b.v + b.size3D;
    b.q = b.t + b.size3D;
    b.px = b.q + b.size3D;
    b.ts = b.px + b.size2D;
    float *pnt = b.ts + b.size2D;
    if (b.ehso4)
    {
      b.so4 = pnt;
      pnt += b.size2D;
    }
    if (b.usgs)
    {
      b.sm = pnt;
      b.icet = b.sm + 4 * b.size2D;
      b.soilt = b.icet + 4 * b.size2D;
      b.snowd = b.icet + 4 * b.size2D;
    }
    rcmdate df(b.rdate);
    int ss = df.datediffh(b.date0)/b.dt;
    bcf.seekg (ss*readsize, std::ios::beg);
    initbc = true;
  }
  size_t pos = bcf.tellg( );
  if (pos+readsize > bcsize)
  {
    if (b.nfiles > 0)
    {
      bcf.close();
      // Add one month
      b.rdate = nextmonth(b.rdate);
      // Open BC file
      char fname[PATH_MAX];
      sprintf(fname, "%s%s%s_ICBC%d", outdir, separator, b.name, b.rdate);
      bcf.open(fname, std::ios::binary);
      if (! bcf.good()) return -1;

      bcf.seekg (0, std::ios::end);
      bcsize = bcf.tellg();
      bcf.seekg (readsize, std::ios::beg);

      b.nfiles--;
    }
    else
    {
      delete [] storage;
      storage = 0;
      bcf.close();
      return 1;
    }
  }
  if (b.nsteps < 0)
  {
    delete [] storage;
    storage = 0;
    bcf.close();
    return 1;
  }
  bcf.read(storage, readsize);
  b.nsteps --;
  if (doseq)
  {
    char *p1 = storage;
    char *p2 = storage;
    for (int i = 0; i < b.n3D; i ++)
    {
      p1 += sizeof(int);
      memcpy(p2, p1, b.size3D*sizeof(float));
      p2 += b.size3D*sizeof(float);
      p1 += b.size3D*sizeof(float)+sizeof(int);
    }
    for (int i = 0; i < b.n2D; i ++)
    {
      p1 += sizeof(int);
      memcpy(p2, p1, b.size2D*sizeof(float));
      p2 += b.size2D*sizeof(float);
      p1 += b.size2D*sizeof(float)+sizeof(int);
    }
  }
  b.date1 = intvalfrombuf(storage);
  vectorfrombuf(storage+b.size2D*sizeof(float), (float *) b.buffer, b.nvals);
  return 0;
}

rcmdate::rcmdate(int date)
{
  struct tm ref;
  memset(&ref, 0, sizeof(struct tm));
  unsigned int base = date;
  basey = base/1000000;
  base = base-basey*1000000;
  ref.tm_year = basey - 1900;
  basem = base/10000;
  base = base-basem*10000;
  ref.tm_mon = basem-1;
  based = base/100;
  base = base-based*100;
  ref.tm_mday = based-1;
  baseh = base;
  ref.tm_hour = baseh;
  tdate = mktime(&ref);
}

static bool isleap(int year)
{
  if((year % 4 == 0 && year % 100 != 0) || year % 400 == 0)
    return true;
  else
    return false;
}

static int mdays(int year, int mon)
{
  if (mon != 2) return mlen[mon];
  return (mlen[2] + (isleap(year) ? 1 : 0));
}

static int idatei(int y, int m, int d, int h)
{
  return y*1000000+m*10000+d*100+h;
}

int rcmdate::idate( )
{
  return idatei(basey, basem, based, baseh);
}

int rcmdate::idatendh( )
{
  return basey*1000000+basem*10000;
}

int rcmdate::hour_adder(int nh)
{
  int nmd;
  int nsum  = nh / 23;
  int nlast = nh % 23;
  for (int n = 0; n < nsum; n ++)
  {
    baseh = baseh + 23;
    if (baseh > 23)
    {
      based = based + 1;
      baseh = baseh - 24;
      nmd = mdays(basey, basem);
      if (based > nmd)
      {
        basem = basem + 1;
        based = based - nmd;
        if (basem > 12)
        {
          basey = basey + 1;
          basem = basem - 12;
        }
      }
    }
  }
  baseh = baseh+nlast;
  if (baseh > 23)
  {
    based = based + 1;
    baseh = baseh - 24;
    nmd = mdays(basey, basem);
    if (based > nmd)
    {
      basem = basem + 1;
      based = based - nmd;
      if (basem > 12)
      {
        basey = basey + 1;
        basem = basem - 12;
      }
    }
  }
  return idatei(basey, basem, based, baseh);
}

unsigned long rcmdate::datediffh(int date1)
{
  struct tm ref;
  memset(&ref, 0, sizeof(struct tm));
  unsigned int base = date1;
  int basey1 = base/1000000;
  base = base-basey1*1000000;
  ref.tm_year = basey1 - 1900;
  int basem1 = base/10000;
  base = base-basem1*10000;
  ref.tm_mon = basem1-1;
  int based1 = base/100;
  base = base-based1*100;
  ref.tm_mday = based1-1;
  ref.tm_hour = base;
  time_t tdate1 = mktime(&ref);
  return ((unsigned long) (difftime(tdate1,tdate)/3600.0));
}

unsigned long rcmdate::unixtime( )
{
  struct tm ref;
  memset(&ref, 0, sizeof(struct tm));
  ref.tm_year = 70;
  time_t tzero = mktime(&ref);
  return ((unsigned long) (difftime(tdate,tzero)/3600.0));
}

rcmio::~rcmio( )
{
  if (storage) delete [ ] storage;
}
