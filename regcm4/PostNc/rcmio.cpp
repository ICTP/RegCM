/***************************************************************************
 *   Copyright (C) 2008-2009 Graziano Giuliani                             *
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

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

using namespace rcm;

static const int aword = 4;
static const char separator[2] = "/";

bool fexist(char *fname)
{
  struct stat finfo;
  int retval = stat(fname, &finfo);
  if (retval == 0) return true;
  return false;
}

typedef struct {
  size_t len;
  char *data;
} t_fortran_record;

header_data::header_data( )
{
  hsigf = ht = htsd = veg2d = landuse = 0;
  xlat = xlon = xmap = dmap = coriol = mask = 0;
}

header_data::~header_data( )
{
  free_space( );
}

void header_data::free_space( )
{
  if (hsigf) delete [ ] hsigf;
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
  hsigf = ht = htsd = veg2d = landuse = 0;
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

atmodata::atmodata(int nx, int ny, int nz, int mdate0, float dto)
{
  date0 = mdate0;
  dt = dto;
  n3D = 6;
  n2D = 5;
  size3D = nx*ny*nz;
  size2D = nx*ny;
  nvals = n3D*size3D + n2D*size2D;
  datasize = nvals*sizeof(float);
  buffer = new char[datasize];
  u = (float *) buffer;
  v = u + size3D;
  omega = v + size3D;
  t = omega + size3D;
  qv = t + size3D;
  qc = qv + size3D;
  psa = qc + size3D;
  tpr = psa + size2D;
  tgb = tpr + size2D;
  swt = tgb + size2D;
  rno = swt + size2D;
}

atmodata::~atmodata( )
{
  delete [] buffer;
}

raddata::raddata(int nx, int ny, int nz, int mdate0, float dtr)
{
  date0 = mdate0;
  dt = dtr;
  n3D = 4;
  n2D = 9;
  size3D = nx*ny*nz;
  size2D = nx*ny;
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
}

raddata::~raddata( )
{
  delete [] buffer;
}

chedata::chedata(int nx, int ny, int nz, int mdate0, float dtr)
{
  date0 = mdate0;
  dt = dtr;
  ntr = 10;
  n3D = 13;
  n2D = 72;
  size3D = nx*ny*nz;
  size2D = nx*ny;
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
}

chedata::~chedata( )
{
  delete [] buffer;
}

srfdata::srfdata(int nx, int ny, int mdate0, float dto)
{
  date0 = mdate0;
  dt = dto;
  n2D = 27;
  size2D = nx*ny;
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

rcmio::rcmio(char *directory, bool lbig, bool ldirect)
{
  // big or little endian swapping ?
  if ((littlearch() && lbig) || (bigarch() && ! lbig)) doswap = true;
  else doswap = false;

  // direct access or sequential file ?
  doseq = ! ldirect;
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
  storage = 0;
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
  if (chk != header.len)
  {
    delete [] header.data;
    throw "ERROR IN READING HEADER DATA";
  }

  // Calculate and store half sigma levels
  if (h.hsigf) delete [] h.hsigf;
  h.hsigf = new float[h.kz];
  float *fsigf = new float[h.kz+1];
  for (int i = 0; i < h.kz+1; i ++)
  {
    fsigf[h.kz-i] = rvalfrombuf(buf); buf = buf + sizeof(float);
  }
  for (int i = 0; i < h.kz; i ++)
    h.hsigf[i] = fsigf[i]+0.5*(fsigf[i+1]-fsigf[i]);
  delete [] fsigf;
  h.ds = rvalfrombuf(buf); buf = buf + sizeof(float);
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
  h.idirect = intvalfrombuf(buf); buf = buf + sizeof(int);

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

  // Last check
  if ((buf-header.data) != header.len)
  {
    h.free_space();
    delete [] header.data;
    throw "READ error in file: size mismatch";
  }

  delete [] header.data;

  // Check what output is present for that header
  sprintf(fname, "%s%sATM.%d", outdir, separator, h.mdate0);
  has_atm = fexist(fname);
  sprintf(fname, "%s%sSRF.%d", outdir, separator, h.mdate0);
  has_srf = fexist(fname);
  sprintf(fname, "%s%sSUB.%d", outdir, separator, h.mdate0);
  has_sub = fexist(fname);
  sprintf(fname, "%s%sCHE.%d", outdir, separator, h.mdate0);
  has_che = fexist(fname);
  sprintf(fname, "%s%sRAD.%d", outdir, separator, h.mdate0);
  has_rad = fexist(fname);
}

int rcmio::atmo_read_tstep(atmodata &a)
{
  if (! has_atm) return 1;
  if (! initatm)
  {
    char fname[PATH_MAX];
    sprintf(fname, "%s%sATM.%d", outdir, separator, a.date0);
    atmf.open(fname, std::ios::binary);
    if (! atmf.good()) return -1;
    if (doseq)
      storage = new char[a.datasize+(a.n3D+a.n2D)*2*sizeof(int)];
    else
      storage = new char[a.datasize];
    initatm = true;
    atmf.seekg (0, std::ios::end);
    atmsize = atmf.tellg();
    atmf.seekg (0, std::ios::beg);
  }
  if (atmf.tellg( ) == atmsize)
  {
    delete [] storage;
    storage = 0;
    atmf.close();
    return 1;
  }
  if (doseq)
  {
    atmf.read(storage, a.datasize+(a.n3D+a.n2D)*2*sizeof(int));
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
  else
    atmf.read(storage, a.datasize);
  vectorfrombuf(storage, (float *) a.buffer, a.nvals);
  return 0;
}

int rcmio::srf_read_tstep(srfdata &s)
{
  if (! has_srf) return 1;
  if (! initsrf)
  {
    char fname[PATH_MAX];
    sprintf(fname, "%s%sSRF.%d", outdir, separator, s.date0);
    srff.open(fname, std::ios::binary);
    if (! srff.good()) return -1;
    if (doseq)
      storage = new char[s.datasize+s.n2D*2*sizeof(int)];
    else
      storage = new char[s.datasize];
    initsrf = true;
    srff.seekg (0, std::ios::end);
    srfsize = srff.tellg();
    srff.seekg (0, std::ios::beg);
  }
  if (srff.tellg( ) == srfsize)
  {
    delete [] storage;
    storage = 0;
    srff.close();
    return 1;
  }
  if (doseq)
  {
    srff.read(storage, s.datasize+s.n2D*2*sizeof(int));
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
  else
    srff.read(storage, s.datasize);
  vectorfrombuf(storage, (float *) s.buffer, s.nvals);
  return 0;
}

int rcmio::rad_read_tstep(raddata &r)
{
  if (! has_rad) return 1;
  if (! initrad)
  {
    char fname[PATH_MAX];
    sprintf(fname, "%s%sRAD.%d", outdir, separator, r.date0);
    radf.open(fname, std::ios::binary);
    if (! radf.good()) return -1;
    if (doseq)
      storage = new char[r.datasize+(r.n3D+r.n2D)*2*sizeof(int)];
    else
      storage = new char[r.datasize];
    initrad = true;
    radf.seekg (0, std::ios::end);
    radsize = radf.tellg();
    radf.seekg (0, std::ios::beg);
  }
  if (radf.tellg( ) == radsize)
  {
    delete [] storage;
    storage = 0;
    radf.close();
    return 1;
  }
  if (doseq)
  {
    radf.read(storage, r.datasize+(r.n3D+r.n2D)*2*sizeof(int));
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
  else
    radf.read(storage, r.datasize);
  vectorfrombuf(storage, (float *) r.buffer, r.nvals);
  return 0;
}

int rcmio::che_read_tstep(chedata &c)
{
  if (! has_che) return 1;
  if (! initche)
  {
    char fname[PATH_MAX];
    sprintf(fname, "%s%sCHE.%d", outdir, separator, c.date0);
    chef.open(fname, std::ios::binary);
    if (! chef.good()) return -1;
    if (doseq)
      storage = new char[c.datasize+(c.n3D+c.n2D)*2*sizeof(int)];
    else
      storage = new char[c.datasize];
    initche = true;
    chef.seekg (0, std::ios::end);
    chesize = chef.tellg();
    chef.seekg (0, std::ios::beg);
  }
  if (chef.tellg( ) == chesize)
  {
    delete [] storage;
    storage = 0;
    chef.close();
    return 1;
  }
  if (doseq)
  {
    chef.read(storage, c.datasize+(c.n3D+c.n2D)*2*sizeof(int));
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
  else
    chef.read(storage, c.datasize);
  vectorfrombuf(storage, (float *) c.buffer, c.nvals);
  return 0;
}

rcmio::~rcmio( )
{
  if (storage) delete [ ] storage;
}
