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
  size_t size3D = nx*ny*nz;
  size_t size2D = nx*ny;
  nvals = 6*size3D + 5*size2D;
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

rcmio::rcmio(char *directory, bool lbig, bool ldirect)
{
  // big or little endian swapping ?
  if ((littlearch() && lbig) || (bigarch() && ! lbig)) doswap = true;
  else doswap = false;

  // direct access or sequential file ?
  doseq = ! ldirect;
  strncpy(outdir, directory, PATH_MAX);

  // All to false
  initatmo = false;
  has_atmo = false;
  initchem = false;
  has_chem = false;
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
  has_atmo = fexist(fname);
  sprintf(fname, "%s%sSRF.%d", outdir, separator, h.mdate0);
  has_srf = fexist(fname);
  sprintf(fname, "%s%sSUB.%d", outdir, separator, h.mdate0);
  has_sub = fexist(fname);
  sprintf(fname, "%s%sCHE.%d", outdir, separator, h.mdate0);
  has_chem = fexist(fname);
  sprintf(fname, "%s%sRAD.%d", outdir, separator, h.mdate0);
  has_rad = fexist(fname);
}

int rcmio::atmo_read_tstep(atmodata &a)
{
  if (! has_atmo) return 1;
  if (! initatmo)
  {
    char fname[PATH_MAX];
    sprintf(fname, "%s%sATM.%d", outdir, separator, a.date0);
    atmof.open(fname, std::ios::binary);
    if (! atmof.good()) return -1;
    storage = new char[a.datasize];
    initatmo = true;
    atmof.seekg (0, std::ios::end);
    atmsize = atmof.tellg();
    atmof.seekg (0, std::ios::beg);
  }
  if (atmof.tellg( ) == atmsize)
  {
    delete [] storage;
    storage = 0;
    atmof.close();
    return 1;
  }
  atmof.read(storage, a.datasize);
  vectorfrombuf(storage, (float *) a.buffer, a.nvals);
  return 0;
}

rcmio::~rcmio( )
{
  if (storage) delete [ ] storage;
}
