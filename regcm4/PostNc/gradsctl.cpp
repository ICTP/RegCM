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

#include <gradsctl.h>
#include <cmath>
#include <sstream>
#include <iostream>

using namespace rcm;

static void extrema(float *vals, int nvals, float *max, float *min)
{
  float fmax = *vals;
  float fmin = *vals;
  for (int i = 1; i < nvals; i ++)
  {
    if (*(vals+i) > fmax) fmax = *(vals+i);
    else if (*(vals+i) < fmin) fmin = *(vals+i);
  }
  *max = fmax;
  *min = fmin;
  return;
}

gradsvar::gradsvar( )
{
  vline = "unknown=>unknown 0 t,z,y,x unknown var (unkn)";
}

gradsvar::gradsvar(char *ncname, char *gname, char *desc, int nlevs, bool tv)
{
  set(ncname, gname, desc, nlevs, tv);
}

void gradsvar::set(char *ncname, char *gname, char *desc, int nlevs, bool tv)
{
  std::stringstream a(std::stringstream::out);
  a << ncname << "=>" << gname << " " << nlevs << " ";
  if (tv)
    a << "t,";
  if (nlevs > 0)
    a << "z,";
  a << "y,x " << desc;
  vline = a.str();
}

std::ostream& rcm::operator<< (std::ostream& os, const gradsvar &g)
{
  os << g.vline;
  return os;
}

gradsctl::gradsctl(char *ctlname, char *dname)
{
  ctlf.open(ctlname);
  ctlf << "dset ^" << dname << std::endl
       << "dtype netcdf" << std::endl;
}

void gradsctl::head(char *title, float missval)
{
  ctlf << "title " << title << std::endl
       << "undef " << missval << "_FillValue" << std::endl;
}

void gradsctl::set_grid(domain_data &d)
{
  float ci, cj;

  cj = ((float) d.nx )/2.0f;
  ci = ((float) d.ny )/2.0f;
  float maxlat, maxlon, minlat, minlon;
  extrema(d.xlat, d.nx*d.ny, &maxlat, &minlat);
  extrema(d.xlon, d.nx*d.ny, &maxlon, &minlon);
  float rlatinc , rloninc;
  rlatinc = rloninc = d.ds*0.001f/111.0f/2.0f;
  int nlat = 2 + rintf(fabs(maxlat-minlat)/rlatinc);
  int nlon = 2 + rintf(fabs((maxlon-minlon)/rlatinc));
  float dx = d.ds / 111000.0f;
  float dy = d.ds /111000.0f*0.95238f;
  if (strcmp(d.proj, "LAMCON") == 0)
  {
    ctlf << "pdef " << d.ny << " " << d.nx << " "
         << "lcc " << d.clat << " " << d.clon << " "
         << ci << " " << cj << " " << d.trlat1 << " " << d.trlat2 << " "
         << d.clon << " " << d.ds << " " << d.ds << std::endl
         << "xdef " << nlon << " linear " << minlon - rloninc << " "
         << rloninc << std::endl
         << "ydef " << nlat << " linear " << minlat - rlatinc << " "
         << rlatinc << std::endl;
  }
  else if (strcmp(d.proj, "POLSTR") == 0)
  {
    std::cout << "WARNING: Approximate projection. Not supported in GrADS"
              << std::endl;
    int cij = (d.nx/2)*d.ny+(d.ny/2);
    ctlf << "pdef " << d.ny << " " << d.nx << " "
         << "ops " << d.clat << " " << d.clon << " 0 0 "
         << ci-1 << " " << cj-1 << " " << d.ds << " " << d.ds << std::endl
         << "xdef " << nlon << " linear " << minlon - rloninc << " "
         << rloninc << std::endl
         << "ydef " << nlat << " linear " << minlat - rlatinc << " "
         << rlatinc << std::endl;;
  }
  else if (strcmp(d.proj, "NORMER") == 0)
  {
    ctlf << "xdef " << d.ny << " linear " << minlon << " "
         << (maxlon-minlon)/(d.ny+1) << std::endl
         << "ydef " << d.nx << " levels " << std::endl << "  ";
    int c = 1;
    for (int j = 0; j < d.nx; j ++, c ++)
    {
      ctlf << d.xlat[j*d.ny] << " ";
      if (c > 8 && j != d.nx-1) 
      {
        ctlf << std::endl << "  ";
        c = 0;
      }
    }
    ctlf << std::endl;    
  }
  else if (strcmp(d.proj, "ROTMER") == 0)
  {
    std::cout << "WARNING: Approximate projection. Not supported in GrADS"
              << std::endl;
    ctlf << "pdef " << d.ny << " " << d.nx << " "
         << "eta.u " << d.xplon << " " << d.xplat << " "
         << dx << " " << dy << std::endl
         << "xdef " << nlon << " linear " << minlon - rloninc << " "
         << rloninc << std::endl
         << "ydef " << nlat << " linear " << minlat - rlatinc << " "
         << rlatinc << std::endl;
  }
  else
    throw "Unknown projection used.";
  ntimes = 0;
  levline = "zdef 1 levels 0.0";
  timeline = "linear 00Z31dec1999 1yr";
}

void gradsctl::set_levs(float *lev, int nz)
{
  std::stringstream a(std::stringstream::out);
  a << "zdef " << nz << " levels ";
  for (int i = 0; i < nz; i ++)
    a << lev[i] << " ";
  levline = a.str( );
  return;
}

void gradsctl::set_time(char *timestr)
{
  timeline = "linear ";
  timeline += timestr;
}

void gradsctl::add_time( )
{
  ntimes++;
}

void gradsctl::add_var(gradsvar &g)
{
  vars.push_back(g);
}

void gradsctl::finalize( )
{
  ctlf << levline << std::endl;
  ctlf << "tdef " << ntimes << " " << timeline << std::endl;
  ctlf << "vars " << vars.size( ) << std::endl;
  int ncount = vars.size( );
  for (int i = 0; i < ncount; i ++)
  {
    ctlf << vars.front( ) << std::endl;
    vars.pop_front( );
  }
  ctlf << "endvars" << std::endl;
  return;
}

char *gradsctl::gradstime(int year, int month, int day, int hour, int incr)
{
  static char gratim[32];
  static char months[12][4] = { "jan", "feb", "mar", "apr", "may", "jun",
                                "jul", "aug", "sep", "oct", "nov", "dec" };

  snprintf(gratim, 32, "%02dZ%02d%s%04d %dhr",
            hour, day, months[month-1], year, incr);
  return gratim;
}

gradsctl::~gradsctl( )
{
  ctlf.close();
}
