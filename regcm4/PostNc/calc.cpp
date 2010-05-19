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

#include <calc.h>
#include <cmath>
#include <cstring>

using namespace rcm;

static const float misval = -1e34;
static const float tzero = 273.15;
static const float rgas = 287.0058;
static const float rgti = 1.0/9.80665;
static const float cpd = 1005.46;
static const float rovcp = rgas/cpd;
static const float ep2 = 0.62197;
static const float svp1 = 6.112;
static const float svp2 = 17.67;
static const float svp3 = 29.65;
static const float svp4 = svp1;
static const float svp5 = 22.514;
static const float svp6 = 6150.0;

inline float rhfromptq(float p, float t, float q)
{
  float r, satvp, qs;
  if (t > tzero)
    satvp = svp1*expf(svp2*(t-tzero)/(t-svp3));
  else
    satvp = svp4*expf(svp5-svp6/t);
  qs =  ep2*satvp/(p-satvp);
  r = (q/qs);
  if (r < 0.0) r = 0.0;
  if (r > 1.0) r = 1.0;
  return r;
}

srfcalc::srfcalc(header_data &h)
{
  nh = h.nx*h.ny;
  r2 = new float[nh];
}

srfcalc::~srfcalc()
{
  delete [] r2;
}

void srfcalc::calcrh(float *sp, float *t2, float *q2)
{
  for (int i = 0; i < nh; i ++)
  {
    if (sp[i] > 0.0)
      r2[i] = rhfromptq(sp[i], t2[i], q2[i]);
    else
      r2[i] = misval;
  }
  return;
}

void srfcalc::do_calc(srfdata &s, t_srf_deriv &d)
{
  calcrh(s.psb, s.t2m, s.q2m);
  d.r2 = r2;
  return;
}

subcalc::subcalc(subdom_data &s)
{
  nh = s.nx*s.ny;
  r2 = new float[nh];
}

subcalc::~subcalc()
{
  delete [] r2;
}

void subcalc::calcrh(float *sp, float *t2, float *q2)
{
  for (int i = 0; i < nh; i ++)
  {
    if (sp[i] > 0.0)
      r2[i] = rhfromptq(sp[i], t2[i], q2[i]);
    else
      r2[i] = misval;
  }
  return;
}

void subcalc::do_calc(subdata &s, t_srf_deriv &d)
{
  calcrh(s.psb, s.t2m, s.q2m);
  d.r2 = r2;
  return;
}

atmcalc::atmcalc(header_data &h)
{
  nk = h.nz;
  nx = h.nx;
  ny = h.ny;
  nh = nx*ny;
  ptop = h.ptop;
  hsigf = h.hsigf;
  hsigm = h.hsigm;
  zs = h.ht;
  ds = h.ds;
  ds2r = 1.0f/(2.0f*ds);
  xmap = h.xmap;
  dmap = h.dmap;
  p  = new float[nh*nk];
  rh = new float[nh*nk];
  td = new float[nh*nk];
  pt = new float[nh*nk];
  ht = new float[nh*nk];
  vr = new float[nh*nk];
  dv = new float[nh*nk];
}

atmcalc::~atmcalc( )
{
  delete [] p;
  delete [] rh;
  delete [] td;
  delete [] pt;
  delete [] ht;
  delete [] vr;
  delete [] dv;
}

void atmcalc::calcp(float *sp)
{
  for(int k = 0; k < nk; k ++)
   for (int i = 0; i < nh; i ++)
     p[k*nh+i] = (sp[i]-ptop)*hsigm[k]+ptop;
  return;
}

void atmcalc::calcrh(float *t, float *q)
{
  for (int i = 0; i < nh*nk; i ++)
    rh[i] = rhfromptq(p[i], t[i], q[i]);
  return;
}

void atmcalc::calctd(float *t)
{
  float rx, tx, dpd;

  for (int i = 0; i < nh*nk; i ++)
  {
    rx = 1.0f - rh[i];
    tx = t[i] - tzero;
    dpd = ((14.55f+0.144f*tx) * rx) +
          (2.0f * powf(((2.5f+0.007f*tx)*rx), 3.0f) +
          (15.9f+0.117f*tx) * powf(rx,14.0f));
    td[i] = t[i] - dpd;
  }
  return;
}

void atmcalc::calcpt(float *t)
{
  for (int i = 0; i < nh*nk; i ++)
    pt[i] = t[i]*powf((1000.0f/p[i]), rovcp);
  return;
}

void atmcalc::calcht(float *ps, float *t)
{
  for (int i = 0; i < nh; i ++)
    ht[i] = zs[i] + rgas*rgti*t[i]*logf(ps[i]/p[i]);
  float tbar;
  for(int k = 1; k < nk; k ++)
  {
    for (int i = 0; i < nh; i ++)
    {
      tbar = 0.5f*(t[k*nh+i]+t[(k-1)*nh+i]);
      ht[k*nh+i] = ht[(k-1)*nh+i] + 
                   rgas*rgti*tbar*logf(p[(k-1)*nh+i]/p[k*nh+i]);
    }
  }
  return;
}

void atmcalc::calcdv(float *u, float *v)
{
  float u1, u2, u3, u4, v1, v2, v3, v4;
  for(int k = 0; k < nk; k ++)
  {
    for (int j = 0; j < ny-1; j ++)
    {
      for (int i = 0; i < nx-1; i ++)
      {
        u1 = u[k*nh+i*ny+j]/dmap[i*ny+j];
        u2 = u[k*nh+i*ny+j+1]/dmap[i*ny+j+1];
        u3 = u[k*nh+(i+1)*ny+j]/dmap[(i+1)*ny+j];
        u4 = u[k*nh+(i+1)*ny+j+1]/dmap[(i+1)*ny+j+1];
        v1 = v[k*nh+i*ny+j]/dmap[i*ny+j];
        v2 = v[k*nh+i*ny+j+1]/dmap[i*ny+j+1];
        v3 = v[k*nh+(i+1)*ny+j]/dmap[(i+1)*ny+j];
        v4 = v[k*nh+(i+1)*ny+j+1]/dmap[(i+1)*ny+j+1];
        vr[k*nh+i*ny+j] = xmap[i*ny+j]*xmap[i*ny+j]*
                            ds2r*((v4-v2+v3-v1)-(u2-u1+u4-u3));
        dv[k*nh+i*ny+j] = xmap[i*ny+j]*xmap[i*ny+j]*
                            ds2r*((u3-u1+u4-u2)+(v2-v1+v4-v3));
      }
    }
  }
  return;
}

void atmcalc::do_calc(atmodata &a, t_atm_deriv &d)
{
  calcp(a.psa);
  calcrh(a.t, a.qv);
  calctd(a.t);
  calcpt(a.t);
  calcht(a.psa, a.t);
  calcdv(a.u, a.v);
  d.p = p;
  d.rh = rh;
  d.td = td;
  d.pt = pt;
  d.ht = ht;
  d.vr = vr;
  d.dv = dv;
}
