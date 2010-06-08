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

#ifndef __CALC__H__
#define __CALC__H__

#include <rcmio.h>

namespace rcm {

  typedef struct {
    float *r2;
  } t_srf_deriv;

  class srfcalc {
    public:
      srfcalc(header_data &h);
      ~srfcalc( );
      void do_calc(srfdata &s, t_srf_deriv &d);
    private:
      void calcrh(float *sp, float *t2, float *q2);
      int nh;
      float *r2;
  };

  class subcalc {
    public:
      subcalc(subdom_data &s);
      ~subcalc( );
      void do_calc(subdata &s, t_srf_deriv &d);
    private:
      void calcrh(float *sp, float *t2, float *q2);
      int nh;
      float *r2;
  };

  typedef struct {
    float *p;
    float *rh;
    float *td;
    float *tp;
    float *hg;
    float *vr;
    float *dv;
  } t_atm_deriv;

  class atmcalc {
    public:
      atmcalc(header_data &h);
      ~atmcalc( );
      void do_calc(atmodata &a, t_atm_deriv &d);
    private:
      void calcp(float *sp);
      void calcrh(float *t, float *q);
      void calctd(float *t);
      void calctp(float *t);
      void calchg(float *sp, float *t);
      void calcdv(float *u, float *v);
      float ptop;
      float ds;
      float ds2r;
      int nk;
      int nx;
      int ny;
      int nh;
      float *hsigf;
      float *hsigm;
      float *zs;
      float *xmap;
      float *dmap;
      float *p;
      float *rh;
      float *td;
      float *tp;
      float *hg;
      float *vr;
      float *dv;
  };

  class presslevs {
    public:
      presslevs(float *plev, int npl, float ptop);
      ~presslevs( );
      void setup_dims(int nx, int ny, int nz, float *sigma);
      void intlin(float *fp, float *f, float *pstar);
      void intlog(float *fp, float *f, float *pstar);
      void slpres(float *h, float *t, float *pstar, float *ht, float *tg,
                  float *slp1, float *slp2);
    private:
      bool dimsok;
      int nx;
      int ny;
      int nz;
      int np;
      int n2D;
      int n3D;
      int km;
      int kbc;
      int start_tl;
      float *plevs;
      float ptop;
      float *sig;
      float sigmax;
      float sigmin;
      float *pmsl;
  };

}

#endif
