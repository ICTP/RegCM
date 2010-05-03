/***************************************************************************
 *   Copyright (C) 2008-2009 by Graziano Giuliani                          *
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

#ifndef __RCMNC__H__
#define __RCMNC__H__

#include <netcdf.hh>
#include <rcmio.h>
#include <calc.h>

namespace rcm
{
  class domNc {
    public:
      domNc(char *fname, char *experiment);
      ~domNc();
      void write(domain_data &d);
    private:
      NcFile *f;
  };

  class bcNc {
    public:
      bcNc(char *fname, char *experiment, domain_data &d);
      ~bcNc( );
      void put_rec(bcdata &b);
    private:
      NcFile *f;
      NcDim *iy;
      NcDim *jx;
      NcDim *kz;
      NcDim *tt;
      NcDim *soil;
      NcVar *timevar;
      NcVar *uvar;
      NcVar *vvar;
      NcVar *tvar;
      NcVar *qvvar;
      NcVar *psvar;
      NcVar *tsvar;
      NcVar *so4var;
      NcVar *smvar;
      NcVar *itvar;
      NcVar *stvar;
      NcVar *snam;
      unsigned long reference_time;
      unsigned int rcount;
      unsigned int tcount;
      bool notinit;
      bool secondt;
      int dt;
  };

  class rcmNc {
    public:
      rcmNc(char *fname, char *experiment, header_data &h, bool full);
      ~rcmNc();
      void increment_time( ) { tcount ++; }
      NcFile *f;
      NcDim *iy;
      NcDim *jx;
      NcDim *kz;
      NcDim *tt;
      NcVar *timevar;
      unsigned long reference_time;
      unsigned int rcount;
      unsigned int tcount;
  };

  class rcmNcAtmo : public rcmNc {
    public:
      rcmNcAtmo(char *fname, char *experiment, header_data &h);
      void put_rec(atmodata &a, t_atm_deriv &d);
    private:
      NcVar *psvar;
      NcVar *tprvar;
      NcVar *tgbvar;
      NcVar *swtvar;
      NcVar *rnovar;
      NcVar *uvar;
      NcVar *vvar;
      NcVar *ovar;
      NcVar *tvar;
      NcVar *rhvar;
      NcVar *tdvar;
      NcVar *tpvar;
      NcVar *pvar;
      NcVar *htvar;
      NcVar *dvvar;
      NcVar *vrvar;
      NcVar *qvvar;
      NcVar *qcvar;
  };

  class rcmNcSrf : public rcmNc {
    public:
      rcmNcSrf(char *fname, char *experiment, header_data &h);
      void put_rec(srfdata &s, t_srf_deriv &d);
    private:
      NcVar *u10mvar;
      NcVar *v10mvar;
      NcVar *uvdragvar;
      NcVar *tgvar;
      NcVar *tlefvar;
      NcVar *t2mvar;
      NcVar *q2mvar;
      NcVar *smwvar;
      NcVar *tprvar;
      NcVar *evpvar;
      NcVar *runoffvar;
      NcVar *scvvar;
      NcVar *senavar;
      NcVar *flwvar;
      NcVar *fswvar;
      NcVar *flwdvar;
      NcVar *sinavar;
      NcVar *prcvvar;
      NcVar *psbvar;
      NcVar *zpblvar;
      NcVar *tgmaxvar;
      NcVar *tgminvar;
      NcVar *t2maxvar;
      NcVar *t2minvar;
      NcVar *w10maxvar;
      NcVar *ps_minvar;
      NcVar *r2mvar;
  };

  class rcmNcRad : public rcmNc {
    public:
      rcmNcRad(char *fname, char *experiment, header_data &h);
      void put_rec(raddata &r);
    private:
      NcVar *psvar;
      NcVar *cldvar;
      NcVar *clwpvar;
      NcVar *qrsvar;
      NcVar *qrlvar;
      NcVar *frsavar;
      NcVar *frlavar;
      NcVar *clrstvar;
      NcVar *clrssvar;
      NcVar *clrltvar;
      NcVar *clrlsvar;
      NcVar *solinvar;
      NcVar *sabtpvar;
      NcVar *firtpvar;
  };

  class rcmNcChe : public rcmNc {
    public:
      rcmNcChe(char *fname, char *experiment, header_data &h);
      void put_rec(chedata &r);
    private:
      NcDim *trc;
      NcVar *psvar;
      NcVar *trac3Dvar;
      NcVar *colbvar;
      NcVar *wdlscvar;
      NcVar *wdcvcvar;
      NcVar *sdrdpvar;
      NcVar *xgascvar;
      NcVar *xaqucvar;
      NcVar *emissvar;
      NcVar *aext8var;
      NcVar *assa8var;
      NcVar *agfu8var;
      NcVar *acstoarfvar;
      NcVar *acstsrrfvar;
      long nx, ny;
  };

  class rcmNcSub : public rcmNc {
    public:
      rcmNcSub(char *fname, char *experiment, header_data &h, subdom_data &s);
      void put_rec(subdata &s, t_srf_deriv &d);
    private:
      NcVar *u10mvar;
      NcVar *v10mvar;
      NcVar *uvdragvar;
      NcVar *tgvar;
      NcVar *tlefvar;
      NcVar *t2mvar;
      NcVar *q2mvar;
      NcVar *r2mvar;
      NcVar *smwvar;
      NcVar *tprvar;
      NcVar *evpvar;
      NcVar *runoffvar;
      NcVar *scvvar;
      NcVar *senavar;
      NcVar *prcvvar;
      NcVar *psbvar;
  };

}

#endif
