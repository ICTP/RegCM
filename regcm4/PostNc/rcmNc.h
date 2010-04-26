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
  class rcmNc {
    public:
      rcmNc(char *fname, char *experiment, header_data &h, bool full);
      ~rcmNc();
      NcFile *f;
      NcDim *iy;
      NcDim *jx;
      NcDim *kz;
      NcDim *tt;
      NcVar *timevar;
      unsigned long reference_time;
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
      unsigned int count;
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
      unsigned int count;
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
      unsigned int count;
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
      unsigned int count;
  };

  class rcmNcSub : public rcmNc {
    public:
      rcmNcSub(char *fname, char *experiment, header_data &h, subdom_data &s);
      void put_rec(subdata &s);
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
      NcVar *prcvvar;
      NcVar *psbvar;
      unsigned int count;
  };

}

#endif
