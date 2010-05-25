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

#ifndef __RCMNC__H__
#define __RCMNC__H__

#include <string>
#include <list>

#include <netcdf.hh>
#include <rcmio.h>
#include <calc.h>
#include <gradsctl.h>

namespace rcm
{
  class regcmvar {
    public:
     void addvar(std::string name);
     bool isthere(std::string name);
    private:
     std::list <std::string> rcmname;
  };

  class regcmout {
    public:
      std::string fname;
      std::string experiment;
      regcmvar vl;
      gradsctl ctl;
  };

  class domNc {
    public:
      domNc(regcmout &fnc);
      ~domNc();
      void write(domain_data &d);
    private:
      gradsctl *ctl;
      NcFile *f;
  };

  class bcNc {
    public:
      bcNc(regcmout &fnc, domain_data &d);
      ~bcNc( );
      void put_rec(bcdata &b);
    private:
      gradsctl *ctl;
      NcFile *f;
      NcDim *iy;
      NcDim *jx;
      NcDim *kz;
      NcDim *tt;
      NcDim *soil;
      NcVar *timevar;
      NcVar *psvar;
      bool varmask[10];
      NcVar *uvar;
      NcVar *vvar;
      NcVar *tvar;
      NcVar *qvvar;
      NcVar *tsvar;
      // Active if so4
      NcVar *so4var;
      // Active for USGS (Are used?)
      NcVar *smvar;
      NcVar *itvar;
      NcVar *stvar;
      NcVar *snam;
      unsigned long reference_time;
      unsigned int rcount;
      unsigned int tcount;
      bool notinit;
  };

  class rcmNc {
    public:
      rcmNc(regcmout &fnc, header_data &h, bool full);
      ~rcmNc();
      void increment_time( ) { tcount ++; }
      gradsctl *ctl;
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
      rcmNcAtmo(regcmout &fnc, header_data &h);
      void put_rec(atmodata &a, t_atm_deriv &d);
    private:
      NcVar *psvar;
      bool varmask[17];
      NcVar *uvar;
      NcVar *vvar;
      NcVar *ovar;
      NcVar *tvar;
      NcVar *rhvar;
      NcVar *tdvar;
      NcVar *tpvar;
      NcVar *pvar;
      NcVar *hgtvar;
      NcVar *dvvar;
      NcVar *vrvar;
      NcVar *qvvar;
      NcVar *qcvar;
      NcVar *tprvar;
      NcVar *tgbvar;
      NcVar *swtvar;
      NcVar *rnovar;
  };

  class rcmNcSrf : public rcmNc {
    public:
      rcmNcSrf(regcmout &fnc, header_data &h);
      void put_rec(srfdata &s, t_srf_deriv &d);
    private:
      NcVar *psbvar;
      NcVar *tbnd;
      bool varmask[26];
      NcVar *u10mvar;
      NcVar *v10mvar;
      NcVar *uvdragvar;
      NcVar *tgvar;
      NcVar *tlefvar;
      NcVar *t2mvar;
      NcVar *r2mvar;
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
      NcVar *zpblvar;
      NcVar *tgmaxvar;
      NcVar *tgminvar;
      NcVar *t2maxvar;
      NcVar *t2minvar;
      NcVar *w10maxvar;
      NcVar *ps_minvar;
      unsigned long last_time;
  };

  class rcmNcRad : public rcmNc {
    public:
      rcmNcRad(regcmout &fnc, header_data &h);
      void put_rec(raddata &r);
    private:
      NcVar *psvar;
      bool varmask[13];
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
      rcmNcChe(regcmout &fnc, header_data &h);
      void put_rec(chedata &r);
    private:
      NcDim *trc;
      NcVar *psvar;
      bool varmask[13];
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
      rcmNcSub(regcmout &fnc, header_data &h, subdom_data &s);
      void put_rec(subdata &s, t_srf_deriv &d);
    private:
      NcVar *psbvar;
      bool varmask[15];
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
  };

}

#endif
