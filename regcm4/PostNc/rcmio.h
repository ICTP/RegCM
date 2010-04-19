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

#ifndef __RCMIO__H__
#define __RCMIO__H__

#include <cstring>
#include <fstream>
#include <rcminp.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

namespace rcm
{

  static const char separator[2] = "/";

  // Header data
  class header_data {
    public:
      header_data(rcminp &inp);
      ~header_data( );
      const char *lbcs( );
      const char *cums( );
      const char *cumsclos( );
      const char *pbls( );
      const char *moists( );
      const char *ocnflxs( );
      const char *pgfs( );
      const char *emisss( );
      const char *lakes( );
      const char *chems( );
      void free_space( );
      // Base dates for run
      int mdate0;
      int idate0;
      int idate1;
      int idate2;
      // Model configuration
      bool ifrest;
      unsigned short ibltyp;
      unsigned short icup;
      unsigned short igcc;
      unsigned short imoist;
      unsigned short iboudy;
      unsigned short idirect;
      unsigned short iocnflx;
      unsigned short ipgf;
      unsigned short iemiss;
      unsigned short lakemod;
      unsigned short ichem;
      float radfrq;
      float abatm;
      float abemh;
      float dt;
      float ibdyfrq;
      // Dimensions
      unsigned int iy, jx, kz;
      unsigned int nx, ny, nz;
      // Geolocation
      char proj[7];
      float clat, clon, ds, xplat, xplon, trlat1, trlat2;
      // Timesteps for output
      float dto, dtb, dtr, dtc;
      // Vertical sigma levels ant top pressure
      float ptop;
      float *hsigf;
      float *ht;
      float *htsd;
      float *veg2d;
      float *landuse;
      float *xlat;
      float *xlon;
      float *xmap;
      float *dmap;
      float *coriol;
      float *mask;
  };

  class subdom_data {
    public:
      subdom_data( );
      ~subdom_data( );
      unsigned int nsg;
      unsigned int nx, ny;
      float ds, clat, clon, xplat, xplon, trlat1, trlat2;
      char proj[6];
      float *ht;
      float *htsd;
      float *landuse;
      float *xlat;
      float *xlon;
      float *xmap;
      float *coriol;
      float *mask;
  };

  class atmodata {
    public:
      atmodata(header_data &h);
      ~atmodata( );
      float *u;
      float *v;
      float *omega;
      float *t;
      float *qv;
      float *qc;
      float *psa;
      float *tpr;
      float *tgb;
      float *swt;
      float *rno;
      int date0;
      float dt;
      size_t datasize;
      int n3D;
      int n2D;
      size_t size3D;
      size_t size2D;
      int nvals;
      char *buffer;
  };

  class raddata {
    public:
      raddata(header_data &h);
      ~raddata( );
      float *cld;
      float *clwp;
      float *qrs;
      float *qrl;
      float *frsa;
      float *frla;
      float *clrst;
      float *clrss;
      float *clrlt;
      float *clrls;
      float *solin;
      float *sabtp;
      float *firtp;
      float *psa;
      int date0;
      float dt;
      size_t datasize;
      int n2D;
      int n3D;
      size_t size2D;
      size_t size3D;
      int nvals;
      char *buffer;
  };

  class chedata {
    public:
      chedata(header_data &h);
      ~chedata( );
      float *trac3D;
      float *trac2D;
      float *aext8;
      float *assa8;
      float *agfu8;
      float *acstoarf;
      float *acstsrrf;
      float *psa;
      int date0;
      float dt;
      size_t datasize;
      int n2D;
      int n3D;
      size_t size2D;
      size_t size3D;
      int nvals;
      int ntr;
      char *buffer;
  };

  class subdata {
    public:
      subdata(header_data &h, subdom_data &s);
      ~subdata( );
      float *u10m;
      float *v10m;
      float *uvdrag;
      float *tg;
      float *tlef;
      float *t2m;
      float *q2m;
      float *smw;
      float *tpr;
      float *evp;
      float *runoff;
      float *scv;
      float *sena;
      float *prcv;
      float *psb;
      int date0;
      float dt;
      size_t datasize;
      int n2D;
      size_t size2D;
      int nvals;
      char *buffer;
  };

  class srfdata {
    public:
      srfdata(header_data &h);
      ~srfdata( );
      float *u10m;
      float *v10m;
      float *uvdrag;
      float *tg;
      float *tlef;
      float *t2m;
      float *q2m;
      float *smw;
      float *tpr;
      float *evp;
      float *runoff;
      float *scv;
      float *sena;
      float *flw;
      float *fsw;
      float *flwd;
      float *sina;
      float *prcv;
      float *psb;
      float *zpbl;
      float *tgmax;
      float *tgmin;
      float *t2max;
      float *t2min;
      float *w10max;
      float *ps_min;
      int date0;
      float dt;
      size_t datasize;
      int n2D;
      size_t size2D;
      int nvals;
      char *buffer;
  };

  class rcmio
  {
    public:
      rcmio(char *dirname, bool lbig, bool ldirect);
      ~rcmio( );
      void read_header(header_data &h);
      void read_subdom(header_data &h, subdom_data &s);
      int atmo_read_tstep(atmodata &a);
      int srf_read_tstep(srfdata &a);
      int rad_read_tstep(raddata &a);
      int che_read_tstep(chedata &a);
      int sub_read_tstep(subdata &a);
      bool has_atm;
      bool has_che;
      bool has_srf;
      bool has_rad;
      bool has_sub;
    private:
      inline void swap2(char *p)
      {
        unsigned short *x;
        x = (unsigned short*) p;
        *x = (*x>>8) |
             (*x<<8);
      }
      inline void swap4(char *p)
      {
        unsigned int *x;
        x = (unsigned int *) p;
        *x = (*x>>24) |
            ((*x<<8) & 0x00FF0000) |
            ((*x>>8) & 0x0000FF00) |
             (*x<<24);
      }
      inline void swap8(char *p)
      {
        unsigned long long *x;
        x = (unsigned long long *) p;
        *x =  (*x>>56) |
             ((*x<<40) & 0x00FF000000000000) |
             ((*x<<24) & 0x0000FF0000000000) |
             ((*x<<8)  & 0x000000FF00000000) |
             ((*x>>8)  & 0x00000000FF000000) |
             ((*x>>24) & 0x0000000000FF0000) |
             ((*x>>40) & 0x000000000000FF00) |
              (*x<<56);
      }
      inline int intvalfrombuf(char *p)
      {
        if (doswap) swap4(p);
        return *((int *) p);
      }
      inline float rvalfrombuf(char *p)
      {
        if (doswap) swap4(p);
        return *((float *) p);
      }
      inline void vectorfrombuf(char *b, float *v, int nvals)
      {
        char *p = b;
        if (doseq) p += sizeof(int);
        for (int i = 0; i < nvals; i ++)
        {
          v[i] = rvalfrombuf(p); p += sizeof(float);
        }
        return;
      }
      inline bool littlearch( )
      {
        unsigned char test[2] = { 1, 0 };
        if( *(short *) test == 1 ) return true;
        return false;
      }
      inline bool bigarch( )
      {
        return (! littlearch());
      }
      // Directory with data and IO flags
      char outdir[PATH_MAX];
      bool doswap;
      bool doseq;
      // ATM pieces
      bool initatm;
      std::ifstream atmf;
      size_t atmsize;
      // CHE pieces
      bool initche;
      std::ifstream chef;
      size_t chesize;
      // SRF pieces
      bool initsrf;
      std::ifstream srff;
      size_t srfsize;
      // RAD pieces
      bool initrad;
      std::ifstream radf;
      size_t radsize;
      // SUB pieces
      bool initsub;
      std::ifstream subf;
      size_t subsize;
      // just I/O space
      char *storage;
      size_t readsize;
  };
}

#endif
