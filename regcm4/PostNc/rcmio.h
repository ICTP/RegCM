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

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

namespace rcm
{
  // Header data
  class header_data {
    public:
      header_data( );
      ~header_data( );
      const char *lbcs( );
      const char *cums( );
      const char *pbls( );
      const char *moists( );
      void free_space( );
      // Base date for run
      int mdate0;
      // Model configuration
      unsigned short ibltyp;
      unsigned short icup;
      unsigned short imoist;
      unsigned short iboudy;
      unsigned short idirect;
      // Dimensions
      unsigned int iy, jx, kz;
      unsigned int nx, ny, nz;
      // Geolocation
      char proj[7];
      float clat, clon, ds, xplat, xplon;
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

  class atmodata {
    public:
      atmodata(int nx, int ny, int nz, int mdate0, float dto);
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
      raddata(int nx, int ny, int nz, int mdate0, float dts);
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

  class srfdata {
    public:
      srfdata(int nx, int ny, int mdate0, float dts);
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
      int atmo_read_tstep(atmodata &a);
      int srf_read_tstep(srfdata &a);
      int rad_read_tstep(raddata &a);
      bool has_atmo;
      bool has_chem;
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
      bool initatmo;
      std::ifstream atmof;
      size_t atmsize;
      // CHE pieces
      bool initchem;
      std::ifstream chemf;
      size_t chemsize;
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
      char *storage;
  };
}

#endif
