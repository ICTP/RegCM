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

namespace rcm
{
  class rcmNc {
    public:
      rcmNc(char *fname, char *experiment, header_data &h);
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
      void put_rec(atmodata &a);
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
      NcVar *qvvar;
      NcVar *qcvar;
      unsigned int count;
  };
}

#endif
