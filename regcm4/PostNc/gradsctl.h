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

#ifndef __GRADSCTL__H__
#define __GRADSCTL__H__

#include <fstream>
#include <list>
#include <rcmio.h>

namespace rcm {

  class gradsvar {
    public:
      gradsvar(char *ncname, char *gname, char *desc, int nlevs, bool tv);
      gradsvar( );
      void set(char *ncname, char *gname, char *desc, int nlevs, bool tv);
    private:
      friend std::ostream& operator<< (std::ostream& os,
                                       const gradsvar &g);
      std::string vline;
  };

  class gradsctl {
    public:
      gradsctl( );
      ~gradsctl( );
      void open(char *ctlname, char *dname);
      void head(char *title, float misval);
      void set_grid(domain_data &d);
      void addentry(char *entry);
      void set_grid(header_data &d);
      void set_levs(float *lev, int nz);
      void set_time(char *timestr);
      void add_var(gradsvar &g);
      void add_time( );
      void finalize( );
      char *gradstime(int year, int month, int day, int hour, int incr);
      bool doit;
    private:
      std::ofstream ctlf;
      std::list <gradsvar> vars;
      std::string timeline;
      std::string levline;
      int ntimes;
  };

}

#endif
