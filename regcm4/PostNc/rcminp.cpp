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

#include <iostream>
#include <fstream>
#include <cstring>
#include <cctype>
#include <rcminp.h>

using namespace rcm;

char *strstrip(char *s)
{
  size_t size = strlen(s);
  if (!size)
    return s;
  char *end = s + size - 1;
  while (end >= s && isspace(*end))
    end--;
  *(end + 1) = '\0';
  while (*s && isspace(*s))
    s++;
  return s;
}

rcminp::rcminp(char *fname)
{
  std::ifstream rcinp;
  rcinp.open(fname);
  if (! rcinp.good()) throw "Error opening RegCM input file";

  char buf[256];
  char tok1[256];
  char tok2[256];
  while (! rcinp.eof())
  {
    rcinp.getline(buf, 256);
    if (sscanf(buf, "%[ A-z0-9_]=%[ A-z0-9.', -]\n", tok1, tok2) < 2)
      continue;
    // Strip last comma from tok2
    char *p = strrchr(tok2, ',');
    if (p != NULL) *p = 0;
    items.insert(std::pair<std::string,std::string>(strstrip(tok1),
                                                    strstrip(tok2)));
  }
  rcinp.close();
}

const char *rcminp::value(const char *key)
{
  std::map<std::string, std::string>::iterator iter = items.find(key);
  return iter->second.c_str();
}

#ifdef TESTME

int main(int argc, char *argv[])
{
  rcminp a(argv[1]);
  std::cout << a.value("idate1") << std::endl;
  return 0;
}

#endif
