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

#include <iostream>
#include <fstream>
#include <cstring>
#include <cctype>
#include <rcminp.h>

using namespace rcm;

char *rcminp::strstrip(char *s)
{
  size_t size = strlen(s);
  if (!size)
    return s;
  for (size_t i = 0; i < size; i ++)
    *s = (char) tolower(*s);
  char *end = s + size - 1;
  while (end >= s && isspace(*end))
    end--;
  *(end + 1) = '\0';
  while (*s && isspace(*s))
    s++;
  char *p = strrchr(s, ',');
  if (p != NULL) *p = 0;
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
  std::string s1;
  std::string s2;
  while (! rcinp.eof())
  {
    rcinp.getline(buf, 256);
    if (sscanf(buf, "%[ A-z0-9_]=%[ A-z0-9.',+_/ -]\n", tok1, tok2) < 2)
      continue;
    s1 = strstrip(tok1);
    s2 = strstrip(tok2);
    items.insert(std::make_pair(s1,s2));
  }
  rcinp.close();
}

const char *rcminp::valuec(const char *key)
{
  std::map<std::string, std::string>::iterator iter = items.find(key);
  if (iter == items.end()) throw "Item not found";
  size_t found;
  while ( (found = iter->second.find("'")) != std::string::npos)
    iter->second.erase(found, 1);
  return iter->second.c_str();
}

int rcminp::valuei(const char *key)
{
  std::map<std::string, std::string>::iterator iter = items.find(key);
  if (iter == items.end()) throw "Item not found";
  int vl;
  if (sscanf(iter->second.c_str(), "%d", &vl) != 1)
    throw "rcminp::value_int : cannot parse to integer";
  return vl;
}

float rcminp::valuef(const char *key)
{
  std::map<std::string, std::string>::iterator iter = items.find(key);
  if (iter == items.end()) throw "Item not found";
  float vl;
  if (sscanf(iter->second.c_str(), "%f", &vl) != 1)
    throw "rcminp::value_real : cannot parse to real";
  return vl;
}

bool rcminp::valueb(const char *key)
{
  std::map<std::string, std::string>::iterator iter = items.find(key);
  if (iter == items.end()) throw "Item not found";
  if (strstr(iter->second.c_str(), "false") != NULL) return false;
  return true;
}

#ifdef TESTME

int main(int argc, char *argv[])
{
  rcminp a(argv[1]);
  std::cout << a.valuei("idate1") << std::endl;
  return 0;
}

#endif
