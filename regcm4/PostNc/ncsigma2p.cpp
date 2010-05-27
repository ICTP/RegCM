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
#include <vector>

#include <getopt.h>

#include <netcdf.hh>
#include <calc.h>

using namespace rcm;

static void help(char *pname);
static const char version[] = SVN_REV;
static void levlist_to_vec(std::string &ll, std::vector<float> &lf);
static void copygatt(NcFile *f, NcAtt *src);
static void copyvatt(NcVar *v, NcAtt *src);

void help(char *pname)
{
  std::cerr << std::endl
      << "                 RegCM V4 ICTP NetCDF Postprocessor." << std::endl
      << std::endl
      << "This simple program interpolates NetCDF files from RegCM V4 "
      << "Sigma levels" << std::endl
      << "onto pressure levels."
      << std::endl << "Usage:" << std::endl << std::endl << "     " << pname
      << " [options] -p lev1[,lev2...] Infile.nc"
      << std::endl << std::endl
      << "where options can be in:" << std::endl << std::endl
      << "   --plev/-p l1[,l2..] : Comma separated list of pressure levels"
      << std::endl
      << "   --help/-h           : Print this help" << std::endl
      << "   --version/-V        : Print versioning information"
      << std::endl << std::endl;
   return;
}

void levlist_to_vec(std::string &ll, std::vector<float> &lf)
{
  float val;
  std::string::size_type lp = ll.find_first_not_of(",",0);
  std::string::size_type pos = ll.find_first_of(",",lp);
  while (std::string::npos != pos || std::string::npos != lp)
  {
    sscanf(ll.substr(lp,pos-lp).c_str(), "%f", &val);
    lf.push_back(val);
    lp = ll.find_first_not_of(",",pos);
    pos = ll.find_first_of(",", lp);
  }
  return;
}

void copygatt(NcFile *f, NcAtt *src)
{
  int nvals = src->num_vals( );
  switch (src->type( ))
  {
    case ncByte:
      {
        ncbyte *vals = new ncbyte [nvals];
        for (int j = 0; j < nvals; j ++)
          vals[j] = src->as_ncbyte(j);
        f->add_att(src->name( ), nvals, vals);
        delete [ ] vals;
      }
      break;
    case ncChar:
      char *ss = src->as_string(0);
      f->add_att(src->name( ), ss);
      delete [ ] ss;
      break;
    case ncShort:
      {
        short *vals = new short [nvals];
        for (int j = 0; j < nvals; j ++)
          vals[j] = src->as_short(j);
        f->add_att(src->name( ), nvals, vals);
        delete [ ] vals;
      }
      break;
    case ncInt:
      {
        int *vals = new int [nvals];
        for (int j = 0; j < nvals; j ++)
          vals[j] = src->as_long(j);
        f->add_att(src->name( ), nvals, vals);
        delete [ ] vals;
      }
      break;
    case ncFloat:
      {
        float *vals = new float [nvals];
        for (int j = 0; j < nvals; j ++)
          vals[j] = src->as_float(j);
        f->add_att(src->name( ), nvals, vals);
        delete [ ] vals;
      }
      break;
    case ncDouble:
      {
        double *vals = new double [nvals];
        for (int j = 0; j < nvals; j ++)
          vals[j] = src->as_double(j);
        f->add_att(src->name( ), nvals, vals);
        delete [ ] vals;
      }
      break;
    default:
      throw "NetCDF error in copy global attributes";
  }
}

void copyvatt(NcVar *v, NcAtt *src)
{
  int nvals = src->num_vals( );
  switch (src->type( ))
  {
    case ncByte:
      {
        ncbyte *vals = new ncbyte [nvals];
        for (int j = 0; j < nvals; j ++)
          vals[j] = src->as_ncbyte(j);
        v->add_att(src->name( ), nvals, vals);
        delete [ ] vals;
      }
      break;
    case ncChar:
      char *ss = src->as_string(0);
      v->add_att(src->name( ), ss);
      delete [ ] ss;
      break;
    case ncShort:
      {
        short *vals = new short [nvals];
        for (int j = 0; j < nvals; j ++)
          vals[j] = src->as_short(j);
        v->add_att(src->name( ), nvals, vals);
        delete [ ] vals;
      }
      break;
    case ncInt:
      {
        int *vals = new int [nvals];
        for (int j = 0; j < nvals; j ++)
          vals[j] = src->as_long(j);
        v->add_att(src->name( ), nvals, vals);
        delete [ ] vals;
      }
      break;
    case ncFloat:
      {
        float *vals = new float [nvals];
        for (int j = 0; j < nvals; j ++)
          vals[j] = src->as_float(j);
        v->add_att(src->name( ), nvals, vals);
        delete [ ] vals;
      }
      break;
    case ncDouble:
      {
        double *vals = new double [nvals];
        for (int j = 0; j < nvals; j ++)
          vals[j] = src->as_double(j);
        v->add_att(src->name( ), nvals, vals);
        delete [ ] vals;
      }
      break;
    default:
      throw "NetCDF error in copy variable attributes";
  }
}

int main(int argc, char *argv[])
{
  std::string plevs = "1000.,925.,850.,700.,500.,400.,300.,250.,200.,150.,100";
  char *pname = basename(argv[0]);

  while (1)
  {
    static struct option long_options[] = {
      { "plev", required_argument, 0, 'p'},
      { "help", no_argument, 0, 'h'},
      { "version", no_argument, 0, 'V'},
      { 0, 0, 0, 0 }
    };
    int optind, c = 0;
    c = getopt_long (argc, argv, "p:hV",
                     long_options, &optind);
    if (c == -1) break;
    switch (c)
    {
      case 'p':
        plevs = optarg;
      case 'h':
        help(pname);
        return 0;
        break;
      case 'V':
        std::cerr << "This is " << pname << " version " << version
                  << std::endl;
        return 0;
      case '?':
        break;
      default:
        std::cerr << "Unknown switch '" << (char) c << "' discarded."
                  << std::endl;
        break;
    }
  }

  if (argc - optind != 1)
  {
    std::cerr << std::endl << "Howdy there, wrong number of arguments."
              << std::endl;
    help(pname);
    return -1;
  }

  try
  {
    float *fplev;
    int np;
    std::vector<float> fl;
    levlist_to_vec(plevs, fl);
    np = fl.size( );
    if (np <= 0) throw "Unable to understande pressure level list...";
    fplev = new float[np];
    for (int i = 0; i < np; i ++)
      fplev[i] = fl[i];
    fl.clear( );
    
    std::string iname = argv[1];
    std::string oname = iname.substr(0, iname.find_last_of('.'))+"_plevs.nc";
    NcFile *fin, *fout;
    fin  = new NcFile(iname.c_str(), NcFile::ReadOnly);
    fout = new NcFile(oname.c_str(), NcFile::Replace);

    std::vector <NcDim *> dims;
    NcDim *pdim, *tdim, *udim, *zdim;

    udim = fin->rec_dim( );
    zdim = fin->get_dim("kz");

    tdim = fout->add_dim(udim->name());
    dims.push_back(tdim);
    pdim = fout->add_dim("np", np);
    dims.push_back(pdim);
    int ndims = fin->num_dims( );
    {
      for (int i = 0; i < ndims; i ++)
      {
        NcDim *dim = fin->get_dim(i);
        if (dim != udim && dim != zdim)
        {
          NcDim *xdim = fout->add_dim(dim->name( ), dim->size( ));
          dims.push_back(xdim);
        }
      }
    }

    int ngatts = fin->num_atts( );
    for (int i = 0; i < ngatts; i ++)
    {
      NcAtt *att = fin->get_att(i);
      copygatt(fout, att);
      delete att;
    }
    NcAtt *hatt = fin->get_att("history");
    NcValues *vhatt = hatt->values( );
    char *cohist = vhatt->as_string(0);
    std::string ohist = cohist;
    delete [ ] cohist;
    delete vhatt;
    delete hatt;
    ohist += "\n";
    for (int i = 0; i < argc-1; i ++)
      ohist = ohist+argv[i]+" ";
    ohist = ohist+argv[argc-1];
    fout->add_att("history", ohist.c_str());

    int nvars = fin->num_vars( );
    bool *is3D = new bool[nvars];
    bool *isvt = new bool[nvars];
    bool *istemp = new bool[nvars];
    NcDim **xdims;
    xdims = new NcDim*[ndims];
    const NcDim **xxd = (const NcDim **) xdims;
    int levid = -1;
    int timeid = -1;
    int prid = -1;
    int psid = -1;
    for (int i = 0; i < nvars; i ++)
    {
      is3D[i] = false;
      isvt[i] = false;
      istemp[i] = false;
      NcVar *var = fin->get_var(i);
      char *name = (char *) var->name( );
      if (strcmp(name, "level") == 0)
      {
        NcVar *pl = fout->add_var("plev", ncFloat, pdim);
        pl->add_att("standard_name", "pressure");
        pl->add_att("long_name", "Pressure level");
        pl->add_att("positive", "down");
        pl->add_att("units", "hPa");
        pl->add_att("axis", "Z");
        levid = i;
        continue;
      }
      else if (strcmp(name, "p") == 0)
      {
        // Do not define a 3D pressure in output: we are going on pressure
        // levels as Z dimension.
        prid = i;
        continue;
      }
      else if (strcmp(name, "time") == 0) timeid = i;
      else if (strcmp(name, "ps") == 0) psid = i;
      else if (strcmp(name, "t") == 0)  istemp[i] = true;
      else if (strcmp(name, "td") == 0) istemp[i] = true;
      else if (strcmp(name, "tp") == 0) istemp[i] = true;
      int nvdims = var->num_dims( );
      int nvatts = var->num_atts( );
      for (int j = 0; j < nvdims; j ++)
      {
        NcDim *d = var->get_dim(j);
        if (d == zdim)
        {
          if (nvdims > 1) is3D[i] = true;
          xdims[j] = pdim;
        }
        else
        {
          if (d == udim) isvt[i] = true;
          for (unsigned int k = 0; k < dims.size( ); k++)
          {
            if (strcmp(d->name( ), dims[k]->name( )) == 0)
            {
              xdims[j] = dims[k];
              break;
            }
          }
        }
      }
      NcVar *vout = fout->add_var(var->name( ), var->type( ), nvdims, xxd);
      for (int j = 0; j < nvatts; j ++)
      {
        NcAtt *a = var->get_att(j);
        copyvatt(vout, a);
        delete a;
      }
    }
    delete [ ] xdims;

    fout->sync( );

    long actsize = 0;
    float *vals = 0;

    float ptop = fin->get_var("ptop")->as_float(0);
    int nz = zdim->size();
    float *sigma = new float[nz];
    NcVar *sigmalev = fin->get_var("level");
    sigmalev->get(sigma, nz);
    int nx = fin->get_dim("jx")->size( );
    int ny = fin->get_dim("iy")->size( );
    long size2D = nx*ny;
    long size3D = size2D*nz;
    long size3Dp = size2D*np;
    float *ps = new float[size2D];
    float *varsig = new float[size3D];
    float *varp = new float[size3Dp];
    presslevs pint(fplev, np, ptop);
    pint.setup_dims(nx, ny, nz, sigma);

    // Write out time independent data
    for (int i = 0, io = 0; i < nvars; i ++, io ++)
    {
      if (isvt[i] || is3D[i]) continue; // Assuming no static 3D+ needed
      if (i == prid)
      {
        io --;
        continue;          // Skip undefined var
      }
      NcVar *vout = fout->get_var(io);
      if (i == levid)
      {
        vout->put(fplev, np);
        continue;
      }
      NcVar *vin = fin->get_var(i);
      long *xs = vin->edges( );
      long totsize = vin->num_vals( );
      if (totsize > actsize)
      {
        if (vals) delete [ ] vals;
        vals = new float[totsize];
        actsize = totsize;
      }

      vin->get(vals, xs);
      vout->put(vals, xs);
      delete [ ] xs;
    }

    std::cout << "Static 2D vars written..." << std::endl;

    // Time dependent data
    int ntimes = udim->size( );
    for (long it = 0; it < ntimes; it ++)
    {
      // Add a timestep
      NcVar *vint = fin->get_var(timeid);
      NcVar *vt = fout->get_var(timeid);
      NcValues *xt = vint->get_rec(it);
      double xtime = xt->as_double(0);
      vt->put_rec(&xtime, it);
      delete xt;

      // Need a copy of surface pressure
      NcVar *psvar = fin->get_var("ps");
      psvar->set_cur(it);
      psvar->get(ps, 1, ny, nx);
      NcVar *psout = fout->get_var("ps");
      psout->set_cur(it);
      psout->put(ps, 1, ny, nx);

      for (int i = 0, io = 0; i < nvars; i ++, io ++)
      {
        if (i == prid)
        {
          io --;    // Do not increment outpointer
          continue; // Do not want pressure in output on plevs
        }
        if (i == psid) continue; // Already written
        if (!isvt[i])  continue; // Not function of time
        if (!is3D[i])            // 2D+ var, not vertical info
        {
          NcVar *vin2D = fin->get_var(i);
          NcVar *vout2D = fout->get_var(io);
          vin2D->set_cur(it);
          vout2D->set_cur(it);
          long *xs = vin2D->edges( );
          xs[0] = 1;
          long totsize = vin2D->num_vals( )/ntimes;
          if (totsize > actsize)
          {
            if (vals) delete [ ] vals;
            vals = new float[totsize];
            actsize = totsize;
          }
          vin2D->get(vals, xs);
          vout2D->put_rec(vals, it);
          delete [ ] xs;
          continue;
        }
        // Time dependent and 3D+
        NcVar *vin3D = fin->get_var(i);
        NcVar *vout3D = fout->get_var(io);
        vin3D->set_cur(it);
        vout3D->set_cur(it);
        long *xs = vin3D->edges( );
        long *xd = vout3D->edges( );
        long totsize = vin3D->num_vals( )/ntimes;
        int nblocks = totsize/size3D;
        xs[0] = 1;
        xd[0] = 1;
        if (nblocks > 1) 
        { 
          xs[1] = 1;
          xd[1] = 1;
        }
        for (int j = 0; j < nblocks; j ++)
        {
          vin3D->get(varsig, xs);
          if (istemp[i])
            pint.intlog(varp, varsig, ps);
          else
            pint.intlin(varp, varsig, ps);
          vout3D->put(varp, xd);
          xs[1] ++;
          xd[1] ++;
        }
        delete [ ] xs;
        delete [ ] xd;
      }
      std::cout << "Done timestep." << std::endl;
    }

    fin->close( );
    fout->close( );

    delete [ ] vals;
    delete [ ] varsig;
    delete [ ] varp;
    delete [ ] is3D;
    delete [ ] isvt;
    delete [ ] istemp;
    delete [ ] fplev;
    delete [ ] sigma;
    delete [ ] ps;
    delete fin;
    delete fout;
  }
  catch (const char *e)
  {
    std::cerr << "Error : " << e << std::endl;
    return -1;
  }
  catch (...)
  {
    return -1;
  }

  std::cout << "Successfully interpolated." << std::endl;
  return 0;
}
