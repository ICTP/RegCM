// Simple program to fill holes in the Era Interim SST dataset.
// It does use a simple distance weight interpolation to set a temperature on
// land points from neighbour sea points.

#include <netcdf.hh>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>

using namespace std;

static int nlat = 121;
static int nlon = 240;
static int ntot = nlon*nlat;
static int nsize = ntot*sizeof(short);
static int ntimes = 29824;
static short missval = -32767;

int main(int argc, char *argv[])
{
  short *iobuf1 = new short[ntot];
  short *iobuf2 = new short[ntot];

//  NcFile *sstwf = new NcFile("sstERAIN.1989-2009.nc",NcFile::Write);
  NcFile *sstwf = new NcFile("peppe.nc",NcFile::Write);
  NcFile *sktwf = new NcFile("tskinERAIN.1989-2009.nc",NcFile::ReadOnly);
  NcVar *ncvtmp1 = sstwf->get_var("sst");
  NcVar *ncvtmp2 = sktwf->get_var("skt");

  for (int i = 0; i < ntimes; i ++)
  {
    // Read data
    ncvtmp1->set_cur(i,0,0);
    ncvtmp2->set_cur(i,0,0);
    ncvtmp1->get(iobuf1,1,nlat,nlon);
    ncvtmp2->get(iobuf2,1,nlat,nlon);

    for ( int j = 0; j < ntot; j ++)
      if ( iobuf1[j] == missval ) iobuf1[j] = iobuf2[j];

    // Write back data
    ncvtmp1->put(iobuf1,1,nlat,nlon);
    sstwf->sync();
  }

  sstwf->close();
  sktwf->close();

  delete [] iobuf1;
  delete [] iobuf2;

  return 0;
}
