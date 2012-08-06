//
// Read a SRF file in C++
//
#include <netcdf>
#include <string>
#include <iostream>
#include <iomanip>

using namespace netCDF;
using namespace netCDF::exceptions;

int main(int argc, char *argv[])
{
  NcFile * srf_file;
  try
  {
    std::cout << "Opening file " << argv[1] << std::endl;
    srf_file = new NcFile(argv[1], NcFile::read);
  }
  catch (NcException& e)
  {
    std::cerr << "Error opening File " << argv[1] << std::endl;
    e.what();
  }

  NcDim jxdim , iydim , kzdim , m10dim , m2dim , sldim , timedim;
  size_t jx , iy , kz , m10 , m2 , sl , ntime;
  try
  {
    jxdim = srf_file->getDim("jx");
    iydim = srf_file->getDim("iy");
    kzdim = srf_file->getDim("kz");
    m10dim = srf_file->getDim("m10");
    m2dim = srf_file->getDim("m2");
    sldim = srf_file->getDim("soil_layer");
    timedim = srf_file->getDim("time");
  }
  catch (NcException& e)
  {
    std::cerr << "Not all SRF dimensions present!" << std::endl;
    e.what();
  }

  jx    = jxdim.getSize( );
  iy    = iydim.getSize( );
  kz    = kzdim.getSize( );
  m10   = m10dim.getSize( );
  m2    = m2dim.getSize( );
  sl    = sldim.getSize( );
  ntime = timedim.getSize( );

  std::cout << "DIM JX   = " << jx << std::endl;
  std::cout << "DIM IY   = " << iy << std::endl;
  std::cout << "DIM KZ   = " << kz << std::endl;
  std::cout << "DIM M10  = " << m10 << std::endl;
  std::cout << "DIM M2   = " << m2 << std::endl;
  std::cout << "DIM SL   = " << sl << std::endl;
  std::cout << "DIM TIME = " << ntime << std::endl;

  NcGroupAtt projatt , dxatt , clatatt , clonatt;
  try
  {
    projatt = srf_file->getAtt("projection");
    dxatt = srf_file->getAtt("grid_size_in_meters");
    clatatt = srf_file->getAtt("latitude_of_projection_origin");
    clonatt = srf_file->getAtt("longitude_of_projection_origin");
  }
  catch (NcException& e)
  {
    std::cerr << "File doens't have projection info." << std::endl;
    e.what();
  }

  std::string proj;
  double dx , clat , clon;

  projatt.getValues(proj);
  dxatt.getValues(&dx);
  clatatt.getValues(&clat);
  clonatt.getValues(&clon);
  std::cout << "PROJECTION CODE : " << proj << std::endl;
  std::cout << "GRID SIZE       : " << dx << std::endl;
  std::cout << "PROJECT CLAT    : " << clat << std::endl;
  std::cout << "PROJECT CLON    : " << clon << std::endl;

  NcVar t2mvar;
  NcVarAtt t2munitsatt;
  try
  {
    t2mvar = srf_file->getVar("t2m");
    t2munitsatt = t2mvar.getAtt("units");
  }
  catch (NcException& e)
  {
    std::cerr << "File doens't have t2m variable." << std::endl;
    e.what();
  }

  float * t2m = new float[jx*iy*m2*ntime];
  std::string units;
  t2mvar.getVar(t2m);
  t2munitsatt.getValues(units);

  double meanval = 0;
  double count = 0;

  for ( int it = 0; it < ntime; it ++)
    for ( int il = 0; il < m2; il ++ )
      for ( int isn = 0; isn < iy; isn ++ )
        for ( int iwe = 0; iwe < jx; iwe ++ )
        {
          meanval = meanval+t2m[iwe+isn*iwe+il*isn*iwe+it*il*isn*iwe];
          count = count + 1.0;
        }

  std::cout << "The mean value of T2M is " << meanval/count 
            << " " << units << std::endl;

  delete srf_file;
  srf_file = 0;
  std::cout << "Done." << std::endl;
}
