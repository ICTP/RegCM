# This file contains variables to be included in all other Makefiles

NETCDF_LIB=.

# libraries and data
NETCDFLIB = -L$(NETCDF_LIB) -lnetcdf
NETCDFINC = -I$(NETCDF_INC)
REGCM_DATA_DIR = ./Input/DATA
REGCM_BASE_DIR = myregcmdir

# general compile options
CPPFLAGS = -DDIAG
USCORING = -fno-underscoring
F90FLAGS =  -O0 -g -fbounds-check -ffpe-trap=zero -fconvert=big-endian -Wall -pedantic
F90 = gfortran
LD =$(F90)
