# This file contains variables to be included in all other Makefiles

# libraries and data
NETCDFLIB = -L$(NETCDF_LIB) -lnetcdf
NETCDFINC = -I$(NETCDF_INC)
REGCM_DATA_DIR = ./Input/DATA
REGCM_BASE_DIR = myregcmdir

# general compile options
CPPFLAGS = -DDIAG -DMPP1
USCORING = -fno-underscoring
FFLAGS =  -O0 -g -fbounds-check -ffpe-trap=zero -fconvert=big-endian -Wall -pedantic
FC = mpif90
F90 = mpif90
LD =$(FC)
