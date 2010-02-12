#!/bin/sh
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of RegCM model.
#
#    RegCM model is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    RegCM model is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if [ $# -lt 2 ]
then
  echo "Welcome to job launcher."
  echo "I have not enough informations. Typical usage is:"
  echo
  echo "   $0 myregcm.in"
  echo
  echo " or"
  echo
  echo "   $0 myregcm.in mpi MPI_ARGUMENTS"
  echo
  echo "where MPI_ARGUMENTS are the arguments to be passed to mpirun,"
  echo "like \"-np 16\" or \"-hostfile filename\""
  echo
  echo "Good bye !"
  exit 1
fi

echo "This is Regional Climatic Model job launcher."

echo "Performing basic environment checks."

if [ ! -f $1 ]
then
  echo "Namelist file $1 not found."
  echo "Error. Input namelist file is not present in this directory."
  echo "Please read documentation on how to write it."
  exit 1
fi

ln -sf $1 regcm.in
idate0=`cat regcm.in | grep idate0 | awk 'BEGIN{FS="="}{print $2}'`
idate1=`cat regcm.in | grep idate1 | awk 'BEGIN{FS="="}{print $2}'`
idate2=`cat regcm.in | grep idate2 | awk 'BEGIN{FS="="}{print $2}'`

# Check I/O directories

if [ ! -d output ]
then
  mkdir output
fi

if [ ! -d ../Postproc ]
then
  mkdir ../Postproc
fi

# Check if a postproc namelist and param already exist and backup it

if [ -f ../Postproc/postproc.in ]
then
  mv -f postproc.in postproc.in.bak
fi
if [ -f ../Postproc/postproc.param ]
then
  mv -f postproc.param postproc.param.bak
fi

echo "Starting simulation on `date` for $USER on `hostname`"

echo "Linking Model executable here."
ln -sf ../Main/regcm .

command="false"
if [ $# -gt 1 ]
then
  if [ "$2" = "mpi" ]
  then
    shift
    shift
    command="mpirun $@ $PWD/regcm"
  fi
else
 command="$PWD/regcm < regcm.in"
fi 

$command
if [ $? -ne 0 ]
then
  echo "Model error. Check it before trying to run again."
  exit 1
fi

echo "Simulation completed on `date` for $USER on `hostname`"

echo "Swapping new input namelist file in place for new run start."
if [ -f regcm.in.new ]
then
  mv regcm.in regcm.in.$idate0
  mv regcm.in.new regcm.in
fi

exit 0
