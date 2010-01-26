#!/bin/bash -f                                     
# a simple script to compile the regcm code for a specific regcm.param2
#
# Regcm distribution can be obtained from svn server hosted on gforge.escience-lab.org portal 
#
# just do: svn checkout --username anonsvn https://svn.gforge.escience-lab.org/svn/regcm
#
#SC 18/09/2009
#
# define the location of the standard Regcm distribution:

REGCM_DIR=/sp6/usersissa/sisci001/scratch/regcm/branches/regcm4/Main

if [ ! -d  "$REGCM_DIR" ]
then
   # Not a directory 
   echo -e "$REGCM_DIR not a directory..." 1>&2
   exit 1
fi
# location of the benchmark given as input 
# at the moment 29/09/2009  we should have four benchmarks:
# africatest-4 
# RegCM_192x208  
# RegCM_48x52 
# RegCM_bi 

if test $# -ne 2
then
   # Not enough arguments
   echo -e "Usage:\n$0 name_of_the_benchmark (full path) NPROC  " 1>&2
   exit 1
fi

MYDIR=$1                    
NPROC=$2
if test ! -d "$MYDIR" 
then
   # Not a directory 
   echo -e "$MYDIR not a directory..." 1>&2
   exit 1
fi

cd $MYDIR

# copy parameter files in the Main dir of the RegCM distribution to compile statically regcm
# 
# Please: be sure to backup previous versions.. ( to be done)

# this is also needed 
if test -f parame
then
   cp $REGCM_DIR/Main/parame  $REGCM_DIR/Main/parame.backup
   echo "parame file present: copying it to $REGCM_DIR/Main/parame.backup" 1>&2
   cp parame $REGCM_DIR/Main/parame
fi

echo "copying regcm.param and regcm.param2 to $REGCM_DIR" 1>&2
cp regcm.param $REGCM_DIR/Main
cp regcm.param2 $REGCM_DIR/Main

sed -i "s/PARAMNPROC/$NPROC/g" $REGCM_DIR/Main/regcm.param2

if [ $? -eq 0 ]         # Test exit status of  previous command.
 then 
 echo "files copied succesfully" 
 else
 echo "problems in copying files: stop" 
 exit
fi
# do compilation: 
# to be fixed accordingly to the new rules we want to adopt..

cd $REGCM_DIR/Main


#be sure parallel version is selected 
make clean                                        
make                                              
if [ $? -eq 0 ]         # Test exit status of  previous command.
then
  echo " make completed succesfull " 
else  
  echo " make failed: exiting.., check  why "
  exit 
fi


# go back in benchmark dir and copy the executable back 
cd $MYDIR                         
NPROC=$(printf "%02d" $NPROC)
cp $REGCM_DIR/Main/regcm Run/regcm-${NPROC}pe.x                               
if [ $? -eq 0 ]         # Test exit status of  previous command.
then
  echo " RegCM copied and ready to be executed in Run directory  " 
else
  echo " copy failed: exiting..  "
fi


# now run the executable uisng parallel enviroment available on the machine..
#  mpirun -np ... ./regcm <./regcm.in                                
