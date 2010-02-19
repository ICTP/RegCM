#!/bin/bash

if [ $# -ne 1 ]
then
  echo "Need output directory name."
  exit 1
fi

outdir=$1

for ctlfile in $outdir/*.ctl
do
  varlist=`cat $ctlfile | awk 'BEGIN{found="0"} \
                         {if ($1 == "vars") found="1"; \
                          if ($1 == "endvars") found="0"; \
                          if (found=="1") print $1}'`

  for var in $varlist
  do
    grads -blc "run grads/plotfinal.gs $ctlfile $var"
  done
done
