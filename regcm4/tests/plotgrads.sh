#!/bin/bash

if [ $# -ne 1 ]
then
  echo "Need input ctl filename"
  exit 1
fi

ctlfile=$1
varlist=`cat $ctlfile | awk 'BEGIN{found="0"} \
                       {if ($1 == "vars") found="1"; \
                        if ($1 == "endvars") found="0"; \
                        if (found=="1") print $1}'`

for var in $varlist
do
  grads -blc "run grads/plotfinal.gs $ctlfile $var"
done
