#!/bin/bash

if [ $# -ne 2 ]
then
  echo "Need input ctl filename and variable to plot"
  exit 1
fi

grads -blc "run grads/plotfinal.gs $1 $2"
