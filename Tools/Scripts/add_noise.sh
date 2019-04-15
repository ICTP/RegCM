#!/bin/bash

# This can be used to add a random noise to an ICBC file
# The range here is -0.25,+0.25
# This method is deprecated. Use the perturbparam in the namelist!

export GSL_RNG_TYPE=mt19937
export GSL_RNG_SEED=`date +%s`
for file in $@
do
  ncap2 -O -s 't=t+0.5*(0.5-gsl_rng_uniform_pos(t))' \
          $file `basename $file.nc`_perturb.nc
done
