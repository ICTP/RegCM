#!/bin/bash

# This can be used to add a random noise to an ICBC file
# The range here is -0.25,+0.25

export GSL_RNG_TYPE=mt19937
for file in $@
do
  ncap2 -O -s 't=t+0.5*(0.5-gsl_rng_uniform(t))' \
          $file `basename $file.nc`_perturb.nc
done
