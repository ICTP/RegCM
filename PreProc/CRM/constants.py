""" Parses Share/mod_constants.F90 and loads variables into this namespace """
import os
from math import sqrt

# set the location of mod_constants
script_path = os.path.dirname(os.path.abspath(__file__))
mod_constants_path = os.path.abspath("{}/../../Share/mod_constants.F90".format(script_path))

# attempt to open mod_constants.F90
try:
    # read mod_constants.F90
    with open(mod_constants_path,'r') as fin: 
        # read mod_constants and remove blank lines, indentations and newline characters
        mod_constants = [ line.rstrip().lstrip() for line in fin.readlines() if line.rstrip().lstrip() != '' ]
except:
    raise(RuntimeError,"Could not open `{}` for reading.".format(mod_constants_path))

# extract variable definition lines
var_defs = [ line.split('::')[-1].lstrip().split('!')[0].rstrip() for line in mod_constants if line.lstrip()[0] != '!' ]

# merge any line continuations
icontinue = [ i for i in range(len(var_defs)) if var_defs[i][-1] == '&' ]
# loop over continuation lines from end to beginning (so pop() doesn't mess up indices)
for i in icontinue[::-1]:
    # join this line with the line after; also remove continuation characters
    var_defs[i] = (var_defs[i] + var_defs[i+1]).replace('&',' ')
    # remove the line after now that it has been merged with this line
    var_defs.pop(i+1)

# remove any non equation lines and remove fortran precision specifiers
var_defs = [ line.replace('_rkx','').replace('_rk4','') for line in var_defs if '=' in line ]

# evaluate these lines as Python code to load the variables into memory
for line in var_defs:
    exec(line)

