#!/usr/bin/python
import math
import sys
import os.path

file1='../Config/mod_preproc_param.f90'
file2='../Config/mod_regcm_param.F90'
fileout='../PreProc/ICBC/icbc_v4.x'

ssttyp=''
dattyp=''
aertyp=''

# Find the sst type in mod_preproc_param

if os.path.isfile(file1):
    f=open(file1,'r')
    content=f.readlines()
    f.close()
    for line in content:
            if line.find('ssttyp') > -1:
                riga = line.split('=')
                ssttyp=riga[1].strip()
                #print sstyp
else :
    print file1,"not found!!!"
    sys.exit()

# Find the data type and aer type in mod_regcm_param

if os.path.isfile(file2):
    f=open(file2,'r')
    content=f.readlines()
    f.close()
    for line in content:
            if line.find('dattyp') > -1:
                riga = line.split('=')
                dattyp=riga[1].strip()
            if line.find('aertyp') > -1:
                riga = line.split('=')
                aertyp=riga[1].strip()
else :
    print file2,"not found!!!"
    sys.exit()

dattyp=dattyp.strip("'")
ssttyp=ssttyp.strip("'")
aertyp=aertyp.strip("'")

print ssttyp
print dattyp
print aertyp

# Write icbc_v4.x accordingly

f = open(fileout,"w")
f.write('#!/bin/csh \n')
f.write('#cleaning files.. \n')
f.write("foreach FILE (RCM_SST.ctl RCM_SST.dat SST.RCM)\n")
f.write("if ( -f $FILE ) /bin/rm $FILE\n")
f.write("end\n")

default="./sst_1deg \n"

if ((dattyp=='FVGCM') or  (dattyp == 'FNEST' and (ssttyp== 'FV_RF' or ssttyp == 'FV_A2') ) ):
    default="./sst_fvgcm\n"
if ((dattyp=='EH5OM') or (dattyp=='FNEST' and (ssttyp=='EH5RF' or ssttyp=='EH5A2' or ssttyp=='EH5B1' or ssttyp=='EHA1B'))):
    default="./sst_eh50m\n"
if ((dattyp=='EIN15' or dattyp=='ERAIN') and (ssttyp=='ERSST' or ssttyp=='ERSKT')):
    default="./sst_ersst\n" 

f.write(default) 

# chemistry enabled ? 
if (aertyp[3:5]!='00'):
    f.write('./aereosol \n') 
      
if (dattyp=='ERAHI'): 
    f.writelines('echo: not yet ported ERAHI \nexit \n')

f.writelines('./icbc \n') 
f.close() 
