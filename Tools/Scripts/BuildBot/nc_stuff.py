#!/usr/bin/python

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP RegCM.
#    
#    Use of this source code is governed by an MIT-style license that can
#    be found in the LICENSE file or at
#
#         https://opensource.org/licenses/MIT.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Module that contains NetCDF related functions
import subprocess,os

# compares two RegCM outputs in NetCDF format
# requirements: subprocess
def compare_nc_file(filename,refname,varname):

    #print filename
    #print refname
    
    try :
        p_1 = subprocess.Popen("ncdiff -v "+varname+" "+filename+" "+refname+" temp.nc",stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
        while not p_1.poll():
            if p_1.returncode is not None: break
            time.sleep(10)
            print "   ...ncdiff: still alive..."
        
    except OSError :
        print "Could not run ncdiff!"
        output,error = p_1.communicate()
        return output

    if p_1.wait() == 0 :
        try :
            p_2 = subprocess.Popen("ncwa -y rms temp.nc rms.nc",stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
            while not p_2.poll():
                if p_2.returncode is not None: break
                time.sleep(10)
                print "   ...ncwa: still alive..."
        except OSError :
           print "Could not run ncwa!"
           output,error = p_2.communicate()
           return output  
    else :
        print "Step 1 failed!"
        output,error = p_1.communicate()
        return output

    if p_2.wait() == 0 :
        os.remove("temp.nc")
        try :
            p_3 = subprocess.Popen('ncks -H -s "%g\n" -v '+varname+' rms.nc',stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
            while not p_3.poll():
                if p_3.returncode is not None: break
                time.sleep(10)
                print "   ...ncks: still alive..."
        except OSError :
            print "Could not run ncks!"
            output,error = p_3.communicate()
            return output
    else :
        print "Step 2 failed!"
        output,error = p_2.communicate()
        return output    
    
    if p_3.wait() != 0:
        print "Step 3 failed!"
        output,error = p_3.communicate()
        return output+error
    else:
        os.remove("rms.nc")
        output,error = p_3.communicate()
        
    return output
