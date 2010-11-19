#!/usr/bin/python

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP RegCM.
#
#    ICTP RegCM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Module that contains NetCDF related functions
import subprocess

# compares two RegCM outputs in NetCDF format
# requirements: subprocess
def compare_nc_file(filename,refname,varname):

    #print filename
    #print refname
    
    try :
        p_1 = subprocess.Popen("ncdiff -v "+varname+" "+filename+" "+refname+" temp.nc",stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
    except OSError :
        print "Could not run ncdiff!"
        output,error = p_1.communicate()
        return output

    if p_1.wait() == 0 :
        try :
            p_2 = subprocess.Popen("ncwa -y rms temp.nc rms.nc",stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
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
