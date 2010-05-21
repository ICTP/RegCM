* This is a script for displaying moisture convergence 
* Written by Michael Maxwell 
*
* rh    = Relative Humidity in % 
* t     = Temp at *set level in degrees Kelvin
* tc    = Temp in degrees C
* td    = Dewpoint at *set level in degrees C
* e     = Vapor pressure	
* mixr  = Mixing ratio
* u     = U-wind in m/s
* v     = V-wind in m/s
* mconv = moisture convergence/divergence. convergence is positive and divergence is negative.

'tc = (t-273.16)'
'td = tc-((14.55+0.114*tc)*(1-0.01*rh) + pow((2.5+0.007*tc)*(1-0.01*rh),3) + (15.9+0.117*tc)*pow((1-0.01*rh),14))'
'vapr = 6.112*exp((17.67*td)/(td+243.5))'
'e = vapr*1.001+(lev-100)/900*0.0034'
'mixr = 0.62197*(e/(lev-e))*1000'
'mconv = (-1)*hdivg(u*mixr,v*mixr)*1e4'
'd mconv'