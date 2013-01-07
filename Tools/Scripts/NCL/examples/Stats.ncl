;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
;
;    This file is part of ICTP RegCM.
;
;    ICTP RegCM is free software: you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation, either version 3 of the License, or
;    (at your option) any later version.
;
;    ICTP RegCM is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
;
;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

; Credit to Noah Diffenbaugh http://es.ucsc.edu/~msnyder/NCL/ncldata.html

;********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;********************************************

begin

; user-specied buffer on each side of model domain
; four variables are used to provide maximum flexibility,
; though all values may be the same
; depending on user intention.
east_buffer = 8
west_buffer = 8
north_buffer = 8
south_buffer = 8

diri = "/home/graziano/test/output/"
fili = "CAS_ATM.1990010100.nc"

dirout = "./"
filout = "CAS_ATM.1990010100_stats.nc"

f = addfile(diri+fili,"r")

fout = addfile(dirout+filout,"c")

; Read in variables
ts_in     = f->ta(:,17,:,:)
lat2d_in  = f->xlat
lon2d_in  = f->xlon

; Get dimensions from input variable
dimvar_in = dimsizes(ts_in)
ntim_in   = dimvar_in(0)
jlat_in   = dimvar_in(1)
ilon_in   = dimvar_in(2)

;trim buffer (buffer size is set by the user above)
ts = ts_in(:,south_buffer:jlat_in-north_buffer-1,west_buffer:ilon_in-north_buffer-1)
lat2d = lat2d_in(south_buffer:jlat_in-north_buffer-1,west_buffer:ilon_in-north_buffer-1)
lon2d = lon2d_in(south_buffer:jlat_in-north_buffer-1,west_buffer:ilon_in-north_buffer-1)

;get dimensions of trimmed variable
dimvar = dimsizes(ts)
ntim = dimvar(0)
jlat = dimvar(1)
ilon = dimvar(2)

;use ncl functions to compute simple statistics on the data
;whole domain for all time planes
d_avg = avg(ts)						;domain average over all time planes
d_max = max(ts)						;domain max over all time planes
d_min = min(ts)						;domain min over all time planes
d_stddev = stddev(ts)					;domain standard deviation over all time planes
d_var = variance(ts)					;domain variance over all time planes

;whole domain at each time plane
;allocate memory for variables
t_avg = new((/ntim/),typeof(ts))
t_max = new((/ntim/),typeof(ts))
t_min = new((/ntim/),typeof(ts))
t_stddev = new((/ntim/),typeof(ts))
t_var = new((/ntim/),typeof(ts))

do i=0,ntim-1
	t_avg(i) = avg(ts(i,:,:))				;domain average for each time plane
	t_max(i) = max(ts(i,:,:))				;domain max for each time plane
	t_min(i) = min(ts(i,:,:))				;domain min for each time plane
	t_stddev(i) = stddev(ts(i,:,:))			;domain standard deviation for each time plane
	t_var(i) = variance(ts(i,:,:))			;domain variance for each time plane
end do

;iy x jx for all time planes
time_avg = dim_avg(ts(iy|:,jx|:,time|:))		;iy x jx average at each grid point
time_med = dim_median(ts(iy|:,jx|:,time|:))		;iy x jx median at each grid point
time_min = dim_min(ts(iy|:,jx|:,time|:))		;iy x jx min at each grid point 
time_max = dim_max(ts(iy|:,jx|:,time|:))		;iy x jx max at each grid point
time_stddev = dim_stddev(ts(iy|:,jx|:,time|:))	;iy x jx standard deviation at each grid point
time_var = dim_variance(ts(iy|:,jx|:,time|:))		;iy x jx variance at each grid point

;assign some meta-data so it's easier to tell what the variables are in the netcdf file
d_avg@long_name = "domain average over all time planes"
d_max@long_name = "domain max over all time planes"
d_min@long_name	= "domain min over all time planes"
d_stddev@long_name = "domain standard deviation over all time planes"
d_var@long_name = "domain variance over all time planes"

t_avg@long_name = "domain average for each time plane"
t_max@long_name = "domain max for each time plane"
t_min@long_name = "domain min for each time plane"
t_stddev@long_name = "domain standard deviation for each time plane"
t_var@long_name = "domain variance for each time plane"

time_avg@long_name = "lat x lon average at each grid point"
time_med@long_name = "lat x lon median at each grid point"
time_min@long_name = "lat x lon min at each grid point"
time_max@long_name = "lat x lon max at each grid point"
time_stddev@long_name = "lat x lon standard deviation at each grid point"
time_var@long_name = "lat x lon variance at each grid point"

;assign coordinate variables
;scalar
d_avg!0 = "avg"
d_max!0 = "avg"
d_min!0 = "avg"
d_stddev!0 = "avg"
d_var!0 = "avg"

;time
t_avg!0 = "time"
t_max!0 = "time"
t_min!0 = "time"
t_stddev!0 = "time"
t_var!0 = "time"

;lat x lon
time_avg!0 = "iy"
time_avg!1 = "jx"

time_med!0 = "iy"
time_med!1 = "jx"

time_min!0 = "iy"
time_min!1 = "jx"

time_max!0 = "iy"
time_max!1 = "jx"

time_stddev!0 = "iy"
time_stddev!1 = "jx"

time_var!0 = "iy"
time_var!1 = "jx"

;output variables to file
fout->d_avg = d_avg
fout->d_max = d_max
fout->d_min = d_min
fout->d_stddev = d_stddev
fout->d_var = d_var

fout->t_avg = t_avg
fout->t_max = t_max
fout->t_min = t_min
fout->t_stddev = t_stddev
fout->t_var = t_var

fout->time_avg = time_avg
fout->time_med = time_med
fout->time_max = time_max
fout->time_min = time_min
fout->time_stddev = time_stddev
fout->time_var = time_var

fout->xlat = lat2d
fout->xlon = lon2d

end
