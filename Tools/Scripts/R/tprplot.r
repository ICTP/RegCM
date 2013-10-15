#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
#
#    This file is part of ICTP RegCM.
#
#    ICTP RegCM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY# without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

library("ncdf")
library("udunits2")
library("fields")
library("chron")

surface = open.ncdf("/scratch/ggiulian/test/output/verysmall_SRF.1999050100.nc")
lonmat  = get.var.ncdf(nc=surface,varid="xlon")
latmat  = get.var.ncdf(nc=surface,varid="xlat")

timearr   = get.var.ncdf(nc=surface,varid="time")
timeunits = att.get.ncdf(nc=surface,varid="time",attname="units")

tpr       = get.var.ncdf(nc=surface,varid="pr")
tprunits  = att.get.ncdf(nc=surface,varid="pr",attname="units")

ts = ud.convert(timearr,timeunits$value,"seconds since 1970-01-01 00:00:00 GMT")
a = as.POSIXct(ts,origin="1970-01-01 00:00:00 GMT",tz="GMT")

inds = (1:dim(timearr))
targettime = chron(dates("1999-05-05",format="y-m-d"), times("21:00:00"))
tind = inds[as.chron(a) == targettime]

image.plot(lonmat, latmat,tpr[,,tind])

close.ncdf(surface)
