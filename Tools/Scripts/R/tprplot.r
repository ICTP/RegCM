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

library("ncdf")
library("PCICt")
library("fields")

surface = open.ncdf("/scratch/ggiulian/test/output/verysmall_SRF.1999050100.nc")
lonmat  = get.var.ncdf(nc=surface,varid="xlon")
latmat  = get.var.ncdf(nc=surface,varid="xlat")

timearr   = get.var.ncdf(nc=surface,varid="time")
timeunits = att.get.ncdf(nc=surface,varid="time",attname="units")
calendar  = att.get.ncdf(nc=surface,varid="time",attname="calendar")

tpr       = get.var.ncdf(nc=surface,varid="pr")
tprunits  = att.get.ncdf(nc=surface,varid="pr",attname="units")

seconds_per_hour <- 3600.0
origin <- as.PCICt(substring(timeunits$value,13,31),cal=calendar$value)
times <- origin + timearr*seconds_per_hour

inds = (1:dim(timearr))
targettime = as.PCICt("1999-05-06 12:00:00",cal=calendar$value)
tind = inds[times == targettime]

image.plot(lonmat, latmat,tpr[,,tind])

close.ncdf(surface)
