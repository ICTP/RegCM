C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine zenith(dec,alat,fjd,rsun,coszrs,frac,imax)
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c ** used by radiation package in ccm to coszrs  *****not yet used
c                 here but could be called from albedo
c
c     this routine computes the zenith angle at each of 'imax'
c     equally-spaced longitudes at latitude 'alat', for a diurnally
c     varying sun.
c
c     input:  dec    - declination of sun in radians. computed in
c                      'solar'.
c             alat   - latitude in radians.
c             fjd    - fractional julian date in days.  the code conven-
c                      tion (see 'compjd') is that fjd=0 at noon  at
c                      greenwich meridian (1200 gmt), so the hour angle
c                      at longitude 'lon' is ha=lon+twopi*fjd.
c             rsun   - distance to sun, normalized to the mean distance.
c                      computed in 'solar'.
c             imax   - number of equally spaced longitude values.
c                      it is assumed that i=1 and i=imax+1 are both
c                      at the greenwich meridian (lon=0).
c     output: coszrs - cos(z)/ where 'z' is the zenith angle at
c                      longitude 'i', latitude 'alat', and time 'fjd'.
c                           cos(z)=sin(alat)*sin(dec)+
c                                  cos(alat)*cos(dec)*cos(ha)
c                      it is assumed that the
c                      annual mean solar constant is used elsewhere
c                      in determining the solar flux.
c                      the 1/r**2 decrease of the solar flux appears
c                      in subroutine radiatn as eccf
c             frac   - not used in diurnal mode: set to 1.  in the
c                      average insolation mode, 'frac' is the daylight
c                      fraction at the point (see 'zenith'); to lowest
c                      order it should be independent of longitude.
c
c     ***************
c
      implicit none
      real(kind=8)  dec,alat,fjd,rsun
      integer imax
      real(kind=8)  coszrs(imax),frac(imax)
c
      real(kind=8)  twopi
      data twopi /6.283185308/
      real(kind=8)  ss,cc,rs,eccf,dlon,tpifjd,ha,cosz
      integer i
c
c***********************************************************************
c
      ss=dsin(alat)*dsin(dec)
      cc=dcos(alat)*dcos(dec)
      rs=1./rsun**2
      eccf=rs
      rs=1.0
      dlon=twopi/imax
      tpifjd=twopi*fjd
      do 20 i=1,imax
         frac(i)=1.
         ha=(i-1)*dlon+tpifjd
c       -if cosz is negative, the sun is below the horizon.
         cosz=dmax1(0.d0,ss+cc*dcos(ha))
         coszrs(i)=cosz
20    continue

      return
      end
