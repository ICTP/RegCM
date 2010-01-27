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
C
      subroutine zenith(dec,alat,fjd,rsun,coszrs,frac,imax)
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C ** used by radiation package in ccm to coszrs  *****not yet used
C                 here but could be called from albedo
C
C     this routine computes the zenith angle at each of 'imax'
C     equally-spaced longitudes at latitude 'alat', for a diurnally
C     varying sun.
C
C     input:  dec    - declination of sun in radians. computed in
C                      'solar'.
C             alat   - latitude in radians.
C             fjd    - fractional julian date in days.  the code conven-
C                      tion (see 'compjd') is that fjd=0 at noon  at
C                      greenwich meridian (1200 gmt), so the hour angle
C                      at longitude 'lon' is ha=lon+twopi*fjd.
C             rsun   - distance to sun, normalized to the mean distance.
C                      computed in 'solar'.
C             imax   - number of equally spaced longitude values.
C                      it is assumed that i=1 and i=imax+1 are both
C                      at the greenwich meridian (lon=0).
C     output: coszrs - cos(z)/ where 'z' is the zenith angle at
C                      longitude 'i', latitude 'alat', and time 'fjd'.
C                           cos(z)=sin(alat)*sin(dec)+
C                                  cos(alat)*cos(dec)*cos(ha)
C                      it is assumed that the
C                      annual mean solar constant is used elsewhere
C                      in determining the solar flux.
C                      the 1/r**2 decrease of the solar flux appears
C                      in subroutine radiatn as eccf
C             frac   - not used in diurnal mode: set to 1.  in the
C                      average insolation mode, 'frac' is the daylight
C                      fraction at the point (see 'zenith'); to lowest
C                      order it should be independent of longitude.
C
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
      implicit none
C
      real(kind=8) dec , alat , fjd , rsun
      integer imax
      real(kind=8) coszrs(imax) , frac(imax)
C
      real(kind=8) twopi
      data twopi/6.283185308/
      real(kind=8) ss , cc , rs , dlon , tpifjd , ha , cosz
      integer i
C
C***********************************************************************
C
      ss = dsin(alat)*dsin(dec)
      cc = dcos(alat)*dcos(dec)
      rs = 1.0/rsun**2
      rs = 1.0
      dlon = twopi/imax
      tpifjd = twopi*fjd
      do i = 1 , imax
        frac(i) = 1.0
        ha = (i-1)*dlon + tpifjd
C       if cosz is negative, the sun is below the horizon.
        cosz = dmax1(0.D0,ss+cc*dcos(ha))
        coszrs(i) = cosz
      end do
C
      end subroutine zenith
