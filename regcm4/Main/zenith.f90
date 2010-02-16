!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      subroutine zenith(dec,alat,fjd,coszrs,frac,imax)
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! ** used by radiation package in ccm to coszrs  *****not yet used
!                 here but could be called from albedo
!
!     this routine computes the zenith angle at each of 'imax'
!     equally-spaced longitudes at latitude 'alat', for a diurnally
!     varying sun.
!
!     input:  dec    - declination of sun in radians. computed in
!                      'solar'.
!             alat   - latitude in radians.
!             fjd    - fractional julian date in days.  the code conven-
!                      tion (see 'compjd') is that fjd=0 at noon  at
!                      greenwich meridian (1200 gmt), so the hour angle
!                      at longitude 'lon' is ha=lon+twopi*fjd.
!             imax   - number of equally spaced longitude values.
!                      it is assumed that i=1 and i=imax+1 are both
!                      at the greenwich meridian (lon=0).
!     output: coszrs - cos(z)/ where 'z' is the zenith angle at
!                      longitude 'i', latitude 'alat', and time 'fjd'.
!                           cos(z)=sin(alat)*sin(dec)+
!                                  cos(alat)*cos(dec)*cos(ha)
!                      it is assumed that the
!                      annual mean solar constant is used elsewhere
!                      in determining the solar flux.
!                      the 1/r**2 decrease of the solar flux appears
!                      in subroutine radiatn as eccf
!             frac   - not used in diurnal mode: set to 1.  in the
!                      average insolation mode, 'frac' is the daylight
!                      fraction at the point (see 'zenith'); to lowest
!                      order it should be independent of longitude.
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      implicit none
!
! Dummy arguments
!
      integer , intent (in) :: imax
      real(kind=8) , intent (in) :: alat , dec , fjd
      real(kind=8) , intent (out) , dimension(imax) :: coszrs , frac
!
! Local variables
!
      real(8) :: cc , cosz , dlon , ha , ss , tpifjd , twopi
      integer :: i
!
      data twopi / 6.28318530717958647692 /
!
!***********************************************************************
!
      ss = dsin(alat)*dsin(dec)
      cc = dcos(alat)*dcos(dec)
      dlon = twopi/imax
      tpifjd = twopi*fjd
      do i = 1 , imax
        frac(i) = 1.0
        ha = (i-1)*dlon + tpifjd
!       if cosz is negative, the sun is below the horizon.
        cosz = dmax1(0.D0,ss+cc*dcos(ha))
        coszrs(i) = cosz
      end do
!
      end subroutine zenith
