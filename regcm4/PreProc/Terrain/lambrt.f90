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

      subroutine lambrt(xlon,xlat,smap,coriol,iy,jx,clon,clat,ds,idot,  &
                      & xn,truelatl,truelath)
      implicit none
!
! Dummy arguments
!
      real(4) :: clat , clon , ds , truelath , truelatl , xn
      integer :: idot , iy , jx
      real(4) , dimension(iy,jx) :: coriol , smap , xlat , xlon
      intent (in) clat , clon , ds , idot , iy , jx , truelath ,        &
                & truelatl
      intent (out) coriol , smap , xlon
      intent (inout) xlat , xn
!
! Local variables
!
      real(8) :: er , cell , cell1 , cell2 , cntri , cntrj , d2r , flp ,&
               & flpp , omega2 , pi , pole , psi1 , psix , psx , r ,    &
               & xsign , truelat1 , truelat2 , x , xcntr , y , ycntr
      integer :: i , j
!
!     CLAT IS THE CENTRAL LATITUDE OF THE COARSE DOMAIN.
!     CLON IS THE CENTRAL LONGITUDE OF THE COARSE DOMAIN.
!
!     THIS ROUTINE CALCULATES MESO MAP(LAT,LONG,CORIOLIS,MAP SCALE)
!     FOR LAMBERT CONFORMAL PROJECTION
!
!     IX IS THE I DIMENSION FOR THIS DOMAIN.
!     JX IS THE J DIMENSION FOR THIS DOMAIN.
!     IDOT IS ICROSS ( = 1)  OR IDOT ( = 0).
!
!---------------------------------------------------------------------
!
!
!     XN IS CONE FACTOR FOR THE PROJECTION (FROM PROGRAM TERRAIN).
!     PSI1 IS THE COLATITUDE OF TRUELAT 1, IN RADIANS.
!     PI IS PI.
!
!---------------------------------------------------------------------
!
      er = 6.371229E6
      if ( clat<0. ) then
        xsign = -1.       ! SOUTH HEMESPHERE
      else
        xsign = 1.        ! NORTH HEMESPHERE
      end if
      pole = xsign*90.0
      pi = 3.1415926535897932384626433832795029D+00
      d2r = pi/180.
 
      truelat1 = truelath
      truelat2 = truelatl
      if ( abs(truelat1-truelat2)>1.E-1 ) then
        xn = (dlog10(cos(truelat1*d2r))-dlog10(cos(truelat2*d2r)))      &
           & /(dlog10(tan((45.0-xsign*truelat1/2.0)*d2r))                &
           & -dlog10(tan((45.0-xsign*truelat2/2.0)*d2r)))
      else
        xn = xsign*sin(truelat1*d2r)
      end if
!     XN=0.716
 
      psi1 = 90. - xsign*truelat1
      if ( clat<0. ) psi1 = -psi1
!
      psi1 = psi1*d2r
      omega2 = 7.2921159D-05*2.0D0
      cntrj = (jx+idot)/2.
      cntri = (iy+idot)/2.
!
      psx = (pole-clat)*d2r
      cell = er*sin(psi1)/xn
      write (*,*) 'PSX,PSI1 = ' , psx , psi1
      cell2 = (tan(psx/2.))/(tan(psi1/2.))
      r = cell*(cell2)**xn
      xcntr = 0.
      ycntr = -r
!
      do j = 1 , jx
        x = xcntr + (j-cntrj)*ds
        do i = 1 , iy
          y = ycntr + (i-cntri)*ds
          r = sqrt(x*x+y*y)
          if ( y==0. ) then
            if ( x>=0. ) then
              flp = 90.*d2r
            else
              flp = -90.*d2r
            end if
          else if ( clat<0.0 ) then
            flp = atan2(x,y)
          else
            flp = atan2(x,-y)
          end if
          flpp = flp/xn/d2r + clon
!         IF (FLPP.GT.180.0) FLPP = FLPP-360.0
!         IF (FLPP.LT.-180.0) FLPP = FLPP+360.0
          xlon(i,j) = flpp
          if ( clat<0.0 ) r = -r
          cell = r*xn/(er*sin(psi1))
          cell1 = tan(psi1/2.)*cell**(1./xn)
          cell2 = atan(cell1)
          psx = 2.*cell2/d2r
          xlat(i,j) = pole - psx
          psix = psx*d2r
          smap(i,j) = (sin(psi1)/sin(psix))                             &
                    & *((tan(psix/2.)/tan(psi1/2.))**xn)
        end do
      end do
      if ( idot==1 ) then
        do i = 1 , iy
          do j = 1 , jx
            coriol(i,j) = omega2*sin(xlat(i,j)*d2r)
          end do
        end do
      end if
      end subroutine lambrt
