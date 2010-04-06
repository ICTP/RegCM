!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine mappol(xlon,xlat,xmap,coriol,iy,jx,clon,clat,delx,idot)
      implicit none
!
! Dummy arguments
!
      real(4) :: clat , clon , delx
      integer :: idot , iy , jx
      real(4) , dimension(iy,jx) :: coriol , xlat , xlon , xmap
      intent (in) clat , clon , delx , idot , iy , jx
      intent (out) coriol , xlon , xmap
      intent (inout) xlat
!
! Local variables
!
      real(4) :: aa , cell , cell1 , cell2 , cntri , cntrj , d2r , flp ,&
               & flpp , pi , pole , psi1 , psix , psx , r , x , xcntr , &
               & xn , xomeg2 , xomega , y , ycntr
      integer :: i , ii1 , j , jj1
!
!     XN IS CONE FACTOR FOR THE PROJECTION (FROM PROGRAM TERRAIN).
!     PSI1 IS THE COLATITUDE OF TRUELAT 1, IN RADIANS.
!     PI IS PI.
!
!---------------------------------------------------------------------
 
!     COMPUTE LATS, LONS, AND MAP-SCALE FACTORS FOR
!     LAMBERT CONFORMAL MAP CENTERED AT CLON,CLAT. TRUE AT 30.N AND
 
!     60.N. IY IS NUMBER OF N-S POINTS.  JX IS NUMBER OF E-W POINTS.
!     CLON, CLAT IS LAT, LON OF CENTER OF GRID (DEGREES EAST, NORTH).
!     DELX IS GRID SPACING IN METERS.
!     ALWAYS FOR CROSS GRID.
 
!
      aa = 6.371229E6
      xn = 1.0
      pi = atan(1.)*4.
      d2r = pi/180.
!
      pole = 90.
      psi1 = 30.
      psi1 = psi1*d2r
      if ( clat<0.0 ) then
        pole = -90.0
        psi1 = -30.
        psi1 = psi1*d2r
      end if
      cntrj = float(jx+idot)/2.
      cntri = float(iy+idot)/2.
!
      psx = (pole-clat)*d2r
      cell = aa*sin(psx)/xn
      cell2 = (1.+cos(psi1))/(1.+cos(psx))
      r = cell*(cell2)**xn
      xcntr = 0.
      ycntr = -r
!
      ii1 = iy
      jj1 = jx
      do i = 1 , ii1
        y = ycntr + (i-cntri)*delx
        do j = 1 , jj1
          x = xcntr + (j-cntrj)*delx
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
          cell = r/aa
          cell1 = cell/(1.0+cos(psi1))
          cell2 = atan(cell1)
          psx = 2.*cell2/d2r
          xlat(i,j) = pole - psx
          psix = psx*d2r
          xmap(i,j) = ((1.0+cos(psi1))/(1.0+cos(psix)))**xn
        end do
      end do
 
      if ( idot==1 ) then
        xomega = 7.2722E-5            ! ANG. ROT OF EARTH IN S**-1
        d2r = atan(1.)/45.          ! CONVERT DEGREES TO RADIANS
        xomeg2 = 2.*xomega
        do i = 1 , iy
          do j = 1 , jx
            coriol(i,j) = xomeg2*sin(xlat(i,j)*d2r)
          end do
        end do
      end if
      end subroutine mappol
