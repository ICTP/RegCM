      subroutine rotmer(xlon,xlat,xmap,coriol,iy,jx,clon,clat,pollon,   &
                      & pollat,ds,idot)
      implicit none
!
! Dummy arguments
!
      real(4) :: clat , clon , ds , pollat , pollon
      integer :: idot , iy , jx
      real(4) , dimension(iy,jx) :: coriol , xlat , xlon , xmap
      intent (in) clat , clon , ds , idot , iy , jx
      intent (out) coriol , xlon , xmap
      intent (inout) xlat
!
! Local variables
!
      real(4) :: a , cntri , cntrj , d2r , ddeg , fai , r2d , x , xoff ,&
               & xomega , xomega2 , xr , y , yoff , yr
      integer :: i , j
!
!---------------------------------------------------------------------
!     COMPUTE LATS, LONS, MAP-SCALE FACTORS, AND CORIOLIS PARAMETER FOR
!     ROTATED POLE MAP CENTERED AT CLON,CLAT. ORIGIN OF ROTATED GRID IS
!     GIVEN BY POLLON AND POLLAT.
!     IMX,JMX,KSV,KSH AND LSIG2 MUST CORRESPOND TO IY,JMAX,KSTRPV,KSTRPH
!     AND NVERT2 IN THE MASTER INPUT FILE; OTHERWISE THE PROGRAM WILL
!     ABORT:LMX MUST BE THE MAXIMUM NUMBER OF LEVELS
!     (PRESSURE--OR--SIGMA) IMAXN AND IMXC ARE NESTED AND NON-EXPANDED
!     GRID DIMENSIONS. IMXC IS EQUAL TO IMX IF YOU ARE NOT USING THE
!     EXPANDED GRID. SAME FOR J.
!
!
      xomega = 7.2722E-5                       ! ANG. ROT OF EARTH IN S**-1
      d2r = atan(1.)/45.                     ! CONVERT DEGREES TO RADIANS
      r2d = 1./d2r                         ! CONVERT RADIANS TO DEGREES
      a = 6371229.                             ! RADIUS OF EARTH IN METERS
!-----CENTER OF GRID
      cntrj = (jx+idot)/2.
      cntri = (iy+idot)/2.
 
      ddeg = ds*r2d/a                        ! GRID SPACING IN DEGREES
      xoff = clon - pollon
      yoff = clat - pollat
!-----CALCULATE X AND Y POSITIONS OF GRID
      do i = 1 , iy
        do j = 1 , jx
          xr = xoff + (j-cntrj)*ddeg
          yr = yoff + (i-cntri)*ddeg
!-----NOW CALCULATE LAT AND LON OF THIS POINT
!-----    ROTATE COORDINATES BACK TO NONRATED POLE
          call rot2nrot(xr,yr,pollon,pollat,x,y)
          xlon(i,j) = x
          xlat(i,j) = y
          fai = d2r*yr
          xmap(i,j) = 1.0/cos(fai)
        end do
      end do
 
      if ( idot==1 ) then
        xomega2 = 2.*xomega
        do i = 1 , iy
          do j = 1 , jx
            coriol(i,j) = xomega2*sin(xlat(i,j)*d2r)
          end do
        end do
      end if
 
      end subroutine rotmer
