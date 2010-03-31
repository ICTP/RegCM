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

      subroutine clm3grid1(nlon,nlat,nlev,ntim,glon1,glon2,glat1,glat2, &
                         & xlonmin,xlonmax,xlatmin,xlatmax,glon,glat,   &
                         & istart,icount)
      implicit none
!
! Dummy arguments
!
      real(4) :: glat1 , glat2 , glon1 , glon2 , xlatmax , xlatmin ,    &
            & xlonmax , xlonmin
      integer :: nlat , nlev , nlon , ntim
      real(4) , dimension(nlat) :: glat
      real(4) , dimension(nlon) :: glon
      integer , dimension(4) :: icount , istart
      intent (in) glat1 , glat2 , nlat , nlev , nlon , ntim
      intent (out) glat , glon , icount , istart
      intent (inout) glon1 , glon2 , xlatmax , xlatmin , xlonmax ,      &
                   & xlonmin
!
! Local variables
!
      integer :: corrlatn , corrlats , i , ilatmax , ilatmin , ilonmax ,&
               & ilonmin , j
      real(4) :: dlat , dlon
 
!     dlon = 360./nlon
!     dlat = 180./nlat

!     ABT added to get mksrf dependent resolution
      dlon = (glon2+abs(glon1)+0.5)/nlon
      dlat = (glat2+abs(glat1)+0.5)/nlat
!     ABT correction terms in case the grid is not from 90S to 90N 
      corrlatn = 90 - nint(glat2)
      corrlats = -90 - nint(glat1)
 
      if ( glon1>=0. ) then
        glon1 = glon1 - 180.
        glon2 = glon2 - 180.
      end if
      do i = 1 , nlon
        glon(i) = glon1 + dlon*float(i-1)
      end do
      do j = 1 , nlat
        glat(j) = glat1 + dlat*float(j-1)
      end do
 
      xlonmin = max(xlonmin-dlon,glon1)
      xlonmax = min(xlonmax+dlon,glon2)
      xlatmin = max(xlatmin-dlat,glat1)
      xlatmax = min(xlatmax+dlat,glat2)
      ilonmin = max(min(nint((glon2+xlonmin)/dlon)-1,nlon),1)
      ilonmax = max(min(nint((glon2+xlonmax)/dlon)+1,nlon),1)
!     abt ilatmin = max(min(nint((glat2+xlatmin)/dlat)-1,nlat),1)
!     abt ilatmax = max(min(nint((glat2+xlatmax)/dlat)+1,nlat),1)
      ilatmin = max(min(nint((glat2+xlatmin+corrlatn+corrlats)/dlat)-1, &
              & nlat),1)   ! ABT added corrlat terms
      ilatmax = max(min(nint((glat2+xlatmax+corrlatn+corrlats)/dlat)+1, &
              & nlat),1)   ! ABT added corrlat terms
      istart(1) = ilonmin
      icount(1) = ilonmax - ilonmin + 1
      istart(2) = ilatmin
      icount(2) = ilatmax - ilatmin + 1
      istart(3) = 1
      icount(3) = nlev
      istart(4) = 1
      icount(4) = ntim
 
      end subroutine clm3grid1
