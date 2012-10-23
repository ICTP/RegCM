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

module mod_uvrot

  use mod_intkinds
  use mod_realkinds
  use mod_constants

  contains

  subroutine uvrot4(u,v,dlon,dlat,clon,clat,gridfc,jx,iy,ll,pollon, &
                    pollat,lgtype)
  implicit none
!
  real(rk8) :: clat , clon , gridfc , pollat , pollon
  integer(ik4) :: iy , jx , ll
  character(6) :: lgtype
  real(rk4) , dimension(jx,iy) :: dlat , dlon
  real(rk4) , dimension(jx,iy,ll) :: u , v
  intent (in) clat , clon , dlat , dlon , gridfc , iy , jx ,        &
              lgtype , ll , pollat , pollon
  intent (inout) u , v
!
  real(rk8) :: cosdel , d , polcphi , pollam , polphi , polsphi ,    &
          sindel , us , vs , x , xc , xs , zarg1 , zarg2 , znorm ,  &
          zphi , zrla , zrlap
  integer(ik4) :: i , j , l
!
!     CHANGE U AND V FROM TRUE (N,E) TO MAP VALUES (X,Y)
!
!     FOR ROTATED MERCATOR PROJECTION
!UVUSVS   -   ROTATES THE TWO WINDCOMPONENTS U AND V AT POINT
!     DLON,DLAT TO THE WINDCOMPONENTS US AND VS IN A
!     ROTATED POLE GRID WHERE THE ORIGIN IS LOCATED
!     AT POLLON,POLLAT
!**   CALL  :   CALL UVUSVS(U,V,US,VS,DLON,DLAT,POLLON,POLLAT)
!**   AUTHOR:   D.MAJEWSKI
!
  if ( lgtype == 'NORMER' ) then
    return
  end if
  if ( lgtype == 'ROTMER' ) then
    if ( pollat > d_zero ) then
      pollam = pollon + 180.0D0
      polphi = 90.0D0 - pollat
    else
      polphi = 90.0D0 + pollat
      pollam = pollon
    end if
    if ( pollam > 180.0D0 ) pollam = pollam - 360.0D0
 
    polcphi = cos(degrad*polphi)
    polsphi = sin(degrad*polphi)
 
    do j = 1 , iy
      do i = 1 , jx
        zphi = dble(dlat(i,j))*degrad
        zrla = dble(dlon(i,j))*degrad
        if ( dlat(i,j) > 89.999999 ) zrla = d_zero
        zrlap = pollam*degrad - zrla
        zarg1 = polcphi*sin(zrlap)
        zarg2 = polsphi*cos(zphi) - polcphi*sin(zphi)*cos(zrlap)
        znorm = d_one/sqrt(zarg1**2+zarg2**2)
        sindel = zarg1*znorm
        cosdel = zarg2*znorm
        do l = 1 , ll
          us = dble(u(i,j,l))*cosdel - dble(v(i,j,l))*sindel
          vs = dble(u(i,j,l))*sindel + dble(v(i,j,l))*cosdel
          u(i,j,l) = real(us)
          v(i,j,l) = real(vs)
        end do
      end do
    end do
  else
    do j = 1 , iy
      do i = 1 , jx
        if ( (clon >= d_zero .and. dlon(i,j) >= 0.0) .or.   &
             (clon < d_zero  .and. dlon(i,j) < 0.0) ) then
          x = (clon-dble(dlon(i,j)))*degrad*gridfc
        else if ( clon >= d_zero ) then
          if ( abs(clon-(dble(dlon(i,j))+360.0D0)) < &
               abs(clon-dble(dlon(i,j))) ) then
            x = (clon-(dble(dlon(i,j))+360.0D0))*degrad*gridfc
          else
            x = (clon-dble(dlon(i,j)))*degrad*gridfc
          end if
        else if ( abs(clon-(dble(dlon(i,j))-360.0D0)) < &
                  abs(clon-dble(dlon(i,j))) ) then
          x = (clon-(dble(dlon(i,j))-360.0D0))*degrad*gridfc
        else
          x = (clon-dble(dlon(i,j)))*degrad*gridfc
        end if
        xs = sin(x)
        xc = cos(x)
        if ( clat >= d_zero ) then
          do l = 1 , ll
            d = dble(v(i,j,l))*xs + dble(u(i,j,l))*xc
            v(i,j,l) = real(dble(v(i,j,l))*xc - dble(u(i,j,l))*xs)
            u(i,j,l) = real(d)
          end do
        else
          do l = 1 , ll
            d = -dble(v(i,j,l))*xs + dble(u(i,j,l))*xc
            v(i,j,l) = real(dble(v(i,j,l))*xc + dble(u(i,j,l))*xs)
            u(i,j,l) = real(d)
          end do
        end if
      end do
    end do
  end if
  end subroutine uvrot4
!
!-----------------------------------------------------------------------
!
  subroutine uvrot4nx(u,v,dlon,dlat,clon,clat,gridfc,jx,iy,ll, &
                      pollon,pollat,lgtype)
  implicit none
!
  real(rk8) :: clat , clon , gridfc , pollat , pollon
  integer(ik4) :: iy , jx , ll
  character(6) :: lgtype
  real(rk4) , dimension(jx,iy) :: dlat , dlon
  real(rk4) , dimension(jx,iy,ll) :: u , v
  intent (in) clat , clon , dlat , dlon , gridfc , iy , jx ,  &
              lgtype , ll , pollat , pollon
  intent (inout) u , v
!
  real(rk8) :: cosdel , d , polcphi , pollam , polphi , polsphi ,   &
          sindel , us , vs , x , xc , xs , zarg1 , zarg2 , znorm , &
          zphi , zrla , zrlap
  integer(ik4) :: i , j , l
!
!     CHANGE U AND V FROM TRUE (N,E) TO MAP VALUES (X,Y)
!
!     FOR ROTATED MERCATOR PROJECTION
!UVUSVS   -   ROTATES THE TWO WINDCOMPONENTS U AND V AT POINT
!     DLON,DLAT TO THE WINDCOMPONENTS US AND VS IN A
!     ROTATED POLE GRID WHERE THE ORIGIN IS LOCATED
!     AT POLLON,POLLAT
!**   CALL  :   CALL UVUSVS(U,V,US,VS,DLON,DLAT,POLLON,POLLAT)
!**   AUTHOR:   D.MAJEWSKI
!
  if ( lgtype == 'ROTMER' ) then
    if ( pollat > d_zero ) then
      pollam = pollon + 180.0D0
      polphi = 90.0D0 - pollat
    else
      polphi = 90.0D0 + pollat
      pollam = pollon
    end if
    if ( pollam > 180.0D0 ) pollam = pollam - 360.0D0
 
    polcphi = cos(degrad*polphi)
    polsphi = sin(degrad*polphi)
 
    do j = 1 , iy
      do i = 1 , jx
        zphi = dble(dlat(i,j))*degrad
        zrla = dble(dlon(i,j))*degrad
        if ( dlat(i,j) > 89.999999 ) zrla = d_zero
        zrlap = pollam*degrad - zrla
        zarg1 = polcphi*sin(zrlap)
        zarg2 = polsphi*cos(zphi) - polcphi*sin(zphi)*cos(zrlap)
        znorm = d_one/sqrt(zarg1**2+zarg2**2)
        sindel = zarg1*znorm
        cosdel = zarg2*znorm
        do l = 1 , ll
          us = dble(u(i,j,l))*cosdel + dble(v(i,j,l))*sindel
          vs = -dble(u(i,j,l))*sindel + dble(v(i,j,l))*cosdel
          u(i,j,l) = real(us)
          v(i,j,l) = real(vs)
        end do
      end do
    end do
  else
    do j = 1 , iy
      do i = 1 , jx
        if ( (clon >= d_zero .and. dlon(i,j) >= 0.) .or.  &
             (clon < d_zero .and. dlon(i,j) < 0.) ) then
          x = (clon-dble(dlon(i,j)))*degrad*gridfc
        else if ( clon >= d_zero ) then
          if ( abs(clon-(dble(dlon(i,j))+360.0D0)) <  &
               abs(clon-dble(dlon(i,j))) ) then
            x = (clon-(dble(dlon(i,j))+360.0D0))*degrad*gridfc
          else
            x = (clon-dble(dlon(i,j)))*degrad*gridfc
          end if
        else if ( abs(clon-(dble(dlon(i,j))-360.0D0)) < &
                  abs(clon-dble(dlon(i,j))) ) then
          x = (clon-(dble(dlon(i,j))-360.0D0))*degrad*gridfc
        else
          x = (clon-dble(dlon(i,j)))*degrad*gridfc
        end if
        xs = sin(x)
        xc = cos(x)
        if ( clat >= d_zero ) then
          do l = 1 , ll
            d = dble(u(i,j,l))*xc - dble(v(i,j,l))*xs
            v(i,j,l) = real(dble(u(i,j,l))*xs + dble(v(i,j,l))*xc)
            u(i,j,l) = real(d)
          end do
        else
          do l = 1 , ll
            d = dble(u(i,j,l))*xc + dble(v(i,j,l))*xs
            v(i,j,l) = real(dble(v(i,j,l))*xc - dble(u(i,j,l))*xs)
            u(i,j,l) = real(d)
          end do
        end if
      end do
    end do
  end if
  end subroutine uvrot4nx
!
end module mod_uvrot
