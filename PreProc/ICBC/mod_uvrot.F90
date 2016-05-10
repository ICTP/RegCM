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
  !
  ! Change U and V from true (N,E) to map values (X,Y)
  !
  subroutine uvrot4(u,v,dlon,dlat,clon,clat,gridfc,jx,iy,ll,pollon, &
                    pollat,lgtype)
    implicit none
    real(rkx) , intent(in) :: clat , clon , gridfc , pollat , pollon
    integer(ik4) , intent(in) :: iy , jx , ll
    character(len=6) , intent(in) :: lgtype
    real(rkx) , dimension(jx,iy) , intent(in) :: dlat , dlon
    real(rkx) , dimension(jx,iy,ll) , intent(inout) :: u , v
    real(rkx) :: cosdel , d , polcphi , pollam , polphi , polsphi ,  &
            sindel , us , vs , x , xc , xs , zarg1 , zarg2 , znorm , &
            zphi , zrla , zrlap
    integer(ik4) :: i , j , l

    if ( lgtype == 'NORMER' ) then
      return
    end if
    if ( lgtype == 'ROTMER' ) then
      if ( pollat > d_zero ) then
        pollam = pollon + deg180
        polphi = deg90 - pollat
      else
        polphi = deg90 + pollat
        pollam = pollon
      end if
      if ( pollam > deg180 ) pollam = pollam - deg360

      polcphi = cos(degrad*polphi)
      polsphi = sin(degrad*polphi)

      do j = 1 , iy
        do i = 1 , jx
          zphi = dlat(i,j)*degrad
          zrla = dlon(i,j)*degrad
          if ( dlat(i,j) > 89.999999_rkx ) zrla = d_zero
          zrlap = pollam*degrad - zrla
          zarg1 = polcphi*sin(zrlap)
          zarg2 = polsphi*cos(zphi) - polcphi*sin(zphi)*cos(zrlap)
          znorm = d_one/sqrt(zarg1**2+zarg2**2)
          sindel = zarg1*znorm
          cosdel = zarg2*znorm
          do l = 1 , ll
            us = u(i,j,l)*cosdel - v(i,j,l)*sindel
            vs = u(i,j,l)*sindel + v(i,j,l)*cosdel
            u(i,j,l) = us
            v(i,j,l) = vs
          end do
        end do
      end do
    else
      do j = 1 , iy
        do i = 1 , jx
          if ( (clon >= d_zero .and. dlon(i,j) >= deg00) .or.   &
               (clon < d_zero  .and. dlon(i,j) < deg00) ) then
            x = (clon-dlon(i,j))*degrad*gridfc
          else if ( clon >= d_zero ) then
            if ( abs(clon-(dlon(i,j)+deg360)) < abs(clon-dlon(i,j)) ) then
              x = (clon-(dlon(i,j)+deg360))*degrad*gridfc
            else
              x = (clon-dlon(i,j))*degrad*gridfc
            end if
          else if ( abs(clon-(dlon(i,j)-deg360)) < abs(clon-dlon(i,j)) ) then
            x = (clon-(dlon(i,j)-deg360))*degrad*gridfc
          else
            x = (clon-dlon(i,j))*degrad*gridfc
          end if
          xs = sin(x)
          xc = cos(x)
          if ( clat >= d_zero ) then
            do l = 1 , ll
              d = v(i,j,l)*xs + u(i,j,l)*xc
              v(i,j,l) = v(i,j,l)*xc - u(i,j,l)*xs
              u(i,j,l) = d
            end do
          else
            do l = 1 , ll
              d = -v(i,j,l)*xs + u(i,j,l)*xc
              v(i,j,l) = v(i,j,l)*xc + u(i,j,l)*xs
              u(i,j,l) = d
            end do
          end if
        end do
      end do
    end if
  end subroutine uvrot4
  !
  ! Change U and V from map values (X,Y) to true (N,E)
  !
  subroutine uvrot4nx(u,v,dlon,dlat,clon,clat,gridfc,jx,iy,ll, &
                      pollon,pollat,lgtype)
    implicit none

    real(rkx) , intent(in) :: clat , clon , gridfc , pollat , pollon
    integer(ik4) , intent(in) :: iy , jx , ll
    character(len=6) , intent(in) :: lgtype
    real(rkx) , dimension(jx,iy) , intent(in) :: dlat , dlon
    real(rkx) , dimension(jx,iy,ll) , intent(inout) :: u , v
    real(rkx) :: cosdel , d , polcphi , pollam , polphi , polsphi ,  &
            sindel , us , vs , x , xc , xs , zarg1 , zarg2 , znorm , &
            zphi , zrla , zrlap
    integer(ik4) :: i , j , l

    if ( lgtype == 'ROTMER' ) then
      if ( pollat > d_zero ) then
        pollam = pollon + deg180
        polphi = deg90 - pollat
      else
        polphi = deg90 + pollat
        pollam = pollon
      end if
      if ( pollam > deg180 ) pollam = pollam - deg360

      polcphi = cos(degrad*polphi)
      polsphi = sin(degrad*polphi)

      do j = 1 , iy
        do i = 1 , jx
          zphi = dlat(i,j)*degrad
          zrla = dlon(i,j)*degrad
          if ( dlat(i,j) > 89.999999_rkx ) zrla = d_zero
          zrlap = pollam*degrad - zrla
          zarg1 = polcphi*sin(zrlap)
          zarg2 = polsphi*cos(zphi) - polcphi*sin(zphi)*cos(zrlap)
          znorm = d_one/sqrt(zarg1**2+zarg2**2)
          sindel = zarg1*znorm
          cosdel = zarg2*znorm
          do l = 1 , ll
            us = u(i,j,l)*cosdel + v(i,j,l)*sindel
            vs = -u(i,j,l)*sindel + v(i,j,l)*cosdel
            u(i,j,l) = us
            v(i,j,l) = vs
          end do
        end do
      end do
    else
      do j = 1 , iy
        do i = 1 , jx
          if ( (clon >= d_zero .and. dlon(i,j) >= deg00) .or.  &
               (clon < d_zero .and. dlon(i,j) < deg00) ) then
            x = (clon-dlon(i,j))*degrad*gridfc
          else if ( clon >= d_zero ) then
            if ( abs(clon-(dlon(i,j)+deg360)) < abs(clon-dlon(i,j)) ) then
              x = (clon-(dlon(i,j)+deg360))*degrad*gridfc
            else
              x = (clon-dlon(i,j))*degrad*gridfc
            end if
          else if ( abs(clon-(dlon(i,j)-deg360)) < abs(clon-dlon(i,j)) ) then
            x = (clon-(dlon(i,j)-deg360))*degrad*gridfc
          else
            x = (clon-dlon(i,j))*degrad*gridfc
          end if
          xs = sin(x)
          xc = cos(x)
          if ( clat >= d_zero ) then
            do l = 1 , ll
              d = u(i,j,l)*xc - v(i,j,l)*xs
              v(i,j,l) = u(i,j,l)*xs + v(i,j,l)*xc
              u(i,j,l) = d
            end do
          else
            do l = 1 , ll
              d = u(i,j,l)*xc + v(i,j,l)*xs
              v(i,j,l) = v(i,j,l)*xc - u(i,j,l)*xs
              u(i,j,l) = d
            end do
          end if
        end do
      end do
    end if
  end subroutine uvrot4nx
!
end module mod_uvrot
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
