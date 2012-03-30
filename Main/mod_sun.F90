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
 
module mod_sun
!
! Sun zenith and declination
!
  use mod_constants
  use mod_dynparam
  use mod_runparams
  use mod_atm_interface
  use mod_mpmessage
  use mod_service
!
  private
!
  public :: solar1 , zenitm
!
  contains
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! This subroutine computes the solar declination angle
! from the julian date.
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine solar1

    implicit none
!
    real(dp) :: decdeg
    real(dp) :: theta
!
!----------------------------------------------------------------------
!
    character (len=64) :: subroutine_name='solar1'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
    calday = yeardayfrac(idatex)
    theta = twopi*calday/dayspy
!
!   Solar declination in radians:
!
    declin = 0.006918D0 - 0.399912D0*dcos(theta) + &
             0.070257D0*dsin(theta) -              &
             0.006758D0*dcos(2.0D0*theta) +        &
             0.000907D0*dsin(2.0D0*theta) -        &
             0.002697D0*dcos(3.0D0*theta) +        &
             0.001480D0*dsin(3.0D0*theta)
    decdeg = declin/degrad
!
    write (aline, 99001) calday, decdeg
    call say
99001 format (11x,'*** Day ',f12.4,' solar declination angle = ',f12.8,&
        &   ' degrees.')
!
    call time_end(subroutine_name,idindx)
!
  end subroutine solar1
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! This subroutine calculates the cosine of the solar zenith angle
! for all longitude points of the RegCM domain. It needs as inputs
! the longitude and latitude of the points, the initial date of the
! simulation and the gmt. All these quantities are specified
! in the initialization procedure of RegCM
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine zenitm(coszrs,jstart,jend,istart,iend)
!
    implicit none
!
    integer, intent (in) :: jstart , jend , istart , iend
    real(dp) , pointer , intent (out), dimension(:,:) :: coszrs
!
    integer :: i , j
    real(dp) :: omga , tlocap , xt24
    character (len=64) :: subroutine_name='zenitm'
    real(dp) :: xxlat
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
!***********************************************************************
!
    xt24 = dble(idatex%second_of_day)/secph
    do i = istart , iend
      do j = jstart , jend
        tlocap = xt24 + mddom%xlon(j,i)/15.0D0
        tlocap = dmod(tlocap+houpd,houpd)
        omga = 15.0D0*(tlocap-12.0D0)*degrad
        xxlat = mddom%xlat(j,i)*degrad
!       coszrs = cosine of solar zenith angle
        coszrs(j,i) = dsin(declin)*dsin(xxlat) +           &
                      dcos(declin)*dcos(xxlat)*dcos(omga)
        coszrs(j,i) = dmax1(0.0D0,coszrs(j,i))
        coszrs(j,i) = dmin1(1.0D0,coszrs(j,i))
      end do
    end do
    call time_end(subroutine_name,idindx)
  end subroutine zenitm
!
end module mod_sun
