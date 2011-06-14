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
  use mod_main
  use mod_message
  use mod_service
!
#ifdef CLM
  use mod_clm
  use clm_time_manager , only : get_curr_calday
  use clm_varsur ,  only : numdays
  use shr_orb_mod , only : shr_orb_cosz , shr_orb_decl , &
                           shr_orb_params
#endif
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
! xtime  : forecast time in minutes.
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine solar1(xtime)

    implicit none
!
    real(8) , intent(in) :: xtime
!
    real(8) :: calday , decdeg
#ifdef CLM
    real(8) :: mvelp , obliq
    integer :: iyear_ad
    logical :: log_print
#else
    real(8) :: theta
#endif
!
!----------------------------------------------------------------------
!
    character (len=50) :: subroutine_name='solar1'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
#ifdef CLM
    log_print = .false.

    iyear_ad = idate1/1000000
    numdays = dayspy

!   Get eccen,obliq,mvelp,obliqr,lambm0,mvelpp
    call shr_orb_params(iyear_ad,r2ceccen,obliq,mvelp,r2cobliqr,      &
                      & r2clambm0,r2cmvelpp,log_print)
!
    calday = dble(julday) + dble(nnnnnn)*xdfbdy + &
                           (xtime/minph+gmt)/houpd

!   Get declin,eccf
    call shr_orb_decl(calday,r2ceccen,r2cmvelpp,r2clambm0,r2cobliqr,  &
                    & declin,r2ceccf)

!   convert declin to degrees
    declin = declin
    decdeg = declin/degrad
#else
    calday = dble(julday) + dble(nnnnnn)*xdfbdy + &
                     (xtime/minph+gmt)/houpd
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
!
    decdeg = declin/degrad
!
#endif
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
  subroutine zenitm(coszrs,ivmx,j)
!
    implicit none
!
    integer, intent (in) :: ivmx , j
    real(8) , intent (inout), dimension(iy) :: coszrs
!
    integer :: i
#ifdef CLM
    real(8) :: cldy , declinp1 , xxlon
#else
    real(8) :: omga , tlocap , xt24
#endif
    character (len=50) :: subroutine_name='zenitm'
    real(8) :: xxlat
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
!***********************************************************************
!
#ifdef CLM
    cldy = get_curr_calday()
!   cldy = dble(julday) + (nnnnnn)/4. + (xtime/60.+gmt)/24.
    call shr_orb_decl(cldy,r2ceccen,r2cmvelpp,r2clambm0,              &
         &            r2cobliqr,declinp1,r2ceccf)
    do i = 1 , ivmx
      xxlat = mddom%xlat(i,j)*degrad
      xxlon = mddom%xlon(i,j)*degrad
      coszrs(i) = shr_orb_cosz(cldy,xxlat,xxlon,declinp1)
      coszrs(i) = dmax1(0.0D0,coszrs(i))
      coszrs(i) = dmin1(1.0D0,coszrs(i))
    end do
#else
    xt24 = dmod(idatex%hour*minph+xtime,minpd)
    do i = 1 , ivmx
      tlocap = xt24/minph + mddom%xlon(i,j)/15.0D0
      tlocap = dmod(tlocap+houpd,houpd)
      omga = 15.0D0*(tlocap-12.0D0)*degrad
      xxlat = mddom%xlat(i,j)*degrad
!     coszrs = cosine of solar zenith angle
      coszrs(i) = dsin(declin)*dsin(xxlat) +           &
                    dcos(declin)*dcos(xxlat)*dcos(omga)
      coszrs(i) = dmax1(0.0D0,coszrs(i))
      coszrs(i) = dmin1(1.0D0,coszrs(i))
    end do
#endif
!
    call time_end(subroutine_name,idindx)
  end subroutine zenitm
!
end module mod_sun
