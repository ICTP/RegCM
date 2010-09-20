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
      use mod_main
      use mod_message
      use mod_date
      use service_mod
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the solar declination angle from the   c
!     julian date.                                                    c
!                                                                     c
!     xtime  : forecast time in minutes.                              c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine solar1(xtime)

      implicit none
!
      real(8) , intent(in) :: xtime
!
      real(8) :: calday , decdeg , delta
#ifdef CLM
      real(8) mvelp , obliq
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

!     Get eccen,obliq,mvelp,obliqr,lambm0,mvelpp
      call shr_orb_params(iyear_ad,r2ceccen,obliq,mvelp,r2cobliqr,      &
                        & r2clambm0,r2cmvelpp,log_print)
!
      calday = dble(julday) + (nnnnnn-nstrt0)/4. + (xtime/60.+gmt)/24.

!     Get delta,eccf
      call shr_orb_decl(calday,r2ceccen,r2cmvelpp,r2clambm0,r2cobliqr,  &
                      & delta,r2ceccf)

!     convert delta to degrees
      declin = delta
      decdeg = declin/degrad
#else
      calday = dble(julday) + (nnnnnn-nstrt0)/4. + (xtime/60.+gmt)/24.
      theta = 2.*mathpi*calday/dayspy
!
!     Solar declination in radians:
!
      delta = .006918 - .399912*dcos(theta) + .070257*dsin(theta)       &
            & - .006758*dcos(2.*theta) + .000907*dsin(2.*theta)         &
            & - .002697*dcos(3.*theta) + .001480*dsin(3.*theta)
!
      declin = delta
      decdeg = declin/degrad
!
#endif
      write (aline, 99001) calday, decdeg
      call say
99001 format (11x,'*** Day ',f12.8,' solar declination angle = ',f12.8,&
          &   ' degrees.')
!
      call time_end(subroutine_name,idindx)
      end subroutine solar1
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!   this subroutine calculates the cosine of the solar zenith angle
!   for all longitude points of the mm42d domain. it needs as inputs
!   the longitude and latitude of the points, the initial date of the
!   simulation and the gmt. all these quantities are specified
!   in the initialization procedure of RegCM
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      subroutine zenitm(coszrs,ivmx,jslc)
!
      implicit none
!
      integer, intent (in) :: ivmx , jslc
      real(8) , intent (inout), dimension(iy) :: coszrs
!
      integer :: ill
#ifdef CLM
      integer :: jj
      real(8) :: cldy , declinp1
#else
      real(8) :: omega , tlocap , xt24 , xxlat
#endif
      character (len=50) :: subroutine_name='zenitm'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
!***********************************************************************
!
#ifdef CLM
      cldy = get_curr_calday()
!      cldy = dble(julday) + (nnnnnn-nstrt0)/4. + (xtime/60.+gmt)/24.
      call shr_orb_decl(cldy,r2ceccen,r2cmvelpp,r2clambm0,              &
           &            r2cobliqr,declinp1,r2ceccf)
      jj = (jxp*myid) + jslc
      do ill = 1 , ivmx
        coszrs(ill) = shr_orb_cosz(cldy,r2cxlat_all(jj,ill),            &
            &                      r2cxlon_all(jj,ill),declinp1)
        coszrs(ill) = dmax1(0.d0,coszrs(ill))
      end do
#else
      xt24 = dmod(lhour*60.+xtime,1440.D0)
      do ill = 1 , ivmx
        tlocap = xt24/60. + xlong(ill,jslc)/15.
        tlocap = dmod(tlocap+24.,24.D0)
        omega = 15.*(tlocap-12.)*degrad
        xxlat = xlat(ill,jslc)*degrad
!       coszrs = cosine of solar zenith angle
        coszrs(ill) = dsin(declin)*dsin(xxlat) + dcos(declin)           &
                    & *dcos(xxlat)*dcos(omega)
        coszrs(ill) = dmax1(0.D0,coszrs(ill))
      end do
#endif
!
      call time_end(subroutine_name,idindx)
      end subroutine zenitm
!
      end module mod_sun
