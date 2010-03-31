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
!
      subroutine zenitm(coszrs,ivmx,jslc)
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
#ifdef CLM
      use mod_regcm_param , only : myid , jxp
!     use clm_time_manager , only : get_curr_calday
      use mod_date , only : declin , julday , gmt , nnnnnn , nstrt0 ,   &
                            xtime
      use shr_orb_mod , only : shr_orb_cosz , shr_orb_decl
      use mod_clm
#else
      use mod_regcm_param
      use mod_main
      use mod_date , only : declin , lhour , xtime
      use mod_constants , only : degrad
#endif
      implicit none
!
! Arguments
!
      integer, intent (in) :: ivmx , jslc
      real(kind=8) , intent (inout), dimension(iy) :: coszrs
!
! Local variables
!
      integer :: ill
#ifdef CLM
      integer :: jj
      real(kind=8) :: cldy , declinp1
#else
      real(kind=8) :: omega , tlocap , xt24 , xxlat
#endif
!
!***********************************************************************
!
#ifdef CLM
!     cldy = get_curr_calday()
      cldy = dble(julday) + (nnnnnn-nstrt0)/4. + (xtime/60.+gmt)/24.
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
      end subroutine zenitm
