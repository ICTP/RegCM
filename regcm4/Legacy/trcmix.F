C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine trcmix(pmid, clat, coslat, n2o, ch4, cfc11, cfc12)
c-------------------------------------------------------------
c Specify zonal mean mass mixing ratios of CH4, N2O, CFC11 and
c CFC12
c          Code: J.T.Kiehl November 21, 1994
C------------------------------Parameters-------------------------------
c
c $Id: trcmix.F,v 1.2 1995/02/17 21:28:58 jhack Exp $
c $Author: jhack $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
c
#     include "parameter.inc"
#     include "include/parrad.cb"
#     include "include/tracer.cb"
C
c------------------------------input------------------------------------
      real(kind=8)  pmid(plond,plev),       ! model pressures
     $     clat(plond),                   ! current latitude in radians
     $     coslat(plond)                  ! cosine of latitude
c------------------------------output-----------------------------------
      real(kind=8)  n2o(plond,plev),        ! nitrous oxide mass mixing ratio
     $     ch4(plond,plev),        ! methane mass mixing ratio
     $     cfc11(plond,plev),      ! cfc11 mass mixing ratio
     $     cfc12(plond,plev)       ! cfc12 mass mixing ratio
c------------------------------local------------------------------------
      integer i,                   ! longitude loop index
     $        k                    ! level index
      real(kind=8)  dlat                    ! latitude in degrees
      real(kind=8)  xn2o,                   ! pressure scale height for n2o
     $     xch4,                   ! pressure scale height for ch4
     $     xcfc11,                 ! pressure scale height for cfc11
     $     xcfc12                  ! pressure scale height for cfc12
c     real(kind=8)  ch40,                   ! tropospheric mass mixing ratio for ch4
c    $     n2o0,                   ! tropospheric mass mixing ratio for n2o
c    $     cfc110,                 ! tropospheric mass mixing ratio for cfc11
c    $     cfc120                  ! tropospheric mass mixing ratio for cfc12
      real(kind=8)  ptrop(plond),                  ! pressure level of tropopause
     $     pratio                  ! pressure divided by ptrop
c
      do i = 1,plon
c
C set tropospheric mass mixing ratios
C
c     ch40   = rmwch4 * ch4vmr
c     n2o0   = rmwn2o * n2ovmr
c     cfc110 = rmwf11 * f11vmr
c     cfc120 = rmwf12 * f12vmr

c set stratospheric scale height factor for gases
      dlat = dabs(57.2958 * clat(i))
      if(dlat.le.45.0) then
        xn2o = 0.3478 + 0.00116 * dlat
        xch4 = 0.2353
        xcfc11 = 0.7273 + 0.00606 * dlat
        xcfc12 = 0.4000 + 0.00222 * dlat
      else
        xn2o = 0.4000 + 0.013333 * (dlat - 45)
        xch4 = 0.2353 + 0.0225489 * (dlat - 45)
        xcfc11 = 1.00 + 0.013333 * (dlat - 45)
        xcfc12 = 0.50 + 0.024444 * (dlat - 45)
      end if
c
c pressure of tropopause
      ptrop(i) = 250.0e2 - 150.0e2*coslat(i)**2.0
c
      end do
c
      do k = 1,plev
         do i = 1,plon
            if(pmid(i,k).ge.ptrop(i)) then
              ch4(i,k) = ch40
              n2o(i,k) = n2o0
              cfc11(i,k) = cfc110
              cfc12(i,k) = cfc120
            else
              pratio = pmid(i,k)/ptrop(i)
              ch4(i,k) = ch40 * (pratio)**xch4
              n2o(i,k) = n2o0 * (pratio)**xn2o
              cfc11(i,k) = cfc110 * (pratio)**xcfc11
              cfc12(i,k) = cfc120 * (pratio)**xcfc12
            end if
cfil
c   test 0. trace mixing ratio
c              ch4(i,k) = 1.e-30
c              n2o(i,k) = 1.e-30
c              cfc11(i,k) = 1.e-30
c              cfc12(i,k) = 1.e-30
cfil
         end do
      end do
      return
      end subroutine trcmix
