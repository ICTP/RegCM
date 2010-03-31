!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine co2

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     model for plant carbon uptake and dead carbon decomposition
!
!            dtbat = surface output interval
!            cari = rate of carbon uptake by plants
!       pbp1d(n,i) = accumulated primary biomass
!      resp1d(n,i) = carbon in soil and in veg (kg c / m**2 / sec)
!              ra = leaf aerodynamic resistance factor
!             rap = leaf boundary layer resistance to co2
!           resps = soil respiration
!             rmp = physical mesophyllic resistance
!              rs = stomatal resistance
!             rsp = stomatal resistance to co2
!              rt = total mechanical resistance to co2
!
!     0.33 pptv co2 = 0.30 pptm dry matter in plants
!
!     decomposable soil carbon assumed incremented by npp-soil resp.
!
!     resistances for co2 are larger than those for h2o due to
!                    difference in molecular weight
!
      use mod_regcm_param
      use mod_param1 , only : dtbat
      use mod_bats , only : fdry , pbp1d , resp1d , ldoc1d , sigf ,     &
                   & rlai , gwmx0, solis , tlef1d , ssw1d , tg1d ,      &
                   & xlai , xlsai
      use mod_ictp01
      implicit none
!
! Local variables
!
      real(8) , dimension(nnsg,iym1) :: cari
      integer :: n , i
      real(8) :: rap , resps , rmp , rsp , rt
      real(8), external :: carbon
 
      rmp = 800.0
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
              rsp = rs(n,i)*1.7
              rap = ra(n,i)*1.5
              rt = rsp + rap + rmp
              cari(n,i) = sigf(n,i)*xlsai(n,i)*fdry(n,i)                &
                         & *carbon(solis(i)*rlai(n,i),tlef1d(n,i),rt,   &
                         & tg1d(n,i),xlai(n,i),xlsai(n,i))
              pbp1d(n,i) = pbp1d(n,i) + cari(n,i)*dtbat
            end if
          end if
        end do
      end do
 
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
              if ( pbp1d(n,i).lt.0 ) pbp1d(n,i) = 0.
              resps = 0.7E-7*resp1d(n,i)*dexp(0.1*(tg1d(n,i)-300.))     &
                    & *dmin1(1.D0,ssw1d(n,i)/(0.6*gwmx0(n,i)))
              resp1d(n,i) = resp1d(n,i) + (cari(n,i)-resps)*dtbat
              if ( resp1d(n,i).lt.0.0 ) resp1d(n,i) = 0.
            end if
          end if
        end do
      end do
      end subroutine co2
