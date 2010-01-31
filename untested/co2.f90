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
!            c(4) = number of seconds in half-hour
!            cari = rate of carbon uptake by plants
!       pbp1d(n,np) = accumulated primary biomass
!      resp1d(n,np) = carbon in soil and in veg (kg c / m**2 / sec)
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
      use regcm_param
      use bats
      use ictp01
      implicit none
!
! Local variables
!
      real(8) , dimension(nnsg,nbmax) :: cari
      integer :: n , np
      real(8) :: rap , resps , rmp , rsp , rt
      real(8), external :: carbon
 
      rmp = 800.0
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
            if ( sigf(n,np).gt.0.001 ) then
              rsp = rs(n,np)*1.7
              rap = ra(n,np)*1.5
              rt = rsp + rap + rmp
              cari(n,np) = sigf(n,np)*xlsai(n,np)*fdry(n,np)            &
                         & *carbon(solis(np)*rlai(n,np),tlef1d(n,np),rt,&
                         & tg1d(n,np),xlai(n,np),xlsai(n,np))
              pbp1d(n,np) = pbp1d(n,np) + cari(n,np)*c(4)
            end if
          end if
        end do
      end do
 
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
            if ( sigf(n,np).gt.0.001 ) then
              if ( pbp1d(n,np).lt.0 ) pbp1d(n,np) = 0.
              resps = 0.7E-7*resp1d(n,np)*dexp(0.1*(tg1d(n,np)-300.))   &
                    & *dmin1(1.D0,ssw1d(n,np)/(0.6*gwmx0(n,np)))
              resp1d(n,np) = resp1d(n,np) + (cari(n,np)-resps)*c(4)
              if ( resp1d(n,np).lt.0.0 ) resp1d(n,np) = 0.
            end if
          end if
        end do
      end do
      end subroutine co2
