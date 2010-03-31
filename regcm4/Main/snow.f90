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
 
      subroutine snow
!
!     update snow cover and snow age
!
!     three-part if block:
!       if snow cover < 0, then snow cover and snow age = 0
!       if antarctica, snow age = 0 (katabatic winds keep snow fresh)
!       if elsewhere, snow age follows given formulae
!
!        ps = snow precipitation rate
!     evaps = moisture flux from ground to atmosphere
!        sm = snow melt rate
!     sdrop = snow fallen from vegetation
!
!     aging of snow consists of three factors:
!           age1: snow crystal growth
!           age2: surface melting
!           age3: accumulation  of other particles, soot, etc., which
!                      is small in southern hemisphere
!
      use mod_regcm_param
      use mod_param1 , only : dtbat
      use mod_bats
      use mod_constants , only : tzero
      implicit none
!
! Local variables
!
      real(8) :: age1 , age2 , age3 , arg , arg2 , dela , dela0 , dels ,&
               & sge , tage
      integer :: n , i
      real(8) , dimension(nnsg,iym1) :: sold
!
      age3 = 0.3
 
!=======================================================================
!     1.   partition soil evaporation and precipitation
!     between water and snow
!=======================================================================
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
 
            evapw(n,i) = fevpg(n,i)
            evaps(n,i) = scvk(n,i)*evapw(n,i)
            if ( ldoc1d(n,i).gt.1.5 ) evaps(n,i) = fevpg(n,i)
            evapw(n,i) = (1.-scvk(n,i))*evapw(n,i)
!
!           ******                tm  is temperature of precipitation
            if ( tm(n,i).ge.tzero ) then
              pw(n,i) = prcp1d(n,i)*(1.-sigf(n,i))
              ps(n,i) = 0.0
            else
!             ******                snowing
              pw(n,i) = 0.0
              ps(n,i) = prcp1d(n,i)*(1.-sigf(n,i))
            end if
          end if
        end do
      end do
!
!=======================================================================
!     2.   update snow cover
!=======================================================================
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            sold(n,i) = scv1d(n,i)
            scv1d(n,i) = scv1d(n,i) + dtbat                             &
                        & *(ps(n,i)-evaps(n,i)-sm(n,i)) + sdrop(n,i)
            scv1d(n,i) = dmax1(scv1d(n,i),0.D0)
            sag1d(n,i) = dmax1(sag1d(n,i),0.D0)
 
!           ******           snow cover except for antarctica
!=======================================================================
!           3.   increment non-dimensional "age" of snow;
!           10 mm snow restores surface to that of new snow.
!=======================================================================
            if ( scv1d(n,i).gt.0. ) then
              arg = 5.E3*(1./tzero-1./tg1d(n,i))
              age1 = dexp(arg)
              arg2 = dmin1(0.D0,10.*arg)
              age2 = dexp(arg2)
              tage = age1 + age2 + age3
              dela0 = 1.E-6*dtbat
              dela = dela0*tage
              dels = 0.1*dmax1(0.D0,scv1d(n,i)-sold(n,i))
              sge = (sag1d(n,i)+dela)*(1.0-dels)
              sag1d(n,i) = dmax1(0.D0,sge)
            end if
 
!           ******           antarctica
            if ( scv1d(n,i).gt.800. ) sag1d(n,i) = 0.
          end if
        end do
      end do
 
      end subroutine snow
