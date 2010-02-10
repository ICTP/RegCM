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
      use mod_bats
      implicit none
!
! Local variables
!
      real(8) :: age1 , age2 , age3 , arg , arg2 , dela , dela0 , dels ,&
               & sge , tage
      integer :: n , np
      real(8) , dimension(nnsg,nbmax) :: sold
!
      age3 = 0.3
 
!=======================================================================
!     1.   partition soil evaporation and precipitation
!     between water and snow
!=======================================================================
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
 
            evapw(n,np) = fevpg(n,np)
            evaps(n,np) = scvk(n,np)*evapw(n,np)
            if ( ldoc1d(n,np).gt.1.5 ) evaps(n,np) = fevpg(n,np)
            evapw(n,np) = (1.-scvk(n,np))*evapw(n,np)
!
!           ******                tm  is temperature of precipitation
            if ( tm(n,np).ge.c(67) ) then
              pw(n,np) = prcp1d(n,np)*(1.-sigf(n,np))
              ps(n,np) = 0.0
            else
!             ******                snowing
              pw(n,np) = 0.0
              ps(n,np) = prcp1d(n,np)*(1.-sigf(n,np))
            end if
          end if
        end do
      end do
!
!=======================================================================
!     2.   update snow cover
!=======================================================================
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
            sold(n,np) = scv1d(n,np)
            scv1d(n,np) = scv1d(n,np) + c(4)                            &
                        & *(ps(n,np)-evaps(n,np)-sm(n,np)) + sdrop(n,np)
            scv1d(n,np) = dmax1(scv1d(n,np),0.D0)
            sag1d(n,np) = dmax1(sag1d(n,np),0.D0)
 
!           ******           snow cover except for antarctica
!=======================================================================
!           3.   increment non-dimensional "age" of snow;
!           10 mm snow restores surface to that of new snow.
!=======================================================================
            if ( scv1d(n,np).gt.0. ) then
              arg = 5.E3*(1./c(67)-1./tg1d(n,np))
              age1 = dexp(arg)
              arg2 = dmin1(0.D0,10.*arg)
              age2 = dexp(arg2)
              tage = age1 + age2 + age3
              dela0 = 1.E-6*c(4)
              dela = dela0*tage
              dels = 0.1*dmax1(0.D0,scv1d(n,np)-sold(n,np))
              sge = (sag1d(n,np)+dela)*(1.0-dels)
              sag1d(n,np) = dmax1(0.D0,sge)
            end if
 
!           ******           antarctica
            if ( scv1d(n,np).gt.800. ) sag1d(n,np) = 0.
          end if
        end do
      end do
 
      end subroutine snow
