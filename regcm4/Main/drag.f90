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
 
      subroutine drag
 
!     *** determines surface transfer coeffs. at anemom. level from
!     *** lowest model level based on monin-obukov theory using
!     *** deardorff parameterization in terms of bulk richardson no.
 
!     ****  a.  calculates neutral drag coefficient (cdrn) as a fn of
!     ****             underlying surface
 
!     ****  b.  modifies cdrn as fn of bulk rich. no. of surface layer
 
      use mod_regcm_param
      use mod_bats , only : ribd , us1d , vs1d , vspda , cdrn , cdrx ,  &
                     & cdr , aarea , z1 , clead , ts1d , lveg , tg1d ,  &
                     & sigf , displa , ldoc1d ,  taf1d
      use mod_constants , only : gti , zoce , vonkar , wtur
      implicit none
!
! Local variables
!
      real(8) , dimension(nnsg,ixm1) :: cdrmin , rib , ribl , ribn
      real(8) :: dthdz , u1 , u2 , zatild
      integer :: n , i
!
!=======================================================================
!     1.   get neutral drag coefficient
!=======================================================================
      call dragdn
 
      do i = 2 , ixm1
        do n = 1 , nnsg
!=======================================================================
!         2.   compute stability as bulk rich. no. = rin/rid =
!         ri(numerator)/ri(denominator)
!=======================================================================
          if ( lveg(n,i).ne.0 ) then
            zatild = (z1(n,i)-displa(lveg(n,i)))*sigf(n,i) + z1(n,i)    &
                   & *(1.-sigf(n,i))
          else
            zatild = z1(n,i)
          end if
          ribn(n,i) = zatild*gti*(ts1d(n,i)-sigf(n,i)*taf1d(n,i)-       &
                     & (1.-sigf(n,i))*tg1d(n,i))/ts1d(n,i)
!=======================================================================
!         2.1  compute the bulk richardson number;
!         first get avg winds to use for ri number by summing the
!         squares of horiz., vertical, and convective velocities
!=======================================================================
          if ( ribn(n,i).le.0. ) then
            dthdz = (1.-sigf(n,i))*tg1d(n,i) + sigf(n,i)*taf1d(n,i)     &
                  & - ts1d(n,i)
            u1 = wtur + 2.*dsqrt(dthdz)
            ribd(n,i) = us1d(i)**2 + vs1d(i)**2 + u1**2
          else
            u2 = wtur
            ribd(n,i) = us1d(i)**2 + vs1d(i)**2 + u2**2
          end if
          vspda(n,i) = dsqrt(ribd(n,i))
          if ( vspda(n,i).lt.1. ) then
            vspda(n,i) = 1.
            ribd(n,i) = 1.
          end if
          rib(n,i) = ribn(n,i)/ribd(n,i)
!=======================================================================
!         3.   obtain drag coefficient as product of neutral value
!         and stability correction
!=======================================================================
!         ****   -0.4 < rib < 0.2   (deardorff, jgr, 1968, 2549-2557)
          if ( rib(n,i).lt.0. ) then
            cdr(n,i) = cdrn(n,i)*(1.0+24.5*dsqrt(-cdrn(n,i)*rib(n,i)))
          else
            cdr(n,i) = cdrn(n,i)/(1.0+11.5*rib(n,i))
          end if
!         3.1  apply lower limit to drag coefficient value
          cdrmin(n,i) = dmax1(0.25*cdrn(n,i),6.D-4)
          if ( cdr(n,i).lt.cdrmin(n,i) ) cdr(n,i) = cdrmin(n,i)
          cdrx(n,i) = cdr(n,i)
 
        end do
      end do
 
!=======================================================================
!     4.   obtain drag coefficient over sea ice as weighted average
!     over ice and leads
!     warning! the lat test below (4.1-4.3) is model dependent!
!=======================================================================
 
!     4.1  test if northern or southern hemisphere
      do i = 2 , ixm1
        do n = 1 , nnsg
                                                    ! check each point
!cc   if(lat(i).eq.    1) aarea(i) = 0.005  ! ccm specific code
!cc   if(lat(i).eq.    2) aarea(i) = 0.01
!cc   if(lat(i).ge.nlat2) aarea(i) = 0.04   !  4.2  antarctic
          if ( ldoc1d(n,i).gt.1.5 ) aarea(n,i) = 0.02
                                                    !  4.3  arctic
        end do
      end do
 
!     4.4  neutral cd over lead water
      do i = 2 , ixm1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.1.5 ) then       !  check each point
            cdrn(n,i) = (vonkar/dlog(z1(n,i)/zoce))**2
 
!           4.5  drag coefficient over leads
            ribl(n,i) = (1.-271.5/ts1d(n,i))*z1(n,i)*gti/ribd(n,i)
            if ( ribl(n,i).ge.0 ) then
              clead(n,i) = cdrn(n,i)/(1.+11.5*ribl(n,i))
            else
              clead(n,i) = cdrn(n,i)*(1.+24.5*dsqrt(-cdrn(n,i)*         &
                    & ribl(n,i)))
            end if
 
!           4.6  calculate weighted avg of ice and lead drag
!           coefficients
            cdrx(n,i) = (1.-aarea(n,i))*cdr(n,i) + aarea(n,i)*clead(n,i)
          end if
        end do
      end do
 
      end subroutine drag
