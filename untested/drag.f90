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
                     & cdr , aarea , z1 , clead , npts , lveg , tg1d ,  &
                     & sigf , displa , ts1d , ldoc1d ,  taf1d
      use mod_constants , only : gti , zoce , vonkar , wtur
      implicit none
!
! Local variables
!
      real(8) , dimension(nnsg,nbmax) :: cdrmin , rib , ribl , ribn
      real(8) :: dthdz , u1 , u2 , zatild
      integer :: n , np
!
!=======================================================================
!     1.   get neutral drag coefficient
!=======================================================================
      call dragdn
 
      do np = np1 , npts
        do n = 1 , nnsg
!=======================================================================
!         2.   compute stability as bulk rich. no. = rin/rid =
!         ri(numerator)/ri(denominator)
!=======================================================================
          if ( lveg(n,np).ne.0 ) then
            zatild = (z1(n,np)-displa(lveg(n,np)))*sigf(n,np) + z1(n,np)&
                   & *(1.-sigf(n,np))
          else
            zatild = z1(n,np)
          end if
          ribn(n,np) = zatild*gti*(ts1d(n,np)-sigf(n,np)*taf1d(n,np)-   &
                     & (1.-sigf(n,np))*tg1d(n,np))/ts1d(n,np)
!=======================================================================
!         2.1  compute the bulk richardson number;
!         first get avg winds to use for ri number by summing the
!         squares of horiz., vertical, and convective velocities
!=======================================================================
          if ( ribn(n,np).le.0. ) then
            dthdz = (1.-sigf(n,np))*tg1d(n,np) + sigf(n,np)*taf1d(n,np) &
                  & - ts1d(n,np)
            u1 = wtur + 2.*dsqrt(dthdz)
            ribd(n,np) = us1d(np)**2 + vs1d(np)**2 + u1**2
          else
            u2 = wtur
            ribd(n,np) = us1d(np)**2 + vs1d(np)**2 + u2**2
          end if
          vspda(n,np) = dsqrt(ribd(n,np))
          if ( vspda(n,np).lt.1. ) then
            vspda(n,np) = 1.
            ribd(n,np) = 1.
          end if
          rib(n,np) = ribn(n,np)/ribd(n,np)
!=======================================================================
!         3.   obtain drag coefficient as product of neutral value
!         and stability correction
!=======================================================================
!         ****   -0.4 < rib < 0.2   (deardorff, jgr, 1968, 2549-2557)
          if ( rib(n,np).lt.0. ) then
            cdr(n,np) = cdrn(n,np)                                      &
                      & *(1.0+24.5*dsqrt(-cdrn(n,np)*rib(n,np)))
          else
            cdr(n,np) = cdrn(n,np)/(1.0+11.5*rib(n,np))
          end if
!         3.1  apply lower limit to drag coefficient value
          cdrmin(n,np) = dmax1(0.25*cdrn(n,np),6.D-4)
          if ( cdr(n,np).lt.cdrmin(n,np) ) cdr(n,np) = cdrmin(n,np)
          cdrx(n,np) = cdr(n,np)
 
        end do
      end do
 
!=======================================================================
!     4.   obtain drag coefficient over sea ice as weighted average
!     over ice and leads
!     warning! the lat test below (4.1-4.3) is model dependent!
!=======================================================================
 
!     4.1  test if northern or southern hemisphere
      do np = np1 , npts
        do n = 1 , nnsg
                                                    ! check each point
!cc   if(lat(np).eq.    1) aarea(np) = 0.005  ! ccm specific code
!cc   if(lat(np).eq.    2) aarea(np) = 0.01
!cc   if(lat(np).ge.nlat2) aarea(np) = 0.04   !  4.2  antarctic
          if ( ldoc1d(n,np).gt.1.5 ) aarea(n,np) = 0.02
                                                    !  4.3  arctic
        end do
      end do
 
!     4.4  neutral cd over lead water
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.1.5 ) then       !  check each point
            cdrn(n,np) = (vonkar/dlog(z1(n,np)/zoce))**2
 
!           4.5  drag coefficient over leads
            ribl(n,np) = (1.-271.5/ts1d(n,np))*z1(n,np)*gti/ribd(n,np)
            if ( ribl(n,np).ge.0 ) then
              clead(n,np) = cdrn(n,np)/(1.+11.5*ribl(n,np))
            else
              clead(n,np) = cdrn(n,np)                                  &
                          & *(1.+24.5*dsqrt(-cdrn(n,np)*ribl(n,np)))
            end if
 
!           4.6  calculate weighted avg of ice and lead drag
!           coefficients
            cdrx(n,np) = (1.-aarea(n,np))*cdr(n,np) + aarea(n,np)       &
                       & *clead(n,np)
          end if
        end do
      end do
 
      end subroutine drag
