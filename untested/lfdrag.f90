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
 
      subroutine lfdrag

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     recalculate stability dependent drag coefficient for vegetation,
!     given the neutral drag coefficient.
!
      use regcm_param
      use mod_bats
      use ictp01
      implicit none
!
! Local variables
!
      real(8) , dimension(nnsg,nbmax) :: cdrmin , dlstaf , rib , rib1
      real(8) :: dthdz , ribi , sqrtf , tkb , u1 , u2 , zatild
      integer :: n , np
!
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
            if ( sigf(n,np).gt.0.001 ) then
              tkb = wta0(n,np)*ts1d(n,np) + wtl0(n,np)*tlef1d(n,np)     &
                  & + wtg0(n,np)*tg1d(n,np)
              dlstaf(n,np) = ts1d(n,np) - sigf(n,np)                    &
                           & *tkb - (1.-sigf(n,np))*tg1d(n,np)
              if ( dlstaf(n,np).le.0 ) then
                dthdz = (1.-sigf(n,np))*tg1d(n,np) + sigf(n,np)         &
                      & *tkb - ts1d(n,np)
                u1 = c(90) + 2.*dsqrt(dthdz)
                ribd(n,np) = us1d(np)**2 + vs1d(np)**2 + u1**2
              else
                u2 = c(90)
                ribd(n,np) = us1d(np)**2 + vs1d(np)**2 + u2**2
              end if
              vspda(n,np) = dsqrt(ribd(n,np))
              if ( vspda(n,np).lt.1. ) then
                vspda(n,np) = 1.
                ribd(n,np) = 1.
              end if
              zatild = (z1(n,np)-displa(lveg(n,np)))*sigf(n,np)         &
                     & + z1(n,np)*(1.-sigf(n,np))
              rib1(n,np) = c(54)*zatild/(ribd(n,np)*ts1d(n,np))
              rib(n,np) = rib1(n,np)*dlstaf(n,np)
              if ( rib(n,np).lt.0. ) then
                cdr(n,np) = cdrn(n,np)                                  &
                          & *(1.+24.5*dsqrt(-cdrn(n,np)*rib(n,np)))
                sqrtf = dmin1(dsqrt(-cdrn(n,np)/rib(n,np)),             &
                      & 11.5D0/12.25D0)
                cdrd(n,np) = cdrn(n,np)*12.25*wtl0(n,np)*rib1(n,np)     &
                           & *sigf(n,np)*sqrtf
              else
                ribi = 1./(1.+11.5*rib(n,np))
                cdr(n,np) = cdrn(n,np)*ribi
                cdrd(n,np) = cdr(n,np)*ribi*11.5*rib1(n,np)*wtl0(n,np)  &
                           & *sigf(n,np)
                cdrmin(n,np) = dmax1(0.25*cdrn(n,np),6.D-4)
              end if
              if ( (rib(n,np).ge.0.) ) then
                if ( (cdr(n,np).lt.cdrmin(n,np)) ) then
                  cdr(n,np) = cdrmin(n,np)
                  cdrd(n,np) = 0.
                end if
              end if
            end if
          end if
        end do
      end do
 
      end subroutine lfdrag
