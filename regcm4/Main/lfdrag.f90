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
      use mod_regcm_param
      use mod_bats
      use mod_ictp01
      use mod_constants , only : gti , wtur
      implicit none
!
! Local variables
!
      real(8) , dimension(nnsg,iym1) :: cdrmin , dlstaf , rib , rib1
      real(8) :: dthdz , ribi , sqrtf , tkb , u1 , u2 , zatild
      integer :: n , i
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
              tkb = wta0(n,i)*ts1d(n,i) + wtl0(n,i)*tlef1d(n,i)         &
                  & + wtg0(n,i)*tg1d(n,i)
              dlstaf(n,i) = ts1d(n,i) - sigf(n,i)                       &
                           & *tkb - (1.-sigf(n,i))*tg1d(n,i)
              if ( dlstaf(n,i).le.0 ) then
                dthdz = (1.-sigf(n,i))*tg1d(n,i) + sigf(n,i)            &
                      & *tkb - ts1d(n,i)
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
              zatild = (z1(n,i)-displa(lveg(n,i)))*sigf(n,i)            &
                     & + z1(n,i)*(1.-sigf(n,i))
              rib1(n,i) = gti*zatild/(ribd(n,i)*ts1d(n,i))
              rib(n,i) = rib1(n,i)*dlstaf(n,i)
              if ( rib(n,i).lt.0. ) then
                cdr(n,i) = cdrn(n,i)*(1.+24.5*dsqrt(-cdrn(n,i)*         &
                      & rib(n,i)))
                sqrtf = dmin1(dsqrt(-cdrn(n,i)/rib(n,i)),11.5D0/12.25D0)
                cdrd(n,i) = cdrn(n,i)*12.25*wtl0(n,i)*rib1(n,i)         &
                           & *sigf(n,i)*sqrtf
              else
                ribi = 1./(1.+11.5*rib(n,i))
                cdr(n,i) = cdrn(n,i)*ribi
                cdrd(n,i) = cdr(n,i)*ribi*11.5*rib1(n,i)*wtl0(n,i)      &
                           & *sigf(n,i)
                cdrmin(n,i) = dmax1(0.25*cdrn(n,i),6.D-4)
              end if
              if ( (rib(n,i).ge.0.) ) then
                if ( (cdr(n,i).lt.cdrmin(n,i)) ) then
                  cdr(n,i) = cdrmin(n,i)
                  cdrd(n,i) = 0.
                end if
              end if
            end if
          end if
        end do
      end do
 
      end subroutine lfdrag
