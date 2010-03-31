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
 
      subroutine vadv(ften,fa,j,ind)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the vertical flux-divergence terms.    c
!                                                                     c
!     ften   : is the tendency of variable 'f'.                       c
!                                                                     c
!     fa     : is p*f.                                                c
!                                                                     c
!     qdot   : is the vertical sigma-velocity                         c
!                                                                     c
!     f      : is the working space used to store the interlated      c
!              values.                                                c
!                                                                     c
!     psa    : is p* used to interpolate the temperature.             c
!                                                                     c
!     j      : j'th slice of variable fa.                             c
!                                                                     c
!     ind = 1 : for t.                                                c
!           2 : for qv.                                               c
!           3 : for qc and qr.                                        c
!           4 : for u and v.                                          c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mod_regcm_param
      use mod_param3 , only : dsigma , ptop , qcon , twt , a , sigma
      use mod_main
      use mod_cvaria
      implicit none
!
! Dummy arguments
!
      integer :: ind , j
      real(8) , dimension(iy,kz) :: fa , ften
      intent (in) fa , ind , j
      intent (inout) ften
!
! Local variables
!
      real(8) :: f1 , f2
      real(8) , dimension(iy,kz) :: fg
      integer :: i , k
!
!----------------------------------------------------------------------
!
      if ( ind.eq.1 ) then
!
!-----vertical advection terms for:
!.....interpolate ta to full sigma levels:
!
        do i = 2 , iym2
          fg(i,1) = 0.
        end do
        do k = 2 , kz
          do i = 2 , iym2
            fg(i,k) = twt(k,1)*fa(i,k)                                  &
                    & *((psa(i,j)*sigma(k)+ptop)/(psa(i,j)*a(k)+ptop))  &
                    & **0.287 + twt(k,2)*fa(i,k-1)                      &
                    & *((psa(i,j)*sigma(k)+ptop)/(psa(i,j)*a(k-1)+ptop))&
                    & **0.287
          end do
        end do
!......k = 1
        do i = 2 , iym2
          ften(i,1) = ften(i,1) - qdot(i,2,j)*fg(i,2)/dsigma(1)
        end do
!......k = 2,kzm1
        do k = 2 , kzm1
          do i = 2 , iym2
            ften(i,k) = ften(i,k)                                       &
                      & - (qdot(i,k+1,j)*fg(i,k+1)-qdot(i,k,j)*fg(i,k)) &
                      & /dsigma(k)
          end do
        end do
!,.....k = kz
        do i = 2 , iym2
          ften(i,kz) = ften(i,kz) + qdot(i,kz,j)*fg(i,kz)/dsigma(kz)
        end do
!
      else if ( ind.eq.2 ) then
!
!-----vertical advection term for qv:
!.....interpolate qv to full sigma levels:
!
        do i = 2 , iym2
          fg(i,1) = 0.
        end do
        do k = 2 , kz
          do i = 2 , iym2
! modif !!
            if ( fa(i,k).gt.1.E-15 .and. fa(i,k-1).gt.1.E-15 ) then
              fg(i,k) = fa(i,k)*(fa(i,k-1)/fa(i,k))**qcon(k)
            else
              fg(i,k) = 0.
            end if
          end do
        end do
!......k = 1
        do i = 2 , iym2
          ften(i,1) = ften(i,1) - qdot(i,2,j)*fg(i,2)/dsigma(1)
        end do
!......k = 2,kzm1
        do k = 2 , kzm1
          do i = 2 , iym2
            ften(i,k) = ften(i,k)                                       &
                      & - (qdot(i,k+1,j)*fg(i,k+1)-qdot(i,k,j)*fg(i,k)) &
                      & /dsigma(k)
          end do
        end do
!,.....k = kz
        do i = 2 , iym2
          ften(i,kz) = ften(i,kz) + qdot(i,kz,j)*fg(i,kz)/dsigma(kz)
        end do
!
      else if ( ind.eq.3 ) then
!
!-----vertical advection terms for qc and qr:
!
!......k = 1
        do i = 2 , iym2
          if ( qdot(i,2,j).ge.0. ) then
            f2 = fa(i,1)
          else
            f2 = fa(i,2)
          end if
          ften(i,1) = ften(i,1) - qdot(i,2,j)*f2/dsigma(1)
        end do
!......k = 2,kzm1
        do k = 2 , kzm1
          do i = 2 , iym2
            if ( qdot(i,k+1,j).ge.0. ) then
              f2 = fa(i,k)
            else
              f2 = fa(i,k+1)
            end if
            if ( qdot(i,k,j).ge.0. ) then
              f1 = fa(i,k-1)
            else
              f1 = fa(i,k)
            end if
            ften(i,k) = ften(i,k) - (qdot(i,k+1,j)*f2-qdot(i,k,j)*f1)   &
                      & /dsigma(k)
          end do
        end do
!......k = kz
        do i = 2 , iym2
          if ( qdot(i,kz,j).ge.0. ) then
            f1 = fa(i,kzm1)
          else
            f1 = fa(i,kz)
          end if
          ften(i,kz) = ften(i,kz) + qdot(i,kz,j)*f1/dsigma(kz)
        end do
!
      else if ( ind.eq.4 ) then
!
!-----vertical advection terms for u and v:
!.....interpolate ua or va to full sigma levels:
!
        do i = 2 , iym1
          fg(i,1) = 0.
        end do
        do k = 2 , kz
          do i = 2 , iym1
            fg(i,k) = 0.5*(fa(i,k)+fa(i,k-1))/msfd(i,j)
          end do
        end do
!......k = 1
        do i = 2 , iym1
          ften(i,1) = ften(i,1)                                         &
                    & - (qdot(i-1,2,j-1)+qdot(i,2,j-1)+qdot(i,2,j)      &
                    & +qdot(i-1,2,j))*fg(i,2)/(4.*dsigma(1))
        end do
!......k = 2,kzm1
        do k = 2 , kzm1
          do i = 2 , iym1
            ften(i,k) = ften(i,k)                                       &
                      & - ((qdot(i,k+1,j-1)+qdot(i-1,k+1,j-1)+qdot(i,   &
                      & k+1,j)+qdot(i-1,k+1,j))*fg(i,k+1)               &
                      & -(qdot(i,k,j-1)+qdot(i-1,k,j-1)+qdot(i,k,j)     &
                      & +qdot(i-1,k,j))*fg(i,k))/(4.*dsigma(k))
          end do
        end do
!......k = kz
        do i = 2 , iym1
          ften(i,kz) = ften(i,kz)                                       &
                     & + (qdot(i,kz,j-1)+qdot(i-1,kz,j-1)+qdot(i,kz,j)  &
                     & +qdot(i-1,kz,j))*fg(i,kz)/(4.*dsigma(kz))
        end do
!
 
 
      else if ( ind.eq.5 ) then
 
        do k = 2 , kz
          do i = 2 , iym2
            fg(i,k) = twt(k,1)*fa(i,k) + twt(k,2)*fa(i,k-1)
          end do
        end do
 
!......k = 1
        do i = 2 , iym2
          ften(i,1) = ften(i,1) - qdot(i,2,j)*fg(i,2)/dsigma(1)
        end do
!......k = 2,kzm1
        do k = 2 , kzm1
          do i = 2 , iym2
            ften(i,k) = ften(i,k)                                       &
                      & - (qdot(i,k+1,j)*fg(i,k+1)-qdot(i,k,j)*fg(i,k)) &
                      & /dsigma(k)
          end do
        end do
!,.....k = kz
        do i = 2 , iym2
          ften(i,kz) = ften(i,kz) + qdot(i,kz,j)*fg(i,kz)/dsigma(kz)
        end do
 
      else
      end if
!
      end subroutine vadv
