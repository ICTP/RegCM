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
 
      module mod_advection
!
! Horizontal and vertical advection.
!
      use mod_runparams
      use mod_main
      use mod_service
      private
 
      public :: hadv_x , hadv_d
      public :: vadv

      real(8) , parameter :: c287 = 0.287D+00
!
!     "relaxed" upstream scheme factors
!
      real(8) , parameter :: fact1 = 0.60D0
!hy   real(8) , parameter :: fact1 = 0.75D0
      real(8) , parameter :: fact2 = d_one - fact1
      real(8) , parameter :: falow = 1.0D-15
!
      contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!  HADV_X                                                             c
!                                                                     c
!     This subroutines computes the horizontal flux-divergence terms. c
!     second-order difference is used.                                c
!                                                                     c
!     ften   : is the tendency for variable 'f'.                      c
!                                                                     c
!     dxx    : is the horizontal distance.                            c
!              = dx4  for ind=1.                                      c
!              = dx   for ind=2.                                      c
!              = dx16 for ind=3.                                      c
!                                                                     c
!     j      : is the j'th slice of f anf ften.                       c
!                                                                     c
!     ind = 1 : for t and qv.                                         c
!         = 2 : for qc and qr.                                        c
!         = 3 : for u and v.                                          c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine hadv_x(ften,var,dxx,j,ind)
!
      implicit none
!
      real(8) ,intent (in) :: dxx
      integer ,intent (in) :: ind , j
      real(8) ,intent (inout), dimension(iy,kz) :: ften
#ifdef MPP1
      real(8) ,intent (in) , dimension(iy,kz,0:jxp+1) :: var
#else
      real(8) ,intent (in) , dimension(iy,kz,jx) :: var
#endif
!
      integer :: jm1 , jp1
      real(8) :: fx1 , fx2 , fy1 , fy2 , uavg1 , uavg2 ,&
               & vavg1 , vavg2
      integer :: i , k
      character (len=50) :: subroutine_name='hadv_x'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
      jm1 = j - 1
      jp1 = j + 1
!
#if defined(BAND) && (!defined(MPP1))
      if (jm1 == 0) jm1 = jx
      if (jp1 == jx+1) jp1 = 1
#endif
!
!----------------------------------------------------------------------
!
! ua, va : are p*u and p*v.
! msfx   : is the map scale factor at cross points.
!
!----------------------------------------------------------------------
!
      if ( ind == 1 ) then
!
!-----for t and qv:
!
        do k = 1 , kz
          do i = 2 , iym2
            ften(i,k) = ften(i,k) -                    &
                ((atm1%u(i+1,k,jp1)+atm1%u(i,k,jp1)) * &
                 (var(i,k,jp1)+var(i,k,j)) -           &
                 (atm1%u(i+1,k,j)+atm1%u(i,k,j)) *     &
                 (var(i,k,j)+var(i,k,jm1)) +           &
                 (atm1%v(i+1,k,jp1)+atm1%v(i+1,k,j)) * &
                 (var(i+1,k,j)+var(i,k,j)) -           &
                 (atm1%v(i,k,jp1)+atm1%v(i,k,j)) *     &
                 (var(i-1,k,j)+var(i,k,j))) /          &
                 (dxx*mddom%msfx(i,j)*mddom%msfx(i,j))
          end do
        end do
!
      else if ( ind == 2 ) then
!
!-----for qc and qr:
!       up-wind values of qc and qr are used.
!
        do k = 1 , kz
          do i = 2 , iym2
            uavg2 = d_half*(atm1%u(i+1,k,jp1)+atm1%u(i,k,jp1))
            uavg1 = d_half*(atm1%u(i+1,k,j)+atm1%u(i,k,j))
            if ( uavg2 >= d_zero ) then
              fx2 = fact1*var(i,k,j) + fact2*var(i,k,jp1)
            else
              fx2 = fact1*var(i,k,jp1) + fact2*var(i,k,j)
            end if
            if ( uavg1 >= d_zero ) then
              fx1 = fact1*var(i,k,jm1) + fact2*var(i,k,j)
            else
              fx1 = fact1*var(i,k,j) + fact2*var(i,k,jm1)
            end if
            vavg2 = d_half*(atm1%v(i+1,k,jp1)+atm1%v(i+1,k,j))
            vavg1 = d_half*(atm1%v(i,k,jp1)+atm1%v(i,k,j))
            if ( vavg2 >= d_zero ) then
              fy2 = fact1*var(i,k,j) + fact2*var(i+1,k,j)
            else
              fy2 = fact1*var(i+1,k,j) + fact2*var(i,k,j)
            end if
            if ( vavg1 >= d_zero ) then
              fy1 = fact1*var(i-1,k,j) + fact2*var(i,k,j)
            else
              fy1 = fact1*var(i,k,j) + fact2*var(i-1,k,j)
            end if
            ften(i,k) = ften(i,k)                                       &
                      & - (uavg2*fx2-uavg1*fx1+vavg2*fy2-vavg1*fy1)     &
                      & /(dxx*mddom%msfx(i,j)*mddom%msfx(i,j))
          end do
        end do
      else
        write(*,*) 'The advection scheme ',ind, &
                 & ' you required is not available.'
        stop
      end if
      call time_end(subroutine_name,idindx)
      end subroutine hadv_x
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine hadv_d(ften,var,dxx,j,ind)
!
      implicit none
!
      real(8) ,intent (in) :: dxx
      integer ,intent (in) :: ind , j
      real(8) ,intent (inout), dimension(iy,kz) :: ften
#ifdef MPP1
      real(8) ,intent (in) , dimension(iy,kz,0:jxp+1) :: var
#else
      real(8) ,intent (in) , dimension(iy,kz,jx) :: var
#endif
!
      integer :: jm1 , jp1
#ifndef BAND
      integer :: jdm1 , jdp1
#endif
      integer :: i , k , idx , idxm1 , idxp1
      real(8) :: ucmona , ucmonb , ucmonc , vcmona , vcmonb , vcmonc
!
      character (len=50) :: subroutine_name='hadv_d'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
      jm1 = j - 1
      jp1 = j + 1
!
!----------------------------------------------------------------------
!
! ua, va : are p*u and p*v.
! msfd   : is the map scale factor at dot points.
!
!----------------------------------------------------------------------
!
      if ( ind == 3 ) then
!
!-----for u and v:
!
#ifdef BAND
#if defined(BAND) && (!defined(MPP1))
        if (jm1 == 0) jm1 = jx
        if (jp1 == jx+1) jp1 = 1
#endif
!
        do k = 1 , kz
          do i = 2 , iym1
            idx = i
            idxp1 = i + 1
            idxp1 = min0(idxp1,iym1)
            idxm1 = i - 1
            idxm1 = max0(idxm1,2)
            ucmona = atm1%u(idxp1,k,j) +  &
                     d_two*atm1%u(idx,k,j) + atm1%u(idxm1,k,j)
            vcmona = atm1%v(idx,k,jp1) +  &
                     d_two*atm1%v(idx,k,j) + atm1%v(idx,k,jm1)
            ucmonb = atm1%u(idxp1,k,jp1) + d_two*atm1%u(idx,k,jp1)  &
                   & + atm1%u(idxm1,k,jp1) + ucmona
            vcmonb = atm1%v(idxp1,k,jp1) + d_two*atm1%v(idxp1,k,j)  &
                   & + atm1%v(idxp1,k,jm1) + vcmona
            ucmonc = atm1%u(idxp1,k,jm1) + d_two*atm1%u(idx,k,jm1)  &
                   & + atm1%u(idxm1,k,jm1) + ucmona
            vcmonc = atm1%v(idxm1,k,jp1) + d_two*atm1%v(idxm1,k,j)  &
                   & + atm1%v(idxm1,k,jm1) + vcmona
            ften(i,k) = ften(i,k)                            &
                      & - ((var(i,k,jp1)+var(i,k,j))*ucmonb- &
                      &    (var(i,k,j)+var(i,k,jm1))*ucmonc+ &
                      &    (var(i+1,k,j)+var(i,k,j))*vcmonb- &
                      &    (var(i,k,j)+var(i-1,k,j))*vcmonc) &
                      & /(dxx*mddom%msfd(i,j)*mddom%msfd(i,j))
          end do
        end do
!
!----------------------------------------------------------------------
#else
        jdp1 = j + 1
        jdm1 = j - 1
#ifdef MPP1
        if ( myid == 0 ) jdm1 = max0(jdm1,2)
        if ( myid == nproc-1 ) jdp1 = min0(jdp1,jendl-1)
#else
        jdp1 = min0(jdp1,jxm1)
        jdm1 = max0(jdm1,2)
#endif
!
        do k = 1 , kz
          do i = 2 , iym1
            idx = i
            idxp1 = i + 1
            idxp1 = min0(idxp1,iym1)
            idxm1 = i - 1
            idxm1 = max0(idxm1,2)
            ucmona = atm1%u(idxp1,k,j) + d_two*atm1%u(idx,k,j)        &
                   & + atm1%u(idxm1,k,j)
            vcmona = atm1%v(idx,k,jdp1) + d_two*atm1%v(idx,k,j)       &
                   & + atm1%v(idx,k,jdm1)
            ucmonb = atm1%u(idxp1,k,jdp1) + d_two*atm1%u(idx,k,jdp1)  &
                   & + atm1%u(idxm1,k,jdp1) + ucmona
            vcmonb = atm1%v(idxp1,k,jdp1) + d_two*atm1%v(idxp1,k,j)   &
                   & + atm1%v(idxp1,k,jdm1) + vcmona
            ucmonc = atm1%u(idxp1,k,jdm1) + d_two*atm1%u(idx,k,jdm1)  &
                   & + atm1%u(idxm1,k,jdm1) + ucmona
            vcmonc = atm1%v(idxm1,k,jdp1) + d_two*atm1%v(idxm1,k,j)   &
                   & + atm1%v(idxm1,k,jdm1) + vcmona
            ften(i,k) = ften(i,k) -                        &
                      & ((var(i,k,jp1)+var(i,k,j))*ucmonb- &
                      &  (var(i,k,j)+var(i,k,jm1))*ucmonc+ &
                      &  (var(i+1,k,j)+var(i,k,j))*vcmonb- &
                      &  (var(i,k,j)+var(i-1,k,j))*vcmonc) &
                      & /(dxx*mddom%msfd(i,j)*mddom%msfd(i,j))
          end do
        end do
#endif
!
      else
        write(*,*) 'The advection scheme ',ind, &
                 & ' you required is not available.'
        stop
      end if
      call time_end(subroutine_name,idindx)
      end subroutine hadv_d
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
! VADV                                                                c
!                                                                     c
!     This subroutine computes the vertical flux-divergence terms.    c
!                                                                     c
!     ften   : is the tendency of variable 'f'.                       c
!                                                                     c
!     fa     : is p*f.                                                c
!                                                                     c
!     j      : j'th slice of variable fa.                             c
!                                                                     c
!     ind = 1 : for t.                                                c
!           2 : for qv.                                               c
!           3 : for qc and qr.                                        c
!           4 : for u and v.                                          c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine vadv(ften,qdot,fa,j,ind)
!
      implicit none
!
      integer :: ind , j
      real(8) , dimension(iy,kz) :: fa , ften
#ifdef MPP1
      real(8) , dimension(iy,kzp1,0:jxp+1) , intent(in) :: qdot
#else
      real(8) , dimension(iy,kzp1,jx) , intent(in) :: qdot
#endif
      intent (in) fa , ind , j
      intent (inout) ften
!
      real(8) :: f1 , f2
      real(8) , dimension(iy,kz) :: fg
      integer :: i , k
!
      character (len=50) :: subroutine_name='vadv'
      integer :: idindx=0
!
#ifdef BAND
      integer :: jm1
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      call time_begin(subroutine_name,idindx)
#ifdef MPP1
      jm1 = j-1
#else
      jm1 = j-1
      if (jm1 == 0) jm1=jx
#endif
!
!     qdot   : is the vertical sigma-velocity
!     fg     : is the working space used to store the interlated
!              values.
!     psa    : is p* used to interpolate the temperature.
!
      if ( ind == 1 ) then
!
!-----vertical advection terms for:
!.....interpolate ta to full sigma levels:
!
        do i = 2 , iym2
          fg(i,1) = d_zero
        end do
        do k = 2 , kz
          do i = 2 , iym2
            fg(i,k) = twt(k,1)*fa(i,k) * &
                    & ((sps1%ps(i,j)*sigma(k)+r8pt)/ &
                    &  (sps1%ps(i,j)*a(k)+r8pt))**c287 + &
                    & twt(k,2)*fa(i,k-1) * &
                    & ((sps1%ps(i,j)*sigma(k)+r8pt)/ &
                    &  (sps1%ps(i,j)*a(k-1)+r8pt))**c287
          end do
        end do
!......k = 1
        do i = 2 , iym2
          ften(i,1) = ften(i,1) - qdot(i,2,j)*fg(i,2)/dsigma(1)
        end do
!......k = 2,kzm1
        do k = 2 , kzm1
          do i = 2 , iym2
            ften(i,k) = ften(i,k)                                     &
                    & - (qdot(i,k+1,j)*fg(i,k+1)-qdot(i,k,j)*fg(i,k)) &
                    & /dsigma(k)
          end do
        end do
!,.....k = kz
        do i = 2 , iym2
          ften(i,kz) = ften(i,kz) + qdot(i,kz,j)*fg(i,kz)/dsigma(kz)
        end do
!
      else if ( ind == 2 ) then
!
!-----vertical advection term for qv:
!.....interpolate qv to full sigma levels:
!
        do i = 2 , iym2
          fg(i,1) = d_zero
        end do
        do k = 2 , kz
          do i = 2 , iym2
! modif !!
            if ( fa(i,k) > falow .and. fa(i,k-1) > falow ) then
              fg(i,k) = fa(i,k)*(fa(i,k-1)/fa(i,k))**qcon(k)
            else
              fg(i,k) = d_zero
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
            ften(i,k) = ften(i,k)                                     &
                    & - (qdot(i,k+1,j)*fg(i,k+1)-qdot(i,k,j)*fg(i,k)) &
                    & /dsigma(k)
          end do
        end do
!,.....k = kz
        do i = 2 , iym2
          ften(i,kz) = ften(i,kz) + qdot(i,kz,j)*fg(i,kz)/dsigma(kz)
        end do
!
      else if ( ind == 3 ) then
!
!-----vertical advection terms for qc and qr:
!
!......k = 1
        do i = 2 , iym2
          if ( qdot(i,2,j) >= d_zero ) then
            f2 = fa(i,1)
          else
            f2 = fa(i,2)
          end if
          ften(i,1) = ften(i,1) - qdot(i,2,j)*f2/dsigma(1)
        end do
!......k = 2,kzm1
        do k = 2 , kzm1
          do i = 2 , iym2
            if ( qdot(i,k+1,j) >= d_zero ) then
              f2 = fa(i,k)
            else
              f2 = fa(i,k+1)
            end if
            if ( qdot(i,k,j) >= d_zero ) then
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
          if ( qdot(i,kz,j) >= d_zero ) then
            f1 = fa(i,kzm1)
          else
            f1 = fa(i,kz)
          end if
          ften(i,kz) = ften(i,kz) + qdot(i,kz,j)*f1/dsigma(kz)
        end do
!
      else if ( ind == 4 ) then
!
!-----vertical advection terms for u and v:
!.....interpolate ua or va to full sigma levels:
!
        do i = 2 , iym1
          fg(i,1) = d_zero
        end do
        do k = 2 , kz
          do i = 2 , iym1
            fg(i,k) = d_half*(fa(i,k)+fa(i,k-1))/mddom%msfd(i,j)
          end do
        end do
!......k = 1
        do i = 2 , iym1
          ften(i,1) = ften(i,1)                                         &
                    & - (qdot(i-1,2,jm1)+qdot(i,2,jm1)+qdot(i,2,j)      &
                    & +qdot(i-1,2,j))*fg(i,2)/(d_four*dsigma(1))
        end do
!......k = 2,kzm1
        do k = 2 , kzm1
          do i = 2 , iym1
            ften(i,k) = ften(i,k)                                       &
                      & - ((qdot(i,k+1,jm1)+qdot(i-1,k+1,jm1)+          &
                      &     qdot(i,k+1,j)+qdot(i-1,k+1,j))*fg(i,k+1)    &
                      & -(qdot(i,k,jm1)+qdot(i-1,k,jm1)+qdot(i,k,j)     &
                      & +qdot(i-1,k,j))*fg(i,k))/(d_four*dsigma(k))
          end do
        end do
!......k = kz
        do i = 2 , iym1
          ften(i,kz) = ften(i,kz)                                       &
                     & + (qdot(i,kz,jm1)+qdot(i-1,kz,jm1)+qdot(i,kz,j)  &
                     & +qdot(i-1,kz,j))*fg(i,kz)/(d_four*dsigma(kz))
        end do
!
 
 
      else if ( ind == 5 ) then
 
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
            ften(i,k) = ften(i,k)                                     &
                    & - (qdot(i,k+1,j)*fg(i,k+1)-qdot(i,k,j)*fg(i,k)) &
                    & /dsigma(k)
          end do
        end do
!,.....k = kz
        do i = 2 , iym2
          ften(i,kz) = ften(i,kz) + qdot(i,kz,j)*fg(i,kz)/dsigma(kz)
        end do
 
      end if
!
#else
!----------------------------------------------------------------------
      call time_begin(subroutine_name,idindx)
!
      if ( ind == 1 ) then
!
!-----vertical advection terms for:
!.....interpolate ta to full sigma levels:
!
        do i = 2 , iym2
          fg(i,1) = d_zero
        end do
        do k = 2 , kz
          do i = 2 , iym2
            fg(i,k) = twt(k,1)*fa(i,k) * &
                    & ((sps1%ps(i,j)*sigma(k)+r8pt)/ &
                    &  (sps1%ps(i,j)*a(k)+r8pt))**c287 + &
                    & twt(k,2)*fa(i,k-1) * &
                    & ((sps1%ps(i,j)*sigma(k)+r8pt)/ &
                    &  (sps1%ps(i,j)*a(k-1)+r8pt))**c287
          end do
        end do
!......k = 1
        do i = 2 , iym2
          ften(i,1) = ften(i,1) - qdot(i,2,j)*fg(i,2)/dsigma(1)
        end do
!......k = 2,kzm1
        do k = 2 , kzm1
          do i = 2 , iym2
            ften(i,k) = ften(i,k)                                     &
                    & - (qdot(i,k+1,j)*fg(i,k+1)-qdot(i,k,j)*fg(i,k)) &
                    & /dsigma(k)
          end do
        end do
!,.....k = kz
        do i = 2 , iym2
          ften(i,kz) = ften(i,kz) + qdot(i,kz,j)*fg(i,kz)/dsigma(kz)
        end do
!
      else if ( ind == 2 ) then
!
!-----vertical advection term for qv:
!.....interpolate qv to full sigma levels:
!
        do i = 2 , iym2
          fg(i,1) = d_zero
        end do
        do k = 2 , kz
          do i = 2 , iym2
! modif !!
            if ( fa(i,k) > falow .and. fa(i,k-1) > falow ) then
              fg(i,k) = fa(i,k)*(fa(i,k-1)/fa(i,k))**qcon(k)
            else
              fg(i,k) = d_zero
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
            ften(i,k) = ften(i,k)                                     &
                    & - (qdot(i,k+1,j)*fg(i,k+1)-qdot(i,k,j)*fg(i,k)) &
                    & /dsigma(k)
          end do
        end do
!,.....k = kz
        do i = 2 , iym2
          ften(i,kz) = ften(i,kz) + qdot(i,kz,j)*fg(i,kz)/dsigma(kz)
        end do
!
      else if ( ind == 3 ) then
!
!-----vertical advection terms for qc and qr:
!
!......k = 1
        do i = 2 , iym2
          if ( qdot(i,2,j) >= d_zero ) then
            f2 = fa(i,1)
          else
            f2 = fa(i,2)
          end if
          ften(i,1) = ften(i,1) - qdot(i,2,j)*f2/dsigma(1)
        end do
!......k = 2,kzm1
        do k = 2 , kzm1
          do i = 2 , iym2
            if ( qdot(i,k+1,j) >= d_zero ) then
              f2 = fa(i,k)
            else
              f2 = fa(i,k+1)
            end if
            if ( qdot(i,k,j) >= d_zero ) then
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
          if ( qdot(i,kz,j) >= d_zero ) then
            f1 = fa(i,kzm1)
          else
            f1 = fa(i,kz)
          end if
          ften(i,kz) = ften(i,kz) + qdot(i,kz,j)*f1/dsigma(kz)
        end do
!
      else if ( ind == 4 ) then
!
!-----vertical advection terms for u and v:
!.....interpolate ua or va to full sigma levels:
!
        do i = 2 , iym1
          fg(i,1) = d_zero
        end do
        do k = 2 , kz
          do i = 2 , iym1
            fg(i,k) = d_half*(fa(i,k)+fa(i,k-1))/mddom%msfd(i,j)
          end do
        end do
!......k = 1
        do i = 2 , iym1
          ften(i,1) = ften(i,1)                                         &
                    & - (qdot(i-1,2,j-1)+qdot(i,2,j-1)+qdot(i,2,j)      &
                    & +qdot(i-1,2,j))*fg(i,2)/(d_four*dsigma(1))
        end do
!......k = 2,kzm1
        do k = 2 , kzm1
          do i = 2 , iym1
            ften(i,k) = ften(i,k)                                       &
                      & - ((qdot(i,k+1,j-1)+qdot(i-1,k+1,j-1)+          &
                            qdot(i,k+1,j)  +qdot(i-1,k+1,j))*fg(i,k+1)  &
                      & -(qdot(i,k,j-1)+qdot(i-1,k,j-1)+qdot(i,k,j)     &
                      & +qdot(i-1,k,j))*fg(i,k))/(d_four*dsigma(k))
          end do
        end do
!......k = kz
        do i = 2 , iym1
          ften(i,kz) = ften(i,kz)                                       &
                     & + (qdot(i,kz,j-1)+qdot(i-1,kz,j-1)+qdot(i,kz,j)  &
                     & +qdot(i-1,kz,j))*fg(i,kz)/(d_four*dsigma(kz))
        end do
!
 
 
      else if ( ind == 5 ) then
 
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
            ften(i,k) = ften(i,k)                                     &
                    & - (qdot(i,k+1,j)*fg(i,k+1)-qdot(i,k,j)*fg(i,k)) &
                    & /dsigma(k)
          end do
        end do
!,.....k = kz
        do i = 2 , iym2
          ften(i,kz) = ften(i,kz) + qdot(i,kz,j)*fg(i,kz)/dsigma(kz)
        end do
 
      end if
!
#endif
      call time_end(subroutine_name,idindx)
      end subroutine vadv
!
      end module mod_advection
