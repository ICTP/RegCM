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
      use mod_dynparam
      use mod_runparams
      use mod_cvaria
      use mod_main

      private
 
      public :: hadv_t , hadv_u , hadv_v , hadvqv , hadvqc , hadvch
      public :: vadv

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
      subroutine hadv_t(ften,dxx,j,ind)
!
      implicit none
!
      real(8) ,intent (in) :: dxx
      integer ,intent (in) :: ind , j
      real(8) ,intent (inout), dimension(iy,kz) :: ften
!
      integer :: jm1 , jp1
      real(8) :: fact1 , fact2 , fx1 , fx2 , fy1 , fy2 , uavg1 , uavg2 ,&
               & vavg1 , vavg2
      integer :: i , k
!
      jm1 = j - 1
      jp1 = j + 1
!
#if defined(BAND) && (!defined(MPP1))
      if(jm1.eq.0) jm1 = jx
      if(jp1.eq.jx+1) jp1 = 1
#endif
!
!----------------------------------------------------------------------
!
! ua, va : are p*u and p*v.
! msfx   : is the map scale factor at cross points.
!
!----------------------------------------------------------------------
!
      if ( ind.eq.1 ) then
!
!-----for t and qv:
!
        do k = 1 , kz
          do i = 2 , iym2
            ften(i,k) = ften(i,k)                                       &
                      & - ((ua(i+1,k,jp1)+ua(i,k,jp1))*(t(i,k,jp1)      &
                      & +t(i,k,j))-(ua(i+1,k,j)+ua(i,k,j))              &
                      & *(t(i,k,j)+t(i,k,jm1))                          &
                      & +(va(i+1,k,jp1)+va(i+1,k,j))                    &
                      & *(t(i+1,k,j)+t(i,k,j))-(va(i,k,jp1)+va(i,k,j))  &
                      & *(t(i-1,k,j)+t(i,k,j)))                         &
                      & /(dxx*msfx(i,j)*msfx(i,j))
          end do
        end do
!
      else if ( ind.eq.2 ) then
!
!       implement a "relaxed" upstream scheme
!
!hy     fact1=0.75
        fact1 = 0.60
        fact2 = 1. - fact1
!
!-----for qc and qr:
!       up-wind values of qc and qr are used.
!
        do k = 1 , kz
          do i = 2 , iym2
            uavg2 = 0.5*(ua(i+1,k,jp1)+ua(i,k,jp1))
            uavg1 = 0.5*(ua(i+1,k,j)+ua(i,k,j))
            if ( uavg2.ge.0. ) then
              fx2 = fact1*t(i,k,j) + fact2*t(i,k,jp1)
            else
              fx2 = fact1*t(i,k,jp1) + fact2*t(i,k,j)
            end if
            if ( uavg1.ge.0. ) then
              fx1 = fact1*t(i,k,jm1) + fact2*t(i,k,j)
            else
              fx1 = fact1*t(i,k,j) + fact2*t(i,k,jm1)
            end if
            vavg2 = 0.5*(va(i+1,k,jp1)+va(i+1,k,j))
            vavg1 = 0.5*(va(i,k,jp1)+va(i,k,j))
            if ( vavg2.ge.0. ) then
              fy2 = fact1*t(i,k,j) + fact2*t(i+1,k,j)
            else
              fy2 = fact1*t(i+1,k,j) + fact2*t(i,k,j)
            end if
            if ( vavg1.ge.0. ) then
              fy1 = fact1*t(i-1,k,j) + fact2*t(i,k,j)
            else
              fy1 = fact1*t(i,k,j) + fact2*t(i-1,k,j)
            end if
            ften(i,k) = ften(i,k)                                       &
                      & - (uavg2*fx2-uavg1*fx1+vavg2*fy2-vavg1*fy1)     &
                      & /(dxx*msfx(i,j)*msfx(i,j))
          end do
        end do
      else
        write(*,*) 'The T advection scheme ',ind, &
                 & ' you required is not available.'
        stop
      end if
      end subroutine hadv_t
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine hadv_u(ften,dxx,j,ind)
!
      implicit none
!
      real(8) ,intent (in) :: dxx
      integer ,intent (in) :: ind , j
      real(8) ,intent (inout), dimension(iy,kz) :: ften
!
      integer :: jm1 , jp1
#ifndef BAND
      integer :: jdm1 , jdp1
#endif
      integer :: i , k , idx , idxm1 , idxp1
      real(8) :: ucmona , ucmonb , ucmonc , vcmona , vcmonb , vcmonc
!
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
      if ( ind.eq.3 ) then
!
!-----for u and v:
!
#ifdef BAND
#if defined(BAND) && (!defined(MPP1))
        if(jm1.eq.0) jm1 = jx
        if(jp1.eq.jx+1) jp1 = 1
#endif
!
        do k = 1 , kz
          do i = 2 , iym1
            idx = i
            idxp1 = i + 1
            idxp1 = min0(idxp1,iym1)
            idxm1 = i - 1
            idxm1 = max0(idxm1,2)
            ucmona = ua(idxp1,k,j) + 2.*ua(idx,k,j) + ua(idxm1,k,j)
            vcmona = va(idx,k,jp1) + 2.*va(idx,k,j) + va(idx,k,jm1)
            ucmonb = ua(idxp1,k,jp1) + 2.*ua(idx,k,jp1)             &
                   & + ua(idxm1,k,jp1) + ucmona
            vcmonb = va(idxp1,k,jp1) + 2.*va(idxp1,k,j)             &
                   & + va(idxp1,k,jm1) + vcmona
            ucmonc = ua(idxp1,k,jm1) + 2.*ua(idx,k,jm1)             &
                   & + ua(idxm1,k,jm1) + ucmona
            vcmonc = va(idxm1,k,jp1) + 2.*va(idxm1,k,j)             &
                   & + va(idxm1,k,jm1) + vcmona
            ften(i,k) = ften(i,k)                                       &
                      & - ((u(i,k,jp1)+u(i,k,j))*ucmonb-(u(i,k,j)       &
                      & +u(i,k,jm1))*ucmonc+(u(i+1,k,j)+u(i,k,j))       &
                      & *vcmonb-(u(i,k,j)+u(i-1,k,j))*vcmonc)           &
                      & /(dxx*msfd(i,j)*msfd(i,j))
          end do
        end do
!
!----------------------------------------------------------------------
#else
        jdp1 = j + 1
        jdm1 = j - 1
#ifdef MPP1
        if ( myid.eq.0 ) jdm1 = max0(jdm1,2)
        if ( myid.eq.nproc-1 ) jdp1 = min0(jdp1,jendl-1)
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
            ucmona = ua(idxp1,k,j) + 2.*ua(idx,k,j)                 &
                   & + ua(idxm1,k,j)
            vcmona = va(idx,k,jdp1) + 2.*va(idx,k,j)                 &
                   & + va(idx,k,jdm1)
            ucmonb = ua(idxp1,k,jdp1) + 2.*ua(idx,k,jdp1)             &
                   & + ua(idxm1,k,jdp1) + ucmona
            vcmonb = va(idxp1,k,jdp1) + 2.*va(idxp1,k,j)             &
                   & + va(idxp1,k,jdm1) + vcmona
            ucmonc = ua(idxp1,k,jdm1) + 2.*ua(idx,k,jdm1)             &
                   & + ua(idxm1,k,jdm1) + ucmona
            vcmonc = va(idxm1,k,jdp1) + 2.*va(idxm1,k,j)             &
                   & + va(idxm1,k,jdm1) + vcmona
            ften(i,k) = ften(i,k)                                       &
                      & - ((u(i,k,jp1)+u(i,k,j))*ucmonb-(u(i,k,j)       &
                      & +u(i,k,jm1))*ucmonc+(u(i+1,k,j)+u(i,k,j))       &
                      & *vcmonb-(u(i,k,j)+u(i-1,k,j))*vcmonc)           &
                      & /(dxx*msfd(i,j)*msfd(i,j))
          end do
        end do
#endif
!
      else
        write(*,*) 'The U advection scheme ',ind, &
                 & ' you required is not available.'
        stop
      end if
      end subroutine hadv_u
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine hadv_v(ften,dxx,j,ind)
!
      implicit none
!
      real(8) ,intent (in) :: dxx
      integer ,intent (in) :: ind , j
      real(8) ,intent (inout), dimension(iy,kz) :: ften
!
      integer :: jm1 , jp1
#ifndef BAND
      integer :: jdm1 , jdp1
#endif
      integer :: i , k , idx , idxm1 , idxp1
      real(8) :: ucmona , ucmonb , ucmonc , vcmona , vcmonb , vcmonc
!
      jm1 = j - 1
      jp1 = j + 1
!----------------------------------------------------------------------
!
      if ( ind.eq.3 ) then
!
!-----for u and v:
!
#ifdef BAND
#if defined(BAND) && (!defined(MPP1))
        if(jm1.eq.0) jm1 = jx
        if(jp1.eq.jx+1) jp1 = 1
#endif
!
        do k = 1 , kz
          do i = 2 , iym1
            idx = i
            idxp1 = i + 1
            idxp1 = min0(idxp1,iym1)
            idxm1 = i - 1
            idxm1 = max0(idxm1,2)
            ucmona = ua(idxp1,k,j) + 2.*ua(idx,k,j) + ua(idxm1,k,j)
            vcmona = va(idx,k,jp1) + 2.*va(idx,k,j) + va(idx,k,jm1)
            ucmonb = ua(idxp1,k,jp1) + 2.*ua(idx,k,jp1)             &
                   & + ua(idxm1,k,jp1) + ucmona
            vcmonb = va(idxp1,k,jp1) + 2.*va(idxp1,k,j)             &
                   & + va(idxp1,k,jm1) + vcmona
            ucmonc = ua(idxp1,k,jm1) + 2.*ua(idx,k,jm1)             &
                   & + ua(idxm1,k,jm1) + ucmona
            vcmonc = va(idxm1,k,jp1) + 2.*va(idxm1,k,j)             &
                   & + va(idxm1,k,jm1) + vcmona
            ften(i,k) = ften(i,k)                                       &
                      & - ((v(i,k,jp1)+v(i,k,j))*ucmonb-(v(i,k,j)       &
                      & +v(i,k,jm1))*ucmonc+(v(i+1,k,j)+v(i,k,j))       &
                      & *vcmonb-(v(i,k,j)+v(i-1,k,j))*vcmonc)           &
                      & /(dxx*msfd(i,j)*msfd(i,j))
          end do
        end do
!
!----------------------------------------------------------------------
#else
        jdp1 = j + 1
        jdm1 = j - 1
#ifdef MPP1
        if ( myid.eq.0 ) jdm1 = max0(jdm1,2)
        if ( myid.eq.nproc-1 ) jdp1 = min0(jdp1,jendl-1)
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
            ucmona = ua(idxp1,k,j) + 2.*ua(idx,k,j)                 &
                   & + ua(idxm1,k,j)
            vcmona = va(idx,k,jdp1) + 2.*va(idx,k,j)                 &
                   & + va(idx,k,jdm1)
            ucmonb = ua(idxp1,k,jdp1) + 2.*ua(idx,k,jdp1)             &
                   & + ua(idxm1,k,jdp1) + ucmona
            vcmonb = va(idxp1,k,jdp1) + 2.*va(idxp1,k,j)             &
                   & + va(idxp1,k,jdm1) + vcmona
            ucmonc = ua(idxp1,k,jdm1) + 2.*ua(idx,k,jdm1)             &
                   & + ua(idxm1,k,jdm1) + ucmona
            vcmonc = va(idxm1,k,jdp1) + 2.*va(idxm1,k,j)             &
                   & + va(idxm1,k,jdm1) + vcmona
            ften(i,k) = ften(i,k)                                       &
                      & - ((v(i,k,jp1)+v(i,k,j))*ucmonb-(v(i,k,j)       &
                      & +v(i,k,jm1))*ucmonc+(v(i+1,k,j)+v(i,k,j))       &
                      & *vcmonb-(v(i,k,j)+v(i-1,k,j))*vcmonc)           &
                      & /(dxx*msfd(i,j)*msfd(i,j))
          end do
        end do
#endif
!
      else
        write(*,*) 'The V advection scheme ',ind, &
                 & ' you required is not available.'
        stop
      end if
      end subroutine hadv_v
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine hadvqv(ften,dxx,j,ind)
!
      implicit none
!
      real(8) ,intent (in) :: dxx
      integer ,intent (in) :: ind , j
      real(8) ,intent (inout), dimension(iy,kz) :: ften
!
      integer :: jm1 , jp1
      real(8) :: fact1 , fact2 , fx1 , fx2 , fy1 , fy2 , uavg1 , uavg2 ,&
               & vavg1 , vavg2
      integer :: i , k
!
      jm1 = j - 1
      jp1 = j + 1
#if defined(BAND) && (!defined(MPP1))
      if(jm1.eq.0) jm1 = jx
      if(jp1.eq.jx+1) jp1 = 1
#endif
!----------------------------------------------------------------------
!
      if ( ind.eq.1 ) then
!
!-----for t and qv:
!
        do k = 1 , kz
          do i = 2 , iym2
            ften(i,k) = ften(i,k)                                       &
                      & - ((ua(i+1,k,jp1)+ua(i,k,jp1))*(qv(i,k,jp1)     &
                      & +qv(i,k,j))-(ua(i+1,k,j)+ua(i,k,j))             &
                      & *(qv(i,k,j)+qv(i,k,jm1))                        &
                      & +(va(i+1,k,jp1)+va(i+1,k,j))                    &
                      & *(qv(i+1,k,j)+qv(i,k,j))-(va(i,k,jp1)+va(i,k,j))&
                      & *(qv(i-1,k,j)+qv(i,k,j)))                       &
                      & /(dxx*msfx(i,j)*msfx(i,j))
          end do
        end do
      else if ( ind.eq.2 ) then
!
!       implement a "relaxed" upstream scheme
!
!hy     fact1=0.75
        fact1 = 0.60
        fact2 = 1. - fact1
!
!-----for qc and qr:
!       up-wind values of qc and qr are used.
!
        do k = 1 , kz
          do i = 2 , iym2
            uavg2 = 0.5*(ua(i+1,k,jp1)+ua(i,k,jp1))
            uavg1 = 0.5*(ua(i+1,k,j)+ua(i,k,j))
            if ( uavg2.ge.0. ) then
              fx2 = fact1*qv(i,k,j) + fact2*qv(i,k,jp1)
            else
              fx2 = fact1*qv(i,k,jp1) + fact2*qv(i,k,j)
            end if
            if ( uavg1.ge.0. ) then
              fx1 = fact1*qv(i,k,jm1) + fact2*qv(i,k,j)
            else
              fx1 = fact1*qv(i,k,j) + fact2*qv(i,k,jm1)
            end if
            vavg2 = 0.5*(va(i+1,k,jp1)+va(i+1,k,j))
            vavg1 = 0.5*(va(i,k,jp1)+va(i,k,j))
            if ( vavg2.ge.0. ) then
              fy2 = fact1*qv(i,k,j) + fact2*qv(i+1,k,j)
            else
              fy2 = fact1*qv(i+1,k,j) + fact2*qv(i,k,j)
            end if
            if ( vavg1.ge.0. ) then
              fy1 = fact1*qv(i-1,k,j) + fact2*qv(i,k,j)
            else
              fy1 = fact1*qv(i,k,j) + fact2*qv(i-1,k,j)
            end if
            ften(i,k) = ften(i,k)                                       &
                      & - (uavg2*fx2-uavg1*fx1+vavg2*fy2-vavg1*fy1)     &
                      & /(dxx*msfx(i,j)*msfx(i,j))
          end do
        end do
      else
        write(*,*) 'The QV advection scheme ',ind, &
                 & ' you required is not available.'
        stop
      end if
      end subroutine hadvqv
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine hadvqc(ften,dxx,j,ind)
!
      implicit none
!
      real(8) ,intent (in) :: dxx
      integer ,intent (in) :: ind , j
      real(8) ,intent (inout), dimension(iy,kz) :: ften
!
      integer :: jm1 , jp1
      real(8) :: fact1 , fact2 , fx1 , fx2 , fy1 , fy2 , uavg1 , uavg2 ,&
               & vavg1 , vavg2
      integer :: i , k
!
      jm1 = j - 1
      jp1 = j + 1
#if defined(BAND) && (!defined(MPP1))
      if(jm1.eq.0) jm1 = jx
      if(jp1.eq.jx+1) jp1 = 1
#endif
!----------------------------------------------------------------------
!
      if ( ind.eq.1 ) then
!
!-----for t and qv:
!
        do k = 1 , kz
          do i = 2 , iym2
            ften(i,k) = ften(i,k)                                       &
                      & - ((ua(i+1,k,jp1)+ua(i,k,jp1))*(qc(i,k,jp1)     &
                      & +qc(i,k,j))-(ua(i+1,k,j)+ua(i,k,j))             &
                      & *(qc(i,k,j)+qc(i,k,jm1))                        &
                      & +(va(i+1,k,jp1)+va(i+1,k,j))                    &
                      & *(qc(i+1,k,j)+qc(i,k,j))-(va(i,k,jp1)+va(i,k,j))&
                      & *(qc(i-1,k,j)+qc(i,k,j)))                       &
                      & /(dxx*msfx(i,j)*msfx(i,j))
          end do
        end do
      else if ( ind.eq.2 ) then
!
!       implement a "relaxed" upstream scheme
!
!hy     fact1=0.75
        fact1 = 0.60
        fact2 = 1. - fact1
!
!-----for qc and qr:
!       up-wind values of qc and qr are used.
!
        do k = 1 , kz
          do i = 2 , iym2
            uavg2 = 0.5*(ua(i+1,k,jp1)+ua(i,k,jp1))
            uavg1 = 0.5*(ua(i+1,k,j)+ua(i,k,j))
            if ( uavg2.ge.0. ) then
              fx2 = fact1*qc(i,k,j) + fact2*qc(i,k,jp1)
            else
              fx2 = fact1*qc(i,k,jp1) + fact2*qc(i,k,j)
            end if
            if ( uavg1.ge.0. ) then
              fx1 = fact1*qc(i,k,jm1) + fact2*qc(i,k,j)
            else
              fx1 = fact1*qc(i,k,j) + fact2*qc(i,k,jm1)
            end if
            vavg2 = 0.5*(va(i+1,k,jp1)+va(i+1,k,j))
            vavg1 = 0.5*(va(i,k,jp1)+va(i,k,j))
            if ( vavg2.ge.0. ) then
              fy2 = fact1*qc(i,k,j) + fact2*qc(i+1,k,j)
            else
              fy2 = fact1*qc(i+1,k,j) + fact2*qc(i,k,j)
            end if
            if ( vavg1.ge.0. ) then
              fy1 = fact1*qc(i-1,k,j) + fact2*qc(i,k,j)
            else
              fy1 = fact1*qc(i,k,j) + fact2*qc(i-1,k,j)
            end if
            ften(i,k) = ften(i,k)                                       &
                      & - (uavg2*fx2-uavg1*fx1+vavg2*fy2-vavg1*fy1)     &
                      & /(dxx*msfx(i,j)*msfx(i,j))
          end do
        end do
      else
        write(*,*) 'The QC advection scheme ',ind, &
                 & ' you required is not available.'
        stop
      end if
      end subroutine hadvqc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine hadvch(ften,dxx,n,j,ind)
!
      implicit none
!
      real(8) ,intent (in) :: dxx
      integer ,intent (in) :: ind , j , n
      real(8) ,intent (inout), dimension(iy,kz) :: ften
!
      integer :: jm1 , jp1
      integer :: i , k
      real(8) :: fact1 , fact2 , fx1 , fx2 , fy1 , fy2 , uavg1 , uavg2 ,&
               & vavg1 , vavg2
!
      jm1 = j - 1
      jp1 = j + 1
#if defined(BAND) && (!defined(MPP1))
      if(jm1.eq.0) jm1 = jx
      if(jp1.eq.jx+1) jp1 = 1
#endif
!----------------------------------------------------------------------
!
      if ( ind.eq.2 ) then
!
!       implement a "relaxed" upstream scheme
!
!hy     fact1=0.75
        fact1 = 0.60
        fact2 = 1. - fact1
!
!-----for qc and qr:
!       up-wind values of qc and qr are used.
!
        do k = 1 , kz
          do i = 2 , iym2
            uavg2 = 0.5*(ua(i+1,k,jp1)+ua(i,k,jp1))
            uavg1 = 0.5*(ua(i+1,k,j)+ua(i,k,j))
            if ( uavg2.ge.0. ) then
              fx2 = fact1*chi(i,k,j,n) + fact2*chi(i,k,jp1,n)
            else
              fx2 = fact1*chi(i,k,jp1,n) + fact2*chi(i,k,j,n)
            end if
            if ( uavg1.ge.0. ) then
              fx1 = fact1*chi(i,k,jm1,n) + fact2*chi(i,k,j,n)
            else
              fx1 = fact1*chi(i,k,j,n) + fact2*chi(i,k,jm1,n)
            end if
            vavg2 = 0.5*(va(i+1,k,jp1)+va(i+1,k,j))
            vavg1 = 0.5*(va(i,k,jp1)+va(i,k,j))
            if ( vavg2.ge.0. ) then
              fy2 = fact1*chi(i,k,j,n) + fact2*chi(i+1,k,j,n)
            else
              fy2 = fact1*chi(i+1,k,j,n) + fact2*chi(i,k,j,n)
            end if
            if ( vavg1.ge.0. ) then
              fy1 = fact1*chi(i-1,k,j,n) + fact2*chi(i,k,j,n)
            else
              fy1 = fact1*chi(i,k,j,n) + fact2*chi(i-1,k,j,n)
            end if
            ften(i,k) = ften(i,k)                                       &
                      & - (uavg2*fx2-uavg1*fx1+vavg2*fy2-vavg1*fy1)     &
                      & /(dxx*msfx(i,j)*msfx(i,j))
          end do
        end do
      else
        write(*,*) 'The CH advection scheme ',ind, &
                 & ' you required is not available.'
        stop
      end if
      end subroutine hadvch
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
      subroutine vadv(ften,fa,j,ind)
!
      implicit none
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
#ifdef BAND
      integer :: jm1
!----------------------------------------------------------------------
#ifdef MPP1
      jm1 = j-1
#else
      jm1 = j-1
      if(jm1.eq.0) jm1=jx
#endif
!
!     qdot   : is the vertical sigma-velocity
!     fg     : is the working space used to store the interlated
!              values.
!     psa    : is p* used to interpolate the temperature.
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
                    & *((psa(i,j)*sigma(k)+r8pt)/(psa(i,j)*a(k)+r8pt))  &
                    & **0.287 + twt(k,2)*fa(i,k-1)                      &
                    & *((psa(i,j)*sigma(k)+r8pt)/(psa(i,j)*a(k-1)+r8pt))&
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
                    & - (qdot(i-1,2,jm1)+qdot(i,2,jm1)+qdot(i,2,j)      &
                    & +qdot(i-1,2,j))*fg(i,2)/(4.*dsigma(1))
        end do
!......k = 2,kzm1
        do k = 2 , kzm1
          do i = 2 , iym1
            ften(i,k) = ften(i,k)                                       &
                      & - ((qdot(i,k+1,jm1)+qdot(i-1,k+1,jm1)+qdot(i,   &
                      & k+1,j)+qdot(i-1,k+1,j))*fg(i,k+1)               &
                      & -(qdot(i,k,jm1)+qdot(i-1,k,jm1)+qdot(i,k,j)     &
                      & +qdot(i-1,k,j))*fg(i,k))/(4.*dsigma(k))
          end do
        end do
!......k = kz
        do i = 2 , iym1
          ften(i,kz) = ften(i,kz)                                       &
                     & + (qdot(i,kz,jm1)+qdot(i-1,kz,jm1)+qdot(i,kz,j)  &
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
#else
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
                    & *((psa(i,j)*sigma(k)+r8pt)/(psa(i,j)*a(k)+r8pt))  &
                    & **0.287 + twt(k,2)*fa(i,k-1)                      &
                    & *((psa(i,j)*sigma(k)+r8pt)/(psa(i,j)*a(k-1)+r8pt))&
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
#endif
      end subroutine vadv
!
      end module mod_advection
