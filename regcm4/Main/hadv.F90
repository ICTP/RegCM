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
 
      subroutine hadv_t(ften,dxx,j,ind)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the horizontal flux-divergence terms.  c
!     second-order difference is used.                                c
!                                                                     c
!     ften   : is the tendency for variable 'f'.                      c
!                                                                     c
!     ua, va : are p*u and p*v.                                       c
!                                                                     c
!     msfx   : is the map scale factor.                               c
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

      use mod_regcm_param
      use mod_cvaria
      use mod_main
      implicit none
!
! Dummy arguments
!
      real(8) :: dxx
      integer :: ind , j
      real(8) , dimension(iy,kz) :: ften
      intent (in) dxx , ind , j
      intent (inout) ften
!
! Local variables
!
      real(8) :: fact1 , fact2 , fx1 , fx2 , fy1 , fy2 , uavg1 , uavg2 ,&
               & ucmona , ucmonb , ucmonc , vavg1 , vavg2 , vcmona ,    &
               & vcmonb , vcmonc
      integer :: i , idx , idxm1 , idxp1 , jdx , jdxm1 , jdxp1 , k
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
                      & - ((ua(i+1,k,j+1)+ua(i,k,j+1))*(t(i,k,j+1)      &
                      & +t(i,k,j))-(ua(i+1,k,j)+ua(i,k,j))              &
                      & *(t(i,k,j)+t(i,k,j-1))                          &
                      & +(va(i+1,k,j+1)+va(i+1,k,j))                    &
                      & *(t(i+1,k,j)+t(i,k,j))-(va(i,k,j+1)+va(i,k,j))  &
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
            uavg2 = 0.5*(ua(i+1,k,j+1)+ua(i,k,j+1))
            uavg1 = 0.5*(ua(i+1,k,j)+ua(i,k,j))
            if ( uavg2.ge.0. ) then
              fx2 = fact1*t(i,k,j) + fact2*t(i,k,j+1)
            else
              fx2 = fact1*t(i,k,j+1) + fact2*t(i,k,j)
            end if
            if ( uavg1.ge.0. ) then
              fx1 = fact1*t(i,k,j-1) + fact2*t(i,k,j)
            else
              fx1 = fact1*t(i,k,j) + fact2*t(i,k,j-1)
            end if
            vavg2 = 0.5*(va(i+1,k,j+1)+va(i+1,k,j))
            vavg1 = 0.5*(va(i,k,j+1)+va(i,k,j))
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
!
      else if ( ind.eq.3 ) then
!
!-----for u and v:
!
        jdx = j
        jdxp1 = j + 1
        jdxm1 = j - 1
#ifdef MPP1
        if ( myid.eq.0 ) jdxm1 = max0(jdxm1,2)
        if ( myid.eq.nproc-1 ) jdxp1 = min0(jdxp1,jendl-1)
#else
        jdxp1 = min0(jdxp1,jxm1)
        jdxm1 = max0(jdxm1,2)
#endif
!
        do k = 1 , kz
          do i = 2 , iym1
            idx = i
            idxp1 = i + 1
            idxp1 = min0(idxp1,iym1)
            idxm1 = i - 1
            idxm1 = max0(idxm1,2)
            ucmona = ua(idxp1,k,jdx) + 2.*ua(idx,k,jdx)                 &
                   & + ua(idxm1,k,jdx)
            vcmona = va(idx,k,jdxp1) + 2.*va(idx,k,jdx)                 &
                   & + va(idx,k,jdxm1)
            ucmonb = ua(idxp1,k,jdxp1) + 2.*ua(idx,k,jdxp1)             &
                   & + ua(idxm1,k,jdxp1) + ucmona
            vcmonb = va(idxp1,k,jdxp1) + 2.*va(idxp1,k,jdx)             &
                   & + va(idxp1,k,jdxm1) + vcmona
            ucmonc = ua(idxp1,k,jdxm1) + 2.*ua(idx,k,jdxm1)             &
                   & + ua(idxm1,k,jdxm1) + ucmona
            vcmonc = va(idxm1,k,jdxp1) + 2.*va(idxm1,k,jdx)             &
                   & + va(idxm1,k,jdxm1) + vcmona
            ften(i,k) = ften(i,k)                                       &
                      & - ((t(i,k,j+1)+t(i,k,j))*ucmonb-(t(i,k,j)       &
                      & +t(i,k,j-1))*ucmonc+(t(i+1,k,j)+t(i,k,j))       &
                      & *vcmonb-(t(i,k,j)+t(i-1,k,j))*vcmonc)           &
                      & /(dxx*msfx(i,j)*msfx(i,j))
          end do
        end do
!
      else
      end if
!
      end subroutine hadv_t
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine hadv_u(ften,dxx,j,ind)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the horizontal flux-divergence terms.  c
!     second-order difference is used.                                c
!                                                                     c
!     ften   : is the tendency for variable 'f'.                      c
!                                                                     c
!     ua, va : are p*u and p*v.                                       c
!                                                                     c
!     msfd   : is the map scale factor.                               c
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

      use mod_regcm_param
      use mod_cvaria
      use mod_main
      implicit none
!
! Dummy arguments
!
      real(8) :: dxx
      integer :: ind , j
      real(8) , dimension(iy,kz) :: ften
      intent (in) dxx , ind , j
      intent (inout) ften
!
! Local variables
!
      real(8) :: fact1 , fact2 , fx1 , fx2 , fy1 , fy2 , uavg1 , uavg2 ,&
               & ucmona , ucmonb , ucmonc , vavg1 , vavg2 , vcmona ,    &
               & vcmonb , vcmonc
      integer :: i , idx , idxm1 , idxp1 , jdx , jdxm1 , jdxp1 , k
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
                      & - ((ua(i+1,k,j+1)+ua(i,k,j+1))*(u(i,k,j+1)      &
                      & +u(i,k,j))-(ua(i+1,k,j)+ua(i,k,j))              &
                      & *(u(i,k,j)+u(i,k,j-1))                          &
                      & +(va(i+1,k,j+1)+va(i+1,k,j))                    &
                      & *(u(i+1,k,j)+u(i,k,j))-(va(i,k,j+1)+va(i,k,j))  &
                      & *(u(i-1,k,j)+u(i,k,j)))                         &
                      & /(dxx*msfd(i,j)*msfd(i,j))
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
            uavg2 = 0.5*(ua(i+1,k,j+1)+ua(i,k,j+1))
            uavg1 = 0.5*(ua(i+1,k,j)+ua(i,k,j))
            if ( uavg2.ge.0. ) then
              fx2 = fact1*u(i,k,j) + fact2*u(i,k,j+1)
            else
              fx2 = fact1*u(i,k,j+1) + fact2*u(i,k,j)
            end if
            if ( uavg1.ge.0. ) then
              fx1 = fact1*u(i,k,j-1) + fact2*u(i,k,j)
            else
              fx1 = fact1*u(i,k,j) + fact2*u(i,k,j-1)
            end if
            vavg2 = 0.5*(va(i+1,k,j+1)+va(i+1,k,j))
            vavg1 = 0.5*(va(i,k,j+1)+va(i,k,j))
            if ( vavg2.ge.0. ) then
              fy2 = fact1*u(i,k,j) + fact2*u(i+1,k,j)
            else
              fy2 = fact1*u(i+1,k,j) + fact2*u(i,k,j)
            end if
            if ( vavg1.ge.0. ) then
              fy1 = fact1*u(i-1,k,j) + fact2*u(i,k,j)
            else
              fy1 = fact1*u(i,k,j) + fact2*u(i-1,k,j)
            end if
            ften(i,k) = ften(i,k)                                       &
                      & - (uavg2*fx2-uavg1*fx1+vavg2*fy2-vavg1*fy1)     &
                      & /(dxx*msfd(i,j)*msfd(i,j))
          end do
        end do
!
      else if ( ind.eq.3 ) then
!
!-----for u and v:
!
        jdx = j
        jdxp1 = j + 1
        jdxm1 = j - 1
#ifdef MPP1
        if ( myid.eq.0 ) jdxm1 = max0(jdxm1,2)
        if ( myid.eq.nproc-1 ) jdxp1 = min0(jdxp1,jendl-1)
#else
        jdxp1 = min0(jdxp1,jxm1)
        jdxm1 = max0(jdxm1,2)
#endif
!
        do k = 1 , kz
          do i = 2 , iym1
            idx = i
            idxp1 = i + 1
            idxp1 = min0(idxp1,iym1)
            idxm1 = i - 1
            idxm1 = max0(idxm1,2)
            ucmona = ua(idxp1,k,jdx) + 2.*ua(idx,k,jdx)                 &
                   & + ua(idxm1,k,jdx)
            vcmona = va(idx,k,jdxp1) + 2.*va(idx,k,jdx)                 &
                   & + va(idx,k,jdxm1)
            ucmonb = ua(idxp1,k,jdxp1) + 2.*ua(idx,k,jdxp1)             &
                   & + ua(idxm1,k,jdxp1) + ucmona
            vcmonb = va(idxp1,k,jdxp1) + 2.*va(idxp1,k,jdx)             &
                   & + va(idxp1,k,jdxm1) + vcmona
            ucmonc = ua(idxp1,k,jdxm1) + 2.*ua(idx,k,jdxm1)             &
                   & + ua(idxm1,k,jdxm1) + ucmona
            vcmonc = va(idxm1,k,jdxp1) + 2.*va(idxm1,k,jdx)             &
                   & + va(idxm1,k,jdxm1) + vcmona
            ften(i,k) = ften(i,k)                                       &
                      & - ((u(i,k,j+1)+u(i,k,j))*ucmonb-(u(i,k,j)       &
                      & +u(i,k,j-1))*ucmonc+(u(i+1,k,j)+u(i,k,j))       &
                      & *vcmonb-(u(i,k,j)+u(i-1,k,j))*vcmonc)           &
                      & /(dxx*msfd(i,j)*msfd(i,j))
          end do
        end do
!
      else
      end if
!
      end subroutine hadv_u
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine hadv_v(ften,dxx,j,ind)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the horizontal flux-divergence terms.  c
!     second-order difference is used.                                c
!                                                                     c
!     ften   : is the tendency for variable 'f'.                      c
!                                                                     c
!     ua, va : are p*u and p*v.                                       c
!                                                                     c
!     msfd   : is the map scale factor.                               c
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

      use mod_regcm_param
      use mod_cvaria
      use mod_main
      implicit none
!
! Dummy arguments
!
      real(8) :: dxx
      integer :: ind , j
      real(8) , dimension(iy,kz) :: ften
      intent (in) dxx , ind , j
      intent (inout) ften
!
! Local variables
!
      real(8) :: fact1 , fact2 , fx1 , fx2 , fy1 , fy2 , uavg1 , uavg2 ,&
               & ucmona , ucmonb , ucmonc , vavg1 , vavg2 , vcmona ,    &
               & vcmonb , vcmonc
      integer :: i , idx , idxm1 , idxp1 , jdx , jdxm1 , jdxp1 , k
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
                      & - ((ua(i+1,k,j+1)+ua(i,k,j+1))*(v(i,k,j+1)      &
                      & +v(i,k,j))-(ua(i+1,k,j)+ua(i,k,j))              &
                      & *(v(i,k,j)+v(i,k,j-1))                          &
                      & +(va(i+1,k,j+1)+va(i+1,k,j))                    &
                      & *(v(i+1,k,j)+v(i,k,j))-(va(i,k,j+1)+va(i,k,j))  &
                      & *(v(i-1,k,j)+v(i,k,j)))                         &
                      & /(dxx*msfd(i,j)*msfd(i,j))
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
            uavg2 = 0.5*(ua(i+1,k,j+1)+ua(i,k,j+1))
            uavg1 = 0.5*(ua(i+1,k,j)+ua(i,k,j))
            if ( uavg2.ge.0. ) then
              fx2 = fact1*v(i,k,j) + fact2*v(i,k,j+1)
            else
              fx2 = fact1*v(i,k,j+1) + fact2*v(i,k,j)
            end if
            if ( uavg1.ge.0. ) then
              fx1 = fact1*v(i,k,j-1) + fact2*v(i,k,j)
            else
              fx1 = fact1*v(i,k,j) + fact2*v(i,k,j-1)
            end if
            vavg2 = 0.5*(va(i+1,k,j+1)+va(i+1,k,j))
            vavg1 = 0.5*(va(i,k,j+1)+va(i,k,j))
            if ( vavg2.ge.0. ) then
              fy2 = fact1*v(i,k,j) + fact2*v(i+1,k,j)
            else
              fy2 = fact1*v(i+1,k,j) + fact2*v(i,k,j)
            end if
            if ( vavg1.ge.0. ) then
              fy1 = fact1*v(i-1,k,j) + fact2*v(i,k,j)
            else
              fy1 = fact1*v(i,k,j) + fact2*v(i-1,k,j)
            end if
            ften(i,k) = ften(i,k)                                       &
                      & - (uavg2*fx2-uavg1*fx1+vavg2*fy2-vavg1*fy1)     &
                      & /(dxx*msfd(i,j)*msfd(i,j))
          end do
        end do
!
      else if ( ind.eq.3 ) then
!
!-----for u and v:
!
        jdx = j
        jdxp1 = j + 1
        jdxm1 = j - 1
#ifdef MPP1
        if ( myid.eq.0 ) jdxm1 = max0(jdxm1,2)
        if ( myid.eq.nproc-1 ) jdxp1 = min0(jdxp1,jendl-1)
#else
        jdxp1 = min0(jdxp1,jxm1)
        jdxm1 = max0(jdxm1,2)
#endif
!
        do k = 1 , kz
          do i = 2 , iym1
            idx = i
            idxp1 = i + 1
            idxp1 = min0(idxp1,iym1)
            idxm1 = i - 1
            idxm1 = max0(idxm1,2)
            ucmona = ua(idxp1,k,jdx) + 2.*ua(idx,k,jdx)                 &
                   & + ua(idxm1,k,jdx)
            vcmona = va(idx,k,jdxp1) + 2.*va(idx,k,jdx)                 &
                   & + va(idx,k,jdxm1)
            ucmonb = ua(idxp1,k,jdxp1) + 2.*ua(idx,k,jdxp1)             &
                   & + ua(idxm1,k,jdxp1) + ucmona
            vcmonb = va(idxp1,k,jdxp1) + 2.*va(idxp1,k,jdx)             &
                   & + va(idxp1,k,jdxm1) + vcmona
            ucmonc = ua(idxp1,k,jdxm1) + 2.*ua(idx,k,jdxm1)             &
                   & + ua(idxm1,k,jdxm1) + ucmona
            vcmonc = va(idxm1,k,jdxp1) + 2.*va(idxm1,k,jdx)             &
                   & + va(idxm1,k,jdxm1) + vcmona
            ften(i,k) = ften(i,k)                                       &
                      & - ((v(i,k,j+1)+v(i,k,j))*ucmonb-(v(i,k,j)       &
                      & +v(i,k,j-1))*ucmonc+(v(i+1,k,j)+v(i,k,j))       &
                      & *vcmonb-(v(i,k,j)+v(i-1,k,j))*vcmonc)           &
                      & /(dxx*msfd(i,j)*msfd(i,j))
          end do
        end do
!
      else
      end if
!
      end subroutine hadv_v
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine hadvqv(ften,dxx,j,ind)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the horizontal flux-divergence terms.  c
!     second-order difference is used.                                c
!                                                                     c
!     ften   : is the tendency for variable 'f'.                      c
!                                                                     c
!     ua, va : are p*u and p*v.                                       c
!                                                                     c
!     msfx   : is the map scale factor.                               c
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

      use mod_regcm_param
      use mod_cvaria
      use mod_main
      implicit none
!
! Dummy arguments
!
      real(8) :: dxx
      integer :: ind , j
      real(8) , dimension(iy,kz) :: ften
      intent (in) dxx , ind , j
      intent (inout) ften
!
! Local variables
!
      real(8) :: fact1 , fact2 , fx1 , fx2 , fy1 , fy2 , uavg1 , uavg2 ,&
               & ucmona , ucmonb , ucmonc , vavg1 , vavg2 , vcmona ,    &
               & vcmonb , vcmonc
      integer :: i , idx , idxm1 , idxp1 , jdx , jdxm1 , jdxp1 , k
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
                      & - ((ua(i+1,k,j+1)+ua(i,k,j+1))*(qv(i,k,j+1)     &
                      & +qv(i,k,j))-(ua(i+1,k,j)+ua(i,k,j))             &
                      & *(qv(i,k,j)+qv(i,k,j-1))                        &
                      & +(va(i+1,k,j+1)+va(i+1,k,j))                    &
                      & *(qv(i+1,k,j)+qv(i,k,j))-(va(i,k,j+1)+va(i,k,j))&
                      & *(qv(i-1,k,j)+qv(i,k,j)))                       &
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
            uavg2 = 0.5*(ua(i+1,k,j+1)+ua(i,k,j+1))
            uavg1 = 0.5*(ua(i+1,k,j)+ua(i,k,j))
            if ( uavg2.ge.0. ) then
              fx2 = fact1*qv(i,k,j) + fact2*qv(i,k,j+1)
            else
              fx2 = fact1*qv(i,k,j+1) + fact2*qv(i,k,j)
            end if
            if ( uavg1.ge.0. ) then
              fx1 = fact1*qv(i,k,j-1) + fact2*qv(i,k,j)
            else
              fx1 = fact1*qv(i,k,j) + fact2*qv(i,k,j-1)
            end if
            vavg2 = 0.5*(va(i+1,k,j+1)+va(i+1,k,j))
            vavg1 = 0.5*(va(i,k,j+1)+va(i,k,j))
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
!
      else if ( ind.eq.3 ) then
!
!-----for u and v:
!
        jdx = j
        jdxp1 = j + 1
        jdxm1 = j - 1
#ifdef MPP1
        if ( myid.eq.0 ) jdxm1 = max0(jdxm1,2)
        if ( myid.eq.nproc-1 ) jdxp1 = min0(jdxp1,jendl-1)
#else
        jdxp1 = min0(jdxp1,jxm1)
        jdxm1 = max0(jdxm1,2)
#endif
!
        do k = 1 , kz
          do i = 2 , iym1
            idx = i
            idxp1 = i + 1
            idxp1 = min0(idxp1,iym1)
            idxm1 = i - 1
            idxm1 = max0(idxm1,2)
            ucmona = ua(idxp1,k,jdx) + 2.*ua(idx,k,jdx)                 &
                   & + ua(idxm1,k,jdx)
            vcmona = va(idx,k,jdxp1) + 2.*va(idx,k,jdx)                 &
                   & + va(idx,k,jdxm1)
            ucmonb = ua(idxp1,k,jdxp1) + 2.*ua(idx,k,jdxp1)             &
                   & + ua(idxm1,k,jdxp1) + ucmona
            vcmonb = va(idxp1,k,jdxp1) + 2.*va(idxp1,k,jdx)             &
                   & + va(idxp1,k,jdxm1) + vcmona
            ucmonc = ua(idxp1,k,jdxm1) + 2.*ua(idx,k,jdxm1)             &
                   & + ua(idxm1,k,jdxm1) + ucmona
            vcmonc = va(idxm1,k,jdxp1) + 2.*va(idxm1,k,jdx)             &
                   & + va(idxm1,k,jdxm1) + vcmona
            ften(i,k) = ften(i,k)                                       &
                      & - ((qv(i,k,j+1)+qv(i,k,j))*ucmonb-(qv(i,k,j)    &
                      & +qv(i,k,j-1))*ucmonc+(qv(i+1,k,j)+qv(i,k,j))    &
                      & *vcmonb-(qv(i,k,j)+qv(i-1,k,j))*vcmonc)         &
                      & /(dxx*msfx(i,j)*msfx(i,j))
          end do
        end do
!
      else
      end if
!
      end subroutine hadvqv
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine hadvqc(ften,dxx,j,ind)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the horizontal flux-divergence terms.  c
!     second-order difference is used.                                c
!                                                                     c
!     ften   : is the tendency for variable 'f'.                      c
!                                                                     c
!     ua, va : are p*u and p*v.                                       c
!                                                                     c
!     msfx   : is the map scale factor.                               c
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

      use mod_regcm_param
      use mod_cvaria
      use mod_main
      implicit none
!
! Dummy arguments
!
      real(8) :: dxx
      integer :: ind , j
      real(8) , dimension(iy,kz) :: ften
      intent (in) dxx , ind , j
      intent (inout) ften
!
! Local variables
!
      real(8) :: fact1 , fact2 , fx1 , fx2 , fy1 , fy2 , uavg1 , uavg2 ,&
               & ucmona , ucmonb , ucmonc , vavg1 , vavg2 , vcmona ,    &
               & vcmonb , vcmonc
      integer :: i , idx , idxm1 , idxp1 , jdx , jdxm1 , jdxp1 , k
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
                      & - ((ua(i+1,k,j+1)+ua(i,k,j+1))*(qc(i,k,j+1)     &
                      & +qc(i,k,j))-(ua(i+1,k,j)+ua(i,k,j))             &
                      & *(qc(i,k,j)+qc(i,k,j-1))                        &
                      & +(va(i+1,k,j+1)+va(i+1,k,j))                    &
                      & *(qc(i+1,k,j)+qc(i,k,j))-(va(i,k,j+1)+va(i,k,j))&
                      & *(qc(i-1,k,j)+qc(i,k,j)))                       &
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
            uavg2 = 0.5*(ua(i+1,k,j+1)+ua(i,k,j+1))
            uavg1 = 0.5*(ua(i+1,k,j)+ua(i,k,j))
            if ( uavg2.ge.0. ) then
              fx2 = fact1*qc(i,k,j) + fact2*qc(i,k,j+1)
            else
              fx2 = fact1*qc(i,k,j+1) + fact2*qc(i,k,j)
            end if
            if ( uavg1.ge.0. ) then
              fx1 = fact1*qc(i,k,j-1) + fact2*qc(i,k,j)
            else
              fx1 = fact1*qc(i,k,j) + fact2*qc(i,k,j-1)
            end if
            vavg2 = 0.5*(va(i+1,k,j+1)+va(i+1,k,j))
            vavg1 = 0.5*(va(i,k,j+1)+va(i,k,j))
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
!
      else if ( ind.eq.3 ) then
!
!-----for u and v:
!
        jdx = j
        jdxp1 = j + 1
        jdxm1 = j - 1
#ifdef MPP1
        if ( myid.eq.0 ) jdxm1 = max0(jdxm1,2)
        if ( myid.eq.nproc-1 ) jdxp1 = min0(jdxp1,jendl-1)
#else
        jdxp1 = min0(jdxp1,jxm1)
        jdxm1 = max0(jdxm1,2)
#endif
!
        do k = 1 , kz
          do i = 2 , iym1
            idx = i
            idxp1 = i + 1
            idxp1 = min0(idxp1,iym1)
            idxm1 = i - 1
            idxm1 = max0(idxm1,2)
            ucmona = ua(idxp1,k,jdx) + 2.*ua(idx,k,jdx)                 &
                   & + ua(idxm1,k,jdx)
            vcmona = va(idx,k,jdxp1) + 2.*va(idx,k,jdx)                 &
                   & + va(idx,k,jdxm1)
            ucmonb = ua(idxp1,k,jdxp1) + 2.*ua(idx,k,jdxp1)             &
                   & + ua(idxm1,k,jdxp1) + ucmona
            vcmonb = va(idxp1,k,jdxp1) + 2.*va(idxp1,k,jdx)             &
                   & + va(idxp1,k,jdxm1) + vcmona
            ucmonc = ua(idxp1,k,jdxm1) + 2.*ua(idx,k,jdxm1)             &
                   & + ua(idxm1,k,jdxm1) + ucmona
            vcmonc = va(idxm1,k,jdxp1) + 2.*va(idxm1,k,jdx)             &
                   & + va(idxm1,k,jdxm1) + vcmona
            ften(i,k) = ften(i,k)                                       &
                      & - ((qc(i,k,j+1)+qc(i,k,j))*ucmonb-(qc(i,k,j)    &
                      & +qc(i,k,j-1))*ucmonc+(qc(i+1,k,j)+qc(i,k,j))    &
                      & *vcmonb-(qc(i,k,j)+qc(i-1,k,j))*vcmonc)         &
                      & /(dxx*msfx(i,j)*msfx(i,j))
          end do
        end do
!
      else
      end if
!
      end subroutine hadvqc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine hadvch(ften,dxx,n,j,ind)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the horizontal flux-divergence terms.  c
!     second-order difference is used.                                c
!                                                                     c
!     ften   : is the tendency for variable 'f'.                      c
!                                                                     c
!     ua, va : are p*u and p*v.                                       c
!                                                                     c
!     msfx   : is the map scale factor.                               c
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

      use mod_regcm_param
      use mod_cvaria
      use mod_main
      implicit none
!
! Dummy arguments
!
      real(8) :: dxx
      integer :: ind , j , n
      real(8) , dimension(iy,kz) :: ften
      intent (in) dxx , ind , j , n
      intent (inout) ften
!
! Local variables
!
      real(8) :: fact1 , fact2 , fx1 , fx2 , fy1 , fy2 , uavg1 , uavg2 ,&
               & ucmona , ucmonb , ucmonc , vavg1 , vavg2 , vcmona ,    &
               & vcmonb , vcmonc
      integer :: i , idx , idxm1 , idxp1 , jdx , jdxm1 , jdxp1 , k
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
                      & - ((ua(i+1,k,j+1)+ua(i,k,j+1))*(chi(i,k,j+1,n)  &
                      & +chi(i,k,j,n))-(ua(i+1,k,j)+ua(i,k,j))          &
                      & *(chi(i,k,j,n)+chi(i,k,j-1,n))                  &
                      & +(va(i+1,k,j+1)+va(i+1,k,j))                    &
                      & *(chi(i+1,k,j,n)+chi(i,k,j,n))                  &
                      & -(va(i,k,j+1)+va(i,k,j))                        &
                      & *(chi(i-1,k,j,n)+chi(i,k,j,n)))                 &
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
            uavg2 = 0.5*(ua(i+1,k,j+1)+ua(i,k,j+1))
            uavg1 = 0.5*(ua(i+1,k,j)+ua(i,k,j))
            if ( uavg2.ge.0. ) then
              fx2 = fact1*chi(i,k,j,n) + fact2*chi(i,k,j+1,n)
            else
              fx2 = fact1*chi(i,k,j+1,n) + fact2*chi(i,k,j,n)
            end if
            if ( uavg1.ge.0. ) then
              fx1 = fact1*chi(i,k,j-1,n) + fact2*chi(i,k,j,n)
            else
              fx1 = fact1*chi(i,k,j,n) + fact2*chi(i,k,j-1,n)
            end if
            vavg2 = 0.5*(va(i+1,k,j+1)+va(i+1,k,j))
            vavg1 = 0.5*(va(i,k,j+1)+va(i,k,j))
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
!
      else if ( ind.eq.3 ) then
!
!-----for u and v:
!
        jdx = j
        jdxp1 = j + 1
        jdxm1 = j - 1
#ifdef MPP1
        if ( myid.eq.0 ) jdxm1 = max0(jdxm1,2)
        if ( myid.eq.nproc-1 ) jdxp1 = min0(jdxp1,jendl-1)
#else
        jdxp1 = min0(jdxp1,jxm1)
        jdxm1 = max0(jdxm1,2)
#endif
!
        do k = 1 , kz
          do i = 2 , iym1
            idx = i
            idxp1 = i + 1
            idxp1 = min0(idxp1,iym1)
            idxm1 = i - 1
            idxm1 = max0(idxm1,2)
            ucmona = ua(idxp1,k,jdx) + 2.*ua(idx,k,jdx)                 &
                   & + ua(idxm1,k,jdx)
            vcmona = va(idx,k,jdxp1) + 2.*va(idx,k,jdx)                 &
                   & + va(idx,k,jdxm1)
            ucmonb = ua(idxp1,k,jdxp1) + 2.*ua(idx,k,jdxp1)             &
                   & + ua(idxm1,k,jdxp1) + ucmona
            vcmonb = va(idxp1,k,jdxp1) + 2.*va(idxp1,k,jdx)             &
                   & + va(idxp1,k,jdxm1) + vcmona
            ucmonc = ua(idxp1,k,jdxm1) + 2.*ua(idx,k,jdxm1)             &
                   & + ua(idxm1,k,jdxm1) + ucmona
            vcmonc = va(idxm1,k,jdxp1) + 2.*va(idxm1,k,jdx)             &
                   & + va(idxm1,k,jdxm1) + vcmona
            ften(i,k) = ften(i,k)                                       &
                      & - ((chi(i,k,j+1,n)+chi(i,k,j,n))*ucmonb-        &
                      & (chi(i,k,j,n)+chi(i,k,j-1,n))                   &
                      & *ucmonc+(chi(i+1,k,j,n)+chi(i,k,j,n))           &
                      & *vcmonb-(chi(i,k,j,n)+chi(i-1,k,j,n))*vcmonc)   &
                      & /(dxx*msfx(i,j)*msfx(i,j))
          end do
        end do
!
      else
      end if
!
      end subroutine hadvch
