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
 
      subroutine sponge_p(ip,wg,ften,j)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine applies sponge boundary condition to the        c
!     tendency term - ften.                                           c
!                                                                     c
!     ip   : is the number of slices affected by sponge boundary.     c
!                                                                     c
!     wg   : are the weightings.                                      c
!                                                                     c
!     ften : is the tendency calculated from the model.               c
!                                                                     c
!     pebt, pwbt, pnbt, psbt : are the large-scale or observed        c
!            tendencies at east, west, north, and south boundaries.   c
!                                                                     c
!     ie = iy, je = jx for dot-point variables.                       c
!     ie = iym1, je = jxm1 for cross-point variables.                 c
!                                                                     c
!     j    : is the j'th slice of the tendency to be adjusted.        c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mod_dynparam
      use mod_bdycod
      implicit none
!
! Dummy arguments
!
      integer :: ip , j
      real(8) , dimension(iy) :: ften
      real(8) , dimension(ip) :: wg
      intent (in) ip , j , wg
      intent (inout) ften
!
! Local variables
!
      integer :: i , ibeg , iend , ii , jj , jsls
#ifdef MPP1
      integer :: jwb , jeb
#endif
!
#ifdef BAND
!----------------------------------------------------------------------
!
!-----interior j slices:
      do i = 2 , ip
         ii = iy - i
!.......south boundary:
         ften(i) = wg(i)*ften(i) + (1.-wg(i))*psbt(i,j)
!.......north boundary:
         ften(ii) = wg(i)*ften(ii) + (1.-wg(i))*pnbt(i,j)
      end do

#else
!----------------------------------------------------------------------
!
#ifdef MPP1
      jsls = j + myid*jxp
      jj = jx - jsls
      if ( jj.le.ip ) jsls = jj
      jwb = jsls
      if ( jwb.gt.jxp ) jwb = mod(jwb,jxp)
      if ( jwb.eq.0 ) jwb = jxp
      if ( myid.eq.nproc-1 ) then
        jeb = jsls
      else
        jeb = jsls + 1
      end if
      if ( jeb.gt.jxp ) jeb = mod(jeb,jxp)
      if ( jeb.eq.0 ) jeb = jxp
#else
      jsls = j
      jj = jx - jsls
      if ( jj.le.ip ) jsls = jj
#endif
!
      if ( jsls.gt.ip ) then
!-----interior j slices:
        do i = 2 , ip
          ii = iy - i
!.......south boundary:
          ften(i) = wg(i)*ften(i) + (1.-wg(i))*psbt(i,j)
!.......north boundary:
          ften(ii) = wg(i)*ften(ii) + (1.-wg(i))*pnbt(i,j)
        end do
!
      else if ( jsls.le.ip ) then
        ibeg = 2
        iend = iym1 - 1
        if ( jsls.gt.2 ) then
          do i = 2 , jsls - 1
            ii = iy - i
!........south boundary:
            ften(i) = wg(i)*ften(i) + (1.-wg(i))*psbt(i,j)
!........north boundary:
            ften(ii) = wg(i)*ften(ii) + (1.-wg(i))*pnbt(i,j)
          end do
          ibeg = jsls
          iend = iy - jsls
        end if
!
        if ( jj.gt.ip ) then
!------west-boundary slice:
          do i = ibeg , iend
#ifdef MPP1
            if ( jsls.le.ip ) ften(i) = wg(jsls)*ften(i) + (1.-wg(jsls))&
                                      & *pwbt(i,jwb)
#else
            ften(i) = wg(jsls)*ften(i) + (1.-wg(jsls))*pwbt(i,jsls)
#endif
          end do
        else if ( jj.le.ip ) then
!------east-boundary slice:
          do i = ibeg , iend
#ifdef MPP1
            if ( jsls.le.ip ) ften(i) = wg(jsls)*ften(i) + (1.-wg(jsls))&
                                      & *pebt(i,jeb)
#else
            ften(i) = wg(jsls)*ften(i) + (1.-wg(jsls))*pebt(i,jsls)
#endif
          end do
        else
        end if
!
      else
      end if

#endif
      end subroutine sponge_p
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine sponge_t(ip,wg,ften,j)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine applies sponge boundary condition to the        c
!     tendency term - ften.                                           c
!                                                                     c
!     ip   : is the number of slices affected by sponge boundary.     c
!                                                                     c
!     wg   : are the weightings.                                      c
!                                                                     c
!     ften : is the tendency calculated from the model.               c
!                                                                     c
!     tebt, twbt, tnbt, tsbt : are the large-scale or observed        c
!            tendencies at east, west, north, and south boundaries.   c
!                                                                     c
!     ie = iy, je = jx for dot-point variables.                       c
!     ie = iym1, je = jxm1 for cross-point variables.                 c
!                                                                     c
!     j    : is the j'th slice of the tendency to be adjusted.        c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_dynparam
      use mod_bdycod
      implicit none
!
! Dummy arguments
!
      integer :: ip , j
      real(8) , dimension(iy,kz) :: ften
      real(8) , dimension(ip) :: wg
      intent (in) ip , j , wg
      intent (inout) ften
!
! Local variables
!
      integer :: i , ibeg , iend , ii , jj , jsls , k
#ifdef MPP1
      integer :: jwb , jeb
#endif
!
#ifdef BAND
!----------------------------------------------------------------------
!
!-----interior j slices:
      do i = 2 , ip
         ii = iy - i
         do k = 1 , kz
!.......south boundary:
            ften(i,k) = wg(i)*ften(i,k) + (1.-wg(i))*tsbt(i,k,j)
!.......north boundary:
            ften(ii,k) = wg(i)*ften(ii,k) + (1.-wg(i))*tnbt(i,k,j)
         end do
      end do

#else
!----------------------------------------------------------------------
!
#ifdef MPP1
      jsls = j + myid*jxp
      jj = jx - jsls
      if ( jj.le.ip ) jsls = jj
      jwb = jsls
      if ( jwb.gt.jxp ) jwb = mod(jwb,jxp)
      if ( jwb.eq.0 ) jwb = jxp
      if ( myid.eq.nproc-1 ) then
        jeb = jsls
      else
        jeb = jsls + 1
      end if
      if ( jeb.gt.jxp ) jeb = mod(jeb,jxp)
      if ( jeb.eq.0 ) jeb = jxp
#else
      jsls = j
      jj = jx - jsls
      if ( jj.le.ip ) jsls = jj
#endif
!
      if ( jsls.gt.ip ) then
!-----interior j slices:
        do i = 2 , ip
          ii = iy - i
          do k = 1 , kz
!.......south boundary:
            ften(i,k) = wg(i)*ften(i,k) + (1.-wg(i))*tsbt(i,k,j)
!.......north boundary:
            ften(ii,k) = wg(i)*ften(ii,k) + (1.-wg(i))*tnbt(i,k,j)
          end do
        end do
!
      else if ( jsls.le.ip ) then
        ibeg = 2
        iend = iym1 - 1
        if ( jsls.gt.2 ) then
          do i = 2 , jsls - 1
            ii = iy - i
            do k = 1 , kz
!........south boundary:
              ften(i,k) = wg(i)*ften(i,k) + (1.-wg(i))*tsbt(i,k,j)
!........north boundary:
              ften(ii,k) = wg(i)*ften(ii,k) + (1.-wg(i))*tnbt(i,k,j)
            end do
          end do
          ibeg = jsls
          iend = iy - jsls
        end if
!
        if ( jj.gt.ip ) then
!------west-boundary slice:
          do k = 1 , kz
            do i = ibeg , iend
#ifdef MPP1
              if ( jsls.le.ip ) ften(i,k) = wg(jsls)*ften(i,k)          &
                 & + (1.-wg(jsls))*twbt(i,k,jwb)
#else
              ften(i,k) = wg(jsls)*ften(i,k) + (1.-wg(jsls))            &
                        & *twbt(i,k,jsls)
#endif
            end do
          end do
        else if ( jj.le.ip ) then
!------east-boundary slice:
          do k = 1 , kz
            do i = ibeg , iend
#ifdef MPP1
              if ( jsls.le.ip ) ften(i,k) = wg(jsls)*ften(i,k)          &
                 & + (1.-wg(jsls))*tebt(i,k,jeb)
#else
              ften(i,k) = wg(jsls)*ften(i,k) + (1.-wg(jsls))            &
                        & *tebt(i,k,jsls)
#endif
            end do
          end do
        else
        end if
!
      else
      end if
#endif

      end subroutine sponge_t
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine spongeqv(ip,wg,ften,j)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine applies sponge boundary condition to the        c
!     tendency term - ften.                                           c
!                                                                     c
!     ip   : is the number of slices affected by sponge boundary.     c
!                                                                     c
!     wg   : are the weightings.                                      c
!                                                                     c
!     ften : is the tendency calculated from the model.               c
!                                                                     c
!     qebt, qwbt, qnbt, qsbt : are the large-scale or observed        c
!            tendencies at east, west, north, and south boundaries.   c
!                                                                     c
!     ie = iy, je = jx for dot-point variables.                       c
!     ie = iym1, je = jxm1 for cross-point variables.                 c
!                                                                     c
!     j    : is the j'th slice of the tendency to be adjusted.        c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_dynparam
      use mod_bdycod
      implicit none
!
! Dummy arguments
!
      integer :: ip , j
      real(8) , dimension(iy,kz) :: ften
      real(8) , dimension(ip) :: wg
      intent (in) ip , j , wg
      intent (inout) ften
!
! Local variables
!
      integer :: i , ibeg , iend , ii , jj , jsls , k
#ifdef MPP1
      integer :: jwb , jeb
#endif
!
#ifdef BAND
!----------------------------------------------------------------------
!
!-----interior j slices:
      do i = 2 , ip
         ii = iy - i
         do k = 1 , kz
!.......south boundary:
            ften(i,k) = wg(i)*ften(i,k) + (1.-wg(i))*qsbt(i,k,j)
!.......north boundary:
            ften(ii,k) = wg(i)*ften(ii,k) + (1.-wg(i))*qnbt(i,k,j)
         end do
      end do
#else
!----------------------------------------------------------------------
!
#ifdef MPP1
      jsls = j + myid*jxp
      jj = jx - jsls
      if ( jj.le.ip ) jsls = jj
      jwb = jsls
      if ( jwb.gt.jxp ) jwb = mod(jwb,jxp)
      if ( jwb.eq.0 ) jwb = jxp
      if ( myid.eq.nproc-1 ) then
        jeb = jsls
      else
        jeb = jsls + 1
      end if
      if ( jeb.gt.jxp ) jeb = mod(jeb,jxp)
      if ( jeb.eq.0 ) jeb = jxp
#else
      jsls = j
      jj = jx - jsls
      if ( jj.le.ip ) jsls = jj
#endif
!
      if ( jsls.gt.ip ) then
!-----interior j slices:
        do i = 2 , ip
          ii = iy - i
          do k = 1 , kz
!.......south boundary:
            ften(i,k) = wg(i)*ften(i,k) + (1.-wg(i))*qsbt(i,k,j)
!.......north boundary:
            ften(ii,k) = wg(i)*ften(ii,k) + (1.-wg(i))*qnbt(i,k,j)
          end do
        end do
!
      else if ( jsls.le.ip ) then
        ibeg = 2
        iend = iym1 - 1
        if ( jsls.gt.2 ) then
          do i = 2 , jsls - 1
            ii = iy - i
            do k = 1 , kz
!........south boundary:
              ften(i,k) = wg(i)*ften(i,k) + (1.-wg(i))*qsbt(i,k,j)
!........north boundary:
              ften(ii,k) = wg(i)*ften(ii,k) + (1.-wg(i))*qnbt(i,k,j)
            end do
          end do
          ibeg = jsls
          iend = iy - jsls
        end if
!
        if ( jj.gt.ip ) then
!------west-boundary slice:
          do k = 1 , kz
            do i = ibeg , iend
#ifdef MPP1
              if ( jsls.le.ip ) ften(i,k) = wg(jsls)*ften(i,k)          &
                 & + (1.-wg(jsls))*qwbt(i,k,jwb)
#else
              ften(i,k) = wg(jsls)*ften(i,k) + (1.-wg(jsls))            &
                        & *qwbt(i,k,jsls)
#endif
            end do
          end do
        else if ( jj.le.ip ) then
!------east-boundary slice:
          do k = 1 , kz
            do i = ibeg , iend
#ifdef MPP1
              if ( jsls.le.ip ) ften(i,k) = wg(jsls)*ften(i,k)          &
                 & + (1.-wg(jsls))*qebt(i,k,jeb)
#else
              ften(i,k) = wg(jsls)*ften(i,k) + (1.-wg(jsls))            &
                        & *qebt(i,k,jsls)
#endif
            end do
          end do
        else
        end if
!
      else
      end if
#endif
      end subroutine spongeqv
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine sponge_u(ip,wg,ften,j)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine applies sponge boundary condition to the        c
!     tendency term - ften.                                           c
!                                                                     c
!     ip   : is the number of slices affected by sponge boundary.     c
!                                                                     c
!     wg   : are the weightings.                                      c
!                                                                     c
!     ften : is the tendency calculated from the model.               c
!                                                                     c
!     uebt, uwbt, unbt, usbt : are the large-scale or observed        c
!            tendencies at east, west, north, and south boundaries.   c
!                                                                     c
!     ie = iy, je = jx for dot-point variables.                       c
!     ie = iym1, je = jxm1 for cross-point variables.                 c
!                                                                     c
!     j    : is the j'th slice of the tendency to be adjusted.        c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_dynparam
      use mod_bdycod
      implicit none
!
! Dummy arguments
!
      integer :: ip , j
      real(8) , dimension(iy,kz) :: ften
      real(8) , dimension(ip) :: wg
      intent (in) ip , j , wg
      intent (inout) ften
!
! Local variables
!
      integer :: i , ibeg , iend , ii , jj , jsls , k
#ifdef MPP1
      integer :: jew
#endif

!
#ifdef BAND
!----------------------------------------------------------------------
!
!-----interior j slices:
      do i = 2 , ip
         ii = iy - i + 1
         do k = 1 , kz
!.......south boundary:
            ften(i,k) = wg(i)*ften(i,k) + (1.-wg(i))*usbt(i,k,j)
!.......north boundary:
            ften(ii,k) = wg(i)*ften(ii,k) + (1.-wg(i))*unbt(i,k,j)
         end do
      end do
#else
!----------------------------------------------------------------------
!
#ifdef MPP1
      jsls = j + myid*jxp
      jj = jxp1 - jsls
      if ( jj.le.ip ) jsls = jj
      jew = jsls
      if ( jew.gt.jxp ) jew = mod(jsls,jxp)
      if ( jew.eq.0 ) jew = jxp
#else
      jsls = j
      jj = jxp1 - jsls
      if ( jj.le.ip ) jsls = jj
#endif
!
      if ( jsls.gt.ip ) then
!-----interior j slices:
        do i = 2 , ip
          ii = iy - i + 1
          do k = 1 , kz
!.......south boundary:
            ften(i,k) = wg(i)*ften(i,k) + (1.-wg(i))*usbt(i,k,j)
!.......north boundary:
            ften(ii,k) = wg(i)*ften(ii,k) + (1.-wg(i))*unbt(i,k,j)
          end do
        end do
!
      else if ( jsls.le.ip ) then
        ibeg = 2
        iend = iym1
        if ( jsls.gt.2 ) then
          do i = 2 , jsls - 1
            ii = iy - i + 1
            do k = 1 , kz
!........south boundary:
              ften(i,k) = wg(i)*ften(i,k) + (1.-wg(i))*usbt(i,k,j)
!........north boundary:
              ften(ii,k) = wg(i)*ften(ii,k) + (1.-wg(i))*unbt(i,k,j)
            end do
          end do
          ibeg = jsls
          iend = iy - jsls + 1
        end if
!
        if ( jj.gt.ip ) then
!------west-boundary slice:
          do k = 1 , kz
            do i = ibeg , iend
#ifdef MPP1
              if ( jsls.le.ip ) ften(i,k) = wg(jsls)*ften(i,k)          &
                 & + (1.-wg(jsls))*uwbt(i,k,jew)
#else
              ften(i,k) = wg(jsls)*ften(i,k) + (1.-wg(jsls))            &
                        & *uwbt(i,k,jsls)
#endif
            end do
          end do
        else if ( jj.le.ip ) then
!------east-boundary slice:
          do k = 1 , kz
            do i = ibeg , iend
#ifdef MPP1
              if ( jsls.le.ip ) ften(i,k) = wg(jsls)*ften(i,k)          &
                 & + (1.-wg(jsls))*uebt(i,k,jew)
#else
              ften(i,k) = wg(jsls)*ften(i,k) + (1.-wg(jsls))            &
                        & *uebt(i,k,jsls)
#endif
            end do
          end do
        else
        end if
!
      else
      end if
#endif
      end subroutine sponge_u
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine sponge_v(ip,wg,ften,j)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine applies sponge boundary condition to the        c
!     tendency term - ften.                                           c
!                                                                     c
!     ip   : is the number of slices affected by sponge boundary.     c
!                                                                     c
!     wg   : are the weightings.                                      c
!                                                                     c
!     ften : is the tendency calculated from the model.               c
!                                                                     c
!     vebt, vwbt, vnbt, vsbt : are the large-scale or observed        c
!            tendencies at east, west, north, and south boundaries.   c
!                                                                     c
!     ie = iy, je = jx for dot-point variables.                       c
!     ie = iym1, je = jxm1 for cross-point variables.                 c
!                                                                     c
!     j    : is the j'th slice of the tendency to be adjusted.        c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_dynparam
      use mod_bdycod
      implicit none
!
! Dummy arguments
!
      integer :: ip , j
      real(8) , dimension(iy,kz) :: ften
      real(8) , dimension(ip) :: wg
      intent (in) ip , j , wg
      intent (inout) ften
!
! Local variables
!
      integer :: i , ibeg , iend , ii , jj , jsls , k
#ifdef MPP1
      integer :: jew
#endif
!
#ifdef BAND
!----------------------------------------------------------------------
!
!-----interior j slices:
      do i = 2 , ip
         ii = iy - i + 1
         do k = 1 , kz
!.......south boundary:
            ften(i,k) = wg(i)*ften(i,k) + (1.-wg(i))*vsbt(i,k,j)
!.......north boundary:
            ften(ii,k) = wg(i)*ften(ii,k) + (1.-wg(i))*vnbt(i,k,j)
         end do
      end do
#else
!----------------------------------------------------------------------
!
#ifdef MPP1
      jsls = j + myid*jxp
      jj = jxp1 - jsls
      if ( jj.le.ip ) jsls = jj
      jew = jsls
      if ( jew.gt.jxp ) jew = mod(jsls,jxp)
      if ( jew.eq.0 ) jew = jxp
#else
      jsls = j
      jj = jxp1 - jsls
      if ( jj.le.ip ) jsls = jj
#endif
!
      if ( jsls.gt.ip ) then
!-----interior j slices:
        do i = 2 , ip
          ii = iy - i + 1
          do k = 1 , kz
!.......south boundary:
            ften(i,k) = wg(i)*ften(i,k) + (1.-wg(i))*vsbt(i,k,j)
!.......north boundary:
            ften(ii,k) = wg(i)*ften(ii,k) + (1.-wg(i))*vnbt(i,k,j)
          end do
        end do
!
      else if ( jsls.le.ip ) then
        ibeg = 2
        iend = iym1
        if ( jsls.gt.2 ) then
          do i = 2 , jsls - 1
            ii = iy - i + 1
            do k = 1 , kz
!........south boundary:
              ften(i,k) = wg(i)*ften(i,k) + (1.-wg(i))*vsbt(i,k,j)
!........north boundary:
              ften(ii,k) = wg(i)*ften(ii,k) + (1.-wg(i))*vnbt(i,k,j)
            end do
          end do
          ibeg = jsls
          iend = iy - jsls + 1
        end if
!
        if ( jj.gt.ip ) then
!------west-boundary slice:
          do k = 1 , kz
            do i = ibeg , iend
#ifdef MPP1
              if ( jsls.le.ip ) ften(i,k) = wg(jsls)*ften(i,k)          &
                 & + (1.-wg(jsls))*vwbt(i,k,jew)
#else
              ften(i,k) = wg(jsls)*ften(i,k) + (1.-wg(jsls))            &
                        & *vwbt(i,k,jsls)
#endif
            end do
          end do
        else if ( jj.le.ip ) then
!------east-boundary slice:
          do k = 1 , kz
            do i = ibeg , iend
#ifdef MPP1
              if ( jsls.le.ip ) ften(i,k) = wg(jsls)*ften(i,k)          &
                 & + (1.-wg(jsls))*vebt(i,k,jew)
#else
              ften(i,k) = wg(jsls)*ften(i,k) + (1.-wg(jsls))            &
                        & *vebt(i,k,jsls)
#endif
            end do
          end do
        else
        end if
!
      else
      end if
#endif
      end subroutine sponge_v
