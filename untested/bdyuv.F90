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
 
      subroutine bdyuv(ib,dtb)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine sets the boundary values of u and v according   c
!     to the boundary conditions specified.                           c
!                                                                     c
!     ua, va, and psa : variables needed                              c
!                                                                     c
!     ib = 0 : fixed                                                  c
!        = 1 : relaxation, linear technique                           c
!        = 2 : time dependent                                         c
!        = 3 : time dependent and inflow/outflow dependent            c
!        = 4 : sponge                                                 c
!        = 5 : relaxation, exponential technique                      c
!                                                                     c
!     dtb    : elapsed time from the initial boundary values.         c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_regcm_param
      use mod_bdycod
      use mod_cvaria
      use mod_main
#ifdef MPP1
      use mod_mppio
      use mpi
#endif
      implicit none
!
! Dummy arguments
!
      real(8) :: dtb
      integer :: ib
      intent (in) dtb , ib
!
! Local variables
!
      integer :: i , j , k
#ifdef MPP1
      integer :: ierr
      integer , dimension(mpi_status_size) :: status
#endif
!
#ifdef MPP1
!----------------------------------------------------------------------
!-----compute the p* at dot points:
!
!=======================================================================
      call mpi_sendrecv(psa(1,jxp),ix,mpi_double_precision,ieast,1,     &
                      & psa(1,0),ix,mpi_double_precision,iwest,1,       &
                      & mpi_comm_world,status,ierr)
#endif
!=======================================================================
!
!-----interior points:
!
#ifdef MPP1
      do j = jbegin , jendx
        do i = 2 , ixm1
          pdota(i,j) = 0.25*(psa(i,j)+psa(i-1,j)+psa(i,j-1)+psa(i-1,j-1)&
                     & )
        end do
      end do
#else
      do j = 2 , jxm1
        do i = 2 , ixm1
          pdota(i,j) = 0.25*(psa(i,j)+psa(i-1,j)+psa(i,j-1)+psa(i-1,j-1)&
                     & )
        end do
      end do
#endif
!
!-----east and west boundaries:
!
      do i = 2 , ixm1
#ifdef MPP1
        if ( myid.eq.0 ) pdota(i,1) = 0.5*(psa(i,1)+psa(i-1,1))
        if ( myid.eq.nproc-1 ) pdota(i,jendl)                           &
           & = 0.5*(psa(i,jendx)+psa(i-1,jendx))
#else
        pdota(i,1) = 0.5*(psa(i,1)+psa(i-1,1))
        pdota(i,jx) = 0.5*(psa(i,jxm1)+psa(i-1,jxm1))
#endif
      end do
!
!-----north and south boundaries:
!
#ifdef MPP1
      do j = jbegin , jendx
        pdota(1,j) = 0.5*(psa(1,j)+psa(1,j-1))
        pdota(ix,j) = 0.5*(psa(ixm1,j)+psa(ixm1,j-1))
      end do
#else
      do j = 2 , jxm1
        pdota(1,j) = 0.5*(psa(1,j)+psa(1,j-1))
        pdota(ix,j) = 0.5*(psa(ixm1,j)+psa(ixm1,j-1))
      end do
#endif
!
!-----corner points:
!
#ifdef MPP1
      if ( myid.eq.0 ) then
        pdota(1,1) = psa(1,1)
        pdota(ix,1) = psa(ixm1,1)
      end if
      if ( myid.eq.nproc-1 ) then
        pdota(1,jendl) = psa(1,jendx)
        pdota(ix,jendl) = psa(ixm1,jendx)
      end if
#else
      pdota(1,1) = psa(1,1)
      pdota(ix,1) = psa(ixm1,1)
      pdota(1,jx) = psa(1,jxm1)
      pdota(ix,jx) = psa(ixm1,jxm1)
#endif
!=======================================================================
!
!-----interior silces:
!
      do k = 1 , kx
!
!.....for j = 2 and j = jlx :
!
        do i = 2 , ixm1
#ifdef MPP1
          if ( myid.eq.0 ) then
            uj2(i,k) = ua(i,k,2)/pdota(i,2)
            vj2(i,k) = va(i,k,2)/pdota(i,2)
          end if
          if ( myid.eq.nproc-1 ) then
            ujlx(i,k) = ua(i,k,jendx)/pdota(i,jendx)
            vjlx(i,k) = va(i,k,jendx)/pdota(i,jendx)
          end if
#else
          uj2(i,k) = ua(i,k,2)/pdota(i,2)
          vj2(i,k) = va(i,k,2)/pdota(i,2)
          ujlx(i,k) = ua(i,k,jxm1)/pdota(i,jxm1)
          vjlx(i,k) = va(i,k,jxm1)/pdota(i,jxm1)
#endif
        end do
!
!.....for i = 2 and i = ixm1 :
!
#ifdef MPP1
        do j = jbegin , jendx
          ui2(k,j) = ua(2,k,j)/pdota(2,j)
          vi2(k,j) = va(2,k,j)/pdota(2,j)
          uilx(k,j) = ua(ixm1,k,j)/pdota(ixm1,j)
          vilx(k,j) = va(ixm1,k,j)/pdota(ixm1,j)
        end do
#else
        do j = 2 , jxm1
          ui2(k,j) = ua(2,k,j)/pdota(2,j)
          vi2(k,j) = va(2,k,j)/pdota(2,j)
          uilx(k,j) = ua(ixm1,k,j)/pdota(ixm1,j)
          vilx(k,j) = va(ixm1,k,j)/pdota(ixm1,j)
        end do
#endif
!
      end do
!
!----------------------------------------------------------------------
!-----boundary silces:
!
      if ( ib.eq.0 ) then
!
!-----fixed boundary conditions:
!
        do k = 1 , kx
!
!.....west (j = 1) and east (j = jx) boundaries:
!
          do i = 1 , ix
#ifdef MPP1
            if ( myid.eq.0 ) then
              uj1(i,k) = uwb(i,k,1)/pdota(i,1)
              vj1(i,k) = vwb(i,k,1)/pdota(i,1)
            end if
            if ( myid.eq.nproc-1 ) then
              ujl(i,k) = ueb(i,k,1)/pdota(i,jendl)
              vjl(i,k) = veb(i,k,1)/pdota(i,jendl)
            end if
#else
            uj1(i,k) = uwb(i,k,1)/pdota(i,1)
            vj1(i,k) = vwb(i,k,1)/pdota(i,1)
            ujl(i,k) = ueb(i,k,1)/pdota(i,jx)
            vjl(i,k) = veb(i,k,1)/pdota(i,jx)
#endif
          end do
!
!.....south (i = 1) and north (i = ix) boundaries:
!
#ifdef MPP1
          do j = 1 , jendl
            ui1(k,j) = usb(1,k,j)/pdota(1,j)
            vi1(k,j) = vsb(1,k,j)/pdota(1,j)
            uil(k,j) = unb(1,k,j)/pdota(ix,j)
            vil(k,j) = vnb(1,k,j)/pdota(ix,j)
          end do
#else
          do j = 1 , jx
            ui1(k,j) = usb(1,k,j)/pdota(1,j)
            vi1(k,j) = vsb(1,k,j)/pdota(1,j)
            uil(k,j) = unb(1,k,j)/pdota(ix,j)
            vil(k,j) = vnb(1,k,j)/pdota(ix,j)
          end do
#endif
        end do
        go to 100
!
      end if
!
!-----time-dependent boundary conditions:
!
      do k = 1 , kx
!
!.....west (j = 1) and east (j = jx) boundaries:
!
        do i = 1 , ix
#ifdef MPP1
          if ( myid.eq.0 ) then
            uj1(i,k) = (uwb(i,k,1)+dtb*uwbt(i,k,1))/pdota(i,1)
            vj1(i,k) = (vwb(i,k,1)+dtb*vwbt(i,k,1))/pdota(i,1)
          end if
          if ( myid.eq.nproc-1 ) then
            ujl(i,k) = (ueb(i,k,1)+dtb*uebt(i,k,1))/pdota(i,jendl)
            vjl(i,k) = (veb(i,k,1)+dtb*vebt(i,k,1))/pdota(i,jendl)
          end if
#else
          uj1(i,k) = (uwb(i,k,1)+dtb*uwbt(i,k,1))/pdota(i,1)
          vj1(i,k) = (vwb(i,k,1)+dtb*vwbt(i,k,1))/pdota(i,1)
          ujl(i,k) = (ueb(i,k,1)+dtb*uebt(i,k,1))/pdota(i,jx)
          vjl(i,k) = (veb(i,k,1)+dtb*vebt(i,k,1))/pdota(i,jx)
#endif
        end do
!
!.....south (i = 1) and north (i = ix) boundaries:
!
#ifdef MPP1
        do j = 1 , jendl
          ui1(k,j) = (usb(1,k,j)+dtb*usbt(1,k,j))/pdota(1,j)
          vi1(k,j) = (vsb(1,k,j)+dtb*vsbt(1,k,j))/pdota(1,j)
          uil(k,j) = (unb(1,k,j)+dtb*unbt(1,k,j))/pdota(ix,j)
          vil(k,j) = (vnb(1,k,j)+dtb*vnbt(1,k,j))/pdota(ix,j)
        end do
#else
        do j = 1 , jx
          ui1(k,j) = (usb(1,k,j)+dtb*usbt(1,k,j))/pdota(1,j)
          vi1(k,j) = (vsb(1,k,j)+dtb*vsbt(1,k,j))/pdota(1,j)
          uil(k,j) = (unb(1,k,j)+dtb*unbt(1,k,j))/pdota(ix,j)
          vil(k,j) = (vnb(1,k,j)+dtb*vnbt(1,k,j))/pdota(ix,j)
        end do
#endif
!
      end do
!
!-----fill up the interior silces:
!
 100  continue

#ifdef MPP1
      do k = 1 , kx
        if ( myid.eq.0 ) then
          uj2(1,k) = ui1(k,2)
          uj2(ix,k) = uil(k,2)
          ui2(k,1) = uj1(2,k)
          uilx(k,1) = uj1(ixm1,k)
          vj2(1,k) = vi1(k,2)
          vj2(ix,k) = vil(k,2)
          vi2(k,1) = vj1(2,k)
          vilx(k,1) = vj1(ixm1,k)
        end if
        if ( myid.eq.nproc-1 ) then
          ujlx(1,k) = ui1(k,jendx)
          ujlx(ix,k) = uil(k,jendx)
          ui2(k,jendl) = ujl(2,k)
          uilx(k,jendl) = ujl(ixm1,k)
          vjlx(1,k) = vi1(k,jendx)
          vjlx(ix,k) = vil(k,jendx)
          vi2(k,jendl) = vjl(2,k)
          vilx(k,jendl) = vjl(ixm1,k)
        end if
      end do
      if ( myid.ne.nproc-1 ) then
        do k = 1 , kx
          var1snd(k,1) = ui1(k,jxp)
          var1snd(k,2) = vi1(k,jxp)
          var1snd(k,3) = ui2(k,jxp)
          var1snd(k,4) = vi2(k,jxp)
          var1snd(k,5) = uilx(k,jxp)
          var1snd(k,6) = vilx(k,jxp)
          var1snd(k,7) = uil(k,jxp)
          var1snd(k,8) = vil(k,jxp)
        end do
      end if
      call mpi_sendrecv(var1snd(1,1),kx*8,mpi_double_precision,ieast,1, &
                      & var1rcv(1,1),kx*8,mpi_double_precision,iwest,1, &
                      & mpi_comm_world,status,ierr)
      if ( myid.ne.0 ) then
        do k = 1 , kx
          ui1(k,0) = var1rcv(k,1)
          vi1(k,0) = var1rcv(k,2)
          ui2(k,0) = var1rcv(k,3)
          vi2(k,0) = var1rcv(k,4)
          uilx(k,0) = var1rcv(k,5)
          vilx(k,0) = var1rcv(k,6)
          uil(k,0) = var1rcv(k,7)
          vil(k,0) = var1rcv(k,8)
        end do
      end if
!
      if ( myid.ne.0 ) then
        do k = 1 , kx
          var1snd(k,1) = ui1(k,1)
          var1snd(k,2) = vi1(k,1)
          var1snd(k,3) = ui2(k,1)
          var1snd(k,4) = vi2(k,1)
          var1snd(k,5) = uilx(k,1)
          var1snd(k,6) = vilx(k,1)
          var1snd(k,7) = uil(k,1)
          var1snd(k,8) = vil(k,1)
        end do
      end if
      call mpi_sendrecv(var1snd(1,1),kx*8,mpi_double_precision,iwest,2, &
                      & var1rcv(1,1),kx*8,mpi_double_precision,ieast,2, &
                      & mpi_comm_world,status,ierr)
      if ( myid.ne.nproc-1 ) then
        do k = 1 , kx
          ui1(k,jxp+1) = var1rcv(k,1)
          vi1(k,jxp+1) = var1rcv(k,2)
          ui2(k,jxp+1) = var1rcv(k,3)
          vi2(k,jxp+1) = var1rcv(k,4)
          uilx(k,jxp+1) = var1rcv(k,5)
          vilx(k,jxp+1) = var1rcv(k,6)
          uil(k,jxp+1) = var1rcv(k,7)
          vil(k,jxp+1) = var1rcv(k,8)
        end do
      end if
#else
      do k = 1 , kx
        uj2(1,k) = ui1(k,2)
        uj2(ix,k) = uil(k,2)
        ui2(k,1) = uj1(2,k)
        uilx(k,1) = uj1(ixm1,k)
        vj2(1,k) = vi1(k,2)
        vj2(ix,k) = vil(k,2)
        vi2(k,1) = vj1(2,k)
        vilx(k,1) = vj1(ixm1,k)
        ujlx(1,k) = ui1(k,jxm1)
        ujlx(ix,k) = uil(k,jxm1)
        ui2(k,jx) = ujl(2,k)
        uilx(k,jx) = ujl(ixm1,k)
        vjlx(1,k) = vi1(k,jxm1)
        vjlx(ix,k) = vil(k,jxm1)
        vi2(k,jx) = vjl(2,k)
        vilx(k,jx) = vjl(ixm1,k)
      end do
#endif
!
      end subroutine bdyuv
