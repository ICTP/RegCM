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
 
      subroutine bdyval(xt,iexec)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine sets the boundary values for p*, p*u, p*v,      c
!     p*t, p*qv, p*qc, and p*qr.                                      c
!                                                                     c
!     ---the boundary values of p*u and p*v are extrapolated from     c
!        the interior points.                                         c
!                                                                     c
!     ---the boundary values of p* and p*t are specified.             c
!                                                                     c
!     ---the boundary values of p*qv, p*qc, and p*qr depend on        c
!        inflow/outflow conditions, if iboudy = 3 or 4.               c
!                                                                     c
!     xt     : is the time in minutes the variables xxa represent.    c
!                                                                     c
!     iexec  : = 1 ; represents this subroutine is called for the     c
!                    first time in this forecast run.                 c
!              > 1 ; represents subsequent calls to this subroutine.  c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use regcm_param
      use param1
      use param2
      use main
      use mainchem
      use bdycod
      use date
      implicit none
!
! Dummy arguments
!
      integer :: iexec
      real(8) :: xt
      intent (in) iexec , xt
!
! Local variables
!
      real(8) :: chix , chix1 , chix2 , dtb , qcx , qcx2 , qvx , qvx1 , &
               & qvx2 , uavg , vavg
      integer :: i , itr , j , k
!
!*********************************************************************
!*****fill up the boundary value for xxb variables from xxa variables:
!     if this subroutine is called for the first time, this part
!     shall be skipped.
!
      if ( iexec.ne.1 ) then
#ifdef MPP1
!
!-----for p*:
!
        do i = 1 , ilx
          if ( myid.eq.0 ) psb(i,1) = psa(i,1)
          if ( myid.eq.nproc-1 ) psb(i,jendx) = psa(i,jendx)
        end do
        do j = jbegin , jendm
          psb(1,j) = psa(1,j)
          psb(ilx,j) = psa(ilx,j)
        end do
!
!-----for p*u and p*v:
!
        do k = 1 , kx
          do i = 1 , ix
            if ( myid.eq.0 ) then
              ub(i,k,1) = ua(i,k,1)/msfd(i,1)
              vb(i,k,1) = va(i,k,1)/msfd(i,1)
            end if
            if ( myid.eq.nproc-1 ) then
              ub(i,k,jendl) = ua(i,k,jendl)/msfd(i,jendl)
              vb(i,k,jendl) = va(i,k,jendl)/msfd(i,jendl)
            end if
          end do
          do j = jbegin , jendx
            ub(1,k,j) = ua(1,k,j)/msfd(1,j)
            ub(ix,k,j) = ua(ix,k,j)/msfd(ix,j)
            vb(1,k,j) = va(1,k,j)/msfd(1,j)
            vb(ix,k,j) = va(ix,k,j)/msfd(ix,j)
          end do
        end do
!
!-----for p*t:
!
        do k = 1 , kx
          do i = 1 , ilx
            if ( myid.eq.0 ) tb(i,k,1) = ta(i,k,1)
            if ( myid.eq.nproc-1 ) tb(i,k,jendx) = ta(i,k,jendx)
          end do
          do j = jbegin , jendm
            tb(1,k,j) = ta(1,k,j)
            tb(ilx,k,j) = ta(ilx,k,j)
          end do
        end do
!
!-----for p*qv:
!
        do k = 1 , kx
          do i = 1 , ilx
            if ( myid.eq.0 ) qvb(i,k,1) = qva(i,k,1)
            if ( myid.eq.nproc-1 ) qvb(i,k,jendx) = qva(i,k,jendx)
          end do
          do j = jbegin , jendm
            qvb(1,k,j) = qva(1,k,j)
            qvb(ilx,k,j) = qva(ilx,k,j)
          end do
        end do
!
!chem2
!
        if ( ichem.eq.1 ) then
!-----for p*chi (tracers)
          do itr = 1 , ntr
            do k = 1 , kx
              do i = 1 , ilx
                if ( myid.eq.0 ) chib(i,k,1,itr) = chia(i,k,1,itr)
                if ( myid.eq.nproc-1 ) chib(i,k,jendx,itr)              &
                   & = chia(i,k,jendx,itr)
              end do
              do j = jbegin , jendm
                chib(1,k,j,itr) = chia(1,k,j,itr)
                chib(ilx,k,j,itr) = chia(ilx,k,j,itr)
              end do
            end do
          end do
        end if
!chem2_
!
!-----for p*qc:
!
        do k = 1 , kx
          do i = 1 , ilx
            if ( myid.eq.0 ) qcb(i,k,1) = qca(i,k,1)
            if ( myid.eq.nproc-1 ) qcb(i,k,jendx) = qca(i,k,jendx)
          end do
          do j = jbegin , jendm
            qcb(1,k,j) = qca(1,k,j)
            qcb(ilx,k,j) = qca(ilx,k,j)
          end do
        end do
!
      end if      !end if(iexec.ne.1) test
!**********************************************************************
!*****compute the boundary values for xxa variables:
!
!-----compute the time interval for boundary tendency:
!
      dtb = xt*60.
      if ( dabs(xt).lt.0.00001 .and. ldatez.gt.idate0 )                 &
         & dtb = ibdyfrq*60.*60.
!
!-----set boundary values for p*:
!-----set boundary conditions for p*u and p*v:
!
      if ( .not.(iexec.eq.1 .and. ifrest) ) then
!
        if ( iboudy.eq.0 ) then
!.....fixed boundary conditions:
          do i = 1 , ilx
            if ( myid.eq.0 ) psa(i,1) = pwb(i,1)
            if ( myid.eq.nproc-1 ) psa(i,jendx) = peb(i,1)
          end do
          do j = jbegin , jendm
            psa(1,j) = pss(1,j)
            psa(ilx,j) = pnb(1,j)
          end do
!
          do k = 1 , kx
            do i = 1 , ix
              if ( myid.eq.0 ) then
                ua(i,k,1) = uwb(i,k,1)
                va(i,k,1) = vwb(i,k,1)
              end if
              if ( myid.eq.nproc-1 ) then
                ua(i,k,jendl) = ueb(i,k,1)
                va(i,k,jendl) = veb(i,k,1)
              end if
            end do
            do j = jbegin , jendx
              ua(1,k,j) = usb(1,k,j)
              ua(ix,k,j) = unb(1,k,j)
              va(1,k,j) = vsb(1,k,j)
              va(ix,k,j) = vnb(1,k,j)
            end do
          end do
        end if
!
!.....time-dependent boundary conditions:
!
        do i = 1 , ilx
          if ( myid.eq.0 ) psa(i,1) = pwb(i,1) + dtb*pwbt(i,1)
          if ( myid.eq.nproc-1 ) psa(i,jendx) = peb(i,1) + dtb*pebt(i,1)
        end do
        do j = jbegin , jendm
          psa(1,j) = pss(1,j) + dtb*psbt(1,j)
          psa(ilx,j) = pnb(1,j) + dtb*pnbt(1,j)
        end do
!
        do k = 1 , kx
          do i = 1 , ix
            if ( myid.eq.0 ) then
              ua(i,k,1) = uwb(i,k,1) + dtb*uwbt(i,k,1)
              va(i,k,1) = vwb(i,k,1) + dtb*vwbt(i,k,1)
            end if
            if ( myid.eq.nproc-1 ) then
              ua(i,k,jendl) = ueb(i,k,1) + dtb*uebt(i,k,1)
              va(i,k,jendl) = veb(i,k,1) + dtb*vebt(i,k,1)
            end if
          end do
          do j = jbegin , jendx
            ua(1,k,j) = usb(1,k,j) + dtb*usbt(1,k,j)
            ua(ix,k,j) = unb(1,k,j) + dtb*unbt(1,k,j)
            va(1,k,j) = vsb(1,k,j) + dtb*vsbt(1,k,j)
            va(ix,k,j) = vnb(1,k,j) + dtb*vnbt(1,k,j)
          end do
        end do
      end if
!
!-----get boundary values of u and v:
!
      call bdyuv(iboudy,dtb)
!
      if ( iexec.eq.1 .and. ifrest ) return
!
!-----set boundary values for p*t:
!-----set boundary values for p*qv:
!
      if ( iboudy.eq.0 ) then
!.....fixed boundary conditions:
        do k = 1 , kx
          do i = 1 , ilx
            if ( myid.eq.0 ) ta(i,k,1) = twb(i,k,1)
            if ( myid.eq.nproc-1 ) ta(i,k,jendx) = teb(i,k,1)
          end do
          do j = jbegin , jendm
            ta(1,k,j) = tsb(1,k,j)
            ta(ilx,k,j) = tnb(1,k,j)
          end do
        end do
        do k = 1 , kx
          do i = 1 , ilx
            if ( myid.eq.0 ) qva(i,k,1) = qwb(i,k,1)
            if ( myid.eq.nproc-1 ) qva(i,k,jendx) = qeb(i,k,1)
          end do
          do j = jbegin , jendm
            qva(1,k,j) = qsb(1,k,j)
            qva(ilx,k,j) = qnb(1,k,j)
          end do
        end do
      end if
!
!.....time-dependent boundary conditions:
!
      do k = 1 , kx
        do i = 1 , ilx
          if ( myid.eq.0 ) ta(i,k,1) = twb(i,k,1) + dtb*twbt(i,k,1)
          if ( myid.eq.nproc-1 ) ta(i,k,jendx) = teb(i,k,1)             &
             & + dtb*tebt(i,k,1)
        end do
        do j = jbegin , jendm
          ta(1,k,j) = tsb(1,k,j) + dtb*tsbt(1,k,j)
          ta(ilx,k,j) = tnb(1,k,j) + dtb*tnbt(1,k,j)
        end do
      end do
      do k = 1 , kx
        do i = 1 , ilx
          if ( myid.eq.0 ) qva(i,k,1) = qwb(i,k,1) + dtb*qwbt(i,k,1)
          if ( myid.eq.nproc-1 ) qva(i,k,jendx) = qeb(i,k,1)            &
             & + dtb*qebt(i,k,1)
        end do
        do j = jbegin , jendm
          qva(1,k,j) = qsb(1,k,j) + dtb*qsbt(1,k,j)
          qva(ilx,k,j) = qnb(1,k,j) + dtb*qnbt(1,k,j)
        end do
      end do
!
      if ( iboudy.eq.3 .or. iboudy.eq.4 ) then
!
!-----determine boundary values depends on inflow/outflow:
!
        do k = 1 , kx
!
!.....west boundary:
!
          if ( myid.eq.0 ) then
            do i = 1 , ilx
              qvx1 = qva(i,k,1)/psa(i,1)
              qvx2 = qva(i,k,2)/psa(i,2)
              uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
              if ( uavg.ge.0. ) then
                qvx = qvx1
              else
                qvx = qvx2
              end if
              qva(i,k,1) = qvx*psa(i,1)
            end do
          end if
!
!.....east boundary:
!
          if ( myid.eq.nproc-1 ) then
            do i = 1 , ilx
              qvx1 = qva(i,k,jendx)/psa(i,jendx)
              qvx2 = qva(i,k,jendm)/psa(i,jendm)
              uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
              if ( uavg.lt.0. ) then
                qvx = qvx1
              else
                qvx = qvx2
              end if
              qva(i,k,jendx) = qvx*psa(i,jendx)
            end do
          end if
!
!.....south boundary:
!
          do j = jbegin , jendm
            qvx1 = qva(1,k,j)/psa(1,j)
            qvx2 = qva(2,k,j)/psa(2,j)
            vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
            if ( vavg.ge.0. ) then
              qvx = qvx1
            else
              qvx = qvx2
            end if
            qva(1,k,j) = qvx*psa(1,j)
          end do
!
!.....north boundary:
!
          do j = jbegin , jendm
            qvx1 = qva(ilx,k,j)/psa(ilx,j)
            qvx2 = qva(ilxm,k,j)/psa(ilxm,j)
            vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
            if ( vavg.lt.0. ) then
              qvx = qvx1
            else
              qvx = qvx2
            end if
            qva(ilx,k,j) = qvx*psa(ilx,j)
          end do
!
        end do
      end if      !end if(iboudy.eq.3.or.4) test
!
!-----set boundary values for p*qc and p*qr:
!     *** note ***
!     for large domain, we assume the boundary tendencies are not
!     available.
!
!
!-----if the boundary values and tendencies are not available,
!     determine boundary values depends on inflow/outflow:
!     inflow  : set it equal to zero.
!     outflow : get from interior point.
!
      do k = 1 , kx
!
!.....west boundary:
!
        if ( myid.eq.0 ) then
          do i = 1 , ilx
            qcx2 = qca(i,k,2)/psa(i,2)
            uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
            if ( uavg.ge.0. ) then
              qcx = 0.
            else
              qcx = qcx2
            end if
            qca(i,k,1) = qcx*psa(i,1)
          end do
        end if
!
!.....east boundary:
!
        if ( myid.eq.nproc-1 ) then
          do i = 1 , ilx
            qcx2 = qca(i,k,jendm)/psa(i,jendm)
            uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
            if ( uavg.lt.0. ) then
              qcx = 0.
            else
              qcx = qcx2
            end if
            qca(i,k,jendx) = qcx*psa(i,jendx)
          end do
        end if
!
!.....south boundary:
!
        do j = jbegin , jendm
          qcx2 = qca(2,k,j)/psa(2,j)
          vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
          if ( vavg.ge.0. ) then
            qcx = 0.
          else
            qcx = qcx2
          end if
          qca(1,k,j) = qcx*psa(1,j)
        end do
!
!.....north boundary:
!
        do j = jbegin , jendm
          qcx2 = qca(ilxm,k,j)/psa(ilxm,j)
          vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
          if ( vavg.lt.0. ) then
            qcx = 0.
          else
            qcx = qcx2
          end if
          qca(ilx,k,j) = qcx*psa(ilx,j)
        end do
!
      end do
 
      if ( ichem.eq.1 ) then
!chem2
 
!----add tracer bc's
!
        do itr = 1 , ntr
          do k = 1 , kx
!
!.....west  boundary:
!
 
            if ( myid.eq.0 ) then
              do i = 1 , ilx
                chix1 = chia(i,k,1,itr)/psa(i,1)
                chix2 = chia(i,k,2,itr)/psa(i,2)
                uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
                if ( uavg.ge.0. ) then
                  chix = chix1
                else
                  chix = chix2
                end if
                chia(i,k,1,itr) = chix*psa(i,1)
              end do
            end if
!
!.....east  boundary:
!
            if ( myid.eq.nproc-1 ) then
              do i = 1 , ilx
                chix1 = chia(i,k,jendx,itr)/psa(i,jendx)
                chix2 = chia(i,k,jendm,itr)/psa(i,jendm)
                uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
                if ( uavg.lt.0. ) then
                  chix = chix1
                else
                  chix = chix2
                end if
                chia(i,k,jendx,itr) = chix*psa(i,jendx)
              end do
            end if
!
!.....south boundary:
!
            do j = jbegin , jendm
              chix1 = chia(1,k,j,itr)/psa(1,j)
              chix2 = chia(2,k,j,itr)/psa(2,j)
              vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
              if ( vavg.ge.0. ) then
                chix = chix1
              else
                chix = chix2
              end if
              chia(1,k,j,itr) = chix*psa(1,j)
            end do
!
!.....north boundary:
!
            do j = jbegin , jendm
              chix1 = chia(ilx,k,j,itr)/psa(ilx,j)
              chix2 = chia(ilxm,k,j,itr)/psa(ilxm,j)
              vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
              if ( vavg.lt.0. ) then
                chix = chix1
              else
                chix = chix2
              end if
              chia(ilx,k,j,itr) = chix*psa(ilx,j)
            end do
          end do
        end do
!chem2_
      end if
#else
!
!-----for p*:
!
        do i = 1 , ilx
          psb(i,1) = psa(i,1)
          psb(i,jlx) = psa(i,jlx)
        end do
        do j = 2 , jlxm
          psb(1,j) = psa(1,j)
          psb(ilx,j) = psa(ilx,j)
        end do
!
!-----for p*u and p*v:
!
        do k = 1 , kx
          do i = 1 , ix
            ub(i,k,1) = ua(i,k,1)/msfd(i,1)
            vb(i,k,1) = va(i,k,1)/msfd(i,1)
            ub(i,k,jx) = ua(i,k,jx)/msfd(i,jx)
            vb(i,k,jx) = va(i,k,jx)/msfd(i,jx)
          end do
          do j = 2 , jlx
            ub(1,k,j) = ua(1,k,j)/msfd(1,j)
            ub(ix,k,j) = ua(ix,k,j)/msfd(ix,j)
            vb(1,k,j) = va(1,k,j)/msfd(1,j)
            vb(ix,k,j) = va(ix,k,j)/msfd(ix,j)
          end do
        end do
!
!-----for p*t:
!
        do k = 1 , kx
          do i = 1 , ilx
            tb(i,k,1) = ta(i,k,1)
            tb(i,k,jlx) = ta(i,k,jlx)
          end do
          do j = 2 , jlxm
            tb(1,k,j) = ta(1,k,j)
            tb(ilx,k,j) = ta(ilx,k,j)
          end do
        end do
!
!-----for p*qv:
!
        do k = 1 , kx
          do i = 1 , ilx
            qvb(i,k,1) = qva(i,k,1)
            qvb(i,k,jlx) = qva(i,k,jlx)
          end do
          do j = 2 , jlxm
            qvb(1,k,j) = qva(1,k,j)
            qvb(ilx,k,j) = qva(ilx,k,j)
          end do
        end do
!
!chem2
!
        if ( ichem.eq.1 ) then
!-----for p*chi (tracers)
          do itr = 1 , ntr
            do k = 1 , kx
              do i = 1 , ilx
                chib(i,k,1,itr) = chia(i,k,1,itr)
                chib(i,k,jlx,itr) = chia(i,k,jlx,itr)
              end do
              do j = 2 , jlxm
                chib(1,k,j,itr) = chia(1,k,j,itr)
                chib(ilx,k,j,itr) = chia(ilx,k,j,itr)
              end do
            end do
          end do
        end if
!chem2_
!
!-----for p*qc:
!
        do k = 1 , kx
          do i = 1 , ilx
            qcb(i,k,1) = qca(i,k,1)
            qcb(i,k,jlx) = qca(i,k,jlx)
          end do
          do j = 2 , jlxm
            qcb(1,k,j) = qca(1,k,j)
            qcb(ilx,k,j) = qca(ilx,k,j)
          end do
        end do
!
      end if      !end if(iexec.ne.1) test
!**********************************************************************
!*****compute the boundary values for xxa variables:
!
!-----compute the time interval for boundary tendency:
!
      dtb = xt*60.
      if ( dabs(xt).lt.0.00001 .and. ldatez.gt.idate0 )                 &
         & dtb = ibdyfrq*60.*60.
!
!-----set boundary values for p*:
!-----set boundary conditions for p*u and p*v:
!
      if ( .not.(iexec.eq.1 .and. ifrest) ) then
!
        if ( iboudy.eq.0 ) then
!.....fixed boundary conditions:
          do i = 1 , ilx
            psa(i,1) = pwb(i,1)
            psa(i,jlx) = peb(i,1)
          end do
          do j = 2 , jlxm
            psa(1,j) = pss(1,j)
            psa(ilx,j) = pnb(1,j)
          end do
!
          do k = 1 , kx
            do i = 1 , ix
              ua(i,k,1) = uwb(i,k,1)
              va(i,k,1) = vwb(i,k,1)
              ua(i,k,jx) = ueb(i,k,1)
              va(i,k,jx) = veb(i,k,1)
            end do
            do j = 2 , jlx
              ua(1,k,j) = usb(1,k,j)
              ua(ix,k,j) = unb(1,k,j)
              va(1,k,j) = vsb(1,k,j)
              va(ix,k,j) = vnb(1,k,j)
            end do
          end do
        end if
!
!.....time-dependent boundary conditions:
!
        do i = 1 , ilx
          psa(i,1) = pwb(i,1) + dtb*pwbt(i,1)
          psa(i,jlx) = peb(i,1) + dtb*pebt(i,1)
        end do
        do j = 2 , jlxm
          psa(1,j) = pss(1,j) + dtb*psbt(1,j)
          psa(ilx,j) = pnb(1,j) + dtb*pnbt(1,j)
        end do
!
        do k = 1 , kx
          do i = 1 , ix
            ua(i,k,1) = uwb(i,k,1) + dtb*uwbt(i,k,1)
            va(i,k,1) = vwb(i,k,1) + dtb*vwbt(i,k,1)
            ua(i,k,jx) = ueb(i,k,1) + dtb*uebt(i,k,1)
            va(i,k,jx) = veb(i,k,1) + dtb*vebt(i,k,1)
          end do
          do j = 2 , jlx
            ua(1,k,j) = usb(1,k,j) + dtb*usbt(1,k,j)
            ua(ix,k,j) = unb(1,k,j) + dtb*unbt(1,k,j)
            va(1,k,j) = vsb(1,k,j) + dtb*vsbt(1,k,j)
            va(ix,k,j) = vnb(1,k,j) + dtb*vnbt(1,k,j)
          end do
        end do
      end if
!
!-----get boundary values of u and v:
!
      call bdyuv(iboudy,dtb)
!
      if ( iexec.eq.1 .and. ifrest ) return
!
!-----set boundary values for p*t:
!-----set boundary values for p*qv:
!
      if ( iboudy.eq.0 ) then
!.....fixed boundary conditions:
        do k = 1 , kx
          do i = 1 , ilx
            ta(i,k,1) = twb(i,k,1)
            ta(i,k,jlx) = teb(i,k,1)
          end do
          do j = 2 , jlxm
            ta(1,k,j) = tsb(1,k,j)
            ta(ilx,k,j) = tnb(1,k,j)
          end do
        end do
        do k = 1 , kx
          do i = 1 , ilx
            qva(i,k,1) = qwb(i,k,1)
            qva(i,k,jlx) = qeb(i,k,1)
          end do
          do j = 2 , jlxm
            qva(1,k,j) = qsb(1,k,j)
            qva(ilx,k,j) = qnb(1,k,j)
          end do
        end do
      end if
!
!.....time-dependent boundary conditions:
!
      do k = 1 , kx
        do i = 1 , ilx
          ta(i,k,1) = twb(i,k,1) + dtb*twbt(i,k,1)
          ta(i,k,jlx) = teb(i,k,1) + dtb*tebt(i,k,1)
        end do
        do j = 2 , jlxm
          ta(1,k,j) = tsb(1,k,j) + dtb*tsbt(1,k,j)
          ta(ilx,k,j) = tnb(1,k,j) + dtb*tnbt(1,k,j)
        end do
      end do
      do k = 1 , kx
        do i = 1 , ilx
          qva(i,k,1) = qwb(i,k,1) + dtb*qwbt(i,k,1)
          qva(i,k,jlx) = qeb(i,k,1) + dtb*qebt(i,k,1)
        end do
        do j = 2 , jlxm
          qva(1,k,j) = qsb(1,k,j) + dtb*qsbt(1,k,j)
          qva(ilx,k,j) = qnb(1,k,j) + dtb*qnbt(1,k,j)
        end do
      end do
!
      if ( iboudy.eq.3 .or. iboudy.eq.4 ) then
!
!-----determine boundary values depends on inflow/outflow:
!
        do k = 1 , kx
!
!.....west boundary:
!
          do i = 1 , ilx
            qvx1 = qva(i,k,1)/psa(i,1)
            qvx2 = qva(i,k,2)/psa(i,2)
            uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
            if ( uavg.ge.0. ) then
              qvx = qvx1
            else
              qvx = qvx2
            end if
            qva(i,k,1) = qvx*psa(i,1)
          end do
!
!.....east boundary:
!
          do i = 1 , ilx
            qvx1 = qva(i,k,jlx)/psa(i,jlx)
            qvx2 = qva(i,k,jlxm)/psa(i,jlxm)
            uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
            if ( uavg.lt.0. ) then
              qvx = qvx1
            else
              qvx = qvx2
            end if
            qva(i,k,jlx) = qvx*psa(i,jlx)
          end do
!
!.....south boundary:
!
          do j = 2 , jlxm
            qvx1 = qva(1,k,j)/psa(1,j)
            qvx2 = qva(2,k,j)/psa(2,j)
            vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
            if ( vavg.ge.0. ) then
              qvx = qvx1
            else
              qvx = qvx2
            end if
            qva(1,k,j) = qvx*psa(1,j)
          end do
!
!.....north boundary:
!
          do j = 2 , jlxm
            qvx1 = qva(ilx,k,j)/psa(ilx,j)
            qvx2 = qva(ilxm,k,j)/psa(ilxm,j)
            vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
            if ( vavg.lt.0. ) then
              qvx = qvx1
            else
              qvx = qvx2
            end if
            qva(ilx,k,j) = qvx*psa(ilx,j)
          end do
!
        end do
      end if      !end if(iboudy.eq.3.or.4) test
!
!-----set boundary values for p*qc and p*qr:
!     *** note ***
!     for large domain, we assume the boundary tendencies are not
!     available.
!
!
!-----if the boundary values and tendencies are not available,
!     determine boundary values depends on inflow/outflow:
!     inflow  : set it equal to zero.
!     outflow : get from interior point.
!
      do k = 1 , kx
!
!.....west boundary:
!
        do i = 1 , ilx
          qcx2 = qca(i,k,2)/psa(i,2)
          uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
          if ( uavg.ge.0. ) then
            qcx = 0.
          else
            qcx = qcx2
          end if
          qca(i,k,1) = qcx*psa(i,1)
        end do
!
!.....east boundary:
!
        do i = 1 , ilx
          qcx2 = qca(i,k,jlxm)/psa(i,jlxm)
          uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
          if ( uavg.lt.0. ) then
            qcx = 0.
          else
            qcx = qcx2
          end if
          qca(i,k,jlx) = qcx*psa(i,jlx)
        end do
!
!.....south boundary:
!
        do j = 2 , jlxm
          qcx2 = qca(2,k,j)/psa(2,j)
          vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
          if ( vavg.ge.0. ) then
            qcx = 0.
          else
            qcx = qcx2
          end if
          qca(1,k,j) = qcx*psa(1,j)
        end do
!
!.....north boundary:
!
        do j = 2 , jlxm
          qcx2 = qca(ilxm,k,j)/psa(ilxm,j)
          vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
          if ( vavg.lt.0. ) then
            qcx = 0.
          else
            qcx = qcx2
          end if
          qca(ilx,k,j) = qcx*psa(ilx,j)
        end do
!
      end do
 
      if ( ichem.eq.1 ) then
!chem2
 
!----add tracer bc's
!
        do itr = 1 , ntr
          do k = 1 , kx
!
!.....west  boundary:
!
 
            do i = 1 , ilx
              chix1 = chia(i,k,1,itr)/psa(i,1)
              chix2 = chia(i,k,2,itr)/psa(i,2)
              uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
              if ( uavg.ge.0. ) then
                chix = chix1
              else
                chix = chix2
              end if
              chia(i,k,1,itr) = chix*psa(i,1)
            end do
!
!.....east  boundary:
!
            do i = 1 , ilx
              chix1 = chia(i,k,jlx,itr)/psa(i,jlx)
              chix2 = chia(i,k,jlxm,itr)/psa(i,jlxm)
              uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
              if ( uavg.lt.0. ) then
                chix = chix1
              else
                chix = chix2
              end if
              chia(i,k,jlx,itr) = chix*psa(i,jlx)
            end do
!
!.....south boundary:
!
            do j = 2 , jlxm
              chix1 = chia(1,k,j,itr)/psa(1,j)
              chix2 = chia(2,k,j,itr)/psa(2,j)
              vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
              if ( vavg.ge.0. ) then
                chix = chix1
              else
                chix = chix2
              end if
              chia(1,k,j,itr) = chix*psa(1,j)
            end do
!
!.....north boundary:
!
            do j = 2 , jlxm
              chix1 = chia(ilx,k,j,itr)/psa(ilx,j)
              chix2 = chia(ilxm,k,j,itr)/psa(ilxm,j)
              vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
              if ( vavg.lt.0. ) then
                chix = chix1
              else
                chix = chix2
              end if
              chia(ilx,k,j,itr) = chix*psa(ilx,j)
            end do
          end do
        end do
!chem2_
      end if
#endif
!
      end subroutine bdyval
