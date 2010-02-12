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
 
      subroutine tend(iexec)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the tendencies of the prognostic       c
!     variables p*, u, v, and t.                                      c
!                                                                     c
!     p*u, p*v, p*t ,p*qv, and p*qc stored in main common block.      c
!                                                                     c
!     all the two-dimension arrays stored in main common block.       c
!                                                                     c
!     east/west boundary conditions stored in common block /bdycod/ . c
!                                                                     c
!     north/south boundary conditions stored in common block          c
!              /bdycod/.                                              c
!                                                                     c
!     all the integers stored in common block /param1/.               c
!                                                                     c
!     all the constants stored in common block /param1/.              c
!                                                                     c
!     iexec  : = 1 ; represents this subroutine is called for the     c
!                    first time in this forecast run.                 c
!              > 1 ; represents subsequent calls to this subroutine.  c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      use mod_regcm_param
      use mod_param1
      use mod_param2
      use mod_param3
      use mod_main
      use mod_mainchem
      use mod_bdycod
      use mod_cvaria
      use mod_pmoist
      use mod_rad
      use mod_bats
      use mod_trachem
      use mod_date
      use mod_message
      use mod_aerosol , only : aermm
#ifdef MPP1
      use mod_slice
      use mpi
#endif
      implicit none
!
! Dummy arguments
!
      integer :: iexec
      intent (inout) iexec
!
! Local variables
!
      real(8) :: alam , cell , chias , chibs , dto2 , dudx , dudy ,     &
               & dvdx , dvdy , p00pg , pgfaa1 , psabar , psasum ,       &
               & pt2bar , pt2tot , ptnbar , ptntot , qcas , qcbs ,      &
               & qvas , qvbs , rovcpm , rtbar , sigpsa , t00pg , tv ,   &
               & tv1 , tv2 , tv3 , tv4 , tva , tvavg , tvb , tvc ,      &
               & xday , xmsf , xtm1
      real(8) , dimension(ix,kx) :: divl
      integer :: i , icons , iptn , itr , j , k , lev , n
#ifdef MPP1
      integer :: ierr , icons_mpi , numrec
      real(8) , dimension(ix,kx,jxp) :: ttld
      real(8) , dimension(ix,kx*16+4) :: bdyewrcv , bdyewsnd
      real(8) , dimension(nspgx,kx*16+4) :: bdynsrcv , bdynssnd
      real(8) , dimension(ix,4,jxp) :: ps4
      real(8) , dimension(ix,4,jx) :: ps_4 
      real(8) , dimension(ix,kx*(ntr+5)*2) :: var2rcv , var2snd
      real(8) , dimension(ix,kx*(ntr+11)+1) :: tvar1rcv , tvar1snd
#else
      real(8) , dimension(ix,kx,jx) :: ttld
#endif
!
!----for pressure gradient force calculations
!
      t00pg = 287.
      p00pg = 101.325
      alam = 6.5E-3
      pgfaa1 = alam*r/g
!
!----------------------------------------------------------------------
!-----fill up the boundary slices:
!
!     if (iexec .eq. 1) then
      if ( .not.ifrest .and. iexec.eq.1 ) then
        call bdyval(xtime,iexec)
        iexec = 2
      else
        iexec = 2
      end if
!
!----------------------------------------------------------------------
!*****for large domain, subroutine tend just needed to go through once.
!
!-----multiply ua and va by inverse of mapscale factor at dot point:
!
#ifdef MPP1
      do j = 1 , jendl
#else
      do j = 1 , jx
#endif
        do k = 1 , kx
          do i = 1 , ix
            ua(i,k,j) = ua(i,k,j)*msfd(i,j)
            va(i,k,j) = va(i,k,j)*msfd(i,j)
          end do
        end do
      end do
#ifdef MPP1
      call mpi_sendrecv(psa(1,jxp),ix,mpi_double_precision,ieast,1,     &
                      & psa(1,0),ix,mpi_double_precision,iwest,1,       &
                      & mpi_comm_world,mpi_status_ignore,ierr)
      call mpi_sendrecv(psa(1,1),ix,mpi_double_precision,iwest,2,       &
                      & psa(1,jxp+1),ix,mpi_double_precision,ieast,2,   &
                      & mpi_comm_world,mpi_status_ignore,ierr)
#endif
!
!-----decouple u, v, t, qv, and qc
!
#ifdef MPP1
      do j = 1 , jendl
        if ( myid.eq.0 .and. j.eq.1 ) then
#else
      do j = 1 , jx
        if ( j.eq.1 ) then
#endif
!-----------lateral slices:
!-----------west boundary:
          do k = 1 , kx
            do i = 1 , ix
              u(i,k,j) = uj1(i,k)
              v(i,k,j) = vj1(i,k)
            end do
          end do
          if ( iboudy.eq.3 .or. iboudy.eq.4 ) then
!..............inflow/outflow dependence:
            do k = 1 , kx
              do i = 1 , ix
                if ( u(i,k,j).lt.0. ) then
                  v(i,k,j) = vj2(i,k)
                  u(i,k,j) = uj2(i,k)
                end if
              end do
            end do
          end if
#ifdef MPP1
        else if ( myid.eq.0 .and. j.eq.2 ) then
#else
        else if ( j.eq.2 ) then
#endif
          do k = 1 , kx
            do i = 1 , ix
              u(i,k,j) = uj2(i,k)
              v(i,k,j) = vj2(i,k)
            end do
          end do
          if ( iboudy.eq.3 .or. iboudy.eq.4 ) then
!..............inflow/outflow dependence:
            do k = 1 , kx
!.................south boundary:
              if ( v(1,k,j).lt.0. ) then
                v(1,k,j) = v(2,k,j)
                u(1,k,j) = u(2,k,j)
              end if
!.................north boundary:
              if ( v(ix,k,j).ge.0. ) then
                v(ix,k,j) = v(ixm1,k,j)
                u(ix,k,j) = u(ixm1,k,j)
              end if
            end do
          end if
#ifdef MPP1
        else if ( myid.eq.nproc-1 .and. j.eq.jendl-1 ) then
#else
        else if ( j.eq.jxm1 ) then
#endif
          do k = 1 , kx
            do i = 1 , ix
              u(i,k,j) = ujlx(i,k)
              v(i,k,j) = vjlx(i,k)
            end do
          end do
          if ( iboudy.eq.3 .or. iboudy.eq.4 ) then
!..............inflow/outflow dependence:
            do k = 1 , kx
!.................south boundary:
              if ( v(1,k,j).lt.0. ) then
                v(1,k,j) = v(2,k,j)
                u(1,k,j) = u(2,k,j)
              end if
!.................north boundary:
              if ( v(ix,k,j).ge.0. ) then
                v(ix,k,j) = v(ixm1,k,j)
                u(ix,k,j) = u(ixm1,k,j)
              end if
            end do
          end if
#ifdef MPP1
        else if ( myid.eq.nproc-1 .and. j.eq.jendl ) then
#else
        else if ( j.eq.jx ) then
#endif
!-----------east boundary:
!
!...........no inflow/outflow dependence:
!
          do k = 1 , kx
            do i = 1 , ix
              u(i,k,j) = ujl(i,k)
              v(i,k,j) = vjl(i,k)
            end do
          end do
          if ( iboudy.eq.3 .or. iboudy.eq.4 ) then
!..............inflow/outflow dependence:
            do k = 1 , kx
              do i = 1 , ix
                if ( u(i,k,j).ge.0. ) then
                  v(i,k,j) = vjlx(i,k)
                  u(i,k,j) = ujlx(i,k)
                end if
              end do
            end do
          end if
        else
!
!-----interior slice:
!-----interior points:
!
          do k = 1 , kx
            do i = 3 , ixm2
              psabar = 0.25*(psa(i,j)+psa(i,j-1)+psa(i-1,j)+psa(i-1,j-1)&
                     & )
              xmsf = msfd(i,j)
              u(i,k,j) = ua(i,k,j)/(psabar*xmsf)
              v(i,k,j) = va(i,k,j)/(psabar*xmsf)
            end do
          end do
!
!-----------north/south boundary points:
!...........no inflow/outflow dependence:
!
          do k = 1 , kx
!..............for i=2 and i=ixm1:
            u(2,k,j) = ui2(k,j)
            u(ixm1,k,j) = uilx(k,j)
            v(2,k,j) = vi2(k,j)
            v(ixm1,k,j) = vilx(k,j)
!..............for i=1 and i=ix:
            u(1,k,j) = ui1(k,j)
            u(ix,k,j) = uil(k,j)
            v(1,k,j) = vi1(k,j)
            v(ix,k,j) = vil(k,j)
          end do
          if ( iboudy.eq.3 .or. iboudy.eq.4 ) then
!..............inflow/outflow dependence:
            do k = 1 , kx
!.................south boundary:
              if ( v(1,k,j).lt.0. ) then
                v(1,k,j) = v(2,k,j)
                u(1,k,j) = u(2,k,j)
              end if
!.................north boundary:
              if ( v(ix,k,j).ge.0. ) then
                v(ix,k,j) = v(ixm1,k,j)
                u(ix,k,j) = u(ixm1,k,j)
              end if
            end do
          end if
        end if
      end do
!
#ifdef MPP1
      do j = 1 , jendx
#else
      do j = 1 , jxm1
#endif
        do k = 1 , kx
          do i = 1 , ixm1
            t(i,k,j) = ta(i,k,j)/psa(i,j)
            qv(i,k,j) = qva(i,k,j)/psa(i,j)
            qc(i,k,j) = qca(i,k,j)/psa(i,j)
          end do
        end do
      end do
!chem2
      if ( ichem.eq.1 ) then
!
!-----call special tracer decoupling routine for multiple (ntr) species
!
        do n = 1 , ntr
#ifdef MPP1
          do j = 1 , jendx
#else
          do j = 1 , jxm1
#endif
            do k = 1 , kx
              do i = 1 , ixm1
                chi(i,k,j,n) = chia(i,k,j,n)/psa(i,j)
              end do
            end do
          end do
        end do
      end if
!chem2_
!
!=======================================================================
#ifdef MPP1
      if ( myid.ne.nproc-1 ) then
        do i = 1 , ix
          tvar1snd(i,1) = psb(i,jxp)
        end do
        do k = 1 , kx
          do i = 1 , ix
            tvar1snd(i,1+k) = tb(i,k,jxp)
            tvar1snd(i,1+kx+k) = qvb(i,k,jxp)
            tvar1snd(i,1+kx*2+k) = ub(i,k,jxp)
            tvar1snd(i,1+kx*3+k) = vb(i,k,jxp)
            tvar1snd(i,1+kx*4+k) = u(i,k,jxp)
            tvar1snd(i,1+kx*5+k) = v(i,k,jxp)
            tvar1snd(i,1+kx*6+k) = t(i,k,jxp)
            tvar1snd(i,1+kx*7+k) = qv(i,k,jxp)
            tvar1snd(i,1+kx*8+k) = qc(i,k,jxp)
            tvar1snd(i,1+kx*9+k) = ua(i,k,jxp)
            tvar1snd(i,1+kx*10+k) = va(i,k,jxp)
          end do
        end do
        if ( ichem.eq.1 ) then
          do n = 1 , ntr
            do k = 1 , kx
              do i = 1 , ix
                tvar1snd(i,kx*11+1+(n-1)*kx+k) = chi(i,k,jxp,n)
              end do
            end do
          end do
        end if
      end if
      numrec = kx*11 + 1
      if ( ichem.eq.1 ) numrec = kx*(ntr+11) + 1
      call mpi_sendrecv(tvar1snd(1,1),ix*numrec,mpi_double_precision,   &
                      & ieast,1,tvar1rcv(1,1),ix*numrec,                &
                      & mpi_double_precision,iwest,1,mpi_comm_world,    &
                      & mpi_status_ignore,ierr)
      if ( myid.ne.0 ) then
        do i = 1 , ix
          psb(i,0) = tvar1rcv(i,1)
        end do
        do k = 1 , kx
          do i = 1 , ix
            tb(i,k,0) = tvar1rcv(i,1+k)
            qvb(i,k,0) = tvar1rcv(i,1+kx+k)
            ub(i,k,0) = tvar1rcv(i,1+kx*2+k)
            vb(i,k,0) = tvar1rcv(i,1+kx*3+k)
            u(i,k,0) = tvar1rcv(i,1+kx*4+k)
            v(i,k,0) = tvar1rcv(i,1+kx*5+k)
            t(i,k,0) = tvar1rcv(i,1+kx*6+k)
            qv(i,k,0) = tvar1rcv(i,1+kx*7+k)
            qc(i,k,0) = tvar1rcv(i,1+kx*8+k)
            ua(i,k,0) = tvar1rcv(i,1+kx*9+k)
            va(i,k,0) = tvar1rcv(i,1+kx*10+k)
          end do
        end do
        if ( ichem.eq.1 ) then
          do n = 1 , ntr
            do k = 1 , kx
              do i = 1 , ix
                chi(i,k,0,n) = tvar1rcv(i,kx*11+1+(n-1)*kx+k)
              end do
            end do
          end do
        end if
      end if
!
      if ( myid.ne.0 ) then
        do i = 1 , ix
          tvar1snd(i,1) = psb(i,1)
        end do
        do k = 1 , kx
          do i = 1 , ix
            tvar1snd(i,1+k) = tb(i,k,1)
            tvar1snd(i,1+kx+k) = qvb(i,k,1)
            tvar1snd(i,1+kx*2+k) = ub(i,k,1)
            tvar1snd(i,1+kx*3+k) = vb(i,k,1)
            tvar1snd(i,1+kx*4+k) = u(i,k,1)
            tvar1snd(i,1+kx*5+k) = v(i,k,1)
            tvar1snd(i,1+kx*6+k) = t(i,k,1)
            tvar1snd(i,1+kx*7+k) = qv(i,k,1)
            tvar1snd(i,1+kx*8+k) = qc(i,k,1)
            tvar1snd(i,1+kx*9+k) = ua(i,k,1)
            tvar1snd(i,1+kx*10+k) = va(i,k,1)
          end do
        end do
        if ( ichem.eq.1 ) then
          do n = 1 , ntr
            do k = 1 , kx
              do i = 1 , ix
                tvar1snd(i,kx*11+1+(n-1)*kx+k) = chi(i,k,1,n)
              end do
            end do
          end do
        end if
      end if
      numrec = kx*11 + 1
      if ( ichem.eq.1 ) numrec = kx*(ntr+11) + 1
      call mpi_sendrecv(tvar1snd(1,1),ix*numrec,mpi_double_precision,   &
                      & iwest,2,tvar1rcv(1,1),ix*numrec,                &
                      & mpi_double_precision,ieast,2,mpi_comm_world,    &
                      & mpi_status_ignore,ierr)
      if ( myid.ne.nproc-1 ) then
        do i = 1 , ix
          psb(i,jxp+1) = tvar1rcv(i,1)
        end do
        do k = 1 , kx
          do i = 1 , ix
            tb(i,k,jxp+1) = tvar1rcv(i,1+k)
            qvb(i,k,jxp+1) = tvar1rcv(i,1+kx+k)
            ub(i,k,jxp+1) = tvar1rcv(i,1+kx*2+k)
            vb(i,k,jxp+1) = tvar1rcv(i,1+kx*3+k)
            u(i,k,jxp+1) = tvar1rcv(i,1+kx*4+k)
            v(i,k,jxp+1) = tvar1rcv(i,1+kx*5+k)
            t(i,k,jxp+1) = tvar1rcv(i,1+kx*6+k)
            qv(i,k,jxp+1) = tvar1rcv(i,1+kx*7+k)
            qc(i,k,jxp+1) = tvar1rcv(i,1+kx*8+k)
            ua(i,k,jxp+1) = tvar1rcv(i,1+kx*9+k)
            va(i,k,jxp+1) = tvar1rcv(i,1+kx*10+k)
          end do
        end do
        if ( ichem.eq.1 ) then
          do n = 1 , ntr
            do k = 1 , kx
              do i = 1 , ix
                chi(i,k,jxp+1,n) = tvar1rcv(i,kx*11+1+(n-1)*kx+k)
              end do
            end do
          end do
        end if
      end if
#endif
!=======================================================================
!
!-----interior points:
!
#ifdef MPP1
      do j = jbegin , jendx
#else
      do j = 2 , jxm1
#endif
        do i = 2 , ixm1
          pdotb(i,j) = 0.25*(psb(i,j)+psb(i-1,j)+psb(i,j-1)+psb(i-1,j-1)&
                     & )
        end do
      end do
!
!-----east and west boundaries:
!
      do i = 2 , ixm1
#ifdef MPP1
        if ( myid.eq.0 ) pdotb(i,1) = 0.5*(psb(i,1)+psb(i-1,1))
        if ( myid.eq.nproc-1 ) pdotb(i,jendl)                           &
           & = 0.5*(psb(i,jendx)+psb(i-1,jendx))
#else
        pdotb(i,1) = 0.5*(psb(i,1)+psb(i-1,1))
        pdotb(i,jx) = 0.5*(psb(i,jxm1)+psb(i-1,jxm1))
#endif
      end do
!
!-----north and south boundaries:
!
#ifdef MPP1
      do j = jbegin , jendx
#else
      do j = 2 , jxm1
#endif
        pdotb(1,j) = 0.5*(psb(1,j)+psb(1,j-1))
        pdotb(ix,j) = 0.5*(psb(ixm1,j)+psb(ixm1,j-1))
      end do
!
!-----corner points:
!
#ifdef MPP1
      if ( myid.eq.0 ) then
        pdotb(1,1) = psb(1,1)
        pdotb(ix,1) = psb(ixm1,1)
      end if
      if ( myid.eq.nproc-1 ) then
        pdotb(1,jendl) = psb(1,jendx)
        pdotb(ix,jendl) = psb(ixm1,jendx)
      end if
#else
      pdotb(1,1) = psb(1,1)
      pdotb(ix,1) = psb(ixm1,1)
      pdotb(1,jx) = psb(1,jxm1)
      pdotb(ix,jx) = psb(ixm1,jxm1)
#endif
!
!=======================================================================
!
      call slice3d
!
!=======================================================================
!
#ifdef MPP1
      if ( myid.ne.nproc-1 ) then
        do k = 1 , kx
          do i = 1 , ix
            var2snd(i,+k) = ubd3d(i,k,jxp-1)
            var2snd(i,kx+k) = ubd3d(i,k,jxp)
            var2snd(i,kx*2+k) = vbd3d(i,k,jxp-1)
            var2snd(i,kx*3+k) = vbd3d(i,k,jxp)
            var2snd(i,kx*4+k) = tb3d(i,k,jxp-1)
            var2snd(i,kx*5+k) = tb3d(i,k,jxp)
            var2snd(i,kx*6+k) = qvb3d(i,k,jxp-1)
            var2snd(i,kx*7+k) = qvb3d(i,k,jxp)
            var2snd(i,kx*8+k) = qcb3d(i,k,jxp-1)
            var2snd(i,kx*9+k) = qcb3d(i,k,jxp)
          end do
        end do
        if ( ichem.eq.1 ) then
          do n = 1 , ntr
            do k = 1 , kx
              do i = 1 , ix
                var2snd(i,kx*10+(n-1)*2*kx+k) = chib3d(i,k,jxp-1,n)
                var2snd(i,kx*10+(n-1)*2*kx+kx+k) = chib3d(i,k,jxp,n)
              end do
            end do
          end do
        end if
      end if
      numrec = kx*5*2
      if ( ichem.eq.1 ) numrec = kx*(ntr+5)*2
      call mpi_sendrecv(var2snd(1,1),ix*numrec,mpi_double_precision,    &
                      & ieast,1,var2rcv(1,1),ix*numrec,                 &
                      & mpi_double_precision,iwest,1,mpi_comm_world,    &
                      & mpi_status_ignore,ierr)
      if ( myid.ne.0 ) then
        do k = 1 , kx
          do i = 1 , ix
            ubd3d(i,k,-1) = var2rcv(i,+k)
            ubd3d(i,k,0) = var2rcv(i,kx+k)
            vbd3d(i,k,-1) = var2rcv(i,kx*2+k)
            vbd3d(i,k,0) = var2rcv(i,kx*3+k)
            tb3d(i,k,-1) = var2rcv(i,kx*4+k)
            tb3d(i,k,0) = var2rcv(i,kx*5+k)
            qvb3d(i,k,-1) = var2rcv(i,kx*6+k)
            qvb3d(i,k,0) = var2rcv(i,kx*7+k)
            qcb3d(i,k,-1) = var2rcv(i,kx*8+k)
            qcb3d(i,k,0) = var2rcv(i,kx*9+k)
          end do
        end do
        if ( ichem.eq.1 ) then
          do n = 1 , ntr
            do k = 1 , kx
              do i = 1 , ix
                chib3d(i,k,-1,n) = var2rcv(i,kx*10+(n-1)*2*kx+k)
                chib3d(i,k,0,n) = var2rcv(i,kx*10+(n-1)*2*kx+kx+k)
              end do
            end do
          end do
        end if
      end if
!
      if ( myid.ne.0 ) then
        do k = 1 , kx
          do i = 1 , ix
            var2snd(i,+k) = ubd3d(i,k,1)
            var2snd(i,kx+k) = ubd3d(i,k,2)
            var2snd(i,kx*2+k) = vbd3d(i,k,1)
            var2snd(i,kx*3+k) = vbd3d(i,k,2)
            var2snd(i,kx*4+k) = tb3d(i,k,1)
            var2snd(i,kx*5+k) = tb3d(i,k,2)
            var2snd(i,kx*6+k) = qvb3d(i,k,1)
            var2snd(i,kx*7+k) = qvb3d(i,k,2)
            var2snd(i,kx*8+k) = qcb3d(i,k,1)
            var2snd(i,kx*9+k) = qcb3d(i,k,2)
          end do
        end do
        if ( ichem.eq.1 ) then
          do n = 1 , ntr
            do k = 1 , kx
              do i = 1 , ix
                var2snd(i,kx*10+(n-1)*2*kx+k) = chib3d(i,k,1,n)
                var2snd(i,kx*10+(n-1)*2*kx+kx+k) = chib3d(i,k,2,n)
              end do
            end do
          end do
        end if
      end if
      numrec = kx*5*2
      if ( ichem.eq.1 ) numrec = kx*(ntr+5)*2
      call mpi_sendrecv(var2snd(1,1),ix*numrec,mpi_double_precision,    &
                      & iwest,2,var2rcv(1,1),ix*numrec,                 &
                      & mpi_double_precision,ieast,2,mpi_comm_world,    &
                      & mpi_status_ignore,ierr)
      if ( myid.ne.nproc-1 ) then
        do k = 1 , kx
          do i = 1 , ix
            ubd3d(i,k,jxp+1) = var2rcv(i,+k)
            ubd3d(i,k,jxp+2) = var2rcv(i,kx+k)
            vbd3d(i,k,jxp+1) = var2rcv(i,kx*2+k)
            vbd3d(i,k,jxp+2) = var2rcv(i,kx*3+k)
            tb3d(i,k,jxp+1) = var2rcv(i,kx*4+k)
            tb3d(i,k,jxp+2) = var2rcv(i,kx*5+k)
            qvb3d(i,k,jxp+1) = var2rcv(i,kx*6+k)
            qvb3d(i,k,jxp+2) = var2rcv(i,kx*7+k)
            qcb3d(i,k,jxp+1) = var2rcv(i,kx*8+k)
            qcb3d(i,k,jxp+2) = var2rcv(i,kx*9+k)
          end do
        end do
        if ( ichem.eq.1 ) then
          do n = 1 , ntr
            do k = 1 , kx
              do i = 1 , ix
                chib3d(i,k,jxp+1,n) = var2rcv(i,kx*10+(n-1)*2*kx+k)
                chib3d(i,k,jxp+2,n) = var2rcv(i,kx*10+(n-1)*2*kx+kx+k)
              end do
            end do
          end do
        end if
      end if
#endif
!
!**********************************************************************
!***** "j" loop begins here:
!
#ifdef MPP1
      do j = 1 , jendx
#else
      do j = 1 , jxm1
#endif
!
        icon(j) = 0
#ifdef MPP1
        if ( (myid.eq.0 .and. j.eq.1) .or.                              &
           & (myid.eq.nproc-1 .and. j.eq.jendx) ) then
#else
        if ( j.eq.1 .or. j.eq.jxm1 ) then
#endif
          do k = 1 , kxp1
            do i = 1 , ixm1
              qdot(i,k,j) = 0.
            end do
          end do
        else
!
!----------------------------------------------------------------------
!**p**compute the pressure tendency:
!
          do i = 2 , ixm2
            pten(i,j) = 0.
          end do
          do k = 1 , kx
            do i = 2 , ixm2
              divl(i,k) = (ua(i+1,k,j+1)+ua(i,k,j+1)-ua(i+1,k,j)        &
                        & -ua(i,k,j))                                   &
                        & + (va(i+1,k,j+1)+va(i+1,k,j)-va(i,k,j+1)      &
                        & -va(i,k,j))
              pten(i,j) = pten(i,j) - divl(i,k)*dsigma(k)               &
                        & /(dx2*msfx(i,j)*msfx(i,j))
            end do
          end do
!
!..p..compute vertical sigma-velocity (qdot):
!
          do k = 1 , kxp1
            do i = 1 , ixm1
              qdot(i,k,j) = 0.
            end do
          end do
          do k = 2 , kx
            do i = 2 , ixm2
              qdot(i,k,j) = qdot(i,k-1,j)                               &
                          & - (pten(i,j)+divl(i,k-1)/(dx2*msfx(i,j)     &
                          & *msfx(i,j)))*dsigma(k-1)/psa(i,j)
            end do
          end do
        end if
      end do
#ifdef MPP1
      call mpi_sendrecv(qdot(1,1,jxp),ix*kxp1,mpi_double_precision,     &
                      & ieast,1,qdot(1,1,0),ix*kxp1,                    &
                      & mpi_double_precision,iwest,1,mpi_comm_world,    &
                      & mpi_status_ignore,ierr)
      call mpi_sendrecv(qdot(1,1,1),ix*kxp1,mpi_double_precision,iwest, &
                      & 2,qdot(1,1,jxp+1),ix*kxp1,mpi_double_precision, &
                      & ieast,2,mpi_comm_world,mpi_status_ignore,ierr)
#endif
!
!..p..compute omega
!
#ifdef MPP1
      do j = 1 , jendx
        if ( (myid.eq.0 .and. j.eq.1) .or.                              &
           & (myid.eq.nproc-1 .and. j.eq.jendx) ) then
#else
      do j = 1 , jxm1
        if ( j.eq.1 .or. j.eq.jxm1 ) then
#endif
          do k = 1 , kx
            do i = 1 , ixm1
              omega(i,k,j) = 0.
            end do
          end do
        else
          do k = 1 , kx
            do i = 2 , ixm2
              omega(i,k,j) = 0.5*psa(i,j)*(qdot(i,k+1,j)+qdot(i,k,j))   &
                           & + a(k)                                     &
                           & *(pten(i,j)+((u(i,k,j)+u(i+1,k,j)+u(i+1,k, &
                           & j+1)+u(i,k,j+1))*(psa(i,j+1)-psa(i,j-1))   &
                           & +(v(i,k,j)+v(i+1,k,j)+v(i+1,k,j+1)         &
                           & +v(i,k,j+1))*(psa(i+1,j)-psa(i-1,j)))      &
                           & /(dx8*msfx(i,j)))
            end do
          end do
        end if
      end do
#ifdef MPP1
      if ( nspgx.ge.jxp ) then
        do i = 1 , ix
          bdyewsnd(i,1) = peb(i,1)
          bdyewsnd(i,2) = pwb(i,jxp)
          bdyewsnd(i,3) = pebt(i,1)
          bdyewsnd(i,4) = pwbt(i,jxp)
        end do
        do k = 1 , kx
          do i = 1 , ix
            bdyewsnd(i,4+k) = teb(i,k,1)
            bdyewsnd(i,4+kx+k) = twb(i,k,jxp)
            bdyewsnd(i,4+kx*2+k) = tebt(i,k,1)
            bdyewsnd(i,4+kx*3+k) = twbt(i,k,jxp)
            bdyewsnd(i,4+kx*4+k) = qeb(i,k,1)
            bdyewsnd(i,4+kx*5+k) = qwb(i,k,jxp)
            bdyewsnd(i,4+kx*6+k) = qebt(i,k,1)
            bdyewsnd(i,4+kx*7+k) = qwbt(i,k,jxp)
            bdyewsnd(i,4+kx*8+k) = ueb(i,k,1)
            bdyewsnd(i,4+kx*9+k) = uwb(i,k,jxp)
            bdyewsnd(i,4+kx*10+k) = uebt(i,k,1)
            bdyewsnd(i,4+kx*11+k) = uwbt(i,k,jxp)
            bdyewsnd(i,4+kx*12+k) = veb(i,k,1)
            bdyewsnd(i,4+kx*13+k) = vwb(i,k,jxp)
            bdyewsnd(i,4+kx*14+k) = vebt(i,k,1)
            bdyewsnd(i,4+kx*15+k) = vwbt(i,k,jxp)
          end do
        end do
        call mpi_sendrecv(bdyewsnd(1,1),ix*(kx*16+4),                   &
                        & mpi_double_precision,ieast,1,bdyewrcv(1,1),   &
                        & ix*(kx*16+4),mpi_double_precision,iwest,1,    &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        do i = 1 , ix
          if ( myid.eq.nproc-1 ) then
            peb(i,jendl) = bdyewrcv(i,1)
            pebt(i,jendl) = bdyewrcv(i,3)
          else
            peb(i,jxp+1) = bdyewrcv(i,1)
            pebt(i,jxp+1) = bdyewrcv(i,3)
          end if
          pwb(i,0) = bdyewrcv(i,2)
          pwbt(i,0) = bdyewrcv(i,4)
        end do
        do k = 1 , kx
          do i = 1 , ix
            if ( myid.eq.nproc-1 ) then
              teb(i,k,jendl) = bdyewrcv(i,4+k)
              tebt(i,k,jendl) = bdyewrcv(i,4+kx*2+k)
              qeb(i,k,jendl) = bdyewrcv(i,4+kx*4+k)
              qebt(i,k,jendl) = bdyewrcv(i,4+kx*6+k)
            else
              teb(i,k,jxp+1) = bdyewrcv(i,4+k)
              tebt(i,k,jxp+1) = bdyewrcv(i,4+kx*2+k)
              qeb(i,k,jxp+1) = bdyewrcv(i,4+kx*4+k)
              qebt(i,k,jxp+1) = bdyewrcv(i,4+kx*6+k)
            end if
            ueb(i,k,jxp+1) = bdyewrcv(i,4+kx*8+k)
            uebt(i,k,jxp+1) = bdyewrcv(i,4+kx*10+k)
            veb(i,k,jxp+1) = bdyewrcv(i,4+kx*12+k)
            vebt(i,k,jxp+1) = bdyewrcv(i,4+kx*14+k)
            twb(i,k,0) = bdyewrcv(i,4+kx+k)
            twbt(i,k,0) = bdyewrcv(i,4+kx*3+k)
            qwb(i,k,0) = bdyewrcv(i,4+kx*5+k)
            qwbt(i,k,0) = bdyewrcv(i,4+kx*7+k)
            uwb(i,k,0) = bdyewrcv(i,4+kx*9+k)
            uwbt(i,k,0) = bdyewrcv(i,4+kx*11+k)
            vwb(i,k,0) = bdyewrcv(i,4+kx*13+k)
            vwbt(i,k,0) = bdyewrcv(i,4+kx*15+k)
          end do
        end do
        do i = 1 , ix
          if ( myid.eq.nproc-1 ) then
            bdyewsnd(i,1) = peb(i,jendx)
            bdyewsnd(i,3) = pebt(i,jendx)
          else
            bdyewsnd(i,1) = peb(i,jxp)
            bdyewsnd(i,3) = pebt(i,jxp)
          end if
          bdyewsnd(i,2) = pwb(i,1)
          bdyewsnd(i,4) = pwbt(i,1)
        end do
        do k = 1 , kx
          do i = 1 , ix
            if ( myid.eq.nproc-1 ) then
              bdyewsnd(i,4+k) = teb(i,k,jendx)
              bdyewsnd(i,4+kx*2+k) = tebt(i,k,jendx)
              bdyewsnd(i,4+kx*4+k) = qeb(i,k,jendx)
              bdyewsnd(i,4+kx*6+k) = qebt(i,k,jendx)
            else
              bdyewsnd(i,4+k) = teb(i,k,jxp)
              bdyewsnd(i,4+kx*2+k) = tebt(i,k,jxp)
              bdyewsnd(i,4+kx*4+k) = qeb(i,k,jxp)
              bdyewsnd(i,4+kx*6+k) = qebt(i,k,jxp)
            end if
            bdyewsnd(i,4+kx*8+k) = ueb(i,k,jxp)
            bdyewsnd(i,4+kx*10+k) = uebt(i,k,jxp)
            bdyewsnd(i,4+kx*12+k) = veb(i,k,jxp)
            bdyewsnd(i,4+kx*14+k) = vebt(i,k,jxp)
            bdyewsnd(i,4+kx+k) = twb(i,k,1)
            bdyewsnd(i,4+kx*3+k) = twbt(i,k,1)
            bdyewsnd(i,4+kx*5+k) = qwb(i,k,1)
            bdyewsnd(i,4+kx*7+k) = qwbt(i,k,1)
            bdyewsnd(i,4+kx*9+k) = uwb(i,k,1)
            bdyewsnd(i,4+kx*11+k) = uwbt(i,k,1)
            bdyewsnd(i,4+kx*13+k) = vwb(i,k,1)
            bdyewsnd(i,4+kx*15+k) = vwbt(i,k,1)
          end do
        end do
        call mpi_sendrecv(bdyewsnd(1,1),ix*(kx*16+4),                   &
                        & mpi_double_precision,iwest,2,bdyewrcv(1,1),   &
                        & ix*(kx*16+4),mpi_double_precision,ieast,2,    &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        do i = 1 , ix
          peb(i,0) = bdyewrcv(i,1)
          pebt(i,0) = bdyewrcv(i,3)
          pwb(i,jxp+1) = bdyewrcv(i,2)
          pwbt(i,jxp+1) = bdyewrcv(i,4)
        end do
        do k = 1 , kx
          do i = 1 , ix
            teb(i,k,0) = bdyewrcv(i,4+k)
            twb(i,k,jxp+1) = bdyewrcv(i,4+kx+k)
            tebt(i,k,0) = bdyewrcv(i,4+kx*2+k)
            twbt(i,k,jxp+1) = bdyewrcv(i,4+kx*3+k)
            qeb(i,k,0) = bdyewrcv(i,4+kx*4+k)
            qwb(i,k,jxp+1) = bdyewrcv(i,4+kx*5+k)
            qebt(i,k,0) = bdyewrcv(i,4+kx*6+k)
            qwbt(i,k,jxp+1) = bdyewrcv(i,4+kx*7+k)
            ueb(i,k,0) = bdyewrcv(i,4+kx*8+k)
            uwb(i,k,jxp+1) = bdyewrcv(i,4+kx*9+k)
            uebt(i,k,0) = bdyewrcv(i,4+kx*10+k)
            uwbt(i,k,jxp+1) = bdyewrcv(i,4+kx*11+k)
            veb(i,k,0) = bdyewrcv(i,4+kx*12+k)
            vwb(i,k,jxp+1) = bdyewrcv(i,4+kx*13+k)
            vebt(i,k,0) = bdyewrcv(i,4+kx*14+k)
            vwbt(i,k,jxp+1) = bdyewrcv(i,4+kx*15+k)
          end do
        end do
      end if
!
      if ( myid.ne.nproc-1 ) then
        do i = 1 , nspgx
          bdynssnd(i,1) = pnb(i,jxp)
          bdynssnd(i,2) = pnbt(i,jxp)
          bdynssnd(i,3) = pss(i,jxp)
          bdynssnd(i,4) = psbt(i,jxp)
        end do
        do k = 1 , kx
          do i = 1 , nspgx
            bdynssnd(i,4+k) = tnb(i,k,jxp)
            bdynssnd(i,4+kx+k) = tnbt(i,k,jxp)
            bdynssnd(i,4+kx*2+k) = tsb(i,k,jxp)
            bdynssnd(i,4+kx*3+k) = tsbt(i,k,jxp)
            bdynssnd(i,4+kx*4+k) = qnb(i,k,jxp)
            bdynssnd(i,4+kx*5+k) = qnbt(i,k,jxp)
            bdynssnd(i,4+kx*6+k) = qsb(i,k,jxp)
            bdynssnd(i,4+kx*7+k) = qsbt(i,k,jxp)
            bdynssnd(i,4+kx*8+k) = unb(i,k,jxp)
            bdynssnd(i,4+kx*9+k) = unbt(i,k,jxp)
            bdynssnd(i,4+kx*10+k) = usb(i,k,jxp)
            bdynssnd(i,4+kx*11+k) = usbt(i,k,jxp)
            bdynssnd(i,4+kx*12+k) = vnb(i,k,jxp)
            bdynssnd(i,4+kx*13+k) = vnbt(i,k,jxp)
            bdynssnd(i,4+kx*14+k) = vsb(i,k,jxp)
            bdynssnd(i,4+kx*15+k) = vsbt(i,k,jxp)
          end do
        end do
      end if
      call mpi_sendrecv(bdynssnd(1,1),nspgx*(kx*16+4),                  &
                      & mpi_double_precision,ieast,1,bdynsrcv(1,1),     &
                      & nspgx*(kx*16+4),mpi_double_precision,iwest,1,   &
                      & mpi_comm_world,mpi_status_ignore,ierr)
      if ( myid.ne.0 ) then
        do i = 1 , nspgx
          pnb(i,0) = bdynsrcv(i,1)
          pnbt(i,0) = bdynsrcv(i,2)
          pss(i,0) = bdynsrcv(i,3)
          psbt(i,0) = bdynsrcv(i,4)
        end do
        do k = 1 , kx
          do i = 1 , nspgx
            tnb(i,k,0) = bdynsrcv(i,4+k)
            tnbt(i,k,0) = bdynsrcv(i,4+kx+k)
            tsb(i,k,0) = bdynsrcv(i,4+kx*2+k)
            tsbt(i,k,0) = bdynsrcv(i,4+kx*3+k)
            qnb(i,k,0) = bdynsrcv(i,4+kx*4+k)
            qnbt(i,k,0) = bdynsrcv(i,4+kx*5+k)
            qsb(i,k,0) = bdynsrcv(i,4+kx*6+k)
            qsbt(i,k,0) = bdynsrcv(i,4+kx*7+k)
            unb(i,k,0) = bdynsrcv(i,4+kx*8+k)
            unbt(i,k,0) = bdynsrcv(i,4+kx*9+k)
            usb(i,k,0) = bdynsrcv(i,4+kx*10+k)
            usbt(i,k,0) = bdynsrcv(i,4+kx*11+k)
            vnb(i,k,0) = bdynsrcv(i,4+kx*12+k)
            vnbt(i,k,0) = bdynsrcv(i,4+kx*13+k)
            vsb(i,k,0) = bdynsrcv(i,4+kx*14+k)
            vsbt(i,k,0) = bdynsrcv(i,4+kx*15+k)
          end do
        end do
      end if
!
      if ( myid.ne.0 ) then
        do i = 1 , nspgx
          bdynssnd(i,1) = pnb(i,1)
          bdynssnd(i,2) = pnbt(i,1)
          bdynssnd(i,3) = pss(i,1)
          bdynssnd(i,4) = psbt(i,1)
        end do
        do k = 1 , kx
          do i = 1 , nspgx
            bdynssnd(i,4+k) = tnb(i,k,1)
            bdynssnd(i,4+kx+k) = tnbt(i,k,1)
            bdynssnd(i,4+kx*2+k) = tsb(i,k,1)
            bdynssnd(i,4+kx*3+k) = tsbt(i,k,1)
            bdynssnd(i,4+kx*4+k) = qnb(i,k,1)
            bdynssnd(i,4+kx*5+k) = qnbt(i,k,1)
            bdynssnd(i,4+kx*6+k) = qsb(i,k,1)
            bdynssnd(i,4+kx*7+k) = qsbt(i,k,1)
            bdynssnd(i,4+kx*8+k) = unb(i,k,1)
            bdynssnd(i,4+kx*9+k) = unbt(i,k,1)
            bdynssnd(i,4+kx*10+k) = usb(i,k,1)
            bdynssnd(i,4+kx*11+k) = usbt(i,k,1)
            bdynssnd(i,4+kx*12+k) = vnb(i,k,1)
            bdynssnd(i,4+kx*13+k) = vnbt(i,k,1)
            bdynssnd(i,4+kx*14+k) = vsb(i,k,1)
            bdynssnd(i,4+kx*15+k) = vsbt(i,k,1)
          end do
        end do
      end if
      call mpi_sendrecv(bdynssnd(1,1),nspgx*(kx*16+4),                  &
                      & mpi_double_precision,iwest,2,bdynsrcv(1,1),     &
                      & nspgx*(kx*16+4),mpi_double_precision,ieast,2,   &
                      & mpi_comm_world,mpi_status_ignore,ierr)
      if ( myid.ne.nproc-1 ) then
        do i = 1 , nspgx
          pnb(i,jxp+1) = bdynsrcv(i,1)
          pnbt(i,jxp+1) = bdynsrcv(i,2)
          pss(i,jxp+1) = bdynsrcv(i,3)
          psbt(i,jxp+1) = bdynsrcv(i,4)
        end do
        do k = 1 , kx
          do i = 1 , nspgx
            tnb(i,k,jxp+1) = bdynsrcv(i,4+k)
            tnbt(i,k,jxp+1) = bdynsrcv(i,4+kx+k)
            tsb(i,k,jxp+1) = bdynsrcv(i,4+kx*2+k)
            tsbt(i,k,jxp+1) = bdynsrcv(i,4+kx*3+k)
            qnb(i,k,jxp+1) = bdynsrcv(i,4+kx*4+k)
            qnbt(i,k,jxp+1) = bdynsrcv(i,4+kx*5+k)
            qsb(i,k,jxp+1) = bdynsrcv(i,4+kx*6+k)
            qsbt(i,k,jxp+1) = bdynsrcv(i,4+kx*7+k)
            unb(i,k,jxp+1) = bdynsrcv(i,4+kx*8+k)
            unbt(i,k,jxp+1) = bdynsrcv(i,4+kx*9+k)
            usb(i,k,jxp+1) = bdynsrcv(i,4+kx*10+k)
            usbt(i,k,jxp+1) = bdynsrcv(i,4+kx*11+k)
            vnb(i,k,jxp+1) = bdynsrcv(i,4+kx*12+k)
            vnbt(i,k,jxp+1) = bdynsrcv(i,4+kx*13+k)
            vsb(i,k,jxp+1) = bdynsrcv(i,4+kx*14+k)
            vsbt(i,k,jxp+1) = bdynsrcv(i,4+kx*15+k)
          end do
        end do
      end if
#endif

#ifdef MPP1
      do j = jbegin , jendx
        if ( myid.ne.nproc-1 .or. j.ne.jendx ) then
#else
      do j = 2 , jxm1
        if ( j.ne.jxm1 ) then
#endif
!EES      omega change: I broke up the write onto two lines
!         and commented out the if statment
!
          if ( iboudy.eq.4 ) then
!..p..apply sponge boundary conditions to pten:
            call sponge_p(ispgx,wgtx,pten(1,j),j)
!....apply  the nudging boundary conditions:
          else if ( iboudy.eq.1 .or. iboudy.eq.5 ) then
            xtm1 = xtime - dtmin
            if ( dabs(xtime).lt.0.00001 .and. ldatez.gt.idate0 )        &
               & xtm1 = -dtmin
            call nudge_p(ispgx,fnudge,gnudge,xtm1,pten(1,j),c203,j,     &
                       & iboudy)
          else
          end if
        end if     !end if(j.ne.jxm1) test
      end do
#ifdef MPP1
      do j = 1 , jendx
        if ( myid.eq.0 .and. j.eq.1 ) then
#else
      do j = 1 , jxm1
        if ( j.eq.1 ) then
#endif
          do i = 1 , ixm1
            psc(i,j) = psb(i,j) + dt*pwbt(i,j)
            psd(i,j) = psa(i,j)
          end do
#ifdef MPP1
        else if ( myid.eq.nproc-1 .and. j.eq.jendx ) then
#else
        else if ( j.eq.jxm1 ) then
#endif
          do i = 1 , ixm1
            psd(i,j) = psa(i,j)
          end do
        else
!
!..p..forecast pressure:
!
          do i = 2 , ixm2
            psc(i,j) = psb(i,j) + pten(i,j)*dt
          end do
!
!..p..weighted p* (psd)
!
          do i = 2 , ixm2
            psd(i,j) = psa(i,j)
          end do
!
          psc(1,j) = psb(1,j) + dt*psbt(1,j)
          psc(ixm1,j) = psb(ixm1,j) + dt*pnbt(1,j)
          psd(1,j) = psa(1,j)
          psd(ixm1,j) = psa(ixm1,j)
        end if
      end do
#ifdef MPP1
      call mpi_sendrecv(psd(1,jxp),ix,mpi_double_precision,ieast,1,     &
                      & psd(1,0),ix,mpi_double_precision,iwest,1,       &
                      & mpi_comm_world,mpi_status_ignore,ierr)
#endif
!
!-----compute bleck (1977) noise parameters:
!
#ifdef MPP1
      do j = 1 , jendl
        do i = 1 , ix
          ps4(i,1,j) = pten(i,j)
          ps4(i,2,j) = psc(i,j)
          ps4(i,3,j) = psb(i,j)
          ps4(i,4,j) = psa(i,j)
        end do
      end do
      call mpi_gather(ps4(1,1,1),ix*4*jxp,mpi_double_precision,         &
                    & ps_4(1,1,1),ix*4*jxp,mpi_double_precision,0,      &
                    & mpi_comm_world,ierr)
#endif
      iptn = 0
      ptntot = 0.
      pt2tot = 0.
#ifdef MPP1
      if ( myid.eq.0 ) then
        do j = 2 , jxm2
          if ( jyear.ne.jyear0 .or. ktau.ne.0 ) then
            do i = 2 , ixm2
              iptn = iptn + 1
              ptntot = ptntot + dabs(ps_4(i,1,j))
              pt2tot = pt2tot +                                         &
                     & dabs((ps_4(i,2,j)+ps_4(i,3,j)-2.*ps_4(i,4,j))    &
                     & /(0.25*dt*dt))
            end do
          end if
        end do
      end if
      call mpi_bcast(iptn,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(ptntot,1,mpi_double_precision,0,mpi_comm_world,    &
                   & ierr)
      call mpi_bcast(pt2tot,1,mpi_double_precision,0,mpi_comm_world,    &
                   & ierr)
#else
      do j = 2 , jxm1
        if ( j.ne.jxm1 ) then
          if ( jyear.ne.jyear0 .or. ktau.ne.0 ) then
            do i = 2 , ixm2
              iptn = iptn + 1
              ptntot = ptntot + dabs(pten(i,j))
              pt2tot = pt2tot + dabs((psc(i,j)+psb(i,j)-2.*psa(i,j))    &
                     & /(0.25*dt*dt))
            end do
          end if
        end if
      end do
#endif
!
#ifdef MPP1
      do j = jbegin , jendx
#else
      do j = 2 , jxm1
#endif

!
!------compute the horizontal diffusion coefficient and stored in xkc:
!       the values are calculated at cross points, but they also used
!       for dot-point variables.
 
        do k = 1 , kx
          do i = 2 , ixm1
            dudx = ub(i,k,j+1) + ub(i+1,k,j+1) - ub(i,k,j) - ub(i+1,k,j)
            dvdx = vb(i,k,j+1) + vb(i+1,k,j+1) - vb(i,k,j) - vb(i+1,k,j)
            dudy = ub(i+1,k,j) + ub(i+1,k,j+1) - ub(i,k,j) - ub(i,k,j+1)
            dvdy = vb(i+1,k,j) + vb(i+1,k,j+1) - vb(i,k,j) - vb(i,k,j+1)
!fil        cell=(xkhz*hgfact(i,j)/5.+c200*dsqrt((dudx-dvdy)*(dudx-dvdy)
            cell = (xkhz*hgfact(i,j)                                    &
                 & +c200*dsqrt((dudx-dvdy)*(dudx-dvdy)+(dvdx+dudy)      &
                 & *(dvdx+dudy)))
            xkc(i,k,j) = dmin1(cell,xkhmax)
 
          end do
        end do
      end do
!
#ifdef MPP1
      do j = jbegin , jendx
        if ( myid.ne.nproc-1 .or. j.ne.jendx ) then
#else
      do j = 2 , jxm1
        if ( j.ne.jxm1 ) then
#endif
!
!---------------------------------------------------------------------
!**t**compute the temperature tendency:
!
          do k = 1 , kx
            do i = 2 , ixm2
              tten(i,k,j) = 0.
              qvten(i,k,j) = 0.
              qcten(i,k,j) = 0.
            end do
          end do
 
!
!..t..compute the horizontal advection term:
!
!         call hadv_T(tten(1,1,j),ua,va,t,msfx,dx4,j,1)
          call hadv_t(tten(1,1,j),dx4,j,1)
!
!..t..compute the vertical advection term:
!
          call vadv(tten(1,1,j),ta(1,1,j),j,1)
!
!..t..compute the adiabatic term:
!
          do k = 1 , kx
            do i = 2 , ixm2
              rovcpm = r/(cp*(1.+0.8*(qv(i,k,j))))
              tv = t(i,k,j)*(1.+ep1*(qv(i,k,j)))
              tten(i,k,j) = tten(i,k,j) + (omega(i,k,j)*rovcpm*tv)      &
                          & /(ptop/psa(i,j)+a(k))
            end do
          end do
!
!..t..compute the diffusion term for t and store in difft:
!
          do k = 1 , kx
            do i = 1 , ixm1
              difft(i,k,j) = 0.
              diffq(i,k,j) = 0.
            end do
          end do
!
          call diffut_t(difft(1,1,j),xkc(1,1,j),c203,j)
!
!**q**compute the moisture tendencies:
!
!....icup = 1 : kuo-anthes cumulus parameterizaion scheme
!....icup = 2 : grell cumulus paramterization scheme
!....icup = 3 : betts-miller (1986)
!....icup = 4 : emanuel (1991)
!
          if ( icup.ne.1 ) then
            call hadvqv(qvten(1,1,j),dx4,j,1)
            call vadv(qvten(1,1,j),qva(1,1,j),j,2)
          end if
 
          if ( icup.eq.1 ) then
            call cupara(j)
          else if ( icup.eq.2 ) then
            call cuparan(tten(1,1,j),qvten(1,1,j),j)
          else if ( icup.eq.3 ) then
            write (aline,*)                                             &
                         &'ICTP RegCM team thinks the Betts-Miller code'&
                        & ,                                             &
                         &' is not ready for Regional Climate Run yet.'
            call say
            call fatal(__FILE__,__LINE__,                               &
                      &'BETTS MILLER CUMULUS OPTION NOT ALLOWED')
            call bmpara(tten(1,1,j),qvten(1,1,j),j)
          else if ( icup.eq.4 ) then
            call cupemandrv(j)
          else
          end if
 
          if ( ipptls.eq.1 ) then
            call hadvqc(qcten(1,1,j),dx4,j,1)
!           call hadvQC(qcten(1,1,j),dx,j,2)
!fix        call vadv(qcten(1,1,j),qca(1,1,j),j,3)
            call vadv(qcten(1,1,j),qca(1,1,j),j,5)
            call pcp(j)
            call cldfrac(j)
 
!           need also to set diffq to 0 here before calling diffut
            do k = 1 , kx
              do i = 1 , ixm1
                diffq(i,k,j) = 0.
              end do
            end do
 
!-----compute the diffusion terms:
!           the diffusion term for qv is stored in diffq. before
!           completing qvten computation, do not use diffq for other
 
!           purpose.
            call diffutqv(diffq(1,1,j),xkc(1,1,j),c203,j)
            call diffutqc(qcten(1,1,j),xkc(1,1,j),c203,j)
          end if
!
!chem2    compute the tracers tendencies
          if ( ichem.eq.1 ) then
            call zenitm(coszrs,ix,j)
            call tractend2(j)
          end if
!chem2_
!
        end if           !end if(j.ne.jxm1) test
!----------------------------------------------------------------------
!*****compute the pbl fluxes:
!       the diffusion and pbl tendencies of t and qv are stored in
!       difft and diffq.
!
        do k = 1 , kx
          do i = 2 , ixm1
            uten(i,k,j) = 0.
            vten(i,k,j) = 0.
          end do
        end do
 
 
!       ****** calculate solar zenith angle
        if ( (jyear.eq.jyear0 .and. ktau.eq.0) .or. mod(ktau+1,nbatst)  &
           & .eq.0 .or. mod(ktau+1,ntrad).eq.0 ) then
          call zenitm(coszrs,ix,j)
          call slice1D(j)
        end if
 
!       ****** calculate albedo
        if ( (jyear.eq.jyear0 .and. ktau.eq.0) .or. mod(ktau+1,ntrad)   &
           & .eq.0 ) call albedov(j,iemiss)
 
!       ****** call ccm3 radiative transfer package
        if ( (jyear.eq.jyear0 .and. ktau.eq.0) .or. mod(ktau+1,ntrad)   &
           & .eq.0 ) call colmod3(j)
 
!       ****** call vector bats for surface physics calculations
        if ( (jyear.eq.jyear0 .and. ktau.eq.0) .or. mod(ktau+1,nbatst)  &
           & .eq.0 ) then
          dtbat = dt/2.*nbatst
          if ( jyear.eq.jyear0 .and. ktau.eq.0 ) dtbat = dt
          call vecbats(j)
          if ( iocnflx.eq.2 ) call zengocndrv(j)
                                      ! Zeng ocean flux model
!         ****** accumulate quantities for energy and moisture budgets
          call interf(2,j)
        end if
 
      end do
      if ( icup.eq.1 ) then
        dto2 = dt/2
        call htdiff(dto2,dxsq,akht1)
      end if
!     call medium resolution pbl
      if ( ibltyp.eq.1 ) call holtbl
!
#ifdef MPP1
      do j = jbegin , jendx
#else
      do j = 2 , jxm1
#endif
!       add ccm radiative transfer package-calculated heating rates to
!       temperature tendency
        do k = 1 , kx
          do i = 2 , ixm2
            tten(i,k,j) = tten(i,k,j) + psb(i,j)*heatrt(i,k,j)
                                                       !heating rate in deg/sec
          end do
        end do
!
#ifdef MPP1
        if ( myid.ne.nproc-1 .or. j.ne.jendx ) then
#else
        if ( j.ne.jxm1 ) then
#endif
!
!..tq.add horizontal diffusion and pbl tendencies for t and qv to tten
!         and qvten for calculating condensational term in subroutine
!         "condtq".
!
          do k = 1 , kx
            do i = 2 , ixm2
              tten(i,k,j) = tten(i,k,j) + difft(i,k,j)
            end do
          end do
!
          do k = 1 , kx
            do i = 2 , ixm2
              qvten(i,k,j) = qvten(i,k,j) + diffq(i,k,j)
            end do
          end do
!
!..tq.compute the condensation and precipitation terms for explicit
!         moisture scheme:
!
          call condtq(j)
!
!..tq.subtract horizontal diffusion and pbl tendencies from tten and
!         qvten for appling the sponge boundary conditions on t and qv:
!
          if ( iboudy.eq.4 ) then
            do k = 1 , kx
              do i = 2 , ixm2
                tten(i,k,j) = tten(i,k,j) - difft(i,k,j)
              end do
            end do
            call sponge_t(ispgx,wgtx,tten(1,1,j),j)
            do k = 1 , kx
              do i = 2 , ixm2
                tten(i,k,j) = tten(i,k,j) + difft(i,k,j)
              end do
            end do
            do k = 1 , kx
              do i = 2 , ixm2
                qvten(i,k,j) = qvten(i,k,j) - diffq(i,k,j)
              end do
            end do
            call spongeqv(ispgx,wgtx,qvten(1,1,j),j)
            do k = 1 , kx
              do i = 2 , ixm2
                qvten(i,k,j) = qvten(i,k,j) + diffq(i,k,j)
              end do
            end do
          end if
!
!..tq.apply the nudging boundary conditions:
!
          if ( iboudy.eq.1 .or. iboudy.eq.5 ) then
            xtm1 = xtime - dtmin
            if ( dabs(xtime).lt.0.00001 .and. ldatez.gt.idate0 )        &
               & xtm1 = -dtmin
            call nudge_t(ispgx,fnudge,gnudge,xtm1,tten(1,1,j),c203,j,   &
                       & iboudy)
            call nudgeqv(ispgx,fnudge,gnudge,xtm1,qvten(1,1,j),c203,j,  &
                       & iboudy)
          end if
!
!..tq.forecast t, qv, and qc at tau+1:
!
          do k = 1 , kx
            do i = 2 , ixm2
              qvc(i,k,j) = qvb(i,k,j) + dt*qvten(i,k,j)
            end do
          end do
!
          do k = 1 , kx
            do i = 2 , ixm2
              qcc(i,k,j) = qcb(i,k,j) + dt*qcten(i,k,j)
            end do
          end do
!
          do k = 1 , kx
            do i = 2 , ixm2
              tc(i,k,j) = tb(i,k,j) + dt*tten(i,k,j)
            end do
          end do
!
!chem2
!         forecast tracer chi at at tau+1:
          if ( ichem.eq.1 ) then
!
            do itr = 1 , ntr
              do k = 1 , kx
                do i = 2 , ixm2
                  chic(i,k,j,itr) = chib(i,k,j,itr)                     &
                                  & + dt*chiten(i,k,j,itr)
                end do
              end do
            end do
          end if
!chem2_
        end if       !end if(j.ne.jxm1),else test
      end do
!
#ifdef MPP1
      do j = 1 , jendx
        if ( myid.eq.0 .and. j.eq.1 ) then
#else
      do j = 1 , jxm1
        if ( j.eq.1 ) then
#endif
          if ( ipgf.eq.1 ) then
            do k = 1 , kx
              do i = 1 , ixm1
                td(i,k,j) = ta(i,k,j)*(1.+ep1*(qv(i,k,j)))
                ttld(i,k,j) = td(i,k,j) - psa(i,j)                      &
                            & *t00pg*((a(k)*psa(i,j)+ptop)/p00pg)       &
                            & **pgfaa1
              end do
            end do
          else if ( ipgf.eq.0 ) then
            do k = 1 , kx
              do i = 1 , ixm1
                td(i,k,j) = ta(i,k,j)*(1.+ep1*(qv(i,k,j)))
              end do
            end do
          else
          end if
!
#ifdef MPP1
        else if ( myid.eq.nproc-1 .and. j.eq.jendx ) then
#else
        else if ( j.eq.jxm1 ) then
#endif
!
!-----set td and psd at j=jlx equal to ta and psa:
!
          if ( ipgf.eq.1 ) then
            do k = 1 , kx
              do i = 1 , ixm1
                td(i,k,j) = ta(i,k,j)*(1.+ep1*(qv(i,k,j)))
                ttld(i,k,j) = td(i,k,j) - psa(i,j)                      &
                            & *t00pg*((a(k)*psa(i,j)+ptop)/p00pg)       &
                            & **pgfaa1
              end do
            end do
          else if ( ipgf.eq.0 ) then
            do k = 1 , kx
              do i = 1 , ixm1
                td(i,k,j) = ta(i,k,j)*(1.+ep1*(qv(i,k,j)))
              end do
            end do
          else
          end if
!
!
!..t..compute weighted p*t (td) for use in ssi:
!
        else if ( ipgf.eq.1 ) then
!
          do k = 1 , kx
            do i = 2 , ixm2
              tvc = tc(i,k,j)*(1.+ep1*(qvc(i,k,j))/psc(i,j))
              tva = ta(i,k,j)*(1.+ep1*(qv(i,k,j)))
              tvb = tb(i,k,j)*(1.+ep1*(qvb(i,k,j))/psb(i,j))
              td(i,k,j) = alpha*(tvc+tvb) + beta*tva
              ttld(i,k,j) = td(i,k,j) - psd(i,j)                        &
                          & *t00pg*((a(k)*psd(i,j)+ptop)/p00pg)**pgfaa1
            end do
          end do
          do k = 1 , kx
            td(1,k,j) = ta(1,k,j)*(1.+ep1*(qv(1,k,j)))
            ttld(1,k,j) = td(1,k,j) - psa(1,j)                          &
                        & *t00pg*((a(k)*psa(1,j)+ptop)/p00pg)**pgfaa1
            td(ixm1,k,j) = ta(ixm1,k,j)*(1.+ep1*(qv(ixm1,k,j)))
            ttld(ixm1,k,j) = td(ixm1,k,j) - psa(ixm1,j)                    &
                          & *t00pg*((a(k)*psa(ixm1,j)+ptop)/p00pg)       &
                          & **pgfaa1
          end do
!
        else if ( ipgf.eq.0 ) then
!
          do k = 1 , kx
            do i = 2 , ixm2
              tvc = tc(i,k,j)*(1.+ep1*(qvc(i,k,j))/psc(i,j))
              tva = ta(i,k,j)*(1.+ep1*(qv(i,k,j)))
              tvb = tb(i,k,j)*(1.+ep1*(qvb(i,k,j))/psb(i,j))
              td(i,k,j) = alpha*(tvc+tvb) + beta*tva
            end do
          end do
          do k = 1 , kx
            td(1,k,j) = ta(1,k,j)*(1.+ep1*(qv(1,k,j)))
            td(ixm1,k,j) = ta(ixm1,k,j)*(1.+ep1*(qv(ixm1,k,j)))
          end do
 
        else         !end if(j.ne.jxm1),else test
        end if
!
      end do
!----------------------------------------------------------------------
!**uv*compute the u and v tendencies:
#ifdef MPP1
      do j = jbegin , jendx
#else
      do j = 2 , jxm1
#endif
!
!..uv.compute the diffusion terms:
!       put diffusion and pbl tendencies of u and v in difuu and difuv.
!
        do k = 1 , kx
          do i = 2 , ixm1
            difuu(i,k,j) = uten(i,k,j)
            difuv(i,k,j) = vten(i,k,j)
          end do
        end do
!
        call diffu_u(difuu(1,1,j),xkc(1,1,j),c203,j,1)
        call diffu_v(difuv(1,1,j),xkc(1,1,j),c203,j,1)
!
!..uv.compute the horizontal advection terms for u and v:
!
        do k = 1 , kx
          do i = 2 , ixm1
            uten(i,k,j) = 0.
            vten(i,k,j) = 0.
          end do
        end do
!
        call hadv_u(uten(1,1,j),dx16,j,3)
        call hadv_v(vten(1,1,j),dx16,j,3)
!
!..uv.compute coriolis terms:
!
        do k = 1 , kx
          do i = 2 , ixm1
            uten(i,k,j) = uten(i,k,j) + f(i,j)*va(i,k,j)/msfd(i,j)
            vten(i,k,j) = vten(i,k,j) - f(i,j)*ua(i,k,j)/msfd(i,j)
          end do
        end do
      end do
!
#ifdef MPP1
      do j = jbegin , jendx
#else
      do j = 2 , jxm1
#endif
!
!..uv.compute pressure gradient terms:
!
        if ( ipgf.eq.1 ) then
          do k = 1 , kx
            do i = 2 , ixm1
              psasum = psd(i,j) + psd(i-1,j) + psd(i,j-1) + psd(i-1,j-1)
              sigpsa = psasum
              tv1 = t(i-1,k,j-1)*(1.+ep1*(qv(i-1,k,j-1)))
              tv2 = t(i,k,j-1)*(1.+ep1*(qv(i,k,j-1)))
              tv3 = t(i-1,k,j)*(1.+ep1*(qv(i-1,k,j)))
              tv4 = t(i,k,j)*(1.+ep1*(qv(i,k,j)))
              rtbar = tv1 + tv2 + tv3 + tv4 -                           &
                    & 4.*t00pg*((a(k)*psasum/4.+ptop)/p00pg)**pgfaa1
              rtbar = r*rtbar*sigpsa/16.
              uten(i,k,j) = uten(i,k,j)                                 &
                          & - rtbar*(dlog(0.5*(psd(i,j)+psd(i-1,j))*a(k)&
                          & +ptop)                                      &
                          & -dlog(0.5*(psd(i,j-1)+psd(i-1,j-1))*a(k)    &
                          & +ptop))/(dx*msfd(i,j))
              vten(i,k,j) = vten(i,k,j)                                 &
                          & - rtbar*(dlog(0.5*(psd(i,j)+psd(i,j-1))*a(k)&
                          & +ptop)                                      &
                          & -dlog(0.5*(psd(i-1,j-1)+psd(i-1,j))*a(k)    &
                          & +ptop))/(dx*msfd(i,j))
            end do
          end do
        else if ( ipgf.eq.0 ) then
          do k = 1 , kx
            do i = 2 , ixm1
              psasum = psd(i,j) + psd(i-1,j) + psd(i,j-1) + psd(i-1,j-1)
              sigpsa = psasum
              tv1 = t(i-1,k,j-1)*(1.+ep1*(qv(i-1,k,j-1)))
              tv2 = t(i,k,j-1)*(1.+ep1*(qv(i,k,j-1)))
              tv3 = t(i-1,k,j)*(1.+ep1*(qv(i-1,k,j)))
              tv4 = t(i,k,j)*(1.+ep1*(qv(i,k,j)))
              rtbar = r*(tv1+tv2+tv3+tv4)*sigpsa/16.
              uten(i,k,j) = uten(i,k,j)                                 &
                          & - rtbar*(dlog(0.5*(psd(i,j)+psd(i-1,j))*a(k)&
                          & +ptop)                                      &
                          & -dlog(0.5*(psd(i,j-1)+psd(i-1,j-1))*a(k)    &
                          & +ptop))/(dx*msfd(i,j))
              vten(i,k,j) = vten(i,k,j)                                 &
                          & - rtbar*(dlog(0.5*(psd(i,j)+psd(i,j-1))*a(k)&
                          & +ptop)                                      &
                          & -dlog(0.5*(psd(i-1,j-1)+psd(i-1,j))*a(k)    &
                          & +ptop))/(dx*msfd(i,j))
            end do
          end do
        else
        end if
      end do
!
#ifdef MPP1
      do j = 1 , jendx
#else
      do j = 1 , jxm1
#endif
!
!..uv.compute geopotential height at half-k levels, cross points:
!
        if ( ipgf.eq.1 ) then
 
          do i = 1 , ixm1
            tv = (ttld(i,kx,j)/psd(i,j))/(1.+qc(i,kx,j)/(1.+qv(i,kx,j)))
            phi(i,kx,j) = ht(i,j)                                       &
                        & + r*t00pg/pgfaa1*((psd(i,j)+ptop)/p00pg)      &
                        & **pgfaa1
            phi(i,kx,j) = phi(i,kx,j)                                   &
                        & - r*tv*dlog((a(kx)+ptop/psd(i,j))/(1.+        &
                        & ptop/psd(i,j)))
          end do
 
          do k = 1 , kxm1
            lev = kx - k
            do i = 1 , ixm1
              tvavg = ((ttld(i,lev,j)*dsigma(lev)+ttld(i,lev+1,j)*dsigma&
                    & (lev+1))/(psd(i,j)*(dsigma(lev)+dsigma(lev+1))))  &
                    & /(1.+qc(i,lev,j)/(1.+qv(i,lev,j)))
              phi(i,lev,j) = phi(i,lev+1,j)                             &
                           & - r*tvavg*dlog((a(lev)+ptop/psd(i,j))      &
                           & /(a(lev+1)+ptop/psd(i,j)))
            end do
          end do
 
        else if ( ipgf.eq.0 ) then
 
          do i = 1 , ixm1
            tv = (td(i,kx,j)/psd(i,j))/(1.+qc(i,kx,j)/(1.+qv(i,kx,j)))
            phi(i,kx,j) = ht(i,j)                                       &
                        & - r*tv*dlog((a(kx)+ptop/psd(i,j))/(1.+ptop/psd&
                        & (i,j)))
          end do
 
          do k = 1 , kxm1
            lev = kx - k
            do i = 1 , ixm1
              tvavg = ((td(i,lev,j)*dsigma(lev)+td(i,lev+1,j)*dsigma(lev&
                    & +1))/(psd(i,j)*(dsigma(lev)+dsigma(lev+1))))      &
                    & /(1.+qc(i,lev,j)/(1.+qv(i,lev,j)))
              phi(i,lev,j) = phi(i,lev+1,j)                             &
                           & - r*tvavg*dlog((a(lev)+ptop/psd(i,j))      &
                           & /(a(lev+1)+ptop/psd(i,j)))
            end do
          end do
 
        else   ! ipgf if block
        end if
      end do
#ifdef MPP1
      call mpi_sendrecv(phi(1,1,jxp),ix*kx,mpi_double_precision,ieast,1,&
                      & phi(1,1,0),ix*kx,mpi_double_precision,iwest,1,  &
                      & mpi_comm_world,mpi_status_ignore,ierr)
      do j = jbegin , jendx
#else
      do j = 2 , jxm1
#endif
!
!..uv.compute the geopotential gradient terms:
!
        do k = 1 , kx
          do i = 2 , ixm1
            uten(i,k,j) = uten(i,k,j)                                   &
                        & - (psd(i-1,j-1)+psd(i,j-1)+psd(i-1,j)+psd(i,j)&
                        & )                                             &
                        & *(phi(i,k,j)+phi(i-1,k,j)-phi(i,k,j-1)-phi(i-1&
                        & ,k,j-1))/(dx8*msfd(i,j))
            vten(i,k,j) = vten(i,k,j)                                   &
                        & - (psd(i-1,j-1)+psd(i,j-1)+psd(i-1,j)+psd(i,j)&
                        & )                                             &
                        & *(phi(i,k,j)+phi(i,k,j-1)-phi(i-1,k,j)-phi(i-1&
                        & ,k,j-1))/(dx8*msfd(i,j))
          end do
        end do
      end do
!
#ifdef MPP1
      do j = jbegin , jendx
#else
      do j = 2 , jxm1
#endif
!
!..uv.compute teh vertical advection terms:
!
        call vadv(uten(1,1,j),ua(1,1,j),j,4)
        call vadv(vten(1,1,j),va(1,1,j),j,4)
!
!..uv.apply the sponge boundary condition on u and v:
!
        if ( iboudy.eq.4 ) then
          call sponge_u(ispgd,wgtd,uten(1,1,j),j)
          call sponge_v(ispgd,wgtd,vten(1,1,j),j)
        end if
!
!..uv.apply the nudging boundary conditions:
!
        if ( iboudy.eq.1 .or. iboudy.eq.5 ) then
          call nudge_u(ispgd,fnudge,gnudge,xtm1,uten(1,1,j),c203,j,     &
                     & iboudy)
          call nudge_v(ispgd,fnudge,gnudge,xtm1,vten(1,1,j),c203,j,     &
                     & iboudy)
        end if
!
!..uv.add the diffusion and pbl tendencies to uten and vten:
!
        do k = 1 , kx
          do i = 2 , ixm1
            uten(i,k,j) = uten(i,k,j) + difuu(i,k,j)
            vten(i,k,j) = vten(i,k,j) + difuv(i,k,j)
          end do
        end do
!
!..uv.forecast p*u and p*v at tau+1:
!
        do k = 1 , kx
          do i = 2 , ixm1
            uc(i,k,j) = ub(i,k,j) + dt*uten(i,k,j)
            vc(i,k,j) = vb(i,k,j) + dt*vten(i,k,j)
          end do
        end do
!
!*****end of j loop.
!**********************************************************************
      end do
!
!---------------------------------------------------------------------
!-----store the xxa variables in xxb and xxc in xxa:
!     perform time smoothing operations.
!
#ifdef MPP1
      do j = jbegin , jendx
#else
      do j = 2 , jxm1
#endif
        do k = 1 , kx
          do i = 2 , ixm1
            ub(i,k,j) = omuhf*ua(i,k,j)/msfd(i,j)                       &
                      & + gnuhf*(ub(i,k,j)+uc(i,k,j))
            vb(i,k,j) = omuhf*va(i,k,j)/msfd(i,j)                       &
                      & + gnuhf*(vb(i,k,j)+vc(i,k,j))
            ua(i,k,j) = uc(i,k,j)
            va(i,k,j) = vc(i,k,j)
          end do
        end do
      end do
!
#ifdef MPP1
      do j = jbegin , jendm
#else
      do j = 2 , jxm2
#endif
        do k = 1 , kx
          do i = 2 , ixm2
            tb(i,k,j) = omuhf*ta(i,k,j) + gnuhf*(tb(i,k,j)+tc(i,k,j))
            ta(i,k,j) = tc(i,k,j)
          end do
          do i = 2 , ixm2
            qvbs = omuhf*qva(i,k,j) + gnuhf*(qvb(i,k,j)+qvc(i,k,j))
            qvas = qvc(i,k,j)
            qvb(i,k,j) = dmax1(qvbs,1.D-99)
            qva(i,k,j) = dmax1(qvas,1.D-99)
          end do
          do i = 2 , ixm2
            qcbs = omu*qca(i,k,j) + gnu*(qcb(i,k,j)+qcc(i,k,j))
            qcb(i,k,j) = dmax1(qcbs,0.D0)
          end do
          do i = 2 , ixm2
            qcas = qcc(i,k,j)
            qca(i,k,j) = dmax1(qcas,0.D0)
          end do
!chem2
          if ( ichem.eq.1 ) then
            do itr = 1 , ntr
              do i = 2 , ixm2
                chibs = omu*chia(i,k,j,itr)                             &
                      & + gnu*(chib(i,k,j,itr)+chic(i,k,j,itr))
                chib(i,k,j,itr) = dmax1(chibs,0.D0)
                chias = chic(i,k,j,itr)
                chia(i,k,j,itr) = dmax1(chias,0.D0)
              end do
            end do
          end if
!chem2_
        end do
        do i = 2 , ixm2
          psb(i,j) = omuhf*psa(i,j) + gnuhf*(psb(i,j)+psc(i,j))
          psa(i,j) = psc(i,j)
        end do
      end do
      if ( ehso4 ) then
        do k = 1 , kx
#ifdef MPP1
          do j = 1 , jendx
#else
          do j = 1 , jxm1
#endif
            do i = 1 , ixm1
              aermm(i,k,j) = so4(i,k,j)
            end do
          end do
        end do
      end if
!
!----------------------------------------------------------------------
!-----increment elapsed forecast time:
!
      ktau = ktau + 1
      xtime = xtime + dtmin
      ntime = ntime + nint(dtmin*60.)
      if ( dabs(xtime-ibdyfrq*60.).lt.0.00001 ) then
        lhour = lhour + ibdyfrq
        if ( lhour.eq.24 ) then
          call finddate(nnnnnn,ldatez)
          ldatez = mdatez(nnnnnn+1)
          lyear = ldatez/1000000
          lmonth = (ldatez-lyear*1000000)/10000
          lday = (ldatez-lyear*1000000-lmonth*10000)/100
          lhour = mod(ldatez,100)
        else
          ldatez = ldatez + ibdyfrq
        end if
        nnnnnn = nnnnnn + 1
        xtime = 0.0
        if ( mod(ldatez,1000000).eq.10100 .and. xtime.lt.0.0001 ) then
          jyear = ldatez/1000000
          ktau = 0
          ntime = 0
        end if
      end if
      if ( jyear.ne.jyear0 .or. ktau.ne.0 ) dt = dt2
!
!-----compute the amounts advected through the lateral boundaries:
!     *** note *** we must calculate the amounts advected through
!     the lateral boundaries before updating the values
!     at boundary slices.
!
#ifdef DIAG
      call conadv
      if ( ichem.eq.1 ) call tracdiag(xkc)
#endif
 
!-----fill up the boundary values for xxb and xxa variables:
!
      call bdyval(xtime,iexec)
!
!-----compute the nonconvective precipitation:
!
!???  call nconvp(psa,psb,ta,tb,qva,qvb,qca,qcb)
!
!chem2_
!     do cumulus transport of tracers
      if ( ichem.eq.1 .and. ichcumtra.eq.1 ) call cumtran
 
!chem2_
 
!-----trace the mass conservation of dry air and water substance:
!
#ifdef DIAG
      call conmas
#endif
!
!
!---- budgets for tracers
      if ( ichem.eq.1 ) call tracbud
!
!-----print out noise parameter:
!
      if ( jyear.ne.jyear0 .or. ktau.gt.1 ) then
        ptnbar = ptntot/dble(iptn)
        pt2bar = pt2tot/dble(iptn)
        icons = 0
#ifdef MPP1
        icons_mpi = 0
        do j = jbegin , jendm
#else
        do j = 2 , jxm2
#endif
          icons = icons + icon(j)
        end do
#ifdef MPP1
        icons_mpi = 0
        call mpi_allreduce(icons,icons_mpi,1,mpi_integer,mpi_sum,       &
                         & mpi_comm_world,ierr)
#endif
        xday = ((nnnnnn-nstrt0)*ibdyfrq*60.+xtime-dtmin)/1440.
#ifdef MPP1
        if ( myid.eq.0 ) then
          if ( mod(ktau,50).eq.0 ) print 99001 , xday , ktau , ptnbar , &
             & pt2bar , icons_mpi
        end if
#else
        if ( mod(ktau,50).eq.0 ) print 99001 , xday , ktau , ptnbar ,   &
                                     & pt2bar , icons
#endif
99001     format (5x,'at day = ',f9.4,', ktau = ',i10,                  &
                 &' :  1st, 2nd time deriv of ps = ',2E12.5,            &
                 &',  no. of points w/convection = ',i7)
      end if
!
!----------------------------------------------------------------------
!
!-----recalculate solar declination angle if forecast time larger than
!     24 hours:
!
      if ( dabs(xtime).lt.0.00001 .and. ldatez.ne.idate1 ) then
        call solar1(xtime)
        dectim = anint(1440.+dectim)
#ifdef MPP1
        if ( myid.eq.0 ) write (*,*) ' dectim = ' , dectim
#else
        write (*,*) ' dectim = ' , dectim
#endif
      end if
!
      end subroutine tend
