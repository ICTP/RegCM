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

      module mod_bdycod
!
! Storage and subroutines for input of boundary values and tendencies
! of p*u, p*v, p*t, p*qv and p*, and outmost 2 slices of u and v for
! large domain.
!
      use mod_dynparam
      use mod_runparams
      use mod_main
      use mod_mainchem
      use mod_bats
      use mod_message 
      use mod_ncio
      use mod_date
      use mod_cvaria
      use mod_service
#ifdef MPP1
      use mod_mppio
#endif
!
      private
!
      public :: allocate_mod_bdycon , bdyin , bdyval
#ifndef BAND
      public :: uj1 , uj2 , ujlx , ujl
      public :: vj1 , vj2 , vjlx , vjl
#endif
      public :: ui1 , ui2 , uilx , uil
      public :: vi1 , vi2 , vilx , vil
      public :: ub0 , vb0 , qb0 , tb0 , ps0 , ts0 , so0
#ifndef BAND
      public :: peb , pebt , pwb , pwbt
      public :: ueb , uebt , veb , vebt , uwb , uwbt , vwb , vwbt
#endif
      public :: pnb , pnbt , pss , psbt
      public :: unb , unbt , vnb , vnbt , usb , usbt , vsb , vsbt
#ifndef BAND
      public :: teb , tebt , qeb , qebt , twb , twbt , qwb , qwbt
#endif
      public :: tnb , tnbt , qnb , qnbt , tsb , tsbt , qsb , qsbt
      public :: ts1 ! FOR DCSST
!
#ifndef BAND
      real(8) , allocatable , dimension(:,:) :: uj1 , uj2 , ujl , ujlx ,&
                      &  vj1 , vj2 , vjl , vjlx
#endif
      real(8) , allocatable , dimension(:,:) :: ps0 , ps1
      real(8) , allocatable , dimension(:,:,:) :: qb0 , qb1 , so0 ,     &
                      & so1 , tb0 , tb1 , ub0 , ub1 , vb0 , vb1
      real(8) , allocatable , dimension(:,:) :: ts0 , ts1
!
#ifndef BAND
      real(8) , allocatable , dimension(:,:) :: peb , pebt , pwb , pwbt
#endif
      real(8) , allocatable , dimension(:,:) :: pnb , pnbt , psbt , pss
#ifndef BAND
      real(8) , allocatable , dimension(:,:,:) :: qeb , qebt , qwb ,    &
           & qwbt , teb , tebt , twb , twbt , ueb , uebt , uwb , uwbt , &
           & veb , vebt , vwb , vwbt
#endif
      real(8) , allocatable , dimension(:,:,:) :: qnb , qnbt , qsb ,    &
           & qsbt , tnb , tnbt , tsb , tsbt
      real(8) , allocatable , dimension(:,:) :: ui1 , ui2 , uil , uilx ,&
           & vi1 , vi2 , vil , vilx
      real(8) , allocatable , dimension(:,:,:) :: unb , unbt , usb ,    &
           & usbt , vnb , vnbt , vsb , vsbt
!
      contains
!
        subroutine allocate_mod_bdycon
        implicit none
        character (len=50) :: subroutine_name='allocate_mod_bdycon'
        integer :: idindx=0
!
        call time_begin(subroutine_name,idindx)
#ifdef MPP1
        allocate(ps0(iy,0:jxp+1))
        allocate(ps1(iy,0:jxp+1))
        allocate(ts0(iy,jxp))
        allocate(ts1(iy,jxp))
!
        allocate(qb0(iy,kz,jxp))
        allocate(qb1(iy,kz,jxp))
        allocate(so0(iy,kz,jxp))
        allocate(so1(iy,kz,jxp))
        allocate(tb0(iy,kz,jxp))
        allocate(tb1(iy,kz,jxp))
        allocate(ub0(iy,kz,jxp))
        allocate(ub1(iy,kz,jxp))
        allocate(vb0(iy,kz,jxp))
        allocate(vb1(iy,kz,jxp))
#ifndef BAND
        allocate(peb(iy,0:jxp+1))
        allocate(pebt(iy,0:jxp+1))
        allocate(pwb(iy,0:jxp+1))
        allocate(pwbt(iy,0:jxp+1))
#endif
        allocate(pnb(nspgx,0:jxp+1))
        allocate(pnbt(nspgx,0:jxp+1))
        allocate(psbt(nspgx,0:jxp+1))
        allocate(pss(nspgx,0:jxp+1))
#ifndef BAND
        allocate(qeb(iy,kz,0:jxp+1))
        allocate(qebt(iy,kz,0:jxp+1))
        allocate(qwb(iy,kz,0:jxp+1))
        allocate(qwbt(iy,kz,0:jxp+1))
        allocate(teb(iy,kz,0:jxp+1))
        allocate(tebt(iy,kz,0:jxp+1))
        allocate(twb(iy,kz,0:jxp+1))
        allocate(twbt(iy,kz,0:jxp+1))
        allocate(ueb(iy,kz,0:jxp+1))
        allocate(uebt(iy,kz,0:jxp+1))
        allocate(uwb(iy,kz,0:jxp+1))
        allocate(uwbt(iy,kz,0:jxp+1))
        allocate(veb(iy,kz,0:jxp+1))
        allocate(vebt(iy,kz,0:jxp+1))
        allocate(vwb(iy,kz,0:jxp+1))
        allocate(vwbt(iy,kz,0:jxp+1))
#endif
        allocate(qnb(nspgx,kz,0:jxp+1))
        allocate(qnbt(nspgx,kz,0:jxp+1))
        allocate(qsb(nspgx,kz,0:jxp+1))
        allocate(qsbt(nspgx,kz,0:jxp+1))
        allocate(tnb(nspgx,kz,0:jxp+1))
        allocate(tnbt(nspgx,kz,0:jxp+1))
        allocate(tsb(nspgx,kz,0:jxp+1))
        allocate(tsbt(nspgx,kz,0:jxp+1))
!
        allocate(ui1(kz,0:jxp+1))
        allocate(ui2(kz,0:jxp+1))
        allocate(uil(kz,0:jxp+1))
        allocate(uilx(kz,0:jxp+1))
        allocate(vi1(kz,0:jxp+1))
        allocate(vi2(kz,0:jxp+1))
        allocate(vil(kz,0:jxp+1))
        allocate(vilx(kz,0:jxp+1))
!
        allocate(unb(nspgd,kz,0:jxp+1))
        allocate(unbt(nspgd,kz,0:jxp+1))
        allocate(usb(nspgd,kz,0:jxp+1))
        allocate(usbt(nspgd,kz,0:jxp+1))
        allocate(vnb(nspgd,kz,0:jxp+1))
        allocate(vnbt(nspgd,kz,0:jxp+1))
        allocate(vsb(nspgd,kz,0:jxp+1))
        allocate(vsbt(nspgd,kz,0:jxp+1))
#else
        allocate(ps0(iy,jx))
        allocate(ps1(iy,jx))
        allocate(ts0(iy,jx))
        allocate(ts1(iy,jx))
        allocate(qb0(iy,kz,jx))
        allocate(qb1(iy,kz,jx))
        allocate(so0(iy,kz,jx))
        allocate(so1(iy,kz,jx))
        allocate(tb0(iy,kz,jx))
        allocate(tb1(iy,kz,jx))
        allocate(ub0(iy,kz,jx))
        allocate(ub1(iy,kz,jx))
        allocate(vb0(iy,kz,jx))
        allocate(vb1(iy,kz,jx))
#ifndef BAND
        allocate(peb(iy,nspgx))
        allocate(pebt(iy,nspgx))
        allocate(pwb(iy,nspgx))
        allocate(pwbt(iy,nspgx))
#endif
        allocate(pnb(nspgx,jx))
        allocate(pnbt(nspgx,jx))
        allocate(psbt(nspgx,jx))
        allocate(pss(nspgx,jx))
#ifndef BAND
        allocate(qeb(iy,kz,nspgx))
        allocate(qebt(iy,kz,nspgx))
        allocate(qwb(iy,kz,nspgx))
        allocate(qwbt(iy,kz,nspgx))
        allocate(teb(iy,kz,nspgx))
        allocate(tebt(iy,kz,nspgx))
        allocate(twb(iy,kz,nspgx))
        allocate(twbt(iy,kz,nspgx))
#endif
        allocate(qnb(nspgx,kz,jx))
        allocate(qnbt(nspgx,kz,jx))
        allocate(qsb(nspgx,kz,jx))
        allocate(qsbt(nspgx,kz,jx))
        allocate(tnb(nspgx,kz,jx))
        allocate(tnbt(nspgx,kz,jx))
        allocate(tsb(nspgx,kz,jx))
        allocate(tsbt(nspgx,kz,jx))
#ifndef BAND
        allocate(ueb(iy,kz,nspgd))
        allocate(uebt(iy,kz,nspgd))
        allocate(uwb(iy,kz,nspgd))
        allocate(uwbt(iy,kz,nspgd))
        allocate(veb(iy,kz,nspgd))
        allocate(vebt(iy,kz,nspgd))
        allocate(vwb(iy,kz,nspgd))
        allocate(vwbt(iy,kz,nspgd))
#endif
        allocate(ui1(kz,jx))
        allocate(ui2(kz,jx))
        allocate(uil(kz,jx))
        allocate(uilx(kz,jx))
        allocate(vi1(kz,jx))
        allocate(vi2(kz,jx))
        allocate(vil(kz,jx))
        allocate(vilx(kz,jx))
        allocate(unb(nspgd,kz,jx))
        allocate(unbt(nspgd,kz,jx))
        allocate(usb(nspgd,kz,jx))
        allocate(usbt(nspgd,kz,jx))
        allocate(vnb(nspgd,kz,jx))
        allocate(vnbt(nspgd,kz,jx))
        allocate(vsb(nspgd,kz,jx))
        allocate(vsbt(nspgd,kz,jx))
#endif 
#ifndef BAND
        allocate(uj1(iy,kz))
        allocate(uj2(iy,kz))
        allocate(ujl(iy,kz))
        allocate(ujlx(iy,kz))
        allocate(vj1(iy,kz))
        allocate(vj2(iy,kz))
        allocate(vjl(iy,kz))
        allocate(vjlx(iy,kz))
#endif
        ps0 = 0.0D0
        ps1 = 0.0D0
        ts0 = 0.0D0
        ts1 = 0.0D0
        qb0 = 0.0D0
        qb1 = 0.0D0
        so0 = 0.0D0
        so1 = 0.0D0
        tb0 = 0.0D0
        tb1 = 0.0D0
        ub0 = 0.0D0
        ub1 = 0.0D0
        vb0 = 0.0D0
        vb1 = 0.0D0
#ifndef BAND
        peb = 0.0D0
        pebt = 0.0D0
        pwb = 0.0D0
        pwbt = 0.0D0
#endif
        pnb = 0.0D0
        pnbt = 0.0D0
        psbt = 0.0D0
        pss = 0.0D0
#ifndef BAND
        qeb = 0.0D0
        qebt = 0.0D0
        qwb = 0.0D0
        qwbt = 0.0D0
        teb = 0.0D0
        tebt = 0.0D0
        twb = 0.0D0
        twbt = 0.0D0
#endif
        qnb = 0.0D0
        qnbt = 0.0D0
        qsb = 0.0D0
        qsbt = 0.0D0
        tnb = 0.0D0
        tnbt = 0.0D0
        tsb = 0.0D0
        tsbt = 0.0D0
#ifndef BAND
        ueb = 0.0D0
        uebt = 0.0D0
        uwb = 0.0D0
        uwbt = 0.0D0
        veb = 0.0D0
        vebt = 0.0D0
        vwb = 0.0D0
        vwbt = 0.0D0
#endif
        ui1 = 0.0D0
        ui2 = 0.0D0
        uil = 0.0D0
        uilx = 0.0D0
        vi1 = 0.0D0
        vi2 = 0.0D0
        vil = 0.0D0
        vilx = 0.0D0
        unb = 0.0D0
        unbt = 0.0D0
        usb = 0.0D0
        usbt = 0.0D0
        vnb = 0.0D0
        vnbt = 0.0D0
        vsb = 0.0D0
        vsbt = 0.0D0
#ifndef BAND
        uj1 = 0.0D0
        uj2 = 0.0D0
        ujl = 0.0D0
        ujlx = 0.0D0
        vj1 = 0.0D0
        vj2 = 0.0D0
        vjl = 0.0D0
        vjlx = 0.0D0
#endif
        call time_end(subroutine_name,idindx)
        end subroutine allocate_mod_bdycon
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine reads in the boundary conditions.               c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine bdyin
!
#ifdef MPP1
      use mod_mppio
#ifndef IBM
      use mpi
#else
      include 'mpif.h'
#endif
#endif
      implicit none
!
      real(8) :: dtbdys
      integer :: i , j , k , nn , nnb , mmrec
#ifdef MPP1
      integer :: ierr , ndeb , ndwb , nxeb , nxwb
#ifndef BAND
      integer :: nkk
#endif
      real(8) , dimension(iy,jxp) :: psdot , tdum
#else
      real(8) , dimension(iy,jx) :: psdot , tdum
#endif
      integer :: n
      character (len=50) :: subroutine_name='bdyin'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
      if ( dabs(xtime).gt.0.0001 ) return
!
#ifdef MPP1
      dtbdys = ibdyfrq*60.*60.
      if ( myid.eq.0 ) then
        if ( ehso4 ) then
          do k = 1 , kz
            do j = 1 , jendl
              do i = 1 , iy
                sulf%so4(i,k,j) = so0(i,k,j)
              end do
            end do
          end do
        end if
        call addhours(ndate1, ibdyfrq)
        mmrec = icbc_search(ndate1)
        if (mmrec < 0) then
          call open_icbc(imonfirst(ndate1))
        end if
        call read_icbc(ndate1,ps1_io,ts1_io,ub1_io,vb1_io, &
                   &     tb1_io,qb1_io,so1_io)
        ps1_io = ps1_io/10.0
        do j = 1 , jx
          do k = 1 , kz
            do i = 1 , iy
              sav_0(i,k,j) = ub1_io(i,k,j)
              sav_0(i,kz+k,j) = vb1_io(i,k,j)
              sav_0(i,kz*2+k,j) = qb1_io(i,k,j)
              sav_0(i,kz*3+k,j) = tb1_io(i,k,j)
            end do
          end do
          do i = 1 , iy
            sav_0(i,kz*4+1,j) = ps1_io(i,j)
            sav_0(i,kz*4+2,j) = ts1_io(i,j)
          end do
        end do
        if ( ehso4 ) then
          do j = 1 , jx
            do k = 1 , kz
              do i = 1 , iy
                sav_0s(i,k,j) = so1_io(i,k,j)
              end do
            end do
          end do
        end if
      end if
!
      call mpi_bcast(ndate1,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_scatter(sav_0,iy*(kz*4+2)*jxp,mpi_real8,        &
                     & sav0, iy*(kz*4+2)*jxp,mpi_real8,        &
                     & 0,mpi_comm_world,ierr)
      do j = 1 , jendl
        do k = 1 , kz
          do i = 1 , iy
            ub1(i,k,j) = sav0(i,k,j)
            vb1(i,k,j) = sav0(i,kz+k,j)
            qb1(i,k,j) = sav0(i,kz*2+k,j)
            tb1(i,k,j) = sav0(i,kz*3+k,j)
          end do
        end do
        do i = 1 , iy
          ps1(i,j) = sav0(i,kz*4+1,j)
          ts1(i,j) = sav0(i,kz*4+2,j)
        end do
      end do
      if ( ehso4 ) then
        call mpi_scatter(sav_0s,iy*kz*jxp,mpi_real8,           &
                       & sav0s,iy*kz*jxp,mpi_real8,            &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendl
          do k = 1 , kz
            do i = 1 , iy
              so1(i,k,j) = sav0s(i,k,j)
            end do
          end do
        end do
      end if
!     Convert surface pressure to pstar
      do j = 1 , jendl
        do i = 1 , iy
          ps1(i,j) = ps1(i,j) - r8pt
        end do
      end do
!=======================================================================
!
!     this routine determines p(.) from p(x) by a 4-point
!     interpolation. on the x-grid, a p(x) point outside the grid
!     domain is assumed to satisfy p(0,j)=p(1,j); p(iy,j)=p(iym1,j);
!     and similarly for the i's.
!
      call mpi_sendrecv(ps1(1,jxp),iy,mpi_real8,ieast,1,              &
                      & ps1(1,0),iy,mpi_real8,iwest,1,                &
                      & mpi_comm_world,mpi_status_ignore,ierr)
      do j = jbegin , jendx
        do i = 2 , iym1
          psdot(i,j) = 0.25*(ps1(i,j)+ps1(i-1,j) +                    &
                     &       ps1(i,j-1)+ps1(i-1,j-1))
        end do
      end do
#ifdef BAND
      do j = jbegin , jendx
        psdot(1,j) = 0.5*(ps1(1,j)+ps1(1,j-1))
        psdot(iy,j) = 0.5*(ps1(iym1,j)+ps1(iym1,j-1))
      end do
#else
!
      do i = 2 , iym1
        if ( myid.eq.0 ) psdot(i,1) = 0.5*(ps1(i,1)+ps1(i-1,1))
        if ( myid.eq.nproc-1 )                                        &
           & psdot(i,jendl) = 0.5*(ps1(i,jendx)+ps1(i-1,jendx))
      end do
!
      do j = jbegin , jendx
        psdot(1,j) = 0.5*(ps1(1,j)+ps1(1,j-1))
        psdot(iy,j) = 0.5*(ps1(iym1,j)+ps1(iym1,j-1))
      end do
!
      if ( myid.eq.0 ) then
        psdot(1,1) = ps1(1,1)
        psdot(iy,1) = ps1(iym1,1)
      end if
      if ( myid.eq.nproc-1 ) then
        psdot(1,jendl) = ps1(1,jendx)
        psdot(iy,jendl) = ps1(iym1,jendx)
      end if
!
#endif
!=======================================================================
!       Couple pressure u,v,t,q
      do k = 1 , kz
        do j = 1 , jendl
          do i = 1 , iy
            ub1(i,k,j) = ub1(i,k,j)*psdot(i,j)
            vb1(i,k,j) = vb1(i,k,j)*psdot(i,j)
            tb1(i,k,j) = tb1(i,k,j)*ps1(i,j)
            qb1(i,k,j) = qb1(i,k,j)*ps1(i,j)
          end do
        end do
      end do
 
      mdate = ndate0
      nnnchk = nnnchk + 1
 
!     print*,'read in datasets at :',ndate1
!
!-----compute boundary conditions for p*:
!
#ifdef BAND
      nxwb=0
      nxeb=0
#else
      if ( nspgx.le.jxp ) then
        nxwb = nspgx
      else
        nkk = nspgx/jxp
        if ( nspgx.eq.nkk*jxp ) then
          nxwb = jxp
        else
          nxwb = nspgx - nkk*jxp
        end if
      end if
      if ( nxwb+myid*jxp.gt.nspgx ) then
        nxwb = 0
      else if ( nxwb+myid*jxp.lt.nspgx ) then
        nxwb = jxp
      end if
!
      if ( nspgx.le.jxp-1 ) then
        nxeb = nspgx
      else
        nkk = (nspgx-jxp+1)/jxp
        if ( (nspgx-jxp+1).eq.nkk*jxp ) then
          nxeb = jxp
        else
          nxeb = (nspgx-jxp+1) - nkk*jxp
        end if
      end if
      if ( jxm1-(myid*jxp+jxp-nxeb).gt.nspgx ) then
        nxeb = 0
      else if ( jxm1-(myid*jxp+jxp-nxeb).lt.nspgx ) then
        nxeb = min(jendx,jxp)
      end if
      do nn = 1 , nxwb
        do i = 1 , iym1
          pwb(i,nn) = ps0(i,nn)
          pwbt(i,nn) = (ps1(i,nn)-ps0(i,nn))/dtbdys
        end do
      end do
      do nn = 1 , nxeb
        nnb = min(jendx,jxp) - nn + 1
        do i = 1 , iym1
          peb(i,nn) = ps0(i,nnb)
          pebt(i,nn) = (ps1(i,nnb)-ps0(i,nnb))/dtbdys
        end do
      end do
#endif
      do nn = 1 , nspgx
        nnb = iym1 - nn + 1
        do j = 1 , jendx
          pnb(nn,j) = ps0(nnb,j)
          pss(nn,j) = ps0(nn,j)
          pnbt(nn,j) = (ps1(nnb,j)-ps0(nnb,j))/dtbdys
          psbt(nn,j) = (ps1(nn,j)-ps0(nn,j))/dtbdys
        end do
      end do
!
!-----compute boundary conditions for p*u and p*v:
!
#ifdef BAND
      ndwb = 0
      ndeb = 0
#else
      if ( nspgd.le.jxp ) then
        ndwb = nspgd
      else
        nkk = nspgd/jxp
        if ( nspgd.eq.nkk*jxp ) then
          ndwb = jxp
        else
          ndwb = nspgd - nkk*jxp
        end if
      end if
      if ( ndwb+myid*jxp.gt.nspgd ) then
        ndwb = 0
      else if ( ndwb+myid*jxp.lt.nspgd ) then
        ndwb = jxp
      end if
!
      if ( nspgd.le.jendl ) then
        ndeb = nspgd
      else
        nkk = nspgd/jxp
        if ( nspgd.eq.nkk*jxp ) then
          ndeb = jxp
        else
          ndeb = nspgd - nkk*jxp
        end if
      end if
      if ( jx-(myid*jxp+jxp-ndeb).gt.nspgd ) then
        ndeb = 0
      else if ( jx-(myid*jxp+jxp-ndeb).lt.nspgd ) then
        ndeb = jxp
      end if
      do nn = 1 , ndwb
        do k = 1 , kz
          do i = 1 , iy
            uwb(i,k,nn) = ub0(i,k,nn)
            vwb(i,k,nn) = vb0(i,k,nn)
            uwbt(i,k,nn) = (ub1(i,k,nn)-ub0(i,k,nn))/dtbdys
            vwbt(i,k,nn) = (vb1(i,k,nn)-vb0(i,k,nn))/dtbdys
          end do
        end do
      end do
      do nn = 1 , ndeb
        nnb = min(jendl,jxp) - nn + 1
        do k = 1 , kz
          do i = 1 , iy
            ueb(i,k,nn) = ub0(i,k,nnb)
            veb(i,k,nn) = vb0(i,k,nnb)
            uebt(i,k,nn) = (ub1(i,k,nnb)-ub0(i,k,nnb))/dtbdys
            vebt(i,k,nn) = (vb1(i,k,nnb)-vb0(i,k,nnb))/dtbdys
          end do
        end do
      end do
#endif
      do nn = 1 , nspgd
        nnb = iy - nn + 1
        do k = 1 , kz
          do j = 1 , jendl
            unb(nn,k,j) = ub0(nnb,k,j)
            usb(nn,k,j) = ub0(nn,k,j)
            vnb(nn,k,j) = vb0(nnb,k,j)
            vsb(nn,k,j) = vb0(nn,k,j)
            unbt(nn,k,j) = (ub1(nnb,k,j)-ub0(nnb,k,j))/dtbdys
            usbt(nn,k,j) = (ub1(nn,k,j)-ub0(nn,k,j))/dtbdys
            vnbt(nn,k,j) = (vb1(nnb,k,j)-vb0(nnb,k,j))/dtbdys
            vsbt(nn,k,j) = (vb1(nn,k,j)-vb0(nn,k,j))/dtbdys
          end do
        end do
      end do
!
!-----compute boundary conditions for p*t and p*qv:
!
#ifndef BAND
      do nn = 1 , nxwb
        do k = 1 , kz
          do i = 1 , iym1
            twb(i,k,nn) = tb0(i,k,nn)
            qwb(i,k,nn) = qb0(i,k,nn)
            twbt(i,k,nn) = (tb1(i,k,nn)-tb0(i,k,nn))/dtbdys
            qwbt(i,k,nn) = (qb1(i,k,nn)-qb0(i,k,nn))/dtbdys
          end do
        end do
      end do
      do nn = 1 , nxeb
        nnb = min(jendx,jxp) - nn + 1
        do k = 1 , kz
          do i = 1 , iym1
            teb(i,k,nn) = tb0(i,k,nnb)
            qeb(i,k,nn) = qb0(i,k,nnb)
            tebt(i,k,nn) = (tb1(i,k,nnb)-tb0(i,k,nnb))/dtbdys
            qebt(i,k,nn) = (qb1(i,k,nnb)-qb0(i,k,nnb))/dtbdys
          end do
        end do
      end do
#endif
      do nn = 1 , nspgx
        nnb = iym1 - nn + 1
        do k = 1 , kz
          do j = 1 , jendx
            tnb(nn,k,j) = tb0(nnb,k,j)
            tsb(nn,k,j) = tb0(nn,k,j)
            qnb(nn,k,j) = qb0(nnb,k,j)
            qsb(nn,k,j) = qb0(nn,k,j)
            tnbt(nn,k,j) = (tb1(nnb,k,j)-tb0(nnb,k,j))/dtbdys
            tsbt(nn,k,j) = (tb1(nn,k,j)-tb0(nn,k,j))/dtbdys
            qnbt(nn,k,j) = (qb1(nnb,k,j)-qb0(nnb,k,j))/dtbdys
            qsbt(nn,k,j) = (qb1(nn,k,j)-qb0(nn,k,j))/dtbdys
          end do
        end do
      end do
      if ( myid.eq.0 ) print * , 'BCs are ready from ' , ndate0 ,     &
                            &'  to ' , ndate1
      idatex = ndate0
      ndate0 = ndate1
      do j = 1 , jendx
        do i = 1 , iym1
          tdum(i,j) = ts1(i,j)
        end do
      end do
      do k = 1 , kz
        do j = 1 , jendl
          do i = 1 , iy
            ub0(i,k,j) = ub1(i,k,j)
            vb0(i,k,j) = vb1(i,k,j)
            qb0(i,k,j) = qb1(i,k,j)
            tb0(i,k,j) = tb1(i,k,j)
          end do
        end do
      end do
      do j = 1 , jendl
        do i = 1 , iy
          ps0(i,j) = ps1(i,j)
          ts0(i,j) = ts1(i,j)
        end do
      end do
      if ( ehso4 ) then
        do k = 1 , kz
          do j = 1 , jendl
            do i = 1 , iy
              so0(i,k,j) = so1(i,k,j)
            end do
          end do
        end do
      end if
!bxqOCT2001_
 
      nnbase = nnnnnn
      call split_idate(mdate, nyear, nmonth, nday, nhour)

!-----------------------------------------------------------------------
      if ( ldatez.lt.ndate1 ) then
 
        do j = 1 , jendx
          do i = 1 , iym1
#ifdef CLM
! manuaully setting ocld2d subgrid to 1 (regcm_clm does not support subgridding)
            if ( ocld2d(1,i,j).le.0.5 .or. ocld2d(1,i,j).gt.1.5 ) then
#else
            if ( veg2d(i,j).le.0.00001 ) then
#endif
              if (idcsst == 1) then
                sts1%tg(i,j) = tdum(i,j) + dtskin(i,j)
                sts2%tg(i,j) = tdum(i,j) + dtskin(i,j)
              else
                sts1%tg(i,j) = tdum(i,j)
                sts2%tg(i,j) = tdum(i,j)
              end if
              if (iseaice == 1) then
                if ( tdum(i,j).le.271.38 ) then
                   sts1%tg(i,j) = 271.38
                   sts2%tg(i,j) = 271.38
                   tdum(i,j) = 271.38
                  do n = 1, nnsg
                    ocld2d(n,i,j) = 2.0
                    sice2d(n,i,j) = 1000.0
                    scv2d(n,i,j) = 0.0
                  end do
                else
                  sts1%tg(i,j) = tdum(i,j)
                  sts2%tg(i,j) = tdum(i,j)
                  do n = 1, nnsg
                    ocld2d(n,i,j) = 0.0
                    sice2d(n,i,j) = 0.0
                    scv2d(n,i,j) = 0.0
                  end do
                end if
              end if
            end if
          end do
        end do
      end if
#else
      dtbdys = ibdyfrq*60.*60.
      if ( ehso4 ) then
        do k = 1 , kz
          do j = 1 , jx
            do i = 1 , iy
              sulf%so4(i,k,j) = so0(i,k,j)
            end do
          end do
        end do
      end if
      call addhours(ndate1, ibdyfrq)
      mmrec = icbc_search(ndate1)
      if (mmrec < 0) then
        call open_icbc(imonfirst(ndate1))
      end if
      call read_icbc(ndate1,ps1,ts1,ub1,vb1,tb1,qb1,so1)

!     Convert surface pressure to pstar
      do j = 1 , jx
        do i = 1 , iy
          ps1(i,j) = ps1(i,j)/10.0 - r8pt
        end do
      end do
!=====================================================================
!
!   this routine determines p(.) from p(x) by a 4-point
!   interpolation. on the x-grid, a p(x) point outside the grid
!   domain is assumed to satisfy p(0,j)=p(1,j); p(iy,j)=p(iym1,j);
!   and similarly for the i's.

#ifdef BAND
      do j = 2 , jx
        do i = 2 , iym1
          psdot(i,j) = 0.25*(ps1(i,j)+ps1(i-1,j) +               &
                     &       ps1(i,j-1)+ps1(i-1,j-1))
        end do
      end do
!
      do i = 2 , iym1
        psdot(i,1) = 0.25*(ps1(i,1)+ps1(i-1,1) +                 &
                   &       ps1(i,jx)+ps1(i-1,jx))
      end do
!
      do j = 2 , jx
        psdot(1,j) = 0.5*(ps1(1,j)+ps1(1,j-1))
        psdot(iy,j) = 0.5*(ps1(iym1,j)+ps1(iym1,j-1))
      end do
!
      psdot(1,1) = 0.5*(ps1(1,1)+ps1(1,jx))
      psdot(iy,1) = 0.5*(ps1(iym1,1)+ps1(iym1,jx))
!
#else
      do j = 2 , jxm1
        do i = 2 , iym1
          psdot(i,j) = 0.25*(ps1(i,j)+ps1(i-1,j) +                    &
                     &       ps1(i,j-1)+ps1(i-1,j-1))
        end do
      end do
!
      do i = 2 , iym1
        psdot(i,1) = 0.5*(ps1(i,1)+ps1(i-1,1))
        psdot(i,jx) = 0.5*(ps1(i,jxm1)+ps1(i-1,jxm1))
      end do
!
      do j = 2 , jxm1
        psdot(1,j) = 0.5*(ps1(1,j)+ps1(1,j-1))
        psdot(iy,j) = 0.5*(ps1(iym1,j)+ps1(iym1,j-1))
      end do
!
      psdot(1,1) = ps1(1,1)
      psdot(iy,1) = ps1(iym1,1)
      psdot(1,jx) = ps1(1,jxm1)
      psdot(iy,jx) = ps1(iym1,jxm1)
!
#endif
!=======================================================================
!     Couple pressure u,v,t,q
      do k = 1 , kz
        do j = 1 , jx
          do i = 1 , iy
            ub1(i,k,j) = ub1(i,k,j)*psdot(i,j)
            vb1(i,k,j) = vb1(i,k,j)*psdot(i,j)
            tb1(i,k,j) = tb1(i,k,j)*ps1(i,j)
            qb1(i,k,j) = qb1(i,k,j)*ps1(i,j)
          end do
        end do
      end do
 
      mdate = ndate0
      nnnchk = nnnchk + 1
 
!     print*,'read in datasets at :',ndate1
!
!-----compute boundary conditions for p*:
!
 
#ifndef BAND
      do nn = 1 , nspgx
        do i = 1 , iym1
          pwb(i,nn) = ps0(i,nn)
          pwbt(i,nn) = (ps1(i,nn)-ps0(i,nn))/dtbdys
        end do
      end do
      do nn = 1 , nspgx
        nnb = jxm1 - nn + 1
        do i = 1 , iym1
          peb(i,nn) = ps0(i,nnb)
          pebt(i,nn) = (ps1(i,nnb)-ps0(i,nnb))/dtbdys
        end do
      end do
#endif
      do nn = 1 , nspgx
        nnb = iym1 - nn + 1
#ifdef BAND
        do j = 1 , jx
#else
        do j = 1 , jxm1
#endif
          pnb(nn,j) = ps0(nnb,j)
          pss(nn,j) = ps0(nn,j)
          pnbt(nn,j) = (ps1(nnb,j)-ps0(nnb,j))/dtbdys
          psbt(nn,j) = (ps1(nn,j)-ps0(nn,j))/dtbdys
        end do
      end do
!
!-----compute boundary conditions for p*u and p*v:
!
#ifdef BAND
#else
      do nn = 1 , nspgd
        do k = 1 , kz
          do i = 1 , iy
            uwb(i,k,nn) = ub0(i,k,nn)
            vwb(i,k,nn) = vb0(i,k,nn)
            uwbt(i,k,nn) = (ub1(i,k,nn)-ub0(i,k,nn))/dtbdys
            vwbt(i,k,nn) = (vb1(i,k,nn)-vb0(i,k,nn))/dtbdys
          end do
        end do
      end do
      do nn = 1 , nspgd
        nnb = jx - nn + 1
        do k = 1 , kz
          do i = 1 , iy
            ueb(i,k,nn) = ub0(i,k,nnb)
            veb(i,k,nn) = vb0(i,k,nnb)
            uebt(i,k,nn) = (ub1(i,k,nnb)-ub0(i,k,nnb))/dtbdys
            vebt(i,k,nn) = (vb1(i,k,nnb)-vb0(i,k,nnb))/dtbdys
          end do
        end do
      end do
#endif
      do nn = 1 , nspgd
        nnb = iy - nn + 1
        do k = 1 , kz
          do j = 1 , jx
            unb(nn,k,j) = ub0(nnb,k,j)
            usb(nn,k,j) = ub0(nn,k,j)
            vnb(nn,k,j) = vb0(nnb,k,j)
            vsb(nn,k,j) = vb0(nn,k,j)
            unbt(nn,k,j) = (ub1(nnb,k,j)-ub0(nnb,k,j))/dtbdys
            usbt(nn,k,j) = (ub1(nn,k,j)-ub0(nn,k,j))/dtbdys
            vnbt(nn,k,j) = (vb1(nnb,k,j)-vb0(nnb,k,j))/dtbdys
            vsbt(nn,k,j) = (vb1(nn,k,j)-vb0(nn,k,j))/dtbdys
          end do
        end do
      end do
!
!-----compute boundary conditions for p*t and p*qv:
!
#ifndef BAND
      do nn = 1 , nspgx
        do k = 1 , kz
          do i = 1 , iym1
            twb(i,k,nn) = tb0(i,k,nn)
            qwb(i,k,nn) = qb0(i,k,nn)
            twbt(i,k,nn) = (tb1(i,k,nn)-tb0(i,k,nn))/dtbdys
            qwbt(i,k,nn) = (qb1(i,k,nn)-qb0(i,k,nn))/dtbdys
          end do
        end do
      end do
      do nn = 1 , nspgx
        nnb = jxm1 - nn + 1
        do k = 1 , kz
          do i = 1 , iym1
            teb(i,k,nn) = tb0(i,k,nnb)
            qeb(i,k,nn) = qb0(i,k,nnb)
            tebt(i,k,nn) = (tb1(i,k,nnb)-tb0(i,k,nnb))/dtbdys
            qebt(i,k,nn) = (qb1(i,k,nnb)-qb0(i,k,nnb))/dtbdys
          end do
        end do
      end do
#endif
      do nn = 1 , nspgx
        nnb = iym1 - nn + 1
        do k = 1 , kz
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            tnb(nn,k,j) = tb0(nnb,k,j)
            tsb(nn,k,j) = tb0(nn,k,j)
            qnb(nn,k,j) = qb0(nnb,k,j)
            qsb(nn,k,j) = qb0(nn,k,j)
            tnbt(nn,k,j) = (tb1(nnb,k,j)-tb0(nnb,k,j))/dtbdys
            tsbt(nn,k,j) = (tb1(nn,k,j)-tb0(nn,k,j))/dtbdys
            qnbt(nn,k,j) = (qb1(nnb,k,j)-qb0(nnb,k,j))/dtbdys
            qsbt(nn,k,j) = (qb1(nn,k,j)-qb0(nn,k,j))/dtbdys
          end do
        end do
      end do
      print * , 'BCs are ready from ' , ndate0 , '  to ' , ndate1
      idatex = ndate0
      ndate0 = ndate1
#ifdef BAND
      do j = 1 , jx
#else
      do j = 1 , jxm1
#endif
        do i = 1 , iym1
          tdum(i,j) = ts1(i,j)
        end do
      end do
      do k = 1 , kz
        do j = 1 , jx
          do i = 1 , iy
            ub0(i,k,j) = ub1(i,k,j)
            vb0(i,k,j) = vb1(i,k,j)
            qb0(i,k,j) = qb1(i,k,j)
            tb0(i,k,j) = tb1(i,k,j)
          end do
        end do
      end do
      do j = 1 , jx
        do i = 1 , iy
          ps0(i,j) = ps1(i,j)
          ts0(i,j) = ts1(i,j)
        end do
      end do
      if ( ehso4 ) then
        do k = 1 , kz
          do j = 1 , jx
            do i = 1 , iy
              so0(i,k,j) = so1(i,k,j)
            end do
          end do
        end do
      end if
!bxqOCT2001_
 
      nnbase = nnnnnn
      call split_idate(mdate, nyear, nmonth, nday, nhour)
 
!-----------------------------------------------------------------------
      if ( ldatez.lt.ndate1 ) then
 
#ifdef BAND
        do j = 1 , jx
#else
        do j = 1 , jxm1
#endif
          do i = 1 , iym1
            if ( veg2d(i,j).le.0.00001 ) then
              if (idcsst == 1) then
                sts1%tg(i,j) = tdum(i,j) + dtskin(i,j)
                sts2%tg(i,j) = tdum(i,j) + dtskin(i,j)
              else
                sts1%tg(i,j) = tdum(i,j)
                sts2%tg(i,j) = tdum(i,j)
              end if
              if (iseaice == 1) then
                if ( tdum(i,j).le.271.38 ) then
                   sts1%tg(i,j) = 271.38
                   sts2%tg(i,j) = 271.38
                   tdum(i,j) = 271.38
                  do n = 1, nnsg
                    ocld2d(n,i,j) = 2.0
                    sice2d(n,i,j) = 1000.0
                    scv2d(n,i,j) = 0.0
                  end do
                else
                  sts1%tg(i,j) = tdum(i,j)
                  sts2%tg(i,j) = tdum(i,j)
                  do n = 1, nnsg
                    ocld2d(n,i,j) = 0.0
                    sice2d(n,i,j) = 0.0
                    scv2d(n,i,j) = 0.0
                  end do
                end if
              end if
            end if
          end do
        end do
      end if
#endif
      call time_end(subroutine_name,idindx)
      end subroutine bdyin
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine sets the boundary values of u and v according   c
!     to the boundary conditions specified.                           c
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
!
      subroutine bdyuv(ib,dtb)
!
#ifdef MPP1
#ifndef IBM
      use mpi
#else 
      include 'mpif.h'
#endif 
#endif
      implicit none
!
      real(8) :: dtb
      integer :: ib
      intent (in) dtb , ib
!
      integer :: i , j , k
#ifdef MPP1
      integer :: ierr
#endif
      character (len=50) :: subroutine_name='bdyuv'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
#ifdef MPP1
!----------------------------------------------------------------------
!-----compute the p* at dot points:
!
!=======================================================================
      call mpi_sendrecv(sps1%ps(1,jxp),iy,mpi_real8,ieast,1,      &
                      & sps1%ps(1,0),iy,mpi_real8,iwest,1,        &
                      & mpi_comm_world,mpi_status_ignore,ierr)
#endif
!=======================================================================
!
!-----interior points:
!
#ifdef MPP1
      do j = jbegin , jendx
        do i = 2 , iym1
          sps1%pdot(i,j) = 0.25*(sps1%ps(i,j)+sps1%ps(i-1,j)+     &
                                 sps1%ps(i,j-1)+sps1%ps(i-1,j-1))
        end do
      end do
#else
#ifdef BAND
      do j = 2 , jx
#else
      do j = 2 , jxm1
#endif
        do i = 2 , iym1
          sps1%pdot(i,j) = 0.25*(sps1%ps(i,j)+sps1%ps(i-1,j)+     &
                                 sps1%ps(i,j-1)+sps1%ps(i-1,j-1))
        end do
      end do
#ifdef BAND
      do i = 2 , iym1
        sps1%pdot(i,1) = 0.25*(sps1%ps(i,1)+sps1%ps(i-1,1)+   &
                               sps1%ps(i,jx)+sps1%ps(i-1,jx))
      enddo
#endif
#endif
!
!-----east and west boundaries:
!
#ifndef BAND
      do i = 2 , iym1
#ifdef MPP1
        if ( myid.eq.0 )  &
          sps1%pdot(i,1) = 0.5*(sps1%ps(i,1)+sps1%ps(i-1,1))
        if ( myid.eq.nproc-1 ) &
          sps1%pdot(i,jendl) = 0.5*(sps1%ps(i,jendx)+sps1%ps(i-1,jendx))
#else
        sps1%pdot(i,1) = 0.5*(sps1%ps(i,1)+sps1%ps(i-1,1))
        sps1%pdot(i,jx) = 0.5*(sps1%ps(i,jxm1)+sps1%ps(i-1,jxm1))
#endif
      end do
#endif
!
!-----north and south boundaries:
!
#ifdef MPP1
      do j = jbegin , jendx
        sps1%pdot(1,j) = 0.5*(sps1%ps(1,j)+sps1%ps(1,j-1))
        sps1%pdot(iy,j) = 0.5*(sps1%ps(iym1,j)+sps1%ps(iym1,j-1))
      end do
#else
#ifdef BAND
      do j = 2 , jx
        sps1%pdot(1,j) = 0.5*(sps1%ps(1,j)+sps1%ps(1,j-1))
        sps1%pdot(iy,j) = 0.5*(sps1%ps(iym1,j)+sps1%ps(iym1,j-1))
      end do
      sps1%pdot(1,1) = 0.5*(sps1%ps(1,1)+sps1%ps(1,jx))
      sps1%pdot(iy,1) = 0.5*(sps1%ps(iym1,1)+sps1%ps(iym1,jx))
#else
      do j = 2 , jxm1
        sps1%pdot(1,j) = 0.5*(sps1%ps(1,j)+sps1%ps(1,j-1))
        sps1%pdot(iy,j) = 0.5*(sps1%ps(iym1,j)+sps1%ps(iym1,j-1))
      end do
#endif
#endif
!
!-----corner points:
!
#ifndef BAND
#ifdef MPP1
      if ( myid.eq.0 ) then
        sps1%pdot(1,1) = sps1%ps(1,1)
        sps1%pdot(iy,1) = sps1%ps(iym1,1)
      end if
      if ( myid.eq.nproc-1 ) then
        sps1%pdot(1,jendl) = sps1%ps(1,jendx)
        sps1%pdot(iy,jendl) = sps1%ps(iym1,jendx)
      end if
#else
      sps1%pdot(1,1) = sps1%ps(1,1)
      sps1%pdot(iy,1) = sps1%ps(iym1,1)
      sps1%pdot(1,jx) = sps1%ps(1,jxm1)
      sps1%pdot(iy,jx) = sps1%ps(iym1,jxm1)
#endif
#endif
!=======================================================================
!
!-----interior silces:
!
      do k = 1 , kz
#ifndef BAND
!
!.....for j = 2 and j = jlx :
!
        do i = 2 , iym1
#ifdef MPP1
          if ( myid.eq.0 ) then
            uj2(i,k) = atm1%u(i,k,2)/sps1%pdot(i,2)
            vj2(i,k) = atm1%v(i,k,2)/sps1%pdot(i,2)
          end if
          if ( myid.eq.nproc-1 ) then
            ujlx(i,k) = atm1%u(i,k,jendx)/sps1%pdot(i,jendx)
            vjlx(i,k) = atm1%v(i,k,jendx)/sps1%pdot(i,jendx)
          end if
#else
          uj2(i,k) = atm1%u(i,k,2)/sps1%pdot(i,2)
          vj2(i,k) = atm1%v(i,k,2)/sps1%pdot(i,2)
          ujlx(i,k) = atm1%u(i,k,jxm1)/sps1%pdot(i,jxm1)
          vjlx(i,k) = atm1%v(i,k,jxm1)/sps1%pdot(i,jxm1)
#endif
        end do
#endif
!
!.....for i = 2 and i = iym1 :
!
#ifdef MPP1
#ifdef BAND
        do j = 1 , jendl
#else
        do j = jbegin , jendx
#endif
          ui2(k,j) = atm1%u(2,k,j)/sps1%pdot(2,j)
          vi2(k,j) = atm1%v(2,k,j)/sps1%pdot(2,j)
          uilx(k,j) = atm1%u(iym1,k,j)/sps1%pdot(iym1,j)
          vilx(k,j) = atm1%v(iym1,k,j)/sps1%pdot(iym1,j)
        end do
#else
#ifdef BAND
        do j = 1 , jx
#else
        do j = 2 , jxm1
#endif
          ui2(k,j) = atm1%u(2,k,j)/sps1%pdot(2,j)
          vi2(k,j) = atm1%v(2,k,j)/sps1%pdot(2,j)
          uilx(k,j) = atm1%u(iym1,k,j)/sps1%pdot(iym1,j)
          vilx(k,j) = atm1%v(iym1,k,j)/sps1%pdot(iym1,j)
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
        do k = 1 , kz
#ifndef BAND
!
!.....west (j = 1) and east (j = jx) boundaries:
!
          do i = 1 , iy
#ifdef MPP1
            if ( myid.eq.0 ) then
              uj1(i,k) = uwb(i,k,1)/sps1%pdot(i,1)
              vj1(i,k) = vwb(i,k,1)/sps1%pdot(i,1)
            end if
            if ( myid.eq.nproc-1 ) then
              ujl(i,k) = ueb(i,k,1)/sps1%pdot(i,jendl)
              vjl(i,k) = veb(i,k,1)/sps1%pdot(i,jendl)
            end if
#else
            uj1(i,k) = uwb(i,k,1)/sps1%pdot(i,1)
            vj1(i,k) = vwb(i,k,1)/sps1%pdot(i,1)
            ujl(i,k) = ueb(i,k,1)/sps1%pdot(i,jx)
            vjl(i,k) = veb(i,k,1)/sps1%pdot(i,jx)
#endif
          end do
#endif
!
!.....south (i = 1) and north (i = iy) boundaries:
!
#ifdef MPP1
          do j = 1 , jendl
            ui1(k,j) = usb(1,k,j)/sps1%pdot(1,j)
            vi1(k,j) = vsb(1,k,j)/sps1%pdot(1,j)
            uil(k,j) = unb(1,k,j)/sps1%pdot(iy,j)
            vil(k,j) = vnb(1,k,j)/sps1%pdot(iy,j)
          end do
#else
          do j = 1 , jx
            ui1(k,j) = usb(1,k,j)/sps1%pdot(1,j)
            vi1(k,j) = vsb(1,k,j)/sps1%pdot(1,j)
            uil(k,j) = unb(1,k,j)/sps1%pdot(iy,j)
            vil(k,j) = vnb(1,k,j)/sps1%pdot(iy,j)
          end do
#endif
        end do
      else
!
!-----time-dependent boundary conditions:
!
        do k = 1 , kz
#ifndef BAND
!
!.....west (j = 1) and east (j = jx) boundaries:
!
          do i = 1 , iy
#ifdef MPP1
            if ( myid.eq.0 ) then
              uj1(i,k) = (uwb(i,k,1)+dtb*uwbt(i,k,1))/sps1%pdot(i,1)
              vj1(i,k) = (vwb(i,k,1)+dtb*vwbt(i,k,1))/sps1%pdot(i,1)
            end if
            if ( myid.eq.nproc-1 ) then
              ujl(i,k) = (ueb(i,k,1)+dtb*uebt(i,k,1))/sps1%pdot(i,jendl)
              vjl(i,k) = (veb(i,k,1)+dtb*vebt(i,k,1))/sps1%pdot(i,jendl)
            end if
#else
            uj1(i,k) = (uwb(i,k,1)+dtb*uwbt(i,k,1))/sps1%pdot(i,1)
            vj1(i,k) = (vwb(i,k,1)+dtb*vwbt(i,k,1))/sps1%pdot(i,1)
            ujl(i,k) = (ueb(i,k,1)+dtb*uebt(i,k,1))/sps1%pdot(i,jx)
            vjl(i,k) = (veb(i,k,1)+dtb*vebt(i,k,1))/sps1%pdot(i,jx)
#endif
          end do
#endif
!
!.....south (i = 1) and north (i = iy) boundaries:
!
#ifdef MPP1
          do j = 1 , jendl
            ui1(k,j) = (usb(1,k,j)+dtb*usbt(1,k,j))/sps1%pdot(1,j)
            vi1(k,j) = (vsb(1,k,j)+dtb*vsbt(1,k,j))/sps1%pdot(1,j)
            uil(k,j) = (unb(1,k,j)+dtb*unbt(1,k,j))/sps1%pdot(iy,j)
            vil(k,j) = (vnb(1,k,j)+dtb*vnbt(1,k,j))/sps1%pdot(iy,j)
          end do
#else
          do j = 1 , jx
            ui1(k,j) = (usb(1,k,j)+dtb*usbt(1,k,j))/sps1%pdot(1,j)
            vi1(k,j) = (vsb(1,k,j)+dtb*vsbt(1,k,j))/sps1%pdot(1,j)
            uil(k,j) = (unb(1,k,j)+dtb*unbt(1,k,j))/sps1%pdot(iy,j)
            vil(k,j) = (vnb(1,k,j)+dtb*vnbt(1,k,j))/sps1%pdot(iy,j)
          end do
#endif
!
        end do
!
      end if
!
!-----fill up the interior silces:
!
#ifdef MPP1

#ifndef BAND
      do k = 1 , kz
        if ( myid.eq.0 ) then
          uj2(1,k) = ui1(k,2)
          uj2(iy,k) = uil(k,2)
          ui2(k,1) = uj1(2,k)
          uilx(k,1) = uj1(iym1,k)
          vj2(1,k) = vi1(k,2)
          vj2(iy,k) = vil(k,2)
          vi2(k,1) = vj1(2,k)
          vilx(k,1) = vj1(iym1,k)
        end if
        if ( myid.eq.nproc-1 ) then
          ujlx(1,k) = ui1(k,jendx)
          ujlx(iy,k) = uil(k,jendx)
          ui2(k,jendl) = ujl(2,k)
          uilx(k,jendl) = ujl(iym1,k)
          vjlx(1,k) = vi1(k,jendx)
          vjlx(iy,k) = vil(k,jendx)
          vi2(k,jendl) = vjl(2,k)
          vilx(k,jendl) = vjl(iym1,k)
        end if
      end do
#endif
#ifndef BAND
      if ( myid.ne.nproc-1 ) then
#endif
        do k = 1 , kz
          var1snd(k,1) = ui1(k,jxp)
          var1snd(k,2) = vi1(k,jxp)
          var1snd(k,3) = ui2(k,jxp)
          var1snd(k,4) = vi2(k,jxp)
          var1snd(k,5) = uilx(k,jxp)
          var1snd(k,6) = vilx(k,jxp)
          var1snd(k,7) = uil(k,jxp)
          var1snd(k,8) = vil(k,jxp)
        end do
#ifndef BAND
      end if
#endif
      call mpi_sendrecv(var1snd(1,1),kz*8,mpi_real8,ieast,1,            &
                      & var1rcv(1,1),kz*8,mpi_real8,iwest,1,            &
                      & mpi_comm_world,mpi_status_ignore,ierr)
#ifndef BAND
      if ( myid.ne.0 ) then
#endif
        do k = 1 , kz
          ui1(k,0) = var1rcv(k,1)
          vi1(k,0) = var1rcv(k,2)
          ui2(k,0) = var1rcv(k,3)
          vi2(k,0) = var1rcv(k,4)
          uilx(k,0) = var1rcv(k,5)
          vilx(k,0) = var1rcv(k,6)
          uil(k,0) = var1rcv(k,7)
          vil(k,0) = var1rcv(k,8)
        end do
#ifndef BAND
      end if
#endif
!
#ifndef BAND
      if ( myid.ne.0 ) then
#endif
        do k = 1 , kz
          var1snd(k,1) = ui1(k,1)
          var1snd(k,2) = vi1(k,1)
          var1snd(k,3) = ui2(k,1)
          var1snd(k,4) = vi2(k,1)
          var1snd(k,5) = uilx(k,1)
          var1snd(k,6) = vilx(k,1)
          var1snd(k,7) = uil(k,1)
          var1snd(k,8) = vil(k,1)
        end do
#ifndef BAND
      end if
#endif
      call mpi_sendrecv(var1snd(1,1),kz*8,mpi_real8,iwest,2,            &
                      & var1rcv(1,1),kz*8,mpi_real8,ieast,2,            &
                      & mpi_comm_world,mpi_status_ignore,ierr)
#ifndef BAND
      if ( myid.ne.nproc-1 ) then
#endif
        do k = 1 , kz
          ui1(k,jxp+1) = var1rcv(k,1)
          vi1(k,jxp+1) = var1rcv(k,2)
          ui2(k,jxp+1) = var1rcv(k,3)
          vi2(k,jxp+1) = var1rcv(k,4)
          uilx(k,jxp+1) = var1rcv(k,5)
          vilx(k,jxp+1) = var1rcv(k,6)
          uil(k,jxp+1) = var1rcv(k,7)
          vil(k,jxp+1) = var1rcv(k,8)
        end do
#ifndef BAND
      end if
#endif

#else

#ifndef BAND
      do k = 1 , kz
        uj2(1,k) = ui1(k,2)
        uj2(iy,k) = uil(k,2)
        ui2(k,1) = uj1(2,k)
        uilx(k,1) = uj1(iym1,k)
        vj2(1,k) = vi1(k,2)
        vj2(iy,k) = vil(k,2)
        vi2(k,1) = vj1(2,k)
        vilx(k,1) = vj1(iym1,k)
        ujlx(1,k) = ui1(k,jxm1)
        ujlx(iy,k) = uil(k,jxm1)
        ui2(k,jx) = ujl(2,k)
        uilx(k,jx) = ujl(iym1,k)
        vjlx(1,k) = vi1(k,jxm1)
        vjlx(iy,k) = vil(k,jxm1)
        vi2(k,jx) = vjl(2,k)
        vilx(k,jx) = vjl(iym1,k)
      end do
#endif

#endif
!
      call time_end(subroutine_name,idindx)
      end subroutine bdyuv
!
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
!
      subroutine bdyval(xt,iexec)
!
      implicit none
!
      integer :: iexec
      real(8) :: xt
      intent (in) iexec , xt
!
      real(8) :: chix , chix1 , chix2 , dtb , qcx , qcx2 , qvx , qvx1 , &
               & qvx2 , vavg
      integer :: itr , j , k
#ifndef BAND
      integer :: i
      real(8) :: uavg
#else
#ifndef MPP1
      integer :: jp1
#endif
#endif
      character (len=50) :: subroutine_name='bdyval'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
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
#ifndef BAND
        do i = 1 , iym1
          if ( myid.eq.0 ) sps2%ps(i,1) = sps1%ps(i,1)
          if ( myid.eq.nproc-1 ) sps2%ps(i,jendx) = sps1%ps(i,jendx)
        end do
#endif
        do j = jbegin , jendm
          sps2%ps(1,j) = sps1%ps(1,j)
          sps2%ps(iym1,j) = sps1%ps(iym1,j)
        end do
!
!-----for p*u and p*v:
!
        do k = 1 , kz
#ifndef BAND
          do i = 1 , iy
            if ( myid.eq.0 ) then
              atm2%u(i,k,1) = atm1%u(i,k,1)/mddom%msfd(i,1)
              atm2%v(i,k,1) = atm1%v(i,k,1)/mddom%msfd(i,1)
            end if
            if ( myid.eq.nproc-1 ) then
              atm2%u(i,k,jendl) = atm1%u(i,k,jendl)/mddom%msfd(i,jendl)
              atm2%v(i,k,jendl) = atm1%v(i,k,jendl)/mddom%msfd(i,jendl)
            end if
          end do
#endif
          do j = jbegin , jendx
            atm2%u(1,k,j) = atm1%u(1,k,j)/mddom%msfd(1,j)
            atm2%u(iy,k,j) = atm1%u(iy,k,j)/mddom%msfd(iy,j)
            atm2%v(1,k,j) = atm1%v(1,k,j)/mddom%msfd(1,j)
            atm2%v(iy,k,j) = atm1%v(iy,k,j)/mddom%msfd(iy,j)
          end do
        end do
!
!-----for p*t:
!
        do k = 1 , kz
#ifndef BAND
          do i = 1 , iym1
            if ( myid.eq.0 ) atm2%t(i,k,1) = atm1%t(i,k,1)
            if ( myid.eq.nproc-1 ) atm2%t(i,k,jendx) = atm1%t(i,k,jendx)
          end do
#endif
          do j = jbegin , jendm
            atm2%t(1,k,j) = atm1%t(1,k,j)
            atm2%t(iym1,k,j) = atm1%t(iym1,k,j)
          end do
        end do
!
!-----for p*qv:
!
        do k = 1 , kz
#ifndef BAND
          do i = 1 , iym1
            if ( myid.eq.0 ) atm2%qv(i,k,1) = atm1%qv(i,k,1)
            if ( myid.eq.nproc-1 )  &
              atm2%qv(i,k,jendx) = atm1%qv(i,k,jendx)
          end do
#endif
          do j = jbegin , jendm
            atm2%qv(1,k,j) = atm1%qv(1,k,j)
            atm2%qv(iym1,k,j) = atm1%qv(iym1,k,j)
          end do
        end do
!
!chem2
!
        if ( ichem.eq.1 ) then
!-----for p*chi (tracers)
          do itr = 1 , ntr
            do k = 1 , kz
#ifndef BAND
              do i = 1 , iym1
                if ( myid.eq.0 ) chib(i,k,1,itr) = chia(i,k,1,itr)
                if ( myid.eq.nproc-1 ) chib(i,k,jendx,itr)              &
                   & = chia(i,k,jendx,itr)
              end do
#endif
              do j = jbegin , jendm
                chib(1,k,j,itr) = chia(1,k,j,itr)
                chib(iym1,k,j,itr) = chia(iym1,k,j,itr)
              end do
            end do
          end do
        end if
!chem2_
!
!-----for p*qc:
!
        do k = 1 , kz
#ifndef BAND
          do i = 1 , iym1
            if ( myid.eq.0 ) atm2%qc(i,k,1) = atm1%qc(i,k,1)
            if ( myid.eq.nproc-1 ) &
              atm2%qc(i,k,jendx) = atm1%qc(i,k,jendx)
          end do
#endif
          do j = jbegin , jendm
            atm2%qc(1,k,j) = atm1%qc(1,k,j)
            atm2%qc(iym1,k,j) = atm1%qc(iym1,k,j)
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
#ifndef BAND
          do i = 1 , iym1
            if ( myid.eq.0 ) sps1%ps(i,1) = pwb(i,1)
            if ( myid.eq.nproc-1 ) sps1%ps(i,jendx) = peb(i,1)
          end do
#endif
          do j = jbegin , jendm
            sps1%ps(1,j) = pss(1,j)
            sps1%ps(iym1,j) = pnb(1,j)
          end do
!
          do k = 1 , kz
#ifndef BAND
            do i = 1 , iy
              if ( myid.eq.0 ) then
                atm1%u(i,k,1) = uwb(i,k,1)
                atm1%v(i,k,1) = vwb(i,k,1)
              end if
              if ( myid.eq.nproc-1 ) then
                atm1%u(i,k,jendl) = ueb(i,k,1)
                atm1%v(i,k,jendl) = veb(i,k,1)
              end if
            end do
#endif
            do j = jbegin , jendx
              atm1%u(1,k,j) = usb(1,k,j)
              atm1%u(iy,k,j) = unb(1,k,j)
              atm1%v(1,k,j) = vsb(1,k,j)
              atm1%v(iy,k,j) = vnb(1,k,j)
            end do
          end do
        end if
!
!.....time-dependent boundary conditions:
!
#ifndef BAND
        do i = 1 , iym1
          if ( myid.eq.0 ) sps1%ps(i,1) = pwb(i,1) + dtb*pwbt(i,1)
          if ( myid.eq.nproc-1 )  &
            sps1%ps(i,jendx) = peb(i,1) + dtb*pebt(i,1)
        end do
#endif
        do j = jbegin , jendm
          sps1%ps(1,j) = pss(1,j) + dtb*psbt(1,j)
          sps1%ps(iym1,j) = pnb(1,j) + dtb*pnbt(1,j)
        end do
!
        do k = 1 , kz
#ifndef BAND
          do i = 1 , iy
            if ( myid.eq.0 ) then
              atm1%u(i,k,1) = uwb(i,k,1) + dtb*uwbt(i,k,1)
              atm1%v(i,k,1) = vwb(i,k,1) + dtb*vwbt(i,k,1)
            end if
            if ( myid.eq.nproc-1 ) then
              atm1%u(i,k,jendl) = ueb(i,k,1) + dtb*uebt(i,k,1)
              atm1%v(i,k,jendl) = veb(i,k,1) + dtb*vebt(i,k,1)
            end if
          end do
#endif
          do j = jbegin , jendx
            atm1%u(1,k,j) = usb(1,k,j) + dtb*usbt(1,k,j)
            atm1%u(iy,k,j) = unb(1,k,j) + dtb*unbt(1,k,j)
            atm1%v(1,k,j) = vsb(1,k,j) + dtb*vsbt(1,k,j)
            atm1%v(iy,k,j) = vnb(1,k,j) + dtb*vnbt(1,k,j)
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
        do k = 1 , kz
#ifndef BAND
          do i = 1 , iym1
            if ( myid.eq.0 ) atm1%t(i,k,1) = twb(i,k,1)
            if ( myid.eq.nproc-1 ) atm1%t(i,k,jendx) = teb(i,k,1)
          end do
#endif
          do j = jbegin , jendm
            atm1%t(1,k,j) = tsb(1,k,j)
            atm1%t(iym1,k,j) = tnb(1,k,j)
          end do
        end do
        do k = 1 , kz
#ifndef BAND
          do i = 1 , iym1
            if ( myid.eq.0 ) atm1%qv(i,k,1) = qwb(i,k,1)
            if ( myid.eq.nproc-1 ) atm1%qv(i,k,jendx) = qeb(i,k,1)
          end do
#endif
          do j = jbegin , jendm
            atm1%qv(1,k,j) = qsb(1,k,j)
            atm1%qv(iym1,k,j) = qnb(1,k,j)
          end do
        end do
      end if
!
!.....time-dependent boundary conditions:
!
      do k = 1 , kz
#ifndef BAND
        do i = 1 , iym1
          if ( myid.eq.0 ) atm1%t(i,k,1) = twb(i,k,1) + dtb*twbt(i,k,1)
          if ( myid.eq.nproc-1 ) atm1%t(i,k,jendx) = teb(i,k,1)   &
             & + dtb*tebt(i,k,1)
        end do
#endif
        do j = jbegin , jendm
          atm1%t(1,k,j) = tsb(1,k,j) + dtb*tsbt(1,k,j)
          atm1%t(iym1,k,j) = tnb(1,k,j) + dtb*tnbt(1,k,j)
        end do
      end do
      do k = 1 , kz
#ifndef BAND
        do i = 1 , iym1
          if ( myid.eq.0 ) atm1%qv(i,k,1) = qwb(i,k,1) + dtb*qwbt(i,k,1)
          if ( myid.eq.nproc-1 ) atm1%qv(i,k,jendx) = qeb(i,k,1)        &
             & + dtb*qebt(i,k,1)
        end do
#endif
        do j = jbegin , jendm
          atm1%qv(1,k,j) = qsb(1,k,j) + dtb*qsbt(1,k,j)
          atm1%qv(iym1,k,j) = qnb(1,k,j) + dtb*qnbt(1,k,j)
        end do
      end do
!
      if ( iboudy.eq.3 .or. iboudy.eq.4 ) then
!
!-----determine boundary values depends on inflow/outflow:
!
        do k = 1 , kz
#ifndef BAND
!
!.....west boundary:
!
          if ( myid.eq.0 ) then
            do i = 1 , iym1
              qvx1 = atm1%qv(i,k,1)/sps1%ps(i,1)
              qvx2 = atm1%qv(i,k,2)/sps1%ps(i,2)
              uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
              if ( uavg.ge.0. ) then
                qvx = qvx1
              else
                qvx = qvx2
              end if
              atm1%qv(i,k,1) = qvx*sps1%ps(i,1)
            end do
          end if
!
!.....east boundary:
!
          if ( myid.eq.nproc-1 ) then
            do i = 1 , iym1
              qvx1 = atm1%qv(i,k,jendx)/sps1%ps(i,jendx)
              qvx2 = atm1%qv(i,k,jendm)/sps1%ps(i,jendm)
              uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
              if ( uavg.lt.0. ) then
                qvx = qvx1
              else
                qvx = qvx2
              end if
              atm1%qv(i,k,jendx) = qvx*sps1%ps(i,jendx)
            end do
          end if
#endif
!
!.....south boundary:
!
          do j = jbegin , jendm
            qvx1 = atm1%qv(1,k,j)/sps1%ps(1,j)
            qvx2 = atm1%qv(2,k,j)/sps1%ps(2,j)
            vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
            if ( vavg.ge.0. ) then
              qvx = qvx1
            else
              qvx = qvx2
            end if
            atm1%qv(1,k,j) = qvx*sps1%ps(1,j)
          end do
!
!.....north boundary:
!
          do j = jbegin , jendm
            qvx1 = atm1%qv(iym1,k,j)/sps1%ps(iym1,j)
            qvx2 = atm1%qv(iym2,k,j)/sps1%ps(iym2,j)
            vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
            if ( vavg.lt.0. ) then
              qvx = qvx1
            else
              qvx = qvx2
            end if
            atm1%qv(iym1,k,j) = qvx*sps1%ps(iym1,j)
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
      do k = 1 , kz
#ifndef BAND
!
!.....west boundary:
!
        if ( myid.eq.0 ) then
          do i = 1 , iym1
            qcx2 = atm1%qc(i,k,2)/sps1%ps(i,2)
            uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
            if ( uavg.ge.0. ) then
              qcx = 0.
            else
              qcx = qcx2
            end if
            atm1%qc(i,k,1) = qcx*sps1%ps(i,1)
          end do
        end if
!
!.....east boundary:
!
        if ( myid.eq.nproc-1 ) then
          do i = 1 , iym1
            qcx2 = atm1%qc(i,k,jendm)/sps1%ps(i,jendm)
            uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
            if ( uavg.lt.0. ) then
              qcx = 0.
            else
              qcx = qcx2
            end if
            atm1%qc(i,k,jendx) = qcx*sps1%ps(i,jendx)
          end do
        end if
#endif
!
!.....south boundary:
!
        do j = jbegin , jendm
          qcx2 = atm1%qc(2,k,j)/sps1%ps(2,j)
          vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
          if ( vavg.ge.0. ) then
            qcx = 0.
          else
            qcx = qcx2
          end if
          atm1%qc(1,k,j) = qcx*sps1%ps(1,j)
        end do
!
!.....north boundary:
!
        do j = jbegin , jendm
          qcx2 = atm1%qc(iym2,k,j)/sps1%ps(iym2,j)
          vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
          if ( vavg.lt.0. ) then
            qcx = 0.
          else
            qcx = qcx2
          end if
          atm1%qc(iym1,k,j) = qcx*sps1%ps(iym1,j)
        end do
!
      end do
 
      if ( ichem.eq.1 ) then
!chem2
 
!----add tracer bc's
!
        do itr = 1 , ntr
          do k = 1 , kz
#ifndef BAND
!
!.....west  boundary:
!
            if ( myid.eq.0 ) then
              do i = 1 , iym1
                chix1 = chia(i,k,1,itr)/sps1%ps(i,1)
                chix2 = chia(i,k,2,itr)/sps1%ps(i,2)
                uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
                if ( uavg.ge.0. ) then
                  chix = chix1
                else
                  chix = chix2
                end if
                chia(i,k,1,itr) = chix*sps1%ps(i,1)
              end do
            end if
!
!.....east  boundary:
!
            if ( myid.eq.nproc-1 ) then
              do i = 1 , iym1
                chix1 = chia(i,k,jendx,itr)/sps1%ps(i,jendx)
                chix2 = chia(i,k,jendm,itr)/sps1%ps(i,jendm)
                uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
                if ( uavg.lt.0. ) then
                  chix = chix1
                else
                  chix = chix2
                end if
                chia(i,k,jendx,itr) = chix*sps1%ps(i,jendx)
              end do
            end if
#endif
!
!.....south boundary:
!
            do j = jbegin , jendm
              chix1 = chia(1,k,j,itr)/sps1%ps(1,j)
              chix2 = chia(2,k,j,itr)/sps1%ps(2,j)
              vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
              if ( vavg.ge.0. ) then
                chix = chix1
              else
                chix = chix2
              end if
              chia(1,k,j,itr) = chix*sps1%ps(1,j)
            end do
!
!.....north boundary:
!
            do j = jbegin , jendm
              chix1 = chia(iym1,k,j,itr)/sps1%ps(iym1,j)
              chix2 = chia(iym2,k,j,itr)/sps1%ps(iym2,j)
              vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
              if ( vavg.lt.0. ) then
                chix = chix1
              else
                chix = chix2
              end if
              chia(iym1,k,j,itr) = chix*sps1%ps(iym1,j)
            end do
          end do
        end do
!chem2_
      end if
#else
!
!-----for p*:
!
#ifndef BAND
        do i = 1 , iym1
          sps2%ps(i,1) = sps1%ps(i,1)
          sps2%ps(i,jxm1) = sps1%ps(i,jxm1)
        end do
#endif
#ifdef BAND
        do j = 1 , jx
#else
        do j = 2 , jxm2
#endif
          sps2%ps(1,j) = sps1%ps(1,j)
          sps2%ps(iym1,j) = sps1%ps(iym1,j)
        end do
!
!-----for p*u and p*v:
!
        do k = 1 , kz
#ifndef BAND
          do i = 1 , iy
            atm2%u(i,k,1) = atm1%u(i,k,1)/mddom%msfd(i,1)
            atm2%v(i,k,1) = atm1%v(i,k,1)/mddom%msfd(i,1)
            atm2%u(i,k,jx) = atm1%u(i,k,jx)/mddom%msfd(i,jx)
            atm2%v(i,k,jx) = atm1%v(i,k,jx)/mddom%msfd(i,jx)
          end do
#endif
#ifdef BAND
          do j = 1 , jx
#else
          do j = 2 , jxm1
#endif
            atm2%u(1,k,j) = atm1%u(1,k,j)/mddom%msfd(1,j)
            atm2%u(iy,k,j) = atm1%u(iy,k,j)/mddom%msfd(iy,j)
            atm2%v(1,k,j) = atm1%v(1,k,j)/mddom%msfd(1,j)
            atm2%v(iy,k,j) = atm1%v(iy,k,j)/mddom%msfd(iy,j)
          end do
        end do
!
!-----for p*t:
!
        do k = 1 , kz
#ifndef BAND
          do i = 1 , iym1
            atm2%t(i,k,1) = atm1%t(i,k,1)
            atm2%t(i,k,jxm1) = atm1%t(i,k,jxm1)
          end do
#endif
#ifdef BAND
          do j = 1 , jx
#else
          do j = 2 , jxm2
#endif
            atm2%t(1,k,j) = atm1%t(1,k,j)
            atm2%t(iym1,k,j) = atm1%t(iym1,k,j)
          end do
        end do
!
!-----for p*qv:
!
        do k = 1 , kz
#ifndef BAND
          do i = 1 , iym1
            atm2%qv(i,k,1) = atm1%qv(i,k,1)
            atm2%qv(i,k,jxm1) = atm1%qv(i,k,jxm1)
          end do
#endif
#ifdef BAND
          do j = 1 , jx
#else
          do j = 2 , jxm2
#endif
            atm2%qv(1,k,j) = atm1%qv(1,k,j)
            atm2%qv(iym1,k,j) = atm1%qv(iym1,k,j)
          end do
        end do
!
!chem2
!
        if ( ichem.eq.1 ) then
!-----for p*chi (tracers)
          do itr = 1 , ntr
            do k = 1 , kz
#ifndef BAND
              do i = 1 , iym1
                chib(i,k,1,itr) = chia(i,k,1,itr)
                chib(i,k,jxm1,itr) = chia(i,k,jxm1,itr)
              end do
#endif
#ifdef BAND
              do j = 1 , jx
#else
              do j = 2 , jxm2
#endif
                chib(1,k,j,itr) = chia(1,k,j,itr)
                chib(iym1,k,j,itr) = chia(iym1,k,j,itr)
              end do
            end do
          end do
        end if
!chem2_
!
!-----for p*qc:
!
        do k = 1 , kz
#ifndef BAND
          do i = 1 , iym1
            atm2%qc(i,k,1) = atm1%qc(i,k,1)
            atm2%qc(i,k,jxm1) = atm1%qc(i,k,jxm1)
          end do
#endif
#ifdef BAND
          do j = 1 , jx
#else
          do j = 2 , jxm2
#endif
            atm2%qc(1,k,j) = atm1%qc(1,k,j)
            atm2%qc(iym1,k,j) = atm1%qc(iym1,k,j)
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
#ifndef BAND
          do i = 1 , iym1
            sps1%ps(i,1) = pwb(i,1)
            sps1%ps(i,jxm1) = peb(i,1)
          end do
#endif
#ifdef BAND
          do j = 1 , jx
#else
          do j = 2 , jxm2
#endif
            sps1%ps(1,j) = pss(1,j)
            sps1%ps(iym1,j) = pnb(1,j)
          end do
!
          do k = 1 , kz
#ifndef BAND
            do i = 1 , iy
              atm1%u(i,k,1) = uwb(i,k,1)
              atm1%v(i,k,1) = vwb(i,k,1)
              atm1%u(i,k,jx) = ueb(i,k,1)
              atm1%v(i,k,jx) = veb(i,k,1)
            end do
#endif
#ifdef BAND
            do j = 1 , jx
#else
            do j = 2 , jxm1
#endif
              atm1%u(1,k,j) = usb(1,k,j)
              atm1%u(iy,k,j) = unb(1,k,j)
              atm1%v(1,k,j) = vsb(1,k,j)
              atm1%v(iy,k,j) = vnb(1,k,j)
            end do
          end do
        end if
!
!.....time-dependent boundary conditions:
!
#ifndef BAND
        do i = 1 , iym1
          sps1%ps(i,1) = pwb(i,1) + dtb*pwbt(i,1)
          sps1%ps(i,jxm1) = peb(i,1) + dtb*pebt(i,1)
        end do
#endif
#ifdef BAND
        do j = 1 , jx
#else
        do j = 2 , jxm2
#endif
          sps1%ps(1,j) = pss(1,j) + dtb*psbt(1,j)
          sps1%ps(iym1,j) = pnb(1,j) + dtb*pnbt(1,j)
        end do
!
        do k = 1 , kz
#ifndef BAND
          do i = 1 , iy
            atm1%u(i,k,1) = uwb(i,k,1) + dtb*uwbt(i,k,1)
            atm1%v(i,k,1) = vwb(i,k,1) + dtb*vwbt(i,k,1)
            atm1%u(i,k,jx) = ueb(i,k,1) + dtb*uebt(i,k,1)
            atm1%v(i,k,jx) = veb(i,k,1) + dtb*vebt(i,k,1)
          end do
#endif
#ifdef BAND
          do j = 1 , jx
#else
          do j = 2 , jxm1
#endif
            atm1%u(1,k,j) = usb(1,k,j) + dtb*usbt(1,k,j)
            atm1%u(iy,k,j) = unb(1,k,j) + dtb*unbt(1,k,j)
            atm1%v(1,k,j) = vsb(1,k,j) + dtb*vsbt(1,k,j)
            atm1%v(iy,k,j) = vnb(1,k,j) + dtb*vnbt(1,k,j)
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
        do k = 1 , kz
#ifndef BAND
          do i = 1 , iym1
            atm1%t(i,k,1) = twb(i,k,1)
            atm1%t(i,k,jxm1) = teb(i,k,1)
          end do
#endif
#ifdef BAND
          do j = 1 , jx
#else
          do j = 2 , jxm2
#endif
            atm1%t(1,k,j) = tsb(1,k,j)
            atm1%t(iym1,k,j) = tnb(1,k,j)
          end do
        end do
        do k = 1 , kz
#ifndef BAND
          do i = 1 , iym1
            atm1%qv(i,k,1) = qwb(i,k,1)
            atm1%qv(i,k,jxm1) = qeb(i,k,1)
          end do
#endif
#ifdef BAND
          do j = 1 , jx
#else
          do j = 2 , jxm2
#endif
            atm1%qv(1,k,j) = qsb(1,k,j)
            atm1%qv(iym1,k,j) = qnb(1,k,j)
          end do
        end do
      end if
!
!.....time-dependent boundary conditions:
!
      do k = 1 , kz
#ifndef BAND
        do i = 1 , iym1
          atm1%t(i,k,1) = twb(i,k,1) + dtb*twbt(i,k,1)
          atm1%t(i,k,jxm1) = teb(i,k,1) + dtb*tebt(i,k,1)
        end do
#endif
#ifdef BAND
        do j = 1 , jx
#else
        do j = 2 , jxm2
#endif
          atm1%t(1,k,j) = tsb(1,k,j) + dtb*tsbt(1,k,j)
          atm1%t(iym1,k,j) = tnb(1,k,j) + dtb*tnbt(1,k,j)
        end do
      end do
      do k = 1 , kz
#ifndef BAND
        do i = 1 , iym1
          atm1%qv(i,k,1) = qwb(i,k,1) + dtb*qwbt(i,k,1)
          atm1%qv(i,k,jxm1) = qeb(i,k,1) + dtb*qebt(i,k,1)
        end do
#endif
#ifdef BAND
        do j = 1 , jx
#else
        do j = 2 , jxm2
#endif
          atm1%qv(1,k,j) = qsb(1,k,j) + dtb*qsbt(1,k,j)
          atm1%qv(iym1,k,j) = qnb(1,k,j) + dtb*qnbt(1,k,j)
        end do
      end do
!
      if ( iboudy.eq.3 .or. iboudy.eq.4 ) then
!
!-----determine boundary values depends on inflow/outflow:
!
        do k = 1 , kz
#ifndef BAND
!
!.....west boundary:
!
          do i = 1 , iym1
            qvx1 = atm1%qv(i,k,1)/sps1%ps(i,1)
            qvx2 = atm1%qv(i,k,2)/sps1%ps(i,2)
            uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
            if ( uavg.ge.0. ) then
              qvx = qvx1
            else
              qvx = qvx2
            end if
            atm1%qv(i,k,1) = qvx*sps1%ps(i,1)
          end do
!
!.....east boundary:
!
          do i = 1 , iym1
            qvx1 = atm1%qv(i,k,jxm1)/sps1%ps(i,jxm1)
            qvx2 = atm1%qv(i,k,jxm2)/sps1%ps(i,jxm2)
            uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
            if ( uavg.lt.0. ) then
              qvx = qvx1
            else
              qvx = qvx2
            end if
            atm1%qv(i,k,jxm1) = qvx*sps1%ps(i,jxm1)
          end do
#endif
!
!.....south boundary:
!
#ifdef BAND
          do j = 1 , jx
            jp1 = j+1
            if(jp1.eq.jx+1) jp1=1
            qvx1 = atm1%qv(1,k,j)/sps1%ps(1,j)
            qvx2 = atm1%qv(2,k,j)/sps1%ps(2,j)
            vavg = vi1(k,j) + vi1(k,jp1) + vi2(k,j) + vi2(k,jp1)
#else
          do j = 2 , jxm2
            qvx1 = atm1%qv(1,k,j)/sps1%ps(1,j)
            qvx2 = atm1%qv(2,k,j)/sps1%ps(2,j)
            vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
#endif
            if ( vavg.ge.0. ) then
              qvx = qvx1
            else
              qvx = qvx2
            end if
            atm1%qv(1,k,j) = qvx*sps1%ps(1,j)
          end do
!
!.....north boundary:
!
#ifdef BAND
          do j = 1 , jx
            jp1 = j+1
            if(jp1.eq.jx+1) jp1=1
            qvx1 = atm1%qv(iym1,k,j)/sps1%ps(iym1,j)
            qvx2 = atm1%qv(iym2,k,j)/sps1%ps(iym2,j)
            vavg = vilx(k,j) + vilx(k,jp1) + vil(k,j) + vil(k,jp1)
#else
          do j = 2 , jxm2
            qvx1 = atm1%qv(iym1,k,j)/sps1%ps(iym1,j)
            qvx2 = atm1%qv(iym2,k,j)/sps1%ps(iym2,j)
            vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
#endif
            if ( vavg.lt.0. ) then
              qvx = qvx1
            else
              qvx = qvx2
            end if
            atm1%qv(iym1,k,j) = qvx*sps1%ps(iym1,j)
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
      do k = 1 , kz
#ifndef BAND
!
!.....west boundary:
!
        do i = 1 , iym1
          qcx2 = atm1%qc(i,k,2)/sps1%ps(i,2)
          uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
          if ( uavg.ge.0. ) then
            qcx = 0.
          else
            qcx = qcx2
          end if
          atm1%qc(i,k,1) = qcx*sps1%ps(i,1)
        end do
!
!.....east boundary:
!
        do i = 1 , iym1
          qcx2 = atm1%qc(i,k,jxm2)/sps1%ps(i,jxm2)
          uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
          if ( uavg.lt.0. ) then
            qcx = 0.
          else
            qcx = qcx2
          end if
          atm1%qc(i,k,jxm1) = qcx*sps1%ps(i,jxm1)
        end do
#endif
!
!.....south boundary:
!
#ifdef BAND
        do j = 1 , jx
          jp1 = j+1
          if(jp1.eq.jx+1) jp1=1
          qcx2 = atm1%qc(2,k,j)/sps1%ps(2,j)
          vavg = vi1(k,j) + vi1(k,jp1) + vi2(k,j) + vi2(k,jp1)
#else
        do j = 2 , jxm2
          qcx2 = atm1%qc(2,k,j)/sps1%ps(2,j)
          vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
#endif
          if ( vavg.ge.0. ) then
            qcx = 0.
          else
            qcx = qcx2
          end if
          atm1%qc(1,k,j) = qcx*sps1%ps(1,j)
        end do
!
!.....north boundary:
!
#ifdef BAND
        do j = 1 , jx
          jp1 = j+1
          if(jp1.eq.jx+1) jp1=1
          qcx2 = atm1%qc(iym2,k,j)/sps1%ps(iym2,j)
          vavg = vilx(k,j) + vilx(k,jp1) + vil(k,j) + vil(k,jp1)
#else
        do j = 2 , jxm2
          qcx2 = atm1%qc(iym2,k,j)/sps1%ps(iym2,j)
          vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
#endif
          if ( vavg.lt.0. ) then
            qcx = 0.
          else
            qcx = qcx2
          end if
          atm1%qc(iym1,k,j) = qcx*sps1%ps(iym1,j)
        end do
!
      end do
 
      if ( ichem.eq.1 ) then
!chem2
 
!----add tracer bc's
!
        do itr = 1 , ntr
          do k = 1 , kz
#ifndef BAND
!
!.....west  boundary:
!
 
            do i = 1 , iym1
!             FAB force to zero inflow conditions
!             chix1 = chia(i,k,1,itr)/sps1%ps(i,1)
              chix1 = 0. 
              chix2 = chia(i,k,2,itr)/sps1%ps(i,2)
              uavg = uj1(i,k) + uj1(i+1,k) + uj2(i,k) + uj2(i+1,k)
              if ( uavg.ge.0. ) then
                chix = chix1
              else
                chix = chix2
              end if
              chia(i,k,1,itr) = chix*sps1%ps(i,1)
            end do
!
!.....east  boundary:
!
            do i = 1 , iym1
!             chix1 = chia(i,k,jxm1,itr)/sps1%ps(i,jxm1)
              chix1 = 0.
              chix2 = chia(i,k,jxm2,itr)/sps1%ps(i,jxm2)
              uavg = ujlx(i,k) + ujlx(i+1,k) + ujl(i,k) + ujl(i+1,k)
              if ( uavg.lt.0. ) then
                chix = chix1
              else
                chix = chix2
              end if
              chia(i,k,jxm1,itr) = chix*sps1%ps(i,jxm1)
            end do
#endif
!
!.....south boundary:
!
#ifdef BAND
            do j = 1 , jx
              jp1 = j+1
              if(jp1.eq.jx+1) jp1=1
              chix1 = 0.
              chix2 = chia(2,k,j,itr)/sps1%ps(2,j)
              vavg = vi1(k,j) + vi1(k,jp1) + vi2(k,j) + vi2(k,jp1)
#else
            do j = 2 , jxm2
!             chix1 = chia(1,k,j,itr)/sps1%ps(1,j)
              chix1 = 0.
              chix2 = chia(2,k,j,itr)/sps1%ps(2,j)
              vavg = vi1(k,j) + vi1(k,j+1) + vi2(k,j) + vi2(k,j+1)
#endif
              if ( vavg.ge.0. ) then
                chix = chix1
              else
                chix = chix2
              end if
              chia(1,k,j,itr) = chix*sps1%ps(1,j)
            end do
!
!.....north boundary:
!
#ifdef BAND
            do j = 1 , jx
              jp1 = j+1
              if(jp1.eq.jx+1) jp1=1
!             chix1 = chia(iym1,k,j,itr)/sps1%ps(iym1,j)
              chix1 = 0.
              chix2 = chia(iym2,k,j,itr)/sps1%ps(iym2,j)
              vavg = vilx(k,j) + vilx(k,jp1) + vil(k,j) + vil(k,jp1)
#else
            do j = 2 , jxm2
!             chix1 = chia(iym1,k,j,itr)/sps1%ps(iym1,j)
              chix1 = 0.
              chix2 = chia(iym2,k,j,itr)/sps1%ps(iym2,j)
              vavg = vilx(k,j) + vilx(k,j+1) + vil(k,j) + vil(k,j+1)
#endif
              if ( vavg.lt.0. ) then
                chix = chix1
              else
                chix = chix2
              end if
              chia(iym1,k,j,itr) = chix*sps1%ps(iym1,j)
            end do
          end do
        end do
!chem2_
      end if
#endif
!
      call time_end(subroutine_name,idindx)
      end subroutine bdyval
!
      end module mod_bdycod
