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
 
      module mod_tendency

      use mod_runparams
      use mod_main
      use mod_mainchem
      use mod_bdycod
      use mod_cvaria
      use mod_pmoist
      use mod_precip
      use mod_rad
      use mod_bats
      use mod_lake
      use mod_vecbats
      use mod_holtbl
      use mod_trachem
      use mod_colmod3
      use mod_cu_grell
      use mod_cu_kuo
      use mod_cu_bm
      use mod_cu_em
      use mod_date
      use mod_message
      use mod_aerosol
      use mod_zengocn
      use mod_sun
      use mod_slice
      use mod_cldfrac
      use mod_cumtran
      use mod_condtq
      use mod_diffusion
      use mod_advection
      use mod_nudge
      use mod_che_tend
      use mod_diagnosis
      use mod_service
#ifdef MPP1
#ifdef CLM
      use mod_clm
      use mod_mtrxclm
      use clm_varsur
#endif
#endif

      private

      public :: allocate_mod_tend , tend

      real(8) , allocatable , dimension(:,:) :: divl
      real(8) , allocatable , dimension(:,:,:) :: ttld , xkc , td
      real(8) , allocatable , dimension(:,:) :: bdyewrcv , bdyewsnd
      real(8) , allocatable , dimension(:,:) :: bdynsrcv , bdynssnd
      real(8) , allocatable , dimension(:,:,:) :: ps4
      real(8) , allocatable , dimension(:,:,:) :: ps_4 
      real(8) , allocatable , dimension(:,:) :: var2rcv , var2snd
      real(8) , allocatable , dimension(:,:) :: tvar1rcv , tvar1snd

      contains

      subroutine allocate_mod_tend(lmpi,lband)
        implicit none
        logical , intent(in) :: lmpi , lband

        allocate(divl(iy,kz))
        if (lmpi) then
          allocate(ttld(iy,kz,jxp))
          allocate(xkc(iy,kz,jxp))
          allocate(td(iy,kz,jxp))
          if ( .not. lband ) then
            allocate(bdyewrcv(iy,kz*16+4))
            allocate(bdyewsnd(iy,kz*16+4))
          end if
          allocate(bdynsrcv(nspgx,kz*16+4))
          allocate(bdynssnd(nspgx,kz*16+4))
          allocate(ps4(iy,4,jxp))
          allocate(ps_4(iy,4,jx))
          allocate(var2rcv(iy,kz*(ntr+5)*2))
          allocate(var2snd(iy,kz*(ntr+5)*2))
          allocate(tvar1rcv(iy,kz*11+1+ntr*kz*2))
          allocate(tvar1snd(iy,kz*11+1+ntr*kz*2))
        else
          allocate(ttld(iy,kz,jx))
          allocate(xkc(iy,kz,jx))
          allocate(td(iy,kz,jx))
        end if
      end subroutine allocate_mod_tend

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
#ifdef MPP1
#ifndef IBM
      use mpi
#else
      include 'mpif.h'
#endif
#endif
      implicit none
!
      integer :: iexec
      intent (inout) iexec
!
      real(8) :: cell , chias , chibs , dudx , dudy ,               &
               & dvdx , dvdy , psabar , psasum ,                    &
               & pt2bar , pt2tot , ptnbar , ptntot , qcas , qcbs ,  &
               & qvas , qvbs , rovcpm , rtbar , sigpsa , tv ,       &
               & tv1 , tv2 , tv3 , tv4 , tva , tvavg , tvb , tvc ,  &
               & xday , xmsf , xtm1
      integer :: i , icons , iptn , itr , j , k , lev , n
      integer :: jm1, jp1
#ifdef MPP1
      integer :: ierr , icons_mpi , numrec
#endif
      character (len=50) :: subroutine_name='tend'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!-----fill up the boundary slices:
!
!     if (iexec == 1) then
      if ( .not. ifrest .and. iexec == 1 ) then
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
        do k = 1 , kz
          do i = 1 , iy
            atm1%u(i,k,j) = atm1%u(i,k,j)*mddom%msfd(i,j)
            atm1%v(i,k,j) = atm1%v(i,k,j)*mddom%msfd(i,j)
          end do
        end do
      end do
#ifdef MPP1
      call mpi_sendrecv(sps1%ps(1,jxp),iy,mpi_real8,ieast,1,            &
                      & sps1%ps(1,0),iy,mpi_real8,iwest,1,              &
                      & mpi_comm_world,mpi_status_ignore,ierr)
      call mpi_sendrecv(sps1%ps(1,1),iy,mpi_real8,iwest,2,              &
                      & sps1%ps(1,jxp+1),iy,mpi_real8,ieast,2,          &
                      & mpi_comm_world,mpi_status_ignore,ierr)
#endif
!
!-----decouple u, v, t, qv, and qc
!
#ifndef BAND
#ifdef MPP1
      do j = 1 , jendl
        if ( myid == 0 .and. j == 1 ) then
#else
      do j = 1 , jx
        if ( j == 1 ) then
#endif
!-----------lateral slices:
!-----------west boundary:
          do k = 1 , kz
            do i = 1 , iy
              atmx%u(i,k,j) = uj1(i,k)
              atmx%v(i,k,j) = vj1(i,k)
            end do
          end do
          if ( iboudy == 3 .or. iboudy == 4 ) then
!..............inflow/outflow dependence:
            do k = 1 , kz
              do i = 1 , iy
                if ( atmx%u(i,k,j) < d_zero ) then
                  atmx%v(i,k,j) = vj2(i,k)
                  atmx%u(i,k,j) = uj2(i,k)
                end if
              end do
            end do
          end if
#ifdef MPP1
        else if ( myid == 0 .and. j == 2 ) then
#else
        else if ( j == 2 ) then
#endif
          do k = 1 , kz
            do i = 1 , iy
              atmx%u(i,k,j) = uj2(i,k)
              atmx%v(i,k,j) = vj2(i,k)
            end do
          end do
          if ( iboudy == 3 .or. iboudy == 4 ) then
!..............inflow/outflow dependence:
            do k = 1 , kz
!.................south boundary:
              if ( atmx%v(1,k,j) < d_zero ) then
                atmx%v(1,k,j) = atmx%v(2,k,j)
                atmx%u(1,k,j) = atmx%u(2,k,j)
              end if
!.................north boundary:
              if ( atmx%v(iy,k,j) >= d_zero ) then
                atmx%v(iy,k,j) = atmx%v(iym1,k,j)
                atmx%u(iy,k,j) = atmx%u(iym1,k,j)
              end if
            end do
          end if
#ifdef MPP1
        else if ( myid == nproc-1 .and. j == jendl-1 ) then
#else
        else if ( j == jxm1 ) then
#endif
          do k = 1 , kz
            do i = 1 , iy
              atmx%u(i,k,j) = ujlx(i,k)
              atmx%v(i,k,j) = vjlx(i,k)
            end do
          end do
          if ( iboudy == 3 .or. iboudy == 4 ) then
!..............inflow/outflow dependence:
            do k = 1 , kz
!.................south boundary:
              if ( atmx%v(1,k,j) < d_zero ) then
                atmx%v(1,k,j) = atmx%v(2,k,j)
                atmx%u(1,k,j) = atmx%u(2,k,j)
              end if
!.................north boundary:
              if ( atmx%v(iy,k,j) >= d_zero ) then
                atmx%v(iy,k,j) = atmx%v(iym1,k,j)
                atmx%u(iy,k,j) = atmx%u(iym1,k,j)
              end if
            end do
          end if
#ifdef MPP1
        else if ( myid == nproc-1 .and. j == jendl ) then
#else
        else if ( j == jx ) then
#endif
!-----------east boundary:
!
!...........no inflow/outflow dependence:
!
          do k = 1 , kz
            do i = 1 , iy
              atmx%u(i,k,j) = ujl(i,k)
              atmx%v(i,k,j) = vjl(i,k)
            end do
          end do
          if ( iboudy == 3 .or. iboudy == 4 ) then
!..............inflow/outflow dependence:
            do k = 1 , kz
              do i = 1 , iy
                if ( atmx%u(i,k,j) >= d_zero ) then
                  atmx%v(i,k,j) = vjlx(i,k)
                  atmx%u(i,k,j) = ujlx(i,k)
                end if
              end do
            end do
          end if
        end if
      end do
#endif

#ifdef MPP1
      do j = 1 , jendl
#else
      do j = 1 , jx
#endif
        jm1 = j-1
#if defined(BAND) && (!defined(MPP1))
        if(jm1 == 0) jm1=jx
#endif
!
!-----interior slice:
!-----interior points:
!
#ifndef BAND
#ifdef MPP1
        if((.not.(myid == 0 .and. (j == 1 .or. j == 2))) .and. &
           (.not.(myid == nproc-1 .and. &
           (j == jendl-1 .or. j == jendl))) ) then
#else
        if(.not.(j == 1 .or. j == 2 .or. j == jxm1 .or. j == jx) ) then
#endif
#endif
          do k = 1 , kz
            do i = 3 , iym2
              psabar=(sps1%ps(i,j)+sps1%ps(i,jm1)+ &
                      sps1%ps(i-1,j)+sps1%ps(i-1,jm1))/d_four
              xmsf = mddom%msfd(i,j)
              atmx%u(i,k,j) = atm1%u(i,k,j)/(psabar*xmsf)
              atmx%v(i,k,j) = atm1%v(i,k,j)/(psabar*xmsf)
            end do
          end do
!
!-----------north/south boundary points:
!...........no inflow/outflow dependence:
!
          do k = 1 , kz
!..............for i=2 and i=iym1:
            atmx%u(2,k,j) = ui2(k,j)
            atmx%u(iym1,k,j) = uilx(k,j)
            atmx%v(2,k,j) = vi2(k,j)
            atmx%v(iym1,k,j) = vilx(k,j)
!..............for i=1 and i=ix:
            atmx%u(1,k,j) = ui1(k,j)
            atmx%u(iy,k,j) = uil(k,j)
            atmx%v(1,k,j) = vi1(k,j)
            atmx%v(iy,k,j) = vil(k,j)
          end do
          if ( iboudy == 3 .or. iboudy == 4 ) then
!..............inflow/outflow dependence:
            do k = 1 , kz
!.................south boundary:
              if ( atmx%v(1,k,j) < d_zero ) then
                atmx%v(1,k,j) = atmx%v(2,k,j)
                atmx%u(1,k,j) = atmx%u(2,k,j)
              end if
!.................north boundary:
              if ( atmx%v(iy,k,j) >= d_zero ) then
                atmx%v(iy,k,j) = atmx%v(iym1,k,j)
                atmx%u(iy,k,j) = atmx%u(iym1,k,j)
              end if
            end do
          end if
#ifndef BAND
        end if
#endif
      end do
!
#ifdef MPP1
      do j = 1 , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 1 , jxm1
#endif
#endif
        do k = 1 , kz
          do i = 1 , iym1
            atmx%t(i,k,j) = atm1%t(i,k,j)/sps1%ps(i,j)
            atmx%qv(i,k,j) = atm1%qv(i,k,j)/sps1%ps(i,j)
            atmx%qc(i,k,j) = atm1%qc(i,k,j)/sps1%ps(i,j)
          end do
        end do
      end do
!chem2
      if ( ichem == 1 ) then
!
!-----call special tracer decoupling routine for multiple (ntr) species
!
        do n = 1 , ntr
#ifdef MPP1
          do j = 1 , jendx
#else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
#endif
            do k = 1 , kz
              do i = 1 , iym1
                chi(i,k,j,n) = chia(i,k,j,n)/sps1%ps(i,j)
              end do
            end do
          end do
        end do
      end if
!
!=======================================================================
#ifdef MPP1
#ifndef BAND
      if ( myid /= nproc-1 ) then
#endif
        do i = 1 , iy
          tvar1snd(i,1) = sps2%ps(i,jxp)
        end do
        do k = 1 , kz
          do i = 1 , iy
            tvar1snd(i,1+k) = atm2%t(i,k,jxp)
            tvar1snd(i,1+kz+k) = atm2%qv(i,k,jxp)
            tvar1snd(i,1+kz*2+k) = atm2%u(i,k,jxp)
            tvar1snd(i,1+kz*3+k) = atm2%v(i,k,jxp)
            tvar1snd(i,1+kz*4+k) = atmx%u(i,k,jxp)
            tvar1snd(i,1+kz*5+k) = atmx%v(i,k,jxp)
            tvar1snd(i,1+kz*6+k) = atmx%t(i,k,jxp)
            tvar1snd(i,1+kz*7+k) = atmx%qv(i,k,jxp)
            tvar1snd(i,1+kz*8+k) = atmx%qc(i,k,jxp)
            tvar1snd(i,1+kz*9+k) = atm1%u(i,k,jxp)
            tvar1snd(i,1+kz*10+k) = atm1%v(i,k,jxp)
          end do
        end do
        if ( ichem == 1 ) then
          do n = 1 , ntr
            do k = 1 , kz
              do i = 1 , iy
                tvar1snd(i,kz*11+1+(n-1)*2*kz+k) = chi(i,k,jxp,n)
                tvar1snd(i,kz*11+1+(n-1)*2*kz+kz+k) = chib(i,k,jxp,n)
              end do
            end do
          end do
        end if
#ifndef BAND
      end if
#endif
      numrec = kz*11 + 1
      if ( ichem == 1 ) numrec = kz*11 + 1 + ntr*2 *kz
      call mpi_sendrecv(tvar1snd,iy*numrec,mpi_real8,ieast,1, &
                      & tvar1rcv,iy*numrec,mpi_real8,iwest,1, &
                      & mpi_comm_world,mpi_status_ignore,ierr)
#ifndef BAND
      if ( myid /= 0 ) then
#endif
        do i = 1 , iy
          sps2%ps(i,0) = tvar1rcv(i,1)
        end do
        do k = 1 , kz
          do i = 1 , iy
            atm2%t(i,k,0) = tvar1rcv(i,1+k)
            atm2%qv(i,k,0) = tvar1rcv(i,1+kz+k)
            atm2%u(i,k,0) = tvar1rcv(i,1+kz*2+k)
            atm2%v(i,k,0) = tvar1rcv(i,1+kz*3+k)
            atmx%u(i,k,0) = tvar1rcv(i,1+kz*4+k)
            atmx%v(i,k,0) = tvar1rcv(i,1+kz*5+k)
            atmx%t(i,k,0) = tvar1rcv(i,1+kz*6+k)
            atmx%qv(i,k,0) = tvar1rcv(i,1+kz*7+k)
            atmx%qc(i,k,0) = tvar1rcv(i,1+kz*8+k)
            atm1%u(i,k,0) = tvar1rcv(i,1+kz*9+k)
            atm1%v(i,k,0) = tvar1rcv(i,1+kz*10+k)
          end do
        end do
        if ( ichem == 1 ) then
          do n = 1 , ntr
            do k = 1 , kz
              do i = 1 , iy
                chi(i,k,0,n) = tvar1rcv(i,kz*11+1+(n-1)*2*kz+k)
                chib(i,k,0,n) = tvar1rcv(i,kz*11+1+(n-1)*2*kz+kz+k)
              end do
            end do
          end do
        end if
#ifndef BAND
      end if
      if ( myid /= 0 ) then
#endif
!
        do i = 1 , iy
          tvar1snd(i,1) = sps2%ps(i,1)
        end do
        do k = 1 , kz
          do i = 1 , iy
            tvar1snd(i,1+k) = atm2%t(i,k,1)
            tvar1snd(i,1+kz+k) = atm2%qv(i,k,1)
            tvar1snd(i,1+kz*2+k) = atm2%u(i,k,1)
            tvar1snd(i,1+kz*3+k) = atm2%v(i,k,1)
            tvar1snd(i,1+kz*4+k) = atmx%u(i,k,1)
            tvar1snd(i,1+kz*5+k) = atmx%v(i,k,1)
            tvar1snd(i,1+kz*6+k) = atmx%t(i,k,1)
            tvar1snd(i,1+kz*7+k) = atmx%qv(i,k,1)
            tvar1snd(i,1+kz*8+k) = atmx%qc(i,k,1)
            tvar1snd(i,1+kz*9+k) = atm1%u(i,k,1)
            tvar1snd(i,1+kz*10+k) = atm1%v(i,k,1)
          end do
        end do
        if ( ichem == 1 ) then
          do n = 1 , ntr
            do k = 1 , kz
              do i = 1 , iy
                tvar1snd(i,kz*11+1+(n-1)*kz*2+k) = chi(i,k,1,n)
                tvar1snd(i,kz*11+1+(n-1)*kz*2+kz+k) = chib(i,k,1,n)
              end do
            end do
          end do
        end if
#ifndef BAND
      end if
#endif
      numrec = kz*11 + 1
      if ( ichem == 1 ) numrec = kz*11 + 1 + ntr*kz*2
      call mpi_sendrecv(tvar1snd,iy*numrec,mpi_real8,iwest,2, &
                      & tvar1rcv,iy*numrec,mpi_real8,ieast,2, &
                      & mpi_comm_world,mpi_status_ignore,ierr)
#ifndef BAND
      if ( myid /= nproc-1 ) then
#endif
        do i = 1 , iy
          sps2%ps(i,jxp+1) = tvar1rcv(i,1)
        end do
        do k = 1 , kz
          do i = 1 , iy
            atm2%t(i,k,jxp+1) = tvar1rcv(i,1+k)
            atm2%qv(i,k,jxp+1) = tvar1rcv(i,1+kz+k)
            atm2%u(i,k,jxp+1) = tvar1rcv(i,1+kz*2+k)
            atm2%v(i,k,jxp+1) = tvar1rcv(i,1+kz*3+k)
            atmx%u(i,k,jxp+1) = tvar1rcv(i,1+kz*4+k)
            atmx%v(i,k,jxp+1) = tvar1rcv(i,1+kz*5+k)
            atmx%t(i,k,jxp+1) = tvar1rcv(i,1+kz*6+k)
            atmx%qv(i,k,jxp+1) = tvar1rcv(i,1+kz*7+k)
            atmx%qc(i,k,jxp+1) = tvar1rcv(i,1+kz*8+k)
            atm1%u(i,k,jxp+1) = tvar1rcv(i,1+kz*9+k)
            atm1%v(i,k,jxp+1) = tvar1rcv(i,1+kz*10+k)
          end do
        end do
        if ( ichem == 1 ) then
          do n = 1 , ntr
            do k = 1 , kz
              do i = 1 , iy
                chi(i,k,jxp+1,n) = tvar1rcv(i,kz*11+1+(n-1)*kz*2+k)
                chib(i,k,jxp+1,n) = tvar1rcv(i,kz*11+1+(n-1)*kz*2+kz+k)
              end do
            end do
          end do
        end if
#ifndef BAND
      end if
#endif
#endif
!=======================================================================
!
!=======================================================================
!
!-----interior points:
!
#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif
         jm1 = j-1
#if defined(BAND) && (!defined(MPP1))
         if(jm1 == 0) jm1=jx
#endif
        do i = 2 , iym1
          sps2%pdot(i,j)=(sps2%ps(i,j)+sps2%ps(i-1,j)+ &
                          sps2%ps(i,jm1)+sps2%ps(i-1,jm1))/d_four
        end do
      end do
!
!-----east and west boundaries:
!
#ifndef BAND
      do i = 2 , iym1
#ifdef MPP1
        if ( myid == 0 )  &
          sps2%pdot(i,1) = (sps2%ps(i,1)+sps2%ps(i-1,1))/d_two
        if ( myid == nproc-1 ) &
          sps2%pdot(i,jendl) = (sps2%ps(i,jendx)+ &
                                sps2%ps(i-1,jendx))/d_two
#else
        sps2%pdot(i,1) = (sps2%ps(i,1)+sps2%ps(i-1,1))/d_two
        sps2%pdot(i,jx) = (sps2%ps(i,jxm1)+sps2%ps(i-1,jxm1))/d_two
#endif
      end do
#endif
!
!-----north and south boundaries:
!
#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif
         jm1 = j-1
#if defined(BAND) && (!defined(MPP1))
         if(jm1 == 0) jm1=jx
#endif
        sps2%pdot(1,j)  = (sps2%ps(1,j)+sps2%ps(1,jm1))/d_two
        sps2%pdot(iy,j) = (sps2%ps(iym1,j)+sps2%ps(iym1,jm1))/d_two
      end do
!
!-----corner points:
!
#ifndef BAND
#ifdef MPP1
      if ( myid == 0 ) then
        sps2%pdot(1,1) = sps2%ps(1,1)
        sps2%pdot(iy,1) = sps2%ps(iym1,1)
      end if
      if ( myid == nproc-1 ) then
        sps2%pdot(1,jendl) = sps2%ps(1,jendx)
        sps2%pdot(iy,jendl) = sps2%ps(iym1,jendx)
      end if
#else
      sps2%pdot(1,1) = sps2%ps(1,1)
      sps2%pdot(iy,1) = sps2%ps(iym1,1)
      sps2%pdot(1,jx) = sps2%ps(1,jxm1)
      sps2%pdot(iy,jx) = sps2%ps(iym1,jxm1)
#endif
#endif
!
!=======================================================================
!
      call slice

#ifdef CLM
      if ( init_grid ) then
        call initclm
        init_grid = .false.
      end if
#endif
!
!=======================================================================
!
#ifdef MPP1
#ifndef BAND
      if ( myid /= nproc-1 ) then
#endif
        do k = 1 , kz
          do i = 1 , iy
            var2snd(i,+k) = ubd3d(i,k,jxp-1)
            var2snd(i,kz+k) = ubd3d(i,k,jxp)
            var2snd(i,kz*2+k) = vbd3d(i,k,jxp-1)
            var2snd(i,kz*3+k) = vbd3d(i,k,jxp)
            var2snd(i,kz*4+k) = tb3d(i,k,jxp-1)
            var2snd(i,kz*5+k) = tb3d(i,k,jxp)
            var2snd(i,kz*6+k) = qvb3d(i,k,jxp-1)
            var2snd(i,kz*7+k) = qvb3d(i,k,jxp)
            var2snd(i,kz*8+k) = qcb3d(i,k,jxp-1)
            var2snd(i,kz*9+k) = qcb3d(i,k,jxp)
          end do
        end do
        if ( ichem == 1 ) then
          do n = 1 , ntr
            do k = 1 , kz
              do i = 1 , iy
                var2snd(i,kz*10+(n-1)*2*kz+k) = chib3d(i,k,jxp-1,n)
                var2snd(i,kz*10+(n-1)*2*kz+kz+k) = chib3d(i,k,jxp,n)
              end do
            end do
          end do
        end if
#ifndef BAND
      end if
#endif
      numrec = kz*5*2
      if ( ichem == 1 ) numrec = kz*(ntr+5)*2
      call mpi_sendrecv(var2snd(1,1),iy*numrec,mpi_real8,               &
                      & ieast,1,var2rcv(1,1),iy*numrec,                 &
                      & mpi_real8,iwest,1,mpi_comm_world,               &
                      & mpi_status_ignore,ierr)
#ifndef BAND
      if ( myid /= 0 ) then
#endif
        do k = 1 , kz
          do i = 1 , iy
            ubd3d(i,k,-1) = var2rcv(i,+k)
            ubd3d(i,k,0) = var2rcv(i,kz+k)
            vbd3d(i,k,-1) = var2rcv(i,kz*2+k)
            vbd3d(i,k,0) = var2rcv(i,kz*3+k)
            tb3d(i,k,-1) = var2rcv(i,kz*4+k)
            tb3d(i,k,0) = var2rcv(i,kz*5+k)
            qvb3d(i,k,-1) = var2rcv(i,kz*6+k)
            qvb3d(i,k,0) = var2rcv(i,kz*7+k)
            qcb3d(i,k,-1) = var2rcv(i,kz*8+k)
            qcb3d(i,k,0) = var2rcv(i,kz*9+k)
          end do
        end do
        if ( ichem == 1 ) then
          do n = 1 , ntr
            do k = 1 , kz
              do i = 1 , iy
                chib3d(i,k,-1,n) = var2rcv(i,kz*10+(n-1)*2*kz+k)
                chib3d(i,k,0,n) = var2rcv(i,kz*10+(n-1)*2*kz+kz+k)
              end do
            end do
          end do
        end if
#ifndef BAND
      end if
!
      if ( myid /= 0 ) then
#endif
        do k = 1 , kz
          do i = 1 , iy
            var2snd(i,+k) = ubd3d(i,k,1)
            var2snd(i,kz+k) = ubd3d(i,k,2)
            var2snd(i,kz*2+k) = vbd3d(i,k,1)
            var2snd(i,kz*3+k) = vbd3d(i,k,2)
            var2snd(i,kz*4+k) = tb3d(i,k,1)
            var2snd(i,kz*5+k) = tb3d(i,k,2)
            var2snd(i,kz*6+k) = qvb3d(i,k,1)
            var2snd(i,kz*7+k) = qvb3d(i,k,2)
            var2snd(i,kz*8+k) = qcb3d(i,k,1)
            var2snd(i,kz*9+k) = qcb3d(i,k,2)
          end do
        end do
        if ( ichem == 1 ) then
          do n = 1 , ntr
            do k = 1 , kz
              do i = 1 , iy
                var2snd(i,kz*10+(n-1)*2*kz+k) = chib3d(i,k,1,n)
                var2snd(i,kz*10+(n-1)*2*kz+kz+k) = chib3d(i,k,2,n)
              end do
            end do
          end do
        end if
#ifndef BAND
      end if
#endif
      numrec = kz*5*2
      if ( ichem == 1 ) numrec = kz*(ntr+5)*2
      call mpi_sendrecv(var2snd(1,1),iy*numrec,mpi_real8,               &
                      & iwest,2,var2rcv(1,1),iy*numrec,                 &
                      & mpi_real8,ieast,2,mpi_comm_world,               &
                      & mpi_status_ignore,ierr)
#ifndef BAND
      if ( myid /= nproc-1 ) then
#endif
        do k = 1 , kz
          do i = 1 , iy
            ubd3d(i,k,jxp+1) = var2rcv(i,+k)
            ubd3d(i,k,jxp+2) = var2rcv(i,kz+k)
            vbd3d(i,k,jxp+1) = var2rcv(i,kz*2+k)
            vbd3d(i,k,jxp+2) = var2rcv(i,kz*3+k)
            tb3d(i,k,jxp+1) = var2rcv(i,kz*4+k)
            tb3d(i,k,jxp+2) = var2rcv(i,kz*5+k)
            qvb3d(i,k,jxp+1) = var2rcv(i,kz*6+k)
            qvb3d(i,k,jxp+2) = var2rcv(i,kz*7+k)
            qcb3d(i,k,jxp+1) = var2rcv(i,kz*8+k)
            qcb3d(i,k,jxp+2) = var2rcv(i,kz*9+k)
          end do
        end do
        if ( ichem == 1 ) then
          do n = 1 , ntr
            do k = 1 , kz
              do i = 1 , iy
                chib3d(i,k,jxp+1,n) = var2rcv(i,kz*10+(n-1)*2*kz+k)
                chib3d(i,k,jxp+2,n) = var2rcv(i,kz*10+(n-1)*2*kz+kz+k)
              end do
            end do
          end do
        end if
#ifndef BAND
      end if
#endif
#endif
!
!**********************************************************************
!***** "j" loop begins here:
!
#ifndef BAND
#ifdef MPP1
      do j = 1 , jendx
        if ( (myid == 0 .and. j == 1) .or.                              &
           & (myid == nproc-1 .and. j == jendx) ) then
#else
      do j = 1 , jxm1
        if ( j == 1 .or. j == jxm1 ) then
#endif
          icon(j) = 0
          do k = 1 , kzp1
            do i = 1 , iym1
              qdot(i,k,j) = d_zero
            end do
          end do
          do k = 1 , kz
            do i = 1 , iym1
              omega(i,k,j) = d_zero
            end do
          end do
        end if
      end do
#endif

#ifdef MPP1
      do j = 1 , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 1 , jxm1
#endif
#endif
        jp1 = j+1
#if defined(BAND) && (!defined(MPP1))
        if(jp1 == jx+1) jp1 = 1
#endif
!
#ifndef BAND
#ifdef MPP1
        if((.not.(myid == 0 .and. j == 1)) .and. &
           (.not.(myid == nproc-1 .and. j == jendx)) ) then
#else
        if( .not.(j == 1 .or. j == jxm1) ) then
#endif
#endif
        icon(j) = 0
!
!----------------------------------------------------------------------
!**p**compute the pressure tendency:
!
        do i = 2 , iym2
           pten(i,j) = d_zero
        end do
        do k = 1 , kz
           do i = 2 , iym2
              divl(i,k) = (atm1%u(i+1,k,jp1)+atm1%u(i,k,jp1)- &
                           atm1%u(i+1,k,j)-atm1%u(i,k,j))+    &
                          (atm1%v(i+1,k,jp1)+atm1%v(i+1,k,j)- &
                           atm1%v(i,k,jp1)-atm1%v(i,k,j))
              pten(i,j) = pten(i,j) - divl(i,k)*dsigma(k)     &
                        & /(dx2*mddom%msfx(i,j)*mddom%msfx(i,j))
           end do
        end do
!
!..p..compute vertical sigma-velocity (qdot):
!
        do k = 1 , kzp1
           do i = 1 , iym1
              qdot(i,k,j) = d_zero
           end do
        end do
        do k = 2 , kz
           do i = 2 , iym2
              qdot(i,k,j) = qdot(i,k-1,j)                               &
                       & - (pten(i,j)+divl(i,k-1)/(dx2*mddom%msfx(i,j) &
                       & *mddom%msfx(i,j)))*dsigma(k-1)/sps1%ps(i,j)
           end do
        end do
#ifndef BAND
        end if
#endif
      end do
#ifdef MPP1
      call mpi_sendrecv(qdot(1,1,jxp),iy*kzp1,mpi_real8,                &
                      & ieast,1,qdot(1,1,0),iy*kzp1,                    &
                      & mpi_real8,iwest,1,mpi_comm_world,               &
                      & mpi_status_ignore,ierr)
      call mpi_sendrecv(qdot(1,1,1),iy*kzp1,mpi_real8,iwest,            &
                      & 2,qdot(1,1,jxp+1),iy*kzp1,mpi_real8,            &
                      & ieast,2,mpi_comm_world,mpi_status_ignore,ierr)
#endif
!
!..p..compute omega
!
#ifdef MPP1
      do j = 1 , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 1 , jxm1
#endif
#endif
        jp1 = j+1
        jm1 = j-1
#if defined(BAND) && (!defined(MPP1))
        if(jp1 == jx+1) jp1 = 1
        if(jm1 == 0) jm1=jx
#endif
#ifndef BAND
#ifdef MPP1
        if((.not.(myid == 0 .and. j == 1)) .and. &
           (.not.(myid == nproc-1 .and. j == jendx)) ) then
#else
        if( .not.(j == 1 .or. j == jxm1) ) then
#endif
#endif
        do k = 1 , kz
           do i = 2 , iym2
              omega(i,k,j) = (sps1%ps(i,j)/d_two)* &
                       & (qdot(i,k+1,j)+qdot(i,k,j))+a(k)*(pten(i,j)+   &
                       & ((atmx%u(i,k,j)+atmx%u(i+1,k,j)+               &
                       &   atmx%u(i+1,k,jp1)+atmx%u(i,k,jp1))*          &
                       & (sps1%ps(i,jp1)-sps1%ps(i,jm1))+               &
                       & (atmx%v(i,k,j)+atmx%v(i+1,k,j)+                &
                       &  atmx%v(i+1,k,jp1)+atmx%v(i,k,jp1))*           &
                       & (sps1%ps(i+1,j)-sps1%ps(i-1,j)))/              &
                       & (dx8*mddom%msfx(i,j)))
           end do
        end do
#ifndef BAND
        end if
#endif
      end do

#ifdef MPP1
#ifndef BAND
      if ( nspgx >= jxp ) then
        do i = 1 , iy
          bdyewsnd(i,1) = peb(i,1)
          bdyewsnd(i,2) = pwb(i,jxp)
          bdyewsnd(i,3) = pebt(i,1)
          bdyewsnd(i,4) = pwbt(i,jxp)
        end do
        do k = 1 , kz
          do i = 1 , iy
            bdyewsnd(i,4+k) = teb(i,k,1)
            bdyewsnd(i,4+kz+k) = twb(i,k,jxp)
            bdyewsnd(i,4+kz*2+k) = tebt(i,k,1)
            bdyewsnd(i,4+kz*3+k) = twbt(i,k,jxp)
            bdyewsnd(i,4+kz*4+k) = qeb(i,k,1)
            bdyewsnd(i,4+kz*5+k) = qwb(i,k,jxp)
            bdyewsnd(i,4+kz*6+k) = qebt(i,k,1)
            bdyewsnd(i,4+kz*7+k) = qwbt(i,k,jxp)
            bdyewsnd(i,4+kz*8+k) = ueb(i,k,1)
            bdyewsnd(i,4+kz*9+k) = uwb(i,k,jxp)
            bdyewsnd(i,4+kz*10+k) = uebt(i,k,1)
            bdyewsnd(i,4+kz*11+k) = uwbt(i,k,jxp)
            bdyewsnd(i,4+kz*12+k) = veb(i,k,1)
            bdyewsnd(i,4+kz*13+k) = vwb(i,k,jxp)
            bdyewsnd(i,4+kz*14+k) = vebt(i,k,1)
            bdyewsnd(i,4+kz*15+k) = vwbt(i,k,jxp)
          end do
        end do
        call mpi_sendrecv(bdyewsnd(1,1),iy*(kz*16+4),                &
                        & mpi_real8,ieast,1,bdyewrcv(1,1),           &
                        & iy*(kz*16+4),mpi_real8,iwest,1,            &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        do i = 1 , iy
          if ( myid == nproc-1 ) then
            peb(i,jendl) = bdyewrcv(i,1)
            pebt(i,jendl) = bdyewrcv(i,3)
          else
            peb(i,jxp+1) = bdyewrcv(i,1)
            pebt(i,jxp+1) = bdyewrcv(i,3)
          end if
          pwb(i,0) = bdyewrcv(i,2)
          pwbt(i,0) = bdyewrcv(i,4)
        end do
        do k = 1 , kz
          do i = 1 , iy
            if ( myid == nproc-1 ) then
              teb(i,k,jendl) = bdyewrcv(i,4+k)
              tebt(i,k,jendl) = bdyewrcv(i,4+kz*2+k)
              qeb(i,k,jendl) = bdyewrcv(i,4+kz*4+k)
              qebt(i,k,jendl) = bdyewrcv(i,4+kz*6+k)
            else
              teb(i,k,jxp+1) = bdyewrcv(i,4+k)
              tebt(i,k,jxp+1) = bdyewrcv(i,4+kz*2+k)
              qeb(i,k,jxp+1) = bdyewrcv(i,4+kz*4+k)
              qebt(i,k,jxp+1) = bdyewrcv(i,4+kz*6+k)
            end if
            ueb(i,k,jxp+1) = bdyewrcv(i,4+kz*8+k)
            uebt(i,k,jxp+1) = bdyewrcv(i,4+kz*10+k)
            veb(i,k,jxp+1) = bdyewrcv(i,4+kz*12+k)
            vebt(i,k,jxp+1) = bdyewrcv(i,4+kz*14+k)
            twb(i,k,0) = bdyewrcv(i,4+kz+k)
            twbt(i,k,0) = bdyewrcv(i,4+kz*3+k)
            qwb(i,k,0) = bdyewrcv(i,4+kz*5+k)
            qwbt(i,k,0) = bdyewrcv(i,4+kz*7+k)
            uwb(i,k,0) = bdyewrcv(i,4+kz*9+k)
            uwbt(i,k,0) = bdyewrcv(i,4+kz*11+k)
            vwb(i,k,0) = bdyewrcv(i,4+kz*13+k)
            vwbt(i,k,0) = bdyewrcv(i,4+kz*15+k)
          end do
        end do
        do i = 1 , iy
          if ( myid == nproc-1 ) then
            bdyewsnd(i,1) = peb(i,jendx)
            bdyewsnd(i,3) = pebt(i,jendx)
          else
            bdyewsnd(i,1) = peb(i,jxp)
            bdyewsnd(i,3) = pebt(i,jxp)
          end if
          bdyewsnd(i,2) = pwb(i,1)
          bdyewsnd(i,4) = pwbt(i,1)
        end do
        do k = 1 , kz
          do i = 1 , iy
            if ( myid == nproc-1 ) then
              bdyewsnd(i,4+k) = teb(i,k,jendx)
              bdyewsnd(i,4+kz*2+k) = tebt(i,k,jendx)
              bdyewsnd(i,4+kz*4+k) = qeb(i,k,jendx)
              bdyewsnd(i,4+kz*6+k) = qebt(i,k,jendx)
            else
              bdyewsnd(i,4+k) = teb(i,k,jxp)
              bdyewsnd(i,4+kz*2+k) = tebt(i,k,jxp)
              bdyewsnd(i,4+kz*4+k) = qeb(i,k,jxp)
              bdyewsnd(i,4+kz*6+k) = qebt(i,k,jxp)
            end if
            bdyewsnd(i,4+kz*8+k) = ueb(i,k,jxp)
            bdyewsnd(i,4+kz*10+k) = uebt(i,k,jxp)
            bdyewsnd(i,4+kz*12+k) = veb(i,k,jxp)
            bdyewsnd(i,4+kz*14+k) = vebt(i,k,jxp)
            bdyewsnd(i,4+kz+k) = twb(i,k,1)
            bdyewsnd(i,4+kz*3+k) = twbt(i,k,1)
            bdyewsnd(i,4+kz*5+k) = qwb(i,k,1)
            bdyewsnd(i,4+kz*7+k) = qwbt(i,k,1)
            bdyewsnd(i,4+kz*9+k) = uwb(i,k,1)
            bdyewsnd(i,4+kz*11+k) = uwbt(i,k,1)
            bdyewsnd(i,4+kz*13+k) = vwb(i,k,1)
            bdyewsnd(i,4+kz*15+k) = vwbt(i,k,1)
          end do
        end do
        call mpi_sendrecv(bdyewsnd(1,1),iy*(kz*16+4),                &
                        & mpi_real8,iwest,2,bdyewrcv(1,1),           &
                        & iy*(kz*16+4),mpi_real8,ieast,2,            &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        do i = 1 , iy
          peb(i,0) = bdyewrcv(i,1)
          pebt(i,0) = bdyewrcv(i,3)
          pwb(i,jxp+1) = bdyewrcv(i,2)
          pwbt(i,jxp+1) = bdyewrcv(i,4)
        end do
        do k = 1 , kz
          do i = 1 , iy
            teb(i,k,0) = bdyewrcv(i,4+k)
            twb(i,k,jxp+1) = bdyewrcv(i,4+kz+k)
            tebt(i,k,0) = bdyewrcv(i,4+kz*2+k)
            twbt(i,k,jxp+1) = bdyewrcv(i,4+kz*3+k)
            qeb(i,k,0) = bdyewrcv(i,4+kz*4+k)
            qwb(i,k,jxp+1) = bdyewrcv(i,4+kz*5+k)
            qebt(i,k,0) = bdyewrcv(i,4+kz*6+k)
            qwbt(i,k,jxp+1) = bdyewrcv(i,4+kz*7+k)
            ueb(i,k,0) = bdyewrcv(i,4+kz*8+k)
            uwb(i,k,jxp+1) = bdyewrcv(i,4+kz*9+k)
            uebt(i,k,0) = bdyewrcv(i,4+kz*10+k)
            uwbt(i,k,jxp+1) = bdyewrcv(i,4+kz*11+k)
            veb(i,k,0) = bdyewrcv(i,4+kz*12+k)
            vwb(i,k,jxp+1) = bdyewrcv(i,4+kz*13+k)
            vebt(i,k,0) = bdyewrcv(i,4+kz*14+k)
            vwbt(i,k,jxp+1) = bdyewrcv(i,4+kz*15+k)
          end do
        end do
      end if
#endif
!
#ifndef BAND
      if ( myid /= nproc-1 ) then
#endif
        do i = 1 , nspgx
          bdynssnd(i,1) = pnb(i,jxp)
          bdynssnd(i,2) = pnbt(i,jxp)
          bdynssnd(i,3) = pss(i,jxp)
          bdynssnd(i,4) = psbt(i,jxp)
        end do
        do k = 1 , kz
          do i = 1 , nspgx
            bdynssnd(i,4+k) = tnb(i,k,jxp)
            bdynssnd(i,4+kz+k) = tnbt(i,k,jxp)
            bdynssnd(i,4+kz*2+k) = tsb(i,k,jxp)
            bdynssnd(i,4+kz*3+k) = tsbt(i,k,jxp)
            bdynssnd(i,4+kz*4+k) = qnb(i,k,jxp)
            bdynssnd(i,4+kz*5+k) = qnbt(i,k,jxp)
            bdynssnd(i,4+kz*6+k) = qsb(i,k,jxp)
            bdynssnd(i,4+kz*7+k) = qsbt(i,k,jxp)
            bdynssnd(i,4+kz*8+k) = unb(i,k,jxp)
            bdynssnd(i,4+kz*9+k) = unbt(i,k,jxp)
            bdynssnd(i,4+kz*10+k) = usb(i,k,jxp)
            bdynssnd(i,4+kz*11+k) = usbt(i,k,jxp)
            bdynssnd(i,4+kz*12+k) = vnb(i,k,jxp)
            bdynssnd(i,4+kz*13+k) = vnbt(i,k,jxp)
            bdynssnd(i,4+kz*14+k) = vsb(i,k,jxp)
            bdynssnd(i,4+kz*15+k) = vsbt(i,k,jxp)
          end do
        end do
#ifndef BAND
      end if
#endif
      call mpi_sendrecv(bdynssnd(1,1),nspgx*(kz*16+4),               &
                      & mpi_real8,ieast,1,bdynsrcv(1,1),             &
                      & nspgx*(kz*16+4),mpi_real8,iwest,1,           &
                      & mpi_comm_world,mpi_status_ignore,ierr)
#ifndef BAND
      if ( myid /= 0 ) then
#endif
        do i = 1 , nspgx
          pnb(i,0) = bdynsrcv(i,1)
          pnbt(i,0) = bdynsrcv(i,2)
          pss(i,0) = bdynsrcv(i,3)
          psbt(i,0) = bdynsrcv(i,4)
        end do
        do k = 1 , kz
          do i = 1 , nspgx
            tnb(i,k,0) = bdynsrcv(i,4+k)
            tnbt(i,k,0) = bdynsrcv(i,4+kz+k)
            tsb(i,k,0) = bdynsrcv(i,4+kz*2+k)
            tsbt(i,k,0) = bdynsrcv(i,4+kz*3+k)
            qnb(i,k,0) = bdynsrcv(i,4+kz*4+k)
            qnbt(i,k,0) = bdynsrcv(i,4+kz*5+k)
            qsb(i,k,0) = bdynsrcv(i,4+kz*6+k)
            qsbt(i,k,0) = bdynsrcv(i,4+kz*7+k)
            unb(i,k,0) = bdynsrcv(i,4+kz*8+k)
            unbt(i,k,0) = bdynsrcv(i,4+kz*9+k)
            usb(i,k,0) = bdynsrcv(i,4+kz*10+k)
            usbt(i,k,0) = bdynsrcv(i,4+kz*11+k)
            vnb(i,k,0) = bdynsrcv(i,4+kz*12+k)
            vnbt(i,k,0) = bdynsrcv(i,4+kz*13+k)
            vsb(i,k,0) = bdynsrcv(i,4+kz*14+k)
            vsbt(i,k,0) = bdynsrcv(i,4+kz*15+k)
          end do
        end do
#ifndef BAND
      end if
      if ( myid /= 0 ) then
#endif
        do i = 1 , nspgx
          bdynssnd(i,1) = pnb(i,1)
          bdynssnd(i,2) = pnbt(i,1)
          bdynssnd(i,3) = pss(i,1)
          bdynssnd(i,4) = psbt(i,1)
        end do
        do k = 1 , kz
          do i = 1 , nspgx
            bdynssnd(i,4+k) = tnb(i,k,1)
            bdynssnd(i,4+kz+k) = tnbt(i,k,1)
            bdynssnd(i,4+kz*2+k) = tsb(i,k,1)
            bdynssnd(i,4+kz*3+k) = tsbt(i,k,1)
            bdynssnd(i,4+kz*4+k) = qnb(i,k,1)
            bdynssnd(i,4+kz*5+k) = qnbt(i,k,1)
            bdynssnd(i,4+kz*6+k) = qsb(i,k,1)
            bdynssnd(i,4+kz*7+k) = qsbt(i,k,1)
            bdynssnd(i,4+kz*8+k) = unb(i,k,1)
            bdynssnd(i,4+kz*9+k) = unbt(i,k,1)
            bdynssnd(i,4+kz*10+k) = usb(i,k,1)
            bdynssnd(i,4+kz*11+k) = usbt(i,k,1)
            bdynssnd(i,4+kz*12+k) = vnb(i,k,1)
            bdynssnd(i,4+kz*13+k) = vnbt(i,k,1)
            bdynssnd(i,4+kz*14+k) = vsb(i,k,1)
            bdynssnd(i,4+kz*15+k) = vsbt(i,k,1)
          end do
        end do
#ifndef BAND
      end if
#endif
      call mpi_sendrecv(bdynssnd(1,1),nspgx*(kz*16+4),               &
                      & mpi_real8,iwest,2,bdynsrcv(1,1),             &
                      & nspgx*(kz*16+4),mpi_real8,ieast,2,           &
                      & mpi_comm_world,mpi_status_ignore,ierr)
#ifndef BAND
      if ( myid /= nproc-1 ) then
#endif
        do i = 1 , nspgx
          pnb(i,jxp+1) = bdynsrcv(i,1)
          pnbt(i,jxp+1) = bdynsrcv(i,2)
          pss(i,jxp+1) = bdynsrcv(i,3)
          psbt(i,jxp+1) = bdynsrcv(i,4)
        end do
        do k = 1 , kz
          do i = 1 , nspgx
            tnb(i,k,jxp+1) = bdynsrcv(i,4+k)
            tnbt(i,k,jxp+1) = bdynsrcv(i,4+kz+k)
            tsb(i,k,jxp+1) = bdynsrcv(i,4+kz*2+k)
            tsbt(i,k,jxp+1) = bdynsrcv(i,4+kz*3+k)
            qnb(i,k,jxp+1) = bdynsrcv(i,4+kz*4+k)
            qnbt(i,k,jxp+1) = bdynsrcv(i,4+kz*5+k)
            qsb(i,k,jxp+1) = bdynsrcv(i,4+kz*6+k)
            qsbt(i,k,jxp+1) = bdynsrcv(i,4+kz*7+k)
            unb(i,k,jxp+1) = bdynsrcv(i,4+kz*8+k)
            unbt(i,k,jxp+1) = bdynsrcv(i,4+kz*9+k)
            usb(i,k,jxp+1) = bdynsrcv(i,4+kz*10+k)
            usbt(i,k,jxp+1) = bdynsrcv(i,4+kz*11+k)
            vnb(i,k,jxp+1) = bdynsrcv(i,4+kz*12+k)
            vnbt(i,k,jxp+1) = bdynsrcv(i,4+kz*13+k)
            vsb(i,k,jxp+1) = bdynsrcv(i,4+kz*14+k)
            vsbt(i,k,jxp+1) = bdynsrcv(i,4+kz*15+k)
          end do
        end do
#ifndef BAND
      end if
#endif
#endif

#ifdef MPP1
      do j = jbegin , jendx
#ifndef BAND
        if ( myid /= nproc-1 .or. j /= jendx ) then
#endif
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
        if ( j /= jxm1 ) then
#endif
#endif
!EES      omega change: I broke up the write onto two lines
!         and commented out the if statment
!
          if ( iboudy == 4 ) then
!..p..apply sponge boundary conditions to pten:
            call sponge_p(ispgx,wgtx,pten(1,j),j)
!....apply  the nudging boundary conditions:
          else if ( iboudy == 1 .or. iboudy == 5 ) then
            xtm1 = xtime - dtmin
            if ( dabs(xtime) < 0.00001D0 .and. ldatez > idate0 ) &
              xtm1 = -dtmin
            call nudge_p(ispgx,fnudge,gnudge,xtm1,pten(:,j),j,iboudy)
          end if
#ifndef BAND
        end if     !end if(j /= jxm1) test
#endif
      end do

#ifndef BAND
#ifdef MPP1
      do j = 1 , jendx
        if ( myid == 0 .and. j == 1 ) then
#else
      do j = 1 , jxm1
        if ( j == 1 ) then
#endif
          do i = 1 , iym1
            psc(i,j) = sps2%ps(i,j) + dt*pwbt(i,j)
            psd(i,j) = sps1%ps(i,j)
          end do
#ifdef MPP1
        else if ( myid == nproc-1 .and. j == jendx ) then
#else
        else if ( j == jxm1 ) then
#endif
          do i = 1 , iym1
            psd(i,j) = sps1%ps(i,j)
          end do
        end if
      end do
#endif

#ifdef MPP1
      do j = 1 , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 1 , jxm1
#endif
#endif
#ifndef BAND
#ifdef MPP1
        if((.not.(myid == 0 .and. j == 1)) .and. &
           (.not.(myid == nproc-1 .and. j == jendx)) ) then
#else
        if( .not.(j == 1 .or. j == jxm1) ) then
#endif
#endif
!
!..p..forecast pressure:
!
         do i = 2 , iym2
           psc(i,j) = sps2%ps(i,j) + pten(i,j)*dt
         end do
!
!..p..weighted p* (psd)
!
         do i = 2 , iym2
           psd(i,j) = sps1%ps(i,j)
         end do
!
         psc(1,j) = sps2%ps(1,j) + dt*psbt(1,j)
         psc(iym1,j) = sps2%ps(iym1,j) + dt*pnbt(1,j)
         psd(1,j) = sps1%ps(1,j)
         psd(iym1,j) = sps1%ps(iym1,j)
#ifndef BAND
        end if
#endif
      end do
#ifdef MPP1
      call mpi_sendrecv(psd(1,jxp),iy,mpi_real8,ieast,1,                &
                      & psd(1,0),iy,mpi_real8,iwest,1,                  &
                      & mpi_comm_world,mpi_status_ignore,ierr)
#endif
!
!-----compute bleck (1977) noise parameters:
!
#ifdef MPP1
      do j = 1 , jendl
        do i = 1 , iy
          ps4(i,1,j) = pten(i,j)
          ps4(i,2,j) = psc(i,j)
          ps4(i,3,j) = sps2%ps(i,j)
          ps4(i,4,j) = sps1%ps(i,j)
        end do
      end do
      call mpi_gather(ps4, iy*4*jxp,mpi_real8, &
                      ps_4,iy*4*jxp,mpi_real8, &
                      0,mpi_comm_world,ierr)
#endif
      iptn = 0
      ptntot = d_zero
      pt2tot = d_zero
#ifdef MPP1
      if ( myid == 0 ) then
#ifdef BAND
        do j = 1 , jx
#else
        do j = 2 , jxm2
#endif
          if ( jyear /= jyear0 .or. ktau /= 0 ) then
            do i = 2 , iym2
              iptn = iptn + 1
              ptntot = ptntot + dabs(ps_4(i,1,j))
              pt2tot = pt2tot +                       &
                     & dabs((ps_4(i,2,j)+ps_4(i,3,j)- &
                             d_two*ps_4(i,4,j))/((dt*dt)/d_four))
            end do
          end if
        end do
      end if
      call mpi_bcast(iptn,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(ptntot,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(pt2tot,1,mpi_real8,0,mpi_comm_world,ierr)
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#ifndef BAND
        if ( j /= jxm1 ) then
#endif
          if ( jyear /= jyear0 .or. ktau /= 0 ) then
            do i = 2 , iym2
              iptn = iptn + 1
              ptntot = ptntot + dabs(pten(i,j))
              pt2tot = pt2tot + dabs((psc(i,j)+sps2%ps(i,j)- &
                      d_two*sps1%ps(i,j))/((dt*dt)/d_four))
            end do
          end if
#ifndef BAND
        end if
#endif
      end do
#endif
!
#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif
        jp1 = j+1
#if defined(BAND) && (!defined(MPP1))
        if(jp1 == jx+1) jp1 = 1
#endif
!
!------compute the horizontal diffusion coefficient and stored in xkc:
!       the values are calculated at cross points, but they also used
!       for dot-point variables.
 
        do k = 1 , kz
          do i = 2 , iym1
            dudx = atm2%u(i,k,jp1) + atm2%u(i+1,k,jp1) - &
                   atm2%u(i,k,j)   - atm2%u(i+1,k,j)
            dvdx = atm2%v(i,k,jp1) + atm2%v(i+1,k,jp1) - &
                   atm2%v(i,k,j)   - atm2%v(i+1,k,j)
            dudy = atm2%u(i+1,k,j) + atm2%u(i+1,k,jp1) - &
                   atm2%u(i,k,j)   - atm2%u(i,k,jp1)
            dvdy = atm2%v(i+1,k,j) + atm2%v(i+1,k,jp1) - &
                   atm2%v(i,k,j)   - atm2%v(i,k,jp1)
!fil        cell=(xkhz*domfc%hgfact(i,j)/5.+c200*dsqrt((dudx-dvdy)*(dudx-dvdy)
            cell = (xkhz*domfc%hgfact(i,j)                         &
                 & +c200*dsqrt((dudx-dvdy)*(dudx-dvdy)+(dvdx+dudy) &
                 & *(dvdx+dudy)))
            xkc(i,k,j) = dmin1(cell,xkhmax)
 
          end do
        end do
      end do
!

#ifdef MPP1
      do j = jbegin , jendx
#ifndef BAND
        if ( myid /= nproc-1 .or. j /= jendx ) then
#endif
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
        if ( j /= jxm1 ) then
#endif
#endif
!
!---------------------------------------------------------------------
!**t**compute the temperature tendency:
!
          do k = 1 , kz
            do i = 2 , iym2
              aten%t(i,k,j) = d_zero
              aten%qv(i,k,j) = d_zero
              aten%qc(i,k,j) = d_zero
            end do
          end do
 
!
!..t..compute the horizontal advection term:
!
          call hadv_x(aten%t(:,:,j),atmx%t,dx4,j,1)
!
!..t..compute the vertical advection term:
!
          call vadv(aten%t(1,1,j),qdot,atm1%t(1,1,j),j,1)
!
!..t..compute the adiabatic term:
!
          do k = 1 , kz
            do i = 2 , iym2
              rovcpm = rgas/(cpd*(d_one+0.8D0*(atmx%qv(i,k,j))))
              tv = atmx%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
              aten%t(i,k,j) = aten%t(i,k,j) + (omega(i,k,j)*rovcpm*tv) &
                          & /(r8pt/sps1%ps(i,j)+a(k))
            end do
          end do
!
!..t..compute the diffusion term for t and store in difft:
!
          do k = 1 , kz
            do i = 1 , iym1
              difft(i,k,j) = d_zero
              diffq(i,k,j) = d_zero
            end do
          end do
!
          call diffu_x(difft(:,:,j),tb3d,sps2%ps,xkc(:,:,j),j)
!
!**q**compute the moisture tendencies:
!
!....icup = 1 : kuo-anthes cumulus parameterizaion scheme
!....icup = 2 : grell cumulus paramterization scheme
!....icup = 3 : betts-miller (1986)
!....icup = 4 : emanuel (1991)
!....icup = 99: grell over land, emanuel over ocean
!....icup = 98: emanuel over land, grell over ocean
!
          if ( icup /= 1 ) then
            call hadv_x(aten%qv(:,:,j),atmx%qv,dx4,j,1)
            call vadv(aten%qv(:,:,j),qdot,atm1%qv(:,:,j),j,2)
          end if
 
          if ( icup == 1 ) then
            call cupara(j)
          end if
          if ( icup == 2 .or. icup == 99 .or. icup == 98 ) then
            call cuparan(j)
          end if
          if ( icup == 3 ) then
            write (aline,*)                                             &
                & 'ICTP RegCM team thinks the Betts-Miller code',       &
                & ' is not ready for Regional Climate Run yet.'
            call say
            call fatal(__FILE__,__LINE__,                               &
                      &'BETTS MILLER CUMULUS OPTION NOT ALLOWED')
            call bmpara(aten%t(:,:,j),aten%qv(:,:,j),j)
          end if
          if ( icup == 4 .or. icup == 99 .or. icup == 98 ) then
            call cupemandrv(j)
          end if

          if ( ipptls == 1 ) then
            call hadv_x(aten%qc(:,:,j),atmx%qc,dx4,j,1)
            call vadv(aten%qc(:,:,j),qdot,atm1%qc(:,:,j),j,5)
            call pcp(j , 2 , iym2 , kz)
            call cldfrac(j)
 
!           need also to set diffq to 0 here before calling diffut
            do k = 1 , kz
              do i = 1 , iym1
                diffq(i,k,j) = d_zero
              end do
            end do
 
!-----compute the diffusion terms:
!           the diffusion term for qv is stored in diffq. before
!           completing aten%qv computation, do not use diffq for other

!           purpose.
            call diffu_x(diffq(:,:,j),qvb3d,sps2%ps,xkc(:,:,j),j)
            call diffu_x(aten%qc(:,:,j),qcb3d,sps2%ps,xkc(:,:,j),j)
          end if
!
!chem2    compute the tracers tendencies
          if ( ichem == 1 ) then
            call zenitm(coszrs,iy,j)
            call tractend2(j,xkc)
          end if
!chem2_
!
#ifndef BAND
        end if           !end if(j /= jxm1) test
#endif
!----------------------------------------------------------------------
!*****compute the pbl fluxes:
!       the diffusion and pbl tendencies of t and qv are stored in
!       difft and diffq.
!
        do k = 1 , kz
          do i = 2 , iym1
            aten%u(i,k,j) = d_zero
            aten%v(i,k,j) = d_zero
          end do
        end do
 
 
!       ****** calculate solar zenith angle
        if ( (jyear == jyear0 .and. ktau == 0) .or. &
             mod(ktau+1,nbatst) == 0 .or. mod(ktau+1,ntrad) == 0 ) then
          call zenitm(coszrs,iy,j)
          call slice1D(j)
        end if
 
!       ****** calculate albedo
        if ( (jyear == jyear0 .and. ktau == 0) .or. &
             mod(ktau+1,ntrad) == 0 ) then
#ifdef CLM
          call albedoclm(j,iemiss)
#else
          call albedov(j,iemiss)
#endif
        end if
 
!       ****** call ccm3 radiative transfer package
        if ( (jyear == jyear0 .and. ktau == 0) .or. &
            mod(ktau+1,ntrad) == 0 ) call colmod3(j)
 
#ifndef CLM
!       ****** call vector bats for surface physics calculations
        if ( (jyear == jyear0 .and. ktau == 0) .or. &
           &  mod(ktau+1,nbatst) == 0 ) then
          dtbat = dto2*dble(nbatst)
          if ( jyear == jyear0 .and. ktau == 0 ) dtbat = dt
          call interf(1 , j , kz , 2 , iym1 , nnsg)
          call vecbats
!         Zeng ocean flux model
          if ( iocnflx == 2 ) call zengocndrv(j , nnsg , 2 , iym1 , kz)
!         Hostetler lake model for every BATS timestep at lake points
          if ( lakemod == 1 ) then
            call lakedrv(j)
          endif
!         ****** accumulate quantities for energy and moisture budgets
          call interf(2 , j , kz , 2 , iym1 , nnsg)
        end if
#endif
 
      end do

#ifdef CLM
      if ( ( jyear == jyear0 .and. ktau == 0 ) .or. &
          & mod(ktau+1,ntrad) == 0 ) then
        r2cdoalb = .true.
      else
        r2cdoalb = .false.
      end if
      if ( (jyear == jyear0 .and. ktau == 0 ) .or. &
         & mod(ktau+1,nbatst) == 0 ) then
        ! Timestep used is the same as for bats
        if ( jyear == jyear0 .and. ktau == 0 ) then
          r2cnstep = 0
        else
          r2cnstep = (ktau+1)/nbatst
        end if
        dtbat = dto2*nbatst
        ! CLM j loop is in mtrxclm
        call mtrxclm
      end if
#endif

      if ( icup == 1 ) then
        call htdiff(dto2,dxsq,akht1)
      end if
!     call medium resolution pbl
      if ( ibltyp == 1 ) call holtbl
!
#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif
!       add ccm radiative transfer package-calculated heating rates to
!       temperature tendency
        do k = 1 , kz
          do i = 2 , iym2
            ! heating rate in deg/sec
            aten%t(i,k,j) = aten%t(i,k,j) + sps2%ps(i,j)*heatrt(i,k,j)
          end do
        end do
!
#ifndef BAND
#ifdef MPP1
        if ( myid /= nproc-1 .or. j /= jendx ) then
#else
        if ( j /= jxm1 ) then
#endif
#endif
!
!..tq.add horizontal diffusion and pbl tendencies for t and qv to aten%t
!         and aten%qv for calculating condensational term in subroutine
!         "condtq".
!
          do k = 1 , kz
            do i = 2 , iym2
              aten%t(i,k,j) = aten%t(i,k,j) + difft(i,k,j)
            end do
          end do
!
          do k = 1 , kz
            do i = 2 , iym2
              aten%qv(i,k,j) = aten%qv(i,k,j) + diffq(i,k,j)
            end do
          end do
!
!..tq.compute the condensation and precipitation terms for explicit
!         moisture scheme:
!
          call condtq(j)
!
!..tq.subtract horizontal diffusion and pbl tendencies from aten%t and
!         aten%qv for appling the sponge boundary conditions on t and qv:
!
          if ( iboudy == 4 ) then
            do k = 1 , kz
              do i = 2 , iym2
                aten%t(i,k,j) = aten%t(i,k,j) - difft(i,k,j)
              end do
            end do
            call sponge_t(ispgx,wgtx,aten%t(1,1,j),j)
            do k = 1 , kz
              do i = 2 , iym2
                aten%t(i,k,j) = aten%t(i,k,j) + difft(i,k,j)
              end do
            end do
            do k = 1 , kz
              do i = 2 , iym2
                aten%qv(i,k,j) = aten%qv(i,k,j) - diffq(i,k,j)
              end do
            end do
            call spongeqv(ispgx,wgtx,aten%qv(1,1,j),j)
            do k = 1 , kz
              do i = 2 , iym2
                aten%qv(i,k,j) = aten%qv(i,k,j) + diffq(i,k,j)
              end do
            end do
          end if
!
!..tq.apply the nudging boundary conditions:
!
          if ( iboudy == 1 .or. iboudy == 5 ) then
            xtm1 = xtime - dtmin
            if ( dabs(xtime) < 0.00001D0 .and. ldatez > idate0 )  &
              xtm1 = -dtmin
            call nudge_t(ispgx,fnudge,gnudge,xtm1,aten%t(:,:,j),j,   &
                       & iboudy)
            call nudgeqv(ispgx,fnudge,gnudge,xtm1,aten%qv(:,:,j),j,  &
                       & iboudy)
          end if
!
!..tq.forecast t, qv, and qc at tau+1:
!
          do k = 1 , kz
            do i = 2 , iym2
              atmc%qv(i,k,j) = atm2%qv(i,k,j) + dt*aten%qv(i,k,j)
            end do
          end do
!
          do k = 1 , kz
            do i = 2 , iym2
              atmc%qc(i,k,j) = atm2%qc(i,k,j) + dt*aten%qc(i,k,j)
            end do
          end do
!
          do k = 1 , kz
            do i = 2 , iym2
              atmc%t(i,k,j) = atm2%t(i,k,j) + dt*aten%t(i,k,j)
            end do
          end do
!
!chem2
!         forecast tracer chi at at tau+1:
          if ( ichem == 1 ) then
!
            do itr = 1 , ntr
              do k = 1 , kz
                do i = 2 , iym2
                  chic(i,k,j,itr) = chib(i,k,j,itr)                     &
                                  & + dt*chiten(i,k,j,itr)
                end do
              end do
            end do
          end if
!chem2_
#ifndef BAND
        end if       !end if(j /= jxm1),else test
#endif
      end do
!
#ifndef BAND
#ifdef MPP1
      do j = 1 , jendx
        if ( myid == 0 .and. j == 1 ) then
#else
      do j = 1 , jxm1
        if ( j == 1 ) then
#endif
          if ( ipgf == 1 ) then
            do k = 1 , kz
              do i = 1 , iym1
                td(i,k,j) = atm1%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
                ttld(i,k,j) = td(i,k,j) - sps1%ps(i,j)                  &
                            & *t00pg*((a(k)*sps1%ps(i,j)+r8pt)/p00pg)   &
                            & **pgfaa1
              end do
            end do
          else if ( ipgf == 0 ) then
            do k = 1 , kz
              do i = 1 , iym1
                td(i,k,j) = atm1%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
              end do
            end do
          end if
!
#ifdef MPP1
        else if ( myid == nproc-1 .and. j == jendx ) then
#else
        else if ( j == jxm1 ) then
#endif
!
!-----set td and psd at j=jlx equal to ta and sps1%ps:
!
          if ( ipgf == 1 ) then
            do k = 1 , kz
              do i = 1 , iym1
                td(i,k,j) = atm1%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
                ttld(i,k,j) = td(i,k,j) - sps1%ps(i,j)                  &
                            & *t00pg*((a(k)*sps1%ps(i,j)+r8pt)/p00pg)   &
                            & **pgfaa1
              end do
            end do
          else if ( ipgf == 0 ) then
            do k = 1 , kz
              do i = 1 , iym1
                td(i,k,j) = atm1%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
              end do
            end do
          end if
!
!
!..t..compute weighted p*t (td) for use in ssi:
!
        end if
      end do
#endif

#ifdef MPP1
      do j = 1 , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 1 , jxm1
#endif
#endif
#ifndef BAND
#ifdef MPP1
        if((.not.(myid == 0 .and. j == 1)) .and. &
           (.not.(myid == nproc-1 .and. j == jendx)) ) then
#else
        if( .not.(j == 1 .or. j == jxm1) ) then
#endif
#endif
!
!..t..compute weighted p*t (td) for use in ssi:
!
        if ( ipgf == 1 ) then
!
          do k = 1 , kz
            do i = 2 , iym2
              tvc = atmc%t(i,k,j)*(d_one+ep1*(atmc%qv(i,k,j))/psc(i,j))
              tva = atm1%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
              tvb = atm2%t(i,k,j)*(d_one+ep1* &
                                   (atm2%qv(i,k,j))/sps2%ps(i,j))
              td(i,k,j) = alpha*(tvc+tvb) + beta*tva
              ttld(i,k,j) = td(i,k,j) - psd(i,j)                        &
                          & *t00pg*((a(k)*psd(i,j)+r8pt)/p00pg)**pgfaa1
            end do
          end do
          do k = 1 , kz
            td(1,k,j) = atm1%t(1,k,j)*(d_one+ep1*(atmx%qv(1,k,j)))
            ttld(1,k,j) = td(1,k,j) - sps1%ps(1,j)                      &
                      & *t00pg*((a(k)*sps1%ps(1,j)+r8pt)/p00pg)**pgfaa1
            td(iym1,k,j) = atm1%t(iym1,k,j)* &
                          (d_one+ep1*(atmx%qv(iym1,k,j)))
            ttld(iym1,k,j) = td(iym1,k,j) - sps1%ps(iym1,j)             &
                          & *t00pg*((a(k)*sps1%ps(iym1,j)+r8pt)/p00pg)  &
                          & **pgfaa1
          end do
!
        else if ( ipgf == 0 ) then
!
          do k = 1 , kz
            do i = 2 , iym2
              tvc = atmc%t(i,k,j)*(d_one+ep1*(atmc%qv(i,k,j))/psc(i,j))
              tva = atm1%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
              tvb = atm2%t(i,k,j)*(d_one+ep1* &
                   (atm2%qv(i,k,j))/sps2%ps(i,j))
              td(i,k,j) = alpha*(tvc+tvb) + beta*tva
            end do
          end do
          do k = 1 , kz
            td(1,k,j) = atm1%t(1,k,j)*(d_one+ep1*(atmx%qv(1,k,j)))
            td(iym1,k,j) = atm1%t(iym1,k,j)* &
                           (d_one+ep1*(atmx%qv(iym1,k,j)))
          end do
 
        end if
#ifndef BAND
        end if
#endif
      end do
!
!----------------------------------------------------------------------
!**uv*compute the u and v tendencies:
#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif
!
!..uv.compute the diffusion terms:
!       put diffusion and pbl tendencies of u and v in difuu and difuv.
!
        do k = 1 , kz
          do i = 2 , iym1
            difuu(i,k,j) = aten%u(i,k,j)
            difuv(i,k,j) = aten%v(i,k,j)
          end do
        end do
!
        call diffu_d(difuu(:,:,j),ubd3d,sps2%pdot,mddom%msfd, &
                     xkc(:,:,j),j,1)
        call diffu_d(difuv(:,:,j),vbd3d,sps2%pdot,mddom%msfd, &
                     xkc(:,:,j),j,1)
!
!..uv.compute the horizontal advection terms for u and v:
!
        do k = 1 , kz
          do i = 2 , iym1
            aten%u(i,k,j) = d_zero
            aten%v(i,k,j) = d_zero
          end do
        end do
!
        call hadv_d(aten%u(:,:,j),atmx%u,dx16,j,3)
        call hadv_d(aten%v(:,:,j),atmx%v,dx16,j,3)
!
!..uv.compute coriolis terms:
!
        do k = 1 , kz
          do i = 2 , iym1
            aten%u(i,k,j) = aten%u(i,k,j) + &
                         mddom%f(i,j)*atm1%v(i,k,j)/mddom%msfd(i,j)
            aten%v(i,k,j) = aten%v(i,k,j) - &
                         mddom%f(i,j)*atm1%u(i,k,j)/mddom%msfd(i,j)
          end do
        end do
      end do
!
#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif
        jm1 = j-1
#if defined(BAND) && (!defined(MPP1))
        if(jm1 == 0) jm1 = jx
#endif
!
!..uv.compute pressure gradient terms:
!
        if ( ipgf == 1 ) then
          do k = 1 , kz
            do i = 2 , iym1
              psasum = psd(i,j) + psd(i-1,j) + psd(i,jm1) + psd(i-1,jm1)
              sigpsa = psasum
              tv1 = atmx%t(i-1,k,jm1)*(d_one+ep1*(atmx%qv(i-1,k,jm1)))
              tv2 = atmx%t(i,k,jm1)*(d_one+ep1*(atmx%qv(i,k,jm1)))
              tv3 = atmx%t(i-1,k,j)*(d_one+ep1*(atmx%qv(i-1,k,j)))
              tv4 = atmx%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
              rtbar = tv1 + tv2 + tv3 + tv4 - d_four*t00pg*             &
                    & ((a(k)*(psasum/d_four)+r8pt)/p00pg)**pgfaa1
              rtbar = rgas*rtbar*sigpsa/16.0D0
              aten%u(i,k,j) = aten%u(i,k,j) - rtbar * &
                     (dlog((psd(i,j)+psd(i-1,j))/d_two*a(k)+r8pt) -     &
                      dlog((psd(i,jm1)+psd(i-1,jm1))/d_two*a(k)+r8pt))/ &
                      (dx*mddom%msfd(i,j))
              aten%v(i,k,j) = aten%v(i,k,j) - rtbar * &
                     (dlog((psd(i,j)+psd(i,jm1))/d_two*a(k)+r8pt) -     &
                      dlog((psd(i-1,jm1)+psd(i-1,j))/d_two*a(k)+r8pt))/ &
                      (dx*mddom%msfd(i,j))
            end do
          end do
        else if ( ipgf == 0 ) then
          do k = 1 , kz
            do i = 2 , iym1
              psasum = psd(i,j) + psd(i-1,j) + psd(i,jm1) + psd(i-1,jm1)
              sigpsa = psasum
              tv1 = atmx%t(i-1,k,jm1)*(d_one+ep1*(atmx%qv(i-1,k,jm1)))
              tv2 = atmx%t(i,k,jm1)*(d_one+ep1*(atmx%qv(i,k,jm1)))
              tv3 = atmx%t(i-1,k,j)*(d_one+ep1*(atmx%qv(i-1,k,j)))
              tv4 = atmx%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
              rtbar = rgas*(tv1+tv2+tv3+tv4)*sigpsa/16.0D0
              aten%u(i,k,j) = aten%u(i,k,j) - rtbar * &
                      (dlog((psd(i,j)+psd(i-1,j))/d_two*a(k)+r8pt) -    &
                       dlog((psd(i,jm1)+psd(i-1,jm1))/d_two*a(k)+r8pt))/&
                       (dx*mddom%msfd(i,j))
              aten%v(i,k,j) = aten%v(i,k,j) - rtbar *                   &
                      (dlog((psd(i,j)+psd(i,jm1))/d_two*a(k)+r8pt) -    &
                       dlog((psd(i-1,jm1)+psd(i-1,j))/d_two*a(k)+r8pt))/&
                       (dx*mddom%msfd(i,j))
            end do
          end do
        else   ! ipgf if block
        end if
      end do
!
#ifdef MPP1
      do j = 1 , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 1 , jxm1
#endif
#endif
!
!..uv.compute geopotential height at half-k levels, cross points:
!
        if ( ipgf == 1 ) then
 
          do i = 1 , iym1
            tv = (ttld(i,kz,j)/psd(i,j))/(d_one+atmx%qc(i,kz,j)/ &
                                         (d_one+atmx%qv(i,kz,j)))
            phi(i,kz,j) = mddom%ht(i,j)                                &
                        & + rgas*t00pg/pgfaa1*((psd(i,j)+r8pt)/p00pg)   &
                        & **pgfaa1
            phi(i,kz,j) = phi(i,kz,j) - rgas * &
                    tv*dlog((a(kz)+r8pt/psd(i,j))/(d_one+r8pt/psd(i,j)))
          end do
 
          do k = 1 , kzm1
            lev = kz - k
            do i = 1 , iym1
              tvavg = ((ttld(i,lev,j)*dsigma(lev)+ttld(i,lev+1,j)* &
                    & dsigma(lev+1))/(psd(i,j)*(dsigma(lev)+       &
                    & dsigma(lev+1))))/(d_one+atmx%qc(i,lev,j)/    &
                    & (d_one+atmx%qv(i,lev,j)))
              phi(i,lev,j) = phi(i,lev+1,j) - rgas *    &
                     tvavg*dlog((a(lev)+r8pt/psd(i,j))/ &
                               (a(lev+1)+r8pt/psd(i,j)))
            end do
          end do
 
        else if ( ipgf == 0 ) then
 
          do i = 1 , iym1
            tv = (td(i,kz,j)/psd(i,j))/(d_one+atmx%qc(i,kz,j)/  &
                 (d_one+atmx%qv(i,kz,j)))
            phi(i,kz,j) = mddom%ht(i,j) - rgas * &
                 tv*dlog((a(kz)+r8pt/psd(i,j))/(d_one+r8pt/psd(i,j)))
          end do
 
          do k = 1 , kzm1
            lev = kz - k
            do i = 1 , iym1
              tvavg = ((td(i,lev,j)*dsigma(lev)+td(i,lev+1,j)*   &
                    & dsigma(lev+1))/(psd(i,j)*(dsigma(lev)+     &
                    & dsigma(lev+1))))/(d_one+atmx%qc(i,lev,j)/  &
                    & (d_one+atmx%qv(i,lev,j)))
              phi(i,lev,j) = phi(i,lev+1,j) - rgas *    &
                     tvavg*dlog((a(lev)+r8pt/psd(i,j))  &
                           & /(a(lev+1)+r8pt/psd(i,j)))
            end do
          end do
 
        else   ! ipgf if block
        end if
      end do
#ifdef MPP1
      call mpi_sendrecv(phi(1,1,jxp),iy*kz,mpi_real8,ieast,1,           &
                      & phi(1,1,0),iy*kz,mpi_real8,iwest,1,             &
                      & mpi_comm_world,mpi_status_ignore,ierr)
#endif
#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif
        jm1 = j-1
#if defined(BAND) && (!defined(MPP1))
        if(jm1 == 0) jm1 = jx
#endif
!
!..uv.compute the geopotential gradient terms:
!
        do k = 1 , kz
          do i = 2 , iym1
            aten%u(i,k,j) = aten%u(i,k,j)                               &
                        & -(psd(i-1,jm1)+psd(i,jm1)+psd(i-1,j)+psd(i,j))&
                        & *(phi(i,k,j)+phi(i-1,k,j)-phi(i,k,jm1)-       &
                        & phi(i-1,k,jm1))/(dx8*mddom%msfd(i,j))
            aten%v(i,k,j) = aten%v(i,k,j)                               &
                        & -(psd(i-1,jm1)+psd(i,jm1)+psd(i-1,j)+psd(i,j))&
                        & *(phi(i,k,j)+phi(i,k,jm1)-phi(i-1,k,j)-       &
                        & phi(i-1,k,jm1))/(dx8*mddom%msfd(i,j))
          end do
        end do
      end do
!
#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif
!
!..uv.compute teh vertical advection terms:
!
        call vadv(aten%u(1,1,j),qdot,atm1%u(1,1,j),j,4)
        call vadv(aten%v(1,1,j),qdot,atm1%v(1,1,j),j,4)
!
!..uv.apply the sponge boundary condition on u and v:
!
        if ( iboudy == 4 ) then
          call sponge_u(ispgd,wgtd,aten%u(1,1,j),j)
          call sponge_v(ispgd,wgtd,aten%v(1,1,j),j)
        end if
!
!..uv.apply the nudging boundary conditions:
!
        if ( iboudy == 1 .or. iboudy == 5 ) then
          call nudge_u(ispgd,fnudge,gnudge,xtm1,aten%u(:,:,j),j,     &
                     & iboudy)
          call nudge_v(ispgd,fnudge,gnudge,xtm1,aten%v(:,:,j),j,     &
                     & iboudy)
        end if
!
!..uv.add the diffusion and pbl tendencies to aten%u and aten%v:
!
        do k = 1 , kz
          do i = 2 , iym1
            aten%u(i,k,j) = aten%u(i,k,j) + difuu(i,k,j)
            aten%v(i,k,j) = aten%v(i,k,j) + difuv(i,k,j)
          end do
        end do
!
!..uv.forecast p*u and p*v at tau+1:
!
        do k = 1 , kz
          do i = 2 , iym1
            atmc%u(i,k,j) = atm2%u(i,k,j) + dt*aten%u(i,k,j)
            atmc%v(i,k,j) = atm2%v(i,k,j) + dt*aten%v(i,k,j)
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
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif
        do k = 1 , kz
          do i = 2 , iym1
            atm2%u(i,k,j) = omuhf*atm1%u(i,k,j)/mddom%msfd(i,j)  &
                      & + gnuhf*(atm2%u(i,k,j)+atmc%u(i,k,j))
            atm2%v(i,k,j) = omuhf*atm1%v(i,k,j)/mddom%msfd(i,j)  &
                      & + gnuhf*(atm2%v(i,k,j)+atmc%v(i,k,j))
            atm1%u(i,k,j) = atmc%u(i,k,j)
            atm1%v(i,k,j) = atmc%v(i,k,j)
          end do
        end do
      end do
!
#ifdef MPP1
      do j = jbegin , jendm
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm2
#endif
#endif
        do k = 1 , kz
          do i = 2 , iym2
            atm2%t(i,k,j) = omuhf*atm1%t(i,k,j) + &
                            gnuhf*(atm2%t(i,k,j)+atmc%t(i,k,j))
            atm1%t(i,k,j) = atmc%t(i,k,j)
            qvas = atmc%qv(i,k,j)
            if ( qvas < dlowval ) qvas = d_zero
            qvbs = omuhf*atm1%qv(i,k,j) + &
                   gnuhf*(atm2%qv(i,k,j)+atmc%qv(i,k,j))
            if ( qvbs < dlowval ) qvbs = d_zero
            atm2%qv(i,k,j) = qvbs
            atm1%qv(i,k,j) = qvas
            qcas = atmc%qc(i,k,j)
            if ( qcas < dlowval ) qcas = d_zero
            qcbs = omu*atm1%qc(i,k,j) + &
                   gnu*(atm2%qc(i,k,j)+atmc%qc(i,k,j))
            if ( qcbs < dlowval ) qcbs = d_zero
            atm2%qc(i,k,j) = qcbs
            atm1%qc(i,k,j) = qcas
          end do
!chem2
          if ( ichem == 1 ) then
            do itr = 1 , ntr
              do i = 2 , iym2
                chias = chic(i,k,j,itr)
                if ( chias < dlowval ) chias = d_zero
                chibs = omu*chia(i,k,j,itr)                             &
                      & + gnu*(chib(i,k,j,itr)+chic(i,k,j,itr))
                if ( chibs < dlowval ) chibs = d_zero
                chib(i,k,j,itr) = chibs
                chia(i,k,j,itr) = chias
              end do
            end do
          end if
!chem2_
        end do
        do i = 2 , iym2
          sps2%ps(i,j) = omuhf*sps1%ps(i,j) + &
                         gnuhf*(sps2%ps(i,j)+psc(i,j))
          sps1%ps(i,j) = psc(i,j)
        end do
      end do
      if ( ehso4 ) then
        do k = 1 , kz
#ifdef MPP1
          do j = 1 , jendx
#else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
#endif
            do i = 1 , iym1
              aermm(i,k,j) = sulf%so4(i,k,j)
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
      ntime = ntime + idnint(dtmin*minph)
      if ( dabs(xtime-ibdyfrq*minph) < 0.00001D0 ) then
        call addhours(ldatez, ibdyfrq)
        call split_idate(ldatez, lyear, lmonth, lday, lhour)
        nnnnnn = nnnnnn + 1
        xtime = d_zero
        if ( lfirstjanatmidnight(ldatez) .and. xtime < 0.0001D0 ) then
          jyear = lyear
          ktau = 0
          ntime = 0
        end if
      end if
      if ( jyear /= jyear0 .or. ktau /= 0 ) dt = dt2
      dto2 = dt/d_two
!
!-----compute the amounts advected through the lateral boundaries:
!     *** note *** we must calculate the amounts advected through
!     the lateral boundaries before updating the values
!     at boundary slices.
!
#ifndef BAND
      if (debug_level > 2) then
        call conadv
        if ( ichem == 1 ) call tracdiag(xkc)
      end if
#endif
 
!-----fill up the boundary values for xxb and xxa variables:
!
      call bdyval(xtime,iexec)
!
!-----compute the nonconvective precipitation:
!
!chem2_
!     do cumulus transport of tracers
      if ( ichem == 1 .and. ichcumtra == 1 ) call cumtran
 
!chem2_
 
!-----trace the mass conservation of dry air and water substance:
!
#ifndef BAND
      if (debug_level > 2) call conmas
#endif
!
!
!---- budgets for tracers
      if ( ichem == 1 ) call tracbud
!
!-----print out noise parameter:
!
      if ( jyear /= jyear0 .or. ktau > 1 ) then
        ptnbar = ptntot/dble(iptn)
        pt2bar = pt2tot/dble(iptn)
        icons = 0
#ifdef MPP1
        icons_mpi = 0
        do j = jbegin , jendm
#else
#ifdef BAND
        do j = 1 , jx
#else
        do j = 2 , jxm2
#endif
#endif
          icons = icons + icon(j)
        end do
#ifdef MPP1
        icons_mpi = 0
        call mpi_allreduce(icons,icons_mpi,1,mpi_integer,mpi_sum,       &
                         & mpi_comm_world,ierr)
#endif
        xday = ((nnnnnn-nstrt0)*ibdyfrq*minph+xtime-dtmin)/minpd
        ! Added a check for nan...
        if ((ptnbar/=ptnbar) .or. &
           ((ptnbar>d_zero).eqv.(ptnbar<=d_zero))) then
#ifdef MPP1
          if ( myid == 0 ) then
#endif
          write (*,*) 'WHUUUUBBBASAAAGASDDWD!!!!!!!!!!!!!!!!'
          write (*,*) 'No more atmosphere here....'
          write (*,*) 'CFL violation detected, so model STOP'
          write (*,*) '#####################################'
          write (*,*) '#            DECREASE DT !!!!       #'
          write (*,*) '#####################################'
          call fatal(__FILE__,__LINE__,'CFL VIOLATION')
#ifdef MPP1
        end if
#endif
        end if
#ifdef MPP1
        if ( myid == 0 ) then
          if ( mod(ktau,50) == 0 ) print 99001 , xday , ktau , ptnbar , &
             & pt2bar , icons_mpi
        end if
#else
        if ( mod(ktau,50) == 0 ) print 99001 , xday , ktau , ptnbar ,   &
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
      if ( dabs(xtime) < 0.00001D0 .and. ldatez /= idate1 ) then
        call solar1(xtime)
        dectim = dnint(minpd+dectim)
#ifdef MPP1
        if ( myid == 0 ) write (*,*) ' dectim = ' , dectim
#else
        write (*,*) ' dectim = ' , dectim
#endif
      end if
!
      call time_end(subroutine_name,idindx)
      end subroutine tend
!
      end module mod_tendency
