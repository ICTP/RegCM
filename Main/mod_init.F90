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
 
      module mod_init
!
! RegCM Init module
!
      use mod_constants
      use mod_dynparam
      use mod_o3blk
      use mod_runparams
      use mod_bats
      use mod_lake
      use mod_vecbats
      use mod_pmoist
      use mod_main
      use mod_mainchem
      use mod_bdycod
      use mod_rad
      use mod_trachem
      use mod_message
      use mod_date
      use mod_radiation
      use mod_sun
      use mod_ncio
      use mod_savefile
      use mod_diagnosis
      use mod_cu_bm
#ifdef MPP1
      use mod_mppio
#ifdef CLM
      use mod_clm
      use clm_varsur , only : init_tgb , init_grid
#endif
#endif
!
      private
!
      public :: init
!
      real(8) , parameter :: tlp = 50.D0
      real(8) , parameter :: ts00 = 288.D0
!
      contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine reads in the initial and boundary conditions.   c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine init
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
      integer :: i , ibin , im1h , ip1h , ist ,itr , j , k , n , &
                 icbc_date
      real(8) :: hg1 , hg2 , hg3 , hg4 , hgmax
      integer :: jp1 , jm1
#ifdef MPP1
      real(8) , dimension(iy,jxp) :: psdot
      integer :: allrec , ierr , l
#else
      real(8) , dimension(iy,jx) :: psdot
#endif
      logical :: existing
!
      existing = .false.
#ifdef MPP1
#ifndef BAND
      peb  = 0.0
      pwb  = 0.0
      pebt = 0.0
      pwbt = 0.0
      teb  = 0.0
      twb  = 0.0
      tebt = 0.0
      twbt = 0.0
      qeb  = 0.0
      qwb  = 0.0
      qebt = 0.0
      qwbt = 0.0
      ueb  = 0.0
      uwb  = 0.0
      uebt = 0.0
      uwbt = 0.0
      veb  = 0.0
      vwb  = 0.0
      vebt = 0.0
      vwbt = 0.0
#endif
      pnb  = 0.0
      pss  = 0.0
      pnbt = 0.0
      psbt = 0.0
      tnb  = 0.0
      tsb  = 0.0
      tnbt = 0.0
      tsbt = 0.0
      qnb  = 0.0
      qsb  = 0.0
      qnbt = 0.0
      qsbt = 0.0
      unb  = 0.0
      usb  = 0.0
      unbt = 0.0
      usbt = 0.0
      vnb  = 0.0
      vsb  = 0.0
      vnbt = 0.0
      vsbt = 0.0
#endif
      tgmx_o = -1.E30
      t2mx_o = -1.E30
      tgmn_o =  1.E30
      t2mn_o =  1.E30
      w10x_o = -1.E30
      psmn_o =  1.E30

      ndate0 = idate1
      ndate1 = ndate0
      if (ndate0.eq.globidate1 .or.                                   &
         (((ndate0/10000)*100+1)*100 .eq.                             &
         ((globidate1/10000)*100+1)*100 ) ) then
        icbc_date = globidate1
      else
          icbc_date = ((ndate0/10000)*100+1)*100
      end if
#ifdef MPP1
      if ( myid.eq.0 ) then
        call open_icbc(icbc_date)
      end if
#else
      call open_icbc(icbc_date)
#endif
!
      if ( .not.ifrest ) then
!-----for initial run--not using restart
!
!------set rainwater and cloud water equal to zero initially.
!
        atm1%qc = 0.0
        atm2%qc = 0.0
!
!chem2
        if ( ichem.eq.1 ) then
!qhy      tchie, tchitb(replace tchidp:deposition)
!         initialize removal terms
          remlsc = 0.0
          remcvc = 0.0
          rxsg   = 0.0
          rxsaq1 = 0.0
          rxsaq2 = 0.0
          remdrd = 0.0
          wdlsc  = 0.0
        end if
!chem2_
!------set the variables related to blackadar pbl equal to 0 initially.
!
        if ( ibltyp.ne.0 ) then
          sfsta%hfx = 0.0
          sfsta%qfx = 0.0
        end if
!
        if ( icup.eq.1 ) then
          rsheat = 0.0
          rswat  = 0.0
        end if
!
!------read in the initial conditions for large domain:
!       the initial conditions are the output from PREPROC/ICBC.
!
#ifdef MPP1
#ifdef CLM
        if ( .not. allocated(init_tgb) ) allocate(init_tgb(iy,jx))
#endif
        if ( myid.eq.0 ) then
          call read_icbc(ndate0,ps0_io,ts0_io,ub0_io,vb0_io, &
                   &     tb0_io,qb0_io,so0_io)
          write (6,*) 'READY IC DATA for ', ndate0
          ps0_io = ps0_io/10.0
          do j = 1 , jx
            do k = 1 , kz
              do i = 1 , iy
                sav_0(i,k,j) = ub0_io(i,k,j)
                sav_0(i,kz+k,j) = vb0_io(i,k,j)
                sav_0(i,kz*2+k,j) = qb0_io(i,k,j)
                sav_0(i,kz*3+k,j) = tb0_io(i,k,j)
              end do
            end do
            do i = 1 , iy
              sav_0(i,kz*4+1,j) = ps0_io(i,j)
              sav_0(i,kz*4+2,j) = ts0_io(i,j)
            end do
          end do
          if ( ehso4 ) then
            do j = 1 , jx
              do k = 1 , kz
                do i = 1 , iy
                  sav_0s(i,k,j) = so0_io(i,k,j)
                end do
              end do
            end do
          end if
        end if
!
!       Start transmission of data to other processors
!
        call mpi_scatter(sav_0,iy*(kz*4+2)*jxp,mpi_real8,        &
                       & sav0, iy*(kz*4+2)*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        if ( ehso4 )                                                    &
          & call mpi_scatter(sav_0s,iy*kz*jxp,mpi_real8,         &
          &                  sav0s, iy*kz*jxp,mpi_real8,         &
          &                  0,mpi_comm_world,ierr)
        do j = 1 , jendl
          do k = 1 , kz
            do i = 1 , iy
              ub0(i,k,j) = sav0(i,k,j)
              vb0(i,k,j) = sav0(i,kz+k,j)
              qb0(i,k,j) = sav0(i,kz*2+k,j)
              tb0(i,k,j) = sav0(i,kz*3+k,j)
            end do
          end do
          do i = 1 , iy
            ps0(i,j) = sav0(i,kz*4+1,j)
            ts0(i,j) = sav0(i,kz*4+2,j)
          end do
          if ( ehso4 ) then
            do k = 1 , kz
              do i = 1 , iy
                so0(i,k,j) = sav0s(i,k,j)
              end do
            end do
          end if
        end do
!
!       Convert surface pressure to pstar
!
        do j = 1 , jendl
          do i = 1 , iy
            ps0(i,j) = ps0(i,j) - r8pt
          end do
        end do
!=======================================================================
!
!       this routine determines p(.) from p(x) by a 4-point
!       interpolation. on the x-grid, a p(x) point outside the grid
!       domain is assumed to satisfy p(0,j)=p(1,j); p(iy,j)=p(iym1,j);
!       and similarly for the i's.
!
        call mpi_sendrecv(ps0(1,jxp),iy,mpi_real8,ieast,1,              &
                        & ps0(1,0),iy,mpi_real8,iwest,1,                &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        do j = jbegin , jendx
          do i = 2 , iym1
            psdot(i,j) = 0.25*(ps0(i,j)   + ps0(i-1,j) +                &
                         &     ps0(i,j-1) + ps0(i-1,j-1))
          end do
        end do
!
#ifndef BAND
        do i = 2 , iym1
          if ( myid.eq.0 ) psdot(i,1) = 0.5*(ps0(i,1)+ps0(i-1,1))
          if ( myid.eq.nproc-1 ) psdot(i,jendl)                         &
             & = 0.5*(ps0(i,jendx)+ps0(i-1,jendx))
        end do
#endif
!
        do j = jbegin , jendx
          psdot(1,j) = 0.5*(ps0(1,j)+ps0(1,j-1))
          psdot(iy,j) = 0.5*(ps0(iym1,j)+ps0(iym1,j-1))
        end do
!
#ifndef BAND
        if ( myid.eq.0 ) then
          psdot(1,1) = ps0(1,1)
          psdot(iy,1) = ps0(iym1,1)
        end if
        if ( myid.eq.nproc-1 ) then
          psdot(1,jendl) = ps0(1,jendx)
          psdot(iy,jendl) = ps0(iym1,jendx)
        end if
#endif
!
!=======================================================================
!       Couple pressure u,v,t,q
        do k = 1 , kz
          do j = 1 , jendl
            do i = 1 , iy
              ub0(i,k,j) = ub0(i,k,j)*psdot(i,j)
              vb0(i,k,j) = vb0(i,k,j)*psdot(i,j)
              tb0(i,k,j) = tb0(i,k,j)*ps0(i,j)
              qb0(i,k,j) = qb0(i,k,j)*ps0(i,j)
            end do
         end do
        end do
!
!       Initialize variables and convert to double precision
!
        do k = 1 , kz
          do j = 1 , jendl
            do i = 1 , iy
              atm1%u(i,k,j) = ub0(i,k,j)
              atm2%u(i,k,j) = ub0(i,k,j)
              atm1%v(i,k,j) = vb0(i,k,j)
              atm2%v(i,k,j) = vb0(i,k,j)
              atm1%qv(i,k,j) = qb0(i,k,j)
              atm2%qv(i,k,j) = qb0(i,k,j)
              atm1%t(i,k,j) = tb0(i,k,j)
              atm2%t(i,k,j) = tb0(i,k,j)
            end do
          end do
        end do
        do j = 1 , jendl
          do i = 1 , iy
            sps1%ps(i,j) = ps0(i,j)
            sps2%ps(i,j) = ps0(i,j)
            sts1%tg(i,j) = ts0(i,j)
            sts2%tg(i,j) = ts0(i,j)
          end do
        end do
        if (iseaice == 1) then
          do j = 1 , jendx
            do i = 1 , iym1
              if ( mddom%satbrt(i,j).gt.13.5 .and. &
                   mddom%satbrt(i,j).lt.15.5 ) then
                if ( ts0(i,j).le.271.38 ) then
                  sts1%tg(i,j) = 271.38
                  sts2%tg(i,j) = 271.38
                  ts0(i,j) = 271.38
                  do n = 1, nnsg
                    ocld2d(n,i,j) = 2.0D0
                  end do
                else
                  do n = 1, nnsg
                    if ( satbrt1(n,i,j).gt.13.5 .and. &
                         satbrt1(n,i,j).lt.15.5 ) then
                      ocld2d(n,i,j) = 0.0D0
                    else
                      ocld2d(n,i,j) = 1.0D0
                    end if
                  end do
                end if
              else
                do n = 1, nnsg
                  if ( satbrt1(n,i,j).gt.13.5 .and. &
                       satbrt1(n,i,j).lt.15.5 ) then
                    ocld2d(n,i,j) = 0.0D0
                  else
                    ocld2d(n,i,j) = 1.0D0
                  end if
                end do
              end if
            end do
          end do
        else
          do j = 1 , jendx
            do i = 1 , iym1
              do n = 1, nnsg
                if ( satbrt1(n,i,j).gt.13.5 .and. &
                     satbrt1(n,i,j).lt.15.5 ) then
                  ocld2d(n,i,j) = 0.0D0
                else
                  ocld2d(n,i,j) = 1.0D0
                end if
              end do
            end do
          end do
        end if
        if (icup == 3) then
          do k = 1 , kz
            do j = 1 , jendl
              do i = 1 , iy
                tbase(i,k,j) = ts00 + &
                               tlp*dlog((sps1%ps(i,j)*a(k)+r8pt)/100.)
              end do
            end do
          end do
        end if
        if ( ehso4 ) then
          do k = 1 , kz
            do j = 1 , jendl
              do i = 1 , iy
                sulf%so4(i,k,j) = so0(i,k,j)
              end do
            end do
          end do
        end if
!
        do j = 1 , jendx
          do i = 1 , iym1
            sts1%tg(i,j) = atm1%t(i,kz,j)/sps1%ps(i,j)
            sts2%tg(i,j) = atm2%t(i,kz,j)/sps2%ps(i,j)
            sfsta%tgbb(i,j) = atm2%t(i,kz,j)/sps2%ps(i,j)
            sfsta%zpbl(i,j) = 500.  ! For Zeng Ocean Flux Scheme
          end do
        end do
        do j = 1 , jendx
          do i = 1 , iym1
            do k = 1 , nnsg
              snowc(k,i,j) = 0.
            end do
          end do
        end do
        if ( ichem.eq.1 ) then
          ssw2da    = 0.0
          sdeltk2d  = 0.0
          sdelqk2d  = 0.0
          sfracv2d  = 0.5
          sfracb2d  = 0.5
          sfracs2d  = 0.0
          svegfrac2d = 0.0
        end if
#else
        call read_icbc(ndate0,ps0,ts0,ub0,vb0,tb0,qb0,so0)
        write (6,*) 'READY IC DATA for ', ndate0
!
!       Convert surface pressure to pstar
!
        ps0 = ps0/10.0 - r8pt

!=======================================================================
!
!       this routine determines p(.) from p(x) by a 4-point
!       interpolation. on the x-grid, a p(x) point outside the grid
!       domain is assumed to satisfy p(0,j)=p(1,j);
!       p(iy,j)=p(iym1,j); and similarly for the i's.
!
#ifdef BAND
        do j = 1 , jx
#else
        do j = 2 , jxm1
#endif
          jm1 = j-1
#if defined(BAND) && (!defined(MPP1))
          if(jm1.eq.0) jm1=jx
#endif
          do i = 2 , iym1
            psdot(i,j) = 0.25*(ps0(i,j)+ps0(i-1,j)+                     &
                       &       ps0(i,jm1)+ps0(i-1,jm1))
          end do
        end do
!
#ifndef BAND
        do i = 2 , iym1
          psdot(i,1) = 0.5*(ps0(i,1)+ps0(i-1,1))
          psdot(i,jx) = 0.5*(ps0(i,jxm1)+ps0(i-1,jxm1))
        end do
#endif
!
#ifdef BAND
        do j = 1 , jx
#else
        do j = 2 , jxm1
#endif
          jm1 = j-1
#if defined(BAND) && (!defined(MPP1))
          if(jm1.eq.0) jm1=jx
#endif
          psdot(1,j) = 0.5*(ps0(1,j)+ps0(1,jm1))
          psdot(iy,j) = 0.5*(ps0(iym1,j)+ps0(iym1,jm1))
        end do
!
#ifndef BAND
        psdot(1,1) = ps0(1,1)
        psdot(iy,1) = ps0(iym1,1)
        psdot(1,jx) = ps0(1,jxm1)
        psdot(iy,jx) = ps0(iym1,jxm1)
#endif
!
!=======================================================================
!       Couple pressure u,v,t,q
!
        do k = 1 , kz
          do j = 1 , jx
            do i = 1 , iy
              ub0(i,k,j) = ub0(i,k,j)*psdot(i,j)
              vb0(i,k,j) = vb0(i,k,j)*psdot(i,j)
              tb0(i,k,j) = tb0(i,k,j)*ps0(i,j)
              qb0(i,k,j) = qb0(i,k,j)*ps0(i,j)
            end do
          end do
        end do
!
        mdate = ndate0
!
!       Initialize variables and convert to double precision
!
        do k = 1 , kz
          do j = 1 , jx
            do i = 1 , iy
              atm1%u(i,k,j) = ub0(i,k,j)
              atm2%u(i,k,j) = ub0(i,k,j)
              atm1%v(i,k,j) = vb0(i,k,j)
              atm2%v(i,k,j) = vb0(i,k,j)
              atm1%qv(i,k,j) = qb0(i,k,j)
              atm2%qv(i,k,j) = qb0(i,k,j)
              atm1%t(i,k,j) = tb0(i,k,j)
              atm2%t(i,k,j) = tb0(i,k,j)
            end do
          end do
        end do
        do j = 1 , jx
          do i = 1 , iy
            sps1%ps(i,j) = ps0(i,j)
            sps2%ps(i,j) = ps0(i,j)
            sts1%tg(i,j) = ts0(i,j)
            sts2%tg(i,j) = ts0(i,j)
          end do
        end do
        if (iseaice == 1) then
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do i = 1 , iym1
              if ( mddom%satbrt(i,j).gt.13.5 .and. &
                   mddom%satbrt(i,j).lt.15.5 ) then
                if ( ts0(i,j).le.271.38 ) then
                  sts1%tg(i,j) = 271.38
                  sts2%tg(i,j) = 271.38
                  ts0(i,j) = 271.38
                  do n = 1, nnsg
                    ocld2d(n,i,j) = 2.0D0
                  end do
                else
                  do n = 1, nnsg
                    if ( satbrt1(n,i,j).gt.13.5 .and. &
                         satbrt1(n,i,j).lt.15.5 ) then
                      ocld2d(n,i,j) = 0.0D0
                    else
                      ocld2d(n,i,j) = 1.0D0
                    end if
                  end do
                end if
              else
                do n = 1, nnsg
                  if ( satbrt1(n,i,j).gt.13.5 .and. &
                       satbrt1(n,i,j).lt.15.5 ) then
                    ocld2d(n,i,j) = 0.0D0
                  else
                    ocld2d(n,i,j) = 1.0D0
                  end if
                end do
              end if
            end do
          end do
        else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do i = 1 , iym1
              do n = 1, nnsg
                if ( satbrt1(n,i,j).gt.13.5 .and. &
                     satbrt1(n,i,j).lt.15.5 ) then
                  ocld2d(n,i,j) = 0.0D0
                else
                  ocld2d(n,i,j) = 1.0D0
                end if
              end do
            end do
          end do
        end if
        if (icup == 3) then
          do k = 1 , kz
            do j = 1 , jx
              do i = 1 , iy
                tbase(i,k,j) = ts00 + &
                           tlp*dlog((sps1%ps(i,j)*a(k)+r8pt)/100.)
              end do
            end do
          end do
        end if
        if ( ehso4 ) then
          do k = 1 , kz
            do j = 1 , jx
              do i = 1 , iy
                sulf%so4(i,k,j) = so0(i,k,j)
              end do
            end do
          end do
        end if
!
#ifdef BAND
        do j = 1 , jx
#else
        do j = 1 , jxm1
#endif
          do i = 1 , iym1
            sts1%tg(i,j) = atm1%t(i,kz,j)/sps1%ps(i,j)
            sts2%tg(i,j) = atm2%t(i,kz,j)/sps2%ps(i,j)
            sfsta%tgbb(i,j) = atm2%t(i,kz,j)/sps2%ps(i,j)
            sfsta%zpbl(i,j) = 500.
                       ! For Zeng Ocean Flux Scheme
          end do
        end do
#ifdef BAND
        do j = 1 , jx
#else
        do j = 1 , jxm1
#endif
          do i = 1 , iym1
            do k = 1 , nnsg
              snowc(k,i,j) = 0.
            end do
          end do
        end do
        if ( ichem.eq.1 ) then
          ssw2da = 0.0
          sdeltk2d = 0.0
          sdelqk2d = 0.0
          sfracv2d = 0.5
          sfracb2d = 0.5
          sfracs2d = 0.0
          svegfrac2d = 0.0
        end if
#endif
#ifndef BAND
        if (debug_level > 2) call initdiag
#endif
!
!chem2
        if ( ichem.eq.1 ) then
!-----set tracer concs to 1 (kg/kg) initially. Must convert this to p*
!-----mixing ratio to compute tendencies:
!US       mass test zero concs init input for advection
!qhy      initial chia is 10ppt
!hy       set the initial tracer concentration 10ppt (1.e-11), 9/4/98
 
          do itr = 1 , ntr
            do k = 1 , kz
#ifdef MPP1
              do j = 1 , jendx
                do i = 1 , iym1
                  chia(i,k,j,itr) = sps1%ps(i,j)*0.0D0
                  chib(i,k,j,itr) = sps2%ps(i,j)*0.0D0
!                 chia(i,k,j,itr)=sps1%ps(i,j)*1.e-11
!                 chib(i,k,j,itr)=sps2%ps(i,j)*1.e-11
                end do
              end do
#else
#ifdef BAND
              do j = 1 , jx
#else
              do j = 1 , jxm1
#endif
                do i = 1 , iym1
                  chia(i,k,j,itr) = sps1%ps(i,j)*0.0D0
                  chib(i,k,j,itr) = sps2%ps(i,j)*0.0D0
!                 chia(i,k,j,itr)=sps1%ps(i,j)*1.e-11
!                 chib(i,k,j,itr)=sps2%ps(i,j)*1.e-11
                end do
              end do
#endif
            end do
          end do
 
        end if
!chem2_
!
!------set rainc and rainnc equal to 0. initially
!
        sfsta%rainc  = 0.0
        sfsta%rainnc = 0.0
 
        if ( icup==4 .or. icup==99 .or. icup==98) then
          cbmf2d = 0.0
        end if
!
      else ! ifrest=.true.
!
!-----when ifrest=.true., read in the data saved from previous run
!       for large domain
!
#ifdef MPP1
        if ( myid.eq.0 ) then
          call read_savefile_part1(ndate0)
!
          print * , 'ozone profiles restart'
          do k = 1 , kzp1
            write (6,99004) o3prof_io(3,3,k)
          end do
          print 99005 , xtime , ktau , jyear
!
        end if
!
        if ( lakemod.eq.1 ) then
          call lakescatter
        endif

        if ( myid.eq.0 ) then
          do j = 1 , jx
            do k = 1 , kz
              do i = 1 , iy
                sav_0(i,k,j) = ub0_io(i,k,j)
                sav_0(i,kz+k,j) = vb0_io(i,k,j)
                sav_0(i,kz*2+k,j) = qb0_io(i,k,j)
                sav_0(i,kz*3+k,j) = tb0_io(i,k,j)
              end do
            end do
            do i = 1 , iy
              sav_0(i,kz*4+1,j) = ps0_io(i,j)
              sav_0(i,kz*4+2,j) = ts0_io(i,j)
            end do
            if ( ehso4 ) then
              do k = 1 , kz
                do i = 1 , iy
                  sav_0s(i,k,j) = so0_io(i,k,j)
                end do
              end do
            end if
          end do
        end if
        call mpi_scatter(sav_0,iy*(kz*4+2)*jxp,mpi_real8,        &
                       & sav0, iy*(kz*4+2)*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        if ( ehso4 )                                                    &
          &  call mpi_scatter(sav_0s,iy*kz*jxp,mpi_real8,        &
          &                   sav0s, iy*kz*jxp,mpi_real8,        &
          &                   0,mpi_comm_world,ierr)
        do j = 1 , jendl
          do k = 1 , kz
            do i = 1 , iy
              ub0(i,k,j) = sav0(i,k,j)
              vb0(i,k,j) = sav0(i,kz+k,j)
              qb0(i,k,j) = sav0(i,kz*2+k,j)
              tb0(i,k,j) = sav0(i,kz*3+k,j)
            end do
          end do
          do i = 1 , iy
            ps0(i,j) = sav0(i,kz*4+1,j)
            ts0(i,j) = sav0(i,kz*4+2,j)
          end do
          if ( ehso4 ) then
            do k = 1 , kz
              do i = 1 , iy
                so0(i,k,j) = sav0s(i,k,j)
              end do
            end do
          end if
        end do

        if ( myid.eq.0 ) then
          do j = 1 , jx
            do k = 1 , kz
              do i = 1 , iy
                sav_0(i,k,j) = atm1_io%u(i,k,j)
                sav_0(i,kz+k,j) = atm2_io%u(i,k,j)
                sav_0(i,kz*2+k,j) = atm1_io%v(i,k,j)
                sav_0(i,kz*3+k,j) = atm2_io%v(i,k,j)
              end do
            end do
            do i = 1 , iy
              sav_0(i,kz*4+1,j) = psa_io(i,j)
              sav_0(i,kz*4+2,j) = psb_io(i,j)
            end do
          end do
        end if
        call mpi_scatter(sav_0,iy*(kz*4+2)*jxp,mpi_real8,        &
                       & sav0, iy*(kz*4+2)*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendl
          do k = 1 , kz
            do i = 1 , iy
              atm1%u(i,k,j) = sav0(i,k,j)
              atm2%u(i,k,j) = sav0(i,kz+k,j)
              atm1%v(i,k,j) = sav0(i,kz*2+k,j)
              atm2%v(i,k,j) = sav0(i,kz*3+k,j)
            end do
          end do
          do i = 1 , iy
            sps1%ps(i,j) = sav0(i,kz*4+1,j)
            sps2%ps(i,j) = sav0(i,kz*4+2,j)
          end do
        end do
        if ( myid.eq.0 ) then
          do j = 1 , jx
            do k = 1 , kz
              do i = 1 , iy
                sav_0(i,k,j) = atm1_io%t(i,k,j)
                sav_0(i,kz+k,j) = atm2_io%t(i,k,j)
                sav_0(i,kz*2+k,j) = atm1_io%qv(i,k,j)
                sav_0(i,kz*3+k,j) = atm2_io%qv(i,k,j)
              end do
            end do
            do i = 1 , iy
              sav_0(i,kz*4+1,j) = tga_io(i,j)
              sav_0(i,kz*4+2,j) = tgb_io(i,j)
            end do
          end do
        end if
        call mpi_scatter(sav_0,iy*(kz*4+2)*jxp,mpi_real8,        &
                       & sav0, iy*(kz*4+2)*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendl
          do k = 1 , kz
            do i = 1 , iy
              atm1%t(i,k,j) = sav0(i,k,j)
              atm2%t(i,k,j) = sav0(i,kz+k,j)
              atm1%qv(i,k,j) = sav0(i,kz*2+k,j)
              atm2%qv(i,k,j) = sav0(i,kz*3+k,j)
            end do
          end do
          do i = 1 , iy
            sts1%tg(i,j) = sav0(i,kz*4+1,j)
            sts2%tg(i,j) = sav0(i,kz*4+2,j)
          end do
        end do
        if ( myid.eq.0 ) then
          do j = 1 , jx
            do k = 1 , kz
              do i = 1 , iy
                sav_0(i,k,j) = atm1_io%qc(i,k,j)
                sav_0(i,kz+k,j) = atm2_io%qc(i,k,j)
                sav_0(i,kz*2+k,j) = fcc_io(i,k,j)
              end do
            end do
            do i = 1 , iy
              sav_0(i,kz*4+1,j) = rainc_io(i,j)
              sav_0(i,kz*4+2,j) = rainnc_io(i,j)
            end do
          end do
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do k = 1 , kz
              do i = 1 , iym1
                sav_0(i,kz*3+k,j) = heatrt_io(i,k,j)
              end do
            end do
          end do
        end if
        call mpi_scatter(sav_0,iy*(kz*4+2)*jxp,mpi_real8,        &
                       & sav0, iy*(kz*4+2)*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendl
          do k = 1 , kz
            do i = 1 , iy
              atm1%qc(i,k,j) = sav0(i,k,j)
              atm2%qc(i,k,j) = sav0(i,kz+k,j)
              fcc(i,k,j) = sav0(i,kz*2+k,j)
            end do
          end do
          do i = 1 , iy
            sfsta%rainc(i,j) = sav0(i,kz*4+1,j)
            sfsta%rainnc(i,j) = sav0(i,kz*4+2,j)
          end do
        end do
        do j = 1 , jendx
          do k = 1 , kz
            do i = 1 , iym1
              heatrt(i,k,j) = sav0(i,kz*3+k,j)
            end do
          end do
        end do
        if ( myid.eq.0 ) then
          do j = 1 , jx
            do i = 1 , iy
              sav_0a(i,1,j) = hfx_io(i,j)
              sav_0a(i,2,j) = qfx_io(i,j)
              sav_0a(i,3,j) = uvdrag_io(i,j)
              sav_0a(i,4,j) = tgbb_io(i,j)
            end do
            do n = 1 , nnsg
              do i = 1 , iy
                sav_0a(i,4+n,j) = snowc_io(n,i,j)
              end do
            end do
          end do
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do k = 1 , kzp1
              do i = 1 , iym1
                sav_0a(i,nnsg+4+k,j) = o3prof_io(i,k,j)
              end do
            end do
          end do
        end if
        allrec = kzp1 + 5 + nnsg
        call mpi_scatter(sav_0a,iy*allrec*jxp,mpi_real8,         &
                       & sav0a, iy*allrec*jxp,mpi_real8,         &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendl
          do i = 1 , iy
            sfsta%hfx(i,j) = sav0a(i,1,j)
            sfsta%qfx(i,j) = sav0a(i,2,j)
            sfsta%uvdrag(i,j) = sav0a(i,3,j)
            sfsta%tgbb(i,j) = sav0a(i,4,j)
          end do
          do n = 1 , nnsg
            do i = 1 , iy
              snowc(n,i,j) = sav0a(i,4+n,j)
            end do
          end do
        end do
        do j = 1 , jendx
          do k = 1 , kzp1
            do i = 1 , iym1
              o3prof(i,k,j) = sav0a(i,nnsg+4+k,j)
            end do
          end do
        end do
        if ( iocnflx.eq.2 )                                        &
          & call mpi_scatter(zpbl_io,iy*jxp,mpi_real8,             &
          &                  sfsta%zpbl,   iy*jxp,mpi_real8,       &
          &                  0,mpi_comm_world,ierr)
        if ( icup.eq.1 ) then
          if ( myid.eq.0 ) then
            do j = 1 , jx
              do k = 1 , kz
                do i = 1 , iy
                  sav_0c(i,k,j) = rsheat_io(i,k,j)
                  sav_0c(i,kz+k,j) = rswat_io(i,k,j)
                end do
              end do
            end do
          end if
          call mpi_scatter(sav_0c,iy*kz*2*jxp,mpi_real8,         &
                         & sav0c, iy*kz*2*jxp,mpi_real8,         &
                         & 0,mpi_comm_world,ierr)
          do j = 1 , jendl
            do k = 1 , kz
              do i = 1 , iy
                rsheat(i,k,j) = sav0c(i,k,j)
                rswat(i,k,j) = sav0c(i,kz+k,j)
              end do
            end do
          end do
        end if
        if ( icup.eq.3 ) then
          if ( myid.eq.0 ) then
            do j = 1 , jx
              do k = 1 , kz
                do i = 1 , iy
                  sav_0b(i,k,j) = tbase_io(i,k,j)
                end do
              end do
              do i = 1 , iy
                sav_0b(i,kzp1,j) = cldefi_io(i,j)
              end do
            end do
          end if
          call mpi_scatter(sav_0b,iy*(kzp1)*jxp,mpi_real8,       &
                         & sav0b, iy*(kzp1)*jxp,mpi_real8,       &
                         & 0,mpi_comm_world,ierr)
          do j = 1 , jendl
            do k = 1 , kz
              do i = 1 , iy
                tbase(i,k,j) = sav0b(i,k,j)
              end do
            end do
            do i = 1 , iy
              cldefi(i,j) = sav0b(i,kzp1,j)
            end do
          end do
        end if
        if ( icup==4 .or. icup==99 .or. icup==98 ) then
          call mpi_scatter(cbmf2d_io,iy*jxp,mpi_real8,             &
                         & cbmf2d,   iy*jxp,mpi_real8,             &
                         & 0,mpi_comm_world,ierr)
        end if
        if ( myid.eq.0 ) then
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do l = 1 , 4
              do k = 1 , kz
                do i = 1 , iym1
                  sav_1(i,(l-1)*kz+k,j) = absnxt_io(i,k,l,j)
                end do
              end do
            end do
          end do
          allrec = kz*4
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do l = 1 , kzp1
              do k = 1 , kzp1
                do i = 1 , iym1
                  sav_1(i,allrec+(l-1)*(kzp1)+k,j) = abstot_io(i,k,l,j)
                end do
              end do
            end do
          end do
          allrec = allrec + (kzp1)*(kz+1)
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do k = 1 , kzp1
              do i = 1 , iym1
                sav_1(i,allrec+k,j) = emstot_io(i,k,j)
              end do
            end do
          end do
          allrec = allrec + kzp1
        end if
        allrec = kz*4 + (kzp1*kzp2)
        call mpi_scatter(sav_1,iym1*allrec*jxp,mpi_real8,        &
                       & sav1, iym1*allrec*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendx
          do l = 1 , 4
            do k = 1 , kz
              do i = 1 , iym1
                absnxt(i,k,l,j) = sav1(i,(l-1)*kz+k,j)
              end do
            end do
          end do
        end do
        allrec = kz*4
        do j = 1 , jendx
          do l = 1 , kzp1
            do k = 1 , kzp1
              do i = 1 , iym1
                abstot(i,k,l,j) = sav1(i,allrec+(l-1)*(kzp1)+k,j)
              end do
            end do
          end do
        end do
        allrec = allrec + (kzp1)*(kz+1)
        do j = 1 , jendx
          do k = 1 , kzp1
            do i = 1 , iym1
              emstot(i,k,j) = sav1(i,allrec+k,j)
            end do
          end do
        end do
        if ( myid.eq.0 ) then
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do n = 1 , nnsg
              do i = 1 , iym1
                sav_2(i,n,j) = taf2d_io(n,i,j)
                sav_2(i,nnsg+n,j) = tlef2d_io(n,i,j)
                sav_2(i,nnsg*2+n,j) = ssw2d_io(n,i,j)
                sav_2(i,nnsg*3+n,j) = srw2d_io(n,i,j)
                sav_2(i,nnsg*4+n,j) = col2d_io(n,i,j)
              end do
            end do
            do i = 1 , iym1
              sav_2(i,nnsg*5+1,j) = sol2d_io(i,j)
              sav_2(i,nnsg*5+2,j) = solvd2d_io(i,j)
              sav_2(i,nnsg*5+3,j) = solvs2d_io(i,j)
              sav_2(i,nnsg*5+4,j) = flw2d_io(i,j)
            end do
          end do
        end if
        allrec = nnsg*5 + 4
        call mpi_scatter(sav_2,iym1*allrec*jxp,mpi_real8,        &
                       & sav2, iym1*allrec*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendx
          do n = 1 , nnsg
            do i = 1 , iym1
              taf2d(n,i,j) = sav2(i,n,j)
              tlef2d(n,i,j) = sav2(i,nnsg+n,j)
              ssw2d(n,i,j) = sav2(i,nnsg*2+n,j)
              srw2d(n,i,j) = sav2(i,nnsg*3+n,j)
              col2d(n,i,j) = sav2(i,nnsg*4+n,j)
            end do
          end do
          do i = 1 , iym1
            sol2d(i,j) = sav2(i,nnsg*5+1,j)
            solvd2d(i,j) = sav2(i,nnsg*5+2,j)
            solvs2d(i,j) = sav2(i,nnsg*5+3,j)
            flw2d(i,j) = sav2(i,nnsg*5+4,j)
          end do
        end do
        if ( myid.eq.0 ) then
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do n = 1 , nnsg
              do i = 1 , iym1
                sav_2(i,n,j) = tgb2d_io(n,i,j)
                sav_2(i,nnsg+n,j) = swt2d_io(n,i,j)
                sav_2(i,nnsg*2+n,j) = scv2d_io(n,i,j)
                sav_2(i,nnsg*3+n,j) = gwet2d_io(n,i,j)
                sav_2(i,nnsg*4+n,j) = tg2d_io(n,i,j)
              end do
            end do
            do i = 1 , iym1
              sav_2(i,nnsg*5+1,j) = flwd2d_io(i,j)
              sav_2(i,nnsg*5+2,j) = fsw2d_io(i,j)
              sav_2(i,nnsg*5+3,j) = sabv2d_io(i,j)
              sav_2(i,nnsg*5+4,j) = sinc2d_io(i,j)
            end do
          end do
        end if
        allrec = nnsg*5 + 4
        call mpi_scatter(sav_2,iym1*allrec*jxp,mpi_real8,        &
                       & sav2, iym1*allrec*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendx
          do n = 1 , nnsg
            do i = 1 , iym1
              tgb2d(n,i,j) = sav2(i,n,j)
              swt2d(n,i,j) = sav2(i,nnsg+n,j)
              scv2d(n,i,j) = sav2(i,nnsg*2+n,j)
              gwet2d(n,i,j) = sav2(i,nnsg*3+n,j)
              tg2d(n,i,j) = sav2(i,nnsg*4+n,j)
            end do
          end do
          do i = 1 , iym1
            flwd2d(i,j) = sav2(i,nnsg*5+1,j)
            fsw2d(i,j) = sav2(i,nnsg*5+2,j)
            sabv2d(i,j) = sav2(i,nnsg*5+3,j)
            sinc2d(i,j) = sav2(i,nnsg*5+4,j)
          end do
        end do
        if ( myid.eq.0 ) then
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do n = 1 , nnsg
              do i = 1 , iym1
                sav_2(i,n,j) = veg2d1_io(n,i,j)
                sav_2(i,nnsg+n,j) = sag2d_io(n,i,j)
                sav_2(i,nnsg*2+n,j) = sice2d_io(n,i,j)
                sav_2(i,nnsg*3+n,j) = dew2d_io(n,i,j)
                sav_2(i,nnsg*4+n,j) = ocld2d_io(n,i,j)
              end do
            end do
            do i = 1 , iym1
              sav_2(i,nnsg*5+1,j) = pptnc_io(i,j)
              sav_2(i,nnsg*5+2,j) = pptc_io(i,j)
              sav_2(i,nnsg*5+3,j) = prca2d_io(i,j)
              sav_2(i,nnsg*5+4,j) = prnca2d_io(i,j)
            end do
          end do
        end if
        allrec = nnsg*5 + 4
        call mpi_scatter(sav_2,iym1*allrec*jxp,mpi_real8,        &
                       & sav2, iym1*allrec*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendx
          do n = 1 , nnsg
            do i = 1 , iym1
              veg2d1(n,i,j) = sav2(i,n,j)
              sag2d(n,i,j) = sav2(i,nnsg+n,j)
              sice2d(n,i,j) = sav2(i,nnsg*2+n,j)
              dew2d(n,i,j) = sav2(i,nnsg*3+n,j)
              ocld2d(n,i,j) = sav2(i,nnsg*4+n,j)
            end do
          end do
          do i = 1 , iym1
            pptnc(i,j) = sav2(i,nnsg*5+1,j)
            pptc(i,j) = sav2(i,nnsg*5+2,j)
            prca2d(i,j) = sav2(i,nnsg*5+3,j)
            prnca2d(i,j) = sav2(i,nnsg*5+4,j)
          end do
        end do
        if ( myid.eq.0 ) then
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do n = 1 , nnsg
              do i = 1 , iym1
                sav_2a(i,n,j)        = ircp2d_io(n,i,j)
                sav_2a(i,nnsg+n,j)   = text2d_io(n,i,j)
              end do
            end do
            do i = 1 , iym1
              sav_2a(i,nnsg*2+1,j) = veg2d_io(i,j)
            end do
          end do
        end if
        allrec = nnsg*2 + 1
        call mpi_scatter(sav_2a,iym1*allrec*jxp,mpi_real8,       &
                       & sav2a, iym1*allrec*jxp,mpi_real8,       &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendx
          do n = 1 , nnsg
            do i = 1 , iym1
              ircp2d(n,i,j) = sav2a(i,n,j)
              text2d(n,i,j) = sav2a(i,nnsg+n,j)
            end do
          end do
          do i = 1 , iym1
            veg2d(i,j) = sav2a(i,nnsg*2+1,j)
          end do
        end do
        if ( ichem.eq.1 ) then
          if ( myid.eq.0 ) then
            do j = 1 , jx
              do n = 1 , ntr
                do k = 1 , kz
                  do i = 1 , iy
                    sav_4(i,(n-1)*kz+k,j) = chia_io(i,k,j,n)
                    sav_4(i,ntr*kz+(n-1)*kz+k,j) = chib_io(i,k,j,n)
                    sav_4(i,ntr*kz*2+(n-1)*kz+k,j) = remlsc_io(i,k,j,n)
                    sav_4(i,ntr*kz*3+(n-1)*kz+k,j) = remcvc_io(i,k,j,n)
                  end do
                end do
              end do
            end do
            allrec = 4*ntr*kz
            do j = 1 , jx
              do n = 1 , ntr
                do i = 1 , iy
                  sav_4(i,allrec+n,j) = remdrd_io(i,j,n)
                end do
              end do
            end do
            allrec = allrec + ntr
          end if
          allrec = ntr*(kz*4+1)
          call mpi_scatter(sav_4,iy*allrec*jxp,mpi_real8,        &
                         & sav4, iy*allrec*jxp,mpi_real8,        &
                         & 0,mpi_comm_world,ierr)
          do j = 1 , jendl
            do n = 1 , ntr
              do k = 1 , kz
                do i = 1 , iy
                  chia(i,k,j,n) = sav4(i,(n-1)*kz+k,j)
                  chib(i,k,j,n) = sav4(i,ntr*kz+(n-1)*kz+k,j)
                  remlsc(i,k,j,n) = sav4(i,ntr*kz*2+(n-1)*kz+k,j)
                  remcvc(i,k,j,n) = sav4(i,ntr*kz*3+(n-1)*kz+k,j)
                end do
              end do
            end do
          end do
          allrec = 4*ntr*kz
          do j = 1 , jendl
            do n = 1 , ntr
              do i = 1 , iy
                remdrd(i,j,n) = sav4(i,allrec+n,j)
              end do
            end do
          end do
        end if
        if ( myid.eq.0 ) then
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do i = 1 , iym1
              sav_4a(i,1,j) = ssw2da_io(i,j)
              sav_4a(i,2,j) = sdeltk2d_io(i,j)
              sav_4a(i,3,j) = sdelqk2d_io(i,j)
              sav_4a(i,4,j) = sfracv2d_io(i,j)
              sav_4a(i,5,j) = sfracb2d_io(i,j)
              sav_4a(i,6,j) = sfracs2d_io(i,j)
              sav_4a(i,7,j) = svegfrac2d_io(i,j)
            end do
          end do
        end if
        call mpi_scatter(sav_4a,iym1*7*jxp,mpi_real8,                 &
                       & sav4a, iym1*7*jxp,mpi_real8,                 &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendx
          do i = 1 , iym1
            ssw2da(i,j) = sav4a(i,1,j)
            sdeltk2d(i,j) = sav4a(i,2,j)
            sdelqk2d(i,j) = sav4a(i,3,j)
            sfracv2d(i,j) = sav4a(i,4,j)
            sfracb2d(i,j) = sav4a(i,5,j)
            sfracs2d(i,j) = sav4a(i,6,j)
            svegfrac2d(i,j) = sav4a(i,7,j)
          end do
        end do
#ifdef CLM
        if ( myid.eq.0 ) then
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do i = 1 , iym1
              sav_clmout(i,1,j)  = sols2d_io(i,j)
              sav_clmout(i,2,j)  = soll2d_io(i,j)
              sav_clmout(i,3,j)  = solsd2d_io(i,j)
              sav_clmout(i,4,j)  = solld2d_io(i,j)
              sav_clmout(i,5,j)  = aldirs2d_io(i,j)
              sav_clmout(i,6,j)  = aldirl2d_io(i,j)
              sav_clmout(i,7,j)  = aldifs2d_io(i,j)
              sav_clmout(i,8,j)  = aldifl2d_io(i,j)
              sav_clmout(i,9,j)  = coszrs2d_io(i,j)
            end do
          end do
        end if
        call mpi_scatter(sav_clmout,iym1*9*jxp,mpi_real8,             &
                       & sav_clmin, iym1*9*jxp,mpi_real8,             &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendx
          do i = 1 , iym1
            sols2d(i,j)   = sav_clmin(i,1,j)
            soll2d(i,j)   = sav_clmin(i,2,j)
            solsd2d(i,j)  = sav_clmin(i,3,j)
            solld2d(i,j)  = sav_clmin(i,4,j)
            aldirs2d(i,j) = sav_clmin(i,5,j)
            aldirl2d(i,j) = sav_clmin(i,6,j)
            aldifs2d(i,j) = sav_clmin(i,7,j)
            aldifl2d(i,j) = sav_clmin(i,8,j)
            coszrs2d(i,j) = sav_clmin(i,9,j)
          end do
        end do
#endif
        call mpi_bcast(mdate0,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(jyear0,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(ktau,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(jyear,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(xtime,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(ldatez,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(lyear,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(lmonth,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(lday,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(lhour,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(ntime,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(jyearr,1,mpi_integer,0,mpi_comm_world,ierr)

#ifndef BAND
        if (debug_level > 2) call mpidiag
#endif
        dt = dt2 ! First timestep successfully read in
#else
        call read_savefile_part1(ndate0)
!
        print * , 'ozone profiles restart'
        do k = 1 , kzp1
          write (6,99004) o3prof(3,3,k)
        end do
        print 99005 , xtime , ktau , jyear
!
        dt = dt2 ! First timestep successfully read in
#endif
!
!-----end of initial/restart if test
!
      end if
!
!     Move from param.F to fix the reatart problem found by
!     Zhang DongFeng
!
      if ( ipptls.eq.1 ) then
#ifdef MPP1
        do j = 1 , jendx
          do i = 1 , iym1
            if ( mddom%satbrt(i,j).gt.13.9 .and. &
                 mddom%satbrt(i,j).lt.15.1 ) then
              qck1(i,j) = qck1oce  ! OCEAN
              cgul(i,j) = guloce   ! OCEAN
              rh0(i,j) = rh0oce    ! OCEAN
            else
              qck1(i,j) = qck1land ! LAND
              cgul(i,j) = gulland  ! LAND
              rh0(i,j) = rh0land   ! LAND
            end if
          end do
        end do
#else
#ifdef BAND
        do j = 1 , jx
#else
        do j = 1 , jxm1
#endif
          do i = 1 , iym1
            if ( mddom%satbrt(i,j).gt.13.9 .and.  &
                 mddom%satbrt(i,j).lt.15.1 ) then
              qck1(i,j) = qck1oce  ! OCEAN
              cgul(i,j) = guloce   ! OCEAN
              rh0(i,j) = rh0oce    ! OCEAN
            else
              qck1(i,j) = qck1land ! LAND
              cgul(i,j) = gulland  ! LAND
              rh0(i,j) = rh0land   ! LAND
            end if
          end do
        end do
#endif
      end if
!chem2
      if ( ichem.eq.1 ) then
        iso2 = 0
        iso4 = 0
        ibchl = 0
        ibchb = 0
        iochl = 0
        iochb = 0
        ibin = 0
        do itr = 1 , ntr
          if ( chtrname(itr).eq.'SO2' ) iso2 = itr
          if ( chtrname(itr).eq.'SO4' ) iso4 = itr
          if ( chtrname(itr).eq.'BC_HL' ) ibchl = itr
          if ( chtrname(itr).eq.'BC_HB' ) ibchb = itr
          if ( chtrname(itr).eq.'OC_HL' ) iochl = itr
          if ( chtrname(itr).eq.'OC_HB' ) iochb = itr
          if ( chtrname(itr).eq.'DUST' ) then
            ibin = ibin + 1
            idust(ibin) = itr
          end if
        end do
      end if
!chem2_
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ****** initialize and define constants for vector bats
 
      if ( jyear.eq.jyear0 .and. ktau.eq.0 ) call initb

      if ( iemiss.eq.1 ) then
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
            do n = 1 , nnsg
              ist = nint(veg2d1(n,i,j))
              if ( ist.eq.0 ) then
                emiss2d(n,i,j) = 0.955D0
              else if ( ist.eq.8 ) then
                emiss2d(n,i,j) = 0.76D0
              else if ( ist.eq.11 ) then
                emiss2d(n,i,j) = 0.85D0
              else if ( ist.eq.12 ) then
                emiss2d(n,i,j) = 0.97D0
              else
                emiss2d(n,i,j) = 0.99D0 - (albvgs(ist)+albvgl(ist))     &
                               & *0.1D0
              end if
!             emiss2d(n,i,j) = 1.0d0
            end do
          end do
        end do
      end if
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!-----read in the boundary conditions for large domain:
!
!-----compute the solar declination angle:
!
      call solar1(xtime)
#ifdef CLM
      if ( .not. ifrest ) then
        init_grid = .true.
      else
        init_grid = .false.
      end if
#endif
      call inirad
!
!-----calculating topographical correction to diffusion coefficient
#ifdef MPP1
      do j = 1 , jendl
#else
      do j = 1 , jx
#endif
        do i = 1 , iy
          domfc%hgfact(i,j) = 1.
        end do
      end do
#ifdef BAND
#ifdef MPP1
      do j = jbegin , jendm
#else
      do j = 1 , jx
#endif
        jm1 = j-1
        jp1 = j+1
#if defined(BAND) && (!defined(MPP1))
        if(jm1.eq.0) jm1 = jx
        if(jp1.eq.jx+1) jp1 = 1
#endif
#else
#ifdef MPP1
      do j = jbegin , jendm
        if ( myid.eq.0 ) then
          jm1 = max0(j-1,2)
        else
          jm1 = j - 1
        end if
        if ( myid.eq.nproc-1 ) then
          jp1 = min0(j+1,jxp-2)
        else
          jp1 = j + 1
        end if
#else
      do j = 2 , jxm2
        jm1 = max0(j-1,2)
        jp1 = min0(j+1,jxm2)
#endif
#endif
        do i = 2 , iym2
          im1h = max0(i-1,2)
          ip1h = min0(i+1,iym2)
          hg1 = dabs((mddom%ht(i,j)-mddom%ht(im1h,j))/dx)
          hg2 = dabs((mddom%ht(i,j)-mddom%ht(ip1h,j))/dx)
          hg3 = dabs((mddom%ht(i,j)-mddom%ht(i,jm1))/dx)
          hg4 = dabs((mddom%ht(i,j)-mddom%ht(i,jp1))/dx)
          hgmax = dmax1(hg1,hg2,hg3,hg4)*rgti
          domfc%hgfact(i,j) = 1./(1.+(hgmax/0.001)**2.)
        end do
      end do
!
!-----set up output time:
!
      dectim = anint(xtime+dectim)
      write (aline, *) 'dectim = ' , dectim
      call say

99004 format (1x,7E12.4)
99005 format (' ***** restart file for large domain at time = ',f8.0,   &
             &' minutes, ktau = ',i7,' in year = ',i4,' read in')
!
      end subroutine init
!
!     compute ozone mixing ratio distribution
!
      subroutine inirad
 
      implicit none
!
      integer :: i , j , k
!
      if ( jyear.eq.jyear0 .and. ktau.eq.0 ) then
        do k = 1 , kz
#ifdef MPP1
          do j = 1 , jendl
#else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
#endif
            do i = 1 , iym1
              heatrt(i,k,j) = 0.
              o3prof(i,k,j) = 0.
            end do
          end do
        end do
#ifdef MPP1
        do j = 1 , jendl
#else
#ifdef BAND
        do j = 1 , jx
#else
        do j = 1 , jxm1
#endif
#endif
          do i = 1 , iym1
            o3prof(i,kzp1,j) = 0.
          end do
        end do
        call o3data
#ifdef MPP1
        if ( myid.eq.0 ) then
#endif
          write (6,*) 'ozone profiles'
          do k = 1 , kzp1
            write (6,99001) o3prof(3,k,2)
          end do
#ifdef MPP1
        end if
#endif
      end if
99001 format (1x,7E12.4)
 
      end subroutine inirad
!
      end module mod_init
