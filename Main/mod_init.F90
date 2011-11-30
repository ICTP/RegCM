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
  use mod_runparams
  use mod_mppparam
  use mod_lm_interface
  use mod_atm_interface
  use mod_che_interface
  use mod_cu_interface
  use mod_rad_interface
  use mod_pbl_interface
  use rrtmg_sw_init
  use rrtmg_lw_init
  use mod_pbl_interface
  use mod_precip
  use mod_bdycod
  use mod_mpmessage
  use mod_sun
  use mod_ncio
  use mod_savefile
  use mod_diagnosis
  use mod_mppio
  use mod_constants
#ifdef CLM
  use mod_clm
  use mod_lm_interface
  use clm_varsur , only : init_tgb , init_grid , numdays
#endif
!
  private
!
  public :: init
!
  real(8) , parameter :: tlp = 50.0D0
  real(8) , parameter :: ts00 = 288.0D0
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
#ifndef IBM
  use mpi
#else
  include 'mpif.h'
#endif
  implicit none
!
  integer :: i , ibin , im1h , ip1h , ist ,itr , j , k , n
  type (rcm_time_and_date) :: icbc_date
  real(8) :: hg1 , hg2 , hg3 , hg4 , hgmax
  integer :: jp1 , jm1
  real(8) , dimension(iy,jxp) :: psdot
  character(len=32) :: appdat
  integer :: mmrec , allrec , ierr , l

  tgmx_o = -1.E30
  t2mx_o = -1.E30
  tgmn_o =  1.E30
  t2mn_o =  1.E30
  w10x_o = -1.E30
  psmn_o =  1.E30

  bdydate1 = idate1
  bdydate2 = idate1
  if ( myid == 0 ) then
    if ( bdydate1 == globidate1 ) then
      icbc_date = bdydate1
    else
      icbc_date = monfirst(bdydate1)
    end if
    call open_icbc(icbc_date)
  end if
!
  if ( .not.ifrest ) then
!-----for initial run--not using restart
!
!------set rainwater and cloud water equal to zero initially.
!
    atm1%qc = d_zero
    atm2%qc = d_zero
!
!chem2
    if ( ichem == 1 ) then
!qhy      tchie, tchitb(replace tchidp:deposition)
!         initialize removal terms
      remlsc = d_zero
      remcvc = d_zero
      rxsg   = d_zero
      rxsaq1 = d_zero
      rxsaq2 = d_zero
      remdrd = d_zero
      wdlsc  = d_zero
    end if
!chem2_
!------set the variables related to blackadar pbl equal to 0 initially.
!
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      atm1%tke = tkemin
      atm2%tke = tkemin
    end if
!
    if ( icup == 1 ) then
      rsheat = d_zero
      rswat  = d_zero
    end if
!
!------read in the initial conditions for large domain:
!       the initial conditions are the output from PREPROC/ICBC.
!
#ifdef CLM
    if ( .not. allocated(init_tgb) ) allocate(init_tgb(iy,jx))
#endif
    if ( myid == 0 ) then
      mmrec = icbc_search(bdydate1)
      if (mmrec < 0) then
        appdat = tochar(bdydate2)
        call fatal(__FILE__,__LINE__, &
                   'ICBC for '//appdat//' not found')
      end if
      call read_icbc(ps0_io,ts0_io,ub0_io,vb0_io,tb0_io,qb0_io)
      appdat = tochar(bdydate1)
      write (6,*) 'READY IC DATA for ', appdat
      ps0_io = ps0_io*d_r10
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
    end if
!
!       Start transmission of data to other processors
!
    call mpi_scatter(sav_0,iy*(kz*4+2)*jxp,mpi_real8,        &
                     sav0, iy*(kz*4+2)*jxp,mpi_real8,        &
                     0,mycomm,ierr)
    do j = 1 , jendl
      do k = 1 , kz
        do i = 1 , iy
          xub%b0(i,k,j) = sav0(i,k,j)
          xvb%b0(i,k,j) = sav0(i,kz+k,j)
          xqb%b0(i,k,j) = sav0(i,kz*2+k,j)
          xtb%b0(i,k,j) = sav0(i,kz*3+k,j)
        end do
      end do
      do i = 1 , iy
        xpsb%b0(i,j) = sav0(i,kz*4+1,j)
        ts0(i,j) = sav0(i,kz*4+2,j)
      end do
    end do
!
!       Convert surface pressure to pstar
!
    do j = 1 , jendl
      do i = 1 , iy
        xpsb%b0(i,j) = xpsb%b0(i,j) - ptop
      end do
    end do
!=======================================================================
!
!       this routine determines p(.) from p(x) by a 4-point
!       interpolation. on the x-grid, a p(x) point outside the grid
!       domain is assumed to satisfy p(0,j)=p(1,j); p(iy,j)=p(iym1,j);
!       and similarly for the i's.
!
    call mpi_sendrecv(xpsb%b0(:,jxp),iy,mpi_real8,ieast,1,              &
                      xpsb%b0(:,0),iy,mpi_real8,iwest,1,                &
                      mycomm,mpi_status_ignore,ierr)
    do j = jbegin , jendx
      do i = 2 , iym1
        psdot(i,j) = (xpsb%b0(i,j)   + xpsb%b0(i-1,j) +   &
                      xpsb%b0(i,j-1) + xpsb%b0(i-1,j-1))*d_rfour
      end do
    end do
!
#ifndef BAND
    do i = 2 , iym1
      if ( myid == 0 ) psdot(i,1) = (xpsb%b0(i,1)+xpsb%b0(i-1,1))*d_half
      if ( myid == nproc-1 ) then
        psdot(i,jendl) = (xpsb%b0(i,jendx)+xpsb%b0(i-1,jendx))*d_half
      end if
    end do
#endif
!
    do j = jbegin , jendx
      psdot(1,j) = (xpsb%b0(1,j)+xpsb%b0(1,j-1))*d_half
      psdot(iy,j) = (xpsb%b0(iym1,j)+xpsb%b0(iym1,j-1))*d_half
    end do
!
#ifndef BAND
    if ( myid == 0 ) then
      psdot(1,1) = xpsb%b0(1,1)
      psdot(iy,1) = xpsb%b0(iym1,1)
    end if
    if ( myid == nproc-1 ) then
      psdot(1,jendl) = xpsb%b0(1,jendx)
      psdot(iy,jendl) = xpsb%b0(iym1,jendx)
    end if
#endif
!
!=======================================================================
!       Couple pressure u,v,t,q
    do k = 1 , kz
      do j = 1 , jendl
        do i = 1 , iy
          xub%b0(i,k,j) = xub%b0(i,k,j)*psdot(i,j)
          xvb%b0(i,k,j) = xvb%b0(i,k,j)*psdot(i,j)
          xqb%b0(i,k,j) = xqb%b0(i,k,j)*xpsb%b0(i,j)
          xtb%b0(i,k,j) = xtb%b0(i,k,j)*xpsb%b0(i,j)
        end do
     end do
    end do
!
!       Initialize variables and convert to double precision
!
    do k = 1 , kz
      do j = 1 , jendl
        do i = 1 , iy
          atm1%u(i,k,j) = xub%b0(i,k,j)
          atm2%u(i,k,j) = xub%b0(i,k,j)
          atm1%v(i,k,j) = xvb%b0(i,k,j)
          atm2%v(i,k,j) = xvb%b0(i,k,j)
          atm1%qv(i,k,j) = xqb%b0(i,k,j)
          atm2%qv(i,k,j) = xqb%b0(i,k,j)
          atm1%t(i,k,j) = xtb%b0(i,k,j)
          atm2%t(i,k,j) = xtb%b0(i,k,j)
        end do
      end do
    end do
    do j = 1 , jendl
      do i = 1 , iy
        sps1%ps(j,i) = xpsb%b0(i,j)
        sps2%ps(j,i) = xpsb%b0(i,j)
        sts1%tg(j,i) = ts0(i,j)
        sts2%tg(j,i) = ts0(i,j)
      end do
    end do
    if ( iseaice == 1 ) then
      do j = 1 , jendx
        do i = 1 , iym1
          if ( isocean(mddom%lndcat(j,i)) ) then
            if ( ts0(i,j) <= icetemp ) then
              sts1%tg(j,i) = icetemp
              sts2%tg(j,i) = icetemp
              ts0(i,j) = icetemp
              ldmsk(i,j) = 2
              do n = 1, nnsg
                ocld2d(n,j,i) = 2
              end do
            end if
          end if
        end do
      end do
    end if
#ifndef CLM
    if ( lakemod == 1 ) then
      do j = 1 , jendx
        do i = 1 , iym1
          if ( islake(mddom%lndcat(j,i)) ) then
            if ( ts0(i,j) <= icetemp ) then
              sts1%tg(j,i) = icetemp
              sts2%tg(j,i) = icetemp
              ts0(i,j) = icetemp
              ldmsk(i,j) = 2
              do n = 1, nnsg
                ocld2d(n,j,i) = 2
              end do
            end if
          end if
        end do
      end do
    end if
#endif
    if (icup == 3) then
      do k = 1 , kz
        do j = 1 , jendl
          do i = 1 , iy
            tbase(i,k,j) = ts00 + &
                      tlp*dlog((sps1%ps(j,i)*a(k)+ptop)*d_r100)
          end do
        end do
      end do
    end if
!
    zpbl(:,:) = 500.0D0  ! For Zeng Ocean Flux Scheme
!
    do j = 1 , jendx
      do i = 1 , iym1
        sts1%tg(j,i) = atm1%t(i,kz,j)/sps1%ps(j,i)
        sts2%tg(j,i) = atm2%t(i,kz,j)/sps2%ps(j,i)
        sfsta%tgbb(j,i) = atm2%t(i,kz,j)/sps2%ps(j,i)
      end do
    end do
    if ( ichem == 1 ) then
      ssw2da    = d_zero
      sdeltk2d  = d_zero
      sdelqk2d  = d_zero
      sfracv2d  = d_half
      sfracb2d  = d_half
      sfracs2d  = d_zero
      svegfrac2d = d_zero
    end if
#ifndef BAND
    if (debug_level > 2) call initdiag
#endif
!
!chem2
    if ( ichem == 1 ) then
!-----set tracer concs to 1 (kg/kg) initially. Must convert this to p*
!-----mixing ratio to compute tendencies:
!US       mass test zero concs init input for advection
!qhy      initial chia is 10ppt
!hy       set the initial tracer concentration 10ppt (1.e-11), 9/4/98
 
      do itr = 1 , ntr
        do k = 1 , kz
          do j = 1 , jendx
            do i = 1 , iym1
              chia(i,k,j,itr) = sps1%ps(j,i)*d_zero
              chib(i,k,j,itr) = sps2%ps(j,i)*d_zero
!                 chia(i,k,j,itr)=sps1%ps(j,i)*1.e-11
!                 chib(i,k,j,itr)=sps2%ps(j,i)*1.e-11
            end do
          end do
        end do
      end do
 
    end if
!chem2_
!
!------set rainc and rainnc equal to 0. initially
!
    sfsta%rainc  = d_zero
    sfsta%rainnc = d_zero
 
    if ( icup==4 .or. icup==99 .or. icup==98) then
      cbmf2d = d_zero
    end if
!
  else ! ifrest=.true.
!
!-----when ifrest=.true., read in the data saved from previous run
!       for large domain
!
    call read_savefile_part1(bdydate1)
!
    mtau = mtau + ktau
!
    if ( myid == 0 ) then
      print * , 'ozone profiles restart'
      do k = 1 , kzp1
        write (6,'(1x,7E12.4)') o3prof_io(3,k,3)
      end do
      appdat = tochar(idatex)
      print 99001 , nbdytime, ktau, appdat
    end if
!
#ifndef CLM
    if ( lakemod == 1 ) then
      call lakescatter
    endif
#endif

    if ( myid == 0 ) then
      do j = 1 , jx
        do k = 1 , kz
          do i = 1 , iy
            sav_0(i,k,j)      = ub0_io(i,k,j)
            sav_0(i,kz+k,j)   = vb0_io(i,k,j)
            sav_0(i,kz*2+k,j) = qb0_io(i,k,j)
            sav_0(i,kz*3+k,j) = tb0_io(i,k,j)
          end do
        end do
        do i = 1 , iy
          sav_0(i,kz*4+1,j) = ps0_io(i,j)
          sav_0(i,kz*4+2,j) = ts0_io(i,j)
        end do
      end do
    end if
    call mpi_scatter(sav_0,iy*(kz*4+2)*jxp,mpi_real8,  &
                     sav0, iy*(kz*4+2)*jxp,mpi_real8,  &
                     0,mycomm,ierr)
    do j = 1 , jendl
      do k = 1 , kz
        do i = 1 , iy
          xub%b0(i,k,j) = sav0(i,k,j)
          xvb%b0(i,k,j) = sav0(i,kz+k,j)
          xqb%b0(i,k,j) = sav0(i,kz*2+k,j)
          xtb%b0(i,k,j) = sav0(i,kz*3+k,j)
        end do
      end do
      do i = 1 , iy
        xpsb%b0(i,j) = sav0(i,kz*4+1,j)
        ts0(i,j) = sav0(i,kz*4+2,j)
      end do
    end do

    if ( myid == 0 ) then
      do j = 1 , jx
        do k = 1 , kz
          do i = 1 , iy
            sav_0(i,k,j)      = atm1_io%u(i,k,j)
            sav_0(i,kz+k,j)   = atm2_io%u(i,k,j)
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
    call mpi_scatter(sav_0,iy*(kz*4+2)*jxp,mpi_real8,  &
                     sav0, iy*(kz*4+2)*jxp,mpi_real8,  &
                     0,mycomm,ierr)
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
        sps1%ps(j,i) = sav0(i,kz*4+1,j)
        sps2%ps(j,i) = sav0(i,kz*4+2,j)
      end do
    end do
    if ( myid == 0 ) then
      do j = 1 , jx
        do k = 1 , kz
          do i = 1 , iy
            sav_0(i,k,j)      = atm1_io%t(i,k,j)
            sav_0(i,kz+k,j)   = atm2_io%t(i,k,j)
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
    call mpi_scatter(sav_0,iy*(kz*4+2)*jxp,mpi_real8,   &
                     sav0, iy*(kz*4+2)*jxp,mpi_real8,   &
                     0,mycomm,ierr)
    do j = 1 , jendl
      do k = 1 , kz
        do i = 1 , iy
          atm1%t(i,k,j)  = sav0(i,k,j)
          atm2%t(i,k,j)  = sav0(i,kz+k,j)
          atm1%qv(i,k,j) = sav0(i,kz*2+k,j)
          atm2%qv(i,k,j) = sav0(i,kz*3+k,j)
        end do
      end do
      do i = 1 , iy
        sts1%tg(j,i) = sav0(i,kz*4+1,j)
        sts2%tg(j,i) = sav0(i,kz*4+2,j)
      end do
    end do
    if ( myid == 0 ) then
      do j = 1 , jx
        do k = 1 , kz
          do i = 1 , iy
            sav_0(i,k,j)      = atm1_io%qc(i,k,j)
            sav_0(i,kz+k,j)   = atm2_io%qc(i,k,j)
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
    call mpi_scatter(sav_0,iy*(kz*4+2)*jxp,mpi_real8,  &
                     sav0, iy*(kz*4+2)*jxp,mpi_real8,  &
                     0,mycomm,ierr)
    do j = 1 , jendl
      do k = 1 , kz
        do i = 1 , iy
          atm1%qc(i,k,j) = sav0(i,k,j)
          atm2%qc(i,k,j) = sav0(i,kz+k,j)
          fcc(i,k,j)     = sav0(i,kz*2+k,j)
        end do
      end do
      do i = 1 , iy
        sfsta%rainc(j,i)  = sav0(i,kz*4+1,j)
        sfsta%rainnc(j,i) = sav0(i,kz*4+2,j)
      end do
    end do
    do j = 1 , jendx
      do k = 1 , kz
        do i = 1 , iym1
          heatrt(j,i,k) = sav0(i,kz*3+k,j)
        end do
      end do
    end do
!
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
!    
!     Begin scatter of the UW variables read in from the restart file
!
      if ( myid == 0 ) then
        do j = 1 , jxm1
          do k = 1 , kzp1
            do i = 1 , iy
               sav_0b(i,k,j) = atm1_io%tke(i,k,j)
            end do
          end do
        end do
      end if
      call mpi_scatter(sav_0b,iy*kzp1*jxp,mpi_real8,   &
                       sav0b, iy*kzp1*jxp,mpi_real8,   &
                       0,mycomm,ierr)
      do j = 1 , jendx
        do k = 1 , kzp1
          do i = 1 , iy
            atm1%tke(i,k,j) = sav0b(i,k,j)
          end do
        end do
      end do
      if ( myid == 0 ) then
        do j = 1 , jxm1
          do k = 1 , kzp1
            do i = 1 , iy
               sav_0b(i,k,j) = atm2_io%tke(i,k,j)
            end do
          end do
        end do
      end if
      call mpi_scatter(sav_0b,iy*kzp1*jxp,mpi_real8,   &
                       sav0b, iy*kzp1*jxp,mpi_real8,   &
                       0,mycomm,ierr)
      do j = 1 , jendx
        do k = 1 , kzp1
          do i = 1 , iy
            atm2%tke(i,k,j) = sav0b(i,k,j)
          end do
        end do
      end do

    end if ! ibltyp == 2 .or. ibltyp == 99
!
    if ( myid == 0 ) then
#ifdef BAND
      do j = 1 , jx
#else
      do j = 1 , jxm1
#endif
        do i = 1 , iym1
          var2d_0(i,j) = kpbl_io(i,j)
        end do
      end do
    end if
    call mpi_scatter(var2d_0,iy*jxp,mpi_integer, &
                     var2d0, iy*jxp,mpi_integer, &
                     0,mycomm,ierr)
    do j = 1 , jendx
      do i = 1 , iym1
        kpbl(j,i) = var2d0(i,j)
      end do
    end do
!
    if ( myid == 0 ) then
      do j = 1 , jx
        do i = 1 , iy
          sav_0a(i,1,j) = hfx_io(i,j)
          sav_0a(i,2,j) = qfx_io(i,j)
          sav_0a(i,3,j) = uvdrag_io(i,j)
          sav_0a(i,4,j) = tgbb_io(i,j)
        end do
      end do
#ifdef BAND
      do j = 1 , jx
#else
      do j = 1 , jxm1
#endif
        do k = 1 , kzp1
          do i = 1 , iym1
            sav_0a(i,4+k,j) = o3prof_io(i,k,j)
          end do
        end do
      end do
    end if
    allrec = 4 + kzp1
    call mpi_scatter(sav_0a,iy*allrec*jxp,mpi_real8,  &
                     sav0a, iy*allrec*jxp,mpi_real8,  &
                     0,mycomm,ierr)
    do j = 1 , jendl
      do i = 1 , iy
        sfsta%hfx(j,i)    = sav0a(i,1,j)
        sfsta%qfx(j,i)    = sav0a(i,2,j)
        sfsta%uvdrag(j,i) = sav0a(i,3,j)
        sfsta%tgbb(j,i)   = sav0a(i,4,j)
      end do
    end do
    do j = 1 , jendx
      do k = 1 , kzp1
        do i = 1 , iym1
          o3prof(j,i,k) = sav0a(i,4+k,j)
        end do
      end do
    end do
    if ( iocnflx == 2 ) then
      call mpi_scatter(zpbl_io, iy*jxp,mpi_real8, &
                       swapv,   iy*jxp,mpi_real8, &
                       0,mycomm,ierr)
      zpbl = transpose(swapv)
    end if
    if ( icup == 1 ) then
      if ( myid == 0 ) then
        do j = 1 , jx
          do k = 1 , kz
            do i = 1 , iy
              sav_0c(i,k,j)    = rsheat_io(i,k,j)
              sav_0c(i,kz+k,j) = rswat_io(i,k,j)
            end do
          end do
        end do
      end if
      call mpi_scatter(sav_0c,iy*kz*2*jxp,mpi_real8,  &
                       sav0c, iy*kz*2*jxp,mpi_real8,  &
                       0,mycomm,ierr)
      do j = 1 , jendl
        do k = 1 , kz
          do i = 1 , iy
            rsheat(j,i,k) = sav0c(i,k,j)
            rswat(j,i,k)  = sav0c(i,kz+k,j)
          end do
        end do
      end do
    end if
    if ( icup == 3 ) then
      if ( myid == 0 ) then
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
      call mpi_scatter(sav_0b,iy*(kzp1)*jxp,mpi_real8,  &
                       sav0b, iy*(kzp1)*jxp,mpi_real8,  &
                       0,mycomm,ierr)
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
      call mpi_scatter(cbmf2d_io,iy*jxp,mpi_real8,  &
                       cbmf2d,   iy*jxp,mpi_real8,  &
                       0,mycomm,ierr)
    end if
    if ( myid == 0 ) then
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
    call mpi_scatter(sav_1,iym1*allrec*jxp,mpi_real8,   &
                     sav1, iym1*allrec*jxp,mpi_real8,   &
                     0,mycomm,ierr)
    do j = 1 , jendx
      do l = 1 , 4
        do k = 1 , kz
          do i = 1 , iym1
            absnxt(j,i,k,l) = sav1(i,(l-1)*kz+k,j)
          end do
        end do
      end do
    end do
    allrec = kz*4
    do j = 1 , jendx
      do l = 1 , kzp1
        do k = 1 , kzp1
          do i = 1 , iym1
            abstot(j,i,k,l) = sav1(i,allrec+(l-1)*(kzp1)+k,j)
          end do
        end do
      end do
    end do
    allrec = allrec + (kzp1)*(kz+1)
    do j = 1 , jendx
      do k = 1 , kzp1
        do i = 1 , iym1
          emstot(j,i,k) = sav1(i,allrec+k,j)
        end do
      end do
    end do
    if ( myid == 0 ) then
#ifdef BAND
      do j = 1 , jx
#else
      do j = 1 , jxm1
#endif
        do n = 1 , nnsg
          do i = 1 , iym1
            sav_2(i,n,j) = taf_io(n,i,j)
            sav_2(i,nnsg+n,j) = tlef_io(n,i,j)
            sav_2(i,nnsg*2+n,j) = ssw_io(n,i,j)
            sav_2(i,nnsg*3+n,j) = rsw_io(n,i,j)
          end do
        end do
        do i = 1 , iym1
          sav_2(i,nnsg*5+1,j) = solis_io(i,j)
          sav_2(i,nnsg*5+2,j) = solvd_io(i,j)
          sav_2(i,nnsg*5+3,j) = solvs_io(i,j)
          sav_2(i,nnsg*5+4,j) = flw_io(i,j)
        end do
      end do
    end if
    allrec = nnsg*5 + 4
    call mpi_scatter(sav_2,iym1*allrec*jxp,mpi_real8,   &
                     sav2, iym1*allrec*jxp,mpi_real8,   &
                     0,mycomm,ierr)
    do j = 1 , jendx
      do n = 1 , nnsg
        do i = 1 , iym1
          taf(n,j,i) = sav2(i,n,j)
          tlef(n,j,i) = sav2(i,nnsg+n,j)
          ssw(n,j,i) = sav2(i,nnsg*2+n,j)
          rsw(n,j,i) = sav2(i,nnsg*3+n,j)
        end do
      end do
      do i = 1 , iym1
        solis(j,i) = sav2(i,nnsg*5+1,j)
        solvd(j,i) = sav2(i,nnsg*5+2,j)
        solvs(j,i) = sav2(i,nnsg*5+3,j)
        flw(j,i)   = sav2(i,nnsg*5+4,j)
      end do
    end do
    if ( myid == 0 ) then
#ifdef BAND
      do j = 1 , jx
#else
      do j = 1 , jxm1
#endif
        do n = 1 , nnsg
          do i = 1 , iym1
            sav_2(i,n,j) = tgbrd_io(n,i,j)
            sav_2(i,nnsg+n,j) = tsw_io(n,i,j)
            sav_2(i,nnsg*2+n,j) = sncv_io(n,i,j)
            sav_2(i,nnsg*3+n,j) = gwet_io(n,i,j)
            sav_2(i,nnsg*4+n,j) = tgrd_io(n,i,j)
          end do
        end do
        do i = 1 , iym1
          sav_2(i,nnsg*5+1,j) = flwd_io(i,j)
          sav_2(i,nnsg*5+2,j) = fsw_io(i,j)
          sav_2(i,nnsg*5+3,j) = sabveg_io(i,j)
          sav_2(i,nnsg*5+4,j) = sinc_io(i,j)
        end do
      end do
    end if
    allrec = nnsg*5 + 4
    call mpi_scatter(sav_2,iym1*allrec*jxp,mpi_real8,   &
                     sav2, iym1*allrec*jxp,mpi_real8,   &
                     0,mycomm,ierr)
    do j = 1 , jendx
      do n = 1 , nnsg
        do i = 1 , iym1
          tgbrd(n,j,i) = sav2(i,n,j)
          tsw(n,j,i) = sav2(i,nnsg+n,j)
          sncv(n,j,i) = sav2(i,nnsg*2+n,j)
          gwet(n,j,i) = sav2(i,nnsg*3+n,j)
          tgrd(n,j,i) = sav2(i,nnsg*4+n,j)
        end do
      end do
      do i = 1 , iym1
        flwd(j,i) = sav2(i,nnsg*5+1,j)
        fsw(j,i) = sav2(i,nnsg*5+2,j)
        sabveg(j,i) = sav2(i,nnsg*5+3,j)
        sinc(j,i) = sav2(i,nnsg*5+4,j)
      end do
    end do
    if ( myid == 0 ) then
#ifdef BAND
      do j = 1 , jx
#else
      do j = 1 , jxm1
#endif
        do n = 1 , nnsg
          do i = 1 , iym1
            sav_2(i,nnsg+n,j)   = snag_io(n,i,j)
            sav_2(i,nnsg*2+n,j) = sfice_io(n,i,j)
            sav_2(i,nnsg*3+n,j) = ldew_io(n,i,j)
            sav_2(i,nnsg*4+n,j) = emiss_io(n,i,j)
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
    call mpi_scatter(sav_2,iym1*allrec*jxp,mpi_real8,   &
                     sav2, iym1*allrec*jxp,mpi_real8,   &
                     0,mycomm,ierr)
    do j = 1 , jendx
      do n = 1 , nnsg
        do i = 1 , iym1
          snag(n,j,i)  = sav2(i,nnsg+n,j)
          sfice(n,j,i) = sav2(i,nnsg*2+n,j)
          ldew(n,j,i)  = sav2(i,nnsg*3+n,j)
          emiss(n,j,i) = sav2(i,nnsg*4+n,j)
        end do
      end do
      do i = 1 , iym1
        pptnc(j,i) = sav2(i,nnsg*5+1,j)
        pptc(j,i) = sav2(i,nnsg*5+2,j)
        prca(j,i) = sav2(i,nnsg*5+3,j)
        prnca(j,i) = sav2(i,nnsg*5+4,j)
      end do
    end do
    if ( myid == 0 ) then
#ifdef BAND
      do j = 1 , jx
#else
      do j = 1 , jxm1
#endif
        do n = 1 , nnsg
          do i = 1 , iym1
            sav_2a(i,n,j)      = veg2d1_io(n,i,j)
            sav_2a(i,nnsg+n,j) = ocld2d_io(n,i,j)
          end do
        end do
        do i = 1 , iym1
          sav_2a(i,nnsg*2+1,j) = veg2d_io(i,j)
          sav_2a(i,nnsg*2+2,j) = ldmsk_io(i,j)
        end do
      end do
    end if
    allrec = nnsg*2 + 2
    call mpi_scatter(sav_2a,iym1*allrec*jxp,mpi_integer,  &
                     sav2a, iym1*allrec*jxp,mpi_integer,  &
                     0,mycomm,ierr)
    do j = 1 , jendx
      do n = 1 , nnsg
        do i = 1 , iym1
          veg2d1(n,j,i) = sav2a(i,n,j)
          ocld2d(n,j,i) = sav2a(i,nnsg+n,j)
        end do
      end do
      do i = 1 , iym1
        veg2d(j,i) = sav2a(i,nnsg*2+1,j)
        ldmsk(j,i) = sav2a(i,nnsg*2+2,j)
      end do
    end do
    if ( ichem == 1 ) then
      if ( myid == 0 ) then
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
      call mpi_scatter(sav_4,iy*allrec*jxp,mpi_real8,   &
                       sav4, iy*allrec*jxp,mpi_real8,   &
                       0,mycomm,ierr)
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
      if ( myid == 0 ) then
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
      call mpi_scatter(sav_4a,iym1*7*jxp,mpi_real8,   &
                       sav4a, iym1*7*jxp,mpi_real8,   &
                       0,mycomm,ierr)
      do j = 1 , jendx
        do i = 1 , iym1
          ssw2da(j,i) = sav4a(i,1,j)
          sdeltk2d(j,i) = sav4a(i,2,j)
          sdelqk2d(j,i) = sav4a(i,3,j)
          sfracv2d(j,i) = sav4a(i,4,j)
          sfracb2d(j,i) = sav4a(i,5,j)
          sfracs2d(j,i) = sav4a(i,6,j)
          svegfrac2d(j,i) = sav4a(i,7,j)
        end do
      end do
    end if
#ifdef CLM
    if ( myid == 0 ) then
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
        end do
      end do
    end if
    call mpi_scatter(sav_clmout,iym1*8*jxp,mpi_real8,   &
                     sav_clmin, iym1*8*jxp,mpi_real8,   &
                     0,mycomm,ierr)
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
      end do
    end do
    call mpi_scatter(lndcat2d_io,iy*jxp,mpi_real8, &
                     swapv,      iy*jxp,mpi_real8, &
                     0,mycomm,ierr)
    lndcat2d = transpose(swapv)
#endif
    call mpi_bcast(ktau,1,mpi_integer,0,mycomm,ierr)
    call mpi_bcast(mtau,1,mpi_integer,0,mycomm,ierr)
    call mpi_bcast(nbdytime,1,mpi_integer,0,mycomm,ierr)
    call date_bcast(idatex,0,mycomm,ierr)
    xbctime = dble(nbdytime)
#ifndef BAND
    if (debug_level > 2) call mpidiag
#endif
    dt = dt2    ! First timestep successfully read in
    dtcum = dt2
    dtpbl = dt2
    dtche = dt2
    rdtpbl = d_one/dt2
    dttke = dt2
!
!-----end of initial/restart if test
!
  end if
!chem2
  if ( ichem == 1 ) then
    iso2  = 0
    iso4  = 0
    ibchl = 0
    ibchb = 0
    iochl = 0
    iochb = 0
    ibin  = 0
    do itr = 1 , ntr
      if ( chtrname(itr) == 'SO2' ) iso2 = itr
      if ( chtrname(itr) == 'SO4' ) then
        iso4 = itr
        ichso4 = itr
      end if
      if ( chtrname(itr) == 'BC_HL' ) then
        ibchl = itr
        ichbc = itr
      end if
      if ( chtrname(itr) == 'BC_HB' ) then
        ibchb = itr
        ichoc = itr
      end if
      if ( chtrname(itr) == 'OC_HL' ) iochl = itr
      if ( chtrname(itr) == 'OC_HB' ) iochb = itr
      if ( chtrname(itr) == 'DUST' ) then
        ibin = ibin + 1
        idust(ibin) = itr
      end if
    end do
  end if
!chem2_
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ****** initialize and define constants for vector bats
 
  if ( ktau == 0 ) call initb(jbegin,jendx,2,iym1)

  if ( iemiss == 1 .and. .not. ifrest ) then
    do j = 1 , jendx
      do i = 1 , iym1
        do n = 1 , nnsg
          ist = veg2d1(n,j,i)
          if ( ist == 14 .or. ist == 15 ) then
            emiss(n,j,i) = 0.955D0
          else if ( ist == 8 ) then
            emiss(n,j,i) = 0.76D0
          else if ( ist == 11 ) then
            emiss(n,j,i) = 0.85D0
          else if ( ist == 12 ) then
            emiss(n,j,i) = 0.97D0
          else
            emiss(n,j,i) = 0.99D0 - &
                    (albvgs(ist)+albvgl(ist))*0.1D0
          end if
!         emiss(n,j,i) = d_one
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
#ifdef CLM
  numdays = dayspy
#endif
  if (myid == 0) then
    write (6,*) 'Calculate solar declination angle at ',toint10(idatex)
  end if
  call solar1
#ifdef CLM
  init_grid = .true.
#endif
  call inirad
!
!-----calculating topographical correction to diffusion coefficient
  do i = 1 , iy
    do j = 1 , jendl
      hgfact(j,i) = d_one
    end do
  end do

  do i = 2 , iym2
    im1h = max0(i-1,2)
    ip1h = min0(i+1,iym2)
#ifdef BAND
    do j = jbegin , jendm
      jm1 = j-1
      jp1 = j+1
#else
    do j = jbegin , jendm
      if ( myid == 0 ) then
        jm1 = max0(j-1,2)
      else
        jm1 = j - 1
      end if
      if ( myid == nproc-1 ) then
        jp1 = min0(j+1,jxp-2)
      else
        jp1 = j + 1
      end if
#endif
      hg1 = dabs((mddom%ht(j,i)-mddom%ht(j,im1h))/dx)
      hg2 = dabs((mddom%ht(j,i)-mddom%ht(j,ip1h))/dx)
      hg3 = dabs((mddom%ht(j,i)-mddom%ht(jm1,i))/dx)
      hg4 = dabs((mddom%ht(j,i)-mddom%ht(jp1,i))/dx)
      hgmax = dmax1(hg1,hg2,hg3,hg4)*regrav
      hgfact(j,i) = d_one/(d_one+(hgmax/0.001D0)**d_two)
    end do
  end do
!
!-----set up output time:
!
#ifdef CLM
  if ( ifrest ) then
    ! CLM modifies landuse table. Get the modified one from
    ! restart file
    mddom%lndcat(:,:) = lndcat2d(:,:)
    do n = 1 , nnsg
      lndcat1(n,:,:) = lndcat2d(:,:)
    end do
  end if
#endif

99001 format (' ***** restart file for large domain at time = ', i8,   &
          ' seconds, ktau = ',i7,' date = ',a,' read in')
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
  if ( ktau == 0 ) then
    do k = 1 , kz
      do j = 1 , jendl
        do i = 1 , iym1
          heatrt(j,i,k) = d_zero
          o3prof(j,i,k) = d_zero
        end do
      end do
    end do
    do j = 1 , jendl
      do i = 1 , iym1
        o3prof(j,i,kzp1) = d_zero
      end do
    end do
    call o3data(1,jendx,1,iym1)
    if ( myid == 0 ) then
      write (6,*) 'ozone profiles'
      do k = 1 , kzp1
        write (6,'(1x,7E12.4)') o3prof(2,3,k)
      end do
    end if
    ! RRTM_SW gas / abs constant initialisation
    if ( irrtm == 1 ) then
      call rrtmg_sw_ini(cpd)
      call rrtmg_lw_ini(cpd)
    end if
  end if
 
  end subroutine inirad
!
end module mod_init
