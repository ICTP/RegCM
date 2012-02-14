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
  real(dp) , parameter :: tlp = 50.0D0
  real(dp) , parameter :: ts00 = 288.0D0
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
  integer :: i , j , k , n , ist
  type (rcm_time_and_date) :: icbc_date
  real(dp) :: hg1 , hg2 , hg3 , hg4 , hgmax
  character(len=32) :: appdat
  integer :: mmrec , ierr

  tgmx_o(:,:) = -1.E30
  tgmn_o(:,:) =  1.E30
  t2mx_o(:,:) = -1.E30
  t2mn_o(:,:) =  1.E30
  tavg_o(:,:) = 0.0
  pcpx_o(:,:) = -1.E30
  w10x_o(:,:) = -1.E30
  pcpa_o(:,:) = 0.0
  sund_o(:,:) = 0.0
  psmn_o (:,:)=  1.E30

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
  !
  ! for initial run--not using restart
  !
  if ( .not. ifrest ) then
    !
    ! set rainwater and cloud water equal to zero initially.
    !
    atm1%qc(:,:,:) = d_zero
    atm2%qc(:,:,:) = d_zero
!
    if ( ichem == 1 ) then
      !qhy tchie, tchitb(replace tchidp:deposition)
      !    initialize removal terms
      remlsc = d_zero
      remcvc = d_zero
      rxsg   = d_zero
      rxsaq1 = d_zero
      rxsaq2 = d_zero
      remdrd = d_zero
      wdlsc  = d_zero
    end if
    !
    ! set the variables related to blackadar pbl equal to 0 initially.
    !
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      atm1%tke(:,:,:) = tkemin
      atm2%tke(:,:,:) = tkemin
    end if
    !
    if ( icup == 1 ) then
      rsheat(:,:,:) = d_zero
      rswat(:,:,:)  = d_zero
    end if
    !
    ! Read in the initial conditions for large domain:
    ! the initial conditions are the output from PREPROC/ICBC.
    !
#ifdef CLM
    if ( .not. allocated(init_tgb) ) allocate(init_tgb(iy,jx))
#endif
    if ( myid == 0 ) then
      mmrec = icbc_search(bdydate1)
      if (mmrec < 0) then
        appdat = tochar(bdydate2)
        call fatal(__FILE__,__LINE__,'ICBC for '//appdat//' not found')
      end if
      call read_icbc(ps0_io,ts0_io,ub0_io,vb0_io,tb0_io,qb0_io)
      appdat = tochar(bdydate1)
      write (6,*) 'READY IC DATA for ', appdat
      !
      ! Convert surface pressure to pstar
      !
      ps0_io(:,:) = ps0_io(:,:)*d_r10 - ptop
    end if
!
    call deco1_scatter(ub0_io,xub%b0,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(vb0_io,xvb%b0,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(tb0_io,xtb%b0,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(qb0_io,xqb%b0,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(ps0_io,xpsb%b0,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(ts0_io,ts0,jcross1,jcross2,icross1,icross2)
    !
    ! this piece of code determines p(.) from p(x) by a 4-point
    ! interpolation. on the x-grid, a p(x) point outside the grid
    ! domain is assumed to satisfy p(j,0)=p(j,1); p(j,iy)=p(j,iym1);
    ! and similarly for the i's.
    !
    call deco1_exchange_left(xpsb%b0,1,icross1,icross2)
    call deco1_exchange_right(xpsb%b0,1,icross1,icross2)
    call psc2psd(xpsb%b0,psdot)
    !
    ! Couple pressure u,v,t,q
    !
    do k = 1 , kz
      do i = ide1 , ide2
        do j = jde1 , jde2
          xub%b0(j,i,k) = xub%b0(j,i,k)*psdot(j,i)
          xvb%b0(j,i,k) = xvb%b0(j,i,k)*psdot(j,i)
        end do
      end do
    end do
    call deco1_exchange_left(xub%b0,1,ide1,ide2,1,kz)
    call deco1_exchange_right(xub%b0,1,ide1,ide2,1,kz)
    call deco1_exchange_left(xvb%b0,1,ide1,ide2,1,kz)
    call deco1_exchange_right(xvb%b0,1,ide1,ide2,1,kz)
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          xtb%b0(j,i,k) = xtb%b0(j,i,k)*xpsb%b0(j,i)
          xqb%b0(j,i,k) = xqb%b0(j,i,k)*xpsb%b0(j,i)
        end do
      end do
    end do
    call deco1_exchange_left(xtb%b0,1,ice1,ice2,1,kz)
    call deco1_exchange_right(xtb%b0,1,ice1,ice2,1,kz)
    call deco1_exchange_left(xqb%b0,1,ice1,ice2,1,kz)
    call deco1_exchange_right(xqb%b0,1,ice1,ice2,1,kz)
    !
    ! Initialize model atmospheric status variables
    !
    do k = 1 , kz
      do i = ide1 , ide2
        do j = jde1 , jde2
          atm1%u(j,i,k) = xub%b0(j,i,k)
          atm1%v(j,i,k) = xvb%b0(j,i,k)
          atm2%u(j,i,k) = xub%b0(j,i,k)
          atm2%v(j,i,k) = xvb%b0(j,i,k)
        end do
      end do
    end do
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          atm1%t(j,i,k) = xtb%b0(j,i,k)
          atm2%t(j,i,k) = xtb%b0(j,i,k)
          atm1%qv(j,i,k) = xqb%b0(j,i,k)
          atm2%qv(j,i,k) = xqb%b0(j,i,k)
        end do
      end do
    end do
    do i = ice1 , ice2
      do j = jce1 , jce2
        sfs%psa(j,i) = xpsb%b0(j,i)
        sfs%psb(j,i) = xpsb%b0(j,i)
        sfs%tga(j,i) = ts0(j,i)
        sfs%tgb(j,i) = ts0(j,i)
      end do
    end do
    if ( iseaice == 1 ) then
      do i = ice1 , ice2
        do j = jce1 , jce2
          if ( isocean(mddom%lndcat(j,i)) ) then
            if ( ts0(j,i) <= icetemp ) then
              sfs%tga(j,i) = icetemp
              sfs%tgb(j,i) = icetemp
              ts0(j,i) = icetemp
              ldmsk(j,i) = 2
              do n = 1, nnsg
                ocld(n,j,i) = 2
              end do
            end if
          end if
        end do
      end do
    end if
#ifndef CLM
    if ( lakemod == 1 ) then
      do i = ice1 , ice2
        do j = jce1 , jce2
          if ( islake(mddom%lndcat(j,i)) ) then
            if ( ts0(j,i) <= icetemp ) then
              sfs%tga(j,i) = icetemp
              sfs%tgb(j,i) = icetemp
              ts0(j,i) = icetemp
              ldmsk(j,i) = 2
              do n = 1, nnsg
                ocld(n,j,i) = 2
              end do
            end if
          end if
        end do
      end do
    end if
#endif
    if (icup == 3) then
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            tbase(j,i,k) = ts00 + &
                      tlp*dlog((sfs%psa(j,i)*a(k)+ptop)*d_r100)
          end do
        end do
      end do
    end if
!
    zpbl(:,:) = 500.0D0  ! For Zeng Ocean Flux Scheme
    do i = ice1 , ice2
      do j = jce1 , jce2
        sfs%tga(j,i) = atm1%t(j,i,kz)/sfs%psa(j,i)
        sfs%tgb(j,i) = atm2%t(j,i,kz)/sfs%psb(j,i)
        sfs%tgbb(j,i) = atm2%t(j,i,kz)/sfs%psb(j,i)
      end do
    end do
    if ( ichem == 1 ) then
      ssw2da(:,:)    = d_zero
      sdeltk2d(:,:)  = d_zero
      sdelqk2d(:,:)  = d_zero
      sfracv2d(:,:)  = d_half
      sfracb2d(:,:)  = d_half
      sfracs2d(:,:)  = d_zero
      svegfrac2d(:,:) = d_zero
    end if
#ifndef BAND
    if (debug_level > 2) call initdiag
#endif
!
    sfs%rainc(:,:)  = d_zero
    sfs%rainnc(:,:) = d_zero
 
    if ( icup==4 .or. icup==99 .or. icup==98) then
      cbmf2d(:,:) = d_zero
    end if
!
  else ! ifrest=.true.
!
!-----when ifrest=.true., read in the data saved from previous run
!       for large domain
!
    call read_savefile(bdydate1)
!
    call mpi_bcast(ktau,1,mpi_integer8,0,mycomm,ierr)
    call mpi_bcast(mtau,1,mpi_integer8,0,mycomm,ierr)
    call mpi_bcast(nbdytime,1,mpi_integer8,0,mycomm,ierr)
    call date_bcast(idatex,0,mycomm,ierr)

    mtau = mtau + ktau
    xbctime = dble(nbdytime)
!
    call deco1_scatter(ub1_io,xub%b1,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(vb1_io,xvb%b1,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(tb1_io,xtb%b1,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(qb1_io,xqb%b1,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(ps1_io,xpsb%b1,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(ts1_io,ts1,jcross1,jcross2,icross1,icross2)

    call deco1_exchange_left(xub%b1,1,idot1,idot2,1,kz)
    call deco1_exchange_right(xub%b1,1,idot1,idot2,1,kz)
    call deco1_exchange_left(xvb%b1,1,idot1,idot2,1,kz)
    call deco1_exchange_right(xvb%b1,1,idot1,idot2,1,kz)
    call deco1_exchange_left(xtb%b1,1,icross1,icross2,1,kz)
    call deco1_exchange_right(xtb%b1,1,icross1,icross2,1,kz)
    call deco1_exchange_left(xqb%b1,1,icross1,icross2,1,kz)
    call deco1_exchange_right(xqb%b1,1,icross1,icross2,1,kz)
    call deco1_exchange_left(xpsb%b1,1,icross1,icross2)
    call deco1_exchange_right(xpsb%b1,1,icross1,icross2)

    call deco1_scatter(atm1_io%u,atm1%u,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(atm1_io%v,atm1%v,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(atm1_io%t,atm1%t,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(atm1_io%qv,atm1%qv,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(atm1_io%qc,atm1%qc,jcross1,jcross2,icross1,icross2,1,kz)

    call deco1_scatter(atm2_io%u,atm2%u,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(atm2_io%v,atm2%v,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(atm2_io%t,atm2%t,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(atm2_io%qv,atm2%qv,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(atm2_io%qc,atm2%qc,jcross1,jcross2,icross1,icross2,1,kz)

    call deco1_scatter(sfs_io%psa,sfs%psa,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(sfs_io%psb,sfs%psb,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(sfs_io%tga,sfs%tga,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(sfs_io%tgb,sfs%tgb,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(sfs_io%hfx,sfs%hfx,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(sfs_io%qfx,sfs%qfx,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(sfs_io%rainc,sfs%rainc,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(sfs_io%rainnc,sfs%rainnc,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(sfs_io%uvdrag,sfs%uvdrag,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(sfs_io%tgbb,sfs%tgbb,jcross1,jcross2,icross1,icross2)

    call deco1_exchange_left(sfs%psa,1,icross1,icross2)
    call deco1_exchange_right(sfs%psa,1,icross1,icross2)
    call deco1_exchange_left(sfs%psb,1,icross1,icross2)
    call deco1_exchange_right(sfs%psb,1,icross1,icross2)

    if ( ipptls == 1 ) then
      call deco1_scatter(fcc_io,fcc,jcross1,jcross2,icross1,icross2,1,kz)
    end if
    call deco1_scatter(heatrt_io,heatrt,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(o3prof_io,o3prof,jcross1,jcross2,icross1,icross2,1,kzp1)

    ! Scatter of the UW variables read in from the restart file
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      call deco1_scatter(atm1_io%tke,atm1%tke, &
                         jcross1,jcross2,icross1,icross2,1,kz)
      call deco1_scatter(atm2_io%tke,atm2%tke, &
                         jcross1,jcross2,icross1,icross2,1,kz)
    end if ! ibltyp == 2 .or. ibltyp == 99
!
    if ( iocnflx == 2 ) then
      call deco1_scatter(zpbl_io,zpbl,jcross1,jcross2,icross1,icross2)
    end if
    if ( icup == 1 ) then
      call deco1_scatter(rsheat_io,rsheat,jcross1,jcross2,icross1,icross2,1,kz)
      call deco1_scatter(rswat_io,rswat,jcross1,jcross2,icross1,icross2,1,kz)
    end if
    if ( icup == 3 ) then
      call deco1_scatter(tbase_io,tbase,jcross1,jcross2,icross1,icross2,1,kz)
      call deco1_scatter(cldefi_io,cldefi,jcross1,jcross2,icross1,icross2)
    end if
    if ( icup==4 .or. icup==99 .or. icup==98 ) then
      call deco1_scatter(cbmf2d_io,cbmf2d,jcross1,jcross2,icross1,icross2)
    end if

    if ( irrtm == 0 ) then 
      call deco1_scatter(gasabsnxt_io,gasabsnxt, &
                         jcross1,jcross2,icross1,icross2,1,kz,1,4)
      call deco1_scatter(gasabstot_io,gasabstot, &
                         jcross1,jcross2,icross1,icross2,1,kzp1,1,kzp1)
      call deco1_scatter(gasemstot_io,gasemstot, &
                         jcross1,jcross2,icross1,icross2,1,kzp1)
    end if ! irrtm test

    call subgrid_deco1_scatter(tlef_io,tlef,jcross1,jcross2,icross1,icross2)
    call subgrid_deco1_scatter(ssw_io,ssw,jcross1,jcross2,icross1,icross2)
    call subgrid_deco1_scatter(rsw_io,rsw,jcross1,jcross2,icross1,icross2)
    call subgrid_deco1_scatter(tgrd_io,tgrd,jcross1,jcross2,icross1,icross2)
    call subgrid_deco1_scatter(tgbrd_io,tgbrd,jcross1,jcross2,icross1,icross2)
    call subgrid_deco1_scatter(sncv_io,sncv,jcross1,jcross2,icross1,icross2)
    call subgrid_deco1_scatter(gwet_io,gwet,jcross1,jcross2,icross1,icross2)
    call subgrid_deco1_scatter(snag_io,snag,jcross1,jcross2,icross1,icross2)
    call subgrid_deco1_scatter(sfice_io,sfice,jcross1,jcross2,icross1,icross2)
    call subgrid_deco1_scatter(ldew_io,ldew,jcross1,jcross2,icross1,icross2)
    call subgrid_deco1_scatter(taf_io,taf,jcross1,jcross2,icross1,icross2)
    call subgrid_deco1_scatter(tsw_io,tsw,jcross1,jcross2,icross1,icross2)
    call subgrid_deco1_scatter(emiss_io,emiss,jcross1,jcross2,icross1,icross2)
    call subgrid_deco1_scatter(ocld_io,ocld,jcross1,jcross2,icross1,icross2)

    call deco1_scatter(kpbl_io,kpbl,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(solis_io,solis,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(solvd_io,solvd,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(solvs_io,solvs,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(sabveg_io,sabveg,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(flw_io,flw,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(flwd_io,flwd,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(fsw_io,fsw,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(sinc_io,sinc,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(pptnc_io,pptnc,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(pptc_io,pptc,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(ldmsk_io,ldmsk,jcross1,jcross2,icross1,icross2)

    if ( iseaice == 1 .or. lakemod == 1 ) then
      do i = ice1 , ice2
        do j = jce1 , jce2
          do n = 1 , nnsg
            if ( ocld(n,j,i) == 2 ) iveg1(n,j,i) = 12
          end do
        end do
      end do
    end if

#ifndef CLM
    if ( lakemod == 1 ) then
      call subgrid_deco1_scatter(eta_io,eta,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_scatter(hi_io,hi,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_scatter(aveice_io,aveice, &
                                 jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_scatter(hsnow_io,hsnow,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_scatter(tlak_io,tlak, &
                                 jcross1,jcross2,icross1,icross2,1,ndpmax)
    endif
#else
    call deco1_scatter(sols2d_io,sols2d,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(soll2d_io,soll2d,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(solsd2d_io,solsd2d,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(solld2d_io,solld2d,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(aldirs2d_io,aldirs2d,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(aldirl2d_io,aldirl2d,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(aldifs2d_io,aldifs2d,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(aldifl2d_io,aldifl2d,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(lndcat2d_io,lndcat2d,jcross1,jcross2,icross1,icross2)
#endif
!
    if ( ichem == 1 ) then
      call deco1_scatter(chia_io,chia, &
                         jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
      call deco1_scatter(chib_io,chib, &
                         jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
      call deco1_scatter(remlsc_io,remlsc, &
                         jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
      call deco1_scatter(remcvc_io,remcvc, &
                         jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
      call deco1_scatter(remdrd_io,remdrd, &
                         jcross1,jcross2,icross1,icross2,1,ntr)
      call deco1_scatter(ssw2da_io,ssw2da,jcross1,jcross2,icross1,icross2)
      call deco1_scatter(sdeltk2d_io,sdeltk2d,jcross1,jcross2,icross1,icross2)
      call deco1_scatter(sdelqk2d_io,sdelqk2d,jcross1,jcross2,icross1,icross2)
      call deco1_scatter(sfracv2d_io,sfracv2d,jcross1,jcross2,icross1,icross2)
      call deco1_scatter(sfracb2d_io,sfracb2d,jcross1,jcross2,icross1,icross2)
      call deco1_scatter(sfracs2d_io,sfracs2d,jcross1,jcross2,icross1,icross2)
      call deco1_scatter(svegfrac2d_io,svegfrac2d, &
                         jcross1,jcross2,icross1,icross2)
    end if

    call deco1_scatter(fbat_io,fbat, &
                       jout1,jout2,iout1,iout2,numbat-numsts+1,numbat)

    call deco1_scatter(dstor_io,dstor,jdot1,jdot2,idot1,idot2,1,nsplit)
    call deco1_scatter(hstor_io,hstor,jdot1,jdot2,idot1,idot2,1,nsplit)
!
    call deco1_scatter(sue_io,sue,jdot1,jdot2,1,kz)
    call deco1_scatter(sui_io,sui,jdot1,jdot2,1,kz)
    call deco1_scatter(nue_io,nue,jdot1,jdot2,1,kz)
    call deco1_scatter(nui_io,nui,jdot1,jdot2,1,kz)
    call deco1_scatter(sve_io,sve,jdot1,jdot2,1,kz)
    call deco1_scatter(svi_io,svi,jdot1,jdot2,1,kz)
    call deco1_scatter(nve_io,nve,jdot1,jdot2,1,kz)
    call deco1_scatter(nvi_io,nvi,jdot1,jdot2,1,kz)
!
    call deco1_exchange_left(sue,1,1,kz)
    call deco1_exchange_right(sue,1,1,kz)
    call deco1_exchange_left(sui,1,1,kz)
    call deco1_exchange_right(sui,1,1,kz)
    call deco1_exchange_left(nue,1,1,kz)
    call deco1_exchange_right(nue,1,1,kz)
    call deco1_exchange_left(nui,1,1,kz)
    call deco1_exchange_right(nui,1,1,kz)
    call deco1_exchange_left(sve,1,1,kz)
    call deco1_exchange_right(sve,1,1,kz)
    call deco1_exchange_left(svi,1,1,kz)
    call deco1_exchange_right(svi,1,1,kz)
    call deco1_exchange_left(nve,1,1,kz)
    call deco1_exchange_right(nve,1,1,kz)
    call deco1_exchange_left(nvi,1,1,kz)
    call deco1_exchange_right(nvi,1,1,kz)
#ifndef BAND
    call mpi_bcast(eui,nidot*kz,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(eue,nidot*kz,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(evi,nidot*kz,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(eve,nidot*kz,mpi_real8,0,mycomm,ierr)
    if (debug_level > 2) call mpidiag
#endif

    dt = dt2    ! First timestep successfully read in
    dtcum = dt2
    dtche = dt2
    dtpbl = dt2
    rdtpbl = d_one/dt2
    dttke = dt2

    if ( myid == 0 ) then
      print * , 'ozone profiles restart'
      do k = 1 , kzp1
        write (6,'(1x,7E12.4)') o3prof_io(3,3,k)
      end do
      appdat = tochar(idatex)
      print 99001 , nbdytime, ktau, appdat
    end if
!
!-----end of initial/restart if test
!
  end if
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ****** initialize and define constants for vector bats

  !
  ! The following allows to change landuse on restart.
  !
  do i = ice1 , ice2
    do j = jce1 , jce2
      iveg(j,i) = idnint(lndcat(j,i)+0.1D0)
      do n = 1 , nnsg
        iveg1(n,j,i) = idnint(lndcat1(n,j,i)+0.1D0)
      end do
    end do
  end do
  
  if ( ktau == 0 ) then
    call initb(jci1,jci2,ici1,ici2)
    if ( lakemod == 1 ) then
      call subgrid_deco1_gather(idep,idep_io,jcross1,jcross2,icross1,icross2)
    end if
  end if

  if ( iemiss == 1 .and. .not. ifrest ) then
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          ist = iveg1(n,j,i)
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
  ! calculating topographical correction to diffusion coefficient
  !
  do i = ide1 , ide2
    do j = jde1 , jde2
      hgfact(j,i) = d_one
    end do
  end do

  do i = idi1 , idi2
    do j = jdi1 , jdi2
      hg1 = dabs((mddom%ht(j,i)-mddom%ht(j,i-1))/dx)
      hg2 = dabs((mddom%ht(j,i)-mddom%ht(j,i+1))/dx)
      hg3 = dabs((mddom%ht(j,i)-mddom%ht(j-1,i))/dx)
      hg4 = dabs((mddom%ht(j,i)-mddom%ht(j+1,i))/dx)
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
  integer :: k
!
  if ( ktau == 0 ) then
    heatrt(:,:,:) = d_zero
    o3prof(:,:,:) = d_zero
    call o3data(jci1,jci2,ici1,ici2)
    if ( myid == 0 ) then
      write (6,*) 'ozone profiles'
      do k = 1 , kzp1
        write (6,'(1x,7E12.4)') o3prof(3,3,k)
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
