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
  use mod_slice
  use mod_constants
#ifdef CLM
  use mod_clm
  use mod_lm_interface
  use clm_varsur , only : init_tgb , init_grid
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
  real(dp) :: hg1 , hg2 , hg3 , hg4 , hgmax
  character(len=32) :: appdat
  integer :: ierr
  character (len=64) :: subroutine_name='init'
  integer :: idindx = 0
!
  call time_begin(subroutine_name,idindx)
  !
  ! Reset the accumulation arrays
  !
  tgmx_o(:,:) = -1.E30
  tgmn_o(:,:) =  1.E30
  t2mx_o(:,:) = -1.E30
  t2mn_o(:,:) =  1.E30
  tavg_o(:,:) = 0.0
  pcpx_o(:,:) = -1.E30
  w10x_o(:,:) = -1.E30
  pcpa_o(:,:) = 0.0
  sund_o(:,:) = 0.0
  sunt_o(:,:) = 0.0
  psmn_o (:,:)=  1.E30
  !
  ! For an initial run -- not a restart
  !
  if ( .not. ifrest ) then
    !
    ! Initialize model atmospheric status variables
    ! Data are from the ICBC input at first timestep.
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
          atm1%qx(j,i,k,iqv) = xqb%b0(j,i,k)
          atm2%t(j,i,k) = xtb%b0(j,i,k)
          atm2%qx(j,i,k,iqv) = xqb%b0(j,i,k)
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
    !
    ! If we have activated SeaIce scheme, on ocean point we consider
    ! the temperature as the signal to cover with ice the sea, changing
    ! the tipe of soil to permanent ice. The landmask is:
    !
    !    0 -> Ocean
    !    1 -> Land
    !    2 -> Sea Ice
    !
    ! We have the grid (ldmsk) and subgrid (ldmsk1) versions of this
    !
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
                ldmsk1(n,j,i) = 2
              end do
            end if
          end if
        end do
      end do
    end if
    !
    ! Repeat the above for lake points if lake model is activated
    !
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
                ldmsk1(n,j,i) = 2
              end do
            end if
          end if
        end do
      end do
    end if
#endif
    !
    ! Initialize the tbase for BM cumulus scheme
    !
    if (icup == 3) then
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            tbase(j,i,k) = ts00 + tlp*dlog((sfs%psa(j,i)*hsigma(k)+ptop)*d_r100)
          end do
        end do
      end do
    end if
    !
    ! Initialize PBL Hgt
    !
    zpbl(:,:) = 500.0D0
    !
    ! Inizialize the surface atmospheric temperature
    !
    do i = ice1 , ice2
      do j = jce1 , jce2
        sfs%tgbb(j,i) = atm2%t(j,i,kz)/sfs%psb(j,i)
      end do
    end do
    !
    ! Initialize surface parameters for aerosol scheme
    !
    if ( ichem == 1 ) then
      sfracv2d(:,:)  = d_half
      sfracb2d(:,:)  = d_half
    end if
    !
    ! Set the TKE variables for UW PBL to a default value
    !
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      atm1%tke(:,:,:) = tkemin
      atm2%tke(:,:,:) = tkemin
    end if
    !
    ! Init the diurnal cycle SST scheme
    !
    if ( idcsst == 1 ) then
      tdeltas(:,:) = tground2(:,:) - 0.001D0
      deltas(:,:) = 0.001D0
    end if
    !
    ! Inizialize Ozone profiles
    !
    call o3data
    if ( myid == 0 ) then
      write (6,*) 'ozone profiles'
      do k = 1 , kzp1
        write (6,'(1x,7E12.4)') o3prof(3,3,k)
      end do
    end if
    !
    ! Diagnostic init
    !
    if ( .not. lband .and. debug_level > 2 ) then
      call initdiag
    end if
    !
    ! End of initial run case
    !
  else
    !
    ! When restarting, read in the data saved from previous run
    !
    call read_savefile(idate1)
    !
    ! Comunicate the data to other processors
    !
    xbctime = d_zero
    nbdytime = 0
    call mpi_bcast(ktau,1,mpi_integer8,0,mycomm,ierr)
    call date_bcast(idatex,0,mycomm,ierr)
!
    mtau = mtau + ktau
!
    call deco1_scatter(atm1_io%u,atm1%u,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(atm1_io%v,atm1%v,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(atm1_io%t,atm1%t,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(atm1_io%qx,atm1%qx, &
                       jcross1,jcross2,icross1,icross2,1,kz,1,nqx)

    call deco1_scatter(atm2_io%u,atm2%u,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(atm2_io%v,atm2%v,jdot1,jdot2,idot1,idot2,1,kz)
    call deco1_scatter(atm2_io%t,atm2%t,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(atm2_io%qx,atm2%qx, &
                       jcross1,jcross2,icross1,icross2,1,kz,1,nqx)

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

    call exchange(sfs%psa,1,jce1,jce2,ice1,ice2)
    call exchange(sfs%psb,1,jce1,jce2,ice1,ice2)

    if ( ipptls == 1 ) then
      call deco1_scatter(fcc_io,fcc,jcross1,jcross2,icross1,icross2,1,kz)
    end if
    call deco1_scatter(heatrt_io,heatrt,jcross1,jcross2,icross1,icross2,1,kz)
    call deco1_scatter(o3prof_io,o3prof,jcross1,jcross2,icross1,icross2,1,kzp1)
    if ( myid == 0 ) then
      print * , 'ozone profiles restart'
      do k = 1 , kzp1
        write (6,'(1x,7E12.4)') o3prof(3,3,k)
      end do
    end if

    ! Scatter of the UW variables read in from the restart file
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      call deco1_scatter(atm1_io%tke,atm1%tke, &
                         jcross1,jcross2,icross1,icross2,1,kzp1)
      call deco1_scatter(atm2_io%tke,atm2%tke, &
                         jcross1,jcross2,icross1,icross2,1,kzp1)
      call deco1_scatter(kpbl_io,kpbl,jcross1,jcross2,icross1,icross2)
    end if
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
    call subgrid_deco1_scatter(ldmsk1_io,ldmsk1,jcross1,jcross2,icross1,icross2)

    call deco1_scatter(solis_io,solis,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(solvd_io,solvd,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(solvs_io,solvs,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(sabveg_io,sabveg,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(flw_io,flw,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(flwd_io,flwd,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(fsw_io,fsw,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(sinc_io,sinc,jcross1,jcross2,icross1,icross2)
    call deco1_scatter(ldmsk_io,ldmsk,jcross1,jcross2,icross1,icross2)

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
    !
    ! CLM modifies landuse table. Get the modified one from restart file
    !
    mddom%lndcat(jci1:jci2,ici1:ici2) = lndcat2d(jci1:jci2,ici1:ici2)
    do n = 1 , nnsg
      lndcat1(n,jci1:jci2,ici1:ici2) = lndcat2d(jci1:jci2,ici1:ici2)
    end do
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

    if ( idcsst == 1 ) then
      call deco1_scatter(dtskin_io,dtskin,jcross1,jcross2,icross1,icross2)
      call deco1_scatter(deltas_io,deltas,jcross1,jcross2,icross1,icross2)
      call deco1_scatter(tdeltas_io,tdeltas,jcross1,jcross2,icross1,icross2)
    end if
    !
    ! Update ground temperature on Ocean/Lakes
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( iswater(mddom%lndcat(j,i)) ) then
          if (idcsst == 1) then
            sfs%tga(j,i) = ts1(j,i) + dtskin(j,i)
            sfs%tgb(j,i) = ts1(j,i) + dtskin(j,i)
          else
            sfs%tga(j,i) = ts1(j,i)
            sfs%tgb(j,i) = ts1(j,i)
          end if
          if ( iseaice == 1 ) then
            if ( lakemod == 1 .and. islake(mddom%lndcat(j,i)) ) cycle
            if ( ts1(j,i) <= icetemp ) then
              sfs%tga(j,i) = icetemp
              sfs%tgb(j,i) = icetemp
              ts1(j,i) = icetemp
              ldmsk(j,i) = 2
              do n = 1, nnsg
                ldmsk1(n,j,i) = 2
                sfice(n,j,i) = d_1000
                sncv(n,j,i) = d_zero
              end do
            else
              sfs%tga(j,i) = ts1(j,i)
              sfs%tgb(j,i) = ts1(j,i)
              ldmsk(j,i) = 0
              do n = 1, nnsg
                ldmsk1(n,j,i) = 0
                sfice(n,j,i) = d_zero
                sncv(n,j,i)  = d_zero
              end do
            end if
          end if
        end if
      end do
    end do

    call deco1_scatter(dstor_io,dstor,jdot1,jdot2,idot1,idot2,1,nsplit)
    call deco1_scatter(hstor_io,hstor,jdot1,jdot2,idot1,idot2,1,nsplit)
!
    if ( .not. lband .and. debug_level > 2 ) then
      call mpidiag
    end if
    !
    ! Setup all timeseps for a restart
    !
    dt = dt2
    dtcum = dt2
    dtche = dt2
    dtpbl = dt2
    rdtpbl = d_one/dt2
    dttke = dt2
    !
    ! Report success
    !
    if ( myid == 0 ) then
      appdat = tochar(idatex)
      print 99001 , appdat
    end if
    !
    ! End of restart phase
    !
  end if
  !
  ! The following allows to change landuse on restart.
  !
  do i = ici1 , ici2
    do j = jci1 , jci2
      iveg(j,i) = idnint(lndcat(j,i))
      do n = 1 , nnsg
        iveg1(n,j,i) = idnint(lndcat1(n,j,i))
      end do
    end do
  end do
  !
  ! Initialize the BATS variable (Used also by CLM)
  !
  if ( ktau == 0 ) then
    call initb
    !
#ifndef CLM
    if ( lakemod == 1 ) then
      call subgrid_deco1_gather(idep,idep_io,jcross1,jcross2,icross1,icross2)
    end if
#endif
  end if
#ifdef CLM
  call mkslice
  call initclm(ifrest,idate1,idate2,dx,dtrad,dtsrf)
#endif
  !
  ! Calculate emission coefficients
  !
  if ( iemiss == 1 .and. ktau == 0 ) then
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
            emiss(n,j,i) = 0.99D0-(albvgs(ist)+albvgl(ist))*0.1D0
          end if
        end do
      end do
    end do
  end if
  !  
  ! Calculate topographical correction to diffusion coefficient
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
  ! pressure of tropopause
  !
  do i = ici1 , ici2
    do j = jci1 , jci2
      ptrop(j,i) = 250.0D2 - 150.0D2*dcos(mddom%xlat(j,i)*degrad)**d_two
    end do
  end do
  !
  ! RRTM_SW gas / abs constant initialisation
  !
  if ( irrtm == 1 ) then
    call rrtmg_sw_ini(cpd)
    call rrtmg_lw_ini(cpd)
  end if
!
! chemistry initialisation
!
  if ( ichem == 1 ) then
    call start_chem(ifrest,idate1,intbdy,dtbdys)
  end if
!
  call time_end(subroutine_name,idindx)
!
! Formats for printout
!
99001 format ('Successfully read restart file at time = ', a)
!
  end subroutine init
!
end module mod_init
