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
    call grid_distribute(atm1_io%u,atm1%u,jde1,jde2,ide1,ide2,1,kz)
    call grid_distribute(atm1_io%v,atm1%v,jde1,jde2,ide1,ide2,1,kz)
    call grid_distribute(atm1_io%t,atm1%t,jce1,jce2,ice1,ice2,1,kz)
    call grid_distribute(atm1_io%qx,atm1%qx,jce1,jce2,ice1,ice2,1,kz,1,nqx)

    call grid_distribute(atm2_io%u,atm2%u,jde1,jde2,ide1,ide2,1,kz)
    call grid_distribute(atm2_io%v,atm2%v,jde1,jde2,ide1,ide2,1,kz)
    call grid_distribute(atm2_io%t,atm2%t,jce1,jce2,ice1,ice2,1,kz)
    call grid_distribute(atm2_io%qx,atm2%qx,jce1,jce2,ice1,ice2,1,kz,1,nqx)

    call grid_distribute(sfs_io%psa,sfs%psa,jce1,jce2,ice1,ice2)
    call grid_distribute(sfs_io%psb,sfs%psb,jce1,jce2,ice1,ice2)
    call grid_distribute(sfs_io%tga,sfs%tga,jce1,jce2,ice1,ice2)
    call grid_distribute(sfs_io%tgb,sfs%tgb,jce1,jce2,ice1,ice2)
    call grid_distribute(sfs_io%hfx,sfs%hfx,jce1,jce2,ice1,ice2)
    call grid_distribute(sfs_io%qfx,sfs%qfx,jce1,jce2,ice1,ice2)
    call grid_distribute(sfs_io%rainc,sfs%rainc,jce1,jce2,ice1,ice2)
    call grid_distribute(sfs_io%rainnc,sfs%rainnc,jce1,jce2,ice1,ice2)
    call grid_distribute(sfs_io%uvdrag,sfs%uvdrag,jce1,jce2,ice1,ice2)
    call grid_distribute(sfs_io%tgbb,sfs%tgbb,jce1,jce2,ice1,ice2)

    call exchange(sfs%psa,1,jce1,jce2,ice1,ice2)
    call exchange(sfs%psb,1,jce1,jce2,ice1,ice2)

    if ( ipptls == 1 ) then
      call grid_distribute(fcc_io,fcc,jce1,jce2,ice1,ice2,1,kz)
    end if
    call grid_distribute(heatrt_io,heatrt,jce1,jce2,ice1,ice2,1,kz)
    call grid_distribute(o3prof_io,o3prof,jce1,jce2,ice1,ice2,1,kzp1)
    if ( myid == 0 ) then
      print * , 'ozone profiles restart'
      do k = 1 , kzp1
        write (6,'(1x,7E12.4)') o3prof(3,3,k)
      end do
    end if

    ! Scatter of the UW variables read in from the restart file
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      call grid_distribute(atm1_io%tke,atm1%tke,jce1,jce2,ice1,ice2,1,kzp1)
      call grid_distribute(atm2_io%tke,atm2%tke,jce1,jce2,ice1,ice2,1,kzp1)
      call grid_distribute(kpbl_io,kpbl,jce1,jce2,ice1,ice2)
    end if
!
    if ( iocnflx == 2 ) then
      call grid_distribute(zpbl_io,zpbl,jce1,jce2,ice1,ice2)
    end if
    if ( icup == 1 ) then
      call grid_distribute(rsheat_io,rsheat,jce1,jce2,ice1,ice2,1,kz)
      call grid_distribute(rswat_io,rswat,jce1,jce2,ice1,ice2,1,kz)
    end if
    if ( icup == 3 ) then
      call grid_distribute(tbase_io,tbase,jce1,jce2,ice1,ice2,1,kz)
      call grid_distribute(cldefi_io,cldefi,jce1,jce2,ice1,ice2)
    end if
    if ( icup==4 .or. icup==99 .or. icup==98 ) then
      call grid_distribute(cbmf2d_io,cbmf2d,jce1,jce2,ice1,ice2)
    end if

    if ( irrtm == 0 ) then 
      call grid_distribute(gasabsnxt_io,gasabsnxt,jce1,jce2,ice1,ice2,1,kz,1,4)
      call grid_distribute(gasabstot_io,gasabstot, &
                           jce1,jce2,ice1,ice2,1,kzp1,1,kzp1)
      call grid_distribute(gasemstot_io,gasemstot,jce1,jce2,ice1,ice2,1,kzp1)
    end if ! irrtm test

    call subgrid_distribute(tlef_io,tlef,jce1,jce2,ice1,ice2)
    call subgrid_distribute(ssw_io,ssw,jce1,jce2,ice1,ice2)
    call subgrid_distribute(rsw_io,rsw,jce1,jce2,ice1,ice2)
    call subgrid_distribute(tgrd_io,tgrd,jce1,jce2,ice1,ice2)
    call subgrid_distribute(tgbrd_io,tgbrd,jce1,jce2,ice1,ice2)
    call subgrid_distribute(sncv_io,sncv,jce1,jce2,ice1,ice2)
    call subgrid_distribute(gwet_io,gwet,jce1,jce2,ice1,ice2)
    call subgrid_distribute(snag_io,snag,jce1,jce2,ice1,ice2)
    call subgrid_distribute(sfice_io,sfice,jce1,jce2,ice1,ice2)
    call subgrid_distribute(ldew_io,ldew,jce1,jce2,ice1,ice2)
    call subgrid_distribute(taf_io,taf,jce1,jce2,ice1,ice2)
    call subgrid_distribute(tsw_io,tsw,jce1,jce2,ice1,ice2)
    call subgrid_distribute(emiss_io,emiss,jce1,jce2,ice1,ice2)
    call subgrid_distribute(ldmsk1_io,ldmsk1,jce1,jce2,ice1,ice2)

    call grid_distribute(solis_io,solis,jce1,jce2,ice1,ice2)
    call grid_distribute(solvd_io,solvd,jce1,jce2,ice1,ice2)
    call grid_distribute(solvs_io,solvs,jce1,jce2,ice1,ice2)
    call grid_distribute(sabveg_io,sabveg,jce1,jce2,ice1,ice2)
    call grid_distribute(flw_io,flw,jce1,jce2,ice1,ice2)
    call grid_distribute(flwd_io,flwd,jce1,jce2,ice1,ice2)
    call grid_distribute(fsw_io,fsw,jce1,jce2,ice1,ice2)
    call grid_distribute(sinc_io,sinc,jce1,jce2,ice1,ice2)
    call grid_distribute(ldmsk_io,ldmsk,jce1,jce2,ice1,ice2)

#ifndef CLM
    if ( lakemod == 1 ) then
      call subgrid_distribute(eta_io,eta,jce1,jce2,ice1,ice2)
      call subgrid_distribute(hi_io,hi,jce1,jce2,ice1,ice2)
      call subgrid_distribute(aveice_io,aveice,jce1,jce2,ice1,ice2)
      call subgrid_distribute(hsnow_io,hsnow,jce1,jce2,ice1,ice2)
      call subgrid_distribute(tlak_io,tlak,jce1,jce2,ice1,ice2,1,ndpmax)
    endif
#else
    call grid_distribute(sols2d_io,sols2d,jce1,jce2,ice1,ice2)
    call grid_distribute(soll2d_io,soll2d,jce1,jce2,ice1,ice2)
    call grid_distribute(solsd2d_io,solsd2d,jce1,jce2,ice1,ice2)
    call grid_distribute(solld2d_io,solld2d,jce1,jce2,ice1,ice2)
    call grid_distribute(aldirs2d_io,aldirs2d,jce1,jce2,ice1,ice2)
    call grid_distribute(aldirl2d_io,aldirl2d,jce1,jce2,ice1,ice2)
    call grid_distribute(aldifs2d_io,aldifs2d,jce1,jce2,ice1,ice2)
    call grid_distribute(aldifl2d_io,aldifl2d,jce1,jce2,ice1,ice2)
    call grid_distribute(lndcat2d_io,lndcat2d,jce1,jce2,ice1,ice2)
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
      call grid_distribute(chia_io,chia,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      call grid_distribute(chib_io,chib,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      call grid_distribute(remlsc_io,remlsc,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      call grid_distribute(remcvc_io,remcvc,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      call grid_distribute(remdrd_io,remdrd,jce1,jce2,ice1,ice2,1,kz)
      call grid_distribute(ssw2da_io,ssw2da,jce1,jce2,ice1,ice2)
      call grid_distribute(sdeltk2d_io,sdeltk2d,jce1,jce2,ice1,ice2)
      call grid_distribute(sdelqk2d_io,sdelqk2d,jce1,jce2,ice1,ice2)
      call grid_distribute(sfracv2d_io,sfracv2d,jce1,jce2,ice1,ice2)
      call grid_distribute(sfracb2d_io,sfracb2d,jce1,jce2,ice1,ice2)
      call grid_distribute(sfracs2d_io,sfracs2d,jce1,jce2,ice1,ice2)
      call grid_distribute(svegfrac2d_io,svegfrac2d,jce1,jce2,ice1,ice2)
    end if

    if ( idcsst == 1 ) then
      call grid_distribute(dtskin_io,dtskin,jce1,jce2,ice1,ice2)
      call grid_distribute(deltas_io,deltas,jce1,jce2,ice1,ice2)
      call grid_distribute(tdeltas_io,tdeltas,jce1,jce2,ice1,ice2)
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
!
    call grid_distribute(dstor_io,dstor,jde1,jde2,ide1,ide2,1,nsplit)
    call grid_distribute(hstor_io,hstor,jde1,jde2,ide1,ide2,1,nsplit)
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
      call subgrid_collect(idep,idep_io,jce1,jce2,ice1,ice2)
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
