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
  use mod_slabocean
  use rrtmg_sw_init
  use rrtmg_lw_init
  use mod_pbl_interface
  use mod_precip
  use mod_bdycod
  use mod_mpmessage
  use mod_sun
  use mod_ncio
  use mod_savefile
  use mod_slice
  use mod_constants
  use mod_outvars
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
  real(rk8) , parameter :: tlp = 50.0D0
  real(rk8) , parameter :: ts00 = 288.0D0
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
  implicit none
!
  integer(ik4) :: i , j , k , n , ist
  real(rk8) :: hg1 , hg2 , hg3 , hg4 , hgmax
  real(rk8) :: sfice_temp
  character(len=32) :: appdat
#ifdef DEBUG
  character(len=dbgslen) :: subroutine_name = 'init'
  integer(ik4) , save :: idindx = 0
  call time_begin(subroutine_name,idindx)
#endif
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
    sfice_temp = icetemp
    if ( iseaice == 1 ) then
      do i = ice1 , ice2
        do j = jce1 , jce2
          if ( isocean(mddom%lndcat(j,i)) ) then
            if ( ts0(j,i) <= icetemp ) then
              sfs%tga(j,i) = sfice_temp
              sfs%tgb(j,i) = sfice_temp
              ts0(j,i) = icetemp
              ldmsk(j,i) = 2
              do n = 1, nnsg
                ldmsk1(n,j,i) = 2
                sfice(n,j,i) = d_10
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
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( islake(mddom%lndcat(j,i)) ) then
            if ( ts0(j,i) <= icetemp ) then
              sfs%tga(j,i) = icetemp
              sfs%tgb(j,i) = icetemp
              ts0(j,i) = icetemp
              ldmsk(j,i) = 2
              do n = 1, nnsg
                ldmsk1(n,j,i) = 2
                sfice(n,j,i) = d_10
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
    if ( icup == 3 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
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
    do i = ici1 , ici2
      do j = jci1 , jci2
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
      dtskin(:,:) = d_zero
      deltas(:,:) = 0.001D0
      sst(:,:) = ts0(jci1:jci2,ici1:ici2)
      tdeltas(:,:) = sst(:,:) - deltas(:,:)
    end if
    !
    ! Inizialize Ozone profiles
    !
    call o3data
    if ( myid == italk ) then
      call vprntv(o3prof(3,3,:),kzp1,'Ozone profile at (3,3)')
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
    call bcast(ktau)
    call bcast(idatex)
!
    mtau = mtau + ktau
!
    call grid_distribute(atm1_u_io,atm1%u,jde1,jde2,ide1,ide2,1,kz)
    call grid_distribute(atm1_v_io,atm1%v,jde1,jde2,ide1,ide2,1,kz)
    call grid_distribute(atm1_t_io,atm1%t,jce1,jce2,ice1,ice2,1,kz)
    call grid_distribute(atm1_qx_io,atm1%qx,jce1,jce2,ice1,ice2,1,kz,1,nqx)

    call grid_distribute(atm2_u_io,atm2%u,jde1,jde2,ide1,ide2,1,kz)
    call grid_distribute(atm2_v_io,atm2%v,jde1,jde2,ide1,ide2,1,kz)
    call grid_distribute(atm2_t_io,atm2%t,jce1,jce2,ice1,ice2,1,kz)
    call grid_distribute(atm2_qx_io,atm2%qx,jce1,jce2,ice1,ice2,1,kz,1,nqx)

    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      call grid_distribute(atm1_tke_io,atm1%tke,jce1,jce2,ice1,ice2,1,kzp1)
      call grid_distribute(atm2_tke_io,atm2%tke,jce1,jce2,ice1,ice2,1,kzp1)
      call grid_distribute(kpbl_io,kpbl,jci1,jci2,ici1,ici2)
    end if

    call grid_distribute(psa_io,sfs%psa,jce1,jce2,ice1,ice2)
    call grid_distribute(psb_io,sfs%psb,jce1,jce2,ice1,ice2)
    call grid_distribute(tga_io,sfs%tga,jce1,jce2,ice1,ice2)
    call grid_distribute(tgb_io,sfs%tgb,jce1,jce2,ice1,ice2)

    call grid_distribute(hfx_io,sfs%hfx,jci1,jci2,ici1,ici2)
    call grid_distribute(qfx_io,sfs%qfx,jci1,jci2,ici1,ici2)
    call grid_distribute(rainc_io,sfs%rainc,jci1,jci2,ici1,ici2)
    call grid_distribute(rainnc_io,sfs%rainnc,jci1,jci2,ici1,ici2)
    call grid_distribute(tgbb_io,sfs%tgbb,jci1,jci2,ici1,ici2)
    call grid_distribute(uvdrag_io,sfs%uvdrag,jci1,jci2,ici1,ici2)

    call exchange(sfs%psa,1,jce1,jce2,ice1,ice2)
    call exchange(sfs%psb,1,jce1,jce2,ice1,ice2)

    if ( ipptls > 0 ) then
      call grid_distribute(fcc_io,fcc,jci1,jci2,ici1,ici2,1,kz)
      if ( ipptls == 2 ) then
        call grid_distribute(snownc_io,sfs%snownc,jci1,jci2,ici1,ici2)
      end if
    end if
    call grid_distribute(heatrt_io,heatrt,jci1,jci2,ici1,ici2,1,kz)
    call grid_distribute(o3prof_io,o3prof,jci1,jci2,ici1,ici2,1,kzp1)
!
    if ( myid == italk ) then
      call vprntv(o3prof(3,3,:),kzp1,'Ozone profiles restart')
    end if
!
    if ( iocnflx == 2 ) then
      call grid_distribute(zpbl_io,zpbl,jci1,jci2,ici1,ici2)
    end if
    if ( icup == 1 ) then
      call grid_distribute(rsheat_io,rsheat,jci1,jci2,ici1,ici2,1,kz)
      call grid_distribute(rswat_io,rswat,jci1,jci2,ici1,ici2,1,kz)
    end if
    if ( icup == 3 ) then
      call grid_distribute(tbase_io,tbase,jci1,jci2,ici1,ici2,1,kz)
      call grid_distribute(cldefi_io,cldefi,jci1,jci2,ici1,ici2)
    end if
    if ( icup == 4 .or. icup == 99 .or. icup == 98 .or. icup == 97 ) then
      call grid_distribute(cbmf2d_io,cbmf2d,jci1,jci2,ici1,ici2)
    end if

    if ( irrtm == 0 ) then 
      call grid_distribute(gasabsnxt_io,gasabsnxt,jci1,jci2,ici1,ici2,1,kz,1,4)
      call grid_distribute(gasabstot_io,gasabstot, &
                           jci1,jci2,ici1,ici2,1,kzp1,1,kzp1)
      call grid_distribute(gasemstot_io,gasemstot,jci1,jci2,ici1,ici2,1,kzp1)
    end if

    call subgrid_distribute(tlef_io,tlef,jci1,jci2,ici1,ici2)
    call subgrid_distribute(ssw_io,ssw,jci1,jci2,ici1,ici2)
    call subgrid_distribute(rsw_io,rsw,jci1,jci2,ici1,ici2)
    call subgrid_distribute(tgrd_io,tgrd,jci1,jci2,ici1,ici2)
    call subgrid_distribute(tgbrd_io,tgbrd,jci1,jci2,ici1,ici2)
    call subgrid_distribute(sncv_io,sncv,jci1,jci2,ici1,ici2)
    call subgrid_distribute(gwet_io,gwet,jci1,jci2,ici1,ici2)
    call subgrid_distribute(snag_io,snag,jci1,jci2,ici1,ici2)
    call subgrid_distribute(sfice_io,sfice,jci1,jci2,ici1,ici2)
    call subgrid_distribute(ldew_io,ldew,jci1,jci2,ici1,ici2)
    call subgrid_distribute(taf_io,taf,jci1,jci2,ici1,ici2)
    call subgrid_distribute(tsw_io,tsw,jci1,jci2,ici1,ici2)
    call subgrid_distribute(emiss_io,emiss,jci1,jci2,ici1,ici2)
    call subgrid_distribute(ldmsk1_io,ldmsk1,jci1,jci2,ici1,ici2)

    call grid_distribute(solis_io,solis,jci1,jci2,ici1,ici2)
    call grid_distribute(solvd_io,solvd,jci1,jci2,ici1,ici2)
    call grid_distribute(solvs_io,solvs,jci1,jci2,ici1,ici2)
    call grid_distribute(sabveg_io,sabveg,jci1,jci2,ici1,ici2)
    call grid_distribute(flw_io,flw,jci1,jci2,ici1,ici2)
    call grid_distribute(flwd_io,flwd,jci1,jci2,ici1,ici2)
    call grid_distribute(fsw_io,fsw,jci1,jci2,ici1,ici2)
    call grid_distribute(sinc_io,sinc,jci1,jci2,ici1,ici2)
    call grid_distribute(ldmsk_io,ldmsk,jci1,jci2,ici1,ici2)

#ifndef CLM
    if ( lakemod == 1 ) then
      call subgrid_distribute(eta_io,xlake,jci1,jci2,ici1,ici2)
      call lake_fillvar(var_eta,xlake,1)
      call subgrid_distribute(hi_io,xlake,jci1,jci2,ici1,ici2)
      call lake_fillvar(var_hi,xlake,1)
      call subgrid_distribute(aveice_io,xlake,jci1,jci2,ici1,ici2)
      call lake_fillvar(var_aveice,xlake,1)
      call subgrid_distribute(hsnow_io,xlake,jci1,jci2,ici1,ici2)
      call lake_fillvar(var_hsnow,xlake,1)
      call subgrid_distribute(tlak_io,tlake,jci1,jci2,ici1,ici2,1,ndpmax)
      call lake_fillvar(var_tlak,tlake,1)
    endif
#else
    call grid_distribute(sols2d_io,sols2d,jci1,jci2,ici1,ici2)
    call grid_distribute(soll2d_io,soll2d,jci1,jci2,ici1,ici2)
    call grid_distribute(solsd2d_io,solsd2d,jci1,jci2,ici1,ici2)
    call grid_distribute(solld2d_io,solld2d,jci1,jci2,ici1,ici2)
    call grid_distribute(lndcat2d_io,lndcat2d,jci1,jci2,ici1,ici2)
    !
    ! CLM modifies landuse table. Get the modified one from restart file
    !
    if ( imask == 2 ) then
      mddom%lndcat(jci1:jci2,ici1:ici2) = lndcat2d(jci1:jci2,ici1:ici2)
      do n = 1 , nnsg
        lndcat1(n,jci1:jci2,ici1:ici2) = lndcat2d(jci1:jci2,ici1:ici2)
      end do
    end if
#endif
!
    if ( idcsst == 1 ) then
      call grid_distribute(dtskin_io,dtskin,jci1,jci2,ici1,ici2)
      call grid_distribute(deltas_io,deltas,jci1,jci2,ici1,ici2)
      call grid_distribute(tdeltas_io,tdeltas,jci1,jci2,ici1,ici2)
    end if

    call grid_distribute(dstor_io,dstor,jde1,jde2,ide1,ide2,1,nsplit)
    call grid_distribute(hstor_io,hstor,jde1,jde2,ide1,ide2,1,nsplit)

    if ( ichem == 1 ) then
      call grid_distribute(chia_io,chia,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      call grid_distribute(chib_io,chib,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      call grid_distribute(remlsc_io,remlsc,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      call grid_distribute(remcvc_io,remcvc,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      call grid_distribute(remdrd_io,remdrd,jce1,jce2,ice1,ice2,1,ntr)
      if ( igaschem == 1 .and. ichsolver > 0 ) then
        call grid_distribute(chemall_io,chemall,jci1,jci2,ici1,ici2, &
                              1,kz,1,totsp)
        call grid_distribute(taucldsp_io,taucldsp,jci1,jci2,ici1,ici2, &
                              0,kz,1,nspi)
      end if

      call grid_distribute(ssw2da_io,ssw2da,jci1,jci2,ici1,ici2)
      call grid_distribute(sdeltk2d_io,sdeltk2d,jci1,jci2,ici1,ici2)
      call grid_distribute(sdelqk2d_io,sdelqk2d,jci1,jci2,ici1,ici2)
      call grid_distribute(sfracv2d_io,sfracv2d,jci1,jci2,ici1,ici2)
      call grid_distribute(sfracb2d_io,sfracb2d,jci1,jci2,ici1,ici2)
      call grid_distribute(sfracs2d_io,sfracs2d,jci1,jci2,ici1,ici2)
      call grid_distribute(svegfrac2d_io,svegfrac2d,jci1,jci2,ici1,ici2)
    end if

    if ( islab_ocean == 1 .and. do_restore_sst ) then
      call grid_distribute(qflux_restore_sst_io,qflux_restore_sst, &
        jci1,jci2,ici1,ici2,1,12)
      call bcast(stepcount)
    end if
    !
    ! Restore ground temperature on Ocean/Lakes
    !
    if ( islab_ocean == 0 ) then
      sfice_temp = icetemp
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( iswater(mddom%lndcat(j,i)) ) then
            if ( ldmsk(j,i) == 0 ) then
              if ( idcsst == 1 ) then
                sst(j,i) = ts1(j,i)
                sfs%tga(j,i) = sst(j,i) + dtskin(j,i)
                sfs%tgb(j,i) = sst(j,i) + dtskin(j,i)
              else
                ! Update temperature where NO ice
                sfs%tga(j,i) = ts1(j,i)
                sfs%tgb(j,i) = ts1(j,i)
              end if
            end if
            if ( iseaice == 1 ) then
              if ( lakemod == 1 .and. islake(mddom%lndcat(j,i)) ) cycle
              if ( ts1(j,i) <= icetemp .and. ldmsk(j,i) == 0 ) then
                sfs%tga(j,i) = sfice_temp
                sfs%tgb(j,i) = sfice_temp
                ts1(j,i) = icetemp
                ldmsk(j,i) = 2
                do n = 1, nnsg
                  ldmsk1(n,j,i) = 2
                  sfice(n,j,i) = d_10
                end do
              else if ( ts1(j,i) > icetemp .and. ldmsk(j,i) == 2 ) then
                sfs%tga(j,i) = ts1(j,i)
                sfs%tgb(j,i) = ts1(j,i)
                ! Decrease the surface ice to melt it
                do n = 1, nnsg
                  sfice(n,j,i) = sfice(n,j,i)*d_r10
                end do
              end if
            end if
          end if
        end do
      end do
    end if
    !
    ! Setup all timeseps for a restart
    !
    dtbat = dt*dble(ntsrf)
    dt = dt2
    rdt = d_one/dt
    dtsq = dt*dt
    dtcb = dt*dt*dt
    !
    ! Report success
    !
    if ( myid == italk ) then
      appdat = tochar(idatex)
      write(stdout,*) 'Successfully read restart file at time = ', appdat
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
  end if
#ifdef CLM
  call mkslice
  call initclm(ifrest,idate1,idate2,dx,dtrad,dtsrf,igaschem,iaerosol,chtrname)
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
      hgfact(j,i) = d_one/(d_one+(hgmax/0.001D0)**2)
    end do
  end do
  !
  ! pressure of tropopause
  !
  do i = ici1 , ici2
    do j = jci1 , jci2
      ptrop(j,i) = 250.0D2 - 150.0D2*dcos(mddom%xlat(j,i)*degrad)**2
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
    call start_chem
  end if
!
#ifdef DEBUG
  call time_end(subroutine_name,idindx)
#endif
  end subroutine init
!
end module mod_init
