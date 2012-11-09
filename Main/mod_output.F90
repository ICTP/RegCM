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
 
module mod_output

  use mod_runparams
  use mod_header
  use mod_mpmessage
  use mod_mppparam
  use mod_service
  use mod_atm_interface
  use mod_che_interface
  use mod_che_output
  use mod_lm_interface
  use mod_rad_interface
  use mod_cu_interface
  use mod_pbl_interface
  use mod_ncio
  use mod_ncout
  use mod_bdycod
  use mod_precip
  use mod_split
  use mod_savefile
  use mod_mppio
#ifdef CLM
  use mod_clm
#endif

  private

  integer(ik4) :: iolak
  logical :: lskipsrf , lskiprad , lskipche

  public :: output , mkfile

  data iolak /1/
  data lskipsrf /.false./
  data lskiprad /.false./
  data lskipche /.false./

  contains

  subroutine output

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine handles all of the output                       c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  implicit none
!
  logical :: ldoatm , ldosrf , ldorad , ldoche , ldosav , ldotmp
  logical :: lstartup
#ifdef DEBUG
  character(len=dbgslen) :: subroutine_name = 'output'
  integer(ik4) , save :: idindx = 0
  call time_begin(subroutine_name,idindx)
#endif
!
  lstartup = .false.
  if ( ktau == 0 .or. doing_restart ) then
    call newoutfiles(idatex)
    if ( myid == iocpu ) then
      call mkfile
    end if
    lstartup = .true.
    if ( doing_restart ) then
      doing_restart = .false.
#ifdef DEBUG
      call time_end(subroutine_name,idindx) 
#endif
      return
    end if
  end if
!
  ldoatm = .false.
  ldosrf = .false.
  ldorad = .false.
  ldoche = .false.
  ldosav = .false.
  ldotmp = .false.

  if ( ktau > 0 ) then
    if ( mod(ktau,ksav) == 0 ) then
      ldotmp = .true.
    end if
    if ( ( idatex == idate2 .or. &
         (lfdomonth(idatex) .and. lmidnight(idatex))) ) then
      ldosav = .true.
      ldotmp = .false.
    end if
    if ( mod(ktau,katm) == 0 ) then
      ldoatm = .true.
    end if
    if ( mod(ktau,ksrf) == 0 ) then
      ldosrf = .true.
    end if
    if ( mod(ktau,krad) == 0 ) then
      ldorad = .true.
    end if
    if ( mod(ktau,kche) == 0 ) then
      ldoche = .true.
    end if
  end if

  if ( lskipsrf ) then
    lskipsrf = .false.
    ldosrf = .true.
  end if
  if ( lskiprad ) then
    lskiprad = .false.
    ldorad = .true.
  end if
  if ( lskipche ) then
    lskipche = .false.
    ldoche = .true.
  end if
!
  if ( ktau == 0 ) then
    if ( debug_level > 2 ) then
      ldoatm = .true.
      lskipsrf = .true.
      lskiprad = .true.
      lskipche = .true.
    end if
  end if
!
!-----output for dataflow analyses:
!
  if ( ifatm ) then
    if ( ldoatm ) then
!
!=======================================================================
!     gather  ua,va,ta,qxa,rainc,rainnc,tgbrd,tsw,olcd2d
!=======================================================================
!
      call grid_collect(atm1%u,atm1_io%u,jde1,jde2,ide1,ide2,1,kz)
      call grid_collect(atm1%v,atm1_io%v,jde1,jde2,ide1,ide2,1,kz)
      call grid_collect(atm1%t,atm1_io%t,jce1,jce2,ice1,ice2,1,kz)
      call grid_collect(atm1%qx,atm1_io%qx,jce1,jce2,ice1,ice2,1,kz,1,nqx)
      call grid_collect(omega,omega_io,jci1,jci2,ici1,ici2,1,kz)
      call grid_collect(sfs%psa,sfs_io%psa,jce1,jce2,ice1,ice2)
      call grid_collect(sfs%rainc,sfs_io%rainc,jci1,jci2,ici1,ici2)
      call grid_collect(sfs%rainnc,sfs_io%rainnc,jci1,jci2,ici1,ici2)
      call subgrid_collect(tgbrd,tgbrd_io,jci1,jci2,ici1,ici2)
      call subgrid_collect(tsw,tsw_io,jci1,jci2,ici1,ici2)
      call subgrid_collect(ldmsk1,ldmsk1_io,jci1,jci2,ici1,ici2)
!
!=======================================================================
!     gather UW Scheme variables
!=======================================================================
!
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call grid_collect(atm1%tke,atm1_io%tke,jce1,jce2,ice1,ice2,1,kzp1)
        call grid_collect(uwstateb%kth,tcmstate_io%kth, &
                          jci1,jci2,ici1,ici2,1,kzp1)
        call grid_collect(uwstateb%kzm,tcmstate_io%kzm, &
                          jci1,jci2,ici1,ici2,1,kzp1)
      end if
!
      if ( myid == iocpu ) then
        call outatm
      end if
      sfs%rainc(:,:)  = d_zero
      sfs%rainnc(:,:) = d_zero
    end if
  end if
 
  ! Call surface output
 
  if ( ifsrf ) then
    if ( ldosrf ) then
#ifndef CLM
      if ( lakemod == 1 .and. iflak .and. iolak == klak ) then
        call subgrid_collect(eta,eta_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(hi,hi_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(aveice,aveice_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(evpr,evl_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(hsnow,hsnow_io,jci1,jci2,ici1,ici2)
#ifndef CLM
        call subgrid_collect(tlak,tlak_io,jci1,jci2,ici1,ici2,1,ndpmax)
        call subgrid_collect(iveg1,iveg1_io,jci1,jci2,ici1,ici2)
        call grid_collect(iveg,iveg_io,jci1,jci2,ici1,ici2)
#endif
      end if
#endif
      if ( iseaice == 1 .or. lakemod == 1 ) then
        call grid_collect(ldmsk,ldmsk_io,jci1,jci2,ici1,ici2)
      end if
      call grid_collect(fbat,fbat_io,jci1,jci2,ici1,ici2,1,numbat)
      if ( myid == iocpu ) then
        call outsrf
      end if
      fbat(:,:,sunt_o) = 0.0
      if ( ifsts .and. mod(ktau+kstsoff,ksts) == 0 .and. &
           ktau > kstsoff+2 ) then
        fbat(:,:,tgmx_o) = -1.E30
        fbat(:,:,tgmn_o) =  1.E30
        fbat(:,:,t2mx_o) = -1.E30
        fbat(:,:,t2mn_o) =  1.E30
        fbat(:,:,tavg_o) = 0.0
        fbat(:,:,pcpx_o) = -1.E30
        fbat(:,:,w10x_o) = -1.E30
        fbat(:,:,pcpa_o) = 0.0
        fbat(:,:,sund_o) = 0.0
        fbat(:,:,psmn_o) = 1.E30
      end if
      if ( lakemod == 1 .and. iflak .and. iolak == klak) then
        iolak = 1
      else
        iolak = iolak + 1
      end if

      if ( ifsub .and. nsg > 1 ) then
        call subgrid_collect(fsub,fsub_io,jci1,jci2,ici1,ici2,1,numsub)
#ifndef CLM
        if ( lakemod == 1 ) then
          call subgrid_collect(tlak,tlak_io,jci1,jci2,ici1,ici2,1,4)
          call subgrid_collect(iveg1,iveg1_io,jci1,jci2,ici1,ici2)
        end if
#endif
        if ( myid == iocpu ) then
          call outsub
        end if
      end if
    end if
  end if
! 
!     Call radiation output
! 
  if ( ifrad ) then
    if ( ldorad ) then
      call grid_collect(frad2d,frad2d_io,jci1,jci2,ici1,ici2,1,nrad2d)
      call grid_collect(frad3d,frad3d_io,jci1,jci2,ici1,ici2,1,kz,1,nrad3d)
      call grid_collect(sfs%psa,sfs_io%psa,jci1,jci2,ici1,ici2)
      if ( myid == iocpu ) then
        call outrad
      end if
    end if
  end if
!
!     Call chem output
! 
  if ( ifchem ) then
    if ( ldoche ) then
      ! call output for chemistry/aerosol
      ! aerosol optical properties be passed in t interface since
      ! they are declared and calculated in RAD module but outputed
      ! in CHE module   
      call output_chem(idatex,aerext,aerssa,aerasp,aertarf, &
                       aersrrf,aertalwrf,aersrlwrf )
   end if
 end if
!
!-----output for restart:
!
  if ( ifsave ) then
    if ( ldosav .or. ldotmp ) then
      call grid_collect(atm1%u,atm1_io%u,jde1,jde2,ide1,ide2,1,kz)
      call grid_collect(atm1%v,atm1_io%v,jde1,jde2,ide1,ide2,1,kz)
      call grid_collect(atm1%t,atm1_io%t,jce1,jce2,ice1,ice2,1,kz)
      call grid_collect(atm1%qx,atm1_io%qx,jce1,jce2,ice1,ice2,1,kz,1,nqx)

      call grid_collect(atm2%u,atm2_io%u,jde1,jde2,ide1,ide2,1,kz)
      call grid_collect(atm2%v,atm2_io%v,jde1,jde2,ide1,ide2,1,kz)
      call grid_collect(atm2%t,atm2_io%t,jce1,jce2,ice1,ice2,1,kz)
      call grid_collect(atm2%qx,atm2_io%qx,jce1,jce2,ice1,ice2,1,kz,1,nqx)

      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call grid_collect(atm1%tke,atm1_io%tke,jce1,jce2,ice1,ice2,1,kzp1)
        call grid_collect(atm2%tke,atm2_io%tke,jce1,jce2,ice1,ice2,1,kzp1)
        call grid_collect(kpbl,kpbl_io,jci1,jci2,ici1,ici2)
      end if

      call grid_collect(sfs%psa,sfs_io%psa,jce1,jce2,ice1,ice2)
      call grid_collect(sfs%psb,sfs_io%psb,jce1,jce2,ice1,ice2)
      call grid_collect(sfs%tga,sfs_io%tga,jce1,jce2,ice1,ice2)
      call grid_collect(sfs%tgb,sfs_io%tgb,jce1,jce2,ice1,ice2)

      call grid_collect(sfs%hfx,sfs_io%hfx,jci1,jci2,ici1,ici2)
      call grid_collect(sfs%qfx,sfs_io%qfx,jci1,jci2,ici1,ici2)
      call grid_collect(sfs%rainc,sfs_io%rainc,jci1,jci2,ici1,ici2)
      call grid_collect(sfs%rainnc,sfs_io%rainnc,jci1,jci2,ici1,ici2)
      call grid_collect(sfs%tgbb,sfs_io%tgbb,jci1,jci2,ici1,ici2)
      call grid_collect(sfs%uvdrag,sfs_io%uvdrag,jci1,jci2,ici1,ici2)

      if ( ipptls == 1 ) then
        call grid_collect(fcc,fcc_io,jci1,jci2,ici1,ici2,1,kz)
      end if
      call grid_collect(heatrt,heatrt_io,jci1,jci2,ici1,ici2,1,kz)
      call grid_collect(o3prof,o3prof_io,jci1,jci2,ici1,ici2,1,kzp1)

      if ( iocnflx == 2 ) then
        call grid_collect(zpbl,zpbl_io,jci1,jci2,ici1,ici2)
      end if
      if ( icup == 1 ) then
        call grid_collect(rsheat,rsheat_io,jci1,jci2,ici1,ici2,1,kz)
        call grid_collect(rswat,rswat_io,jci1,jci2,ici1,ici2,1,kz)
      end if
      if ( icup == 3 ) then
        call grid_collect(tbase,tbase_io,jci1,jci2,ici1,ici2,1,kz)
        call grid_collect(cldefi,cldefi_io,jci1,jci2,ici1,ici2)
      end if
      if ( icup==4 .or. icup==99 .or. icup==98 ) then
        call grid_collect(cbmf2d,cbmf2d_io,jci1,jci2,ici1,ici2)
      end if

      if ( irrtm == 0 ) then
        call grid_collect(gasabsnxt,gasabsnxt_io,jci1,jci2,ici1,ici2,1,kz,1,4)
        call grid_collect(gasabstot,gasabstot_io, &
                          jci1,jci2,ici1,ici2,1,kzp1,1,kzp1)
        call grid_collect(gasemstot,gasemstot_io,jci1,jci2,ici1,ici2,1,kzp1)
      end if

      call subgrid_collect(tlef,tlef_io,jci1,jci2,ici1,ici2)
      call subgrid_collect(ssw,ssw_io,jci1,jci2,ici1,ici2)
      call subgrid_collect(rsw,rsw_io,jci1,jci2,ici1,ici2)
      call subgrid_collect(tgrd,tgrd_io,jci1,jci2,ici1,ici2)
      call subgrid_collect(tgbrd,tgbrd_io,jci1,jci2,ici1,ici2)
      call subgrid_collect(sncv,sncv_io,jci1,jci2,ici1,ici2)
      call subgrid_collect(gwet,gwet_io,jci1,jci2,ici1,ici2)
      call subgrid_collect(snag,snag_io,jci1,jci2,ici1,ici2)
      call subgrid_collect(sfice,sfice_io,jci1,jci2,ici1,ici2)
      call subgrid_collect(ldew,ldew_io,jci1,jci2,ici1,ici2)
      call subgrid_collect(taf,taf_io,jci1,jci2,ici1,ici2)
      call subgrid_collect(tsw,tsw_io,jci1,jci2,ici1,ici2)
      call subgrid_collect(emiss,emiss_io,jci1,jci2,ici1,ici2)
      call subgrid_collect(ldmsk1,ldmsk1_io,jci1,jci2,ici1,ici2)

      call grid_collect(solis,solis_io,jci1,jci2,ici1,ici2)
      call grid_collect(solvd,solvd_io,jci1,jci2,ici1,ici2)
      call grid_collect(solvs,solvs_io,jci1,jci2,ici1,ici2)
      call grid_collect(sabveg,sabveg_io,jci1,jci2,ici1,ici2)
      call grid_collect(flw,flw_io,jci1,jci2,ici1,ici2)
      call grid_collect(flwd,flwd_io,jci1,jci2,ici1,ici2)
      call grid_collect(fsw,fsw_io,jci1,jci2,ici1,ici2)
      call grid_collect(sinc,sinc_io,jci1,jci2,ici1,ici2)
      call grid_collect(ldmsk,ldmsk_io,jci1,jci2,ici1,ici2)

#ifndef CLM
      if ( lakemod == 1 ) then
        call subgrid_collect(eta,eta_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(hi,hi_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(aveice,aveice_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(hsnow,hsnow_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(tlak,tlak_io,jci1,jci2,ici1,ici2,1,ndpmax)
      end if
#else
      call grid_collect(sols2d,sols2d_io,jci1,jci2,ici1,ici2)
      call grid_collect(soll2d,soll2d_io,jci1,jci2,ici1,ici2)
      call grid_collect(solsd2d,solsd2d_io,jci1,jci2,ici1,ici2)
      call grid_collect(solld2d,solld2d_io,jci1,jci2,ici1,ici2)
      call grid_collect(aldirs2d,aldirs2d_io,jci1,jci2,ici1,ici2)
      call grid_collect(aldirl2d,aldirl2d_io,jci1,jci2,ici1,ici2)
      call grid_collect(aldifs2d,aldifs2d_io,jci1,jci2,ici1,ici2)
      call grid_collect(aldifl2d,aldifl2d_io,jci1,jci2,ici1,ici2)
      call grid_collect(lndcat2d,lndcat2d_io,jci1,jci2,ici1,ici2)
#endif
      if ( idcsst == 1 ) then
        call grid_collect(dtskin,dtskin_io,jci1,jci2,ici1,ici2)
        call grid_collect(deltas,deltas_io,jci1,jci2,ici1,ici2)
        call grid_collect(tdeltas,tdeltas_io,jci1,jci2,ici1,ici2)
      end if

      call grid_collect(dstor,dstor_io,jde1,jde2,ide1,ide2,1,nsplit)
      call grid_collect(hstor,hstor_io,jde1,jde2,ide1,ide2,1,nsplit)

      if ( ichem == 1 ) then
        call grid_collect(chia,chia_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        call grid_collect(chib,chib_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        call grid_collect(remlsc,remlsc_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        call grid_collect(remcvc,remcvc_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        call grid_collect(remdrd,remdrd_io,jce1,jce2,ice1,ice2,1,ntr)

        call grid_collect(ssw2da,ssw2da_io,jci1,jci2,ici1,ici2)
        call grid_collect(sdeltk2d,sdeltk2d_io,jci1,jci2,ici1,ici2)
        call grid_collect(sdelqk2d,sdelqk2d_io,jci1,jci2,ici1,ici2)
        call grid_collect(sfracv2d,sfracv2d_io,jci1,jci2,ici1,ici2)
        call grid_collect(sfracb2d,sfracb2d_io,jci1,jci2,ici1,ici2)
        call grid_collect(sfracs2d,sfracs2d_io,jci1,jci2,ici1,ici2)
        call grid_collect(svegfrac2d,svegfrac2d_io,jci1,jci2,ici1,ici2)
      end if

      if ( ldosav ) then
        call write_savefile(idatex, .false.)
      else
        call write_savefile(idatex, .true.)
      end if
    end if
  end if

  if ( lfdomonth(idatex) .and. lmidnight(idatex) ) then
    if ( .not. lstartup .and. idatex /= idate2 ) then
      call newoutfiles(idatex)
    end if
  end if
  if ( myid == iocpu ) then
    if ( lfdomonth(idatex) .and. lmidnight(idatex) ) then
      if ( .not. lstartup .and. idatex /= idate2 ) then
        call mkfile
        call checktime(myid)
      end if
    end if
  end if
#ifdef DEBUG
  call time_end(subroutine_name,idindx) 
#endif
  end subroutine output
!
  subroutine mkfile
    implicit none
    if ( myid /= iocpu ) return
    print * , ' '
    print * , '******* OPENING NEW OUTPUT FILES : ' , tochar(idatex)
    print * , ' '
    if ( ifatm ) then
      call prepare_common_out(idatex,'ATM')
    end if
    if ( ifsrf ) then
      call prepare_common_out(idatex,'SRF')
      if (lakemod == 1 .and. iflak) then
        call prepare_common_out(idatex,'LAK')
      end if
      call prepare_common_out(idatex,'STS')
    end if
    if ( nsg > 1 .and. ifsub ) then
      call prepare_common_out(idatex,'SUB')
    end if
    if ( ifrad ) then
      call prepare_common_out(idatex,'RAD')
    end if
    if ( ichem == 1 ) then
      if ( ifchem ) then
        call prepare_chem_out(idatex,ifrest)
      end if
    end if
  end subroutine mkfile
!
  subroutine outatm
    implicit none
    call writerec_atm(atm1_io%u,atm1_io%v,omega_io,atm1_io%t,atm1_io%qx, &
                      atm1_io%tke,tcmstate_io%kth,tcmstate_io%kzm, &
                      sfs_io%psa,sfs_io%rainc,sfs_io%rainnc,tgbrd_io, &
                      tsw_io,ldmsk1_io,idatex)
    print *, 'ATM variables written at ' , tochar(idatex)
  end subroutine outatm
!
  subroutine outsrf
    implicit none
    call writerec_srf(fbat_io,ldmsk_io,idatex)
    print *, 'SRF variables written at ' , tochar(idatex)

    if ( ifsts .and. mod(ktau+kstsoff,ksts) == 0 .and. ktau > kstsoff+2 ) then
      call writerec_sts(fbat_io,idatex)
      print *, 'STS variables written at ' , tochar(idatex)
    end if
#ifndef CLM
    if ( lakemod == 1 .and. iflak .and. iolak == klak ) then
      call writerec_lak(fbat_io,evl_io,aveice_io, &
                        hsnow_io,tlak_io,iveg1_io,iveg_io,idatex)
      print *, 'LAK variables written at ' , tochar(idatex)
    end if
#endif
  end subroutine outsrf
!
  subroutine outsub
    implicit none
#ifndef CLM
    call writerec_sub(fsub_io,iveg1_io,tlak_io,idatex)
#endif
    print *, 'SUB variables written at ' , tochar(idatex)
  end subroutine outsub
!
  subroutine outrad
    implicit none
    call writerec_rad(4,12,frad3d_io,frad2d_io,sfs_io%psa,idatex)
    print * , 'RAD variables written at ' , tochar(idatex)
  end subroutine outrad
!
end module mod_output
