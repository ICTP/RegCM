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
  use mod_bdycod
  use mod_precip
  use mod_split
  use mod_savefile
  use mod_mppio
#ifdef CLM
  use mod_clm
#endif

  private

  integer :: iolak
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
  use mpi
  implicit none
!
  integer :: j2c , i2c , j2o
  logical :: ldoatm , ldosrf , ldorad , ldoche , ldosav , ldotmp
  logical :: lstartup
  character (len=64) :: subroutine_name='output'
  integer :: idindx=0
!
!
  call time_begin(subroutine_name,idindx)
!
#ifdef BAND
  j2c = jx
  j2o  = jx
#else
  j2c = jxm1
  j2o  = jxm2
#endif
  i2c = iym1
!
!----------------------------------------------------------------------
!
  lstartup = .false.
  if ( ktau == 0 .or. doing_restart ) then
    if ( myid == 0 ) then
      call mkfile
    end if
    lstartup = .true.
    if ( doing_restart ) then
      doing_restart = .false.
      call time_end(subroutine_name,idindx) 
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
!     gather  ua,va,ta,qva,qca,rainc,rainnc,tgbrd,tsw,olcd2d
!=======================================================================
!
      call deco1_gather(atm1%u,atm1_io%u,jdot1,jdot2,idot1,idot2,1,kz)
      call deco1_gather(atm1%v,atm1_io%v,jdot1,jdot2,idot1,idot2,1,kz)
      call deco1_gather(atm1%t,atm1_io%t,jcross1,jcross2,icross1,icross2,1,kz)
      call deco1_gather(atm1%qv,atm1_io%qv,jcross1,jcross2,icross1,icross2,1,kz)
      call deco1_gather(atm1%qc,atm1_io%qc,jcross1,jcross2,icross1,icross2,1,kz)
      call deco1_gather(omega,omega_io,jcross1,jcross2,icross1,icross2,1,kz)
      call deco1_gather(sfs%psa,sfs_io%psa,jcross1,jcross2,icross1,icross2)
      call deco1_gather(sfs%rainc,sfs_io%rainc,jcross1,jcross2,icross1,icross2)
      call deco1_gather(sfs%rainnc,sfs_io%rainnc, &
                        jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_gather(tgbrd,tgbrd_io,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_gather(tsw,tsw_io,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_gather(ldmsk1,ldmsk1_io,jcross1,jcross2,icross1,icross2)
!
!=======================================================================
!     gather UW Scheme variables
!=======================================================================
!
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call deco1_gather(atm1%tke,atm1_io%tke, &
                          jcross1,jcross2,icross1,icross2,1,kz)
        call deco1_gather(uwstateb%kth,tcmstate_io%kth, &
                          jcross1,jcross2,icross1,icross2,1,kz)
        call deco1_gather(uwstateb%kzm,tcmstate_io%kzm, &
                          jcross1,jcross2,icross1,icross2,1,kz)
      end if
!
      if ( myid == 0 ) then
        call outatm
      end if
      sfs%rainc(:,:)  = d_zero
      sfs%rainnc(:,:) = d_zero
    end if
  end if
 
!     Call surface output
 
  if ( ifsrf ) then
    if ( ldosrf ) then
#ifndef CLM
      if ( lakemod == 1 .and. iflak .and. iolak == klak ) then
        call subgrid_deco1_gather(eta,eta_io,jcross1,jcross2,icross1,icross2)
        call subgrid_deco1_gather(hi,hi_io,jcross1,jcross2,icross1,icross2)
        call subgrid_deco1_gather(aveice,aveice_io, &
                                  jcross1,jcross2,icross1,icross2)
        call subgrid_deco1_gather(hsnow,hsnow_io, &
                                  jcross1,jcross2,icross1,icross2)
        call subgrid_deco1_gather(tlak,tlak_io, &
                                  jcross1,jcross2,icross1,icross2,1,ndpmax)
      end if
#endif
      if ( iseaice == 1 .or. lakemod == 1 ) then
        call deco1_gather(ldmsk,ldmsk_io,jcross1,jcross2,icross1,icross2)
      end if
      call deco1_gather(fbat,fbat_io,jout1,jout2,iout1,iout2,1,numbat)
      if ( myid == 0 ) then
        call outsrf
      end if
      sunt_o(:,:) = 0.0
      if ( ifsts .and. mod(ktau+kstsoff,ksts) == 0 .and. &
           ktau > kstsoff+2 ) then
        tgmx_o(:,:) = -1.E30
        tgmn_o(:,:) =  1.E30
        t2mx_o(:,:) = -1.E30
        t2mn_o(:,:) =  1.E30
        tavg_o(:,:) = 0.0
        w10x_o(:,:) = -1.E30
        pcpa_o(:,:) = 0.0
        pcpx_o(:,:) = -1.E30
        sund_o(:,:) = 0.0
        psmn_o(:,:) =  1.E30
      end if
      if ( lakemod == 1 .and. iflak .and. iolak == klak) then
        iolak = 1
      else
        iolak = iolak + 1
      end if

      if ( ifsub .and. nsg > 1 ) then

        call subgrid_deco1_gather(fsub,fsub_io,jout1,jout2,iout1,iout2,1,numsub)
        if ( myid == 0 ) then
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
      call deco1_gather(frad2d,frad2d_io,jout1,jout2,iout1,iout2,1,nrad2d)
      call deco1_gather(frad3d,frad3d_io,jout1,jout2,iout1,iout2,1,kz,1,nrad3d)
      call deco1_gather(sfs%psa,sfs_io%psa,jout1,jout2,iout1,iout2)
      if ( myid == 0 ) then
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
      call deco1_gather(atm1%u,atm1_io%u,jdot1,jdot2,idot1,idot2,1,kz)
      call deco1_gather(atm1%v,atm1_io%v,jdot1,jdot2,idot1,idot2,1,kz)
      call deco1_gather(atm1%t,atm1_io%t,jcross1,jcross2,icross1,icross2,1,kz)
      call deco1_gather(atm1%qv,atm1_io%qv,jcross1,jcross2,icross1,icross2,1,kz)
      call deco1_gather(atm1%qc,atm1_io%qc,jcross1,jcross2,icross1,icross2,1,kz)

      call deco1_gather(atm2%u,atm2_io%u,jdot1,jdot2,idot1,idot2,1,kz)
      call deco1_gather(atm2%v,atm2_io%v,jdot1,jdot2,idot1,idot2,1,kz)
      call deco1_gather(atm2%t,atm2_io%t,jcross1,jcross2,icross1,icross2,1,kz)
      call deco1_gather(atm2%qv,atm2_io%qv,jcross1,jcross2,icross1,icross2,1,kz)
      call deco1_gather(atm2%qc,atm2_io%qc,jcross1,jcross2,icross1,icross2,1,kz)

      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call deco1_gather(atm1%tke,atm1_io%tke, &
                          jcross1,jcross2,icross1,icross2,1,kz)
        call deco1_gather(atm2%tke,atm2_io%tke, &
                          jcross1,jcross2,icross1,icross2,1,kz)
        call deco1_gather(kpbl,kpbl_io,jcross1,jcross2,icross1,icross2)
      end if

      call deco1_gather(sfs%psa,sfs_io%psa,jcross1,jcross2,icross1,icross2)
      call deco1_gather(sfs%psb,sfs_io%psb,jcross1,jcross2,icross1,icross2)
      call deco1_gather(sfs%tga,sfs_io%tga,jcross1,jcross2,icross1,icross2)
      call deco1_gather(sfs%tgb,sfs_io%tgb,jcross1,jcross2,icross1,icross2)
      call deco1_gather(sfs%hfx,sfs_io%hfx,jcross1,jcross2,icross1,icross2)
      call deco1_gather(sfs%qfx,sfs_io%qfx,jcross1,jcross2,icross1,icross2)
      call deco1_gather(sfs%rainc,sfs_io%rainc,jcross1,jcross2,icross1,icross2)
      call deco1_gather(sfs%rainnc,sfs_io%rainnc, &
                        jcross1,jcross2,icross1,icross2)
      call deco1_gather(sfs%tgbb,sfs_io%tgbb,jcross1,jcross2,icross1,icross2)
      call deco1_gather(sfs%uvdrag,sfs_io%uvdrag, &
                        jcross1,jcross2,icross1,icross2)

      if ( ipptls == 1 ) then
        call deco1_gather(fcc,fcc_io,jcross1,jcross2,icross1,icross2,1,kz)
      end if
      call deco1_gather(heatrt,heatrt_io,jcross1,jcross2,icross1,icross2,1,kz)
      call deco1_gather(o3prof,o3prof_io,jcross1,jcross2,icross1,icross2,1,kzp1)

      if ( iocnflx == 2 ) then
        call deco1_gather(zpbl,zpbl_io,jcross1,jcross2,icross1,icross2)
      end if
      if ( icup == 1 ) then
        call deco1_gather(rsheat,rsheat_io,jcross1,jcross2,icross1,icross2,1,kz)
        call deco1_gather(rswat,rswat_io,jcross1,jcross2,icross1,icross2,1,kz)
      end if
      if ( icup == 3 ) then
        call deco1_gather(tbase,tbase_io,jcross1,jcross2,icross1,icross2,1,kz)
        call deco1_gather(cldefi,cldefi_io,jcross1,jcross2,icross1,icross2)
      end if
      if ( icup==4 .or. icup==99 .or. icup==98 ) then
        call deco1_gather(cbmf2d,cbmf2d_io,jcross1,jcross2,icross1,icross2)
      end if

      if ( irrtm == 0 ) then
        call deco1_gather(gasabsnxt,gasabsnxt_io, &
                          jcross1,jcross2,icross1,icross2,1,kz,1,4)
        call deco1_gather(gasabstot,gasabstot_io, &
                          jcross1,jcross2,icross1,icross2,1,kzp1,1,kzp1)
        call deco1_gather(gasemstot,gasemstot_io, &
                          jcross1,jcross2,icross1,icross2,1,kzp1)
      end if

      call subgrid_deco1_gather(tlef,tlef_io,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_gather(ssw,ssw_io,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_gather(rsw,rsw_io,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_gather(tgrd,tgrd_io,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_gather(tgbrd,tgbrd_io,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_gather(sncv,sncv_io,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_gather(gwet,gwet_io,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_gather(snag,snag_io,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_gather(sfice,sfice_io,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_gather(ldew,ldew_io,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_gather(taf,taf_io,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_gather(tsw,tsw_io,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_gather(emiss,emiss_io,jcross1,jcross2,icross1,icross2)
      call subgrid_deco1_gather(ldmsk1,ldmsk1_io,jcross1,jcross2,icross1,icross2)

      call deco1_gather(solis,solis_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(solvd,solvd_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(solvs,solvs_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(sabveg,sabveg_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(flw,flw_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(flwd,flwd_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(fsw,fsw_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(sinc,sinc_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(ldmsk,ldmsk_io,jcross1,jcross2,icross1,icross2)
#ifndef CLM
      if ( lakemod == 1 ) then
        call subgrid_deco1_gather(eta,eta_io,jcross1,jcross2,icross1,icross2)
        call subgrid_deco1_gather(hi,hi_io,jcross1,jcross2,icross1,icross2)
        call subgrid_deco1_gather(aveice,aveice_io, &
                                  jcross1,jcross2,icross1,icross2)
        call subgrid_deco1_gather(hsnow,hsnow_io, &
                                  jcross1,jcross2,icross1,icross2)
        call subgrid_deco1_gather(tlak,tlak_io, &
                                  jcross1,jcross2,icross1,icross2,1,ndpmax)
      end if
#else
      call deco1_gather(sols2d,sols2d_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(soll2d,soll2d_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(solsd2d,solsd2d_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(solld2d,solld2d_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(aldirs2d,aldirs2d_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(aldirl2d,aldirl2d_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(aldifs2d,aldifs2d_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(aldifl2d,aldifl2d_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(lndcat2d,lndcat2d_io,jcross1,jcross2,icross1,icross2)
#endif
 
      if ( ichem == 1 ) then
        call deco1_gather(chia,chia_io, &
                          jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
        call deco1_gather(chib,chib_io, &
                          jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
        call deco1_gather(remlsc,remlsc_io, &
                          jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
        call deco1_gather(remcvc,remcvc_io, &
                          jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
        call deco1_gather(remdrd,remdrd_io, &
                          jcross1,jcross2,icross1,icross2,1,ntr)
        call deco1_gather(ssw2da,ssw2da_io, &
                          jcross1,jcross2,icross1,icross2)
        call deco1_gather(sdeltk2d,sdeltk2d_io, &
                          jcross1,jcross2,icross1,icross2)
        call deco1_gather(sdelqk2d,sdelqk2d_io, &
                          jcross1,jcross2,icross1,icross2)
        call deco1_gather(sfracv2d,sfracv2d_io, &
                          jcross1,jcross2,icross1,icross2)
        call deco1_gather(sfracb2d,sfracb2d_io, &
                          jcross1,jcross2,icross1,icross2)
        call deco1_gather(sfracs2d,sfracs2d_io, &
                          jcross1,jcross2,icross1,icross2)
        call deco1_gather(svegfrac2d,svegfrac2d_io, &
                          jcross1,jcross2,icross1,icross2)
      end if

      call deco1_gather(dstor,dstor_io,jdot1,jdot2,idot1,idot2,1,nsplit)
      call deco1_gather(hstor,hstor_io,jdot1,jdot2,idot1,idot2,1,nsplit)

      if ( ldosav ) then
        call write_savefile(idatex, .false.)
      else
        call write_savefile(idatex, .true.)
      end if
    end if
  end if

  if ( myid == 0 ) then
    if ( lfdomonth(idatex) .and. lmidnight(idatex) ) then
      if ( .not. lstartup .and. idatex /= idate2 ) then
        call mkfile
      end if
    end if
  end if
!
  call time_end(subroutine_name,idindx) 

  end subroutine output
!
  subroutine mkfile
 
  implicit none
!
  if (myid /= 0) return

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
!      call prepare_common_out(idatex,'CHE')
!      if (iaerosol == 1)
!       call prepare_opt_out(idatex)
      call prepare_chem_out(idatex,ifrest)
    end if
  end if

  end subroutine mkfile
!
  subroutine outatm

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine writes the model output to tape or disk for use c
!     in dataflow analyses.                                           c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  implicit none
  integer :: jjx , iiy
#ifdef BAND
  jjx = jx
  iiy = iym1
#else
  jjx = jxm1
  iiy = iym1
#endif

  call writerec_atm(jx,iy,jjx,iiy,kz,nnsg,atm1_io%u,atm1_io%v,     &
          omega_io,atm1_io%t,atm1_io%qv,atm1_io%qc,atm1_io%tke,    &
          tcmstate_io%kth,tcmstate_io%kzm,sfs_io%psa,sfs_io%rainc, &
          sfs_io%rainnc,tgbrd_io,tsw_io,ldmsk1_io,idatex)
 
  print *, 'ATM variables written at ' , tochar(idatex)
 
  end subroutine outatm
!
  subroutine outsrf

  implicit none
!
  integer :: i , j
!
#ifdef BAND
  i = iym3
  j = jx
#else
  i = iym3
  j = jxm3
#endif

  call writerec_srf(j,i,fbat_io,ldmsk_io,idatex)
  print *, 'SRF variables written at ' , tochar(idatex)

  if ( ifsts .and. mod(ktau+kstsoff,ksts) == 0 .and. ktau > kstsoff+2 ) then
    call writerec_sts(j,i,numbat,fbat_io,idatex)
    print *, 'STS variables written at ' , tochar(idatex)
  end if

#ifndef CLM
  if ( lakemod == 1 .and. iflak .and. iolak == klak ) then
    call writerec_lak(j,i,numbat,fbat_io,evl_io,aveice_io, &
                      hsnow_io,tlak_io,idatex)
    print *, 'LAK variables written at ' , tochar(idatex)
  end if
#endif

  end subroutine outsrf
!
  subroutine outsub

  implicit none

  integer :: i , j

#ifdef BAND
  i = iym3sg
  j = jxsg
#else
  i = iym3sg
  j = jxm3sg
#endif

  call writerec_sub(j,i,nsg,numsub,fsub_io,idatex)

  print *, 'SUB variables written at ' , tochar(idatex)

  end subroutine outsub
!
  subroutine outrad
!
  implicit none
!
  integer :: i , j , imax , jmax , jmin
 
!   character (len=64) :: subroutine_name='outrad'
!   integer :: idindx=0
! 
!  call time_begin(subroutine_name,idindx)
#ifdef BAND
  imax = iym2
  jmin = 1
  jmax = jx
#else
  imax = iym2
  jmin = 2
  jmax = jxm2
#endif

  do i = 2 , imax
    do j = jmin , jmax
      radpsa_io(j,i) = real((sfs_io%psa(j,i)+ptop)*d_10)
    end do
  end do

  call writerec_rad(jxm3,iym3,kz,4,12,frad3d_io,frad2d_io,radpsa_io,idatex)

  print * , 'RAD variables written at ' , tochar(idatex)
!   call time_end(subroutine_name,idindx)

  end subroutine outrad
!
  subroutine outche

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine writes the model chem                           c
!                                                                     c
!     iutl : is the output unit number for large-domain variables.    c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  implicit none
!
  integer :: ni , itr , nj , nk

!      character (len=64) :: subroutine_name='outche'
!      integer :: idindx=0
!
!      call time_begin(subroutine_name,idindx)
#ifdef BAND
  ni = iym3
  nj = jx
  nk = kz
  itr = ntr
#else
  ni = iym3
  nj = jxm3
  nk = kz
  itr = ntr
#endif

     call writerec_che(nj, ni, nk, itr, chia_io, &
            aerext_io, aerssa_io, aerasp_io, dtrace_io,  &
            wdlsc_io, wdcvc_io, ddsfc_io, wxsg_io,       &
            wxaq_io, cemtrac_io, aertarf_io, aersrrf_io, &
            aertalwrf_io, aersrlwrf_io, sfs_io%psa, idatex)

  print *, 'CHE variables written at ' , tochar(idatex)

  end subroutine outche
!
end module mod_output
