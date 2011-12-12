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

module mod_mppio

  use mod_runparams
  use mod_atm_interface , only : atmstate , allocate_atmstate
  use mod_atm_interface , only : domain , allocate_domain
  use mod_lm_interface
  use mod_cu_interface
  use mod_che_interface
  use mod_rad_interface
  use mod_pbl_interface , only : ibltyp , tcm_state , allocate_tcm_state
  use mod_memutil
  use mod_mpmessage
!
  integer , pointer , dimension(:,:,:) :: ocld2d_io , veg2d1_io
  integer , pointer , dimension(:,:) :: veg2d_io , ldmsk_io
  integer , pointer , dimension(:,:) :: kpbl_io

  real(8) , pointer , dimension(:,:,:) :: ldew_io , gwet_io , rno_io , &
         snag_io , sncv_io , sfice_io , rsw_io , ssw_io , tsw_io ,     &
         taf_io , text2d_io , tgrd_io , tgbrd_io , tlef_io , emiss_io

  real(8) , pointer , dimension(:,:,:) :: ht1_io , lndcat1_io ,     &
                                          xlat1_io , xlon1_io
!
  real(8) , pointer , dimension(:,:) :: flw_io , flwd_io , fsw_io , &
         sabveg_io , sinc_io , solis_io ,  solvd_io , solvs_io
!
  real(4) , pointer , dimension(:,:,:) :: fbat_io
  real(4) , pointer , dimension(:,:,:,:) :: fsub_io
  real(4) , pointer , dimension(:,:,:) :: fsavsts_io
  real(4) , pointer , dimension(:,:,:) :: fsavsts
  real(4) , pointer , dimension(:,:,:) :: frad2d_io
  real(4) , pointer , dimension(:,:,:,:) :: frad3d_io
  real(4) , pointer , dimension(:,:) :: radpsa_io

  real(8) , pointer , dimension(:,:) :: cbmf2d_io
  real(8) , pointer , dimension(:,:,:) :: fcc_io , rsheat_io ,  &
                                    rswat_io

  real(8) , pointer , dimension(:,:,:,:) :: gasabsnxt_io
  real(8) , pointer , dimension(:,:,:,:) :: gasabstot_io
  real(8) , pointer , dimension(:,:,:) :: gasemstot_io

  real(8) , pointer , dimension(:,:,:) :: heatrt_io
  real(8) , pointer , dimension(:,:,:) :: o3prof_io

  real(8) , pointer , dimension(:,:,:) :: dstor_io , hstor_io

  real(8) , pointer , dimension(:,:,:) :: aerasp_io ,           &
                              aerext_io , aerssa_io
  real(8) , pointer , dimension(:,:) :: aersrrf_io , aertarf_io,&
                            aertalwrf_io , aersrlwrf_io

  real(8) , pointer , dimension(:,:) :: ps0_io , ps1_io , ts0_io ,  &
                           ts1_io
  real(8) , pointer , dimension(:,:,:) :: qb0_io , qb1_io , &
            tb0_io , tb1_io , ub0_io , ub1_io , vb0_io , vb1_io
  real(8) , pointer , dimension(:,:) :: ui1_io , ui2_io , uilx_io , &
                                 uil_io , vi1_io , vi2_io ,         &
                                 vilx_io , vil_io

  real(8) , pointer , dimension(:,:) :: pptc_io , pptnc_io ,    &
                                 prca2d_io , prnca2d_io

  real(8) , pointer , dimension(:,:) :: cldefi_io , hfx_io , &
                                 psa_io , psb_io , qfx_io ,  &
                                 rainc_io , rainnc_io ,    &
                                 tga_io , tgbb_io ,     &
                                 tgb_io , uvdrag_io , &
                                 zpbl_io
  real(8) , pointer , dimension(:,:,:) :: omega_io , tbase_io

  type(atmstate) :: atm1_io , atm2_io
  type(domain) :: mddom_io
  type(tcm_state) :: tcmstate_io

  real(8) , pointer , dimension(:,:,:) :: inisrf0
  real(8) , pointer , dimension(:,:,:) :: inisrf_0

  real(8) , pointer , dimension(:,:) :: var1snd , var1rcv
  integer , pointer , dimension(:,:) :: var2d0
  integer , pointer , dimension(:,:) :: var2d_0
  integer , pointer , dimension(:,:,:) :: var2d1
  integer , pointer , dimension(:,:,:) :: var2d_1
  real(8) , pointer , dimension(:,:) :: swapv
 
  real(8) , pointer , dimension(:,:,:) :: atm0
  real(8) , pointer , dimension(:,:,:) :: atm_0
  real(8) , pointer , dimension(:,:,:) :: uw0
  real(8) , pointer , dimension(:,:,:) :: uw_0
  real(4) , pointer , dimension(:,:,:) :: bat0
  real(4) , pointer , dimension(:,:,:) :: bat_0
  real(4) , pointer , dimension(:,:,:) :: rad0
  real(4) , pointer , dimension(:,:,:) :: rad_0
  real(4) , pointer , dimension(:,:,:,:) :: sub0
  real(4) , pointer , dimension(:,:,:,:) :: sub_0
#ifdef CLM
  real(8) , pointer , dimension(:,:) :: sols2d_io , soll2d_io ,     &
                      solsd2d_io , solld2d_io , aldifl2d_io ,       &
                      aldirs2d_io , aldirl2d_io , aldifs2d_io
  real(8) , pointer , dimension(:,:) :: lndcat2d_io
#endif
!
  real(8) , pointer , dimension(:,:,:) :: sav0
  real(8) , pointer , dimension(:,:,:) :: sav_0
  real(8) , pointer , dimension(:,:,:) :: sav0a
  real(8) , pointer , dimension(:,:,:) :: sav_0a
  real(8) , pointer , dimension(:,:,:) :: sav0b
  real(8) , pointer , dimension(:,:,:) :: sav_0b
  real(8) , pointer , dimension(:,:,:) :: sav0c
  real(8) , pointer , dimension(:,:,:) :: sav_0c
  real(8) , pointer , dimension(:,:,:) :: sav0s
  real(8) , pointer , dimension(:,:,:) :: sav_0s
  real(8) , pointer , dimension(:,:,:) :: sav0d
  real(8) , pointer , dimension(:,:,:) :: sav_0d
  real(8) , pointer , dimension(:,:,:) :: sav1
  real(8) , pointer , dimension(:,:,:) :: sav_1
  real(8) , pointer , dimension(:,:,:) :: sav2
  real(8) , pointer , dimension(:,:,:) :: sav_2
  integer , pointer , dimension(:,:,:) :: sav2a
  integer , pointer , dimension(:,:,:) :: sav_2a
  real(8) , pointer , dimension(:,:,:) :: sav4
  real(8) , pointer , dimension(:,:,:) :: sav_4
  real(8) , pointer , dimension(:,:,:) :: sav4a
  real(8) , pointer , dimension(:,:,:) :: sav_4a
  real(8) , pointer , dimension(:,:,:) :: sav6
  real(8) , pointer , dimension(:,:,:) :: sav_6
#ifdef CLM
  real(8) , pointer , dimension(:,:,:) :: sav_clmout
  real(8) , pointer , dimension(:,:,:) :: sav_clmin
#endif

!---------- DATA init section--------------------------------------------

  contains 
!
!     This routines allocate all the arrays contained in the module
!
  subroutine allocate_mod_mppio(lband)
    implicit none
    logical , intent(in) :: lband
    integer :: mjj , mojj

    if ( lband ) then
      mjj = jx
      mojj = jx
    else
      mjj = jxm1
      mojj = jxm2
    end if

    call getmem2d(var1snd,1,kz,1,8,'mppio:var1snd')
    call getmem2d(var1rcv,1,kz,1,8,'mppio:var1rcv')
    call getmem2d(var2d0,1,iy,1,jxp,'mppio:var2d0')
    call getmem2d(swapv,1,iy,1,jxp,'mppio:swap')
    call getmem3d(var2d1,1,iy,1,nnsg,1,jxp,'mppio:var2d1')
    call getmem3d(inisrf0,1,iy,1,nnsg*4+7,1,jxp,'mppio:inisrf0')
    call getmem3d(atm0,1,iy,1,kz*6+3+nnsg*3,1,jxp,'mppio:atm0')
    call getmem3d(bat0,1,iym2,1,numbat,1,jxp,'mppio:bat0')
    call getmem3d(rad0,1,iym2,1,nrad3d*kz+nrad2d,1,jxp,'mppio:rad0')
    call getmem4d(sub0,1,iym2,1,nnsg,1,numsub,1,jxp,'mppio:sub0')
    call getmem3d(fsavsts,1,iym2,1,numsts,1,jxp,'mod_mppio:fsavsts')
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      call getmem3d(uw0,1,iy,1,kz*3,1,jxp,'mppio:uw0')
    end if

    if (myid == 0) then
      call allocate_domain(mddom_io,.false.)
      call allocate_atmstate(atm1_io,ibltyp,.false.,0,0)
      call allocate_atmstate(atm2_io,ibltyp,.false.,0,0)
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call allocate_tcm_state(tcmstate_io,.false.)
      end if

      call getmem2d(kpbl_io,1,iy,1,jx,'mppio:kpbl_io')

      call getmem2d(var2d_0,1,iy,1,jx,'mppio:var2d_0')
      call getmem3d(var2d_1,1,iy,1,nnsg,1,jx,'mppio:var2d_1')
      call getmem3d(inisrf_0,1,iy,1,nnsg*4+7,1,jx,'mppio:inisrf_0')
      call getmem3d(atm_0,1,iy,1,kz*6+3+nnsg*3,1,jx,'mppio:atm_0')
      call getmem3d(bat_0,1,iym2,1,numbat,1,jx,'mppio:bat_0')
      call getmem3d(rad_0,1,iym2,1,nrad3d*kz+nrad2d,1,jx,'mppio:rad_0')
      call getmem4d(sub_0,1,iym2,1,nnsg,1,numsub,1,jx,'mppio:sub_0')
      call getmem3d(fsavsts_io,1,iym2,1,numsts,1,jx,'mod_mppio:fsavsts_io')
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call getmem3d(uw_0,1,iy,1,kz*3,1,jx,'mppio:uw_0')
      end if

      call getmem3d(ldew_io,1,nnsg,1,iym1,1,mjj,'mppio:ldew_io')
      call getmem3d(gwet_io,1,nnsg,1,iym1,1,mjj,'mppio:gwet_io')
      call getmem3d(rno_io,1,nnsg,1,iym1,1,mjj,'mppio:rno_io')
      call getmem3d(snag_io,1,nnsg,1,iym1,1,mjj,'mppio:snag_io')
      call getmem3d(sncv_io,1,nnsg,1,iym1,1,mjj,'mppio:sncv_io')
      call getmem3d(sfice_io,1,nnsg,1,iym1,1,mjj,'mppio:sfice_io')
      call getmem3d(rsw_io,1,nnsg,1,iym1,1,mjj,'mppio:rsw_io')
      call getmem3d(ssw_io,1,nnsg,1,iym1,1,mjj,'mppio:ssw_io')
      call getmem3d(tsw_io,1,nnsg,1,iym1,1,mjj,'mppio:tsw_io')
      call getmem3d(taf_io,1,nnsg,1,iym1,1,mjj,'mppio:taf_io')
      call getmem3d(text2d_io,1,nnsg,1,iym1,1,mjj,'mppio:text2d_io')
      call getmem3d(tgrd_io,1,nnsg,1,iym1,1,mjj,'mppio:tgrd_io')
      call getmem3d(tgbrd_io,1,nnsg,1,iym1,1,mjj,'mppio:tgbrd_io')
      call getmem3d(tlef_io,1,nnsg,1,iym1,1,mjj,'mppio:tlef_io')
      call getmem3d(emiss_io,1,nnsg,1,iym1,1,mjj,'mppio:emiss_io')
      call getmem3d(veg2d1_io,1,nnsg,1,iym1,1,mjj,'mppio:veg2d1_io')
      call getmem3d(ocld2d_io,1,nnsg,1,iym1,1,mjj,'mppio:ocld2d_io')
      call getmem2d(veg2d_io,1,iym1,1,mjj,'mppio:veg2d_io')
      call getmem2d(ldmsk_io,1,iym1,1,mjj,'mppio:ldmsk_io')
      call getmem3d(ht1_io,1,nnsg,1,iy,1,jx,'mppio:ht1_io')
      call getmem3d(lndcat1_io,1,nnsg,1,iy,1,jx,'mppio:lndcat1_io')
      call getmem3d(xlat1_io,1,nnsg,1,iy,1,jx,'mppio:xlat1_io')
      call getmem3d(xlon1_io,1,nnsg,1,iy,1,jx,'mppio:xlon1_io')
      call getmem2d(flw_io,1,iym1,1,mjj,'mppio:flw_io')
      call getmem2d(flwd_io,1,iym1,1,mjj,'mppio:flwd_io')
      call getmem2d(fsw_io,1,iym1,1,mjj,'mppio:fsw_io')
      call getmem2d(sabveg_io,1,iym1,1,mjj,'mppio:sabveg_io')
      call getmem2d(sinc_io,1,iym1,1,mjj,'mppio:sinc_io')
      call getmem2d(solis_io,1,iym1,1,mjj,'mppio:solis_io')
      call getmem2d(solvd_io,1,iym1,1,mjj,'mppio:solvd_io')
      call getmem2d(solvs_io,1,iym1,1,mjj,'mppio:solvs_io')

      call getmem3d(fbat_io,1,mojj,1,iym2,1,numbat,'mppio:fbat_io')
      call getmem4d(fsub_io,1,nnsg,1,mojj,1,iym2,1,numsub,'mppio:fsub_io')
      call getmem3d(frad2d_io,1,mojj,1,iym2,1,nrad2d,'mppio:frad2d_io')
      call getmem4d(frad3d_io,1,mojj,1,iym2,1,kz,1,nrad3d,'mppio:frad3d_io')
      call getmem2d(radpsa_io,1,mojj,1,iym2,'mppio:radpsa_io')

      call getmem2d(cbmf2d_io,1,iy,1,jx,'mppio:cbmf2d_io')
      call getmem3d(fcc_io,1,iy,1,kz,1,jx,'mppio:fcc_io')
      call getmem3d(rsheat_io,1,iy,1,kz,1,jx,'mppio:rsheat_io')
      call getmem3d(rswat_io,1,iy,1,kz,1,jx,'mppio:rswat_io')
      call getmem3d(dstor_io,1,iy,1,jx,1,nsplit,'mppio:dstor_io')
      call getmem3d(hstor_io,1,iy,1,jx,1,nsplit,'mppio:hstor_io')

      call getmem4d(gasabsnxt_io,1,iym1,1,kz,1,4,1,mjj,'mppio:gasabsnxt_io')
      call getmem4d(gasabstot_io,1,iym1,1,kzp1,1,kzp1,1,mjj, &
                    'mppio:gasabstot_io')
      call getmem3d(gasemstot_io,1,iym1,1,kzp1,1,mjj,'mppio:gasemstot_io')
      call getmem3d(heatrt_io,1,iym1,1,kz,1,mjj,'mppio:heatrt_io')
      call getmem3d(o3prof_io,1,iym1,1,kzp1,1,mjj,'mppio:o3prof_io')
      call getmem3d(aerasp_io,1,iym1,1,kz,1,mjj,'mppio:aerasp_io')
      call getmem3d(aerext_io,1,iym1,1,kz,1,mjj,'mppio:aerext_io')
      call getmem3d(aerssa_io,1,iym1,1,kz,1,mjj,'mppio:aerssa_io')
      call getmem2d(aersrrf_io,1,iym1,1,mjj,'mppio:aersrrf_io')
      call getmem2d(aertarf_io,1,iym1,1,mjj,'mppio:aertarf_io')
      call getmem2d(aertalwrf_io,1,iym1,1,mjj,'mppio:aertalwrf_io')
      call getmem2d(aersrlwrf_io,1,iym1,1,mjj,'mppio:aersrlwrf_io')

      call getmem2d(ps0_io,1,iy,1,jx,'mppio:ps0_io')
      call getmem2d(ps1_io,1,iy,1,jx,'mppio:ps1_io')
      call getmem2d(ts0_io,1,iy,1,jx,'mppio:ts0_io')
      call getmem2d(ts1_io,1,iy,1,jx,'mppio:ts1_io')
      call getmem3d(qb0_io,1,iy,1,kz,1,jx,'mppio:qb0_io')
      call getmem3d(qb1_io,1,iy,1,kz,1,jx,'mppio:qb1_io')
      call getmem3d(tb0_io,1,iy,1,kz,1,jx,'mppio:tb0_io')
      call getmem3d(tb1_io,1,iy,1,kz,1,jx,'mppio:tb1_io')
      call getmem3d(ub0_io,1,iy,1,kz,1,jx,'mppio:ub0_io')
      call getmem3d(ub1_io,1,iy,1,kz,1,jx,'mppio:ub1_io')
      call getmem3d(vb0_io,1,iy,1,kz,1,jx,'mppio:vb0_io')
      call getmem3d(vb1_io,1,iy,1,kz,1,jx,'mppio:vb1_io')
      call getmem2d(ui1_io,1,kz,1,jx,'mppio:ui1_io')
      call getmem2d(ui2_io,1,kz,1,jx,'mppio:ui2_io')
      call getmem2d(uilx_io,1,kz,1,jx,'mppio:uilx_io')
      call getmem2d(uil_io,1,kz,1,jx,'mppio:uil_io')
      call getmem2d(vi1_io,1,kz,1,jx,'mppio:vi1_io')
      call getmem2d(vi2_io,1,kz,1,jx,'mppio:vi2_io')
      call getmem2d(vilx_io,1,kz,1,jx,'mppio:vilx_io')
      call getmem2d(vil_io,1,kz,1,jx,'mppio:vil_io')
      call getmem2d(pptc_io,1,iym1,1,mjj,'mppio:pptc_io')
      call getmem2d(pptnc_io,1,iym1,1,mjj,'mppio:pptnc_io')
      call getmem2d(prca2d_io,1,iym1,1,mjj,'mppio:prca2d_io')
      call getmem2d(prnca2d_io,1,iym1,1,mjj,'mppio:prnca2d_io')
      call getmem2d(hfx_io,1,iy,1,jx,'mppio:hfx_io')
      call getmem2d(psa_io,1,iy,1,jx,'mppio:psa_io')
      call getmem2d(psb_io,1,iy,1,jx,'mppio:psb_io')
      call getmem2d(qfx_io,1,iy,1,jx,'mppio:qfx_io')
      call getmem2d(rainc_io,1,iy,1,jx,'mppio:rainc_io')
      call getmem2d(rainnc_io,1,iy,1,jx,'mppio:rainnc_io')
      call getmem2d(tga_io,1,iy,1,jx,'mppio:tga_io')
      call getmem2d(tgbb_io,1,iy,1,jx,'mppio:tgbb_io')
      call getmem2d(tgb_io,1,iy,1,jx,'mppio:tgb_io')
      call getmem2d(uvdrag_io,1,iy,1,jx,'mppio:uvdrag_io')
      call getmem2d(zpbl_io,1,iy,1,jx,'mppio:zpbl_io')
      call getmem3d(omega_io,1,iy,1,kz,1,jx,'mppio:omega_io')
      if (icup == 3) then
        call getmem2d(cldefi_io,1,iy,1,jx,'mppio:cldefi_io')
        call getmem3d(tbase_io,1,iy,1,kz,1,jx,'mppio:tbase_io')
      end if
#ifdef CLM
      call getmem2d(sols2d_io,1,iym1,1,mjj,'mppio:sols2d_io')
      call getmem2d(soll2d_io,1,iym1,1,mjj,'mppio:soll2d_io')
      call getmem2d(solsd2d_io,1,iym1,1,mjj,'mppio:solsd2d_io')
      call getmem2d(solld2d_io,1,iym1,1,mjj,'mppio:solld2d_io')
      call getmem2d(aldifl2d_io,1,iym1,1,mjj,'mppio:aldifl2d_io')
      call getmem2d(aldifs2d_io,1,iym1,1,mjj,'mppio:aldifs2d_io')
      call getmem2d(aldirl2d_io,1,iym1,1,mjj,'mppio:aldirl2d_io')
      call getmem2d(aldirs2d_io,1,iym1,1,mjj,'mppio:aldirs2d_io')
      call getmem2d(lndcat2d_io,1,iy,1,jx,'mppio:lndcat2d_io')
#endif
    endif
    if (myid == 0) then
      call getmem3d(sav_0,1,iy,1,kz*4+2,1,jx,'mppio:sav_0')
      call getmem3d(sav_0a,1,iy,1,kzp1+4,1,jx,'mppio:sav_0a')
      call getmem3d(sav_0b,1,iy,1,kzp1,1,jx,'mppio:sav_0b')
      call getmem3d(sav_0c,1,iy,1,kz*2,1,jx,'mppio:sav_0c')
      call getmem3d(sav_0s,1,iy,1,kz,1,jx,'mppio:sav_0s')
      call getmem3d(sav_0d,1,iy,1,nsplit*2,1,jx,'mppio:sav_0d')
      call getmem3d(sav_1,1,iym1,1,kz*4+(kzp1*kzp2),1,jx,'mppio:sav_1')
      call getmem3d(sav_2,1,iym1,1,nnsg*5+4,1,jx,'mppio:sav_2')
      call getmem3d(sav_2a,1,iym1,1,nnsg*2+2,1,jx,'mppio:sav_2a')
      if ( ichem == 1 ) then
        call getmem3d(sav_4,1,iy,1,ntr*(kz*4+1),1,jx,'mppio:sav_4')
        call getmem3d(sav_4a,1,iym1,1,7,1,jx,'mppio:sav_4a')
      end if
      call getmem3d(sav_6,1,kz,1,8,1,jx,'mppio:sav_6')
#ifdef CLM
      call getmem3d(sav_clmout,1,iym1,1,8,1,jx,'mppio:sav_clmout')
#endif
    end if

    call getmem3d(sav0,1,iy,1,kz*4+2,1,jxp,'mppio:sav0')
    call getmem3d(sav0a,1,iy,1,kzp1+4,1,jxp,'mppio:sav0a')
    call getmem3d(sav0b,1,iy,1,kzp1,1,jxp,'mppio:sav0b')
    call getmem3d(sav0c,1,iy,1,kz*2,1,jxp,'mppio:sav0c')
    call getmem3d(sav0s,1,iy,1,kz,1,jxp,'mppio:sav0s')
    call getmem3d(sav0d,1,iy,1,nsplit*2,1,jxp,'mppio:sav0d')
    call getmem3d(sav1,1,iym1,1,kz*4+(kzp1*kzp2),1,jxp,'mppio:sav1')
    call getmem3d(sav2,1,iym1,1,nnsg*5+4,1,jxp,'mppio:sav2')
    call getmem3d(sav2a,1,iym1,1,nnsg*2+2,1,jxp,'mppio:sav2a')
    if ( ichem == 1 ) then
      call getmem3d(sav4,1,iy,1,ntr*(kz*4+1),1,jxp,'mppio:sav4')
      call getmem3d(sav4a,1,iym1,1,7,1,jxp,'mppio:sav4a')
    end if
    call getmem3d(sav6,1,kz,1,8,1,jxp,'mppio:sav6')
#ifdef CLM
    call getmem3d(sav_clmin,1,iym1,1,8,1,jxp,'mppio:sav_clmin')
#endif

  end subroutine allocate_mod_mppio
!
  subroutine free_mpp_initspace
    implicit none
    call relmem3d(inisrf0)
    if (ichem == 1) then
      call relmem4d(src0)
      call relmem3d(src1)
    end if
    if (myid == 0) then
      call relmem3d(inisrf_0)
      if (ichem == 1) then
        call relmem4d(src_0)
        call relmem3d(src_1)
      end if
    end if
  end subroutine free_mpp_initspace
!
end module mod_mppio
