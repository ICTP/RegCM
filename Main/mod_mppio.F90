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
  use mod_message
!
  integer , pointer , dimension(:,:,:) :: ocld2d_io , veg2d1_io
  integer , pointer , dimension(:,:) :: veg2d_io , ldmsk_io
  integer , pointer , dimension(:,:) :: kpbl_io

  real(8) , pointer , dimension(:,:,:) :: col2d_io , dew2d_io ,     &
         evpa2d_io , gwet2d_io , ircp2d_io , rno2d_io , &
         rnos2d_io , sag2d_io , scv2d_io , sena2d_io , sice2d_io ,  &
         srw2d_io , ssw2d_io , swt2d_io , taf2d_io , text2d_io ,    &
         tg2d_io , tgb2d_io , tlef2d_io , emiss2d_io

  real(8) , pointer , dimension(:,:,:) :: ht1_io , lndcat1_io ,     &
                                          xlat1_io , xlon1_io
!
  real(8) , pointer , dimension(:,:) :: flw2d_io , flwd2d_io ,      &
                                   fsw2d_io , sabv2d_io ,           &
                                   sinc2d_io , sol2d_io ,           &
                                   solvd2d_io , solvs2d_io
!
  real(4) , pointer , dimension(:,:,:) :: fbat_io
  real(4) , pointer , dimension(:,:,:,:) :: fsub_io
  real(4) , pointer , dimension(:,:,:) :: frad2d_io
  real(4) , pointer , dimension(:,:,:,:) :: frad3d_io
  real(4) , pointer , dimension(:,:) :: radpsa_io

  real(8) , pointer , dimension(:,:) :: cbmf2d_io
  real(8) , pointer , dimension(:,:,:) :: fcc_io , rsheat_io ,  &
                                    rswat_io

  real(8) , pointer , dimension(:,:,:,:) :: absnxt_io
  real(8) , pointer , dimension(:,:,:,:) :: abstot_io
  real(8) , pointer , dimension(:,:,:) :: emstot_io

  real(8) , pointer , dimension(:,:,:) :: heatrt_io
  real(8) , pointer , dimension(:,:,:) :: o3prof_io

  real(8) , pointer , dimension(:,:,:) :: dstor_io , hstor_io

  real(8) , pointer , dimension(:,:,:) :: aerasp_io ,           &
                              aerext_io , aerssa_io
  real(8) , pointer , dimension(:,:) :: aersrrf_io , aertarf_io,&
                            aertalwrf_io , aersrlwrf_io

  real(8) , pointer , dimension(:,:) :: ps0_io , ps1_io , ts0_io ,  &
                           ts1_io
  real(8) , pointer , dimension(:,:,:) :: qb0_io , qb1_io , so0_io ,&
                                    so1_io , tb0_io , tb1_io ,      &
                                    ub0_io , ub1_io , vb0_io ,      &
                                    vb1_io
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

    call getmem2d(var1snd,1,kz,1,8,'mod_mppio:var1snd')
    call getmem2d(var1rcv,1,kz,1,8,'mod_mppio:var1rcv')
    call getmem2d(var2d0,1,iy,1,jxp,'mod_mppio:var2d0')
    call getmem3d(var2d1,1,iy,1,nnsg,1,jxp,'mod_mppio:var2d1')
    call getmem3d(inisrf0,1,iy,1,nnsg*4+7,1,jxp,'mod_mppio:inisrf0')
    call getmem3d(atm0,1,iy,1,kz*6+3+nnsg*3,1,jxp,'mod_mppio:atm0')
    call getmem3d(bat0,1,iym2,1,numbat,1,jxp,'mod_mppio:bat0')
    call getmem3d(rad0,1,iym2,1,nrad3d*kz+nrad2d,1,jxp,'mod_mppio:rad0')
    call getmem4d(sub0,1,iym2,1,nnsg,1,numsub,1,jxp,'mod_mppio:sub0')
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      call getmem3d(uw0,1,iy,1,kz*3,1,jxp,'mod_mppio:uw0')
    end if

    if (myid == 0) then
      call allocate_domain(mddom_io,.false.)
      call allocate_atmstate(atm1_io,ibltyp,.false.,0,0)
      call allocate_atmstate(atm2_io,ibltyp,.false.,0,0)
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call allocate_tcm_state(tcmstate_io,.false.)
      end if

      call getmem2d(kpbl_io,1,iy,1,jx,'mod_mppio:kpbl_io')

      call getmem2d(var2d_0,1,iy,1,jx,'mod_mppio:var2d_0')
      call getmem3d(var2d_1,1,iy,1,nnsg,1,jx,'mod_mppio:var2d_1')
      call getmem3d(inisrf_0,1,iy,1,nnsg*4+7,1,jx,'mod_mppio:inisrf_0')
      call getmem3d(atm_0,1,iy,1,kz*6+3+nnsg*3,1,jx,'mod_mppio:atm_0')
      call getmem3d(bat_0,1,iym2,1,numbat,1,jx,'mod_mppio:bat_0')
      call getmem3d(rad_0,1,iym2,1,nrad3d*kz+nrad2d,1,jx,'mod_mppio:rad_0')
      call getmem4d(sub_0,1,iym2,1,nnsg,1,numsub,1,jx,'mod_mppio:sub_0')
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call getmem3d(uw_0,1,iy,1,kz*3,1,jx,'mod_mppio:uw_0')
      end if

      call getmem3d(col2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:col2d_io')
      call getmem3d(dew2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:dew2d_io')
      call getmem3d(evpa2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:evpa2d_io')
      call getmem3d(gwet2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:gwet2d_io')
      call getmem3d(ircp2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:ircp2d_io')
      call getmem3d(rno2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:rno2d_io')
      call getmem3d(rnos2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:rnos2d_io')
      call getmem3d(sag2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:sag2d_io')
      call getmem3d(scv2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:scv2d_io')
      call getmem3d(sena2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:sena2d_io')
      call getmem3d(sice2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:sice2d_io')
      call getmem3d(srw2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:srw2d_io')
      call getmem3d(ssw2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:ssw2d_io')
      call getmem3d(swt2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:swt2d_io')
      call getmem3d(taf2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:taf2d_io')
      call getmem3d(text2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:text2d_io')
      call getmem3d(tg2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:tg2d_io')
      call getmem3d(tgb2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:tgb2d_io')
      call getmem3d(tlef2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:tlef2d_io')
      call getmem3d(emiss2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:emiss2d_io')
      call getmem3d(veg2d1_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:veg2d1_io')
      call getmem3d(ocld2d_io,1,nnsg,1,iym1,1,mjj,'mod_mppio:ocld2d_io')
      call getmem2d(veg2d_io,1,iym1,1,mjj,'mod_mppio:veg2d_io')
      call getmem2d(ldmsk_io,1,iym1,1,mjj,'mod_mppio:ldmsk_io')
      call getmem3d(ht1_io,1,nnsg,1,iy,1,jx,'mod_mppio:ht1_io')
      call getmem3d(lndcat1_io,1,nnsg,1,iy,1,jx,'mod_mppio:lndcat1_io')
      call getmem3d(xlat1_io,1,nnsg,1,iy,1,jx,'mod_mppio:xlat1_io')
      call getmem3d(xlon1_io,1,nnsg,1,iy,1,jx,'mod_mppio:xlon1_io')
      call getmem2d(flw2d_io,1,iym1,1,mjj,'mod_mppio:flw2d_io')
      call getmem2d(flwd2d_io,1,iym1,1,mjj,'mod_mppio:flwd2d_io')
      call getmem2d(fsw2d_io,1,iym1,1,mjj,'mod_mppio:fsw2d_io')
      call getmem2d(sabv2d_io,1,iym1,1,mjj,'mod_mppio:sabv2d_io')
      call getmem2d(sinc2d_io,1,iym1,1,mjj,'mod_mppio:sinc2d_io')
      call getmem2d(sol2d_io,1,iym1,1,mjj,'mod_mppio:sol2d_io')
      call getmem2d(solvd2d_io,1,iym1,1,mjj,'mod_mppio:solvd2d_io')
      call getmem2d(solvs2d_io,1,iym1,1,mjj,'mod_mppio:solvs2d_io')

      call getmem3d(fbat_io,1,mojj,1,iym2,1,numbat,'mod_mppio:fbat_io')
      call getmem4d(fsub_io,1,nnsg,1,mojj,1,iym2,1,numsub,'mod_mppio:fsub_io')
      call getmem3d(frad2d_io,1,mojj,1,iym2,1,nrad2d,'mod_mppio:frad2d_io')
      call getmem4d(frad3d_io,1,mojj,1,iym2,1,kz,1,nrad3d,'mod_mppio:frad3d_io')
      call getmem2d(radpsa_io,1,mojj,1,iym2,'mod_mppio:radpsa_io')

      call getmem2d(cbmf2d_io,1,iy,1,jx,'mod_mppio:cbmf2d_io')
      call getmem3d(fcc_io,1,iy,1,kz,1,jx,'mod_mppio:fcc_io')
      call getmem3d(rsheat_io,1,iy,1,kz,1,jx,'mod_mppio:rsheat_io')
      call getmem3d(rswat_io,1,iy,1,kz,1,jx,'mod_mppio:rswat_io')
      call getmem3d(dstor_io,1,iy,1,jx,1,nsplit,'mod_mppio:dstor_io')
      call getmem3d(hstor_io,1,iy,1,jx,1,nsplit,'mod_mppio:hstor_io')

      call getmem4d(absnxt_io,1,iym1,1,kz,1,4,1,mjj,'mod_mppio:absnxt_io')
      call getmem4d(abstot_io,1,iym1,1,kzp1,1,kzp1,1,mjj,'mod_mppio:abstot_io')
      call getmem3d(emstot_io,1,iym1,1,kzp1,1,mjj,'mod_mppio:emstot_io')
      call getmem3d(heatrt_io,1,iym1,1,kz,1,mjj,'mod_mppio:heatrt_io')
      call getmem3d(o3prof_io,1,iym1,1,kzp1,1,mjj,'mod_mppio:o3prof_io')
      call getmem3d(aerasp_io,1,iym1,1,kz,1,mjj,'mod_mppio:aerasp_io')
      call getmem3d(aerext_io,1,iym1,1,kz,1,mjj,'mod_mppio:aerext_io')
      call getmem3d(aerssa_io,1,iym1,1,kz,1,mjj,'mod_mppio:aerssa_io')
      call getmem2d(aersrrf_io,1,iym1,1,mjj,'mod_mppio:aersrrf_io')
      call getmem2d(aertarf_io,1,iym1,1,mjj,'mod_mppio:aertarf_io')
      call getmem2d(aertalwrf_io,1,iym1,1,mjj,'mod_mppio:aertalwrf_io')
      call getmem2d(aersrlwrf_io,1,iym1,1,mjj,'mod_mppio:aersrlwrf_io')

      call getmem2d(ps0_io,1,iy,1,jx,'mod_mppio:ps0_io')
      call getmem2d(ps1_io,1,iy,1,jx,'mod_mppio:ps1_io')
      call getmem2d(ts0_io,1,iy,1,jx,'mod_mppio:ts0_io')
      call getmem2d(ts1_io,1,iy,1,jx,'mod_mppio:ts1_io')
      call getmem3d(qb0_io,1,iy,1,kz,1,jx,'mod_mppio:qb0_io')
      call getmem3d(qb1_io,1,iy,1,kz,1,jx,'mod_mppio:qb1_io')
      call getmem3d(so0_io,1,iy,1,kz,1,jx,'mod_mppio:so0_io')
      call getmem3d(so1_io,1,iy,1,kz,1,jx,'mod_mppio:so1_io')
      call getmem3d(tb0_io,1,iy,1,kz,1,jx,'mod_mppio:tb0_io')
      call getmem3d(tb1_io,1,iy,1,kz,1,jx,'mod_mppio:tb1_io')
      call getmem3d(ub0_io,1,iy,1,kz,1,jx,'mod_mppio:ub0_io')
      call getmem3d(ub1_io,1,iy,1,kz,1,jx,'mod_mppio:ub1_io')
      call getmem3d(vb0_io,1,iy,1,kz,1,jx,'mod_mppio:vb0_io')
      call getmem3d(vb1_io,1,iy,1,kz,1,jx,'mod_mppio:vb1_io')
      call getmem2d(ui1_io,1,kz,1,jx,'mod_mppio:ui1_io')
      call getmem2d(ui2_io,1,kz,1,jx,'mod_mppio:ui2_io')
      call getmem2d(uilx_io,1,kz,1,jx,'mod_mppio:uilx_io')
      call getmem2d(uil_io,1,kz,1,jx,'mod_mppio:uil_io')
      call getmem2d(vi1_io,1,kz,1,jx,'mod_mppio:vi1_io')
      call getmem2d(vi2_io,1,kz,1,jx,'mod_mppio:vi2_io')
      call getmem2d(vilx_io,1,kz,1,jx,'mod_mppio:vilx_io')
      call getmem2d(vil_io,1,kz,1,jx,'mod_mppio:vil_io')
      call getmem2d(pptc_io,1,iym1,1,mjj,'mod_mppio:pptc_io')
      call getmem2d(pptnc_io,1,iym1,1,mjj,'mod_mppio:pptnc_io')
      call getmem2d(prca2d_io,1,iym1,1,mjj,'mod_mppio:prca2d_io')
      call getmem2d(prnca2d_io,1,iym1,1,mjj,'mod_mppio:prnca2d_io')
      call getmem2d(hfx_io,1,iy,1,jx,'mod_mppio:hfx_io')
      call getmem2d(psa_io,1,iy,1,jx,'mod_mppio:psa_io')
      call getmem2d(psb_io,1,iy,1,jx,'mod_mppio:psb_io')
      call getmem2d(qfx_io,1,iy,1,jx,'mod_mppio:qfx_io')
      call getmem2d(rainc_io,1,iy,1,jx,'mod_mppio:rainc_io')
      call getmem2d(rainnc_io,1,iy,1,jx,'mod_mppio:rainnc_io')
      call getmem2d(tga_io,1,iy,1,jx,'mod_mppio:tga_io')
      call getmem2d(tgbb_io,1,iy,1,jx,'mod_mppio:tgbb_io')
      call getmem2d(tgb_io,1,iy,1,jx,'mod_mppio:tgb_io')
      call getmem2d(uvdrag_io,1,iy,1,jx,'mod_mppio:uvdrag_io')
      call getmem2d(zpbl_io,1,iy,1,jx,'mod_mppio:zpbl_io')
      call getmem3d(omega_io,1,iy,1,kz,1,jx,'mod_mppio:omega_io')
      if (icup == 3) then
        call getmem2d(cldefi_io,1,iy,1,jx,'mod_mppio:cldefi_io')
        call getmem3d(tbase_io,1,iy,1,kz,1,jx,'mod_mppio:tbase_io')
      end if
#ifdef CLM
      call getmem2d(sols2d_io,1,iym1,1,mjj,'mod_mppio:sols2d_io')
      call getmem2d(soll2d_io,1,iym1,1,mjj,'mod_mppio:soll2d_io')
      call getmem2d(solsd2d_io,1,iym1,1,mjj,'mod_mppio:solsd2d_io')
      call getmem2d(solld2d_io,1,iym1,1,mjj,'mod_mppio:solld2d_io')
      call getmem2d(aldifl2d_io,1,iym1,1,mjj,'mod_mppio:aldifl2d_io')
      call getmem2d(aldifs2d_io,1,iym1,1,mjj,'mod_mppio:aldifs2d_io')
      call getmem2d(aldirl2d_io,1,iym1,1,mjj,'mod_mppio:aldirl2d_io')
      call getmem2d(aldirs2d_io,1,iym1,1,mjj,'mod_mppio:aldirs2d_io')
      call getmem2d(lndcat2d_io,1,iy,1,jx,'mod_mppio:lndcat2d_io')
#endif
    endif
    if (myid == 0) then
      call getmem3d(sav_0,1,iy,1,kz*4+2,1,jx,'mod_mppio:sav_0')
      call getmem3d(sav_0a,1,iy,1,kzp1+4,1,jx,'mod_mppio:sav_0a')
      call getmem3d(sav_0b,1,iy,1,kzp1,1,jx,'mod_mppio:sav_0b')
      call getmem3d(sav_0c,1,iy,1,kz*2,1,jx,'mod_mppio:sav_0c')
      call getmem3d(sav_0s,1,iy,1,kz,1,jx,'mod_mppio:sav_0s')
      call getmem3d(sav_0d,1,iy,1,nsplit*2,1,jx,'mod_mppio:sav_0d')
      call getmem3d(sav_1,1,iym1,1,kz*4+(kzp1*kzp2),1,jx,'mod_mppio:sav_1')
      call getmem3d(sav_2,1,iym1,1,nnsg*5+4,1,jx,'mod_mppio:sav_2')
      call getmem3d(sav_2a,1,iym1,1,nnsg*2+2,1,jx,'mod_mppio:sav_2a')
      if ( ichem == 1 ) then
        call getmem3d(sav_4,1,iy,1,ntr*(kz*4+1),1,jx,'mod_mppio:sav_4')
        call getmem3d(sav_4a,1,iym1,1,7,1,jx,'mod_mppio:sav_4a')
      end if
      call getmem3d(sav_6,1,kz,1,8,1,jx,'mod_mppio:sav_6')
#ifdef CLM
      call getmem3d(sav_clmout,1,iym1,1,8,1,jx,'mod_mppio:sav_clmout')
#endif
    end if

    call getmem3d(sav0,1,iy,1,kz*4+2,1,jxp,'mod_mppio:sav0')
    call getmem3d(sav0a,1,iy,1,kzp1+4,1,jxp,'mod_mppio:sav0a')
    call getmem3d(sav0b,1,iy,1,kzp1,1,jxp,'mod_mppio:sav0b')
    call getmem3d(sav0c,1,iy,1,kz*2,1,jxp,'mod_mppio:sav0c')
    call getmem3d(sav0s,1,iy,1,kz,1,jxp,'mod_mppio:sav0s')
    call getmem3d(sav0d,1,iy,1,nsplit*2,1,jxp,'mod_mppio:sav0d')
    call getmem3d(sav1,1,iym1,1,kz*4+(kzp1*kzp2),1,jxp,'mod_mppio:sav1')
    call getmem3d(sav2,1,iym1,1,nnsg*5+4,1,jxp,'mod_mppio:sav2')
    call getmem3d(sav2a,1,iym1,1,nnsg*2+2,1,jxp,'mod_mppio:sav2a')
    if ( ichem == 1 ) then
      call getmem3d(sav4,1,iy,1,ntr*(kz*4+1),1,jxp,'mod_mppio:sav4')
      call getmem3d(sav4a,1,iym1,1,7,1,jxp,'mod_mppio:sav4a')
    end if
    call getmem3d(sav6,1,kz,1,8,1,jxp,'mod_mppio:sav6')
#ifdef CLM
    call getmem3d(sav_clmin,1,iym1,1,8,1,jxp,'mod_mppio:sav_clmin')
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
