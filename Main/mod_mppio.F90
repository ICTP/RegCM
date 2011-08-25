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
  use mod_outrad , only : nrad2d , nrad3d
  use mod_tcm_interface , only : tcm_state , allocate_tcm_state
  use mod_memutil
  use mod_message
  use mod_che_mppio
  use mod_bats_mppio
  use mod_cu_common
!
  real(8) , pointer , dimension(:,:,:,:) :: spacesubm1
  real(8) , pointer , dimension(:,:,:,:) :: spacesub
  real(8) , pointer , dimension(:,:,:) :: spacebat
  real(8) , pointer , dimension(:,:,:) :: space2d
  real(8) , pointer , dimension(:,:,:,:) :: space3d
  real(8) , pointer , dimension(:,:,:) :: spacev
  real(8) , pointer , dimension(:,:,:) :: spacesurf
  private :: spacesubm1 , spacesub , spacebat
  private :: space2d , space3d , spacev , spacesurf
#ifdef CLM
  real(8) , pointer , dimension(:,:,:) :: spaceclm
  private :: spaceclm
#endif
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
      call allocate_atmstate(atm1_io,.false.,0,0)
      call allocate_atmstate(atm2_io,.false.,0,0)
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

      if (lband) then
        call getmem4d(spacesubm1,1,nnsg,1,iym1,1,jx,1,20,'mod_mppio:spacesubm1')
      else
        call getmem4d(spacesubm1,1,nnsg,1,iym1, &
                                 1,jxm1,1,20,'mod_mppio:spacesubm1')
      end if
      col2d_io  => spacesubm1(:,:,:,1)
      dew2d_io  => spacesubm1(:,:,:,2)
      evpa2d_io => spacesubm1(:,:,:,3)
      gwet2d_io => spacesubm1(:,:,:,4)
      ircp2d_io => spacesubm1(:,:,:,5)
      rno2d_io  => spacesubm1(:,:,:,6)
      rnos2d_io => spacesubm1(:,:,:,7)
      sag2d_io  => spacesubm1(:,:,:,8)
      scv2d_io  => spacesubm1(:,:,:,9)
      sena2d_io => spacesubm1(:,:,:,10)
      sice2d_io => spacesubm1(:,:,:,11)
      srw2d_io  => spacesubm1(:,:,:,12)
      ssw2d_io  => spacesubm1(:,:,:,13)
      swt2d_io  => spacesubm1(:,:,:,14)
      taf2d_io  => spacesubm1(:,:,:,15)
      text2d_io => spacesubm1(:,:,:,16)
      tg2d_io   => spacesubm1(:,:,:,17)
      tgb2d_io  => spacesubm1(:,:,:,18)
      tlef2d_io => spacesubm1(:,:,:,19)
      emiss2d_io => spacesubm1(:,:,:,20)
      if (lband) then
        call getmem3d(veg2d1_io,1,nnsg,1,iym1,1,jx,'mod_mppio:veg2d1_io')
        call getmem3d(ocld2d_io,1,nnsg,1,iym1,1,jx,'mod_mppio:ocld2d_io')
        call getmem2d(veg2d_io,1,iym1,1,jx,'mod_mppio:veg2d_io')
        call getmem2d(ldmsk_io,1,iym1,1,jx,'mod_mppio:ldmsk_io')
      else
        call getmem3d(veg2d1_io,1,nnsg,1,iym1,1,jxm1,'mod_mppio:veg2d1_io')
        call getmem3d(ocld2d_io,1,nnsg,1,iym1,1,jxm1,'mod_mppio:ocld2d_io')
        call getmem2d(veg2d_io,1,iym1,1,jxm1,'mod_mppio:veg2d_io')
        call getmem2d(ldmsk_io,1,iym1,1,jxm1,'mod_mppio:ldmsk_io')
      end if
      call getmem4d(spacesub,1,nnsg,1,iy,1,jx,1,4,'mod_mppio:spacesub')
      ht1_io     => spacesub(:,:,:,1)
      lndcat1_io => spacesub(:,:,:,2)
      xlat1_io   => spacesub(:,:,:,3)
      xlon1_io   => spacesub(:,:,:,4)
      if (lband) then
        call getmem3d(spacebat,1,iym1,1,jx,1,8,'mod_mppio:spacebat')
      else
        call getmem3d(spacebat,1,iym1,1,jxm1,1,8,'mod_mppio:spacebat')
      end if
      flw2d_io      => spacebat(:,:,1)
      flwd2d_io     => spacebat(:,:,2)
      fsw2d_io      => spacebat(:,:,3)
      sabv2d_io     => spacebat(:,:,4)
      sinc2d_io     => spacebat(:,:,5)
      sol2d_io      => spacebat(:,:,6)
      solvd2d_io    => spacebat(:,:,7)
      solvs2d_io    => spacebat(:,:,8)
      if (lband) then
        call getmem3d(fbat_io,1,jx,1,iym2,1,numbat,'mod_mppio:fbat_io')
        call getmem4d(fsub_io,1,nnsg,1,jx,1,iym2,1,numsub,'mod_mppio:fsub_io')
        call getmem3d(frad2d_io,1,jx,1,iym2,1,nrad2d,'mod_mppio:frad2d_io')
        call getmem4d(frad3d_io,1,jx,1,iym2,1,kz,1,nrad3d,'mod_mppio:frad3d_io')
        call getmem2d(radpsa_io,1,jx,1,iym2,'mod_mppio:radpsa_io')
      else
        call getmem3d(fbat_io,1,jxm2,1,iym2,1,numbat,'mod_mppio:fbat_io')
        call getmem4d(fsub_io,1,nnsg,1,jxm2,1,iym2,1,numsub,'mod_mppio:fsub_io')
        call getmem3d(frad2d_io,1,jxm2,1,iym2,1,nrad2d,'mod_mppio:frad2d_io')
        call getmem4d(frad3d_io,1,jxm2,1,iym2, &
                                1,kz,1,nrad3d,'mod_mppio:frad3d_io')
        call getmem2d(radpsa_io,1,jxm2,1,iym2,'mod_mppio:radpsa_io')
      end if

      call getmem2d(cbmf2d_io,1,iy,1,jx,'mod_mppio:cbmf2d_io')
      call getmem3d(fcc_io,1,iy,1,kz,1,jx,'mod_mppio:fcc_io')
      call getmem3d(rsheat_io,1,iy,1,kz,1,jx,'mod_mppio:rsheat_io')
      call getmem3d(rswat_io,1,iy,1,kz,1,jx,'mod_mppio:rswat_io')
      call getmem3d(dstor_io,1,iy,1,jx,1,nsplit,'mod_mppio:dstor_io')
      call getmem3d(hstor_io,1,iy,1,jx,1,nsplit,'mod_mppio:hstor_io')
      if (lband) then
        call getmem4d(absnxt_io,1,iym1,1,kz,1,4,1,jx,'mod_mppio:absnxt_io')
        call getmem4d(abstot_io,1,iym1,1,kzp1,1,kzp1,1,jx,'mod_mppio:abstot_io')
        call getmem3d(emstot_io,1,iym1,1,kzp1,1,jx,'mod_mppio:emstot_io')
        call getmem3d(heatrt_io,1,iym1,1,kz,1,jx,'mod_mppio:heatrt_io')
        call getmem3d(o3prof_io,1,iym1,1,kzp1,1,jx,'mod_mppio:o3prof_io')
        call getmem3d(aerasp_io,1,iym1,1,kz,1,jx,'mod_mppio:aerasp_io')
        call getmem3d(aerext_io,1,iym1,1,kz,1,jx,'mod_mppio:aerext_io')
        call getmem3d(aerssa_io,1,iym1,1,kz,1,jx,'mod_mppio:aerssa_io')
        call getmem2d(aersrrf_io,1,iym1,1,jx,'mod_mppio:aersrrf_io')
        call getmem2d(aertarf_io,1,iym1,1,jx,'mod_mppio:aertarf_io')
        call getmem2d(aertalwrf_io,1,iym1,1,jx,'mod_mppio:aertalwrf_io')
        call getmem2d(aersrlwrf_io,1,iym1,1,jx,'mod_mppio:aersrlwrf_io')
      else
        call getmem4d(absnxt_io,1,iym1,1,kz,1,4,1,jxm1,'mod_mppio:absnxt_io')
        call getmem4d(abstot_io,1,iym1,1,kzp1, &
                                1,kzp1,1,jxm1,'mod_mppio:abstot_io')
        call getmem3d(emstot_io,1,iym1,1,kzp1,1,jxm1,'mod_mppio:emstot_io')
        call getmem3d(heatrt_io,1,iym1,1,kz,1,jxm1,'mod_mppio:heatrt_io')
        call getmem3d(o3prof_io,1,iym1,1,kzp1,1,jxm1,'mod_mppio:o3prof_io')
        call getmem3d(aerasp_io,1,iym1,1,kz,1,jxm1,'mod_mppio:aerasp_io')
        call getmem3d(aerext_io,1,iym1,1,kz,1,jxm1,'mod_mppio:aerext_io')
        call getmem3d(aerssa_io,1,iym1,1,kz,1,jxm1,'mod_mppio:aerssa_io')
        call getmem2d(aersrrf_io,1,iym1,1,jxm1,'mod_mppio:aersrrf_io')
        call getmem2d(aertarf_io,1,iym1,1,jxm1,'mod_mppio:aertarf_io')
        call getmem2d(aertalwrf_io,1,iym1,1,jxm1,'mod_mppio:aertalwrf_io')
        call getmem2d(aersrlwrf_io,1,iym1,1,jxm1,'mod_mppio:aersrlwrf_io')
      end if
      call getmem3d(space2d,1,iy,1,jx,1,4,'mod_mppio:space2d')
      ps0_io => space2d(:,:,1)
      ps1_io => space2d(:,:,2)
      ts0_io => space2d(:,:,3)
      ts1_io => space2d(:,:,4)
      call getmem4d(space3d,1,iy,1,kz,1,jx,1,10,'mod_mppio:space3d')
      qb0_io => space3d(:,:,:,1)
      qb1_io => space3d(:,:,:,2)
      so0_io => space3d(:,:,:,3)
      so1_io => space3d(:,:,:,4)
      tb0_io => space3d(:,:,:,5)
      tb1_io => space3d(:,:,:,6)
      ub0_io => space3d(:,:,:,7)
      ub1_io => space3d(:,:,:,8)
      vb0_io => space3d(:,:,:,9)
      vb1_io => space3d(:,:,:,10)
      call getmem3d(spacev,1,kz,1,jx,1,8,'mod_mppio:spacev')
      ui1_io  => spacev(:,:,1)
      ui2_io  => spacev(:,:,2)
      uilx_io => spacev(:,:,3)
      uil_io  => spacev(:,:,4)
      vi1_io  => spacev(:,:,5)
      vi2_io  => spacev(:,:,6)
      vilx_io => spacev(:,:,7)
      vil_io  => spacev(:,:,8)
      if (lband) then
        call getmem2d(pptc_io,1,iym1,1,jx,'mod_mppio:pptc_io')
        call getmem2d(pptnc_io,1,iym1,1,jx,'mod_mppio:pptnc_io')
        call getmem2d(prca2d_io,1,iym1,1,jx,'mod_mppio:prca2d_io')
        call getmem2d(prnca2d_io,1,iym1,1,jx,'mod_mppio:prnca2d_io')
      else
        call getmem2d(pptc_io,1,iym1,1,jxm1,'mod_mppio:pptc_io')
        call getmem2d(pptnc_io,1,iym1,1,jxm1,'mod_mppio:pptnc_io')
        call getmem2d(prca2d_io,1,iym1,1,jxm1,'mod_mppio:prca2d_io')
        call getmem2d(prnca2d_io,1,iym1,1,jxm1,'mod_mppio:prnca2d_io')
      end if
      if (icup == 3) then
        call getmem3d(spacesurf,1,iy,1,jx,1,11,'mod_mppio:spacesurf')
      else
        call getmem3d(spacesurf,1,iy,1,jx,1,12,'mod_mppio:spacesurf')
      end if
      hfx_io    => spacesurf(:,:,1)
      psa_io    => spacesurf(:,:,2)
      psb_io    => spacesurf(:,:,3)
      qfx_io    => spacesurf(:,:,4)
      rainc_io  => spacesurf(:,:,5)
      rainnc_io => spacesurf(:,:,6)
      tga_io    => spacesurf(:,:,7)
      tgbb_io   => spacesurf(:,:,8)
      tgb_io    => spacesurf(:,:,9)
      uvdrag_io => spacesurf(:,:,10)
      zpbl_io   => spacesurf(:,:,11)
      if (icup == 3) cldefi_io => spacesurf(:,:,12)
      call getmem3d(omega_io,1,iy,1,kz,1,jx,'mod_mppio:omega_io')
      if (icup == 3) then
        call getmem3d(tbase_io,1,iy,1,kz,1,jx,'mod_mppio:tbase_io')
      end if
#ifdef CLM
      if (lband) then
        call getmem3d(spaceclm,1,iym1,1,jx,1,8,'mod_mppio:spaceclm')
      else
        call getmem3d(spaceclm,1,iym1,1,jxm1,1,8,'mod_mppio:spaceclm')
      end if
      sols2d_io   => spaceclm(:,:,1)
      soll2d_io   => spaceclm(:,:,2)
      solsd2d_io  => spaceclm(:,:,3)
      solld2d_io  => spaceclm(:,:,4)
      aldifl2d_io => spaceclm(:,:,5)
      aldirs2d_io => spaceclm(:,:,6)
      aldirl2d_io => spaceclm(:,:,7)
      aldifs2d_io => spaceclm(:,:,8)
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
