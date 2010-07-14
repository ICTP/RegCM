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

      use mod_dynparam

      implicit none

#ifdef MPP1
!
      real(8) , allocatable , target , dimension(:,:,:,:) :: spacesubm1
      real(8) , allocatable , target , dimension(:,:,:,:) :: spacesub
      real(8) , allocatable , target , dimension(:,:,:) :: spacebat
      real(8) , allocatable , target , dimension(:,:,:) :: space2d
      real(8) , allocatable , target , dimension(:,:,:,:) :: space3d
      real(8) , allocatable , target , dimension(:,:,:) :: spacev
      real(8) , allocatable , target , dimension(:,:,:) :: spacesurf
      real(8) , allocatable , target , dimension(:,:,:,:) :: spaceair
      private :: spacesubm1 , spacesub , spacebat
      private :: space2d , space3d , spacev , spacesurf , spaceair
#ifdef CLM
      real(8) , allocatable , target , dimension(:,:,:) :: spaceclm
      private :: spaceclm
#endif
      real(8) , pointer , dimension(:,:,:) :: col2d_io
      real(8) , pointer , dimension(:,:,:) :: dew2d_io
      real(8) , pointer , dimension(:,:,:) :: evpa2d_io
      real(8) , pointer , dimension(:,:,:) :: gwet2d_io 
      real(8) , pointer , dimension(:,:,:) :: ircp2d_io
      ! ocean/land mask
      real(8) , pointer , dimension(:,:,:) :: ocld2d_io
      real(8) , pointer , dimension(:,:,:) :: rno2d_io
      real(8) , pointer , dimension(:,:,:) :: rnos2d_io
      real(8) , pointer , dimension(:,:,:) :: sag2d_io
      real(8) , pointer , dimension(:,:,:) :: scv2d_io
      real(8) , pointer , dimension(:,:,:) :: sena2d_io
      real(8) , pointer , dimension(:,:,:) :: sice2d_io
      real(8) , pointer , dimension(:,:,:) :: srw2d_io
      real(8) , pointer , dimension(:,:,:) :: ssw2d_io
      ! total soil wetness of column
      real(8) , pointer , dimension(:,:,:) :: swt2d_io
      ! 2-m air temperature
      real(8) , pointer , dimension(:,:,:) :: taf2d_io
      real(8) , pointer , dimension(:,:,:) :: text2d_io
      real(8) , pointer , dimension(:,:,:) :: tg2d_io
      real(8) , pointer , dimension(:,:,:) :: tgb2d_io
      real(8) , pointer , dimension(:,:,:) :: tlef2d_io
      real(8) , pointer , dimension(:,:,:) :: veg2d1_io

      real(8) , pointer , dimension(:,:,:) :: ht1_io , satbrt1_io ,     &
                                         &    snowc_io
!
      ! Downward Longwave radiation
      real(8) , pointer , dimension(:,:) :: flwd2d_io
      ! Net Surface Longwave
      real(8) , pointer , dimension(:,:) :: flw2d_io
      ! surface absorbed radiation
      real(8) , pointer , dimension(:,:) :: fsw2d_io
      real(8) , pointer , dimension(:,:) :: sabv2d_io
      real(8) , pointer , dimension(:,:) :: sdelqk2d_io
      real(8) , pointer , dimension(:,:) :: sdeltk2d_io
      real(8) , pointer , dimension(:,:) :: sfracb2d_io
      real(8) , pointer , dimension(:,:) :: sfracs2d_io
      real(8) , pointer , dimension(:,:) :: sfracv2d_io
      ! total incident solar radiation
      real(8) , pointer , dimension(:,:) :: sinc2d_io
      real(8) , pointer , dimension(:,:) :: sol2d_io
      real(8) , pointer , dimension(:,:) :: solvd2d_io
      real(8) , pointer , dimension(:,:) :: solvs2d_io
      real(8) , pointer , dimension(:,:) :: ssw2da_io
      real(8) , pointer , dimension(:,:) :: svegfrac2d_io
      real(8) , pointer , dimension(:,:) :: veg2d_io
      
      real(4) , allocatable , dimension(:,:,:) :: fbat_io
      real(4) , allocatable , dimension(:,:,:,:) :: fsub_io

      real(4) , allocatable , dimension(:,:,:) :: frad2d_io
      real(4) , allocatable , dimension(:,:,:,:) :: frad3d_io

      real(8) , allocatable , dimension(:,:) :: cbmf2d_io
      real(8) , allocatable , dimension(:,:,:) :: fcc_io , rsheat_io ,  &
                                      & rswat_io

      real(8) , allocatable , dimension(:,:,:,:) :: absnxt_io
      real(8) , allocatable , dimension(:,:,:,:) :: abstot_io
      real(8) , allocatable , dimension(:,:,:) :: emstot_io

      ! residual heat
      real(8) , allocatable , dimension(:,:,:) :: heatrt_io
      ! ozone profile
      real(8) , allocatable , dimension(:,:,:) :: o3prof_io

      real(8) , allocatable , dimension(:,:,:) :: dstor_io , hstor_io

      real(8) , allocatable , dimension(:,:,:) :: aerasp_io ,           &
                                & aerext_io , aerssa_io
      real(8) , allocatable , dimension(:,:) :: aersrrf_io , aertarf_io,&
                              & aertalwrf_io , aersrlwrf_io
      real(8) , allocatable , dimension(:,:,:) :: cemtrac_io ,          &
                                & cemtr_io , wxaq_io , wxsg_io
      real(8) , allocatable , dimension(:,:,:,:) :: rxsaq1_io ,         &
                                & rxsaq2_io , rxsg_io

      real(8) , allocatable , dimension(:,:,:,:) :: remcvc_io
      real(8) , allocatable , dimension(:,:,:,:) :: remlsc_io
      real(8) , allocatable , dimension(:,:,:) :: remdrd_io

      real(8) , pointer , dimension(:,:) :: ps0_io , ps1_io , ts0_io ,  &
                            &  ts1_io
      real(8) , pointer , dimension(:,:,:) :: qb0_io , qb1_io , so0_io ,&
                                      & so1_io , tb0_io , tb1_io ,      &
                                      & ub0_io , ub1_io , vb0_io ,      &
                                      & vb1_io
      real(8) , pointer , dimension(:,:) :: ui1_io , ui2_io , uilx_io , &
                                   & uil_io , vi1_io , vi2_io ,         &
                                   & vilx_io , vil_io

      real(8) , allocatable , dimension(:,:,:,:) :: chemsrc_io
      real(8) , allocatable , dimension(:,:,:) :: ddsfc_io , dtrace_io ,&
                                   & wdcvc_io , wdlsc_io
      real(8) , allocatable , dimension(:,:) :: pptc_io , pptnc_io ,    &
                                   & prca2d_io , prnca2d_io

      real(8) , allocatable , dimension(:,:,:,:) :: chia_io , chib_io

      real(8) , pointer , dimension(:,:) :: cldefi_io
      real(8) , pointer , dimension(:,:) :: f_io
      real(8) , pointer , dimension(:,:) :: hfx_io
      real(8) , pointer , dimension(:,:) :: htsd_io
      real(8) , pointer , dimension(:,:) :: ht_io
      real(8) , pointer , dimension(:,:) :: msfd_io
      real(8) , pointer , dimension(:,:) :: msfx_io
      real(8) , pointer , dimension(:,:) :: psa_io
      real(8) , pointer , dimension(:,:) :: psb_io
      real(8) , pointer , dimension(:,:) :: qfx_io
      real(8) , pointer , dimension(:,:) :: rainc_io
      real(8) , pointer , dimension(:,:) :: rainnc_io
      real(8) , pointer , dimension(:,:) :: satbrt_io
      real(8) , pointer , dimension(:,:) :: tga_io
      ! ground blackbody temperature
      real(8) , pointer , dimension(:,:) :: tgbb_io
      real(8) , pointer , dimension(:,:) :: tgb_io
      real(8) , pointer , dimension(:,:) :: uvdrag_io
      real(8) , pointer , dimension(:,:) :: xlat_io
      real(8) , pointer , dimension(:,:) :: xlong_io
      real(8) , pointer , dimension(:,:) :: zpbl_io

      real(8) , pointer , dimension(:,:,:) :: omega_io , qca_io ,       &
                              & qcb_io , qva_io , qvb_io , ta_io ,      &
                                      & tbase_io , tb_io , ua_io ,      &
                                      & ub_io , va_io , vb_io

      real(8) , allocatable , dimension(:,:,:) :: inisrf0
      real(8) , allocatable , dimension(:,:,:) :: inisrf_0

      real(8) , allocatable , dimension(:,:) :: var1snd , var1rcv
 
      real(8) , allocatable , dimension(:,:,:) :: atm0
      real(8) , allocatable , dimension(:,:,:) :: atm_0
      real(4) , allocatable , dimension(:,:,:) :: bat0
      real(4) , allocatable , dimension(:,:,:) :: bat_0
      real(8) , allocatable , dimension(:,:,:) :: out0
      real(8) , allocatable , dimension(:,:,:) :: out_0
      real(4) , allocatable , dimension(:,:,:) :: rad0
      real(4) , allocatable , dimension(:,:,:) :: rad_0
      real(4) , allocatable , dimension(:,:,:,:) :: sub0
      real(4) , allocatable , dimension(:,:,:,:) :: sub_0
      real(8) , allocatable , dimension(:,:,:) :: chem0
      real(8) , allocatable , dimension(:,:,:) :: chem_0

#ifdef CLM
      ! Direct Shortwave radiation
      real(8) , pointer , dimension(:,:) :: sols2d_io
      ! Direct Longwave radiation
      real(8) , pointer , dimension(:,:) :: soll2d_io
      ! Diffuse Shortwave radiation
      real(8) , pointer , dimension(:,:) :: solsd2d_io
      ! Diffuse Longwave radiation
      real(8) , pointer , dimension(:,:) :: solld2d_io
      ! Albedo of Direct Shortwave radiation
      real(8) , pointer , dimension(:,:) :: aldirs2d_io
      ! Albedo of Direct Longwave radiation
      real(8) , pointer , dimension(:,:) :: aldirl2d_io
      ! Albedo of Diffuse Shortwave radiation
      real(8) , pointer , dimension(:,:) :: aldifs2d_io
      ! Albedo of Diffuse Longwave radiation
      real(8) , pointer , dimension(:,:) :: aldifl2d_io
      ! cosine of solar zenith angle
      real(8) , pointer , dimension(:,:) :: coszrs2d_io
      ! 2-m mixing ratio
      real(8) , pointer , dimension(:,:) :: q2m_io
#endif
#endif

      contains 
!
!     This routines allocate all the arrays contained in the module
!
      subroutine allocate_mod_mppio
      use mod_dust , only : nats
      implicit none

#ifdef MPP1
#ifdef BAND
      allocate(spacesubm1(nnsg,iym1,jx,21))
#else
      allocate(spacesubm1(nnsg,iym1,jxm1,21))
#endif
      col2d_io  => spacesubm1(:,:,:,1)
      dew2d_io  => spacesubm1(:,:,:,2)
      evpa2d_io => spacesubm1(:,:,:,3)
      gwet2d_io => spacesubm1(:,:,:,4)
      ircp2d_io => spacesubm1(:,:,:,5)
      ocld2d_io => spacesubm1(:,:,:,6)
      rno2d_io  => spacesubm1(:,:,:,7)
      rnos2d_io => spacesubm1(:,:,:,8)
      sag2d_io  => spacesubm1(:,:,:,9)
      scv2d_io  => spacesubm1(:,:,:,10)
      sena2d_io => spacesubm1(:,:,:,11)
      sice2d_io => spacesubm1(:,:,:,12)
      srw2d_io  => spacesubm1(:,:,:,13)
      ssw2d_io  => spacesubm1(:,:,:,14)
      swt2d_io  => spacesubm1(:,:,:,15)
      taf2d_io  => spacesubm1(:,:,:,16)
      text2d_io => spacesubm1(:,:,:,17)
      tg2d_io   => spacesubm1(:,:,:,18)
      tgb2d_io  => spacesubm1(:,:,:,19)
      tlef2d_io => spacesubm1(:,:,:,20)
      veg2d1_io => spacesubm1(:,:,:,21)
      allocate(spacesub(nnsg,iy,jx,3))
      ht1_io     => spacesub(:,:,:,1)
      satbrt1_io => spacesub(:,:,:,2)
      snowc_io   => spacesub(:,:,:,3)
#ifdef BAND
      allocate(spacebat(iym1,jx,16))
#else
      allocate(spacebat(iym1,jxm1,16))
#endif
      flw2d_io      => spacebat(:,:,1)
      flwd2d_io     => spacebat(:,:,2)
      fsw2d_io      => spacebat(:,:,3)
      sabv2d_io     => spacebat(:,:,4)
      sdelqk2d_io   => spacebat(:,:,5)
      sdeltk2d_io   => spacebat(:,:,6)
      sfracb2d_io   => spacebat(:,:,7)
      sfracs2d_io   => spacebat(:,:,8)
      sfracv2d_io   => spacebat(:,:,9)
      sinc2d_io     => spacebat(:,:,10)
      sol2d_io      => spacebat(:,:,11)
      solvd2d_io    => spacebat(:,:,12)
      solvs2d_io    => spacebat(:,:,13)
      ssw2da_io     => spacebat(:,:,14)
      svegfrac2d_io => spacebat(:,:,15)
      veg2d_io      => spacebat(:,:,16)
#ifdef BAND
      allocate(fbat_io(jx,iym2,numbat))
      allocate(fsub_io(nnsg,jx,iym2,numsub)) 
      allocate(frad2d_io(jx,iym2,nrad2d))
      allocate(frad3d_io(jx,iym2,kz,nrad3d))
#else
      allocate(fbat_io(jxm2,iym2,numbat))
      allocate(fsub_io(nnsg,jxm2,iym2,numsub)) 
      allocate(frad2d_io(jxm2,iym2,nrad2d))
      allocate(frad3d_io(jxm2,iym2,kz,nrad3d))
#endif
      allocate(cbmf2d_io(iy,jx))
      allocate(fcc_io(iy,kz,jx))
      allocate(rsheat_io(iy,kz,jx))
      allocate(rswat_io(iy,kz,jx))
      allocate(dstor_io(iy,jx,nsplit))
      allocate(hstor_io(iy,jx,nsplit))
#ifdef BAND
      allocate(absnxt_io(iym1,kz,4,jx))
      allocate(abstot_io(iym1,kzp1,kz + 1,jx))
      allocate(emstot_io(iym1,kzp1,jx))
      allocate(heatrt_io(iym1,kz,jx))
      allocate(o3prof_io(iym1,kzp1,jx))
      allocate(aerasp_io(iym1,kz,jx))
      allocate(aerext_io(iym1,kz,jx))
      allocate(aerssa_io(iym1,kz,jx))
      allocate(aersrrf_io(iym1,jx))
      allocate(aertarf_io(iym1,jx))
      allocate(aertalwrf_io(iym1,jx))
      allocate(aersrlwrf_io(iym1,jx))
#else
      allocate(absnxt_io(iym1,kz,4,jxm1))
      allocate(abstot_io(iym1,kzp1,kz + 1,jxm1))
      allocate(emstot_io(iym1,kzp1,jxm1))
      allocate(heatrt_io(iym1,kz,jxm1))
      allocate(o3prof_io(iym1,kzp1,jxm1))
      allocate(aerasp_io(iym1,kz,jxm1))
      allocate(aerext_io(iym1,kz,jxm1))
      allocate(aerssa_io(iym1,kz,jxm1))
      allocate(aersrrf_io(iym1,jxm1))
      allocate(aertarf_io(iym1,jxm1))
      allocate(aertalwrf_io(iym1,jxm1))
      allocate(aersrlwrf_io(iym1,jxm1))
#endif
      allocate(cemtrac_io(iy,jx,ntr))
      allocate(cemtr_io(iy,jx,ntr))
      allocate(wxaq_io(iy,jx,ntr))
      allocate(wxsg_io(iy,jx,ntr))
      allocate(rxsaq1_io(iy,kz,jx,ntr))
      allocate(rxsaq2_io(iy,kz,jx,ntr))
      allocate(rxsg_io(iy,kz,jx,ntr))
      allocate(remcvc_io(iy,kz,jx,ntr))
      allocate(remlsc_io(iy,kz,jx,ntr))
      allocate(remdrd_io(iy,jx,ntr))
      allocate(space2d(iy,jx,4))
      ps0_io => space2d(:,:,1)
      ps1_io => space2d(:,:,2)
      ts0_io => space2d(:,:,3)
      ts1_io => space2d(:,:,4)
      allocate(space3d(iy,kz,jx,10))
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
      allocate(spacev(kz,jx,8))
      ui1_io  => spacev(:,:,1)
      ui2_io  => spacev(:,:,2)
      uilx_io => spacev(:,:,3)
      uil_io  => spacev(:,:,4)
      vi1_io  => spacev(:,:,5)
      vi2_io  => spacev(:,:,6)
      vilx_io => spacev(:,:,7)
      vil_io  => spacev(:,:,8)
      allocate(chemsrc_io(iy,jx,nats,ntr))
      allocate(ddsfc_io(iy,jx,ntr))
      allocate(dtrace_io(iy,jx,ntr))
      allocate(wdcvc_io(iy,jx,ntr))
      allocate(wdlsc_io(iy,jx,ntr))
#ifdef BAND
      allocate(pptc_io(iym1,jx))
      allocate(pptnc_io(iym1,jx))
      allocate(prca2d_io(iym1,jx))
      allocate(prnca2d_io(iym1,jx))
#else
      allocate(pptc_io(iym1,jxm1))
      allocate(pptnc_io(iym1,jxm1))
      allocate(prca2d_io(iym1,jxm1))
      allocate(prnca2d_io(iym1,jxm1))
#endif
      allocate(chia_io(iy,kz,jx,ntr))
      allocate(chib_io(iy,kz,jx,ntr))
      allocate(spacesurf(iy,jx,20))
      cldefi_io => spacesurf(:,:,1)
      f_io      => spacesurf(:,:,2)
      hfx_io    => spacesurf(:,:,3)
      htsd_io   => spacesurf(:,:,4)
      ht_io     => spacesurf(:,:,5)
      msfd_io   => spacesurf(:,:,6)
      msfx_io   => spacesurf(:,:,7)
      psa_io    => spacesurf(:,:,8)
      psb_io    => spacesurf(:,:,9)
      qfx_io    => spacesurf(:,:,10)
      rainc_io  => spacesurf(:,:,11)
      rainnc_io => spacesurf(:,:,12)
      satbrt_io => spacesurf(:,:,13)
      tga_io    => spacesurf(:,:,14)
      tgbb_io   => spacesurf(:,:,15)
      tgb_io    => spacesurf(:,:,16)
      uvdrag_io => spacesurf(:,:,17)
      xlat_io   => spacesurf(:,:,18)
      xlong_io  => spacesurf(:,:,19)
      zpbl_io   => spacesurf(:,:,20)
      allocate(spaceair(iy,kz,jx,12))
      omega_io => spaceair(:,:,:,1)
      qca_io   => spaceair(:,:,:,2)
      qcb_io   => spaceair(:,:,:,3)
      qva_io   => spaceair(:,:,:,4)
      qvb_io   => spaceair(:,:,:,5)
      ta_io    => spaceair(:,:,:,6)
      tbase_io => spaceair(:,:,:,7)
      tb_io    => spaceair(:,:,:,8)
      ua_io    => spaceair(:,:,:,9)
      ub_io    => spaceair(:,:,:,10)
      va_io    => spaceair(:,:,:,11)
      vb_io    => spaceair(:,:,:,12)
      allocate(inisrf0(iy,nnsg*3+8,jxp))
      allocate(inisrf_0(iy,nnsg*3+8,jx))
      allocate(var1snd(kz,8))
      allocate(var1rcv(kz,8))
      allocate(atm0(iy,kz*6+3+nnsg*4,jxp))
      allocate(atm_0(iy,kz*6+3+nnsg*4,jx))
      allocate(bat0(iym2,numbat,jxp))
      allocate(bat_0(iym2,numbat,jx))
      allocate(out0(iy,3,jxp))
      allocate(out_0(iy,3,jx))
      allocate(rad0(iym2,nrad3d*kz+nrad2d,jxp))
      allocate(rad_0(iym2,nrad3d*kz+nrad2d,jx))
      allocate(sub0(iym2,nnsg,numsub,jxp))
      allocate(sub_0(iym2,nnsg,numsub,jx))
      allocate(chem0(iy,ntr*kz+kz*3+ntr*7+5,jxp))
      allocate(chem_0(iy,ntr*kz+kz*3+ntr*7+5,jx))
#ifdef CLM
#ifdef BAND
      allocate(spaceclm(iym1,jx,9))
#else
      allocate(spaceclm(iym1,jxm1,10))
#endif
      sols2d_io   => spaceclm(:,:,1)
      soll2d_io   => spaceclm(:,:,2)
      solsd2d_io  => spaceclm(:,:,3)
      solld2d_io  => spaceclm(:,:,4)
      aldifl2d_io => spaceclm(:,:,5)
      aldirs2d_io => spaceclm(:,:,6)
      aldirl2d_io => spaceclm(:,:,7)
      aldifs2d_io => spaceclm(:,:,8)
      coszrs2d_io => spaceclm(:,:,9)
      q2m_io      => spaceclm(:,:,10)
#endif       
#endif       
      end subroutine allocate_mod_mppio
      
      end module mod_mppio
