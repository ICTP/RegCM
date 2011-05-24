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
      use mod_main , only : atmstate , allocate_atmstate
      use mod_main , only : domain , allocate_domain
      use mod_outrad , only : nrad2d , nrad3d
      use mod_message
!
      real(8) , allocatable , target , dimension(:,:,:,:) :: spacesubm1
      real(8) , allocatable , target , dimension(:,:,:,:) :: spacesub
      real(8) , allocatable , target , dimension(:,:,:) :: spacebat
      real(8) , allocatable , target , dimension(:,:,:) :: space2d
      real(8) , allocatable , target , dimension(:,:,:,:) :: space3d
      real(8) , allocatable , target , dimension(:,:,:) :: spacev
      real(8) , allocatable , target , dimension(:,:,:) :: spacesurf
      private :: spacesubm1 , spacesub , spacebat
      private :: space2d , space3d , spacev , spacesurf
#ifdef CLM
      real(8) , allocatable , target , dimension(:,:,:) :: spaceclm
      private :: spaceclm
#endif
      integer , allocatable , dimension(:,:,:) :: ocld2d_io , veg2d1_io
      integer , allocatable , dimension(:,:) :: veg2d_io , ldmsk_io

      real(8) , pointer , dimension(:,:,:) :: col2d_io , dew2d_io ,     &
           & evpa2d_io , gwet2d_io , ircp2d_io , rno2d_io , &
           & rnos2d_io , sag2d_io , scv2d_io , sena2d_io , sice2d_io ,  &
           & srw2d_io , ssw2d_io , swt2d_io , taf2d_io , text2d_io ,    &
           & tg2d_io , tgb2d_io , tlef2d_io , emiss2d_io

      integer , allocatable , dimension(:,:,:) :: idep2d_io
      real(8) , allocatable , dimension(:,:,:) :: dhlake1_io
      real(8) , allocatable , dimension(:,:,:) :: eta2d_io
      real(8) , allocatable , dimension(:,:,:) :: hi2d_io
      real(8) , allocatable , dimension(:,:,:) :: aveice2d_io
      real(8) , allocatable , dimension(:,:,:) :: hsnow2d_io
      real(8) , allocatable , dimension(:,:,:) :: evl2d_io
      real(8) , allocatable , dimension(:,:,:,:) :: tlak3d_io

      real(8) , pointer , dimension(:,:,:) :: ht1_io , satbrt1_io ,     &
                                         &    xlat1_io , xlon1_io
!
      real(8) , pointer , dimension(:,:) :: flw2d_io , flwd2d_io ,      &
                                     & fsw2d_io , sabv2d_io ,           &
                                     & sinc2d_io , sol2d_io ,           &
                                     & solvd2d_io , solvs2d_io
!
      real(8) , pointer , dimension(:,:) :: ssw2da_io , sdeltk2d_io ,  &
                           & sdelqk2d_io , sfracv2d_io , sfracb2d_io , &
                           & sfracs2d_io , svegfrac2d_io
!      
      real(4) , allocatable , dimension(:,:,:) :: fbat_io
      real(4) , allocatable , dimension(:,:,:,:) :: fsub_io
      real(4) , allocatable , dimension(:,:,:) :: frad2d_io
      real(4) , allocatable , dimension(:,:,:,:) :: frad3d_io
      real(4) , allocatable , dimension(:,:) :: radpsa_io

      real(8) , allocatable , dimension(:,:) :: cbmf2d_io
      real(8) , allocatable , dimension(:,:,:) :: fcc_io , rsheat_io ,  &
                                      & rswat_io

      real(8) , allocatable , dimension(:,:,:,:) :: absnxt_io
      real(8) , allocatable , dimension(:,:,:,:) :: abstot_io
      real(8) , allocatable , dimension(:,:,:) :: emstot_io

      real(8) , allocatable , dimension(:,:,:) :: heatrt_io
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

      real(8) , pointer , dimension(:,:) :: cldefi_io , hfx_io , &
                                   & psa_io , psb_io ,        &
                                   & qfx_io , rainc_io , rainnc_io ,    &
                                   & tga_io , tgbb_io ,     &
                                   & tgb_io , uvdrag_io , &
                                   & zpbl_io
      real(8) , pointer , dimension(:,:,:) :: omega_io , tbase_io

      type(atmstate) :: atm1_io , atm2_io
      type(domain) :: mddom_io

      real(8) , allocatable , dimension(:,:,:) :: inisrf0
      real(8) , allocatable , dimension(:,:,:) :: inisrf_0

      real(8) , allocatable , dimension(:,:) :: var1snd , var1rcv
      integer , allocatable , dimension(:,:) :: var2d0
      integer , allocatable , dimension(:,:) :: var2d_0
      integer , allocatable , dimension(:,:,:) :: var2d1
      integer , allocatable , dimension(:,:,:) :: var2d_1
 
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
      real(8) , allocatable , dimension(:,:,:) :: dustsotex_io
#ifdef CLM
      real(8) , pointer , dimension(:,:) :: sols2d_io , soll2d_io ,     &
                   &      solsd2d_io , solld2d_io , aldifl2d_io ,       &
                   &      aldirs2d_io , aldirl2d_io , aldifs2d_io
      real(8) , allocatable , dimension(:,:) :: satbrt2d_io
#endif
      real(8) , allocatable , dimension(:,:,:,:) :: src0
      real(8) , allocatable , dimension(:,:,:,:) :: src_0
      real(8) , allocatable , dimension(:,:,:) :: src1
      real(8) , allocatable , dimension(:,:,:) :: src_1
!
      real(8) , allocatable , dimension(:,:,:) :: sav0
      real(8) , allocatable , dimension(:,:,:) :: sav_0
      real(8) , allocatable , dimension(:,:,:) :: sav0a
      real(8) , allocatable , dimension(:,:,:) :: sav_0a
      real(8) , allocatable , dimension(:,:,:) :: sav0b
      real(8) , allocatable , dimension(:,:,:) :: sav_0b
      real(8) , allocatable , dimension(:,:,:) :: sav0c
      real(8) , allocatable , dimension(:,:,:) :: sav_0c
      real(8) , allocatable , dimension(:,:,:) :: sav0s
      real(8) , allocatable , dimension(:,:,:) :: sav_0s
      real(8) , allocatable , dimension(:,:,:) :: sav0d
      real(8) , allocatable , dimension(:,:,:) :: sav_0d
      real(8) , allocatable , dimension(:,:,:) :: sav1
      real(8) , allocatable , dimension(:,:,:) :: sav_1
      real(8) , allocatable , dimension(:,:,:) :: sav2
      real(8) , allocatable , dimension(:,:,:) :: sav_2
      integer , allocatable , dimension(:,:,:) :: sav2a
      integer , allocatable , dimension(:,:,:) :: sav_2a
      real(8) , allocatable , dimension(:,:,:) :: sav4
      real(8) , allocatable , dimension(:,:,:) :: sav_4
      real(8) , allocatable , dimension(:,:,:) :: sav4a
      real(8) , allocatable , dimension(:,:,:) :: sav_4a
      real(8) , allocatable , dimension(:,:,:) :: sav6
      real(8) , allocatable , dimension(:,:,:) :: sav_6
#ifdef CLM
      real(8) , allocatable , dimension(:,:,:) :: sav_clmout
      real(8) , allocatable , dimension(:,:,:) :: sav_clmin
#endif

!---------- DATA init section--------------------------------------------

      contains 
!
!     This routines allocate all the arrays contained in the module
!
      subroutine allocate_mod_mppio(lband)
        implicit none
        logical , intent(in) :: lband
        integer :: ierr            ! control variable for allocation
        character(len=50) :: myname = 'allocate_mod_mppio'

        allocate(var1snd(kz,8),stat=ierr)
        call check_alloc(ierr,myname,'var1snd',size(var1snd))
        var1snd = d_zero
        allocate(var1rcv(kz,8),stat=ierr)
        call check_alloc(ierr,myname,'var1rcv',size(var1rcv))
        var1rcv = d_zero
        allocate(var2d0(iy,jxp),stat=ierr)
        var2d0 = -1
        allocate(var2d1(iy,nnsg,jxp),stat=ierr)
        var2d1 = -1
        allocate(inisrf0(iy,nnsg*4+7,jxp),stat=ierr)
        call check_alloc(ierr,myname,'inisrf0',size(inisrf0))
        inisrf0 = d_zero
        allocate(atm0(iy,kz*6+3+nnsg*3,jxp),stat=ierr)
        call check_alloc(ierr,myname,'atm0',size(atm0))
        atm0 = d_zero
        allocate(bat0(iym2,numbat,jxp),stat=ierr)
        call check_alloc(ierr,myname,'bat0',size(bat0))
        bat0 = 0.0
        allocate(out0(iy,3,jxp),stat=ierr)
        call check_alloc(ierr,myname,'out0',size(out0))
        out0 = d_zero
        allocate(rad0(iym2,nrad3d*kz+nrad2d,jxp),stat=ierr)
        call check_alloc(ierr,myname,'rad0',size(rad0))
        rad0 = 0.0
        allocate(sub0(iym2,nnsg,numsub,jxp),stat=ierr)
        call check_alloc(ierr,myname,'sub0',size(sub0))
        sub0 = 0.0
        if (ichem == 1) then
          allocate(chem0(iy,ntr*kz+kz*3+ntr*7+5,jxp),stat=ierr)
          call check_alloc(ierr,myname,'chem0',size(chem0))
          chem0 = d_zero
          allocate(src0(iy,mpy,ntr,jxp),stat=ierr)
          call check_alloc(ierr,myname,'src0',size(src0))
          src0 = d_zero
          allocate(src1(iy,nats,jxp),stat=ierr)
          call check_alloc(ierr,myname,'src1',size(src1))
          src1 = d_zero
        end if

        if (myid == 0) then

          call allocate_domain(mddom_io,.false.)
          call allocate_atmstate(atm1_io,.false.,0,0)
          call allocate_atmstate(atm2_io,.false.,0,0)

          allocate(inisrf_0(iy,nnsg*4+7,jx),stat=ierr)
          call check_alloc(ierr,myname,'inisrf_0',size(inisrf_0))
          inisrf_0 = d_zero
          allocate(var2d_0(iy,jx),stat=ierr)
          var2d_0 = -1
          allocate(var2d_1(iy,nnsg,jx),stat=ierr)
          var2d_1 = -1
          allocate(atm_0(iy,kz*6+3+nnsg*3,jx),stat=ierr)
          call check_alloc(ierr,myname,'atm_0',size(atm_0))
          atm_0 = d_zero
          allocate(bat_0(iym2,numbat,jx),stat=ierr)
          call check_alloc(ierr,myname,'bat_0',size(bat_0))
          bat_0 = 0.0
          allocate(out_0(iy,3,jx),stat=ierr)
          call check_alloc(ierr,myname,'out_0',size(out_0))
          out_0 = d_zero
          allocate(rad_0(iym2,nrad3d*kz+nrad2d,jx),stat=ierr)
          call check_alloc(ierr,myname,'rad_0',size(rad_0))
          rad_0 = 0.0
          allocate(sub_0(iym2,nnsg,numsub,jx),stat=ierr)
          call check_alloc(ierr,myname,'sub_0',size(sub_0))
          sub_0 = 0.0
          if (ichem == 1) then
            allocate(chem_0(iy,ntr*kz+kz*3+ntr*7+5,jx),stat=ierr)
            call check_alloc(ierr,myname,'chem_0',size(chem_0))
            chem_0 = d_zero
            allocate(src_0(iy,mpy,ntr,jx),stat=ierr)
            call check_alloc(ierr,myname,'src_0',size(src_0))
            src_0 = d_zero
            allocate(src_1(iy,nats,jx),stat=ierr)
            call check_alloc(ierr,myname,'src_1',size(src_1))
            src_1 = d_zero
          end if
          if (lband) then
            allocate(spacesubm1(nnsg,iym1,jx,20),stat=ierr)
          else
            allocate(spacesubm1(nnsg,iym1,jxm1,20),stat=ierr)
          end if
          call check_alloc(ierr,myname,'spacesubm1',size(spacesubm1))
          spacesubm1 = d_zero
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
            allocate(veg2d1_io(nnsg,iym1,jx),stat=ierr)
            allocate(ocld2d_io(nnsg,iym1,jx),stat=ierr)
            allocate(veg2d_io(iym1,jx),stat=ierr)
            allocate(ldmsk_io(iym1,jx),stat=ierr)
          else
            allocate(veg2d1_io(nnsg,iym1,jxm1),stat=ierr)
            allocate(ocld2d_io(nnsg,iym1,jxm1),stat=ierr)
            allocate(veg2d_io(iym1,jxm1),stat=ierr)
            allocate(ldmsk_io(iym1,jxm1),stat=ierr)
          end if
          veg2d1_io = -1
          ocld2d_io = -1
          veg2d_io = -1
          ldmsk_io = -1
          allocate(spacesub(nnsg,iy,jx,4),stat=ierr)
          call check_alloc(ierr,myname,'spacesub',size(spacesub))
          spacesub = d_zero
          ht1_io     => spacesub(:,:,:,1)
          satbrt1_io => spacesub(:,:,:,2)
          xlat1_io   => spacesub(:,:,:,3)
          xlon1_io   => spacesub(:,:,:,4)
          if (lakemod == 1) then
            allocate(dhlake1_io(nnsg,iy,jx),stat=ierr)
            call check_alloc(ierr,myname,'dhlake1_io',size(dhlake1_io))
            allocate(idep2d_io(nnsg,iym1,jx),stat=ierr)
            call check_alloc(ierr,myname,'idep2d_io',size(idep2d_io))
            allocate(eta2d_io(nnsg,iym1,jx),stat=ierr)
            call check_alloc(ierr,myname,'eta2d_io',size(eta2d_io))
            allocate(hi2d_io(nnsg,iym1,jx),stat=ierr)
            call check_alloc(ierr,myname,'hi2d_io',size(hi2d_io))
            allocate(aveice2d_io(nnsg,iym1,jx),stat=ierr)
            call check_alloc(ierr,myname,'aveice2d_io', size(aveice2d_io))
            allocate(hsnow2d_io(nnsg,iym1,jx),stat=ierr)
            call check_alloc(ierr,myname,'hsnow2d_io',size(hsnow2d_io))
            allocate(evl2d_io(nnsg,iym1,jx),stat=ierr)
            call check_alloc(ierr,myname,'evl2d_io',size(evl2d_io))
            allocate(tlak3d_io(ndpmax,nnsg,iym1,jx),stat=ierr)
            call check_alloc(ierr,myname,'tlak3d_io',size(tlak3d_io))
          endif
          if (lband) then
            if ( ichem == 1 ) then
              allocate(spacebat(iym1,jx,15),stat=ierr)
            else
              allocate(spacebat(iym1,jx,8),stat=ierr)
            end if
          else
            if ( ichem == 1 ) then
              allocate(spacebat(iym1,jxm1,15),stat=ierr)
            else
              allocate(spacebat(iym1,jxm1,8),stat=ierr)
            end if
          end if
          call check_alloc(ierr,myname,'spacebat',size(spacebat))
          spacebat(:,:,:) = d_zero
          flw2d_io      => spacebat(:,:,1)
          flwd2d_io     => spacebat(:,:,2)
          fsw2d_io      => spacebat(:,:,3)
          sabv2d_io     => spacebat(:,:,4)
          sinc2d_io     => spacebat(:,:,5)
          sol2d_io      => spacebat(:,:,6)
          solvd2d_io    => spacebat(:,:,7)
          solvs2d_io    => spacebat(:,:,8)
          if ( ichem == 1 ) then
            ssw2da_io     => spacebat(:,:,9)
            sdelqk2d_io   => spacebat(:,:,10)
            sdeltk2d_io   => spacebat(:,:,11)
            sfracb2d_io   => spacebat(:,:,12)
            sfracs2d_io   => spacebat(:,:,13)
            sfracv2d_io   => spacebat(:,:,14)
            svegfrac2d_io => spacebat(:,:,15)
          end if
          if (lband) then
            allocate(fbat_io(jx,iym2,numbat),stat=ierr)
            call check_alloc(ierr,myname,'fbat_io',size(fbat_io))
            allocate(fsub_io(nnsg,jx,iym2,numsub),stat=ierr) 
            call check_alloc(ierr,myname,'fsub_io',size(fsub_io))
            allocate(frad2d_io(jx,iym2,nrad2d),stat=ierr)
            call check_alloc(ierr,myname,'frad2d_io',size(frad2d_io))
            allocate(frad3d_io(jx,iym2,kz,nrad3d),stat=ierr)
            call check_alloc(ierr,myname,'frad3d_io',size(frad3d_io))
            allocate(radpsa_io(jx,iym2),stat=ierr)
            call check_alloc(ierr,myname,'radpsa_io',size(radpsa_io))
          else
            allocate(fbat_io(jxm2,iym2,numbat),stat=ierr)
            call check_alloc(ierr,myname,'fbat_io',size(fbat_io))
            allocate(fsub_io(nnsg,jxm2,iym2,numsub),stat=ierr) 
            call check_alloc(ierr,myname,'fsub_io',size(fsub_io))
            allocate(frad2d_io(jxm2,iym2,nrad2d),stat=ierr)
            call check_alloc(ierr,myname,'frad2d_io',size(frad2d_io))
            allocate(frad3d_io(jxm2,iym2,kz,nrad3d),stat=ierr)
            call check_alloc(ierr,myname,'frad3d_io',size(frad3d_io))
            allocate(radpsa_io(jxm2,iym2),stat=ierr)
            call check_alloc(ierr,myname,'radpsa_io',size(radpsa_io))
          end if
          fbat_io = 0.0
          fsub_io = 0.0
          frad2d_io = 0.0
          frad3d_io = 0.0
          radpsa_io = 0.0
          allocate(cbmf2d_io(iy,jx),stat=ierr)
          call check_alloc(ierr,myname,'cbmf2d_io',size(cbmf2d_io))
          cbmf2d_io = d_zero
          allocate(fcc_io(iy,kz,jx),stat=ierr)
          call check_alloc(ierr,myname,'fcc_io',size(fcc_io))
          fcc_io = d_zero
          allocate(rsheat_io(iy,kz,jx),stat=ierr)
          call check_alloc(ierr,myname,'rsheat_io',size(rsheat_io))
          rsheat_io = d_zero
          allocate(rswat_io(iy,kz,jx),stat=ierr)
          call check_alloc(ierr,myname,'rswat_io',size(rswat_io))
          rswat_io = d_zero
          allocate(dstor_io(iy,jx,nsplit),stat=ierr)
          call check_alloc(ierr,myname,'dstor_io',size(dstor_io))
          dstor_io = d_zero
          allocate(hstor_io(iy,jx,nsplit),stat=ierr)
          call check_alloc(ierr,myname,'hstor_io',size(hstor_io))
          hstor_io = d_zero
          if (lband) then
            allocate(absnxt_io(iym1,kz,4,jx),stat=ierr)
            call check_alloc(ierr,myname,'absnxt_io',size(absnxt_io))
            allocate(abstot_io(iym1,kzp1,kz + 1,jx),stat=ierr)
            call check_alloc(ierr,myname,'abstot_io',size(abstot_io))
            allocate(emstot_io(iym1,kzp1,jx),stat=ierr)
            call check_alloc(ierr,myname,'emstot_io',size(emstot_io))
            allocate(heatrt_io(iym1,kz,jx),stat=ierr)
            call check_alloc(ierr,myname,'heatrt_io',size(heatrt_io))
            allocate(o3prof_io(iym1,kzp1,jx),stat=ierr)
            call check_alloc(ierr,myname,'o3prof_io',size(o3prof_io))
            allocate(aerasp_io(iym1,kz,jx),stat=ierr)
            call check_alloc(ierr,myname,'aerasp_io',size(aerasp_io))
            allocate(aerext_io(iym1,kz,jx),stat=ierr)
            call check_alloc(ierr,myname,'aerext_io',size(aerext_io))
            allocate(aerssa_io(iym1,kz,jx),stat=ierr)
            call check_alloc(ierr,myname,'aerssa_io',size(aerssa_io))
            allocate(aersrrf_io(iym1,jx),stat=ierr)
            call check_alloc(ierr,myname,'aersrrf_io',size(aersrrf_io))
            allocate(aertarf_io(iym1,jx),stat=ierr)
            call check_alloc(ierr,myname,'aertarf_io',size(aertarf_io))
            allocate(aertalwrf_io(iym1,jx),stat=ierr)
            call check_alloc(ierr,myname,'aertalwrf_io',size(aertalwrf_io))
            allocate(aersrlwrf_io(iym1,jx),stat=ierr)
            call check_alloc(ierr,myname,'aersrlwrf_io',size(aersrlwrf_io))
          else
            allocate(absnxt_io(iym1,kz,4,jxm1),stat=ierr)
            call check_alloc(ierr,myname,'absnxt_io',size(absnxt_io))
            allocate(abstot_io(iym1,kzp1,kz + 1,jxm1),stat=ierr)
            call check_alloc(ierr,myname,'abstot_io',size(abstot_io))
            allocate(emstot_io(iym1,kzp1,jxm1),stat=ierr)
            call check_alloc(ierr,myname,'emstot_io',size(emstot_io))
            allocate(heatrt_io(iym1,kz,jxm1),stat=ierr)
            call check_alloc(ierr,myname,'heatrt_io',size(heatrt_io))
            allocate(o3prof_io(iym1,kzp1,jxm1),stat=ierr)
            call check_alloc(ierr,myname,'o3prof_io',size(o3prof_io))
            allocate(aerasp_io(iym1,kz,jxm1),stat=ierr)
            call check_alloc(ierr,myname,'aerasp_io',size(aerasp_io))
            allocate(aerext_io(iym1,kz,jxm1),stat=ierr)
            call check_alloc(ierr,myname,'aerext_io',size(aerext_io))
            allocate(aerssa_io(iym1,kz,jxm1),stat=ierr)
            call check_alloc(ierr,myname,'aerssa_io',size(aerssa_io))
            allocate(aersrrf_io(iym1,jxm1),stat=ierr)
            call check_alloc(ierr,myname,'aersrrf_io',size(aersrrf_io))
            allocate(aertarf_io(iym1,jxm1),stat=ierr)
            call check_alloc(ierr,myname,'aertarf_io',size(aertarf_io))
            allocate(aertalwrf_io(iym1,jxm1),stat=ierr)
            call check_alloc(ierr,myname,'aertalwrf_io',size(aertalwrf_io))
            allocate(aersrlwrf_io(iym1,jxm1),stat=ierr)
            call check_alloc(ierr,myname,'aersrlwrf_io',size(aersrlwrf_io))
          end if
          absnxt_io = d_zero
          abstot_io = d_zero
          emstot_io = d_zero
          heatrt_io = d_zero
          o3prof_io = d_zero
          aerasp_io = d_zero
          aerext_io = d_zero
          aerssa_io = d_zero
          aersrrf_io = d_zero
          aertarf_io = d_zero
          aertalwrf_io = d_zero
          aersrlwrf_io = d_zero
          if ( ichem == 1 ) then
            allocate(cemtrac_io(iy,jx,ntr),stat=ierr)
            call check_alloc(ierr,myname,'cemtrac_io',size(cemtrac_io))
            cemtrac_io = d_zero
            allocate(cemtr_io(iy,jx,ntr),stat=ierr)
            call check_alloc(ierr,myname,'cemtr_io',size(cemtr_io))
            cemtr_io = d_zero
            allocate(wxaq_io(iy,jx,ntr),stat=ierr)
            call check_alloc(ierr,myname,'wxaq_io',size(wxaq_io))
            wxaq_io = d_zero
            allocate(wxsg_io(iy,jx,ntr),stat=ierr)
            call check_alloc(ierr,myname,'wxsg_io',size(wxsg_io))
            wxsg_io = d_zero
            allocate(rxsaq1_io(iy,kz,jx,ntr),stat=ierr)
            call check_alloc(ierr,myname,'rxsaq1_io',size(rxsaq1_io))
            rxsaq1_io = d_zero
            allocate(rxsaq2_io(iy,kz,jx,ntr),stat=ierr)
            call check_alloc(ierr,myname,'rxsaq2_io',size(rxsaq2_io))
            rxsaq2_io = d_zero
            allocate(rxsg_io(iy,kz,jx,ntr),stat=ierr)
            call check_alloc(ierr,myname,'rxsg_io',size(rxsg_io))
            rxsg_io = d_zero
            allocate(remcvc_io(iy,kz,jx,ntr),stat=ierr)
            call check_alloc(ierr,myname,'remcvc_io',size(remcvc_io))
            remcvc_io = d_zero
            allocate(remlsc_io(iy,kz,jx,ntr),stat=ierr)
            call check_alloc(ierr,myname,'remlsc_io',size(remlsc_io))
            remlsc_io = d_zero
            allocate(remdrd_io(iy,jx,ntr),stat=ierr)
            call check_alloc(ierr,myname,'remdrd_io',size(remdrd_io))
            remdrd_io = d_zero
          end if
          allocate(space2d(iy,jx,4),stat=ierr)
          call check_alloc(ierr,myname,'space2d',size(space2d))
          space2d = d_zero
          ps0_io => space2d(:,:,1)
          ps1_io => space2d(:,:,2)
          ts0_io => space2d(:,:,3)
          ts1_io => space2d(:,:,4)
          allocate(space3d(iy,kz,jx,10),stat=ierr)
          call check_alloc(ierr,myname,'space3d',size(space3d))
          space3d = d_zero
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
          allocate(spacev(kz,jx,8),stat=ierr)
          call check_alloc(ierr,myname,'spacev',size(spacev))
          spacev = d_zero
          ui1_io  => spacev(:,:,1)
          ui2_io  => spacev(:,:,2)
          uilx_io => spacev(:,:,3)
          uil_io  => spacev(:,:,4)
          vi1_io  => spacev(:,:,5)
          vi2_io  => spacev(:,:,6)
          vilx_io => spacev(:,:,7)
          vil_io  => spacev(:,:,8)
          if ( ichem == 1 ) then
            allocate(chemsrc_io(iy,jx,mpy,ntr),stat=ierr)
            call check_alloc(ierr,myname,'chemsrc_io',size(chemsrc_io))
            chemsrc_io = d_zero
            allocate(ddsfc_io(iy,jx,ntr),stat=ierr)
            call check_alloc(ierr,myname,'ddsfc_io',size(ddsfc_io))
            ddsfc_io = d_zero
            allocate(dtrace_io(iy,jx,ntr),stat=ierr)
            call check_alloc(ierr,myname,'dtrace_io',size(dtrace_io))
            dtrace_io = d_zero
            allocate(wdcvc_io(iy,jx,ntr),stat=ierr)
            call check_alloc(ierr,myname,'wdcvc_io',size(wdcvc_io))
            wdcvc_io = d_zero
            allocate(wdlsc_io(iy,jx,ntr),stat=ierr)
            call check_alloc(ierr,myname,'wdlsc_io',size(wdlsc_io))
            wdlsc_io = d_zero
          end if
          allocate(dustsotex_io(iy,jx,nats),stat=ierr)
          call check_alloc(ierr,myname,'dustsotex_io',size(dustsotex_io))
          dustsotex_io = d_zero
          if (lband) then
            allocate(pptc_io(iym1,jx),stat=ierr)
            call check_alloc(ierr,myname,'pptc_io',size(pptc_io))
            allocate(pptnc_io(iym1,jx),stat=ierr)
            call check_alloc(ierr,myname,'pptnc_io',size(pptnc_io))
            allocate(prca2d_io(iym1,jx),stat=ierr)
            call check_alloc(ierr,myname,'prca2d_io',size(prca2d_io))
            allocate(prnca2d_io(iym1,jx),stat=ierr)
            call check_alloc(ierr,myname,'prnca2d_io',size(prnca2d_io))
          else
            allocate(pptc_io(iym1,jxm1),stat=ierr)
            call check_alloc(ierr,myname,'pptc_io',size(pptc_io))
            allocate(pptnc_io(iym1,jxm1),stat=ierr)
            call check_alloc(ierr,myname,'pptnc_io',size(pptnc_io))
            allocate(prca2d_io(iym1,jxm1),stat=ierr)
            call check_alloc(ierr,myname,'prca2d_io',size(prca2d_io))
            allocate(prnca2d_io(iym1,jxm1),stat=ierr)
            call check_alloc(ierr,myname,'prnca2d_io',size(prnca2d_io))
          end if
          pptc_io = d_zero
          pptnc_io = d_zero
          prca2d_io = d_zero
          prnca2d_io = d_zero
          if ( ichem == 1 ) then
            allocate(chia_io(iy,kz,jx,ntr),stat=ierr)
            call check_alloc(ierr,myname,'chia_io',size(chia_io))
            chia_io = d_zero
            allocate(chib_io(iy,kz,jx,ntr),stat=ierr)
            call check_alloc(ierr,myname,'chib_io',size(chib_io))
            chib_io = d_zero
          end if
          if (icup == 3) then
            allocate(spacesurf(iy,jx,11),stat=ierr)
          else
            allocate(spacesurf(iy,jx,12),stat=ierr)
          end if
          call check_alloc(ierr,myname,'spacesurf',size(spacesurf))
          spacesurf = d_zero
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
          allocate(omega_io(iy,kz,jx),stat=ierr)
          call check_alloc(ierr,myname,'omega_io',size(omega_io))
          omega_io = d_zero
          if (icup == 3) then
            allocate(tbase_io(iy,kz,jx),stat=ierr)
            call check_alloc(ierr,myname,'tbase_io',size(tbase_io))
            tbase_io = d_zero
          end if
#ifdef CLM
          if (lband) then
            allocate(spaceclm(iym1,jx,8))
          else
            allocate(spaceclm(iym1,jxm1,8))
          end if
          call check_alloc(ierr,myname,'spaceclm',size(spaceclm))
          spaceclm = d_zero
          sols2d_io   => spaceclm(:,:,1)
          soll2d_io   => spaceclm(:,:,2)
          solsd2d_io  => spaceclm(:,:,3)
          solld2d_io  => spaceclm(:,:,4)
          aldifl2d_io => spaceclm(:,:,5)
          aldirs2d_io => spaceclm(:,:,6)
          aldirl2d_io => spaceclm(:,:,7)
          aldifs2d_io => spaceclm(:,:,8)
          allocate(satbrt2d_io(iy,jx))
          call check_alloc(ierr,myname,'satbrt2d_io',size(satbrt2d_io))
          satbrt2d_io = d_zero
#endif
        endif
        if (myid == 0) then
          allocate(sav_0(iy,kz*4+2,jx),stat=ierr)
          call check_alloc(ierr,myname,'sav_0',size(sav_0))
          sav_0 = d_zero
          allocate(sav_0a(iy,kzp1+4,jx) ,stat=ierr)
          call check_alloc(ierr,myname,'sav_0a',size(sav_0a))
          sav_0a = d_zero
          allocate(sav_0b(iy,kzp1,jx),stat=ierr)
          call check_alloc(ierr,myname,'sav_0b',size(sav_0b))
          sav_0b = d_zero
          allocate(sav_0c(iy,kz*2,jx),stat=ierr)
          call check_alloc(ierr,myname,'sav_0c',size(sav_0c))
          sav_0c = d_zero
          allocate(sav_0s(iy,kz,jx),stat=ierr)
          call check_alloc(ierr,myname,'sav_0s',size(sav_0s))
          sav_0s = d_zero
          allocate(sav_0d(iy,nsplit*2,jx),stat=ierr)
          call check_alloc(ierr,myname,'sav_0d',size(sav_0d))
          sav_0d = d_zero
          allocate(sav_1(iym1,kz*4+(kzp1*kzp2),jx),stat=ierr)
          call check_alloc(ierr,myname,'sav_1',size(sav_1))
          sav_1 = d_zero
          allocate(sav_2(iym1,nnsg*5+4,jx),stat=ierr)
          call check_alloc(ierr,myname,'sav_2',size(sav_2))
          sav_2 = d_zero
          allocate(sav_2a(iym1,nnsg*2+2,jx),stat=ierr)
          sav_2a = -1
          if ( ichem == 1 ) then
            allocate(sav_4(iy,ntr*(kz*4+1),jx),stat=ierr)
            call check_alloc(ierr,myname,'sav_4',size(sav_4))
            sav_4 = d_zero
            allocate(sav_4a(iym1,7,jx),stat=ierr)
            call check_alloc(ierr,myname,'sav_4a',size(sav_4a))
            sav_4a = d_zero
          end if
          allocate(sav_6(kz,8,jx),stat=ierr)
          call check_alloc(ierr,myname,'sav_6',size(sav_6))
          sav_6 = d_zero
#ifdef CLM
          allocate(sav_clmout(iym1,8,jx),stat=ierr)
          call check_alloc(ierr,myname,'sav_clmout',size(sav_clmout))
          sav_clmout = d_zero
#endif
        end if
        allocate(sav0(iy,kz*4+2,jxp),stat=ierr)
        call check_alloc(ierr,myname,'sav0',size(sav0))
        sav0 = d_zero
        allocate(sav0a(iy,kzp1+4,jxp) ,stat=ierr)
        call check_alloc(ierr,myname,'sav0a',size(sav0a))
        sav0a = d_zero
        allocate(sav0b(iy,kzp1,jxp),stat=ierr)
        call check_alloc(ierr,myname,'sav0b',size(sav0b))
        sav0b = d_zero
        allocate(sav0c(iy,kz*2,jxp),stat=ierr)
        call check_alloc(ierr,myname,'sav0c',size(sav0c))
        sav0c = d_zero
        allocate(sav0s(iy,kz,jxp),stat=ierr)
        call check_alloc(ierr,myname,'sav0s',size(sav0s))
        sav0s = d_zero
        allocate(sav0d(iy,nsplit*2,jxp),stat=ierr)
        call check_alloc(ierr,myname,'sav0d',size(sav0d))
        sav0d = d_zero
        allocate(sav1(iym1,kz*4+(kzp1*kzp2),jxp),stat=ierr)
        call check_alloc(ierr,myname,'sav1',size(sav1))
        sav1 = d_zero
        allocate(sav2(iym1,nnsg*5+4,jxp),stat=ierr)
        call check_alloc(ierr,myname,'sav2',size(sav2))
        sav2 = d_zero
        allocate(sav2a(iym1,nnsg*2+2,jxp),stat=ierr)
        sav2a = -1
        if ( ichem == 1 ) then
          allocate(sav4(iy,ntr*(kz*4+1),jxp),stat=ierr)
          call check_alloc(ierr,myname,'sav4',size(sav4))
          sav4 = d_zero
          allocate(sav4a(iym1,7,jxp),stat=ierr)
          call check_alloc(ierr,myname,'sav4a',size(sav4a))
          sav4a = d_zero
        end if
        allocate(sav6(kz,8,jxp),stat=ierr)
        call check_alloc(ierr,myname,'sav6',size(sav6))
        sav6 = d_zero
#ifdef CLM
        allocate(sav_clmin(iym1,8,jxp),stat=ierr)
        call check_alloc(ierr,myname,'sav_clmin',size(sav_clmin))
        sav_clmin = d_zero
#endif

        write(aline,*) 'allocate_mod_mppio'
        call say
        call report_alloc('allocate_mod_mppio')
        write(aline,*) 'allocate_mod_mppio'
        call say
      end subroutine allocate_mod_mppio
!
      subroutine free_mpp_initspace
        implicit none
        deallocate(inisrf0)
        if (ichem == 1) then
          deallocate(src0)
          deallocate(src1)
        end if
        if (myid == 0) then
          deallocate(inisrf_0)
          if (ichem == 1) then
            deallocate(src_0)
            deallocate(src_1)
          end if
        end if
      end subroutine free_mpp_initspace
!
      end module mod_mppio
