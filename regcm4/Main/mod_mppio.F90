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

      use mod_regcm_param

      implicit none

#ifdef MPP1
!
      real(8) , dimension(nnsg,iym1,jxm1) :: col2d_io , dew2d_io ,      &
           & evpa2d_io , gwet2d_io , ircp2d_io , ocld2d_io , rno2d_io , &
           & rnos2d_io , sag2d_io , scv2d_io , sena2d_io , sice2d_io ,  &
           & srw2d_io , ssw2d_io , swt2d_io , taf2d_io , text2d_io ,    &
           & tg2d_io , tgb2d_io , tlef2d_io , veg2d1_io
      real(8) , dimension(nnsg,iy,jx) :: ht1_io , satbrt1_io
!
      real(8) , dimension(iym1,jxm1) :: flw2d_io , flwd2d_io ,          &
                          & fsw2d_io , sabv2d_io , sdelqk2d_io ,        &
                                     & sdeltk2d_io , sfracb2d_io ,      &
                                     & sfracs2d_io , sfracv2d_io ,      &
                                     & sinc2d_io , sol2d_io ,           &
                                     & solvd2d_io , solvs2d_io ,        &
                                     & ssw2da_io , svegfrac2d_io ,      &
                                     & veg2d_io

      
      real(kind=4) , dimension(jxm2,iym2,numbat) :: fbat_io
      real(4) , dimension(nnsg,jxm2,iym2,numsub) :: fsub_io

      real(4) , dimension(jxm2,iym2,nrad2d) :: frad2d_io
      real(4) , dimension(jxm2,iym2,kz,nrad3d) :: frad3d_io

      real(8) , dimension(iy,jx) :: cbmf2d_io
      real(8) , dimension(iy,kz,jx) :: fcc_io , rsheat_io , rswat_io

      real(8) , dimension(iym1,kz,4,jxm1) :: absnxt_io
      real(8) , dimension(iym1,kzp1,kz + 1,jxm1) :: abstot_io
      real(8) , dimension(iym1,kzp1,jxm1) :: emstot_io

      real(8) , dimension(iym1,kz,jxm1) :: heatrt_io
      real(8) , dimension(iym1,kzp1,jxm1) :: o3prof_io

      real(8) , dimension(iy,jx,nsplit) :: dstor_io , hstor_io

      real(8) , dimension(iym1,kz,jxm1) :: aerasp_io , aerext_io ,      &
                                   & aerssa_io
      real(8) , dimension(iym1,jxm1) :: aersrrf_io , aertarf_io
      real(8) , dimension(iy,jx,ntr) :: cemtrac_io , cemtr_io ,         &
                                   & wxaq_io , wxsg_io
      real(8) , dimension(iy,jx) :: dustsotex_io
      real(8) , dimension(iy,kz,jx,ntr) :: rxsaq1_io , rxsaq2_io ,      &
                                   & rxsg_io

      real(8) , dimension(iy,kz,jx,ntr) :: remcvc_io , remlsc_io
      real(8) , dimension(iy,jx,ntr) :: remdrd_io

      real(8) , dimension(iy,jx) :: ps0_io , ps1_io , ts0_io , ts1_io
      real(8) , dimension(iy,kz,jx) :: qb0_io , qb1_io , so0_io ,       &
                                      & so1_io , tb0_io , tb1_io ,      &
                                      & ub0_io , ub1_io , vb0_io ,      &
                                      & vb1_io
      real(8) , dimension(kz,jx) :: ui1_io , ui2_io , uilx_io ,         &
                                   & uil_io , vi1_io , vi2_io ,         &
                                   & vilx_io , vil_io

      real(8) , dimension(iy,jx,12,ntr) :: chemsrc_io
      real(8) , dimension(iy,jx,ntr) :: ddsfc_io , dtrace_io ,          &
                                   & wdcvc_io , wdlsc_io
      real(8) , dimension(iym1,jxm1) :: pptc_io , pptnc_io , prca2d_io ,&
                                    & prnca2d_io

      real(8) , dimension(iy,kz,jx,ntr) :: chia_io , chib_io
      real(8) , dimension(iy,jx) :: cldefi_io , f_io , hfx_io ,         &
                                   & htsd_io , ht_io , msfd_io ,        &
                                   & msfx_io , psa_io , psb_io ,        &
                                   & qfx_io , rainc_io , rainnc_io ,    &
                                   & satbrt_io , tga_io , tgbb_io ,     &
                                   & tgb_io , uvdrag_io , xlat_io ,     &
                                   & xlong_io , zpbl_io
      real(8) , dimension(iy,kz,jx) :: omega_io , qca_io , qcb_io ,     &
                                      & qva_io , qvb_io , ta_io ,       &
                                      & tbase_io , tb_io , ua_io ,      &
                                      & ub_io , va_io , vb_io
      real(8) , dimension(nnsg,iy,jx) :: snowc_io

      real(8) ,allocatable, dimension(:,:,:) :: inisrf0
      real(8) , dimension(iy,nnsg*3+8,jx) :: inisrf_0

      real(8) , dimension(kz,8) :: var1snd , var1rcv

#ifdef CLM
      real(8) , dimension(iym1,jxm1) :: sols2d_io , soll2d_io ,         &
                   &      solsd2d_io , solld2d_io , aldifl2d_io ,       &
                   &      aldirs2d_io , aldirl2d_io , aldifs2d_io ,     &
                   &      coszrs2d_io
#endif

#endif

contains 

	subroutine allocate_mppio
!
!     This routines allocate all the arrays contained in the module
!	
#ifdef MPP1
       allocate (inisrf0(iy,nnsg*3+8,jxp)) 	
#endif       
       end subroutine allocate_mppio
      
      end module mod_mppio
