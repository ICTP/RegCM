!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_bats

      use mod_regcm_param
      use mod_bats_param

      implicit none
!
      real(8) , dimension(nnsg,iym1) :: delq1d , delt1d , drag1d ,      &
           & emiss_1d , evpr1d , gwet1d , ircp1d , ldew1d , ldoc1d ,    &
           & p1d , pbp1d , prcp1d , q2m_1d , qg1d , qs1d , resp1d ,     &
           & rhs1d , rno1d , rnos1d , rsw1d , sag1d , scv1d , sent1d ,  &
           & sice1d , ssw1d , t2m_1d , taf1d , tg1d , tgb1d , tlef1d ,  &
           & ts1d , tsw1d , u10m1d , v10m1d , veg1d , z1d
      real(8) , dimension(iym1) :: flw1d , fsw1d , us1d , vs1d
!
      real(8) , dimension(nnsg,iym1) :: bfc , bsw , evmx0 , fdry ,      &
           & fwet , gwmx0 , gwmx1 , gwmx2 , porsl , relfc , rnet ,      &
           & texrat , vegt , wiltr , wt , xkmx
      real(8) , dimension(iym1) :: czen , sola , vpdd
      real(8) :: difrat
!
      integer :: ilat , ihis , mhis , ncase
!
      real(8) , dimension(nnsg,iym1) :: aarea , cdr , cdrn , cdrx ,     &
           & cf , cgrnd , cgrndl , cgrnds , clead , densi , efpr , eg , &
           & etr , etrrun , evaps , evapw , fevpg , flnet , flneto ,    &
           & fseng , htvp , ps , pw , qice , qsatl , rhosw , ribd ,     &
           & rlai , rpp , scrat , scvk , sdrop , seasb , sigf , sm ,    &
           & tm , uaf , vspda , wata , watr , watt , watu , wta , xlai ,&
           & xlsai , xrun , z1 , z1log
      real(8) , dimension(iym1) :: ems
      integer , dimension(nnsg,iym1) :: imelt , lveg
!
      real(8) , dimension(nnsg,iym1) :: cn1 , df , rgr , wta0 , wtaq0 , &
           & wtg , wtg0 , wtg2 , wtga , wtgaq , wtgl , wtglq , wtgq ,   &
           & wtgq0 , wtl0 , wtlh , wtlq , wtlq0 , wtshi , wtsqi
!
      real(8) , dimension(iym1) :: albdif , albdir , albvl , albvld ,   &
                                & albvs , albvsd , aldifl , aldifs ,    &
                                & aldirl , aldirs , emiss1d , fracd ,   &
                                & sabveg , solis , solvd , solvs
      real(8) , dimension(iy) :: coszrs
!
#ifdef MPP1
      real(8) , dimension(iym1,jxp) :: flw2d , flwa2d , flwd2d ,        &
                                    & flwda2d , fsw2d , fswa2d , pptc , &
                                    & pptnc , prca2d , prnca2d ,        &
                                    & sabv2d , sdelqk2d , sdeltk2d ,    &
                                    & sfracb2d , sfracs2d , sfracv2d ,  &
                                    & sina2d , sinc2d , sol2d ,         &
                                    & solvd2d , solvs2d , ssw2da ,      &
                                    & svegfrac2d , svga2d , veg2d
#else
      real(8) , dimension(iym1,jxm1) :: flw2d , flwa2d , flwd2d ,       &
                                    & flwda2d , fsw2d , fswa2d , pptc , &
                                    & pptnc , prca2d , prnca2d ,        &
                                    & sabv2d , sdelqk2d , sdeltk2d ,    &
                                    & sfracb2d , sfracs2d , sfracv2d ,  &
                                    & sina2d , sinc2d , sol2d ,         &
                                    & solvd2d , solvs2d , ssw2da ,      &
                                    & svegfrac2d , svga2d , veg2d
#endif
!
#ifdef MPP1
      real(8) , dimension(nnsg,iym1,jxp) :: col2d , dew2d , emiss2d ,   &
           & evpa2d , gwet2d , ircp2d , ocld2d , rno2d , rnos2d ,       &
           & sag2d , scv2d , sena2d , sice2d , srw2d , ssw2d , swt2d ,  &
           & taf2d , text2d , tg2d , tgb2d , tlef2d , veg2d1
      real(8) , dimension(nnsg,iy,jxp) :: ht1 , satbrt1
#else
      real(8) , dimension(nnsg,iym1,jxm1) :: col2d , dew2d , emiss2d ,  &
           & evpa2d , gwet2d , ircp2d , ocld2d , rno2d , rnos2d ,       &
           & sag2d , scv2d , sena2d , sice2d , srw2d , ssw2d , swt2d ,  &
           & taf2d , text2d , tg2d , tgb2d , tlef2d , veg2d1
      real(8) , dimension(nnsg,iy,jx) :: ht1 , satbrt1
#endif
!
#ifdef MPP1
      real(4) , pointer , dimension(:,:) :: drag_o , evpa_o , flwa_o ,  &
                                     & flwd_o , fswa_o , prcv_o ,       &
                                     & psmn_o , ps_o , q2m_o , rnos_o , &
                                     & rsw_o , scv_o , sena_o , sina_o ,&
                                     & ssw_o , t2mn_o , t2mx_o , t2m_o ,&
                                     & tgmn_o , tgmx_o , tg_o , tlef_o ,&
                                     & tpr_o , u10m_o , v10m_o ,        &
                                     & w10x_o , zpbl_o
#else
      real(4) , pointer , dimension(:,:) :: drag_o , evpa_o , flwa_o ,  &
                                      & flwd_o , fswa_o , prcv_o ,      &
                                      & psmn_o , ps_o , q2m_o , rnos_o ,&
                                      & rsw_o , scv_o , sena_o ,        &
                                      & sina_o , ssw_o , t2mn_o ,       &
                                      & t2mx_o , t2m_o , tgmn_o ,       &
                                      & tgmx_o , tg_o , tlef_o , tpr_o ,&
                                      & u10m_o , v10m_o , w10x_o ,      &
                                      & zpbl_o
#endif
!
#ifdef MPP1
      real(4) , pointer , dimension(:,:,:) :: drag_s , evpa_s , prcv_s ,&
           & ps_s , q2m_s , rnos_s , rsw_s , scv_s , sena_s , ssw_s ,   &
           & t2m_s , tg_s , tlef_s , tpr_s , u10m_s , v10m_s
#else
      real(4) , pointer , dimension(:,:,:) :: drag_s , evpa_s , prcv_s ,&
           & ps_s , q2m_s , rnos_s , rsw_s , scv_s , sena_s , ssw_s ,   &
           & t2m_s , tg_s , tlef_s , tpr_s , u10m_s , v10m_s
#endif
!
#ifdef MPP1
      real(4) , target , dimension(jxp,iym2,numbat) :: fbat
#else
      real(4) , target , dimension(jxm2,iym2,numbat) :: fbat
#endif

#ifdef MPP1
      real(4) , target , dimension(nnsg,jxp,iym2,numsub) :: fsub
#else
      real(4) , target , dimension(nnsg,jxm2,iym2,numsub) :: fsub
#endif
!
      real(8) , dimension(nnsg,iym1) :: p1d0 , qs1d0 , ts1d0
!
#ifdef MPP1
#ifdef CLM
      ! Direct solar rad incident on surface (<0.7)
      real(8) , dimension(iym1,jxp):: sols2d
      ! Direct solar rad incident on surface (>=0.7)
      real(8) , dimension(iym1,jxp):: soll2d
      ! Diffuse solar rad incident on surface (<0.7)
      real(8) , dimension(iym1,jxp):: solsd2d
      ! Diffuse solar rad incident on surface (>=0.7)
      real(8) , dimension(iym1,jxp):: solld2d
      real(8) , dimension(iym1,jxp):: aldirs2d
      real(8) , dimension(iym1,jxp):: aldirl2d
      real(8) , dimension(iym1,jxp):: aldifs2d
      real(8) , dimension(iym1,jxp):: aldifl2d
      real(8) , dimension(iym1,jxp):: coszrs2d
      real(8) , dimension(iym1,jxp):: rs2d
      real(8) , dimension(iym1,jxp):: ra2d
      real(8) , dimension(iym1,jxp):: q2d   ! 2 meter specific humidity
#endif
#endif

#ifdef DCSST
      ! dtskin is difference between skin tem and bulk sst
#ifdef MPP1
      real(8) , dimension(iy,jxp) :: deltas , tdeltas , dtskin
      logical , dimension(iy,jxp) :: firstcall
#else
      real(8) , dimension(iy,jx) :: deltas , tdeltas , dtskin
      logical , dimension(iy,jx) :: firstcall
#endif

#endif

      contains

#ifdef DCSST
        subroutine inidcsst
          implicit none
          firstcall(:,:) = .false.
        end subroutine inidcsst
#endif

        subroutine inibat
          implicit none
          u10m_o => fbat(:,:,1)
          v10m_o => fbat(:,:,2)
          drag_o => fbat(:,:,3)
          tg_o   => fbat(:,:,4)
          tlef_o => fbat(:,:,5)
          t2m_o  => fbat(:,:,6)
          q2m_o  => fbat(:,:,7)
          ssw_o  => fbat(:,:,8)
          rsw_o  => fbat(:,:,9)
          tpr_o  => fbat(:,:,10)
          evpa_o => fbat(:,:,11)
          rnos_o => fbat(:,:,12)
          scv_o  => fbat(:,:,13)
          sena_o => fbat(:,:,14)
          flwa_o => fbat(:,:,15)
          fswa_o => fbat(:,:,16)
          flwd_o => fbat(:,:,17)
          sina_o => fbat(:,:,18)
          prcv_o => fbat(:,:,19)
          ps_o   => fbat(:,:,20)
          zpbl_o => fbat(:,:,21)
          tgmx_o => fbat(:,:,22)
          tgmn_o => fbat(:,:,23)
          t2mx_o => fbat(:,:,24)
          t2mn_o => fbat(:,:,25)
          w10x_o => fbat(:,:,26)
          psmn_o => fbat(:,:,27)
        end subroutine inibat

        subroutine inisub
        implicit none
          u10m_s => fsub(:,:,:,1)
          v10m_s => fsub(:,:,:,2)
          drag_s => fsub(:,:,:,3)
          tg_s   => fsub(:,:,:,4)
          tlef_s => fsub(:,:,:,5)
          t2m_s  => fsub(:,:,:,6)
          q2m_s  => fsub(:,:,:,7)
          ssw_s  => fsub(:,:,:,8)
          rsw_s  => fsub(:,:,:,9)
          tpr_s  => fsub(:,:,:,10)
          evpa_s => fsub(:,:,:,11)
          rnos_s => fsub(:,:,:,12)
          scv_s  => fsub(:,:,:,13)
          sena_s => fsub(:,:,:,14)
          prcv_s => fsub(:,:,:,15)
          ps_s   => fsub(:,:,:,16)
        end subroutine inisub

      end module mod_bats
