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

      implicit none
!
      real(8) , dimension(nnsg,ixm1) :: delq1d , delt1d , drag1d ,      &
           & emiss_1d , evpr1d , gwet1d , ircp1d , ldew1d , ldoc1d ,    &
           & p1d , pbp1d , prcp1d , q2m_1d , qg1d , qs1d , resp1d ,     &
           & rhs1d , rno1d , rnos1d , rsw1d , sag1d , scv1d , sent1d ,  &
           & sice1d , ssw1d , t2m_1d , taf1d , tg1d , tgb1d , tlef1d ,  &
           & ts1d , tsw1d , u10m1d , v10m1d , veg1d , z1d
      real(8) , dimension(ixm1) :: flw1d , fsw1d , us1d , vs1d
!
      real(8) , dimension(nnsg,ixm1) :: bfc , bsw , evmx0 , fdry ,      &
           & fwet , gwmx0 , gwmx1 , gwmx2 , porsl , relfc , rnet ,      &
           & texrat , vegt , wiltr , wt , xkmx
      real(8) , dimension(ixm1) :: czen , sola , vpdd
      real(8) :: difrat
!
      real(8) , dimension(20) :: albvgl , albvgs , crough , deprv ,     &
                               & deptv , depuv , displa , fc , freza ,  &
                               & frezu , rough , rsmin , sai , seasf ,  &
                               & sqrtdi , vegc , xla , xlai0
      real(8) , dimension(12) :: bee , skrat , xmofc , xmohyd , xmopor ,&
                               & xmosuc , xmowil
      integer , dimension(20) :: iexsol , kolsol
      integer :: ihis , lat , mhis , ncase
!
      real(8) , dimension(nnsg,ixm1) :: aarea , cdr , cdrn , cdrx ,     &
           & cf , cgrnd , cgrndl , cgrnds , clead , densi , efpr , eg , &
           & etr , etrrun , evaps , evapw , fevpg , flnet , flneto ,    &
           & fseng , htvp , ps , pw , qice , qsatl , rhosw , ribd ,     &
           & rlai , rpp , scrat , scvk , sdrop , seasb , sigf , sm ,    &
           & tm , uaf , vspda , wata , watr , watt , watu , wta , xlai ,&
           & xlsai , xrun , z1 , z1log
      real(8) , dimension(ixm1) :: ems
      integer , dimension(nnsg,ixm1) :: imelt , lveg
!
      real(8) , dimension(nnsg,ixm1) :: cn1 , df , rgr , wta0 , wtaq0 , &
           & wtg , wtg0 , wtg2 , wtga , wtgaq , wtgl , wtglq , wtgq ,   &
           & wtgq0 , wtl0 , wtlh , wtlq , wtlq0 , wtshi , wtsqi
!
      real(8) , dimension(ixm1) :: albdif , albdir , albvl , albvld ,   &
                                & albvs , albvsd , aldifl , aldifs ,    &
                                & aldirl , aldirs , emiss1d , fracd ,   &
                                & sabveg , solis , solvd , solvs
      real(8) , dimension(ix) :: coszrs
!
#ifdef MPP1
      real(8) , dimension(ixm1,jxp) :: flw2d , flwa2d , flwd2d ,        &
                                    & flwda2d , fsw2d , fswa2d , pptc , &
                                    & pptnc , prca2d , prnca2d ,        &
                                    & sabv2d , sdelqk2d , sdeltk2d ,    &
                                    & sfracb2d , sfracs2d , sfracv2d ,  &
                                    & sina2d , sinc2d , sol2d ,         &
                                    & solvd2d , solvs2d , ssw2da ,      &
                                    & svegfrac2d , svga2d , veg2d
#else
      real(8) , dimension(ixm1,jxm1) :: flw2d , flwa2d , flwd2d ,       &
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
      real(8) , dimension(nnsg,ixm1,jxp) :: col2d , dew2d , emiss2d ,   &
           & evpa2d , gwet2d , ircp2d , ocld2d , rno2d , rnos2d ,       &
           & sag2d , scv2d , sena2d , sice2d , srw2d , ssw2d , swt2d ,  &
           & taf2d , text2d , tg2d , tgb2d , tlef2d , veg2d1
      real(8) , dimension(nnsg,ix,jxp) :: ht1 , satbrt1
#else
      real(8) , dimension(nnsg,ixm1,jxm1) :: col2d , dew2d , emiss2d ,  &
           & evpa2d , gwet2d , ircp2d , ocld2d , rno2d , rnos2d ,       &
           & sag2d , scv2d , sena2d , sice2d , srw2d , ssw2d , swt2d ,  &
           & taf2d , text2d , tg2d , tgb2d , tlef2d , veg2d1
      real(8) , dimension(nnsg,ix,jx) :: ht1 , satbrt1
#endif
!
#ifdef MPP1
      real(4) , dimension(jxp,ixm2) :: drag_o , evpa_o , flwa_o ,       &
                                     & flwd_o , fswa_o , prcv_o ,       &
                                     & psmn_o , ps_o , q2m_o , rnos_o , &
                                     & rsw_o , scv_o , sena_o , sina_o ,&
                                     & ssw_o , t2mn_o , t2mx_o , t2m_o ,&
                                     & tgmn_o , tgmx_o , tg_o , tlef_o ,&
                                     & tpr_o , u10m_o , v10m_o ,        &
                                     & w10x_o , zpbl_o
#else
      real(4) , dimension(jxm2,ixm2) :: drag_o , evpa_o , flwa_o ,      &
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
      real(4) , dimension(nnsg,jxp,ixm2) :: drag_s , evpa_s , prcv_s ,  &
           & ps_s , q2m_s , rnos_s , rsw_s , scv_s , sena_s , ssw_s ,   &
           & t2m_s , tg_s , tlef_s , tpr_s , u10m_s , v10m_s
#else
      real(4) , dimension(nnsg,jxm2,ixm2) :: drag_s , evpa_s , prcv_s , &
           & ps_s , q2m_s , rnos_s , rsw_s , scv_s , sena_s , ssw_s ,   &
           & t2m_s , tg_s , tlef_s , tpr_s , u10m_s , v10m_s
#endif
!
#ifdef MPP1
      real(4) , dimension(jxp,ixm2,numbat) :: fbat
#else
      real(4) , dimension(jxm2,ixm2,numbat) :: fbat
#endif

#ifdef MPP1
      real(4) , dimension(nnsg,jxp,ixm2,numsub) :: fsub
#else
      real(4) , dimension(nnsg,jxm2,ixm2,numsub) :: fsub
#endif
!
      real(8) , dimension(nnsg,ixm1) :: p1d0 , qs1d0 , ts1d0
!
!clm35
#ifdef MPP1
      ! Direct solar rad incident on surface (<0.7)
      real(8) , dimension(ixm1,jxp):: sols2d
      ! Direct solar rad incident on surface (>=0.7)
      real(8) , dimension(ixm1,jxp):: soll2d
      ! Diffuse solar rad incident on surface (<0.7)
      real(8) , dimension(ixm1,jxp):: solsd2d
      ! Diffuse solar rad incident on surface (>=0.7)
      real(8) , dimension(ixm1,jxp):: solld2d
      real(8) , dimension(ixm1,jxp):: aldirs2d
      real(8) , dimension(ixm1,jxp):: aldirl2d
      real(8) , dimension(ixm1,jxp):: aldifs2d
      real(8) , dimension(ixm1,jxp):: aldifl2d
      real(8) , dimension(ixm1,jxp):: coszrs2d
      real(8) , dimension(ixm1,jxp):: rs2d
      real(8) , dimension(ixm1,jxp):: ra2d
      real(8) , dimension(ixm1,jxp):: q2d   ! 2 meter specific humidity
#else
      ! Direct solar rad incident on surface (<0.7)
      real(8) , dimension(ixm1,jxm1):: sols2d
      ! Direct solar rad incident on surface (>=0.7)
      real(8) , dimension(ixm1,jxm1):: soll2d
      ! Diffuse solar rad incident on surface (<0.7)
      real(8) , dimension(ixm1,jxm1):: solsd2d
      ! Diffuse solar rad incident on surface (>=0.7)
      real(8) , dimension(ixm1,jxm1):: solld2d
      real(8) , dimension(ixm1,jxm1):: aldirs2d
      real(8) , dimension(ixm1,jxm1):: aldirl2d
      real(8) , dimension(ixm1,jxm1):: aldifs2d
      real(8) , dimension(ixm1,jxm1):: aldifl2d
      real(8) , dimension(ixm1,jxm1):: coszrs2d
      real(8) , dimension(ixm1,jxm1):: rs2d
      real(8) , dimension(ixm1,jxm1):: ra2d
      real(8) , dimension(ixm1,jxm1):: q2d ! 2 meter specific humidity
#endif
!clm35

!
!------------------ DATA SECTION ----------------------------------------
!
!*    vegc is maximum fractional cover of vegetation
      data vegc /0.85D0 , 0.8D0 , 0.8D0 , 0.8D0 , 0.8D0  , 0.9D0 ,      &
              &  0.8D0  , 0.0D0 , 0.6D0 , 0.8D0 , 0.35D0 , 0.0D0 ,      &
              &  0.8D0  , 2*0.0D0 , 5*0.8D0 /
!*    seasf is the difference between vegc and fractional cover at 269k
      data seasf /0.6D0 , 0.1D0 , 0.1D0 , 0.3D0 , 0.3D0 , 0.5D0 ,       &
              &   0.3D0 , 0.0D0 , 0.2D0 , 0.6D0 , 0.1D0 , 0.0D0 ,       &
              &   0.4D0 , 0.0D0 , 0.0D0 , 0.2D0 , 0.3D0 , 0.2D0 ,       &
              &   2*0.4D0 /
!*    rough is an aerodynamic roughness length (m) =approx 0.1*veg
!*    height also used snow masking depth in subrout albedo
      data rough /0.08D0 , 0.05D0 , 2*1.0D0 , 0.8D0 , 2.0D0  , 0.1D0  , &
              &   0.05D0 , 0.04D0 , 0.06D0 ,  0.1D0 , 0.01D0 , 0.03D0 , &
              &   2*0.0004D0 , 2*0.1D0 , 0.8D0 , 2*0.3D0 /
!     ******      displacement height (meter)
!     ******      if great parts of veg. are covered by snow, use
!     displa=0 ******      because mod_then the new displa-theory is not
!     valid
      data displa/0.0D0 , 0.0D0 , 9.0D0 , 9.0D0 , 0.0D0 , 18.0D0 ,      &
              &   14*0.0D0/
!     ******      min stomatl resistance (s/m)
!cc   data rsmin/153.0,4*200.0,150.0,14*200.0/   ! shuttleworth
!     data rsmin/120.0,4*200.0,150.0,14*200.0/   ! bats1e numbers
!Sara
!     data rsmin/45.,60.,2*80.,120.,2*60.,200.,80.,45.,150.,200.,45.
!     &          ,2*200.,80.,120.,100.,2*120./
      data rsmin/20*200.0D0/
!Sara_
!     ******      max leaf area index (ratio unit cover per unit ground)
!     ORIGINAL
!     data xla/6.,2.,5*6.,0.,3*6.,0.,6.,2*0.,5*6./
!     Laura 21/04/08
      data xla/4.0D0 , 2.0D0 , 4*6.0D0 , 3.0D0 , 0.0D0 , 2.0D0 , 4.0D0 ,&
            &  1.0D0 , 0.0D0 , 4.0D0 , 2*0.0D0 , 4.0D0 , 4.0D0 , 5.0D0 ,&
            &  4.0D0 , 1.0D0 /
!     ******      min leaf area index **lai depends on temp as veg cover
!     ORIGINAL
!     data xlai0/0.5,0.5,5.,2*1.,5.,0.5,0.,3*0.5,0.,0.5,2*0.,5.,1.,3.,
!     &2*0.5/
!     Laura 21/04/08
      data xlai0/2*0.5D0 , 5.0D0 , 2*1.0D0 , 5.0D0 , 1.0D0 , 0.0D0 ,    &
            &    0.5D0 , 2.0D0 , 0.5D0 , 0.0D0 , 2.0D0 , 2*0.0D0 ,      &
            &    3.0D0 , 1.0D0 , 3.0D0 , 0.5D0 , 1.0D0/

!     ******      stem area index (projected area of non-transpiring
!     sfcs)
      data sai/0.5D0 , 4.0D0 , 5*2.0D0 , 2*0.5D0 , 11*2.0D0/
!     ******      inverse square root of leaf dimension - used for
!     ******      calculating fluxes from foliage
      data sqrtdi/10.0D0 , 19*5.0D0/
!     ******      fc = light dependence of stomatal resistance
      data fc/0.02D0 , 0.02D0 , 4*0.06D0 , 11*0.02D0 , 0.06D0 ,         &
          &   2*0.02D0/

!     ******      depuv is depth of upper soil layer (mm)
!     ******      deprv is depth of root zone (mm)
!     ******      deptv is depth of total soil (mm)
      data depuv/20*100.0D0/
      data deprv/2*1000.0D0 , 2*1500.0D0 , 2000.0D0 , 1500.0D0 ,        &
           &     11*1000.D0 , 2000.D0 , 2*2000.D0/
      data deptv/20*3000.D0/
!     ******      iexsol is soil texture type (see subr soilbc)
!     ORIGINAL
!     data iexsol/6,6,6,6,7,8,6,3,6,6,5,12,6,6,6,6,5,6,6,6/
!     Laura  04/04/08 changed soil texture for desert: 3->1
      data iexsol/6 , 6 , 6 , 6 , 7 , 8 , 6 , 1 , 6 , 6 , 5 , 12 , 6 ,  &
              &   6 , 6 , 6 , 5 , 6 , 6 , 6/
!     ******      kolsol is soil color type (see subr. albedo)
!     Dec. 15, 2008
!     data kolsol/5,3,4,4,4,4,4,1,3,3,2,1,5,5,5,4,3,4,4,4/
      data kolsol/6 , 4 , 5 , 5 , 5 , 5 , 5 , 1 , 4 , 4 , 2 , 1 , 6 ,   &
             &    6 , 6 , 5 , 4 , 5 , 5 , 5/
!     Dec. 15, 2008_
!     ******      xmopor is fraction of soil that is voids
      data xmopor/0.33D0 , 0.36D0 , 0.39D0 , 0.42D0 , 0.45D0 , 0.48D0 , &
             &    0.51D0 , 0.54D0 , 0.57D0 , 0.6D0 , 0.63D0 , 0.66D0/
!     ******      xmosuc is the minimum soil suction (mm)
      data xmosuc/3*30.0D0 , 9*200.D0/
!     ******      xmohyd is the max. hydraulic conductivity (mm/s)
      data xmohyd/0.20D-0 , 0.80D-1 , 0.32D-1 , 0.13D-1 , 0.89D-2 ,     &
         & 0.63D-2 , 0.45D-2 , 0.32D-2 , 0.22D-2 , 0.16D-2 , 0.11D-2 ,  &
         & 0.80D-3/
!     ******      xmowilt is fraction of water content at which
!     permanent wilting occurs
      data xmowil/0.095D0 , 0.128D0 , 0.161D0 , 0.266D0 , 0.3D0 ,       &
         &        0.332D0 , 0.378D0 , 0.419D0 , 0.455D0 , 0.487D0 ,     &
         &        0.516D0 , 0.542D0 /
      data xmofc/0.404D0 , 0.477D0 , 0.547D0 , 0.614D0 , 0.653D0 ,      &
         &       0.688D0 , 0.728D0 , 0.763D0 , 0.794D0 , 0.820D0 ,      &
         &       0.845D0 , 0.866D0/
!     ******      bee is the clapp and hornbereger "b" parameter
      data bee/3.5D0 , 4.0D0 , 4.5D0 , 5.0D0 , 5.5D0 , 6.0D0 , 6.8D0 ,  &
         &     7.6D0 , 8.4D0 , 9.2D0 , 10.0D0 , 10.8D0/
!     ******      bskrat is ratio of soil thermal conduc. to that of
!     loam - a function of texture
      data skrat/1.7D0 , 1.5D0 , 1.3D0 , 1.2D0 , 1.1D0 , 1.0D0 , .95D0 ,&
         &       .90D0 , .85D0 , .80D0 , .75D0 , .7D0/

!     Dec. 15, 2008
!     ******      albvgs is vegetation albedo for wavelengths < 0.7
!     microns data
!     albvgs/.1,.1,.05,.05,.08,.04,.08,.2,.1,.08,.17,.8,.06,2*.07,
!     &.05,.08,.06,2*0.06/
      data albvgs/.1D0 , .1D0 , .04D0 , .04D0 , .06D0 , .04D0 , .08D0 , &
          &       .2D0 , .1D0 , .08D0 , .17D0 , .8D0 , .06D0 , 2*.07D0 ,&
          &       .05D0 , .08D0 , .05D0 , 2*0.06D0/
!     ******      albvgl is vegetation albedo for wavelengths > 0.7
!     microns data
!     albvgl/.3,.3,.23,.23,.28,.20,.30,.4,.3,.28,.34,.6,.18,2*.2,
!     &.23,.28,.24,2*.18/
      data albvgl/.3D0 , .3D0 , .20D0 , .20D0 , .26D0 , .20D0 , .30D0 , &
         &        .4D0 , .3D0 , .28D0 , .34D0 , .6D0 , .18D0 , 2*.2D0 , &
         &        .23D0 , .28D0 , .23D0 , 2*.18D0/
!     Dec. 15, 2008_

      contains

        subroutine fillbat
          implicit none
          fbat(:,:,1)  = u10m_o(:,:)
          fbat(:,:,2)  = v10m_o(:,:)
          fbat(:,:,3)  = drag_o(:,:)
          fbat(:,:,4)  = tg_o(:,:)
          fbat(:,:,5)  = tlef_o(:,:)
          fbat(:,:,6)  = t2m_o(:,:)
          fbat(:,:,7)  = q2m_o(:,:)
          fbat(:,:,8)  = ssw_o(:,:)
          fbat(:,:,9)  = rsw_o(:,:)
          fbat(:,:,10) = tpr_o(:,:)
          fbat(:,:,11) = evpa_o(:,:)
          fbat(:,:,12) = rnos_o(:,:)
          fbat(:,:,13) = scv_o(:,:)
          fbat(:,:,14) = sena_o(:,:)
          fbat(:,:,15) = flwa_o(:,:)
          fbat(:,:,16) = fswa_o(:,:)
          fbat(:,:,17) = flwd_o(:,:)
          fbat(:,:,18) = sina_o(:,:)
          fbat(:,:,19) = prcv_o(:,:)
          fbat(:,:,20) = ps_o(:,:)
          fbat(:,:,21) = zpbl_o(:,:)
          fbat(:,:,22) = tgmx_o(:,:)
          fbat(:,:,23) = tgmn_o(:,:)
          fbat(:,:,24) = t2mx_o(:,:)
          fbat(:,:,25) = t2mn_o(:,:)
          fbat(:,:,26) = w10x_o(:,:)
          fbat(:,:,27) = psmn_o(:,:)
        end subroutine fillbat

        subroutine fillsub
        implicit none
          fsub(:,:,:,1)  = u10m_s(:,:,:)
          fsub(:,:,:,2)  = v10m_s(:,:,:)
          fsub(:,:,:,3)  = drag_s(:,:,:)
          fsub(:,:,:,4)  = tg_s(:,:,:)
          fsub(:,:,:,5)  = tlef_s(:,:,:)
          fsub(:,:,:,6)  = t2m_s(:,:,:)
          fsub(:,:,:,7)  = q2m_s(:,:,:)
          fsub(:,:,:,8)  = ssw_s(:,:,:)
          fsub(:,:,:,9)  = rsw_s(:,:,:)
          fsub(:,:,:,10) = tpr_s(:,:,:)
          fsub(:,:,:,11) = evpa_s(:,:,:)
          fsub(:,:,:,12) = rnos_s(:,:,:)
          fsub(:,:,:,13) = scv_s(:,:,:)
          fsub(:,:,:,14) = sena_s(:,:,:)
          fsub(:,:,:,15) = prcv_s(:,:,:)
          fsub(:,:,:,16) = ps_s(:,:,:)
        end subroutine fillsub

      end module mod_bats
