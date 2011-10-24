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

module mod_bats_common
!
! Storage for Surface (BATS and shared by CLM) variables
!
  use mod_realkinds
  use mod_memutil
  use mod_dynparam
  use mod_bats_param
  use mod_bats_internal
!
  real(dp) :: xdtsec ! Atmosferic Model dt in seconds
  real(dp) :: dtbat  ! BATS1e internal timestep
  real(dp) :: dtlake ! Lake model internal timestep

  logical :: lemiss , lchem , ldcsst , llake , lseaice , ldesseas

  integer(8) :: kbats  ! Step frequency in calling BATS1e LSM

  integer :: iocnrough , iocnflx

  real(dp) , pointer , dimension(:,:) :: delq1d , delt1d , aldifs1d , &
         drag1d , evpr1d , gwet1d , ircp1d , ldew1d , aldirs1d ,      &
         p1d , q2m1d , rno1d , rnos1d , rsw1d , sag1d , scv1d , sent1d ,   &
         sice1d , ssw1d , t2m1d , tg1d , tgb1d , tlef1d ,    &
         tsw1d , u10m1d , v10m1d , veg1d
!
  real(dp) :: rdnnsg
  real(sp) :: rrnnsg
!
  real(dp) , pointer , dimension(:) :: flw1d , fsw1d
  real(dp) , pointer , dimension(:) :: czen , sola , vpdd
  real(dp) , pointer , dimension(:) :: ems , us1d , vs1d
  real(dp) , pointer , dimension(:) :: albdif , albdir , albvl ,      &
         albvs , albvsd , aldifl , aldifs , aldirl , aldirs , emiss ,&
         fracd , sabveg , solis , solvd , solvs , albvld
!
  real(dp) , pointer , dimension(:) :: coszrs
!
  real(dp) , pointer , dimension(:,:) :: flw2d , flwa2d , flwd2d ,     &
         flwda2d , fsw2d , fswa2d , pptc , pptnc , prca2d , prnca2d , &
         sabv2d , sina2d , sinc2d , sol2d , solvd2d , solvs2d , svga2d
!
  real(dp) , pointer , dimension(:,:) :: ssw2da , sdeltk2d , &
        sdelqk2d , sfracv2d , sfracb2d , sfracs2d , svegfrac2d
!
  integer , pointer , dimension(:,:,:) :: ocld2d , veg2d1
  integer , pointer , dimension(:,:) :: veg2d , ldmsk
!
  real(dp) , pointer , dimension(:,:,:) :: col2d , dew2d ,        &
       emiss2d , evpa2d , gwet2d , ircp2d , rno2d , rnos2d ,     &
       sag2d , scv2d , sena2d , sice2d , srw2d , ssw2d , swt2d , &
       taf2d , tg2d , tgb2d , tlef2d
!
  real(dp) , pointer , dimension(:,:,:) :: ht1 , lndcat1 , xlat1 , xlon1
!
  real(sp) , pointer , dimension(:,:,:) :: fbat
!
  real(sp) , pointer , dimension(:,:) :: drag_o , evpa_o , flwa_o ,  &
        flwd_o , fswa_o , prcv_o , psmn_o , ps_o , q2m_o , rnos_o , &
        rsw_o , scv_o , sena_o , sina_o , ssw_o , t2mn_o , t2mx_o , &
        t2m_o , tgmn_o , tgmx_o , tg_o , tlef_o , tpr_o , u10m_o ,  &
        v10m_o , w10x_o , zpbl_o , aldirs_o , aldifs_o
!
  real(sp) , pointer , dimension(:,:,:,:) :: fsub
!
  real(sp) , pointer , dimension(:,:,:) :: drag_s , evpa_s , prcv_s ,&
         ps_s , q2m_s , rnos_s , rsw_s , scv_s , sena_s , ssw_s ,   &
         t2m_s , tg_s , tlef_s , tpr_s , u10m_s , v10m_s
!
  ! dtskin is difference between skin temp and bulk sst
  real(dp) , pointer , dimension(:,:) :: deltas , tdeltas , dtskin
  logical , pointer , dimension(:,:) :: firstcall
!
  data lchem  /.false./
  data lemiss /.false./
  data ldcsst /.false./
  data llake  /.false./
  data ldesseas /.false./

  real(dp) , pointer , dimension(:,:) :: xlat          ! mddom%xlat
  real(dp) , pointer , dimension(:,:) :: xlon          ! mddom%xlon
  real(dp) , pointer , dimension(:,:) :: lndcat        ! mddom%lndcat
  real(dp) , pointer , dimension(:,:) :: ht            ! mddom%ht
  real(dp) , pointer , dimension(:,:) :: htf           ! mddom_io%ht
  real(dp) , pointer , dimension(:,:) :: tground1      ! sts1%tg
  real(dp) , pointer , dimension(:,:) :: tground2      ! sts2%tg
  real(dp) , pointer , dimension(:,:,:) :: uatm , vatm ! atms%ubx3d , atms%vbx3d
  real(dp) , pointer , dimension(:,:,:) :: tatm        ! atms%tb3d
  real(dp) , pointer , dimension(:,:,:) :: thatm       ! atms%thx3d
  real(dp) , pointer , dimension(:,:,:) :: qvatm       ! atms%qvb3d
  real(dp) , pointer , dimension(:,:) :: zpbl          ! sfsta%zpbl
  real(dp) , pointer , dimension(:,:) :: hfx           ! sfsta%hfx
  real(dp) , pointer , dimension(:,:) :: qfx           ! sfsta%qfx
  real(dp) , pointer , dimension(:,:) :: uvdrag        ! sfsta%uvdrag
  real(dp) , pointer , dimension(:,:) :: tgbb          ! sfsta%tgbb
  real(dp) , pointer , dimension(:,:) :: sfps          ! sps2%ps
  real(dp) , pointer , dimension(:,:,:) :: hgt         ! za
  real(dp) , pointer , dimension(:,:) :: ts            ! ts1
  real(dp) , pointer , dimension(:,:) :: tsf           ! ts0_io
  real(dp) , pointer , dimension(:,:) :: rho           ! rhox2d
  integer , pointer , dimension(:,:) :: lmask         ! CLM landmask

  contains

    subroutine allocate_mod_bats_common(ichem,idcsst)
    implicit none
    integer , intent(in) :: ichem , idcsst

    rrnnsg = 1.0/real(nnsg)
    rdnnsg = d_one/dble(nnsg)

    call getmem2d(veg2d,1,iym1,1,jxp,'bats:veg2d')
    call getmem2d(ldmsk,1,iym1,1,jxp,'bats:ldmsk')
    call getmem2d(flw2d,1,iym1,1,jxp,'bats:flw2d')
    call getmem2d(flwa2d,1,iym1,1,jxp,'bats:flwa2d')
    call getmem2d(flwd2d,1,iym1,1,jxp,'bats:flwd2d')
    call getmem2d(flwda2d,1,iym1,1,jxp,'bats:flwda2d')
    call getmem2d(fsw2d,1,iym1,1,jxp,'bats:fsw2d')
    call getmem2d(fswa2d,1,iym1,1,jxp,'bats:fswa2d')
    call getmem2d(pptc,1,iym1,1,jxp,'bats:pptc')
    call getmem2d(pptnc,1,iym1,1,jxp,'bats:pptnc')
    call getmem2d(prca2d,1,iym1,1,jxp,'bats:prca2d')
    call getmem2d(prnca2d,1,iym1,1,jxp,'bats:prnca2d')
    call getmem2d(sabv2d,1,iym1,1,jxp,'bats:sabv2d')
    call getmem2d(sina2d,1,iym1,1,jxp,'bats:sina2d')
    call getmem2d(sinc2d,1,iym1,1,jxp,'bats:sinc2d')
    call getmem2d(sol2d,1,iym1,1,jxp,'bats:sol2d')
    call getmem2d(solvd2d,1,iym1,1,jxp,'bats:solvd2d')
    call getmem2d(solvs2d,1,iym1,1,jxp,'bats:solvs2d')
    call getmem2d(svga2d,1,iym1,1,jxp,'bats:svga2d')
    if ( ichem == 1 ) then
      call getmem2d(ssw2da,1,iym1,1,jxp,'bats:ssw2da')
      call getmem2d(sdeltk2d,1,iym1,1,jxp,'bats:sdeltk2d')
      call getmem2d(sdelqk2d,1,iym1,1,jxp,'bats:sdelqk2d')
      call getmem2d(sfracv2d,1,iym1,1,jxp,'bats:sfracv2d')
      call getmem2d(sfracb2d,1,iym1,1,jxp,'bats:sfracb2d')
      call getmem2d(sfracs2d,1,iym1,1,jxp,'bats:sfracs2d')
      call getmem2d(svegfrac2d,1,iym1,1,jxp,'bats:svegfrac2d')
    end if
    call getmem3d(ocld2d,1,nnsg,1,iym1,1,jxp,'bats:ocld2d')
    call getmem3d(veg2d1,1,nnsg,1,iym1,1,jxp,'bats:veg2d1')
    call getmem3d(col2d,1,nnsg,1,iym1,1,jxp,'bats:col2d')
    call getmem3d(dew2d,1,nnsg,1,iym1,1,jxp,'bats:dew2d')
    call getmem3d(emiss2d,1,nnsg,1,iym1,1,jxp,'bats:emiss2d')
    call getmem3d(evpa2d,1,nnsg,1,iym1,1,jxp,'bats:evpa2d')
    call getmem3d(gwet2d,1,nnsg,1,iym1,1,jxp,'bats:gwet2d')
    call getmem3d(ircp2d,1,nnsg,1,iym1,1,jxp,'bats:ircp2d')
    call getmem3d(rno2d,1,nnsg,1,iym1,1,jxp,'bats:rno2d')
    call getmem3d(rnos2d,1,nnsg,1,iym1,1,jxp,'bats:rnos2d')
    call getmem3d(sag2d,1,nnsg,1,iym1,1,jxp,'bats:sag2d')
    call getmem3d(scv2d,1,nnsg,1,iym1,1,jxp,'bats:scv2d')
    call getmem3d(sena2d,1,nnsg,1,iym1,1,jxp,'bats:sena2d')
    call getmem3d(sice2d,1,nnsg,1,iym1,1,jxp,'bats:sice2d')
    call getmem3d(srw2d,1,nnsg,1,iym1,1,jxp,'bats:srw2d')
    call getmem3d(ssw2d,1,nnsg,1,iym1,1,jxp,'bats:ssw2d')
    call getmem3d(swt2d,1,nnsg,1,iym1,1,jxp,'bats:swt2d')
    call getmem3d(taf2d,1,nnsg,1,iym1,1,jxp,'bats:taf2d')
    call getmem3d(tg2d,1,nnsg,1,iym1,1,jxp,'bats:tg2d')
    call getmem3d(tgb2d,1,nnsg,1,iym1,1,jxp,'bats:tgb2d')
    call getmem3d(tlef2d,1,nnsg,1,iym1,1,jxp,'bats:tlef2d')

    call getmem3d(ht1,1,nnsg,1,iy,1,jxp,'bats:ht1')
    call getmem3d(lndcat1,1,nnsg,1,iy,1,jxp,'bats:lndcat1')
    call getmem3d(xlat1,1,nnsg,1,iy,1,jxp,'bats:xlat1')
    call getmem3d(xlon1,1,nnsg,1,iy,1,jxp,'bats:xlon1')

    if (idcsst == 1) then
      call getmem2d(deltas,1,iy,1,jxp,'bats:deltas')
      call getmem2d(tdeltas,1,iy,1,jxp,'bats:tdeltas')
      call getmem2d(dtskin,1,iy,1,jxp,'bats:dtskin')
      call getmem2d(firstcall,1,iy,1,jxp,'bats:firstcall')
    end if

    call getmem2d(delq1d,1,nnsg,1,iym1,'bats:delq1d')
    call getmem2d(delt1d,1,nnsg,1,iym1,'bats:delt1d')
    call getmem2d(drag1d,1,nnsg,1,iym1,'bats:drag1d')
    call getmem2d(evpr1d,1,nnsg,1,iym1,'bats:evpr1d')
    call getmem2d(gwet1d,1,nnsg,1,iym1,'bats:gwet1d')
    call getmem2d(ircp1d,1,nnsg,1,iym1,'bats:ircp1d')
    call getmem2d(ldew1d,1,nnsg,1,iym1,'bats:ldew1d')
    call getmem2d(p1d,1,nnsg,1,iym1,'bats:lp1d')
    call getmem2d(q2m1d,1,nnsg,1,iym1,'bats:q2m1d')
    call getmem2d(rno1d,1,nnsg,1,iym1,'bats:rno1d')
    call getmem2d(rnos1d,1,nnsg,1,iym1,'bats:rnos1d')
    call getmem2d(rsw1d,1,nnsg,1,iym1,'bats:rsw1d')
    call getmem2d(sag1d,1,nnsg,1,iym1,'bats:sag1d')
    call getmem2d(scv1d,1,nnsg,1,iym1,'bats:scv1d')
    call getmem2d(sent1d,1,nnsg,1,iym1,'bats:sent1d')
    call getmem2d(sice1d,1,nnsg,1,iym1,'bats:sice1d')
    call getmem2d(ssw1d,1,nnsg,1,iym1,'bats:ssw1d')
    call getmem2d(t2m1d,1,nnsg,1,iym1,'bats:t2m1d')
    call getmem2d(tg1d,1,nnsg,1,iym1,'bats:tg1d')
    call getmem2d(tgb1d,1,nnsg,1,iym1,'bats:tgb1d')
    call getmem2d(tlef1d,1,nnsg,1,iym1,'bats:tlef1d')
    call getmem2d(tsw1d,1,nnsg,1,iym1,'bats:tsw1d')
    call getmem2d(u10m1d,1,nnsg,1,iym1,'bats:u10m1d')
    call getmem2d(v10m1d,1,nnsg,1,iym1,'bats:v10m1d')
    call getmem2d(veg1d,1,nnsg,1,iym1,'bats:veg1d')
    call getmem2d(aldirs1d,1,nnsg,1,iym1,'bats:aldirs1d')
    call getmem2d(aldifs1d,1,nnsg,1,iym1,'bats:aldifs1d')
    call getmem1d(coszrs,1,iy,'bats:coszrs')
    call getmem1d(flw1d,1,iym1,'bats:flw1d')
    call getmem1d(fsw1d,1,iym1,'bats:fsw1d')
    call getmem1d(us1d,1,iym1,'bats:us1d')
    call getmem1d(vs1d,1,iym1,'bats:vs1d')
    call getmem1d(czen,1,iym1,'bats:czen')
    call getmem1d(sola,1,iym1,'bats:sola')
    call getmem1d(vpdd,1,iym1,'bats:vpdd')
    call getmem1d(ems,1,iym1,'bats:ems')
    call getmem1d(albdif,1,iym1,'bats:albdif')
    call getmem1d(albdir,1,iym1,'bats:albdir')
    call getmem1d(albvl,1,iym1,'bats:albvl')
    call getmem1d(albvld,1,iym1,'bats:albvld')
    call getmem1d(albvs,1,iym1,'bats:albvs')
    call getmem1d(albvsd,1,iym1,'bats:albvsd')
    call getmem1d(aldifl,1,iym1,'bats:aldifl')
    call getmem1d(aldifs,1,iym1,'bats:aldifs')
    call getmem1d(aldirl,1,iym1,'bats:aldirl')
    call getmem1d(aldirs,1,iym1,'bats:aldirs')
    call getmem1d(emiss,1,iym1,'bats:emiss')
    call getmem1d(fracd,1,iym1,'bats:fracd')
    call getmem1d(sabveg,1,iym1,'bats:sabveg')
    call getmem1d(solis,1,iym1,'bats:solis')
    call getmem1d(solvd,1,iym1,'bats:solvd')
    call getmem1d(solvs,1,iym1,'bats:solvs')

    call getmem3d(fbat,1,jxp,1,iym2,1,numbat,'bats:fbat')
    ps_o   => fbat(:,:,1)
    u10m_o => fbat(:,:,2)
    v10m_o => fbat(:,:,3)
    drag_o => fbat(:,:,4)
    tg_o   => fbat(:,:,5)
    tlef_o => fbat(:,:,6)
    t2m_o  => fbat(:,:,7)
    q2m_o  => fbat(:,:,8)
    ssw_o  => fbat(:,:,9)
    rsw_o  => fbat(:,:,10)
    tpr_o  => fbat(:,:,11)
    evpa_o => fbat(:,:,12)
    rnos_o => fbat(:,:,13)
    scv_o  => fbat(:,:,14)
    sena_o => fbat(:,:,15)
    flwa_o => fbat(:,:,16)
    fswa_o => fbat(:,:,17)
    flwd_o => fbat(:,:,18)
    sina_o => fbat(:,:,19)
    prcv_o => fbat(:,:,20)
    zpbl_o => fbat(:,:,21)
    tgmx_o => fbat(:,:,22)
    tgmn_o => fbat(:,:,23)
    t2mx_o => fbat(:,:,24)
    t2mn_o => fbat(:,:,25)
    w10x_o => fbat(:,:,26)
    psmn_o => fbat(:,:,27)
    aldirs_o => fbat(:,:,28)
    aldifs_o => fbat(:,:,29)

    call getmem4d(fsub,1,nnsg,1,jxp,1,iym2,1,numsub,'bats:fsub')
    ps_s   => fsub(:,:,:,1)
    u10m_s => fsub(:,:,:,2)
    v10m_s => fsub(:,:,:,3)
    drag_s => fsub(:,:,:,4)
    tg_s   => fsub(:,:,:,5)
    tlef_s => fsub(:,:,:,6)
    t2m_s  => fsub(:,:,:,7)
    q2m_s  => fsub(:,:,:,8)
    ssw_s  => fsub(:,:,:,9)
    rsw_s  => fsub(:,:,:,10)
    tpr_s  => fsub(:,:,:,11)
    evpa_s => fsub(:,:,:,12)
    rnos_s => fsub(:,:,:,13)
    scv_s  => fsub(:,:,:,14)
    sena_s => fsub(:,:,:,15)
    prcv_s => fsub(:,:,:,16)
!
    call allocate_mod_bats_internal
!
  end subroutine allocate_mod_bats_common
!
end module mod_bats_common
