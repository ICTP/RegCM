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
  integer(8) :: ntcpl  ! Number of time step to call ROMS update 
  integer(8) :: ntsrf2 ! Number of time step to call BATs 

  integer :: iocnrough , iocnflx

  real(dp) , pointer , dimension(:,:,:) :: delq , delt , albdifs , &
         drag , evpr , gwet , ircp , ldew , albdirs ,      &
         sfcp , q2m , trnof , srnof , rsw , snag , sncv , sent ,   &
         sfice , ssw , t2m , tgrd , tgbrd , tlef ,    &
         tsw , u10m , v10m , lncl
!
  real(dp) :: rdnnsg
  real(sp) :: rrnnsg
!
  real(dp) , pointer , dimension(:,:) :: flw , fsw , fracd , &
         solis , czen , aemiss
  real(dp) , pointer , dimension(:,:) :: albdif , albdir , albvl , &
         albvs , albvsd , aldifl , aldifs , aldirl , aldirs ,    &
         sabveg
!
  real(dp) , pointer , dimension(:,:) :: coszrs
!
  real(dp) , pointer , dimension(:,:) :: flw2d , flwd2d ,    &
         fsw2d , pptc , pptnc , prca2d , prnca2d , &
         sinc , solvd , solvs
!
  real(dp) , pointer , dimension(:,:) :: ssw2da , sdeltk2d , &
        sdelqk2d , sfracv2d , sfracb2d , sfracs2d , svegfrac2d
!
  integer , pointer , dimension(:,:,:) :: ocld2d , veg2d1
  integer , pointer , dimension(:,:) :: veg2d , ldmsk
!
  real(dp) , pointer , dimension(:,:,:) :: dew2d , runoff , &
       srfrno , emiss , evpa2d , sena2d
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
  real(dp) , pointer , dimension(:,:,:) :: uatm        ! atms%ubx3d
  real(dp) , pointer , dimension(:,:,:) :: vatm        ! atms%vbx3d
  real(dp) , pointer , dimension(:,:,:) :: tatm        ! atms%tb3d
  real(dp) , pointer , dimension(:,:,:) :: thatm       ! atms%thx3d
  real(dp) , pointer , dimension(:,:,:) :: qvatm       ! atms%qvb3d
  real(dp) , pointer , dimension(:,:) :: hpbl          ! zpbl
  real(dp) , pointer , dimension(:,:) :: hfx           ! sfsta%hfx
  real(dp) , pointer , dimension(:,:) :: qfx           ! sfsta%qfx
  real(dp) , pointer , dimension(:,:) :: uvdrag        ! sfsta%uvdrag
  real(dp) , pointer , dimension(:,:) :: tgbb          ! sfsta%tgbb
  real(dp) , pointer , dimension(:,:) :: sfps          ! sps2%ps
  real(dp) , pointer , dimension(:,:,:) :: hgt         ! za
  real(dp) , pointer , dimension(:,:) :: ts            ! ts1
  real(dp) , pointer , dimension(:,:) :: tsf           ! ts0_io
  real(dp) , pointer , dimension(:,:) :: rho           ! rhox2d
  integer , pointer , dimension(:,:) :: lmask          ! CLM landmask

  contains

    subroutine allocate_mod_bats_common(ichem,idcsst)
    implicit none
    integer , intent(in) :: ichem , idcsst

    rrnnsg = 1.0/real(nnsg)
    rdnnsg = d_one/dble(nnsg)

    call getmem2d(veg2d,1,jxp,1,iym1,'bats:veg2d')
    call getmem2d(ldmsk,1,jxp,1,iym1,'bats:ldmsk')
    call getmem2d(flw2d,1,jxp,1,iym1,'bats:flw2d')
    call getmem2d(flwd2d,1,jxp,1,iym1,'bats:flwd2d')
    call getmem2d(fsw2d,1,jxp,1,iym1,'bats:fsw2d')
    call getmem2d(pptc,1,jxp,1,iym1,'bats:pptc')
    call getmem2d(pptnc,1,jxp,1,iym1,'bats:pptnc')
    call getmem2d(prca2d,1,jxp,1,iym1,'bats:prca2d')
    call getmem2d(prnca2d,1,jxp,1,iym1,'bats:prnca2d')
    call getmem2d(sinc,1,jxp,1,iym1,'bats:sinc')
    call getmem2d(solvd,1,jxp,1,iym1,'bats:solvd')
    call getmem2d(solvs,1,jxp,1,iym1,'bats:solvs')
    if ( ichem == 1 ) then
      call getmem2d(ssw2da,1,jxp,1,iym1,'bats:ssw2da')
      call getmem2d(sdeltk2d,1,jxp,1,iym1,'bats:sdeltk2d')
      call getmem2d(sdelqk2d,1,jxp,1,iym1,'bats:sdelqk2d')
      call getmem2d(sfracv2d,1,jxp,1,iym1,'bats:sfracv2d')
      call getmem2d(sfracb2d,1,jxp,1,iym1,'bats:sfracb2d')
      call getmem2d(sfracs2d,1,jxp,1,iym1,'bats:sfracs2d')
      call getmem2d(svegfrac2d,1,jxp,1,iym1,'bats:svegfrac2d')
    end if
    call getmem3d(ocld2d,1,nnsg,1,jxp,1,iym1,'bats:ocld2d')
    call getmem3d(veg2d1,1,nnsg,1,jxp,1,iym1,'bats:veg2d1')
    call getmem3d(dew2d,1,nnsg,1,jxp,1,iym1,'bats:dew2d')
    call getmem3d(emiss,1,nnsg,1,jxp,1,iym1,'bats:emiss')
    call getmem3d(evpa2d,1,nnsg,1,jxp,1,iym1,'bats:evpa2d')
    call getmem3d(gwet,1,nnsg,1,jxp,1,iym1,'bats:gwet')
    call getmem3d(ircp,1,nnsg,1,jxp,1,iym1,'bats:ircp')
    call getmem3d(runoff,1,nnsg,1,jxp,1,iym1,'bats:runoff')
    call getmem3d(srfrno,1,nnsg,1,jxp,1,iym1,'bats:srfrno')
    call getmem3d(sena2d,1,nnsg,1,jxp,1,iym1,'bats:sena2d')

    call getmem3d(ht1,1,nnsg,1,jxp,1,iy,'bats:ht1')
    call getmem3d(lndcat1,1,nnsg,1,jxp,1,iy,'bats:lndcat1')
    call getmem3d(xlat1,1,nnsg,1,jxp,1,iy,'bats:xlat1')
    call getmem3d(xlon1,1,nnsg,1,jxp,1,iy,'bats:xlon1')

    if (idcsst == 1) then
      call getmem2d(deltas,1,iy,1,jxp,'bats:deltas')
      call getmem2d(tdeltas,1,iy,1,jxp,'bats:tdeltas')
      call getmem2d(dtskin,1,iy,1,jxp,'bats:dtskin')
      call getmem2d(firstcall,1,iy,1,jxp,'bats:firstcall')
    end if

    call getmem3d(delq,1,nnsg,1,jxp,1,iym1,'bats:delq')
    call getmem3d(delt,1,nnsg,1,jxp,1,iym1,'bats:delt')
    call getmem3d(drag,1,nnsg,1,jxp,1,iym1,'bats:drag')
    call getmem3d(evpr,1,nnsg,1,jxp,1,iym1,'bats:evpr')
    call getmem3d(ldew,1,nnsg,1,jxp,1,iym1,'bats:ldew')
    call getmem3d(sfcp,1,nnsg,1,jxp,1,iym1,'bats:sfcp')
    call getmem3d(q2m,1,nnsg,1,jxp,1,iym1,'bats:q2m')
    call getmem3d(trnof,1,nnsg,1,jxp,1,iym1,'bats:trnof')
    call getmem3d(srnof,1,nnsg,1,jxp,1,iym1,'bats:srnof')
    call getmem3d(rsw,1,nnsg,1,jxp,1,iym1,'bats:rsw')
    call getmem3d(snag,1,nnsg,1,jxp,1,iym1,'bats:snag')
    call getmem3d(sncv,1,nnsg,1,jxp,1,iym1,'bats:sncv')
    call getmem3d(sent,1,nnsg,1,jxp,1,iym1,'bats:sent')
    call getmem3d(sfice,1,nnsg,1,jxp,1,iym1,'bats:sfice')
    call getmem3d(ssw,1,nnsg,1,jxp,1,iym1,'bats:ssw')
    call getmem3d(t2m,1,nnsg,1,jxp,1,iym1,'bats:t2m')
    call getmem3d(tgrd,1,nnsg,1,jxp,1,iym1,'bats:tgrd')
    call getmem3d(tgbrd,1,nnsg,1,jxp,1,iym1,'bats:tgbrd')
    call getmem3d(tlef,1,nnsg,1,jxp,1,iym1,'bats:tlef')
    call getmem3d(tsw,1,nnsg,1,jxp,1,iym1,'bats:tsw')
    call getmem3d(u10m,1,nnsg,1,jxp,1,iym1,'bats:u10m')
    call getmem3d(v10m,1,nnsg,1,jxp,1,iym1,'bats:v10m')
    call getmem3d(lncl,1,nnsg,1,jxp,1,iym1,'bats:lncl')
    call getmem3d(albdirs,1,nnsg,1,jxp,1,iym1,'bats:albdirs')
    call getmem3d(albdifs,1,nnsg,1,jxp,1,iym1,'bats:albdifs')

    call getmem2d(aemiss,1,jxp,1,iym1,'bats_internal:aemiss')
    call getmem2d(flw,1,jxp,1,iym1,'bats_internal:flw')
    call getmem2d(fsw,1,jxp,1,iym1,'bats_internal:fsw')
    call getmem2d(fracd,1,jxp,1,iym1,'bats_internal:fracd')
    call getmem2d(solis,1,jxp,1,iym1,'bats_internal:solis')
    call getmem2d(coszrs,1,jxp,1,iy,'bats:coszrs')
    call getmem2d(sabveg,1,jxp,1,iym1,'bats:sabveg')
    call getmem2d(albdif,1,jxp,1,iym1,'bats:albdif')
    call getmem2d(albdir,1,jxp,1,iym1,'bats:albdir')
    call getmem2d(albvl,1,jxp,1,iym1,'bats:albvl')
    call getmem2d(albvs,1,jxp,1,iym1,'bats:albvs')
    call getmem2d(aldifl,1,jxp,1,iym1,'bats:aldifl')
    call getmem2d(aldifs,1,jxp,1,iym1,'bats:aldifs')
    call getmem2d(aldirl,1,jxp,1,iym1,'bats:aldirl')
    call getmem2d(aldirs,1,jxp,1,iym1,'bats:aldirs')

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
