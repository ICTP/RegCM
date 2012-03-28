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

  ! How many soil model steps for a day
  real(sp) :: fdaysrf

  integer :: iocnrough , iocnflx, iocncpl

  real(dp) , pointer , dimension(:,:,:) :: delq , delt ,  taf , &
         drag , evpr , gwet , ldew , q2m , sfcp , trnof , &
         srnof , rsw , snag , sncv , sent , sfice , ssw , &
         t2m , tgrd , tgbrd , tlef , tsw , u10m , v10m , lncl
!
  real(dp) :: rdnnsg
  real(sp) :: rrnnsg
!
  real(dp) , pointer , dimension(:,:) :: flw , fsw , fracd , &
         solis , czen , aemiss
  real(dp) , pointer , dimension(:,:) :: albvl , albvs , albvsd , &
         aldifl , aldifs , aldirl , aldirs , sabveg
!
  real(dp) , pointer , dimension(:,:) :: coszrs
!
  real(dp) , pointer , dimension(:,:) :: flwd , pptc , pptnc , &
         prca , prnca , sinc , solvd , solvs , totpr
!
  real(dp) , pointer , dimension(:,:) :: ssw2da , sdeltk2d , &
        sdelqk2d , sfracv2d , sfracb2d , sfracs2d , svegfrac2d
!
  integer , pointer , dimension(:,:,:) :: ldmsk1 , iveg1
  integer , pointer , dimension(:,:) :: iveg , ldmsk
!
  real(dp) , pointer , dimension(:,:,:) :: runoff , emiss , evpa , sena , &
        srfrna
!
  real(dp) , pointer , dimension(:,:,:) :: ht1 , lndcat1 , xlat1 , xlon1
!
  real(sp) , pointer , dimension(:,:,:) :: fbat
!
  real(sp) , pointer , dimension(:,:) :: drag_o , evpa_o , flwa_o , &
        flwd_o , fswa_o , prcv_o , psmn_o , ps_o , q2m_o , rnos_o , &
        rsw_o , scv_o , sena_o , sina_o , ssw_o , t2mn_o , t2mx_o , &
        t2m_o , tgmn_o , tgmx_o , tg_o , tlef_o , tpr_o , u10m_o ,  &
        v10m_o , w10x_o , zpbl_o , aldirs_o , aldifs_o , pcpx_o ,   &
        pcpa_o , tavg_o , sunt_o , sund_o
!
  real(sp) , pointer , dimension(:,:,:,:) :: fsub
!
  real(sp) , pointer , dimension(:,:,:) :: drag_s , evpa_s , prcv_s ,&
         ps_s , q2m_s , rnos_s , rsw_s , scv_s , sena_s , ssw_s ,   &
         t2m_s , tg_s , tlef_s , tpr_s , u10m_s , v10m_s
!
  ! dtskin is difference between skin temp and bulk sst
  real(dp) , pointer , dimension(:,:) :: deltas , tdeltas , dtskin
  ! Lake model
  real(dp) , pointer , dimension(:,:,:) :: dhlake1
  integer , pointer , dimension(:,:,:) :: idep
  real(dp) , pointer , dimension(:,:,:) :: eta
  real(dp) , pointer , dimension(:,:,:) :: hi
  real(dp) , pointer , dimension(:,:,:) :: aveice
  real(dp) , pointer , dimension(:,:,:) :: hsnow
  real(dp) , pointer , dimension(:,:,:,:) :: tlak
!
  data lchem  /.false./
  data lemiss /.false./
  data ldcsst /.false./
  data llake  /.false./
  data ldesseas /.false./

  real(dp) , pointer , dimension(:,:) :: xlat          ! mddom%xlat
  real(dp) , pointer , dimension(:,:) :: lndcat        ! mddom%lndcat
  real(dp) , pointer , dimension(:,:) :: ht            ! mddom%ht
  real(dp) , pointer , dimension(:,:,:) :: uatm        ! atms%ubx3d
  real(dp) , pointer , dimension(:,:,:) :: vatm        ! atms%vbx3d
  real(dp) , pointer , dimension(:,:,:) :: tatm        ! atms%tb3d
  real(dp) , pointer , dimension(:,:,:) :: thatm       ! atms%thx3d
  real(dp) , pointer , dimension(:,:,:) :: qvatm       ! atms%qvb3d
  real(dp) , pointer , dimension(:,:) :: hpbl          ! zpbl
  real(dp) , pointer , dimension(:,:) :: hfx           ! sfs%hfx
  real(dp) , pointer , dimension(:,:) :: qfx           ! sfs%qfx
  real(dp) , pointer , dimension(:,:) :: uvdrag        ! sfs%uvdrag
  real(dp) , pointer , dimension(:,:) :: tgbb          ! sfs%tgbb
  real(dp) , pointer , dimension(:,:) :: tground1      ! sfs%tga
  real(dp) , pointer , dimension(:,:) :: tground2      ! sfs%tgb
  real(dp) , pointer , dimension(:,:) :: sfps          ! sfs%psb
  real(dp) , pointer , dimension(:,:,:) :: hgt         ! za
  real(dp) , pointer , dimension(:,:) :: ts            ! ts1
  real(dp) , pointer , dimension(:,:) :: tsf           ! ts0_io
  real(dp) , pointer , dimension(:,:) :: rho           ! rhox2d
  integer , pointer , dimension(:,:) :: lmask          ! CLM landmask

  contains

    subroutine allocate_mod_bats_common(ichem,idcsst,lakemod)
    implicit none
    integer , intent(in) :: ichem , idcsst , lakemod

    rrnnsg = 1.0/real(nnsg)
    rdnnsg = d_one/dble(nnsg)

    call getmem2d(iveg,jci1,jci2,ici1,ici2,'bats:iveg')
    call getmem2d(pptc,jci1,jci2,ici1,ici2,'bats:pptc')
    call getmem2d(pptnc,jci1,jci2,ici1,ici2,'bats:pptnc')
    call getmem2d(totpr,jci1,jci2,ici1,ici2,'bats:totpr')
    call getmem2d(prca,jci1,jci2,ici1,ici2,'bats:prca')
    call getmem2d(prnca,jci1,jci2,ici1,ici2,'bats:prnca')

    call getmem2d(sinc,jce1,jce2,ice1,ice2,'bats:sinc')
    call getmem2d(ldmsk,jce1,jce2,ice1,ice2,'bats:ldmsk')
    call getmem2d(flwd,jce1,jce2,ice1,ice2,'bats:flwd')
    call getmem2d(solvd,jce1,jce2,ice1,ice2,'bats:solvd')
    call getmem2d(solvs,jce1,jce2,ice1,ice2,'bats:solvs')
    if ( ichem == 1 ) then
      call getmem2d(ssw2da,jce1,jce2,ice1,ice2,'bats:ssw2da')
      call getmem2d(sdeltk2d,jce1,jce2,ice1,ice2,'bats:sdeltk2d')
      call getmem2d(sdelqk2d,jce1,jce2,ice1,ice2,'bats:sdelqk2d')
      call getmem2d(sfracv2d,jce1,jce2,ice1,ice2,'bats:sfracv2d')
      call getmem2d(sfracb2d,jce1,jce2,ice1,ice2,'bats:sfracb2d')
      call getmem2d(sfracs2d,jce1,jce2,ice1,ice2,'bats:sfracs2d')
      call getmem2d(svegfrac2d,jce1,jce2,ice1,ice2,'bats:svegfrac2d')
    end if
    call getmem3d(iveg1,1,nnsg,jci1,jci2,ici1,ici2,'bats:iveg1')
    call getmem3d(runoff,1,nnsg,jci1,jci2,ici1,ici2,'bats:runoff')
    call getmem3d(evpa,1,nnsg,jci1,jci2,ici1,ici2,'bats:evpa')
    call getmem3d(srfrna,1,nnsg,jci1,jci2,ici1,ici2,'bats:srfrna')

    call getmem3d(ldmsk1,1,nnsg,jce1,jce2,ice1,ice2,'bats:ldmsk1')
    call getmem3d(emiss,1,nnsg,jce1,jce2,ice1,ice2,'bats:emiss')
    call getmem3d(gwet,1,nnsg,jce1,jce2,ice1,ice2,'bats:gwet')
    call getmem3d(sena,1,nnsg,jce1,jce2,ice1,ice2,'bats:sena')
    call getmem3d(rsw,1,nnsg,jce1,jce2,ice1,ice2,'bats:rsw')
    call getmem3d(snag,1,nnsg,jce1,jce2,ice1,ice2,'bats:snag')
    call getmem3d(sncv,1,nnsg,jce1,jce2,ice1,ice2,'bats:sncv')
    call getmem3d(sent,1,nnsg,jce1,jce2,ice1,ice2,'bats:sent')
    call getmem3d(sfice,1,nnsg,jce1,jce2,ice1,ice2,'bats:sfice')
    call getmem3d(ssw,1,nnsg,jce1,jce2,ice1,ice2,'bats:ssw')
    call getmem3d(tgrd,1,nnsg,jce1,jce2,ice1,ice2,'bats:tgrd')
    call getmem3d(tgbrd,1,nnsg,jce1,jce2,ice1,ice2,'bats:tgbrd')
    call getmem3d(tlef,1,nnsg,jce1,jce2,ice1,ice2,'bats:tlef')
    call getmem3d(tsw,1,nnsg,jce1,jce2,ice1,ice2,'bats:tsw')
    call getmem3d(taf,1,nnsg,jce1,jce2,ice1,ice2,'bats:taf')
    call getmem3d(ldew,1,nnsg,jce1,jce2,ice1,ice2,'bats:ldew')
    call getmem2d(aemiss,jce1,jce2,ice1,ice2,'bats:aemiss')

    call getmem3d(ht1,1,nnsg,jde1,jde2,ide1,ide2,'bats:ht1')
    call getmem3d(lndcat1,1,nnsg,jde1,jde2,ide1,ide2,'bats:lndcat1')
    call getmem3d(xlat1,1,nnsg,jde1,jde2,ide1,ide2,'bats:xlat1')
    call getmem3d(xlon1,1,nnsg,jde1,jde2,ide1,ide2,'bats:xlon1')

    if (idcsst == 1) then
      call getmem2d(deltas,jce1,jce2,ice1,ice2,'bats:deltas')
      call getmem2d(tdeltas,jce1,jce2,ice1,ice2,'bats:tdeltas')
      call getmem2d(dtskin,jce1,jce2,ice1,ice2,'bats:dtskin')
    end if

    call getmem3d(delq,1,nnsg,jci1,jci2,ici1,ici2,'bats:delq')
    call getmem3d(delt,1,nnsg,jci1,jci2,ici1,ici2,'bats:delt')
    call getmem3d(drag,1,nnsg,jci1,jci2,ici1,ici2,'bats:drag')
    call getmem3d(evpr,1,nnsg,jci1,jci2,ici1,ici2,'bats:evpr')
    call getmem3d(sfcp,1,nnsg,jci1,jci2,ici1,ici2,'bats:sfcp')
    call getmem3d(q2m,1,nnsg,jci1,jci2,ici1,ici2,'bats:q2m')
    call getmem3d(trnof,1,nnsg,jci1,jci2,ici1,ici2,'bats:trnof')
    call getmem3d(srnof,1,nnsg,jci1,jci2,ici1,ici2,'bats:srnof')
    call getmem3d(t2m,1,nnsg,jci1,jci2,ici1,ici2,'bats:t2m')
    call getmem3d(u10m,1,nnsg,jci1,jci2,ici1,ici2,'bats:u10m')
    call getmem3d(v10m,1,nnsg,jci1,jci2,ici1,ici2,'bats:v10m')
    call getmem3d(lncl,1,nnsg,jci1,jci2,ici1,ici2,'bats:lncl')

    call getmem2d(flw,jce1,jce2,ice1,ice2,'bats:flw')
    call getmem2d(fsw,jce1,jce2,ice1,ice2,'bats:fsw')
    call getmem2d(fracd,jce1,jce2,ice1,ice2,'bats:fracd')
    call getmem2d(solis,jce1,jce2,ice1,ice2,'bats:solis')
    call getmem2d(sabveg,jce1,jce2,ice1,ice2,'bats:sabveg')

    call getmem2d(coszrs,jci1,jci2,ici1,ici2,'bats:coszrs')
    call getmem2d(albvl,jci1,jci2,ici1,ici2,'bats:albvl')
    call getmem2d(albvs,jci1,jci2,ici1,ici2,'bats:albvs')
    call getmem2d(aldifl,jci1,jci2,ici1,ici2,'bats:aldifl')
    call getmem2d(aldifs,jci1,jci2,ici1,ici2,'bats:aldifs')
    call getmem2d(aldirl,jci1,jci2,ici1,ici2,'bats:aldirl')
    call getmem2d(aldirs,jci1,jci2,ici1,ici2,'bats:aldirs')

    call getmem3d(fbat,jce1,jce2,ici1,ici2,1,numbat,'bats:fbat')
    ps_o(1:,2:) => fbat(:,:,1)
    u10m_o(1:,2:) => fbat(:,:,2)
    v10m_o(1:,2:) => fbat(:,:,3)
    drag_o(1:,2:) => fbat(:,:,4)
    tg_o(1:,2:)   => fbat(:,:,5)
    tlef_o(1:,2:) => fbat(:,:,6)
    t2m_o(1:,2:)  => fbat(:,:,7)
    q2m_o(1:,2:)  => fbat(:,:,8)
    ssw_o(1:,2:)  => fbat(:,:,9)
    rsw_o(1:,2:)  => fbat(:,:,10)
    tpr_o(1:,2:)  => fbat(:,:,11)
    evpa_o(1:,2:) => fbat(:,:,12)
    rnos_o(1:,2:) => fbat(:,:,13)
    scv_o(1:,2:)  => fbat(:,:,14)
    sena_o(1:,2:) => fbat(:,:,15)
    flwa_o(1:,2:) => fbat(:,:,16)
    fswa_o(1:,2:) => fbat(:,:,17)
    flwd_o(1:,2:) => fbat(:,:,18)
    sina_o(1:,2:) => fbat(:,:,19)
    prcv_o(1:,2:) => fbat(:,:,20)
    zpbl_o(1:,2:) => fbat(:,:,21)
    aldirs_o(1:,2:) => fbat(:,:,22)
    aldifs_o(1:,2:) => fbat(:,:,23)
    sunt_o(1:,2:) => fbat(:,:,24)
    tgmx_o(1:,2:) => fbat(:,:,25)
    tgmn_o(1:,2:) => fbat(:,:,26)
    t2mx_o(1:,2:) => fbat(:,:,27)
    t2mn_o(1:,2:) => fbat(:,:,28)
    tavg_o(1:,2:) => fbat(:,:,29)
    w10x_o(1:,2:) => fbat(:,:,30)
    pcpx_o(1:,2:) => fbat(:,:,31)
    pcpa_o(1:,2:) => fbat(:,:,32)
    sund_o(1:,2:) => fbat(:,:,33)
    psmn_o(1:,2:) => fbat(:,:,34)

    call getmem4d(fsub,1,nnsg,jce1,jce2,ici1,ici2,1,numsub,'bats:fsub')
    ps_s(1:,1:,2:)   => fsub(:,:,:,1)
    u10m_s(1:,1:,2:) => fsub(:,:,:,2)
    v10m_s(1:,1:,2:) => fsub(:,:,:,3)
    drag_s(1:,1:,2:) => fsub(:,:,:,4)
    tg_s(1:,1:,2:)   => fsub(:,:,:,5)
    tlef_s(1:,1:,2:) => fsub(:,:,:,6)
    t2m_s(1:,1:,2:)  => fsub(:,:,:,7)
    q2m_s(1:,1:,2:)  => fsub(:,:,:,8)
    ssw_s(1:,1:,2:)  => fsub(:,:,:,9)
    rsw_s(1:,1:,2:)  => fsub(:,:,:,10)
    tpr_s(1:,1:,2:)  => fsub(:,:,:,11)
    evpa_s(1:,1:,2:) => fsub(:,:,:,12)
    rnos_s(1:,1:,2:) => fsub(:,:,:,13)
    scv_s(1:,1:,2:)  => fsub(:,:,:,14)
    sena_s(1:,1:,2:) => fsub(:,:,:,15)
    prcv_s(1:,1:,2:) => fsub(:,:,:,16)
!
    if ( lakemod == 1 ) then
      call getmem3d(dhlake1,1,nnsg,jce1,jce2,ice1,ice2,'bats:dhlake1')
      call getmem3d(idep,1,nnsg,jce1,jce2,ice1,ice2,'bats:idep')
      call getmem3d(eta,1,nnsg,jce1,jce2,ice1,ice2,'bats:eta')
      call getmem3d(hi,1,nnsg,jce1,jce2,ice1,ice2,'bats:hi')
      call getmem3d(aveice,1,nnsg,jce1,jce2,ice1,ice2,'bats:aveice')
      call getmem3d(hsnow,1,nnsg,jce1,jce2,ice1,ice2,'bats:hsnow')
      call getmem4d(tlak,1,nnsg,jce1,jce2,ice1,ice2,1,ndpmax,'bats:tlak')
    end if

    call allocate_mod_bats_internal
!
  end subroutine allocate_mod_bats_common
!
end module mod_bats_common
