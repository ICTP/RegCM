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
  use mod_intkinds
  use mod_realkinds
  use mod_memutil
  use mod_dynparam
  use mod_bats_param
  use mod_bats_internal
!
  real(rk8) :: xdtsec ! Atmosferic Model dt in seconds
  real(rk8) :: dtbat  ! BATS1e internal timestep
  real(rk8) :: dtlake ! Lake model internal timestep

  logical :: lemiss , lchem , ldcsst , llake , lseaice , ldesseas

  integer(ik4) , parameter :: ps_o     = 1
  integer(ik4) , parameter :: u10m_o   = 2
  integer(ik4) , parameter :: v10m_o   = 3
  integer(ik4) , parameter :: drag_o   = 4
  integer(ik4) , parameter :: tg_o     = 5
  integer(ik4) , parameter :: tlef_o   = 6
  integer(ik4) , parameter :: t2m_o    = 7
  integer(ik4) , parameter :: q2m_o    = 8
  integer(ik4) , parameter :: ssw_o    = 9
  integer(ik4) , parameter :: rsw_o    = 10
  integer(ik4) , parameter :: tpr_o    = 11
  integer(ik4) , parameter :: evpa_o   = 12
  integer(ik4) , parameter :: rnos_o   = 13
  integer(ik4) , parameter :: scv_o    = 14
  integer(ik4) , parameter :: sena_o   = 15
  integer(ik4) , parameter :: flwa_o   = 16
  integer(ik4) , parameter :: fswa_o   = 17
  integer(ik4) , parameter :: flwd_o   = 18
  integer(ik4) , parameter :: sina_o   = 19
  integer(ik4) , parameter :: prcv_o   = 20
  integer(ik4) , parameter :: zpbl_o   = 21
  integer(ik4) , parameter :: aldirs_o = 22
  integer(ik4) , parameter :: aldifs_o = 23
  integer(ik4) , parameter :: sunt_o   = 24
  integer(ik4) , parameter :: tgmx_o   = 25
  integer(ik4) , parameter :: tgmn_o   = 26
  integer(ik4) , parameter :: t2mx_o   = 27
  integer(ik4) , parameter :: t2mn_o   = 28
  integer(ik4) , parameter :: tavg_o   = 29
  integer(ik4) , parameter :: w10x_o   = 30
  integer(ik4) , parameter :: pcpx_o   = 31
  integer(ik4) , parameter :: pcpa_o   = 32
  integer(ik4) , parameter :: sund_o   = 33
  integer(ik4) , parameter :: psmn_o   = 34

  integer(ik4) , parameter :: ps_s   = 1
  integer(ik4) , parameter :: u10m_s = 2
  integer(ik4) , parameter :: v10m_s = 3
  integer(ik4) , parameter :: drag_s = 4
  integer(ik4) , parameter :: tg_s   = 5
  integer(ik4) , parameter :: tlef_s = 6
  integer(ik4) , parameter :: t2m_s  = 7
  integer(ik4) , parameter :: q2m_s  = 8
  integer(ik4) , parameter :: ssw_s  = 9
  integer(ik4) , parameter :: rsw_s  = 10
  integer(ik4) , parameter :: tpr_s  = 11
  integer(ik4) , parameter :: evpa_s = 12
  integer(ik4) , parameter :: rnos_s = 13
  integer(ik4) , parameter :: scv_s  = 14
  integer(ik4) , parameter :: sena_s = 15
  integer(ik4) , parameter :: prcv_s = 16

  integer(8) :: kbats  ! Step frequency in calling BATS1e LSM
  integer(8) :: ntcpl  ! Number of time step to call ROMS update 
  integer(8) :: ntsrf2 ! Number of time step to call BATs 

  ! How many soil model steps for a day
  real(rk4) :: fdaysrf

  real(rk8) , pointer , dimension(:,:,:) :: delq , delt ,  taf , &
         drag , evpr , gwet , ldew , q2m , sfcp , trnof , &
         srnof , rsw , snag , sncv , sent , sfice , ssw , &
         t2m , tgrd , tgbrd , tlef , tsw , u10m , v10m , lncl, &
         taux, tauy
!
  real(rk8) :: rdnnsg
  real(rk4) :: rrnnsg
!
  real(rk8) , pointer , dimension(:,:) :: flw , fsw , fracd , &
         solis , czen , aemiss
  real(rk8) , pointer , dimension(:,:) :: albvl , albvs , albvsd , &
         aldifl , aldifs , aldirl , aldirs , sabveg
!
  real(rk8) , pointer , dimension(:,:) :: coszrs
!
  real(rk8) , pointer , dimension(:,:) :: flwd , pptc , pptnc , &
         prca , prnca , sinc , solvd , solvs , totpr
!
  real(rk8) , pointer , dimension(:,:) :: ssw2da , sdeltk2d , &
        sdelqk2d , sfracv2d , sfracb2d , sfracs2d , svegfrac2d
!
  integer(ik4) , pointer , dimension(:,:,:) :: ldmsk1 , iveg1
  integer(ik4) , pointer , dimension(:,:) :: iveg , ldmsk
!
  real(rk8) , pointer , dimension(:,:,:) :: runoff , emiss , evpa , sena , &
        srfrna
!
  real(rk8) , pointer , dimension(:,:,:) :: ht1 , lndcat1 , xlat1 , xlon1
!
  real(rk4) , pointer , dimension(:,:,:) :: fbat
!
  real(rk4) , pointer , dimension(:,:,:,:) :: fsub
!
  ! dtskin is difference between skin temp and bulk sst
  real(rk8) , pointer , dimension(:,:) :: deltas , tdeltas , dtskin
  ! Lake model
  real(rk8) , pointer , dimension(:,:,:) :: dhlake1
  integer(ik4) , pointer , dimension(:,:,:) :: idep
  real(rk8) , pointer , dimension(:,:,:) :: eta
  real(rk8) , pointer , dimension(:,:,:) :: hi
  real(rk8) , pointer , dimension(:,:,:) :: aveice
  real(rk8) , pointer , dimension(:,:,:) :: hsnow
  real(rk8) , pointer , dimension(:,:,:,:) :: tlak
!
  data lchem  /.false./
  data lemiss /.false./
  data ldcsst /.false./
  data llake  /.false./
  data ldesseas /.false./

  real(rk8) , pointer , dimension(:,:) :: xlat          ! mddom%xlat
  real(rk8) , pointer , dimension(:,:) :: xlon          ! mddom%xlon
  real(rk8) , pointer , dimension(:,:) :: lndcat        ! mddom%lndcat
  real(rk8) , pointer , dimension(:,:) :: ht            ! mddom%ht
  real(rk8) , pointer , dimension(:,:,:) :: uatm        ! atms%ubx3d
  real(rk8) , pointer , dimension(:,:,:) :: vatm        ! atms%vbx3d
  real(rk8) , pointer , dimension(:,:,:) :: tatm        ! atms%tb3d
  real(rk8) , pointer , dimension(:,:,:) :: thatm       ! atms%thx3d
  real(rk8) , pointer , dimension(:,:,:,:) :: qxatm     ! atms%qxb3d
  real(rk8) , pointer , dimension(:,:) :: hpbl          ! zpbl
  real(rk8) , pointer , dimension(:,:) :: hfx           ! sfs%hfx
  real(rk8) , pointer , dimension(:,:) :: qfx           ! sfs%qfx
  real(rk8) , pointer , dimension(:,:) :: uvdrag        ! sfs%uvdrag
  real(rk8) , pointer , dimension(:,:) :: tgbb          ! sfs%tgbb
  real(rk8) , pointer , dimension(:,:) :: tground1      ! sfs%tga
  real(rk8) , pointer , dimension(:,:) :: tground2      ! sfs%tgb
  real(rk8) , pointer , dimension(:,:) :: sfps          ! sfs%psb
  real(rk8) , pointer , dimension(:,:,:) :: hgt         ! za
  real(rk8) , pointer , dimension(:,:) :: tsf           ! ts0_io
  real(rk8) , pointer , dimension(:,:) :: rhox          ! rhox2d
  integer(ik4) , pointer , dimension(:,:) :: lmask          ! CLM landmask

  contains

    subroutine allocate_mod_bats_common(ichem,idcsst,lakemod)
      implicit none
      integer(ik4) , intent(in) :: ichem , idcsst , lakemod

      rrnnsg = 1.0/real(nnsg)
      rdnnsg = d_one/dble(nnsg)

      call getmem2d(iveg,jci1,jci2,ici1,ici2,'bats:iveg')
      call getmem2d(pptc,jci1,jci2,ici1,ici2,'bats:pptc')
      call getmem2d(pptnc,jci1,jci2,ici1,ici2,'bats:pptnc')
      call getmem2d(totpr,jci1,jci2,ici1,ici2,'bats:totpr')
      call getmem2d(prca,jci1,jci2,ici1,ici2,'bats:prca')
      call getmem2d(prnca,jci1,jci2,ici1,ici2,'bats:prnca')

      call getmem2d(sinc,jci1,jci2,ici1,ici2,'bats:sinc')
      call getmem2d(ldmsk,jci1,jci2,ici1,ici2,'bats:ldmsk')
      call getmem2d(flwd,jci1,jci2,ici1,ici2,'bats:flwd')
      call getmem2d(solvd,jci1,jci2,ici1,ici2,'bats:solvd')
      call getmem2d(solvs,jci1,jci2,ici1,ici2,'bats:solvs')
      if ( ichem == 1 ) then
        call getmem2d(ssw2da,jci1,jci2,ici1,ici2,'bats:ssw2da')
        call getmem2d(sdeltk2d,jci1,jci2,ici1,ici2,'bats:sdeltk2d')
        call getmem2d(sdelqk2d,jci1,jci2,ici1,ici2,'bats:sdelqk2d')
        call getmem2d(sfracv2d,jci1,jci2,ici1,ici2,'bats:sfracv2d')
        call getmem2d(sfracb2d,jci1,jci2,ici1,ici2,'bats:sfracb2d')
        call getmem2d(sfracs2d,jci1,jci2,ici1,ici2,'bats:sfracs2d')
        call getmem2d(svegfrac2d,jci1,jci2,ici1,ici2,'bats:svegfrac2d')
      end if
      call getmem3d(iveg1,1,nnsg,jci1,jci2,ici1,ici2,'bats:iveg1')
      call getmem3d(runoff,1,nnsg,jci1,jci2,ici1,ici2,'bats:runoff')
      call getmem3d(evpa,1,nnsg,jci1,jci2,ici1,ici2,'bats:evpa')
      call getmem3d(srfrna,1,nnsg,jci1,jci2,ici1,ici2,'bats:srfrna')

      call getmem3d(ldmsk1,1,nnsg,jci1,jci2,ici1,ici2,'bats:ldmsk1')
      call getmem3d(emiss,1,nnsg,jci1,jci2,ici1,ici2,'bats:emiss')
      call getmem3d(gwet,1,nnsg,jci1,jci2,ici1,ici2,'bats:gwet')
      call getmem3d(sena,1,nnsg,jci1,jci2,ici1,ici2,'bats:sena')
      call getmem3d(rsw,1,nnsg,jci1,jci2,ici1,ici2,'bats:rsw')
      call getmem3d(snag,1,nnsg,jci1,jci2,ici1,ici2,'bats:snag')
      call getmem3d(sncv,1,nnsg,jci1,jci2,ici1,ici2,'bats:sncv')
      call getmem3d(sent,1,nnsg,jci1,jci2,ici1,ici2,'bats:sent')
      call getmem3d(sfice,1,nnsg,jci1,jci2,ici1,ici2,'bats:sfice')
      call getmem3d(ssw,1,nnsg,jci1,jci2,ici1,ici2,'bats:ssw')
      call getmem3d(tgrd,1,nnsg,jci1,jci2,ici1,ici2,'bats:tgrd')
      call getmem3d(tgbrd,1,nnsg,jci1,jci2,ici1,ici2,'bats:tgbrd')
      call getmem3d(tlef,1,nnsg,jci1,jci2,ici1,ici2,'bats:tlef')
      call getmem3d(tsw,1,nnsg,jci1,jci2,ici1,ici2,'bats:tsw')
      call getmem3d(taf,1,nnsg,jci1,jci2,ici1,ici2,'bats:taf')
      call getmem3d(ldew,1,nnsg,jci1,jci2,ici1,ici2,'bats:ldew')
      call getmem2d(aemiss,jci1,jci2,ici1,ici2,'bats:aemiss')

      call getmem3d(ht1,1,nnsg,jde1,jde2,ide1,ide2,'bats:ht1')
      call getmem3d(lndcat1,1,nnsg,jde1,jde2,ide1,ide2,'bats:lndcat1')
      call getmem3d(xlat1,1,nnsg,jde1,jde2,ide1,ide2,'bats:xlat1')
      call getmem3d(xlon1,1,nnsg,jde1,jde2,ide1,ide2,'bats:xlon1')

      if (idcsst == 1) then
        call getmem2d(deltas,jci1,jci2,ici1,ici2,'bats:deltas')
        call getmem2d(tdeltas,jci1,jci2,ici1,ici2,'bats:tdeltas')
        call getmem2d(dtskin,jci1,jci2,ici1,ici2,'bats:dtskin')
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
      call getmem3d(taux,1,nnsg,jci1,jci2,ici1,ici2,'bats:taux')
      call getmem3d(tauy,1,nnsg,jci1,jci2,ici1,ici2,'bats:tauy')
      call getmem3d(lncl,1,nnsg,jci1,jci2,ici1,ici2,'bats:lncl')

      call getmem2d(flw,jci1,jci2,ici1,ici2,'bats:flw')
      call getmem2d(fsw,jci1,jci2,ici1,ici2,'bats:fsw')
      call getmem2d(fracd,jci1,jci2,ici1,ici2,'bats:fracd')
      call getmem2d(solis,jci1,jci2,ici1,ici2,'bats:solis')
      call getmem2d(sabveg,jci1,jci2,ici1,ici2,'bats:sabveg')

      call getmem2d(coszrs,jci1,jci2,ici1,ici2,'bats:coszrs')
      call getmem2d(albvl,jci1,jci2,ici1,ici2,'bats:albvl')
      call getmem2d(albvs,jci1,jci2,ici1,ici2,'bats:albvs')
      call getmem2d(aldifl,jci1,jci2,ici1,ici2,'bats:aldifl')
      call getmem2d(aldifs,jci1,jci2,ici1,ici2,'bats:aldifs')
      call getmem2d(aldirl,jci1,jci2,ici1,ici2,'bats:aldirl')
      call getmem2d(aldirs,jci1,jci2,ici1,ici2,'bats:aldirs')

      call getmem3d(fbat,jci1,jci2,ici1,ici2,1,numbat,'bats:fbat')
      call getmem4d(fsub,1,nnsg,jci1,jci2,ici1,ici2,1,numsub,'bats:fsub')
!
      if ( lakemod == 1 ) then
        call getmem3d(dhlake1,1,nnsg,jci1,jci2,ici1,ici2,'bats:dhlake1')
        call getmem3d(idep,1,nnsg,jci1,jci2,ici1,ici2,'bats:idep')
        call getmem3d(eta,1,nnsg,jci1,jci2,ici1,ici2,'bats:eta')
        call getmem3d(hi,1,nnsg,jci1,jci2,ici1,ici2,'bats:hi')
        call getmem3d(aveice,1,nnsg,jci1,jci2,ici1,ici2,'bats:aveice')
        call getmem3d(hsnow,1,nnsg,jci1,jci2,ici1,ici2,'bats:hsnow')
        call getmem4d(tlak,1,nnsg,jci1,jci2,ici1,ici2,1,ndpmax,'bats:tlak')
      end if

      call allocate_mod_bats_internal
!
    end subroutine allocate_mod_bats_common
!
end module mod_bats_common
