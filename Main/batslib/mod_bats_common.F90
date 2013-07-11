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
  real(rk8) :: dtbat  ! BATS1e internal timestep
  real(rk8) :: dtlake ! Lake model internal timestep

  logical :: ldcsst , llake , lseaice , ldesseas

  integer(ik8) :: ntcpl  ! Number of time step to call ROMS update 
  integer(ik8) :: ntsrf2 ! Number of time step to call BATs 

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
  real(rk8) , pointer , dimension(:,:,:) :: ht1 , lndcat1 , &
    mask1 , xlat1 , xlon1 , emiss
!
  ! Coupling variables
  real(rk8) , pointer , dimension(:,:,:) :: dailyrnf
  integer(ik4) , pointer , dimension(:,:) :: cplmsk
  real(rk8) :: runoffcount = 0.0D0
  ! dtskin is difference between skin temp and bulk sst
  real(rk8) , pointer , dimension(:,:) :: deltas , tdeltas , dtskin
  real(rk8) , pointer , dimension(:,:) :: sst
  ! Lake model
  real(rk8) , pointer , dimension(:,:,:) :: dhlake1
  integer(ik4) , pointer , dimension(:,:,:) :: idep
  real(rk8) , pointer , dimension(:,:,:) :: eta
  real(rk8) , pointer , dimension(:,:,:) :: hi
  real(rk8) , pointer , dimension(:,:,:) :: aveice
  real(rk8) , pointer , dimension(:,:,:) :: hsnow
  real(rk8) , pointer , dimension(:,:,:,:) :: tlak
!
  data ldcsst /.false./
  data llake  /.false./
  data ldesseas /.false./

  real(rk8) , pointer , dimension(:,:) :: xlat          ! mddom%xlat
  real(rk8) , pointer , dimension(:,:) :: xlon          ! mddom%xlon
  real(rk8) , pointer , dimension(:,:) :: lndcat        ! mddom%lndcat
  real(rk8) , pointer , dimension(:,:) :: ht            ! mddom%ht
  real(rk8) , pointer , dimension(:,:) :: snowam        ! mddom%snowam
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
  real(rk8) , pointer , dimension(:,:) :: rhox          ! rhox2d
  integer(ik4) , pointer , dimension(:,:) :: lmask      ! CLM landmask

  contains

    subroutine allocate_mod_bats_common(ichem,idcsst,lakemod,iocncpl)
      implicit none
      integer(ik4) , intent(in) :: ichem , idcsst , lakemod , iocncpl

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
      if ( iocncpl == 1 ) then
        call getmem2d(cplmsk,jci1,jci2,ici1,ici2,'bats:cplmsk')
        ! This is for the RTM component
        call getmem3d(dailyrnf,jci1,jci2,ici1,ici2,1,2,'bats:dailyrnf')
      end if
      cplmsk = 0 
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

      call getmem3d(ldmsk1,1,nnsg,jci1,jci2,ici1,ici2,'bats:ldmsk1')
      call getmem3d(emiss,1,nnsg,jci1,jci2,ici1,ici2,'bats:emiss')
      call getmem3d(gwet,1,nnsg,jci1,jci2,ici1,ici2,'bats:gwet')
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
      call getmem3d(mask1,1,nnsg,jde1,jde2,ide1,ide2,'bats:mask1')

      if (idcsst == 1) then
        call getmem2d(deltas,jci1,jci2,ici1,ici2,'bats:deltas')
        call getmem2d(tdeltas,jci1,jci2,ici1,ici2,'bats:tdeltas')
        call getmem2d(dtskin,jci1,jci2,ici1,ici2,'bats:dtskin')
        call getmem2d(sst,jci1,jci2,ici1,ici2,'bats:sst')
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
