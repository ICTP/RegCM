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
  use mod_bats_internal , only : allocate_mod_bats_internal
!
  real(rk8) :: dtbat  ! BATS1e internal timestep
  real(rk8) :: dtlake ! Lake model internal timestep

  logical :: ldcsst , llake , lseaice , ldesseas

  integer(ik8) :: ntcpl  ! Number of time step to call ROMS update 
  integer(ik8) :: ntsrf2 ! Number of time step to call BATs 
  real(rk8) :: rdnnsg
  real(rk4) :: rrnnsg

  ! How many soil model steps for a day
  real(rk4) :: fdaysrf

  real(rk8) , pointer , dimension(:,:,:) :: delt ,  taf , &
         drag , evpr , gwet , ldew , q2m , sfcp , trnof ,        &
         srnof , rsw , snag , sncv , sent , sfice , ssw ,        &
         t2m , tgrd , tgbrd , tlef , tsw , u10m , v10m , lncl ,  &
         taux , tauy , sfcemiss , cdrx , prcp
!
  real(rk8) , pointer , dimension(:,:) :: ssw2da , &
        sfracv2d , sfracb2d , sfracs2d , svegfrac2d
!
  integer(ik4) , pointer , dimension(:,:,:) :: lakemsk
  integer(ik4) , pointer , dimension(:,:) :: landmsk
!
  ! Coupling variables
  real(rk8) , pointer , dimension(:,:,:) :: dailyrnf
  integer(ik4) , pointer , dimension(:,:) :: cplmsk
  real(rk8) :: runoffcount = 0.0D0
  ! dtskin is difference between skin temp and bulk sst
  real(rk8) , pointer , dimension(:,:) :: deltas , tdeltas , dtskin
  real(rk8) , pointer , dimension(:,:) :: sst
  ! Lake model
  real(rk8) , pointer , dimension(:,:,:,:) :: tlake
  real(rk8) , pointer , dimension(:,:,:) :: xlake
!
  data ldcsst /.false./
  data llake  /.false./
  data ldesseas /.false./

  real(rk8) , pointer , dimension(:,:) :: xlat          ! mddom%xlat
  real(rk8) , pointer , dimension(:,:) :: xlon          ! mddom%xlon
  real(rk8) , pointer , dimension(:,:) :: lndcat        ! mddom%lndcat
  real(rk8) , pointer , dimension(:,:) :: ht            ! mddom%ht
  real(rk8) , pointer , dimension(:,:) :: snowam        ! mddom%snowam
  real(rk8) , pointer , dimension(:,:,:) :: ht1         ! mdsub%ht
  real(rk8) , pointer , dimension(:,:,:) :: lndcat1     ! mdsub%lndcat
  real(rk8) , pointer , dimension(:,:,:) :: mask1       ! mdsub%mask
  real(rk8) , pointer , dimension(:,:,:) :: xlat1       ! mdsub%xlat
  real(rk8) , pointer , dimension(:,:,:) :: xlon1       ! mdsub%xlon
  real(rk8) , pointer , dimension(:,:,:) :: dhlake1     ! mdsub%dhlake
  integer(ik4) , pointer , dimension(:,:,:) :: ldmsk1   ! mdsub%ldmsk
  integer(ik4) , pointer , dimension(:,:,:) :: iveg1    ! mdsub%iveg
  real(rk8) , pointer , dimension(:,:) :: uatm          ! atms%ubx3d(:,:,kz)
  real(rk8) , pointer , dimension(:,:) :: vatm          ! atms%vbx3d(:,:,kz)
  real(rk8) , pointer , dimension(:,:) :: tatm          ! atms%tb3d(:,:,kz)
  real(rk8) , pointer , dimension(:,:) :: thatm         ! atms%thx3d(:,:,kz)
  real(rk8) , pointer , dimension(:,:) :: qvatm         ! atms%qxb3d(:,:,kz,iqv)
  real(rk8) , pointer , dimension(:,:) :: hgt           ! za(:,:,kz)
  real(rk8) , pointer , dimension(:,:) :: hpbl          ! zpbl
  real(rk8) , pointer , dimension(:,:) :: hfx           ! sfs%hfx
  real(rk8) , pointer , dimension(:,:) :: qfx           ! sfs%qfx
  real(rk8) , pointer , dimension(:,:) :: tground1      ! sfs%tga
  real(rk8) , pointer , dimension(:,:) :: tground2      ! sfs%tgb
  real(rk8) , pointer , dimension(:,:) :: sfps          ! sfs%psb
  real(rk8) , pointer , dimension(:,:) :: uvdrag        ! sfs%uvdrag
  real(rk8) , pointer , dimension(:,:) :: tgbb          ! sfs%tgbb
  real(rk8) , pointer , dimension(:,:) :: rhox          ! rhox2d
  real(rk8) , pointer , dimension(:,:) :: rswf          ! fsw
  real(rk8) , pointer , dimension(:,:) :: rlwf          ! flw
  real(rk8) , pointer , dimension(:,:) :: dwrlwf        ! flwd
  real(rk8) , pointer , dimension(:,:) :: zencos        ! coszrs
  real(rk8) , pointer , dimension(:,:) :: ncprate       ! pptnc 
  real(rk8) , pointer , dimension(:,:) :: cprate        ! cprate 
  real(rk8) , pointer , dimension(:,:) :: vegswab       ! sabveg 
  real(rk8) , pointer , dimension(:,:) :: lwalb         ! albvl
  real(rk8) , pointer , dimension(:,:) :: swalb         ! albvs
  real(rk8) , pointer , dimension(:,:) :: swdiralb      ! aldirs
  real(rk8) , pointer , dimension(:,:) :: swdifalb      ! aldifs
  real(rk8) , pointer , dimension(:,:) :: lwdiralb      ! aldirl
  real(rk8) , pointer , dimension(:,:) :: lwdifalb      ! aldifl
  real(rk8) , pointer , dimension(:,:) :: swdir         ! solvs
  real(rk8) , pointer , dimension(:,:) :: swdif         ! solvsd
  real(rk8) , pointer , dimension(:,:) :: lwdir         ! solvl
  real(rk8) , pointer , dimension(:,:) :: lwdif         ! solvld
  real(rk8) , pointer , dimension(:,:) :: solinc        ! sinc
  real(rk8) , pointer , dimension(:,:) :: solar         ! solis
  real(rk8) , pointer , dimension(:,:) :: emissivity    ! emiss
  real(rk8) , pointer , dimension(:,:) :: deltaq        ! sdelq
  real(rk8) , pointer , dimension(:,:) :: deltat        ! sdelt
  integer(ik4) , pointer , dimension(:,:) :: lmask      ! CLM landmask

  integer(ik4) :: nlakep = 0
  integer(ik4) :: totlakep = 0

  contains

    subroutine allocate_mod_bats_common(ichem,idcsst,lakemod,iocncpl)
      implicit none
      integer(ik4) , intent(in) :: ichem , idcsst , lakemod , iocncpl

      rrnnsg = 1.0/real(nnsg)
      rdnnsg = d_one/dble(nnsg)

      if ( iocncpl == 1 ) then
        call getmem2d(cplmsk,jci1,jci2,ici1,ici2,'bats:cplmsk')
        cplmsk(:,:) = 0 
        ! This is for the RTM component
        call getmem3d(dailyrnf,jci1,jci2,ici1,ici2,1,2,'bats:dailyrnf')
      end if
      if ( ichem == 1 ) then
        call getmem2d(ssw2da,jci1,jci2,ici1,ici2,'bats:ssw2da')
        call getmem2d(sfracv2d,jci1,jci2,ici1,ici2,'bats:sfracv2d')
        call getmem2d(sfracb2d,jci1,jci2,ici1,ici2,'bats:sfracb2d')
        call getmem2d(sfracs2d,jci1,jci2,ici1,ici2,'bats:sfracs2d')
        call getmem2d(svegfrac2d,jci1,jci2,ici1,ici2,'bats:svegfrac2d')
      end if

      call getmem3d(delt,1,nnsg,jci1,jci2,ici1,ici2,'bats:delt')
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

      if (idcsst == 1) then
        call getmem2d(deltas,jci1,jci2,ici1,ici2,'bats:deltas')
        call getmem2d(tdeltas,jci1,jci2,ici1,ici2,'bats:tdeltas')
        call getmem2d(dtskin,jci1,jci2,ici1,ici2,'bats:dtskin')
        call getmem2d(sst,jci1,jci2,ici1,ici2,'bats:sst')
      end if

      call getmem3d(drag,1,nnsg,jci1,jci2,ici1,ici2,'bats:drag')
      call getmem3d(prcp,1,nnsg,jci1,jci2,ici1,ici2,'bats:prcp')
      call getmem3d(cdrx,1,nnsg,jci1,jci2,ici1,ici2,'bats:cdrx')
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
      call getmem3d(sfcemiss,1,nnsg,jci1,jci2,ici1,ici2,'bats:sfcemiss')
      call getmem3d(lncl,1,nnsg,jci1,jci2,ici1,ici2,'bats:lncl')

      if ( lakemod == 1 ) then
        call getmem3d(lakemsk,1,nnsg,jci1,jci2,ici1,ici2,'bats:lakemsk')
        call getmem3d(xlake,1,nnsg,jci1,jci2,ici1,ici2,'bats:xlake')
        call getmem4d(tlake,1,nnsg,jci1,jci2,ici1,ici2,1,ndpmax,'bats:tlake')
      end if

      call allocate_mod_bats_internal

    end subroutine allocate_mod_bats_common

end module mod_bats_common
