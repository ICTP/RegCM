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
  use mod_bats_internal , only : allocate_mod_bats_internal , dtlake , dtbat
  use mod_bats_internal , only : llake

  public :: dtlake , dtbat , llake

  logical :: ldcsst , lseaice , ldesseas

  integer(ik8) :: ntcpl  ! Number of time step to call ROMS update 
  integer(ik8) :: ntsrf2 ! Number of time step to call BATs 
  real(rk8) :: rdnnsg
  real(rk4) :: rrnnsg

  ! How many soil model steps for a day
  real(rk4) :: fdaysrf

  ! TO be pushed in savefile
  real(rk8) , pointer , dimension(:,:,:) :: ldew1
  real(rk8) , pointer , dimension(:,:,:) :: gwet1
  real(rk8) , pointer , dimension(:,:,:) :: snag1
  real(rk8) , pointer , dimension(:,:,:) :: sncv1
  real(rk8) , pointer , dimension(:,:,:) :: sfice1
  real(rk8) , pointer , dimension(:,:,:) :: rsw1
  real(rk8) , pointer , dimension(:,:,:) :: ssw1
  real(rk8) , pointer , dimension(:,:,:) :: tsw1
  real(rk8) , pointer , dimension(:,:,:) :: taf1
  real(rk8) , pointer , dimension(:,:,:) :: tgrd1
  real(rk8) , pointer , dimension(:,:,:) :: tgbrd1
  real(rk8) , pointer , dimension(:,:,:) :: tlef1
  real(rk8) , pointer , dimension(:,:,:) :: emiss1

  real(rk8) , pointer , dimension(:,:,:) :: sent1
  real(rk8) , pointer , dimension(:,:,:) :: drag1
  real(rk8) , pointer , dimension(:,:,:) :: evpr1
  real(rk8) , pointer , dimension(:,:,:) :: q2m
  real(rk8) , pointer , dimension(:,:,:) :: ps1
  real(rk8) , pointer , dimension(:,:,:) :: trnof1
  real(rk8) , pointer , dimension(:,:,:) :: srnof1
  real(rk8) , pointer , dimension(:,:,:) :: t2m
  real(rk8) , pointer , dimension(:,:,:) :: u10m
  real(rk8) , pointer , dimension(:,:,:) :: v10m
  real(rk8) , pointer , dimension(:,:,:) :: taux
  real(rk8) , pointer , dimension(:,:,:) :: tauy
  real(rk8) , pointer , dimension(:,:,:) :: prcp1
  logical , pointer , dimension(:,:,:) :: llakmsk1
  integer(ik4) , pointer , dimension(:,:,:) :: lakmsk1
!
  real(rk8) , pointer , dimension(:,:) :: ssw2da
  real(rk8) , pointer , dimension(:,:) :: sfracv2d
  real(rk8) , pointer , dimension(:,:) :: sfracb2d
  real(rk8) , pointer , dimension(:,:) :: sfracs2d
  real(rk8) , pointer , dimension(:,:) :: svegfrac2d
!
  ! Coupling variables
  real(rk8) , pointer , dimension(:,:,:) :: dailyrnf
  integer(ik4) , pointer , dimension(:,:) :: cplmsk
  real(rk8) :: runoffcount = 0.0D0
  ! dtskin is difference between skin temp and bulk sst
  real(rk8) , pointer , dimension(:,:) :: deltas
  real(rk8) , pointer , dimension(:,:) :: tdeltas
  real(rk8) , pointer , dimension(:,:) :: dtskin
  real(rk8) , pointer , dimension(:,:) :: sst
  ! Lake model
  real(rk8) , pointer , dimension(:,:,:,:) :: tlake
  real(rk8) , pointer , dimension(:,:,:) :: xlake
!
  data ldcsst /.false./
  data ldesseas /.false./

  real(rk8) , pointer , dimension(:,:) :: xlat          ! mddom%xlat
  real(rk8) , pointer , dimension(:,:) :: xlon          ! mddom%xlon
  real(rk8) , pointer , dimension(:,:) :: lndcat        ! mddom%lndcat
  real(rk8) , pointer , dimension(:,:) :: ht            ! mddom%ht
  real(rk8) , pointer , dimension(:,:) :: snowam        ! mddom%snowam
  integer(ik4) , pointer , dimension(:,:) :: iveg       ! mddom%iveg
  integer(ik4) , pointer , dimension(:,:) :: ldmsk      ! mddom%ldmsk

  real(rk8) , pointer , dimension(:,:,:) :: ht1         ! mdsub%ht
  real(rk8) , pointer , dimension(:,:,:) :: lndcat1     ! mdsub%lndcat
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

end module mod_bats_common
