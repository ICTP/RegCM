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

      call getmem3d(gwet1,1,nnsg,jci1,jci2,ici1,ici2,'bats:gwet1')
      call getmem3d(rsw1,1,nnsg,jci1,jci2,ici1,ici2,'bats:rsw1')
      call getmem3d(snag1,1,nnsg,jci1,jci2,ici1,ici2,'bats:snag1')
      call getmem3d(sncv1,1,nnsg,jci1,jci2,ici1,ici2,'bats:sncv1')
      call getmem3d(sfice1,1,nnsg,jci1,jci2,ici1,ici2,'bats:sfice1')
      call getmem3d(ssw1,1,nnsg,jci1,jci2,ici1,ici2,'bats:ssw1')
      call getmem3d(tgrd1,1,nnsg,jci1,jci2,ici1,ici2,'bats:tgrd1')
      call getmem3d(tgbrd1,1,nnsg,jci1,jci2,ici1,ici2,'bats:tgbrd1')
      call getmem3d(tlef1,1,nnsg,jci1,jci2,ici1,ici2,'bats:tlef1')
      call getmem3d(tsw1,1,nnsg,jci1,jci2,ici1,ici2,'bats:tsw1')
      call getmem3d(taf1,1,nnsg,jci1,jci2,ici1,ici2,'bats:taf1')
      call getmem3d(ldew1,1,nnsg,jci1,jci2,ici1,ici2,'bats:ldew1')

      if (idcsst == 1) then
        call getmem2d(deltas,jci1,jci2,ici1,ici2,'bats:deltas')
        call getmem2d(tdeltas,jci1,jci2,ici1,ici2,'bats:tdeltas')
        call getmem2d(dtskin,jci1,jci2,ici1,ici2,'bats:dtskin')
        call getmem2d(sst,jci1,jci2,ici1,ici2,'bats:sst')
      end if

      call getmem3d(sent1,1,nnsg,jci1,jci2,ici1,ici2,'bats:sent1')
      call getmem3d(evpr1,1,nnsg,jci1,jci2,ici1,ici2,'bats:evpr1')
      call getmem3d(drag1,1,nnsg,jci1,jci2,ici1,ici2,'bats:drag1')
      call getmem3d(prcp1,1,nnsg,jci1,jci2,ici1,ici2,'bats:prcp1')
      call getmem3d(q2m,1,nnsg,jci1,jci2,ici1,ici2,'bats:q2m')
      call getmem3d(ps1,1,nnsg,jci1,jci2,ici1,ici2,'bats:ps1')
      call getmem3d(trnof1,1,nnsg,jci1,jci2,ici1,ici2,'bats:trnof1')
      call getmem3d(srnof1,1,nnsg,jci1,jci2,ici1,ici2,'bats:srnof1')
      call getmem3d(t2m,1,nnsg,jci1,jci2,ici1,ici2,'bats:t2m')
      call getmem3d(u10m,1,nnsg,jci1,jci2,ici1,ici2,'bats:u10m')
      call getmem3d(v10m,1,nnsg,jci1,jci2,ici1,ici2,'bats:v10m')
      call getmem3d(taux,1,nnsg,jci1,jci2,ici1,ici2,'bats:taux')
      call getmem3d(tauy,1,nnsg,jci1,jci2,ici1,ici2,'bats:tauy')
      call getmem3d(emiss1,1,nnsg,jci1,jci2,ici1,ici2,'bats:emiss1')

      if ( lakemod == 1 ) then
        call getmem3d(lakmsk1,1,nnsg,jci1,jci2,ici1,ici2,'bats:lakmsk1')
        call getmem3d(llakmsk1,1,nnsg,jci1,jci2,ici1,ici2,'bats:llakmsk1')
        call getmem3d(xlake,1,nnsg,jci1,jci2,ici1,ici2,'bats:xlake')
        call getmem4d(tlake,1,nnsg,jci1,jci2,ici1,ici2,1,ndpmax,'bats:tlake')
      end if

      call allocate_mod_bats_internal

    end subroutine allocate_mod_bats_common

end module mod_bats_common
