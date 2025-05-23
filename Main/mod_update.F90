!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_update
  !
  !-----------------------------------------------------------------------
  !     Used module declarations
  !-----------------------------------------------------------------------
  !
  use mod_intkinds, only : ik4
  use mod_realkinds, only : rk8
  use mod_runparams, only : icopcpl
  use mod_regcm_types, only : exp_data, imp_data, exp_data3d
  use mod_memutil

  implicit none

  private

  public :: exp_data
  public :: imp_data
  type(imp_data), public :: importFields
  type(exp_data), public :: exportFields
  type(exp_data3d), public :: exportFields3d

  integer(ik4), pointer, contiguous, dimension(:,:) :: ldmskb => null( )
  integer(ik4), pointer, contiguous, dimension(:,:) :: wetdry => null( )

  real(rk8), parameter :: zeroval = 0.0_rk8
  real(rk8), parameter :: missing_r8 = 1.0e20_rk8
  real(rk8), parameter :: tol = missing_r8/2.0_rk8
  !
  !-----------------------------------------------------------------------
  !     Public subroutines
  !-----------------------------------------------------------------------
  !
  public :: RCM_Get
  public :: RCM_Put
  public :: RCM_Allocate
  !
  !-----------------------------------------------------------------------
  !     Module constants
  !-----------------------------------------------------------------------
  !
  real(rk8), parameter :: beta = 1.25_rk8 ! gustiness coeff
  real(rk8), parameter :: von  = 0.4_rk8  ! von Karman constant
  real(rk8), parameter :: fdg  = 1.00_rk8 ! ratio of thermal to wind von Karman
  real(rk8), parameter :: tdk  = 273.16_rk8
  real(rk8), parameter :: grav = 9.82_rk8 ! accel of earth grav

  contains

  subroutine RCM_Allocate()
    !
    !-----------------------------------------------------------------------
    !     Used module declarations
    !-----------------------------------------------------------------------
    !
    use mod_atm_interface, only : mddom
    use mod_dynparam, only : kz
    use mod_dynparam, only : ice1, ice2, jce1, jce2
    use mod_dynparam, only : ici1, ici2, jci1, jci2

    implicit none
    !
    !-----------------------------------------------------------------------
    !     Local variable declarations
    !-----------------------------------------------------------------------
    !
    integer(ik4) :: i, j, k
    real(rk8), parameter :: initval = 1.0e20_rk8
    real(rk8), parameter :: zeroval = 0.0e20_rk8
    !
    !-----------------------------------------------------------------------
    !     Allocate arrays
    !-----------------------------------------------------------------------
    !
    call getmem2d(exportFields%psfc,jce1,jce2,ice1,ice2,'cpl:psfc')
    call getmem2d(exportFields%tsfc,jce1,jce2,ice1,ice2,'cpl:tsfc')
    call getmem2d(exportFields%qsfc,jce1,jce2,ice1,ice2,'cpl:qsfc')
    call getmem2d(exportFields%swrd,jce1,jce2,ice1,ice2,'cpl:swrd')
    call getmem2d(exportFields%lwrd,jce1,jce2,ice1,ice2,'cpl:lwrd')
    call getmem2d(exportFields%dlwr,jce1,jce2,ice1,ice2,'cpl:dlwr')
    call getmem2d(exportFields%lhfx,jce1,jce2,ice1,ice2,'cpl:lhfx')
    call getmem2d(exportFields%shfx,jce1,jce2,ice1,ice2,'cpl:shfx')
    call getmem2d(exportFields%prec,jce1,jce2,ice1,ice2,'cpl:prec')
    call getmem2d(exportFields%wndu,jce1,jce2,ice1,ice2,'cpl:wndu')
    call getmem2d(exportFields%wndv,jce1,jce2,ice1,ice2,'cpl:wndv')
    call getmem2d(exportFields%rnof,jci1,jci2,ici1,ici2,'cpl:rnof')
    call getmem2d(exportFields%taux,jce1,jce2,ice1,ice2,'cpl:taux')
    call getmem2d(exportFields%tauy,jce1,jce2,ice1,ice2,'cpl:tauy')
    call getmem2d(exportFields%wspd,jce1,jce2,ice1,ice2,'cpl:wspd')
    call getmem2d(exportFields%wdir,jce1,jce2,ice1,ice2,'cpl:wdir')
    call getmem2d(exportFields%ustr,jce1,jce2,ice1,ice2,'cpl:ustr')
    call getmem2d(exportFields%nflx,jce1,jce2,ice1,ice2,'cpl:nflx')
    call getmem2d(exportFields%sflx,jce1,jce2,ice1,ice2,'cpl:sflx')
    call getmem2d(exportFields%snow,jce1,jce2,ice1,ice2,'cpl:snow')
    call getmem2d(exportFields%dswr,jce1,jce2,ice1,ice2,'cpl:dswr')
    call getmem2d(exportFields%rhoa,jce1,jce2,ice1,ice2,'cpl:rhoa')

    call getmem2d(importFields%sst,jce1,jce2,ice1,ice2,'cpl:sst')
    call getmem2d(importFields%sit,jce1,jce2,ice1,ice2,'cpl:sit')
    call getmem2d(importFields%msk,jce1,jce2,ice1,ice2,'cpl:msk')
    call getmem2d(importFields%zo,jce1,jce2,ice1,ice2,'cpl:zo')
    call getmem2d(importFields%ustar,jce1,jce2,ice1,ice2,'cpl:ustar')

    call getmem2d(ldmskb,jci1,jci2,ici1,ici2,'cpl:ldmsk')
    call getmem2d(wetdry,jci1,jci2,ici1,ici2,'cpl:wetdry')

    if ( icopcpl == 1 ) then
      call getmem3d(exportFields3d%u,jce1,jce2,ice1,ice2,1,kz,'cpl:u')
      call getmem3d(exportFields3d%v,jce1,jce2,ice1,ice2,1,kz,'cpl:v')
      call getmem3d(exportFields3d%w,jce1,jce2,ice1,ice2,1,kz,'cpl:w')
      call getmem3d(exportFields3d%t,jce1,jce2,ice1,ice2,1,kz,'cpl:t')
      call getmem3d(exportFields3d%q,jce1,jce2,ice1,ice2,1,kz,'cpl:q')
      call getmem3d(exportFields3d%cldfrc,jce1,jce2,ice1,ice2,1,kz,'cpl:cldfrc')
      call getmem3d(exportFields3d%cldlwc,jce1,jce2,ice1,ice2,1,kz,'cpl:cldlwc')
    end if
    !
    ! Initialize arrays
    !
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      exportFields%psfc(j,i) = initval
      exportFields%tsfc(j,i) = initval
      exportFields%qsfc(j,i) = initval
      exportFields%swrd(j,i) = initval
      exportFields%lwrd(j,i) = initval
      exportFields%dlwr(j,i) = initval
      exportFields%lhfx(j,i) = initval
      exportFields%shfx(j,i) = initval
      exportFields%prec(j,i) = initval
      exportFields%wndu(j,i) = initval
      exportFields%wndv(j,i) = initval
      exportFields%rnof(j,i) = initval
      exportFields%taux(j,i) = initval
      exportFields%tauy(j,i) = initval
      exportFields%wspd(j,i) = initval
      exportFields%wdir(j,i) = initval
      exportFields%ustr(j,i) = initval
      exportFields%nflx(j,i) = initval
      exportFields%sflx(j,i) = initval
      exportFields%snow(j,i) = initval
      exportFields%dswr(j,i) = initval
      exportFields%rhoa(j,i) = initval

      importFields%sst(j,i) = initval
      importFields%sit(j,i) = initval
      importFields%msk(j,i) = initval
      importFields%zo(j,i) = initval
      importFields%ustar(j,i) = initval

      ldmskb(j,i) = mddom%ldmsk(j,i)
      wetdry(j,i) = 0
    end do

    if ( icopcpl == 1 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        exportFields3d%u(j,i,k) = initval
        exportFields3d%v(j,i,k) = initval
        exportFields3d%w(j,i,k) = initval
        exportFields3d%t(j,i,k) = initval
        exportFields3d%q(j,i,k) = initval
        exportFields3d%cldfrc(j,i,k) = initval
        exportFields3d%cldlwc(j,i,k) = initval
      end do
    end if

  end subroutine RCM_Allocate

  subroutine RCM_Get(localPet)
    !-----------------------------------------------------------------------
    !     Used module declarations
    !-----------------------------------------------------------------------
    !
    use mod_constants
    use mod_lm_interface, only : import_data_into_surface

    implicit none
    !
    !-----------------------------------------------------------------------
    !     Imported variable declarations
    !-----------------------------------------------------------------------
    !
    integer, intent(in) :: localPet
    !
    !-----------------------------------------------------------------------
    !     Get information from OCN component
    !-----------------------------------------------------------------------
    !
    call import_data_into_surface(importFields,ldmskb,wetdry,tol)

   end subroutine RCM_Get

   subroutine RCM_Put(localPet)
    !
    !-----------------------------------------------------------------------
    !     Used module declarations
    !-----------------------------------------------------------------------
    !
    use mod_runparams, only : dtsec, alarm_day
    use mod_lm_interface, only : export_data_from_surface
    use mod_atm_interface, only : export_data_from_atm, sfs
    use mod_rad_interface, only : export_data_from_rad

    implicit none
    !
    !-----------------------------------------------------------------------
    !     Imported variable declarations
    !-----------------------------------------------------------------------
    !
    integer, intent(in) :: localPet
    !
    !-----------------------------------------------------------------------
    !     Send information to OCN, WAV and RTM component
    !-----------------------------------------------------------------------
    !
    call export_data_from_surface(exportFields)

    if ( icopcpl == 1 ) then
      !
      !-----------------------------------------------------------------------
      !     Send information to COP component
      !-----------------------------------------------------------------------
      !
      call export_data_from_atm(exportFields3d)
      call export_data_from_rad(exportFields3d)
    end if

    if ( alarm_day%act() ) then
      sfs%dtrnof = zeroval
    end if
  end subroutine RCM_Put

end module mod_update

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
