!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     This file is part of ICTP RegCM.
!
!     ICTP RegCM is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     ICTP RegCM is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
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
  use mod_regcm_types , only : exp_data , imp_data
  use mod_memutil
!
  implicit none
  private

  public :: exp_data
  public :: imp_data
  type(imp_data) , public :: importFields
  type(exp_data) , public :: exportFields
!
  integer(ik4) , pointer :: ldmskb(:,:)
  integer(ik4) , pointer :: wetdry(:,:)
!
  real(rk8) , parameter :: zeroval = 0.0d0
  real(rk8) , parameter :: missing_r8 = 1.0d20
  real(rk8) , parameter :: tol = missing_r8/2.0d0
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
  real(rk8), parameter :: beta = 1.25 ! gustiness coeff
  real(rk8), parameter :: von  = 0.4  ! von Karman constant
  real(rk8), parameter :: fdg  = 1.00 ! ratio of thermal to wind von Karman
  real(rk8), parameter :: tdk  = 273.16
  real(rk8), parameter :: grav = 9.82 ! accel of earth grav
!
  contains
!
  subroutine RCM_Allocate()
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
    use mod_atm_interface , only : mddom
    use mod_dynparam, only : ice1, ice2, jce1, jce2
    use mod_dynparam, only : ici1, ici2, jci1, jci2
!
    implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
    integer :: i, j
    real(rk8), parameter :: initval = 1.0d20
    real(rk8), parameter :: zeroval = 0.0d20
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
    call getmem2d(exportFields%snof,jci1,jci2,ici1,ici2,'cpl:snof')
    call getmem2d(exportFields%taux,jce1,jce2,ice1,ice2,'cpl:taux')
    call getmem2d(exportFields%tauy,jce1,jce2,ice1,ice2,'cpl:tauy')
    call getmem2d(exportFields%wspd,jce1,jce2,ice1,ice2,'cpl:wspd')
    call getmem2d(exportFields%nflx,jce1,jce2,ice1,ice2,'cpl:nflx')
    call getmem2d(exportFields%sflx,jce1,jce2,ice1,ice2,'cpl:sflx')
    call getmem2d(exportFields%snow,jce1,jce2,ice1,ice2,'cpl:snow')
    call getmem2d(exportFields%dswr,jce1,jce2,ice1,ice2,'cpl:dswr')
!
    call getmem2d(importFields%sst,jce1,jce2,ice1,ice2,'cpl:sst')
    call getmem2d(importFields%sit,jce1,jce2,ice1,ice2,'cpl:sit')
    call getmem2d(importFields%msk,jce1,jce2,ice1,ice2,'cpl:msk')
!
    call getmem2d(ldmskb,jci1,jci2,ici1,ici2,'cpl:ldmsk')
    call getmem2d(wetdry,jci1,jci2,ici1,ici2,'cpl:wetdry')
!
!-----------------------------------------------------------------------
!     Initialize arrays 
!-----------------------------------------------------------------------
!
    do i = ici1, ici2
      do j = jci1, jci2
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
        exportFields%rnof(j,i) = zeroval
        exportFields%snof(j,i) = zeroval
        exportFields%taux(j,i) = zeroval
        exportFields%tauy(j,i) = zeroval
        exportFields%wspd(j,i) = zeroval
        exportFields%nflx(j,i) = zeroval
        exportFields%sflx(j,i) = zeroval
        exportFields%snow(j,i) = zeroval
        exportFields%dswr(j,i) = zeroval
!
        importFields%sst(j,i) = initval
        importFields%sit(j,i) = initval
        importFields%msk(j,i) = initval 
!
        ldmskb(j,i) = mddom%ldmsk(j,i)
        wetdry(j,i) = 0
      end do
    end do
!
  end subroutine RCM_Allocate
!
  subroutine RCM_Get(localPet)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
    use mod_constants
    use mod_atm_interface, only : sfs, mddom , mdsub
    use mod_lm_interface, only : import_data_into_surface
!
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
    use mod_lm_interface, only : export_data_from_surface
!
    implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
    integer, intent(in) :: localPet
!
!-----------------------------------------------------------------------
!     Send information to OCN component
!-----------------------------------------------------------------------
!
    call export_data_from_surface(exportFields)

  end subroutine RCM_Put
!
end module mod_update
