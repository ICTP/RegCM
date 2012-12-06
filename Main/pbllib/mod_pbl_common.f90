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

module mod_pbl_common
!
! Storage parameters and constants related to the boundary layer
!
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_memutil
  use mod_runparams , only : ibltyp , ichem
!
  private

  real(rk8) , public , pointer ,  dimension(:,:) :: rhox2d
!
  integer(ik4) , public , pointer , dimension(:,:) :: kpbl
  real(rk8) , public , pointer , dimension(:,:) :: zpbl
!
  type tcm_state
    !
    ! TKE*ps (m^2/s^2 * cb)
    !
    real(rk8) , pointer , dimension(:,:,:) :: tkeps
    !
    ! Coupled TKE Advective Tendency (m^2/s^3 * cb) 
    !
    real(rk8) , pointer , dimension(:,:,:) :: advtke
    !
    ! Vertical momentum diffusivity (m^2/s)
    !
    real(rk8) , pointer , dimension(:,:,:) :: kzm
    !
    ! Vertical scalar diffusivity (m^2/s)
    !
    real(rk8) , pointer , dimension(:,:,:) :: kth
    !
    ! Boundary layer height (m)
    !
    real(rk8) , pointer , dimension(:,:) :: zpbl
    !
    ! Surface layer TKE (m^2/s^2)
    !
    real(rk8) , pointer , dimension(:,:) :: srftke
    !
  end type tcm_state

  public :: tcm_state

  integer(ik4) , public :: kmxpbl

  !
  ! Pointers to point to the TCM's state variable
  !
  type(tcm_state) , public :: uwstatea , uwstateb

  real(rk8) , public :: tkemin
  real(rk8) , public , pointer , dimension(:,:,:,:) :: chiuwten
  real(rk8) , public , pointer , dimension(:,:,:) :: chifxuw
  !
  ! Specific instances of the model's state variables (at the b time step)
  !
  real(rk8) , public , pointer , dimension(:,:,:) :: uten      ! aten%u
  real(rk8) , public , pointer , dimension(:,:,:) :: vten      ! aten%v
  real(rk8) , public , pointer , dimension(:,:,:) :: tten      ! aten%t
  real(rk8) , public , pointer , dimension(:,:,:) :: tketen    ! aten%tke
  real(rk8) , public , pointer , dimension(:,:,:,:) :: qxten   ! aten%qx
  real(rk8) , public , pointer , dimension(:,:,:) :: uuwten    ! uwten%u
  real(rk8) , public , pointer , dimension(:,:,:) :: vuwten    ! uwten%v
  real(rk8) , public , pointer , dimension(:,:,:) :: tuwten    ! uwten%t
  real(rk8) , public , pointer , dimension(:,:,:) :: tkeuwten  ! uwten%tke
  real(rk8) , public , pointer , dimension(:,:,:,:) :: qxuwten ! uwten%qx
  real(rk8) , public , pointer , dimension(:,:,:) :: uxatm     ! atms%ubx3d
  real(rk8) , public , pointer , dimension(:,:,:) :: vxatm     ! atms%vbx3d
  real(rk8) , public , pointer , dimension(:,:,:) :: udatm     ! atms%ubd3d
  real(rk8) , public , pointer , dimension(:,:,:) :: vdatm     ! atms%vbd3d
  real(rk8) , public , pointer , dimension(:,:,:) :: tatm      ! atms%tb3d
  real(rk8) , public , pointer , dimension(:,:,:,:) :: qxatm   ! atms%qx
  real(rk8) , public , pointer , dimension(:,:,:) :: tkests    ! atms%tke
  real(rk8) , public , pointer , dimension(:,:,:) :: thxatm    ! atms%thx3d
  real(rk8) , public , pointer , dimension(:,:,:) :: zq        ! atms%zq
  real(rk8) , public , pointer , dimension(:,:,:) :: za        ! atms%za
  real(rk8) , public , pointer , dimension(:,:,:) :: dzq       ! atms%dzq
  real(rk8) , public , pointer , dimension(:,:,:) :: difft     ! adf%difft
  real(rk8) , public , pointer , dimension(:,:,:,:) :: diffqx  ! adf%diffqx
  real(rk8) , public , pointer , dimension(:,:,:) :: radheatrt ! heatrt
  real(rk8) , public , pointer , dimension(:,:,:,:) :: diagqx  ! holtten%qx
  real(rk8) , public , pointer , dimension(:,:,:,:) :: chmx    ! chib
  real(rk8) , public , pointer , dimension(:,:,:,:) :: chten   ! chiten
  real(rk8) , public , pointer , dimension(:,:,:) :: drmr      ! remdrd
  real(rk8) , public , pointer , dimension(:,:) :: sfcpd       ! psdot
  real(rk8) , public , pointer , dimension(:,:) :: sfcps       ! sfs%psb
  real(rk8) , public , pointer , dimension(:,:) :: tg          ! sfs%tgb
  real(rk8) , public , pointer , dimension(:,:) :: qfx         ! sfs%qfx
  real(rk8) , public , pointer , dimension(:,:) :: hfx         ! sfs%hfx
  real(rk8) , public , pointer , dimension(:,:) :: uvdrag      ! sfs%uvdrag
  real(rk8) , public , pointer , dimension(:,:) :: coriolis    ! mddom%coriol
  real(rk8) , public , pointer , dimension(:,:) :: mapfcx      ! mddom%msfx
  integer(ik4) , public , pointer , dimension(:,:) :: landmsk  ! ldmsk
  real(rk8) , public , pointer , dimension(:) :: hlev          ! a
  real(rk8) , public , pointer , dimension(:) :: flev          ! sigma
  real(rk8) , public , pointer , dimension(:) :: dlev          ! dsigma
  real(rk8) , public , pointer , dimension(:,:) :: depvel      ! chtrdpv

  real(rk8) , public , pointer , dimension(:,:,:) :: dotqdot , ftmp

  public :: allocate_mod_pbl_common , allocate_tcm_state

  contains

  subroutine allocate_tcm_state(tcmstate,lpar)
    implicit none
    type(tcm_state) , intent(out) :: tcmstate
    logical , intent(in) :: lpar
    if ( lpar ) then
      call getmem3d(tcmstate%tkeps,jce1,jce2,ice1,ice2,1,kzp1,'pbl_common:tkeps')
      call getmem3d(tcmstate%advtke,jce1,jce2, &
                                    ice1,ice2,1,kzp1,'pbl_common:advtke')
      call getmem3d(tcmstate%kzm,jci1,jci2,ici1,ici2,1,kzp1,'pbl_common:kzm')
      call getmem3d(tcmstate%kth,jci1,jci2,ici1,ici2,1,kzp1,'pbl_common:kth')
      call getmem2d(tcmstate%zpbl,jci1,jci2,ici1,ici2,'pbl_common:zpbl')
      call getmem2d(tcmstate%srftke,jci1,jci2,ici1,ici2,'pbl_common:srftke')
    else
      call getmem3d(tcmstate%kzm,jcross1,jcross2, &
                                 icross1,icross2,1,kzp1,'pbl_common:kzm')
      call getmem3d(tcmstate%kth,jcross1,jcross2, &
                                 icross1,icross2,1,kzp1,'pbl_common:kth')
    end if
  end subroutine allocate_tcm_state

  subroutine allocate_mod_pbl_common
    implicit none
    call getmem2d(kpbl,jci1,jci2,ici1,ici2,'pbl_common:kpbl')
    call getmem2d(zpbl,jci1,jci2,ici1,ici2,'pbl_common:zpbl')
    !
    ! Allocate the tcm state variables
    !
    if ( ibltyp == 2 .or. ibltyp == 99) then
      call allocate_tcm_state(uwstatea,.true.)
      call allocate_tcm_state(uwstateb,.true.)
      ! To be used in vertical advection scheme
      call getmem3d(dotqdot,jci1,jci2,ici1,ici2,1,kz,'mod_uwtcm:dotqdot')
      call getmem3d(ftmp,jci1,jci2,ici1,ici2,1,kz,'mod_uwtcm:ftmp')
      if ( ichem == 1 ) then
        call getmem4d(chiuwten,jci1,jci2,ici1,ici2, &
                               1,kz,1,ntr,'pbl_common:chiuwten')
        call getmem3d(chifxuw,jci1,jci2,ici1,ici2,1,ntr,'pbl_common:chifxuw')
      end if
    end if
  end subroutine allocate_mod_pbl_common
!
end module mod_pbl_common
