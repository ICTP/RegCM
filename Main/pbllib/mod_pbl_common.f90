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
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_memutil
!
  private

  real(dp) , public , pointer , dimension(:,:,:) :: zq
  real(dp) , public , pointer , dimension(:,:,:) :: za
  real(dp) , public , pointer , dimension(:,:,:) :: dzq
  real(dp) , public , pointer ,  dimension(:,:) :: rhox2d
!
  integer , public , pointer , dimension(:,:) :: kpbl
  real(dp) , public , pointer , dimension(:,:) :: zpbl
!
  type tcm_state
    !
    ! TKE*ps (m^2/s^2 * cb)
    !
    real(dp) , pointer , dimension(:,:,:) :: tkeps
    !
    ! Coupled TKE Advective Tendency (m^2/s^3 * cb) 
    !
    real(dp) , pointer , dimension(:,:,:) :: advtke
    !
    ! Vertical momentum diffusivity (m^2/s)
    !
    real(dp) , pointer , dimension(:,:,:) :: kzm
    !
    ! Vertical scalar diffusivity (m^2/s)
    !
    real(dp) , pointer , dimension(:,:,:) :: kth
    !
    ! Boundary layer height (m)
    !
    real(dp) , pointer , dimension(:,:) :: zpbl
    !
    ! Surface layer TKE (m^2/s^2)
    !
    real(dp) , pointer , dimension(:,:) :: srftke
    !
  end type tcm_state

  public :: tcm_state

  integer , public :: ibltyp
  integer , public :: kmxpbl

  !
  ! Pointers to point to the TCM's state variable
  !
  type(tcm_state) , public :: uwstatea , uwstateb

  real(dp) , public :: dtpbl ! dt
  real(dp) , public :: rdtpbl ! 1/dt
  real(dp) , public :: dttke ! TKE time step
  real(dp) , public :: tkemin
  real(dp) , public , pointer , dimension(:,:,:,:) :: chiuwten! chiuwten
  real(dp) , public , pointer , dimension(:,:,:) :: chifxuw   ! chifxuw

  !
  ! Specific instances of the model's state variables (at the b time step)
  !
  real(dp) , public , pointer , dimension(:,:,:) :: uten      ! aten%u
  real(dp) , public , pointer , dimension(:,:,:) :: vten      ! aten%v
  real(dp) , public , pointer , dimension(:,:,:) :: tten      ! aten%t
  real(dp) , public , pointer , dimension(:,:,:) :: tketen    ! aten%tke
  real(dp) , public , pointer , dimension(:,:,:) :: qvten     ! aten%qv
  real(dp) , public , pointer , dimension(:,:,:) :: qcten     ! aten%qc
  real(dp) , public , pointer , dimension(:,:,:) :: uuwten    ! uwten%u
  real(dp) , public , pointer , dimension(:,:,:) :: vuwten    ! uwten%v
  real(dp) , public , pointer , dimension(:,:,:) :: tuwten    ! uwten%t
  real(dp) , public , pointer , dimension(:,:,:) :: tkeuwten  ! uwten%tke
  real(dp) , public , pointer , dimension(:,:,:) :: qvuwten   ! uwten%qv
  real(dp) , public , pointer , dimension(:,:,:) :: qcuwten   ! uwten%qc
  real(dp) , public , pointer , dimension(:,:,:) :: uatm      ! atms%ubx3d
  real(dp) , public , pointer , dimension(:,:,:) :: vatm      ! atms%vbx3d
  real(dp) , public , pointer , dimension(:,:,:) :: udatm     ! atms%ubd3d
  real(dp) , public , pointer , dimension(:,:,:) :: vdatm     ! atms%vbd3d
  real(dp) , public , pointer , dimension(:,:,:) :: tatm      ! atms%tb3d
  real(dp) , public , pointer , dimension(:,:,:) :: qvatm     ! atms%qv
  real(dp) , public , pointer , dimension(:,:,:) :: qcatm     ! atms%qc
  real(dp) , public , pointer , dimension(:,:,:) :: tkests    ! atms%tke
  real(dp) , public , pointer , dimension(:,:,:) :: thxatm    ! atms%thx3d
  real(dp) , public , pointer , dimension(:,:,:) :: difft     ! adf%difft
  real(dp) , public , pointer , dimension(:,:,:) :: diffq     ! adf%diffq
  real(dp) , public , pointer , dimension(:,:,:) :: radheatrt ! heatrt
  real(dp) , public , pointer , dimension(:,:,:) :: diagqv    ! holtten%qv
  real(dp) , public , pointer , dimension(:,:,:) :: diagqc    ! holtten%qc
  real(dp) , public , pointer , dimension(:,:,:,:) :: chmx    ! chib
  real(dp) , public , pointer , dimension(:,:,:,:) :: chten   ! chiten
  real(dp) , public , pointer , dimension(:,:,:) :: drmr      ! remdrd
  real(dp) , public , pointer , dimension(:,:) :: sfcpd       ! psdot
  real(dp) , public , pointer , dimension(:,:) :: sfcps       ! sfs%psb
  real(dp) , public , pointer , dimension(:,:) :: tg          ! sfs%tgb
  real(dp) , public , pointer , dimension(:,:) :: qfx         ! sfs%qfx
  real(dp) , public , pointer , dimension(:,:) :: hfx         ! sfs%hfx
  real(dp) , public , pointer , dimension(:,:) :: uvdrag      ! sfs%uvdrag
  real(dp) , public , pointer , dimension(:,:) :: coriolis    ! mddom%coriol
  real(dp) , public , pointer , dimension(:,:) :: mapfcx      ! mddom%msfx
  integer , public , pointer , dimension(:,:) :: landmsk      ! ldmsk
  real(dp) , public , pointer , dimension(:) :: hlev          ! a
  real(dp) , public , pointer , dimension(:) :: flev          ! sigma
  real(dp) , public , pointer , dimension(:) :: dlev          ! dsigma
  real(dp) , public :: ptp                                    ! ptop

  real(dp) , public , pointer , dimension(:,:) :: depvel       ! chtrdpv
  character(len=5) , public , pointer , dimension(:) :: chname ! chtrname
  logical , public :: lchem , lchdrydepo

  real(dp) , public , pointer , dimension(:,:,:) :: dotqdot , ftmp

  data lchem /.false./
  data lchdrydepo /.false./

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
      call getmem3d(tcmstate%kzm,jce1,jce2,ice1,ice2,1,kzp1,'pbl_common:kzm')
      call getmem3d(tcmstate%kth,jce1,jce2,ice1,ice2,1,kzp1,'pbl_common:kth')
      call getmem2d(tcmstate%zpbl,jci1,jci2,ici1,ici2,'pbl_common:zpbl')
      call getmem2d(tcmstate%srftke,jci1,jci2,ici1,ici2,'pbl_common:srftke')
    else
      call getmem3d(tcmstate%kzm,jcross1,jcross2, &
                                 icross1,icross2,1,kzp1,'pbl_common:kzm')
      call getmem3d(tcmstate%kth,jcross1,jcross2, &
                                 icross1,icross2,1,kzp1,'pbl_common:kth')
    end if
  end subroutine allocate_tcm_state

  subroutine allocate_mod_pbl_common(ibltyp,ichem)
    implicit none
    integer , intent(in) :: ibltyp
    integer , intent(in) :: ichem
    call getmem3d(zq,jce1,jce2,ice1,ice2,1,kzp1,'pbl_common:zq')
    call getmem3d(za,jce1,jce2,ice1,ice2,1,kz,'pbl_common:za')
    call getmem3d(dzq,jce1,jce2,ice1,ice2,1,kz,'pbl_common:dzq')
    call getmem2d(rhox2d,jce1,jce2,ice1,ice2,'pbl_common:rhox2d')
    call getmem2d(kpbl,jce1,jce2,ice1,ice2,'pbl_common:kpbl')
    call getmem2d(zpbl,jce1,jce2,ice1,ice2,'pbl_common:zpbl')
    !
    ! Allocate the tcm state variables
    !
    if ( ibltyp == 2 .or. ibltyp == 99) then
      call allocate_tcm_state(uwstatea,.true.)
      call allocate_tcm_state(uwstateb,.true.)
      ! To be used in vertical advection scheme
      call getmem3d(dotqdot,jci1,jci2,ici1,ici2,1,kz,'mod_uwtcm:dotqdot')
      call getmem3d(ftmp,jci1,jci2,ici1,ici2,1,kz,'mod_uwtcm:ftmp')
      if(ichem == 1)then
        lchem = .true. 
        call getmem4d(chiuwten,jci1,jci2,ici1,ici2, &
                               1,kz,1,ntr,'pbl_common:chiuwten')
!       call getmem3d(chifxuw,jci1,jci2,ici1,ici2,1,ntr,'pbl_common:chifxuw')
      end if
    end if
  end subroutine allocate_mod_pbl_common
!
end module mod_pbl_common
