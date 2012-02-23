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
  public

  public :: tcm_state

  real(dp) , pointer , dimension(:,:,:) :: zq
  real(dp) , pointer , dimension(:,:,:) :: za
  real(dp) , pointer , dimension(:,:,:) :: dzq
  real(dp) , pointer ,  dimension(:,:) :: rhox2d
!
  public :: kpbl , zpbl
  integer , pointer , dimension(:,:) :: kpbl
  real(dp) , pointer , dimension(:,:) :: zpbl
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

  integer :: ibltyp
  integer :: kmxpbl
  integer :: itcmstart , itcmend
  integer :: jtcmstart , jtcmend

  !
  ! Pointers to point to the TCM's state variable
  !
  type(tcm_state) :: uwstatea , uwstateb

  real(dp) :: dtpbl ! dt
  real(dp) :: rdtpbl ! 1/dt
  real(dp) :: dttke ! TKE time step
  real(dp) :: tkemin
  real(dp) , pointer , dimension(:,:,:,:) :: chiuwten! chiuwten
  real(dp) , pointer , dimension(:,:,:) :: chifxuw   ! chifxuw

  !
  ! Specific instances of the model's state variables (at the b time step)
  !
  real(dp) , pointer , dimension(:,:,:) :: uten      ! aten%u
  real(dp) , pointer , dimension(:,:,:) :: vten      ! aten%v
  real(dp) , pointer , dimension(:,:,:) :: tten      ! aten%t
  real(dp) , pointer , dimension(:,:,:) :: tketen    ! aten%tke
  real(dp) , pointer , dimension(:,:,:) :: qvten     ! aten%qv
  real(dp) , pointer , dimension(:,:,:) :: qcten     ! aten%qc
  real(dp) , pointer , dimension(:,:,:) :: uuwten    ! uwten%u
  real(dp) , pointer , dimension(:,:,:) :: vuwten    ! uwten%v
  real(dp) , pointer , dimension(:,:,:) :: tuwten    ! uwten%t
  real(dp) , pointer , dimension(:,:,:) :: tkeuwten  ! uwten%tke
  real(dp) , pointer , dimension(:,:,:) :: qvuwten   ! uwten%qv
  real(dp) , pointer , dimension(:,:,:) :: qcuwten   ! uwten%qc
  real(dp) , pointer , dimension(:,:,:) :: uatm      ! atms%ubx3d
  real(dp) , pointer , dimension(:,:,:) :: vatm      ! atms%vbx3d
  real(dp) , pointer , dimension(:,:,:) :: udatm     ! atms%ubd3d
  real(dp) , pointer , dimension(:,:,:) :: vdatm     ! atms%vbd3d
  real(dp) , pointer , dimension(:,:,:) :: tatm      ! atms%tb3d
  real(dp) , pointer , dimension(:,:,:) :: qvatm     ! atms%qv
  real(dp) , pointer , dimension(:,:,:) :: qcatm     ! atms%qc
  real(dp) , pointer , dimension(:,:,:) :: tkests    ! atms%tke
  real(dp) , pointer , dimension(:,:,:) :: thxatm    ! atms%thx3d
  real(dp) , pointer , dimension(:,:,:) :: difft     ! adf%difft
  real(dp) , pointer , dimension(:,:,:) :: diffq     ! adf%diffq
  real(dp) , pointer , dimension(:,:,:) :: radheatrt ! heatrt
  real(dp) , pointer , dimension(:,:,:) :: diagqv    ! holtten%qv
  real(dp) , pointer , dimension(:,:,:) :: diagqc    ! holtten%qc
  real(dp) , pointer , dimension(:,:,:,:) :: chmx    ! chib
  real(dp) , pointer , dimension(:,:,:,:) :: chten   ! chiten
  real(dp) , pointer , dimension(:,:,:) :: drmr      ! remdrd
  real(dp) , pointer , dimension(:,:) :: sfcpd       ! psdot
  real(dp) , pointer , dimension(:,:) :: sfcps       ! sfs%psb
  real(dp) , pointer , dimension(:,:) :: tg          ! sfs%tgb
  real(dp) , pointer , dimension(:,:) :: qfx         ! sfs%qfx
  real(dp) , pointer , dimension(:,:) :: hfx         ! sfs%hfx
  real(dp) , pointer , dimension(:,:) :: uvdrag      ! sfs%uvdrag
  real(dp) , pointer , dimension(:,:) :: coriolis    ! mddom%coriol
  real(dp) , pointer , dimension(:,:) :: mapfcx      ! mddom%msfx
  integer , pointer , dimension(:,:) :: landmsk      ! ldmsk
  real(dp) , pointer , dimension(:) :: hlev          ! a
  real(dp) , pointer , dimension(:) :: flev          ! sigma
  real(dp) , pointer , dimension(:) :: dlev          ! dsigma
  real(dp) :: ptp                                    ! ptop

  real(dp) , pointer , dimension(:,:) :: depvel       ! chtrdpv
  character(len=5) , pointer , dimension(:) :: chname ! chtrname
  logical :: lchem , lchdrydepo

  data lchem /.false./
  data lchdrydepo /.false./

  contains

  subroutine allocate_tcm_state(tcmstate,lpar)
    implicit none
    type(tcm_state) , intent(out) :: tcmstate
    logical , intent(in) :: lpar
    if ( lpar ) then
      call getmem3d(tcmstate%tkeps,1,jxp,1,iy,1,kzp1,'pbl_common:tkeps')
      call getmem3d(tcmstate%advtke,1,jxp,1,iy,1,kzp1,'pbl_common:advtke')
      call getmem3d(tcmstate%kzm,1,jxp,1,iy,1,kzp1,'pbl_common:kzm')
      call getmem3d(tcmstate%kth,1,jxp,1,iy,1,kzp1,'pbl_common:kth')
      call getmem2d(tcmstate%zpbl,1,jxp,1,iy,'pbl_common:zpbl')
      call getmem2d(tcmstate%srftke,1,jxp,1,iy,'pbl_common:srftke')
    else
      call getmem3d(tcmstate%kzm,1,jx,1,iy,1,kzp1,'pbl_common:kzm')
      call getmem3d(tcmstate%kth,1,jx,1,iy,1,kzp1,'pbl_common:kth')
    end if
  end subroutine allocate_tcm_state

  subroutine allocate_mod_pbl_common(ibltyp,ichem)
    implicit none
    integer , intent(in) :: ibltyp
    integer , intent(in) :: ichem
    call getmem3d(zq,1,jxp,1,iy,1,kzp1,'pbl_common:zq')
    call getmem3d(za,1,jxp,1,iy,1,kz,'pbl_common:za')
    call getmem3d(dzq,1,jxp,1,iy,1,kz,'pbl_common:dzq')
    call getmem2d(rhox2d,1,jxp,1,iy,'pbl_common:rhox2d')
    call getmem2d(kpbl,1,jxp,1,iy,'pbl_common:kpbl')
    call getmem2d(zpbl,1,jxp,1,iy,'pbl_common:zpbl')

    !
    ! Allocate the tcm state variables
    !
    if ( ibltyp == 2 .or. ibltyp == 99) then
      call allocate_tcm_state(uwstatea,.true.)
      call allocate_tcm_state(uwstateb,.true.)
      if(ichem == 1)then
        call getmem4d(chiuwten,1,jxp,1,iy,1,kz,1,ntr,'pbl_common:chiuwten')
!        call getmem3d(chifxuw,1,jxp,1,iy,1,ntr,'pbl_common:chifxuw')
      end if
    end if
  end subroutine allocate_mod_pbl_common
!
end module mod_pbl_common
