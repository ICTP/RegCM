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

module mod_pbl_interface

  use mod_realkinds
  use mod_service
  use mod_constants
  use mod_dynparam
  use mod_memutil
  use mod_mppparam
  use mod_regcm_types
  use mod_pbl_common , only : ricr , chiuwten , uwstatea , uwstateb , kmxpbl
  use mod_pbl_holtbl , only : holtbl , allocate_mod_pbl_holtbl
  use mod_pbl_uwtcm , only : allocate_tcm_state
  use mod_pbl_uwtcm , only : uwtcm , get_data_from_tcm
  use mod_pbl_uwtcm , only : init_mod_pbl_uwtcm , tkemin
  use mod_runparams , only : ibltyp
  use mod_runparams , only : iqc , iqv , dt , rdt , ichem , hsigma , dsigma

  implicit none

  private

  type(mod_2_pbl) :: m2p
  type(pbl_2_mod) :: p2m

  public :: allocate_pblscheme
  public :: init_pblscheme
  public :: pblscheme

  public :: uwstatea
  public :: uwstateb
  public :: kmxpbl
  public :: ricr
  public :: tkemin

  contains

  subroutine allocate_pblscheme
    implicit none
    if ( ibltyp == 1 ) then
      call getmem2d(ricr,jci1,jci2,ici1,ici2,'pbl_common:ricr')
      call allocate_mod_pbl_holtbl
    end if
    if ( ibltyp == 2 ) then
      call allocate_tcm_state(uwstatea)
      call allocate_tcm_state(uwstateb)
      if ( ichem == 1 ) then
        call getmem4d(chiuwten,jci1,jci2,ici1,ici2, &
                      1,kz,1,ntr,'pbl_common:chiuwten')
      end if
      call init_mod_pbl_uwtcm
    end if
  end subroutine allocate_pblscheme

  subroutine init_pblscheme
    use mod_atm_interface
    use mod_che_interface
    use mod_che_common
    implicit none

    ! INPUT to PBL
    call assignpnt(mddom%coriol,m2p%coriol)
    call assignpnt(sfs%psdotb,m2p%psdot)
    call assignpnt(sfs%psb,m2p%psb)
    call assignpnt(sfs%tgb,m2p%tgb)
    call assignpnt(sfs%qfx,m2p%qfx)
    call assignpnt(sfs%hfx,m2p%hfx)
    call assignpnt(sfs%uvdrag,m2p%uvdrag)
    call assignpnt(atms%ubx3d,m2p%uxatm)
    call assignpnt(atms%vbx3d,m2p%vxatm)
    call assignpnt(atms%ubd3d,m2p%udatm)
    call assignpnt(atms%vbd3d,m2p%vdatm)
    call assignpnt(atms%tb3d,m2p%tatm)
    call assignpnt(atms%tp3d,m2p%tpatm)
    call assignpnt(atms%pb3d,m2p%patm)
    call assignpnt(atms%pf3d,m2p%patmf)
    call assignpnt(atms%qxb3d,m2p%qxatm)
    call assignpnt(atms%chib3d,m2p%chib)
    call assignpnt(atms%th3d,m2p%thatm)
    call assignpnt(atms%za,m2p%za)
    call assignpnt(atms%zq,m2p%zq)
    call assignpnt(atms%dzq,m2p%dzq)
    call assignpnt(atms%rhox2d,m2p%rhox2d)
    call assignpnt(atm2%tke,m2p%tkests)
    call assignpnt(heatrt,m2p%heatrt)
    call assignpnt(drydepv,m2p%drydepv)
    call assignpnt(chifxuw,m2p%chifxuw)
    call assignpnt(ktrop,m2p%ktrop)

    ! OUTPUT FROM PBL
    call assignpnt(aten%t,p2m%tten)
    call assignpnt(aten%u,p2m%uten)
    call assignpnt(aten%v,p2m%vten)
    call assignpnt(aten%qx,p2m%qxten)
    call assignpnt(aten%tke,p2m%tketen)
    call assignpnt(uwten%u,p2m%uuwten)
    call assignpnt(uwten%v,p2m%vuwten)
    call assignpnt(uwten%t,p2m%tuwten)
    call assignpnt(uwten%tke,p2m%tkeuwten)
    call assignpnt(uwten%qx,p2m%qxuwten)
    call assignpnt(adf%t,p2m%difft)
    call assignpnt(adf%qx,p2m%diffqx)
    call assignpnt(chiten,p2m%chiten)
    call assignpnt(remdrd,p2m%remdrd)
    call assignpnt(zpbl,p2m%zpbl)
    call assignpnt(kpbl,p2m%kpbl)

  end subroutine init_pblscheme

  subroutine pblscheme
    use mod_atm_interface
    implicit none
    select case ( ibltyp )
      case (1)
        call holtbl(m2p,p2m)
      case (2)
        call uwtcm(m2p,p2m)
        call uvcross2dot(uwten%u,uwten%v,aten%u,aten%v,1,kz)
        call get_data_from_tcm(p2m,atm1,atm2)
      case default
        return
    end select
  end subroutine pblscheme

end module mod_pbl_interface
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
