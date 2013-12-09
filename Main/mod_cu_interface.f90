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

!
! Link atmospheric model and cumulus schemes
!
module mod_cu_interface
  use mod_realkinds
  use mod_intkinds
  use mod_runparams
  use mod_regcm_types

  use mod_cu_common , only : cuscheme , total_precip_points , cevapu ,     &
      model_cumulus_cloud , init_mod_cumulus
  use mod_cu_tiedtke , only : allocate_mod_cu_tiedtke , tiedtkedrv , q_detr
  use mod_cu_tables , only : init_convect_tables
  use mod_cu_bm , only : allocate_mod_cu_bm , bmpara , lutbl , cldefi ,    &
      tbase
  use mod_cu_em , only : allocate_mod_cu_em , cupemandrv , cbmf2d ,        &
      elcrit2d , epmax2d
  use mod_cu_kuo , only : allocate_mod_cu_kuo , cupara , htdiff , rsheat , &
      rswat , twght , vqflx , k700
  use mod_cu_grell , only : allocate_mod_cu_grell , cuparan , mincld2d ,   &
      shrmax2d , shrmin2d , edtmax2d , edtmin2d ,  edtmaxo2d , edtmaxx2d , &
      edtmino2d , edtminx2d , htmax2d , htmin2d , pbcmax2d , dtauc2d ,     &
      pbcmax2d , kbmax2d

  implicit none

  private

  public :: allocate_cumulus
  public :: init_cuscheme
  public :: cumulus

  public :: htdiff
  public :: lutbl

  public :: cuscheme
  public :: cbmf2d
  public :: cldefi
  public :: rsheat
  public :: rswat
  public :: tbase
  public :: q_detr
  public :: twght
  public :: vqflx
  public :: shrmax2d
  public :: shrmin2d
  public :: edtmax2d
  public :: edtmin2d
  public :: edtmaxo2d
  public :: edtmino2d
  public :: edtmaxx2d
  public :: edtminx2d
  public :: pbcmax2d
  public :: mincld2d
  public :: kbmax2d
  public :: htmax2d
  public :: htmin2d
  public :: dtauc2d
  public :: elcrit2d
  public :: epmax2d
  public :: cevapu
  public :: k700
  public :: total_precip_points

  type(mod_2_cum) :: m2c
  type(cum_2_mod) :: c2m

  contains

  subroutine allocate_cumulus
    implicit none
    if ( icup > 90 ) then
      call getmem2d(cuscheme,jci1,jci2,ici1,ici2,'cumulus:cuscheme')
     end if
    select case ( icup )
      case (1)
        call allocate_mod_cu_kuo
      case (2)
        call allocate_mod_cu_grell
      case (3)
        call allocate_mod_cu_bm
      case (4)
        call allocate_mod_cu_em
      case (5)
        call init_convect_tables
        call allocate_mod_cu_tiedtke
      case (96)
        call init_convect_tables
        call allocate_mod_cu_tiedtke
        call allocate_mod_cu_grell
      case (97)
        call init_convect_tables
        call allocate_mod_cu_tiedtke
        call allocate_mod_cu_em
      case (98)
        call allocate_mod_cu_grell
        call allocate_mod_cu_em
      case (99)
        call allocate_mod_cu_grell
        call allocate_mod_cu_em
      case default
        return
    end select
  end subroutine allocate_cumulus

  subroutine init_cuscheme
    use mod_atm_interface
    use mod_che_interface
    implicit none
    ! INPUT
    call assignpnt(mddom%ht,m2c%ht)
    call assignpnt(ldmsk,m2c%ldmsk)
    call assignpnt(sfs%psb,m2c%psb)
    call assignpnt(atms%za,m2c%zas)
    call assignpnt(atms%zq,m2c%zfs)
    call assignpnt(atms%pb3d,m2c%pas)
    call assignpnt(atms%tb3d,m2c%tas)
    call assignpnt(atms%ubx3d,m2c%uas)
    call assignpnt(atms%vbx3d,m2c%vas)
    call assignpnt(atms%qsb3d,m2c%qsas)
    call assignpnt(atms%qxb3d,m2c%qxas)
    call assignpnt(atms%chib3d,m2c%chias)
    call assignpnt(qdot,m2c%qdot)
    call assignpnt(sfs%qfx,m2c%qfx)
    call assignpnt(sfs%hfx,m2c%hfx)
    call assignpnt(ktrop,m2c%ktrop)
    ! OUTPUT
    call assignpnt(aten%t,c2m%tten)
    call assignpnt(aten%u,c2m%uten)
    call assignpnt(aten%v,c2m%vten)
    call assignpnt(aten%qx,c2m%qxten)
    call assignpnt(chiten,c2m%chiten)
    call assignpnt(sfs%rainc,c2m%rainc)
    call assignpnt(pptc,c2m%pcratec)
    call assignpnt(cldfra,c2m%cldfrc)
    call assignpnt(cldlwc,c2m%cldlwc)
    call assignpnt(icumtop,c2m%kcumtop)
    call assignpnt(icumbot,c2m%kcumbot)
    call assignpnt(convpr,c2m%convpr)
    call init_mod_cumulus
  end subroutine init_cuscheme
!
  subroutine cumulus
    implicit none
    select case ( icup )
      case (1)
        call cupara(m2c,c2m)
      case (2)
        call cuparan(m2c,c2m)
      case (3)
        call bmpara(m2c,c2m)
      case (4)
        call cupemandrv(m2c,c2m)
      case (5)
        call tiedtkedrv(m2c,c2m)
      case (96)
        call tiedtkedrv(m2c,c2m)
        call cuparan(m2c,c2m)
      case (97)
        call tiedtkedrv(m2c,c2m)
        call cupemandrv(m2c,c2m)
      case (98)
        call cuparan(m2c,c2m)
        call cupemandrv(m2c,c2m)
      case (99)
        call cuparan(m2c,c2m)
        call cupemandrv(m2c,c2m)
      case default
        return
    end select
    call model_cumulus_cloud(m2c,c2m)
  end subroutine cumulus

end module mod_cu_interface
