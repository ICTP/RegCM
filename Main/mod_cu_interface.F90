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
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_runparams
  use mod_memutil
  use mod_regcm_types

  use mod_cu_common , only : cuscheme , total_precip_points , cevapu ,     &
      model_cumulus_cloud , init_mod_cumulus , rain_cc
  use mod_cu_tiedtke , only : allocate_mod_cu_tiedtke , tiedtkedrv
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
  use mod_cu_kf

  implicit none

  private

  public :: allocate_cumulus
  public :: init_cumulus
  public :: cumulus

  public :: htdiff
  public :: lutbl

  public :: cuscheme
  public :: cbmf2d
  public :: cldefi
  public :: rsheat
  public :: rswat
  public :: tbase
  public :: rain_cc
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
    use mod_atm_interface
    implicit none
    integer(ik4) :: i , j
    call getmem3d(rain_cc,jci1,jci2,ici1,ici2,1,kz+1,'cumulus:rain_cc')
    call getmem2d(cuscheme,jci1,jci2,ici1,ici2,'cumulus:cuscheme')
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( mddom%lndcat(j,i) > 14.5D0 .and. &
             mddom%lndcat(j,i) < 15.5D0 ) then
          cuscheme(j,i) = icup_ocn
        else
          cuscheme(j,i) = icup_lnd
        end if
      end do
    end do
    if ( any(icup == 1) ) then
      call allocate_mod_cu_kuo
    end if
    if ( any(icup == 2) ) then
      call allocate_mod_cu_grell
    end if
    if ( any(icup == 3) ) then
      call allocate_mod_cu_bm
    end if
    if ( any(icup == 4) ) then
      call allocate_mod_cu_em
    end if
    if ( any(icup == 5) ) then
      if ( iconv /= 4 ) call init_convect_tables
      call allocate_mod_cu_tiedtke
    end if
    if ( any(icup == 6) ) then
      call allocate_mod_cu_kf
      call kf_lutab
    end if
  end subroutine allocate_cumulus

  subroutine init_cumulus
    use mod_atm_interface
    use mod_che_interface
    implicit none
    ! INPUT
    call assignpnt(mddom%ht,m2c%ht)
    call assignpnt(mddom%ldmsk,m2c%ldmsk)
    call assignpnt(sfs%psb,m2c%psb)
    call assignpnt(atms%za,m2c%zas)
    call assignpnt(atms%zq,m2c%zfs)
    call assignpnt(atms%dzq,m2c%dzq)
    call assignpnt(atms%ps2d,m2c%psf)
    call assignpnt(atms%pb3d,m2c%pas)
    call assignpnt(atms%pf3d,m2c%pasf)
    call assignpnt(atms%tb3d,m2c%tas)
    call assignpnt(atms%ubx3d,m2c%uas)
    call assignpnt(atms%vbx3d,m2c%vas)
    call assignpnt(atms%wpx3d,m2c%wpas)
    call assignpnt(atms%qsb3d,m2c%qsas)
    call assignpnt(atms%qxb3d,m2c%qxas)
    call assignpnt(atms%rhob3d,m2c%rhoas)
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
    call assignpnt(q_detr,c2m%q_detr)
    call init_mod_cumulus
  end subroutine init_cumulus

  subroutine cumulus
    implicit none
    if ( icup_lnd == icup_ocn ) then
      select case ( icup_lnd )
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
        case (6)
          call kfdrv(m2c,c2m)
      end select
    else
      select case ( icup_lnd )
        case (2)
          call cuparan(m2c,c2m)
        case (4)
          call cupemandrv(m2c,c2m)
        case (5)
          call tiedtkedrv(m2c,c2m)
        case (6)
          call kfdrv(m2c,c2m)
      end select
      select case ( icup_ocn )
        case (2)
          call cuparan(m2c,c2m)
        case (4)
          call cupemandrv(m2c,c2m)
        case (5)
          call tiedtkedrv(m2c,c2m)
        case (6)
          call kfdrv(m2c,c2m)
      end select
    end if
    call model_cumulus_cloud(m2c,c2m)
  end subroutine cumulus

end module mod_cu_interface
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
