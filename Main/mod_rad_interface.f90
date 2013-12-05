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

module mod_rad_interface
!
  use mod_realkinds
  use mod_regcm_types
  use mod_rad_common
  use mod_rad_aerosol
  use mod_rad_colmod3
  use mod_rrtmg_driver
  use mod_rad_o3blk
  use mod_rad_outrad
  use mod_rad_radiation
  use mod_rad_scenarios
  use mod_rad_tracer

  private

  public :: init_rad , init_rad_clm , radiation

  public :: o3prof , gasabsnxt , gasabstot , gasemstot , taucldsp

  public :: allocate_mod_rad_common
  public :: allocate_mod_rad_aerosol
  public :: allocate_mod_rad_o3blk
  public :: allocate_mod_rad_outrad
  public :: allocate_mod_rad_colmod3
  public :: allocate_mod_rad_radiation
  public :: allocate_mod_rad_rrtmg
  public :: set_scenario , o3data , read_o3data

  type(mod_2_rad) :: m2r
  type(rad_2_mod) :: r2m

  contains

  subroutine init_rad
    use mod_atm_interface
    use mod_che_interface
    implicit none
    ! Set pipings from atm_interface to radiation I/O data types

    ! INPUT
    call assignpnt(atms%tb3d,m2r%tatms)
    call assignpnt(atms%qxb3d,m2r%qxatms)
    call assignpnt(atms%rhb3d,m2r%rhatms)
    call assignpnt(sfs%tgbb,m2r%tg)
    call assignpnt(sfs%psa,m2r%psa)
    call assignpnt(sfs%psb,m2r%psb)
    call assignpnt(mddom%xlat,m2r%xlat)
    call assignpnt(coszrs,m2r%coszrs)
    call assignpnt(aldirs,m2r%aldirs)
    call assignpnt(aldifs,m2r%aldifs)
    call assignpnt(aldirl,m2r%aldirl)
    call assignpnt(aldifl,m2r%aldifl)
    call assignpnt(albvs,m2r%albvs)
    call assignpnt(albvl,m2r%albvl)
    call assignpnt(emiss,m2r%emiss)
    call assignpnt(chia,m2r%chia)
    call assignpnt(ldmsk,m2r%ldmsk)
    call assignpnt(cldfra,m2r%cldfrc)
    call assignpnt(cldlwc,m2r%cldlwc)
    call assignpnt(ptrop,m2r%ptrop)
    ! OUTPUT
    call assignpnt(solis,r2m%solis)
    call assignpnt(sabveg,r2m%sabveg)
    call assignpnt(sinc,r2m%sinc)
    call assignpnt(solvs,r2m%solvs)
    call assignpnt(solvd,r2m%solvd)
    call assignpnt(fsw,r2m%fsw)
    call assignpnt(flw,r2m%flw)
    call assignpnt(flwd,r2m%flwd)
    call assignpnt(heatrt,r2m%heatrt)
  end subroutine init_rad
!
  subroutine init_rad_clm(sols,soll,solsd,solld)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:,:) :: sols
    real(rk8) , pointer , intent(in) , dimension(:,:) :: soll
    real(rk8) , pointer , intent(in) , dimension(:,:) :: solsd
    real(rk8) , pointer , intent(in) , dimension(:,:) :: solld
    call assignpnt(sols,r2m%sols)
    call assignpnt(soll,r2m%soll)
    call assignpnt(solsd,r2m%solsd)
    call assignpnt(solld,r2m%solld)
  end subroutine init_rad_clm

  subroutine radiation(iyear,eccf,loutrad,labsem)
    implicit none
    integer(ik4) , intent(in) :: iyear
    real(rk8) , intent(in) :: eccf
    logical , intent(in) :: loutrad
    logical , intent(in) :: labsem
    if ( irrtm == 1 ) then
      call rrtmg_driver(iyear,eccf,loutrad,m2r,r2m)
    else
      call colmod3(iyear,eccf,loutrad,labsem,m2r,r2m)
    end if
  end subroutine radiation

end module mod_rad_interface
