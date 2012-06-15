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
  use mod_atm_interface , only : atmstate , slice , surfstate , domain
  use mod_rad_common
  use mod_rad_aerosol
  use mod_rad_colmod3
  use mod_rrtmg_driver
  use mod_rad_o3blk
  use mod_rad_outrad
  use mod_rad_radiation
  use mod_rad_scenarios
  use mod_rad_tracer
!
  public
!
  contains 
!
  subroutine init_rad(ichem,ptop,a,sigma,twt,atms,sfs,mddom,sabveg, &
                      solis,coszrs,aldirs,aldifs,aldirl,aldifl,albvs, &
                      albvl,emiss,sinc,solvs,solvd,fsw,flw,flwd,ldmsk1, &
                      chia,chtrname)
    implicit none
    integer , intent(in) :: ichem
    real(dp) , intent(in) :: ptop
    real(dp) , pointer , dimension(:) :: a , sigma
    type(slice) , intent(in) :: atms
    type(surfstate) , intent(in) :: sfs
    type(domain) , intent(in) :: mddom
    real(dp) , pointer , intent(in) , dimension(:,:) :: sabveg
    real(dp) , pointer , intent(in) , dimension(:,:) :: solis
    real(dp) , pointer , intent(in) , dimension(:,:) :: coszrs
    real(dp) , pointer , intent(in) , dimension(:,:) :: aldirs
    real(dp) , pointer , intent(in) , dimension(:,:) :: aldifs
    real(dp) , pointer , intent(in) , dimension(:,:) :: aldirl
    real(dp) , pointer , intent(in) , dimension(:,:) :: aldifl
    real(dp) , pointer , intent(in) , dimension(:,:) :: albvs
    real(dp) , pointer , intent(in) , dimension(:,:) :: albvl
    real(dp) , pointer , intent(in) , dimension(:,:) :: emiss
    real(dp) , pointer , intent(in) , dimension(:,:) :: twt
    real(dp) , pointer , intent(in) , dimension(:,:) :: sinc
    real(dp) , pointer , intent(in) , dimension(:,:) :: solvs
    real(dp) , pointer , intent(in) , dimension(:,:) :: solvd
    real(dp) , pointer , intent(in) , dimension(:,:) :: fsw
    real(dp) , pointer , intent(in) , dimension(:,:) :: flw
    real(dp) , pointer , intent(in) , dimension(:,:) :: flwd
    integer , pointer , intent(in) , dimension(:,:,:) :: ldmsk1
    real(dp) , pointer , intent(in) , dimension(:,:,:,:) :: chia
    character(len=6) , pointer , intent(in) , dimension(:) :: chtrname

    if ( ichem == 1 ) lchem = .true.
    ptp = ptop
    call assignpnt(sigma,flev)
    call assignpnt(a,hlev)
    call assignpnt(twt,twtr)
    call assignpnt(atms%tb3d,tatms)
    call assignpnt(atms%qxb3d,qxatms)
    call assignpnt(atms%rhb3d,rhatms)
    call assignpnt(sfs%tgbb,tground)
    call assignpnt(sfs%psa,psfps)
    call assignpnt(sfs%psb,sfps)
    call assignpnt(mddom%xlat,xlat)
    call assignpnt(sabveg,abveg)
    call assignpnt(solis,solar)
    call assignpnt(coszrs,coszen)
    call assignpnt(aldirs,swdiralb)
    call assignpnt(aldifs,swdifalb)
    call assignpnt(aldirl,lwdiralb)
    call assignpnt(aldifl,lwdifalb)
    call assignpnt(albvs,swalb)
    call assignpnt(albvl,lwalb)
    call assignpnt(emiss,emsvt)
    call assignpnt(sinc,totsol)
    call assignpnt(solvs,soldir)
    call assignpnt(solvd,soldif)
    call assignpnt(fsw,srfabswflx)
    call assignpnt(flw,srflwflxup)
    call assignpnt(flwd,srflwflxdw)
    call assignpnt(ldmsk1,lndocnicemsk)
    call assignpnt(chia,chspmix)
    if ( associated(chtrname) ) tracname => chtrname
  end subroutine init_rad
!
  subroutine init_rad_clm(sols,soll,solsd,solld)
    implicit none
    real(dp) , pointer , intent(in) , dimension(:,:) :: sols
    real(dp) , pointer , intent(in) , dimension(:,:) :: soll
    real(dp) , pointer , intent(in) , dimension(:,:) :: solsd
    real(dp) , pointer , intent(in) , dimension(:,:) :: solld
    call assignpnt(sols,solswdir)
    call assignpnt(soll,sollwdir)
    call assignpnt(solsd,solswdif)
    call assignpnt(solld,sollwdif)
  end subroutine init_rad_clm

end module mod_rad_interface
