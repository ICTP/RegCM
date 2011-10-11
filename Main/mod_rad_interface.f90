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
  use mod_atm_interface , only : atmstate , slice , surfpstate , surfstate , &
                                 domain
  use mod_rad_common
  use mod_rad_aerosol
  use mod_rad_colmod3
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
  subroutine init_rad(ichem,ptop,a,sigma,sps1,sps2,atms,sfsta, &
                      mddom,sabveg,solis,solvs,solvd,coszrs,aldirs,    &
                      aldifs,aldirl,aldifl,albdir,albdif,albvs,albvl,  &
                      emiss,sabv2d,sol2d,sinc2d,solvs2d,solvd2d,     &
                      fsw2d,flw2d,flwd2d,ocld2d,chia,chtrname)
    implicit none
    integer , intent(in) :: ichem
    real(8) , intent(in) :: ptop
    real(8) , pointer , dimension(:) :: a , sigma
    type(surfpstate) , intent(in) :: sps1 , sps2
    type(slice) , intent(in) :: atms
    type(surfstate) , intent(in) :: sfsta
    type(domain) , intent(in) :: mddom
    real(8) , pointer , intent(in) , dimension(:) :: sabveg
    real(8) , pointer , intent(in) , dimension(:) :: solis
    real(8) , pointer , intent(in) , dimension(:) :: solvs
    real(8) , pointer , intent(in) , dimension(:) :: solvd
    real(8) , pointer , intent(in) , dimension(:) :: coszrs
    real(8) , pointer , intent(in) , dimension(:) :: aldirs
    real(8) , pointer , intent(in) , dimension(:) :: aldifs
    real(8) , pointer , intent(in) , dimension(:) :: aldirl
    real(8) , pointer , intent(in) , dimension(:) :: aldifl
    real(8) , pointer , intent(in) , dimension(:) :: albdir
    real(8) , pointer , intent(in) , dimension(:) :: albdif
    real(8) , pointer , intent(in) , dimension(:) :: albvs
    real(8) , pointer , intent(in) , dimension(:) :: albvl
    real(8) , pointer , intent(in) , dimension(:) :: emiss
    real(8) , pointer , intent(in) , dimension(:,:) :: sabv2d
    real(8) , pointer , intent(in) , dimension(:,:) :: sol2d
    real(8) , pointer , intent(in) , dimension(:,:) :: sinc2d
    real(8) , pointer , intent(in) , dimension(:,:) :: solvs2d
    real(8) , pointer , intent(in) , dimension(:,:) :: solvd2d
    real(8) , pointer , intent(in) , dimension(:,:) :: fsw2d
    real(8) , pointer , intent(in) , dimension(:,:) :: flw2d
    real(8) , pointer , intent(in) , dimension(:,:) :: flwd2d
    integer , pointer , intent(in) , dimension(:,:,:) :: ocld2d
    real(8) , pointer , intent(in) , dimension(:,:,:,:) :: chia
    character(len=5) , pointer , intent(in) , dimension(:) :: chtrname

    if ( ichem == 1 ) lchem = .true.
    ptp = ptop
    call assignpnt(sigma,flev)
    call assignpnt(a,hlev)
    call assignpnt(sps1%ps,psfps)
    call assignpnt(sps2%ps,sfps)
    call assignpnt(atms%tb3d,tatms)
    call assignpnt(atms%qvb3d,qvatms)
    call assignpnt(atms%rhb3d,rhatms)
    call assignpnt(sfsta%tgbb,tground)
    call assignpnt(mddom%xlat,xlat)
    call assignpnt(sabveg,abveg)
    call assignpnt(solis,solar)
    call assignpnt(solvs,soldir)
    call assignpnt(solvd,soldif)
    call assignpnt(coszrs,coszen)
    call assignpnt(aldirs,swdiralb)
    call assignpnt(aldifs,swdifalb)
    call assignpnt(aldirl,lwdiralb)
    call assignpnt(aldifl,lwdifalb)
    call assignpnt(albdir,diralb)
    call assignpnt(albdif,difalb)
    call assignpnt(albvs,swalb)
    call assignpnt(albvl,lwalb)
    call assignpnt(emiss,emsvt)
    call assignpnt(sabv2d,abveg2d)
    call assignpnt(sol2d,solar2d)
    call assignpnt(sinc2d,totsol2d)
    call assignpnt(solvs2d,soldir2d)
    call assignpnt(solvd2d,soldif2d)
    call assignpnt(fsw2d,srfabswflx)
    call assignpnt(flw2d,srflwflxup)
    call assignpnt(flwd2d,srflwflxdw)
    call assignpnt(ocld2d,lndocnicemsk)
    call assignpnt(chia,chspmix)
    if ( associated(chtrname) ) tracname => chtrname
  end subroutine init_rad
!
  subroutine init_rad_clm(sols2d,soll2d,solsd2d,solld2d)
    implicit none
    real(8) , pointer , intent(in) , dimension(:,:) :: sols2d
    real(8) , pointer , intent(in) , dimension(:,:) :: soll2d
    real(8) , pointer , intent(in) , dimension(:,:) :: solsd2d
    real(8) , pointer , intent(in) , dimension(:,:) :: solld2d
    call assignpnt(sols2d,solswdir)
    call assignpnt(soll2d,sollwdir)
    call assignpnt(solsd2d,solswdif)
    call assignpnt(solld2d,sollwdif)
  end subroutine init_rad_clm

end module mod_rad_interface
