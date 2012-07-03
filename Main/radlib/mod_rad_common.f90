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

module mod_rad_common

  use mod_constants
  use mod_dynparam
  use mod_realkinds
  use mod_memutil

  public

  integer :: irrtm , irrtm_cldov , irrtm_sw_opcliq , irrtm_sw_opcice

  logical :: lchem ! ichem logical equiv
  real(dp) :: ptp ! ptop
  real(dp) , pointer , dimension(:) :: flev , hlev ! sigma , a
  real(dp),  pointer , dimension(:,:) :: twtr
  real(dp) , pointer , dimension(:,:) :: sfps    ! sfs%psb
  real(dp) , pointer , dimension(:,:) :: psfps   ! sfs%psa

  real(dp) , pointer , dimension(:,:,:) :: tatms    ! atms%tb3d
  real(dp) , pointer , dimension(:,:,:,:) :: qxatms ! atms%qxb3d
  real(dp) , pointer , dimension(:,:,:) :: rhatms   ! atms%rhb3d
  real(dp) , pointer , dimension(:,:) :: tground    ! sfs%tgbb
  real(dp) , pointer , dimension(:,:) :: xlat       ! mddom%xlat

  ! vegetation absorbed radiation (full solar spectrum)
  real(dp) , pointer , dimension(:,:) :: abveg   ! sabveg
  ! Incident solar flux
  real(dp) , pointer , dimension(:,:) :: solar   ! solis
  ! Cosine of zenithal solar angle
  real(dp) , pointer , dimension(:,:) :: coszen  ! coszrs
  ! 0.2-0.7 micro-meter srfc alb to direct radiation
  real(dp) , pointer , dimension(:,:) :: swdiralb ! aldirs
  ! 0.2-0.7 micro-meter srfc alb to diffuse radiation
  real(dp) , pointer , dimension(:,:) :: swdifalb ! aldifs
  ! 0.7-5.0 micro-meter srfc alb to direct radiation
  real(dp) , pointer , dimension(:,:) :: lwdiralb ! aldirl
  ! 0.7-5.0 micro-meter srfc alb to diffuse radiation
  real(dp) , pointer , dimension(:,:) :: lwdifalb ! aldifl
  ! Total Short wave albedo (0.2-0.7 micro-meter)
  real(dp) , pointer , dimension(:,:) :: swalb   ! albvs
  ! Total Long wave albedo (0.7-5.0 micro-meter)
  real(dp) , pointer , dimension(:,:) :: lwalb   ! albvl
  ! Emissivity at surface
  real(dp) , pointer , dimension(:,:) :: emsvt   ! emiss1d
  ! Bidimensional collector storage for above
  real(dp) , pointer , dimension(:,:) :: totsol ! sinc
  real(dp) , pointer , dimension(:,:) :: soldir ! solvs
  real(dp) , pointer , dimension(:,:) :: soldif ! solvd
  real(dp) , pointer , dimension(:,:) :: solswdir ! sols
  real(dp) , pointer , dimension(:,:) :: sollwdir ! soll
  real(dp) , pointer , dimension(:,:) :: solswdif ! solsd
  real(dp) , pointer , dimension(:,:) :: sollwdif ! solld
  real(dp) , pointer , dimension(:,:) :: srfabswflx ! fsw
  real(dp) , pointer , dimension(:,:) :: srflwflxup ! flw
  real(dp) , pointer , dimension(:,:) :: srflwflxdw ! flwd

  ! Land Ocean Ice (1,0,2) mask
  integer , pointer , dimension(:,:,:) :: lndocnicemsk ! ldmsk12d

  real(dp) , pointer , dimension(:,:,:,:) :: chspmix  ! chia

  character(len=6) , pointer , dimension(:) :: tracname ! chtrname

  integer :: ncld ! # of bottom model levels with no clouds

  real(dp) , pointer , dimension(:,:,:) :: cldfra , cldlwc
  real(dp) , pointer , dimension(:,:,:) :: heatrt
  real(dp) , pointer , dimension(:,:,:) :: o3prof
  real(dp) , pointer , dimension(:,:,:) :: aerasp , aerext , aerssa
  real(dp) , pointer , dimension(:,:) :: aersrrf , aertarf
  real(dp) , pointer , dimension(:,:) :: aertalwrf , aersrlwrf
! absnxt  - Nearest layer absorptivities
! abstot  - Non-adjacent layer absorptivites
! emstot  - Total emissivity
  real(dp) , pointer , dimension(:,:,:,:)  :: gasabsnxt
  real(dp) , pointer , dimension(:,:,:,:)  :: gasabstot
  real(dp) , pointer , dimension(:,:,:) :: gasemstot

  real(dp) , pointer , dimension(:,:,:,:) :: taucldsp

  real(dp) , pointer , dimension(:,:) :: ptrop

  integer :: iclimao3
  logical :: doabsems , dolw , dosw
  integer :: ichso4 , ichbc , ichoc

  real(dp) :: chfrovrradfr ! chfrq/rafrq

  integer(8) :: ntabem

  data lchem /.false./

  contains 

  subroutine allocate_mod_rad_common(ichem)
    implicit none
    integer , intent(in) :: ichem
    call getmem3d(cldfra,jci1,jci2,ici1,ici2,1,kz,'rad:cldfra')
    call getmem3d(cldlwc,jci1,jci2,ici1,ici2,1,kz,'rad:cldlwc')
    call getmem3d(heatrt,jce1,jce2,ice1,ice2,1,kz,'rad:heatrt')
    call getmem3d(o3prof,jce1,jce2,ice1,ice2,1,kzp1,'rad:o3prof')
    call getmem2d(ptrop,jci1,jci2,ici1,ici2,'rad:ptrop')
    call getmem4d(taucldsp,jci1,jci2,ici1,ici2,0,kz,1,nspi,'rad:taucldsp')
    if ( irrtm == 0 ) then
      call getmem4d(gasabsnxt,jce1,jce2,ice1,ice2,1,kz,1,4,'rad:gasabsnxt')
      call getmem4d(gasabstot,jce1,jce2,ice1,ice2,1,kzp1,1,kzp1,'rad:gasabstot')
      call getmem3d(gasemstot,jce1,jce2,ice1,ice2,1,kzp1,'rad:gasemstot')
    end if

    if ( ichem == 1 ) then
      call getmem3d(aerasp,jce1,jce2,ice1,ice2,1,kz,'rad:aerasp')
      call getmem3d(aerext,jce1,jce2,ice1,ice2,1,kz,'rad:aerext')
      call getmem3d(aerssa,jce1,jce2,ice1,ice2,1,kz,'rad:aerssa')
      call getmem2d(aersrrf,jce1,jce2,ice1,ice2,'rad:aersrrf')
      call getmem2d(aertalwrf,jce1,jce2,ice1,ice2,'rad:aertalwrf')
      call getmem2d(aersrlwrf,jce1,jce2,ice1,ice2,'rad:aersrlwrf')
      call getmem2d(aertarf,jce1,jce2,ice1,ice2,'rad:aertarf')
    end if
  end subroutine  allocate_mod_rad_common

end module mod_rad_common
