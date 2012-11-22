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

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams
  use mod_memutil

  public

  logical :: lchem ! ichem logical equiv
  real(rk8) :: ptp ! ptop
  real(rk8) , pointer , dimension(:) :: flev , hlev ! sigma , a
  real(rk8),  pointer , dimension(:,:) :: twtr
  real(rk8) , pointer , dimension(:,:) :: sfps    ! sfs%psb
  real(rk8) , pointer , dimension(:,:) :: psfps   ! sfs%psa

  real(rk8) , pointer , dimension(:,:,:) :: tatms    ! atms%tb3d
  real(rk8) , pointer , dimension(:,:,:,:) :: qxatms ! atms%qxb3d
  real(rk8) , pointer , dimension(:,:,:) :: rhatms   ! atms%rhb3d
  real(rk8) , pointer , dimension(:,:) :: tground    ! sfs%tgbb
  real(rk8) , pointer , dimension(:,:) :: xlat       ! mddom%xlat

  ! vegetation absorbed radiation (full solar spectrum)
  real(rk8) , pointer , dimension(:,:) :: abveg   ! sabveg
  ! Incident solar flux
  real(rk8) , pointer , dimension(:,:) :: solar   ! solis
  ! Cosine of zenithal solar angle
  real(rk8) , pointer , dimension(:,:) :: coszen  ! coszrs
  ! 0.2-0.7 micro-meter srfc alb to direct radiation
  real(rk8) , pointer , dimension(:,:) :: swdiralb ! aldirs
  ! 0.2-0.7 micro-meter srfc alb to diffuse radiation
  real(rk8) , pointer , dimension(:,:) :: swdifalb ! aldifs
  ! 0.7-5.0 micro-meter srfc alb to direct radiation
  real(rk8) , pointer , dimension(:,:) :: lwdiralb ! aldirl
  ! 0.7-5.0 micro-meter srfc alb to diffuse radiation
  real(rk8) , pointer , dimension(:,:) :: lwdifalb ! aldifl
  ! Total Short wave albedo (0.2-0.7 micro-meter)
  real(rk8) , pointer , dimension(:,:) :: swalb   ! albvs
  ! Total Long wave albedo (0.7-5.0 micro-meter)
  real(rk8) , pointer , dimension(:,:) :: lwalb   ! albvl
  ! Emissivity at surface
  real(rk8) , pointer , dimension(:,:) :: emsvt   ! emiss1d
  ! Bidimensional collector storage for above
  real(rk8) , pointer , dimension(:,:) :: totsol ! sinc
  real(rk8) , pointer , dimension(:,:) :: soldir ! solvs
  real(rk8) , pointer , dimension(:,:) :: soldif ! solvd
  real(rk8) , pointer , dimension(:,:) :: solswdir ! sols
  real(rk8) , pointer , dimension(:,:) :: sollwdir ! soll
  real(rk8) , pointer , dimension(:,:) :: solswdif ! solsd
  real(rk8) , pointer , dimension(:,:) :: sollwdif ! solld
  real(rk8) , pointer , dimension(:,:) :: srfabswflx ! fsw
  real(rk8) , pointer , dimension(:,:) :: srflwflxup ! flw
  real(rk8) , pointer , dimension(:,:) :: srflwflxdw ! flwd

  ! Land Ocean Ice (1,0,2) mask
  integer(ik4) , pointer , dimension(:,:,:) :: lndocnicemsk ! ldmsk12d

  real(rk8) , pointer , dimension(:,:,:,:) :: chspmix  ! chia

  character(len=6) , pointer , dimension(:) :: tracname ! chtrname

  real(rk8) , pointer , dimension(:,:,:) :: cldfra , cldlwc
  real(rk8) , pointer , dimension(:,:,:) :: heatrt
  real(rk8) , pointer , dimension(:,:,:) :: o3prof
  real(rk8) , pointer , dimension(:,:,:) :: aerasp , aerext , aerssa
  real(rk8) , pointer , dimension(:,:) :: aersrrf , aertarf
  real(rk8) , pointer , dimension(:,:) :: aertalwrf , aersrlwrf
! absnxt  - Nearest layer absorptivities
! abstot  - Non-adjacent layer absorptivites
! emstot  - Total emissivity
  real(rk8) , pointer , dimension(:,:,:,:)  :: gasabsnxt
  real(rk8) , pointer , dimension(:,:,:,:)  :: gasabstot
  real(rk8) , pointer , dimension(:,:,:) :: gasemstot

  real(rk8) , pointer , dimension(:,:,:,:) :: taucldsp

  real(rk8) , pointer , dimension(:,:) :: ptrop
  integer(ik4) , pointer , dimension(:,:) :: ktrop

  logical :: doabsems , dolw , dosw
  integer(ik4) :: ichso4 , ichbc , ichoc

  real(rk8) :: chfrovrradfr ! chfrq/rafrq

  integer(ik8) :: ntabem

  data lchem /.false./

  contains 

  subroutine allocate_mod_rad_common(ichem)
    implicit none
    integer(ik4) , intent(in) :: ichem
    call getmem3d(cldfra,jci1,jci2,ici1,ici2,1,kz,'rad:cldfra')
    call getmem3d(cldlwc,jci1,jci2,ici1,ici2,1,kz,'rad:cldlwc')
    call getmem3d(heatrt,jci1,jci2,ici1,ici2,1,kz,'rad:heatrt')
    call getmem3d(o3prof,jci1,jci2,ici1,ici2,1,kzp1,'rad:o3prof')
    call getmem2d(ptrop,jci1,jci2,ici1,ici2,'rad:ptrop')
    call getmem2d(ktrop,jci1,jci2,ici1,ici2,'rad:ktrop')
    if ( irrtm == 0 ) then
      call getmem4d(gasabsnxt,jci1,jci2,ici1,ici2,1,kz,1,4,'rad:gasabsnxt')
      call getmem4d(gasabstot,jci1,jci2,ici1,ici2,1,kzp1,1,kzp1,'rad:gasabstot')
      call getmem3d(gasemstot,jci1,jci2,ici1,ici2,1,kzp1,'rad:gasemstot')
    end if

    if ( ichem == 1 ) then
      call getmem4d(taucldsp,jci1,jci2,ici1,ici2,0,kz,1,nspi,'rad:taucldsp')
      call getmem3d(aerasp,jci1,jci2,ici1,ici2,1,kz,'rad:aerasp')
      call getmem3d(aerext,jci1,jci2,ici1,ici2,1,kz,'rad:aerext')
      call getmem3d(aerssa,jci1,jci2,ici1,ici2,1,kz,'rad:aerssa')
      call getmem2d(aersrrf,jci1,jci2,ici1,ici2,'rad:aersrrf')
      call getmem2d(aertalwrf,jci1,jci2,ici1,ici2,'rad:aertalwrf')
      call getmem2d(aersrlwrf,jci1,jci2,ici1,ici2,'rad:aersrlwrf')
      call getmem2d(aertarf,jci1,jci2,ici1,ici2,'rad:aertarf')
    end if
  end subroutine  allocate_mod_rad_common

end module mod_rad_common
