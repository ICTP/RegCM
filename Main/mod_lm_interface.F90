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

module mod_lm_interface
!
! Link surface and atmospheric models
!
  use mod_bats_common
  use mod_runparams
  use mod_memutil
  use mod_regcm_types
#ifdef CLM
  use mod_mtrxclm
  use mod_clm
  use mod_bats_mtrxbats
#else
  use mod_bats_param
  use mod_bats_bndry
  use mod_bats_co2
  use mod_bats_drag
  use mod_bats_lake
  use mod_bats_leaftemp
  use mod_bats_mtrxbats
  use mod_bats_zengocn
#endif
!
  public

  contains

  subroutine init_bats(dom,sub,atm,sfs,zpbl,pptc,pptnc,coszrs, &
                  fsw,flw,flwd,sabveg,albvs,albvl,aldirs,aldifs,aldirl, &
                  aldifl,solis,emiss,sinc,solvs,solvsd,solvl,solvld)
    implicit none
    type(domain) , intent(in) :: dom
    type(domain_subgrid) , intent(in) :: sub
    type(slice) , intent(in) :: atm
    type(surfstate) , intent(in) :: sfs
    real(rk8) , pointer , intent(in) , dimension(:,:) :: zpbl
    real(rk8) , pointer , intent(in) , dimension(:,:) :: pptc
    real(rk8) , pointer , intent(in) , dimension(:,:) :: pptnc
    real(rk8) , pointer , intent(in) , dimension(:,:) :: coszrs
    real(rk8) , pointer , intent(in) , dimension(:,:) :: flw
    real(rk8) , pointer , intent(in) , dimension(:,:) :: fsw
    real(rk8) , pointer , intent(in) , dimension(:,:) :: flwd
    real(rk8) , pointer , intent(in) , dimension(:,:) :: sabveg
    real(rk8) , pointer , intent(in) , dimension(:,:) :: albvs
    real(rk8) , pointer , intent(in) , dimension(:,:) :: albvl
    real(rk8) , pointer , intent(in) , dimension(:,:) :: aldirs
    real(rk8) , pointer , intent(in) , dimension(:,:) :: aldifs
    real(rk8) , pointer , intent(in) , dimension(:,:) :: aldirl
    real(rk8) , pointer , intent(in) , dimension(:,:) :: aldifl
    real(rk8) , pointer , intent(in) , dimension(:,:) :: solis
    real(rk8) , pointer , intent(in) , dimension(:,:) :: emiss
    real(rk8) , pointer , intent(in) , dimension(:,:) :: sinc
    real(rk8) , pointer , intent(in) , dimension(:,:) :: solvs
    real(rk8) , pointer , intent(in) , dimension(:,:) :: solvsd
    real(rk8) , pointer , intent(in) , dimension(:,:) :: solvl
    real(rk8) , pointer , intent(in) , dimension(:,:) :: solvld
    ntcpl  = idnint(cpldt/dtsec)
    ntsrf2 = idnint(dtsrf/dtsec)
    if ( idcsst   == 1 ) ldcsst    = .true.
    if ( lakemod  == 1 ) llake     = .true.
    if ( idesseas == 1 ) ldesseas  = .true.
    if ( iseaice  == 1 ) lseaice   = .true.
    call assignpnt(dom%xlat,xlat)
    call assignpnt(dom%xlon,xlon)
    call assignpnt(dom%lndcat,lndcat)
    call assignpnt(dom%ldmsk,landmsk)
    call assignpnt(dom%ht,ht)
    call assignpnt(sub%xlat,xlat1)
    call assignpnt(sub%xlon,xlon1)
    call assignpnt(sub%lndcat,lndcat1)
    call assignpnt(sub%ldmsk,ldmsk1)
    call assignpnt(sub%ht,ht1)
    call assignpnt(sub%dhlake,dhlake1)
    call assignpnt(dom%snowam,snowam)
    call assignpnt(atm%ubx3d,uatm,kz)
    call assignpnt(atm%vbx3d,vatm,kz)
    call assignpnt(atm%tb3d,tatm,kz)
    call assignpnt(atm%qxb3d,qvatm,kz,iqv)
    call assignpnt(atm%thx3d,thatm,kz)
    call assignpnt(atm%rhox2d,rhox)
    call assignpnt(atm%za,hgt,kz)
    call assignpnt(sfs%hfx,hfx)
    call assignpnt(sfs%qfx,qfx)
    call assignpnt(sfs%uvdrag,uvdrag)
    call assignpnt(sfs%tgbb,tgbb)
    call assignpnt(sfs%psb,sfps)
    call assignpnt(sfs%tga,tground1)
    call assignpnt(sfs%tgb,tground2)
    call assignpnt(zpbl,hpbl)
    call assignpnt(pptc,cprate)
    call assignpnt(pptnc,ncprate)
    call assignpnt(coszrs,zencos)
    call assignpnt(fsw,rswf)
    call assignpnt(flw,rlwf)
    call assignpnt(flwd,dwrlwf)
    call assignpnt(sabveg,vegswab)
    call assignpnt(albvl,lwalb)
    call assignpnt(albvs,swalb)
    call assignpnt(aldirs,swdiralb)
    call assignpnt(aldifs,swdifalb)
    call assignpnt(aldirl,lwdiralb)
    call assignpnt(aldifl,lwdifalb)
    call assignpnt(solis,solar)
    call assignpnt(emiss,emissivity)
    call assignpnt(sinc,solinc)
    call assignpnt(solvs,swdir)
    call assignpnt(solvsd,swdif)
    call assignpnt(solvl,lwdir)
    call assignpnt(solvld,lwdif)
  end subroutine init_bats
!
#ifdef CLM
  subroutine init_clm(lm)
    implicit none
     integer , pointer , intent(in) , dimension(:,:) :: lm
    call assignpnt(lm,lmask)
  end subroutine init_clm
#endif
!
end module mod_lm_interface
