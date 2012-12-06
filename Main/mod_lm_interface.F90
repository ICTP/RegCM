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
  use mod_atm_interface , only : slice , surfstate , domain
#ifdef CLM
  use mod_mtrxclm
  use mod_clm
  use mod_bats_mtrxbats
#else
  use mod_bats_param
  use mod_bats_mppio
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

#ifdef CLM
  real(rk8) , dimension(:,:,:) , pointer :: dhlake1_io
#endif

  contains

  subroutine init_bats(dt,ksrf,ichem,iemiss,dom,atm,sfs,zpbl)
    implicit none
    real(rk8) , intent(in) :: dt
    integer(ik8) , intent(in) :: ksrf
    integer(ik4) , intent(in) :: ichem , iemiss
    type(domain) , intent(in) :: dom
    type(slice) , intent(in) :: atm
    type(surfstate) , intent(in) :: sfs
    real(rk8) , pointer , intent(in) , dimension(:,:) :: zpbl
    kbats = ksrf
    ntcpl  = idnint(cpldt/dtsec)
    ntsrf2 = idnint(dtsrf/dtsec)
    if ( ichem    == 1 ) lchem     = .true.
    if ( iemiss   == 1 ) lemiss    = .true.
    if ( idcsst   == 1 ) ldcsst    = .true.
    if ( lakemod  == 1 ) llake     = .true.
    if ( idesseas == 1 ) ldesseas  = .true.
    if ( iseaice  == 1 ) lseaice   = .true.
    call assignpnt(dom%xlat,xlat)
    call assignpnt(dom%xlon,xlon)
    call assignpnt(dom%lndcat,lndcat)
    call assignpnt(dom%ht,ht)
    call assignpnt(dom%snowam,snowam)
    call assignpnt(atm%ubx3d,uatm)
    call assignpnt(atm%vbx3d,vatm)
    call assignpnt(atm%tb3d,tatm)
    call assignpnt(atm%qxb3d,qxatm)
    call assignpnt(atm%thx3d,thatm)
    call assignpnt(atm%rhox2d,rhox)
    call assignpnt(atm%za,hgt)
    call assignpnt(sfs%hfx,hfx)
    call assignpnt(sfs%qfx,qfx)
    call assignpnt(sfs%uvdrag,uvdrag)
    call assignpnt(sfs%tgbb,tgbb)
    call assignpnt(sfs%psb,sfps)
    call assignpnt(sfs%tga,tground1)
    call assignpnt(sfs%tgb,tground2)
    call assignpnt(zpbl,hpbl)
  end subroutine init_bats
!
#ifdef CLM
  subroutine init_clm(dt,ksrf,ichem,iemiss,dom,atm,sfs,zpbl,lm)
    implicit none
    real(rk8)       , intent(in) :: dt
    integer(ik8)    , intent(in) :: ksrf
    integer(ik4)    , intent(in) :: ichem , iemiss
    type(domain)    , intent(in) :: dom
    type(slice)     , intent(in) :: atm
    type(surfstate) , intent(in) :: sfs
    real(rk8)       , pointer , intent(in) , dimension(:,:) :: zpbl
    integer         , pointer , intent(in) , dimension(:,:) :: lm

    call init_bats(dt,ksrf,ichem,iemiss,dom,atm,sfs,zpbl)
    call assignpnt(lm,lmask)
  end subroutine init_clm
#endif
!
end module mod_lm_interface
