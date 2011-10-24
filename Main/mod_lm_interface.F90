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
  use mod_atm_interface , only : slice , surfstate , surfpstate , &
                                 surftstate , domain
#ifdef CLM
  use mod_mtrxclm
  use mod_clm
  use mod_bats_mppio
  use mod_bats_vecbats
#else
  use mod_bats_param
  use mod_bats_mppio
  use mod_bats_bndry
  use mod_bats_co2
  use mod_bats_drag
  use mod_bats_lake
  use mod_bats_leaftemp
  use mod_bats_vecbats
  use mod_bats_zengocn
  use mod_bats_romsocn
#endif
!
  public

  integer :: idcsst , lakemod , idesseas , iseaice

  contains

  subroutine init_bats(dt,ksrf,ichem,iemiss,dom,atm,sfs, &
                       sps,st1,st2,za,ts1,rhox2d)
    implicit none
    real(8) , intent(in) :: dt
    integer(8) , intent(in) :: ksrf
    integer , intent(in) :: ichem , iemiss
    type(domain) , intent(in) :: dom
    type(slice) , intent(in) :: atm
    type(surfstate) , intent(in) :: sfs
    type(surfpstate) , intent(in) :: sps
    type(surftstate) , intent(in) :: st1 , st2
    real(8) , pointer , intent(in) , dimension(:,:,:) :: za
    real(8) , pointer , intent(in) , dimension(:,:) :: ts1 , rhox2d
    xdtsec = dt
    kbats = ksrf
    if ( ichem    == 1 ) lchem     = .true.
    if ( iemiss   == 1 ) lemiss    = .true.
    if ( idcsst   == 1 ) ldcsst    = .true.
    if ( lakemod  == 1 ) llake     = .true.
    if ( idesseas == 1 ) ldesseas  = .true.
    if ( iseaice  == 1 ) lseaice   = .true.
    call assignpnt(dom%xlat,xlat)
    call assignpnt(dom%lndcat,lndcat)
    call assignpnt(dom%ht,ht)
    call assignpnt(atm%ubx3d,uatm)
    call assignpnt(atm%vbx3d,vatm)
    call assignpnt(atm%tb3d,tatm)
    call assignpnt(atm%qvb3d,qvatm)
    call assignpnt(atm%thx3d,thatm)
    call assignpnt(sfs%zpbl,zpbl)
    call assignpnt(sfs%hfx,hfx)
    call assignpnt(sfs%qfx,qfx)
    call assignpnt(sfs%uvdrag,uvdrag)
    call assignpnt(sfs%tgbb,tgbb)
    call assignpnt(sps%ps,sfps)
    call assignpnt(st1%tg,tground1)
    call assignpnt(st2%tg,tground2)
    call assignpnt(za,hgt)
    call assignpnt(ts1,ts)
    call assignpnt(rhox2d,rho)
  end subroutine init_bats
!
#ifdef CLM
  subroutine init_clm(dt,ksrf,ichem,iemiss,dom,dom1,atm,sfs,&
                      sps,st1,st2,za,ts1,ts0,rhox2d,lm)
    implicit none
    real(8) , intent(in) :: dt
    integer(8) , intent(in) :: ksrf
    integer , intent(in) :: ichem , iemiss
    type(domain) , intent(in) :: dom , dom1
    type(slice) , intent(in) :: atm
    type(surfstate) , intent(in) :: sfs
    type(surfpstate) , intent(in) :: sps
    type(surftstate) , intent(in) :: st1 , st2
    real(8) , pointer , intent(in) , dimension(:,:,:) :: za
    real(8) , pointer , intent(in) , dimension(:,:) :: ts0 , ts1 , rhox2d
    integer , target , intent(in) , dimension(:,:) :: lm

    call init_bats(dt,ksrf,ichem,iemiss,dom,atm,sfs,sps,st1,st2,za,ts1,rhox2d)
    call assignpnt(ts0,tsf)
    call assignpnt(dom1%ht,htf)
    call assignpnt(dom%xlon,xlon)
    call assignpnt(lm,lmask)
  end subroutine init_clm
#endif
!
end module mod_lm_interface
