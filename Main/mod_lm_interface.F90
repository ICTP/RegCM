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
  use mod_bats
  use mod_atm_interface , only : slice , surfstate , surfpstate , &
                                 surftstate , domain
!
  contains

  subroutine init_bats(dt,ksrf,ichem,iemiss,idcsst,lakemod,idesseas, &
                       iseaice,dom,atm,sfs,sps,st1,st2,za,ts1,rhox2d)
    implicit none
    real(8) , intent(in) :: dt
    integer(8) , intent(in) :: ksrf
    integer , intent(in) :: ichem , iemiss , idcsst , &
                            lakemod , idesseas , iseaice
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
    xlat     => dom%xlat
    lndcat   => dom%lndcat
    ht       => dom%ht
    uatm     => atm%ubx3d
    vatm     => atm%vbx3d
    tatm     => atm%tb3d
    qvatm    => atm%qvb3d
    thatm    => atm%thx3d
    zpbl     => sfs%zpbl
    hfx      => sfs%hfx
    qfx      => sfs%qfx
    uvdrag   => sfs%uvdrag
    tgbb     => sfs%tgbb
    sfps     => sps%ps
    tground1 => st1%tg
    tground2 => st2%tg
    hgt      => za
    ts       => ts1
    rho      => rhox2d
  end subroutine init_bats
!
#ifdef CLM
  subroutine init_clm(dt,ksrf,ichem,iemiss,idcsst,lakemod,idesseas, &
                      iseaice,dom,dom1,atm,sfs,sps,st1,st2,za,ts1,  &
                      ts0,rhox2d,lm)
    implicit none
    real(8) , intent(in) :: dt
    integer(8) , intent(in) :: ksrf
    integer , intent(in) :: ichem , iemiss , idcsst , &
                            lakemod , idesseas , iseaice
    type(domain) , intent(in) :: dom , dom1
    type(slice) , intent(in) :: atm
    type(surfstate) , intent(in) :: sfs
    type(surfpstate) , intent(in) :: sps
    type(surftstate) , intent(in) :: st1 , st2
    real(8) , pointer , intent(in) , dimension(:,:,:) :: za
    real(8) , pointer , intent(in) , dimension(:,:) :: ts0 , ts1 , rhox2d
    integer , target , intent(in) , dimension(:,:) :: lm

    call init_bats(dt,ksrf,ichem,iemiss,idcsst,lakemod,idesseas, &
                   iseaice,dom,atm,sfs,sps,st1,st2,za,ts1,rhox2d)
    tsf   => ts0
    htf   => dom1%ht
    xlon  => dom%xlon
    lmask => lm
  end subroutine init_clm
#endif
!
end module mod_lm_interface
