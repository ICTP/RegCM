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

module mod_cu_interface
!
! Link atmospheric model and cumulus schemes
!
  use mod_constants
  use mod_cu_common
  use mod_atm_interface , only : atmstate , slice , surfstate , &
                                 surfpstate , domain
!
  contains

  subroutine init_cuscheme(ichem,dtsec,ntsrf,mddom,atm1,atm2,aten,atms, &
                           sfsta,sps1,sps2,za,qdot,pptc,ldmsk,sigma,a,  &
                           dsigma,qcon,cldfra,cldlwc)
    implicit none
    real(8) , intent(in) :: dtsec
    integer(8) , intent(in) :: ntsrf
    integer , intent(in) :: ichem
    type(domain) , intent(in) :: mddom
    type(atmstate) , intent(in) :: atm1 , atm2 , aten
    type(slice) , intent(in) :: atms
    type(surfstate) , intent(in) :: sfsta
    type(surfpstate) , intent(in) :: sps1 , sps2
    real(8) , pointer , intent(in) , dimension(:,:,:) :: za , qdot
    real(8) , pointer , intent(in) , dimension(:,:) :: pptc
    integer , pointer , intent(in) , dimension(:,:) :: ldmsk
    real(8) , pointer , intent(in) , dimension(:) :: sigma , a , dsigma , qcon
    real(8) , pointer , dimension(:,:) :: cldlwc , cldfra
!
    if ( ichem    == 1 ) lchem = .true.
    dtcum  = dtsec
    dtmdl  = dtsec
    dtcum2 = dtsec*d_two
    aprdiv = d_one/dble(ntsrf)

    sfhgt   => mddom%ht
    psfcps  => sps1%ps
    sfcps   => sps2%ps
    hgt     => za
    ptatm   => atm1%t
    puatm   => atm1%u
    pvatm   => atm1%v
    pvqvtm  => atm1%qv
    tatm    => atm2%t
    uatm    => atm2%u
    vatm    => atm2%v
    qvatm   => atm2%qv
    qcatm   => atm2%qc
    tas     => atms%tb3d
    uas     => atms%ubx3d
    vas     => atms%vbx3d
    pas     => atms%pb3d
    qsas    => atms%qsb3d
    qvas    => atms%qvb3d
    tten    => aten%t
    uten    => aten%u
    vten    => aten%v
    qvten   => aten%qv
    qcten   => aten%qc
    rainc   => sfsta%rainc
    qfx     => sfsta%qfx
    svv     => qdot
    lmpcpc  => pptc
    lmask   => ldmsk
    flev    => sigma
    hlev    => a
    dflev   => dsigma
    wlev    => qcon
    rcldfra => cldfra
    rcldlwc => cldlwc

  end subroutine init_cuscheme
!
end module mod_cu_interface
