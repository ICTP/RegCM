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

module mod_che_interface

  use m_realkinds
  use mod_che_common

  public
!
  contains 

  subroutine init_chem(ichem,idirect,dt,chemfrq,dtrad,dsigma,ps1,rh, &
                       icutop,icubot)
    implicit none
    integer , intent(in) :: ichem , idirect
    real(dp) , intent(in) :: dt , chemfrq , dtrad
    real(dp) , pointer , dimension(:) , intent(in) :: dsigma ! dsigma
    real(dp) , pointer , dimension(:,:) , intent(in) :: ps1  ! sps1%ps
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: rh ! rhb3d
    integer , pointer , dimension(:,:) :: icutop , icubot

    if ( ichem == 1 ) lch = .true.
    ichdir = idirect

    chfrq = chemfrq
    rafrq = dtrad
    mdfrq = dt

    chlevs => dsigma
    chps1 => ps1
    chrh => rh
    kcumtop => icutop
    kcumbot => icubot

  end subroutine init_chem
!
end module mod_che_interface
