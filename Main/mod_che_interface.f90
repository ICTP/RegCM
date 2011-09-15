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
!
  use m_realkinds
  use mod_atm_interface , only : surfpstate
  use mod_che_common
  use mod_che_cumtran
  use mod_che_dust
  use mod_che_indices
  use mod_che_mppio
  use mod_che_ncio
  use mod_che_param
  use mod_che_drydep
  use mod_che_emission
  use mod_che_species
!
  public
!
  contains 
!
  subroutine init_chem(idirect,dt,chemfrq,dtrad,dsigma,sps2,icutop,icubot)
    implicit none
    integer , intent(in) :: idirect
    real(dp) , intent(in) :: dt , chemfrq , dtrad
    real(dp) , pointer , dimension(:) , intent(in) :: dsigma ! dsigma
    integer , pointer , dimension(:,:) :: icutop , icubot
    type(surfpstate) , intent(in) :: sps2

    ichdir = idirect

    chfrq = chemfrq
    rafrq = dtrad
    dtche = dt

    call assignpnt(dsigma,chlevs)
    call assignpnt(icutop,kcumtop)
    call assignpnt(icubot,kcumbot)
    call assignpnt(sps2%ps,sfcp)

  end subroutine init_chem
!
end module mod_che_interface
