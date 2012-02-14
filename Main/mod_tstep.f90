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
 
module mod_tstep
!
! Calculate smaller timesteps for model startup
!
  use mod_runparams
!
  private
!
  public :: tstep
!
  contains
!
! This subroutine makes a refined start to the model, i.e. the
! model is started with a series of small time steps.
!
  subroutine tstep(extime,dtinc)
!
    implicit none
!
    real(dp) , intent (in) :: extime
    real(dp) , intent (out) :: dtinc
!
    real(dp) :: deltmn , tscale
    integer :: idtmax
!
!---------------------------------------------------------------------
!
    deltmn = d_zero
    tscale = d_zero
    if ( ktau == 0 ) then
      idtmax = 1
      tscale = d_five*dt
      deltmn = 0.1D0*dt
      deltmx = dt
    end if
    if ( idtmax == 1 ) then
      dt = deltmx*(d_one-dexp(-extime/tscale)) + deltmn
      dtinc = dt
      idtmax = 2
    end if
    dt2 = d_two*dt
    if ( ktau /= 0 ) then
      dt = dt2
    end if
!
  end subroutine tstep
!
end module mod_tstep
