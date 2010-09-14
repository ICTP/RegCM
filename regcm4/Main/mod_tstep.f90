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
      use mod_dynparam
      use mod_runparams
      use mod_date
!
      private
!
      public :: tstep
!
      contains
!
!     This subroutine makes a refined start to the model, i.e. the
!     model is started with a series of small time steps.
!
      subroutine tstep(extime,dtinc,deltmx)

      implicit none
!
      real(8) :: deltmx , dtinc , extime
      intent (in) extime
      intent (out) dtinc
      intent (inout) deltmx
!
      real(8) :: deltmn , tscale
      integer :: idtmax
!
!---------------------------------------------------------------------
!
      deltmn = 0.0
      tscale = 0.0
      if ( jyear.eq.jyear0 .and. ktau.eq.0 ) then
        idtmax = 1
        tscale = 5.*dt
        deltmn = 0.1*dt
        deltmx = dt
      end if
      if ( idtmax.eq.1 ) then
        dt = deltmx*(1.-dexp(-extime/tscale)) + deltmn
        dtinc = dt
        dtmin = dt/60.
        idtmax = 2
      end if
      dt2 = 2.*dt
      if ( jyear.ne.jyear0 .or. ktau.ne.0 ) dt = dt2
!
      end subroutine tstep
!
      end module mod_tstep
