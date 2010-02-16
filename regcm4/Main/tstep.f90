!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine tstep(extime,dtinc,deltmx)

      use mod_regcm_param
      use mod_param1 , only : dt , dt2 , dtmin
      use mod_date , only : jyear , jyear0 , ktau
      implicit none

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!                                                                     c
!     this subroutine makes a refined start to the model, i.e. the    c
!     model is started with a series of small time steps.             c
!                                                                     c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Dummy arguments
!
      real(8) :: deltmx , dtinc , extime
      intent (in) extime
      intent (out) dtinc
      intent (inout) deltmx
!
! Local variables
!
      real(8) :: deltmn , tscale
      integer :: idtmax
!
!---------------------------------------------------------------------
!
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
