C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine tstep(extime,dtinc,deltmx)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c                                                                     c
c     this subroutine makes a refined start to the model, i.e. the    c
c     model is started with a series of small time steps.             c
c                                                                     c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real(kind=8)  extime,dtinc,deltmx
#     include "parameter.inc"
#     include "include/param1.cb"

      integer idtmax
      real(kind=8)  tscale,deltmn
c---------------------------------------------------------------------
c
      if (jyear.eq.jyear0.and.ktau.eq.0) then
       idtmax=1
       tscale=5.*dt
       deltmn=0.1*dt
       deltmx=dt
      end if
      if (idtmax .eq. 1) then
       dt=deltmx*(1.-dexp(-extime/tscale))+deltmn
       dtinc=dt
       dtmin=dt/60.
       idtmax=2
      end if
      dt2=2.*dt
      if (.not.(jyear.eq.jyear0.and.ktau.eq.0)) dt=dt2
c
      return
      end
