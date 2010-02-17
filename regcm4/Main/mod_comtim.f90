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

      module mod_comtim

      implicit none
!
      real(8) :: calday , dtime , twodt
      logical :: doabsems , dolw , dosw
      integer :: mbdate , mbsec , mcdate , mcsec , mdbase , mdcur ,     &
               & msbase , mscur , nelapse , nestep , nnbdat , nnbsec ,  &
               & nndbas , nnsbas , nrstrt , nstep , nstepr , nstop

      end module mod_comtim
