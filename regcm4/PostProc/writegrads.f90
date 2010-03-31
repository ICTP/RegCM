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

      subroutine writegrads(idout,outvar,nx,ny,nk,nrec)
 
      implicit none
!
! Dummy arguments
!
      integer :: idout , nk , nrec , nx , ny
      real(4) , dimension(nx,ny,nk) :: outvar
      intent (in) idout , nk , nx , ny , outvar
      intent (inout) nrec
!
! Local variables
!
      integer :: i , j , k
!
      do k = 1 , nk
        nrec = nrec + 1
        write (idout,rec=nrec) ((outvar(i,j,k),i=1,nx),j=1,ny)
      end do
 
      end subroutine writegrads
