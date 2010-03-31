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
 
      subroutine reorder(fdp,fsp,jx,iy,nz)
 
      implicit none
!
! Dummy arguments
!
      integer :: iy , jx , nz
      real(4) , dimension(nz*nz,jx,iy) :: fdp
      real(4) , dimension(jx*nz,iy*nz) :: fsp
      intent (in) fdp , iy , jx , nz
      intent (out) fsp
!
! Local variables
!
      integer :: i , ii , j , jj , k
!
      do j = 1 , jx*nz
        do i = 1 , iy*nz
          jj = mod(j,nz)
          if ( jj.eq.0 ) jj = nz
          ii = mod(i,nz)
          if ( ii.eq.0 ) ii = nz
          k = (jj-1)*nz + ii
          jj = (j+nz-1)/nz
          ii = (i+nz-1)/nz
          fsp(j,i) = fdp(k,jj,ii)
        end do
      end do
 
      end subroutine reorder
