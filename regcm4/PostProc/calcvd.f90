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

      subroutine calcvd(fld3d,nx,ny,kz,nfld3d,ds,dmap,xmap,nua,nva,nvor,&
                      & ndiv,nx1,ny1)
 
      implicit none
!
! Dummy arguments
!
      real(4) :: ds
      integer :: ndiv , nfld3d , nua , nva , nvor , nx , nx1 , ny ,     &
               & ny1 , kz
      real(4) , dimension(nx,ny) :: dmap , xmap
      real(4) , dimension(nx,ny,kz,nfld3d) :: fld3d
      intent (in) dmap , ds , ndiv , nfld3d , nua , nva , nvor , nx ,   &
                & nx1 , ny , ny1 , kz , xmap
      intent (inout) fld3d
!
! Local variables
!
      real(4) :: ds2r , u1 , u2 , u3 , u4 , v1 , v2 , v3 , v4
      integer :: i , j , k
      real(4) , dimension(nx,ny,kz) :: u , v
!
      ds2r = 1.0/(2.0*ds)
 
      u(:,:,:) = fld3d(:,:,:,nua)
      v(:,:,:) = fld3d(:,:,:,nva)
 
      do k = 1 , kz
        do j = 1 , ny1
          do i = 1 , nx1
            u1 = u(i,j,k)/dmap(i,j)
            u2 = u(i+1,j,k)/dmap(i+1,j)
            u3 = u(i,j+1,k)/dmap(i,j+1)
            u4 = u(i+1,j+1,k)/dmap(i+1,j+1)
            v1 = v(i,j,k)/dmap(i,j)
            v2 = v(i+1,j,k)/dmap(i+1,j)
            v3 = v(i,j+1,k)/dmap(i,j+1)
            v4 = v(i+1,j+1,k)/dmap(i+1,j+1)
            fld3d(i,j,k,nvor) = xmap(i,j)*xmap(i,j)                     &
                              & *ds2r*((v4-v2+v3-v1)-(u2-u1+u4-u3))
            fld3d(i,j,k,ndiv) = xmap(i,j)*xmap(i,j)                     &
                              & *ds2r*((u3-u1+u4-u2)+(v2-v1+v4-v3))
 
          end do
        end do
      end do
      end subroutine calcvd
