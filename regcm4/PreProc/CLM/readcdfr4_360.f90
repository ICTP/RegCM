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

      subroutine readcdfr4_360(idcdf,vnam,lnam,units,nlon1,nlon,nlat1,  &
                             & nlat,nlev1,nlev,ntim1,ntim,nglon,nglat,  &
                             & nglev,ngtim,vals)
 
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: idcdf , nglat , nglev , nglon , ngtim , nlat , nlat1 , &
               & nlev , nlev1 , nlon , nlon1 , ntim , ntim1
      character(64) :: lnam , units , vnam
      real(4) , dimension(nlon,nlat,nlev,ntim) :: vals
      intent (in) nglat , nglev , nglon , ngtim , nlat , nlat1 , nlev , &
                & nlev1 , nlon , nlon1 , ntim , ntim1
      intent (out) vals
!
! Local variables
!
      integer :: i , iflag , ii , ilon5 , invarid , j , jj , k , kk ,   &
               & l , ll , nlat2 , nlev2 , nlon2 , ntim2
      integer , dimension(4) :: icount , istart
      real(4) , allocatable , dimension(:,:,:,:) :: vals1 , vals2
! 
      istart(1) = 1
      icount(1) = nglon
      istart(2) = 1
      icount(2) = nglat
      istart(3) = 1
      icount(3) = nglev
      istart(4) = 1
      icount(4) = ngtim
!     /*get variable and attributes*/
      ilon5 = nglon/2
      allocate(vals1(nglon,nglat,nglev,ngtim))
      allocate(vals2(nglon,nglat,nglev,ngtim))

      iflag = nf90_inq_varid(idcdf,vnam,invarid)
      if (iflag /= nf90_noerr) go to 920

      iflag = nf90_get_var(idcdf,invarid,vals1,istart,icount)
      if (iflag /= nf90_noerr) go to 920

      do l = 1 , ngtim
        do k = 1 , nglev
          do j = 1 , nglat
            do i = 1 , nglon
              vals2(i,j,k,l) = vals1(i,j,k,l)
            end do
          end do
        end do
      end do
      do l = 1 , ngtim
        do k = 1 , nglev
          do j = 1 , nglat
            do i = 1 , nglon
              ii = ilon5 + i
              if ( ii>nglon ) ii = ii - nglon
              vals1(ii,j,k,l) = vals2(i,j,k,l)
            end do
          end do
        end do
      end do
 
      ntim2 = ntim1 + ntim - 1
      nlev2 = nlev1 + nlev - 1
      nlat2 = nlat1 + nlat - 1
      nlon2 = nlon1 + nlon - 1
      do l = ntim1 , ntim2
        do k = nlev1 , nlev2
          do j = nlat1 , nlat2
            do i = nlon1 , nlon2
              ll = l - ntim1 + 1
              kk = k - nlev1 + 1
              jj = j - nlat1 + 1
              ii = i - nlon1 + 1
              vals(ii,jj,kk,ll) = vals1(i,j,k,l)
            end do
          end do
        end do
      end do
 
      deallocate(vals1)
      deallocate(vals2)
 
      iflag = nf90_get_att(idcdf,invarid,'long_name',lnam)
      if (iflag /= nf90_noerr) go to 920

      iflag = nf90_get_att(idcdf,invarid,'units',units)
      if (iflag /= nf90_noerr) go to 920
 
      return

 920  write (6, *) 'ERROR: An error occurred while attempting to ',     &
            &      'read variable ', vnam , ' at time ', ntim1
      write (6, *) nf90_strerror(iflag)

      stop 'READ ERROR'

      end subroutine readcdfr4_360
