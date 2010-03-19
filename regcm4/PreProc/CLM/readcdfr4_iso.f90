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

      subroutine readcdfr4_iso(idcdf,vnam,lnam,units,nlon1,nlon,nlat1,  &
                             & nlat,nlev1,nlev,ntim1,ntim,vals)
 
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: idcdf , nlat , nlat1 , nlev , nlev1 , nlon , nlon1 ,   &
               & ntim , ntim1
      character(64) :: lnam
      character(64) :: units
      character(64) :: vnam
      real(4) , dimension(nlon,nlat,nlev,ntim) :: vals
      intent (in) nlat , nlat1 , nlev , nlev1 , nlon , nlon1 , ntim ,   &
                & ntim1
      intent (out) vals
!
! Local variables
!
      integer , dimension(4) :: icount , istart
      integer , dimension(2) :: icount1 , istart1
      integer :: iflag , invarid
!
      istart1(1) = nlon1
      icount1(1) = nlon
      istart1(2) = nlat1
      icount1(2) = nlat
 
      istart(1) = nlon1
      icount(1) = nlon
      istart(2) = nlat1
      icount(2) = nlat
      istart(3) = nlev1
      icount(3) = nlev
      istart(4) = ntim1
      icount(4) = ntim
 
      iflag = nf90_inq_varid(idcdf,vnam,invarid)
      if (iflag /= nf90_noerr) go to 920

      iflag = nf90_get_var(idcdf,invarid,vals,istart,icount)
      if (iflag /= nf90_noerr) go to 920

      iflag = nf90_get_att(idcdf,invarid,'long_name',lnam)
      if (iflag /= nf90_noerr) go to 920

      iflag = nf90_get_att(idcdf,invarid,'units',units)
      if (iflag /= nf90_noerr) go to 920
 
      return

 920  write (6, *) 'ERROR: An error occurred while attempting to ',     &
            &      'read variable ', vnam , ' at time ', ntim1
      write (6, *) nf90_strerror(iflag)

      stop 'READ ERROR'

      end subroutine readcdfr4_iso
