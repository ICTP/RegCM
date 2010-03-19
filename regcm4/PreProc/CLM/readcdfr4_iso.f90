      subroutine readcdfr4_iso(idcdf,vnam,lnam,units,nlon1,nlon,nlat1,  &
                             & nlat,nlev1,nlev,ntim1,ntim,vals)
 
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: idcdf , nlat , nlat1 , nlev , nlev1 , nlon , nlon1 ,   &
               & ntim , ntim1
      character(30) :: lnam
      character(22) :: units
      character(*) :: vnam
      real , dimension(nlon,nlat,nlev,ntim) :: vals
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
