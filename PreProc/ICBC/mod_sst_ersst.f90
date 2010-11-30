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

      module mod_sst_ersst

      use m_realkinds
      use m_die
      use m_stdio

      contains

      subroutine sst_ersst

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Comments on dataset sources and location:                          !
!                                                                    !
! ERAIN    ERAIN_SST is provided by ERA-Interim project              !
!          6 hourly frequncy, 1.5x1.5 degree resolution              !
!          from 1989010100 to 2009053118.                            !
!          'ERSST' for using the sea surface temperature;            !
!          'ERSKT' for using the skin temperature.                   !
!                                                                    !
!          ML= 1 is   0.0; ML= 2 is   1.5; => ML=240 is 358.5E       !
!          NL= 1 is  90.0; ML= 2 is  88.5; => ML=121 is -90.         !
!                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use mod_sst_grid
      use mod_date
      use mod_interp , only : bilinx
      use netcdf

      implicit none
!
! PARAMETERS
!
      integer , parameter :: ilon = 240 , jlat = 121
      integer , parameter :: idtbc = 6
!
! Local variables
!
      integer :: i , it , j , nday , nhour , nmo , nyear
      real(sp) , dimension(jlat) :: lati
      real(sp) , dimension(ilon) :: loni
      integer :: idate , ierastart , ierrec , nsteps
      logical :: there
      real(sp) , dimension(ilon,jlat) :: sst
      character(256) :: inpfile
!
      if ( ssttyp=='ERSST' ) then
        there = .false.
        if ( (globidate1>=1989010100 .and. globidate1<=2009053118) .or. &
           & (globidate2>=1989010100 .and. globidate2<=2009053118) ) then
          inquire (file=trim(inpglob)//'/SST/sstERAIN.1989-2009.nc',    &
            &      exist=there)
          if ( .not.there ) then
            write(stderr,*) 'sstERAIN.1989-2009.nc is not available' ,  &
                 &' under ',trim(inpglob),'/SST/'
            call die('sst_ersst')
          end if
        end if
        if ( .not.there ) then
          write(stderr,*) 'ERSST Sea Surface Temp is just available' , &
               &' from 1989010100 to 2009053118'
          call die('sst_ersst')
        end if
      else if ( ssttyp=='ERSKT' ) then
        there = .false.
        if ( (globidate1>=1989010100 .and. globidate1<=2009053118) .or. &
           & (globidate2>=1989010100 .and. globidate2<=2009053118) ) then
          inquire (file=trim(inpglob)//'/SST/tskinERAIN.1989-2009.nc',  &
                 & exist=there)
          if ( .not.there ) then
            write(stderr,*) 'tskinERAIN.1989-2009.nc is not available' ,&
                 &' under ',trim(inpglob),'/SST/'
            call die('sst_ersst')
          end if
        end if
        if ( .not.there ) then
          write(stderr,*) 'ERSKT Skin Temperature is just available' ,  &
               &' from 1989010100 to 2009053118'
          call die('sst_ersst')
        end if
      else
        write (stderr,*) 'PLEASE SET right SSTTYP in regcm.in'
        write (stderr,*) 'Supported types are ERSST ERSKT'
        call die('sst_ersst')
      end if

      nsteps = idatediff(globidate2,globidate1)/idtbc + 1

      write (stdout,*) 'GLOBIDATE1 : ' , globidate1
      write (stdout,*) 'GLOBIDATE2 : ' , globidate2
      write (stdout,*) 'NSTEPS     : ' , nsteps

      call open_sstfile(globidate1)
 
!     ******    SET UP LONGITUDES AND LATITUDES FOR SST DATA
      do i = 1 , ilon
        loni(i) = float(i-1)*1.5
      end do
      do j = 1 , jlat
        lati(j) = -90. + 1.5*float(j-1)
      end do
 
      idate = globidate1
      ierastart = 1989010100
      do it = 1 , nsteps

        ierrec = idatediff(idate,ierastart)/idtbc+1

        if ( ssttyp=='ERSST' ) then
          inpfile = trim(inpglob)//'/SST/sstERAIN.1989-2009.nc'
          call sst_erain(ierrec,ilon,jlat,sst,inpfile)
        else if ( ssttyp=='ERSKT' ) then
          inpfile = trim(inpglob)//'/SST/tskinERAIN.1989-2009.nc'
          call skt_erain(ierrec,ilon,jlat,sst,inpfile)
        else
        end if

        call split_idate(idate, nyear, nmo, nday, nhour)
 
        call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
        write(stdout,*) &
            'XLON,XLAT,SST=' , xlon(1,1) , xlat(1,1) , sstmm(1,1)
 
!       ******           WRITE OUT SST DATA ON MM4 GRID
        call writerec(idate,.false.)
        write(stdout,*) 'WRITING OUT MM4 SST DATA:' , nmo , nyear

        call addhours(idate, idtbc)

      end do
 
      end subroutine sst_ersst
!
!-----------------------------------------------------------------------
!
      subroutine skt_erain(it,ilon,jlat,sst,pathaddname)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: it , ilon , jlat
      character(256) :: pathaddname
      intent (in) it , ilon , jlat , pathaddname
      real(sp) , dimension(ilon,jlat) :: sst
      intent (out) :: sst
!
! Local variables
!
      integer :: i , j , n
      logical :: there
      character(5) :: varname
      integer(2) , dimension(ilon,jlat) :: work
      integer :: istatus
!
      integer , dimension(10) , save :: icount , istart
      integer , save :: inet , ivar
      real(dp) , save :: xadd , xscale , xmiss
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers. The array 'sst'
!     will contain the unpacked data.
!
!     DATA ARRAY AND WORK ARRAY
!
      data varname/'skt'/
!
      if ( it==1 ) then
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          call die('skt_erain',trim(pathaddname)//' is not available',1)
        end if
        istatus = nf90_open(pathaddname,nf90_nowrite,inet)
        if ( istatus/=nf90_noerr ) then
          call die('skt_erain','Cannot open input file '// &
                   trim(pathaddname),1)
        end if
        istatus = nf90_inq_varid(inet,varname,ivar)
        if ( istatus/=nf90_noerr ) then
          call die('skt_erain','Cannot open input file '// &
                   trim(pathaddname),1)
        end if
        istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
        istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
        istatus = nf90_get_att(inet,ivar,'_FillValue',xmiss)
        istart(1) = 1
        istart(2) = 1
        icount(1) = 240
        icount(2) = 121
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
! 
      istart(3) = it
      icount(3) = 1
      istatus = nf90_get_var(inet,ivar,work,istart,icount)
      if ( istatus/=nf90_noerr ) then
        call die('skt_erain','Cannot get '//varname//' from file', &
                 istatus,nf90_strerror(istatus),0)
      end if
!
      do j = 1 , jlat
        do i = 1 , ilon
          if (work(i,j)/=xmiss) then
            sst(i,jlat+1-j) = work(i,j)*xscale + xadd
          else
            sst(i,jlat+1-j)=-9999
          end if
        end do
      end do
!
      end subroutine skt_erain
!
!-----------------------------------------------------------------------
!
      subroutine sst_erain(it,ilon,jlat,sst,pathaddname)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: it , ilon , jlat
      character(256) :: pathaddname
      intent (in) it , ilon , jlat , pathaddname
      real(sp) , dimension(ilon,jlat) :: sst
      intent (out) :: sst
!
! Local variables
!
      integer :: i , j , n
      logical :: there
      character(5) :: varname
      integer(2) , dimension(ilon,jlat) :: work
      integer :: istatus
!
      integer , save :: inet , ivar
      real(dp) , save :: xadd , xscale , xmiss
      integer , dimension(10) , save :: icount , istart
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers. The array 'sst'
!     will contain the unpacked data.
!
!     DATA ARRAY AND WORK ARRAY
!
      data varname/'sst'/
!
      if ( it==1 ) then
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          call die('sst_erain',trim(pathaddname)//' is not available',1)
        end if
        istatus = nf90_open(pathaddname,nf90_nowrite,inet)
        if ( istatus/=nf90_noerr ) then
          call die('sst_erain','Cannot open input file '// &
                   trim(pathaddname)//' : '//nf90_strerror(istatus), &
                   istatus)
        end if
        istatus = nf90_inq_varid(inet,varname,ivar)
        if ( istatus/=nf90_noerr ) then
          call die('sst_erain','Cannot find variable '//trim(varname)// &
                   ' in file '//trim(pathaddname)//' : '//              &
                   nf90_strerror(istatus),istatus)
        end if
        istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
        istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
        istatus = nf90_get_att(inet,ivar,'_FillValue',xmiss)
        istart(1) = 1
        istart(2) = 1
        icount(1) = 240
        icount(2) = 121
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
!
      istart(3) = it
      icount(3) = 1
      istatus = nf90_get_var(inet,ivar,work,istart,icount)
      if ( istatus/=nf90_noerr ) then
        call die('sst_erain','Cannot read '//trim(varname)// &
                 ' from file '//trim(pathaddname)//' : '//              &
                 nf90_strerror(istatus),istatus)
      end if
!
      do j = 1 , jlat
        do i = 1 , ilon
          if (work(i,j)/=xmiss) then
            sst(i,jlat+1-j) = work(i,j)*xscale + xadd
          else
            sst(i,jlat+1-j)=-9999
          end if
        end do
      end do
!
      end subroutine sst_erain
!
      end module mod_sst_ersst
