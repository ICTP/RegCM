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

      module mod_rdldtr

      use mod_block

      contains

      subroutine rdldtr(inpter,ntypec,nveg,ntex,aertyp,ibyte)

      implicit none
!
! Dummy arguments
!
      integer :: ntypec , ibyte , nveg , ntex
      character(7) :: aertyp
      character(*) :: inpter
      intent(in) inpter , ntypec , ibyte , nveg , ntex , aertyp
!
! Local variables
!
      integer :: nxmax , nymax
      real(4) :: center , rlat , rlon
      character(2) :: char_2
      character(256) :: filcat , filelev
      character(256) :: filtext
      integer :: ihmax , ilat1 , ilat2 , irec , irect , j , jhmax , k , &
               & lrec
      logical :: there
      character(1) , allocatable , dimension(:,:) :: ch_cat
      character(2) , allocatable , dimension(:) :: ch_htsd , ch_topo
      character(1) , allocatable , dimension(:,:) :: ch_tex
!
      nxmax = 21600/ntypec
      nymax = nxmax/2

      allocate(ch_cat(nxmax,nveg))
      allocate(ch_tex(nxmax,ntex))
      allocate(ch_topo(nxmax))
      allocate(ch_htsd(nxmax))

      if ( ntypec<10 ) then
        write (filelev,99001) trim(inpter) , ntypec
      else if ( ntypec<100 ) then
        write (filelev,99002) trim(inpter) , ntypec
      else
        write (*,*) 'For terrain, ntypec is not set correctly ' , ntypec
        stop 'subroutine RDLDTR'
      end if
      inquire (file=filelev,exist=there)
      if ( .not.there ) then
        print * , 'ERROR OPENING ' , filelev ,                          &
             &' FILE:  FILE DOES NOT EXIST'
        stop '4810 IN SUBROUTINE RDLDTR'
      end if
      open (46,file=filelev,form='unformatted',recl=nxmax*ibyte/2,      &
          & access='direct')
 
      if ( ntypec<10 ) then
        write (filcat,99003) trim(inpter) , ntypec
      else if ( ntypec<100 ) then
        write (filcat,99004) trim(inpter) , ntypec
      else
        write (*,*) 'For landuse, ntypec is not set correctly ' ,     &
                  & ntypec
        stop 'subroutine RDLDTR'
      end if
      inquire (file=filcat,exist=there)
      if ( .not.there ) then
        print * , 'ERROR OPENING ' , filcat ,                         &
             &' FILE:  FILE DOES NOT EXIST'
        stop '4820 IN SUBROUTINE RDLDTR'
      end if
      open (47,file=filcat,form='unformatted',recl=nxmax*nveg*ibyte/4,&
          & access='direct')
      if ( aertyp(7:7)=='1' ) then
        if ( ntypec<10 ) then
          write (filtext,99011) ntypec
        else if ( ntypec<100 ) then
          write (filtext,99012) ntypec
        else
          write (*,*) 'Error input of ntypec: ' , ntypec
        end if
        inquire (file=trim(inpter)//'/SURFACE/'//filtext,exist=there)
        if ( .not.there ) then
          print * , 'ERROR OPENING ' ,                                  &
              & trim(inpter)//'/SURFACE/'//filtext ,                    &
              & ' FILE: FILE DOES NOT EXIST'
          stop '4830 IN SUBROUTINE RDLDTR'
        end if
        open (45,file=trim(inpter)//'/SURFACE/'//filtext,               &
           & form='unformatted',recl=nxmax*ntex*ibyte/4,access='direct')
      end if
 
      lrec = 0
      rewind (48)
      center = xnc
      ilat1 = nint((90.-xmaxlat)/xnc)
      ilat1 = max(min(nymax,ilat1),1)
      ilat2 = nint((90.-xminlat)/xnc) + 1
      ilat2 = max(min(nymax,ilat2),1)
      do irec = ilat1 , ilat2
        rlat = 90. - center*irec + center/2.
        irect = nymax - irec + 1
        read (46,rec=irect) (ch_topo(j),j=1,nxmax)
!       print*,'   ELEVATION READ IN', ilat1,ilat2,irec
        read (46,rec=nymax+irect) (ch_htsd(j),j=1,nxmax)
!       print*,'   ELEVATION STD DEV READ IN', ilat1,ilat2,irec
        read (47,rec=irect) ((ch_cat(j,k),k=1,nveg),j=1,nxmax)
        if ( aertyp(7:7)=='1' ) read (45,rec=irec)                      &
                                    & ((ch_tex(j,k),k=1,ntex),j=1,nxmax)
!       print*,'   LANDUSE READ IN', ilat1,ilat2,irec
!........process slice and store in array stores
        do j = 1 , nxmax
          rlon = j*center - center/2. - 180.
          if ( rlon >  180.0) rlon = rlon - 360.0
          if ( rlon < -180.0) rlon = rlon + 360.0
          if ((lcrosstime .and. &
              (rlon>=xminlon .or. rlon<=xmaxlon)) .or. &
              (rlon>=xminlon .and. rlon<=xmaxlon)) then
            lrec = lrec + 1
            stores(1) = rlat
            stores(2) = rlon
            char_2 = ch_topo(j)
            if ( ichar(char_2(1:1))*256+ichar(char_2(2:2))-1000<-200 )  &
               & then      ! OCEAN/UNDEFINED
              stores(3) = 0.0
            else
              stores(3) = ichar(char_2(1:1))*256 + ichar(char_2(2:2))   &
                        & - 1000
            end if
            char_2 = ch_htsd(j)
            stores(4) = ichar(char_2(1:1))*256 + ichar(char_2(2:2))
            do k = 1 , nveg
              stores(k+4) = float(ichar(ch_cat(j,k)))
            end do
            if ( aertyp(7:7)=='1' ) then
              do k = 1 , ntex
                stores(nveg+4+k) = float(ichar(ch_tex(j,k)))
              end do
            end if
            write (48) stores
            yobs(lrec) = rlat
            xobs(lrec) = rlon
            ht(lrec) = stores(3)
            htsd(lrec) = stores(4)
            ht2(lrec) = stores(4)**2 + stores(3)**2
            if (lcrosstime) then
              if ( mod((xobs(lrec)+360.0D0),360.0D0)<grdlnmn ) &
                 grdlnmn = xobs(lrec)
              if ( mod((xobs(lrec)+360.0D0),360.0D0)>(grdlnma+360.0) ) &
                 grdlnma = xobs(lrec)
            else
              if ( xobs(lrec)<grdlnmn ) grdlnmn = xobs(lrec)
              if ( xobs(lrec)>grdlnma ) grdlnma = xobs(lrec)
            end if
            if ( yobs(lrec)<grdltmn ) grdltmn = yobs(lrec)
            if ( yobs(lrec)>grdltma ) grdltma = yobs(lrec)
          end if
        end do
      end do
 
      close (46)
      close (47)
!
      print 99013 , lrec
      ihmax = (xmaxlat-xminlat)/xnc
      if (xminlon > 0.0 .and. xmaxlon < 0.0) then
        jhmax = ((xmaxlon+360.0)-xminlon)/xnc
      else
        jhmax = (xmaxlon-xminlon)/xnc
      end if
      if ( ihmax>iter .or. jhmax>jter ) print 99014 , iter , ihmax ,    &
         & jter , jhmax
      if ( ihmax*jhmax>iblk ) print 99015 , ihmax*jhmax , iblk
 
      nobs = lrec
 
      rewind (48)

      deallocate(ch_cat)
      deallocate(ch_tex)
      deallocate(ch_topo)
      deallocate(ch_htsd)

99001 format (a,'/SURFACE/GTOPO30_',i1,'MIN.dat')
99002 format (a,'/SURFACE/GTOPO30_',i2,'MIN.dat')
99003 format (a,'/SURFACE/GLCC',i1,'MIN_BATS.dat')
99004 format (a,'/SURFACE/GLCC',i2,'MIN_BATS.dat')
99011 format ('SOILCAT.0',i1)
99012 format ('SOILCAT.',i2)
99013 format (1x,i10,' terrain heights read from land use volume')
99014 format (1x,'***array dimension error***',/,'     iter = ',i5,     &
             &' must be greater than ',i5,10x,'jter = ',i5,             &
             &' must be greater than ',i5)
99015 format (1x,'***array dimension error***',/,'  ihmax*jhmax = ',i8, &
             &' must be greater than ',i8,10x,'iblk = ',i8)
      end subroutine rdldtr

      subroutine rdldtr_nc(inpter,ntypec,nveg,ntex,aertyp)

      use netcdf

      implicit none
!
! Dummy arguments
!
      integer :: ntypec , nveg , ntex
      character(7) :: aertyp
      character(*) :: inpter
      intent(in) inpter , ntypec , nveg , ntex , aertyp
!
! Local variables
!
      integer :: nxmax , nymax
      real(4) :: center , rlat , rlon
      character(256) :: filcat
      character(256) :: filelev , filtext
      integer , dimension(3) :: icount3 , istart3
      integer , dimension(4) :: icount4 , istart4
      integer :: ihmax , ilat1 , ilat2 , irec , istat ,  j , jhmax , k ,&
               & lrec , ncid_cat , ncid_lev , ncid_tex
      integer :: idtopo , idhtsd , idlufrac
      logical :: there
      integer(1) , allocatable , dimension(:,:) :: icat
      integer(1) , allocatable , dimension(:,:) :: itex
      integer(2) , allocatable , dimension(:) :: ihtsd , itopo
!
      nxmax = 21600/ntypec
      nymax = nxmax/2

      allocate(icat(nxmax,nveg))
      allocate(itex(nxmax,ntex))
      allocate(itopo(nxmax))
      allocate(ihtsd(nxmax))

      if ( ntypec<10 ) then
        write (filelev,99001) trim(inpter), ntypec
      else if ( ntypec<100 ) then
        write (filelev,99002) trim(inpter), ntypec
      else
        write (*,*) 'For terrain, ntypec is not set correctly ' , ntypec
        stop 'subroutine RDLDTR_nc'
      end if
      inquire (file=filelev,exist=there)
      if ( .not.there ) then
        print * , 'ERROR OPENING ' , filelev ,                          &
             &' FILE:  FILE DOES NOT EXIST'
        stop '4810 IN SUBROUTINE RDLDTR_nc'
      end if
      istat = nf90_open(filelev,nf90_nowrite,ncid_lev)
      if ( istat /= nf90_noerr ) then
        write (6, *) 'Error opening ', filelev
        write (6, *) nf90_strerror(istat)
        stop
      end if
      istat = nf90_inq_varid(ncid_lev, 'HT', idtopo)
      if ( istat /= nf90_noerr ) then
        write (6, *) 'Var HT not found'
        write (6, *) nf90_strerror(istat)
        stop
      end if
      istat = nf90_inq_varid(ncid_lev, 'HTSD', idhtsd)
      if ( istat /= nf90_noerr ) then
        write (6, *) 'Var HTSD not found'
        write (6, *) nf90_strerror(istat)
        stop
      end if
 
      if ( ntypec<10 ) then
        write (filcat,99003) trim(inpter), ntypec
      else if ( ntypec<100 ) then
        write (filcat,99004) trim(inpter), ntypec
      else
        write (*,*) 'For landuse, ntypec is not set correctly ' ,     &
                  & ntypec
        stop 'subroutine RDLDTR_nc'
      end if
      inquire (file=filcat,exist=there)
      if ( .not.there ) then
        print * , 'ERROR OPENING ' , filcat ,                         &
             &' FILE:  FILE DOES NOT EXIST'
        stop '4820 IN SUBROUTINE RDLDTR_nc'
      end if
      istat = nf90_open(filcat,nf90_nowrite,ncid_cat)
      if ( istat /= nf90_noerr ) then
        write (6, *) 'Error opening ', filcat
        write (6, *) nf90_strerror(istat)
        stop
      end if
      istat = nf90_inq_varid(ncid_cat, 'lufrac', idlufrac)
      if ( istat /= nf90_noerr ) then
        write (6, *) 'Var lufrac not found'
        write (6, *) nf90_strerror(istat)
        stop
      end if
      if ( aertyp(7:7)=='1' ) then
        if ( ntypec<10 ) then
          write (filtext,99009) ntypec
        else if ( ntypec<100 ) then
          write (filtext,99010) ntypec
        else
          write (*,*) 'Error input of ntypec: ' , ntypec
        end if
        inquire (file=filtext,exist=there)
        if ( .not.there ) then
          print * , 'ERROR OPENING ' , filtext ,                        &
               &' FILE: FILE DOES NOT EXIST'
          stop '4830 IN SUBROUTINE RDLDTR_nc'
        end if
        istat = nf90_open(filtext,nf90_nowrite,ncid_tex)
        if ( istat /= nf90_noerr ) then
          write (6, *) 'Error opening ', filtext
          write (6, *) nf90_strerror(istat)
          stop
        end if
      end if
 
      lrec = 0
      rewind (48)
      center = xnc
      ilat1 = nint((90.-xmaxlat)/xnc)
      ilat1 = max(min(nymax,ilat1),1)
      ilat2 = nint((90.-xminlat)/xnc) + 1
      ilat2 = max(min(nymax,ilat2),1)
      istart3(1) = 1
      icount3(1) = nxmax
      istart4(1) = 1
      icount4(1) = nxmax
      istart3(3) = 1
      icount3(3) = 1
      istart4(4) = 1
      icount4(4) = 1
      do irec = ilat1 , ilat2
        rlat = 90. - center*irec + center/2.
        istart3(2) = nymax - irec + 1
        icount3(2) = 1
        istart4(2) = nymax - irec + 1
        icount4(2) = 1
        istat = nf90_get_var(ncid_lev,idtopo,itopo,istart3,icount3)
        if ( istat /= nf90_noerr ) then
          write (6, *) 'Error reading var itopo'
          write (6, *) nf90_strerror(istat)
          stop
        end if
        istat = nf90_get_var(ncid_lev,idhtsd,ihtsd,istart3,icount3)
        if ( istat /= nf90_noerr ) then
          write (6, *) 'Error reading var ihtsd'
          write (6, *) nf90_strerror(istat)
          stop
        end if
 
        istart4(3) = 1
        icount4(3) = nveg
        istat = nf90_get_var(ncid_cat,idlufrac,icat,istart4,icount4)
        if ( istat /= nf90_noerr ) then
          write (6, *) 'Error reading var icat'
          write (6, *) nf90_strerror(istat)
          stop
        end if
        if ( aertyp(7:7)=='1' ) then
          istart3(3) = 1
          icount3(3) = ntex
          istat = nf90_get_var(ncid_tex,4,itex,istart3,icount3)
          if ( istat /= nf90_noerr ) then
            write (6, *) 'Error reading var itex'
            write (6, *) nf90_strerror(istat)
            stop
          end if
        end if
!       print*,'   LANDUSE READ IN', ilat1,ilat2,irec
!........process slice and store in array stores
        do j = 1 , nxmax
          rlon = j*center - center/2. - 180.
          if ( rlon >  180.0) rlon = rlon - 360.0
          if ( rlon < -180.0) rlon = rlon + 360.0
          if ((lcrosstime .and. &
              (rlon>=xminlon .or. rlon<=xmaxlon)) .or. &
              (rlon>=xminlon .and. rlon<=xmaxlon)) then
            lrec = lrec + 1
            stores(1) = rlat
            stores(2) = rlon
            stores(3) = itopo(j)*1.0
            stores(4) = ihtsd(j)*1.0
            do k = 1 , nveg
              stores(k+4) = icat(j,k)*1.0
            end do
            if ( aertyp(7:7)=='1' ) then
              do k = 1 , ntex
                stores(nveg+4+k) = itex(j,k)*1.0
              end do
            end if
            write (48) stores
            yobs(lrec) = rlat
            xobs(lrec) = rlon
            ht(lrec) = stores(3)
            htsd(lrec) = stores(4)
            ht2(lrec) = stores(4)**2 + stores(3)**2
            if (lcrosstime) then
              if ( mod((xobs(lrec)+360.0D0),360.0D0)<grdlnmn ) &
                 grdlnmn = xobs(lrec)
              if ( mod((xobs(lrec)+360.0D0),360.0D0)>(grdlnma+360.0) ) &
                 grdlnma = xobs(lrec)
            else
              if ( xobs(lrec)<grdlnmn ) grdlnmn = xobs(lrec)
              if ( xobs(lrec)>grdlnma ) grdlnma = xobs(lrec)
            end if
            if ( yobs(lrec)<grdltmn ) grdltmn = yobs(lrec)
            if ( yobs(lrec)>grdltma ) grdltma = yobs(lrec)
          end if
        end do
      end do
 
      istat = nf90_close(ncid_lev)
      if ( istat /= nf90_noerr ) then
        write (6, *) 'Error close'
        write (6, *) nf90_strerror(istat)
        stop
      end if
      istat = nf90_close(ncid_cat)
      if ( istat /= nf90_noerr ) then
        write (6, *) 'Error close'
        write (6, *) nf90_strerror(istat)
        stop
      end if
      if ( aertyp(7:7)=='1' ) then
        istat = nf90_close(ncid_tex)
        if ( istat /= nf90_noerr ) then
          write (6, *) 'Error close'
          write (6, *) nf90_strerror(istat)
          stop
        end if
      end if
!
      print 99011 , lrec
      ihmax = (xmaxlat-xminlat)/xnc
      if (lcrosstime) then
        jhmax = ((xmaxlon+360.0)-xminlon)/xnc
      else
        jhmax = (xmaxlon-xminlon)/xnc
      end if
      if ( ihmax>iter .or. jhmax>jter ) print 99012 , iter , ihmax ,    &
         & jter , jhmax
      if ( ihmax*jhmax>iblk ) print 99013 , ihmax*jhmax , iblk
 
      nobs = lrec
 
      rewind (48)

      deallocate(icat)
      deallocate(itex)
      deallocate(itopo)
      deallocate(ihtsd)

99001 format (a,'/SURFACE/GTOPO30_',i1,'min.nc')
99002 format (a,'/SURFACE/GTOPO30_',i2,'min.nc')
99003 format (a,'/SURFACE/GLCC_BATS_',i1,'min.nc')
99004 format (a,'/SURFACE/GLCC_BATS_',i2,'min.nc')
99009 format (a,'/SURFACE/SOILCAT_',i1,'min.nc')
99010 format (a,'/SURFACE/SOILCAT_',i2,'min.nc')
99011 format (1x,i10,' terrain heights read from land use volume')
99012 format (1x,'***array dimension error***',/,'     iter = ',i5,     &
             &' must be greater than ',i5,10x,'jter = ',i5,             &
             &' must be greater than ',i5)
99013 format (1x,'***array dimension error***',/,'  ihmax*jhmax = ',i8, &
             &' must be greater than ',i8,10x,'iblk = ',i8)
      end subroutine rdldtr_nc
!
      end module mod_rdldtr
