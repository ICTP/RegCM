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

      subroutine rdldtr(ntypec,nveg,ntex,lsmtyp,aertyp,ibyte)

      implicit none
!
! Dummy arguments
!
      integer :: ntypec , ibyte , nveg , ntex
      character(4) :: lsmtyp
      character(7) :: aertyp
      intent(in) ntypec , ibyte , nveg , ntex , lsmtyp , aertyp
!
! Local variables
!
      integer :: nxmax , nymax
      real(4) :: center , rlat , rlon
      character(2) :: char_2
      character(36) :: filcat , filelev
      character(7) :: filclay , filsand
      character(10) :: filtext
      character(11) :: filusgs
      integer :: ihmax , ilat1 , ilat2 , irec , irect , j , jhmax , k , &
               & lrec
      logical :: there
      character(1) , allocatable , dimension(:,:) :: ch_cat
      character(2) , allocatable , dimension(:) :: ch_htsd , ch_topo
      character(1) , allocatable , dimension(:,:) :: ch_tex
      integer(2) , allocatable , dimension(:,:) :: iclay , isand
      integer(2) , allocatable , dimension(:,:) :: iusgs
!
      nxmax = 21600/ntypec
      nymax = nxmax/2

      allocate(ch_cat(nxmax,nveg))
      allocate(iusgs(nxmax,nveg))
      allocate(iclay(nxmax,2))
      allocate(isand(nxmax,2))
      allocate(ch_tex(nxmax,ntex))
      allocate(ch_topo(nxmax))
      allocate(ch_htsd(nxmax))

      if ( ntypec<10 ) then
        write (filelev,99001) ntypec
      else if ( ntypec<100 ) then
        write (filelev,99002) ntypec
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
 
      if ( lsmtyp=='BATS' ) then
        if ( ntypec<10 ) then
          write (filcat,99003) ntypec
        else if ( ntypec<100 ) then
          write (filcat,99004) ntypec
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
      else if ( lsmtyp=='USGS' ) then
        if ( ntypec<10 ) then
          write (filusgs,99005) ntypec
          write (filsand,99006) ntypec
          write (filclay,99007) ntypec
        else if ( ntypec<100 ) then
          write (filusgs,99008) ntypec
          write (filsand,99009) ntypec
          write (filclay,99010) ntypec
        else
          write (*,*) 'Error input of ntypec: ' , ntypec
        end if
        if ( ntypec==2 ) then
          inquire (file='../DATA/SURFACE/'//filusgs//'A',exist=there)
          if ( .not.there ) then
            print * , 'ERROR OPENING ' , '../DATA/SURFACE/'//filusgs ,  &
                 &' FILE: FILE DOES NOT EXIST'
            stop '4820 IN SUBROUTINE RDLDTR'
          end if
          inquire (file='../DATA/SURFACE/'//filusgs//'B',exist=there)
          if ( .not.there ) then
            print * , 'ERROR OPENING ' , '../DATA/SURFACE/'//filusgs ,  &
                 &' FILE: FILE DOES NOT EXIST'
            stop '4820 IN SUBROUTINE RDLDTR'
          end if
        else
          inquire (file='../DATA/SURFACE/'//filusgs,exist=there)
          if ( .not.there ) then
            print * , 'ERROR OPENING ' , '../DATA/SURFACE/'//filusgs ,  &
                 &' FILE: FILE DOES NOT EXIST'
            stop '4820 IN SUBROUTINE RDLDTR'
          end if
        end if
        inquire (file='../DATA/SURFACE/'//filsand,exist=there)
        if ( .not.there ) then
          print * , 'ERROR OPENING ' , '../DATA/SURFACE/'//filsand ,    &
               &' FILE: FILE DOES NOT EXIST'
          stop '4820 IN SUBROUTINE RDLDTR'
        end if
        inquire (file='../DATA/SURFACE/'//filclay,exist=there)
        if ( .not.there ) then
          print * , 'ERROR OPENING ' , '../DATA/SURFACE/'//filclay ,    &
               &' FILE: FILE DOES NOT EXIST'
          stop '4820 IN SUBROUTINE RDLDTR'
        end if
        if ( ntypec==2 ) then
          open (41,file='../DATA/SURFACE/'//filusgs//'A',               &
               &form='unformatted',recl=nxmax*ibyte/2,access='direct')
          open (44,file='../DATA/SURFACE/'//filusgs//'B',               &
               &form='unformatted',recl=nxmax*ibyte/2,access='direct')
        else
          open (41,file='../DATA/SURFACE/'//filusgs,form='unformatted', &
              & recl=nxmax*ibyte/2,access='direct')
        end if
        open (42,file='../DATA/SURFACE/'//filsand,form='unformatted',   &
            & recl=nxmax*ibyte/2,access='direct')
        open (43,file='../DATA/SURFACE/'//filclay,form='unformatted',   &
            & recl=nxmax*ibyte/2,access='direct')
      else
      end if
      if ( aertyp(7:7)=='1' ) then
        if ( ntypec<10 ) then
          write (filtext,99011) ntypec
        else if ( ntypec<100 ) then
          write (filtext,99012) ntypec
        else
          write (*,*) 'Error input of ntypec: ' , ntypec
        end if
        inquire (file='../DATA/SURFACE/'//filtext,exist=there)
        if ( .not.there ) then
          print * , 'ERROR OPENING ' , '../DATA/SURFACE/'//filtext ,    &
               &' FILE: FILE DOES NOT EXIST'
          stop '4830 IN SUBROUTINE RDLDTR'
        end if
        open (45,file='../DATA/SURFACE/'//filtext,form='unformatted',   &
            & recl=nxmax*ntex*ibyte/4,access='direct')
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
        if ( lsmtyp=='BATS' ) then
          read (47,rec=irect) ((ch_cat(j,k),k=1,nveg),j=1,nxmax)
        else if ( lsmtyp=='USGS' ) then
          if ( ntypec==2 ) then
            do k = 1 , 13
              read (41,rec=nymax*(k-1)+irec) (iusgs(j,k),j=1,nxmax)
            end do
            do k = 14 , nveg
              read (44,rec=nymax*(k-14)+irec) (iusgs(j,k),j=1,nxmax)
            end do
          else
            do k = 1 , nveg
              read (41,rec=nymax*(k-1)+irec) (iusgs(j,k),j=1,nxmax)
            end do
          end if
          read (42,rec=irec) (isand(j,1),j=1,nxmax)
          read (42,rec=nymax+irec) (isand(j,2),j=1,nxmax)
          read (43,rec=irec) (iclay(j,1),j=1,nxmax)
          read (43,rec=nymax+irec) (iclay(j,2),j=1,nxmax)
        else
        end if
        if ( aertyp(7:7)=='1' ) read (45,rec=irec)                      &
                                    & ((ch_tex(j,k),k=1,ntex),j=1,nxmax)
!       print*,'   LANDUSE READ IN', ilat1,ilat2,irec
!........process slice and store in array stores
        do j = 1 , nxmax
          rlon = j*center - center/2. - 180.
          if ( xminlon<-180. .and. rlon>0. ) rlon = rlon - 360.
          if ( xmaxlon>180. .and. rlon<0. ) rlon = rlon + 360.
          if ( rlon>=xminlon .and. rlon<=xmaxlon ) then
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
            if ( lsmtyp=='BATS' ) then
              do k = 1 , nveg
                stores(k+4) = float(ichar(ch_cat(j,k)))
              end do
            else if ( lsmtyp=='USGS' ) then
              do k = 1 , nveg
                stores(k+4) = float(iusgs(j,k))*50./32767. + 50.
              end do
              stores(nveg+5) = float(isand(j,1))*50./32767. + 50.
              stores(nveg+6) = float(isand(j,2))*50./32767. + 50.
              stores(nveg+7) = float(iclay(j,1))*50./32767. + 50.
              stores(nveg+8) = float(iclay(j,2))*50./32767. + 50.
            else
            end if
            if ( aertyp(7:7)=='1' ) then
              if ( lsmtyp=='BATS' ) then
                do k = 1 , ntex
                  stores(nveg+4+k) = float(ichar(ch_tex(j,k)))
                end do
              else if ( lsmtyp=='USGS' ) then
                do k = 1 , ntex
                  stores(nveg+8+k) = float(ichar(ch_tex(j,k)))
                end do
              else
              end if
            end if
            write (48) stores
!add
            yobs(lrec) = rlat
            xobs(lrec) = rlon
            ht(lrec) = stores(3)
            htsd(lrec) = stores(4)
            ht2(lrec) = stores(4)**2 + stores(3)**2
            if ( grdlnmn<=-180.0 .and. xobs(lrec)>0.0 ) xobs(lrec)      &
               & = xobs(lrec) - 360.
            if ( xobs(lrec)<grdlnmn ) grdlnmn = xobs(lrec)
            if ( yobs(lrec)<grdltmn ) grdltmn = yobs(lrec)
!add_
          end if
        end do
      end do
 
      close (46)
      if ( lsmtyp=='BATS' ) then
        close (47)
      else if ( lsmtyp=='USGS' ) then
        close (41)
        close (42)
        close (43)
        if ( ntypec==2 ) close (44)
      else
      end if
!
      print 99013 , lrec
      ihmax = (xmaxlat-xminlat)/xnc
      jhmax = (xmaxlon-xminlon)/xnc
      if ( ihmax>iter .or. jhmax>jter ) print 99014 , iter , ihmax ,    &
         & jter , jhmax
      if ( ihmax*jhmax>iblk ) print 99015 , ihmax*jhmax , iblk
 
      nobs = lrec
 
      rewind (48)

      deallocate(ch_cat)
      deallocate(iusgs)
      deallocate(iclay)
      deallocate(isand)
      deallocate(ch_tex)
      deallocate(ch_topo)
      deallocate(ch_htsd)

99001 format ('../DATA/SURFACE/GTOPO30_',i1,'MIN.dat')
99002 format ('../DATA/SURFACE/GTOPO30_',i2,'MIN.dat')
99003 format ('../DATA/SURFACE/GLCC',i1,'MIN_BATS.dat')
99004 format ('../DATA/SURFACE/GLCC',i2,'MIN_BATS.dat')
99005 format ('VEG-USGS.0',i1)
99006 format ('SAND.0',i1)
99007 format ('CLAY.0',i1)
99008 format ('VEG-USGS.',i2)
99009 format ('SAND.',i2)
99010 format ('CLAY.',i2)
99011 format ('SOILCAT.0',i1)
99012 format ('SOILCAT.',i2)
99013 format (1x,i10,' terrain heights read from land use volume')
99014 format (1x,'***array dimension error***',/,'     iter = ',i5,     &
             &' must be greater than ',i5,10x,'jter = ',i5,             &
             &' must be greater than ',i5)
99015 format (1x,'***array dimension error***',/,'  ihmax*jhmax = ',i8, &
             &' must be greater than ',i8,10x,'iblk = ',i8)
      end subroutine rdldtr

      subroutine rdldtr_nc(ntypec,nveg,ntex,lsmtyp,aertyp)

      use netcdf

      implicit none
!
! Dummy arguments
!
      integer :: ntypec , nveg , ntex
      character(4) :: lsmtyp
      character(7) :: aertyp
      intent(in) ntypec , nveg , ntex , lsmtyp , aertyp
!
! Local variables
!
      integer :: nxmax , nymax
      real(4) :: center , rlat , rlon
      character(64) :: filcat , filsandclay , filusgs
      character(64) :: filelev , filtext
      integer , dimension(3) :: icount3 , istart3
      integer , dimension(4) :: icount4 , istart4
      integer :: ihmax , ilat1 , ilat2 , irec , istat ,  j , jhmax , k ,&
               & lrec , ncid_cat , ncid_lev , ncid_sandclay , ncid_tex ,&
               & ncid_usgs
      integer :: idtopo , idhtsd , idlufrac
      logical :: there
      integer(1) , allocatable , dimension(:,:) :: icat , iusgs
      integer(1) , allocatable , dimension(:,:) :: iclay , isand
      integer(1) , allocatable , dimension(:,:) :: itex
      integer(2) , allocatable , dimension(:) :: ihtsd , itopo
!
      nxmax = 21600/ntypec
      nymax = nxmax/2

      allocate(icat(nxmax,nveg))
      allocate(iusgs(nxmax,nveg))
      allocate(iclay(nxmax,2))
      allocate(isand(nxmax,2))
      allocate(itex(nxmax,ntex))
      allocate(itopo(nxmax))
      allocate(ihtsd(nxmax))

      if ( ntypec<10 ) then
        write (filelev,99001) ntypec
      else if ( ntypec<100 ) then
        write (filelev,99002) ntypec
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
 
      if ( lsmtyp=='BATS' ) then
        if ( ntypec<10 ) then
          write (filcat,99003) ntypec
        else if ( ntypec<100 ) then
          write (filcat,99004) ntypec
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
      else if ( lsmtyp=='USGS' ) then
        if ( ntypec<10 ) then
          write (filusgs,99005) ntypec
          write (filsandclay,99006) ntypec
        else if ( ntypec<100 ) then
          write (filusgs,99007) ntypec
          write (filsandclay,99008) ntypec
        else
          write (*,*) 'Error input of ntypec: ' , ntypec
        end if
        inquire (file=filusgs,exist=there)
        if ( .not.there ) then
          print * , 'ERROR OPENING ' , filusgs ,                        &
               &' FILE: FILE DOES NOT EXIST'
          stop '4820 IN SUBROUTINE RDLDTR_nc'
        end if
        inquire (file=filsandclay,exist=there)
        if ( .not.there ) then
          print * , 'ERROR OPENING ' , filsandclay ,                    &
               &' FILE: FILE DOES NOT EXIST'
          stop '4820 IN SUBROUTINE RDLDTR_nc'
        end if
        istat = nf90_open(filusgs,nf90_nowrite,ncid_usgs)
        if ( istat /= nf90_noerr ) then
          write (6, *) 'Error opening ', filusgs
          write (6, *) nf90_strerror(istat)
          stop
        end if
        istat = nf90_open(filsandclay,nf90_nowrite,ncid_sandclay)
        if ( istat /= nf90_noerr ) then
          write (6, *) 'Error opening ', filsandclay
          write (6, *) nf90_strerror(istat)
          stop
        end if
      else
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
 
        if ( lsmtyp=='BATS' ) then
          istart4(3) = 1
          icount4(3) = nveg
          istat = nf90_get_var(ncid_cat,idlufrac,icat,istart4,icount4)
          if ( istat /= nf90_noerr ) then
            write (6, *) 'Error reading var icat'
            write (6, *) nf90_strerror(istat)
            stop
          end if
        else if ( lsmtyp=='USGS' ) then
          istart3(3) = 1
          icount3(3) = nveg
          istat = nf90_get_var(ncid_usgs,4,iusgs,istart3,icount3)
          if ( istat /= nf90_noerr ) then
            write (6, *) 'Error reading var iusgs'
            write (6, *) nf90_strerror(istat)
            stop
          end if
          istart3(3) = 1
          icount3(3) = 2
          istat = nf90_get_var(ncid_sandclay,4,isand,istart3,icount3)
          if ( istat /= nf90_noerr ) then
            write (6, *) 'Error reading var isand'
            write (6, *) nf90_strerror(istat)
            stop
          end if
          istat = nf90_get_var(ncid_sandclay,5,iclay,istart3,icount3)
          if ( istat /= nf90_noerr ) then
            write (6, *) 'Error reading var iclay'
            write (6, *) nf90_strerror(istat)
            stop
          end if
        else
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
          if ( xminlon<-180. .and. rlon>0. ) rlon = rlon - 360.
          if ( xmaxlon>180. .and. rlon<0. ) rlon = rlon + 360.
          if ( rlon>=xminlon .and. rlon<=xmaxlon ) then
            lrec = lrec + 1
            stores(1) = rlat
            stores(2) = rlon
            stores(3) = itopo(j)*1.0
            stores(4) = ihtsd(j)*1.0
            if ( lsmtyp=='BATS' ) then
              do k = 1 , nveg
                stores(k+4) = icat(j,k)*1.0
              end do
            else if ( lsmtyp=='USGS' ) then
              do k = 1 , nveg
                stores(k+4) = iusgs(j,k)*1.0
              end do
              stores(nveg+5) = isand(j,1)*1.0
              stores(nveg+6) = isand(j,2)*1.0
              stores(nveg+7) = iclay(j,1)*1.0
              stores(nveg+8) = iclay(j,2)*1.0
            else
            end if
            if ( aertyp(7:7)=='1' ) then
              if ( lsmtyp=='BATS' ) then
                do k = 1 , ntex
                  stores(nveg+4+k) = itex(j,k)*1.0
                end do
              else if ( lsmtyp=='USGS' ) then
                do k = 1 , ntex
                  stores(nveg+8+k) = itex(j,k)*1.0
                end do
              else
              end if
            end if
            write (48) stores
!add
            yobs(lrec) = rlat
            xobs(lrec) = rlon
            ht(lrec) = stores(3)
            htsd(lrec) = stores(4)
            ht2(lrec) = stores(4)**2 + stores(3)**2
            if ( grdlnmn<=-180.0 .and. xobs(lrec)>0.0 ) xobs(lrec)      &
               & = xobs(lrec) - 360.
            if ( xobs(lrec)<grdlnmn ) grdlnmn = xobs(lrec)
            if ( yobs(lrec)<grdltmn ) grdltmn = yobs(lrec)
!add_
          end if
        end do
      end do
 
      istat = nf90_close(ncid_lev)
      if ( istat /= nf90_noerr ) then
        write (6, *) 'Error close'
        write (6, *) nf90_strerror(istat)
        stop
      end if
      if ( lsmtyp=='BATS' ) then
        istat = nf90_close(ncid_cat)
        if ( istat /= nf90_noerr ) then
          write (6, *) 'Error close'
          write (6, *) nf90_strerror(istat)
          stop
        end if
      else if ( lsmtyp=='USGS' ) then
        istat = nf90_close(ncid_usgs)
        if ( istat /= nf90_noerr ) then
          write (6, *) 'Error close'
          write (6, *) nf90_strerror(istat)
          stop
        end if
        istat = nf90_close(ncid_sandclay)
        if ( istat /= nf90_noerr ) then
          write (6, *) 'Error close'
          write (6, *) nf90_strerror(istat)
          stop
        end if
      else
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
      jhmax = (xmaxlon-xminlon)/xnc
      if ( ihmax>iter .or. jhmax>jter ) print 99012 , iter , ihmax ,    &
         & jter , jhmax
      if ( ihmax*jhmax>iblk ) print 99013 , ihmax*jhmax , iblk
 
      nobs = lrec
 
      rewind (48)

      deallocate(icat)
      deallocate(iusgs)
      deallocate(iclay)
      deallocate(isand)
      deallocate(itex)
      deallocate(itopo)
      deallocate(ihtsd)

99001 format ('../DATA/SURFACE/GTOPO30_',i1,'min.nc')
99002 format ('../DATA/SURFACE/GTOPO30_',i2,'min.nc')
99003 format ('../DATA/SURFACE/GLCC_BATS_',i1,'min.nc')
99004 format ('../DATA/SURFACE/GLCC_BATS_',i2,'min.nc')
99005 format ('../DATA/SURFACE/GLCC_USGS_',i1,'min.nc')
99006 format ('../DATA/SURFACE/SAND_CLAY_',i1,'min.nc')
99007 format ('../DATA/SURFACE/GLCC_USGS_',i2,'min.nc')
99008 format ('../DATA/SURFACE/SAND_CLAY_',i2,'min.nc')
99009 format ('../DATA/SURFACE/SOILCAT_',i1,'min.nc')
99010 format ('../DATA/SURFACE/SOILCAT_',i2,'min.nc')
99011 format (1x,i10,' terrain heights read from land use volume')
99012 format (1x,'***array dimension error***',/,'     iter = ',i5,     &
             &' must be greater than ',i5,10x,'jter = ',i5,             &
             &' must be greater than ',i5)
99013 format (1x,'***array dimension error***',/,'  ihmax*jhmax = ',i8, &
             &' must be greater than ',i8,10x,'iblk = ',i8)
      end subroutine rdldtr_nc
!
      end module mod_rdldtr
