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

      subroutine rdldtr_s
      use mod_param
      use mod_a
      use mod_aa
      use mod_block
      use mod_const
      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: nxmax = 21600/ntypec_s , nymax = nxmax/2
!
! Local variables
!
      real(4) :: center , rlat , rlon
      character(2) :: char_2
      character(1) , dimension(nxmax,nveg) :: ch_cat
      character(2) , dimension(nxmax) :: ch_htsd , ch_topo
      character(1) , dimension(nxmax,ntex) :: ch_tex
      character(36) :: filcat , filelev
      character(7) :: filclay , filsand
      character(10) :: filtext
      character(11) :: filusgs
      integer(2) , dimension(nxmax,2) :: iclay , isand
      integer :: ihmax , ilat1 , ilat2 , irec , irect , j , jhmax , k , &
               & lrec
      integer(2) , dimension(nxmax,nveg) :: iusgs
      logical :: there
!
      if ( ntypec_s<10 ) then
        write (filelev,99001) ntypec_s
      else if ( ntypec_s<100 ) then
        write (filelev,99002) ntypec_s
      else
        write (*,*) 'For terrain, ntypec_s is not set correctly ' ,     &
                  & ntypec_s
        stop 'subroutine RDLDTR_s'
      end if
      inquire (file=filelev,exist=there)
      if ( .not.there ) then
        print * , 'ERROR OPENING ' , filelev ,                          &
             &' FILE:  FILE DOES NOT EXIST'
        stop '4810 IN SUBROUTINE RDLDTR_s'
      end if
      open (46,file=filelev,form='unformatted',recl=nxmax*ibyte/2,      &
          & access='direct')
 
      if ( lsmtyp=='BATS' ) then
        if ( ntypec_s<10 ) then
          write (filcat,99003) ntypec_s
        else if ( ntypec_s<100 ) then
          write (filcat,99004) ntypec_s
        else
          write (*,*) 'For landuse, ntypec_s is not set correctly ' ,   &
                    & ntypec_s
          stop 'subroutine RDLDTR_s'
        end if
        inquire (file=filcat,exist=there)
        if ( .not.there ) then
          print * , 'ERROR OPENING ' , filcat ,                         &
               &' FILE:  FILE DOES NOT EXIST'
          stop '4820 IN SUBROUTINE RDLDTR_s'
        end if
        open (47,file=filcat,form='unformatted',recl=nxmax*nveg*ibyte/4,&
            & access='direct')
      else if ( lsmtyp=='USGS' ) then
        if ( ntypec_s<10 ) then
          write (filusgs,99005) ntypec_s
          write (filsand,99006) ntypec_s
          write (filclay,99007) ntypec_s
        else if ( ntypec_s<100 ) then
          write (filusgs,99008) ntypec_s
          write (filsand,99009) ntypec_s
          write (filclay,99010) ntypec_s
        else
          write (*,*) 'Error input of ntypec_s: ' , ntypec_s
        end if
        if ( ntypec_s==2 ) then
          inquire (file='../DATA/SURFACE/'//filusgs//'A',exist=there)
          if ( .not.there ) then
            print * , 'ERROR OPENING ' , '../DATA/SURFACE/'//filusgs ,  &
                 &' FILE: FILE DOES NOT EXIST'
            stop '4820 IN SUBROUTINE RDLDTR_s'
          end if
          inquire (file='../DATA/SURFACE/'//filusgs//'B',exist=there)
          if ( .not.there ) then
            print * , 'ERROR OPENING ' , '../DATA/SURFACE/'//filusgs ,  &
                 &' FILE: FILE DOES NOT EXIST'
            stop '4820 IN SUBROUTINE RDLDTR_s'
          end if
        else
          inquire (file='../DATA/SURFACE/'//filusgs,exist=there)
          if ( .not.there ) then
            print * , 'ERROR OPENING ' , '../DATA/SURFACE/'//filusgs ,  &
                 &' FILE: FILE DOES NOT EXIST'
            stop '4820 IN SUBROUTINE RDLDTR_s'
          end if
        end if
        inquire (file='../DATA/SURFACE/'//filsand,exist=there)
        if ( .not.there ) then
          print * , 'ERROR OPENING ' , '../DATA/SURFACE/'//filsand ,    &
               &' FILE: FILE DOES NOT EXIST'
          stop '4820 IN SUBROUTINE RDLDTR_s'
        end if
        inquire (file='../DATA/SURFACE/'//filclay,exist=there)
        if ( .not.there ) then
          print * , 'ERROR OPENING ' , '../DATA/SURFACE/'//filclay ,    &
               &' FILE: FILE DOES NOT EXIST'
          stop '4820 IN SUBROUTINE RDLDTR_s'
        end if
        if ( ntypec_s==2 ) then
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
        if ( ntypec_s<10 ) then
          write (filtext,99011) ntypec_s
        else if ( ntypec_s<100 ) then
          write (filtext,99012) ntypec_s
        else
          write (*,*) 'Error input of ntypec_s: ' , ntypec_s
        end if
        inquire (file='../DATA/SURFACE/'//filtext,exist=there)
        if ( .not.there ) then
          print * , 'ERROR OPENING ' , '../DATA/SURFACE/'//filtext ,    &
               &' FILE: FILE DOES NOT EXIST'
          stop '4830 IN SUBROUTINE RDLDTR_s'
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
        read (46,rec=irect) (ch_htsd(j),j=1,nxmax)
!       print*,'   ELEVATION STD DEV READ IN', ilat1,ilat2,irec
        if ( lsmtyp=='BATS' ) then
          read (47,rec=irect) ((ch_cat(j,k),k=1,nveg),j=1,nxmax)
        else if ( lsmtyp=='USGS' ) then
          if ( ntypec_s==2 ) then
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
               & then    ! OCEAN/UNDEFINED
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
99015 format (1x,'***array dimension error***',/,'  ihmax*jhmax = ',i5, &
             &' must be greater than ',i5,10x,'iblk = ',i5)
      end subroutine rdldtr_s
