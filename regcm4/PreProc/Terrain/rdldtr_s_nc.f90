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

      subroutine rdldtr_s_nc
      use netcdf
      use mod_regcm_param , only : ibyte , lsmtyp , aertyp , nveg
      use mod_preproc_param , only : ntypec , ntypec_s , ntex
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
      character(34) :: filcat , filsandclay , filusgs
      character(32) :: filelev , filtext
      integer(1) , dimension(nxmax,nveg) :: icat , iusgs
      integer(1) , dimension(nxmax,2) :: iclay , isand
      integer , dimension(2) :: icount , istart
      integer , dimension(3) :: icount3 , istart3
      integer :: ihmax , ilat1 , ilat2 , irec , istat , j , jhmax , k , &
               & lrec , ncid_cat , ncid_lev , ncid_sandclay , ncid_tex ,&
               &  ncid_usgs
      integer(2) , dimension(nxmax) :: ihtsd , itopo
      integer(1) , dimension(nxmax,ntex) :: itex
      logical :: there
!
      if ( ntypec_s<10 ) then
        write (filelev,99001) ntypec_s
      else if ( ntypec_s<100 ) then
        write (filelev,99002) ntypec_s
      else
        write (*,*) 'For terrain, ntypec_s is not set correctly ' ,     &
                  & ntypec_s
        stop 'subroutine RDLDTR_s_nc'
      end if
      inquire (file=filelev,exist=there)
      if ( .not.there ) then
        print * , 'ERROR OPENING ' , filelev ,                          &
             &' FILE:  FILE DOES NOT EXIST'
        stop '4810 IN SUBROUTINE RDLDTR_s_nc'
      end if
      istat = nf90_open(filelev,nf90_nowrite,ncid_lev)
 
      if ( lsmtyp=='BATS' ) then
        if ( ntypec_s<10 ) then
          write (filcat,99003) ntypec_s
        else if ( ntypec_s<100 ) then
          write (filcat,99004) ntypec_s
        else
          write (*,*) 'For landuse, ntypec_s is not set correctly ' ,   &
                    & ntypec_s
          stop 'subroutine RDLDTR_s_nc'
        end if
        inquire (file=filcat,exist=there)
        if ( .not.there ) then
          print * , 'ERROR OPENING ' , filcat ,                         &
               &' FILE:  FILE DOES NOT EXIST'
          stop '4820 IN SUBROUTINE RDLDTR_s_nc'
        end if
        istat = nf90_open(filcat,nf90_nowrite,ncid_cat)
      else if ( lsmtyp=='USGS' ) then
        if ( ntypec_s<10 ) then
          write (filusgs,99005) ntypec_s
          write (filsandclay,99006) ntypec_s
        else if ( ntypec_s<100 ) then
          write (filusgs,99007) ntypec_s
          write (filsandclay,99008) ntypec_s
        else
          write (*,*) 'Error input of ntypec_s: ' , ntypec_s
        end if
        inquire (file=filusgs,exist=there)
        if ( .not.there ) then
          print * , 'ERROR OPENING ' , filusgs ,                        &
               &' FILE: FILE DOES NOT EXIST'
          stop '4820 IN SUBROUTINE RDLDTR_s_nc'
        end if
        inquire (file=filsandclay,exist=there)
        if ( .not.there ) then
          print * , 'ERROR OPENING ' , filsandclay ,                    &
               &' FILE: FILE DOES NOT EXIST'
          stop '4820 IN SUBROUTINE RDLDTR_s_nc'
        end if
        istat = nf90_open(filusgs,nf90_nowrite,ncid_usgs)
        istat = nf90_open(filsandclay,nf90_nowrite,ncid_sandclay)
      else
      end if
      if ( aertyp(7:7)=='1' ) then
        if ( ntypec_s<10 ) then
          write (filtext,99009) ntypec_s
        else if ( ntypec_s<100 ) then
          write (filtext,99010) ntypec_s
        else
          write (*,*) 'Error input of ntypec_s: ' , ntypec_s
        end if
        inquire (file=filtext,exist=there)
        if ( .not.there ) then
          print * , 'ERROR OPENING ' , filtext ,                        &
               &' FILE: FILE DOES NOT EXIST'
          stop '4830 IN SUBROUTINE RDLDTR_s_nc'
        end if
        istat = nf90_open(filtext,nf90_nowrite,ncid_tex)
      end if
 
      lrec = 0
      rewind (48)
      center = xnc
      ilat1 = nint((90.-xmaxlat)/xnc)
      ilat1 = max(min(nymax,ilat1),1)
      ilat2 = nint((90.-xminlat)/xnc) + 1
      ilat2 = max(min(nymax,ilat2),1)
      istart(1) = 1
      icount(1) = nxmax
      istart3(1) = 1
      icount3(1) = nxmax
      do irec = ilat1 , ilat2
        rlat = 90. - center*irec + center/2.
        istart(2) = nymax - irec + 1
        icount(2) = 1
        istat = nf90_get_var(ncid_lev,3,itopo,istart,icount)
        istat = nf90_get_var(ncid_lev,4,ihtsd,istart,icount)
 
        istart3(2) = nymax - irec + 1
        icount3(2) = 1
        if ( lsmtyp=='BATS' ) then
          istart3(3) = 1
          icount3(3) = nveg
          istat = nf90_get_var(ncid_cat,4,icat,istart3,icount3)
        else if ( lsmtyp=='USGS' ) then
          istart3(3) = 1
          icount3(3) = nveg
          istat = nf90_get_var(ncid_usgs,4,iusgs,istart3,icount3)
          istart3(3) = 1
          icount3(3) = 2
          istat = nf90_get_var(ncid_sandclay,4,isand,istart3,icount3)
          istat = nf90_get_var(ncid_sandclay,5,iclay,istart3,icount3)
        else
        end if
        if ( aertyp(7:7)=='1' ) then
          istart3(3) = 1
          icount3(3) = ntex
          istat = nf90_get_var(ncid_tex,4,itex,istart3,icount3)
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
      if ( lsmtyp=='BATS' ) then
        istat = nf90_close(ncid_cat)
      else if ( lsmtyp=='USGS' ) then
        istat = nf90_close(ncid_usgs)
        istat = nf90_close(ncid_sandclay)
      else
      end if
      if ( aertyp(7:7)=='1' ) istat = nf90_close(ncid_tex)
!
      print 99011 , lrec
      ihmax = (xmaxlat-xminlat)/xnc
      jhmax = (xmaxlon-xminlon)/xnc
      if ( ihmax>iter .or. jhmax>jter ) print 99012 , iter , ihmax ,    &
         & jter , jhmax
      if ( ihmax*jhmax>iblk ) print 99013 , ihmax*jhmax , iblk
 
      nobs = lrec
 
      rewind (48)
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
99013 format (1x,'***array dimension error***',/,'  ihmax*jhmax = ',i5, &
             &' must be greater than ',i5,10x,'iblk = ',i5)
      end subroutine rdldtr_s_nc
 
