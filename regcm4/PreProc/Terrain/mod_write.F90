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

      module mod_write

      contains

      subroutine setup(nunit,unctl,iy,jx,ntypec,iproj,ds,clat,clon,     &
                     & igrads,ibyte,filout,filctl)
      use mod_block
      implicit none
!
! Dummy arguments
!
      real(4) :: clat , clon , ds
      character(256) :: filctl , filout
      integer :: ibyte , igrads , iy , jx , ntypec , nunit , unctl
      character(6) :: iproj
      intent (in) clat , clon , ds , ibyte , igrads , iproj , iy , jx , &
                & ntypec , nunit , unctl
      intent(inout) :: filctl , filout
!
      rin = 1.5          ! 1.5 rad of influence-coarse mesh
 
      write (6,*) ' '
      write (6,*) 'Doing Domain Setup with following parameters'
      write (6,*) ' '
      write (6,*) 'ntypec = ' , ntypec
      write (6,*) 'iy     = ' , iy
      write (6,*) 'jx     = ' , jx
      write (6,*) 'ds     = ' , ds
      write (6,*) 'clat   = ' , clat
      write (6,*) 'clon   = ' , clon
      write (6,*) 'rin    = ' , rin
      write (6,*) 'iproj  = ' , iproj
      write (6,*) ' '
!
      call fexist(filout)
      open (nunit,file=filout,status='replace',form='unformatted',      &
          & access='direct',recl=iy*jx*ibyte)
      if ( igrads==1 ) then
        call fexist(filctl)
        open (unctl,file=filctl,status='unknown')
      end if
!
      dsinm = ds*1000.
!
      nnc = nint(60./float(ntypec))
      xnc = float(ntypec)/60.
      print * , '***** Terrain resolution (min): ' , xnc*60.
!
      end subroutine setup

      subroutine fexist(filnam)
      implicit none
!
! Dummy arguments
!
      character(256) :: filnam
      intent (inout) filnam
!
! Local variables
!
      logical :: there
      character(1) :: yesno

 100  continue
      inquire (file=filnam,exist=there)
      if ( there ) then
 150    continue
        print * , ' '
        print * , ' '
        print * , '**************************************************'
        print * , 'FILE ALREADY EXISTS:  ' , trim(filnam)
        print * , 'Do you want to overwrite the existing file? [y/n/q]'
        read (*,*) yesno
        if ( yesno=='y' ) then
          return
        else if ( yesno=='n' ) then
          print * , 'ENTER NEW FILE NAME'
          read (*,'(a)') filnam
          print * , 'USING ', trim(filnam)
          goto 100
        else if ( yesno=='q' ) then
          stop 999
        else
          go to 150
        end if
      end if
 
      end subroutine fexist
!
!
!
      subroutine output(nunitc,iunc,iy,jx,dsinm,clat,clon,plat,plon,    &
                      & iproj,htgrid,htsdgrid,lndout,xlat,xlon,dlat,    &
                      & dlon,xmap,dattyp,dmap,coriol,snowam,igrads,     &
                      & ibigend,kz,sigma,mask,ptop,nsg,truelatl,        &
                      & truelath,grdfac,filout,aertyp,texout,frac_tex,  &
                      & ntex,lcoarse)

      implicit none
!
! Dummy arguments
!
      character(7) :: aertyp
      real(4) :: clat , clon , dsinm , grdfac , plat , plon , ptop ,    &
               & truelath , truelatl
      character(5) :: dattyp
      character(*) :: filout
      integer :: ibigend , igrads , iy , jx , kz , nsg , ntex ,         &
               & nunitc , iunc
      character(6) :: iproj
      real(4) , dimension(iy,jx) :: coriol , dlat , dlon , dmap ,       &
                                  & htgrid , htsdgrid , lndout , mask , &
                                  & snowam , texout , xlat , xlon , xmap
      real(4) , dimension(iy,jx,ntex) :: frac_tex
      real(4) , dimension(kz+1) :: sigma
      logical :: lcoarse
      intent (in) aertyp , clat , clon , coriol , dattyp , dlat , dlon ,&
                & dmap , dsinm , filout , frac_tex , grdfac ,           &
                & htgrid , htsdgrid , ibigend , igrads , iproj , iy ,   &
                & jx , kz , lndout , nsg , ntex , nunitc , plat , plon ,&
                & ptop , snowam , texout , truelath , truelatl , xlat , &
                & xlon , xmap , lcoarse , sigma , iunc
      intent (inout) mask
!
! Local variables
!
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
               & centerj , lat0 , lat1 , lon0 , lon1 , rlatinc , rloninc
      integer :: i , j , k , nx , ny
!
      alatmin = 999999.
      alatmax = -999999.
      alonmin = 999999.
      alonmax = -999999.
      nx = 0
      ny = 0

      do i = 1 , iy
        do j = 1 , jx
          if ( ((lndout(i,j)>20. .or. lndout(i,j)<0.)) ) then
            print * , i , j , lndout(i,j)
            stop 999
          end if
        end do
      end do

      if ( .not. lcoarse .or.                                           &
         & ( dattyp/='FVGCM' .and. dattyp/='NRP2W' .and.  &
         &   dattyp/='GFS11' .and. dattyp/='EH5OM') ) then
        write (nunitc,rec=1) iy , jx , kz , dsinm , clat , clon , plat ,&
                           & plon , grdfac , iproj , (sigma(k),k=1,kz+1)&
                           & , ptop , igrads , ibigend , truelatl ,     &
                           & truelath
      else
        write (*,*) 'please input lon0,lon1,lat0,lat1'
        write (*,*) 'Note: lon0 < lon1, and lat0 < lat1'
        read (*,*) lon0 , lon1 , lat0 , lat1
        write (nunitc,rec=1) iy , jx , kz , dsinm , clat , clon , plat ,&
                           & plon , grdfac , iproj , (sigma(k),k=1,kz+1)&
                           & , ptop , igrads , ibigend , truelatl ,     &
                           & truelath , lon0 , lon1 , lat0 , lat1
      end if

      write (nunitc,rec=2) ((htgrid(i,j),j=1,jx),i=1,iy)
      write (nunitc,rec=3) ((htsdgrid(i,j),j=1,jx),i=1,iy)
      write (nunitc,rec=4) ((lndout(i,j),j=1,jx),i=1,iy)
      write (nunitc,rec=5) ((xlat(i,j),j=1,jx),i=1,iy)
      write (nunitc,rec=6) ((xlon(i,j),j=1,jx),i=1,iy)
      write (nunitc,rec=7) ((dlat(i,j),j=1,jx),i=1,iy)
      write (nunitc,rec=8) ((dlon(i,j),j=1,jx),i=1,iy)
      write (nunitc,rec=9) ((xmap(i,j),j=1,jx),i=1,iy)
      write (nunitc,rec=10) ((dmap(i,j),j=1,jx),i=1,iy)
      write (nunitc,rec=11) ((coriol(i,j),j=1,jx),i=1,iy)
      write (nunitc,rec=12) ((snowam(i,j),j=1,jx),i=1,iy)
      do i = 1 , iy
        do j = 1 , jx
          if ( lndout(i,j)>13.5 .and. lndout(i,j)<15.5 ) then
            mask(i,j) = 0.0
          else
            mask(i,j) = 2.0
          end if
        end do
      end do
      write (nunitc,rec=13) ((mask(i,j),j=1,jx),i=1,iy)
      if ( aertyp(7:7)=='1' ) then
        write (nunitc,rec=14) ((texout(i,j),j=1,jx),i=1,iy)
        do k = 1 , ntex
          write (nunitc,rec=14+k) ((frac_tex(i,j,k),j=1,jx),i=1,iy)
        end do
      end if
      close (nunitc)
 
      if ( igrads==1 ) then
        if ( lcoarse ) then
          write (iunc,99004) trim(filout), '.INFO'
        else if ( nsg>1 ) then
          write (iunc,99005) trim(filout), nsg , '.INFO'
        end if
        write (iunc,99007)
        if ( ibigend==1 ) then
          write (iunc,99008)
        else
          write (iunc,99009)
        end if
        write (iunc,99010)
        if ( iproj=='LAMCON' .or. iproj=='ROTMER' ) then
          do j = 1 , jx
            if ( xlat(1,j)<alatmin ) alatmin = xlat(1,j)
            if ( xlat(iy,j)>alatmax ) alatmax = xlat(iy,j)
          end do
          do i = 1 , iy
            do j = 1 , jx
              if ( clon>=0.0 ) then
                if ( xlon(i,j)>=0.0 ) then
                  alonmin = amin1(alonmin,xlon(i,j))
                  alonmax = amax1(alonmax,xlon(i,j))
                else if ( abs(clon-xlon(i,j))<abs(clon-(xlon(i,j)+360.))&
                        & ) then
                  alonmin = amin1(alonmin,xlon(i,j))
                  alonmax = amax1(alonmax,xlon(i,j))
                else
                  alonmin = amin1(alonmin,xlon(i,j)+360.)
                  alonmax = amax1(alonmax,xlon(i,j)+360.)
                end if
              else if ( xlon(i,j)<0.0 ) then
                alonmin = amin1(alonmin,xlon(i,j))
                alonmax = amax1(alonmax,xlon(i,j))
              else if ( abs(clon-xlon(i,j))<abs(clon-(xlon(i,j)-360.)) )&
                      & then
                alonmin = amin1(alonmin,xlon(i,j))
                alonmax = amax1(alonmax,xlon(i,j))
              else
                alonmin = amin1(alonmin,xlon(i,j)-360.)
                alonmax = amax1(alonmax,xlon(i,j)-360.)
              end if
            end do
          end do
          rlatinc = dsinm*0.001/111./2.
          rloninc = dsinm*0.001/111./2.
          ny = 2 + nint(abs(alatmax-alatmin)/rlatinc)
          nx = 1 + nint(abs((alonmax-alonmin)/rloninc))
 
          centerj = jx/2.
          centeri = iy/2.
        end if
        if ( iproj=='LAMCON' ) then        ! Lambert projection
          write (iunc,99011) jx , iy , clat , clon , centerj , centeri ,&
                         & truelatl , truelath , clon , dsinm , dsinm
          write (iunc,99012) nx + 2 , alonmin - rloninc , rloninc
          write (iunc,99013) ny + 2 , alatmin - rlatinc , rlatinc
        else if ( iproj=='POLSTR' ) then   !
        else if ( iproj=='NORMER' ) then
          write (iunc,99014) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
          write (iunc,99015) iy
          write (iunc,99016) (xlat(i,1),i=1,iy)
        else if ( iproj=='ROTMER' ) then
          write (*,*) 'Note that rotated Mercartor (ROTMER)' ,          &
                     &' projections are not supported by GrADS.'
          write (*,*) '  Although not exact, the eta.u projection' ,    &
                     &' in GrADS is somewhat similar.'
          write (*,*) ' FERRET, however, does support this projection.'
          write (iunc,99017) jx , iy , plon , plat , dsinm/111000. ,    &
                         & dsinm/111000.*.95238
          write (iunc,99012) nx + 2 , alonmin - rloninc , rloninc
          write (iunc,99013) ny + 2 , alatmin - rlatinc , rlatinc
        else
          write (*,*) 'Are you sure your map projection is right ?'
          stop
        end if
        write (iunc,99018) 1 , 1000.
        write (iunc,99019) 1
        if ( aertyp(7:7)=='1' ) then
          write (iunc,99020) 13 + ntex + 1
        else
          write (iunc,99020) 13
        end if
        write (iunc,99021) 'head    ' , 'header information         '
        write (iunc,99021) 'ht      ' , 'surface elevation          '
        write (iunc,99021) 'htsd    ' , 'surface elevation std. dev.'
        write (iunc,99021) 'landuse ' , 'surface landuse type       '
        write (iunc,99021) 'xlat    ' , 'latitude  of cross points  '
        write (iunc,99021) 'xlon    ' , 'longitude of cross points  '
        write (iunc,99021) 'dlat    ' , 'latitude  of dot points    '
        write (iunc,99021) 'dlon    ' , 'longitude of dot points    '
        write (iunc,99021) 'xmap    ' , 'map factors of cross points'
        write (iunc,99021) 'dmap    ' , 'map factors of dot points  '
        write (iunc,99021) 'coriol  ' , 'coriol force               '
        write (iunc,99021) 'snowam  ' , 'initial snow amount        '
        write (iunc,99021) 'mask    ' , 'land/sea mask              '
        if ( aertyp(7:7)=='1' ) then
          write (iunc,99021) 'texture ' , 'soil texture               '
          write (iunc,99021) 'text01  ' , 'Sand              frac.    '
          write (iunc,99021) 'text02  ' , 'Loamy Sand        frac.    '
          write (iunc,99021) 'text03  ' , 'Sandy Loam        frac.    '
          write (iunc,99021) 'text04  ' , 'Silt Loam         frac.    '
          write (iunc,99021) 'text05  ' , 'Silt              frac.    '
          write (iunc,99021) 'text06  ' , 'Loam              frac.    '
          write (iunc,99021) 'text07  ' , 'Sandy Clay Loam   frac.    '
          write (iunc,99021) 'text08  ' , 'Silty Clay Loam   frac.    '
          write (iunc,99021) 'text09  ' , 'Clay Loam         frac.    '
          write (iunc,99021) 'text10  ' , 'Sandy Clay        frac.    '
          write (iunc,99021) 'text11  ' , 'Silty Clay        frac.    '
          write (iunc,99021) 'text12  ' , 'Clay              frac.    '
          write (iunc,99021) 'text13  ' , 'OM                frac.    '
          write (iunc,99021) 'text14  ' , 'Water             frac.    '
          write (iunc,99021) 'text15  ' , 'Bedrock           frac.    '
          write (iunc,99021) 'text16  ' , 'Other             frac.    '
          write (iunc,99021) 'text17  ' , 'No data           frac.    '
        end if
        write (iunc,99022)
        close (iunc)
      end if
99004 format ('dset ^',a,a)
99005 format ('dset ^',a,i0.3,a)
99007 format ('title RegCM domain information')
99008 format ('options big_endian')
99009 format ('options little_endian')
99010 format ('undef -9999.')
99011 format ('pdef ',i4,1x,i4,1x,'lcc',7(1x,f7.2),1x,2(f7.0,1x))
99012 format ('xdef ',i4,' linear ',f7.2,1x,f7.4)
99013 format ('ydef ',i4,' linear ',f7.2,1x,f7.4)
99014 format ('xdef ',i3,' linear ',f9.4,' ',f9.4)
99015 format ('ydef ',i3,' levels')
99016 format (10F7.2)
99017 format ('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
99018 format ('zdef ',i1,' levels ',f7.2)
99019 format ('tdef ',i1,' linear 00z01Jan2001 1mo')
99020 format ('vars ',i2)
99021 format (a8,'0 99 ',a26)
99022 format ('endvars')
!
      end subroutine output

      subroutine write_domain(lsub)
        use netcdf
        use mod_dynparam
        use mod_block
        use mod_maps
        implicit none
        logical , intent (in) :: lsub

        integer :: istatus , i , j
        integer :: incout
        integer , dimension(3) :: idims
        integer , dimension(3) :: istart
        integer , dimension(3) :: icount
        integer , dimension(12) :: ivar
        integer , dimension(2) :: itvar
        integer , dimension(2) :: ivdim
        integer , dimension(8) :: tvals
        character(256) :: fname , history
        character(3) :: cnsg
        real(4) , dimension(2) :: trlat
        real(4) , allocatable , dimension(:) :: yiy
        real(4) , allocatable , dimension(:) :: xjx

        trlat(1) = truelatl
        trlat(2) = truelath

        if (lsub) then
          allocate(yiy(iysg))
          allocate(xjx(jxsg))
          yiy(1) = -(dble(iysg-1)/2.0) * ds
          xjx(1) = -(dble(jxsg-1)/2.0) * ds
          do i = 2 , iysg
            yiy(i) = yiy(i-1)+ds
          end do
          do j = 2 , jxsg
            xjx(j) = xjx(j-1)+ds
          end do
        else
          allocate(yiy(iy))
          allocate(xjx(jx))
          yiy(1) = -(dble(iy-1)/2.0) * ds
          xjx(1) = -(dble(jx-1)/2.0) * ds
          do i = 2 , iy
            yiy(i) = yiy(i-1)+ds
          end do
          do j = 2 , jx
            xjx(j) = xjx(j-1)+ds
          end do
        end if

        if (lsub) then
          write (cnsg, '(i0.3)') nsg
          fname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN'//      &
                    & cnsg//'.nc'
        else
          fname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_create(fname, nf90_clobber.or.nf90_hdf5, incout)
#else
        istatus = nf90_create(fname, nf90_clobber, incout)
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error creating NetCDF output ', trim(fname)
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        istatus = nf90_put_att(incout, nf90_global, 'title',            &
                 & 'ICTP Regional Climatic model V4 domain')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global title'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global, 'institution',      &
                 & 'ICTP')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global institution'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global, 'source',           &
                 & 'RegCM Model simulation Terrain output')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global source'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global, 'Conventions',      &
                 & 'CF-1.4')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global Conventions'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        call date_and_time(values=tvals)
        write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)')   &
             tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,         &
             tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,               &
             ' : Created by RegCM terrain'
        istatus = nf90_put_att(incout, nf90_global, 'history', history)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global history'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global, 'references',       &
                 & 'http://eforge.escience-lab.org/gf/project/regcm')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global references'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global, 'experiment',       &
                 & domname)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global experiment'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global, 'projection', iproj)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global projection'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global,                     &
                 &   'grid_size_in_meters', ds*1000.0)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global gridsize'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global,                     &
                 &   'latitude_of_projection_origin', clat)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global clat'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global,                     &
                 &   'longitude_of_projection_origin', clon)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global clon'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (iproj == 'ROTMER') then
          istatus = nf90_put_att(incout, nf90_global,                   &
                   &   'latitude_of_projection_pole', plat)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global plat'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, nf90_global,                   &
                   &   'longitude_of_projection_pole', plon)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global plon'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        else if (iproj == 'LAMCON') then
          istatus = nf90_put_att(incout, nf90_global,                   &
                   &   'standard_parallel', trlat)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global truelat'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if
        if (ifanal) then
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'data_interpolation', 'Cressman type objective analysis')
        else
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'data_interpolation', 'Overlapping parabolic 16 points')
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global data_interpolation'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (smthbdy) then
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'boundary_smoothing', 'Yes')
        else
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'boundary_smoothing', 'No')
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global boundary_smoothing'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lakadj) then
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'great_lakes_adjustment', 'Yes')
        else
          istatus = nf90_put_att(incout, nf90_global,                   &
           &   'great_lakes_adjustment', 'No')
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global great_lakes_adjustment'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, nf90_global,                     &
                 &   'minimum_h2o_pct_for_water', h2opct)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error adding global minimum_h2o_pct_for_water'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        if (lsub) then
          istatus = nf90_put_att(incout, nf90_global,                   &
                   &   'input_dataset_resolution_in_minutes', ntypec_s)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global ntypec_s'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          if (fudge_lnd_s) then
            istatus = nf90_put_att(incout, nf90_global,                 &
             &     'landuse_fudging', 'Yes')
          else
            istatus = nf90_put_att(incout, nf90_global,                 &
             &     'landuse_fudging', 'No')
          end if
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global landuse_fudging_s'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          if ( aertyp(7:7)=='1' ) then
            if (fudge_tex_s) then
              istatus = nf90_put_att(incout, nf90_global,               &
               &     'texture_fudging', 'Yes')
            else
              istatus = nf90_put_att(incout, nf90_global,               &
               &     'texture_fudging', 'No')
            end if
            if (istatus /= nf90_noerr) then
              write (6,*) 'Error adding global texture_fudging_s'
              write (6,*) nf90_strerror(istatus)
              stop
            end if
          end if
        else
          istatus = nf90_put_att(incout, nf90_global,                   &
                   &   'input_dataset_resolution_in_minutes', ntypec)
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global ntypec_s'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          if (fudge_lnd) then
            istatus = nf90_put_att(incout, nf90_global,                 &
             &     'landuse_fudging', 'Yes')
          else
            istatus = nf90_put_att(incout, nf90_global,                 &
             &     'landuse_fudging', 'No')
          end if
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error adding global landuse_fudging'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          if ( aertyp(7:7)=='1' ) then
            if (fudge_tex) then
              istatus = nf90_put_att(incout, nf90_global,               &
             &       'texture_fudging', 'Yes')
            else
              istatus = nf90_put_att(incout, nf90_global,               &
             &       'texture_fudging', 'No')
            end if
            if (istatus /= nf90_noerr) then
              write (6,*) 'Error adding global texture_fudging'
              write (6,*) nf90_strerror(istatus)
              stop
            end if
          end if
        end if

        if (lsub) then
          istatus = nf90_def_dim(incout, 'IY', iysg, idims(2))
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error creating dimension IY'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_def_dim(incout, 'JX', jxsg, idims(1))
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error creating dimension JX'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        else
          istatus = nf90_def_dim(incout, 'IY', iy, idims(2))
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error creating dimension IY'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_def_dim(incout, 'JX', jx, idims(1))
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error creating dimension JX'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if
        if ( aertyp(7:7)=='1' ) then
          istatus = nf90_def_dim(incout, 'NTEX', ntex, idims(3))
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error creating dimension NVEG'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'IY', nf90_float, idims(2),      &
                            &  ivdim(1), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'IY', nf90_float, idims(2),      &
                            &  ivdim(1))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable iy definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivdim(1), 'standard_name',       &
                            &  'projection_y_coordinate')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable iy standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivdim(1), 'long_name',           &
                            &  'y-coordinate in Cartesian system')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable iy long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivdim(1), 'units', 'km')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable iy units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'JX', nf90_float, idims(1),      &
                            &  ivdim(2), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'JX', nf90_float, idims(1),      &
                            &  ivdim(2))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable jx definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivdim(2), 'standard_name',       &
                            &  'projection_x_coordinate')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable jx standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivdim(2), 'long_name',           &
                            &  'x-coordinate in Cartesian system')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable jx long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivdim(2), 'units', 'km')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable jx units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        ! XLAT
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'xlat', nf90_float, idims(1:2),  &
                            &  ivar(1), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'xlat', nf90_float, idims(1:2),  &
                            &  ivar(1))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlat definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(1), 'standard_name',        &
                            &  'latitude')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlat standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(1), 'long_name',            &
                            &  'Latitude at cross points')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlat long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(1), 'units',                &
                            &  'degrees_north')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlat units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        ! XLAT

        ! XLON
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'xlon', nf90_float, idims(1:2),  &
                            &  ivar(2), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'xlon', nf90_float, idims(1:2),  &
                            &  ivar(2))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlon definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(2), 'standard_name',        &
                            &  'longitude')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlon standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(2), 'long_name',            &
                            &  'Longitude at cross points')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlon long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(2), 'units',                &
                            &  'degrees_east')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlon units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        ! XLON

        ! DLAT
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'dlat', nf90_float, idims(1:2),  &
                            &  ivar(3), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'dlat', nf90_float, idims(1:2),  &
                            &  ivar(3))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlat definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(3), 'standard_name',        &
                            &  'latitude')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlat standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(3), 'long_name',            &
                            &  'Latitude at dot points')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlat long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(3), 'units',                &
                            &  'degrees_north')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlat units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        ! DLAT

        ! DLON
#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'dlon', nf90_float, idims(1:2),  &
                            &  ivar(4), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'dlon', nf90_float, idims(1:2),  &
                            &  ivar(4))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlon definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(4), 'standard_name',        &
                            &  'longitude')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlon standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(4), 'long_name',            &
                            &  'Longitude at dot points')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlon long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(4), 'units',                &
                            &  'degrees_east')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlon units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        ! DLON

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'topo', nf90_float, idims(1:2),  &
                            &  ivar(5), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'topo', nf90_float, idims(1:2),  &
                            &  ivar(5))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ht definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(5), 'standard_name',        &
                            &  'surface_altitude')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ht standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(5), 'long_name',            &
                            &  'Domain surface elevation')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ht long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(5), 'units', 'm')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ht units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(5), 'coordinates',          &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ht coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'htsd', nf90_float, idims(1:2),  &
                            &  ivar(6), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'htsd', nf90_float, idims(1:2),  &
                            &  ivar(6))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable htsd definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(6), 'cell_method',          &
                            &  'area: standard_deviation')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable htsd cell_method attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(6), 'standard_name',        &
                            &  'surface_altitude')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable htsd standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(6), 'long_name',            &
                      &  'Domain surface elevation standard deviation')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable htsd long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(6), 'units', 'm')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable htsd units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(6), 'coordinates',          &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable htsd coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'landuse', nf90_float,idims(1:2),&
                            &  ivar(7), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'landuse', nf90_float,idims(1:2),&
                            &  ivar(7))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse def in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(7), 'legend',               &
                & '1  => Crop/mixed farming'//char(10)//                &
                & '2  => Short grass'//char(10)//                       &
                & '3  => Evergreen needleleaf tree'//char(10)//         &
                & '4  => Deciduous needleleaf tree'//char(10)//         &
                & '5  => Deciduous broadleaf tree'//char(10)//          &
                & '6  => Evergreen broadleaf tree'//char(10)//          &
                & '7  => Tall grass'//char(10)//                        &
                & '8  => Desert'//char(10)//                            &
                & '9  => Tundra'//char(10)//                            &
                & '10 => Irrigated Crop'//char(10)//                    &
                & '11 => Semi-desert'//char(10)//                       &
                & '12 => Ice cap/glacier'//char(10)//                   &
                & '13 => Bog or marsh'//char(10)//                      &
                & '14 => Inland water'//char(10)//                      &
                & '15 => Ocean'//char(10)//                             &
                & '16 => Evergreen shrub'//char(10)//                   &
                & '17 => Deciduous shrub'//char(10)//                   &
                & '18 => Mixed Woodland'//char(10)//                    &
                & '19 => Forest/Field mosaic'//char(10)//               &
                & '20 => Water and Land mixture')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse legend attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(7), 'standard_name',        &
                            &  'land_type')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(7), 'long_name',            &
                      &  'Landuse category as defined in BATS1E')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(7), 'units', '1')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(7), 'coordinates',          &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'xmap', nf90_float, idims(1:2),  &
                            &  ivar(8), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'xmap', nf90_float, idims(1:2),  &
                            &  ivar(8))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xmap definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(8), 'standard_name',        &
                            &  'map_factor')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xmap standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(8), 'long_name',            &
                            &  'Map factor in domain cross points')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xmap long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(8), 'units', '1')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xmap units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(8), 'coordinates',          &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xmap coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'dmap', nf90_float, idims(1:2),  &
                            &  ivar(9), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'dmap', nf90_float, idims(1:2),  &
                            &  ivar(9))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dmap definition in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(9), 'standard_name',        &
                            &  'map_factor')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dmap standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(9), 'long_name',            &
                            &  'Map factor in domain dot points')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dmap long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(9), 'units', '1')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dmap units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(9), 'coordinates',          &
                            &  'dlon dlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dmap coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'coriol', nf90_float, idims(1:2),&
                            &  ivar(10), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'coriol', nf90_float, idims(1:2),&
                            &  ivar(10))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable coriol def in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(10), 'standard_name',       &
                            &  'coriolis_parameter')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable coriol standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(10), 'long_name',           &
                            &  'Coriolis force parameter')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable coriol long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(10), 'units', 's-1')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable coriol units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(10), 'coordinates',         &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable coriol coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'snowam', nf90_float, idims(1:2),&
                            &  ivar(11), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'snowam', nf90_float, idims(1:2),&
                            &  ivar(11))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable snowam def in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(11), 'standard_name',       &
                            &  'snowfall_amount')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable snowam standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(11), 'long_name',           &
                            &  'Snow initial amount in mm')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable snowam long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(11), 'units', 'kg m-2')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable snowam units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(11), 'coordinates',         &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable snowam coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

#ifdef NETCDF4_HDF5
        istatus = nf90_def_var(incout, 'mask', nf90_float, idims(1:2),  &
                            &  ivar(12), deflate_level=9)
#else
        istatus = nf90_def_var(incout, 'mask', nf90_float, idims(1:2),  &
                            &  ivar(12))
#endif
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable mask def in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(12), 'standard_name',       &
                            &  'land_binary_mask')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable mask standard_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(12), 'long_name',           &
                            &  'Land Sea mask')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable mask long_name attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(12), 'units', '1')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable mask units attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_att(incout, ivar(12), 'coordinates',         &
                            &  'xlon xlat')
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable mask coordinates attribute'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        if ( aertyp(7:7)=='1' ) then
#ifdef NETCDF4_HDF5
          istatus = nf90_def_var(incout, 'texture', nf90_float,         &
                              & idims(1:2), itvar(1), deflate_level=9)
#else
          istatus = nf90_def_var(incout, 'texture', nf90_float,         &
                              & idims(1:2), itvar(1))
#endif
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture def in NetCDF output'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(1), 'legend',            &
                & '1  => Sand'//char(10)//                              &
                & '2  => Loamy Sand'//char(10)//                        &
                & '3  => Sandy Loam'//char(10)//                        &
                & '4  => Silt Loam'//char(10)//                         &
                & '5  => Silt'//char(10)//                              &
                & '6  => Loam'//char(10)//                              &
                & '7  => Sandy Clay Loam'//char(10)//                   &
                & '8  => Silty Clay Loam'//char(10)//                   &
                & '9  => Clay Loam'//char(10)//                         &
                & '10 => Sandy Clay'//char(10)//                        &
                & '11 => Silty Clay'//char(10)//                        &
                & '12 => Clay'//char(10)//                              &
                & '13 => OM'//char(10)//                                &
                & '14 => Water'//char(10)//                             &
                & '15 => Bedrock'//char(10)//                           &
                & '16 => Other'//char(10)//                             &
                & '17 => No data')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture legend attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(1), 'standard_name',     &
                              &  'soil_type')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture standard_name attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(1), 'long_name',         &
                              &  'Texture dominant category')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture long_name attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(1), 'units', '1')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture units attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(1), 'coordinates',       &
                              &  'xlon xlat')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture coordinates attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if

#ifdef NETCDF4_HDF5
          istatus = nf90_def_var(incout, 'texture_fraction', nf90_float,&
                              & idims, itvar(2), deflate_level=9)
#else
          istatus = nf90_def_var(incout, 'texture_fraction', nf90_float,&
                              & idims, itvar(2))
#endif
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture_fraction def in NetCDF'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(2), 'standard_name',     &
                              &  'soil_type_fraction')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture_fraction standard_name'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(2), 'long_name',         &
                              &  'Texture category fraction')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture_fraction long_name'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(2), 'units', '1')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture_fraction units'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istatus = nf90_put_att(incout, itvar(2), 'coordinates',       &
                              &  'xlon xlat')
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture coordinates attribute'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
        end if
!
!-----------------------------------------------------------------------
!
        istatus = nf90_enddef(incout)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error End Definitions NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
!
!-----------------------------------------------------------------------
!
        istatus = nf90_put_var(incout, ivdim(1), yiy)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable iy write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        istatus = nf90_put_var(incout, ivdim(2), xjx)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable jx write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        if (lsub) then
          istatus = nf90_put_var(incout, ivar(1), transpose(xlat_s))
        else
          istatus = nf90_put_var(incout, ivar(1), transpose(xlat))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlat write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(2), transpose(xlon_s))
        else
          istatus = nf90_put_var(incout, ivar(2), transpose(xlon))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xlon write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(3), transpose(dlat_s))
        else
          istatus = nf90_put_var(incout, ivar(3), transpose(dlat))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlat write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(4), transpose(dlon_s))
        else
          istatus = nf90_put_var(incout, ivar(4), transpose(dlon))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dlon write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(5), transpose(htgrid_s))
        else
          istatus = nf90_put_var(incout, ivar(5), transpose(htgrid))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable ht write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(6), transpose(htsdgrid_s))
        else
          istatus = nf90_put_var(incout, ivar(6), transpose(htsdgrid))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable htsd write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(7), transpose(lndout_s))
        else
          istatus = nf90_put_var(incout, ivar(7), transpose(lndout))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable landuse write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(8), transpose(xmap_s))
        else
          istatus = nf90_put_var(incout, ivar(8), transpose(xmap))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable xmap write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(9), transpose(dmap_s))
        else
          istatus = nf90_put_var(incout, ivar(9), transpose(dmap))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable dmap write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(10), transpose(coriol_s))
        else
          istatus = nf90_put_var(incout, ivar(10), transpose(coriol))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable coriol write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(11), transpose(snowam_s))
        else
          istatus = nf90_put_var(incout, ivar(11), transpose(snowam))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable snowam write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        if (lsub) then
          istatus = nf90_put_var(incout, ivar(12), transpose(mask_s))
        else
          istatus = nf90_put_var(incout, ivar(12), transpose(mask))
        end if
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error Variable mask write in NetCDF output'
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        if ( aertyp(7:7)=='1' ) then
          if (lsub) then
            istatus = nf90_put_var(incout, itvar(1), transpose(texout_s))
          else
            istatus = nf90_put_var(incout, itvar(1), transpose(texout))
          end if
          if (istatus /= nf90_noerr) then
            write (6,*) 'Error Variable texture write in NetCDF output'
            write (6,*) nf90_strerror(istatus)
            stop
          end if
          istart(1) = 1
          istart(2) = 1
          icount(3) = 1
          icount(2) = iy
          icount(1) = jx
          do i = 1 , ntex
            istart(3) = i
            if (lsub) then
              istatus = nf90_put_var(incout, itvar(2),                  &
                         & transpose(frac_tex_s(:,:,i)),istart,icount)
            else
              istatus = nf90_put_var(incout, itvar(2),                  &
                         & transpose(frac_tex(:,:,i)),istart,icount)
            end if
            if (istatus /= nf90_noerr) then
              write (6,*) 'Error Variable texture_frac write in NetCDF'
              write (6,*) nf90_strerror(istatus)
              stop
            end if
          end do
        end if

        istatus = nf90_close(incout)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error closing NetCDF output ', trim(fname)
          write (6,*) nf90_strerror(istatus)
          stop
        end if

        deallocate(yiy)
        deallocate(xjx)

      end subroutine write_domain

      end module mod_write
