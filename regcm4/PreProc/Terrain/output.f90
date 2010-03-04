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

      subroutine output(nunitc,iy,jx,nsg,dsinm,clat,clon,plat,plon,     &
                      & iproj,htgrid,htsdgrid,lndout,xlat,xlon,dlat,    &
                      & dlon,xmap,dattyp,dmap,coriol,snowam,igrads,     &
                      & ibigend,kz,sigma,mask,ptop,htgrid_s,lndout_s,   &
                      & ibyte,ngrid,truelatl,truelath,grdfac,filout,    &
                      & lsmtyp,sanda,sandb,claya,clayb,frac_lnd,nveg,   &
                      & aertyp,texout,frac_tex,ntex)
      implicit none
!
! Dummy arguments
!
      character(7) :: aertyp
      real(4) :: clat , clon , dsinm , grdfac , plat , plon , ptop ,    &
               & truelath , truelatl
      character(5) :: dattyp
      character(50) :: filout
      integer :: ibigend , ibyte , igrads , iy , jx , kz , ngrid , nsg ,&
               & ntex , nunitc , nveg
      character(6) :: iproj
      character(4) :: lsmtyp
      real(4) , dimension(iy,jx) :: claya , clayb , coriol , dlat ,     &
                                  & dlon , dmap , htgrid , htsdgrid ,   &
                                  & lndout , mask , sanda , sandb ,     &
                                  & snowam , texout , xlat , xlon , xmap
      real(4) , dimension(iy,jx,nveg) :: frac_lnd
      real(4) , dimension(iy,jx,ntex) :: frac_tex
      real(4) , dimension(iy*nsg,jx*nsg) :: htgrid_s , lndout_s
      real(4) , dimension(kz+1) :: sigma
      intent (in) aertyp , clat , claya , clayb , clon , coriol ,       &
                & dattyp , dlat , dlon , dmap , dsinm , filout ,        &
                & frac_lnd , frac_tex , grdfac , htgrid , htsdgrid ,    &
                & ibigend , ibyte , igrads , iproj , iy , jx , kz ,     &
                & lndout , lsmtyp , ngrid , nsg , ntex , nunitc , nveg ,&
                & plat , plon , ptop , sanda , sandb , snowam , texout ,&
                & truelath , truelatl , xlat , xlon , xmap
      intent (inout) htgrid_s , lndout_s , mask , sigma
!
! Local variables
!
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
               & centerj , htave , htgrid_a , lat0 , lat1 , lon0 ,      &
               & lon1 , rlatinc , rloninc
      character(25) :: fsubname
      integer :: i , i0 , iii , j , j0 , jjj , k , m , n , nx , ny
      logical :: there
!
      if ( kz==14 ) then                      ! RegCM2
        sigma(1) = 0.0
        sigma(2) = 0.04
        sigma(3) = 0.10
        sigma(4) = 0.17
        sigma(5) = 0.25
        sigma(6) = 0.35
        sigma(7) = 0.46
        sigma(8) = 0.56
        sigma(9) = 0.67
        sigma(10) = 0.77
        sigma(11) = 0.86
        sigma(12) = 0.93
        sigma(13) = 0.97
        sigma(14) = 0.99
        sigma(15) = 1.0
      else if ( kz==18 ) then                 ! RegCM3, default
        sigma(1) = 0.0
        sigma(2) = 0.05
        sigma(3) = 0.10
        sigma(4) = 0.16
        sigma(5) = 0.23
        sigma(6) = 0.31
        sigma(7) = 0.39
        sigma(8) = 0.47
        sigma(9) = 0.55
        sigma(10) = 0.63
        sigma(11) = 0.71
        sigma(12) = 0.78
        sigma(13) = 0.84
        sigma(14) = 0.89
        sigma(15) = 0.93
        sigma(16) = 0.96
        sigma(17) = 0.98
        sigma(18) = 0.99
        sigma(19) = 1.0
      else if ( kz==23 ) then                 ! MM5V3
        sigma(1) = 0.0
        sigma(2) = 0.05
        sigma(3) = 0.1
        sigma(4) = 0.15
        sigma(5) = 0.2
        sigma(6) = 0.25
        sigma(7) = 0.3
        sigma(8) = 0.35
        sigma(9) = 0.4
        sigma(10) = 0.45
        sigma(11) = 0.5
        sigma(12) = 0.55
        sigma(13) = 0.6
        sigma(14) = 0.65
        sigma(15) = 0.7
        sigma(16) = 0.75
        sigma(17) = 0.8
        sigma(18) = 0.85
        sigma(19) = 0.89
        sigma(20) = 0.93
        sigma(21) = 0.96
        sigma(22) = 0.98
        sigma(23) = 0.99
        sigma(24) = 1.0
      else
        write (*,*) 'You vertical level number is not 14, 18, or 23'
        write (*,*) 'Please set your sigma parameters in OUTPUT'
        stop
      end if
 
      if ( nsg>1 ) then
        if ( nsg<10 ) then
          write (fsubname,99001) nsg
        else
          write (fsubname,99002) nsg
        end if
        write (*,*) fsubname
        inquire (file=fsubname,exist=there)
        if ( .not.there ) then
          write (*,*) 'Subgrid Terrain and Landuse must be available'
          stop
        end if
        open (19,file=fsubname,form='unformatted',                      &
            & recl=iy*jx*nsg*nsg*ibyte,access='direct')
        read (19,rec=2) ((htgrid_s(i,j),j=1,jx*nsg),i=1,iy*nsg)
        read (19,rec=4) ((lndout_s(i,j),j=1,jx*nsg),i=1,iy*nsg)
        do i = 1 , iy*nsg
          do j = 1 , jx*nsg
            if ( (lsmtyp=='BATS' .and. (lndout_s(i,j)>20. .or. lndout_s(&
               & i,j)<0.)) .or.                                         &
               & (lsmtyp=='USGS' .and. (lndout_s(i,j)>25. .or.          &
               & lndout_s(i,j)<0.)) ) then
              print * , i , j , lndout_s(i,j)
              stop 999
            end if
          end do
        end do
        do i = 1 , iy
          do j = 1 , jx
            i0 = (i-1)*nsg
            j0 = (j-1)*nsg
            htave = 0.0
            do m = 1 , nsg
              do n = 1 , nsg
                if ( lsmtyp=='BATS' ) then
                  if ( htgrid(i,j)<0.1 .and.                            &
                     & (lndout(i,j)>14.5 .and. lndout(i,j)<15.5) ) then
                    htgrid_s(i0+m,j0+n) = 0.0
                    lndout_s(i0+m,j0+n) = 15.
                  end if
                else if ( lsmtyp=='USGS' ) then
                  if ( htgrid(i,j)<0.1 .and. lndout(i,j)>24.5 ) then
                    htgrid_s(i0+m,j0+n) = 0.0
                    lndout_s(i0+m,j0+n) = 25.
                  end if
                else
                end if
                htave = htave + htgrid_s(i0+m,j0+n)
              end do
            end do
            htgrid_a = htave/float(nsg*nsg)
            do m = 1 , nsg
              do n = 1 , nsg
                htgrid_s(i0+m,j0+n) = htgrid_s(i0+m,j0+n) - htgrid_a +  &
                                    & htgrid(i,j)
              end do
            end do
          end do
        end do
        write (19,rec=2) ((htgrid_s(i,j),j=1,jx*nsg),i=1,iy*nsg)
        write (19,rec=4) ((lndout_s(i,j),j=1,jx*nsg),i=1,iy*nsg)
        close (19)
      end if
 
      do i = 1 , iy
        do j = 1 , jx
          if ( (lsmtyp=='BATS' .and. (lndout(i,j)>20. .or. lndout(i,j)< &
             & 0.)) .or.                                                &
             & (lsmtyp=='USGS' .and. (lndout(i,j)>25. .or. lndout(i,j)  &
             & <0.)) ) then
            print * , i , j , lndout(i,j)
            stop 999
          end if
        end do
      end do
      if ( ngrid/=nsg .or. (dattyp/='FVGCM' .and. dattyp/='NRP2W' .and. &
         & dattyp/='GFS11' .and. dattyp/='EH5OM') ) then
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
      if ( ngrid==nsg ) then
        open (23,file='../ICBC/icbcWIN.param')
        write (23,'(a)') '      INTEGER III'
        write (23,'(a)') '      INTEGER JJJ'
        if ( dattyp=='NRP2W' ) then
          iii = nint((lon1-lon0)/2.5) + 1
          jjj = nint((lat1-lat0)/2.5) + 1
          write (23,99003) 'III   =' , iii
          write (23,99003) 'JJJ   =' , jjj
        else
          write (23,99003) 'III   =' , 1
          write (23,99003) 'JJJ   =' , 1
        end if
        close (23)
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
          if ( lsmtyp=='BATS' ) then
            if ( lndout(i,j)>13.5 .and. lndout(i,j)<15.5 ) then
              mask(i,j) = 0.0
            else
              mask(i,j) = 2.0
            end if
          else if ( lsmtyp=='USGS' ) then
            if ( lndout(i,j)>24.5 ) then
              mask(i,j) = 0.0
            else
              mask(i,j) = 2.0
            end if
          else
          end if
        end do
      end do
      write (nunitc,rec=13) ((mask(i,j),j=1,jx),i=1,iy)
      if ( lsmtyp=='USGS' ) then
        write (nunitc,rec=14) ((sanda(i,j),j=1,jx),i=1,iy)
        write (nunitc,rec=15) ((sandb(i,j),j=1,jx),i=1,iy)
        write (nunitc,rec=16) ((claya(i,j),j=1,jx),i=1,iy)
        write (nunitc,rec=17) ((clayb(i,j),j=1,jx),i=1,iy)
        do k = 1 , nveg
          write (nunitc,rec=17+k) ((frac_lnd(i,j,k),j=1,jx),i=1,iy)
        end do
      end if
      if ( aertyp(7:7)=='1' ) then
        if ( lsmtyp=='BATS' ) then
          write (nunitc,rec=14) ((texout(i,j),j=1,jx),i=1,iy)
          do k = 1 , ntex
            write (nunitc,rec=14+k) ((frac_tex(i,j,k),j=1,jx),i=1,iy)
          end do
        else if ( lsmtyp=='USGS' ) then
          write (nunitc,rec=18+nveg) ((texout(i,j),j=1,jx),i=1,iy)
          do k = 1 , ntex
            write (nunitc,rec=18+nveg+k)                                &
                 & ((frac_tex(i,j,k),j=1,jx),i=1,iy)
          end do
        else
        end if
      end if
      close (nunitc)
 
      if ( igrads==1 ) then
        if ( ngrid==nsg ) then
          write (31,99004)
        else if ( ngrid<10 ) then
          write (31,99005) filout(13:24)
        else
          write (31,99006) filout(13:25)
        end if
        write (31,99007)
        if ( ibigend==1 ) then
          write (31,99008)
        else
          write (31,99009)
        end if
        write (31,99010)
        if ( iproj=='LAMCON' .or. iproj=='ROTMER' ) then
          alatmin = 999999.
          alatmax = -999999.
          do j = 1 , jx
            if ( xlat(1,j)<alatmin ) alatmin = xlat(1,j)
            if ( xlat(iy,j)>alatmax ) alatmax = xlat(iy,j)
          end do
          alonmin = 999999.
          alonmax = -999999.
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
          write (31,99011) jx , iy , clat , clon , centerj , centeri ,  &
                         & truelatl , truelath , clon , dsinm , dsinm
          write (31,99012) nx + 2 , alonmin - rloninc , rloninc
          write (31,99013) ny + 2 , alatmin - rlatinc , rlatinc
        else if ( iproj=='POLSTR' ) then   !
        else if ( iproj=='NORMER' ) then
          write (31,99014) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
          write (31,99015) iy
          write (31,99016) (xlat(i,1),i=1,iy)
        else if ( iproj=='ROTMER' ) then
          write (*,*) 'Note that rotated Mercartor (ROTMER)' ,          &
                     &' projections are not supported by GrADS.'
          write (*,*) '  Although not exact, the eta.u projection' ,    &
                     &' in GrADS is somewhat similar.'
          write (*,*) ' FERRET, however, does support this projection.'
          write (31,99017) jx , iy , plon , plat , dsinm/111000. ,      &
                         & dsinm/111000.*.95238
          write (31,99012) nx + 2 , alonmin - rloninc , rloninc
          write (31,99013) ny + 2 , alatmin - rlatinc , rlatinc
        else
          write (*,*) 'Are you sure your map projection is right ?'
          stop
        end if
        write (31,99018) 1 , 1000.
        write (31,99019) 1
        if ( lsmtyp=='BATS' ) then
          if ( aertyp(7:7)=='1' ) then
            write (31,99020) 13 + ntex + 1
          else
            write (31,99020) 13
          end if
        else if ( lsmtyp=='USGS' ) then
          if ( aertyp(7:7)=='1' ) then
            write (31,99020) 42 + ntex + 1
          else
            write (31,99020) 42
          end if
        else
        end if
        write (31,99021) 'head    ' , 'header information         '
        write (31,99021) 'ht      ' , 'surface elevation          '
        write (31,99021) 'htsd    ' , 'surface elevation std. dev.'
        write (31,99021) 'landuse ' , 'surface landuse type       '
        write (31,99021) 'xlat    ' , 'latitude  of cross points  '
        write (31,99021) 'xlon    ' , 'longitude of cross points  '
        write (31,99021) 'dlat    ' , 'latitude  of dot points    '
        write (31,99021) 'dlon    ' , 'longitude of dot points    '
        write (31,99021) 'xmap    ' , 'map factors of cross points'
        write (31,99021) 'dmap    ' , 'map factors of dot points  '
        write (31,99021) 'coriol  ' , 'coriol force               '
        write (31,99021) 'snowam  ' , 'initial snow amount        '
        write (31,99021) 'mask    ' , 'land/sea mask              '
        if ( lsmtyp=='USGS' ) then
          write (31,99021) 'sanda   ' , 'sand percentage (0-30cm)   '
          write (31,99021) 'sandb   ' , 'sand percentage (30-100cm) '
          write (31,99021) 'claya   ' , 'clay percentage (0-30cm)   '
          write (31,99021) 'clayb   ' , 'clay percentage (30-100cm) '
          write (31,99021) 'per1    ' , 'percentage landuse type 1  '
          write (31,99021) 'per2    ' , 'percentage landuse type 2  '
          write (31,99021) 'per3    ' , 'percentage landuse type 3  '
          write (31,99021) 'per4    ' , 'percentage landuse type 4  '
          write (31,99021) 'per5    ' , 'percentage landuse type 5  '
          write (31,99021) 'per6    ' , 'percentage landuse type 6  '
          write (31,99021) 'per7    ' , 'percentage landuse type 7  '
          write (31,99021) 'per8    ' , 'percentage landuse type 8  '
          write (31,99021) 'per9    ' , 'percentage landuse type 9  '
          write (31,99021) 'per10   ' , 'percentage landuse type 10 '
          write (31,99021) 'per11   ' , 'percentage landuse type 11 '
          write (31,99021) 'per12   ' , 'percentage landuse type 12 '
          write (31,99021) 'per13   ' , 'percentage landuse type 13 '
          write (31,99021) 'per14   ' , 'percentage landuse type 14 '
          write (31,99021) 'per15   ' , 'percentage landuse type 15 '
          write (31,99021) 'per16   ' , 'percentage landuse type 16 '
          write (31,99021) 'per17   ' , 'percentage landuse type 17 '
          write (31,99021) 'per18   ' , 'percentage landuse type 18 '
          write (31,99021) 'per19   ' , 'percentage landuse type 19 '
          write (31,99021) 'per20   ' , 'percentage landuse type 20 '
          write (31,99021) 'per21   ' , 'percentage landuse type 21 '
          write (31,99021) 'per22   ' , 'percentage landuse type 22 '
          write (31,99021) 'per23   ' , 'percentage landuse type 23 '
          write (31,99021) 'per24   ' , 'percentage landuse type 24 '
          write (31,99021) 'per25   ' , 'percentage landuse type 24 '
        end if
        if ( aertyp(7:7)=='1' ) then
          write (31,99021) 'texture ' , 'soil texture               '
          write (31,99021) 'text01  ' , 'Sand              frac.    '
          write (31,99021) 'text02  ' , 'Loamy Sand        frac.    '
          write (31,99021) 'text03  ' , 'Sandy Loam        frac.    '
          write (31,99021) 'text04  ' , 'Silt Loam         frac.    '
          write (31,99021) 'text05  ' , 'Silt              frac.    '
          write (31,99021) 'text06  ' , 'Loam              frac.    '
          write (31,99021) 'text07  ' , 'Sandy Clay Loam   frac.    '
          write (31,99021) 'text08  ' , 'Silty Clay Loam   frac.    '
          write (31,99021) 'text09  ' , 'Clay Loam         frac.    '
          write (31,99021) 'text10  ' , 'Sandy Clay        frac.    '
          write (31,99021) 'text11  ' , 'Silty Clay        frac.    '
          write (31,99021) 'text12  ' , 'Clay              frac.    '
          write (31,99021) 'text13  ' , 'OM                frac.    '
          write (31,99021) 'text14  ' , 'Water             frac.    '
          write (31,99021) 'text15  ' , 'Bedrock           frac.    '
          write (31,99021) 'text16  ' , 'Other             frac.    '
          write (31,99021) 'text17  ' , 'No data           frac.    '
        end if
        write (31,99022)
        close (31)
      end if
99001 format ('../../Input/DOMAIN',i1,'.INFO')
99002 format ('../../Input/DOMAIN',i2,'.INFO')
99003 format ('      parameter(',a8,i4,')')
99004 format ('dset ^DOMAIN.INFO')
99005 format ('dset ^',a12)
99006 format ('dset ^',a13)
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
