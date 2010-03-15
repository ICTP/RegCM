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

      subroutine gradsctlb(finame,idate,inumber)
      use mod_param
      use mod_grid
      implicit none
!
! Dummy arguments
!
      character(*) :: finame
      integer :: idate , inumber
      intent (in) finame , idate , inumber
!
! Local variables
!
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
               & centerj , rlatinc , rloninc
      character(2) , dimension(31) :: cday
      character(3) , dimension(12) :: cmonth
      integer :: i , j , k , month , nday , nhour , nx , ny , nyear
!
!
      data cday/'01' , '02' , '03' , '04' , '05' , '06' , '07' , '08' , &
          &'09' , '10' , '11' , '12' , '13' , '14' , '15' , '16' ,      &
         & '17' , '18' , '19' , '20' , '21' , '22' , '23' , '24' ,      &
         & '25' , '26' , '27' , '28' , '29' , '30' , '31'/
      data cmonth/'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,       &
         & 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec'/
!
!
!     DOMAIN VARIABLES FOR RCM HORIZONTAL GRID
!
      open (71,file=finame//'.CTL',status='replace')
      write (71,'(a)') 'dset ^'//finame(13:26)
      write (71,'(a)') 'title ICBC fields for RegCM domain'
      if ( ibigend==1 ) then
        write (71,'(a)') 'options big_endian'
      else
        write (71,'(a)') 'options little_endian'
      end if
      write (71,'(a)') 'undef -9999.'
      if ( iproj=='LAMCON' .or. iproj=='ROTMER' ) then
        alatmin = 999999.
        alatmax = -999999.
        do j = 1 , jx
          if ( xlat(j,1)<alatmin ) alatmin = xlat(j,1)
          if ( xlat(j,ix)>alatmax ) alatmax = xlat(j,ix)
        end do
        alonmin = 999999.
        alonmax = -999999.
        do i = 1 , ix
          do j = 1 , jx
            if ( clon>=0.0 ) then
              if ( xlon(j,i)>=0.0 ) then
                alonmin = amin1(alonmin,xlon(j,i))
                alonmax = amax1(alonmax,xlon(j,i))
              else if ( abs(clon-xlon(j,i))<abs(clon-(xlon(j,i)+360.)) )&
                      & then
                alonmin = amin1(alonmin,xlon(j,i))
                alonmax = amax1(alonmax,xlon(j,i))
              else
                alonmin = amin1(alonmin,xlon(j,i)+360.)
                alonmax = amax1(alonmax,xlon(j,i)+360.)
              end if
            else if ( xlon(j,i)<0.0 ) then
              alonmin = amin1(alonmin,xlon(j,i))
              alonmax = amax1(alonmax,xlon(j,i))
            else if ( abs(clon-xlon(j,i))<abs(clon-(xlon(j,i)-360.)) )  &
                    & then
              alonmin = amin1(alonmin,xlon(j,i))
              alonmax = amax1(alonmax,xlon(j,i))
            else
              alonmin = amin1(alonmin,xlon(j,i)-360.)
              alonmax = amax1(alonmax,xlon(j,i)-360.)
            end if
          end do
        end do
        rlatinc = delx*0.001/111./2.
        rloninc = delx*0.001/111./2.
        ny = 2 + nint(abs(alatmax-alatmin)/rlatinc)
        nx = 1 + nint(abs((alonmax-alonmin)/rloninc))
 
        centerj = jx/2.
        centeri = ix/2.
      end if
      if ( iproj=='LAMCON' ) then       ! Lambert projection
        write (71,99001) jx , ix , clat , clon , centerj , centeri ,    &
                       & truelatl , truelath , clon , delx , delx
        write (71,99002) nx + 2 , alonmin - rloninc , rloninc
        write (71,99003) ny + 2 , alatmin - rlatinc , rlatinc
      else if ( iproj=='POLSTR' ) then  !
      else if ( iproj=='NORMER' ) then
        write (71,99004) jx , xlon(1,1) , xlon(2,1) - xlon(1,1)
        write (71,99005) ix
        write (71,99006) (xlat(1,i),i=1,ix)
      else if ( iproj=='ROTMER' ) then
        write (*,*) 'Note that rotated Mercartor (ROTMER)' ,            &
                   &' projections are not supported by GrADS.'
        write (*,*) '  Although not exact, the eta.u projection' ,      &
                   &' in GrADS is somewhat similar.'
        write (*,*) ' FERRET, however, does support this projection.'
        write (71,99007) jx , ix , plon , plat , delx/111000. ,         &
                       & delx/111000.*.95238
        write (71,99002) nx + 2 , alonmin - rloninc , rloninc
        write (71,99003) ny + 2 , alatmin - rlatinc , rlatinc
      else
        write (*,*) 'Are you sure your map projection is correct ?'
        stop
      end if
      write (71,99008) kz , ((1013.25-ptop*10.)*sigma2(k)+ptop*10.,k=kz,&
                     & 1,-1)
      nyear = idate/1000000
      month = (idate-nyear*1000000)/10000
      nday = (idate-nyear*1000000-month*10000)/100
      nhour = mod(idate,100)
      write (71,99009) inumber , nhour , cday(nday) , cmonth(month) ,    &
                     & nyear
      if ( lsmtyp=='USGS' ) then
        write (71,99011) 20
      else
        write (71,99010) 7
      end if
      write (71,'(a)') 'date 0 99 header information'
      write (71,99012) 'u   ' , kz , 'westerly wind    '
      write (71,99012) 'v   ' , kz , 'southerly wind   '
      write (71,99012) 't   ' , kz , 'air temperature  '
      write (71,99012) 'q   ' , kz , 'specific moisture'
      write (71,99013) 'px  ' , 'surface pressure           '
      write (71,99013) 'ts  ' , 'surface air temperature    '
      if ( lsmtyp=='USGS' ) then
        write (71,99013) 'qs1 ' , 'soil moisture level 1      '
        write (71,99013) 'qs2 ' , 'soil moisture level 2      '
        write (71,99013) 'qs3 ' , 'soil moisture level 3      '
        write (71,99013) 'qs4 ' , 'soil moisture level 4      '
        write (71,99013) 'ti1 ' , 'ice  temperature level 1   '
        write (71,99013) 'ti2 ' , 'ice  temperature level 2   '
        write (71,99013) 'ti3 ' , 'ice  temperature level 3   '
        write (71,99013) 'ti4 ' , 'ice  temperature level 4   '
        write (71,99013) 'ts1 ' , 'soil temperature level 1   '
        write (71,99013) 'ts2 ' , 'soil temperature level 2   '
        write (71,99013) 'ts3 ' , 'soil temperature level 3   '
        write (71,99013) 'ts4 ' , 'soil temperature level 4   '
        write (71,99013) 'snd ' , 'snow depth (in metre)      '
      end if
      write (71,'(a)') 'endvars'
      close (71)
99001 format ('pdef ',i4,1x,i4,1x,'lcc',7(1x,f7.2),1x,2(f7.0,1x))
99002 format ('xdef ',i4,' linear ',f7.2,1x,f7.4)
99003 format ('ydef ',i4,' linear ',f7.2,1x,f7.4)
99004 format ('xdef ',i3,' linear ',f9.4,' ',f9.4)
99005 format ('ydef ',i3,' levels')
99006 format (10F7.2)
99007 format ('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
99008 format ('zdef ',i2,' levels ',30F7.2)
99009 format ('tdef ',i4,' linear ',i2,'z',a2,a3,i4,' 6hr')
99010 format ('vars ',i1)
99011 format ('vars ',i2)
99012 format (a4,i2,' 0 ',a17)
99013 format (a4,'0 99 ',a26)
!
      end subroutine gradsctlb
