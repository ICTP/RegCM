!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or
!    modify
!    it under the terms of the GNU General Public License as
!    published by
!    the Free Software Foundation, either version 3 of the
!    License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty
!    of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public
!    License
!    along with ICTP RegCM.  If not, see
!    <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      module mod_wrtoxd

      implicit none
      integer :: noutrec
      integer :: iny , jnx , knz
      real(4) , allocatable , dimension(:,:,:) :: oh4,ho24,o34,no34,h2o24

      data noutrec /0/

      contains

      subroutine init_outoxd(jx,iy,kz)
      implicit none
      integer , intent(in) :: jx , iy , kz
      iny = iy
      jnx = jx
      knz = kz
      allocate(oh4(jnx,iny,knz))
      allocate(ho24(jnx,iny,knz))
      allocate(o34(jnx,iny,knz))
      allocate(no34(jnx,iny,knz))
      allocate(h2o24(jnx,iny,knz))
      end subroutine init_outoxd

      subroutine free_outoxd
      implicit none
      deallocate(oh4)
      deallocate(ho24)
      deallocate(o34)
      deallocate(no34)
      deallocate(h2o24)
      end subroutine free_outoxd

      subroutine writeox(ptop,idate)
      implicit none
!
! Dummy arguments
!
      integer :: idate
      real(4) :: ptop
      intent (in) idate , ptop
!
! Local variables
!
      integer :: i , j , k
!
!     THIS ROUTINE WRITES OUT AN INPUT FILE FOR THE RCM
!
!     PRINT *,'WRITING OUTPUT:  IDATE= ',IDATE

      noutrec = noutrec + 1
      write(*,*)'INSID',noutrec,idate
      write (64,rec=noutrec) idate , jnx , iny , knz
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((oh4(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((ho24(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((o34(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((no34(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((h2o24(i,j,k),i=1,jnx),j=1,iny)
      end do
!
      end subroutine writeox
      subroutine gradsctl(finame,idate,inumber)
      use mod_dynparam
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
      alatmin = 999999.
      alatmax = -999999.
      alonmin = 999999.
      alonmax = -999999.
      nx = 0
      ny = 0
      open (71,file=trim(finame)//'.CTL',status='replace')
      write (71,'(a,a,a,i10)') 'dset ^',trim(domname),'_OXBC',idate
      write (71,'(a)') 'title oxidant fields for RegCM domain'
      if ( ibigend==1 ) then
        write (71,'(a)') 'options big_endian'
      else
        write (71,'(a)') 'options little_endian'
      end if
      write (71,'(a)') 'undef -9999.'
      if ( iproj=='LAMCON' .or. iproj=='ROTMER' ) then
        do j = 1 , jx
          if ( xlat(j,1)<alatmin ) alatmin = xlat(j,1)
          if ( xlat(j,iy)>alatmax ) alatmax = xlat(j,iy)
        end do
        do i = 1 , iy
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
        centeri = iy/2.
      end if
      if ( iproj=='LAMCON' ) then       ! Lambert projection
        write (71,99001) jx , iy , clat , clon , centerj , centeri ,    &
                       & truelatl , truelath , clon , delx , delx
        write (71,99002) nx + 2 , alonmin - rloninc , rloninc
        write (71,99003) ny + 2 , alatmin - rlatinc , rlatinc
      else if ( iproj=='POLSTR' ) then  !
      else if ( iproj=='NORMER' ) then
        write (71,99004) jx , xlon(1,1) , xlon(2,1) - xlon(1,1)
        write (71,99005) iy
        write (71,99006) (xlat(1,i),i=1,iy)
      else if ( iproj=='ROTMER' ) then
        write (*,*) 'Note that rotated Mercartor (ROTMER)' ,            &
                   &' projections are not supported by GrADS.'
        write (*,*) '  Although not exact, the eta.u projection' ,      &
                   &' in GrADS is somewhat similar.'
        write (*,*) ' FERRET, however, does support this projection.'
        write (71,99007) jx , iy , plon , plat , delx/111000. ,         &
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
      write (71,99010) 6 
      write (71,'(a)') 'date 0 99 header information'
      write (71,99012) 'oh   ' , kz , ' moleculs cm-3  '
      write (71,99012) 'ho2  ' , kz , ' moleculs cm-3  '
      write (71,99012) 'o3   ' , kz , 'volum mixing ratio'
      write (71,99012) 'no3  ' , kz , ' moleculs cm-3'
      write (71,99012) 'h2o2 ' , kz , 'volum mixing ratio'
      write (71,'(a)') 'endvars'
      close (71)
99001 format ('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
99002 format ('xdef ',i4,' linear ',f7.2,1x,f7.4)
99003 format ('ydef ',i4,' linear ',f7.2,1x,f7.4)
99004 format ('xdef ',i4,' linear ',f9.4,' ',f9.4)
99005 format ('ydef ',i4,' levels')
99006 format (10F7.2)
99007 format ('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
99008 format ('zdef ',i2,' levels ',30F7.2)
99009 format ('tdef ',i4,' linear ',i2,'z',a2,a3,i4,' 6hr')
99010 format ('vars ',i1)
99011 format ('vars ',i2)
99012 format (a4,1x,i2,' 0 ',a17)
99013 format (a4,i2,' 33,100 ',a17)
99014 format (a4,i2,' 34,100 ',a17)
99015 format (a4,'0 99 ',a26)
!
      end subroutine gradsctl




      end module mod_wrtoxd
