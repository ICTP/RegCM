!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) ainy later version.
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

      subroutine griddef(year,trec,aertyp,kz)

      use mod_emission

      implicit none
!
! Dummy arguments
!
      integer :: trec , year , kz
      character(7) :: aertyp
      intent (in) trec , year , aertyp , kz
!
! Local variables
!
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
            & centerj , xcla , xclo , dsinm , grdfac , xpla , xplo ,    &
            & xpto , rlatinc , rloninc , xtrul , xtruh
      character(3) , dimension(12) :: cmonth
      integer :: i , isbige , ierr , dograd , nyy , j , nxx , k , nl ,  &
               & month , jnx , iny , period , xnspc1a , xnspc1b ,       &
               & xnspc2a , xnspc2b , xnspc3 , xnspc4a , xnspc4b ,       &
               & xnspc5a , xnspc5b
      character(6) :: cprj
      real(4) , dimension(kz+1) :: sigmaf
!
      data cmonth/'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,       &
         & 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec'/
!
!     new names for number of species
!
      namelist /species/ ele_retroa , ele_retrob , ele_poeta ,          &
      & ele_poetb , ele_gfed , ele_edgara , ele_edgarb , ele_mozrta ,   &
      & ele_mozrtb
!
      rewind (30)
      read (30,nml=species)
      if ( iretro ) then
        xnspc1a = nspc1a
        xnspc1b = nspc1b
      else
        xnspc1a = 0
        xnspc1b = 0
      end if
      if ( ipoet ) then
        xnspc2a = nspc2a
        xnspc2b = nspc2b
      else
        xnspc2a = 0
        xnspc2b = 0
      end if
      if ( iedgar ) then
        xnspc4a = nspc4a
        xnspc4b = nspc4b
      else
        xnspc4a = 0
        xnspc4b = 0
      end if
      if ( imozart ) then
        xnspc5a = nspc5a
        xnspc5b = nspc5b
      else
        xnspc5a = 0
        xnspc5b = 0
      end if
      if ( igfed ) then
        xnspc3 = nspc3
      else
        xnspc3 = 0
      end if
!
      alatmin = 999999.
      alatmax = -999999.
      alonmin = 999999.
      alonmax = -999999.
      iny = 0
      jnx = 0
      read (10,rec=1,iostat=ierr) nyy , nxx , nl , dsinm , xcla , xclo ,&
                                & xpla , xplo , grdfac , cprj ,         &
                                & (sigmaf(k),k=1,kz+1) , xpto , dograd ,&
                                & isbige , xtrul , xtruh
      if ( nyy/=ny .or. nxx/=nx ) then
        print * , 'IMPROPER DIMENSION SPECIFICATION'
        print * , '  namelist   : ' , ny , nx
        print * , '  DOMAIN.INFO: ' , nyy , nxx
        stop
      end if
      read (10,rec=5,iostat=ierr) ((xlat(i,j),j=1,nx),i=1,ny)
      read (10,rec=6,iostat=ierr) ((xlon(i,j),j=1,nx),i=1,ny)
      if ( ierr/=0 ) then
        stop 'EOF'
      end if
!
      if ( dograd==1 ) then
!       inquire(file='../../Input/AERO.ctl',exist=there)
!       if(there) isystm=system('/bin/rm ../../Input/AERO.ctl')
!       OPEN(31,file='AERO.ctl',status='new')
        write (31,'(a)')                                                &
                     &'title AEROSOL fields for RegCM domain, kg/m2/sec'
        if ( isbige==1 ) then
          write (31,'(a)') 'options big_endian'
        else
          write (31,'(a)') 'options little_endian'
        end if
        write (31,'(a)') 'undef -9999.'
        if ( cprj=='LAMCON' .or. cprj=='ROTMER' ) then
          do j = 1 , nx
            if ( xlat(1,j)<alatmin ) alatmin = xlat(1,j)
            if ( xlat(ny,j)>alatmax ) alatmax = xlat(ny,j)
          end do
          do i = 1 , ny
            do j = 1 , nx
              if ( xclo>=0.0 ) then
                if ( xlon(i,j)>=0.0 ) then
                  alonmin = amin1(alonmin,xlon(i,j))
                  alonmax = amax1(alonmax,xlon(i,j))
                else if ( abs(xclo-xlon(i,j))<abs(xclo-(xlon(i,j)+360.))&
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
              else if ( abs(xclo-xlon(i,j))<abs(xclo-(xlon(i,j)-360.)) )&
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
          iny = 2 + nint(abs(alatmax-alatmin)/rlatinc)
          jnx = 1 + nint(abs((alonmax-alonmin)/rloninc))
 
          centerj = nx/2.
          centeri = ny/2.
        end if
        if ( cprj=='LAMCON' ) then        ! Lambert projection
          write (31,99001) nx , ny , xcla , xclo , centerj , centeri ,  &
                         & xtrul , xtruh , xclo , dsinm , dsinm
          write (31,99002) jnx + 2 , alonmin - rloninc , rloninc
          write (31,99003) iny + 2 , alatmin - rlatinc , rlatinc
        else if ( cprj=='POLSTR' ) then   !
        else if ( cprj=='NORMER' ) then
          write (31,99004) nx , xlon(1,1) , xlon(1,2) - xlon(1,1)
          write (31,99005) ny
          write (31,99006) (xlat(i,1),i=1,ny)
        else if ( cprj=='ROTMER' ) then
          write (*,*) 'Note that rotated Mercartor (ROTMER)' ,          &
                     &' projections are not supported by GrADS.'
          write (*,*) '  Although not exact, the eta.u projection' ,    &
                     &' in GrADS is somewhat similar.'
          write (*,*) ' FERRET, however, does support this projection.'
          write (31,99007) nx , ny , xplo , xpla , dsinm/111000. ,      &
                         & dsinm/111000.*.95238
          write (31,99002) jnx + 2 , alonmin - rloninc , rloninc
          write (31,99003) iny + 2 , alatmin - rlatinc , rlatinc
        else
          write (*,*) 'Are you sure your map projection is correct ?'
          stop
        end if
!       write(31,300)KZ,((1013.25-PTOP*10.)*(.5*(SIGMAF(K)+SIGMAF(K+1)))
!       &                                          +PTOP*10.,K=KZ,1,-1)
        write (31,99008) 1 , 1000.
 
!       write(31,300) 1,1000.
!       300  format('zdef ',I1,' levels ',f7.2)
        month = 1
        period = trec
        write (31,99009) period , cmonth(month) , 1960
!       write(31,600) 'so2   ','Anthropogenic SO2 emission, EDGAR      
!       '
        select case (aertyp)
        case ('AER01D0','AER01D1')
          if ( year>2000 ) then
            write (31,99012) xnspc1b + xnspc2b + xnspc3 + xnspc4b
          else
            write (31,99012) xnspc1b + xnspc2b + xnspc4b
          end if
          if ( iretro ) then
            do i = 1 , xnspc1b
              write (31,99011) ele_retrob(i) ,                          &
                              &'Biogenic emission, RETRO      '
            end do
          end if
          if ( ipoet ) then
            do i = 1 , xnspc2b
              write (31,99011) ele_poetb(i) ,                           &
                              &'Biogenic emission, POET  '
            end do
          end if
          if ( igfed ) then
            if ( year>2000 ) then
              do i = 1 , xnspc3
                write (31,99011) ele_gfed(i) ,                          &
                                &'Biogenic emission,GFED     '
              end do
            end if
          end if
          if ( iedgar ) then
            do i = 1 , xnspc4b
              write (31,99011) ele_edgarb(i) ,                          &
                              &'Biogenic emission, EDGAR       '
            end do
          end if
!----------------------------------------------------------------
        case ('AER10D0','AER10D1')
          write (31,99012) xnspc1a + xnspc2a + xnspc4a
          if ( iretro ) then
            do i = 1 , xnspc1a
              write (31,99010) ele_retroa(i) ,                          &
                              &'Anthropogenic emission, RETRO  '
            end do
          end if
          if ( ipoet ) then
            do i = 1 , xnspc2a
              write (31,99010) ele_poeta(i) ,                           &
                              &'Anthropogenic emission, POET  '
            end do
          end if
          if ( iedgar ) then
            do i = 1 , xnspc4a
              write (31,99010) ele_edgara(i) ,                          &
                              &'Anthropogenic emission, EDGAR  '
            end do
          end if
 
        case ('AER11D0','AER11D1')
          if ( year>2000 ) then
!           write(31,500) NSPC1A+NSPC4A+NSPC4B+NSPC3+NSPC2
            write (31,99012) xnspc1a + xnspc1b + xnspc2a + xnspc2b +    &
                           & xnspc3 + xnspc4a + xnspc4b + xnspc5a +     &
                           & xnspc5b
          else
            write (31,99012) xnspc1a + xnspc4a + xnspc4b + xnspc1b +    &
                           & xnspc2a + xnspc2b
          end if
          if ( iretro ) then
            do i = 1 , xnspc1a
              write (31,99010) ele_retroa(i) ,                          &
                              &'Anthropogenic emission, RETRO  '
            end do
          end if
          if ( ipoet ) then
            do i = 1 , xnspc2a
              write (31,99010) ele_poeta(i) ,                           &
                              &'Anthropogenic emission, POET  '
            end do
          end if
          if ( iedgar ) then
            do i = 1 , xnspc4a
              write (31,99010) ele_edgara(i) ,                          &
                              &'Anthropogenic emission, EDGAR  '
            end do
          end if
          if ( imozart ) then
            do i = 1 , xnspc5a
              write (31,99010) ele_mozrta(i) ,                          &
                              &'Anthropogenic emission, MOZART  '
            end do
          end if
!------------------START BIOGENIC EMISSION---------------------
          if ( iretro ) then
            do i = 1 , xnspc1b
              write (31,99011) ele_retrob(i) ,                          &
                              &'Biogenic emission, RETRO      '
            end do
          end if
          if ( ipoet ) then
            do i = 1 , xnspc2b
              write (31,99011) ele_poetb(i) ,                           &
                              &'Biogenic emission, POET  '
            end do
          end if
          if ( igfed ) then
            do i = 1 , xnspc3
              write (31,99011) ele_gfed(i) ,                            &
                              &'Biogenic emission,GFED          '
            end do
          end if
          if ( iedgar ) then
            do i = 1 , xnspc4b
              write (31,99011) ele_edgarb(i) ,                          &
                              &'Biogenic emission, EDGAR       '
            end do
          end if
          if ( imozart ) then
            do i = 1 , xnspc5b
              write (31,99011) ele_mozrtb(i) ,                          &
                              &'Biogenic emission, MOZART'
            end do
          end if
        case default
        end select
 
        write (31,'(a)') 'endvars'
        close (31)
      end if
!
99001 format ('pdef ',i4,1x,i4,1x,'lcc',7(1x,f7.2),1x,2(f7.0,1x))
99002 format ('xdef ',i4,' linear ',f7.2,1x,f7.4)
99003 format ('ydef ',i4,' linear ',f7.2,1x,f7.4)
99004 format ('xdef ',i3,' linear ',f9.4,' ',f9.4)
99005 format ('ydef ',i3,' levels')
99006 format (10F7.2)
99007 format ('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
99008 format ('zdef ',i2,' levels ',30F7.2)
99009 format ('tdef ',i4,' linear 00z16',a3,i4,' 1mo')
99010 format ('a_',a10,'0 99 ',a40)
99011 format ('b_',a10,'0 99 ',a40)
99012 format ('vars ',i2)
!
      end subroutine griddef
