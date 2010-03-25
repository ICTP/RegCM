      subroutine griddef(year,trec)
      use mod_regcm_param , only : aertyp , ibyte , kx
      use mod_emission
      implicit none
!
! Dummy arguments
!
      integer :: trec , year
      intent (in) trec , year
!
! Local variables
!
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
            & centerj , xcla , xclo , dsinm , grdfac , xpla , xplo ,    &
            & xpto , rlatinc , rloninc , xtrul , xtruh
      character(3) , dimension(12) :: cmonth
      integer :: i , isbige , ierr , dograd , ixx , j , jxx , k , nl ,  &
               & month , nx , ny , period , xnspc1a , xnspc1b ,         &
               & xnspc2a , xnspc2b , xnspc3 , xnspc4a , xnspc4b ,       &
               & xnspc5a , xnspc5b
      character(6) :: cprj
      real(4) , dimension(kx+1) :: sigmaf
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
      read (10,rec=1,iostat=ierr) ixx , jxx , nl , dsinm , xcla , xclo ,&
                                & xpla , xplo , grdfac , cprj ,         &
                                & (sigmaf(k),k=1,kx+1) , xpto , dograd ,&
                                & isbige , xtrul , xtruh
      if ( ixx/=ix .or. jxx/=jx ) then
        print * , 'IMPROPER DIMENSION SPECIFICATION (AEROSOL.f)'
        print * , '  icbc.param: ' , ix , jx
        print * , '  DOMAIN.INFO: ' , ixx , jxx
        print * , '  Also check ibyte in icbc.param: ibyte= ' , ibyte
        stop 'Dimensions (subroutine gridml)'
      end if
      read (10,rec=5,iostat=ierr) ((xlat(i,j),j=1,jx),i=1,ix)
      read (10,rec=6,iostat=ierr) ((xlon(i,j),j=1,jx),i=1,ix)
      if ( ierr/=0 ) then
        print * , 'END OF FILE REACHED (AEROSOL.f)'
        print * , '  Check ibyte in icbc.param: ibyte= ' , ibyte
        stop 'EOF (subroutine gridml)'
      end if
!
      if ( dograd==1 ) then
!       inquire(file='../../Input/AERO.ctl',exist=there)
!       if(there) isystm=system('/bin/rm ../../Input/AERO.ctl')
        open (31,file='../../Input/AERO_new.ctl')
        open (32,file='../../Main/emission.param')
        write (32,'(T7,A,I3)') 'INTEGER, PARAMETER :: TREC =' , trec
        write (32,'(T7,A)') 'include ''Commons/emission.def'' '
!       OPEN(31,file='AERO.ctl',status='new')
        write (31,'(a)') 'dset ^AERO_new.dat'
        write (31,'(a)')                                                &
                     &'title AEROSOL fields for RegCM domain, kg/m2/sec'
        if ( isbige==1 ) then
          write (31,'(a)') 'options big_endian'
        else
          write (31,'(a)') 'options little_endian'
        end if
        write (31,'(a)') 'undef -9999.'
        if ( cprj=='LAMCON' .or. cprj=='ROTMER' ) then
          alatmin = 999999.
          alatmax = -999999.
          do j = 1 , jx
            if ( xlat(1,j)<alatmin ) alatmin = xlat(1,j)
            if ( xlat(ix,j)>alatmax ) alatmax = xlat(ix,j)
          end do
          alonmin = 999999.
          alonmax = -999999.
          do i = 1 , ix
            do j = 1 , jx
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
          ny = 2 + nint(abs(alatmax-alatmin)/rlatinc)
          nx = 1 + nint(abs((alonmax-alonmin)/rloninc))
 
          centerj = jx/2.
          centeri = ix/2.
        end if
        if ( cprj=='LAMCON' ) then        ! Lambert projection
          write (31,99001) jx , ix , xcla , xclo , centerj , centeri ,  &
                         & xtrul , xtruh , xclo , dsinm , dsinm
          write (31,99002) nx + 2 , alonmin - rloninc , rloninc
          write (31,99003) ny + 2 , alatmin - rlatinc , rlatinc
        else if ( cprj=='POLSTR' ) then   !
        else if ( cprj=='NORMER' ) then
          write (31,99004) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
          write (31,99005) ix
          write (31,99006) (xlat(i,1),i=1,ix)
        else if ( cprj=='ROTMER' ) then
          write (*,*) 'Note that rotated Mercartor (ROTMER)' ,          &
                     &' projections are not supported by GrADS.'
          write (*,*) '  Although not exact, the eta.u projection' ,    &
                     &' in GrADS is somewhat similar.'
          write (*,*) ' FERRET, however, does support this projection.'
          write (31,99007) jx , ix , xplo , xpla , dsinm/111000. ,      &
                         & dsinm/111000.*.95238
          write (31,99002) nx + 2 , alonmin - rloninc , rloninc
          write (31,99003) ny + 2 , alatmin - rlatinc , rlatinc
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
!             write(32,701) ELE_RETROB(I)
              write (32,99018) ele_retrob(i) , i
              write (32,*)
            end do
          end if
          if ( ipoet ) then
            do i = 1 , xnspc2b
              write (31,99011) ele_poetb(i) ,                           &
                              &'Biogenic emission, POET  '
!             write(32,701) ELE_POETB(I)
              write (32,99018) ele_poetb(i) , i + xnspc1b
              write (32,*)
            end do
          end if
          if ( igfed ) then
            if ( year>2000 ) then
              do i = 1 , xnspc3
                write (31,99011) ele_gfed(i) ,                          &
                                &'Biogenic emission,GFED     '
!               write(32,701) ELE_GFED(I)
                write (32,99018) ele_gfed(i) , i + xnspc1b + xnspc2b
                write (32,*)
              end do
            end if
          end if
          if ( iedgar ) then
            do i = 1 , xnspc4b
              write (31,99011) ele_edgarb(i) ,                          &
                              &'Biogenic emission, EDGAR       '
!             write(32,701) ELE_EDGARB(I)
              write (32,99018) ele_edgarb(i) , i + xnspc1b + xnspc3 +   &
                             & xnspc2b
              write (32,*)
            end do
          end if
!----------------------------------------------------------------
        case ('AER10D0','AER10D1')
          write (31,99012) xnspc1a + xnspc2a + xnspc4a
          if ( iretro ) then
            do i = 1 , xnspc1a
              write (31,99010) ele_retroa(i) ,                          &
                              &'Anthropogenic emission, RETRO  '
!             write(32,700) ELE_RETROA(I)
              write (32,99017) ele_retroa(i) , i
              write (32,*)
            end do
          end if
          if ( ipoet ) then
            do i = 1 , xnspc2a
              write (31,99010) ele_poeta(i) ,                           &
                              &'Anthropogenic emission, POET  '
!             write(32,700) ELE_POETA(I)
              write (32,99017) ele_poeta(i) , i + xnspc1a
              write (32,*)
            end do
          end if
          if ( iedgar ) then
            do i = 1 , xnspc4a
              write (31,99010) ele_edgara(i) ,                          &
                              &'Anthropogenic emission, EDGAR  '
!             write(32,700) ELE_EDGARA(I)
              write (32,99017) ele_edgara(i) , i + xnspc2a + xnspc1a
              write (32,*)
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
!             write(32,700) ELE_RETROA(I)
              write (32,99017) ele_retroa(i) , i
              write (32,*)
            end do
          end if
          if ( ipoet ) then
            do i = 1 , xnspc2a
              write (31,99010) ele_poeta(i) ,                           &
                              &'Anthropogenic emission, POET  '
!             write(32,700) ELE_POETA(I)
              write (32,99017) ele_poeta(i) , i + xnspc1a
              write (32,*)
            end do
          end if
          if ( iedgar ) then
            do i = 1 , xnspc4a
              write (31,99010) ele_edgara(i) ,                          &
                              &'Anthropogenic emission, EDGAR  '
!             write(32,700) ELE_EDGARA(I)
              write (32,99017) ele_edgara(i) , i + xnspc2a + xnspc1a
              write (32,*)
            end do
          end if
          if ( imozart ) then
            do i = 1 , xnspc5a
              write (31,99010) ele_mozrta(i) ,                          &
                              &'Anthropogenic emission, MOZART  '
!             write(32,700) ELE_MOZRTA(I)
              write (32,99017) ele_mozrta(i) , i + xnspc4a + xnspc2a +  &
                             & xnspc1a
            end do
          end if
!------------------START BIOGENIC EMISSION---------------------
          if ( iretro ) then
            do i = 1 , xnspc1b
              write (31,99011) ele_retrob(i) ,                          &
                              &'Biogenic emission, RETRO      '
!             write(32,701) ELE_RETROB(I)
              write (32,99018) ele_retrob(i) , i + xnspc1a + xnspc2a +  &
                             & xnspc4a + xnspc5a
              write (32,*)
            end do
          end if
          if ( ipoet ) then
            do i = 1 , xnspc2b
              write (31,99011) ele_poetb(i) ,                           &
                              &'Biogenic emission, POET  '
!             write(32,701) ELE_POETB(I)
              write (32,99018) ele_poetb(i) , i + xnspc1a + xnspc2a +   &
                             & xnspc4a + xnspc5a + xnspc1b
              write (32,*)
            end do
          end if
          if ( igfed ) then
            do i = 1 , xnspc3
              write (31,99011) ele_gfed(i) ,                            &
                              &'Biogenic emission,GFED          '
!             write(32,701) ELE_GFED(I)
              write (32,99018) ele_gfed(i) , i + xnspc1a + xnspc2a +    &
                             & xnspc4a + +xnspc5a + xnspc1b + xnspc2b
              write (32,*)
            end do
          end if
          if ( iedgar ) then
            do i = 1 , xnspc4b
              write (31,99011) ele_edgarb(i) ,                          &
                              &'Biogenic emission, EDGAR       '
!             write(32,701) ELE_EDGARB(I)
              write (32,99018) ele_edgarb(i) , i + xnspc1a + xnspc2a +  &
                             & xnspc4a + xnspc1b + xnspc5a + xnspc2b +  &
                             & xnspc3
              write (32,*)
            end do
          end if
          if ( imozart ) then
            do i = 1 , xnspc5b
              write (31,99011) ele_mozrtb(i) ,                          &
                              &'Biogenic emission, MOZART'
!             write(32,701) ELE_MOZRTB(I)
              write (32,99018) ele_mozrtb(i) , i + xnspc1a + xnspc2a +  &
                             & xnspc4a + xnspc5a + xnspc1b + xnspc2b +  &
                             & xnspc3 + xnspc4b
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
99013 format (a4,i2,' 0 ',a17)
99014 format (a5,i2,' 0 ',a17)
99015 format (8x,'INTEGER',2x,'ia_',a10)
99016 format (8x,'INTEGER',2x,'ib_',a10)
99017 format (8x,'PARAMETER',2x,'(ia_',a10,'=',i2,')')
99018 format (8x,'PARAMETER',2x,'(ib_',a10,'=',i2,')')
!
      end subroutine griddef
