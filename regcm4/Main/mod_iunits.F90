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

      module mod_iunits
      implicit none
!
      integer , parameter :: nfmax = 999

      character(256) :: ffout
      character(256) :: ffin
!
      integer :: iutbat , iutbc , iutchem , iutchsrc , iutdat , iutin , &
               & iutin1 , iutlak , iutopt , iutrad , iutrs , iutsav ,   &
               & iutsub , mindisp
!
      integer :: iin , iout , lcount , numpts
!
      integer :: nrcout
      integer :: nrcbat
      integer :: nrcsub
      integer :: nrcchem
      integer :: nrcrad

      contains

      subroutine outheadname
        use mod_dynparam
        implicit none
        ffout = trim(dirout)//pthsep//trim(domname)//'_OUT_HEAD'
      end subroutine outheadname

      subroutine indomain
        use mod_dynparam
        implicit none
        ffin = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
      end subroutine indomain

      subroutine insubdom
        use mod_dynparam
        implicit none
        character(3) :: sbstring
        write (sbstring,'(i0.3)') nsg
        ffin = trim(dirter)//pthsep//trim(domname)//'_DOMAIN'//         &
              & sbstring//'.nc'
      end subroutine insubdom

      subroutine inaero
        use mod_dynparam
        implicit none
        ffin = trim(dirglob)//pthsep//trim(domname)//'_AERO.dat'
      end subroutine inaero

      subroutine outname(ctype,idate)
        use mod_dynparam
        implicit none
        character(*) , intent(in) :: ctype
        integer , intent(in) :: idate
        character(32) :: fbname
        write (fbname,99001) trim(ctype) , idate
        ffout = trim(dirout)//pthsep//trim(domname)//'_'//trim(fbname)
99001 format (a,'.',i10)
      end subroutine outname

      subroutine inname(ctype,idate)
        use mod_dynparam
        implicit none
        character(*) , intent(in) :: ctype
        integer , intent(in) :: idate
        character(32) :: fbname
        write (fbname,99001) trim(ctype) , idate
        ffin = trim(dirglob)//pthsep//trim(domname)//'_'//trim(fbname)
99001 format (a,'.',i10)
      end subroutine inname

      subroutine mkfile
 
      use mod_dynparam
      use mod_param1
      use mod_param2
      use mod_date
      use mod_bats
      use mod_mainchem
      use mod_outrad
      implicit none
!
! Local variables
!
      integer , dimension(nfmax) :: idate1d , imo , iyr
      integer :: idatepp , n , nmo
!
      print * , ' '
      print * , '******* OPENING NEW OUTPUT FILES:' , idatex
      if ( iftape ) then
        close (iutdat)
        call outname('ATM',idatex)
        if ( iotyp.eq.1 ) then
#ifdef BAND
          open (iutdat,file=ffout,status='replace',form='unformatted',  &
              & recl=iym2*jx*ibyte,access='direct')
#else
          open (iutdat,file=ffout,status='replace',form='unformatted',  &
              & recl=iym2*jxm2*ibyte,access='direct')
#endif
          nrcout = 0
        else if ( iotyp.eq.2 ) then
          open (iutdat,file=ffout,status='replace',form='unformatted')
        else
        end if
        call gradsout(trim(ffout))
        print * , 'OPENING NEW OUT FILE: ',trim(ffout)
      end if
 
      if ( ifbat ) then
        close (iutbat)
        call outname('SRF',idatex)
        if ( iotyp.eq.1 ) then
#ifdef BAND
          open (iutbat,file=ffout,status='replace',form='unformatted',  &
              & recl=iym2*jx*ibyte,access='direct')
#else
          open (iutbat,file=ffout,status='replace',form='unformatted',  &
              & recl=iym2*jxm2*ibyte,access='direct')
#endif
          nrcbat = 0
        else if ( iotyp.eq.2 ) then
          open (iutbat,file=ffout,status='replace',form='unformatted')
        else
        end if
        call gradsbat(trim(ffout))
        print * , 'OPENING NEW SRF FILE: ',trim(ffout)
      end if
 
      if ( nsg.gt.1 .and. ifsub ) then
        close (iutsub)
        call outname('SUB',idatex)
        if ( iotyp.eq.1 ) then
#ifdef BAND
          open (iutsub,file=ffout,status='replace',form='unformatted',  &
              & recl=iym2*jx*nnsg*ibyte,access='direct')
#else
          open (iutsub,file=ffout,status='replace',form='unformatted',  &
              & recl=iym2*jxm2*nnsg*ibyte,access='direct')
#endif
          nrcsub = 0
        else if ( iotyp.eq.2 ) then
          open (iutsub,file=ffout,status='replace',form='unformatted')
        else
        end if
        call gradssub(trim(ffout))
        print * , 'OPENING NEW SUB FILE: ',trim(ffout)
      end if
 
      if ( ifrad ) then
        close (iutrad)
        call outname('RAD',idatex)
        if ( iotyp.eq.1 ) then
#ifdef BAND
          open (iutrad,file=ffout,status='replace',form='unformatted',  &
              & recl=iym2*jx*ibyte,access='direct')
#else
          open (iutrad,file=ffout,status='replace',form='unformatted',  &
              & recl=iym2*jxm2*ibyte,access='direct')
#endif
          nrcrad = 0
        else if ( iotyp.eq.2 ) then
          open (iutrad,file=ffout,status='replace',form='unformatted')
        else
        end if
        call gradsrad(trim(ffout))
        print * , 'OPENING NEW RAD FILE: ',trim(ffout)
      end if
 
      if ( ichem.eq.1 ) then
        if ( ifchem ) then
          close (iutchem)
          call outname('CHE',idatex)
          if ( iotyp.eq.1 ) then
#ifdef BAND
            open (iutchem,file=ffout,status='replace',access='direct',  &
                & form='unformatted',recl=iym2*jx*ibyte)
#else
            open (iutchem,file=ffout,status='replace',access='direct',  &
                & form='unformatted',recl=iym2*jxm2*ibyte)
#endif
            nrcchem = 0
          else if ( iotyp.eq.2 ) then
            open (iutchem,file=ffout,status='replace',                  &
               &  form='unformatted')
          else
          end if
          call gradschem(trim(ffout))
          print * , 'OPENING NEW CHEM FILE: ',trim(ffout)
        end if
      end if

      if ( lakemod.eq.1 ) then
        close (iutlak)
        call outname('LAK',idatex)
        open (iutlak,file=ffout,status='replace',form='unformatted')
        print * , 'OPENING NEW LAK FILE: ',trim(ffout)
      end if
 
      if ( jyear.eq.jyear0 .and. ktau.eq.0 ) then
        nmo = (idate2/1000000-idate0/1000000)                           &
            & *12 + (idate2/10000-(idate2/1000000)*100)                 &
            & - (idate0/10000-(idate0/1000000)*100)
        nmo = min(nmo,nfmax)
        idatepp = idate0
        iyr(1) = idatepp/1000000
        imo(1) = (idatepp-iyr(1)*1000000)/10000
        idate1d(1) = idatepp
        do n = 2 , nmo
          idatepp = (idatepp/10000)*10000 + 10100
          iyr(n) = idatepp/1000000
          imo(n) = (idatepp-iyr(n)*1000000)/10000
          if ( imo(n).gt.12 ) then
            iyr(n) = iyr(n) + 1
            imo(n) = 1
            idatepp = iyr(n)*1000000 + imo(n)*10000 + 100
          end if
          idate1d(n) = idatepp
        end do
        if ( nmo.le.1 ) then
          nmo = 1
          idate1d(2) = idate0
        end if
 
      end if
 
      end subroutine mkfile

      end module mod_iunits
