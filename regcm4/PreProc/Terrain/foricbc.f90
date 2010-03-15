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

      subroutine foricbc(ix,jx,kz,nsg,idate1,idate2,ibyte,dattyp,ssttyp,&
                       & ehso4,lsmtyp,aertyp)
      implicit none
!
! Dummy arguments
!
      character(7) :: aertyp
      character(5) :: dattyp , ssttyp
      logical :: ehso4
      integer :: ibyte , idate1 , idate2 , ix , jx , kz , nsg
      character(4) :: lsmtyp
      intent (in) aertyp , dattyp , ehso4 , ibyte , idate1 , idate2 ,   &
                & ix , jx , kz , lsmtyp , nsg , ssttyp
!
! Local variables
!
      integer :: ix_o , jx_o , kz_o
!
      open (23,file='../ICBC/icbc.param')
      write (23,'(a)') '      INTEGER JX'
      write (23,'(a)') '      INTEGER IX'
      write (23,'(a)') '      INTEGER KZ'
      write (23,'(a)') '      INTEGER NSG'
      write (23,'(a)') '      INTEGER JX_O'
      write (23,'(a)') '      INTEGER IX_O'
      write (23,'(a)') '      INTEGER KZ_O'
      write (23,'(a)') '      INTEGER NP'
      write (23,'(a)') '      INTEGER IDATE1'
      write (23,'(a)') '      INTEGER IDATE2'
      write (23,'(a)') '      INTEGER IBYTE'
      write (23,'(a)') '      CHARACTER*5 DATTYP'
      write (23,'(a)') '      CHARACTER*5 SSTTYP'
      write (23,'(a)') '      LOGICAL     EHSO4 '
      write (23,'(a)') '      CHARACTER*4 LSMTYP'
      write (23,'(a)') '      CHARACTER*7 AERTYP'
      write (23,99001) 'JX     =' , jx
      write (23,99001) 'IX     =' , ix
      write (23,99001) 'KZ     =' , kz
      write (23,99001) 'NSG    =' , nsg
      if ( dattyp=='FNEST' ) then
        write (*,*) 'You are preparing the NEST run from previous run'
        write (*,*) 'Please input the original size of ARRAYs'
        write (*,*) 'JX_O, IX_O, KZ_O'
        write (*,*) 'You need cut 2 rows for both JX_O and IX_O'
        read (*,*) jx_o , ix_o , kz_o
        write (23,99001) 'JX_O   =' , jx_o
        write (23,99001) 'IX_O   =' , ix_o
        write (23,99001) 'KZ_O   =' , kz_o
        write (23,99001) 'NP     =' , 15
      else
        write (23,99001) 'JX_O   =' , 1
        write (23,99001) 'IX_O   =' , 1
        write (23,99001) 'KZ_O   =' , 14
        write (23,99001) 'NP     =' , 15
      end if
      write (23,99002) 'IDATE1 =' , idate1
      write (23,99002) 'IDATE2 =' , idate2
      write (23,99002) 'IBYTE  =' , ibyte
      if ( dattyp=='ECMWF' ) then
        write (23,'(a)') "      parameter(DATTYP='ECMWF')"
      else if ( dattyp=='ERA40' ) then
        write (23,'(a)') "      parameter(DATTYP='ERA40')"
      else if ( dattyp=='ERAIN' .or. dattyp=='EIN15' ) then
        write (23,'(a)') "      parameter(DATTYP='EIN15')"
      else if ( dattyp=='EIN75' ) then
        write (23,'(a)') "      parameter(DATTYP='EIN75')"
      else if ( dattyp=='EIN25' ) then
        write (23,'(a)') "      parameter(DATTYP='EIN25')"
      else if ( dattyp=='ERAHI' ) then
        write (23,'(a)') "      parameter(DATTYP='ERAHI')"
      else if ( dattyp=='NNRP1' ) then
        write (23,'(a)') "      parameter(DATTYP='NNRP1')"
      else if ( dattyp=='NNRP2' ) then
        write (23,'(a)') "      parameter(DATTYP='NNRP2')"
      else if ( dattyp=='NRP2W' ) then
        write (23,'(a)') "      parameter(DATTYP='NRP2W')"
      else if ( dattyp=='GFS11' ) then
        write (23,'(a)') "      parameter(DATTYP='GFS11')"
      else if ( dattyp=='EH5OM' ) then
        write (23,'(a)') "      parameter(DATTYP='EH5OM')"
      else if ( dattyp=='FVGCM' ) then
        write (23,'(a)') "      parameter(DATTYP='FVGCM')"
      else if ( dattyp=='FNEST' ) then
        write (23,'(a)') "      parameter(DATTYP='FNEST')"
      else
        print * , 'BOUNDARY CONDITION DATASET DOES NOT EXIST'
        stop 'subroutine foricbc'
      end if
      if ( ssttyp=='OISST' .or. ssttyp=='OI_NC' ) then
        write (23,'(a)') "      parameter(SSTTYP='OISST')"
      else if ( ssttyp=='OI2ST' ) then
        write (23,'(a)') "      parameter(SSTTYP='OI2ST')"
      else if ( ssttyp=='OI_WK' ) then
        write (23,'(a)') "      parameter(SSTTYP='OI_WK')"
      else if ( ssttyp=='OI2WK' ) then
        write (23,'(a)') "      parameter(SSTTYP='OI2WK')"
      else if ( ssttyp=='GISST' ) then
        write (23,'(a)') "      parameter(SSTTYP='GISST')"
      else if ( ssttyp=='FV_RF' ) then
        write (23,'(a)') "      parameter(SSTTYP='FV_RF')"
      else if ( ssttyp=='FV_A2' ) then
        write (23,'(a)') "      parameter(SSTTYP='FV_A2')"
      else if ( ssttyp=='FV_B2' ) then
        write (23,'(a)') "      parameter(SSTTYP='FV_B2')"
      else if ( ssttyp=='EH5RF' ) then
        write (23,'(a)') "      parameter(SSTTYP='EH5RF')"
      else if ( ssttyp=='EH5A2' ) then
        write (23,'(a)') "      parameter(SSTTYP='EH5A2')"
      else if ( ssttyp=='EH5B1' ) then
        write (23,'(a)') "      parameter(SSTTYP='EH5B1')"
      else if ( ssttyp=='EHA1B' ) then
        write (23,'(a)') "      parameter(SSTTYP='EHA1B')"
      else if ( ssttyp=='ERSST' ) then
        write (23,'(a)') "      parameter(SSTTYP='ERSST')"
      else if ( ssttyp=='ERSKT' ) then
        write (23,'(a)') "      parameter(SSTTYP='ERSKT')"
      else
        print * , 'SST DATASET DOES NOT EXIST'
        stop 'subroutine foricbc'
      end if
      if ( dattyp=='EH5OM' .and. ehso4 ) then
        write (23,'(a)') "      parameter(EHSO4 =.true. )"
      else
        write (23,'(a)') "      parameter(EHSO4 =.false.)"
      end if
      if ( lsmtyp=='BATS' ) then
        write (23,'(a)') "      parameter(LSMTYP='BATS')"
      else if ( lsmtyp=='USGS' ) then
        write (23,'(a)') "      parameter(LSMTYP='USGS')"
      else
        print * , 'LANDUSE LEGEND DOES NOT EXIST'
        stop 'subroutine foricbc'
      end if
      if ( aertyp=='AER00D0' ) then
        write (23,'(a)') "      parameter(AERTYP='AER00D0')"
      else if ( aertyp=='AER01D0' ) then
        write (23,'(a)') "      parameter(AERTYP='AER01D0')"
      else if ( aertyp=='AER10D0' ) then
        write (23,'(a)') "      parameter(AERTYP='AER10D0')"
      else if ( aertyp=='AER11D0' ) then
        write (23,'(a)') "      parameter(AERTYP='AER11D0')"
      else if ( aertyp=='AER00D1' ) then
        write (23,'(a)') "      parameter(AERTYP='AER00D1')"
      else if ( aertyp=='AER01D1' ) then
        write (23,'(a)') "      parameter(AERTYP='AER01D1')"
      else if ( aertyp=='AER10D1' ) then
        write (23,'(a)') "      parameter(AERTYP='AER10D1')"
      else if ( aertyp=='AER11D1' ) then
        write (23,'(a)') "      parameter(AERTYP='AER11D1')"
      else
        print * , 'AEROSOL TYPE DOES NOT EXIST'
        stop 'subroutine foricbc'
      end if
      close (23)
      open (23,file='../ICBC/icbc.x')
      write (23,'(a)') '#!/bin/csh -f'
      write (23,'(a)') 'make clean'
      write (23,'(a)') 'foreach FILE (RCM_SST.ctl RCM_SST.dat SST.RCM)'
      write (23,'(a)') 'if ( -f $FILE ) /bin/rm $FILE'
      write (23,'(a)') 'end'
      if ( dattyp=='FVGCM' .or.                                         &
         & (dattyp=='FNEST' .and. (ssttyp=='FV_RF' .or. ssttyp=='FV_A2')&
         & ) ) then
        write (23,'(a)') 'make SST_FVGCM'
        write (23,'(a)') './SST_FVGCM'
        write (23,'(a)') '/bin/rm -f SST_FVGCM*.o SST_FVGCM'
      else if ( dattyp=='EH5OM' .or.                                    &
               &(dattyp=='FNEST' .and. (ssttyp=='EH5RF' .or.            &
               &ssttyp=='EH5A2' .or. ssttyp=='EH5B1' .or.               &
               &ssttyp=='EHA1B')) ) then
        write (23,'(a)') 'make SST_EH5OM'
        write (23,'(a)') './SST_EH5OM'
        write (23,'(a)') '/bin/rm -f SST_EH5OM*.o SST_EH5OM'
      else if ( (dattyp=='EIN15' .or. dattyp=='ERAIN') .and.            &
               &(ssttyp=='ERSST' .or. ssttyp=='ERSKT') ) then
        write (23,'(a)') 'make SST_ERSST'
        write (23,'(a)') './SST_ERSST'
        write (23,'(a)') '/bin/rm -f SST_ERSST*.o SST_ERSST'
      else
        write (23,'(a)') 'make SST'
        write (23,'(a)') './SST'
        write (23,'(a)') '/bin/rm -f SST_1DEG*.o SST'
      end if
      if ( aertyp(4:5)/='00' ) then
        write (23,'(a)') 'make AEROSOL'
        write (23,'(a)') './AEROSOL'
        write (23,'(a)') '/bin/rm -f AEROSOL.o AEROSOL'
      end if
      if ( dattyp=='ERAHI' ) then
        write (23,'(a)') 'cp ../../Commons/tools/srcERAHI/*.f .'
        write (23,'(a)') 'make ERAHI_HT'
        write (23,'(a)') './ERAHI_HT'
        write (23,'(a)') 'make ERAHI_PS'
        write (23,'(a)') './ERAHI_PS'
        write (23,'(a)') 'make ERAHI_T'
        write (23,'(a)') './ERAHI_T'
        write (23,'(a)') 'make ERAHI_Q'
        write (23,'(a)') './ERAHI_Q'
        write (23,'(a)') 'make ERAHI_U'
        write (23,'(a)') './ERAHI_U'
        write (23,'(a)') 'make ERAHI_V'
        write (23,'(a)') './ERAHI_V'
      end if
      write (23,'(a)') 'make ICBC'
      write (23,'(a)') './ICBC'
      write (23,'(a)') '/bin/rm -f ICBC.o ICBC SST.RCM'
      close (23)
      if ( dattyp=='ERAHI' ) then
        open (23,file='../ICBC/ERAHI.param')
        write (23,'(a)') '      INTEGER IDATE1'
        write (23,'(a)') '      INTEGER IDATE2'
        write (23,'(a)') '      INTEGER IBYTE'
        write (23,99002) 'IDATE1 =' , idate1
        write (23,99002) 'IDATE2 =' , idate2
        write (23,99002) 'IBYTE  =' , ibyte
        close (23)
      end if
99001 format ('      parameter(',a8,i4,')')
99002 format ('      parameter(',a8,i10,')')
 
      end subroutine foricbc
