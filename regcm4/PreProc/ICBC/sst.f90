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

      program sst

      use mod_dynparam

      implicit none
      integer :: ierr
      character(256) :: namelistfile , prgname
!
!     Read input global namelist
!
      call getarg(0, prgname)
      call getarg(1, namelistfile)
      call initparam(namelistfile, ierr)
      if ( ierr/=0 ) then
        write ( 6, * ) 'Parameter initialization not completed'
        write ( 6, * ) 'Usage : '
        write ( 6, * ) '          ', trim(prgname), ' regcm.in'
        write ( 6, * ) ' '
        write ( 6, * ) 'Check argument and namelist syntax'
        stop
      end if

      if ( ssttyp=='GISST' .or. ssttyp=='OISST' .or.                    &
       &   ssttyp=='OI_NC' .or. ssttyp=='OI2ST' .or.                    &
       &   ssttyp=='OI_WK' .or. ssttyp=='OI2WK' ) then
        call sst_1deg
      else if ( ssttyp=='EH5RF' .or. ssttyp=='EH5A2' .or.               &
       &        ssttyp=='EH5B1' .or. ssttyp=='EHA1B' ) then
        call sst_eh5om
      else if ( ssttyp=='ERSST' .or. ssttyp=='ERSKT' ) then
        call sst_ersst
      else if ( ssttyp=='FV_RF' .or. ssttyp=='FV_A2' .or.               &
       &        ssttyp=='FV_B2' ) then
        call sst_fvgcm
      else
        print *, 'Unknown SSTTYP ', ssttyp , ' specified in ',          &
              & trim(namelistfile)
        stop
      end if

      print *, 'Successfully generated SST'

      end program sst
