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

      subroutine setup(nunit,unctl,iy,jx,ntypec,iproj,ds,clat,clon,     &
                     & igrads,ibyte,filout,filctl)
      use mod_interfaces
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
      intent (inout) :: filctl , filout
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
      open (nunit,file=filout,status='unknown',form='unformatted',      &
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
