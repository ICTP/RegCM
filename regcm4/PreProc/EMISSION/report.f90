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

      program report
      use mod_dynparam
      implicit none

      integer :: ierr
      character(256) :: namelistfile, prgname
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

      call printspecies(ntr, nbin)

      end program report

      subroutine printspecies(ntr, nbin)
      use mod_emission
      implicit none
!
! Dummy arguments
!
      integer , intent(in) :: ntr , nbin
!
! PARAMETERS
!
      integer , parameter :: maxtr = 20
      integer , parameter :: maxbin = 10
!
! Local variables
!
      character(8) , dimension(maxtr) :: chtrname
      real , dimension(maxbin) :: dustbsiz
!
!     This program produce a report have all information about
!     emission species and transported species
!
      namelist /species/ ele_retroa , ele_retrob , ele_poeta ,          &
      & ele_poetb , ele_gfed , ele_edgara , ele_edgarb , ele_mozrta ,   &
      & ele_mozrtb
      namelist /transport/ chtrname
      namelist /dust_bins/ dustbsiz

      chtrname = '       '
      dustbsiz = 0.0
      open (10,file='namelist.report',status='old')
      read (10,nml=species)
      read (10,nml=transport)
      read (10,nml=dust_bins)
 
      write (*,*) ele_retroa
      write (*,*) ele_edgara
      write (*,*) chtrname(1:ntr)
      write (*,*) dustbsiz(1:nbin)
      end subroutine printspecies
