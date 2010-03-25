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

      program report
      use mod_regcm_param , only : ntr ,  nbin
      use mod_emission
      implicit none
!
! Local variables
!
      character(8) , dimension(ntr) :: chtrname
      real , dimension(nbin) :: dustbsiz
!
!     This program produce a report have all information about
!     emission species and transported species
!
      namelist /species/ ele_retroa , ele_retrob , ele_poeta ,          &
      & ele_poetb , ele_gfed , ele_edgara , ele_edgarb , ele_mozrta ,   &
      & ele_mozrtb
      namelist /transport/ chtrname
      namelist /dust_bins/ dustbsiz

      open (10,file='namelist.report',status='old')
      read (10,nml=species)
      read (10,nml=transport)
      read (10,nml=dust_bins)
 
      write (*,*) ele_retroa
      write (*,*) ele_edgara
      write (*,*) chtrname
      write (*,*) dustbsiz

      end program report
