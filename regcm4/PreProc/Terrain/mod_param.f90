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

      module mod_param
      implicit none
!
! PARAMETER definitions
!
      character(6) , parameter :: iproj = 'LAMCON'
      integer , parameter :: iy = 34
      integer , parameter :: jx = 48
      integer , parameter :: kz = 18
      integer , parameter :: nsg = 1
      real(4) , parameter :: ds = 60.0
      real(4) , parameter :: ptop = 5.0
      real(4) , parameter :: clat = 45.39
      real(4) , parameter :: clon = 13.48
      real(4) , parameter :: plat = clat
      real(4) , parameter :: plon = clon
      real(4) , parameter :: truelatl = 30.
      real(4) , parameter :: truelath = 60.
      real(4) , parameter :: h2opct = 75.
      integer , parameter :: ntypec = 10
      integer , parameter :: ntypec_s = 10
      logical , parameter :: ifanal = .true.
      logical , parameter :: smthbdy = .false.
      logical , parameter :: lakadj = .false.
!
!     itype_in is I/O format
!
!     1 => direct access binary
!     2 => netcdf
      integer , parameter :: itype_in  = 1
!     integer , parameter :: itype_out  = 1

      integer , parameter :: igrads = 1
      integer , parameter :: ibigend = 1
      integer , parameter :: ibyte = 4
      logical , parameter :: fudge_lnd = .false.
      logical , parameter :: fudge_lnd_s = .false.
      logical , parameter :: fudge_tex = .false.
      logical , parameter :: fudge_tex_s = .false.
      character(50) , parameter :: filout = '../../Input/DOMAIN.INFO'
      character(50) , parameter :: filctl = '../../Input/DOMAIN.CTL'
      integer , parameter :: idate1 = 1990060100
      integer , parameter :: idate2 = 1990070100
!
! One in:
!        ECMWF,ERA40,ERAIN,EIN75,EIN15,EIM25,ERAHI,NNRP1,NNRP2,
!        NRP2W,GFS11,FVGCM,FNEST,EH5OM
!
      character(5) , parameter :: dattyp = 'EIN15'

! One in:
!        GISST,OISST,OI2ST,OI_WK,OI2WK,FV_RF,FV_A2,FV_B2,EH5RF,
!        EH5A2,EH5B1,EHA1B,ERSST,ERSKT
!
      character(5) , parameter :: ssttyp = 'OI_WK'

      logical , parameter :: ehso4 = .false.
      character(4) , parameter :: lsmtyp = 'BATS'
      integer , parameter :: nveg = 20
      character(7) , parameter :: aertyp = 'AER00D0'
      integer , parameter :: ntex = 17
      integer , parameter :: nproc = 16
      end module mod_param
