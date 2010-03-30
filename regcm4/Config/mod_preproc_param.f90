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

      module mod_preproc_param
      implicit none
!
! PARAMETER definitions
!

!####################### GRID Geolocation ##############################

! Projection
!
! One in : 'LAMCON', Lambert conformal
!          'POLSTR', Polar stereographic
!          'NORMER', Normal  Mercator (ROTMER w/ plat = clat
!          'ROTMER', Rotated Mercator

      character(6) , parameter :: iproj = 'LAMCON'

! Grid point horizontal resolution in km

      real(4) , parameter :: ds = 60.0

! Pressure of model top in cbar

      real(4) , parameter :: ptop = 5.0

! Central latitude  of model domain in degrees, north hem. is positive

      real(4) , parameter :: clat = 45.39

! Central longitude of model domain in degrees, west is negative

      real(4) , parameter :: clon = 13.48

! Pole latitude (only for rotated Mercator Proj, else set = clat)

      real(4) , parameter :: plat = clat

! Pole longitude (only for rotated Mercator Proj, else set = clon)

      real(4) , parameter :: plon = clon

! Lambert true latitude (low latitude side)

      real(4) , parameter :: truelatl = 30.

! Lambert true latitude (high latitude side)

      real(4) , parameter :: truelath = 60.

! Surface minimum H2O percent to be considered water

      real(4) , parameter :: h2opct = 75.

! Resolution of the global terrain and landuse data be used
!
!     Use 60, for  1  degree resolution
!         30, for 30 minutes resolution
!         10, for 10 minutes resolution
!          5, for  5 minutes resolution
!          3, for  3 minutes resolution
!          2, for  2 minutes resolution

      integer , parameter :: ntypec = 10

! Same for subgrid (Used only if nsg > 1)

      integer , parameter :: ntypec_s = 10

! Interpolation Control flag.
!
!     true  -> Perform cressman-type objective analysis
!     false -> 16-point overlapping parabolic interpolation

      logical , parameter :: ifanal = .true.

! Smoothing Control flag
!
!     true  -> Perform extra smoothing in boundaries

      logical , parameter :: smthbdy = .false.

! Great Lakes levels adjustment Control Flag (Set true only if US EC)
!
!     true  -> Adjust Great Lakes Levels according to obs

      logical , parameter :: lakadj = .false.

! I/O format
!
!     1 => direct access binary
!     2 => netcdf

      integer , parameter :: itype_in  = 1
!     integer , parameter :: itype_out  = 1

! Fudging for landuse and texture for grid and subgrid

      logical , parameter :: fudge_lnd = .false.
      logical , parameter :: fudge_lnd_s = .false.
      logical , parameter :: fudge_tex = .false.
      logical , parameter :: fudge_tex_s = .false.

! Terrain output files

      character(50) , parameter :: filout = '../../Input/DOMAIN.INFO'
      character(50) , parameter :: filctl = '../../Input/DOMAIN.CTL'

! Global Begin and End date for Input Pre processing

      integer , parameter :: idate1 = 1990060100 ! BEGIN
      integer , parameter :: idate2 = 1990070100 ! END

! Type of Sea Surface Temperature used
!
! One in: GISST,OISST,OI2ST,OI_WK,OI2WK,FV_RF,FV_A2,FV_B2,EH5RF,
!         EH5A2,EH5B1,EHA1B,ERSST,ERSKT
!
      character(5) , parameter :: ssttyp = 'OI_WK'

! Number of Soil texture categories, leave it to 17

      integer , parameter :: ntex = 17

      end module mod_preproc_param
