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
!
      module mod_regcm_param

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: ipunit = 255
!
!################### GRID DIMENSION ####################################
!

! Point in Y (latitude) direction

      integer :: iy

! Point in X (longitude) direction

      integer :: jx

! Point in vertical

      integer :: kz

! Projection
!
! One in : 'LAMCON', Lambert conformal
!          'POLSTR', Polar stereographic
!          'NORMER', Normal  Mercator (ROTMER w/ plat = clat
!          'ROTMER', Rotated Mercator

      character(6) :: iproj

! Grid point horizontal resolution in km

      real(4) :: ds

! Pressure of model top in cbar

      real(4) :: ptop

! Central latitude  of model domain in degrees, north hem. is positive

      real(4) :: clat

! Central longitude of model domain in degrees, west is negative

      real(4) :: clon

! Pole latitude (only for rotated Mercator Proj, else set = clat)

      real(4) :: plat

! Pole longitude (only for rotated Mercator Proj, else set = clon)

      real(4) :: plon

! Lambert true latitude (low latitude side)

      real(4) :: truelatl

! Lambert true latitude (high latitude side)

      real(4) :: truelath

!###################### I/O control flag ###############################

! Create GrADS CTL files

      integer :: igrads

! Machine endianess. LEAVE IT UNTOUCHED IF WANT TO EXCHANGE I/O FILES

      integer :: ibigend

! Number of bytes in reclen. Usually 4

      integer :: ibyte

!####################### MPI parameters ################################

! Number of processor used

      integer :: nproc

!####################### MPI parameters ################################

! Sub grid decomposition

      integer :: nsg

! Set amount of printout (still unused, sorry)

      integer :: debug_level

! Buffer Zone Depth
! nspgx-1,nspgd-1 represent the number of cross/dot point slices
! on the boundary sponge or relaxation boundary conditions.

      integer :: nspgx
      integer :: nspgd

! Number od split exp modes

      integer :: nsplit

! Number of lake points for lake model

      integer :: lkpts

! Type of global analysis datasets used in Pre processing
!
! One in: ECMWF,ERA40,ERAIN,EIN75,EIN15,EIM25,ERAHI,NNRP1,NNRP2,
!         NRP2W,GFS11,FVGCM,FNEST,EH5OM
!

      character(5) :: dattyp

! Type of Sea Surface Temperature used
!
! One in: GISST,OISST,OI2ST,OI_WK,OI2WK,FV_RF,FV_A2,FV_B2,EH5RF,
!         EH5A2,EH5B1,EHA1B,ERSST,ERSKT
!
      character(5) :: ssttyp

! SO4 Control Flag

      logical :: ehso4

! Land Surface Legend type
!
! One in : BATS,USGS

      character(4) :: lsmtyp
      integer :: nveg

! Aerosol dataset used
!
! One in : AER00D0 -> Neither aerosol, nor dust used
!          AER01D0 -> Biomass, SO2 + BC + OC, no dust
!          AER10D0 -> Anthropogenic, SO2 + BC + OC, no dust
!          AER11D0 -> Anthropogenic+Biomass, SO2 + BC + OC, no dust
!          AER00D1 -> No aerosol, with dust
!          AER01D1 -> Biomass, SO2 + BC + OC, with dust
!          AER10D1 -> Anthropogenic, SO2 + BC + OC, with dust
!          AER11D1 -> Anthropogenic+Biomass, SO2 + BC + OC, with dust

      character(7) :: aertyp

! Tracer parameters: number of tracers and bins number for dust

      integer :: ntr
      integer :: nbin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of configureation. Below this point things are
!    calculated from above or should be considered as fixed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: iym1
      integer :: iym2
      integer :: iym3
      integer :: jxp1
      integer :: jxm1
      integer :: jxm2
      integer :: kzm1
      integer :: kzm2
      integer :: kzp1
      integer :: kzp2
      integer :: kzp3
      integer :: kzp4
      integer :: iysg
      integer :: jxsg
      integer :: iym1sg
      integer :: jxm1sg
      integer :: iym2sg
      integer :: jxm2sg
      integer :: nnsg
      integer :: nspgv
      integer :: nspgp

#ifdef MPP1
      integer :: myid
      integer :: jxp
      integer :: jxpsg
      integer :: iwest , ieast , isouth , inorth
      integer :: jbegin , ibegin
      integer :: jendl , iendl
      integer :: jendx , iendx
      integer :: jendm , iendm
#endif

! Surface minimum H2O percent to be considered water

      real(4) :: h2opct

! Resolution of the global terrain and landuse data be used
!
!     Use 60, for  1  degree resolution
!         30, for 30 minutes resolution
!         10, for 10 minutes resolution
!          5, for  5 minutes resolution
!          3, for  3 minutes resolution
!          2, for  2 minutes resolution

      integer :: ntypec

! Same for subgrid (Used only if nsg > 1)

      integer :: ntypec_s

! Interpolation Control flag.
!
!     true  -> Perform cressman-type objective analysis
!     false -> 16-point overlapping parabolic interpolation

      logical :: ifanal

! Smoothing Control flag
!
!     true  -> Perform extra smoothing in boundaries

      logical :: smthbdy

! Great Lakes levels adjustment Control Flag (Set true only if US EC)
!
!     true  -> Adjust Great Lakes Levels according to obs

      logical :: lakadj

! I/O format
!
!     1 => direct access binary
!     2 => netcdf

      integer :: itype_in
!     integer :: itype_out

! Fudging for landuse and texture for grid and subgrid

      logical :: fudge_lnd
      logical :: fudge_lnd_s
      logical :: fudge_tex
      logical :: fudge_tex_s

! Number of Soil texture categories, leave it to 17

      integer :: ntex

! Terrain output files

      character(64) :: terfilout
      character(64) :: terfilctl

! Global Begin and End date for Input Pre processing

      integer :: globidate1 ! BEGIN
      integer :: globidate2 ! END

! Fixed dimensions

      integer , dimension(289276) :: mdatez

      integer , parameter :: numbat = 21 + 6
      integer , parameter :: numsub = 16
      integer , parameter :: nrad2d = 22
      integer , parameter :: nrad3d = 5

      contains

      subroutine initparam(filename)
        implicit none
        character (len=*) , intent(in) :: filename

        namelist /geoparam/ iproj , ds , ptop , clat , clon , plat ,    &
                     & plon , truelatl, truelath
        namelist /terrainparam/ itype_in , ntypec , ntypec_s , ifanal , &
                     & smthbdy , lakadj , fudge_lnd , fudge_lnd_s ,     &
                     & fudge_tex , fudge_tex_s , ntex , h2opct ,        &
                     & terfilout , terfilctl
        namelist /dimparam/ iy , jx , kz , nsg , nproc
        namelist /ioparam/ igrads , ibigend , ibyte
        namelist /debugparam/ debug_level
        namelist /boundaryparam/ nspgx , nspgd
        namelist /modesparam/ nsplit
        namelist /lakemodparam/ lkpts
        namelist /globdatparam/ dattyp , ssttyp , ehso4 , globidate1 ,  &
                     & globidate2
        namelist /lsmparam/ lsmtyp
        namelist /aerosolparam/ aertyp , ntr, nbin

        open(ipunit, file=filename, status='old', &
                     action='read', err=100)
!
        read(ipunit, geoparam, err=100)
        read(ipunit, terrainparam, err=100)
        read(ipunit, dimparam, err=100)

!       Setup all convenience dimensions

        iym1 = iy - 1
        iym2 = iy - 2
        iym3 = iy - 3
        jxp1 = jx + 1
        jxm1 = jx - 1
        jxm2 = jx - 2
        kzm1 = kz - 1
        kzm2 = kz - 2
        kzp1 = kz + 1
        kzp2 = kz + 2
        kzp3 = kz + 3
        kzp4 = kz + 4
        iysg = iy * nsg
        jxsg = jx * nsg
        iym1sg = (iy-1) * nsg
        jxm1sg = (jx-1) * nsg
        iym2sg = (iy-2) * nsg
        jxm2sg = (jx-2) * nsg
        nnsg = nsg*nsg

        read(ipunit, ioparam, err=100)
        read(ipunit, debugparam, err=100)
        read(ipunit, boundaryparam, err=100)

        nspgv = (nspgd+nspgx)*8 + 8
        nspgp = nspgx*4

        read(ipunit, modesparam, err=100)
        read(ipunit, lakemodparam, err=100)
        read(ipunit, globdatparam, err=100)
        read(ipunit, lsmparam, err=100)

        if (lsmtyp == 'BATS') then
          nveg = 20
        else if (lsmtyp == 'USGS') then
          nveg = 25
        else
          write ( 6 , * ) 'Unknown LSM data type. Use BATS or USGS'
          stop
        end if

        read(ipunit, aerosolparam, err=100)

        return

  100   write ( 6, * ) 'Cannot read namelist file'
        stop

      end subroutine

      end module mod_regcm_param
