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

module mod_cbmz_boxvars
!
  use mod_intkinds
  use mod_realkinds
!
  public
!
! boxvars.EXT      December, 2006
!
!   INCLUDE and COMMON file for the BOX MODEL (boxmain.f, boxproc1.f)
!    with the RADICAL BALANCE-BACK EULER solver for chemistry (quadchem)
!      (chemmain.f cheminit.f chemrates chemsolve.f, linslv.f, jval2.f)
!
!  Variables for the CHEMISTRY SOLVER are not contained here.
!
! Program written by Sanford Sillman
! History:
!  12/06. Program written by Sandy Sillman based on commq7.EXT
! -------------------------------------------------------------

!    BASIC PARAMETERS - also needed by chemistry
!
!    time:     time step (sec)
!    xr:       working species concentrations, molec/cm3
!    (xremit:  emitted component of xr (=direct emissions or entrainment))
!    (xrwdep:  wet deposition sum )
!

!FAB : reduce a bi these dimensions to save memory

!  integer(ik4) , parameter :: nna = 99
!  integer(ik4) , parameter :: nnb = 2010


integer(ik4) , parameter :: nna = 1
integer(ik4) , parameter :: nnb = 2010

  real(rkx) :: time         ! time step (sec)
  real(rkx) :: xr(nna,nnb) ! concentration molec/cm3
!FAB add
  real(rkx) :: xrin(nna,nnb), xrout(nna,nnb)
  real(rkx) :: xremit(nnb) ! emitted concentration molec/cm3
  real(rkx) :: xrwdep(nnb) ! summed wet deposition

! BOXSET - box meteo+chem settings to be used for chemistry
!    HVSET - box meteo+chem  settings for HV rates

  real(rkx) :: temp(nna)   ! Temp K
  real(rkx) :: dens(nna)   ! Density molec cm-3
  real(rkx) :: acqua(nna)  ! LWC grams/cm-3
  real(rkx) :: rainfr(nna) ! rainout fraction during time step

  real(rkx) :: xhour       ! Hour,  EST, decimal hrs
  real(rkx) :: altmid(nna) ! altitude at midpoint, km
  real(rkx) :: xlat(nna)   ! Latitude, decimal degrees
  real(rkx) :: xlon(nna)   ! Longitude, decimal degrees
  integer(ik4) :: idate        ! Date YYMMDD (YY=100 for 2000)

! UNIT INDICES

  integer(ik4) :: iit , itr , itime , nhkem
!  integer(ik4) :: isol                 ! input starting hour
  integer(ik4) :: iplu , jcloud
  integer(ik4) :: lin                  ! Basic input (fort.41)
  integer(ik4) :: lcin                 ! Chem concentration input (fort.87)
  integer(ik4) :: lout                 ! Basic output (fort.43)
  integer(ik4) :: iprt , iprt2 , iprt3 ! I/O indices and flags
  integer(ik4) :: nwrit , iwrit(nnb) , nhc(10)
  integer(ik4) :: nchem , nchem1 , nchem2
  logical :: lprt
  integer(ik4) :: kmax                 ! max for vector loop, if used in box
!    (vector dimension is in chemistry subroutine, always = 1 in box)

! CONCENTRATIONS AND EMISSIONS - note working XR in BASIC.
! xgen:       Emission rate molec/cm2-sec

  real(rkx) :: dr(24,nnb)     ! stored concentration array
  real(rkx) :: oxr(nna,nnb)   ! Prior concentration
  real(rkx) :: gen(nnb)       !  Emissions
  real(rkx) :: xgen(nnb)      !  Specific Emissions
  real(rkx) :: updr(nnb)      !  Upwind concentration
  real(rkx) :: topdr(nnb)     ! Concentin entrainment layer
  real(rkx) :: upxr(nnb)      !  Upwind concentration
  real(rkx) :: topxr(nnb)     ! Concentin entrainment layer
  real(rkx) :: depo(nnb)      ! Dry dep velocity cm/s
  real(rkx) :: xdepo(nnb)     ! Dry dep velocity cm/s
  character(len=4) :: ttchem(nnb)    ! Chemistry species name
!   real(rkx) :: rrec(3,nnb)    ! Record reaction rates
  real(rkx) :: xh2o           ! H2O molec/cm3
  real(rkx) :: xnacl          ! NaCl concentration M/lit

! DYNAMICS AND TIMING SETUP
  real(rkx) :: sol , dilu , alt1 , altm , tem , dtem , tlap
  real(rkx) :: fhc(5) , xhm , xhrsun
  real(rkx) :: wind , wx , side , wmax , wmin
  real(rkx) :: alt      ! Layer thickness, cm
  real(rkx) :: zdzdt    ! Entrainment parameter
  real(rkx) :: zmax     ! Maximum layer height, cm
  real(rkx) :: zmin     ! Minimum layer height, cm
  real(rkx) :: altbase
  real(rkx) :: frz(24) , fgen(24) , frw(24) , fngen(24) , frs(24)
  real(rkx) :: frgen , frngen

! GENERAL VARIABLES USED IN BOX MODEL PROGRAMS (not in common)
!FAB AGIN local variable shopuld notbe declared here !!
!  real(rkx) :: xhm2 , fr1 , fr0 , frx    ! bsetr fraction setters
!  real(rkx) :: sxh , ct , temps , windd , frwin , fn1 , fn0
!  real(rkx) :: conc
!  real(rkx) :: xn          ! brplac time  calc
!  real(rkx) :: tint        ! bprint time calc
!  real(rkx) :: xh1 , xh2   ! bprint time calc
!  real(rkx) :: gas         ! bprint gas concentration
!  real(rkx) :: ffgen       ! boxdyn parameter
!

!!$  integer(ik4) :: ihour , ihour2 , ix1 , ix0 , ix21 , ix20  ! bsetr hour integers
!!$  integer(ik4) :: ind , j1 , jj2 , ij2 , ixn , n            ! brplac indices
!!$  integer(ik4) :: ni , ic1 , ic2 , ip                       ! bsplac indices
!!$  integer(ik4) :: j2 , nx2 , iw , i , j                     ! bprint indices


!!$  integer(ik4) :: kk , ic                                   ! Standard counters
!!$




  real(rkx) :: zenith                          ! zenith angle from RegCM3
!
  ! O3 column DU
  real(rkx) , parameter :: dobson(1) = 347.8_rkx
  ! Cloud above/below - no longer used
  real(rkx) , parameter :: cldindx(1) = 1.00_rkx
  ! Cloud optical depth - no longer used
  real(rkx) , parameter :: depth(1) = 0.0_rkx
  ! Aerosol optical depth
  real(rkx) , parameter :: aaerx(1) = 0.38_rkx
  ! SO2 column, Dobson units
  real(rkx) , parameter :: so2col = 0.1_rkx
  ! NO2 column, Dobson units
  real(rkx) , parameter :: no2col = 0.1_rkx
  ! Aerosol single scattering albedo
  real(rkx) , parameter :: aerssa = 0.75_rkx
  real(rkx) , parameter :: albedo = 0.1_rkx   ! Albedo= (0.1)
  real(rkx) , parameter :: saersa(1)= 0.0_rkx ! sulf aerosol surf area cm2/cm3
  ! Droplet radius, cm. (normally .001)
  real(rkx) , parameter :: droplet(1) = 1.000e-3_rkx
  ! Gas-aerosol transfer rate
  real(rkx) , parameter :: rgasaq(1) = -100.0_rkx
!
  logical , parameter :: lstsaq = .false.  ! future flag for steady state aq
  logical , parameter :: lexpo = .false.   ! flag for expo. decay solution
  real(rkx) :: deptha     ! Cloud above optical depth
  real(rkx) :: depthb     ! Cloud below optical depth
  real(rkx) :: altabove   ! Cloud above weighted altitude
  real(rkx) :: altbelow   ! Cloud below weighted altitude
!
end module mod_cbmz_boxvars
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
