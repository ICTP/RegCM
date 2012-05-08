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

!  integer , parameter :: nna = 99
!  integer , parameter :: nnb = 2010


integer , parameter :: nna = 1 
integer , parameter :: nnb = 2010

  real(dp) :: time         ! time step (sec)
  real(dp) :: xr(nna,nnb) ! concentration molec/cm3
!FAB add
  real(dp) :: xrin(nna,nnb), xrout(nna,nnb)  
  real(dp) :: xremit(nnb) ! emitted concentration molec/cm3
  real(dp) :: xrwdep(nnb) ! summed wet deposition

! BOXSET - box meteo+chem settings to be used for chemistry
!    HVSET - box meteo+chem  settings for HV rates

  real(dp) :: temp(nna)   ! Temp K
  real(dp) :: dens(nna)   ! Density molec cm-3
  real(dp) :: acqua(nna)  ! LWC grams/cm-3
  real(dp) :: rainfr(nna) ! rainout fraction during time step

  real(dp) :: xhour       ! Hour,  EST, decimal hrs
  real(dp) :: altmid(nna) ! altitude at midpoint, km
  real(dp) :: xlat(nna)   ! Latitude, decimal degrees
  real(dp) :: xlon(nna)   ! Longitude, decimal degrees
  integer :: idate        ! Date YYMMDD (YY=100 for 2000)

! UNIT INDICES

  integer :: iit , itr , itime , nhkem
!  integer :: isol                 ! input starting hour
  integer :: iplu , jcloud
  integer :: lin                  ! Basic input (fort.41)
  integer :: lcin                 ! Chem concentration input (fort.87)
  integer :: lout                 ! Basic output (fort.43)
  integer :: iprt , iprt2 , iprt3 ! I/O indices and flags
  integer :: nwrit , iwrit(nnb) , nhc(10)
  integer :: nchem , nchem1 , nchem2
  logical :: lprt
  integer :: kmax                 ! max for vector loop, if used in box
!    (vector dimension is in chemistry subroutine, always = 1 in box)

! CONCENTRATIONS AND EMISSIONS - note working XR in BASIC.
! xgen:       Emission rate molec/cm2-sec

  real(dp) :: dr(24,nnb)     ! stored concentration array
  real(dp) :: oxr(nna,nnb)   ! Prior concentration
  real(dp) :: gen(nnb)       !  Emissions
  real(dp) :: xgen(nnb)      !  Specific Emissions
  real(dp) :: updr(nnb)      !  Upwind concentration
  real(dp) :: topdr(nnb)     ! Concentin entrainment layer
  real(dp) :: upxr(nnb)      !  Upwind concentration
  real(dp) :: topxr(nnb)     ! Concentin entrainment layer
  real(dp) :: depo(nnb)      ! Dry dep velocity cm/s
  real(dp) :: xdepo(nnb)     ! Dry dep velocity cm/s
  character(len=4) :: ttchem(nnb)    ! Chemistry species name
!   real(dp) :: rrec(3,nnb)    ! Record reaction rates
  real(dp) :: xh2o           ! H2O molec/cm3
  real(dp) :: xnacl          ! NaCl concentration M/lit

! DYNAMICS AND TIMING SETUP
  real(dp) :: sol , dilu , alt1 , altm , tem , dtem , tlap
  real(dp) :: fhc(5) , xhm , xhrsun
  real(dp) :: wind , wx , side , wmax , wmin
  real(dp) :: alt      ! Layer thickness, cm
  real(dp) :: zdzdt    ! Entrainment parameter
  real(dp) :: zmax     ! Maximum layer height, cm         
  real(dp) :: zmin     ! Minimum layer height, cm
  real(dp) :: altbase
  real(dp) :: frz(24) , fgen(24) , frw(24) , fngen(24) , frs(24)
  real(dp) :: frgen , frngen

! GENERAL VARIABLES USED IN BOX MODEL PROGRAMS (not in common)
!FAB AGIN local variable shopuld notbe declared here !!
!  real(dp) :: xhm2 , fr1 , fr0 , frx    ! bsetr fraction setters
!  real(dp) :: sxh , ct , temps , windd , frwin , fn1 , fn0
!  real(dp) :: conc
!  real(dp) :: xn          ! brplac time  calc
!  real(dp) :: tint        ! bprint time calc
!  real(dp) :: xh1 , xh2   ! bprint time calc
!  real(dp) :: gas         ! bprint gas concentration
!  real(dp) :: ffgen       ! boxdyn parameter
!

!!$  integer :: ihour , ihour2 , ix1 , ix0 , ix21 , ix20  ! bsetr hour integers
!!$  integer :: ind , j1 , jj2 , ij2 , ixn , n            ! brplac indices
!!$  integer :: ni , ic1 , ic2 , ip                       ! bsplac indices
!!$  integer :: j2 , nx2 , iw , i , j                     ! bprint indices


!!$  integer :: kk , ic                                   ! Standard counters
!!$




  real(dp) :: zenith                          ! zenith angle from RegCM3
!
  ! O3 column DU
  real(dp) , parameter :: dobson(1) = 347.8D0
  ! Cloud above/below - no longer used
  real(dp) , parameter :: cldindx(1) = 1.00D0
  ! Cloud optical depth - no longer used
  real(dp) , parameter :: depth(1) = 0.0D0
  ! Aerosol optical depth
  real(dp) , parameter :: aaerx(1) = 0.38D0
  ! SO2 column, Dobson units
  real(dp) , parameter :: so2col = 0.1D0
  ! NO2 column, Dobson units
  real(dp) , parameter :: no2col = 0.1D0
  ! Aerosol single scattering albedo
  real(dp) , parameter :: aerssa = 0.75D0
  real(dp) , parameter :: albedo = 0.1D0   ! Albedo= (0.1)
  real(dp) , parameter :: saersa(1)= 0.0D0 ! sulf aerosol surf area cm2/cm3
  ! Droplet radius, cm. (normally .001)
  real(dp) , parameter :: droplet(1) = 1.000D-03
  ! Gas-aerosol transfer rate
  real(dp) , parameter :: rgasaq(1) = -100.0D+00
!
  logical , parameter :: lstsaq = .false.  ! future flag for steady state aq
  logical , parameter :: lexpo = .false.   ! flag for expo. decay solution
  real(dp) :: deptha     ! Cloud above optical depth
  real(dp) :: depthb     ! Cloud below optical depth
  real(dp) :: altabove   ! Cloud above weighted altitude
  real(dp) :: altbelow   ! Cloud below weighted altitude
!
end module mod_cbmz_boxvars
