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

  real(rk8) :: time         ! time step (sec)
  real(rk8) :: xr(nna,nnb) ! concentration molec/cm3
!FAB add
  real(rk8) :: xrin(nna,nnb), xrout(nna,nnb)  
  real(rk8) :: xremit(nnb) ! emitted concentration molec/cm3
  real(rk8) :: xrwdep(nnb) ! summed wet deposition

! BOXSET - box meteo+chem settings to be used for chemistry
!    HVSET - box meteo+chem  settings for HV rates

  real(rk8) :: temp(nna)   ! Temp K
  real(rk8) :: dens(nna)   ! Density molec cm-3
  real(rk8) :: acqua(nna)  ! LWC grams/cm-3
  real(rk8) :: rainfr(nna) ! rainout fraction during time step

  real(rk8) :: xhour       ! Hour,  EST, decimal hrs
  real(rk8) :: altmid(nna) ! altitude at midpoint, km
  real(rk8) :: xlat(nna)   ! Latitude, decimal degrees
  real(rk8) :: xlon(nna)   ! Longitude, decimal degrees
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

  real(rk8) :: dr(24,nnb)     ! stored concentration array
  real(rk8) :: oxr(nna,nnb)   ! Prior concentration
  real(rk8) :: gen(nnb)       !  Emissions
  real(rk8) :: xgen(nnb)      !  Specific Emissions
  real(rk8) :: updr(nnb)      !  Upwind concentration
  real(rk8) :: topdr(nnb)     ! Concentin entrainment layer
  real(rk8) :: upxr(nnb)      !  Upwind concentration
  real(rk8) :: topxr(nnb)     ! Concentin entrainment layer
  real(rk8) :: depo(nnb)      ! Dry dep velocity cm/s
  real(rk8) :: xdepo(nnb)     ! Dry dep velocity cm/s
  character(len=4) :: ttchem(nnb)    ! Chemistry species name
!   real(rk8) :: rrec(3,nnb)    ! Record reaction rates
  real(rk8) :: xh2o           ! H2O molec/cm3
  real(rk8) :: xnacl          ! NaCl concentration M/lit

! DYNAMICS AND TIMING SETUP
  real(rk8) :: sol , dilu , alt1 , altm , tem , dtem , tlap
  real(rk8) :: fhc(5) , xhm , xhrsun
  real(rk8) :: wind , wx , side , wmax , wmin
  real(rk8) :: alt      ! Layer thickness, cm
  real(rk8) :: zdzdt    ! Entrainment parameter
  real(rk8) :: zmax     ! Maximum layer height, cm         
  real(rk8) :: zmin     ! Minimum layer height, cm
  real(rk8) :: altbase
  real(rk8) :: frz(24) , fgen(24) , frw(24) , fngen(24) , frs(24)
  real(rk8) :: frgen , frngen

! GENERAL VARIABLES USED IN BOX MODEL PROGRAMS (not in common)
!FAB AGIN local variable shopuld notbe declared here !!
!  real(rk8) :: xhm2 , fr1 , fr0 , frx    ! bsetr fraction setters
!  real(rk8) :: sxh , ct , temps , windd , frwin , fn1 , fn0
!  real(rk8) :: conc
!  real(rk8) :: xn          ! brplac time  calc
!  real(rk8) :: tint        ! bprint time calc
!  real(rk8) :: xh1 , xh2   ! bprint time calc
!  real(rk8) :: gas         ! bprint gas concentration
!  real(rk8) :: ffgen       ! boxdyn parameter
!

!!$  integer(ik4) :: ihour , ihour2 , ix1 , ix0 , ix21 , ix20  ! bsetr hour integers
!!$  integer(ik4) :: ind , j1 , jj2 , ij2 , ixn , n            ! brplac indices
!!$  integer(ik4) :: ni , ic1 , ic2 , ip                       ! bsplac indices
!!$  integer(ik4) :: j2 , nx2 , iw , i , j                     ! bprint indices


!!$  integer(ik4) :: kk , ic                                   ! Standard counters
!!$




  real(rk8) :: zenith                          ! zenith angle from RegCM3
!
  ! O3 column DU
  real(rk8) , parameter :: dobson(1) = 347.8D0
  ! Cloud above/below - no longer used
  real(rk8) , parameter :: cldindx(1) = 1.00D0
  ! Cloud optical depth - no longer used
  real(rk8) , parameter :: depth(1) = 0.0D0
  ! Aerosol optical depth
  real(rk8) , parameter :: aaerx(1) = 0.38D0
  ! SO2 column, Dobson units
  real(rk8) , parameter :: so2col = 0.1D0
  ! NO2 column, Dobson units
  real(rk8) , parameter :: no2col = 0.1D0
  ! Aerosol single scattering albedo
  real(rk8) , parameter :: aerssa = 0.75D0
  real(rk8) , parameter :: albedo = 0.1D0   ! Albedo= (0.1)
  real(rk8) , parameter :: saersa(1)= 0.0D0 ! sulf aerosol surf area cm2/cm3
  ! Droplet radius, cm. (normally .001)
  real(rk8) , parameter :: droplet(1) = 1.000D-03
  ! Gas-aerosol transfer rate
  real(rk8) , parameter :: rgasaq(1) = -100.0D+00
!
  logical , parameter :: lstsaq = .false.  ! future flag for steady state aq
  logical , parameter :: lexpo = .false.   ! flag for expo. decay solution
  real(rk8) :: deptha     ! Cloud above optical depth
  real(rk8) :: depthb     ! Cloud below optical depth
  real(rk8) :: altabove   ! Cloud above weighted altitude
  real(rk8) :: altbelow   ! Cloud below weighted altitude
!
end module mod_cbmz_boxvars
