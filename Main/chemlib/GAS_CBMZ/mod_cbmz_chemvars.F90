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

module mod_cbmz_chemvars
!
  use mod_intkinds
  use mod_realkinds
  use mod_cbmz_chemmech
!
  public
!
! chemvars.EXT    April, 2007
!
!  for RADICAL BALANCE-BACK EULER solver for chemistry (quadchem)
!      (chemmain.f cheminit.f chemrates chemsolve.f, linslv.f, jval2.f)
!
! NOTE- chemvars.EXT and chemmech.EXT must always go together
!       with chemmech.EXT first - it contains indices. (change 7/07)
!
!  Divided into ESSENTIAL input/output 
!    and OPTIONAL output which is rarely used.
!
!  CRITICAL PARAMETERS AND CHOICES:  - moved to chemmech 7/07
!
!  Other INCLUDE files: 
!   chemmech.EXT =  variables for the chem. mechanism
!     that are set once and do not change as the mechanism is called.
!   chemlocal.EXT = Variables that are shared among the subroutines
!     in the chemistry solver, but not used elsewhere and not saved.
!
! Program written by Sanford Sillman
! History:
!  12/06. Program written by Sandy Sillman based on commq7.EXT
!   7/07  Modified to put critical indices in chemmech.
! -------------------------------------------------------------

! ESSENTIAL index parameters:   moved to chemmech.EXT 7/07

! READ/WRITE AND VECTOR INDICES - moved to chemmech.EXT.
!
! BASIC INPUTS:  CONCENTRATIONS, DEPOSITION ARRAYS, TIME STEP.
! 
! c_xcin :
! c_xcin:       input  species  concentration, molec/cm3
!               including emissions during the time step
! c_xcout:      final species concentrations, molec/cm3
!               for gas species: sum of gas + linked aqueous,
!               in gas units (molec/cm3).
!               for aqueous species: M/liter.  
! c_xcav:       average species concentration during the time step
!               (to be used for wet deposition: 
!               xcwdep(M/cm2) = xcav(M/L)*rainfr*alt(cm)*.01 (dm2/cm2)
! c_xcemit:     emissions during time step, molec/cm3
!               (only used in expo. decay solution, usually zero)
! c_xcwdep:     Wet deposition in solver - Not used.
! c_time:       time step (sec)
!         
   real(rk8) :: c_xcin(c_kvec,c_cdim)   ! concentration molec/cm3
   real(rk8) :: c_xcout(c_kvec,c_cdim)  ! concentration, mol/cm3
   real(rk8) :: c_xcav(c_kvec,c_cdim)   ! concentration, mol/cm3
   real(rk8) :: c_xcemit(c_kvec,c_cdim) ! emissions, molec/cm3
   real(rk8) :: c_time                  ! time step (sec)
   real(rk8) :: c_jval(c_kvec,56)
! 
! (note = separate out TIME?)

!  CHEMISTRY INPUT PARAMETERS (CPARAMS):  
!         Parameters used to calculate rate constants

! c_h2oliq             Liquid Water Content (LWC) grams cm-3 (acqua)
! c_dens               DENSITY molec cm-3
! c_rainfr             Rainout fraction for time step
! c_temp               TEMPERATURE K
! c_saersa             SULFATE AEROSOL SURFACE AREA cm2 cm-3
! c_DROPLET             Droplet radius for gas-aq calc, cm (.001)
!                        (not yet vectorized - hard set at .001)
! 
! FUTURE, NOT CURRENTLY USED:  
! c_rgasaq(kvec,nr)     Gas-aq transfer rate, s-1 (optional, <0 cancels)
! c_lgasaq(kvec)        Flag to calculate (T) gas-aqueous transfer rate
! c_lstsaq               Flag for steady state gas/aq partitioning
!                        (T for steady state)
! c_lexpo(kvec)         Flag for modified backward Euler solution
!                        with exponential decay
! 
   real(rk8) :: c_temp(c_kvec)          ! temperature K
   real(rk8) :: c_dens(c_kvec)          ! density molec cm-3
   real(rk8) :: c_h2oliq(c_kvec)        ! LWC grams cm-3
   real(rk8) :: c_rainfr(c_kvec)        ! rainout fraction
   real(rk8) :: c_saersa(c_kvec)        ! s aerosol surf area cm2/cm3
   real(rk8) :: c_h2ogas(c_kvec)        ! H2O absolute humidity molec/cm3

   real(rk8) :: c_rgasaq(c_kvec,c_rdim) ! Input gas->aq rate
   real(rk8) :: c_DROPLET(c_kvec)       ! droplet radius, cm
   logical :: c_lgasaq(c_kvec)         ! Flag to calc. gas-aq rate
   logical :: c_lstsaq(c_kvec)         ! Flag for steady state gas-aq
   logical :: c_lexpo(c_kvec)          ! Flag for expo decay solution

!  PHOTOLYSIS INPUT PARAMETERS (HVPARAMS):  
!      Parameters used in calculating HV RATES.
!      Note: see also DATA FOR HV PARAMETERIZATION (chemmech)
!            and RETURNED J-VALUE ARRAY (jval, local in hvrates.f)
!
! c_jparam(   22)      Input parameters for j-value calculation
!                      (vectorize in future?)
!    c_jparam( 1)=zenith angle, degrees (calculated in program)
!    c_jparam( 2)=altitude (KM)  (OPTION - kPa)
!    c_jparam( 3)=ozone column (DU)
!    c_jparam( 4)=SO2   column (DU)  (0.1=1 ppb, SO2 0-1km)
!    c_jparam( 5)=NO2   column (DU)  (0.1=1 ppb, NO2 0-1km)
! 
!    c_jparam( 6)=aerosol optical depth
!                        (0.36 = Elterman 1968, normal, 0.76=polluted)
!    c_jparam( 7)=surface albedo   (0.1)
!    c_jparam( 8)= aerosol single scattering albedo
!                             (.75-.99, .75 Logan, .99 Madronich)
!    c_jparam( 9)= cloud optical depth for cloud ABOVE current alt.
!    c_jparam(10)= cloud optical depth for cloud BELOW current alt.
!
!    c_jparam(11)= opt-depth-weighted altitude for clouds above (KM)
!    c_jparam(12)= opt-depth-weighted altitude for clouds below (KM)
!    c_jparam(13)= temperature (K)
! 
!    c_jparam(20)= c_IDATE yymmddd (2000=100) or  day number (1-366)
!
!    c_lat (kk)     Latitude, degrees
!    c_lon (kk)     Longitude, degrees
!    c_hour         Hour at end of simulated time interval (decimal hrs)
!    (c_time also used - time interval - above.)
!    c_IDATE        Date YYMMDD (YY=>100) 
!                   Note, ZEN1 has alt.for YYYYDDD
!        (lat, lon, hour and date used in SOLAR ZENITH ANGLE  calc.)
!
  real(rk8) :: c_lat(c_kvec) ! Latitude, degrees
  real(rk8) :: c_lon(c_kvec) ! Longitude, degrees
  real(rk8) :: c_hour        ! Hour at end , EST (hrs)
  integer(ik4) :: c_idate        ! Date YYMMDD (YY=100 for 2000)
  real(rk8) :: c_jparam(22)  ! J-value input parameters

! OPTIONAL CHEM OUTPUT VALUES:    
!  May be used for subsequent analysis but not necessary.
!  OPTION: Move this from chemvars to chemlocal
!
!  c_rp(c_kvec,c_cdim)  Rate of production of species,molec/cm3/timestep
!  c_rl(c_kvec,c_cdim)  Rate of loss  of species ic, molec/cm3/timestep
!  c_rr(c_kvec,c_rdim)    Rate of reaction, molec/cm3/timestep

  real(rk8) :: c_rp(c_kvec,c_cdim)  ! Production, molec/cm3
  real(rk8) :: c_rl(c_kvec,c_cdim)  ! Loss, molec/cm3
  real(rk8) :: c_rr(c_kvec,c_rdim)  ! Reaction rate m/cm3

! OPTIONAL OUTPUTS RELATING TO NUMERICS: 
!   Currently in chemlocal:  
!    xohtest, xnotest, fohtest, final iter, history, geomavg

  real(rk8) :: c_ohtest   ! test: dOH/OH or dOH/HO2
  real(rk8) :: c_notest   ! test: dNO2/NO2
  integer(ik4) :: c_iter      ! chem. number of iterations

end module mod_cbmz_chemvars
