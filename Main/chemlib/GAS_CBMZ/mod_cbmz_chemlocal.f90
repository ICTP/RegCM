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

module mod_cbmz_chemlocal
!
  use mod_realkinds
  use mod_cbmz_chemmech
!
  public
!
! chemlocal.EXT    April, 2007
!     for RADICAL BALANCE-BACK EULER solver for chemistry (quadchem)
!      (chemmain.f cheminit.f chemrates chemsolve.f, linslv.f, jval2.f)
!
!   for SOLVER VARIABLES 
!   which are shared among the solver subroutines in common blocks
!   but are not saved or passed to the main program.  

!  It also contains declarations for LOCAL VARIABLES
!   which are used in solver subroutines but not retained 
!   or passed to other subroutines  
!
!  Other INCLUDE files: 
!   chemvars.EXT =  input-output variables 
!    that are passed to and from the solver each time it is called.
!   chemmech.EXT =  variables for the chem. mechanism
!     that are set once and do not change as the mechanism is called.
!
! Program written by Sanford Sillman
! History:
!  12/06. Program written by Sandy Sillman based on commq7.EXT
! -------------------------------------------------------------
! LOCAL SPECIES CONCENTRATIONS 
!     (see also input-output variables in chemvars.EXT)
!
! xc:         average species concentration over time step, molec/cm3
!                 (internally converted to mol/lit-atm for aqueous)
!                 (note: in expo decay solution, 
!              this is distinguished from final species concentration)
!
! xcfinr:   Ratio of FINAL (end of time step) to AVERAGE concentration
!             (=1 for back-Euler. Used for expo decay solution.
!
! xclastq : OPTION, xc after last aquasolve,  used in next aquasolve
! 
  real(dp) ::  xc( c_kvec,c_cdim) ! concentration molec/cm3
  real(dp) ::  xcfinr( c_kvec,c_cdim) ! ratio FINAL/AVG cncn
  real(dp) ::  xclastq( c_kvec,c_cdim) ! conc molec/cm3

! CHEMICAL REACTION RATES
!     (see also rate parameters in chemmech.EXT
!       and optional output in chemvars.EXT)

!  ratek     Rate constants (s-1 or cm3 s-1 or equiv. aqueous)
!    (c_rk     Parameters for rate constant calculation in chemmech)
!    (c_rr     Calculated rate in solver: ratek*C1*C2 in chemvars)
!  rateh     Henry's law coefficients Moles/liter-atm
!    (c_rkh    Parameters for Henry's law coefficients in chemmech)
!  rateq     Aqueous equilibrium coefficients Moles/liter
!    (c_rkq    Parameters for aqueous eq. coefficients in chemmech)
!  rateqq    Special equilibrium coefficients Moles/liter
!    (c_rkqq   Parameters for Special eq. coefficients in chemmech)
!  rhdif     Henry's law coeff. modified by droplet diffusion
!  egasaq    Gas-to-aqueous transfer coefficient, s-1
!              (Varies for reaction based on accom, molwt)
!              (This is either input or calculated internally)
!   (c_gasaq   Input gas-aq transfer coefficient in chemvars)
! 2009 addition
! ratero2(kk,nr,2) = rate constants for parameterized RO2-RO2: 1 self 2 cross R
!
  real(dp) :: ratek(c_kvec,c_rdim)     ! Rate constants
  real(dp) :: rateh(c_kvec,c_rdim)     ! Henry's law constants
  real(dp) :: rateq(c_kvec,c_rdim)     ! Aqueous equil. constants
  real(dp) :: rateqq(c_kvec,c_rdim)    ! Special equil. constants
  real(dp) :: rhdif(c_kvec,c_rdim)     ! H. modified by diffusion
  real(dp) :: egasaq(c_kvec,c_rdim)    ! Gas-to-aq transfer s-1
  real(dp) :: ratero2(c_kvec,c_rdim,2) ! Rate constants for RO2-RO2

! LOCAL CHEMICAL SOLVER VARIABLES:  
! PRODUCTION AND LOSS RATES AND INTERIM CONCENTRATIONS
! SINGLE SPECIES AND LINKED PAIR SOLUTIONS
!  (OPTION - iteration counter is in common for written output)
!
! cpm(c_kvec,is1,is2)    Cross production is1->is2 for multi-group
!                    conversion from species 1 to species 2, molec/cm3
!                         (is1, is2 are indices w/in multi group)
! cpro(c_kvec,is1,is2) Cross-production for pair species:
!              conversion from species 1 to species 2, molec/cm3/tstep
!
! geomavg(kk, ic)   Geometric average adjustment factor 
!                      for difficult convergence.
!                    Solution is geom. avg of current and prior iter,
!                      using geomavg factor, adjusted based on history.
!                      (1= current iter; 0.5=pure geometric avg)
! history(kk, ic, iter) Solution for concentration of species ic
!                         in each iteration, molec/cm3
!                         (used in geometic weighting adjustment)
!
! rlm(c_kvec,c_cdim)     Summed loss for multi group, molec/cm3/timestep
!                         (including inter-group conversion)
! rlmulti(c_kvec       ) Combined loss rate for multi group of species
!                         (Omitting inter-group conversion)
! rloss(c_kvec,c_cdim)  Solver loss rate + final (nonstst), molec/cm3
! rloss1(c_kvec,c_cdim)  Interim rloss without cpro pairsol. adjustment
! rlpair(c_kvec,c_cdim)  Net loss rate for pair group`(mol/cm3/timestep)
!                          weighted by pairfac to conserve pair mass
!
! rpm(c_kvec,c_cdim)     Summed production for all of multi group
!                         (including inter-group conversion)
! rpmulti(c_kvec       ) Combined production rate for multi group`
!                         (Omitting inter-group conversion)
! rppair(c_kvec,c_cdim)  Net production rate for pair group`
!                          weighted by pairfac to conserve pair mass
! rpro(c_kvec,c_cdim)   Solver production rate+initial concentration, 
!                                           molec/cm3/timestep
! rpro1(c_kvec,c_cdim)  Interim rpro w/out cpro pair production
!
! rrp(c_kvec,c_cdim)   Continuous production sum thru cascade, molec/cm3
! rrl(c_kvec,c_cdim)   Continuous loss sum thru cascade,  molec/cm3
! 
! xrm(c_kvec,c_cdim)     Initial concentration for each member of  
!                                      multi-species group, molec/cm3
!                           (used in chemsolve and noxsolve)
! xrp(c_kvec, c_cdim)   Solver initial concentration for single species
!    or for each species in group for paired solution, molec/cm3
! xrppair(c_kvec, c_cdim)  Solver initial concentration for pair group
!    sum, weighted by pairfac to conserve pair mass.  molec/cm3
!       (used in chemsolve only, not nec. in common)
! xrr(c_kvec, c_cdim)   Solver interim solution      for single species
!    or for each species in group for paired solution, molec/cm3
! xrrm(c_kvec,c_cdim)    Initial solution      for each member of  
!               multi-species group (before normalization), molec/cm3
! xrrpair(c_kvec, c_cdim)  Solver solution for pair group sum
!              (weigted by pairfac to conserve mass). molec/cm3
!               used in normalization of indiv. species solutions
!       (used in chemsolve only, not nec. in common)
!
! 
  real(dp) :: cpm(c_kvec,c_cdim,c_cdim)   ! Cross-prod for multi
  real(dp) :: cpro(c_kvec,c_cdim,c_cdim)  ! Cross-production 
! 
  real(dp) :: history(c_kvec,c_cdim,400)  ! Cncn for each iter
  real(dp) :: geomavg(c_kvec,c_cdim)      ! Geom. avg factor
!
  real(dp) :: rlm(c_kvec,c_cdim)          ! Loss for multi-spec 
  real(dp) :: rlmulti(c_kvec)             ! Loss group sum for multi
  real(dp) :: rloss(c_kvec,c_cdim)        ! Solver loss rate
  real(dp) :: rloss1(c_kvec,c_cdim)       ! Interim rloss w/out cpro
  real(dp) :: rlpair(c_kvec,c_cdim)       ! Loss rate for pair grp
  real(dp) :: rpm(c_kvec,c_cdim)          ! Production rate for multi
  real(dp) :: rpmulti(c_kvec)             ! Product'n multi group sum
  real(dp) :: rppair(c_kvec,c_cdim)       ! Production for pair grp
  real(dp) :: rpro(c_kvec,c_cdim)         ! Solver production rate
  real(dp) :: rpro1(c_kvec,c_cdim)        ! Interim rpro w/out cpro 
  real(dp) :: rrp(c_kvec,c_cdim)          ! Interim production rate 
  real(dp) :: rrl(c_kvec,c_cdim)          ! Interim loss rate 
! 
  real(dp) :: xrm(c_kvec,c_cdim)          ! Initial conc. for multi
  real(dp) :: xrp(c_kvec,c_cdim)          ! Solver prior concentration
  real(dp) :: xrppair(c_kvec)             ! Solver prior pair group cn.
  real(dp) :: xrr(c_kvec,c_cdim)          ! Solver concentration
  real(dp) :: xrrm(c_kvec,c_cdim)         ! Initial solution for multi
  real(dp) :: xrrpair(c_kvec)             ! Solution for pair group sum

! OPTION - iteration counter is in common for output, may cut from common.
!          -> moved to chemvars as c_iter
!     integer iter                     ! Iteration counter

! LCHEMRATES,  LREAC  
! VARIABLES ASSOCIATED WITH ODD HYDROGEN RADICALS AND ODD NITROGEN
!
! foh(kk)           OH/HO2 fraction. (ohsolv and presolve. not common)
!
! oddhdel(kk,2)     Summed sensitivity of net production of Hx to OH,
!                    summed as PHx* (d ln (pHx)/d ln([OH]))  (1=w RO2 2=w/o RO2)
!                    (molec/cm3)
!                    (This results in a linearized update of Hx: 
!                     PHx = oddhsum + oddhdel(OH-OHp)) 
! oddhloh(kk,2)       Partial sum of net production/loss of Hx, molec/cm3
!                    only from reactions that remove OH
! oddhlho2(kk,2)       Partial sum of net production/loss of Hx, molec/cm3
!                    only from reactions that remove HO2
! oddhsrc(kk,2)       Sumed production of Hx (molec/cm3)  
!                                   from net production reactions only 
!                  based on prior species concentrations and reactions
! oddhsum(kk,2)       Sum of net production/loss of Hx, molec/cm3
!                  based on prior species concentrations and reactions
! 
! oddhsump(kk)      Sum of net production/loss of Hx from prior iteration
! oddhfacp(kk)      Hx concentration adjustment from prior iteration
!
! (oddnsum(kk) - deleted)
!
! senhcat(kk,icat)  Estimated ssensitivity of odd-H species conc. 
!                    to changes in OH, by species category.  
!                     d ln [odd-h spec]/d ln([OH])
!                     (Zero for non-odd-H species)
!                     (Sometimes used also for NOx)
! senshx(kk,2)       Sensitivity of Hx production to OH (w/ and w/o RO2)
!             for individual reaction: d ln (pHx)/d ln([OH]) 
!
! sourcnx(kk)     Summed NOx source from chemistry (molec/cm3/timestep)
!                   (includes net NOx source from PAN reactions)
! sinknx(kk)      Summed NOx sink   from chemistry (molec/cm3/timestep)
!                   (includes net NOx sink   from PAN reactions)
!
! xfohtest          Convergence test for OH/HO2 ratio (nonvectorized)
!                     =((OH/HO2)p-(OH/HO2))/(OH/HO2p
!
! THE FOLLOWING ARE MOVED TO chemvars.EXT: 
! xnotest          Convergence test ratio for NOx (nonvectorized)
!                     =(NOxp-NOx)/NOxp
! xohtest          Convergence test ratio for OH (nonvectorized)
!                     =(OHp-OH)/OHp

!  real(dp) :: foh(c_kvec)        ! OH/HO2.  not common

  real(dp) :: oddhdel(c_kvec,2)  ! Summed Hx sens to Oh, mol/cm3
  real(dp) :: oddhloh(c_kvec,2)  ! Hx loss from OH,  mol/cm3
  real(dp) :: oddhlho2(c_kvec,2) ! Hx loss from HO2,  mol/cm3
  real(dp) :: oddhsrc(c_kvec,2)  ! Summed P Hx, mol/cm3
  real(dp) :: oddhsum(c_kvec,2)  ! Summed net P/L Hx, mol/cm3
  real(dp) :: senshx(c_kvec,2)   ! Hx sens: d ln(pHx)/d ln([OH])
  real(dp) :: senhcat(c_kvec,99) ! Hx sens. by species categ.
!
  real(dp) :: sourcnx(c_kvec)    ! Summed NOx source, mol/cm3
  real(dp) :: sinknx(c_kvec)     ! Summed NOx sink, mol/cm3

  real(dp) :: xfohtest           ! OH/HO2 convergence test 

  real(dp) :: oddhfacp(c_kvec)   ! Factor for OH from prior iteration
  real(dp) :: oddhsump(c_kvec)   ! Net pro/loss of OH from prior it.
!
! SPECIES INDICES AND COUNTERS FOR CHEMISTRY SOLVER. 

! ncsol(nc)              Chem. species number for head of pair group
!                        being solved for.  
!                        NOT IN COMMON.  This is passed to chemsolve.
! 
! ncsolv(nc)        Species number for each species in solver call, 
!                    including subspecies in pair chains,   
!                    nc is number within the solver.
!                    ncsolv(nc) gives the species number (ic)
!                    (Local in chemsolve)
!
! nssolv(nc)        Pair/multisolve group number for each species 
!                    in solver call, including subspec in pair chains.
!                    nc is number within the solver.
!                    nssolv(nc) gives the pair group number (is) 
!                    (Local in chemsolve)
!
! nsol                  Number of pair groups to be solved for
!                        (Used in chemread, chemsolve, not in common)
!
! nsolv                 Number of species to be solved for 
!                       in individual call to solver (chemsolve)
!                        (Used in chemread, chemsolve, not in common)
! ncdim                 Maximum number of total species 
!               to be solved simultaneously in pair or multisolve
!                        (Chemread only. May be unnec.)
! nsdim                 Maximum number of species pair groups 
!               to be solved simultaneously in pair or multisolve
!                   (excluding pair chains).
!       (in quadchem, chemsolve, chemread. Must be in common for now.)
!

!  integer :: nsol           ! Number of pair groups being solved for
!  integer :: nsolv          ! Number of species being solved for
  integer :: nsdim          ! Maximum number of pair groups being solved
  integer :: ncdim          ! Maximum Number of species being solved for
!  integer :: ncsolv(c_cdim) ! Species number (ic) for species in solver
!  integer :: nssolv(c_cdim) ! Pair group number (is) for spec. in solver


!
! LOCAL VARIABLES USED IN CHEMISTRY SUBROUTINES 
!        (not passed between subroutines or included in common block)
! 
! prior (kk)   General vector vbl for prior species concentrations
!                                                    (molec/cm3)
!                 (used for initial concentration or prior estimate)
! alpha (kk)   General vector vbl 
! beta  (kk)   General vector vbl
! gamma (kk)   General vector vbl
!

! stoicx                             stoichiometry sum
! ionsum                              aquasolve ion counter
! ic, ic1, ic2,ic3, icc, ics, ics2, icc1, icc2, icc3, is, iss, iscs
!  ich, icq
!  icx, icx1, icy1, icx2, icy2, icp, icp1, icp2:   
! ica1, ica2, icb1, icb2, nra1,nra2,nrb1,nrb2 
!                                Chem species index numbers
!  neq                          Aqueous equilibrium counters
! nc, nc1, nc2, ncc, ncf, nn, nne           Chem. species counters
! icat1, icat2, icatp, icatp2     indices used for species categories
! nr, nr1, nr2, nrh, nrq, nrqq, nrx, np  :          Reaction counters
! kk                            :          Vectorization counter
! kw       :   Index for vector diagnostic output (see c_kkw)

!  real(dp) :: calpha(c_kvec) ! General vector variable
!  real(dp) :: cbeta(c_kvec)  ! General vector variable
!  real(dp) :: cgamma(c_kvec) ! General vector variable
!  real(dp) :: prior(c_kvec)  ! Prior species conc (molec/cm3)

!  real(dp) :: stoicx            ! stoichiometry sum


!!$!FAB : VERY DANGEROUS HERE, LOCAL COUNTER SHOULD NOT BE DEFINED in GLOBAL MODULE
!!$
!!$  ! Chem index
!!$  integer :: ic , ic1 , ic2 , ic3 , iic , icc , ics , ics2 , icc1 , icc2 , icc3
!!$  ! Chem local index
!!$  integer :: is , iss , iscs , isc2
!!$  ! Chem index
!!$  integer :: ich , icq , icx , icx1 , icx2 , icy1 , icy2 , icp , icp1 , icp2
!!$  ! chem index - pair and multi
!!$  integer :: icr1 , icr2 , icr , isr1 , isr2 , icpair
!!$  ! aquasolve ion counter
!!$  integer :: ionsum
!!$  ! Aquasolve chem index
!!$  integer :: ica1 , ica2 , icb1 , icb2 , nra1 , nra2 , nrb1 , nrb2
!!$  ! Aqueous counters
!!$  integer :: neq
!!$  ! Chem species counters
!!$  integer :: nc , nc1 , nc2 , ncc , ncf , nn , nne
!!$  ! Reaction counters
!!$  integer :: nr , nr1 , nr2 , nrh , nrq , nrqq , nrx , np
!!$  ! indices used for species categories
!!$  integer :: icat1 , icat2 , icatp , icatp2
!!$  ! Vectorization counters
!!$  integer :: kk , kw
!!$  ! General counters
!!$  integer :: i , j , k , ii , ij , iii , n
!!$! 
! CUT LOCAL VARIABLES  - Local in individual subroutines
!
! lloss        Local: Flag for identifying exchange loss reaction
!                         (cheminit)
! lpro         Local: Flag for identifying exchange production  reaction
!                         (cheminit)
!   hvrate(kk, 56)    Interpolated hv rates averaged over time interval.
!                                      (sec-1)  (bhvrates)
!   fhvint             Fraction assigned to each discrete time period
!                       used to interpolate for the full interval
!                       (=1/ihvint) (bhvrates)
!   ihvint             Number of discrete times  used  to interpolate
!                       hv over the simulated time interval  (bhvrates)
! 
!   jval(      56)       J-value output from jval2.f parameterization.
!                                      (sec-1)  (bhvrates)
!   zen(kk)           Zenith angle (degrees) (90=horizon) (bhvrates)
!
!  ncsol(c_cdim)      Species list passed in call to chemsolve

end module mod_cbmz_chemlocal
