!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_cbmz_chemmech

  use mod_intkinds
  use mod_realkinds

  public
!
! chemmech.EXT    April, 2007
!     for RADICAL BALANCE-BACK EULER solver for chemistry (quadchem)
!      (chemmain.f cheminit.f chemrates chemsolve.f, linslv.f, jval2.f)
!
! NOTE- chemvars.EXT and chemmech.EXT must always go together
!       with chemmech.EXT first - it contains indices. (changed 7/07)
!
!   for MECHANISM and NUMERICAL VARIABLES:  reaction and j-value data
!   which are read and/or set up  once at the start of the program
!   and then remain unchanged during calls to the solver.
!
! ----------------------------
!  CRITICAL PARAMETERS AND CHOICES:
!   c_kvec, parameter for VECTORIZATION.
!   c_cdim and c_rdim = maximum number of species and reactions
!                       (Currently 1511,3211. Can be increased.)
!
!   Real vs real(rkx):  (Double precision is standard)
!   If real(rkx), change 'real' to 'double precision'
!    throughout this file and all fortran files.
!    also change exponents (DO vs EO) in BOD and YTN
! ----------------------------
!
!  Other INCLUDE files:
!   chemvars.EXT =  input-output variables
!    that are passed to and from the solver each time it is called.
!   chemlocal.EXT = Variables that are shared among the subroutines
!     in the chemistry solver, but not used elsewhere and not saved.
!
! Program written by Sanford Sillman
! History:
!  12/06. Program written by Sandy Sillman based on commq7.EXT
!   7/07  Modified to put critical indices in chemmech.
! -------------------------------------------------------------

! ESSENTIAL index parameters:
!     c_kvec:  maximum index for vectorization
!     c_icdim: maximum dimension for the number of species
!     c_rdim: maximum dimension for the number of reactions

  ! Maximum dimension for vectorization
  integer(ik4), parameter :: c_kvec = 1

  ! max dimension for number of species
  integer(ik4), parameter :: c_cdim = 1511
  ! max dimension for number of reactions
  integer(ik4), parameter :: c_rdim = 3211
!
! READ/WRITE AND VECTOR INDICES
! c_kkw               Index for diagnostic write output (0 for none)
!                     (May be set to write for specific grid and time)
! c_kmax              Maximum number for vectorization
! c_out               Output file index
! c_rin               Input file index for REACTION.DAT (mechanism)
! c_hvin              Input file index for TUVGRID2 (hv data)
!
  integer(ik4) :: c_kkw ! Index for diagnostic write output (0 for none)
  integer(ik4) :: c_kmax    ! Maximum number for vectorization
!
  integer(ik4), parameter :: c_out = 27  ! Output unit
  integer(ik4), parameter :: c_rin = 25  ! Input unit for REACTION.DAT (mechanism)
  integer(ik4) :: c_hvin ! Input unit for TUVGRID2 (hv data)
!
! REACTION RATE PARAMETERS:
!
! c_rk      Gas-phase rate parameters (for calculating rate:  up to 7)
!            connected with rate calc. options in chemrates
! c_rkh     Henry's law rate parameters
! c_rkq     Aqueous equilibrium parameters (A=B+C, H+ and OH- reactions)
! c_rkqq    Special equilibrium parameters
!
!  c_nrk(c_rdim)    Format index for reaction rate
!                    This ID's read format and reaction rate formula
!  c_nrkh(c_rdim)   Format index for Henry's law coefficient
!  c_nrkq(c_rdim)   Format index for regular equilibrium constant
!  c_nrkqq(c_rdim)  Format index for special equilibrium constant
!
!  c_accom(c_rdim)   Accomodation coefficients
!  c_molwt(c_rdim)   Molecular weights

  real(rkx) :: c_rk(7,c_rdim)   ! gas-phase rate parameters
  real(rkx) :: c_rkh(2,c_rdim)  ! Henry's law parameters
  real(rkx) :: c_rkq(2,c_rdim)  ! Aqueous equilib. params.
  real(rkx) :: c_rkqq(2,c_rdim) ! Special equilib. params.
!
  integer(ik4) :: c_nrk(c_rdim)     ! Format index for reaction rate
  integer(ik4) :: c_nrkh(c_rdim)    ! Format index for Henry's law coef.
  integer(ik4) :: c_nrkq(c_rdim)    ! Format index for equilibrium constant
  integer(ik4) :: c_nrkqq(c_rdim)   ! Format index for special eq constant
!
  real(rkx) :: c_accom(c_rdim)  ! Accomodation coefficients
  real(rkx) :: c_molwt(c_rdim)  ! Molecular weights

! CHEMICAL MECHANISM AND SPECIES INDICES
!  c_reactant(c_rdim,2):   Index for reactant species in A+B=C..
!  c_product(c_rdim,20)    Index for product species in A+B=C...
!  c_nnpro(c_rdim)   : number of reaction products for reaction nr.
!  c_stoich(c_rdim,20)     Stoichiometry for each product
!  c_prodarr(nr,ic) Array of stoichiometry for production of chem ic
!                   from reaction nr - assigned to gas species.
!                   (used in solution for chem pairs and multisolve
!
!  c_henry(c_rdim,2)       Index for Henry's law species: A<=>B
!  c_aqueous(c_rdim,3)     Index for aqueous equil. species: A<=>B+C
!  c_aqspec(c_rdim,3)    Index for special equil. species: A<=>B+C
!
!  c_nreac     Number of chemical reactions
!  c_nreach    Number of Henry's law equilibrium constants
!  c_nreacq    Number of regular equilibrium constants (w/H+, OH-)
!  c_nreaqq    Number of special equilibrium constants
!
!  c_nchem1    Number if input/output (transported) species
!  c_nchem2    Number of total species
!
!  c_icat      Category of chemical species
!  c_lsts      Flag to identify steady state species (T for stst)
!  c_lump(nlump,3)  Lumped species: species n1= sum from n2 to n3
!  c_llump    Flag to identify lumped species accumulator (n1) (T for n1)
!
!  (nk() and iarray() - from old NR solver - cut.)

  integer(ik4) :: c_reactant(c_rdim,2)       ! Reactant species numbers
  integer(ik4) :: c_product(c_rdim,20)       ! Product  species numbers
  integer(ik4) :: c_nnpro(c_rdim)            ! Number of reaction products
  real(rkx) :: c_stoich(c_rdim,20)       ! Stoich. for products
  real(rkx) :: c_prodarr(c_rdim,c_cdim)  ! Stoich. nr for ic
  integer(ik4) :: c_henry(c_rdim,2)          ! Henry's law species numbers
  integer(ik4) :: c_aqueous(c_rdim,3)        ! Aqueous equil. species numbers
  integer(ik4) :: c_aqspec(c_rdim,3)         ! Special equil. species numbers

  integer(ik4) :: c_nreac    ! Number of chemical reactions
  integer(ik4) :: c_nreach   ! Number of Henry's law equilibrium constants
  integer(ik4) :: c_nreacq   ! Number of regular equilibrium constants
  integer(ik4) :: c_nreaqq   ! Number of special equilibrium constants
  integer(ik4) :: c_nchem1   ! Number if input/output (transported) species
  integer(ik4) :: c_nchem2   ! Number of total species
!
  integer(ik4) :: c_icat(c_cdim)   ! Category of chemical species
  logical :: c_lsts(c_cdim)   ! Flag to identify steady state species
  integer(ik4) :: c_lump(c_cdim,3) ! Lumped species:  n1= sum n2 to n3
  logical :: c_llump(c_cdim)  ! Flag to identify lumped species (n1)

! INDICES FOR CHEMICAL REACTION SOLUTION PROCEDURE (CASCADE)

!  c_cascade(ic,n)  Array of cascade species: identified species
!                   in order of solution and single/multisolve
!
! c_nppair(c_cdim,23)  pair pointers: for linked pair chain in cascade
! c_nppair(ic,1) = pointer to directly linked pair species
!                  (to species directly above in pair hierarchy)
! c_nppair(ic,2) = pointer to main primary species for pair group
! c_nppair(ic,3) = number of pair group subspecies for this species
! c_nppair(ic,4+) = chem. index for pair group subspecies of species ic
! c_npmulti(ic,1) = chem index for head of multi group for this spec.
!
! c_nnrchem(c_cdim) : Cascade solution order:
!                     number of reactions linked to species ic
! c_nrchem(c_cdim,c_rdim)  Cascade solution order:
!                     Reaction numbers (nr) linked to species ic
! c_nnrchp(c_cdim)   Cascade solution order:
!                     number of reactions linked to ic as  product
! c_nrchmp(c_cdim,c_rdim)  Cascade solution order:
!                     Reaction numbers (nr) linked to ic as product
! n_nicreac(c_rdim) : Chemical species associated with reaction.
!
! c_exloss(ic, is,i)  Reaction number for 'exchange' loss reaction:
!                    loss to species ic;  counter for number of
!                    exchange species is, number of reactions i
! c_expro(ic,is,i)   Reaction number for 'exchange' product reaction:
!                    produces species ic;  counter for number of
!                    exchange species is, number of reactions i
! c_exspec(ic,is)    Exchange species number for species #is
!                    that is in back-forth reaction with species ic
  integer(ik4) :: c_cascade(c_cdim,c_cdim)  ! Array to call cascade
  integer(ik4) :: c_nppair(c_cdim,23)       ! Pointers to paired species
  integer(ik4) :: c_npmulti(c_cdim,1)       ! Pointer to multi spec. head

  integer(ik4) :: c_nnrchem(c_cdim)         ! Number of reactions linked to spec.
  integer(ik4) :: c_nrchem(c_cdim, c_rdim)  ! Reaction numbers for spec.
  integer(ik4) :: c_nnrchp(c_cdim)          ! Number of product reactions for ic
  integer(ik4) :: c_nrchmp(c_cdim, c_rdim)  ! Productreaction numbers for ic
  integer(ik4) :: c_nicreac(c_rdim)         ! Species ic associated with reaction nr

  integer(ik4) :: c_exloss(c_cdim,20,5)   ! Exchange loss reaction number
  integer(ik4) :: c_expro(c_cdim,20,5)    ! Exchange product reaction number
  integer(ik4) :: c_exspec(c_cdim,20)     ! Exchange species number

!  CHEMISTRY REACTION COUNTERS:
!   c_stoiloss(c_rdim): counter for net loss of reactant per reaction
!                     (used to ID self-reaction or product=reactant)
!   stoicpro, stoicprx - cut, no longer used
!
!   c_oddhx(c_rdim,n):    counter for net prod/loss of odd-H per reaction
!                        n=1 with RO2 n=2 w/o RO2
!   c_pronox(c_rdim):   Counter for net production of NOx per reaction.
!                      (net production only; not loss)
! c_pairfac(c_cdim)  Pair group conservation  of mass factor per spec.
!                    Counts double for C=A+B, for A, B in pair group
! c_multfac(c_cdim)  Multi group conservation  of mass factor per spec.
!                    Counts double for C=A+B, for A, B in pair group
!
  real(rkx) :: c_stoiloss(c_rdim)  ! net reactant loss
  real(rkx) :: c_oddhx(c_rdim,2)   ! counter for net odd-H
  real(rkx) :: c_pronox(c_rdim)    ! counter for net P(NOx)

  real(rkx) :: c_pairfac(c_cdim)  ! Pair group net rp factor
  real(rkx) :: c_multfac(c_cdim)  ! Multi group net rp factor

! INDICES AND POINTERS FOR AQUEOUS EQUILIBRIA
! c_nequil(c_cdim)   Number of aqueous/eq species linked to gas ic
! c_npequil(c_cdim)  Pointer to linked gas-master species
! c_ncequil(c_cdim,6)  List of aqueous/eq species for each gas-master
! c_nrequil(c_cdim,6)  Henry's law number (nrh) for first aqueous spec,
!                      Aqueous equil. number (nrq) for subsequent spec,
!                      for species associated with each gas-master.
! c_ion(c_cdim)       Charge associated with each species
!
  integer(ik4) :: c_nequil(c_cdim)    ! Number of aq species linked to gas ic
  integer(ik4) :: c_npequil(c_cdim)   ! Pointer to linked gas-master species
  integer(ik4) :: c_ncequil(c_cdim,6) ! List of aq species for gas-master
  integer(ik4) :: c_nrequil(c_cdim,6) ! nrh/nrq for linked aqueous spec.
  integer(ik4) :: c_ion(c_cdim)       ! Charge associated with each species
!
! INDICES AND POINTERS FOR AQUEOUS EQUILIBRIA
!
! CHEM MECHANISM NAMES  tchem
! c_tchem(c_cdim)          Chem. species name, a8
! c_treac(6, c_rdim)       Chem reaction species names: A+B->3 products
! c_treach(2, c_rdim)      Henrys law species names: A<=>B
! c_treacq(3, c_rdim)      Equilibrium spec names:  A<=>B+C
! c_treacqq(3, c_rdim)     Special eq spec names:   A<=>B+C

  character(len=8) :: c_tchem(c_cdim)       ! Chem. species name
  character(len=8) :: c_treac(6,c_rdim)     ! Chem reaction species names
  character(len=8) :: c_treach(2,c_rdim)    ! Henrys law species names
  character(len=8) :: c_treacq(3,c_rdim)    ! Equilibrium spec names
  character(len=8) :: c_treaqq(3,c_rdim)    ! Special eq spec names

! CHEM INDICES FOR SPECIAL SPECIES: To identify H2O, H+, OH, etc.

! c_nh2o    Index to identify H2O  in species list
! c_nco2    Index to identify CO2  in species list
! c_nhplus  Index to identify H+   in species list
! c_nohmin  Index to identify OH-  in species list
! c_noh     Index to identify OH   in species list
! c_nho2    Index to identify HO2  in species list
! c_nh2o2   Index to identify H2O2 in species list
! c_no3     Index to identify O3   in species list
! c_nno2    Index to identify NO2  in species list
! c_nno     Index to identify NO   in species list
! c_nno3    Index to identify NO3  in species list
! c_nn2o5   Index to identify N2O5 in species list
! c_nhno3   Index to identify HNO3 in species list
!
  integer(ik4) :: c_nh2o   ! Index to identify H2O  in species list
  integer(ik4) :: c_nco2   ! Index to identify CO2  in species list
  integer(ik4) :: c_nhplus ! Index to identify H+   in species list
  integer(ik4) :: c_nohmin ! Index to identify OH-  in species list
  integer(ik4) :: c_noh    ! Index to identify OH   in species list
  integer(ik4) :: c_nho2   ! Index to identify HO2  in species list
  integer(ik4) :: c_nh2o2  ! Index to identify H2O2 in species list
  integer(ik4) :: c_no3    ! Index to identify O3   in species list
  integer(ik4) :: c_nno2   ! Index to identify NO2  in species list
  integer(ik4) :: c_nno    ! Index to identify NO   in species list
  integer(ik4) :: c_nno3   ! Index to identify NO3  in species list
  integer(ik4) :: c_nn2o5  ! Index to identify N2O5 in species list
  integer(ik4) :: c_nhno3  ! Index to identify HNO3 in species list
  integer(ik4) :: c_nco    ! Index to identify CO   in species list
  integer(ik4) :: c_nch4   ! Index to identify CH4  in species list
  integer(ik4) :: c_nco3   ! Index to identify CO3- in species list

! HVRATES:  Data input for J-value parameterization
!            (jval2.f, developed at Michigan based on Madronich TUV)
!
! Note: jval(56) = Output from J-value parameterization.
!                  Declared internally in hvrates
!  c_nhv(ipar)       Number of parameter intervals for each parameter
!  c_hvmat(ipar,ialt)  Value of parameter interval points
!                   for each of 22 parameters, 40 altitudes
!  c_hvmatb(ipar)    Values for temperature adjustment
!  c_jarray( k,ig,jc)  Output array for 56 jvalues (jc)
!                 for matrix of 510 altitudes and zenith angles (ig)
!                 and 80 conditions (k): 1=base value, >1 adjustments

   integer(ik4) :: c_nhv(22)
   real(rkx) :: c_hvmat(22,40)
   real(rkx) :: c_hvmatb(22)
   real(rkx) :: c_jarray(80,510,56)

!!NUMERICAL SOLUTION  PARAMETERS
!  c_numitr              Maximum number of iterations
!  c_converge            Convergence criteria (1e-02)

  integer(ik4) :: c_numitr               ! Maximum number of iterations
  real(rkx) :: c_converge    ! Convergence criteria (1E-02)

! NOTE  OPTIONAL OUTPUTS RELATING TO NUMERICS:
!    currently in chemvars.EXT
!    xohtest, xnotest, fohtest, final iter, history, geomavg

!     real(rkx) c_ohtest       ! test: dOH/OH or dOH/HO2
!     real(rkx) c_notest       ! test: dNO2/NO2
!     integer c_iter                  ! chem. number of iterations
! NOVEMBER 2007 ADDITION:
!   chemistry mechanism parameters to be used in tracer calculations
!
! STOICHIOMETRY PARAMETERS used in global tracer calculations
!  c_noxchem(13,nr) = stoich parameters for:
!                1-5: net rp Ox; NOx, PAN, HNO3, RNO3;
!                6-13:   NOx->PAN, PAN->NOx;  N<->HNO3;  N<->RNO3, N<=>H+RNO3
!

! REACTION NUMBERS IDENTIFIED FOR SPECIAL REACTIONS.
!   (used in preliminary calculation, oh solution, and global tracers)
!   (2009 addition: include nr for RO2-RO2 reactions for CBMZ special rate calc)
!  c_nro3no   = nr for O3+NO=>NO2
!  c_nrno2x   = nr for NO2+hv=>NO+O3

   real(rkx) :: c_noxchem(13,c_rdim) ! net rp Ox NOx PAN HNO3 RNO3
   integer(ik4) :: c_nro3no              ! nr for O3+NO=>NO2
   integer(ik4) :: c_nrno2x              ! nr for NO2+hv
   integer(ik4) :: c_nro3hv              ! nr for  O3+hv->2OH
   integer(ik4) :: c_nrohno2             ! nr for OH+NO2->HNO3
   integer(ik4) :: c_nrohco              ! nr for OH+CO
   integer(ik4) :: c_nrho2no             ! nr for HO2+NO
   integer(ik4) :: c_nrho22              ! nr for HO2+HO2
   integer(ik4) :: c_nro3no2             ! nr for O3+NO2->NO3
   integer(ik4) :: c_nroho3              ! nr for OH+O3
   integer(ik4) :: c_nrohch4             ! nr for OH+CH4
! 2009 addition
   integer(ik4) :: c_nnrro2              ! Number of special RO2-RO2 reactions
   integer(ik4) :: c_nrro2(c_rdim)       ! nr for special RO2-RO2 reactions

end module mod_cbmz_chemmech
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
