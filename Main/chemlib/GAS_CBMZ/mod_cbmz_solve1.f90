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

module mod_cbmz_solve1
!
  use m_realkinds
  use mod_constants
  use mod_cbmz_chemmech
  use mod_cbmz_chemvars
  use mod_cbmz_chemlocal

  private

  public :: quadchem
!
  contains
!
!  chemsolve.f (quadchem.f)   June, 2009
!   chemsolve_609.f June, 2009
!
!  Programs to solve for species concentrations
!   using the RADICAL BALANCE-BACK EULER solution for photochemistry.
!
!  Full notes on development are in the FILE chemsolnotes
!
!   CURRENT UPDATES:
!
!    June 2009: SETGEOM for NOX CONVERGENCE.  (OPTION).
!    June 2009: Separate Hx sums for Hx with and without RO2 (2009 CHANG
!    changed also in cheminit1.f, chemmech.EXT, chemlocal.EXT )
!    changed LSELF flag - depends on rate but not during 1st iter
!
!    NOTE:  2008 version:  chemmech.EXT before chemvars.EXT.
!
!  CURRENT ISSUES:
!    IFORT compiler.
!    slow and  uncertain aqueous/Br, Cl, NOx solution.  (see PROBLEMS)
!    NOX CONVERGENCE OPTION: geom avg or SETGEOM
!      (currently SETGEOM - experimental.)
!
!  Compile for box model:
!   boxmain.f boxpro.f chemmainbox.f cheminit.f chemrates.f quadchem.f
!    (aquasolve.f) linslv.f jval2.f
!  Include:  chemvars.EXT, chemmech.EXT, chemlocal.EXT, boxvars.EXT
!
!  Input files:
!     REACTION.DAT (REACTION7_GMI_AQHG06)  :  mechanism
!     TUVGRID2                             :  hv data, see jval2.f
!     fort.87  (fort.87hg_gmi_rpa)         :  concentrations+emissions
!     fort.41                              :  Run setup
!       (also REACTION7_GMI_LKS06hv with fort.87gmifull-global - update)
!           (fort.87gmi-global - GMI w/telomers. No REACTION7 yet.)
!
!     Tests:  fort.43test407.
!       Includes 12-hr expo decay (2a vs 2b), + full chem P-L cases.
!     Time test:  See quadchv7.f for comparison between versions.
!
!     Future option: nonsteady  state aqueous
!     Future option: exponential decay  solution (done, currently availa
!     Long-term option: integrate with aerosol solution
!
!  Subroutines:
!     quadchem:  driver program for radical balance solver
!     chemsolve: solution for individual species, pairs or groups.
!     brreac:    calculates rate for specified reaction
!     brpro:  calculates rate, assigns production and loss for reaction
!     excorr:  corrects OH production sums for 'exchange' reactions
!     noxsolve:  solution for NO, NO2 and O3
!     ohsolve:   solution for OH and HO2 (also modifies CO3-)
!     setgeom:   Sets geometric averaging based on past iterations
!     presolve:  Preliminary approximate solution for O3, NOx, OH, NO3.
!     prelump:   Sets initial concentrations for 'lumped' species
!                  and sets initial aqueous=0.
!    midlump:    Re-sets 'prior' partitioning for lumped species
!                  based on interim chemical production rates.
!    postlump:   Sums lumped species and aqueous species into gas master
!                  after calculation is complete.
!
!  Time tests of versions - see quadchv7.f and fort.43test1106,407
!
!  Program  written by Sanford Sillman
!  ---------------------------------------------------------
!  Program history:
!    12/06 Initial program by Sandy Sillman based on boxchemv7.f
!    6/09 Modifications: expo decay, separate Hx w/ and w/o RO2
!  ---------------------------------------------------------
!
! Notes:
! Features of this solution:
!   This part of the program uses sums of gas and aqueous species
!     that are linked through Henry's law and aqueous equilibria.
!   It solves for them as a single "species",
!     based on chemical production/loss,
!     and leaving the gas/aq partitioning intact.
!   A separate program (aquasolve) solves for gas/aqueous partitioning.
!     (the chemistry  and gas-aqueous are solved in iterative sequence)
!
!   This solver uses the 'radical balance' approach.
!    It solves for OH and HO2 based on production/loss of radicals
!    Each radical source/sink is also assigned a sensitivity to OH.
!    OH is solved from a linear interpolation:
!      radical sources-sinks = A + B*OH
!      (This is a fragment of the full back-Euler equation
!         with a sparse matrix)
!
!   Other species are solved individually
!                                 in reactant-to-product sequence.
!
!      (species with mutual interactions are solved individually
!        if they interact slowly - e.g. RO2 and ROOH)
!
!   This version includes solutions for
!    (1) interacting pairs of species (e.g. MCO3-PAN)
!    (2) "chains" of paired species (e.g. Hg0<->HgOH<->HgCl etc)
!         (format:  A<->B, A<->C, C<->C2, etc.)
!    (3) full back-Euler sub-solution for 3 or more interacting species.
!
!   Other features:
!
!    OH-HO2 solution based on Hx sources/sinks WITH or W/O RO2
!      depending on rate of RO2+NO reaction vs timestep (sec-1).
!      If RO2+NO is slow (night), OH-HO2 decoupled from RO2.
!
!    Aqueous CO3- is included as odd-H
!
!    Geometric averaging of present and prior solutions for OH:
!      the geom. averaging parameter is set based on iterative history.
!      (range 0.1-1.0 - it decreases when the iteration oscillates).
!      OPTION:  set minimum value for geometric parameter (0.1 or 0.5)
!              (0.1 for slow, sure convergence.  Max is always 1.)
!
!  PLANNED FUTURE IMPROVEMENTS (June 2009)
!   Vectorized version
!   (CBM-Z version)
!   Test condesned version, add condensed aqueous
!   (Write paper on numerics following EPA 2007 presentation)
!   (write condensed version)
!   Eventually add integrated solution for aerosols
!
!  --------CORRECTED AND ONGOING PROBLEMS-----------------------
!  FUTURE ADD:  Identified below by 'FUTURE ADD'.
!     Exponential decay solution option DONE
!     lself flag for MCO3, AHO2 - not nec. because self reaction is smal
!       but is critical for H2O2, SO5l - add a bypass for better converg
!       MUST BE SKIPPED for NO3, N2O5
!       ADD SKIP   WHEN SMALL.    - rself   - (2009 dont skip for 1st it
!       DONE- NEED TO TEST
!     Non-steady state gas/aqueous partitioning
!     Categories for aqueous/soluble aerosol species
!           (for non-steady-state gas/aq partitioning)
!     NOXSOLVE - why does MULTISOLVE not work;  add O1D;  setgeom for no
!     Vectorize the remaining parts of the program:  hv rates.
!
!     Modify to skip aqueous reactions if LWC=0:
!      assign reactants to direct species, not gas pointer
!      and in chem, do reactions only when species is called
!           - if LWC=0 skip aqueous
!
!  --------CORRECTED AND ONGOING PROBLEMS-----------------------
!
!  ONGOING ISSUES:  (June 2009)
!    IFORT compiler fails unless -C flag is used (xeonsrv)
!    Check convergence criteria - is recent correction OK?
!    NOX CONVERGENCE OPTION: geom avg or SETGEOM?
!      (currently SETGEOM - experimental.)

!    Br-Cl aqueous reactions may be difficult to include, may need multi
!     test cases to see.
!
!   Br nighttime nonconvergence/oscillation in AQHG06,
!    better if HBRG separate in cascade, then BRG-BR2G, rather than trip
!    -> DEBUG TRIPLE IN CASCADE (HBRG-BRG-BR2G)
!    -> NO2-HNO3 OSCILLATION, try setgeom for NOx;  NO.
!     (fort.87rb_hg_rpai w/ aq .3e-6)
!
!   ClOH- has always been difficult due to back-forth reaction with OH.
!     This can cause oscillating solution.
!     Solved using multisolve + geometric averaging for OH.
!
!  STANDARD TESTS:
!  Hour 2 and 12; condensed w/ steady state; AQHG06 w/ acqua=0.3e-6, 298
!      (REACTION.COHC3 w/ fort.87cohc3, fort.87 cohc_u)
!      (REACTION7_GMI_AQHG06 w/ fort.87rb_hg_rpai -> NOTE HBRG CASCADE O
!
!   MOST RECENT CORRECTIONS:
!
!    NOX CONVERGENCE OPTION: SETGEOM instead of geom avg.
!
!      CONVERGENCE OPTION (June 2009): geom avg, lself and 1st few itera
!              (1, 2 or 3 iters, and skip if 1st iter and steady state)
!        NOT USED.
!       (This may help with aqueous test case:
!                               AQHG06 fort.87rb_hg_rpai ,aq=0.3e-6)
!         but hurts slightly with condensed steady state case just below
!
!      BR cascade: HBRG, then BRG-BR2G, rather than triple HBRG-BRG-BR2G
!           (need to investigate further, June 2009)
!
!      Hx WITH OR WITHOUT RO2 (June 2009):
!      Hx sums (oddhsum, etc.) are done for Hx with or without RO2
!      Final solution is weighted based on rate of RO2+NO vs 1/timestep
!      Corrects for night/low NO: OH-HO2 and RO2 become dissassociated.
!      -> Corrects for nighttime  case w/ steady state OH HO2 RO2
!         (REACTION.COHC3 and fort.87cohc3 rural)
!
!     SENHCAT OH ADJUSTMENT (June 2009):
!      senhcat based on rlHO2 vs rpH2O2 - assumes rpH2O2 = HO2+HO2
!      error can occur when there are other sources (TERP+O3)
!      correction: insure rpH2O2<rlHO2
!
!     CONVERGENCE CRITERIA corrected (June 2009)
!
!     PRIOR HX SUM changed to not double-count lumped species (June 2009
!
!     LSELF: flag set for 1st 2 iterations, regardless of reaction rates
!      after 2nd iteration, flag set only if significant self-reaction (
!
!   LSELF  408 version. lself correction  (rself) added.
!     lself flag false unless significant self-reaction.
!     when true, use geometric averaging. (2009: skip for 1st iter)
!
!   NOTE:  408 OPTION, SKIP PAIR ADJUSTMENT if zero production.
!
!      PRODUCT=REACTANT PAIR ADJUSTMENT in brpro/stoicx (2008):
!        subtract from rp, rl if product=reactant,
!        subtract from rppair, rlpair if prod, react in same pair group,
!          but (2008 correction) not if product=reactant.
!
!      PAIR ADJUSTMENT OPTION (minor - in chemsolve):
!      SKIP PAIR ADJUSTMENT in case with zero rpro.  (TEST 2008 OPTION)
!          else small num error DCO3-DPAN-aeDPAN if aero rate=0;
!         since pair structure A-C, C-C2 solves C-C2 first,
!         and solves A-C with fixed C-C2 ratio. trivial error if C2=0.
!         CORRECTED with skip ('TEST 2008').  see dlolcese08/fort.43pair
!
!     * SELF-REACTION CORRECTION: lself-rself in chemsolve (important!)
!       geom avg solution only if significant self-reaction.
!         to prevent error DCO3-DPAN-aeDPAN (lolcese) and NO3-N2O5.
!       Tested for SO5G-HSOG error (aq). (REACTION7_AQHG06, fort.87rb_hg
!
!     NO3-N2O5 no lself: geom avg causes nonconvergence w/NOx at high NO
!             hardwired lself=F, to avoid geom avg.
!             later lself/rself error also fixes this.
!     MULTISOLVE option for aqueous Cl-, Br-
!     Conservation-of-mass adjustment for pair groups.
!     H2O2 still has special geometric solution.
!     CO3 is linked to the radical solution, adjusted based on Hx.
!
!   Corrected problems:
!     OH+Cl=>ClOH caused nonconvergence
!       - fixed with multisolve and geometric averaging.
!     Hg oscillated with aqueous Hg=0:
!       - the problem was a too-high zero protect (if xr<1e-20).
!     EXCORR can cause errors if it links two NON-PAIRED species.
!       - the excorr species must always be solved as pairs.
!
!  Full notes (including TIME TESTS and IMPROVEMENTS)
!     are in quadchv7.f -> chemsolnotes
!
! =============
! --------------------END PROGRAM INFORMATION---------------------
!
! This is the driver program
!   for the radical-balance back-Euler solution for photochemistry.
!
! It calls a preliminary partitioning (prelump),
! Establishes a preliminary solution (presolve)
! And then sets up a standard iteration.
!
! The iteration includes:
!    aquasolve: program to solve for gas-aquoues partitioning
!    chemsolve:solution for individual species or for species pairs
!      called in reactant-to-product order.
!    ohsolve:  radical balance solution for OH and HO2.
!
! Input and output to the program operate through the common blocks
!  (chemvars.EXT and chemmech.EXT)
!
! Input and output are processed by chemmain (file chemmainbox)
!
! Called by:    boxmain.  (This is called for main solution)
! Calls to:
!           chemrates   - set rate constants
!           hvrates     - set hv rate constants
!           prelump     - initial partition of lumped and gs-aqueous sp.
!           presolve    - initial solution for OH
!           midlump     - interim partition of lumped species.
!           aquasolve   - solution for gas-aqueous partitioning
!           chemsolve   - solution for individual species and groups
!           ohsolve     - solution for OH/HO2
!           postlump    - sum lumped species into lump sum,
!                          and aqueous into gas-master.
!
! ---
!
!  Program  written by Sanford Sillman
!  ---------------------------------------------------------
!  Program history:
!    12/06 Initial program by Sandy Sillman based on boxchemv7.f
!  ---------------------------------------------------------
! ---------------------------------------------------------------
!
    subroutine quadchem
!
      implicit none
!
      ! Species list passed to chemsolve
      integer ncsol(c_cdim)
!
      kk=1
!
      ! PRELIMINARY CALLS:  BEFORE SOLUTION
      !   Note, must also call the initial read + process at the start.
      !     chemread, hvread and cheminit
      !
      ! OPTION:  CALL HERE TO SET REACTION RATE AND HV CONSTANTS.
      !   (can be here or in chemmain1.f)
      !
      !  SET REACTION RATE CONSTANTS
      !     call chemrates
      !
      ! SET PHOTOLYSIS RATE CONSTANTS
      !     call hvrates
      !
      ! (end OPTION)

      !
      ! SET H2O CONCENTRATION IN CHEMISTRY ARRAY FROM INPUT VARIABLE
      !
      if (c_nh2o > 0 ) then
       xc(kk,c_nh2o) = c_h2ogas(kk)
      end if
      !
      ! PRELUMP:  SETS 'PRIOR' CONCENTRATIONS (XXO) AND LUMPED SPECIES
      !
      call prelump

      ! PRELIMINARY SOLUTION FOR OH, HO2, O3, NOx.  :
      ! Solution uses the following reactions:
      !       #1:  NO2+hv->NO+O3          #2: NO+O3->NO2
      !       #9:  O3+hv->2OH             #12: NO2+OH->HNO3
      !       #16: OH+CO->HO2             #17: OH+O3->HO2
      !       #18: HO2+NO->OH+NO2         #22: HO2+HO2->H2O2
      !       #46: OH+CH4->HO2
      !
      call presolve

      ! ---------------------------------------------------------------
      !  BEGIN ITERATIVE RUNS
      ! ---------------------------------------------------------------
      runiter: &
      do c_iter = 1 , c_numitr

        ! RE-SET PARTITIONING OF LUMPED SPECIES BASED ON CHEM. PRODUCTION.
        if ( c_iter >= 2 ) call midlump

        ! SET RATES FOR PARAMETERIZED RO2-RO2 REACTIONS, IF ANY
        if ( c_nnrro2 > 0 ) call setro2

        !
        ! CALL AQEOUS CHEMISTRY SOLVER FROM WITHIN ITERATION
        !
        if ( c_h2oliq(1) >= 1.0D-20) call aquasolve

        ! SET SENHCAT FOR FIRST ITERATION. (sensitivity to change in OH).
        ! INITIAL SENHCAT = 1 FOR ODD HYDROGEN; 0 FOR OTHER SPECIES

        if ( c_iter == 1 ) then
          do ii = 1 , 30
            is = 0
            if ( (ii == 8 .or. ii == 9 .or. ii == 10) .or. ii == 3) is = 1
            senhcat(kk,ii) = dble(is)
          end do
        end if

        ! ZERO ACCUMULATED SUMS AT START OF EACH CASCADE LOOP:
        !   oddhsum = summed net source/sink for odd hydrogen radicals
        !               1=w/RO2 2=just OH, HO2, CO3
        !   oddhdel = sum of each net source/sink of oddh
        !                multiplied by sensitivity to OH concentration.
        !   oddhsrc = summed source only
        !   oddhloh = summed loss of OH
        !   oddhlho2 = summed loss of HO2 (removal)
        !   sourcnx = summed NOx source
        !   sinknx = summed NOx sink
        !
        !   rrp, rrl = running sum of production, loss of each species
        !              not set to zero here - zero when spec. is solved.
        !
        !   c_rr  reaction rates.
        !
        ! ALSO,  PRESERVE PRIOR CONCENTRATION HERE
        ! xclastq =PRIOR BUT AFTER AQUASOL)
        do ic = 1 , c_nchem2
          xclastq(kk,ic) = xc(kk,ic)
        end do
        do i = 1 , 2
          oddhsum(kk,i) = d_zero
          oddhsrc(kk,i) = d_zero
          oddhdel(kk,i) = d_zero
          oddhloh(kk,i) = d_zero
          oddhlho2(kk,i) = d_zero
        end do
        sourcnx(kk) = d_zero
        sinknx(kk) = d_zero

        ! -------------------------------------------
        ! CASCADE SOLUTION for individual species or species groups:
        !  SPECIES CALL IS SET BY ORDERED ARRAY c_cascade ->ncsol
        ! -------------------------------------------

        ! LOOP TO CALL SPECIES SOLUTION
        do i = 1 , c_nchem2
          if ( c_cascade(i,1) == 0 ) exit
          do ic = 1 , nsdim
            ncsol(ic) = c_cascade(i,ic)
          end do
          call chemsolve(ncsol)
        end do ! END CASCADE LOOP

        ! -------------------------------------------
        !
        ! -----------------------------
        ! SPECIAL TREATMENT FOR NO3-N2O5 - HNO3
        ! -----------------------------
        !
        ! NO3-N2O5 is solved normally, as paired species
        ! HNO3 is solved first, then NO3-N2O5.
        !  (pair declaration is hard-wired in chemread)
        !
        ! IT IS 'SPECIAL' ONLY BECAUSE IT REQUIRES A SPECIAL CATEGORY
        ! NO3 REACTS WITH MOST SPECIES; SO SPECIES MUST NOT BE ASSIGNED
        ! TO NO3.
        !
        ic1 = c_nno3
        ic2 = c_nn2o5
        ic3 = c_nhno3
        !
        ! ERROR? comment out?
        ! call chemsolve(ncsol)
        do ic = 2 , nsdim
          ncsol(ic) = 0
        end do
        ncsol(1) = ic3
        call chemsolve(ncsol)
        ncsol(1) = ic1
        call chemsolve(ncsol)

        ! -------------------------------------------
        ! SPECIAL TREATMENT FOR NOX-OX AND ODD HYDROGEN FAMILIES
        ! -------------------------------------------
        ! O3-NOx SPECIAL SOLUTION:
        ! SIMULTANEOUS FOR THREE SPECIES, USING OX AND NOX RELATIONSHIPS.
        ! THIS IS INCLUDED IN SUBROUTINE 'NOXSOLVE', INVOKED BY CHEMSOLVE.
        !
        ! TO USE, INVOKE CHEMSOLVE(IC1,IC2,IC3) IN PRECISE ORDER: O3, NO2, NO.
        ! IDENTIFIED BY SPECIAL SPECIES INDEX
        ! THE KEY REACTIONS O3+NO->NO2, NO2+hv->NO+O3 ARE SUMMED IN CPRO.
        !
        do ic = 2 , nsdim
          ncsol(ic) = 0
        end do
        ncsol(1) = c_no3
        ncsol(2) = c_nno2
        ncsol(3) = c_nno
        call chemsolve(ncsol)

        ! ---------------------
        ! ODD HYDROGEN SOLUTION:
        ! ---------------------
        !   CALL SUBROUTINE OHSOLVE(OH,HO2,H2O2) TO SOLVE FOR OH, HO2.
        !   SPECIES IDENTIFIED BY SPECIAL  INDEX

        call ohsolve(c_noh, c_nho2, c_nh2o2)

        !  (computer alternative)
        !     ic9=c_noh
        !     ic10=c_nho2
        !     ic11=c_nh2o2
        !     call ohsolve(ic9,ic10,ic11)

        ! --------------------------------------
        ! END OF ITERATION.  WRITE, TEST AND EXIT.
        ! --------------------------------------

        if ( c_kmax == 1 ) then
          if ( ratek(1,1) >= 1.00D-04 .and. c_ohtest < xfohtest ) then
            c_ohtest = xfohtest
          end if
          if ( c_iter > 3 .and. (c_ohtest < c_converge .and. &
               c_notest > c_converge)) exit runiter
        end if

      end do runiter
      !
      ! -----------------------
      ! END ITERATIVE LOOP
      ! -----------------------

      ! POSTLUMP:  SUM LUMPED SPECIES INTO LUMP-SUM AND CONVERT CONCENTRATIONS
      ! TO LUMP-PARTITION FRACTIONS.  (Also sum Aq-Equil species?)

      call postlump

 1611 format(/,'TEST SENHCAT 1-20,31-40:',/,(10f6.3))
 1701 format(/' TEST FOR OH,NO CONVERGENCE =',i3,2(1pe10.3))
 1801 format(//,' SUMMARY WRITE:  ITER =',i3, '    VECTOR KKW=',i3)
 1802 format(/,' #  IC     XCOUT     XCIN  XCF/AV       RL',       &
                    '        RP        dXC      dR')
 1803 format(i4,a8,2(1pe10.3),0pf7.3,1x,(4(1pe10.3)))
 1804 format('IDATE xhour lat lon temp=',i8,4f8.2)
 1805 format('JPARAMS: zenith altmid dobson so2 no2 aaerx aerssa', &
        ' albedo',/,'cld-below cld-above claltB claltA temp date', &
         /,(8f10.3))
 1806 format(//,' REACTION RATES:  ITER =',i3, '    VECTOR KKW=',i3)
 1807 format(i4,2x,a8,'+',a8,'=>',a8,'+',a8,'+',a8,2(1pe10.3))
 1901 format(//,70('-'),/,'BEGIN ITERATION =',i3)
 1902 format(70('-'),/)
!
    end subroutine quadchem
!
! ----------------------------------------------------
! SOLVER SUBROUTINES - SOLVE FOR INDIVIDUAL SPECIES CONCENTRATIONS
! ----------------------------------------------------
!
!     FUTURE DO:  back-Euler solution for multi-species
!
! This is the main solver for individual species concentrations.
! It includes:
! (1) A solution for an individual species -
!      (regarded as a sum of gas + aqueous linked through equilibria)
!
! (2a) A solution for a pair of rapidly interacting species
!
! (2b) A solution for a 'chain' of linked pairs
!                                       (A<->B, A<->C1, C1<->C2, etc.)
!
! (3) A reverse-Euler solution for a group of 3+  interacting species.
!      ('MULTISOLVE')
!
! (4) Call to special solution for O3-NO2-NO ('NOXSOLVE')
!
! It also establishes rates for reactions associated with the species;
! And calculates summed production, loss of reaction products,
!  and odd hydrogen radical sums.
!
!  FOR A SINGLE SPECIES OR PAIR, MAKE IC2<=0.
!
! THE SOLUTION SEQUENCE:
!  (a) Calculate losses/rates of reactions tht affect the species
!      (including linked gas+aqueous equilibrium species)
!  (b) Solve for species concentrations
!        with separate solutions for:
!         SINGLE SPECIES
!         PAIR AND PAIR CHAIN
!         NOXSOLVE
!         MULTISOLVE
!  (c) Add production of product species and odd-h radicals
!       to the production sums.
!
!  Input: Prior species concentrations, reactions and rate constants.
!  Output:  Updated species concentrations (xc).
!
! Called by:  quadchem.
! Calls to:
!     brpro - sums production and loss from reaction
!     brreac - calculates reaction rates.
!     setgeom - establish geometric mean parameter for iteration
!     noxsolve - solution for O3-NO-NO2
!     LINSLV - matrix iteration for MULTISOLVE solution.
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
! ----------------------------------------------------------------
!
!  PREVIOUS:  FOR TWO SPECIES THE 'EXCHANGE' SPECIES MUST COME FIRST.
!  NO LONGER.  EXCHANGE IS NOT DEPENDENT ON ORDER-OF-SPECIES.
!  (ONLY NOXSOLVE USES RPRO, RLOSS, CPRO FROM CHEMSOLVE
!    IN A SPECIFIC ORDER.  EXCORR IS NOT CALLED FOR NOX.)
! ----------------------------------------------------------------
!
!
    subroutine chemsolve(ncsol)
!
      implicit none
      ! Species list passed to chemsolve
      integer :: ncsol(c_cdim)
      ! Flag for exponential decay solution
      logical :: llexpo(c_kvec)
      ! Flag for completion of expo deca
      logical :: doneexpo(c_kvec)
      ! For back-Euler solution
      real(dp) :: ax(100,100), bx(100), xx(100)
      real(dp) :: xpair
!     net pair group  production w/ cons. mass.
      ! Summed factor for multi group
      real(dp) :: xmulti
!     net multi group  production w/ cons. mass.
      ! Rate of self-reaction
      real(dp) :: rself(c_kvec, 20)
      ! flag for self-reaction
      logical :: lself
!
      kk = 1
      !
      ! PRELIMINARY:  SELF-REACTION FLAG
      !
      lself = .false.

      !
      ! PRELIMINARY:  make ic point to GAS-MASTER and PRIMARY OF PAIRED SPECIE
      ! and count nsol = number of active species groups solved for.
      !
      nsol = 0
      loopdim: &
      do nc = 1 , nsdim
        ic = ncsol(nc)
        if ( ic > 0 ) then
          ncsol(nc) = c_nppair(c_npequil(ic),2)
          nsol = nc
        else
          exit loopdim
        end if
      end do loopdim
      !
      ! ---------------------------------------------
      ! SUMMATION OF RLOSS - CHEMICAL LOSSES AND CROSS-PRODUCTS
      ! ---------------------------------------------
      !
      !  -----------------------
      !  PRELIMINARY:  ESTABLISH LIST OF LINKED SPECIES
      !  -----------------------
      !  nsol = number of pair groups to be solved for (set above)
      !  ncsol(n) = ic of head of each pair group to be solved for (input abov
      !  nsolv = number of species to be solved for
      !          (includes pairs and multi-solve)
      !
      ! ncsolv(i) = ic of each species in order (ncsol(1) and paired, ncsol(2)
      !
      ! nssolv(i) = is of each species = multi-species group for ncsolv
      !
      ! Note:  loop over nsolv automatically solves for multiple-species group
      !         as well as single pair grouping
      nsolv = 0
      do is = 1 , nsol
        nsolv = nsolv + 1 + c_nppair(ncsol(is),3)
      end do
      nc = 0
      do is = 1 , nsol
        ic = ncsol(is)
        if ( ic > 0 ) then
          do i = 1 , (c_nppair(ic,3)+1)
            icc = ic
            if ( i > 1 ) then
              icc = c_nppair(ic,i+2)
            end if
            nc = nc+1
           ncsolv(nc) = icc
            nssolv(nc) = is
          end do
        end if
      end do

      ! PRELIMINARY:  FLAG FOR EXPONENTIAL DECAY SOLUTION
      !   Use expo decay only for single species solution (unpaired, no aq.)
      !   (FUTURE DO.  Not included yet.)

      llexpo(kk) = .false.
      if ( c_lexpo(kk) ) then
        if ( (nsolv == 1) .and. (c_nequil(ncsol(1))== 0) .and. &
             (c_tchem(ncsol(1) ) /= '    H2O2') .and. &
             (.not. c_lsts(ncsol(1))) ) llexpo(kk) = .true.
      end if

      ! PRELIMINARY: CALCULATE REACTION RATES,RP and RL  FOR PRODUCT REACTIONS
      ! ASSOCIATED WITH SPECIES (i.e. that are sources, not sinks)
      !
      do nc = 1 , nsolv
        ics = ncsolv(nc)
        icpair = c_nppair(ics,2)
        if ( c_nnrchp(ics) > 0 ) then
          do i = 1 , c_nnrchp(ics)
            nr = c_nrchmp(ics,i)
            call brpro(nr)
          end do
        end if
      end do
      !
      ! ZERO RPRO, RLOSS, XRP, CPRO.  (XRM, XRRM for multisolve) (
      !      (RPPAIR, RLPAIR, XRPPAIR, XRRPAIR for pair group sums)
      !
      do i = 1 , nsolv
        rpro(kk,i) = d_zero
        rloss(kk,i) = d_zero
        xrp(kk,i) = d_zero
        xrm(kk,i) = d_zero
        xrrm(kk,i) = d_zero
        rpm(kk,i) = d_zero
        rlm(kk,i) = d_zero
        rself(kk,i) = d_zero
        do ii = 1 , nsolv
          cpro(kk,i,ii) = d_zero
          cpm(kk,i,ii) = d_zero
        end do
      end do
      xrppair(kk) = d_zero
      xrrpair(kk) = d_zero
      rpmulti(kk) = d_zero
      rlmulti(kk) = d_zero

      !  -----------------
      !  LOOP TO CALCULATE CHEM. LOSS RATES AND CROSS-PRODUCTS
      !  -----------------
      loopnsolv: &
      do nc = 1, nsolv
        ics = ncsolv(nc)
        iscs = nssolv(nc)
        icpair = c_nppair(ics,2)
        !
        ! ESTABLISH RPRO, RLOSS AND XRP FOR THE 'BASE' SPECIES
        !
        rpro(kk,nc) = rpro(kk,nc) + rrp(kk, ics) + c_xcin(kk,ics)
        rloss(kk,nc) = rloss(kk,nc) + rrl(kk, ics)
        xrp(kk,nc) = xrp(kk,nc) + xc(kk,ics)
        !
        ! ADD ALL AQUEOUS-EQUILIBRIUM SPECIES INTO SUMMED RPRO, RLOSS AND XRP.
        ! CONVERTING AQUEOUS INTO GAS UNITS (AVOGADRL)
        ! Also include RAINOUT as loss.
        !
        if ( c_nequil(ics) > 0 ) then
          do neq = 1 , c_nequil(ics)
            ic = c_ncequil(ics,neq)
            xrp(kk,nc) = xrp(kk,nc) + xc(kk,ic)*c_h2oliq(kk)*avogadrl
            rpro(kk,nc) = rpro(kk,nc) + rrp(kk,ic)
            rloss(kk,nc) = rloss(kk,nc)+rrl(kk,ic) + &
                           c_rainfr(kk) * xc(kk,ic)*c_h2oliq(kk)*avogadrl
          end do
        end if
        !
        ! END EQUILIBRIUM SUMMATION.
        !
        ! ADD INITIAL SPECIES XRP, RPRO, RLOSS TO XRPPAIR, RPPAIR, RLPAIR
        !  IF PAIR GROUP HAS >3 MEMBERS AND NO MULTISOLVE (FOR NORMALIZATION)
        !
        if ( nsolv >= 3 .and. nsol == 1 ) then
          xrppair(kk) = xrppair(kk) + xrp(kk,nc) * c_pairfac(ics)
          rppair(kk,icpair) = rppair(kk,icpair) + rpro(kk,nc) * c_pairfac(ics)
          rlpair(kk,icpair) = rlpair(kk,icpair) + rloss(kk,nc) * c_pairfac(ics)
        end if
        !
        ! ADD INITIAL SPECIES XRP, RPRO, RLOSS TO XRM, RPM, RLM: for MULTISOLVE
        !  (note:  no PAIRFAC, different from PAIR NORMALIZATION)
        !
        if ( nsol > 1 ) then
          xrm(kk,iscs) = xrm(kk,iscs) + xrp(kk,nc)
          rpm(kk,iscs) = rpm(kk,iscs) + rpro(kk,nc)
          rlm(kk,iscs) = rlm(kk,iscs) + rloss(kk,nc)
          rpmulti(kk) = rpmulti(kk) + rpro(kk,nc) * c_multfac(ics)
          rlmulti(kk) = rlmulti(kk) + rloss(kk,nc) * c_multfac(ics)
        end if
        !
        ! ADD LOSS REACTIONS AND CROSS-PRODUCT REACTIONS TO RLOSS,RPRO, CPRO.
        ! ALSO SET SELF-REACTION FLAG (lself)
        ! ALSO SET FLAG FOR SECOND REACTANT ALSO INCLUDED IN SPECIES LIST (nc2)
        !     (pair or multisolve.  Important for NO+O3->NO2.)
        !
        if ( c_nnrchem(ics) > 0 ) then
          !
          ! SET INDICES
          !
          do i = 1 , c_nnrchem(ics)
            nr = c_nrchem(ics,i)
            icr1 = c_reactant(nr,1)
            icr2 = c_reactant(nr,2)
            icr1 = c_npequil(icr1)
            ! index to count loss from 2nd reactant
            nc2 = 0
            if ( icr2 > 0 ) then
              icr2 = c_npequil(icr2)
              if ( icr1 == ics .and. icr2 == ics) lself = .true.
              !  FEB 2005OPTION:
              !     IF REACTANTS ARE IN SAME PAIR GROUP, SET SELF-REACTION FLAG
              !
              !  FUTURE DO, FIX SELF-REACTION.
              !  (In SELF-REACTION, do not set 2nd reactant - it messes up CPRO
              !   unless CPRO algorithm changed. )
              !
              if ( c_nppair(icr1,2) == c_nppair(icr2,2)) lself = .true.
              do ncc = 1 , nsolv
                if ( ncc /= nc ) then
                  if ( ncsolv(ncc) == icr1 .or. ncsolv(ncc) == icr2) then
                    nc2=ncc
                    ics2 = ncsolv(ncc)
                  end if
                end if
              end do
            end if
            !
            ! IF REACTION IS PARAMETERIZED RO2-RO2 SET SELF-REACTION FLAG
            !
            if ( c_nrk(nr) == -13 ) lself = .true.
            !
            ! CALL BRREAC - to establish reaction rate
            !
            call brreac(nr)
            !
            ! ADD LOSS FOR FIRST REACTANT
            !
            rloss(kk,nc) = rloss(kk,nc) + c_rr(kk,nr)
            !  stoiloss OPTION:
            ! IF SELF-REACTION: ADD AGAIN.  (2nd reactant index is zero)
            if ( icr1 == icr2 ) then
              rloss(kk,nc) = rloss(kk,nc) + c_rr(kk,nr)
              rself(kk,nc) = rself(kk,nc) + d_two*c_rr(kk,nr)
            end if
            !
            ! IF PARAMETERIZED RO2: add to RSELF.
            ! (note: rate const includes 2x for
            !
            if ( c_nrk(nr) == -13 ) then
              rself(kk,nc) = rself(kk,nc) + c_rr(kk,nr)
            end if
            !
            ! SECOND SPECIES OPTION:
            ! ADD LOSS OF 2ND REACTANT FOR DIFFERENT SPECIES. (for NO+O3->NO2).
            ! ADD TO RLOSS1 ALSO (which may  have been written before
            !    - RLOSS1 only excludes CPRO and PSEUDO-SOLUTION.
            !
            ! NOTE: omit stoiloss, which adjusts only for main species
            !       appearing in reactant and product.
            !
            if ( nc2 > 0 ) then
              rloss(kk,nc2) = rloss(kk,nc2) + c_rr(kk,nr)
              rloss1(kk,nc2)=rloss1(kk,nc2) + c_rr(kk,nr)
            end if
            !
            ! ADD PRODUCTS TO CPRO or RPRO.
            ! CPRO MATRIX  (kk,reactant,product).
            ! (note prodarr=gas-master product)
            !
            do ncc = 1 , nsolv
              if ( c_prodarr(nr,ncsolv(ncc)) /= 0 ) then
                icc = ncsolv(ncc)
                !
                ! Adjust LOSS if PRODUCT EQUALS REACTANT, skip PRO.
                !
                if ( ncc == nc .or. ncc == nc2 ) then
                  stoicx = c_prodarr(nr,ncsolv(ncc))
                  if ( stoicx > d_one) stoicx = d_one
                  rloss(kk,ncc)=rloss(kk,ncc) - c_rr(kk,nr)*stoicx
                  if ( ncc == nc2 ) then
                    rloss1(kk,ncc) = rloss1(kk,ncc) - c_rr(kk,nr)*stoicx
                  end if
                else
                  !
                  ! ADD TO CPRO - ONLY IF SPECIES ARE DIRECT PAIRS
                  ! ALSO IF MULTI.
                  ! (MULTI CPRO is added to RPRO for pair solution below,
                  !  not here - to avoid doublecounting in RPM.)
                  !
                  if ( (c_nppair(ics,1) == icc .or. &
                        c_nppair(icc,1) == ics) .or. &
                        (nssolv(nc) /= nssolv(ncc)) ) then
                    cpro(kk,nc,ncc) = cpro(kk,nc,ncc) + &
                                      c_rr(kk,nr)*c_prodarr(nr,ncsolv(ncc))
                  else
                    !
                    ! ADD TO RPRO IF NOT DIRECT PAIR
                    !
                    rpro(kk,ncc) = rpro(kk,ncc) + &
                                   c_rr(kk,nr)*c_prodarr(nr,ncsolv(ncc))
                  end if
                end if
                !
                ! CPRO FROM SECOND REACTANT
                ! FEBRUARY 2005 CHANGE:
                ! IF DIRECT PAIR, ADD TO CPRO AND SUBTRACT FROM RPRO
                ! (it should have just been added to RPRO for 1st reactant)
                ! OTHERWISE, SKIP.
                ! (Note:  for MULTI, CPRO is recorded and added to RPRO for
                ! PAIRSOL. This must be done just once, for 1st reactant,
                ! to avoid double change
                !
                if ( nc2 > 0 ) then
                  if ( (c_nppair(ics2,1) == icc .or. &
                        c_nppair(icc,1) == ics2) ) then
                    !
                    ! DIRECT PAIR: ADD 2ND REACTANT CPRO, ADJUST RPRO.
                    ! (NCPRO CUT)
                    !
                    cpro(kk,nc2,ncc) = cpro(kk,nc2,ncc) + &
                                       c_rr(kk,nr)*c_prodarr(nr,ncsolv(ncc))
                    rpro(kk,ncc) = rpro(kk,ncc) - &
                                   c_rr(kk,nr)*c_prodarr(nr,ncsolv(ncc))
                  end if
                end if
              end if
            end do
            !
            ! END LOOP - CPRO MATRIX
            !
            ! CALCULATE NET PRODUCTION AND LOSS FOR PAIR GROUP
            !  (FOR NORMALIZATION OF GROUP w>3 MEMBERS)
            !
            if (nsolv >= 3 .and. nsol == 1) then
              xpair = d_zero
              if ( c_nppair(icr1,2) == icpair ) then
                xpair = xpair - c_pairfac(icr1)
              end if
              if ( icr2 > 0 ) then
                if ( c_nppair(icr2,2) == icpair ) then
                  xpair =  xpair - c_pairfac(icr2)
                end if
              end if
              do ncc = 1 , nsolv
                xpair = xpair + c_prodarr(nr,ncsolv(ncc)) * &
                        c_pairfac(ncsolv(ncc))
              end do
              if ( xpair < d_zero ) then
                rlpair(kk,icpair)  = rlpair(kk,icpair) - xpair*c_rr(kk,nr)
              end if
              if ( xpair > d_zero ) then
                rppair(kk,icpair) =rppair(kk,icpair) + xpair*c_rr(kk,nr)
              end if
            end if
            !
            ! CALCULATE NET PRODUCTION AND LOSS FOR PAIR GROUP
            !  FOR MULTISOLVE (AND NORMALIZATION OF GROUP w>3 MEMBERS)
            !
            if ( nsol > 1 ) then
              isr1 = nssolv(nc)
              isr2 = 0
              if ( nc2 > 0 ) then
                isr2 = nssolv(nc2)
              end if
              !
              ! LOOP THROUGH PAIR GROUPS
              !
              xmulti = d_zero
              do is = 1 , nsol
                ic = ncsol(is)
                !
                ! ESTABLISH NET PRO/LOSS FOR PAIR GROUP
                !
                xpair = d_zero
                if ( is == isr1 ) xpair = xpair - d_one
                if ( is == isr2 ) xpair = xpair - d_one
                if ( is == isr1 ) xmulti = xmulti - c_multfac(ic)
                if ( is == isr2 ) xmulti = xmulti - c_multfac(ic)
                !
                ! ERROR HERE  - THIS WAS c_nppair(icpair,.) THROUGHOUT LOOP
                !             -  SHOULD BE c_nppair(ic,..)
                !
                do iss = 1 , (c_nppair(ic,3)+1)
                 icc = ic
                 if ( iss > 1 ) icc = c_nppair(ic,iss+2)
                 xpair = xpair + c_prodarr(nr,(icc))
                 xmulti = xmulti + c_prodarr(nr,(icc))*c_multfac(icc)
                end do
                !
                ! ADD TO RLM LOSS
                !
                if ( xpair < d_zero ) then
                  rlm(kk,is) = rlm(kk,is) - xpair*c_rr(kk,nr)
                  !
                  ! SUBTRACT CROSS-LOSS FROM CPM.  (add negative XPAIR)
                  ! (NOTE: Just one cross-loss; 2nd cross-loss will be
                  ! counted when loop reaches the other reactan
                  !
                  if ( isr1 > 0 .and. isr1 /= is ) then
                    cpm(kk,isr1,is) = cpm(kk,isr1,is) + xpair*c_rr(kk,nr)
                  end if
                  if ( isr2 > 0 .and. isr2 /= is ) then
                    cpm(kk,isr2,is) = cpm(kk,isr2,is) + xpair*c_rr(kk,nr)
                  end if
                  !
                  ! ADJUSTMENT FOR SELF-REACTION
                  ! (includes 2 reactants in one pair group
                  ! Back-Euler: dR/dx for rate kx^2 is 2kx;
                  !             loss per reaction=2.
                  ! So total loss is 2R;  sensitivity is 4 R.
                  ! Add 2R to RPM, 4R to RLM.)
                  ! HERE:  XPAIR*R already added to RLM.
                  ! Add again to RLM, RPM.
                  !
                  if ( isr1 == is .and. isr2 == is ) then
                    rlm(kk,is) = rlm(kk,is) - xpair*c_rr(kk,nr)
                    rpm(kk,is) = rpm(kk,is) - xpair*c_rr(kk,nr)
                  end if
                end if
                !
                ! ADD TO RPM
                !
                if ( xpair > d_zero ) then
                  rpm(kk,is) = rpm(kk,is) + xpair*c_rr(kk,nr)
                  !
                  ! ADD TO CPM  CROSS-PRO.
                  !
                  if ( isr1 > 0 .and. isr1 /= is ) then
                    cpm(kk,isr1,is) = cpm(kk,isr1,is) + xpair*c_rr(kk,nr)
                  end if
                  if ( isr2 > 0 .and. isr2 /= is ) then
                    cpm(kk,isr2,is) = cpm(kk,isr2,is) + xpair*c_rr(kk,nr)
                  end if
                end if
              end do
              !
              ! END LOOP THROUGH PAIR GROUPS
              !
              ! MULTI GROUP SUM
              !
              if ( xmulti < d_zero ) then
                rlmulti(kk) = rlmulti(kk) - xmulti*c_rr(kk,nr)
              end if
              if ( xmulti > d_zero ) then
                rpmulti(kk) =rpmulti(kk) + xmulti*c_rr(kk,nr)
              end if
            end if
            !
            ! END CALCULATE NET PRODUCTION AND LOSS FOR PAIR GROUP
            !
          end do
        end if
        !
        ! END - ADD LOSS REACTIONS AND CROSS-PRODUCTS
        !
        !
        ! NONSTEADY STATE ADJUSTMENT FOR PAIR GROUP SUM
        !   (ahead of STEADY STATE ADJUSTMENT to use RLOSS=0)
        !   (note RPPAIR already includes XXO)
        !
        if ( nsolv >= 3 .and. nsol == 1 ) then
          if ( .not. c_lsts(ics) .or. rloss(kk,nc) == 0 ) then
            rlpair(kk,icpair) = rlpair(kk,icpair) + xrp(kk,nc)*c_pairfac(ics)
          end if
          if ( rloss(kk,nc) <=  d_zero ) then
            rlpair(kk,icpair) = rlpair(kk,icpair) + 1.0D-08*c_pairfac(ics)
          end if
        end if
        !
        ! NONSTEADY STATE ADJUSTMENT FOR MULTISOLVE SUM
        !   (ahead of STEADY STATE ADJUSTMENT to use RLOSS=0)
        !
        if ( nsol > 1 ) then
          if ( .not. c_lsts(ics) .or. rloss(kk,nc) == 0 ) then
            rlm(kk,iscs) = rlm(kk,iscs) + xrp(kk,nc)
            rlmulti(kk) = rlmulti(kk) + xrp(kk,nc)*c_multfac(ics)
          end if
          if ( rloss(kk,nc) <= d_zero )  then
            rlm(kk,iscs) = rlm(kk,iscs) + 1.0D-08
            rlmulti(kk) = rlmulti(kk) + 1.0D-08*c_multfac(ics)
          end if
        end if
        !
        ! NONSTEADY STATE ADJUSTMENT: IF NONSTEADY STATE OR ZERO RLOSS, ADD XRP
        !  (EQUIVALENT TO RL/XR + 1).
        !  (NOTE:  STEADY-STATE OPTION STILL INCLUDES XXO IN RPRO.
        !   FOR STEADY STATE XXO SHOULD EQUAL EITHER EMISSIONS OR ZERO.)
        !    (XXO SET IN PRELUMP)
        ! WITH ZERO PROTECT FOR RLOSS
        if ( .not. c_lsts(ics) .or. rloss(kk,nc) == 0 ) then
          rloss(kk,nc) = rloss(kk,nc) + xrp(kk,nc)
        end if
        if ( rloss(kk,nc) <= 0. ) rloss(kk,nc) = 1.0D-08
        !
        !  (OPTION: FOR A SELF-REACTION (stoiloss>1, e.g. HO2+HO2, stoiloss=2.)
        !  AN EXACT LINEARIZED NR SOLUTION WOULD REQUIRE
        !  RLOSS=+2*RR*STOILOSS, RPRO=+RR*STOILOSS. (d/dR = twice loss rate.)
        !  TO IMPLEMENT THIS, IF STOILOSS>1, ADD ADDITIONAL RR*STOILOSS
        !  TO BOTH RLOSS AND RPRO, WITHIN LOOP 125 (and 225 for twosolve).
        !  This is cut because IT MAY SCREW HO2+HO2 REACTIONS.
        !  SAVE RLOSS1, RPRO1
        !  =RPRO without CPRO; RLOSS w/o internal solution adjustments.
        rloss1(kk,nc) = rloss(kk,nc)
        rpro1(kk,nc) = rpro(kk,nc)
      end do loopnsolv
      !  ---------------------
      !  END LOOP TO CALCULATE CHEM. LOSS RATES AND CROSS-PRODUCTS
      !  ---------------------
      !
      ! SET SELF-REACTION FLAG:  false if self-reaction is insignificant
      !  (ADDED 408 APRIL 2008)
      !  (note future option:  save history and use to prevent oscill.)
      !  (note - this prevents NO3-N2O5 error, but hard-wired NO3 retaine
      !  LSELF OPTION: 2009 CORRECTION - skip for iter=1)
      !
      if ( lself .and. c_iter > 1 ) then
        lself = .false.
        do nc = 1 , nsolv
          if ( rself(kk,nc) > 0.33D0*rloss1(kk,nc) ) lself = .true.
        end do
      end if
      !
      ! END-  SET SELF-REACTION FLAG
      !
      ! OPTION:  Add external CPRO to RPRO for multi-species.
      ! CURRENT OPTION:
      ! Internal XRR will be solved with RPRO from prior values of multi-spe
      ! RLOSS already includes loss to other multi-species groups.
      ! ALTERNATIVE:
      ! subtract CP to multi-species from RLOSS
      !
      if ( nsol > 1 ) then
        do nc = 1 , nsolv
          ics = ncsolv(nc)
          do ncc = 1 , nsolv
            if ( nssolv(ncc) /= nssolv(nc) ) then
              rpro(kk,nc) = rpro(kk,nc) + cpro(kk,ncc,nc)
            end if
          end do
        end do
      end if
      ! ---------------------------------------------
      ! END - SUMMATION OF RLOSS - CHEMICAL LOSSES AND CROSS-PRODUCTS
      ! ---------------------------------------------
      !
      ! ---------------------------------------------
      ! CALCULATE SPECIES CONCENTRATIONS:  INTERNAL PAIRS
      ! ---------------------------------------------
      !
      !  THE ALGORITHM:
      !   For two interacting species:  A<->B
      !   Back-Euler solution:  xr = xrprior * (Pi + CjiPj/Lj)/(Li-CijCji/Lj)
      !                            = xrprior * (   Pi'       )/(  Li'      )
      !     where Li = prior absol. loss rate + xrprior (non-steady state)
      !           Pi = prior production + c_xcin
      !           Cji = cross production of i from j
      !           Pi', Li' = pseudo-P, L with cross prod'n from paired subspec
      !
      !     (neg. not possible since Li>Cij, Lj>Cji)
      !
      !   For a chain of interacting species:  A<->B, A<->C, C<->C2, etc.
      !     in order from primary to subspecies:  Cij=0 for j<i
      !     CjiPj, etc. is implicit sum over all paired subspecies j
      !     Solve from sub-species to primary species Pi'/Li'
      !        where Pi', Li' includes Cji, j>i
      !     Then solve from primary Xi to subspecies Xj
      !      with  Cij from primary Xi updated and added to Pj
      !
      !   This effectively uses a prior secondary partitioning (C<->C2)
      !    to generate primary partitioning (A<->B, A<->C) and primary XR
      !    It uses primary (A) to solve for secondary (B, C, etc.)
      !      while preserving C<->C2 partitioning (implicit in Pi'/Li')
      !
      !   FOR MULTISOLVE/NOXSOLVE (STANDARD OPTION):
      !    PAIR SOLUTION is normalized to keep original within-pair sum unchan
      !    NOXSOLVE or MULTISOLVE changes pair sums (preserving internal parti
      !     based on interaction between pair groups.
      !   OPTION for NOXSOLVE/MULTSOLVE:
      !    should the PAIR SOLUTION include cross-production from other groups
      !    LOSS includes CP to other group.
      !    Either subtract prior Cij from Li, or add prior Cij to Lj
      !     (ABOVE)
      ! ------------------------------------
      !
      ! ---------
      ! LOOP TO  ADJUST RP, RL FOR SUB-SPECIES CROSS PRODUCTION (PSEUDO-RP, RL
      ! ---------
      !  RPRO, RLOSS ADJUSTED FOR INTERNAL CP PRODUCTION FROM SUB-SPECIES ONLY
      !  (This loop automatically includes each multi-species grouping.)
      !  (Loop order:  from sub-species to primary species;
      !     so that primary species includes subspecies RP with sub-CP adjustm
      ! (RPRO1, RLOSS1 preserves original RPRO, RLOSS w/o adjustments)
      !
      if ( nsolv > 1 ) then
        do nc = (nsolv-1) , 1 , -1
          ics = ncsolv(nc)
          do ncc = (nc+1) , nsolv
            if ( nssolv(nc) == nssolv(ncc) ) then
              icc = ncsolv(ncc)
              rpro(kk,nc) = rpro(kk,nc) + &
                  cpro(kk,ncc,nc) * (rpro(kk,ncc)/rloss(kk,ncc))
              calpha(kk) = cpro(kk,ncc,nc)*(cpro(kk,nc,ncc)/rloss(kk,ncc))
              rloss(kk,nc) = rloss(kk,nc) - calpha(kk)
              if ( rloss(kk,nc) <= d_zero ) then
                write(c_out,121 ) c_tchem(ncsol(1)), ncsol(1),        &
                  nc, ncc, rloss(kk,nc), calpha(kk), cpro(kk,ncc,nc),  &
                  cpro(kk,nc,ncc), rloss(kk,ncc)
              end if
            end if
          end do
        end do
      end if
      ! ---------
      !  END LOOP TO ADJUST RP, RL FOR CP
      ! ---------
      !
      !
      ! ---------
      !  LOOP TO SOLVE SPECIES CONCENTRATIONS- INTERNAL PAIRS
      ! ---------
      !   Order:  primary-to-subspecies chain among paired species
      !   After calculation, add modified CPRO to CP-down-species.
      !
      do nc = 1 , nsolv
        ics = ncsolv(nc)
        is = nssolv(nc)
        icpair = c_nppair(ics,2)
        !
        ! BACKWARD EULER SOLUTION FOR SINGLE SPECIES
        ! The true equation is =RP/(RL/XR); modified to avoid second divide.
        ! (skip from here - move lower)
        !
        !  EXPONENTIAL DECAY OPTION:  SINGLE SPECIES ONLY.
        !
        !  EQUATION:  dc/dt= r-kc, c=Co at t=0. (r=rp/time +emission rate)
        !  Ct = r/k + (Co - r/k) exp(-kt)
        !  Cavg = r/k + ((Co-r/k)/kt)*(1- exp(-kt))    =  (Co-Cf)/kt + r/k
        !       Cf= Co + rt - kt*Cav    --checks
        !
        ! Here, Co = xcin - xcemit
        !       rt = rpro-xcin+xcemit
        !       kt = (rloss-xrp)/xrp
        !
        !  Result saved as:
        !   xc = Cavg (used for reaction rate calculations)
        !   (xrr = Cavg, then turned into Cavg/Cprior 'NEW/OLD' below,,
        !   and used to get xc gas,aq -as in standard back-Euler.)
        !   xcfinr = Cf/Cavg  (used in postlump to get xcout)
        !   NOTE: xcfinr=1 in BACK-EULER solution.
        ! NOTE ERROR CORRECTIONS:
        !     cbeta: form A(1-exp) + B(exp) , not A + (B-A)*exp
        !        (else large-number error and negative value)
        !     ZERO PROTECT for xrr
        !     Use back-Euler for very small values.
        doneexpo(kk) = .false.
        if ( llexpo(kk) ) then
          if ( xrp(kk,nc) > 0  .and. &
              ( rloss(kk,nc)-xrp(kk,nc) > d_zero ) .and. &
              ( rloss(kk,nc)+rpro(kk,nc) > d_one ) ) then
            calpha(kk) = rloss(kk,nc)/xrp(kk,nc) - d_one
            ! cbeta = rt/kt = (rpro-xcin+xcemit)/calpha
            cbeta(kk) = (rpro(kk,nc)-c_xcin(kk,ics) + &
                        c_xcemit(kk,ics))/calpha(kk)
            !    Cf =   r/k + (    Co     - r/k) exp(-kt)
            !    Cf =  cbeta + (xcin-xcemit-cbeta)*exp(-calpha)
            !    NOTE: This must be greater than (xcin-xcemit)
            !
            ! The following line caused LARGE-NUMBER NUMERICAL ERRORS
            !    xcfinr(kk,ic) = cbeta(kk)
            !              +(c_xcin(kk,ic) - c_xcemit(kk,ic) - cbeta(kk))
            !              *exp(0.-calpha(kk))
            !
            ! Corrected:
            !
            xcfinr(kk,ics) = cbeta(kk) * (d_one - exp(d_zero-calpha(kk)) ) + &
                   (c_xcin(kk,ics) - c_xcemit(kk,ics))*exp(d_zero-calpha(kk))
            !    Cav = (    Co     -Cf)/kt    + r/k
            !    Cav = (xcin-xcemit-Cf)/calpha + cbeta
            !
            xrr(kk,nc) = cbeta(kk) + (c_xcin(kk,ics) - c_xcemit(kk,ics) - &
                         xcfinr(kk,ics))/calpha(kk)
            !
            ! Cf/Cav  - with ZERO PROTECT
            if ( xrr(kk,nc) > 0 ) then
              xcfinr(kk,ics) = xcfinr(kk,ics)/xrr(kk,nc)
              !
              ! Done-expo flag. (within ZERO PROTECT loop)
              !  OPTION:  protect against wild values
              doneexpo(kk) = .true.
            end if
          end if
        end if
        !
        ! BACKWARD EULER SOLUTION FOR SINGLE SPECIES
        ! The true equation is =RP/(RL/XR); modified to avoid second divide.
        ! Use Back-Euler solution if option is selected, or as
        ! fallback if expo
        ! Note:  ZERO PROTECT above for rloss.
        !
        if ( .not. doneexpo(kk) ) then
          xrr(kk, nc)= xrp(kk,nc)*rpro(kk,nc)/rloss(kk,nc)
          xcfinr(kk,ics) = d_one
        end if
        !
        ! H2O2 SPECIAL OPTION (for AQUEOUS CO3 NONCONVERGENCE)
        ! MAXIMUM H2O2 = PRIOR H2O2 + ODDHSRC
        ! TO USE THIS, NEED TO SUM ODDHSRC (sources only) ALONG WITH ODDHSU
        !
        if ( c_tchem(ics) == '    H2O2') then
          if ( xrr(kk,nc) > oddhsrc(kk,1) + c_xcin(kk,ics)) then
            xrr(kk, nc) = oddhsrc(kk,1) + c_xcin(kk,ics)
            !
            ! ZERO PROTECT CORRECTION 52004 - cut, should never be zero
            !  --> if it  goes below 1e-30 can be zero!!!
            !           if(xrr(kk, 1) < 1.0e-30) xrr(kk,1)=1.0e-30
          end if
          !
          ! (THIS LINE REPLACED WITH SECTION IMMEDIATELY BELOW.)
          !  xrr(kk, nc) =(xrr(kk,nc) **0.70) * (xrp(kk,nc)) **0.30
        end if
        !
        ! ***DIFFICULT CONVERGENCE OPTION***
        ! GEOMETRIC MEAN (0.7 or 0.5) IF FLAGGED.
        !
        ! CONVERGENCE OPTION (June 2009): geom avg for first 1-3 iterations
        !   (but not for  NO3)   ==>NOT USED.
        !   (test with fort.87cohc3, cohc_u, rb_hg_rpai)
        !   (note: lself correction above is necessary)
        !   FLAG FOR:  H2O2   and SELF-REACTION
        !   (Self-reaction flag for SO5- +SO5- as well as H2O2).
        !
        ! If CO3, H2O2 here - slow convergence. Else, CO3-H2O2 small oscillati
        ! (changed:  CO3 linked to Hx. H2O2 remains here)
        !
        ! SPECIAL FIX FOR NO3, N2O5:  set lself=F;  these must not be solved
        !    with geometric avg - messes up NOx
        !
        !
        if ( c_tchem(ics) == '     NO3' .or. &
             c_tchem(ics) == '    N2O5' ) lself = .false.
        !
        ! CONVERGENCE OPTION - ADD THIS CONTROL if geom avg used for 1st iters
        ! geom avg
        if ( (c_iter > 1.or..not.c_lsts(ics) ) .and. &
              c_tchem(ics) /= '     NO3' .and. &
              c_tchem(ics) /= '    N2O5' ) then
          if ( (lself) .or. c_tchem(ics) == '    H2O2' ) then
            !
            ! 2009 CONVERGENCE OPTIONS  (best: c_iter <= 1 or skip)
            !  (if this is used, also add if control just above and endi
            !    &      .or.(c_iter <= 1                            )
            !    &      .or.(c_iter <= 1.and.c_icat(ncsol(1)) <= 8)
            !    &      .or.(c_iter <= 3.and..not.c_lsts(ics)       )
            ! always include
            ! FAILED OPTIONS
            !    *      .or. c_tchem(ics) == '     CLG'    !  fails
            !    *      .or. c_tchem(ics) == '    HCLG'
            !    *      .or. c_tchem(ics) == '    HCL2'
            !    *      .or. c_tchem(ics) == '     CLO'
            !    *      .or. c_tchem(ics) == '     NO3'    ! fails!! must skip
            !    *      .or. c_tchem(ics) == '    N2O5'    ! fails
            !    *      .or. c_tchem(ics) == '     CO3'
            !    *      .or. c_tchem(ncsol(1)) == '    HBRG'
            ! xrr(kk, nc) =(xrr(kk,nc) **0.85) * (xrp(kk,nc)) **0.15
            ! xrr(kk, nc) =(xrr(kk,nc) **0.70) * (xrp(kk,nc)) **0.30
            xrr(kk,nc) = (xrr(kk,nc) **0.50D0) * (xrp(kk,nc)) **0.50D0
          end if
        end if
        !
        ! DIFFICULT CONVERGENCE OPTION:  SETGEO IF FLAGGED
        !   FLAG FOR H2O2 - EITHER ABOVE (AUTOMATIC) OR HERE
        !   FAILS
        !
        ! DIFFICULT CONVERGENCE OPTION:   (ADD 30905)
        !  PROTECT AGAINST WILD OSCILLATION OF HNO3 FOR FIRST FEW ITERS
        !  (For NO3->HNO3->NOx error with high BR)
        !  (HNO3 icat 16;  may also use for all slow species icat=1)
        !
        if ( c_iter <= 3 .and. c_icat(ics) == 16 ) then
          if ( xrr(kk, nc) > 5.0D0*xrp(kk,nc) ) then
            xrr(kk, nc) =  5.0D0*xrp(kk,nc)
          end if
          if ( xrr(kk, nc) < 0.2D0*xrp(kk,nc) ) then
            xrr(kk, nc) = 0.2D0*xrp(kk,nc)
          end if
        end if
        !
        ! ADD SOLUTION TO PAIR SUM FOR PAIR GROUP
        !
        if ( nsolv >= 3 .and. nssolv(nsolv) == 1 ) then
          xrrpair(kk) = xrrpair(kk) + xrr(kk,nc)*c_pairfac(ics)
        end if
        !
        ! SUM XRRM FOR NORMALIZATION FOR MULTISOLVE
        !
        if ( nsol > 1 ) then
!          xrm(kk,is) =  xrm(kk,is) + xrp(kk,nc)*c_pairfac(ics) ! OLD
          xrrm(kk,is) = xrrm(kk,is) + xrr(kk,nc)
        end if
        !
        !    RESET XRR (TOTAL SPECIES SOLUTION) EQUAL TO NEW/OLD RATIO
        !  xrm(kk,is) = xrm(kk,is) + xrp(kk,nc)        ! original
        !  xrrm(kk,is) = xrrm(kk,is) + xrr(kk,nc)        ! original
        xrr(kk,nc) = xrr(kk,nc) / xrp(kk, nc)
        !
        ! PARTITION THE TOTAL SPECIES BETWEEN GAS-MASTER AND AQUEOUS
        ! USING THE SAME GAS-AQUEOUS RATIOS AS PRIOR.
        ! ALL CONCENTRATIONS (GAS AND AQUEOUS) ARE UPDATED BY THE RATIO XR/XR
        ! (PRIOR TIME STEP, XRIT, PRESERVED HERE)
        do  neq = 1 , (c_nequil(ics)+1)
          ic = ics
          if ( neq > 1 ) ic = c_ncequil(ics,(neq-1))
          xc(kk,ic) = xc(kk,ic) *xrr(kk,nc)
          c_xcout(kk,ic) = c_xcout(kk,ic)
        end do
        !
        ! END - PARTITION BETWEEN GAS AND AQUEOUS
        !
        ! UPDATE RLOSS, CPRO, AND XRP BASED ON NEW SPECIES CONCENTRATION
        ! AND ADD UPDATED CPRO TO RPRO FOR PRODUCTION TO PAIRED SUBSPECIES
        !
        ! (updated CPRO only used here for sub-species;
        ! since primary species formula assumes original CPRO, XRP
        ! but update done for all CPRO, to be consistent with updated XRP.
        ! Updated RLOSS, CPRO, XRP will be used in MULTISOLVE.)
        rloss(kk,nc) = rloss(kk,nc) * xrr(kk,nc)
        xrp(kk,nc) = xrp(kk,nc) * xrr(kk,nc)
        !
        ! UPDATE CPRO
        !
        do ncc = 1 , nsolv
          cpro(kk,nc,ncc) = cpro(kk,nc,ncc) * xrr(kk,nc)
        end do
        !
        ! UPDATE RPRO
        !
        if ( nsolv > nc ) then
          do ncc = (nc+1) , nsolv
            if ( nssolv(nc) == nssolv(ncc) ) then
              rpro(kk,ncc) = rpro(kk,ncc) + cpro(kk,nc,ncc)
            end if
          end do
        end if
        !
        ! GOT HERE - RPRO ABOVE, ADD WRITE?
        !
        !  FUTURE CHANGE - use this instead of UPDATE CPRO, RPRO just above.
        ! UPDATE CPRO AND RPRO for sub-species only
        !        if(nsolv > nc) then
        !          do ncc=(nc+1), nsolv
        !            if(nssolv(nc) == nssolv(ncc)) then
        ! c              do kk=1,c_kmax      ! kk vector loop
        !                cpro(kk,nc,ncc) = cpro(kk,nc,ncc) * xrr(kk,nc)
        !                rpro(kk,ncc) = rpro(kk,ncc) + cpro(kk,nc,ncc)
        ! c              end do               ! kk vector loop
        !            end if              !if(nssolv(nc) == nssolv(ncc))
        !          end do              !do ncc=(nc+1), nsolv
        !        end if           !if(nsolv > nc)
      end do
      ! ---------
      !  END LOOP TO SOLVE SPECIES CONCENTRATIONS- INTERNAL PAIRS
      ! ---------
      ! ---------
      !  NORMALIZATION:  UPDATE SPECIES FOR PAIR GROUP CONSERVATION OF MASS
      !                  OR FOR MULTISOLVE TO PRESERVE ORIGINAL PAIR GROUP SUM
      !  OPTION:  UPDATE SPECIES CONCENTRATIONS FOR PAIR GROUP CONSERVATION OF
      ! ---------
      !   PAIR GROUP SUM IS USED ONLY IF PAIR GROUP HAS >3 MEMBERS
      !     AND NO MULTISOLVE (since multisolve normalizes from prior )
      ! FROM THIS VERSION - FUTURE CHANGE (remove normalization in multisolve)
      if ( nsolv >= 3 .or. nsol > 1 ) then
        !
        !  DO FOR ALL PAIR GROUPS
        !   Indices:
        !     nc=chemsolve species counter number
        !     is = chemsolve group number
        !     icpair = pair group lead species
        !     ics = this species of pair group.
        !      (This loop is equivalent to:
        !           nc = 1, nsolv;   ics = ncsolv(nc) )
        do is = 1 , nsol
!          icc=ncsol(is)
!          icpair=ncsol(is)
          icpair = c_nppair(ncsol(is),2)
          !
          ! MULTISOLVE NORMALIZATION FACTOR:  XRM/XRRM
          ! (note, XRM different from XRPAIR:  XRPAIR includes PAIRFAC)
          !
          if ( nsol > 1 ) then
            calpha(kk) = d_one
            if ( xrrm(kk,is) > 0 ) then
              calpha(kk) = xrm(kk,is)/xrrm(kk,is)
            end if
            !
            ! PAIR GROUP NORMALIZATION FACTOR:
            ! Group sum = RP/RL.  Adj*XRR = XRP*RP/RL, Adj=(XRP/XRR)*(RP/RL)
            !
          else
            calpha(kk) = d_one
            if ( rlpair(kk,icpair) > 0 .and. xrrpair(kk) > 0 ) then
              calpha(kk) = (xrppair(kk)/xrrpair(kk)) * &
                          (rppair(kk,icpair)/rlpair(kk,icpair))
            end if
          end if
          !
          ! DO FOR ALL SPECIES IN PAIR GROUP
          !
          nc = 0
          do iss = 1 , (c_nppair(icpair,3)+1)
            nc = nc+1
            ics = icpair
!            ics = c_nppair(icpair,2)   ! no difference
            if ( iss > 1 ) then
              ics = c_nppair(icpair,iss+2)
            end if
!            icpair = c_nppair(ics,2)    ! should be equal to setting
            !
            !
            ! PAIR GROUP NORMALIZATION.
            ! PAIR ADJUSTMENT OPTION (2008): - SKIP if XRP=RPRO
            ! (zero production)
            do neq = 1 , (c_nequil(ics)+1)
              ic = ics
              if ( neq > 1 ) ic = c_ncequil(ics,(neq-1))
              ! TEST 2008 OPTION:  skip if XRP=RPRO
              if ( xrp(kk,nc) /= rpro(kk,nc) ) then
                xc(kk,ic) = xc(kk,ic) *calpha(kk)
              end if
            end do
          end do
          !
          ! END LOOP - DO FOR ALL SPECIES IN PAIR GROUP
          !
        end do
        !
        ! END LOOP - DO FOR ALL PAIR GROUPS.
        !
      end if
      ! ---------
      !  END UPDATE SPECIES CONCENTRATIONS FOR PAIR GROUP SUM CONSERVATION OF
      ! ---------
      ! ---------------------------------------------
      ! END CALCULATE SPECIES CONCENTRATIONS:  INTERNAL PAIRS
      ! ---------------------------------------------
      !
      ! ---------------------------------------------
      ! CALCULATE SPECIES CONCENTRATIONS: MULTISOLVE AND NOXSOLVE
      ! ---------------------------------------------
      !
      !   The algorithms:
      !    MULTISOLVE does a full back-euler matrix inversion
      !     using LINSLV/RESOLV to solve for >2 rapidly interacting species
      !    It solves for interacting species groups
      !     where each group may consist of a pair or pair chain.
      !    The groups are identified by input to CHEMSOLV (currently ic1, ic2,
      !    MULTISOLVE would preserve the previously calculated ratio
      !     among paired species or species pair chains
      !     and aqueous equilibria.
      !
      !    NOXSOLVE is a special solution for O3-NO-NO2
      !      based on the rapid exchange O3+NO<->NO2
      !      it uses preliminary calculations in CHEMSOLV and then calls NOXSO
      !
      !   FOR MULTISOLVE/NOX-SOLVE (STANDARD OPTION):
      !    PAIR SOLUTION is normalized to keep original within-pair sum unchan
      !    NOXSOLVE or MULTISOLVE changes pair sums (preserving internal parti
      !     based on interaction between pair groups.
      !    (implemented above in INTERNAL PAIRS )
      !
      ! NOTE:  pair group sums were set above in INTERNAL PAIRS
      !   xrm = prior pair group sum
      !   xrrm = updated pair group sum (not done yet)
      !   xrp = updated individual (gas+aq) species sum from PAIR SOLUTION.
      !
      !  ---------
      !  MULTISOLVE CONTROL:  DO ONLY IF MULTISOLVE SPECIES NONZERO
      !  ---------
      !      if(nssolv(nsolv) > 1) then
      if ( nsol > 1 ) then
        !
        ! ---------
        ! LOOP TO NORMALIZE INTERNAL PAIR SOLUTION FOR NOXSOLVE/MULTISOLVE
        ! ---------
        !  CUT - DONE ABOVE
        !  NOTE:  CHEMISTRY PRODUCTION, LOSS SET ABOVE
        !   NOT UPDATED AFTER PAIR SOLUTION.
        ! ---------
        !  CALL NOXSOLVE
        ! ---------
        ! NOXSOLVE solves for (O3, NO2, NO) in that order
        !    using O3+NO<->NO2 - Ox and NOx, and counted sourcnx, sourcox, etc.
        !    O3, NO2, etc. can also be part of a pair (for O3<->O1D)
        ! Solution is returned as XRRM(KK,IS)
        !
        !  (NOXSOLVE OPTION - EITHER CALL NOXSOLVE OR USE MULTISOLVE HERE.
        !    IF MULTISOLVE, CONVERGENCE IS SLOWER AND MORE DIFFICULT)
        !
        if ( ncsol(1) == c_no3 .and. ncsol(2) == c_nno2 ) then
          call noxsolve(ncsol(1),ncsol(2),ncsol(3))
          !
          ! CONVERT SOLUTION INTO NEW/OLD RATIO (XRRM/XRM)
          !
          do is = 1 , nsol
            xrrm(kk,is) = xrrm(kk,is)  / xrm(kk, is)
          end do
        else
          !
          !  ---------
          !  SET UP AND CALL MULTISOLVE
          !  ---------
          !   FUTURE DO:
          !      Control to replace MULTISOLVE with simple solution if CPM=0.
          !      Possible:  call RESOLV for easy convergence -
          !        must then save and replace reduced AX, IPA.
          !
          !  ENTER VECTOR AND MATRIX FOR AX=B BACK-EULER SOLUTION
          !  Equation: Xi = Xio + Pi(X)-Li(X) = Xio +Pi(Xp)-Li(Xp) +dRi/dxj
          !            Xi-Xip= (Xio + Pi(Xp)) - (Xip + Li(Xp)) + d(P-L)/dxj
          !            (I-dR/dxj)(Xj-Xjp) = ((Xio + Pi(Xp)) - (Xip + Li(Xp))
          !            [Xjp - (dR/dxj)*Xjp][Xj/Xjp-1] = (Xio + Pi(Xp)) - (Xip +
          !            AX     *    XX          =   BX
          !  XX = xrrm - 1
          !  BX = rpm - rlm
          !  AA(i,i) = rlm
          !  AA(i,j) = -cpm(j,i) for j /= i
          !
          !  SET ARRAY
          !
          do is = 1 , nsol
            bx(is) = rpm(kk,is) - rlm(kk,is)
            do iss = 1 , nsol
              ax(is,iss) = d_zero - cpm(kk,iss,is)
            end do
            ax(is,is) = rlm(kk,is)
          end do
          !
          ! CALL MULTISOLVE
          !
          call linslv(ax, bx, xx, nsol)
          !
          ! ENTER RESULT INTO XRRM (ratio adjustment)
          ! OPTION:  Limit size of change, also protect against zero.
          !
          do is = 1 , nsol
            xrrm(kk,is) = xx(is) + d_one
            if ( xrrm(kk,is) < 0.001D0 ) xrrm(kk,is) = 0.001D0
          end do
          !
          !  MULTISOLVE NORMALIZATION OPTION
          !
          cgamma(kk) = d_zero
          cbeta(kk) = d_zero
          do is = 1 , nsol
            ic = c_nppair(ncsol(is),2)
            cgamma(kk) = cgamma(kk) + xrm(kk,is)*c_multfac(ic)
            cbeta(kk) = cbeta(kk) + xrm(kk,is)*xrrm(kk,is)*c_multfac(ic)
          end do
          calpha(kk) = d_one
          if ( rlmulti(kk) > 0 .and. cbeta(kk) > 0 ) then
            calpha(kk) = (cgamma(kk)/cbeta(kk))*(rpmulti(kk)/rlmulti(kk))
          end if
          do is = 1 , nsol
            !
            ! OPTION HERE
            !
            xrrm(kk,is) = xrrm(kk,is) * calpha(kk)
          end do
        end if
        !  ---------
        !  END:  CALL NOXSOLVE/MULTISOLVE
        !  ---------
        !  -------
        !  ENTER NOXSOLVE/MULTISOLVE VALUES FOR PAIRS AND GAS/AQUEOUS
        !  -------
        !  Enter for all paired species; partition among gas/aqueous
        do nc = 1 , nsolv
          is = nssolv(nc)
          ics = ncsolv(nc)
          !
          ! UPDATE EACH PAIR SPECIES AND GAS/AQUEOUS SPECIES BY
          ! THE RATIO XRM/XR
          ! (This preserves prior partitioning among pair species and gas/aqueo
          !
          do neq = 1 , (c_nequil(ics)+1)
            ic = ics
            if ( neq > 1 ) ic = c_ncequil(ics,(neq-1))
            xc(kk,ic) = xc(kk,ic) *xrrm(kk,is)
          end do
          !
          !   END - PARTITION BETWEEN GAS AND AQUEOUS
          !
        end do
        !  -------
        !  END - ENTER NOXSOLVE/MULTISOLVE VALUES INTO PAIRS AND GAS/AQUEOUS
        !  -------
      end if
      !  ---------
      !  END MULTISOLVE CONTROL
      !  ---------
      !
      ! ---------------------------------------------
      ! END - CALCULATE SPECIES CONCENTRATIONS: MULTISOLVE AND NOXSOLVE
      ! ---------------------------------------------
      !
      ! ---------------------------------------------
      ! CALCULATE CHEMICAL PRODUCTION AND LOSSES FOR DOWN-CASCADE SPECIES.
      ! ---------------------------------------------
      !
      !  *** END IF-LOOP -1 TO SKIP LOSS AND XR CALC FOR SLOW-SPECIES 'PRE'.
      !          (END OF ONESOLVE  'POST' LOOP)
      !  CUT JANUARY 2005
      !  end if               ! if(ic2 /= -1.and.ic2 /= -3) then
      !  *** IF-LOOP -2 TO SKIP DOWN-CASCADE CALC. FOR SLOW-SPECIES 'POST'.
      !         (LOOP EXECUTED FOR full ONESOLVE AND FOR 'PRE')
      !  CUT JANUARY 2005
      !             if(ic2 /= -2) then
      !  ***
      !

      !  LOOP TO CALCULATE DOWN-CASCADE PRODUCTION
      !  LOOP TO  RECORD RRP, RRL AND ZERO RP, RL.

      !   Note:  brpro adds losses, cross-production to RP, RL (for record).
      !   Subsequent down-cascade production will be include in next iteration

      do nc = 1 , nsolv
        ics = ncsolv(nc)
        if ( c_nnrchem(ics) > 0 ) then
          do i = 1 , c_nnrchem(ics)
            nr = c_nrchem(ics,i)
            call brpro(nr)
          end do
        end if
      end do
      !
      ! LOOP TO CALL EXCORR
      !  EXCORR adjusts down-cascade RP, RL, PRONOX and  ODD-H SENSITIVITY
      !    for back-forth exchange reactions.
      !
      !  This avoids bad solution when down-cascade summed RP, RL are dominate
      !    by  a rapid back-forth reaction (MCO3<->PAN, HNO4<->NO2+HO2)
      !
      !  The exchange reaction impact would be exaggerated down-cascade
      !    by a straight solution, the down-cascade solution would not
      !    account for rapid back-forth exchange
      !      (= negative feedback for changed rates)
      !
      !  To avoid this, down-cascade RP, RL, PRONOX and ODD-H SENSITIVITY
      !   are reduced by an amount that reflects back-forth  response to
      !   change in rates.  Effective RP, RL is  reduced by linked back-reacti
      !
      do nc = 1 , nsolv
        ics = ncsolv(nc)
        if ( c_exspec(ics,1) /= 0 ) call excorr(ics)
      end do
      !
      !  LOOP TO  RECORD RP, RL
      !  NOTE:  RRP, RRL represent running sum through cascade.
      !         RP, RL represent final sum through cascade.
      !  ALSO- these are used in AQUASOLVE.
      !
      ! OPTION:  SUM GAS+AQUEOUS RP, RL INTO GAS-ONLY RP, RL?
      !   no, this may mess up AQUASOLVE.
      !
      ! OPTION:  THIS CAN BE MOVED ABOVE EXCORR
      !    so that it preserves full RP, RL for exchanged species
      !    (but may not record back-forth for  down-cascade)
      !
      do nc = 1 , nsolv
        ics = ncsolv(nc)
        do neq = 1 , (c_nequil(ics)+1)
          ic = ics
          if ( neq > 1 ) then
            ic = c_ncequil(ics,neq-1)
            !
            ! OPTION: RECORD SUM OF GAS+AQUEOUS PRODUCTION IN GAS-MASTER
            !            c_rp(kk,ics) = c_rp(kk,ics) + rrp(kk,ic)
            !            c_rl(kk,ics) = c_rl(kk,ics) + rrl(kk,ic)
            !
          end if
          c_rp(kk, ic) = rrp(kk,ic)
          c_rl(kk,ic) = rrl(kk,ic)
        end do
      end do
      !
      ! LOOP TO ZERO RRP, RRL
      ! ZERO RRP HERE so that running sum for next iteration includes
      ! down-cascade production.
      !
      do nc = 1 , nsolv
        ics = ncsolv(nc)
        icpair = c_nppair(ics,2)
        rppair(kk,icpair) = d_zero
        rlpair(kk,icpair) = d_zero
        do neq = 1 , (c_nequil(ics)+1)
          ic = ics
          if ( neq > 1 ) then
            ic = c_ncequil(ics,neq-1)
          end if
          rrp(kk,ic) = d_zero
          rrl(kk,ic) = d_zero
        end do
      end do
      !
      !  *** END IF-LOOP -2 TO SKIP DOWN-CASCADE CALC FOR SLOW SPECIES 'POST'.
      !        (END OF LOOP FOR full ONESOLVE AND 'PRE')
      !  CUT JAUNARY 2005
      !             end if                  ! if(ic2 /= -2) then
      !  ***
      !
      ! ----------------------------------------------------------------
      ! END CHEMSOLVE.
      ! ----------------------------------------------------------------

  121 format(/,'MAJOR ERROR IN CHEMSOLVE: ',             &
              'RLOSS = 0 FROM CPRO ADJUSTMENT.',/,        &
              ' ic1, nc, ncc=  ', a8,3i4,/,               &
              ' rloss, calpha, cpro-nc, cpro-ncc, rl-ncc=',/, (5(1pe10.3)))
                   rloss(kk,nc) = rloss(kk,nc) + calpha(kk)
      !
    end subroutine chemsolve
!
! -------------------------------------------
!
! BRREAC CALCULATES REACTION RATE FOR THE SPECIFIED REACTION
! AND STORES IT IN RR().
! FOR AQUEOUS REACTIONS, IT CALCULATES RR IN GAS UNITS.
! CALLED BY CHEMSOLVE, NOXSOLVE, EXCORR, etc.
!
! Input:  nr  (reaction number)
! Output:  c_rr (rate of reaction, molec/cm3/step)
!
! Called by:
!     chemsolve
!     noxsolve
!     excorr
!
! Calls to:  none.
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------

    subroutine brreac(nr)
      implicit none
      integer , intent(in) :: nr
!
      kk=1
      !
      ! NOTE - ASSUME THAT ICR1>0 THROUGHOUT.
      !
      icr1 = c_reactant(nr,1)
      icr2 = c_reactant(nr,2)
      !
      !  CALCULATE REACTION RATE.
      !
      c_rr(kk,nr) = c_time*ratek(kk,nr)*xc(kk,icr1)
      if ( icr2 > 0 ) then
        c_rr(kk,nr) = c_rr(kk,nr) * xc(kk,icr2)
      end if
      !
      ! CONVERT AQUEOUS REACTIONS TO GAS UNITS
      ! AQUEOUS REACTION IS IDENTIFIED BY REACTANT NUMBER
      ! DIFFERENT FROM ITS GAS-MASTER SPECIES.
      if ( icr1 /= c_npequil(icr1) ) then
        c_rr(kk,nr) = c_rr(kk,nr)*c_h2oliq(kk)*avogadrl
      end if
    end subroutine brreac
!
! ----------------------------------------
!
! This calculates reaction rates (RR),
!  and adds to species production (RP) and loss (RL)
!  for a given reaction (nr).
!
! It also sums species loss rates (rloss) and radical sums
!   (oddhsum, etc) for the NOx and HOx solution.
!
! It also adjusts pair sums (RPPAIR, RLPAIR)
!
! This is called from chemsolve (solution for individual species)
!  before the calculation of species concentrations,
!  for reactions associated with species production
!
! And is called after the calculation of the species concentration
!  to get the complete sum of species losses and production
!  for all reactants and products of reactions linked to the species.
!
!  (The resulting species production/loss is used in the solution
!    for species concentrations.
!   Note, RP and RL are set to zero when the species concentration
!    is calculated, so that species RP and RL are sums for
!    all reactions over a full cycle of present+past iterations.)
!
! Inputs:  reaction number (nr)
!
! Outputs:  Reaction rate (RR) based on species concentrations;
!           Species production and loss (RP, RL)
!           Chemistry production sums (oddhsum, etc.)
!
! Called by:    chemsolve.
!
! Calls to:  brreac
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
      subroutine brpro(nr)
!
      implicit none
      integer , intent(in) :: nr
      kk=1
      !
      ! PRELIM:  CALCULATE REACTION RATE  (alt: cut and invoke elsewhere?)
      !
      call brreac(nr)
      !
      ! PAIR INDEX:  IDENTIFY PAIR GROUP OF NON-KEY-SPECIES REACTANT ONLY.
      ! REQUIRES TWO REACTANTS, WHERE 1ST REACTANT IS KEY SPECIES (nicreac)
      icpair = 0
      ic1 = c_npequil(c_reactant(nr,1))
      ic2 = 0
      if ( c_reactant(nr,2) > 0 ) then
        ic2 = c_npequil(c_reactant(nr,2))
        if ( c_nicreac(nr) == ic1 ) then
          if ( c_nppair(ic2,2) /= c_nppair(ic1,2) ) then
            icpair = c_nppair(ic2,2)
          end if
        end if
        !
        ! IF 2ND REACTANT IS KEY SPECIES - can't handle, would make adj. below m
        !           if(c_nicreac(nr) == ic2) then
        !             if(c_nppair(ic2,2) /= c_nppair(ic1,2) )
        !    *             i       icpair = c_nppair(ic1,2)
        !           end if                !if(c_nicreac(nr) == ic1) then
      end if
      !
      ! ENTER RATE AS LOSS (RL) FOR REACTANTS:
      !
      ic = c_reactant(nr,1)
      rrl(kk,ic) = rrl(kk,ic) + c_rr(kk,nr)
      ic = c_reactant(nr,2)
      if ( ic > 0 ) then
        rrl(kk,ic) = rrl(kk,ic) + c_rr(kk,nr)
      end if
      !
      ! ENTER RATE AS PRODUCT (RP) FOR PRODUCTS.
      ! NOTE:  ASSUME IC>0 FOR PRODUCTS IDENTIFIED BY NNPRO.
      !
      if ( c_nnpro(nr) > 0 ) then
        do n = 1 , c_nnpro(nr)
          ic = c_product(nr,n)
          rrp(kk,ic) = rrp(kk,ic) + c_rr(kk,nr)*c_stoich(nr,n)
          !
          ! OPTION - IF PRODUCT EQUALS REACTANT, SUBTRACT FROM RP AND RL
          ! (This is important for OHL and CLOH in new representation
          !  of reaction CL-+OHL=>CLOH (+OHL)
          !
          ! ( Prior:  This only matters for R1O2+R2O2->R1O2+products.
          ! Reaction is entered twice, once for R1O2 products, once for R2O2.
          ! So RR is entered as R1O2 loss and product - possible error.
          !  For all other reactions, STOILOSS would correct.
          !  Since it converges anyway, no problem.
          !
          if ( ic == c_reactant(nr,1) .or. ic == c_reactant(nr,2) ) then
            stoicx = c_stoich(nr,n)
            if ( stoicx > d_one ) stoicx = d_one
            rrl(kk,ic) = rrl(kk,ic) - c_rr(kk,nr)*stoicx
            rrp(kk,ic) = rrp(kk,ic) - c_rr(kk,nr)*stoicx
          else
            !
            ! END  PRODUCT=REACTANT OPTION.
            !   2008 CORRECTION:  - PAIR ADJUSTMENT JUST BELOW IS SKIPPED
            !                           IF PRODUCT=REACTANT OPTION IS DONE.
            ! since PRODUCT=REACTANT adjustment does the same thi
            ! PAIR ADDITION:  FOR DOWN-CASCADE 2nd REACTANT ONLY
            !         IF REACTANT AND PRODUCT ARE IN SAME PAIR GROUP,
            !            THEN SUBTRACT REDUNDANT PRO/LOSS  FROM RPPAIR, RLPAIR
            ! (for same pair group as key species, calc in CHEMSOLVE)
            !
            !  PAIR INDICES (ic2, icpair) IDENTIFIED ABOVE
            !
            ! KEY SPECIES IS c_nicreac(nr).
            !   PAIR SPECIES IS icpair = c_nppair(c_npequil(ic1), 2)
            !
            if ( icpair > 0 ) then
              if ( c_nppair(c_npequil(ic),2) == icpair ) then
                stoicx = c_pairfac(c_npequil(ic)) * c_stoich(nr,n)
                if ( stoicx > c_pairfac(ic2) ) stoicx = c_pairfac(ic2)
                rppair(kk,icpair) = rppair(kk,icpair) - c_rr(kk,nr)*stoicx
                rlpair(kk,icpair) = rlpair(kk,icpair) - c_rr(kk,nr)*stoicx
              end if
            end if
            !
            ! END PAIR ADDITION FOR PRODUCT=REACTANT PAIRS
            !
          end if
          !
          ! END PRODUCT=REACTION OPTION CONTROL,
          ! which includes PAIR ADJUSTMENT (2008 CORRECTION)
          !
        end do
      end if
      !
      ! SUM ODD-H NET AND ODD NITROGEN PRODUCTION (ONLY)
      ! INCLUDES IDENTIFICATION AND SUM OF dHX/DOH FOR THE REACTION.
      ! (SENSHX CALC MUST BE IN ITERATIVE LOOP - SENHCAT CHANGES.)
      senshx(kk,1) = d_zero
      senshx(kk,2) = d_zero
      do n = 1 , 2
        ic = c_reactant(nr,n)
        if ( ic > 0 ) then
          ic = c_npequil(ic)
          if ( c_icat(ic) > 0 ) then
            senshx(kk,1) = senshx(kk,1) + senhcat(kk, c_icat(ic))
            if ( c_icat(ic) /= 3 ) then
              senshx(kk,2) = senshx(kk,2) + senhcat(kk, c_icat(ic))
            end if
          end if
        end if
      end do
      !
      ! OPTION: double sensitivity for PARAMETERIZED RO2-RO2 reactions
      !
      if ( c_nrk(nr) == -13 ) then
        senshx(kk,1) = d_two*senshx(kk,1)
      end if
      do i = 1 , 2
        oddhsum(kk,i) = oddhsum(kk,i) + c_rr(kk,nr)*c_oddhx(nr,i)
        oddhdel(kk,i) = oddhdel(kk,i) + c_rr(kk,nr)*c_oddhx(nr,i)*senshx(kk,i)
        if ( c_oddhx(nr,i) > 0 ) then
          oddhsrc(kk,i) = oddhsrc(kk,i) + c_rr(kk,nr)*c_oddhx(nr,i)
        end if
        if ( ic1 == c_noh ) then
          oddhloh(kk,i) = oddhloh(kk,i) - c_rr(kk,nr)*c_oddhx(nr,i)
        end if
        if ( ic1 == c_nho2 ) then
          oddhlho2(kk,i) = oddhlho2(kk,i) - c_rr(kk,nr)*c_oddhx(nr,i)
        end if
        if ( ic2 > 0 ) then
          if ( ic2 == c_noh ) then
            oddhloh(kk,i) = oddhloh(kk,i) - c_rr(kk,nr)*c_oddhx(nr,i)
          end if
          if ( ic2 == c_nho2 ) then
            oddhlho2(kk,i) = oddhlho2(kk,i) - c_rr(kk,nr)*c_oddhx(nr,i)
          end if
        end if
      end do
!
      if ( c_pronox(nr) > 0 ) then
        sourcnx(kk) = sourcnx(kk) + c_rr(kk,nr)*c_pronox(nr)
      end if
      if ( c_pronox(nr) < 0 ) then
        sinknx(kk) =  sinknx(kk) + c_rr(kk,nr)*c_pronox(nr)
      end if
    end subroutine brpro
!
! ----------------------------------------
!
! EXCORR IS A CORRECTION FOR "EXCHANGE" (BACK-FORTH) REACTIONS
! (e.g. PAN<->MCO3, HNO4<->HO2, etc.)
! It goes through the back-forth reactions (expro, exloss)
! associated with the exchange species and "undoes" part of them.
!
! For a SINGLE exchange species (eg HNO4<->NOx, Hx)
! it undoes all except NET FORWARD or BACKWARD.
! 1998 CHANGE:  it undoes based on species lifetime.

! For TWO LINKED SPECIES (PAN-MCO3) it calculates allowable
!  back-forth from PRIOR PAN->MCO3->PAN, and prior MCO3->PAN->MCO3
!  and undoes the rest.
!
! NOTE: SOLUTION FOR TWO LINKED SPECIES is possible
!  only if species are solved simultaneously in same reaction pair
!  or (FUTURE DO) in MULTISOLVE.
!  The 2nd species is identified by EXSPEC - set in QUADINIT.
!
! The UNDO part undoes all RP, RL for down-cascade species,
! including ODDNSUM and ODDHDEL.
!  (But not ODDHSUM, that is unaffected by back-forth exchange.)
!
!  UPDATED FOR AQUEOUS TWOSOLVE.  WATCH OUT FOR PAN, CLOH.
!
!  2005 CHANGE (OPTION):  updates loss (rrl) for exchanged species also.
!   normally, excorr comes after species is solved for, so no matter
!   but it screws up when down-cascade species
!    (e.g. HG3A <=>HG# + SO3- special equilibrium)
!
! Inputs:  species number (ic1) for either individual species
!            or 1st of species pair.
!          Rate of reactions (RR), species production and loss (RP, RL)
!
! Outputs:  Updated RP, RL for subsequent products
!             from back-forth exchange reactions;
!             also updated sums (ODDHDEL)
!
! Called by:    chemsolve (after species solution)
!
! Calls to:  None.
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
!
     subroutine excorr(ic1)
!
      implicit none
      integer , intent(in) :: ic1
!
! LOCAL VARIABLES
!
!  ratex(kk)  Rate of back-forth exchange reactions (molec/cm3/timestep)
!              Summed for all exchange reactions affecting species.
!              Equal to minimum of forward and backward rates
!              (more complex for paired species)
!             Preliminary value, used to set 'undo'.
!
! undo (kk)  Rate of back-forth exchange reaction to be 'undone'
!              Adjusted for each  indiv. reactc_ion(molec/cm3/timestep)
!
      ! Exchange reaction sum  mol/cm3
      real(dp) :: ratex(c_kvec)
      ! Exchange undo rate mol/cm3
      real(dp) :: undo(c_kvec)
      ! Index for number of ex vars
      integer :: nex
!
      kk = 1
!
      if ( ic1 == 0 ) return
      if ( c_exspec(ic1,1) == 0 ) return
      !
      ! MAIN LOOP FOR MULTIPLE EXCHANGED SPECIES
      !   Exchange reactions are possible with different species partners.
      !   Sum up all exchange reactions with the same species
      !          and solve  simultaneously
      do nex = 1 , 20
        if ( c_exspec(ic1,nex) == 0 ) return
        icx = c_exspec(ic1,nex)
        !
        ! ZERO
        !
        cpro(kk,1,2) = d_zero
        cpro(kk,2,1) = d_zero
        ratex(kk) = d_zero
        prior(kk) = d_zero
        !
        ! SUM EXCHANGE REACTIONS  AS CPRO
        ! (multiple reactions are allowed for single exchange)
        !
        do  ii = 1 , 5
          nr = c_exloss(ic1,nex,ii)
          if ( nr > 0 ) then
            cpro(kk,1,2) = cpro(kk,1,2) +  c_rr(kk,nr)
          end if
          np = c_expro(ic1,nex,ii)
          if ( np > 0 ) then
            cpro(kk,2,1) = cpro(kk,2,1) + c_rr(kk,np)
          end if
        end do
        !
        ! ESTABLISH RATEX FOR SINGLE-REACTION CASE (e.g. HNO4<->Nx,Hx)
        ! EQUAL TO MIN(RPRO,RLOSS)
        ! 1998 CHANGE:  MULTIPLIED BY FACTOR FOR EXCHANGE LIFETIME
        ! 2000 CHANGE:  IF BETA=0, DO ALSO FOR CASE WHERE IC2 NONZERO
        !
        if ( icx <= 0 ) then
          ratex(kk) = cpro(kk,1,2)
          if ( ratex(kk) > cpro(kk,2,1) ) ratex(kk)=cpro(kk,2,1)
          if ( cpro(kk,1,2) > 0 ) then
            ratex(kk) = ratex(kk)*(cpro(kk,1,2)/(cpro(kk,1,2)+xc(kk,ic1)))
          end if
        end if
        !
        ! IF IC2>0 AND GAS, TWO LINKED SPECIES.
        ! CALCULATE THE SIZE OF THE REDUNDENT EXCHANGE (RATEX):
        !
        !  (REDUNDENT EXCHANGE DERIVED FROM TWOSOLVE SOLUTION ABOVE.
        !   REDUNDENT EX = CP(2->1) FROM SOLUTION WITH ZERO XR2 SOURCE (RP2=0)
        !   PLUS CP(1->2) FROM MATRIX SOLUTION WITH ZERO INITIAL XR(1) (RP1=0).
        !   FROM TWOSOLVE SOLUTION, ABOVE,
        !       PARTIAL X2=XP2*(CP2*RP1)/(RL1*RL2 - CP1*CP2)
        !       REDUNDANT 2->1EX = CP1'*XP2*CP2*RP1/denom.
        !                        =CP1*CP2*RP1/denom.  SAME FOR 1->2EX.
        ! MODIFIED TO ALLOW AQUEOUS, SPECIAL EQUILIBRIUM SPECIES (2000)
        !
        if ( icx > 0 ) then
          !
          ! ZERO AND SUM RPRO, RLOSS, XRP FOR TWO EXCHANGE SPECIES
          !
          do is = 1 , 2
            rloss(kk,is) = d_zero
            rpro(kk,is)  = d_zero
            xrp(kk,is)   = d_zero
            ics = c_npequil(ic1)
            if ( is == 2 ) ics = c_npequil(icx)
            do neq = 1 , (c_nequil(ics)+1)
              ic = ics
              if ( neq > 1 ) ic = c_ncequil(ics,(neq-1))
              calpha(kk) = d_one
              if ( neq > 1 ) calpha(kk)= c_h2oliq(kk)*avogadrl
              rloss(kk,is) = rloss(kk,is) + rrl(kk,ic)
              rpro(kk,is) = rpro(kk,is) + rrp(kk,ic)
              xrp(kk,is) = xrp(kk,is) + xc(kk,ic)*calpha(kk)
              prior(kk) = prior(kk) + c_xcin(kk,ic)*calpha(kk)
            end do
          end do
          !
          ! POSSIBLE ERROR, MUST SUBTRACT CPRO FROM RPRO - ELSE ITS INCLUDED!
          ! SOLUTION FOR RATEX.
          !   (If denominator is zero, use single-species solution instead)
          cbeta(kk) = (rloss(kk,1)+xrp(kk,1))*(rloss(kk,2)+xrp(kk,2)) - &
                      cpro(kk,1,2)*cpro(kk,2,1)
          if ( cbeta(kk) > 0 ) then
            ratex(kk) = cpro(kk,1,2)*(cpro(kk,2,1)*(rpro(kk,1) + &
                        rpro(kk,2) + prior(kk))/cbeta(kk) )
          else
            ratex(kk) = cpro(kk,1,2)
            if ( ratex(kk) > cpro(kk,2,1) ) ratex(kk) = cpro(kk,2,1)
            if ( cpro(kk,1,2) > d_zero ) then
              ratex(kk) = ratex(kk)*(cpro(kk,1,2)/(cpro(kk,1,2)+xc(kk,ic1)))
            end if
          end if
        end if
        !
        ! PROTECT AGAINST RATEX>EXCHANGE, RPRO FOR TWO-REACTION CASE
        !        do kk=1,c_kmax        ! kk vector loop
        !
        if ( ratex(kk) > cpro(kk,1,2) ) ratex(kk) = cpro(kk,1,2)
        if ( ratex(kk) > cpro(kk,2,1) ) ratex(kk) = cpro(kk,2,1)
        !
        ! UNDO LOOP FOR INDIVIDUAL REACTIONS;
        ! UNDO EXCHANGE PRODUCTION AND LOSS REACTIONS IN PROPORTION TO SIZES
        ! SO THAT TOTAL UNDO IN EACH DIRECTION EQUALS RATEX.
        !
        do ii = 1 , 5
          do i = 1 , 2
            if ( i == 1 ) nr = c_exloss(ic1,nex,ii)
            if ( i == 2 ) nr = c_expro(ic1,nex,ii)
            !
            ! SET UNDO FOR INDIVIDUAL REACTION.
            ! NOTE DIVIDE-BY-ZERO DANGER.  RPRO, RLOSS SHOULD NEVER BE ZERO.
            !
            if ( nr > 0 ) then
              if ( i == 1 ) then
                undo(kk) = d_zero
                if ( cpro(kk,1,2) > 0 ) then
                  undo(kk) = ratex(kk)*c_rr(kk,nr)/cpro(kk,1,2)
                end if
              else
                undo(kk) = d_zero
                if ( cpro(kk,2,1) > 0 ) then
                  undo(kk) = ratex(kk)*c_rr(kk,nr)/cpro(kk,2,1)
                end if
              end if
              !
              ! UNDO LOSS FROM REACTANTS
              ! 2000 MODIFICATION: UNDO FOR REAL REACTANTS, NOT
              ! GAS EQUIVALENTS (npequil comment out)
              do n = 1 , 2
                ic = c_reactant(nr,n)
                if ( ic > 0 ) then
                  ! 2005 CHANGE OPTION:  DON'T SKIP FOR MAIN SPECIES.
                  !   CAUSES ERROR IF MAIN SPECIES IS OFF CASCADE.
                  if ( ic /= ic1 .and. ic /= icx ) then
                    rrl(kk,ic) = rrl(kk,ic) - undo(kk)
                  end if
                end if
              end do
              !
              !  UNDO PRODUCTION.
              !
              if ( c_nnpro(nr) > 0 ) then
                do n = 1 , c_nnpro(nr)
                  ic = c_product(nr,n)
                  if ( ic > 0 ) then
                    ! 2005 CHANGE OPTION:  DON'T SKIP FOR MAIN SPECIES.
                    !   CAUSES ERROR IF MAIN SPECIES IS OFF CASCADE.
                    if ( ic /= ic1 .and. ic /= icx ) then
                      rrp(kk,ic) = rrp(kk,ic) - c_stoich(nr,n)*undo(kk)
                    end if
                  end if
                end do
              end if
              !
              ! UNDO ODD NITROGEN PRODUCTION AND ODD-H SENSITIVITY
              !
              senshx(kk,1) = d_zero
              senshx(kk,2) = d_zero
              do n = 1 , 2
                ic = c_reactant(nr,n)
                if ( ic > 0 ) then
                  ic = c_npequil(ic)
                  if ( c_icat(ic) > 0 ) then
                    senshx(kk,1) = senshx(kk,1) + senhcat(kk,c_icat(ic))
                    if ( c_icat(ic) /= 3 ) then
                      senshx(kk,2) = senshx(kk,2) + senhcat(kk, c_icat(ic))
                    end if
                  end if
                end if
              end do
              do n = 1 , 2
                oddhdel(kk,n) = oddhdel(kk,n) - &
                        c_oddhx(nr,n)*undo(kk)*senshx(kk,n)
              end do
              if ( c_pronox(nr) > 0 ) then
                sourcnx(kk) = sourcnx(kk)-c_pronox(nr)*undo(kk)
              end if
              if ( c_pronox(nr) < 0 ) then
                sinknx(kk) =  sinknx(kk)-c_pronox(nr)*undo(kk)
              end if
            end if
          end do
        end do
      !
      ! END UNDO LOOP FOR INDIVIDUAL REACTIONS
      !
      end do
     !
     ! END MAIN LOOP  FOR  MULTIPLE EXCHANGED SPECIES
     !
     end subroutine excorr
!
! --------------------------------
!
! Special solution for concentrations of O3, NO2 and NO.
!
! This uses reaction rates, sums of RLOSS, gas-aqueous sums
!   and cross-production (CPRO)  from CHEMSOLVE.
!
! FUTURE OPTION:  Flag to keep NOx at pre-set value,
!                  while adjusting NO and NO2
!                 (fails so far.  Set SOURCNX = SINKNX?)
!
! NOTE:  THIS SUBROUTINE MAY NEED TO USE R12: OH+NO2->HNO3 SPECIAL.
!
!
! Inputs:    Species numbers (ic) for O3, NO2, NO.
!            Also uses:
!            Species sums for production, loss and cross-production
!            (xrpm, rlm, rpm, cpm) from MULTISOLVE part of CHEMSOLVE.
!            NOx sums (sourcnx, sinknx) from full program.
!
! Outputs:  xrrm:  concentrations for O3, NO, NO2.
!
!           (used in MULTISOLVE part of CHEMSOLVE
!             to create updated species solution array
!             with gas/aqueous partitioning)
!
! Called by:    chemsolve (as alternative species solution)
!
! Calls to:  None.
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
!

    subroutine noxsolve(ic1,ic2,ic3)
!
      implicit none
      integer , intent(in) :: ic1 , ic2 , ic3
!
! LOCAL VARIABLES
! xnox     NOx concentration (NO+NO2) molec/cm3
! xox      Ox concentration  (O3+NO2) molec/cm3

      ! NOx, molec/cm3
      real(dp) :: xnox(c_kvec)
      !  Ox, molec/cm3
      real(dp) :: xox(c_kvec)

      ! Prior Ox for write
      real(dp) :: oxox
      ! Ox source for write
      real(dp) :: sourcox
      ! Ox sink  for write
      real(dp) :: sinkox

      ! Quadratic parameter for write
      real(dp) :: xk2
      ! Quadratic parameter for write
      real(dp) :: xk2n
      ! Quadratic parameter for write
      real(dp) :: xrln
      ! Quadratic parameter for write
      real(dp) :: xjn
      ! Quadratic parameter for write
      real(dp) :: xjnx
      ! Quadratic parameter for write
      real(dp) :: xrpn
      ! Sum pans for write
      real(dp) :: pansum
      ! Sum hno2 hno4 for write
      real(dp) :: hnosum
!
      kk=1
      !
      ! FIRST:  SOLVE FOR NOX.  USE ODDNSUM = CHEM PRODUCTION OF NOX
      !  (excluding NO-NO2 conversions and NOx sinks).
      ! SOURCNX = NOX SOURCE WITH PRIOR.
      ! SINKNX = NOX SINKS, IGNORING NO-NO2 BACK-FORTH. (ALSO W/PRIOR)
      ! XNOXN = SOURCE/(SINK/XNOX)
      !
      sourcnx(kk) = sourcnx(kk) + c_xcin(kk,ic2) + c_xcin(kk,ic3)
      sinknx(kk) = xrm(kk,2) + xrm(kk,3) - sinknx(kk)
      xnox(kk) = (xrm(kk,2) + xrm(kk,3))*sourcnx(kk)/sinknx(kk)
      !
      ! OPTION:  INSERT THIS LINE FOR PRESET NOx CONCENTRATION
      !       xnox(kk) =   c_xcin(kk,ic2)+  c_xcin(kk,ic3)
      !
      ! SENHCAT FOR NOX.  (PRIOR, dHx/dNOx dNOx/dOH product added to SENH.
      !                    REPLACED, just include dNOx/dOH in HxDELTA sum.)
      ! THIS IS FROM THE NOX EQUATION:  NOXPRO=a+bOH
      ! dlnNOx/dlnOH = -[(OH+NO2) + sencatHO2*otherNOx losses]/ total loss;
      !                                        total loss includes XRP(NOx).
      !
      ! OPTIONS: USE REACTION R12 (OH+NO2->HNO3 EXPLICITLY.
      !
      ! OPTION WITH R-12
      !       senhcat(kk,12) = 0.- (  c_rr(kk,12)
      !    *  + senhcat(kk,10)
      !    *       *(sinknx(kk)-c_rr(kk,12)- xrm(kk,  2)- xrm(kk,  3))
      !    *                  )/sinknx(kk)
      ! OPTION WITHOUT EXPLICIT R-12
      !
      senhcat(kk,12) = d_zero - (senhcat(kk,10)*(sinknx(kk)-xrm(kk,2) - &
                           xrm(kk,3)))/sinknx(kk)
      !
      ! END OPTION
      !
      senhcat(kk,13) = senhcat(kk,12)
      !
      ! STEADY-STATE.  THERE IS NO STEADY-STATE OPTION FOR NOX, ->INFINITE.
      ! BUT ISTS(NOX) CAN BE USED TO SET NOX EQUAL TO CONSTANT VALUE.
      !  (FAILS 1999 - MUST USE HARD-WIRE OPTION, ABOVE.)
      !
      if ( c_lsts(ic3) ) then
        xnox(kk) = xrm(kk,2) + xrm(kk,3)
        senhcat(kk,12) = d_zero
        senhcat(kk,13) = d_zero
      end if
      !
      ! ODD OXYGEN (XOX):  SOLVE USING BACK-EULER: OX = RP/(RL/XRP).
      ! RP AND RL FOR ODD OXYGEN EQUAL THE SUM OF RP FOR O3 AND NO2
      ! WITH CROSS-CONVERSIONS (O3->NO2 AND NO2->O3) REMOVED.
      ! IN CHEMSOLVE(O3,NO2,NO), CROSS-PRODUCTION IS ADDED TO CPRO, NOT RPRO,
      !  BUT IT IS ALSO INCLUDED IN RLOSS.  SO RLOSS IS CORRECTED.
      !
      ! PRELIMINARY:  ADD NO-to-NO2 conversions to RPM(NO2)
      !  These are equal to CPM(NO->NO2) minus CPM(O3->NO2)
      !  They count as Ox production
      !  They also count in the solution for NO,
      !      which separates out O3+NO<->NO2 only.
      !
      ! CHANGE:  WITH RPM INCLUDING CPM:
      ! PRELIMINARY:  REMOVE O3+NO<=>NO2 from RPM.
      !     (Removal is CPM(1->2) and (2->1).  NOT (3->1).
      !       In original version, 3->1 is zero.  Why?? - error)
      !   RPM includes all CPM. Other NO-to-NO2 count for Ox and NO solutions.
      !
      ! ORIGINAL
      !       rpm(kk,2) = rpm(kk,2) + cpm(kk,3,2) - cpm(kk,1,2)
      ! OPTION
      !
      rpm(kk,2) = rpm(kk,2) - cpm(kk,1,2)
      rpm(kk,1) = rpm(kk,1) - cpm(kk,2,1)
      rpm(kk,3) = rpm(kk,3) - cpm(kk,2,1)
      !
      ! TEMPORARY BUG FIX:  O3+NO->NO2 COUNTED DOUBLE IN ORIGINAL ALGORITHM
      !   CORRECT HERE.  NOXBUG - delete when CPM is correct
      ! rpm(kk,2) = rpm(kk,2) - cpm(kk,1,2)
      !
      xox(kk) = (xrm(kk,1) + xrm(kk,2))*(rpm(kk,1)+rpm(kk,2)) / &
                (rlm(kk,1) + rlm(kk,2) - cpm(kk,1,2) - cpm(kk,2,1))
      !
      ! OPTION:  INSERT THIS LINE FOR PRESET Ox CONCENTRATION
      !  NOTE - ALSO CHANGE 'PRESET' OPTION IN presolve
      ! xox(kk) =   c_xcin(kk,ic1)+  c_xcin(kk,ic1)
      ! --------------------------------------
      ! NO:  SOLVE BACKWARD EULER EQUATION FOR NO AS FUNCTION OF OX, NOx.
      !
      ! THE EQUATION:  RP" + r1*NO2 = RL" + r2*NO*O3  BECOMES
      !  RP" + r1*(NOx-NO) = RL"*NO + r2*NO*(Ox-NOx+NO).
      !  WHERE RP", RL" = production, loss of NO (including NO->NO2)
      !      but without R1:  NO2->NO+O3 and R2: NO+O3->NO2.
      !
      !  THIS PROVIDES A QUADRATIC FOR NO, SINCE Ox, NOx VARY SLOWLY.
      !  SUBSTITUTE CPRO(O3->NO2) = r1*NO2prior;
      !             CPRO(NO2->O3)=r2*O3p*NOp, w/ aq sums.
      !
      !  -> RP" + CP21(NOx-NO)/NO2p = RL"*NO + CP12*[(Ox-NOx+NO)/O3p]*[NO/NOp]
      !
      !  -> RP" + CP21(NOx/NO2p) = NO* [RL" + CP21/NO2p + CP12*(Ox-NOx)/(O3p*N
      !                          +NO^2 *[CP12/(O3p*NOp)]
      !
      !  JANUARY 2005: SOLVE THIS QUADRATIC WITH THE FOLLOWING SUBSTITUTIONS:
      !
      !   RP" = rpm3  = p(NO), excluding NO2+hv.  (NO2+hv goes to cpro, not rp
      !   RL" = (rlm3 - cpm(1,2).  rlm3 includes all losses of NO.
      !          subtract cpm12 (not cpm32) because cpm12 is O3+NO->NO2 exactl
      !                   (cpm32 includes NO+HO2=>NO2)
      !  CP21, CP12 = cpm(1,2), cpm(2,1)
      !
      !   NOTE:  DIVIDES ADD COMPUTER TIME BUT PREVENT OVERFLOWS
      ! --------------------------------------
      calpha(kk) =  cpm(kk,1,2)/( xrm(kk,3)* xrm(kk,1))
      cgamma(kk) = rpm(kk,3) + (cpm(kk,2,1)/ xrm(kk,2)) * xnox(kk)
      cbeta(kk) = (rlm(kk,3) - cpm(kk,1,2))/xrm(kk,3) + &
                              cpm(kk,2,1)/xrm(kk,2) +  &
                 (cpm(kk,1,2)/(xrm(kk,1)*xrm(kk,3))) * (xox(kk)-xnox(kk))
      xrrm(kk,3) = 0.5D0 * (dsqrt(cbeta(kk)**d_two + &
                     d_four*calpha(kk)*cgamma(kk)) - cbeta(kk))/calpha(kk)
      xrrm(kk,2) = xnox(kk) - xrrm(kk,3)
      xrrm(kk,1) = xox(kk) - xrrm(kk,2)
      !
      ! ATTEMPTED NOXBUG FIX: IF O3<0.
      !  ASSUME THAT O3+NO2->NO3 IS CAUSE,  it represents 2*Ox loss, 1*NOx los
      !  USE BACK-EULER TO SET O3
      !  FUTURE - ADD BETTER CRITERIA, TAKE WHICHEVER-IS-LOWER ESTIMATE?
      if ( xrrm(kk,1) < 0.001D0*xrm(kk,1) ) then
        xrrm(kk,1) = xrm(kk,1) * (rpm(kk,1)/rlm(kk,1))
      end if
      !
      ! NOXTEST:  IF KMAX=1 ONLY; SAVE TEST RATIO FOR NO2 vs PRIOR NO2
      c_notest = dabs(d_one-xrrm(1,2)/xrm(1,2))
      pansum = d_zero
      hnosum = d_zero
      do ic = 1 , c_nchem2
        if ( c_icat(ic) == 5 ) pansum = pansum + xc(c_kkw,ic)
        if ( c_icat(ic) == 6 ) hnosum = hnosum + xc(c_kkw,ic)
      end do

      oxox = xrm(c_kkw,ic1) + xrm(c_kkw,ic2)
      sourcox = rpm(c_kkw,1) + rpm(c_kkw,2)
      sinkox = rlm(c_kkw,1) + rlm(c_kkw,2) - cpm(c_kkw,2,1) - cpm(c_kkw,1,2)
      !
      ! QUADRATIC ANALYSIS:  xk2 + (xk2n+xrln) = xjnx + xrpn
      !                     (~NO**2)  (~NO)        (~1)
      !                   calpha*NO**2 + cbeta*NO = cgamma  with new NO.
      !
      if ( xrrm(c_kkw,3) == 0 ) xrrm(c_kkw,3) = d_one
      xk2 = (xrrm(c_kkw,3)**d_two) * calpha(c_kkw)
      xk2n = xrrm(c_kkw,3) * (cpm(c_kkw,1,2)/(xrm(c_kkw,1)*xrm(c_kkw,3))) * &
             (xox(c_kkw)-xnox(c_kkw))
      xrln = xrrm(c_kkw,3) * (rlm(c_kkw,3) - cpm(c_kkw,1,2))/xrm(c_kkw,3)
      xjn = xrrm(c_kkw,3) * cpm(c_kkw,2,1)/xrm(c_kkw,2)
      xjnx = (cpm(c_kkw,2,1)/xrm(c_kkw,2)) * xnox(c_kkw)
!     xrpn =  rpm(c_kkw,3)+ rpm(c_kkw,2)
      xrpn =  rpm(c_kkw,3)
      !
      ! ZERO PROTECT: THIS ALGORITHM SHOULD NEVER GENERATE NEGATIVE
      ! UNLESS THERE IS A MATH ERROR OR TINY NUMERICAL NEGATIVE.
      !
      if ( xrrm(kk,1) <= 0 ) then
        xrrm(kk,1) = 0.001D0
      end if
      if ( xrrm(kk,2) <= 0 ) then
        xrrm(kk,2) = 0.001D0
      end if
      if ( xrrm(kk,3) <= 0 ) then
        xrrm(kk,3) = 0.001D0
      end if
      !
      ! NOTE:  SAVE PRIOR AS XNO, XNO2, XO3?
      !
      ! DIFFICULT CONVERGENCE AID:
      ! AVERAGE O3, NO2, NO WITH PRIOR TIME STEP HERE IF DESIRED.
      !  - note 2009 ALTERNATIVE below.
      if ( c_iter > 5 ) then
        xrrm(kk,3) = 0.5D0*(xrm(kk,3)+xrrm(kk,3))
        xrrm(kk,1) = 0.5D0*(xrm(kk,1)+xrrm(kk,1))
        xrrm(kk,2) = 0.5D0*(xrm(kk,2)+xrrm(kk,2))
      end if
      !
      ! 2009 NOX CONVERGENCE OPTION/ALTERNATIVE:  SETGEOM
      !  adjust based on history of NOx.  (or Ox).
      !
      !  Note, alt:  geom adjust xnox, xox above? No.
      !        cannot do separate geom. avg for NOx; O3-NO-NO2 must be all tog
      if ( c_iter == -1 ) then
        call setgeom(ic2)
      end if
      xrrm(kk,3) = (xrrm(kk,3)**geomavg(kk,ic2)) * &
                   (xrm(kk,3)**(d_one-geomavg(kk,ic2)))
      xrrm(kk,1) = (xrrm(kk,1)**geomavg(kk,ic2)) * &
                   (xrm(kk,1)**(d_one-geomavg(kk,ic2)))
      xrrm(kk,2) = (xrrm(kk,2)**geomavg(kk,ic2)) * &
                   (xrm(kk,2)**(d_one-geomavg(kk,ic2)))
      history(kk,ic2,c_iter) = xrrm(kk,2)+xrrm(kk,3)
!
    end subroutine noxsolve
!
! ---------------------------------------------
!
! Radical balance solution for odd-hydrogen radicals:
!   OH and HO2, also with adjustment for H2O2 and CO3 (aq)
!
! This uses the algorithm from Sillman, 1991 (J. Geophys. Res.)
!
! Critical parameters are:
!   oddhsum:  prior sum of radical net production-loss
!   oddhdel:  prior sum weighted by d(lnHx)/d(lnOH)
!             = sensitivity of radical sum to OH.
!
! These form a linear sum (A+B*OH) for net radical production
!   which is used to solve for OH
!
! (Note:  change in Hx concentration is also added to oddhsum
!   for a non-steady-state solution).
!
! Solution is also adjusted by geometric mean with prior iteration
!   with geometric mean parameter set by iteration history (setgeom)
!
! The following solutions are generated:
!
!   OH is solved from the radical balance equation.
!
!   HO2 is solved from OH/HO2 in the production/loss equation for OH.
!
!   Both are partitioned between gas and aqueous species.
!
!   CO3 (aq) is adjusted based on the radical balance.
!
!   H2O2 is solved here by a normal call to CHEMSOLVE
!     (needed only to insure H2O2 is solved in the proper order)
!
! 2005 CHANGE:  USE RRP, RRL, where RP, RL IS SET TO ZERO.
!
! Inputs:    Species numbers (ic) for OH, HO2, H2O2
!            Also uses:
!             oddhsum, oddhdel, oddhsrc, oddhloh, oddhlho2
!            (oddhsump, oddhfacp calculated during prior iteration -not
!
! Outputs:  xc:    concentrations for OH and HO2
!
!
! Called by:    quadchem
!
! Calls to:     chemsolve (to solve H2O2)
!               setgeom
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
!
    subroutine ohsolve(ic1,ic2,ic3)
!
      implicit none
      integer , intent(in) :: ic1 , ic2 , ic3
!
! LOCAL VARIABLES:
!
! LOCAL VARIABLES USED IN SOLUTION
! foh(kk)           OH/HO2 ratio (declared in 'chemlocal.EXT')
! foh1(kk)          OH/HO2 ratio from prior iteration
! foh1a(kk)         Production of HO2 from OH->HO2 conv, molec/cm3
! foh1b(kk)         Production of HO2 from other sources, molec/cm3
!
! ncsol(kk)              Species list passed to solver subroutine
!
! oddhfac(kk)       Calculated factor for updating OH:
!                        (OH = OHp*oddhfac)
!
      ! OH/HO2 from prior iteration
      real(dp) :: foh1(c_kvec)
      ! Prod HO2 from OH, mol/cm3
      real(dp) :: foh1a(c_kvec)
      ! Prod HO2 from other, mol/cm3
      real(dp) :: foh1b(c_kvec)
      ! Factor for updating OH
      real(dp) :: oddhfac(c_kvec)

      ! Factor for updating OH - w/RO
      real(dp) :: oddhfac1(c_kvec)
      ! Factor for updating OH - w/o
      real(dp) :: oddhfac2(c_kvec)
      ! Weighing Factor for oddhfac 1
      real(dp) :: oddhro2f(c_kvec)

      ! Species list passed to chemsolve
      integer ncsol(c_cdim)
      ! Function to return chem index ic
      integer namechem

      ! Sum for written output
      real(dp) :: pansum
      ! Sum for written output
      real(dp) :: hnosum
      ! Sum for written output
      real(dp) :: xrooh
      ! Sum for written output
      real(dp) :: xhno4
      ! Sum for written output
      real(dp) :: xhno3
      ! Sum for written output
      real(dp) :: xco3
!
      kk = 1
      !
      ! PRELIMINARY OPTION:  SOLVE H2O2 BY A NORMAL CALL TO CHEMSOLVE.
      !  NOTE:  H2O2 HAS A SPECIAL TREATMENT WITHIN CHEMSOLVE
      !  THIS INSURES THAT NET PRODUCTION OF H2O2 IS NO GREATER THAN HX SOURCE
      !
      ncsol(1) = ic3
      do ic = 2 , nsdim
        ncsol(ic) = 0
      end do
      call chemsolve(ncsol)
      !
      ! PRELIMINARY:  CALCULATE REMAINING ODD-H REACTION RATES AND HX SUMS.
      ! INVOKE CHEMSOLVE FOR (OH,-3,HO2) AND (OH,-1,HO2).
      ! THE RESULTING XRP, RPRO, RLOSS ARE USED BELOW.
      ! (NOTE SWITCHED ORDER FOR MAY 1995.  CALC RATES, RP, RL BEFORE RLOSS.)
      ! 2005:  PRELIMINARY CALL TO CHEMSOLV REPLACED - ADD IN HERE
      !  Add for do ic=1,3,2  => oh, ho2
      !  Generate rrp, rrl (include products)
      !  Generate rpro, rloss, xrp
      !  Zero rrp, rrl.
      ! -----------
      ! PRELIMINARY CALCULATION: REMAINING REACTIONS AND XR, RP, RL SUM
      ! -----------
      !  SUM INTO RLOSS, RPRO, XRP (w/ ST ST ADJUSTMENT)
      !  AND ZERO RP, RL
      !  ADOPTED FROM CHEMSOLVE
      !
      do is = 1 , 3 , 2
        ics = ic1
        if ( is == 3 ) ics = ic2
        !
        ! REACTIONS
        !
        if ( c_nnrchem(ics) > 0 ) then
          do i = 1 , c_nnrchem(ics)
            nr = c_nrchem(ics,i)
            call brpro(nr)
          end do
        end if
        if ( c_nnrchp(ics) > 0 ) then
          do i = 1 , c_nnrchp(ics)
            nr = c_nrchmp(ics,i)
            call brpro(nr)
          end do
        end if
        !
        ! ZERO RUNNING SUMS
        !
        rpro(kk,is) =  d_zero
        rloss(kk,is) = d_zero
        xrp(kk,is) = d_zero
        !
        ! ESTABLISH RPRO, RLOSS AND XRP FOR THE 'BASE' SPECIES
        !
        rpro(kk,is) = rpro(kk,is) + rrp(kk,ics) + c_xcin(kk,ics)
        rloss(kk,is) = rloss(kk,is) + rrl(kk,ics)
        xrp(kk,is) = xrp(kk,is) +  xc(kk,ics)
        !
        ! ADD ALL AQUEOUS-EQUILIBRIUM SPECIES INTO SUMMED RPRO, RLOSS AND XRP.
        ! CONVERTING AQUEOUS INTO GAS UNITS (AVOGADRL)
        !
        if ( c_nequil(ics) > 0 ) then
          do neq = 1 , c_nequil(ics)
            ic = c_ncequil(ics,neq)
            xrp(kk,is) = xrp(kk,is) + xc(kk,ic)*c_h2oliq(kk)*avogadrl
            rpro(kk,is) = rpro(kk,is) + rrp(kk,ic)
            rloss(kk,is) = rloss(kk,is) + rrl(kk,ic)
          end do
        end if
        !
        ! NONSTEADY STATE ADJ:  IF NONSTEADY STATE OR ZERO RLOSS, ADD XRP to RPR
        !
        if ( .not. c_lsts(ics) .or. rloss(kk,is) == 0 ) then
          rloss(kk,is) = rloss(kk,is) + xrp(kk,is)
        end if
        if ( rloss(kk,is) <= d_zero ) rloss(kk,is) = 1.0D-08
        !
        ! RECORD RP, RL AND ZERO RRP, RRL
        !  NOTE:  RRP, RRL represent running sum through cascade.
        !         RP, RL represent final sum through cascade.
        !  SUBSEQUENT OH CALCULATION uses final sum (RP, RL), not RRP, RRL.
        !
        do neq = 1 , (c_nequil(ics)+1)
          ic = ics
          if ( neq > 1 ) ic = c_ncequil(ics,(neq-1))
          c_rp(kk, ic) = rrp(kk,ic)
          c_rl(kk,ic) = rrl(kk,ic)
          rrp(kk,ic) = d_zero
          rrl(kk,ic) = d_zero
        end do
      end do
      !
      ! 2009 CORRECTION?
      ! RECORD RP, RL AND ZERO RRP, RRL for H2O2?
      ! GOT HERE.
      !
      ! c_rp(h2o2) should have been set in call to chemsolve, and should be
      !  consistent. It seems bad. Or is it the very low NO?
      !  NO+HO2 8e-12  6e4 7e2  (->4e2) 1.8e3  -> 6e-1
      ! HO2+HO2 2e-12 7e2 7e2 1.8e3  -> 2e-3  correct in REACTION RATE
      !
      !  but calpha, cbeta very different!  unit change?
      ! round 2 rl -ho2 is much much lower, rp H2O2 unchanged.
      !
      ! Q:
      !  Is there a unit change, or a problem with rrp, rrl?
      !  Try senhcat=1, does this help?
      ! Try oddhfac1, 2 and use rlHO2/*+rpH2O2 to scale
      !
      ! -----------
      ! END LOOP - PRELIMINARY CALCULATION OF RP, RL
      ! -----------
      ! ------------
      !  PRELIMINARY:  ADD PRIOR OH, HO2, RO2 INTO ODDHSUM, ODDHDEL
      ! ------------
      ! FIRST:   ADD PRIOR OH, HO2, RO2 INTO ODDHSUM, ODDHDEL
      !   (XR DOES NOT INCLUDE AQ SUM.  JUST RPRO, XRP, ETC.)
      !  CHANGE 1195:  DO NOT INCLUDE PRIOR SUM IF STEADY STATE
      !
      !  CHANGE 1195:  PRIOR ODDHDEL.  If RP, RL is small, then there
      !    PRIOR ODDH shows little sensitivity to OH.
      !  This was solved by changing PRIOR ODDHDEL here and SENCAT below.
      !  JULY 1996:  PRIOR ODDHDEL = XR or (RP+RL), whichever is smaller.
      !  JULY 1996:  INCLUDE AQUEOUS IN PRIOR ODDHSUM AND ODDHDEL.
      !  (NOTE:  COPY THIS SECTION TO 'PRIOR' IN OHWRITE.)
      !
      ! 2005 CHANGE:  This is now done AFTER OH, HO2 reactions processed into
      !  NOTE:  This uses FINAL RP, RL rather than running sum RRP, RRL
      !         since RRP, RRL have been set to zero.
      !
      ! 2009 CHANGE: Do not count lumped species accumulators in PRIOR Hx sum.
      !
      do ic = 1 , c_nchem2
        if ( (c_icat(ic) == 3 .or. c_icat(ic) == 9 .or. &
              c_icat(ic) == 10 .or. c_icat(ic) == 8 ) .and. &
             .not. c_lsts(ic) .and. .not. c_llump(ic) )  then
          if ( ic == c_npequil(ic) ) then
            oddhsum(kk,1) = oddhsum(kk,1) + c_xcin(kk,ic) - xc(kk,ic)
            calpha(kk) = d_zero
            if ( c_rl(kk,ic)+xc(kk,ic) > 0 ) then
              calpha(kk) = (c_rp(kk,ic)+c_rl(kk,ic))*xc(kk,ic) / &
                          (c_rp(kk,ic)+c_rl(kk,ic)+xc(kk,ic))
              !
              !   TEMPORARY 2005 BUG CHANGE!!!! 2005 -
              !     SHOULDN'T oddhdel include PRIOR MINUS ANY SENSITIVITY????
              !            calpha(kk) =  (xc(kk,ic)          )*xc(kk,ic)
              !    *           /(c_rp(kk,ic)+c_rl(kk,ic)+xc(kk,ic))
              !
              oddhdel(kk,1) = oddhdel(kk,1) - senhcat(kk,c_icat(ic))*calpha(kk)
            end if
          else
            if ( c_h2oliq(kk) > 0 ) then
              oddhsum(kk,1) = oddhsum(kk,1) + &
                       (c_xcin(kk,ic)-xc(kk,ic)) * c_h2oliq(kk)*avogadrl
              calpha(kk) = d_zero
              if ( c_rl(kk,ic)+xc(kk,ic) > 0 ) then
                calpha(kk) = &
                  (c_rp(kk,ic)+c_rl(kk,ic))*xc(kk,ic)*c_h2oliq(kk)*avogadrl / &
                  (c_rp(kk,ic)+c_rl(kk,ic)+xc(kk,ic)*c_h2oliq(kk)*avogadrl)
                oddhdel(kk,1) = oddhdel(kk,1) - &
                     senhcat(kk,c_icat(ic)) * calpha(kk)
              end if
            end if
          end if
        end if
      end do
      !
      ! END : ADD PRIOR OH, HO2 TO ODDHSUM INCLUDING RO2
      !
      ! REPEAT ADD PRIOR OH, HO2 WITH ODDHSUM-2, WITHOUT RO2
      !
      do ic = 1 , c_nchem2
        if ( (c_icat(ic) == 8 .or. c_icat(ic) == 9 .or. &
              c_icat(ic) == 10 ) .and. &
              .not. c_lsts(ic) .and. .not. c_llump(ic) )  then
          if ( ic == c_npequil(ic) ) then
            oddhsum(kk,2) = oddhsum(kk,2) + c_xcin(kk,ic) - xc(kk,ic)
            calpha(kk) = d_zero
            if ( c_rl(kk,ic)+xc(kk,ic) > 0 ) then
              calpha(kk) = (c_rp(kk,ic)+c_rl(kk,ic))*xc(kk,ic) / &
                          (c_rp(kk,ic)+c_rl(kk,ic)+xc(kk,ic))
              !
              !   TEMPORARY 2005 BUG CHANGE!!!! 2005 -
              !     SHOULDN'T oddhdel include PRIOR MINUS ANY SENSITIVITY????
              !            calpha(kk) =  (xc(kk,ic)          )*xc(kk,ic)
              !    *           /(c_rp(kk,ic)+c_rl(kk,ic)+xc(kk,ic))
              !
              oddhdel(kk,2) = oddhdel(kk,2) - senhcat(kk,c_icat(ic))*calpha(kk)
            end if
          else
            if ( c_h2oliq(kk) > 0 ) then
              oddhsum(kk,2) = oddhsum(kk,2) + &
                        (c_xcin(kk,ic)-xc(kk,ic))*c_h2oliq(kk)*avogadrl
              calpha(kk) = d_zero
              if ( c_rl(kk,ic)+xc(kk,ic) > 0 ) then
                calpha(kk) = &
                  (c_rp(kk,ic)+c_rl(kk,ic))*xc(kk,ic)*c_h2oliq(kk)*avogadrl / &
                  (c_rp(kk,ic)+c_rl(kk,ic)+xc(kk,ic)*c_h2oliq(kk)*avogadrl)
                oddhdel(kk,2) = oddhdel(kk,2) - &
                      senhcat(kk,c_icat(ic)) * calpha(kk)
              end if
            end if
          end if
        end if
      end do
      !
      ! END REPEAT: ADD PRIOR OH, HO2 TO ODDHSUM
      !
      ! (dHxdNOx*dNOxdOH IS CUT.  IT IS COUNTED THROUGH SENHCAT(NOX).)
      ! Add PRIOR OH, HO2 to oddhloh, oddhlho2
      ! 2009 CORRECTION: steady state control
      ! (ISSUE TO CHECK!)
      if ( .not. c_lsts(ic1) ) then
        do i = 1 , 2
          do neq = 1 , (c_nequil(ic1)+1)
            ic = ic1
            if ( neq > 1 ) ic = c_ncequil(ic1,(neq-1))
            calpha(kk) = d_one
            if ( neq > 1 ) calpha(kk) = c_h2oliq(kk)*avogadrl
            oddhloh(kk,i) = oddhloh(kk,i) + xc(kk,ic)*calpha(kk)
          end do
          do neq = 1 , (c_nequil(ic2)+1)
            ic = ic2
            if ( neq > 1 ) ic = c_ncequil(ic2,(neq-1))
            calpha(kk) = d_one
            if ( neq > 1 ) calpha(kk) = c_h2oliq(kk)*avogadrl
            oddhlho2(kk,i) = oddhlho2(kk,i) + xc(kk,ic)*calpha(kk)
          end do
        end do
      end if
      ! ------------
      !  END PRELIMINARY:  ADD PRIOR OH, HO2, RO2 INTO ODDHSUM, ODDHDEL
      ! ------------
      ! ------------
      ! CONVERGENCE TEST INDICES - SET TO PRIOR FOH, XOH.
      ! ------------
      xfohtest = foh(1)
      c_ohtest = xc(1,ic1)
      ! ------------
      !  FOH = OH/HO2 RATIO
      ! ------------
      !
      ! FOH=OH/HO2 RATIO.  THIS IS SPLIT INTO THREE COMPONENTS:
      ! FOH1A = OH SOURCE FROM HO2.  FOH1B=OTHER OH SOURCES
      ! FOH2 = OH SINKS.
      ! WHEN OH-HO2 EXCHANGE DOMINATES, OH/HO2 = FOH1A/FOH2.
      ! WHEN OTHER SOURCES DOMINATE, IT WORKS TO WEIGH FOH1B BY HO2 SINK/SRC.
      !
      ! FROM OH BACK-EULER EQUATION: k2*HO2+Soh=kl*OH; FOH=k2/k1 + Soh/k2*HO2;
      ! -> FOH=FOHp*(FOH1a/FOH2 + FOH1b/(FOH2*xho2src/sink).
      !
      ! NOTE:  USE XRP(KK,1) = OH PRIOR.  XRP(KK,3)=HO2 PRIOR.
      ! THESE WILL INCLUDE AQUEOUS EQUIVALENT SPECIES SUMS.
      !
      !  CHANGE 1195 - FOH1A, THE OH SOURCE FROM HO2, IS SUMMED
      !  FROM A SPECIAL REACTION INDEX (NRFOH).
      !
      ! CHANGE 1996
      !  OH/HO2 IS ADJUSTED FROM PRIOR BASED ON RP/RL FOR OH, HO2
      !  ALSO: DIFFICULT CONVERGENCE AID TO PROTECT AGAINST CRAZY HO2-H2O2
      !  IN FIRST TIME STEP.
      !    *** CHECK THAT THIS WORKS IN REMOTE TROP, OTHER ENVIRONMENTS.***
      !         RP, RL SHOULD INCLUDE PRIOR OH, HO2 ALSO.
      !
      ! 1996 ALTERNATIVE:   (GOES WITH 1996 WRITE, BELOW)
      !
      foh1(kk) = (xrp(kk,1) / xrp(kk,3))
      foh(kk) = foh1(kk) * (rpro(kk,1)/rloss(kk,1))*(rloss(kk,3)/rpro(kk,3))
      !
      ! 1997 - ADD PROTECT AGAINST WILD SWINGS IN FOH.  (10x, then 50/50 avg.)
      !
      if ( foh(kk) > 100.0D0*foh1(kk) ) foh(kk) = 100.0D0*foh1(kk)
      if ( foh(kk) < 0.01D0*foh1(kk) ) foh(kk) = 0.01D0*foh1(kk)
      !
      !  AUTOMATED DIFFICULT CONVERGENCE OPTION   (FEBRUARY 2005) :
      !   SET GEOMETRIC MEAN vs PRIOR FOR CASES WITH DIFFICULT CONVERGENCE.
      !   Set here for foh.   NOTE:  use OH to store FOH history.
      !
      if ( c_iter >= 4 ) then
        call setgeom(ic1)
      end if
      foh(kk) = (foh(kk)**geomavg(kk,ic1))*(foh1(kk)**(d_one-geomavg(kk,ic1)))
      history(kk,ic1,c_iter) = foh(kk)
      !
      ! FOHTEST = CONVERGENCE TEST FOR FOH. (ABOVE,PRIOR FOHTEST=PRIOR FOH.)
      xfohtest = dabs(d_one-foh(1)/(xfohtest+1.0D-08))
      ! -------------------
      ! ODD HYDROGEN BALANCE.
      ! -------------------
      !     (A + B*OHp = ODDHSUM; A+B*OH=0; B*OHPp=ODDHDEL<0 )
      ! WITH ZERO PROTECT.  WARNING!  THIS HIDES A MULTITUDE OF (NIGHT) SINS.
      !
      ! 1996 ADDITION:  OH, HO2 IS ADJUSTED
      ! SO THAT % CHANGE IN OH*HO2 IS CONSTANT.  (OH=OH*sqrt(foh/fohp))
      !
      !  2005 CHANGE:  OH, HO2 adjusted based on Hx losses attributed to OH vs
      !  So that species responsible for losses changes closest to change in H
      !  Equations:
      !   Lho2*Fho2 + Loh*Foh = (Lho2+Loh )*Fh
      !   OH/HO2 = (OHp/HO2p) * Fdroh      (Fdroh = FOH*HO2p/OHp)
      !   => Foh*(Lho2/Fdroh + Loh) = (Lho2+Loh)*Fh
      !    => Foh = (Lho2+Loh)*Fh/(Lho2/Fdroh+Loh)
      !
      ! (2009 OPTION to set dHx/dOH based on iteration history - deleted.
      !   algorithm description is in chemsolnotes)
      !
      ! 2009 CHANGE:  ODDHFAC = weighted sum of two alternative oddh factors,
      !   one with RO2 included and one WITHOUT RO2.
      !   WEIGHTING: RATEK(no+ho2)*NO versus TIME  (r*NO/(r*NO + 1/time) sec-1
      !    (= include RO2 when rapid conversion RO2 to HO2)
      !
      ! 2009 CHANGE: This only works for NON-AQUEOUS.
      !
      !   ODDH WEIGHTING FACTOR
      oddhro2f(kk) = d_one
      oddhro2f(kk) = d_two*(ratek(kk,c_nrho2no)*xc(kk,c_nno))
      oddhro2f(kk) = oddhro2f(kk)/(oddhro2f(kk)+ d_one/c_time)
      !
      oddhfac1(kk) = d_one
      if ( oddhdel(kk,1) /= 0 ) then
        oddhfac1(kk) = d_one - oddhsum(kk,1)/oddhdel(kk,1)
        if ( oddhfac1(kk) <= d_zero ) oddhfac1(kk) = 0.2D0
      end if
      oddhfac2(kk) = d_one
      if ( oddhdel(kk,2) /= 0 ) then
        oddhfac2(kk) = d_one - oddhsum(kk,2)/oddhdel(kk,2)
        if ( oddhfac2(kk) <= d_zero ) oddhfac2(kk) = 0.2D0
      end if
      oddhfac(kk) = oddhro2f(kk)*oddhfac1(kk)+(d_one-oddhro2f(kk))*oddhfac2(kk)
!     oddhfac(kk) = oddhfac1(kk)
!     xc(kk,ic1) = xrp(kk,1)*oddhfac(kk)   ! moved just below
      ! ---------------
      ! 2009 OPTION to adjust dHx/dOH based on iter history was deleted from h
      ! ---------------
      ! --------------------------
      ! APPLY HX FACTOR TO OH, HO2
      ! --------------------------
      xc(kk,ic1) = xrp(kk,1)*oddhfac(kk)
      if ( foh(kk) /= 0 ) then
        !
        ! 2005 OPTION
        !  (2009 NOTE: what about steady state? Included in foh, added for oddhl
        !     OK:  oddhloh represents loss of Hx linked to OH, not source (=OHin
        !     Prior OH counts as loss in nonsteady state case only.)
        !
        calpha(kk) = oddhfac1(kk)*oddhloh(kk,1)+oddhfac2(kk)*oddhloh(kk,2)
        cbeta(kk) = oddhfac1(kk)*oddhlho2(kk,1)+oddhfac2(kk)*oddhlho2(kk,2)
        if ( calpha(kk) > d_zero .and. cbeta(kk) > d_zero ) then
          xc(kk,ic1) = xc(kk,ic1)*(oddhloh(kk,2)+oddhlho2(kk,2)) / &
                                  (oddhloh(kk,2)+oddhlho2(kk,2) / &
                       (foh(kk)*xrp(kk,3)/xrp(kk,1)))
        else
          xc(kk,ic1) = xc(kk,ic1)*dsqrt(foh(kk)*xrp(kk,3)/xrp(kk,1))
        end if
        !
        ! HO2
        !
        xc(kk,ic2) = xc(kk,ic1)/foh(kk)
      end if
      !
      ! XOHTEST = CONVERGENCE TEST FOR OH.
      !  Based on prior GAS+AQ SUM (XRP) vs NEWLY CALCULATED
      !   = ratio of CHANGE (new-prior)/PRIOR OH
      !   => NIGHTTIME OVERLY STRINGENT: night OH is very small;
      !      better to use change vs RPRO, RLOSS. vs exact balance???
      !
      c_ohtest = xrp(1,1)
      if ( dabs(c_ohtest) < dlowval ) c_ohtest = 1.0D-08
      c_ohtest = dabs(d_one-xc(1,ic1)/(c_ohtest))
      !
      ! 2009 CHANGE:  CONVERGENCE TEST = change in OH, HO2 vs RPRO-RLOSS
      !      a1 =  rloss - rpro  (includes xc-xcin if non-steady-state)
      !      a2 = max:  xc, (rpro+rloss)/2. 1e-8
      !     test = abs(a1/a2)
      !     test for OH, HO2, maximum
      !
      ! CONVERGENCE TEST FOR OH
      !
      ! rloss, rpro already include
      calpha(1) = rloss(1,1) - rpro(1,1)
!     if ( .not. c_lsts(ic1) ) calpha(1) = calpha(1) + xc(1,ic1) - c_xcin(1,ic1)
      cbeta(1) = 0.5D0*(rloss(1,1)+rpro(1,1))
      if ( xc(1,ic1) > cbeta(1) ) cbeta(1) = xc(1,ic1)
      if ( cbeta(1) <= d_zero ) cbeta(1) = 1.0D-8
      c_ohtest = dabs(calpha(1)/cbeta(1) )
      !
      ! CONVERGENCE TEST FOR HO2 (note, rpro(1,3) for HO2)
      !
      ! rloss, rpro include xc-xcin
      calpha(1) = rloss(1,3) - rpro(1,3)
!     if ( .not. c_lsts(ic2) ) calpha(1) = calpha(1) + xc(1,ic2) - c_xcin(1,ic2)
      cbeta(1) = 0.5D0*(rloss(1,3)+rpro(1,3))
      if ( xc(1,ic2) > cbeta(1) ) cbeta(1) = xc(1,ic2)
      if ( cbeta(1) <= d_zero ) cbeta(1) = 1.0D-8
      cgamma(1) = dabs(calpha(1)/cbeta(1) )
      if ( cgamma(1) > c_ohtest ) c_ohtest = cgamma(1)
      !
      ! ADDED CONVERGENCE EXIT AT NIGHT:  IF XOH IS LOW AND NOX CONVERGES.
      !  CUT

      ! DIFFICULT CONVERGENCE OPTION:  GEOMETRIC  AVERAGE WITH PRIOR OH, HO2.

      !  AUTOMATED CONVERGENCE OPTION (FEBRUARY 2005):
      !   SET GEOMETRIC MEAN BASED ON HISTORY (vs HARD-WIRED OPTION above)
      !  HARD-WIRED OPTION:  cut call to setgeom
      !
      !   Set here for Hx.   Use HO2 to store Hx history.
      !   For Hx, adjust based on Hx parameter 1.-oddhsum/oddhdel
      !           and history is running product of adjustments
      !
      !   Note OPTION WITHIN AUTOMATED CONV (SETGEO):  delta vs ratio.
      !
      if ( c_iter >= 4 ) then
        call setgeom(ic2)
      end if
      xc(kk,ic1) = (xc(kk,ic1)**geomavg(kk,ic2)) * &
                   (xrp(kk,1)**(d_one-geomavg(kk,ic2)))
      xc(kk,ic2) = (xc(kk,ic2)**geomavg(kk,ic2)) * &
                   (xrp(kk,3)**(d_one-geomavg(kk,ic2)))
      history(kk,ic2,c_iter) = oddhfac(kk)
!     oddhfacp(kk) = oddhfac(kk)**geomavg(kk,ic2)
      if ( c_iter > 1 ) then
        !
        ! 2009 CORRECTION: ics is error
        ! history(kk,ic2,c_iter) =  &
        !            history(kk,ics,c_iter)*history(kk,ics,(c_iter-1))
        history(kk,ic2,c_iter) = history(kk,ic2,c_iter) * &
                                 history(kk,ic2,(c_iter-1))
      end if
      !
      ! POSSIBLE ADJUSTMENT HERE:
      !   if history positive, and sqrt high vs oddhfac - make oddhfac higher?
      ! FEBRUARY 2005 CO3 HX ADJUSTMENT
      !  ADJUST CO3 BY SAME ADJUSTMENT AS ODD-H.
      !  WITH TEST TO MAKE SURE P(CO3) IS SIGNIFICANT
      !
      ics = namechem('     CO3')
      if ( ics > 0 ) then
        calpha(kk) = d_zero
        do neq = 1 , (c_nequil(ics)+1)
          icc = ics
          if ( neq > 1 )icc = c_ncequil(ics,(neq-1))
          calpha(kk) = calpha(kk) + c_rp(kk,icc)
        end do
        do neq = 1 , (c_nequil(ics)+1)
          icc = ics
          if ( neq > 1 )icc = c_ncequil(ics,(neq-1))
          if ( calpha(kk) >= 0.1D0*oddhsrc(kk,1) ) then
!           xc(kk,icc) = xc(kk,icc) * (oddhfac(kk)**0.5)
            xc(kk,icc) = xc(kk,icc) * (oddhfac(kk)**geomavg(kk,ic2) )
          end if
        end do
      end if
      !
      ! PARTITION OH AND HO2 BETWEEN GAS AND AQUEOUS-EQUILIBRIUM SPECIES
      ! AS IN CHEMSOLVE.  WATCH INDICES! FOR OH, XR(IC2) AND XRP(3).
      !
      do is = 1 , 3
        icc = ic1
        if ( is == 2 ) exit
        if ( is == 3 ) icc = ic2
        if ( icc > 0 ) then
          if ( c_nequil(icc) > 0 ) then
            !
            ! AQUEOUS CONCENTRATIONS ARE UPDATED BY THE RATIO XR/XRP.
            !
            do neq = 1 , c_nequil(icc)
              ic = c_ncequil(icc,neq)
              if ( ic > 0 ) then
                xc(kk,ic) = xc(kk,ic) * xc(kk,icc)/xrp(kk,is)
              end if
            end do
            !
            ! GAS-MASTER XR IS REDUCED BY AQUEOUS CONCENTRATIONS
            !  WITH UNIT CONVERSION (AVOGADRL)
            !
            do neq = 1 , c_nequil(icc)
              ic = c_ncequil(icc,neq)
              if ( ic > 0 ) then
                xc(kk,icc) = xc(kk,icc) - xc(kk,ic)*c_h2oliq(kk)*avogadrl
              end if
            end do
          end if
        end if
      end do
      !
      ! ESTABLISH SENHCAT (=d lnXR/d lnOH) FOR ODD HYDROGEN SPECIES.
      !  (SENHCAT IS USED TO RESOLVE HX = A+B*OH.  IN THE SUM,
      !   ODDHSUM = + dHx; ODDHDEL = + dHx*SENHCAT =>B*OH.)
      ! FOR HO2:  dHO2/dOH from HO2 EQ:
      !  kOH+otherSho2=k1HO2+2k2HO2**2 = 2RPh2o2 + other RLho2.
      !   (NOTE:  THIS USES IC INDEX FOR HO2, H2O2.)
      ! (NOTE:  TECHNICALLY, d/dOH = d/dOH + d/dHO2 dHO2/dOH
      !  + d/dNOx dNOx/dOH.  NORMALLY, d/dNOx IS USED ONLY FOR NOX SPECIES.
      !  OLD PROBLEM:  PAN OSCILLATION.  Maybe solve this by adding
      !   dRO2/dOH = dHO2/dOH + dRO2/dNOx dNOx/dOH = dHO2/dOH (1-dNOx/dOH)
      !   where dNOx/dOH<0.  )
      !
      !  Possible error (bmaingtest).  Originally both HO2 and RO2 senhcat
      !  were multiplied by RP/RP+XXO.  This term was first dropped for RO2.
      !  Then after 'bmaingtest' it worked best to drop the term for HO2 and R
      !
      ! 1996 - AQUEOUS CHANGE.  RL, RP IS SCREWED UP BY EXCORR WITH AQUEOUS.
      ! TO FIX - REPLACE RL, RP WITH RLOSS, (RPRO-XRP) = GAS+AQUEOUS SUM.
      !  ALGORITHM:  SENHCAT = RLOSS(HO2)/(RLOSS(HO2) + 2 RP(H2O2))
      !
      ! 2009 CORRECTION: rp(H2O2) includes other sources (ALK+O3), leads to er
      ! Trying an alternative:  oddhlho2/2. = rp(H2O2) from Hx...
      !     but it also includes rp (ROOH).
      !
      ! note - is not the maximum dln HO2/dln OH = 0.5, not 0? from HO2+HO2,
      !   but less if other sources
      !  -> it should be rp(HO2 from OH)/(*+rp(HO2 from other sources))
      !    = rp(HO2)
      !
      ! OPTION TO TRY: rp(HO2)*(1+RO2/HO2)/(+oddhsrc)?
      !    rpHO2(1+oddhsrc2/oddhsrc)/(*+oddhsrc)  ?
      ! ALT: loh/(+(oddhsrc-2O3+hv)
      !
      ! SUM GAS AND AQUEOUS:  RLOSS(HO2), RPRO(H2O2)
      !
      calpha(kk) = c_rl(kk,ic2)
      cbeta(kk) = c_rp(kk,ic3)
      if ( c_nequil(ic2) > 0 ) then
        do neq = 1 , c_nequil(ic2)
          ic = c_ncequil(ic2,neq)
          calpha(kk) = calpha(kk) + c_rl(kk,ic)
        end do
      end if
      if ( c_nequil(ic3) > 0 ) then
        do neq = 1 , c_nequil(ic3)
          ic = c_ncequil(ic3,neq)
          cbeta(kk) = cbeta(kk) + c_rp(kk,ic)
        end do
      end if
      !
      ! 2009 correction: cbeta (ph2o2) < calpha (lho2)
      !    protects against pH2O2 from source other than HO2 (TERP+O3->0.02 H2
      !
      if ( cbeta(kk) > calpha(kk) ) cbeta(kk) = calpha(kk)
      !
      ! MAIN SENHCAT CALCULATION.
      !
      ! 1996 VERSION
      ! 2005:  GEOM MEAN for DIFFICULT CONVERGENCE (HO2-MCO3-PAN oscillation)
      !
      senhcat(kk,10) = calpha(kk) / (calpha(kk) + d_two*cbeta(kk))
      !
      ! OPTION w/o  GEOM MEAN
      !
      ! senhcat(kk,3) = senhcat(kk,10)
      ! OPTION with GEOM MEAN
      !
      senhcat(kk,3) = (senhcat(kk,10)**0.5D0)*(senhcat(kk,3)**0.5D0)
      senhcat(kk,10) = senhcat(kk,3)
      !
      ! OLDER HO2-RO2 SENHCAT ALGORITHMS
      ! PRIOR 1994 VERSION
      ! senhcat(kk,10) =  c_rl(kk,ic2) / &
      !                   (c_rl(kk,ic2)+2.*c_rp(kk,ic3))
      ! senhcat(kk,3) = senhcat(kk,10)
      !
      ! OLDER HO2-RO2 SENHCAT ALGORITHMS
      !       senhcat(kk,10) = (c_rp(kk,ic2)-  c_xcin(kk,ic2)) * c_rl(kk,ic2)
      !    *          /(      c_rp(kk,ic2) * (c_rl(kk,ic2)+2.*c_rp(kk,ic3)) )
      !       senhcat(kk,10) = (c_rp(kk,ic2)            ) * c_rl(kk,ic2)
      !    *                /(     ( c_rp(kk,ic2)+  c_xcin(kk,ic2))
      !    *                        * (c_rl(kk,ic2)+2.*c_rp(kk,ic3)) )
      !       senhcat(kk, 3) =
      !    *    c_rl(kk,ic2)    /(   (c_rl(kk,ic2)+2.*c_rp(kk,ic3)) )
      !
      ! dRO2/dNOx OPTION
      !    *     *(1.-senhcat(kk,12))
      !
      ! SENHCAT FOR RO2, RCO3:  FOR NOW, SET EQUAL TO HO2.
      ! THE true VALUE IS 1-fhp*(1-fhr)
      ! WHERE fhp=k(RO2+HO2)/c_rl(RO2); fhr=k(ROOH+OH)/c_rl(ROOH).
      ! AND FOR RCO3:  1 - (1-fpp)(fhp)(1-fhr)
      ! WHERE fpp=k(RCO3+HO2)/c_rl(RCO3); fhp, fhr as above.
      !
      ! IF ROOH IS STEADY STATE, ITS SENHCAT IS EQUAL TO THAT FOR RO2.
      ! THIS REQUIRES THAT ROOH STEADY STATE IS THE SAME AS H2O2.
      !   (PROBABLY DOESN'T WORK ANYWAY.... FOR ROOH STST REMOVE ROOH.)
      !
      if ( c_lsts(ic3) ) then
        senhcat(kk,4) = senhcat(kk,3)
      end if
      !
      ! ----------------------------------------------------------
      ! END OHSOLVE
      ! ----------------------------------------------------------
    end subroutine ohsolve
!
! ----------------------------------------------------------
!
! This calculates the geometric average factor (geomavg(kk,kc))
!    which is used to adjust concentrations as a  geom. avg. with prior.
!      ( XR = (XR**F)*(XRP**(1-F))
!
! The factor is set based on the HISTORY of the past 3 iterations.
! When XR oscillates, GEOMAVG gets lower towards minimum (0.5 or 0.1)
! When XR keeps increasing or decreasing, GEOMAVG moves towards 1.
!
! OPTION:  How strongly does GEOMAVG react to history?
! OPTION:  Set minimum value for GEOMAVG (geomin)
!           Currently set to 0.1 for slow, sure convergence.
!            (previously 0.5)
!
! Inputs:    Species number (ic)
!            history(kk,ic,c_iter)  - past iterative solutions
!            geomavg(kk,ic) - previous geometric avg factor.
!
!
! Outputs:  geomavg(kk,ic) = new geometric average factor for species.
!
!
! Called by:    quadchem
!
! Calls to:     chemsolve (to solve H2O2)
!               setgeom
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
    subroutine setgeom(ic)
!
      implicit none
      integer , intent(in) :: ic
      ! Minimum value for geo. avg. parameter
      real(dp) :: geomin
!
      kk = 1
      !
      ! LINE TO COMMENT OUT AUTOMATIC SETTING
      !
!     if ( c_iter > 0 ) return
 
      !
      ! RETURN IF ITER < 4
      !
      if ( c_iter < 4 ) return
      !
      ! OPTION:  SET MINIMUM VALUE FOR GEOMETRIC AVERAGE PARAMETER.
      !   Originally set at 0.5, then 0.3.
      !   Currently set at 0.1.  Low value insures slow, sure convergence.
      !
      geomin = 0.1D0
      !
      ! ESTABLISH TREND OVER PAST 3 ITERATIONS.
      !  A = delta(last iter) / delta (previous iter)
      !  OPTION:  A= [new/old ratio(last iter)] / [new/old ratio(previous iter
      !     A<0:  oscillates.
      !     A<-1:  oscillation getting worse
      !     A>0:  monotone increase/decrease
      !     A=> 0:  approaching convergence
      !
      calpha(kk) = d_zero
      if ( history(kk,ic,(c_iter-2)) /= 0 .and. &
           history(kk,ic,(c_iter-3)) /= 0) then
        !
        ! DELTA OPTION
        !
        !  cbeta(kk) = history(kk,ic,(c_iter -1 )) - &
        !             history(kk,ic,(c_iter-2))
        !  cgamma(kk) = history(kk,ic,(c_iter-2)) - &
        !               history(kk,ic,(c_iter-3))
        !
        ! RATIO OPTION
        !
        cbeta(kk) = history(kk,ic,(c_iter-1))/history(kk,ic,(c_iter-2))-d_one
        cgamma(kk) = history(kk,ic,(c_iter-2))/history(kk,ic,(c_iter-3))-d_one
        calpha(kk) = cbeta(kk)/cgamma(kk)
      end if
      if ( calpha(kk) >  d_one ) calpha(kk) =  d_one
      if ( calpha(kk) < -d_one ) calpha(kk) = -d_one
      !
      ! DIFFICULT CONVERGENCE  OPTION:  MULTIPLY SETGEOM BY DAMPENING FACTOR.
      !    Multiply by 1:  moves factor towards 1 or 0.5 instantly.
      !    Multiply by 0.5 (standard):  moves factor more slowly.
      !
      calpha(kk) = calpha(kk) * 0.5D0
      !
      ! ADJUST GEOMETRIC AVERAGING FACTOR:
      !   If A<0 (oscillating), move towards minimum (0.5)(0.1)(geomin)
      !   If A>0 (steady trend),move towards maximum (1.)
      !
      if ( calpha(kk) > 0.2D0 ) then
        geomavg(kk,ic) = geomavg(kk,ic) + calpha(kk)*(d_one-geomavg(kk,ic))
      end if
      if ( calpha(kk) < -0.2D0 ) then
        geomavg(kk,ic) = geomavg(kk,ic) + calpha(kk)*(geomavg(kk,ic)-geomin)
      end if
      !
    end subroutine setgeom
!
! ----------------------------------------
!
!  This sets the initial value for six key species
!                        (OH, HO2,O3,NO,NO2,NO3, also H+)
!   based on NO-NO2-O3 equilibrium for NOx-Ox
!   and simple OH/HO2 ratio and jO3 for odd-h radicals.

!  The subroutine uses the following reactions:
!
!       #1:  NO2+hv->NO+O3          #2: NO+O3->NO2
!       #9:  O3+hv->2OH             #12: NO2+OH->HNO3
!       #16: OH+CO->HO2             #17: OH+O3->HO2
!       #18: HO2+NO->OH+NO2         #22: HO2+HO2->H2O2
!       #3:  NO2+O3->NO3
!  Also:
!       #17 OH+O3
!       #46 OH+CH4->HO2.
!
! The reactions are identified from "family" array:
!   family(1,i) = OH, HO2, H2O2
!   family(2,i) = O3, NO2, NO
!   family(3,i) = NO3, N2O5, HNO3
! These are used to identify reactions based on reactants.
!
! Initial  H+ is set to 1e-5.
!
! Inputs:    Initial species concentrations (c_xcin)
!            Chemical reactions
!            family array to ID reactions
!
! Outputs:   Estimated xc for OH, HO2, O3, NO, NO2, NO3, H+
!
!
! Called by:    quadchem
!
! Calls to:     None.
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
     subroutine presolve
!
       implicit none

       integer :: nspecial(40)
       ! (ic1,ic2,ic3,ic9,ic10 ic4)
       ! (O3  NO2 NO  OH  HO2  NO3)
!
       kk = 1
       !
       ! REACTION SUM:  IDENTIFY SPECIAL REACTIONS FOR HX, NOX
       ! AND SUM 'CRATE' (ALPHA) = OH+C REACTIONS.
       ! AND SUM GAMMA = NO3 LOSS REACTIONS.
       !
       do i = 1 , 40
         nspecial(i) = 0
       end do
       calpha(kk) = d_one
       cgamma(kk) = d_one
       do nr = 1 , c_nreac
         !
         ! IDENTIFY SPECIAL REACTIONS BY REACTANT NUMBER
         !       #1:  NO2+hv->NO+O3          #2: NO+O3->NO2
         !       #9:  O3+hv->2OH             #12: NO2+OH->HNO3
         !       #16: OH+CO->HO2             #17: OH+O3->HO2
         !       #18: HO2+NO->OH+NO2         #22: HO2+HO2->H2O2
         !       #3:  NO2+O3->NO3
         if ( c_reactant(nr,1) == c_nno2 .and. &
              c_reactant(nr,2) == -1) nspecial(1) = nr
         if ( (c_reactant(nr,1) == c_no3 .and. &
               c_reactant(nr,2) == c_nno) .or. &
              (c_reactant(nr,1) == c_nno .and. &
               c_reactant(nr,2) == c_no3) ) nspecial(2) = nr
         if ( c_reactant(nr,1) == c_no3 .and. &
              c_reactant(nr,2) == -1 .and. &
              c_product(nr,1) == c_noh ) nspecial(9) = nr
         if ( (c_reactant(nr,1) == c_nno2 .and. &
               c_reactant(nr,2) == c_noh) .or. &
              (c_reactant(nr,1) == c_noh .and. &
               c_reactant(nr,2) == c_nno2) ) nspecial(12) = nr
         if ( (c_reactant(nr,1) == c_nno .and. &
               c_reactant(nr,2) == c_nho2) .or. &
              (c_reactant(nr,1) == c_nho2 .and. &
               c_reactant(nr,2) == c_nno) ) nspecial(18) = nr
         if ( c_reactant(nr,1) == c_nho2 .and. &
              c_reactant(nr,2) == c_nho2 ) nspecial(22) = nr
         if ( (c_reactant(nr,1) == c_nno2 .and. &
               c_reactant(nr,2) == c_no3) .or. &
              (c_reactant(nr,1) == c_no3 .and. &
               c_reactant(nr,2) == c_nno2) ) nspecial(3) = nr
         !
         ! IDENTIFY OH+CO,HC OR O3 REACTIONS BY CATEGORY. SUM AS 'CRATE'(calpha).
         !
         icat1 = 0
         icat2 = 0
         icr = 0
         if ( c_reactant(nr,1) > 0 ) then
           icat1 = c_icat(c_reactant(nr,1))
         end if
         if ( c_reactant(nr,2) > 0 ) then
           icat2 = c_icat(c_reactant(nr,2))
         end if
         if ( icat1 == 9 .and. (icat2 == 11 .or. icat2 <= 3) ) then
           icr = c_reactant(nr,2)
         end if
         if ( icat2 == 9 .and. (icat1 == 11 .or. icat1 <= 3) ) then
           icr = c_reactant(nr,1)
         end if
         if ( icr > 0) then
           calpha(kk) = calpha(kk)+ratek(kk,nr)*xc(kk,icr)
         end if
         !
         ! IDENTIFY NO3 LOSS REACTIONS (cgamma)
         !
         icr = -1
         if ( c_reactant(nr,1) == c_nno3 ) icr = c_reactant(nr,2)
         if ( c_reactant(nr,2) == c_nno3 ) icr = c_reactant(nr,1)
         if ( icr > 0 ) then
           cgamma(kk) = cgamma(kk)+ratek(kk,nr)*xc(kk,icr)
         end if
         if ( icr == 0 ) then
           cgamma(kk) = cgamma(kk)+ratek(kk,nr)
         end if
       end do
       !
       ! NOX: A PRECISE INITIALIZATION WOULD SET EQUILIBRIUM BTWEEN NO+O3=NO2.
       ! EQUATION:  (NO)**2 + NO (Ox-NOx+j1/k2) - NOx j1/k2 = 0
       ! WHERE j1:  NO2+hv; k2: NO+O3.
       !
       if ( nspecial(1) > 0 .and. nspecial(2) > 0 ) then
         cbeta(kk) = ratek(kk,nspecial(1))/ratek(kk,nspecial(2))
         xc(kk,c_nno) = 0.5D0 * &
              (dsqrt((c_xcin(kk,c_no3)-c_xcin(kk,c_nno)+cbeta(kk))**d_two + &
                    d_four*cbeta(kk)*(c_xcin(kk,c_nno2)+c_xcin(kk,c_nno)))- &
                    (c_xcin(kk,c_no3)-c_xcin(kk,c_nno)+cbeta(kk)))
         xc(kk,c_nno2) = c_xcin(kk,c_nno2)+c_xcin(kk,c_nno)-xc(kk,c_nno)
         xc(kk,c_no3) = c_xcin(kk,c_no3)+  c_xcin(kk,c_nno2)-xc(kk,c_nno2)
         !
         ! ZERO-PROTECT SHOULD APPLY ONLY IF ERROR IN FORMULA.
         !
         if ( xc(kk,c_nno2) <= d_zero ) xc(kk,c_nno2) = 0.1D0
         if ( xc(kk,c_nno) <= d_zero ) xc(kk,c_nno) = 0.1D0
         if ( xc(kk,c_no3) <= d_zero ) xc(kk,c_no3) = 0.1D0
         !
         ! INITIAL NOX REDUCED BASED ON 5-HOUR LIFETIME
         !  COMMENT OUT IF PRESET NOx OPTION)
         !
         if ( .not. c_lsts(c_nno2) .and. .not. c_lsts(c_nno)) then
           xc(kk,c_nno) = xc(kk,c_nno)/(d_one+c_time/18000.0D0)
           xc(kk,c_nno2) = xc(kk,c_nno2)/(d_one+c_time/18000.0D0)
         end if
       end if
       !
       ! NO3:  PRELIMINARY VALUE SET FROM BACK-EULER FORMULA:
       ! NO2+O3 SOURCE, ALL NO3 SINKS. (Added to deal with large-isop. crash.)
       xc(kk,c_nno3) = c_xcin(kk,c_nno3)
       if ( nspecial(3) > 0 ) then
         xc(kk,c_nno3) = xc(kk,c_nno3) + &
                         ratek(kk,nspecial(3))*xc(kk,c_no3)*xc(kk,c_nno2)
       end if
       xc(kk,c_nno3) = xc(kk,c_nno3) / (d_one + cgamma(kk) )
       !
       ! xno, xno2, xo3, xoh, xho2 - ALL CUT.
       ! oxno, oxno2, oxo3, rno, rno1, rno2, xnox - ALL CUT
       !
       ! OH, HO2 INIITIALIZE:REQUIRES PRE-SET REACTION NUMBERS (see above)
       !  PLUS AUTOMATED CRATE = #16, OH+CO, #17 OH+O3,#46 OH+CH4->HO2.
       !  DECEMBER 1994:  CHANGE TO FULL QUADRATIC WITH PRIOR OH, HO2
       !  AND PROTECT AGAINST NIGHTTIME STEADY STATE ZERO.  HO2 prior>1E6.
       !
       foh(kk ) = 0.01D0
       if ( ratek(kk ,nspecial(1)) >= 1.0D-03 ) then
         foh(kk ) = ratek(kk,nspecial(18)) * xc(kk,c_nno)/calpha(kk)
       end if
       !
       ! FULL QUADRATIC SOLVE WITH PRIOR OH, HO2
       !
       cgamma(kk ) = d_two*c_time*ratek(kk ,nspecial(9))*xc(kk ,c_no3) + &
                     c_xcin(kk,c_noh)+c_xcin(kk,c_nho2)
       if ( cgamma(kk ) < 1.01D+04 ) then
         c_xcin(kk,c_nho2) = c_xcin(kk,c_nho2) + 1.0D+04
         c_xcin(kk,c_noh) = c_xcin(kk,c_noh) + 1.0D+02
         cgamma(kk) = cgamma(kk) + 1.01D+04
       end if
       cbeta(kk ) = d_one + d_one/foh(kk ) +  &
                   c_time*ratek(kk,nspecial(12))*xc(kk,c_nno2)
       calpha(kk ) = c_time* ratek(kk,nspecial(22))/(foh(kk)**d_two)
       if ( calpha(kk )*cgamma(kk ) < cbeta(kk)**d_two ) then
         xc(kk,c_noh) = cgamma(kk)/cbeta(kk)
       else
         xc(kk,c_noh) = (dsqrt(cbeta(kk)**d_two + &
                d_four*calpha(kk)*cgamma(kk))-cbeta(kk))/(d_two*calpha(kk))
       end if
       xc(kk,c_nho2) = xc(kk,c_noh)/foh(kk)
       !
       ! INITIAL CONCENTRATION FOR H+
       ! --------------------------
       if ( c_aqueous(1,2) > 0 ) then
         xc(kk,c_aqueous(1,2)) = 1.0D-05
         xc(kk,c_aqueous(1,3)) = 1.0D-09
       end if
       ! --------------------------
       ! END INITIALIZATION ROUTINE
       ! --------------------------
     end subroutine presolve
!
! ---------------------------------------------------------------
!
! This sets initial concentrations for 'lumped' species.
! And sets other initial values  for the start of the chemistry solver
!  ('prechem').
!
! For lumped species:  It converts lumped sum and partition fractions
!    into concentrations for individual species.
!    and sets the input concentration (c_xcin) for lumped species.
!
! For all species:  it enters the input concentration (c_xcin) as xc.
!
! It  sets aqueous species to zero.
!   (Assumes that the saved values and initial concentrations
!     represent gas+aqueous sum for the gas-master species;
!     output aqueous concentrations are for information only.)
!
! It  sets zero initial values for running sums
!   (production, loss - RP, RL)
!
! It sets initial value for geometric average.
! it sets initial value for xcfinr (FINAL/AVG RATIO)
!
! OPTION:  SET XR=0.1 FOR ALL STEADY-STATE SPECIES?
! HOWEVER, IF A SPECIES IS SET IN STEADY STATE VS EMISSIONS,
! (e.g. ISOPRENE) THEN THIS SHOULD BE OVERRIDDEN.

! OPTION:  SET XR=0 FOR ALL UNPARTITIONED LUMPED SPECIES (RCO3-RPAN)
! SO THAT RP() IN MIDLUMP SETS PARTITIONING CORRECTLY
!
! Inputs:    Initial species concentrations (c_xcin)
!
! Outputs:   Working species concentrations (xc)
!            Initial species concentrations for lumped sp. (as c_xcin)
!            Initial zero values  for production, loss (rp, rl)
!            Initial geomavg
!
!
! Called by:    quadchem
!
! Calls to:     None.
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
     subroutine prelump
!
       implicit none
!
       kk = 1
       !
       ! SET H2O CONCENTRATION IN CHEMISTRY ARRAY FROM INPUT VARIABLE
       !
       if ( c_nh2o > 0 ) then
         xc(kk,c_nh2o) = c_h2ogas(kk)
       end if
       !
       !  SET XC=XCIN, ZERO PROTECT,  AND MAYBE XC=0.1 FOR STEADY STATE.
       !  FOR ALL TRANSPORTED SPECIES.
       !
       loopnchem1:&
       do ic = 1 , c_nchem1
         !
         ! possible error with this and other ZERO PROTECT?
         ! 2004 fix - insures no descent to zero->NaN.
         !
         if ( c_xcin(kk,ic) <= 0.000001D0 ) c_xcin(kk,ic) = 0.000001D0
         xc(kk,ic) = c_xcin(kk,ic)
         !
         ! SET XR=0.1 AND XXO=0. FOR STEADY-STATE.
         !  XXO MUST BE ZERO, ELSE TROUBLE WHEN RP, RL=0.
         !  THIS IS THE LINE TO COMMENT OUT, MAYBE.
         if ( c_lsts(ic) ) then
           xc( kk ,ic) = 0.1D0
           c_xcin(kk,ic) = d_zero
           c_xcemit(kk,ic) = d_zero
         end if
         !
         ! ----------------------
         ! LUMPED SPECIES ARE IDENTIFIED FROM ARRAY LUMP(I,J)
         ! READ FROM REACTION.DAT
         ! ----------------------
         !
         do i = 1 , c_cdim
           ics = c_lump(i,1)
           if ( ics == 0 ) exit loopnchem1
           ic1 = c_lump(i,2)
           ic2 = c_lump(i,3)
           !
           ! SUM PARTITION FRACTIONS AND CORRECT IF ZERO.
           ! IF PRIOR PARTITION FRACTIONS ARE NOT SAVED: INITIAL XR
           ! SHOULD BE ZERO.
           ! ORIGINAL OPTION:  SET PARTITION = 100% FIRST SPECIES (IC1)
           ! AND CORRECT AT MIDLUMP BASED ON RP.  PROBLEM:  RCO3-RPAN.
           ! ALTERNATIVE OPTION: LEAVE PARTITIONED XR=0 AND CORRECT AT MIDLUMP.
           ! TO IMPLEMENT, SEE 'LUMPED SPECIES' BELOW.
           !
           calpha(kk) = d_zero
           do iic = ic1 , ic2
             calpha(kk) = calpha(kk) + xc(kk,iic)
           end do
           if ( calpha(kk) == 0 ) then
             xc(kk,ic1) = d_one
             calpha(kk) = d_one
           end if
           !
           ! CONVERT PARTITION FRACTIONS INTO CONCENTRATIONS. ALSO ENTER XXO.
           ! ALSO ZERO-PROTECT AND ENTER XXO.
           !
           do iic = ic1 , ic2
             c_lsts(iic) = c_lsts(ics)
             !
             ! ORIGINAL VERSION - SET LUMPED SPECIES HERE.
             !            xc(kk,ic) = xc(kk,ic)*xc(kk,ics)/calpha(kk)
             !            if(xc(kk,ic) <= 0.1) xc(kk,ic) = 0.1
             ! ALTERNATIVE OPTION:  SET LUMPED SPECIES = 0. FOR FIRST ITERATION;
             !          THEN RESET BASED ON RP IN MIDLUMP.
             xc(kk,iic) = 0.1D0
             !
             ! PUT LUMPED VALUES IN xcin
             !
             c_xcin(kk,iic) = xc(kk,iic)
           end do
           !
           ! SET EMISSIONS TO BE LESS THAN INPUT CONCENTRATION
           !   (Input concentration includes emissions;
           !      emissions are only used for timing in expo decay solution)
           do iic = 1 , c_nchem2
             if ( c_xcemit(kk,iic) > 0.99D0*c_xcin(kk,iic) ) then
               c_xcemit(kk,iic) = 0.99D0*c_xcin(kk,iic)
             end if
           end do
         end do
       end do loopnchem1
       !
       ! SET AQUEOUS CONCENTRATIONS EQUAL TO ZERO
       !
       do nrh = 1 , c_nreach
         icc = c_henry(nrh,1)
         if ( icc > 0 ) then
           if ( c_nequil(icc) > 0 ) then
             do neq = 1 , c_nequil(icc)
               ic = c_ncequil(icc,neq)
               xc(kk,ic) = d_zero
             end do
           end if
         end if
       end do
       !
       ! ZERO REACTION RATES (RR) TO PREVENT CARRY-OVER FROM PREVIOUS RUN.
       !
       do nr = 1 , c_nreac
         c_rr(kk,nr) = d_zero
       end do
       !
       ! ZERO RP, RL, RRP, RRL TO PREVENT CARRY-OVER FROM PREVIOUS RUN
       !         (does not zero RRP, RRL before each iteration, only at start)
       !
       do ic = 1 , c_nchem2
         c_rp(kk,ic) = d_zero
         c_rl(kk,ic) = d_zero
         rrp(kk,ic) = d_zero
         rrl(kk,ic) = d_zero
         rppair(kk,ic) = d_zero
         rlpair(kk,ic) = d_zero
       end do
       !
       ! TEST AND CORRECT RAINOUT PARAMETER (0<=rainfr<1)
       !
       if ( c_rainfr(kk) < d_zero ) c_rainfr(kk) = d_zero
       if ( c_rainfr(kk) > 0.9999D0 ) c_rainfr(kk) = 0.9999D0
       !
       ! SET DIFFICULT CONVERGENCE PARAMETER - INITIAL VALUES
       !   Normally set to 1 for all except odd hydrogen.
       !   Odd hydrogen:  HO2 (=Hxsum) factor set to 0.7
       !                  OH  (=FOH)   factor set to 1.
       do ic = 1 , c_nchem2
        geomavg(kk,ic) = d_one
       end do
       ic = c_nho2
       geomavg(kk,ic) = 0.7D0
       ic = c_noh
       geomavg(kk,ic) = 1.0D0
       !
       ! Set ratio of FINAL to AVERAGE species concentrations (xcfinr)
       !   This ratio is always ONE for back-Euler solution.
       !   It is set to a different value when the EXPO DECAY solution is used.
       !
       do ic = 1 , c_nchem2
         xcfinr(kk,ic) = d_one
       end do
!
     end subroutine prelump
!
! -------------------------------------
!
! This re-sets the initial concentrations (c_xcin) for lumped species
!  based on calculated chemical production and loss.
!
! Initial input for lumped species consists of the lumped sum only.
! The partitioning of the sum into individual species
!  is assumed to be in proportion to chemical production.
!
! As production rates are iteratively calculated,
!   the assumed initial concentration of lumped species is adjusted.
!
! This is called at the start of each iteration (except #1).
! It adjusts the concentration at the start of the time step (c_xcin)
!   ,not the solution for the end of the time step (xc).
!
! 2006 MODIFICATION - SUM LUMPED GAS-PHASE SPECIES INTO XC
!   FOR USE IN REACTION RATES.
!
! 2008 LUMP OPTION: partition initial concentrations based on rp/(1+rl)
!      rather than rp (to avoid major error if lumped species are
!      chemically different - THEY SHOULD NOT BE DIFFERENT.
!
!
! Inputs:    Initial species concentrations (xc)
!
! Outputs:   Modified initial concentrations for lumped species (c_xcin)
!            Sums for lumped species (xc)
!
!
! Called by:    quadchem
!
! Calls to:     None.
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
     subroutine midlump
!
       implicit none
!
       kk = 1
       !
       ! LUMPED SPECIES ARE IDENTIFIED FROM ARRAY LUMP(I,J)
       ! ICS = ACCUMULATOR; IC1-IC2 ARE INDIVIDUAL SPECIES.
       !
       do i = 1 , c_cdim
        ics = c_lump(i,1)
        if ( ics == 0 ) exit
        ic1 = c_lump(i,2)
        ic2 = c_lump(i,3)
        !
        ! SUM LUMPED SPECIES (IC1-IC2) INTO ACCUMULATOR (ICS)
        !   (for use in rate calculations)
        !
        xc(kk,ics) = d_zero
        do ic = ic1 , ic2
          xc(kk,ics) = xc(kk,ics) + xc(kk,ic)
        end do
        !
        ! SUM PARTITION CHEM. PRODUCTION RATES
        !
        calpha(kk) = d_zero
        do ic = ic1 , ic2
          !
          ! 2006 ORIGINAL
          !
          calpha(kk) = calpha(kk) + c_rp(kk,ic)
          !
          ! 2008 LUMP OPTION
          !
          ! calpha(kk) = calpha(kk) + c_rp(kk,ic)/(1.+c_rl(kk,ic))
          !
          ! END OPTION
          !
        end do
        !
        ! MODIFY XXO FOR PARTITIONED SPECIES BASED ON PARTITIONING OF
        !  CHEM PRODUCTION.
        ! OCTOBER 1998 CORRECTION - XXO only.  XR would require POSTLUMP
        !  SUM OF LUMPED SPECIES AND GAS+AQUEOUS.
        !
        do ic = ic1 , ic2
          if ( calpha(kk) > d_zero ) then
            !
            ! 2006 ORIGINAL
            !
            c_xcin(kk,ic) = c_xcin(kk,ics)*c_rp(kk,ic)/calpha(kk)
            !
            ! 2008 LUMP OPTION
            !
            ! c_xcin(kk,ic) = c_xcin(kk,ics) * &
            !           (c_rp(kk,ic)/(d_one+c_rl(kk,ic)))/calpha(kk))
            !
            ! END OPTION
            !
            if ( c_xcin(kk,ic) <= d_zero ) c_xcin(kk,ic) = 0.1D0
          end if
        end do
        !
        ! END LUMPED-SPECIES LOOP
        !
      end do
!
     end subroutine midlump
!
! -------------------------------------
!
! This calculates the final output species concentration (c_xcout)
!  from the working concentration (xc) (average over time step).
!  It also preserves the average species concentrations (c_xcav)
!  and the FINAL/AVERAGE ratio (xcfinr).
!
! It also  calculates sums for all 'lumped' species
!  at the end of the chemistry solution, and saves them in xc.
!
! This also sums aqueous equilibrium species
!  into a combined gas+aqueous value (in gas-phase units, molec/cm3)
!  which is saved as the gas-master species concentration.
!
! Aqueous concentrations are saved as output (in aqueous M/liter),
!   but it is assumed that these are not transported from step to step.
!   The gas-aqueous partitioning would also change with LWC.
!
! Partition fractions for individual lumped species
!   (as fraction of the lumped sum) are also saved for output.
!
! RAINOUT:  The sum of aqueous species into the gas-master
!   includes removal through rainout, based on the parameter RAINFR.
!
! WET DEPOSITION is calculated and saved for each aqueous species,
!   but it can be calculated more easily by using the output
!   concentrations of aqueous species (in M/L) and the rainfall rate.
!    (in gas-phase units).
! Here, this is omitted. Wet deposition can be derived instead
!  from the average concentrations of aqueous species (xcav, M/L)
!  multiplied by the  rain volume (liters) equivalent to rainfr
!  (=c_rainfr * c_h2oliq (kg/L) * alt. thickness (m) * 10 (dm/m))
!
!
! Inputs:    Average species concentrations (xc),
!              FINAL/AVERAGE ratio (xcfinr)
!              final spec ratio (xcf)
!
! Outputs:   Modified species concentrations (xc)
!            Final species concentrations (c_xcout)
!            Wet deposition:  c_xcwdep
!
!
! Called by:    quadchem
!
! Calls to:     None.
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
     subroutine postlump
!
       implicit none
!
       kk = 1
       !
       !   WRITE FINAL SPECIES CONCENTRATIONS TO OUTPUT ARRAY:
       !    final conc. (xcout) = avg concentration (xc) * ratio (xcfinr)
       !     where ratio is always for gas-master species
       !    average concentrations also entered into output array (xcav)
       !
       do ic = 1 , c_nchem2
         ic1 = c_npequil(ic)
         c_xcav(kk,ic) = xc(kk,ic)
         c_xcout(kk,ic) = xc(kk,ic) * xcfinr(kk, ic1)
       end do
       !
       !  AQUEOUS SPECIES:
       !        SUM AQUEOUS SPECIES INTO GAS-MASTER SUM (in GAS UNITS)
       !       (Individual aqueous species, aqueous units, are preserved,
       !          but the gas-aq sum is the main output)
       do ic = 1 , c_nchem2
         !
         !    AQUEOUS SPECIES identified by npequil (pointer, aq-to-gas)
         !
         ic1 = c_npequil(ic)
         if ( ic1 /= ic ) then
           calpha(kk) = c_xcav(kk,ic)*c_h2oliq(kk)*avogadrl
           c_xcav(kk,ic1) = c_xcav(kk,ic1) + calpha(kk)
           c_xcout(kk,ic1) = c_xcout(kk,ic1) + calpha(kk)*xcfinr(kk,ic)
         end if
       end do
       !
       ! ------------------------------------
       ! LUMPED SPECIES ARE IDENTIFIED FROM ARRAY LUMP(I,J)
       ! READ FROM REACTION.DAT
       !
       do i = 1 , c_cdim
         ics = c_lump(i,1)
         if ( ics == 0 ) exit
         ic1 = c_lump(i,2)
         ic2 = c_lump(i,3)
         !
         ! SUM LUMPED SPECIES (IC1-IC2) INTO ACCUMULATOR (ICS)
         ! OPTION - INCLUDE RP, RL, PRODUCTION AND LOSS
         c_xcout(kk,ics) = d_zero
         c_rp(kk,ics) = d_zero
         c_rl(kk,ics) = d_zero
         do ic = ic1 , ic2
           c_xcout(kk,ics) = c_xcout(kk,ics) + c_xcout(kk,ic)
           c_rp(kk,ics) = c_rp(kk,ics) + c_rp(kk,ic)
           c_rl(kk,ics) = c_rl(kk,ics) + c_rl(kk,ic)
         end do
         !
         ! OPTION:  CONVERT LUMPED SPECIES VALUES (IC1-IC2) INTO PARTITION FR.
         !
         do ic = ic1 , ic2
           kk = 1
           if ( c_kkw > 0 ) kk = c_kkw
           if ( c_xcout(kk,ics) > 0 ) then
!            c_xcout(kk,ic) = c_xcout(kk,ic)/c_xcout(kk,ics)
           end if
         end do
       end do
     end subroutine postlump
!
! -------------------------------------
!
! This sets reaction rates for parameterized RO2-RO2 reactions
!  as described in the CBM-Z mechanism (Zaveri, JGR, 1998)
!
! The parameterization represents an ensemble of reactions
!     RO2i + RO2j -> PRODi + PRODj   (rate constant kij)
! These are represented by reactions
!     RO2i -> PRODi
!  with pseudo-1st-order rate constant equal to
!     kii*[RO2i] +  sum over j (kij*[RO2j])
!
!  (the self-reaction RO2i+RO2i is represented with a rate 2*kii
!   because it removes 2 RO2i and produces 2 PRODi)
!
! Rate constants kii are equal to ratero2(kk,ni,1) from chemrates.
! Rate constants kij = 2*sqrt(Ki*Kj) (from Zaveri, 1998)
!  where Ki and Kj are partial rate coefficients from RO2i and RO2j (i /
!  (ratero2(kk,ni,2) = sqrt(Ki))
!
! Rate constants are calculated using [RO2] from prior iteration
!  or initial estimate.
! For steady state and iter=1 only, assume that [RO2]=[HO2]/nnro2
!  (self-reaction only;  zero for cross-reactions)
!  otherwise  1st iter RO2 prior estimate = 0.
!
! The parameterized RO2 reaction also is identified as a self-reaction
!   in chemsolve.
!
! Note OPTION: double sensitivity to OH for PARAMETERIZED RO2
!              (see c_nrk(nr) == -13)
!
! ----------------------
!
! Inputs:    Prior species concentrations (from previous iter) (xc)
!            RO2 reaction list (nrro2)
!            and RO2 partial rate constants (ratero2)
!
! Outputs:   Rate constants (ratek) for parameterized RO2 reactions
!            Sums for lumped species (xc)
!
! Called by:    quadchem
!
! Calls to:     None.
!
! ---------------------------------------------
! History:
!   6/09 Written by Sandy Sillman based on Zaveri, cbmz.f
!
! -------------------------------------------------------------------
!
     subroutine setro2
!
       implicit none

       if ( c_nnrro2 == 0 ) return
       !
       kk = 1
       !
       ! LOOP FOR PARAMETERIZED RO2 REACTIONS
       !
       do i = 1 , c_nnrro2
         nr = c_nrro2(i)
         ic = c_reactant(nr,1)
         !
         ! SELF REACTION: Use [HO2] if stst and 1st iter, otherwise use [RO2]
         !   stst 1st iter, estimate counts self-reaction rate, not summed RO2
         !   and assumes summed RO2 = HO2 or HO2/2
         !
         if ( c_lsts(ic) .and. c_iter == 1 ) then
           calpha(kk) = xc(kk,c_nho2)/(0.5D0*dble(c_nnrro2))
         else
           calpha(kk) = xc(kk,ic)
         end if
         !
         ! PSEUDO-RATE CONSTANT FOR RO2-RO2 SELF-REACTION
         ! (Doubled so that RO2->PROD parameterized reaction represents RO2+RO2
         !
         ratek(kk,nr) = d_two*ratero2(kk,i,1)* calpha(kk)
         !
         ! LOOP TO ADD RO2-RO2 CROSS REACTIONS
         !
         do j = 1 , c_nnrro2
           if ( i /= j ) then
             nr1 = c_nrro2(j)
             ic1 = c_reactant(nr1,1)
             !
             ! RO2i+RO2j RATE CONSTANT = 2*sqrt(RKi*RKj),
             ! ratero2(kk,i,2) = sqrt(RKi
             ratek(kk,nr) = ratek(kk,nr) + d_two*ratero2(kk,i,2) * &
                            ratero2(kk,j,2)*xc(kk,ic1)
           end if
         end do
         !
         ! END LOOP FOR PARAMETERIZED RO2 REACTIONS
         !
       end do
!
     end subroutine setro2
!
! -------------------------------------
!
! THIS  CALCULATES (H+) AND AQUEOUS SPECIES CONCENTRATIONS
!   (MOLES PER LITER)  BASED ON PRIOR (H+).
!
! THE SOLUTION FOR (H+)  MUST BE RUN ITERATIVELY.
! EACH AQUASOLVE CALL REPRESENTS A SINGLE ITERATION.
!
! CHEMISTRY IS ASSUMED TO BE REPRESENTED BY DISSOCIATING SPECIES
! WITH CHAINS OF UP TO THREE DISSOCIATIONS
! (GAS->AQUEOUS, AQUEOUS->SINGLE ION, ION->DOUBLE ION)
! EACH RELEASING H+ AND OH-.

! THIS RESULTS IN A SINGLE FOURTH-ORDER EQUATION FOR H+
!  (0 = THE SUM OF H+, OH-, CATIONS AND ANIONS, ALL FUNCTIONS OF H+)
! WHICH IS SOLVED ITERATIVELY USING NEWTON-RAFSON.
! TYPICALLY IT IS CALLED FROM WITHIN THE GAS-PHASE ITERATION.
!
! INCLUDES LELIEVELD 1991,etc. FOR WATER DROPLET DIFFUSION
! (SEE 'NEWTON-RAFSON' AND 'Lelieveld' FOR DETAILS BELOW.)
!
! INCLUDES GAS-PHASE DIFFUSION MODIFICATION FOR PARTITIONING
!  WHERE GAS->AQ TRANSFER IS SLOWER THAN AQUEOUS LOSS.
!  (with ACCOMODATION coefficients - see ACCOMODATION)
!
! NEW ADD:  GAS/AQUEOUS PARTITIONING BASED ON GAS/AQUEOUS SPEED.
!
! NOTE:  DROPLET RADIUS, GAS DIFFUSION CONSTANT, HARD-WIRED. PI also.
!
! ----------------
! A NOTE ON UNITS:  WITHIN THIS LOOP, THE GAS-MASTER SPECIES
! IS CONVERTED INTO LIQUID UNITS (MOLES GAS PER LIQUID WATER).
!   (=MOLECULES/CM3 /(AVOGADRL*AQUA(KK)).
! THE HENRY'S LAW COEFFICIENTS (RATEH) WERE ALSO CONVERTED
!   (MULTIPLIED BY AVOGADRL*AQUA(KK).)
! THIS ALLOWS FOR EASY CONSERVATION OF THE GAS+AQUEOUS SPECIES SUM.
! GAS-MASTER SPECIES IS CONVERTED BACK TO MOLECULES/CM3 AT THE END.
!
! NOTE:  0.1-1e6 gram/cm3 is typical.
! ----------------
!
! Inputs:    Species concentrations (xc), chemistry, LWC
!
! Outputs:   Gas and aqueous concentrations (xc),
!               in mole/cm3 (gas) and mole/liter (aq)
!
!
! Called by:    quadchem
!
! Calls to:     None.
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------

     subroutine aquasolve
!
       implicit none
!
       character(len=8) :: tsum
       real(dp) :: xsum

       ! 1/[H+]
       real(dp) :: xhinv(c_kvec)
       ! 1/[OH]
       real(dp) :: xohinv(c_kvec)
       ! 1/Kw =1/[H+][OH-}
       real(dp) :: xkwinv(c_kvec)
       ! Test for LWC>0 in loop
       real(dp) :: xtest
       ! Variable for write statement
       real(dp) :: xxx1
       ! Variable for write statement
       real(dp) :: xxx2
       ! Aq conversion fac for write
       real(dp) :: acquacon
!
       kk=1
       !
       ! IF ACQUA=0, RETURN.  ONLY FOR NONVECTORIZED VERSION.
       !
       xtest = d_zero
       !
       if ( c_h2oliq(kk) > xtest ) xtest = c_h2oliq(kk)
       if ( xtest <= 1.000D-25 ) return
       !
       ! NEW PRELIMINARY:  CALCULATE ADJUSTMENT TO HENRY'S LAW COEFFICIENTS
       !  TO REPRESENT DROPLET DIFFUSION LIMITATION
       !  (Lelieveld, J. At. Chem 12, 229, 1991 - see p. 241).
       !
       !  Q = Cavg/Csurf = 3 (coth q /q - 1/q**2)**-1;  q=r*(ka/Da)**0.5
       !
       !   where r= droplet radius (.001cm)  Da=droplet diffusion (2e-5 cm)
       !   ka = pseudo-1st-order aqueous loss rate.
       !  This Q applies only to the species component that originated in the
       !   gas phase and was transported to aqueous (Sg) as opposed to being
       !   produced chemically in the aqueous phase (Pa).
       !
       !  Total Q' = (Pa + QSg)/(Pa+Sg);  Sg=source to aqueous from gas.
       !
       !  (Pa, Sg should be from prior iteration.  Here,
       !   Sg = (rpgas+c_xcin) * prior aq/gas ratio (= rlaq/(rlaq+rlgas).  Slig
       !   other options:  (i)= (Pgas+c_xcin)*(xraq/xrtot).  (ii)=rlaq-rpaq  )
       !
       ! ------------
       ! BEGIN LOOP FOR DIFFUSION-HENRY'S LAW MODIFICATION.
       ! ------------
       do nrh = 1 , c_nreach
         ic = c_henry(nrh,1)
         ich = c_ncequil(ic,1)
         !
         ! DIFFUSION-MODIFIED HENRY'S LAW COEFFICIENT: FIRST SET EQUAL TO ONE
         ! ALSO, SAVE PRIOR RHDIF
         prior(kk) = rhdif(kk,nrh)
         rhdif(kk,nrh) = d_one
         !
         ! FOR ITER>1, USE PRIOR RP, RL TO SET DIFFUSION-MODIFICATION.
         !
         if ( c_iter > 1 .and. c_nequil(ic) > 0 )  then
           !
           ! PRELIMINARY ZERO
           !
           do i = 1 , 3
             rpro(kk,i) = d_zero
             rloss(kk,i) = d_zero
             xrp(kk,i) = d_zero
!            rpro(kk,i) = 0.00001
!            rloss(kk,i) = 0.00001
!            xrp(kk,i) = 0.1
           end do
           !
           ! SUM AQUEOUS CONCENTRATIONS IN GAS-EQUIVALENT UNITS (xrp1)
           ! ALSO SUM AQUEOUS PRODUCTION (rpro1) AND LOSS (rloss1)
           do neq = 1 , c_nequil(ic)
             icq = c_ncequil(ic,neq)
             if ( icq > 0 ) then
               if ( c_h2oliq(kk) > 0 ) then
                 if ( c_iter <= 2 ) then
                   xrp(kk,1) =xrp(kk,1) + xc(kk,icq)*c_h2oliq(kk)*avogadrl
                 else
                   ! xclastq = from prior iteration (=xc end of last aquasolve)
                   !
                   ! OPTION and possible ERROR:
                   !    (1)  Why use xclastq and not current xr?
                   !    (2)  Were rp preserved from prior iteration?
                   !
                   xrp(kk,1) = xrp(kk,1) + xc(kk,icq)*c_h2oliq(kk)*avogadrl
                 end if
                 rpro(kk,1) = rpro(kk,1) + c_rp(kk,icq)
                 rloss(kk,1) = rloss(kk,1) + c_rl(kk,icq)
                 if ( neq == 1 ) xrp(kk,2) = xrp(kk,1)
               end if
             end if
           end do
           !
           ! PSEUDO-FIRST-ORDER AQUEOUS LOSS CONSTANT (calpha)
           !  (NOTE:  if RL and XR=0, initial values above make lifetime long.)
           !
           if ( c_h2oliq(kk) > 0 ) then
             calpha(kk) = 0.00001D0
             if ( xrp(kk,1) > d_zero .and. rloss(kk,1) > d_zero) then
               calpha(kk) = rloss(kk,1)/(xrp(kk,1)*c_time)
             end if
           end if
           !
           ! Lelieveld Q-FACTOR FOR AQUEOUS DIFFUSION  (cbeta)
           !
           if ( c_h2oliq(kk) > d_zero ) then
             cgamma(kk) = c_droplet(kk) * dsqrt(calpha(kk)/dropdif)
             !
             ! PROTECT AGAINST EXTREME q (=droplet ratio)
             !
             if ( cgamma(kk) < 0.01D0 ) then
               cbeta(kk) = d_one
             else
               if ( cgamma(kk) > 100.0D0 ) then
                 cbeta(kk) = 0.001D0
               else
                 cbeta(kk) = (dexp(cgamma(kk)) + dexp(d_zero-cgamma(kk))) / &
                            (dexp(cgamma(kk)) - dexp(d_zero-cgamma(kk)))
                 if ( c_kkw > d_zero ) then
                   xxx1 = cbeta(c_kkw)
                 end if
                 cbeta(kk) = d_three*(cbeta(kk)/cgamma(kk) - &
                            d_one/(cgamma(kk)**d_two) )
                 if ( c_kkw > d_zero ) then
                   xxx2 = cbeta(c_kkw)
                 end if
                 if ( cbeta(kk) > d_one ) cbeta(kk) = d_one
                 if ( cbeta(kk) < d_zero ) cbeta(kk) = d_zero
               end if
             end if
           end if
           !
           ! Lelieveld Q-FACTOR, ADJUSTED FOR AQUEOUS PRODUCTION VS DIFFUSION
           ! FROM Q APPLIED TO GAS DIFFUSION ONLY:
           !    Q' = (Sgas*Q + Saq)/(Sgas+Saq)
           !
           ! (Final option for S=source of aqueous from gas:
           !  S = (rpgas+c_xcin) * aq/gas ratio (= rlaq/(rlaq+rlgas).  Slight overe
           !
           !  alternatives: (i)  aq/gas=1; aq/gas=Ca/Ct where H=Ca/Cg
           !               (ii) S = aqueous rl-rp if >0.
           !
           if ( c_h2oliq(kk) > d_zero ) then
             cgamma(kk) = (c_xcin(kk,ic) + c_rp(kk,ic))
             if ( rloss(kk,1) > d_zero ) then
               cgamma(kk) = (c_xcin(kk,ic) + c_rp(kk,ic)) * &
                            rloss(kk,1)/(rloss(kk,1)+c_rl(kk,ic))
             end if
             !
             ! Options
             ! cgamma(kk) = (  c_xcin(kk,ic) + c_rp(kk,ic))
             ! cgamma(kk) = (  c_xcin(kk,ic) + c_rp(kk,ic)) *rateh(kk,nrh) /
             !              (rateh(kk,nrh)+1.)
             ! cgamma(kk) = rloss(kk,1) - rpro(kk,1)
             ! if(cgamma(kk) < 0) cgamma(kk) == 0.)
             !
             ! Q adjustment:  Q'=rhdif;  Q=cbeta; Sgas=cgamma; Saq=rpro
             !  Apply only if aqueous concentration is not zero.
             !
             if ( xrp(kk,1) > d_zero .and. &
                  (cgamma(kk)+rpro(kk,1) > d_zero) ) then
               rhdif(kk,nrh) = (cbeta(kk)*cgamma(kk) + rpro(kk,1)) / &
                               (cgamma(kk) + rpro(kk,1))
             end if
           end if
         end if
         !
         ! END IF FOR ITER>1, AQUEOUS>0.
         !
         ! ENTER COMBINED HENRY-S LAW-w AQ DIFFUSION (relative units).
         !
         if ( c_h2oliq(kk) > d_zero ) then
           rhdif(kk,nrh) = rhdif(kk,nrh)*rateh(kk,nrh)
         end if
         !
         ! -----------------------------------------------------------
         ! GAS-PHASE DIFFUSION MODIFICATION:
         ! ADJUST GAS/AQUEOUS RATIO TO ACCOUNT FOR CASE WHERE AQUEOUS LOSS RATE
         !  (FROM AQUEOUS AND IONIC EQUILIBRIA SPECIES)
         !  IS FASTER THEN THE NEEDED GAS-AQUEOUS TRANSFER.
         ! -----------------------------------------------------------
         !
         !
         ! USES ACCOMODATION COEFFICIENTS, FORMULAS IN Lelieveld AND INFO
         !  FROM MARY BARTH.
         !
         ! PARTITION COEFFICIENT IS BASED ON STEADY STATE BETWEEN Pg, Pa,
         ! Lg, La, Eg, Ea=Eg/H.  Eg calculated from Lelieveld, Barth.
         ! Lg, La, Eg in s-1;  Pg, Pa in equivalent units.
         !
         ! (note rpro(kk,1) = aqueous production; rloss(kk,1)=aqloss,
         !  c_rp(kk,ic) = gas pro, 
         !  c_rl(kk,ic) = aqueous pro; xc(kk,ic) = gas cnc
         !  xrp(kk,2) = aqueous concentr (w/o ion sum), gas units.
         !
         !  (Pg-(Lg+Eg)Cg+Ea"Ca' = 0.;  Pa'+EgCg-(Ea'+La')Ca' = 0.
         !    where Ca', Pa', La',Ea' are for sum of aq-equil species;
         !    See hand notes in Lelieveld 1991)
         !
         ! AQUEOUS/GAS = (H*(Pa+Pg) + H(Lg/Eg)*Pa)/(Pa+Pg+H(La/Eg)*Pg)
         ! HENRY ADJUSTMENT = (Pa + Pg + (Lg/Eg)Pa)/(Pa+Pg+H(La/Eg)Pg)
         !
         ! *** ONLY IF c_iter>1. AND WITH SOFTENING.
         ! ----------------------------------------------------
         !
         ! ESTABLISH GAS=>AQUEOUS EXCHANGE RATE (Eg).
         !  MOLEC SPEED (VMOLEC, cgamma) = sqrt(8*ru*temp/(pi*c_molwt))
         !    = 3e4 cm/s, speed of sound.  (Barth).
         !
         !   Eg (s-1) = [(DROPLET**2/3DIFGAS) + 4DROPLET/(3*VMOLEC*ACCOM)]**-1
         !   (Lelieveld).
         !   THEN MULTIPLY BY LIQUID WATER CONTENT (acqua)
         !
         ! (STILL WITHIN HENRY'S LAW LOOP)
         ! ----------------------------------------------------
         !
         ! OPTION:  INCLUDE OR OMIT, WITH ITER CONTROL HERE.
         !
         if ( c_iter > 1 ) then
           cgamma(kk) = dsqrt(8.0D0*rumolec*c_temp(kk)/(mathpi*c_molwt(nrh)))
           egasaq(kk,nrh) = c_h2oliq(kk) / &
                ((c_droplet(kk)**d_two/(d_three*difgas)) + &
                 (d_four*c_droplet(kk)/(d_three*cgamma(kk)*c_accom(nrh))))
           cbeta(kk) = rpro(kk,1) + c_rp(kk,ic)
           !
           ! GOT HERE
           ! 2009 CORRECTION
           ! THESE NEXT LINES WERE COMMENTED OUT IN quadchv7.f
           !   WHICH GENERATES A SUCCESSFUL SOLUTION
           !   (w/ these lines commented in GAS DIFF ADJUST calpha = 1.000)
           !
!          if ( cbeta(kk) < 0. ) cbeta(kk) = 0.
!          cgamma(kk) = cbeta(kk)
!          if ( c_rl(kk,ic) > 0 .and. (xc(kk,ic)-xrp(kk,1) ) > 0 .and. &
!               rpro(kk,1) > 0 .and. egasaq(kk,nrh) > 0 ) then
!            cbeta(kk) = cbeta(kk) + &
!                      ( ( c_rl(kk,ic)/(xc(kk,ic)-xrp(kk,1)) ) / &
!                        (c_time*egasaq(kk,nrh)) ) *rpro(kk,1)
!
           !
           ! OLD  CRASH HERE
           !
           if ( rloss(kk,1) > d_zero .and. (xrp(kk,2)) > d_zero .and.    &
                c_rp(kk,ic) > d_zero .and. egasaq(kk,nrh) > d_zero .and. &
                dabs(rhdif(kk,nrh)) < dlowval ) then
             cgamma(kk) = cgamma(kk) + (((rloss(kk,1)/(xrp(kk,2))) / &
                     (c_time*egasaq(kk,nrh))*rhdif(kk,nrh)))*c_rp(kk,ic)
           end if
           calpha(kk) = d_one
           if ( cbeta(kk) > d_zero .and. cgamma(kk) > d_zero ) then
             calpha(kk) = cbeta(kk)/cgamma(kk)
           end if
           rhdif(kk,nrh) = rhdif(kk,nrh)*calpha(kk)
         end if
         !
         ! END GAS-PHASE DIFFUSION MODIFICATION:
         !
         ! MODIFY GAS/AQ RATIO FOR TROUBLESOME CONVERGENCE
         !
         if ( c_iter > 2 ) then
           rhdif(kk,nrh) = (rhdif(kk,nrh)**0.5D0)*(prior(kk)**0.5D0)
         end if
       end do
       ! -----------------------------------------------------------
       ! END LOOP FOR GAS-MASTER SPECIES FOR AQUEOUS DIFFUSION MODIFICATION
       ! -----------------------------------------------------------
       !
       ! PRELIMINARY:  ZERO XR PRIOR (USED IN AQUEOUS SUMS). INITIALIZE H+,OH-.
       ! SUM PRIOR GAS AND AQUEOUS SPECIES INTO THE GAS-MASTER
       ! ALSO SUM PRODUCTION (RP in gas units)
       !        RPRO1 =  SUMMED NET RP FOR THE INDIV. SPECIES
       !        RPRO2 =  RUNNING SUM OF NET RP FOR IONS.
       !        RPRO3 =  RUNNING SUM OF PRIOR ION CONCENTRATION
       !         RPRO2/3, EFFECT OF PRIOR CHEM PRODUCTION OF IONS
       !         WILL BE FACTORED INTO BETA (d/dH slope).
       !
       if ( c_h2oliq(kk) > d_zero ) then
         rpro(kk,1) = d_zero
         rpro(kk,2) = d_zero
         rpro(kk,3) = d_zero
         if ( xc(kk,c_nhplus) <= d_zero ) then
           xc(kk,c_nhplus) = 1.0D-05
           xc(kk,c_nohmin) = rateq(kk,1)/xc(kk,c_nhplus)
         end if
         if ( xc(kk,c_nohmin) <= d_zero ) then
           xc(kk,c_nohmin) = rateq(kk,1)/xc(kk,c_nhplus)
         end if
       end if
       do nrh = 1 , c_nreach
         ic = c_henry(nrh,1)
         if ( c_nequil(ic) > 0 )  then
           if ( c_h2oliq(kk) > 0 ) then
             rpro(kk,1) = c_rp(kk,ic)-c_rl(kk,ic)
           end if
           !
           ! SUM AQUEOUS SPECIES INTO THE GAS-MASTER, SUM GAS+AQUEOUS RPRO1.
           !
           do neq = 1 , c_nequil(ic)
             icq = c_ncequil(ic,neq)
             if ( icq > 0 ) then
               if ( c_h2oliq(kk) > 0 ) then
                 xc(kk,ic) = xc(kk,ic) + xc(kk,icq)*c_h2oliq(kk)*avogadrl
                 rpro(kk,1) = rpro(kk,1) + c_rp(kk,icq)-c_rl(kk,icq)
               end if
             end if
           end do
           !
           ! SUM PRIOR NET CHEM ION PRODUCTION (RPRO2) AND ION SUM (RPRO3)
           ! ION CHEM. INFERRED FROM PRIOR GAS-AQ-ION PARTITIONING AND RPRO1.
           ! ALSO ZERO AQUEOUS CONCENTRATIONS
           !
           do neq = 1 , c_nequil(ic)
             icq = c_ncequil(ic,neq)
             if ( icq > 0 ) then
               if ( c_iter > 1 ) then
                 if ( c_h2oliq(kk) > 0 ) then
                   if ( xc(kk,ic) > 0 ) then
                     rpro(kk,3) = rpro(kk,3) + (rpro(kk,1)/xc(kk,ic)) * &
                                  iabs(c_ion(icq))*c_h2oliq(kk)*avogadrl * &
                                  xc(kk,icq)
                   end if
                   rpro(kk,2) = rpro(kk,2) + iabs(c_ion(icq)) * &
                                 c_h2oliq(kk)*avogadrl*xc(kk,icq)
                 end if
               end if
               !
               ! ZERO AQUEOUS CONCENTRATIONS
               !
               xc(kk,icq) = d_zero
             end if
           end do
         end if
       end do
       !
       ! --------------------------------------------------
       ! CALCULATE AQUEOUS CONCENTRATIONS BASED ON PRIOR H+.
       ! ALSO CALCULATE ACIDSUM, NET SUM OF AQUEOUS IONS (=ALPHA)
       !  AND ACIDDEL, d(ACIDSUM)/dH+  (=BETA).
       ! THESE ARE THE NEWTON-RAFSON PARAMETERS.
       ! -------------------------------------------------------------
       ! FOR A GAS+AQUEOUS GROUP Xg, Xa, X1, X2 w/ constants Kh,K1, K2
       ! where Xa = Kh*Xg  X1*H = Xa*K1  X2*H=X1*K2  and Kw=H*OH:
       !
       ! UNITS:  Kh was converted from (MOLES/LITER)/ATMOSPHERE
       ! to (MOLES/LITER)/(MOLEC/CM3) to (RELATIVE UNITS) in BRATES.
       ! K1, K2 are in MOLES/LITER.  Xg is converted to MOLES/LITER here.
       ! The MOLES/LITER conversion is based on LIQUID WATER (acqua).
       !
       ! DOUBLE CATION:    (Xt=Xg+Xa+X1+X2)
       !  Xg = Xt/ (1 + Kh*(1 + K1/H*(1 + K2/H) ) )
       !  Xa = Xg*Kh   X1=Xa*K1/H  X2=X1*K2/H
       !           ( X1 = Xt*Rh*R1*H / (H**2*(1+Kh) + H*Kh*K1 + Kh*K1*K2)   )
       !  dX1/dH= X1/H - (X1**2/(Xt*Kh*K1*H) (2H(1+Kh)+Kh*K1)
       !           ( X2 = XtKhK1K2/(H**2(1+Kh)+ H*Kh*K1 + Kh*K1*K2)         )
       !  dX2/dH = (X2**2/(Xt*Kh*K1*K2)(2H(1+Kh)+Kh*K1)
       ! SINGLE CATION:    (K2=0)
       !  dX1/dt = -(X1**2/(Xt*Kh*K1*H))*(1+Kh)
       !
       ! DOUBLE ANION:
       ! Xg = Xt/ (1 + Kh*(1 + R1/OH*(1 + R2/OH) ) )
       ! Xa = Xg*Kh    X1=Xa*K1/OH   X2=X1*K2/OH
       !           ( X1=XtKhK1OH/( OH**2(1+Kh) + OHKhK1 + KhK1K2 )
       !           ( dX1/dH = -(OH/H)*dX1/dOH                                )
       ! dX1/dt= -X1/H + (X1**2/Xt*Kh*K1*H)*(2OH*(1+Kh)+Kh*K1) )
       !           ( X2=XtKhK1K2/(OH**2(1+Kh)+OH*K1*Kh + K2*K1*Kh )
       ! dX2/dH =  (OH/H) * (X2**2/Xt*Kh*K1*K2)*(2OH(1+Kh)*Kh*K1)
       !
       ! SINGLE ANION:
       ! dX1/dH =  (OH/H)* (X1**2)/(Xt*Kh*K1)*(1+Kh)
       !        (   dX1/dH = (X1*Xg/Xt*H)*(1+Kh)                             )
       !
       ! -------------------------------------------------------------
       !
       ! SET INITIAL ALPHA (=ion sum) AND BETA (=d/dH)
       ! TO INCLUDE THE IMPACT OF (H+) and (OH-).  (Note:  H+>0, above.)
       !
       ! SPECIAL FUNCTIONS:  1/H=xhinv(kk)   1/OH = xohinv(kk)
       !                     1/Kw = xkwinv(kk)  1/Xt = cgamma(kk)
       !
       if ( c_h2oliq(kk) > 0 ) then
         xkwinv(kk) = d_one/rateq(kk,1)
         xhinv(kk) = xkwinv(kk)*xc(kk,c_nohmin)
         xohinv(kk) = xkwinv(kk)*xc(kk,c_nhplus)
         calpha(kk) = xc(kk,c_nohmin)  - xc(kk,c_nhplus)
         cbeta(kk) =  d_one + xc(kk,c_nohmin)*xhinv(kk)
       end if
       !
       ! BEGIN XR CALCULATION.
       ! LOOP THROUGH HENRY'S LAW TO IDENTIFY GAS-MASTER
       ! AND SOLVE FOR THE ASSOCIATED AQUEOUS GROUP.
       !
       do nrh = 1 , c_nreach
         ic = c_henry(nrh,1)
         ich = c_henry(nrh,2)
         if ( c_nequil(ic) <= 0 ) exit
         !
         ! ZERO COUNTER FOR ION CHARGE
         ! AND INDICES FOR FIRST, SECOND ACID  OR BASE REACTIONS
         !
         ionsum = 0
         ica1 = 0
         ica2 = 0
         icb1 = 0
         icb2 = 0
         nra1 = 0
         nra2 = 0
         nrb1 = 0
         nrb2 = 0
         !
         !  CONVERT GAS-MASTER SUM TO LIQUID-EQUIVALENT UNITS (MOLES/LITER)
         !  WITH ZERO-PROTECT
         !
         if ( c_h2oliq(kk) > 0 ) then
           xc(kk,ic) = xc(kk,ic)/(c_h2oliq(kk)*avogadrl)
         end if
         !
         ! HENRY'S LAW CALCULATION FOR SPECIES WITH NO AQUEOUS EQUILIBRIA
         !  (dimensionless H=Ca/Cg, Cg=Ct/(1+H), Ca=Ct*H/(1+H) )
         !
         if ( c_nequil(ic) == 1 ) then
           if ( c_h2oliq(kk) > 0 ) then
             cgamma(kk) = d_one/(rhdif(kk,nrh) + d_one)
             xc(kk,ich) = (xc(kk,ic)*rhdif(kk,nrh)) * cgamma(kk)
             xc(kk,ic) = xc(kk,ic) * (cgamma(kk)*(c_h2oliq(kk)*avogadrl))
           end if
           !
           ! LOOP FOR AQUEOUS EQUILIBRIA.  (c_nequil(ic) > 1)
           !
         else
           !
           ! PRELIMINARY: IDENTIFICATION OF FIRST AND SECOND ACID-FORMING
           ! OR BASE-FORMING REACTIONS FOR THE AQUEOUS GROUP (ICA,NRA,ICB,NRB).
           ! THIS ESTABLISHES SOLUTION PROCEDURE FOR THE GROUP,
           ! TO BE USED BELOW.
           do neq = 2 , c_nequil(ic)
             icq = c_ncequil(ic,neq)
             nrq = c_nrequil(ic,neq)
             if ( icq > 0 .and. nrq > 0 ) then
               if ( c_aqueous(nrq,3) == c_nhplus ) then
                 if ( ica1 == 0 ) then
                   ica1 = icq
                   nra1 = nrq
                 else
                   ica2 = icq
                   nra2 = nrq
                 end if
               end if
               if ( c_aqueous(nrq,3) == c_aqueous(1,3) ) then
                 if ( icb1 == 0 )  then
                   icb1 = icq
                   nrb1 = nrq
                 else
                   icb2 = icq
                   nrb2 = nrq
                 end if
               end if
             end if
           end do
           !
           ! AQUEOUS IC, NR RECORDED FOR FIRST AQUEOUS REACTION
           !  (USED IF IT IS NOT ACID OR BASE REACTION).
           !
           icq = c_ncequil(ic,2)
           nrq = c_nrequil(ic,2)
           if ( icq <= 0 .and. nrq <= 0 ) then
             write(c_out,91) ic, c_tchem(ic), icq,nrq
           end if
           !
           ! MAIN SOLVER FOR AQUEOUS IONS.
           ! OPTIONS FOR: NEUTRAL ION, SINGLE OR DOUBLE CATION, ANION.
           !
           ! NEUTRAL AQUEOUS EQUILIBRIA:  GAS <-> AQUEOUS <-> NEUTRAL AQUEOUS.
           !   THIS WORKS ONLY FOR A SINGLE CHAIN:  ONE NEUTRAL ION.
           !
           ! Xa = Xg*Kh    X1=Xa*K1
           !
           if ( ica1 == 0 .and. icb1 == 0 ) then
             if ( c_h2oliq(kk) > 0 ) then
               cgamma(kk) = d_one/((rateq(kk,nrq)+d_one)*rhdif(kk,nrh)+d_one)
               xc(kk,icq) = xc(kk,ic)*(rateq(kk,nrq)*rhdif(kk,nrh))*cgamma(kk)
               xc(kk,ich) = (xc(kk,ic)*rhdif(kk,nrh)) * cgamma(kk)
               xc(kk,ic) = xc(kk,ic) * (cgamma(kk)*(c_h2oliq(kk)*avogadrl))
             end if
           end if
           !
           ! SINGLE CATION  - NOTE, SOLVE X1 FIRST, IN CASE Xgas, Xaq => 0.
           !  Xg = Xt/ (1 + Kh*(1 + K1/H            ) )
           !  Xa = Xg*Kh   X1=Xa*K1/H
           !  dX1/dt = -(X1**2/(Xt*Kh*K1*H))*(1+Kh)
           ! NOTE:  SOLVE X1 FIRST, INCASE Xgas, Xaq =>0.
           ! ALSO:  Initially use XR(IC) = TOTAL GAS+AQUEOUS, in aq. units.
           ! Then solve for GAS at the end.
           !
           if ( ica1 > 0 .and. ica2 == 0 ) then
             if ( c_h2oliq(kk) > 0 .and. xc(kk,ic) > 1.0D-40 ) then
               cgamma(kk) = d_one / &
                 ((rateq(kk,nra1)*xhinv(kk)+d_one)*rhdif(kk,nrh)+d_one)
               xc(kk,ica1) = (xc(kk,ic)*xhinv(kk) * &
                           rhdif(kk,nrh)*rateq(kk,nra1))*cgamma(kk)
               xc(kk,ich) = (xc(kk,ic)*rhdif(kk,nrh)) * cgamma(kk)
               !
               ! ION SUM (calpha), d/dH ION SUM (cbeta):
               !
               calpha(kk) = calpha(kk) - xc(kk,ica1)*c_ion(ica1)
               cbeta(kk) = cbeta(kk) - c_ion(ica1)*((d_one + rhdif(kk,nrh)) * &
                           (xc(kk,ica1)/xc(kk,ic))*(xc(kk,ica1) / &
                           (rhdif(kk,nrh)*rateq(kk,nra1))))
               !
               ! GAS SUM, converted to GAS units.
               !
               xc(kk,ic) = xc(kk,ic) * (cgamma(kk)*(c_h2oliq(kk)*avogadrl))
             end if
           end if
           !
           ! DOUBLE CATION
           !  Xg = Xt/ (1 + Kh*(1 + K1/H*(1 + K2/H) ) )
           !  Xa = Xg*Kh   X1=Xa*K1/H  X2=X1*K2/H
           !  dX1/dH= X1/H - (X1**2/(Xt*Kh*K1*H) (2H(1+Kh)+Kh*K1)
           !  dX2/dH = (X2**2/(Xt*Kh*K1*K2)(2H(1+Kh)+Kh*K1)
           !
           if ( ica1 > 0 .and. ica2 > 0 ) then
             if ( c_h2oliq(kk) > 0 .and. xc(kk,ic) > 1.0D-40 ) then
               cgamma(kk) = d_one/(((rateq(kk,nra2)*xhinv(kk) + d_one) * &
                            rateq(kk,nra1)*xhinv(kk)+d_one)*rhdif(kk,nrh)+d_one)
               xc(kk,ica1) = (xc(kk,ic)*xhinv(kk) * &
                     rhdif(kk,nrh)*rateq(kk,nra1))*cgamma(kk)
               xc(kk,ica2) = xc(kk,ica1)*rateq(kk,nra2)*xhinv(kk)
               xc(kk, ich) = (xc(kk,ic)*rhdif(kk,nrh)) * cgamma(kk)
               !
               ! ION SUM (calpha), d/dH ION SUM (cbeta):
               !
               calpha(kk) = calpha(kk) - xc(kk,ica1)*c_ion(ica1) - &
                           xc(kk,ica2)*c_ion(ica2)
               cbeta(kk) = cbeta(kk) + c_ion(ica1)*xc(kk,ica1)*xhinv(kk) - &
                          (d_two*xc(kk,c_nhplus)*(d_one+rhdif(kk,nrh)) + &
                          rhdif(kk,nrh)*rateq(kk,nra1)) * &
                        (c_ion(ica1)*(xc(kk,ica1)/xc(kk,ic)) * &
                        (xc(kk,ica1)/(rhdif(kk,nrh)*rateq(kk,nra1) * &
                         xc(kk,c_nhplus)))+c_ion(ica2) * &
                        (xc(kk,ica2)/xc(kk,ic))*(xc(kk,ica2) / &
                        (rhdif(kk,nrh)*rateq(kk,nra1)*rateq(kk,nra2))))
               !
               ! GAS SUM, converted to GAS units.
               !
               xc(kk,ic) = xc(kk,ic)*(cgamma(kk)*(c_h2oliq(kk)*avogadrl))
             end if
           end if
           !
           ! SINGLE ANION
           ! Xg = Xt/ (1 + Kh*(1 + R1/OH             ) )
           ! Xa = Xg*Kh    X1=Xa*K1/OH
           ! dX1/dH =  (OH/H)* (X1**2)/(Xt*Kh*K1)*(1+Kh)
           !
           if ( icb1 > 0 .and. icb2 == 0 ) then
             if ( c_h2oliq(kk) > 0 .and. xc(kk,ic) > 1.0D-40 ) then
               cgamma(kk) = d_one / &
                 ((rateq(kk,nrb1)*xohinv(kk)+d_one)*rhdif(kk,nrh)+d_one)
               xc(kk,icb1) = (xc(kk,ic)*xohinv(kk) * &
                      rhdif(kk,nrh)*rateq(kk,nrb1))*cgamma(kk)
               xc(kk,ich) = (xc(kk,ic)*rhdif(kk,nrh)) * cgamma(kk)
               !
               ! ION SUM (calpha), d/dH ION SUM (cbeta):
               !
               calpha(kk) = calpha(kk) - xc(kk,icb1)*c_ion(icb1)
               cbeta(kk) = cbeta(kk) + & 
                 c_ion(icb1)*((xc(kk,c_nohmin)*xhinv(kk)) * &
                 (d_one+rhdif(kk,nrh))*(xc(kk,icb1)/xc(kk,ic)) * &
                 (xc(kk,icb1)/(rhdif(kk,nrh)*rateq(kk,nrb1))))
               !
               ! GAS SUM, converted to GAS units.
               !
               xc(kk,ic) = xc(kk,ic)*(cgamma(kk)*(c_h2oliq(kk)*avogadrl))
             end if
           end if
           !
           ! DOUBLE ANION
           ! Xg = Xt/ (1 + Kh*(1 + R1/OH*(1 + R2/OH) ) )
           ! Xa = Xg*Kh    X1=Xa*K1/OH   X2=X1*K2/OH
           ! dX1/dt= -X1/H + (X1**2/Xt*Kh*K1*H)*(2OH*(1+Kh)+Kh*K1) )
           ! dX2/dH =  (OH/H) * (X2**2/Xt*Kh*K1*K2)*(2OH(1+Kh)*Kh*K1)
           !
           if ( icb1 > 0 .and. icb2 > 0 ) then
             if ( c_h2oliq(kk) > 0 .and. xc(kk,ic) > 1.0D-40 ) then
               cgamma(kk) = d_one/(((rateq(kk,nrb2)*xohinv(kk)+d_one) * &
                  rateq(kk,nrb1)*xohinv(kk)+d_one)*rhdif(kk,nrh)+d_one)
               xc(kk,icb1)  = xc(kk,ic) * &
                 ((xohinv(kk)*rhdif(kk,nrh)*rateq(kk,nrb1))*cgamma(kk))
               xc(kk,icb2) = xc(kk,icb1)*(rateq(kk,nrb2)*xohinv(kk))
               xc(kk,ich) = xc(kk,ic)*(rhdif(kk,nrh) * cgamma(kk) )
               !
               ! ION SUM (calpha), d/dH ION SUM (cbeta):
               !
               calpha(kk) = calpha(kk) - xc(kk,icb1)*c_ion(icb1) - &
                           xc(kk,icb2)*c_ion(icb2)
               cbeta(kk) = cbeta(kk) - c_ion(icb1) *xc(kk,icb1)*xhinv(kk)
               cbeta(kk) = cbeta(kk) + &
                 (d_two*xc(kk,c_nohmin)*(d_one+rhdif(kk,nrh)) + &
                  rhdif(kk,nrh)*rateq(kk,nrb1)) * &
                 (c_ion(icb1)*(xc(kk,icb1)/xc(kk,ic))*(xc(kk,icb1) / &
                  (rhdif(kk,nrh)*rateq(kk,nrb1)*xc(kk,c_nhplus))) + &
                  c_ion(icb2)*(xc(kk,icb2)/xc(kk,ic)) * &
                  (xc(kk,c_nohmin)*xhinv(kk))*(xc(kk,icb2) / &
                  (rhdif(kk,nrh)*rateq(kk,nrb1)*rateq(kk,nrb2))))
               !
               ! GAS SUM, converted to GAS units.
               !
               xc(kk,ic) = xc(kk,ic) * (cgamma(kk)*(c_h2oliq(kk)*avogadrl))
             end if
           end if
           !
           ! 2006 ERROR CORRECTION - PROTECT AGAINST ZERO AQUEOUS - CUT
           !     if(xc(1,ich) == 0.) then
           !       if(rhdif(1,nrh) > 1.0e-10) xc(1,ich) = 1.0e-34
           !     end if                        !if(xc(kk,ich) == 0.) then
         end if
       end do
       !
       ! ----------------------------------------------
       ! END OF LOOP TO CALCULATE AQUEOUS CONCENTRATIONS
       ! ----------------------------------------------
       !
       ! -------------------------------
       ! CALCULATE H+ (AND OH-) FROM NEWTON-RAPHSON
       ! -------------------------------
       !
       ! MODIFICATION FOR DIFFICULT CONVERGENCE: H+=geometric mean w/ prior H+
       !   (To prevent oscillation with H+/SO2->HSO3-->SO4=.)
       !
       ! ALSO:  INCREASE BETA BY RATIO:  ION RPRO/PRIOR ION SUM
       ! TO ACCOUNT FOR CHEM PRODUCTION ->H+ FEEDBACK.
       !
       if ( c_h2oliq(kk) > 0 ) then
         if ( rpro(kk,2) > rpro(kk,3) .and. rpro(kk,3) > 0 ) then
           cbeta(kk) = cbeta(kk)*(d_one+rpro(kk,3)/(rpro(kk,2)-rpro(kk,3)))
         end if
         if ( dabs(cbeta(kk)) < dlowval ) calpha(kk) = calpha(kk)/cbeta(kk)
         if ( calpha(kk) < -0.99D0*xc(kk,c_nhplus)) then
           calpha(kk) = -0.9D0*xc(kk,c_nhplus)
         end if
         if ( calpha(kk) > 100.0D0*xc(kk,c_nhplus)) then
           calpha(kk) = 10.0D0*xc(kk,c_nhplus)
         end if
         xc(kk,c_nhplus) = dsqrt(xc(kk,c_nhplus)*(xc(kk,c_nhplus)+calpha(kk)))
         if ( xc(kk,c_nhplus) > 0 ) then
           xc(kk,c_nohmin) = rateq(kk,1)/xc(kk,c_nhplus)
         end if
       end if

   91  format(/,' CHEMISTRY INDEX ERROR IN AQUASOLVE:',/,   &
                ' IC, TCHEM =  ', i5, 2x,a8,/,              &
                ' 2ND AQUEOUS IC, NR <=0 IN AQUEOUS LOOP; = ',2i5)

     end subroutine aquasolve
!
end module mod_cbmz_solve1
