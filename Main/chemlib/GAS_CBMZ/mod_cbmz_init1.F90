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

module mod_cbmz_init1
!
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_mpmessage
  use mod_cbmz_chemmech
  use mod_cbmz_chemlocal
  use mod_cbmz_chemvars
  use mod_cbmz_jval1
!
  private
!
  public :: chemread , namechem , hvread , cheminit
!
  contains
!
!  cheminit.f   April, 2007
!
!    4-2009: error with ifort, but not with -C compile option.
!      indices c_noh, etc. are written incorrectly, possibly related to
!       warning message about real(rk8) :: in COMMON.
!
!     Nov 2007 addition: save net RP stoich for tracers.
!  NOTE CHANGES:  rbchemmech.EXT, quadinit.f
!
!  Program to read input, do initial set-up, and write output
!   for the RADICAL BALANCE-BACK EULER solution for photochemistry.
!
!  Full notes on development in file  quadchem.f (chemsol.f)
!   and older boxchemv7.f
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
!
!     Time test:
!     Future option: nonsteady  state aqueous
!     Future option: exponential decay  solution
!     Long-term option: integrate with aerosol solution
!
!  Subroutines:
!    chemread:  reads input for mechanism: species, reactions, etc.
!    hvread:    Reads data for  j-value calculations
!    cheminit (quadinit) :  Sets up variables relating to
!                chem. mechanism for solver
!                (settings that remain unchanged when solver is called)
!    chemwrite.f: Write complete output for chemistry
!    analyze.f:   Writes output relating to specif chemical species
!
!  FUTURE ADD:  Identified below by 'FUTURE ADD'.
!     Exponential decay solution option
!     Non-steady state gas/aqueous partitioning
!     Categories for aqueous/soluble aerosol species
!           (for non-steady-state gas/aq partitioning)
!     Vectorize the remaining parts of the program:  hv rates.
!
!     Modify to skip aqueous reactions if LWC=0:
!      assign reactants to direct species, not gas pointer
!      and in chem, do reactions only when species is called
!           - if LWC=0 skip aqueous
!
!
!  Program  written by Sanford Sillman
!  ---------------------------------------------------------
!  Program history:
!    12/06 Initial program by Sandy Sillman based on boxchemv7.f
!  ---------------------------------------------------------
!
! This reads the chemistry input file (REACTION.DAT)
!   which includes:
!    species names, steady-state indices, and 'lumping' for I/O
!    aqueous equilibria
!    reaction rates
!    reaction list and stoichiometries
!    order  of solution ('cascade')

! Note:  photolysis rate input is in HVREAD.
!        box model input parameters are in BREDIN.
!
! Inputs:  None.
!
! Outputs:  All chemistry parameters listed above
!
!
! Called by:  boxmain (as part of chemical setup)
!
! Calls to:   cheminit
!
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
!
!
    subroutine chemread
!
      implicit none
      ! Chem index
      integer(ik4) :: ic , ic1 , ic2 , ic3 , icc , ics , icc1 , icc2 , icc3
      ! Chem index
      integer(ik4) :: icp
      ! chem index - pair and multi
      integer(ik4) :: icr1 , icr2
      ! Aqueous counters
      integer(ik4) :: neq
      ! Chem species counters
      integer(ik4) :: nc , nc1 , nc2 , ncf , nn , nne,nsolv,nsol
      ! Reaction counters
      integer(ik4) :: nr , nrh , nrq , nrqq , np
      ! Vectorization counters
      integer(ik4) :: kk
      ! General counters
      integer(ik4) :: i , j , ii , iii , n
!
      ! dummy input character variable
      character(len=8) :: tdum(5)
      ! dummy vbl to identify READ
      character(len=4) :: titl
      ! dummy input integer variable
      integer(ik4) :: ndum(5)
      ! dummy input real variable
      real :: xdum(5)
      !  caspair(ic,2)  Array to identify linked pair in cascade solver;
      ! 2nd species is linked in pair chain to 1st
      ! used to establish pair chain array (nppair)
      ! Array to ID paired species
      integer(ik4) :: caspair(c_cdim,2)
      ! Index for species number in cascade
      integer(ik4) :: ncas
      ! Integer for spec. number in caspair
      integer(ik4) :: ncasp
      !
      !  lcastest: Flag set to true when species is read in cascade.
      ! Flag to check species in cascade
      logical :: lcastest(c_cdim)
      !
      ! Counter for number of lumped species
      integer(ik4) :: nlump
      ! General counter
      integer(ik4) :: jj
      ! Added reaction counter (nr in chemlocal)
      integer(ik4) :: nnr

!
      ! Set vector variable for non-vectorized case
      kk = 1
      !
      ! PRELIMINARY ZERO FOR READ
      !
      do nr = 1 , c_rdim
        c_nnpro(nr) = 0
        do i = 1 , 20
          if ( i <= 2 ) c_reactant(nr,i) = 0
          if ( i <= 6 ) c_treac(i,nr) = '        '
          if ( i <= 2 .and. nr <= 61 ) c_treach(i,nr) = '        '
          if ( i <= 2 .and. nr <= 61 ) c_henry(nr,i) = 0
          if ( i <= 3 .and. nr <= 61 ) c_treacq(i,nr) = '        '
          if ( i <= 3 .and. nr <= 61 ) c_aqueous(nr,i) = 0
          c_product(nr,i) = 0
          c_stoich(nr,i) = d_zero
        end do
        do ic = 1 , c_cdim
          c_prodarr(nr,ic) = d_zero
        end do
      end do

      do i = 1 , c_cdim
        do ii = 1 , c_cdim
          c_cascade(i,ii) = 0
        end do
        do ii = 1 , 3
          c_lump(i,ii) = 0
        end do
        do ii = 1 , 2
          caspair(i,ii) = 0
        end do
        ! PAIR POINTERS to linked species in a pair chain:
        !  (c_nppair(i,1) points to directly linked pair above in chain     )
        !  (c_nppair(i,2) points to primary species for pair chain          )
        !  (c_nppair(i,3) number of pair group subspecies for this species  )
        !  (c_nppair(i,4+) chem. index (ic) for pair subspeice for this sp. )
        !
        ! initially point at self;
        !   counter and pointers to secondary species zero (see nequil)
        !
        c_nppair(i,1) = i
        c_nppair(i,2) = i
        do ii = 3 , 23
          c_nppair(i,ii) = 0
        end do
      end do
      !
      ! ZERO SPECIAL SPECIES INDICES
      !
      c_nh2o   = 0
      c_nco2   = 0
      c_nco    = 0
      c_nhplus = 0
      c_nohmin = 0
      c_noh    = 0
      c_nho2   = 0
      c_nh2o2  = 0
      c_no3    = 0
      c_nno2   = 0
      c_nno    = 0
      c_no3    = 0
      c_nn2o5  = 0
      c_nhno3  = 0
      c_nch4   = 0
      c_nco3   = 0
      !
      ! ZERO SPECIAL INDICES FOR PARAMETERIZED RO2-RO2 REACTIONS (CBMZ)
      !
      c_nnrro2 = 0
      do nr = 1 , c_rdim
        c_nrro2 = 0
      end do
      !
      ! READ LOOP:  Read char*4 until you find 'START'.
      !  Then BREAK and begin first read cycle.
      !
      do i = 1 , 1111
        read(c_rin,11) titl
        if ( titl == 'STAR' ) exit
      end do
      !
      ! READ PRELIMINARY INDICES:  VECTOR KKW, KMAX; NUMITER AND CONVERGE
      !
      read(c_rin,12) c_kkw
      read(c_rin,12) c_kmax
      read(c_rin,12) c_numitr
      if ( c_kkw == 5 ) write(c_out,12) c_kkw, c_kmax, c_numitr
      read(c_rin,13) c_converge
      if ( c_kkw == 5 ) write(c_out,13) c_converge
      !
      ! ------------------------
      ! READ LOOP: NAMES AND CATEGORIES FOR TRANSPORTED (INPUT/OUTPUT) SPECIES
      ! ------------------------
      ic = 0
      do i = 1 , 1111
        read(c_rin,11) titl
        if ( titl == 'STAR' ) exit
      end do
      if ( c_kkw == 5 ) write(c_out,11) titl
      !
      ! OPTION: READ IN GROUPS OF 1 OR 5:  ii=1,1 vs ii=1,5 THROUGHOUT
      !   (main species read here and lumped species read below).
      !   (also can change in summary write)
      !
      loopcdim: &
      do i = 1 , c_cdim
        read(c_rin,29) (tdum(ii) , ndum(ii) , ii=1,5)
        if ( c_kkw == 5 ) write(c_out,29) (tdum(ii) , ndum(ii), ii=1,5)
        !
        ! WHILE NAME NE 'END', ENTER DUMMY READ INTO TCHEM ARRAY.
        ! AT 'END' BREAK AND LEAVE LOOP
        !
        do ii = 1 , 5
          if ( tdum(ii) == '     END' .or. &
               tdum(ii) == '    END ' ) exit loopcdim
          ic = ic+1
          if ( ic == c_cdim ) then
            write(c_out,901) ic
            write(6,901) ic
          end if
          c_tchem(ic) = tdum(ii)
          c_icat(ic) = ndum(ii)
        end do
      end do loopcdim
      !
      ! READ INSERT:  RECORD NUMBER OF INPUT   SPECIES (NCHEM1)
      !
      c_nchem1 = ic
      if ( c_kkw > 0 ) write(c_out,49) c_nchem1
      !
      ! SUMMARY WRITE:  INPUT/OUTPUT SPECIES
      !
      write(c_out,1301)
      ncf = idint(0.2D0*(dble(c_nchem1)-0.01D0))+1
      do nc = 1 , ncf
        nc1 = 5*(nc-1) + 1
        nc2 = nc1 + 4
        if ( nc1 > c_nchem1 ) nc1 = c_nchem1
        if ( nc2 > c_nchem1 ) nc2 = c_nchem1
        write(c_out, 1302) (nn, c_tchem(nn),c_icat(nn) , nn=nc1,nc2 )
      end do
      !
      ! SET CATEGORY AND STEADY STATE INDEX (after WRITE)
      !    ASSIGN LSTS=T IF ICAT<0, MAKE CATEGORY POSITIVE.
      !    ALSO INITIALIZE LUMP FLAG=false
      !
      do ic = 1 , c_nchem1
        ! lump accumulator flag initially F
        c_llump(ic) = .false.
        c_lsts(ic) = .false.
        if ( c_icat(ic) < 0 ) c_lsts(ic) = .true.
        if ( c_icat(ic) < 0 ) c_icat(ic) = 0 - c_icat(ic)
      end do
      !
      ! RECORD INTERIM NUMBER OF TOTAL SPECIES (INPUT + INTERNAL) (NCHEM2)
      !  (This will be increased as other species are added)
      !
      c_nchem2 = c_nchem1
      !
      ! ------------------------
      ! END OF CHEMICAL NAME LOOP.  TOTAL NUMBER OF CHEMICALS=NCHEM2
      !  (to be moved below, after HENRY'S LAW and AQUEOUS EQULIBRIA)
      ! ------------------------
      !
      ! -------------------------------------------------------
      ! READ LUMPED SPECIES
      ! -------------------------------------------------------
      ! LUMP(IC, IC1, IC2): SPECIES FROM IC1 TO IC2 ARE LUMPED
      ! INTO ACCUMULATOR (IC), AND SHARE LSTS, ICAT WITH LUMPED SUM.
      ! -------------------------------------------------------
      !
      !  11/22/06
      ! NEW READ: List of individual lumped species read here.
      !
      ! NOTE: Individual lumped species are automatically added
      !       to the species list.
      ! The name of the lump group must be already included
      !      in the list of 'primary' (i/o) species above.
      !
      nlump = 0
      do i = 1 , 1111
        read(c_rin,11) titl
        if ( titl == 'STAR' ) exit
      end do
      if ( c_kkw == 5 ) write(c_out,11) titl
      loopindex: &
      do i = 1 , 20
        read(c_rin,209) (tdum(ii),ii=1,1)
        if ( c_kkw == 5 ) write(c_out,209) (tdum(ii),ii=1,1)
        if ( tdum(1) == '     END' .or. &
             tdum(1) == '    END ') exit loopindex
        !
        ! ADVANCE COUNTER, SET SPECIES INDICES AT ZERO.
        !
        nlump = nlump + 1
        c_lump(nlump,1) = 0
        c_lump(nlump,2) = 0
        c_lump(nlump,3) = 0
        !
        ! ENTER LUMP SPECIES NUMBER INTO LUMP ARRAY.
        ! IF LUMPED SPECIES  NAME IS UNKNOWN, ERROR AND EXIT LOOP.
        !
        ic = namechem(tdum(1))
        if ( ic == 0 ) then
          write(c_out,219) (tdum(j),j=1,1),ii
          write(6,219) (tdum(j),j=1,1),ii
          c_lump(nlump,1) = 0
          c_lump(nlump,2) = 0
          c_lump(nlump,3) = 0
          nlump = nlump-1
          cycle loopindex
        end if
        c_lump(nlump,1) = ic
        c_llump(ic) = .true.
        !
        ! READ LIST OF SPECIES TO BE INCLUDED IN LUMP, ADD TO SPECIES LIST.
        !
        looplumpspec: &
        do jj = 1 , 111
          read(c_rin,39) (tdum(ii),ii=1,5)
          if ( c_kkw == 5 ) write(c_out,39) (tdum(ii),ii=1,5)
          do ii = 1 , 5
            if ( tdum(ii) == '       x' .or. &
                 tdum(ii) == '        ' .or. &
                 tdum(ii) == '      x ' ) exit looplumpspec
            !
            !  ENTER SPECIES INTO SPECIES LIST AND INTO LUMP ARRAY.
            !     (lump2 = 1st lumped species, lump3 = last lumped species)
            !
            ic = namechem(tdum(ii))
            if ( ic == 0 ) then
              c_nchem2 = c_nchem2+1
              if ( c_nchem2 == c_cdim ) then
                write(c_out,901) c_nchem2
                write(6,901) c_nchem2
              end if
              c_tchem(c_nchem2) = tdum(ii)
              !
              !  ADD TO LUMP ARRAY.
              !
              if ( c_lump(nlump,2) == 0 ) c_lump(nlump, 2) = c_nchem2
              c_lump(nlump, 3) = c_nchem2
            else
              !
              !  TEST IF SPECIES IS ALREADY INCLUDED IN SPECIES LIST -
              !  IF LUMPED SPECIES ALREADY INCLUDED, ERROR AND EXIT LUMP LOOP
              write(c_out,218) (tdum(j),j=1,1)
              write(6,218) (tdum(j),j=1,1)
              c_lump(nlump,1) = 0
              c_lump(nlump,2) = 0
              c_lump(nlump,3) = 0
              nlump = nlump-1
              cycle loopindex
            end if
          end do
        end do looplumpspec
        !
        ! ENTER CATEGORY AND STEADY STATE FOR LUMPED SPECIES.
        ! SET LUMP FLAG= F (T for accumulator only)
        !
        do ic = c_lump(nlump,2) , c_lump(nlump,3)
          c_icat(ic) = c_icat(c_lump(nlump,1))
          c_lsts(ic) = c_lsts(c_lump(nlump,1))
          c_llump(ic) = .false.
        end do
        !
        ! SUMMARY WRITE: LUMPED SPECIES  SUMMARY
        !
        if ( nlump == 1 ) write(c_out, 1310)
        write(c_out,1311) c_tchem(c_lump(nlump,1))
        ncf = idint(0.2D0*(dble(c_lump(nlump,3) - &
                                c_lump(nlump,2)+d_one)-0.01D0))+1
        do nc = 1 , ncf
          nc1 = 5*(nc-1) + c_lump(nlump,2)
          nc2 = nc1 + 4
          if ( nc1 > c_lump(nlump,3) ) nc1 = c_lump(nlump,3)
          if ( nc2 > c_lump(nlump,3) ) nc2 = c_lump(nlump,3)
          write(c_out, 1302) (nn, c_tchem(nn),c_icat(nn) , nn=nc1,nc2 )
        end do
      end do loopindex
      !
      ! -------------------------------------------------------
      ! END LUMPED SPECIES READ
      ! -------------------------------------------------------
      !
      ! ------------------------------------
      ! READ LOOP:  HENRY'S LAW COEFFICIENTS
      ! ------------------------------------
      nrh = 0
      do i = 1 , 1111
        read(c_rin,11) titl
        if ( titl == 'STAR' ) exit
      end do
      if ( c_kkw == 5 ) write(c_out,11) titl
      !
      ! LOOP FOR HENRY'S LAW READS: CONTINUE UNTIL 'END'.
      !
      loophenry: &
      do i = 1 , c_cdim
        !
        ! ADVANCE COUNTER.
        !
        nrh = nrh+1
        ! --------------------------------
        ! READ NAMES OF HENRY'S LAW SPECIES.
        ! --------------------------------
        if ( nrh == 161 .or. nrh == 121 .or. nrh == c_cdim ) then
          write(c_out,903) ic
          write(6,903) ic
        end if
        read(c_rin, 59) tdum(1),tdum(2)
        if ( c_kkw == 5 ) write(c_out, 59) tdum(1),tdum(2)
        !
        ! ENTER NAMES INTO HENRY'S LAW MATRIX
        ! 'END' BREAKS THE LOOP TO END THE HENRY'S LAW READ.
        !
        do iii = 1 , 2
          if ( tdum(iii) == '     END' .or. &
               tdum(iii) == '    END ' ) exit loophenry
          c_treach(iii,nrh) = tdum(iii)
          ic = namechem(tdum(iii))
          !
          ! ADD AQUEOUS HENRY'S LAW SPECIES INTO SPECIES LIST
          ! (if not already ther
          !
          if ( ic == 0 .and. iii == 2 ) then
            ic = c_nchem2+1
            c_nchem2 = ic
            if ( c_nchem2 == c_cdim ) then
              write(c_out,901) c_nchem2
              write(6,901) c_nchem2
            end if
            c_tchem(c_nchem2) = tdum(iii)
          end if
          !
          ! ENTER SPECIES INTO HENRY LAW ARRAY
          !
          if ( ic > 0 ) then
            c_henry(nrh,iii) = ic
          else
            !
            ! TEST FOR POSSIBLE ERROR IN CHEMICAL NAMES.
            !  THE GAS-PHASE SPECIES SHOULD ALREADY BE IN THE SPECIES LIST.
            !
            if( tdum(iii) /= '        ' .and. &
                tdum(iii) /= '       x' .and. &
                tdum(iii) /= '       X' .and. &
                tdum(iii) /= '      x ' .and. &
                tdum(iii) /= '      X ' .and. &
                tdum(iii) /= '    MORE') then
              write(c_out, 69) tdum(iii), nrh,iii,(tdum(j),j=1,2)
              write(6, 69) tdum(iii), nrh,iii,(tdum(j),j=1,2)
            end if
          end if
        end do
        !
        ! DEFAULT ACCOMODATION COEFFICIENT AND MOLECULAR WEIGHT
        !
        c_accom(nrh) = 0.05D0
        c_molwt(nrh) = 30.0D0
        !
        !  READ RATE PARAMETERS
        !
        ! CURRENT OPTION: ACCOMODATION COEFFICIENT AND MOLECULAR WT READ.
        ! FUTURE MODIFICATION: ADD ALTERNATIVE RATE PARAMETER FORMATS IF NEEDED.
        !
        !  READ LINES WITH THE FOLLOWING  FORMAT:
        !    CO2 =     CO2L
        !   1    0 3.400E-02     2400. 5.000E-02       44.
        !
        read(c_rin, 33) c_nrkh(nrh), c_rkh(1,nrh), c_rkh(2,nrh), &
                        c_accom(nrh), c_molwt(nrh)
        if ( c_kkw == 5 ) write(c_out, 33)  &
           c_nrkh(nrh), c_rkh(1,nrh), c_rkh(2,nrh), &
           c_accom(nrh) , c_molwt(nrh)
      end do loophenry
      !
      ! --------------------------------------------------------
      ! END HENRY'S LAW LOOP.
      ! --------------------------------------------------------
      c_nreach = nrh-1
      if ( c_kkw > 0 ) write(c_out,57) c_nreach
      ! ------------------------------------
      ! READ LOOP:  AQUEUS EQUILIBRIUM CONSTANTS WITH H+
      ! ------------------------------------
      !
      nrq = 0
      do i = 1 , 1111
        read(c_rin,11) titl
        if ( titl == 'STAR' ) exit
      end do
      if ( c_kkw == 5 ) write(c_out,11) titl
      !
      ! LOOP FOR EQUILIBRIUM READS: CONTINUE UNTIL 'END'.
      !
      loopeq: &
      do i = 1 , c_cdim
        !
        ! ADVANCE COUNTER.
        !
        nrq = nrq+1
        if ( nrq == 161 .or. nrq == 121 .or. nrq == c_cdim) then
          write(c_out,904) nrq
        end if
        !
        ! --------------------------------
        ! READ NAMES OF H+ EQUILIBRIUM SPECIES:  A = B + C* (*C=H+or OH-)
        !    NOTE, 1st EQUILIBRIUM MUST BE 0<=>H+ + OH-
        ! --------------------------------
        !
        read(c_rin, 59) tdum(1), tdum(2), tdum(3)
        if ( c_kkw == 5 ) write(c_out, 59) tdum(1), tdum(2), tdum(3)
        !
        ! SWITCH 2nd AND 3rd SPECIES if 2nd SPECIES = H+ or OH-
        !
        if (nrq > 1) then
          if (tdum(2) == c_treacq(2,1) .or. &
              tdum(2) == c_treacq(3,1) .or. &
              tdum(2) == '      H+' .or. &
              tdum(2) == '     OH-' .or. &
              tdum(2) == '     H+ ' .or. &
              tdum(2) == '    OH- ') then
            if ( tdum(3) /= '        ' .and. &
                 tdum(3) /= '       x' .and. &
                 tdum(3) /= '       X' .and. &
                 tdum(3) /= '      x ' .and. &
                 tdum(3) /= '      X ' .and. &
                 tdum(3) /= '    MORE') then
              tdum(4) = tdum(2)
              tdum(2) = tdum(3)
              tdum(3) = tdum(4)
              if ( c_kkw > 0 ) write(c_out,88) nrq,tdum(1),tdum(2),tdum(3)
            end if
          end if
        end if
        !
        ! ENTER NAMES INTO AQUEOUS EQUILIBRIUM MATRIX
        ! 'END' BREAKS THE LOOP TO END THE EQUILIBRIUM READ.
        !
        do iii = 1 , 3
          if ( tdum(iii) == '     END' .or. &
               tdum(iii) == '    END ' ) exit loopeq
          c_treacq(iii,nrq) = tdum(iii)
          ic = namechem(tdum(iii))
          !
          ! ENTER NEW CATION/ANION INTO SPECIES LIST IF NOT ALREADY INCLUDED.
          !   (NOTE, FIRST EQUILIBRIUM SHOULD BE 0 <=> H+ + OH-;
          !    SUBSEQUENT SHOULD BE A <=> B + H+ or OH-)

          if ( ic == 0 .and. (iii == 2 .or. (iii == 3 .and. nrq == 1)) ) then
            ic = c_nchem2+1
            c_nchem2 = ic
            if ( c_nchem2 == c_cdim ) then
              write(c_out,901) c_nchem2
              write(6,901) c_nchem2
            end if
            c_tchem(c_nchem2) = tdum(iii)
          end if
          !
          ! ENTER SPECIES INTO AQUEOUS EQUILIBRIUM ARRAY
          !
          if ( ic > 0 ) then
            c_aqueous(nrq,iii) = ic
          else
            !
            ! TEST FOR POSSIBLE ERROR IN CHEMICAL NAMES
            !
            if ( tdum(iii) /= '        ' .and. &
                 tdum(iii) /= '       x' .and. &
                 tdum(iii) /= '       X' .and. &
                 tdum(iii) /= '      x ' .and. &
                 tdum(iii) /= '      X ' .and. &
                 tdum(iii) /= '    MORE')  then
              write(c_out, 89) tdum(iii), nrq,iii,(tdum(j),j=1,3)
              write(   6, 89) tdum(iii), nrq,iii,(tdum(j),j=1,3)
            end if
          end if
        end do
        !
        !  END NAMEREAD LOOP.  READ RATE PARAMETERS
        !
        read(c_rin,33) c_nrkq(nrq), c_rkq(1,nrq), c_rkq(2,nrq)
        if ( c_kkw == 5 ) write(c_out,33) c_nrkq(nrq), &
                       c_rkq(1,nrq), c_rkq(2,nrq)
      end do loopeq

      c_nreacq = nrq-1
      if ( c_kkw > 0 ) write(c_out,79) c_nreacq
      !
      ! ------------------------------------
      ! READ LOOP:  AQUEUS SPECIAL EQUILIBRIUM CONSTANTS WITHOUT H+/OH-
      ! ------------------------------------
      ! (note:  These are special.  A <->B+C; the combined species A are set
      !  in equilibria to be automatically solved when B is called.
      !  Species C MUST have much higher concentrations then either A or B,
      !  or else be zero.)
      !
      nrqq = 0
      do i = 1 , 1111
        read(c_rin,11) titl
        if ( titl == 'STAR' ) exit
      end do
      if ( c_kkw == 5 ) write(c_out,11) titl
      !
      ! LOOP FOR EQUILIBRIUM READS: CONTINUE UNTIL 'END'.
      !
      loopaqeq: &
      do i = 1 , c_cdim
        !
        ! ADVANCE COUNTER.
        !
        nrqq = nrqq+1
        if ( nrqq == 121 .or. nrqq ==  91 .or. nrqq == c_cdim ) then
          write(c_out,905) nrqq
        end if
        !
        ! --------------------------------
        ! READ NAMES OF SPECIAL EQUILIBRIUM SPECIES:  A = B + C* (*C=common spec
        ! --------------------------------
        !
        read(c_rin, 59) tdum(1), tdum(2), tdum(3)
        if ( c_kkw == 5 ) write(c_out, 59) tdum(1), tdum(2), tdum(3)
        !
        ! ENTER NAMES INTO AQUEOUS SPECIAL EQUILIBRIUM MATRIX
        ! 'END' BREAKS THE LOOP TO END THE EQUILIBRIUM READ.
        !
        do iii = 1 , 3
          if ( tdum(iii) == '     END' .or. &
               tdum(iii) == '    END ' ) exit loopaqeq
          c_treaqq(iii,nrqq) = tdum(iii)
          ic = namechem(tdum(iii))
          if ( ic > 0 ) then
            c_aqspec(nrqq,iii)=ic
          else
            !
            ! TEST FOR POSSIBLE ERROR IN CHEMICAL NAMES
            !
            if ( tdum(iii) /= '        ' .and. &
                 tdum(iii) /= '       x' .and. &
                 tdum(iii) /= '       X' .and. &
                 tdum(iii) /= '      x ' .and. &
                 tdum(iii) /= '      X ' .and. &
                 tdum(iii) /= '    MORE') then
              write(c_out, 87) tdum(iii), nrqq,iii,(tdum(j),j=1,3)
              write(   6, 87) tdum(iii), nrqq,iii,(tdum(j),j=1,3)
            end if
          end if
        end do
        !
        !  END NAMEREAD LOOP.  READ RATE PARAMETERS
        !   NOTE FUTURE DO OPTION: Make rate=rate of return reaction
        !     with nn corresponding to nn for regular reactions (   nrk)
        !
        read(c_rin,33) c_nrkqq(nrqq), c_rkqq(1,nrqq), c_rkqq(2,nrqq)
        if ( c_kkw == 5 ) write(c_out,33) c_nrkqq(nrqq), &
                   c_rkqq(1,nrqq), c_rkqq(2,nrqq)
      end do loopaqeq
      c_nreaqq = nrqq-1
      if ( c_kkw > 0 ) write(c_out,99) c_nreaqq
      !
      ! --------------------------------------
      ! READ LOOP:  REACTION RATES; REACTANTS, PRODUCTS AND STOICHIOMETRIES
      ! --------------------------------------
      !
      nr = 0
      do i = 1 , 1111
        read(c_rin,11) titl
        if ( titl == 'STAR' ) exit
      end do
      if ( c_kkw == 5 ) write(c_out,11) titl
      !
      ! LOOP FOR REACTION READS: CONTINUE UNTIL 'END'.
      !
      loopreaction: &
      do i = 1 , c_rdim
        !
        ! ADVANCE REACTION COUNTER
        !
        nr = nr+1
        if ( nr == 1711 .or. nr == 1911 .or. nr == 3211 ) then
          write(c_out,906) nr
        end if
        !
        ! --------------------------------------------------------
        ! READ NAMES OF REACTANTS, PRODUCTS AND STOICHIOMETRIES.
        ! --------------------------------------------------------
        !
        loopreactant: &
        do ii = 1 , 10
          read(c_rin,109) tdum(1), tdum(2), xdum(3), tdum(3), &
                          xdum(4), tdum(4), xdum(5), tdum(5)
          if ( c_kkw == 5 ) write(c_out,109) tdum(1), tdum(2), xdum(3), &
                  tdum(3), xdum(4), tdum(4), xdum(5), tdum(5)
          !
          ! ENTER NAMES INTO REACTANT OR PRODUCT MATRIX. (also TREAC)
          ! 'NEXT' BREAKS THE LOOP TO ADVANCE TO NEXT REACTION
          ! 'END' BREAKS THE LOOP TO END THE REACTION READ.
          !
          do iii = 1 , 5
            if ( ii == 1 .and. &
                 tdum(iii) /= '    NEXT' .and. &
                 tdum(iii) /= '     END' .and. &
                 tdum(iii) /= '   NEXT ' .and. &
                 tdum(iii) /= '    END ' ) then
              c_treac(iii,nr)=tdum(iii)
            end if

            if ( tdum(iii) == '    NEXT' .or. &
                 tdum(iii) == '   NEXT ' ) exit loopreactant
            if ( tdum(iii) == '     END' .or. &
                 tdum(iii) == '    END ') exit loopreaction
            ic = namechem(tdum(iii))
            if ( ic > 0 ) then
              if ( iii <= 2 ) then
                c_reactant(nr,iii) = ic
              else
                do j = 1 , 20
                  if ( c_product(nr,j) == 0 ) exit
                end do
                c_product(nr,j) = ic
                c_stoich(nr,j) = xdum(iii)
                c_nnpro(nr) = j
                c_prodarr(nr,ic) = xdum(iii)
                !
                ! TEST WRITE PRODARR
                ! (=stoichiometry for production of species ic from reaction nr)
                !
                if ( c_kkw > 0 ) write(c_out,127) nr, c_rk(1,nr), tdum
                if ( c_kkw > 0 ) write(c_out,128) &
                               j, ic, c_stoich(nr,j), c_prodarr(nr,ic)
              end if
            else
              !
              ! FLAG FOR HV REACTION:  Identified by 'hv' as Reactiant #2.
              !   Set REACTANT(2) =-1,
              !    and default hv rate index 4 (rate proportional to jNO2)
              !
              if ( iii == 2 .and. &
                  (tdum(iii) == '      hv' .or. &
                   tdum(iii) == '     hv ') ) then
                c_reactant(nr,iii) = -1
                if ( c_nrk(nr) == 0 ) c_nrk(nr) = 4
              end if
              !
              ! TEST FOR POSSIBLE ERROR IN CHEMICAL NAMES
              !
              if ( tdum(iii) /= '      hv' .and. &
                   tdum(iii) /= '        ' .and. &
                   tdum(iii) /= '       x' .and. &
                   tdum(iii) /= '       X' .and. &
                   tdum(iii) /= '     hv ' .and. &
                   tdum(iii) /= '      x ' .and. &
                   tdum(iii) /= '    MORE' ) then
                write(c_out,139) tdum(iii), nr,iii,tdum
                write(   6,139) tdum(iii), nr,iii,tdum
              end if
            end if
          end do
        end do loopreactant
        !
        ! -------------
        ! END READ REACTION NAMES
        ! -------------
        !
        ! --------------------------------------------------------
        !  READ RATE PARAMETERS (moved after NAMES)
        ! --------------------------------------------------------
        read(c_rin,131) ii, c_nrk(nr)
        if ( c_kkw == 5 ) write(c_out,131) ii, c_nrk(nr)
        !
        ! READ OPTIONS AND FLAG FOR READ
        !
        if ( c_nrk(nr) == 0 .or. c_nrk(nr) == -1 .or. c_nrk(nr) == -20 ) then
          read(c_rin,132) c_rk(1,nr), c_rk(2,nr)
          if ( c_kkw == 5 ) write(c_out,132) c_rk(1,nr), c_rk(2,nr)
        end if
        if ( c_nrk(nr) == -2 .or. c_nrk(nr) == -5 .or. c_nrk(nr) == -6) then
          read(c_rin,132) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)
          if ( c_kkw == 5 ) write(c_out,132) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)
        end if
        if ( c_nrk(nr) == -3 ) then
          read(c_rin,133) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr), &
                          c_rk(4,nr), c_rk(5,nr)
          if ( c_kkw == 5 ) write(c_out,133) c_rk(1,nr), c_rk(2,nr), &
                          c_rk(3,nr), c_rk(4,nr), c_rk(5,nr)
        end if
        if ( c_nrk(nr) == -4 ) then
          read(c_rin,134) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr), &
                          c_rk(4,nr), c_rk(5,nr),c_rk(6,nr), c_rk(7,nr)
          if ( c_kkw == 5 ) write(c_out,134) c_rk(1,nr), c_rk(2,nr), &
                   c_rk(3,nr), c_rk(4,nr), c_rk(5,nr),c_rk(6,nr), c_rk(7,nr)
        end if
        if ( c_nrk(nr) == -7 ) then
          read(c_rin,135) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr), &
                          c_rk(4,nr), c_rk(5,nr), c_rk(6,nr)
          if ( c_kkw == 5 ) write(c_out,135) c_rk(1,nr), c_rk(2,nr), &
               c_rk(3,nr), c_rk(4,nr), c_rk(5,nr), c_rk(6,nr)
        end if
        if ( c_nrk(nr) == -8 ) then
          read(c_rin,136) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr), &
                          c_rk(4,nr), c_rk(5,nr)
          if ( c_kkw == 5 ) write(c_out,136) c_rk(1,nr), c_rk(2,nr), &
                          c_rk(3,nr), c_rk(4,nr), c_rk(5,nr)
        end if
        if ( c_nrk(nr) == -9 .or. c_nrk(nr) == -10 ) then
          read(c_rin,136) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr), &
                          c_rk(4,nr), c_rk(5,nr), c_rk(6,nr)
          if ( c_kkw == 5 ) write(c_out,136) c_rk(1,nr), c_rk(2,nr), &
                          c_rk(3,nr), c_rk(4,nr), c_rk(5,nr), c_rk(6,nr)
        end if
        if ( c_nrk(nr) == -11 ) then
          read(c_rin,137) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr), &
                          c_rk(4,nr), c_rk(5,nr)
          if ( c_kkw == 5 ) write(c_out,137) c_rk(1,nr), c_rk(2,nr), &
                          c_rk(3,nr), c_rk(4,nr), c_rk(5,nr)
        end if
        if ( c_nrk(nr) == -12 .or. c_nrk(nr) == -13 ) then
          read(c_rin,136) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr), c_rk(4,nr)
          if ( c_kkw == 5 ) write(c_out,136) c_rk(1,nr), c_rk(2,nr), &
                         c_rk(3,nr), c_rk(4,nr), c_rk(5,nr), c_rk(6,nr)
        end if
        !
        !  Addition for secondary organic aerosols (Luis Olcese)
        !
        if ( c_nrk(nr) == -21 ) then
          read(c_rin,142) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr), c_rk(4,nr)
          if ( c_kkw == 5 ) write(c_out,142) c_rk(1,nr), c_rk(2,nr), &
                         c_rk(3,nr), c_rk(4,nr)
        end if
        if ( c_nrk(nr) == -22 ) then
          read(c_rin,142) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr), c_rk(4,nr)
          if ( c_kkw == 5 ) write(c_out,142) c_rk(1,nr), c_rk(2,nr), &
                         c_rk(3,nr), c_rk(4,nr)
        end if
        !
        ! READ FOR hv REACTIONS
        !
        if ( c_nrk(nr) > 0 .or. c_nrk(nr) <= -30 ) then
          read(c_rin,132) c_rk(1,nr)
          if ( c_kkw == 5 ) write(c_out,132) c_rk(1,nr)
        end if
        !
        ! COUNT NUMBER OF PARAMETERIZED RO2 REACTIONS AND RECORD REACTION NUMBER
        !
        if ( c_nrk(nr) == -13 ) then
          c_nnrro2 = c_nnrro2 + 1
          c_nrro2(c_nnrro2) = nr
        end if
      end do loopreaction
      !
      ! --------------------------------------------------------
      ! END REACTION LOOP.
      ! --------------------------------------------------------
      !
      c_nreac = nr -1
      !
      ! -------------------------------------------------------
      ! READ SPECIAL FAMILY SPECIES (CURRENTLY HX AND NX, IN ORDER)
      ! -------------------------------------------------------
      ! DELETED - REPLACED WITH SPECIAL SPECIES INDICES
      !           AND AUTOMATIC ID FROM SPECIES NAMES
      !
      ! -------------------------------------------------------
      ! READ CASCADE SPECIES LIST.  (INCLUDES CASCADE PAIRS)
      ! -------------------------------------------------------
      !   FUTURE DO:
      !     CHANGE CASCADE READ to be c_cascade(i,1), n;
      !      then if n>1, read c_cascade(i,ii),i=2,n
      !      (maybe add flag for DIFFICULT CONVERGENCE)
      !
      !
      !  THE CASCADE:
      !   IT CONSISTS OF A SEQUENCE OF SPECIES, GIVING THE ORDER OF SOLUTION.
      !   THE CASCADE INCLUDES THREE TYPES OF SOLUTIONS
      !
      !    (1) A SINGLE SPECIES ON A LINE, to be solved by itself.
      !    (2) A PAIR OF INTERACTING SPECIES, or A SERIES OF LINKED PAIRS
      !         These are solved simultaneously by an algorithm for linked pai
      !          and can include:  A<->B1, A<->B2, B1<->C1, etc.
      !    (3) A GROUP OF THREE OR INTERACTING MORE SPECIES,
      !          to be solved by the MULTISOLVE algorithm.
      !
      !     SINGLE SPECIES are entered on a single line:  (SPEC1  x   x)
      !
      !     CASCADE PAIRS are entered either as a PAIR on a single line
      !                                                  (SPEC1 SPEC2 x)
      !        This tells the cascade to solve for the species now.
      !        or with a PAIR declaration                (pair SPEC1 SPEC2)
      !        This tells the cascade to solve SPEC2 elsewhere,
      !           whenever SPEC1 is called
      !       (which can be a single species, part of a pair, or multisolve.)
      !
      !     MULTISOLVE SPECIES are entered as groups of three,
      !       continuing until a blank is entered (   x).
      !
      !  Notes on CASCADE PAIRS:
      !   A primary species may be paired with many secondary species.
      !   The secondary species may then be paired with its own subspecies.
      !   Typically there are rapid reactions between directly paired species.
      !   It is OK for the paired species to have no interactions
      !      (e.g. species that interact only in aqueous phase)
      !
      !   The chain of paired  species is built
      !     with parent primary species first, followed by its (multiple)
      !      secondary species each with its  own  subspecies
      !
      !   The pair chain uses  pointers, which are set below:
      !  c_nppair(ic,1) = ics:  pointer from species ic
      !       to its directly linked primary pair species ics (or self)
      !  c_nppair(ic,2) = ics:  pointer from species ic
      !    to ultimate primary pair species ics (or self)
      !  c_nppair(ics,3) = number of chemical pairs associated w/ species ics.
      !  c_nppair(ics,np),np>=4:  identifies nth pair species or subspecies
      !    for species ics, for total number given by nppair(ics,3).
      !
      !    (old format: pointers as with nequil:
      !      npair(ics) = # of chem pairs  associated with the species
      !      nppair(ic) = ics:  pointer to primary  pair species (or self)
      !     ncpair(ics,np) = ic:  Identifies nth  pair species orsubspec.
      !    )
      !
      ! BEGIN CASCADE READ
      !
      do i = 1 ,1111
        read(c_rin,11) titl
        if(titl == 'STAR') exit
      end do
      !
      ! CASCADE READ LOOP
      !
      ncas = 1
      ncasp = 0
      nc = 0
      loopcascade: &
      do i = 1 , c_cdim
        read(c_rin,309)(tdum(ii),ii=1,3)
        !
        ! LOOP TO READ THREE SPECIES AND CONVERT TO SPECIES NUMBERS
        !   NOTE:  This read loop repeats to read full group
        !          until it reads a break mark (x)
        ! LOOP TO READ GROUP OF THREE SPECIES (320)
        !
        do ii = 1 , 3
          !
          ! IF NAME IS 'END', BREAK THE CASCADE READ AND EXIT
          !     (ALSO ADJUST COUNTER)
          !
          if ( tdum(1) == '     END' .or. &
               tdum(1) == '    END ' ) then
            if ( nc  == 0 ) ncas = ncas - 1
            exit loopcascade
          end if
          !
          ! IF NAME IS  x, BREAK, PROCESS, AND GO TO THE NEXT CASCADE  GROUP.
          !
          if ( tdum(ii) == '       x' .or. &
               tdum(ii) == '        ' .or. &
               tdum(ii) == '     pre' .or. &
               tdum(ii) == '    post' ) exit
          !
          ! IDENTIFY SPECIES NAME
          !
          ic = namechem(tdum(ii))
          !
          ! IF FIRST NAME IS 'pair' , SET FLAG TO ESTABLISH CASCADE PAIRS
          !
          if ( ii == 1 .and. &
              (tdum(ii) == '    pair' .or. &
               tdum(ii) == '    PAIR' .or. &
               tdum(ii) == '   pair ' .or. &
               tdum(ii) == '   PAIR ') ) ic = -1
          !
          ! WARNING FOR NO-NAME IN CASCADE. (no-name is regarded as break)
          !
          if ( ic == 0 ) then
            write(c_out,319) (tdum(j),j=1,3),ii
            write(6,319) (tdum(j),j=1,3),ii
            exit
          end if
          !
          ! ENTER SPECIES IN CASCADE LIST (AND ADVANCE COUNTER)
          !  (NOTE: LATER ADJUSTMENT FOR CASCADE PAIR FLAG)
          !
          nc = nc + 1
          c_cascade(ncas ,nc ) = ic
        end do
        !
        ! PROCESS CASCADE SPECIES GROUP  AFTER READ OF COMPLETE GROUP.
        !
        ! ESTABLISH CASCADE PAIRS BASED ON READ GROUP OF THREE SPECIES.
        !
        ! (1) IF CASCADE LIST CONSISTS OF JUST TWO SPECIES,
        !     THEN IDENTIFY SECOND SPECIES AS PAIRED WITH THE FIRST
        !     AND REMOVE THE SECOND SPECIES  FROM THE CASCADE ARRAY
        !    (paired species will automatically be solved when primary is called
        !
        if ( c_cascade(ncas,1) > 0 .and. &
             c_cascade(ncas,2) > 0 .and. &
             c_cascade(ncas,3) == 0) then
          ncasp = ncasp + 1
          caspair(ncasp,1) = c_cascade(ncas,1)
          caspair(ncasp,2) = c_cascade(ncas,2)
          c_cascade(ncas,2) = 0
        end if
        !
        ! (2) IF CASCADE LIST HAS 'PAIR' FLAG (first species is -1)
        !     THEN IDENTIFY 2nd AND 3rd SPECIES AS PAIR
        !     AND REMOVE BOTH FROM CASCADE ARRAY.
        !      (2nd species should be in cascade or in cascade pair elsewhere)
        !    (with 1st species set to zero, this cascade number will be re-done)
        !
        if ( c_cascade(ncas,1) == -1 .and. &
             c_cascade(ncas,2) > 0 .and. &
             c_cascade(ncas,3) > 0) then
          ncasp = ncasp + 1
          caspair(ncasp,1) = c_cascade(ncas,2)
          caspair(ncasp,2) = c_cascade(ncas,3)
          c_cascade(ncas,1) = 0
          c_cascade(ncas,2) = 0
          c_cascade(ncas,3) = 0
        end if
        !
        !  ADVANCE THE SPECIES COUNTER UNLESS 1st SPECIES IS ZERO
        !    (note: do not advance if the line was just to identify a pair;
        !      or if it is a mistake.)
        !
        if ( c_cascade(ncas,1) > 0 ) ncas = ncas + 1
        nc = 0
      end do loopcascade
      !
      ! -------------------------------------------------------
      ! END CASCADE READ
      ! -------------------------------------------------------
      !
      ! --------------------------
      ! IDENTIFY SPECIAL SPECIES NUMBERS AND CATEGORIES
      ! --------------------------
      !  This replaces 'family'.
      ! LOOP TO IDENTIFY SPECIAL SPECIES
      !
      do ic = 1 , c_nchem2
        if ( c_tchem(ic) == '     H2O' .or. &
             c_tchem(ic) == '    H2O ' ) then
          c_nh2o =ic
          c_icat(ic) = 0
        end if
        if ( c_tchem(ic) == '     CO2' .or. &
             c_tchem(ic) == '    CO2 ' ) then
          c_nco2 =ic
          c_icat(ic) =20
        end if
        if ( c_tchem(ic) == '      CO' .or. &
             c_tchem(ic) == '     CO ' ) then
          c_nco  =ic
        end if
        if ( c_tchem(ic) == '     CH4' .or. &
             c_tchem(ic) == '    CH4 ' ) then
          c_nch4 =ic
        end if
        if ( c_tchem(ic) == '     CO3' .or. &
             c_tchem(ic) == '    CO3 ' ) then
          c_nco3 =ic
          c_icat(ic) = 8
        end if
        if ( c_tchem(ic) == '      H+' .or. &
             c_tchem(ic) == '     H+ ' ) then
          c_nhplus  =ic
          c_icat(ic) =17
        end if
        if ( c_tchem(ic) == '     OH-' .or. &
             c_tchem(ic) == '    OH- ' ) then
          c_nohmin =ic
          c_icat(ic) =18
        end if
        if ( c_tchem(ic) == '      OH' .or. &
             c_tchem(ic) == '     OH ' ) then
          c_noh  =ic
          c_icat(ic) = 9
        end if
        if ( c_tchem(ic) == '     HO2' .or. &
             c_tchem(ic) == '    HO2 ' ) then
          c_nho2 =ic
          c_icat(ic) =10
        end if
        if ( c_tchem(ic) == '    H2O2' .or. &
             c_tchem(ic) == '   H2O2 ' ) then
          c_nh2o2=ic
          c_icat(ic) =19
        end if
        if ( c_tchem(ic) == '      O3' .or. &
             c_tchem(ic) == '     O3 ' ) then
          c_no3  =ic
          c_icat(ic) =11
        end if
        if ( c_tchem(ic) == '     NO2' .or. &
             c_tchem(ic) == '    NO2 ' ) then
          c_nno2 =ic
          c_icat(ic) =12
        end if
        if ( c_tchem(ic) == '      NO' .or. &
             c_tchem(ic) == '     NO ' ) then
          c_nno  =ic
          c_icat(ic) =13
        end if
        if ( c_tchem(ic) == '     NO3' .or. &
             c_tchem(ic) == '    NO3 ' ) then
          c_nno3 =ic
          c_icat(ic) =14
        end if
        if ( c_tchem(ic) == '    N2O5' .or. &
             c_tchem(ic) == '   N2O5 ' ) then
          c_nn2o5=ic
          c_icat(ic) =15
        end if
        if ( c_tchem(ic) == '    HNO3' .or. &
             c_tchem(ic) == '   HNO3 ' ) then
          c_nhno3=ic
          c_icat(ic) =16
        end if
      end do
      !
      ! END  LOOP TO IDENTIFY SPECIAL SPECIES
      !
      ! WARN IF THESE SPECIES ARE MISSING
      !
      if ( c_nh2o   == 0 ) write(c_out, 1381)
      if ( c_nco2   == 0 ) write(c_out, 1382)
      if ( c_nhplus == 0 ) write(c_out, 1383)
      if ( c_nohmin == 0 ) write(c_out, 1384)
      if ( c_noh    == 0 ) write(c_out, 1385)
      if ( c_nho2   == 0 ) write(c_out, 1386)
      if ( c_nh2o2  == 0 ) write(c_out, 1387)
      if ( c_no3    == 0 ) write(c_out, 1388)
      if ( c_nno2   == 0 ) write(c_out, 1389)
      if ( c_nno    == 0 ) write(c_out, 1390)
      if ( c_nno3   == 0 ) write(c_out, 1391)
      if ( c_nn2o5  == 0 ) write(c_out, 1392)
      if ( c_nhno3  == 0 ) write(c_out, 1393)
      !
      !
      ! ----------------------
      ! ASSIGN AQUEOUS SPECIES TO ASSOCIATED GAS-MASTER SPECIES
      ! IN SEQUENCE:  GAS->AQUEOUS->AQUEOUS ION->AQUEOUS DOUBLE ION
      ! IDENTIFIED FROM HENRY'S LAW AND AQUEOUS EQUILIBRIUM SEQUENCE
      ! ----------------------
      !
      ! GAS-AQUEOUS LINK INCLUDES:
      !  c_npequil = pointer to gas-master species for aqueous link
      !  c_nequil = number of aqueous species associated with gas-master.
      !  c_ncequil = ic for each aqueous species assoc. with gas-master.
      !   (up to 3 species in sequence: aqueous, ion1, ion2)
      !
      !   (FUTURE DO ALTERNATIVE:  combine.
      !      c_npequil(ic,1) = pointer from aqueous to gas-master speciees
      !      c_npequil(ic,2) = number of aqueous spec. assoc with gas-master
      !      c_npequil(ic,n),n=3,4,5:  ic for each aq species assoc w/gas
      !   (
      !  nrequil = nrh, nrq for reaction number associated with
      !            corresponding aqueous species in ncequil
      !
      ! ------------------------------------------------------------
      !
      ! PRELIMINARY ZERO
      do ic = 1 , c_nchem2
        c_nequil(ic) = 0
        c_npequil(ic) = ic
        c_ion(ic) = 0
        do i = 1 , 6
          c_ncequil(ic,i) = 0
        end do
      end do
      !
      ! SET ION INDEX FOR H+, OH-
      ! H+
      !
      if ( c_aqueous(1,2) > 0 ) then
        c_ion(c_aqueous(1,2))= 1
      end if
      !
      ! OH-
      !
      if ( c_aqueous(1,3) > 0 ) then
        c_ion(c_aqueous(1,3))=-1
      end if
      !
      ! LOOP FOR HENRY'S LAW COEFFICIENTS
      !
      do nrh = 1 , c_nreach
        if ( c_henry(nrh,1) == 0 ) exit
        !
        ! FOR EACH HENRY'S LAW COEFFICIENT, ASSOCIATE AQUEOUS WITH GAS SPECIES.
        !
        ic = c_henry(nrh,1)
        ic1 = c_henry(nrh,2)
        c_npequil(ic1) = ic
        c_icat(ic1) = c_icat(ic)
        c_lsts(ic1) = c_lsts(ic)
        c_nequil(ic) = c_nequil(ic) + 1
        nne = c_nequil(ic)
        c_ncequil(ic,nne) = ic1
        c_nrequil(ic,nne) = nrh
        !
        ! CHANGE PRODUCT ARRAY TO GAS-MASTER
        ! (add to avoid replacing already-set array)
        do nr = 1 , c_nreac
          c_prodarr(nr,ic) = c_prodarr(nr,ic) + c_prodarr(nr,ic1)
          c_prodarr(nr,ic1) = d_zero
        end do
        !
        ! GO THROUGH AQUEOUS EQUILIBRIA.  ADD EQUILIBRIUM SPECIES
        ! TO GAS-MASTER SPECIES LIST.
        ! NOTE EQUILIBRIA MUST COME IN SEQUENCE:
        !             GAS->AQUEOUS->FIRST ION->DOUBLE ION.
        ! ALSO SET ION INDEX.
        !
        loopreaq: &
        do nrq = 1 , c_nreacq
          if ( c_aqueous(nrq,1) == 0 .and. &
               c_aqueous(nrq,2) == 0 ) exit loopreaq
          !
          !  H+ SPECIAL here: for H+, OH-  skip equilibrium index assignment
          !
          if ( c_aqueous(nrq,1) == 0 ) cycle
          if ( c_npequil(c_aqueous(nrq,1)) == ic ) then
            ic1 = c_aqueous (nrq,2)
            c_npequil(ic1) = ic
            c_icat(ic1) = c_icat(ic)
            c_lsts(ic1) = c_lsts(ic)
            c_nequil(ic) = c_nequil(ic) + 1
            nne = c_nequil(ic)
            c_ncequil(ic,nne) = ic1
            c_nrequil(ic,nne) = nrq
            c_ion(ic1) = c_ion(c_aqueous(nrq,1))
            if ( c_aqueous(nrq,3) > 0 ) then
              c_ion(ic1) = c_ion(ic1) - c_ion(c_aqueous(nrq,3))
              if ( c_kkw > 0 ) write(c_out,659) &
                (c_treacq(j,nrq),j=1,3), c_ion(ic1), c_ion(c_aqueous(nrq,3))
            end if
            !
            ! CHANGE PRODUCT ARRAY TO GAS-MASTER.
            !                (add to avoid replacing already-set array)
            do nr = 1 , c_nreac
              c_prodarr(nr,ic) = c_prodarr(nr,ic) + c_prodarr(nr,ic1)
              c_prodarr(nr,ic1) = d_zero
            end do
          end if
        end do loopreaq
      end do
      !
      ! SUMMARY WRITE: HENRY'S LAW AND AQUEOUS EQUILIBRIUM SPECIES (ncf)
      !
      write(c_out,1321)
      write(c_out,1322) c_tchem(c_aqueous(1,2)),c_tchem(c_aqueous(1,3))
      do nrh = 1 , c_nreach
        ic = c_henry(nrh,1)
        if ( c_nequil(ic) > 0 ) then
          write(c_out, 1322) c_tchem(ic), &
            (c_tchem(c_ncequil(ic,i)),i=1,c_nequil(ic))
        end if
      end do
      !
      ! -----------------------------------------
      !  END ASSIGN-AQUEOUS
      ! -----------------------------------------
      !
      ! -----------------------------------------
      !  CONVERT SPECIAL EQUILIBRIUM SPECIES TO BACK-FORTH REACTIONS
      ! -----------------------------------------
      !  Special equilibrium forward reaction 1e8 sec-1,
      !   backward reaction based on equilibrium coefficient.
      !
      !  (Special equilibrium reactions NOT entered into cascade pair.)
      !
      ! LOOP FOR SPECIAL AQUEOUS EQUILIBRIUM SPECIES
      !
      nr = c_nreac
      do nrqq = 1 , c_nreaqq
        if ( c_aqspec(nrqq,1) /= 0 .and. c_aqspec(nrqq,2) /= 0 ) then
          icc1 = c_aqspec(nrqq,1)
          icc2 = c_aqspec(nrqq,2)
          icc3 = c_aqspec(nrqq,3)
          ic1 = c_npequil(icc1)
          ic2 = c_npequil(icc2)
          ic3 = c_npequil(icc3)
          !
          ! ENTER INTO CASCADE PAIR - CUT.
          ! Note:  1st SPECIAL EQUIL species is secondary paired species
          !        2nd species is primary paired species  (A <=> B + C)
          !        3rd species is totally independent
          ! If this is included, also change 'ncasp =' above.
          !         caspair(nrqq,  1 ) = icc2
          !         caspair(nrqq,  2 ) = icc1
          !
          ! ENTER FORWARD REACTION
          !
          nr = nr+1
          if ( nr == 1711 .or. nr == 1911 .or. nr == 3211 ) then
            write(c_out,906) nr
          end if
          c_treac(1,nr) = c_tchem(icc1)
          c_treac(3,nr) = c_tchem(icc2)
          c_reactant(nr,1) = icc1
          c_product(nr,1) = icc2
          c_stoich(nr,1) = d_one
          c_nnpro(nr) = 1
          c_prodarr(nr,ic2 ) = d_one
          if ( icc3 > 0 ) then
            c_treac(4,nr) = c_tchem(icc3)
            c_product(nr,2) = icc3
            c_stoich(nr,2) = d_one
            c_nnpro(nr) = 2
            c_prodarr(nr,ic3 ) = d_one
          end if
          !
          ! FORWARD RATE CONSTANT:  1e8 s-1.
          ! OPTION - maybe  this causes cncn-> 0.  alt 1e2 s-1
          !
!         c_rk(1,nr) = 1.0D+08
          c_rk(1,nr) = 1.0D+02
          !
          ! ENTER BACKWARD REACTION
          !
          nr = nr+1
          if ( nr == 711 .or. nr == 911 .or. nr == 1211) then
            write(c_out,906) nr
          end if
          c_treac(1,nr) = c_tchem(icc2)
          c_treac(3,nr) = c_tchem(icc1)
          c_reactant(nr,1) = icc2
          c_product(nr,1) = icc1
          c_stoich(nr,1) = d_one
          c_nnpro(nr) = 1
          c_prodarr(nr,ic1) = d_one
          if ( icc3 > 0 ) then
            c_treac(2,nr) = c_tchem(icc3)
            c_reactant(nr,2) = icc3
          end if
          !
          ! BACKWARD RATE CONSTANT
          !  Equil  constant in  moles/liter A exp(-B * (1/temp - 1/298) )
          !  A<->B+C  KA = BC  K in moles/liter
          !  k1a = k2BC, k1/k2=K, k2=k1/K, k2 as A2 exp(-B2/temp) also moles/li
          !  k1=1e8 s-1.  A2=k1/(Aexp(+B/298),  B2=-B
          c_rk(1,nr) = 1.0D+08/(c_rkqq(1,nrqq)*dexp(c_rkqq(2,nrqq)/298.0D0))
          c_rk(2,nr) = d_zero - c_rkqq(2,nrqq)
        end if
      end do
      !
      ! END LOOP FOR SPECIAL EQUILIBRIUM SPECIES
      !
      c_nreac = nr
      !
      ! -----------------------------------------
      !  END CONVERT SPECIAL EQUILIBRIUM SPECIES TO BACK-FORTH REACTIONS
      ! -----------------------------------------
      !
      ! -----------------------------------------
      !  ADD NO3-N2O5 (hard-wired) INTO CASCADE PAIR LIST
      ! -----------------------------------------
      ncasp = ncasp + 1
      caspair(ncasp,1) = c_nno3
      caspair(ncasp,2) = c_nn2o5
      !
      ! -----------------------------------------
      ! GENERATE POINTERS FOR CASCADE PAIRS
      ! -----------------------------------------
      !  Notes:
      !   A primary species may be paired with many secondary pair species.
      !   The secondary pair species may then have its own pair subspecies.
      !   Typically there are rapid reactions between directly linked pairs.
      !   It is OK for the paired species to have no interactions
      !      (e.g. species that interact only in aqueous phase)
      !
      !   The chain of paired  species is built
      !     with parent primary species first,
      !       followed by its secondary pair  species,
      !          each with its  own pair subspecies.
      !
      !   Pointers all point from and to gas-master species
      !
      !   The pair chain pointers are:
      !  c_nppair(ic,1) = ics:  pointer from species ic
      !       to its directly linked primary pair species ics (or self)
      !  c_nppair(ic,2) = ics:  pointer from species ic
      !    to ultimate primary pair species ics (or self)
      !  c_nppair(ics,3) = number of chemical pairs associated w/ species ics.
      !  c_nppair(ics,np),np>=4:  identifies nth pair species or subspecies
      !    for species ics, for total number given by nppair(ics,3).
      !
      !    (old format: pointers as with nequil:
      !      npair(ics) = # of chem pairs  associated with the species
      !      nppair(ic) = ics:  pointer to primary  pair species (or self)
      !     ncpair(ics,np) = ic:  Identifies nth  pair species orsubspec.
      !    )
      !
      !
      !      GENERATE PAIR POINTERS
      !    - if c_nppair(2nd spec) does not equal self, error
      !    - make 2nd species and its subspecies point to primary of 1st spec.
      !    - if 2nd has subspecies, make point to 1st species
      !    - if 1st species points elsewhere, make 2nd and its subspecies
      !       also point there.
      !
      ! LOOP THROUGH CASCADE PAIRS
      !
      do i = 1 , ncasp
        icc1 = caspair(i,1)
        icc2 = caspair(i,2)
        ic1 = 0
        ic2 = 0
        if ( icc1 > 0 ) then
          ic1 = c_npequil(icc1)
        end if
        if ( icc2 > 0 ) then
          ic2 = c_npequil(icc2)
        end if
        ic3 = 0
        if ( ic2 > 0 ) then
          ic3 = c_nppair(ic2,1)
        end if
        !
        ! ERROR CHECK:  IF 2ND SPECIES IS ALREADY PAIRED, MAJOR ERROR
        !
        if ( ic3 /= ic2 .or. icc1 == 0 .or. icc2 == 0 .or. &
            ic1 == 0 .or. ic2 == 0 )  then
          write(c_out,667) i, icc1, icc2, ic1, ic2, c_nppair(ic2,1)
          write(6,667) i, icc1, icc2, ic1, ic2, c_nppair(ic2,1)
          if ( icc1 > 0 ) write (c_out,309) c_tchem(icc1)
          if ( icc2 > 0 ) write (c_out,309) c_tchem(icc2)
          if ( ic1 > 0 )  write (c_out,309) c_tchem(ic1)
          if ( ic2 > 0 )  write (c_out,309) c_tchem(ic2)
          if ( icc1 > 0 ) write (c_out,309) c_tchem(icc1)
          if ( icc1 > 0 ) write (6,309) c_tchem(icc1)
          if ( icc2 > 0 ) write (6,309) c_tchem(icc2)
          if ( ic1 > 0 )  write (6,309) c_tchem(ic1)
          if ( ic2 > 0 )  write (6,309) c_tchem(ic2)
          if ( icc1 > 0 ) write (6,309) c_tchem(icc1)
        else
          !
          ! Assign all species and subspecies to the primary species
          !               of the first of the pair.
          !  Note:  npair(ic2,3) = number of subspecies;
          !         npair(ic2,4) = 1st subsp, etc.
          !  Also assign direct pair species - just for this secondary species.
          !
          c_nppair(ic2,1) = ic1
          ic3 = c_nppair(ic1,2)
          do ii = 1 , (c_nppair(ic2,3)+1)
            ic = ic2
            if ( ii > 1 ) ic = c_nppair(ic2,(ii+2))
            c_nppair(ic ,2) = ic3
            c_nppair(ic3,3) = c_nppair(ic3,3) + 1
            np = c_nppair(ic3,3) + 3
            c_nppair(ic3,np) = ic
            !
            ! Assign CATEGORY and STEADY STATE INDEX to secondary species
            ! also assign to all its subspecies, including aqueous
            ! equil. species.
            !
!           do neq = 1 , (c_nequil(ic )+1)
!             icc = ic
!             if ( neq > 1 ) icc = c_ncequil(ic, (neq-1) )
!             c_nppair(ic ,1) = ic3
!             c_icat(icc) = c_icat(ic3)
!             c_lsts(icc) = c_lsts(ic3)
!           end do
          end do
          !
          ! End Assign loop
          !
        end if
        !
        ! End error check if
        !
      end do
      !
      ! END LOOP THROUGH CASCADE PAIRS
      !
      ! -----------------------------------------
      ! END GENERATE POINTERS FOR CASCADE PAIRS
      ! -----------------------------------------
      !
      ! -------------------------------------------------------
      ! RUN CASCADE - IDENTIFY MISSING SPECIES.
      !   ALSO IDENTIFY MAX. NUMBER OF SPECIES SOLVED TOGETHER (NCDIM)
      !   AND MAXIMUM MULTISOLVE GROUP (NSDIM)
      ! TEST CATEGORY LIST.  VARIABLES W/O CATEGORIES ARE ASSIGNED 'xx'
      !  AND A WARNING ISSUED
      ! -------------------------------------------------------
      !
      !  Set LCASTEST = false initially
      !   Set LCASTEST = true for H+, OH-, CO2, H2O
      !                   - these should never be in cascade
      !
      do ic = 1 , c_nchem2
        lcastest(ic) = .false.
        if ( ic == c_aqueous(1,2) .or. ic == c_aqueous(1,3)) then
          lcastest(ic) = .true.
        end if
        if ( c_tchem(c_npequil(ic )) == '     CO2' ) lcastest(ic) = .true.
        if ( c_tchem(c_npequil(ic )) == '     H2O' ) lcastest(ic) = .true.
      end do
      !
      !   Set LCASTEST = true for LUMPED SPECIES - these are never included.
      !
      do i = 1 , nlump
        icc = c_lump(i,1)
        lcastest(icc) = .true.
      end do
      !
      ! SET COUNTERS (NSOL, NSOLV) AND MAXIMA (NCDIM, NSDIM)
      !
      ncdim = 1
      nsdim = 1
      !
      ! TEST WRITE CASCADE
      !
      if ( c_kkw > 0 ) write(c_out,347)
      !
      !  RUN  CASCADE.  IDENTIFY INCLUDED SPECIES.
      !   COUNT SPECIES IN CASCADE GROUP.   ALSO ASSIGN NMULTI
      !
      do i = 1 , c_nchem2
        if ( c_cascade(i,1) == 0 ) exit
        nsolv = 0
        nsol = 0
        loopcount: &
        do ii = 1 , c_nchem2
          ics = c_cascade(i,ii)
          if ( ics == 0 ) exit loopcount
          ics = c_nppair(c_npequil(ics ),2)
          if ( ii == 1 ) ic1 = ics
          nsol = nsol + 1
          nsolv = nsolv + 1 + c_nppair(ics,3)
          do n = 1 , (c_nppair(ics,3)+1)
            icc = ics
            if ( n > 1 ) then
              icc = c_nppair(ics,n+2)
            end if
            !
            ! TEST WRITE CASCADE
            !
            if ( c_kkw >= 1 .and. icc > 0 ) then
              write(c_out,348) n, c_tchem(icc)
            end if
            do neq = 1 , (c_nequil(icc)+1)
              ic = icc
              if ( neq > 1 ) ic = c_ncequil(icc,(neq-1))
              lcastest(ic) = .true.
              c_npmulti(ic,1) = ic1
            end do
          end do
        end do loopcount
        if ( ncdim < nsolv ) ncdim=nsolv
        if ( nsdim < nsol ) nsdim=nsol
      end do
      !
      !  INCLUDE FAMILY SPECIES, SAME AS FOR CASCADE.
      !
      do i = 1 , 7
        nsolv = 0
        nsol = 0
        do ii = 1 , 3
          ! ADDED 2009
          ics = 0
          if ( i == 1 .and. ii == 1 ) ics = c_nhno3
          if ( i == 2 .and. ii == 1 ) ics = c_nno3
!         if ( i == 3 .and. ii == 1 ) ics = c_nn2o5
          if ( i == 4 .and. ii == 1 ) ics = c_no3
          if ( i == 5 .and. ii == 1 ) ics = c_nno2
          if ( i == 6 .and. ii == 1 ) ics = c_nno
          if ( i == 7 .and. ii == 1 ) ics = c_noh
          if ( i == 7 .and. ii == 2 ) ics = c_nho2
          if ( i == 7 .and. ii == 3 ) ics = c_nh2o2
          if ( ics > 0 ) then
            ics = c_nppair(c_npequil(ics ) , 2)
            nsol = nsol + 1
            nsolv = nsolv + 1 + c_nppair(ics,3)
            do n = 1 , (c_nppair(ics,3)+1)
              icc = ics
              if ( n > 1 ) then
                icc = c_nppair(ics,n+2)
              end if
              !
              ! TEST WRITE CASCADE - SPECIAL ASSIGNED SPECIES
              !
              if ( c_kkw >= 1 ) then
                write(c_out,*) 'AUTOMATIC CASCADE SPECIES:', n,c_tchem(icc)
              end if
              do neq = 1 , (c_nequil(icc)+1)
                ic = icc
                if ( neq > 1 ) ic = c_ncequil(icc,(neq-1))
                lcastest(ic) = .true.
              end do
            end do
          end if
        end do
        if ( ncdim < nsolv ) ncdim = nsolv
        if ( nsdim < nsol ) nsdim = nsol
      end do
      !
      ! TEST WRITE
      !
      if ( c_kkw > 0 ) then
        write(c_out,336) ncdim, nsdim
      end if
      if ( c_kkw == 1 .or. c_kkw == 5) then
        write(c_out,342)
      end if
      !
      ! CATEGORY TEST LOOP AMONG SPECIES.
      !   IDENTIFY SPECIES MISSING FROM CASCADE AND MISSING CATEGORIES.
      !
      do ic = 1 , c_nchem2
        !
        !  Control for diagnostic write
        !
        if ( c_kkw == 1 .or. c_kkw == 5 ) then
          write(c_out,343) ic, c_tchem(ic), c_icat(ic), c_lsts(ic), &
                      c_nequil(ic), c_npequil(ic),  c_nppair(ic,3), &
                      c_nppair(ic,2), c_nppair(ic,1)
          if ( c_npequil(ic) > 0 ) then
            write(c_out,309) c_tchem(c_npequil(ic))
          end if
          if ( c_nppair(ic,2) > 0 ) then
            write(c_out,309) c_tchem(c_nppair(ic,2))
          end if
          if ( c_nppair(ic,1) > 0 ) then
            write(c_out,309) c_tchem(c_nppair(ic,1))
          end if
          if ( c_nequil(ic) > 0 ) then
            write(c_out,344) (c_ncequil(ic,i),i=1,c_nequil(ic))
            do i = 1 , c_nequil(ic)
              if ( c_ncequil(ic,i) > 0 ) then
                write(c_out,309) c_tchem(c_ncequil(ic,i))
              end if
            end do
          end if
          if ( c_nppair(ic,3) > 0 ) then
            write(c_out,345) ( c_nppair(ic,i),i=4,(3+c_nppair(ic,3)) )
            do i = 4 , (3+c_nppair(ic,3))
              if ( c_nppair(ic,i) > 0 ) then
                write(c_out,309) c_tchem(c_nppair(ic,i))
              end if
            end do
          end if
        end if
        !
        ! End control for diagnostic write.
        !
        if ( c_icat(ic) <= 0 ) then
          write(c_out,341) ic,c_tchem(ic), c_icat(ic)
!         c_icat(ic) = 23
        end if
        if ( .not.lcastest(ic) ) then
          write(c_out,352) ic, c_tchem(ic), c_npequil(ic), &
                          c_tchem(c_npequil(ic))
        end if
      end do
      !
      ! END OF CATEGORY TEST LOOP
      !
      ! -------------------------------------------------------
      ! CASCADE TEST:  RUN THROUGH CASCADE TO MAKE SURE THAT
      ! THERE ARE NO SPECIES SOURCES AFTER THE SPECIES HAS BEEN
      ! SOLVED FOR.
      !
      ! NOTE:  IN REVISED VERSION, PRODUCTION FURTHER DOWN THE CHAIN
      !  IS COUNTED IN THE NEXT ITERATION.
      !  (This is retained as a warning for sources of nonconvergence only)
      !
      !  OLD NOTE:  THE 'PRE' LABEL IN THE CASCADE RUNS REACTIONS ASSOCIATED
      !  WITH A SPECIES USING PRIOR VALUES. THIS OBTAINS REACTION PRODUCTS
      !  EVEN BEFORE THE SPECIES HAS BEEN SOLVED FOR (SLOW SPEC. ONLY.)
      !  BUT THERE MUST BE NO FURTHER PRODUCTION AFTER SOLVE.
      !
      do i = 1 , c_nchem2
        lcastest(i) = .false.
      end do
      !
      ! RUN CASCADE.  SET INDEX=1 AS EACH SPECIES IS SOLVED FOR.
      ! ( DO WHILE NCASCADE(1).NE.0)
      ! LOOP:  RUN CASCADE
      !
      do i = 1 , c_cdim
        if ( c_cascade(i,1) == 0 ) exit
        !
        ! LOOP:  TEST REACTIONS CALLED BY THE CASCADE SPECIES
        !
        do ii = 1 , c_cdim
          icc = c_cascade(i,ii)
          if ( icc > 0 ) then
            icc = c_npequil(icc)
            ic = icc
            !
            ! NOTE PROBLEM HERE (write 419): TEST OF REACTANTS AND PRODUCTS
            !   BELONGS IN QUADINIT
            ! FIRST TEST REACTANTS AND PRODUCTS FROM CASCADE REACTIONS.
            !  EXCEPT FOR 'POST' FLAG, PRODUCTS SHOULD BE ALL DOWN-CASCADE.
            !  CASCADE MODIFIED TO REFLECT ALL CHEMICAL PAIRS
            !  ->POINT TO PRIMARY SPEC.
            !
            if ( c_cascade(i,2) /= -2 ) then
              nnr = c_nnrchem(ic)
              if ( nnr > 0 ) then
                do n = 1 , nnr
                  nr = c_nrchem(ic,n)
                  if ( nr > 0 ) then
                    icr1 = c_reactant(nr,1)
                    if ( icr1 > 0 ) then
                      icr1 = c_nppair(c_npequil(icr1),2)
                    end if
                    icr2 = c_reactant(nr,2)
                    if ( icr2 > 0 ) then
                      icr2 = c_nppair(c_npequil(icr2),2)
                    end if
                    if ( c_nnpro(nr) > 0 ) then
                      do nn = 1 , c_nnpro(nr)
                        icp = c_product(nr,nn)
                        if ( icp > 0 ) then
                          icp = c_nppair(c_npequil(icp),2)
                        end if
                        if ( icp == icr1 ) then
                          icr1 = 0
                          icp = 0
                        end if
                        if ( icp == icr2 ) then
                          icr2 = 0
                          icp = 0
                        end if
                        !
                        ! PRODUCT TEST
                        !
                        if ( icp > 0 ) then
                          if ( lcastest(icp) ) then
                            write(c_out,419) c_tchem(icp), nr, c_tchem(ic), &
                                            (c_treac(jj,nr),jj=1,5)
                            write(6,419) c_tchem(icp), nr, c_tchem(ic), &
                                            (c_treac(jj,nr),jj=1,5)
                          end if
                        end if
                      end do
                    end if
                    !
                    ! REACTANT TEST
                    !
                    if ( icr1 > 0 ) then
                      if ( lcastest(icr1) ) then
                        write(c_out,419) c_tchem(icr1), nr, c_tchem(ic), &
                                        (c_treac(jj,nr),jj=1,5)
                        write(6,419) c_tchem(icr1), nr, c_tchem(ic), &
                                        (c_treac(jj,nr),jj=1,5)
                      end if
                    end if
                    if ( icr2 > 0 ) then
                      if ( lcastest(icr2) )then
                        write(c_out,419) c_tchem(icr2), nr, c_tchem(ic), &
                                        (c_treac(jj,nr),jj=1,5)
                        write(6,419) c_tchem(icr2), nr, c_tchem(ic), &
                                        (c_treac(jj,nr),jj=1,5)
                      end if
                    end if
                  end if
                end do
              end if
            end if
            !
            ! LAST, ADD CASCADE REACTIONS TO REACTION LIST, UNLESS 'pre' FLAG
            !
            if ( c_cascade(i,2) /= -1 ) then
              if ( lcastest(ic) ) then
                write(c_out,418) c_tchem(ic), i, (c_cascade(i,jj),jj=1,3)
                write(6,418) c_tchem(ic), i, (c_cascade(i,jj),jj=1,3)
              end if
              lcastest(ic) = .true.
            end if
          end if
        end do
        !
        ! END LOOP:  TEST REACTIONS CALLED BY THE CASCADE SPECIES
        !
      end do
      !
      ! END LOOP:  RUN CASCADE
      !
      !
      ! -------------------------------------------------------
      ! END CASCADE TEST
      ! -------------------------------------------------------
      !
      ! OPTION:  call cheminit here to set initial chemistry parameters
      !
      ! ----------------------------------------------------------------
      !  END READ.
      ! ----------------------------------------------------------------

   11 format(a4)
   12 format(i6)
   13 format(1pe10.3)
   29 format(5(3x,a8,i3,2x))
  901 format(/,'WARNING: NCHEM MAY EXCEED DIMENSION', &
               ' IN COMMON BLOCK.  NCHEM=',i4)
   49 format(' TOTAL NUMBER OF CHEMICALS READ =', i5)
 1301 format(/,'TRANSPORTED (INPUT) SPECIES:  ',  &
               '  (negative=steady state)'   ,/)
 1302 format(5(i3,   a8,i3,2x))
  209 format(3(   a8,3x))
  219 format(/,'MAJOR ERROR: UNKNOWN SPECIES IN LUMPED SPECIES READ', &
             /, 'SPECIES LIST = ',1(a8,4x),'UNKNOWN #',i3)
   39 format(5(3x,a8,5x))
  218 format(/,'MAJOR ERROR: PREVIOUSLY SET SPECIES IN ', &
             /, 'LUMPED LIST SPECIES = ',1(a8   ))
 1310 format(/, 'LUMPED SPECIES: NAME AND SPECIES CONTENTS')
 1311 format(/,   a8)
  903 format(/,'WARNING: NHENRY MAY EXCEED DIMENSION', &
               ' IN COMMON BLOCK.  NHENRY=',i4)
   59 format(3(   a8,3x))
   69 format(/,'WARNING:  UNRECOGNIZED CHEMICAL NAME:',a8,/, &
               ' NRH=',i4,'  NAME # ',i4,'  HENRY:  ',a8,'=',a8)
   33 format(5x,i5, 1pe10.3, 0pf10.0, 1pe10.3, 0pf10.0)
   57 format(' TOTAL NUMBER OF HENRYS-LAW COEFFICIENTS READ =', i5)
  904 format(/,'WARNING: NAQUEOUS MAY EXCEED DIMENSION', &
               ' IN COMMON BLOCK.  NAQUEOUS=',i4)
   88 format(/,'AQUEOUS ORDER SWITCH: NUMBER:',i3,'  REACTION:',&
             a8,'=',a8,'&',a8)
   89 format(/,'WARNING:  UNRECOGNIZED CHEMICAL NAME:',a8,/, &
               ' NRQ=',i4,'  NAME # ',i4,' AQUEOUS:  ',a8,'=',a8,' + ',a8)
   79 format(' TOTAL NUMBER OF AQUEOUS EQUILIBRIUM CONSTANTS =', i5)
  905 format(/,'WARNING: NAQUEOUSQ MAY EXCEED DIMENSION', &
               ' IN COMMON BLOCK.  NAQUEOUSQ=',i4)
   87 format(/,'WARNING:  UNRECOGNIZED CHEMICAL NAME:',a8,/, &
               ' NRQQ=',i4,'  NAME # ',i4,' AQUEOUSQ:  ',a8,'=',a8,' + ',a8)
   99 format(' TOTAL NUMBER OF AQUEOUSQ EQUILIBRIUM CONSTANTS =', i5)
  906 format(/,'WARNING: NREACT MAY EXCEED DIMENSION', &
               ' IN COMMON BLOCK.  NREACT  =',i4)
  109 format(a8,3x,a8,4x,2(f6.3,2x,a8,3x),f6.3,2x,a8)
  127 format(' ID PRODUCTS:  NR, RK1=', i5,1pe10.3,/, &
             '  TREACT= ',a8,'+',a8,'->',a8,'+',a8,'+',a8)
  128 format(' J IC STOIC PRODARR = ', 2i5,2f10.3)
  139 format(/,'WARNING:  UNRECOGNIZED CHEMICAL NAME:',a8,/, &
               ' NR=',i4,'  NAME # ',i4,'  TREACT= ', &
               a8,'+',a8,'->',a8,'+',a8,'+',a8)
  131 format(2i5)
  132 format(10x,1pe10.3,0pf8.0,0pf8.3)
  133 format(10x,0pf6.2,1pe10.3,0pf6.1,1pe10.3,0pf6.1)
  134 format(10x,1pe10.3,0pf8.0,0pf6.2,1pe10.3,0pf6.1,1pe10.3,0pf6.1)
  135 format(10x,1pe10.3,0pf8.0,4(1pe10.3))
  136 format(10x,3(1pe10.3,0pf8.0) )
  137 format(10x,1pe10.3,0pf8.0,1pe10.3,0pf8.3,0pf8.0)
  142 format(10x,1pe10.3,1x,1pe10.3,1x,1pe10.3,1x,1pe10.3)
  309 format(3(a8,2x))
  347 format(/,'CASCADE TEST: LIST OF SPECIES IN ORDER OF SOLUTION')
  348 format(i4,4x,a8)
  336 format(/,' MAXIMUM CASCADE SPECIES, PAIR GROUPS = ', 2i5)
  342 format(/'TEST IC CHEM ICAT ISTS  NEQUIL PRIMARY  NPAIR PRIMARY')
  343 format(i4,2x,a8,2x,i5,l4,2x,8i5)
  344 format(' c_ncequil(aqueous subspecies) = ',10i5)
  345 format(' ncpair (chem. paired subspecies)  = ',10i5)
  319 format(/,'MAJOR ERROR: UNKNOWN SPECIES IN CASCADE SPEC READ.', &
             /, 'SPECIES LIST = ',3(a8,4x),'UNKNOWN #',i3)
 1381 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING', &
               ' OR MISNAMED:  H2O')
 1382 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING', &
               ' OR MISNAMED:  CO2')
 1383 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING', &
               ' OR MISNAMED:   H+')
 1384 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING', &
               ' OR MISNAMED:  OH-')
 1385 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING', &
               ' OR MISNAMED:   OH')
 1386 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING', &
               ' OR MISNAMED:  HO2')
 1387 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING', &
               ' OR MISNAMED: H2O2')
 1388 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING', &
               ' OR MISNAMED:   O3')
 1389 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING', &
               ' OR MISNAMED:  NO2')
 1390 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING', &
               ' OR MISNAMED:   NO')
 1391 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING', &
               ' OR MISNAMED:  NO3')
 1392 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING', &
               ' OR MISNAMED: N2O5')
 1393 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING', &
               ' OR MISNAMED: HNO3')
  659 format(a8,'=',a8,'+',a8,2x, 'ION2=',i3, '  ION3=',i3)
 1321 format(/,'HENRYs LAW AND LINKED AQUEOUS EQUILIBRIUM SPECIES',/)
 1322 format(5(a8,2x))
  667 format(/,' MAJOR ERROR IN CASCADE PAIR SETUP:',/, &
               '  EITHER ZERO OR SECOND CASCADE PAIR ALREADY PAIRED.', &
               /,' PAIR NUMBER, icc1, icc2, ic1, ic2, c_nppair(ic2,1)=',6i5)
  341 format(/,' WARNING:  MISSING SPECIES CATEGORY',/, &
               'SPECIES =',i4,2x,a8,'  CATEGORY=',i4)
  352 format(/,' WARNING:  SPECIES OMITTED FROM CASCADE',/, &
               'SPECIES =',i4,4x,a8,'  NPEQUIL =',i4,2x,a8)
  419 format(/,'WARNING:  CASCADE OUT OF SEQUENCE FOR SPECIES NAME=  ', &
            a8,/,'PRODUCED AT REACTION #',i5,'  CALLED FOR SPECIES= ',  &
            a8,/,' REACTION: ',a8,'+',a8,'=>',a8,'+',a8,'+',a8)
  418 format(/,'WARNING:  SPECIES CALLED TWICE IN CASCADE. SPECIES=',   &
            a8,/,' SECOND CALL AT CASCADE#',i3,' WITH INDICES=',3i3)

    end subroutine chemread
!
! ----------------------------------------------------------------
!
! integer function namechem(titl)
!
! INPUT:  chemical species name, char*8
! OUTPUT: Species index number (integer)  associated with the name,
!           or zero if not found in species list.
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
!
    integer function namechem(titl)
      implicit none
      character(len=8) , intent(in) :: titl
      integer(ik4) :: i
      namechem = 0
      if ( titl == '        ' ) return
      do i = 1 , c_nchem2
        if ( titl == c_tchem(i) ) then
          namechem = i
          exit
        end if
      end do
    end function namechem
!
! ----------------------------------------------------------------
!
! This reads the MADRONICH LOOKUP TABLE (2002 VERSION).
!   Input file:  TUVGRID2  (kept in dhvmad/TUVcode4.1a)
!    (see TUVINFO in /l/kudzu/k-1/sillman/dhvmad)
!                 (also tuvtab2.f, tested in tuvtest2.f)
!
!   OPTION:  Link to jval2.f for altitude in km (standard)
!            or to   jval1.f for altitude in kPa
!
! Result: creates the hv input arrays:
!     c_hvin,c_nhv,c_hvmat, c_hvmatb, c_jarray
!
! Called by:  boxmain (as part of chemical setup)
! Calls to:   readhv (in file jval2.f (KM) or jval1.f (kPA).)
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
    subroutine hvread
      implicit none
      call readhv(c_hvin,c_nhv,c_hvmat,c_hvmatb,c_jarray)
    end subroutine hvread
!
! --------------------------------------------------------------
!
! QUADINIT does a preliminary analysis of the chemical mechanism
!  to link reactions with individual species for the solution procedure
!  count odd hydrogen,
!  set up species pair chains for solution, etc.
!
!  Initial read for aqueous equilibrium was here, now in chemread
!
!  DETAILS:
!   1.  Each reaction is assigned to a species.
!        The numerical solution solves for each species in order
!        using the reactions needed for that species.
!
!   2.  Counts odd hydrogen and odd nitrogen -
!                     net change for each  reaction.
!
!   3.  Identifies back-forth 'exchange' reaction numbers
!         associated with direct exchange of species
!         (e.g. special equilibria, MCO3-PAN, HNO4, etc - for excorr)
!
!   4.  Assigns reactions associated with 'cross-production' (cpro)
!        among paired species in 'cascade' solution (FEB 2005)
!
!   5.  Identifies 'pairfac': conservation factor for pair group species
!         based on back-forth exchange reactions (FEB 2005)
!         (e.g. Br+Br- <=>Br2- : Br2- counts double in pair group.)
!
!  FUTURE ADD:  Add category for soluble aerosol/aqueous species,
!          and for aerosol odd-H.
!      Part of algorithm for non-steady state gas-aq partitioning
!       where initial gas/aqueous partitioning is set by input
!       for gas vs aqueous/aerosol species
!
!  FUTURE ADD:   Modify to skip aqueous reactions if LWC=0:
!      assign reactants to direct species, not gas pointer
!      and in chem, do reactions only when species is called
!           - if LWC=0 skip aqueous
!
! ****NOTE, TREATMENT OF AQUEOUS EQUILIBRIUM SPECIES NEEDS TO
!     BE REWRITTEN AND CHANGED.
!
! ---------------------------------------------
! Called by:  boxmain (as part of chemical setup)
!              or by chemread
!   Called just once, as part of setup.
! Calls to:   none
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
! ---------------------------------------------------------------------
!
    subroutine cheminit
!
      implicit none
      ! Chem index
      integer(ik4) :: ic , ic1 , ic2 , icc
      ! Chem local index
      integer(ik4) :: is
      ! Chem index
      integer(ik4) :: icx , icx1 , icx2 , icy1 , icy2 , icp , icp1 , icp2
      ! Reaction counters
      integer(ik4) :: nr , nrx
      ! indices used for species categories
      integer(ik4) :: icat1 , icat2 , icatp
      ! Vectorization counters
      integer(ik4) :: kk
      ! General counters
      integer(ik4) :: i , j , ii , n
!
! LOCAL VARIABLES
! lloss        Local: Flag for identifying exchange loss reaction
! lpro         Local: Flag for identifying exchange production reaction
!  xoddhx      Counter for net change in odd hydrogen in reaction
!                (intermediate for c_oddhx) - split into RO2; OH-HO2-CO3
!  xpronox     Counter for net production of odd nitrogen in reaction
!                (intermediate for c_pronox)
!
      ! Flag for identifying exchange loss reaction
      logical :: lloss
      ! Flag for id exchange production reaction
      logical :: lpro
      ! Counter for odd hydrogen RO2 only
      real(rk8) :: xoddhx
      ! Counter for odd hydrogen w/o RO2
      real(rk8) :: xoddhx2
      ! Counter for odd nitrogen
      real(rk8) :: xpronox
      ! Counter for product reactions
      integer nrp
      !
      ! -----------------------------------------
      !  ANALYZE CHEMICAL MECHANISM:
      !   FOR EACH REACTION IDENTIFY THE (ONE) SPECIES THAT THE REACTION
      !   IS LINKED WITH IN THE SOLVER.  ALWAYS GAS-PHASE.
      !    ALSO CALCULATE ODD HYDROGEN AND ODD NITROGEN CHANGE
      ! -----------------------------------------
      !
      ! REACTIONS ARE ASSIGNED TO SPECIES IN REACTANT-TO-PRODUCT SEQUENCE
      ! WITH SPECIAL CATEGORIES (NO3, NOx, Hx) AT THE END.
      !
      ! REACTIONS LINKED WITH A (GAS-PHASE) SPECIES INCLUDE:
      !   (1) ALL SPECIES LOSS REACTIONS
      !   (2) SPECIES PRODUCTION REACTIONS ONLY IF PRODUCTION IS FROM A
      !       "SPECIAL" CATEGORY THAT IS NOT INCLUDED
      !                    IN THE REACTION-TO-PRODUCT CASCADE.
      !   (3) REACTIONS FOR ASSOCIATED AQUEOUS SPECIES.
      !   (4) REACTIONS FOR SPECIES LINKED THROUGH SPECIAL AQUEOUS EQUILIBRIA.
      !
      ! THE IDENTIFICATION OF REACTIONS WITH SPECIES INSURES THAT
      ! ALL REACTIONS WILL BE PROCESSED AT THE PROPER PLACE IN THE CASCADE.
      !
      ! (NEW) LOSS AND PRODUCTION REACTIONS ARE ID'D SEPARATELY.
      ! -----------------------------------------
      kk = 1
      !
      !  ZERO THE IMPORTANT ARRAYS
      !
      do ic = 1 , c_nchem2
        c_nnrchem(ic) = 0
        c_nnrchp(ic) = 0
        do i = 1 , 25
          c_nrchem(ic,i) = 0
          c_nrchmp(ic,i) = 0
        end do
      end do
      do nr = 1 , c_nreac
        c_stoiloss(nr) = d_zero
      end do
      !
      ! LOOP FOR EACH INDIVIDUAL REACTION
      !
      loopmechanal: &
      do nr = 1 , c_nreac
        !
        ! ZERO INDEX FOR IC ASSOCIATED WITH REACTION.
        !
        ic = 0
        !
        ! ESTABLISH REACTANT CATEGORIES.
        !     NOTE: REACTION IS ASSIGNED ONLY FOR ICAT>0.  CHECK ERROR.
        !
        icat1 = 0
        icat2 = 0
        if ( c_reactant(nr,1) > 0 ) then
          icat1 = c_icat(c_reactant(nr,1))
        end if
        if ( c_reactant(nr,2) > 0 ) then
          icat2 = c_icat(c_reactant(nr,2))
        end if
        !
        ! FIRST ELIMINATE REDUNDANT REACTANTS:
        !   REACTANTS THAT ARE ALSO PRODUCTS,  AND ALSO H2O
        !
        !  (This is for reactions RCO3+RO2=> products.
        !    Entered as two reactions:
        !       RCO3+RO2=>RCO3+productB
        !       RCO3+RO2=>RO2 +productA
        !    Assign reaction with RO2 products to RO2, not RCO3.)
        !
        if ( c_nnpro(nr) > 0 ) then
          do n = 1 , c_nnpro(nr)
            if ( c_reactant(nr,1) == c_product(nr,n) .and. &
                 c_stoich(nr,n) >= d_one ) then
              icat1 = 0 - icat1
            end if
            if ( c_reactant(nr,2) == c_product(nr,n) .and. &
                 c_stoich(nr,n) >= d_one ) then
              icat2 = 0 - icat2
            end if
          end do
        end if
        !
        ! FIRST ASSOCIATE REACTION WITH REACTANTS
        !
        if ( icat1 < 9 .and. icat1 > 0 ) then
          ic = c_reactant(nr,1)
        else
          if ( icat2 < 9 .and. icat2 > 0 ) then
            ic = c_reactant(nr,2)
          end if
        end if
        !
        ! RESTORE CATEGORY OF REDUNDANT REACTIONS (RCO3+RO2)
        !
        icat1 = iabs(icat1)
        icat2 = iabs(icat2)
        !
        ! (PRIOR:  Look for NO3-N2O5 reactants before non-special products.
        !  Now, Look for nonspecial products first, then NO3-N2O5 reactants.)
        !
        ! IF NEITHER REACTANT IS A NON-SPECIAL SPECIES, THEN CHECK PRODUCTS
        !  (e.g. HNO3, H2O2)
        !
        if ( ic <= 0 ) then
          if ( c_nnpro(nr) > 0 ) then
            do n = 1 , c_nnpro(nr)
              icatp = c_icat(c_npequil(c_product(nr,n)))
              if ( icatp < 9 .and. icatp > 0 ) then
                ic = c_product(nr,n)
                exit
              end if
            end do
          end if
        end if
        !
        ! IF NEITHER REACTANT IS A NON-SPECIAL SPECIES, LOOK FOR NO3-N2O5-HNO3.
        !
        if ( ic <= 0 ) then
          if ( icat1 == 14 .or. icat1 == 15 .or. icat1 == 16 ) then
            ic = c_reactant(nr,1)
          else
            if ( icat2 == 14 .or. icat2 == 15 .or. icat2 == 16 ) then
              ic = c_reactant(nr,2)
            end if
          end if
        end if
        !
        ! IF NOTHING, LOOK FOR NO3/N2O5/HNO3 IN PRODUCTS
        !
        if ( ic <= 0 ) then
          if ( c_nnpro(nr) > 0 ) then
            do n = 1 , c_nnpro(nr)
              icatp = c_icat(c_npequil(c_product(nr,n)))
              if ( icatp == 14 .or. icatp == 15 .or. icatp == 16 ) then
                ic = c_product(nr,n)
                exit
              end if
            end do
          end if
        end if
        !
        ! IF STILL NOTHING, LOOK FOR O3 or NOx IN REACTANTS OR PRODUCTS
        !    O3 (11), NO2 (12) or NO (13)
        ! 2005 CHANGE   O3, NO, NO2 ORDER DOESN'T MATTER,
        !   BUT DECLARE THROUGH REACTANTS NOT PRODUCTS - else brpro too early
        !   (HO2+NO=>NO2+OH must  be assigned to NO2)
        !   MUST BE  11 or 12 or 13 here
        !
        if ( ic <= 0 ) then
          if ( icat1 == 11 .or. icat1 == 12 .or. icat1 == 13 ) then
            ic = c_reactant(nr,1)
          else
            if ( icat2 == 11 .or. icat2 == 12 .or. icat2 == 13 ) then
              ic = c_reactant(nr,2)
            else
              if ( c_nnpro(nr) > 0 ) then
                do n = 1 , c_nnpro(nr)
                  icatp = c_icat(c_npequil(c_product(nr,n)))
                  if ( icatp == 11 .or. icatp == 12 .or. icatp == 13 ) then
                    ic = c_product(nr,n)
                    exit
                  end if
                end do
              end if
            end if
          end if
        end if
        !
        ! H2O2 OPTION:  MAKE H2O2 A SPECIAL SPECIES FOR ASSIGNING REACTIONS:
        ! LOOK FOR H2O2 (19 ) IN REACTANTS OR PRODUCTS.
        !
        if ( ic <= 0 ) then
          if ( icat1 == 19 ) then
            ic = c_reactant(nr,1)
          else
            if ( icat2 == 19 ) then
              ic = c_reactant(nr,2)
            else
              if ( c_nnpro(nr) > 0 ) then
                do n = 1 , c_nnpro(nr)
                  icatp = c_icat(c_npequil(c_product(nr,n)))
                  if ( icatp == 19 ) then
                    ic = c_product(nr,n)
                    exit
                  end if
                end do
              end if
            end if
          end if
        end if
        !
        ! IF STILL NOTHING, THE REMAINDER HAD BETTER HAVE OH OR HO2!
        !  (MAYBE TRY TO USE OH RATHER THAN HO2?)
        !
        if ( ic <= 0 ) then
          if ( icat1 == 9 .or. icat1 == 10 ) then
            ic = c_reactant(nr,1)
          else
            if ( icat2 == 9 .or. icat2 == 10 ) then
              ic = c_reactant(nr,2)
            else
              if ( c_nnpro(nr) > 0 ) then
                do n = 1 , c_nnpro(nr)
                  icatp = c_icat(c_npequil(c_product(nr,n)))
                  if ( icatp == 9 .or. icatp == 10 ) then
                    ic = c_product(nr,n)
                    exit
                  end if
                end do
              end if
            end if
          end if
        end if
        !
        ! IF STILL IC=0, ERROR AND EXIT!! (unless reaction=0)
        !
        if ( ic <= 0 .and. &
            (c_reactant(nr,1) > 0 .or. c_reactant(nr,2) > 0) ) then
          write(c_out,355) nr, c_reactant(nr,1), icat1, &
                c_reactant(nr,2), icat2, (c_product(nr,n),n=1,c_nnpro(nr) )
          write(6,355) nr, c_reactant(nr,1), icat1, c_reactant(nr,2), icat2, &
                           (c_product(nr,n),n=1,c_nnpro(nr) )
          call fatal(__FILE__,__LINE__,'MAJOR ERROR IN CBMZ CODE')
        end if
        !
        ! FINALLY!  ENTER REACTION INTO PROPER SPECIES LIST
        ! NEW - SEPARATE LIST FOR SPECIES-LOSS REACTIONS AND PRODUCT REACTIONS,
        !       AND CALC OF STOILOSS FOR LOSS REACTIONS.
        ! 2000-REACTIONS ASSIGNED TO MAIN SPECIES FOR SPECIAL AQUEOUS LINKS.
        !
        !  FUTURE ADD:   Modify to skip aqueous reactions if LWC=0:
        !      assign reactants to direct species, not gas pointer
        !      and in chem, do reactions only when species is called
        !           - if LWC=0 skip aqueous
        !
        ! 92407 PROBLEM HERE:  EXCHANGE reactions.
        !    They should have the same KEY SPECIES
        !
        !   possibility:  set c_nicreac(nr) here; modify in EXCHANGE analysis be
        !    so that both have the same key species (product of A+B=C),
        !   then do this process below.
        !   ALT: require link for paired species??
        !
        if ( ic > 0 ) then
          icc = c_npequil(ic)
          if ( ic == c_reactant(nr,1) .or. ic == c_reactant(nr,2) ) then
            c_nnrchem(icc) = c_nnrchem(icc) + 1
            c_nrchem(icc,c_nnrchem(icc)) = nr
            c_nicreac(nr) = icc
            c_stoiloss(nr) = d_one
            if ( c_reactant(nr,1 ) == c_reactant(nr,2) ) c_stoiloss(nr) = d_two
            if ( c_nnpro(nr) > 0 ) then
              do n = 1 , c_nnpro(nr)
                if ( ic == c_product(nr,n) ) then
                  c_stoiloss(nr) = c_stoiloss(nr)-c_stoich(nr,n)
                end if
              end do
            end if
          else
            c_nnrchp(icc) = c_nnrchp(icc) + 1
            c_nrchmp(icc,c_nnrchp(icc)) = nr
          end if
        end if
        !
        ! END ANALYSIS OF SPECIES-LINK WITH REACTIONS.
        !
        !  CALCULATE CHANGE IN ODD HYDROGEN (ODDHX)
        !  AND NET PRODUCTION OF ODD NITROGEN (PRONOX, >0).
        !  (CURRENT OPTION - RECORD NET PRONOX, THEN IN PROGRAM
        !   SUM NOx PRODUCTION-ONLY.)
        !  (ODDHX IS NET CHANGE, +/-.  BUT PRONOX IS PRODUCTION ONLY.
        !   NET CHANGE IN NOX CAN BE OBTAINED FROM RP, RL FOR NO AND NO2,
        !   BUT NOX ANALYSIS REQUIRES SEPARATION OF NOX SOURCES AND SINKS.
        !   NOTE - PAN->NOX AND HNO3->NOX OVERPRODUCTION IS FIXED
        !   BY PANCORR AND NO3CORR - DISCOUNTS BACK-FORTH NOX PRODUCTION.)
        !
        xoddhx = d_zero
        xoddhx2 = d_zero
        xpronox = d_zero
        if ( icat1 == 9 .or. icat1 == 10 .or. icat1 == 8) then
          xoddhx2 = xoddhx2 - d_one
        end if
        if ( icat2 == 9 .or. icat2 == 10 .or. icat1 == 8) then
          xoddhx2 = xoddhx2 - d_one
        end if
        if ( icat1 == 3 ) then
          xoddhx = xoddhx - d_one
        end if
        if ( icat2 == 3 ) then
          xoddhx = xoddhx - d_one
        end if
        if ( icat1 == 12 .or. icat1 == 13 ) xpronox = xpronox - d_one
        if ( icat2 == 12 .or. icat2 == 13 ) xpronox = xpronox - d_one
        if ( c_nnpro(nr) > 0 ) then
          do n = 1 , c_nnpro(nr)
            icatp = c_icat(c_npequil(c_product(nr,n)))
            if ( icatp == 9 .or. icatp == 10 .or. icatp == 8) then
              xoddhx2 = xoddhx2 + c_stoich(nr,n)
            end if
            if ( icatp == 3 ) xoddhx = xoddhx + c_stoich(nr,n)
            if ( icatp == 12 .or. icatp == 13 ) then
              xpronox = xpronox + c_stoich(nr,n)
            end if
          end do
        end if
        !
        ! GOT HERE
        !
        c_oddhx(nr,1) = xoddhx + xoddhx2
        c_oddhx(nr,2) = xoddhx2
        c_pronox(nr) = xpronox
        !
        ! ADDITION FOR ODD HYDROGEN:  IDENTIFY THE HO2->OH DIRECT REACTIONS.
        !
        if ( c_nnpro(nr) > 0 ) then
          if ( c_npequil(c_reactant(nr,1)) /= c_nho2 ) then
            if ( c_reactant(nr,2) > 0 ) then
              if ( c_reactant(nr,2) == c_nho2 ) then
                do n = 1 , c_nnpro(nr)
                  if ( c_npequil(c_product(nr,n)) == c_noh ) then
                    exit
                  end if
                end do
              end if
            end if
          end if
        end if

      end do loopmechanal
      !
      ! END REACTION LOOP - END MECHANISM ANALYSIS.
      !
      ! ----------------------------------------------------------
      ! IDENTIFY 'EXCHANGE' REACTIONS WITH THE FORM A+B=>C; C=>A+B.
      ! AND PAIRFAC/MULTIFAC CONSERVATION INDICES.
      !
      ! THESE REACTION NUMBERS ARE RECORDED IN EXPRO(IC,I) AND EXLOSS(IC,I)
      ! WHERE IC IS FOR THE PRODUCT C  (GAS-PHASE).
      !
      ! THE EXCHANGE COEFFICIENTS INVOKE THE SUBROUTINE HNCORR (w/PANCORR)
      ! TO INSURE THAT THE BACK-AND-FORTH PRODUCTION AND LOSS
      ! DO NOT INTERFERE WITH BACK-EULER SOLUTION FOR PRODUCT SPECIES
      !  (USUALLY NOX OR HX).
      ! THIS IS ALSO IMPORTANT FOR AQUEOUS:  HCO3+O2- ->HO2- + CO3-.
      !
      ! NOTE:  THIS ONLY WORKS IF EXCHANGE REACTIONS HAVE THE SAME KEY SPECIES
      ! OR IF THE TWO KEY SPECIES ARE DONE WITH TWOSOLVE (RCO3-PAN).
      !
      ! NOTE:  EXLOSS AND EXPRO ARE ENTERED AS PAIRED REACTIONS;
      !  BUT SAME EXPRO MAY GO TO MORE THAN ONE EXLOSS, OR  VICE VERSUS.
      !  (so same reaction may be listed twice)
      !
      ! ----------------------------------------------------------
      !  ZERO THE IMPORTANT ARRAYS
      !
      do ic = 1 , c_nchem2
        do i = 1 , 20
          c_exspec(ic,i) = 0
          do ii = 1 , 5
            c_expro(ic,i,ii) = 0
            c_exloss(ic,i,ii) = 0
          end do
        end do
      end do
      !
      ! INITIALIZE PAIRFAC (conservation parameter for pair species)
      !   = 1 unless adjusted by EXCHANGE REACTIONS
      !
      do ic = 1 , c_nchem2
        c_pairfac(ic) = d_one
        c_multfac(ic) = d_one
      end do
      !
      ! --------
      ! BEGIN REACTION LOOP FOR EXCHANGE REACTIONS (1200)
      !
      do nrx = 1 , c_nreac
        !
        ! EXLOSS REACTION:
        !  ICX (NICREAC) is key species for reaction (1st non-special reactant).
        !
        icx = c_nicreac(nrx)
        !
        !  Identify up to two REACTANTS ICY1, ICY2;
        !               and up to two PRODUCTS ICX1,ICX2.
        !
        icy1 = c_npequil(c_reactant(nrx,1))
        icy2 = 0
        if ( c_reactant(nrx,2) > 0 ) then
          icy2 = c_npequil(c_reactant(nrx,2))
        end if
        icx1 = 0
        icx2 = 0
        if ( c_product(nrx,1) > 0 ) icx1 = c_npequil(c_product(nrx,1))
        if ( c_product(nrx,2) > 0 ) icx2 = c_npequil(c_product(nrx,2))
        icat1 = 0
        icat2 = 0
        icatp = 0
        if ( icy1 > 0 ) then
          icat1 = c_icat(icy1)
        end if
        if ( icy2 > 0 ) then
          icat2 = c_icat(icy2)
        end if
        if ( icx1 > 0 ) then
          icatp = c_icat(icx1)
        end if
        !
        ! CONTROLS TO SKIP LOOP (1200):
        !
        ! SKIP AND CONTINUE LOOP IF (a) TWO REACTANTS, ONE PRODUCT;
        !  OR IF (b) KEY SPECIES IS A PRODUCT SPECIES.  (IN WHICH CASE
        !  THIS REACTION WOULD BE THE EXCHANGE 'RETURN' - EXPRO, not EXLOSS)
        if ( icy2 > 0 .and. icx2 == 0 ) then
        else
          if ( icx /= icy1 .and. icx /= icy2 .and. icy2 > 0 ) then
          else
            !
            ! OPTIONAL CRITERIA:  EXCHANGE FOR FAST SPECIES ONLY.
            ! FOR SLOW SPECIES, SKIP AND CONTINUE LOOP
            !
            if ( c_icat(icy1) == 1 .or. c_icat(icy1) == 4 ) then
            else
              !
              ! SEARCH THROUGH REACTIONS TO FIND MATCHING EXPRO
              ! REACTION(s) (1300)
              do nrp = 1 , c_nreac
                !
                !  ICP (NICREAC) is key species for EXPRO reaction
                !
                icp = c_nicreac(nrp)
                !
                !  Identify EXPRO reactants IC1,  IC2 and products ICP1, ICP2.
                ic1 = 0
                ic2 = 0
                if ( c_reactant(nrp,1) > 0 ) then
                  ic1 = c_npequil(c_reactant(nrp,1))
                end if
                if ( c_reactant(nrp,2) > 0 ) then
                  ic2 = c_npequil(c_reactant(nrp,2))
                end if
                icp1 = 0
                icp2 = 0
                if ( c_product(nrp,1) > 0 ) then
                  icp1 = c_npequil(c_product(nrp,1))
                end if
                if ( c_product(nrp,2) > 0 ) then
                  icp2 = c_npequil(c_product(nrp,2))
                end if
                !
                ! CONTROLS TO IDENTIFY  EXCHANGE REACTION
                ! (EXLOSS-EXPRO):  A+B=C, C=A+B
                if ( (icy1 == icp1 .and. icy2 == icp2) .or. &
                     (icy1 == icp2 .and. icy2 == icp1) ) then
                  if ( (icx1 == ic1 .and. icx2 == ic2) .or. &
                       (icx1 == ic2 .and. icx2 == ic1) ) then
                    !
                    ! SET PAIRFAC (pair group conservation of mass)
                    !   FOR EXCHANGE REACTIONS A+B=C
                    !    (note A+B=>C must be return reaction)
                    !
                    if ( ic1 /= 0 .and. ic2 /= 0 .and. icp1 /= 0 ) then
                      if ( c_nppair(ic1,2) == c_nppair(ic2,2) ) then
                        if ( c_nppair(icp1,2) == c_nppair(ic1,2) ) then
                          if ( icp2 == 0 ) then
                            c_pairfac(icp1) = c_pairfac(ic1)+c_pairfac(ic2)
                            !
                            ! TEST WRITE - QUADINIT PAIRFAC (STANDARD)
                            !
                            if ( c_kkw > 0 ) then
                              write(c_out,1209) c_tchem(ic1), c_tchem(ic2), &
                                c_tchem(icp1), c_pairfac(ic1), c_pairfac(ic2), &
                                c_pairfac(icp1)
                            end if
                          else
                            if ( c_nppair(icp2,2) /= c_nppair(ic1,2) ) then
                              c_pairfac(icp2) =  c_pairfac(ic1)+c_pairfac(ic2)
                              !
                              ! TEST WRITE - QUADINIT PAIRFAC (STANDARD)
                              !
                              if ( c_kkw > 0 ) then
                                write(c_out,1209) c_tchem(ic1), c_tchem(ic2), &
                                                  c_tchem(icp1)
                              end if
                            end if
                          end if
                        else
                          if ( icp2 > 0 ) then
                            if ( c_nppair(icp2,2) == c_nppair(ic1,2) ) then
                              c_pairfac(icp2) = c_pairfac(ic1)+c_pairfac(ic2)
                              !
                              ! TEST WRITE - QUADINIT PAIRFAC (STANDARD)
                              !
                              if ( c_kkw > 0 ) then
                                write(c_out,1209) c_tchem(ic1), c_tchem(ic2), &
                                                  c_tchem(icp1)
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                    !
                    !
                    ! SET MULTIFAC (multi group conservation of mass)
                    !   FOR EXCHANGE REACTIONS A+B=C
                    !    (note A+B=>C must be return reaction)
                    !
                    if ( ic1 /= 0 .and. ic2 /= 0 .and. icp1 /= 0 ) then
                      if ( c_npmulti(ic1,1) == c_npmulti(ic2,1) ) then
                        if ( c_npmulti(icp1,1) == c_npmulti(ic1,1) ) then
                          if ( icp2 == 0 ) then
                            c_multfac(icp1) = c_multfac(ic1)+c_multfac(ic2)
                            !
                            ! TEST WRITE - QUADINIT MULTIFAC (STANDARD)
                            !
                            if ( c_kkw > 0 ) then
                              write(c_out,1208) c_tchem(ic1), c_tchem(ic2), &
                                c_tchem(icp1), c_multfac(ic1), c_multfac(ic2), &
                                c_multfac(icp1)
                            end if
                          else
                            if ( c_npmulti(icp2,1) /= c_npmulti(ic1,1) ) then
                              c_multfac(icp2) = c_multfac(ic1)+c_multfac(ic2)
                              !
                              ! TEST WRITE - QUADINIT MULTIFAC (STANDARD)
                              !
                              if ( c_kkw > 0 ) then
                                write(c_out,1208) c_tchem(ic1), c_tchem(ic2), &
                                  c_tchem(icp1)
                              end if
                            end if
                          end if
                        else
                          if ( icp2 > 0 ) then
                            if ( c_npmulti(icp2,1) == c_npmulti(ic1,1) ) then
                              c_multfac(icp2) = c_multfac(ic1)+c_multfac(ic2)
                              !
                              ! TEST WRITE - QUADINIT MULTIFAC (STANDARD)
                              !
                              if ( c_kkw > 0 ) then
                                write(c_out,1208) c_tchem(ic1), c_tchem(ic2), &
                                  c_tchem(icp1)
                              end if
                            end if
                          end if
                        end if
                      end if
                    end if
                    !
                    ! SET LOGIC FLAGS TO TEST IF REACTION IS ALREADY LISTED
                    !
                    lloss = .true.
                    lpro  = .true.
                    !
                    ! LOOP TO TEST WHETHER EXCHANGE REACTION IS ALREADY LISTED
                    ! IN REVERSE
                    !   If listed, EXPRO reaction would be  EXLOSS
                    !    for its key species
                    !
                    if ( icp > 0 ) then
                      do i = 1 , 20
                        do ii = 1 , 5
                          !
                          ! IF EXCHANGE RN IS LISTED IN REVERSE, USE LPRO,
                          ! LLOSS to ID reaction as
                          if ( c_exloss(icp,i,ii) == nrp ) then
                            lloss= .false.
                            lpro= .false.
                          end if
                        end do
                      end do
                    end if
                    !
                    ! LOOP TO TEST WHETHER EXLOSS, PRO REACTIONS ARE ALREADY
                    ! IN FORWARD LIST
                    do i = 1 , 20
                      do ii = 1 , 5
                        !
                        ! IF EXCHANGE RN IS LISTED IN REVERSE, USE LPRO,
                        ! LLOSS to ID reaction as
                        !
                        if ( c_exloss(icx,i,ii) == nrx ) lloss = .false.
                        if ( c_expro(icx,i,ii) == nrp )  lpro = .false.
                      end do
                    end do
                    !
                    ! IDENTIFY EXCHANGED SPECIES (icp):
                    !  Record species only if base and exchanged species
                    !  are in same pair group.
                    !   and if exchanged  species is in reaction product list
                    !  If no exchanged species, flag  species number as -1.
                    if ( icp > 0 ) then
                      if ( icp /= icx1 .and. icp /= icx2 ) icp = 0
                      if ( c_nppair(icx,2) /= c_nppair(icp,2) ) icp = 0
                    end if
                    if ( icp == 0 ) icp = -1
                    !
                    ! RECORD EXCHANGE SPECIES ONLY IF MAIN AND EXCHANGE
                    ! SPECIES ARE PAIRED
                    !  AND IF KEY SPECIES OF 2ND REACTION IS PRODUCT OF 1ST
                    !
                    ! (EXCORR must be called for single species also, e.g. HNO4.
                    ! LOOP TO ADD EXCHANGED SPECIES TO EXSPEC LIST
                    !
                    if ( lloss .or. lpro ) then
                      do i = 1 , 20
                        is = i
                        if ( c_exspec(icx,i) == icp ) exit
                        if ( c_exspec(icx,i) == 0 ) then
                          c_exspec(icx,i) = icp
                          exit
                        end if
                        if ( i == 20 ) write(c_out,1206)
                        if ( i == 20 ) write(  6 ,1206)
                      end do
                    end if
                    !
                    ! LOOP TO ADD EXLOSS REACTION
                    !
                    if ( lloss ) then
                      loop1210: &
                      do ii = 1 , 5
                        !
                        ! ENTER THE EXCHANGE REACTION
                        !
                        if ( c_exloss(icx,is,ii) == 0 ) then
                          c_exloss(icx,is,ii) = nrx
                          !
                          ! TEST WRITE - QUADINIT EXCHANGE(STANDARD)
                          !
                          if ( c_kkw > 0 ) then
                            write(c_out,1212) icx ,icp, is, i, nrx,  &
                               c_reactant(nrx,1),c_reactant(nrx,2),  &
                               c_product(nrx,1),c_product(nrx,2)
                          end if
                          if ( c_kkw > 0 .and. c_reactant(nrx,1) > 0 ) then
                            write(c_out,*) c_tchem(c_reactant(nrx,1))
                          end if
                          if ( c_kkw > 0 .and. c_reactant(nrx,2) > 0 ) then
                            write(c_out,*) c_tchem(c_reactant(nrx,2))
                          end if
                          if ( c_kkw > 0 .and. c_product(nrx,1) > 0 ) then
                            write(c_out,*) c_tchem( c_product(nrx,1))
                          end if
                          if ( c_kkw > 0 .and. c_product(nrx,2) > 0 ) then
                            write(c_out,*) c_tchem( c_product(nrx,2))
                          end if
                          !
                          ! HAVING ENTERED THE REACTION, BREAK THE LOOP
                          !
                          exit loop1210
                          ! ERROR WRITE
                          write(c_out,1213)
                          write( 6  ,1213)
                        end if
                      end do loop1210
                    end if
                    !
                    ! END LOOP TO ADD EXLOSS REACTION (1210)
                    !    CONTINUE AFTER BREAK LOOP
                    !
                    ! LOOP TO ADD EXPRO  REACTION
                    !
                    if ( lpro ) then
                      loop1220: &
                      do ii = 1 , 5
                        !
                        ! ENTER THE EXCHANGE REACTION
                        !
                        if ( c_expro(icx,is,ii) == 0 ) then
                          c_expro(icx,is,ii)=nrp
                          !
                          ! TEST WRITE - QUADINIT EXCHANGE(STANDARD)
                          !
                          if ( c_kkw > 0 ) then
                             write(c_out,1222) icx, icp, is, i, nrp,   &
                                 c_reactant(nrp,1), c_reactant(nrp,2), &
                                 c_product(nrp,1), c_product(nrp,2)
                          end if
                          if ( c_kkw > 0 .and. c_reactant(nrp,1) > 0 ) then
                            write(c_out,*) c_tchem(c_reactant(nrp,1))
                          end if
                          if ( c_kkw > 0 .and. c_reactant(nrp,2) > 0 ) then
                            write(c_out,*) c_tchem(c_reactant(nrp,2))
                          end if
                          if ( c_kkw > 0 .and. c_product(nrp,1) > 0 ) then
                            write(c_out,*) c_tchem( c_product(nrp,1))
                          end if
                          if ( c_kkw > 0 .and. c_product(nrp,2) > 0 ) then
                            write(c_out,*) c_tchem( c_product(nrp,2))
                          end if
                          !
                          ! ENTER WARNING FOR EXCHANGE REACTION AMONG
                          ! UNPAIRED SPECIES
                          !
                          if ( icp < 0 ) then
                            if ( (c_reactant(nrp,1) > 11) .or. &
                                 (c_reactant(nrp,2) > 11)) then
                              write(c_out,1224) nrp,  (c_treac(j,nrp),j=1,5)
                            end if
                          end if
                          !
                          ! HAVING ENTERED THE REACTION, BREAK THE LOOP
                          !
                          exit loop1220
                          !
                          ! ERROR WRITE
                          !
                          write(c_out,1223)
                          write( 6  ,1223)
                        end if
                      end do loop1220
                    end if
                    !
                    ! END LOOP TO ADD EXPRO  REACTION (1220)
                    !    CONTINUE AFTER BREAK LOOP
                    !
                    !
                    ! END CONTROLS TO IDENTIFY EXCHANGE REACTIONS
                  end if
                end if
              end do
              !
              ! END LOOP - SEARCH THROUGH REACTIONS TO FIND MATCHING
              ! EXPRO REACTION (1
              ! END - CONTROLS TO SKIP AND CONTINUE LOOP
            end if
          end if
        end if
      end do
      !
      ! END   REACTION LOOP FOR EXCHANGE REACTIONS (1200)
      !
      ! END EXCHANGE ANALYSIS.
      ! ----------------------
      !
      !   IDENTIFY 'CROSS' REACTIONS FOR PAIRED SPECIES IN CASCADE
      ! NNRCPRO, NRCPRO, STOICPRO WITH SAME FORMAT AS NRCHEM, ABOVE.
      !  THESE REPRESENT REACTIONS THAT PRODUCE THE SPECIES FROM 1/2 PAIR.
      !
      ! NNRCPRX, NRCPRX, STOICPRX REPRESENT PRODUCTION REACTIONS ASSOCIATED
      !  WITH THE THIRD PAIR IN A TRIPLET - USED ONLY FOR O3-NO-NO2.
      !
      ! JANUARY 2005:  THIS HAS BEEN DELETED.   IT RELATES TO OLD TWOSOLVE
      !  AND NOXSOLVE - REPLACED BY NEW PAIR AND MULTISOLVE TREATMENTS.
      !
      ! END CROSS-REACTION ANALYSIS.  END MECHANISM ANALYSIS.
      !
      ! ----------------------
      !  NOV 2007:  Identify net RP for NOx, Ox tracers
      ! IDENTIFY REACTION NUMBERS FOR 'SPECIAL' REACTIONS
      !  (used in preliminary chem, oh solver and in NOx-Ox global tracers)
      ! ----------------------
      do i = 1 , 13
        do nr = 1 , c_nreac
          c_noxchem(i,nr) = d_zero
        end do
      end do
      c_nro3no = 0
      c_nrno2x = 0
      c_nro3hv = 0
      c_nrohno2 = 0
      c_nrohco = 0
      c_nrho2no = 0
      c_nrho22 = 0
      c_nro3no2 = 0
      c_nroho3 = 0
      c_nrohch4 = 0
      !
      ! SPECIAL REACTIONS: LOOP FOR EACH INDIVIDUAL REACTION
      !    moved from presolve
      !
      do nr = 1 , c_nreac
        if ( (c_reactant(nr,1) == c_nno .and. &
              c_reactant(nr,2) == c_no3 .and. &
              c_product(nr,1) == c_nno2) .or. &
             (c_reactant(nr,2) == c_nno .and. &
              c_reactant(nr,1) == c_no3 .and. &
              c_product(nr,1) == c_nno2) ) c_nro3no = nr
        if ( (c_product (nr,1) == c_nno .and. &
              c_product (nr,2) == c_no3 .and. &
              c_reactant(nr,1) == c_nno2) .or. &
             (c_product (nr,2) == c_nno .and. &
              c_product (nr,1) == c_no3 .and. &
              c_reactant(nr,1) == c_nno2) ) c_nrno2x = nr
        if ( c_reactant(nr,1) == c_no3 .and. &
             c_reactant(nr,2) == -1 .and. &
             c_product(nr,1) == c_noh) c_nro3hv = nr
        if ( (c_reactant(nr,1) == c_nno2 .and. &
              c_reactant(nr,2) == c_noh) .or. &
             (c_reactant(nr,1) == c_noh .and. &
              c_reactant(nr,2) == c_nno2) ) c_nrohno2 = nr
        if ( (c_reactant(nr,1) == c_nco .and. &
              c_reactant(nr,2) == c_noh) .or. &
             (c_reactant(nr,1) == c_noh .and. &
              c_reactant(nr,2) == c_nco) ) c_nrohco = nr
        if ( (c_reactant(nr,1) == c_nno .and. &
              c_reactant(nr,2) == c_nho2) .or. &
             (c_reactant(nr,1) == c_nho2 .and. &
              c_reactant(nr,2) == c_nno) ) c_nrho2no = nr
        if ( c_reactant(nr,1) == c_nho2 .and. &
             c_reactant(nr,2) == c_nho2 ) c_nrho22 = nr
        if ( (c_reactant(nr,1) == c_nno2 .and. &
              c_reactant(nr,2) == c_no3) .or. &
             (c_reactant(nr,1) == c_no3 .and. &
              c_reactant(nr,2) == c_nno2) ) c_nro3no2 = nr
        if ( (c_reactant(nr,1) == c_no3 .and. &
              c_reactant(nr,2) == c_noh) .or. &
             (c_reactant(nr,1) == c_noh .and. &
              c_reactant(nr,2) == c_no3 ) ) c_nroho3 = nr
        if ( (c_reactant(nr,1) == c_nch4 .and. &
              c_reactant(nr,2) == c_noh) .or. &
             (c_reactant(nr,1) == c_noh .and. &
              c_reactant(nr,2) == c_nch4) ) c_nrohch4 = nr
      end do
      !
      ! TRACER ANALYSIS:
      ! LOOP FOR EACH INDIVIDUAL REACTION
      !
      do nr = 1 , c_nreac
        !
        !  noxchem(i,nr) = net rp O3 NOx PANs HNO3 RNO3;  NOx-PAN, PAN-NOX; HNO3
        !
        ! REACTANT1
        !
        if ( c_icat(c_reactant(nr,1)) == 11 .or. &
             c_icat(c_reactant(nr,1)) == 12 )  then
          c_noxchem(1,nr) = c_noxchem(1,nr) - d_one
        end if
        if ( c_icat(c_reactant(nr,1)) == 6 .or. &
            (c_icat(c_reactant(nr,1)) >= 12 .and. &
             c_icat(c_reactant(nr,1)) <= 15) )  then
          c_noxchem(2,nr) = c_noxchem(2,nr) - d_one
        end if
        if ( c_icat(c_reactant(nr,1)) ==  5 ) then
          c_noxchem(3,nr) = c_noxchem(3,nr) - d_one
        end if
        if ( c_icat(c_reactant(nr,1)) == 16 )  then
          c_noxchem(4,nr) = c_noxchem(4,nr) - d_one
        end if
        if ( c_icat(c_reactant(nr,1)) ==  7 ) then
          c_noxchem(5,nr) = c_noxchem(5,nr) - d_one
        end if
        if ( c_tchem(c_reactant(nr,1)) == '    R4N1' ) then
          c_noxchem(5,nr) = c_noxchem(5,nr) - d_one
        end if
        if ( c_tchem(c_reactant(nr,1)) == '    R3N1' ) then
          c_noxchem(5,nr) = c_noxchem(5,nr) - d_one
        end if
        if ( c_tchem(c_reactant(nr,1)) == '    R7N1' ) then
          c_noxchem(5,nr) = c_noxchem(5,nr) - d_one
        end if
        if ( c_tchem(c_reactant(nr,1)) == '    R6N1' ) then
          c_noxchem(5,nr) = c_noxchem(5,nr) - d_one
        end if
        if ( c_tchem(c_reactant(nr,1)) == '    R5N1' ) then
          c_noxchem(5,nr) = c_noxchem(5,nr) - d_one
        end if
        if ( c_tchem(c_reactant(nr,1)) == '    ISN1' ) then
          c_noxchem(5,nr) = c_noxchem(5,nr) - d_one
        end if
        if ( c_tchem(c_reactant(nr,1)) == '    ISNR' ) then
          c_noxchem(5,nr) = c_noxchem(5,nr) - d_one
        end if
        !
        ! REACTANT2
        !
        if ( c_reactant(nr,2) > 0 ) then
          if ( c_icat(c_reactant(nr,2)) == 11 .or. &
               c_icat(c_reactant(nr,2)) == 12 ) then
            c_noxchem(1,nr) = c_noxchem(1,nr) - d_one
          end if
          if ( c_icat(c_reactant(nr,2)) == 6 .or. &
              (c_icat(c_reactant(nr,2)) >= 12 .and. &
               c_icat(c_reactant(nr,2)) <= 15) ) then
            c_noxchem(2,nr) = c_noxchem(2,nr) - d_one
          end if
          if ( c_icat(c_reactant(nr,2)) == 5 ) then
            c_noxchem(3,nr) = c_noxchem(3,nr) - d_one
          end if
          if ( c_icat(c_reactant(nr,2)) == 16 ) then
            c_noxchem(4,nr) = c_noxchem(4,nr) - d_one
          end if
          if ( c_icat(c_reactant(nr,2)) == 7 ) then
            c_noxchem(5,nr) = c_noxchem(5,nr) - d_one
          end if
          if ( c_tchem(c_reactant(nr,2)) == '    R4N1' ) then
            c_noxchem(5,nr) = c_noxchem(5,nr) - d_one
          end if
          if ( c_tchem(c_reactant(nr,2)) == '    R3N1' ) then
            c_noxchem(5,nr) = c_noxchem(5,nr) - d_one
          end if
          if ( c_tchem(c_reactant(nr,2)) == '    R7N1' ) then
            c_noxchem(5,nr) = c_noxchem(5,nr) - d_one
          end if
          if ( c_tchem(c_reactant(nr,2)) == '    R6N1' ) then
            c_noxchem(5,nr) = c_noxchem(5,nr) - d_one
          end if
          if ( c_tchem(c_reactant(nr,2)) == '    R5N1' ) then
            c_noxchem(5,nr) = c_noxchem(5,nr) - d_one
          end if
          if ( c_tchem(c_reactant(nr,2)) == '    ISN1' ) then
            c_noxchem(5,nr) = c_noxchem(5,nr) - d_one
          end if
          if ( c_tchem(c_reactant(nr,2)) == '    ISNR' ) then
            c_noxchem(5,nr) = c_noxchem(5,nr) - d_one
          end if
        end if
        !
        ! PRODUCTS
        !
        do n = 1 , c_nnpro(nr)
          if ( c_product (nr,n) > 0 ) then
            if ( c_icat(c_product (nr,n)) == 11 .or. &
                 c_icat(c_product (nr,n)) == 12 ) then
              c_noxchem(1,nr) = c_noxchem(1,nr) + c_stoich(nr,n)
            end if
            if ( c_icat(c_product (nr,n)) == 6 .or. &
                (c_icat(c_product (nr,n)) >= 12 .and. &
                 c_icat(c_product (nr,n)) <= 15) ) then
              c_noxchem(2,nr) = c_noxchem(2,nr) + c_stoich(nr,n)
            end if
            if ( c_icat(c_product (nr,n)) ==  5 ) then
              c_noxchem(3,nr) = c_noxchem(3,nr) + c_stoich(nr,n)
            end if
            if ( c_icat(c_product (nr,n)) == 16 ) then
              c_noxchem(4,nr) = c_noxchem(4,nr) + c_stoich(nr,n)
            end if
            if ( c_icat(c_product (nr,n)) ==  7 ) then
              c_noxchem(5,nr) = c_noxchem(5,nr) + c_stoich(nr,n)
            end if
            if (c_tchem(c_product (nr,n)) == '    R4N1' ) then
              c_noxchem(5,nr) = c_noxchem(5,nr) + c_stoich(nr,n)
            end if
            if ( c_tchem(c_product (nr,n)) == '    R3N1' ) then
              c_noxchem(5,nr) = c_noxchem(5,nr) + c_stoich(nr,n)
            end if
            if ( c_tchem(c_product (nr,n)) == '    R7N1' ) then
              c_noxchem(5,nr) = c_noxchem(5,nr) + c_stoich(nr,n)
            end if
            if ( c_tchem(c_product (nr,n)) == '    R6N1' ) then
              c_noxchem(5,nr) = c_noxchem(5,nr) + c_stoich(nr,n)
            end if
            if ( c_tchem(c_product (nr,n)) == '    R5N1' ) then
              c_noxchem(5,nr) = c_noxchem(5,nr) + c_stoich(nr,n)
            end if
            if ( c_tchem(c_product (nr,n)) == '    ISN1' ) then
              c_noxchem(5,nr) = c_noxchem(5,nr) + c_stoich(nr,n)
            end if
            if (c_tchem(c_product (nr,n)) == '    ISNR' ) then
              c_noxchem(5,nr) = c_noxchem(5,nr) + c_stoich(nr,n)
            end if
          end if
        end do
        !
        ! ENTER STOICHIOMETRIES FOR
        !    6 and 7:  NOx=>PAN, PAN=>NOx
        !    8 and 9:  NOx=>HNO3, HNO3=>NOX
        !   10 and 11:  NOx=>RNO3;  RNO3=>NOx
        !   12 and 13:  NOx=>H+RNO3; H+RNO3=>NOx
        !
        if ( c_noxchem(3,nr) > 0 ) c_noxchem(6,nr) = c_noxchem(3,nr)
        if ( d_zero-c_noxchem(2,nr) < c_noxchem(6,nr) ) then
          c_noxchem(6,nr) = 0.-c_noxchem(2,nr)
        end if
        if ( c_noxchem(6,nr) < d_zero ) c_noxchem(6,nr) = d_zero
        if ( c_noxchem(3,nr) < d_zero ) then
          c_noxchem(7,nr) = d_zero-c_noxchem(3,nr)
        end if
        if ( c_noxchem(2,nr) < c_noxchem(7,nr) ) then
          c_noxchem(7,nr) =    c_noxchem(2,nr)
        end if
       if ( c_noxchem(7,nr) < d_zero ) c_noxchem(7,nr) = d_zero
       if ( c_noxchem(4,nr) > d_zero ) c_noxchem(8,nr) = c_noxchem(4,nr)
       if ( d_zero-c_noxchem(2,nr) < c_noxchem(8,nr) ) then
         c_noxchem(8,nr) = d_zero-c_noxchem(2,nr)
       end if
       if ( c_noxchem(8,nr) < d_zero ) c_noxchem(8,nr) = d_zero
       if ( c_noxchem(4,nr) < 0 ) then
         c_noxchem(9,nr) = d_zero - c_noxchem(4,nr)
       end if
       if ( c_noxchem(2,nr) < c_noxchem(9,nr) ) then
         c_noxchem(9,nr) = c_noxchem(2,nr)
       end if
       if ( c_noxchem(9,nr) < d_zero ) c_noxchem(9,nr) = d_zero
       if ( c_noxchem(5,nr) > d_zero ) c_noxchem(10,nr) = c_noxchem(5,nr)
       if ( d_zero-c_noxchem(2,nr) < c_noxchem(10,nr) ) then
         c_noxchem(10,nr) = d_zero-c_noxchem(2,nr)
       end if
       if ( c_noxchem(10,nr) < d_zero ) c_noxchem(10,nr) = d_zero
       if ( c_noxchem(5,nr) < d_zero ) then
         c_noxchem(11,nr) = d_zero - c_noxchem(5,nr)
       end if
       if ( c_noxchem(2,nr) < c_noxchem(11,nr) ) then
         c_noxchem(11,nr) = c_noxchem(2,nr)
       end if
       if ( c_noxchem(11,nr) < d_zero ) c_noxchem(11,nr) = d_zero
       if ( (c_noxchem(4,nr) + c_noxchem(5,nr) > d_zero) ) then
         c_noxchem(12,nr) = c_noxchem(4,nr) + c_noxchem(5,nr)
       end if
       if ( d_zero-c_noxchem(2,nr) < c_noxchem(12,nr) ) then
         c_noxchem(12,nr) = d_zero-c_noxchem(2,nr)
       end if
       if ( c_noxchem(12,nr) < d_zero ) c_noxchem(12,nr) = d_zero
       if ( (c_noxchem(4,nr)+c_noxchem(5,nr) < 0) ) then
         c_noxchem(13,nr) = d_zero - c_noxchem(4,nr) - c_noxchem(5,nr)
       end if
       if ( c_noxchem(2,nr) < c_noxchem(13,nr) ) then
         c_noxchem(13,nr) = c_noxchem(2,nr)
       end if
       if ( c_noxchem(13,nr) < d_zero ) c_noxchem(13,nr) = d_zero
     end do
!
  355 format(//,'MAJOR ERROR:  REACTION NOT IDENTIFIED WITH SPECIES', &
             //,' REACTION NUMBER =',i4,/,                            &
                ' FIRST REACTANT =', i4,'  CATEGORY =',i3,/,          &
                ' SECOND REACTANT=', i4,'  CATEGORY =',i3,/,          &
                ' PRODUCTS =', (10i5)  )
 1209 format(' SET PAIRFAC: species A + B <->C.', &
             ' SPECIES AND PAIRFACS = ',/,3(a8,2x),3(1pe10.3))
 1208 format(' SET MULTIFAC:  species A + B <->C.',&
             ' SPECIES AND MULTIFACS = ',/,3(a8,2x),3(1pe10.3))
 1206 format(/,' MAJOR ERROR IN QUADINIT: ', &
               'EXCHANGED SPECIES INDEX EXCEEDED LIMIT' ,' (i=1,20, 1205)',/)
 1212 format (/,'EXLOSS(ICX,IS,I)=NRX. ','ICX,ICP,IS,I,NRX=',5i5,/, &
                ' REACTANTS=',2i5, '  PRODUCTS=',2i5   )
 1213 format(/,' MAJOR ERROR IN QUADINIT: ','EXLOSS INDEX EXCEEDED LIMIT' , &
               ' (ii=1,5, 1213)',/)
 1222 format (/,' EXPRO(ICX,IS,I)=NRX.  ','ICX,ICP,IS,I,NRX=', 5i5,/, &
                ' REACTANTS=',2i5,'  PRODUCTS=',2i5   )
 1224 format(/,'WARNING: POSSIBLE EXCHANGE REACTION '&
               ,'AMONG UNPAIRED SPECIES: nr =',       &
               /, i5, a8,'+',a8,'=>',a8,'+',a8,'+', a8)
 1223 format(/,' MAJOR ERROR IN QUADINIT: ', 'EXPRO INDEX EXCEEDED LIMIT' , &
               ' (ii=1,5, 1223)',/)

    end subroutine cheminit
!
! ---------------------------------------------
!
! This writes chemical concentrations and reaction rates (gas+aqueous)
!   for the specified vector loop (kw).
!
! Written to the output file (lout)
!
! Called by:  boxmain (as part of output at specified time steps)
! Calls to:   None.
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------

    subroutine chemwrit(kw)
      implicit none
      integer(ik4) , intent(in) :: kw

      ! Chem index
      integer(ik4) :: ic
      ! Aqueous counters
      integer(ik4) :: neq
      ! Chem species counters
      integer(ik4) :: nc , nc1 , nc2 , ncf , nn
      ! Reaction counters
      integer(ik4) :: nr , nrh , nrq , nrqq
      ! Vectorization counters
      integer(ik4) :: kk
      ! General counters
      integer(ik4) :: i , j

      real(rk8) :: calpha(c_kvec) ! General vector variable
      real(rk8) :: cbeta(c_kvec)  ! General vector variable
      real(rk8) :: cgamma(c_kvec) ! General vector variable
      ! 'SUM' name for output
      character(len=8) :: tsum
      ! Gas phase concentration
      real(rk8) :: xcgas
      ! Aqueous   concentration
      real(rk8) :: acquacon
!
      if ( kw <= 0 ) return
      kk = 1
      !
      ! SPECIES CONCENTRATIONS.
      !
      write(c_out,1141) c_hour
      ncf = idint(0.2D0*(dble(c_nchem2)-0.01D0))+1
      do nc = 1 , ncf
        nc1 = 4*(nc-1) + 1
        nc2 = nc1 + 3
        if ( nc1 > c_nchem2 ) nc1 = c_nchem2
        if ( nc2 > c_nchem2 ) nc2 = c_nchem2
        write(c_out,1142) (c_tchem(nn),c_xcout(kw,(nn)),nn=nc1,nc2)
      end do
      !
      !   AQUEOUS SPECIES CONCENTRATIONS
      !
      ! NOTE  c_xcout(kw,ic) for GAS SPECIES IS GAS+AQ SUM, molec/cm3 equiv.
      !       c_xcout(kw,icq) for AQUEOUS SPECIES is moles/liter
      !    .  GAS-ONLY IS CALCULATED HERE.
      ! NOTE!  WRITE error for lumped species such as  MP (in rooh);
      !        c_xcout(kw,ic) is preserved as species fraction, not sum?
      !
      if ( c_h2oliq(kw ) > d_zero ) then
        acquacon = c_h2oliq(kw)*avogadrl
        write(c_out,1301) c_h2oliq(kw )
        write(c_out,1303) acquacon
        if ( c_aqueous(1,2) > 0 ) then
          write(c_out,1302) c_xcout(kw ,c_aqueous(1,2))
        end if
        tsum = '     SUM'
!       do ic = 1 , c_nchem2
        do nrh = 1 , c_nreach
          ic = c_henry(nrh,1)
          if ( c_nequil(ic) > 0 )  then
            xcgas = c_xcout(kw,ic)
            do i = 1 , c_nequil(ic)
              xcgas = xcgas - c_xcout(kw,c_ncequil(ic,i))*acquacon
            end do
            write(c_out,1305) c_tchem(ic), xcgas, (c_tchem(c_ncequil(ic,i)), &
                  c_xcout(kw,c_ncequil(ic,i)),i=1,c_nequil(ic)), &
                  tsum, c_xcout(kw,ic)
            if ( rhdif(kw,nrh) > 0 ) then
              xcgas = (c_xcout(kw,c_ncequil(ic,1))*acquacon)/rhdif(kw,nrh)
            end if
            write(c_out,1306)xcgas, rateh(kw ,nrh), rhdif(kw ,nrh)
          end if
        end do
      end if
      !
      ! --------------
      ! CHEMWRIT  SUMMARY WRITE:
      !    CONCENTRATIONS, XCIN, RP, RL from QUADCHEM
      ! --------------
      write(c_out,1801) c_iter, kw
      write(c_out,1804) c_idate, c_hour, c_lat(1), c_lon(1), c_temp(1)
      write(c_out,1805) (c_jparam(i),i=1,13),c_jparam(20)
      write(c_out,1802)
      !
      ! WRITE GAS-AQUEOUS SUMS
      !   MAKE XC,  RP, RL=GAS+AQ SUM for gas species.
      !    CONVERT AQUEOUS XC TO GAS UNITS
      !     (Note, sum aqueous into gas for xc, but not for xcout.
      !      xcout is made gas-aqueous sum in postlump)
      !
      do j = 1 , c_nchem2
        if ( c_npequil(j) == j ) then
          xrr(kw,1) = d_zero
          rpro(kw,1) = d_zero
          rloss(kw,1) = d_zero
          do neq = 1 , (c_nequil(j)+1)
            ic = j
            if ( neq > 1 ) ic = c_ncequil(j,(neq-1))
            calpha(kw) = d_one
            if ( neq > 1 ) calpha(1) = c_h2oliq(1)*avogadrl
            rloss(kw,1) = rloss(kw,1) + c_rl(kw,ic)
            rpro(kw,1) = rpro(kw,1) + c_rp(kw,ic)
            !
            ! GAS-AQUEOUS SUM + CONVERSION
            !
!           xrr(kw,1) = xrr(kw,1) + c_xcout(kw,ic)*calpha( kw)
            !
            ! AQUEOUS CONVERSION, NO GAS-AQ SUM
            !
            xrr(kw,1) = c_xcout(kw,ic)*calpha(kw)
          end do
          !
          ! CORRECTION- XR IS GAS-AQ SUM
          !
          xrr(kw,1) = c_xcout(kw,j)
          cbeta(1) = xrr(kw,1) - c_xcin(kw,j)
          cgamma(1) = rpro(kw,1) - rloss(kw,1)
          write(c_out,1803) j, c_tchem(j), xrr(kw,1), c_xcin(kw,j), &
                   xcfinr(kw,j), rloss( kw,1), rpro( kw,1), cbeta(1), cgamma(1)
        end if
      end do
      !
      ! ------
      ! END SUMMARY WRITE OPTION
      ! ------
      !
      ! GAS-PHASE REACTION RATES
      !
      write(c_out,1143)
      do nr = 1 , c_nreac
        write(c_out,1144) nr, (c_treac(j,nr),j=1,5), &
                     c_rr(kw ,nr), ratek(kw ,nr)
      end do

      write(c_out,1146)
      if ( c_nreach > 0 ) then
        do i = 1 , c_nreach
          write(c_out,1147) (c_treach(j,i),j=1,2), rateh(kw ,i), rhdif(kw ,i)
        end do
      end if
      !
      ! AQUEOUS EQUILIBRIUM CONSTANTS
      !
      write(c_out,1153)
      do nrq = 1 , c_nreacq
        write(c_out,1154) nrq, (c_treacq(j,nrq),j=1,3), rateq(kw ,nrq)
      end do
      !
      ! SPECIAL EQUILIBRIUM CONSTANTS
      !
      write(c_out,1156)
      do nrqq = 1 , c_nreaqq
        write(c_out,1154) nrqq, (c_treaqq(j,nrqq),j=1,3), rateqq(kw,nrqq)
      end do
!
 1141 format(/,' SPECIES CONCENTRATIONS AT  XHOUR =',f6.2)
 1142 format(4(a8,1pe10.3,2x))
 1301 format(/,'AQUEOUS CHEMISTRY ',/,'WATER (grams/cm3) =',1pe10.3)
 1302 format( ' [H+] (moles per liter) =',1pe10.3)
 1303 format( ' GAS SPECIES and G-A SUM molec/cm3. ',  &
              ' AQUEOUS moles/liter. CONVERSION= ', 1pe10.3)
 1305 format(4(a8,1pe10.3,2x),/,'         0.000E+00  ',3(a8,1pe10.3,2x))
 1306 format('   (',1pe10.3,')                      HENRY, Hw/DIFF=',3(1pe10.3))
 1801 format(//,' SUMMARY WRITE:  ITER =',i3, '    VECTOR KKW=',i3)
 1802 format(/,' #  IC     XCOUT     XCIN  XCF/AV       RL',          &
               '        RP        dXC      dR')
 1803 format(i4,a8,2(1pe10.3),0pf7.3,1x,(4(1pe10.3)))
 1804 format('IDATE xhour lat lon temp=',i8,4f8.2)
 1805 format('JPARAMS: zenith altmid dobson so2 no2 aaerx aerssa',      &
             ' albedo',/,'cld-below cld-above claltB claltA temp date', &
             /,(8f10.3))
 1143 format(/' GAS-PHASE REACTION RATES')
 1144 format(i4,2x,a8,'+',a8,'=>',a8,'+',a8,'+',a8,2(1pe10.3))
 1146 format(/,'HENRYS LAW CONSTANT + MODIFIED CONSTANT', &
               ' WITH DIFFUSION ADJUSTMENT:')
 1147 format(a8,'=',a8,2x,2(1pe10.3))
 1153 format(/,' AQUEOUS EQUILIBRIUM CONSTANTS (H+)')
 1154 format(i4,2x,a8,'<=>',a8,'+',a8,2(1pe10.3))
 1156 format(/,' SPECIAL EQUILIBRIUM CONSTANTS ')

     end subroutine chemwrit
!
! --------------------------------------------------------------
!
! This generates and writes a summary of chemistry
!  for the specified species (TITL) and vector dimension (KW).
!
! Output (written to lout) includes:
!   list of all reactions that produce/remove the species,
!   production and loss rates (molec/cm3 per time step) ,
!   and change in concentration (molec/cm3)
!
! Gas/aqueous treatment:
!  Assumes that the species concentration XR
!   represents the gas+aqueous sum
!   (i.e. called after 'postlump' does the gas-aqueous sum)
!
! Called by:  boxmain (as part of output at specific times)
!               (used to explore chem. of specific species)
! Calls to:   None.
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
!
    subroutine analyze(titl, kw)
!
      implicit none

      ! Name of specified chem species
      character(len=8) , intent(in) :: titl
      integer(ik4) , intent(in) :: kw

      ! Chem index
      integer(ik4) :: ic , icc , ics
      ! Aqueous counters
      integer(ik4) :: neq
      ! Reaction counters
      integer(ik4) :: nr , nrh , nrq
      ! Vectorization counters
      integer(ik4) :: kk
      ! General counters
      integer(ik4) :: i , n

      ! Gas phase concentration
      real(rk8) :: xcgas
      ! Aqueous   concentration
      real(rk8) :: acquacon
      ! Dimensionless Henry coefficient
      real(rk8) :: xcoeff
      ! Droplet diffusion factor
      real(rk8) :: xcoeff2
      ! Chem. reaction rate molec/cm3/step
      real(rk8) :: tpro
      ! Species production  molec/cm3/step
      real(rk8) :: xpro
      ! Stoichiometry for spec. production
      real(rk8) :: stopro
!
      ! Chem. reaction loss molec/cm3/step
      real(rk8) :: tloss
      ! Species loss rate   molec/cm3/step
      real(rk8) :: xloss
      ! Stoichiometry for species loss
      real(rk8) :: stoloss
!
      ! Net production minus loss /cm3/step
      real(rk8) :: tnetpro
      ! Change in species conc molec/cm3
      real(rk8) :: tdelta
!
      ! counter for  aqueous  spec
      integer(ik4) :: neq1
!
      if ( kw <= 0 ) return
      kk = 1
      !
      ics = namechem(titl)
      if ( ics <= 0 ) then
        write(c_out,21) titl, ics
        return
      end if
      !
      ! LOOP FOR ALL SPECIES LINKED THROUGH CHEMICAL PAIRS  - CUT
      !
      ic = ics
      write(c_out,22) titl, ic
      if ( c_nequil(ic) > 0 ) then
        write(c_out,23) c_h2oliq(kw)
      end if
      !
      !  IF ACQUA>0, WRITE ASSOCIATED HENRY'S LAW AND EQUILIBRIA.
      !
      acquacon = c_h2oliq(kw)*avogadrl
      if ( c_nequil(ic) > 0 .and. c_h2oliq(kw) > d_zero ) then
        write(c_out,31) acquacon
      end if
      !
      !  IF ACQUA>0, LOOP FOR AQUEOUS PHASE SPECIES
      !
      if ( c_nequil(ic) > 0 .and. c_h2oliq( kw) > d_zero ) then
        !
        ! CONCENTRATION OF GAS-SPECIES ONLY
        !
        xcgas = c_xcout( kw,ic)
        do n = 1 , c_nequil(ic)
          icc = c_ncequil(ic,n)
          xcgas = xcgas - c_xcout( kw,icc)*acquacon
          write(c_out,10031) icc, c_tchem(icc),c_xcout(kw,icc),  &
                            acquacon, xcgas
        end do
        write(c_out,32) c_xcout( kw,ic), xcgas
        !
        ! AQUEOUS CONCENTRATIONS AND CORRESPONDING HENRY'S LAW COEFFICIENTS
        ! OR EQUILIBRIUM CONSTANTS
        do n = 1 , c_nequil(ic)
          icc = c_ncequil(ic,n)
          nrq = c_nrequil(ic,n)
          nrh = c_nrequil(ic,1)
          if ( n == 1 ) then
            xcoeff = rateh( kw,nrh)*atmos/acquacon
            xcoeff2 = d_zero
            if ( rateh(kw,nrh) > d_zero ) then
              xcoeff2 = rhdif(kw,nrh)/rateh(kw,nrh)
              write(c_out,51) c_tchem(icc),c_xcout( kw,icc)
              write(c_out,52) (c_treach(i,nrh),i=1,2), xcoeff, xcoeff2
            else
!             if ( xcoeff > 0 ) then
!               xcoeff = rateq(kw,nrq)/xcoeff
!             else
!               xcoeff = rateq(kw,nrq)
!             end if
              xcoeff = rateq( kw,nrq)
              write(c_out,51) c_tchem(icc),c_xcout( kw,icc)
              write(c_out,53) (c_treacq(i,nrq),i=1,3),xcoeff
            end if
          end if
        end do
        !
        ! END IF FOR AQUEOUS SPECIES
        !
      end if
      !
      !  SUM AND WRITE CHEMICAL PRODUCTION  (RATE MOL CM-3 PER TIME STEP)
      !
      tpro = d_zero
      write(c_out,201)
      !
      ! DO FOR ALL SPECIES LINKED THROUGH CHEMICAL PAIRS - CUT.
      !
      ic = ics
      !
      ! DO FOR GAS AND ALL ACQUEOUS EQUIVALENT SPECIES
      !
      do neq1 = 1 , (c_nequil(ic)+1)
        if ( neq1 == 1 ) then
          acquacon = d_one
          icc = ic
        else
          acquacon = c_h2oliq( kw)*avogadrl
          neq = neq1-1
          icc = c_ncequil(ic,neq)
        end if
        if ( acquacon == 0 ) cycle
        !
        ! LOOP THROUGH ALL REACTIONS TO FIND LOSSES FOR SPECIES
        ! NOTE THAT RR, REACTION RATE, IS ALREADY IN GAS UNITS FOR  AQ SPECIES.
        !
        do nr = 1 , c_nreac
          stopro  = d_zero
          if (c_reactant(nr,1) == icc) stopro = stopro - d_one
          if (c_reactant(nr,2) == icc) stopro = stopro - d_one
          do i = 1 , 20
            if ( c_product(nr,i) == icc ) stopro = stopro + c_stoich(nr,i)
          end do
          if ( stopro  > d_zero ) then
            xpro = c_rr(kw,nr) * stopro
            tpro = tpro + xpro
            write(c_out,102) nr, (c_treac(i,nr),i=1,5), &
                  xpro, c_rr( kw,nr) , ratek( kw,nr)
          end if
        end do
      end do
      !
      ! SUM AND WRITE CHEMICAL LOSSES
      !
      tloss = d_zero
      write(c_out,101)
      !
      ! DO FOR ALL SPECIES LINKED THROUGH CHEMICAL PAIRS - CUT
      !
      ic = ics
      !
      ! DO FOR GAS AND ALL ACQUEOUS EQUIVALENT SPECIES
      !
      do neq1 = 1 , (c_nequil(ic)+1)
        if ( neq1 == 1 ) then
          acquacon = d_one
          icc = ic
        else
          acquacon = c_h2oliq(kw)*avogadrl
          neq = neq1-1
          icc = c_ncequil(ic,neq)
        end if
        if ( acquacon == 0 ) cycle
        !
        ! LOOP THROUGH ALL REACTIONS TO FIND LOSSES FOR SPECIES
        !
        do nr = 1 , c_nreac
          stoloss = d_zero
          if ( c_reactant(nr,1) == icc ) stoloss = stoloss + d_one
          if ( c_reactant(nr,2) == icc ) stoloss = stoloss + d_one
          do i = 1 , 20
            if ( c_product(nr,i) == icc ) stoloss = stoloss - c_stoich(nr,i)
          end do
          if ( stoloss > d_zero ) then
            xloss = c_rr(kw,nr) * stoloss
            tloss = tloss + xloss
            write(c_out,102) nr,(c_treac(i,nr),i=1,5), &
                      xloss, c_rr( kw,nr) , ratek( kw,nr)
          end if
        end do
      end do
      !
      ! WRITE FINAL SUMMARY
      !
      tnetpro = tpro-tloss
      write(c_out,301) titl, tpro, tloss, tnetpro
      !
      ! FOR ALL SPECIES LINKED THROUGH CHEMICAL PAIRS - CUT
      !
      ic = ics
      tdelta = c_xcout(kw,ic) - c_xcin(kw,ic)
      write(c_out,302) c_tchem(ic), c_xcout(kw,ic), c_xcin(kw,ic), tdelta,  &
                       c_rp(kw,ic), c_rl(kw,ic)

   21 format(/,'WARNING:  UNKNOWN CHEMICAL NAME IN SUBROUTINE', &
              'ANALYZE',/,'NAME = ',a8,'  IC=',i5)
   22 format(/,'CHEMICAL PRODUCTION AND LOSS ANALYSIS FOR:',/,  &
               'SPECIES = ',a8,'  IC=',i5)
   23 format(/,'LIQUID WATER (gr/cm3) =', 1pe10.3)
   31 format(/,'AQUATIC CONVERSION FACTOR:', &
               ' MOLE/LITER per MOLEC/CM3 =', 1pe10.3)
   32 format(/,'  TOTAL SPECIES CONCENTRATION =',1pe10.3,/,  &
               '    GAS SPECIES CONCENTRATION =',1pe10.3)
10031 format(' TEST ICC TCHEM XR*AQUACON subtr fr XGSUM',/,  &
             i5,2x,a8,2x,8(1pe10.3))
   51 format(/,'AQUEOUS SPECIES = ',a8,' MOLES/LITER=',1pe10.3)
   52 format(a8,'=',a8,'  HENRYs LAW COEF=',1pe10.3,  &
             '  DROPLET DIFF FAC=', 1pe10.3)
   53 format(a8,'=',a8,'+',a8,'  EQUILIBRIUM CONSTANT= ',1pe10.3)
  101 format(/,'PHOTOCHEMICAL LOSSES:',/,           &
               '         REACTION                ', &
               ' LOSS        RATE        RATEK')
  102 format(i4,1x,a8,'+',a8,'=>',a8,'+',a8,'+',a8,3(1pe10.3   ))
  201 format(/,'PHOTOCHEMICAL PRODUCTION:',/,          &
               '         REACTION                ',    &
               ' PRODUCTION  RATE        RATEK     ')
  301 format(/,'SPECIES = ',a8,/,                            &
               ' SUMMED PHOTOCHEMICAL PRODUCTION = ',1pe10.3,/,  &
               ' SUMMED PHOTOCHEMICAL LOSS       = ',1pe10.3,/,  &
               ' NET    PHOTOCHEMICAL PRODUCTION = ',1pe10.3)
  302 format(/,' SPEC GAS+AQ CONC. PRIOR CONC. NET CHANGE   ', &
               '(INTERNAL RPRO     RLOSS)',/,a8,3(1pe12.4),2x,2(1pe12.4))
    end subroutine analyze

end module mod_cbmz_init1
