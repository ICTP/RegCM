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

module mod_cbmz_rates1
!
  use m_realkinds
  use mod_constants
  use mod_cbmz_chemvars
  use mod_cbmz_chemmech
  use mod_cbmz_chemlocal
  use mod_cbmz_jval1
!
  private
!
  public :: chemrates , hvrates

  contains

!
!  chemrates.f   April, 2007
!
!  Programs to calculate chemical reaction rates and hv rates
!   for the RADICAL BALANCE-BACK EULER solution for photochemistry.
!
!  Note alternative versions: YYMMDD (standard) or YYYYDDD in ZEN1
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
!  Subroutines in this file:
!    chemrates: sets rate constants for non-photolysis reactions
!      bod:  function to set termolecular rate constants
!      ytn:  function to set alkyl nitrate yield from RO2+NO reactions.
!    hvrates:   sets rate constants for photolysis reactions
!      zen1:      sets solar zenith angle
!
!  FUTURE ADD:  Identified below by 'FUTURE ADD'.
!     Exponential decay solution option (done)
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
!
! This establishes reaction rate constants
!  as a function of input rate parameters, temperature and density.
!
! Final rate constants, temperature and density are vectorized.
!
! Rates are based on k= A*exp(-B/temp) with A, B from REACTION.DAT.
! Special formats are included for special rate constants (3-body, etc.)
! ALK4, ALK7 stoichometries have been converted to special rates.
!
! ---------------------------------------------------------
! NOTE:  HVRATE CONSTANTS ARE MOVED TO SUBROUTINE BHVRATES.
! ---------------------------------------------------------
!
! Inputs:
!   chemistry reaction parameters (RK)
!   temperature
!   density
!
! Outputs:  rate constants (RATEK)
!
! Called by:  either boxmain (as part of chem. solution, or for info)
!             or by quadchem (as part of chemistry solution)
!
! Calls to:   None.
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
!
    subroutine chemrates
!
      implicit none

      ! temperature K
      real(dp) :: tempx
      ! density, molec/cm3
      real(dp) :: denx
      ! Function: 3-bodyreaction rate cm3-sec
      real(dp) :: bod
      ! Function: RNO3 yield from RO2+NO
      real(dp) :: ytn

      kk = 1
      !
      !  SET ALL RATE CONSTANTS (including initial index values for HV rates,)
      !
      ! --------------------------
      ! LOOP TO SET RATE CONSTANTS
      ! --------------------------
      !
      ! counter for parameterized ro2-ro2 reactions
      !
      nr1 = 0
      do j = 1 , c_nreac
        !
        ! Default
        !
        ratek(kk, j) = c_rk(1,j)
        !
        ! Standard 2-body reaction format
        !
        if ( c_nrk(j) == 0 ) then
          ratek(kk,j) = c_rk(1,j)*dexp(-c_rk(2,j)/c_temp(kk))
        end if
        !
        ! Standard 2-body reaction format with input rate at 298 K.
        !
        if ( c_nrk(j) == -1 ) then
          ratek(kk,j) = c_rk(1,j)*dexp(-c_rk(2,j) * &
                     (d_one/c_temp(kk)-0.0033557047D0))
        end if
        !
        ! Standard format *T/300
        !
        if ( c_nrk(j) == -2 ) then
          ratek(kk,j) = c_rk(1,j)*dexp(-c_rk(2,j)/c_temp(kk)) * &
                        (c_temp(kk)/300.0D0)**c_rk(3,j)
        end if
        !
        ! Standard 3-body reaction format - bod function
        !
        if ( c_nrk(j) == -3 ) then
          tempx = c_temp(kk)
          denx = c_dens(kk)
          ratek(kk,j) = bod(c_rk(1,j),c_rk(2,j),c_rk(3,j), &
                            c_rk(4,j),c_rk(5,j),tempx,denx )
        end if
        !
        ! Combined 2-body and  3-body reaction format
        !
        if ( c_nrk(j) == -4 ) then
          tempx = c_temp(kk)
          denx = c_dens(kk)
          ratek(kk,j) = bod(c_rk(3,j),c_rk(4,j),c_rk(5,j), &
                            c_rk(6,j),c_rk(7,j),tempx, denx)
          ratek(kk,j)= ratek(kk,j)*c_rk(1,j)*dexp(-c_rk(2,j)/tempx)
        end if
        !
        ! 2-body reaction format with nitrate yield subtracted
        !
        if ( c_nrk(j) == -5 ) then
          ratek(kk,j) = c_rk(1,j)*dexp(-c_rk(2,j)/c_temp(kk))
          ratek(kk,j) = ratek(kk,j) * &
            (d_one-ytn(c_rk(3,j),c_temp(kk),c_dens(kk)))
        end if
        !
        ! 2-body reaction format for  nitrate yield
        !
        if ( c_nrk(j) == -6 ) then
          ratek(kk,j) = c_rk(1,j)*dexp(-c_rk(2,j)/c_temp(kk))
          ratek(kk,j) = ratek(kk,j)*(ytn(c_rk(3,j),c_temp(kk),c_dens(kk)))
        end if
        !
        ! Density-dependent format A exp(-B/T)* (C + D*dens)/(E + F*dens)
        !    for CO, GLYX+OH
        !
        if ( c_nrk(j) == -7 ) then
          ratek(kk,j) = c_rk(1,j)*dexp(-c_rk(2,j)/c_temp(kk))
          ratek(kk,j) = ratek(kk,j)*(c_rk(3,j)+c_rk(4,j)*c_dens(kk)) / &
                                    (c_rk(5,j)+c_rk(6,j)*c_dens(kk))
        end if
        !
        ! Combined format, interacting two-body reactions with density
        !   (for MO2+MO2 and DMS)
        ! A exp(-B/T) / (1 + (C+E*dens)* exp(-D/T))
        !
        if ( c_nrk(j) == -8 ) then
          ratek(kk,j) = c_rk(1,j)*dexp(-c_rk(2,j)/c_temp(kk))
          ratek(kk,j) = ratek(kk,j) / &
            (d_one+(c_rk(3,j)+c_rk(5,j)*c_dens(kk)) * &
            dexp(-c_rk(4,j)/c_temp(kk)))
        end if
        !
        ! Special HO2+HO2 rate:
        !   (Aexp(-B/t) + C*dens*exp(-D/t)) *(1+ E exp(-F/T)
        !  from Chameides around 1996
        !
        if ( c_nrk(j) == -9 ) then
          ratek(kk,j) = c_rk(1,j)*dexp(-c_rk(2,j)/c_temp(kk))
          ratek(kk,j) = ratek(kk,j) + &
                        c_rk(3,j)*dexp(-c_rk(4,j)/c_temp(kk)) * c_dens(kk)
          ratek(kk,j) = ratek(kk,j) * &
            (d_one+c_rk(3,j)*dexp(-c_rk(4,j)/c_temp(kk))*xc(kk,14))
        end if
        !
        ! Special HNO3+OH rate:
        !   A exp(-B/t) + C exp(-D/t)*E exp(-F/t)/(Cexp(-D/T)+Eexp(-F/T))
        !  from Evans, Fiore 2003
        !
        if ( c_nrk(j) == -10 ) then
          ratek(kk,j) = c_rk(1,j)*dexp(-c_rk(2,j)/c_temp(kk))
          cbeta(kk) = c_rk(3,j)*dexp(-c_rk(4,j)/c_temp(kk))*c_dens(kk)
          cgamma(kk) = c_rk(5,j)*dexp(-c_rk(6,j)/c_temp(kk))
          ratek(kk,j) = ratek(kk,j) + & 
            cbeta(kk)*cgamma(kk)/(cbeta(kk)+cgamma(kk))
        end if
        !
        ! Special C3H8 format (IUPAC 2002):
        !        Aexp(-B/t)/(1+C*T/300^(D*exp(-E/t))
        !
        if ( c_nrk(j) == -11 ) then
          ratek(kk,j) = c_rk(1,j)*dexp(-c_rk(2,j)/c_temp(kk))
          cbeta(kk) = c_rk(4,j)*dexp(-c_rk(5,j)/c_temp(kk))
          ratek(kk,j) = ratek(kk,j) / &
            (d_one+c_rk(3,j)*(c_temp(kk)/300.0D0)**cbeta(kk))
        end if
        !
        ! Sum of two reactions (for ACET+OH, JPL 2002):
        !       A exp(-B/T) + C exp(-D/
        !
        if ( c_nrk(j) == -12 ) then
          ratek(kk,j) = c_rk(1,j)*dexp(-c_rk(2,j)/c_temp(kk)) + &
                        c_rk(3,j)*dexp(-c_rk(4,j)/c_temp(kk))
        end if
        !
        ! Parameterized RO2-RO2 reactions
        !      (CBMZ format: RO2->partial product from RO2-RO2 reaction,
        !        see setro2 in chemsolve.f for exact information.)
        !
        if ( c_nrk(j) == -13 ) then
          nr1 = nr1+1
          ratek(kk,j) = c_rk(1,j)*dexp(-c_rk(2,j)/c_temp(kk))
          ratero2(kk,nr1,1) = ratek(kk,j)
          ratero2(kk,nr1,2) = (c_rk(3,j)*dexp(-c_rk(4,j)/c_temp(kk)))**0.5D0
        end if
        !
        ! REACTION ON SULFATE AEROSOL SURFACES.  c_nrk(j) = -20
        !   ratek = saersa*(cgamma/4)* [8RT/(pi*MW)]^0.5 (Pandis and Seinfeld)
        !   for saersa, aerosol surface area in cm2/cm3
        !   cgamma = reaction efficiency (input)
        !   MW = molecular weight (input)
        !
        if ( c_nrk(j) == -20 ) then
          ratek(kk,j) = c_saersa(kk) * (0.25D0*c_rk(1,j)) * &
             ((8.0D0*rumolec*c_temp(kk))/(mathpi*c_rk(2,j)))**0.5D0
        end if
        !
        ! END SET ALL RATE CONSTANTS LOOP
        !
      end do
      ! --------------------------
      ! END - LOOP TO SET RATE CONSTANTS
      ! --------------------------
      ! ---------------------------
      ! HENRY'S LAW AND AQUEOUS RATE CONSTANTS
      ! ---------------------------
      ! HENRY'S LAW CONSTANT:
      ! UNIT CONVERSION FROM MOLES/LITER-ATMOS TO MOLES/LITER-molecules/cm2
      ! AND COMBINED WITH THE CONVERSION FACTOR FOR MOLES/LITER TO
      !  MOLECULES PER CM2 (for use in AQUASOLVE).
      ! (NOTE:  Dimensionless Henry's law H' = HRT,  R=8.314E-02. L-atm/mol-K.
      !   BARTH VALUE:  0.08206.
      !   H' = HRT = Caq/Cg*L;  H" = HRTL.  L=acqua.  RT=avogadrl/atmos. OK).
      !
      ! (with LIQ WATER 0.3e-6 gr/cm;  H = H'*1.248e5, H'=H*8E-6.)
      !
      !  FUTURE CHANGE: ADD ALT. SPECIAL FORMULAS FOR HENRYS LAW CONSTANTS
      !   AS NEEDED, CONTROLLED BY   c_nrkh()
      do j = 1 , c_nreach
        rateh(kk,j) = c_rkh(1,j)*dexp(d_zero-c_rkh(2,j) * &
                                      (d_one/c_temp(kk)-d_one/298.0D0))
        rateh(kk,j) = rateh(kk,j)*rtcon*c_temp(kk)*c_h2oliq(kk)
        !
        ! (OLDER CONVERSION)
        !        rateh(kk ,j) = rateh(kk, j)
        !    *       *c_h2oliq(kk) * avogadrl/atmos
        !
      end do
      !
      ! H+ AQUEOUS EQUILIBRIUM COEFFICIENTS (MOLES/LITER)
      !
      ! FUTURE CHANGE: ADD ALT. SPECIAL FORMULAS FOR EQUILIBRIUM COEFFICIENTS
      ! AS NEEDED, CONTROLLED BY   c_nrkq()
      !
      do j = 1 , c_nreacq
        rateq(kk,j) = c_rkq(1,j)*dexp(d_zero-c_rkq(2,j) * &
                                      (d_one/c_temp(kk)-d_one/298.0D0))
      end do
      !
      ! (CONVERSION:  AQUEOUS K1 IS MULTIPLIED BY HENRY'S LAW CONSTANT
      !   AQUEOUS K2 IS MULTIPLIED BY K1KH, etc.
      ! DELETED )
      !  SPECIAL AQUEOUS EQUILIBRIUM COEFFICIENTS (MOLES/LITER)
      !  FUTURE CHANGE: ADD ALT. SPECIAL FORMULAS FOR SPECIAL EQ  COEFFICIENTS
      !   AS NEEDED, CONTROLLED BY   c_nrkqq()
      do j = 1 , c_nreaqq
        rateqq(kk,j) = c_rkqq(1,j)*dexp(d_zero-c_rkqq(2,j) * &
                                       (d_one/c_temp(kk)-d_one/298.0D0))
      end do
      !
      ! ---------------------------
      ! END AQUEOUS RATE CONSTANTS
      ! ---------------------------
    end subroutine chemrates
!
! ---------------------------------------------------------------
!
! This calculates photolysis rates for specified location and time
!  based on interpolation in the lookup table TUVGRID2
!
! Average J-values are calculated over the specified time interval
!   Based on an looked-up j-values at specific times.
!
! The J-values are entered into the rate constant array (RATEK)
!  for use in the chemistry solver.
!
! Photolysis reactions in the file REACTION.DAT
!  are identified by a RATE INDEX (NRK) from 1 to 56.
!  The index identifies the photolysis rate number in the lookup table.
!  The final j-value is the looked-up value multiplied by
!   the coefficient RK1 (from REACTION.DAT)
!
! ----
! The input file TUVGRID2 is generated by Sandy Sillman
!   based on full 8-stream TUV photolysis calculation
!   from Sasha Madronich (Madronich and Flocke, 1998). (2002 version)
!   Program to generate TUVGRID2 is 'tuvtab2.f'.
! ----
!
! Input:  (from variables in common block)
!
!   Vectorized meteorological variables -
!     DOBSON  (dobson units)
!     ALBEDO  (wrong name for CLOUD ABOVE/BELOW INDEX)
!     DEPTH   cloud optical depth
!     AAERX   Aerosol optical depth
!     zenith angle, calculated by ZEN1 (entered as c_jparam(1))
!     ===> All entered into NONVECTORIZED  array c_jparam
!           (for subr jvalpro - vectorize in future?)
!
!   Time and location variables (used to get ZENITH ANGLE)
!     altmid  Altitude (km) at midpoint
!     c_lat   Latitude, degrees
!     c_lon   Longitude, degrees
!     c_hour  Hour at the end of the time interval (decimal hrs)
!            (note alternative OPTION: hour at CENTER of time interval;
!             useful if comparing with measured OH at specified time)
!     IDATE   Integer date YYMMDD (YY=100 for 2000)
!    ===> used to calculate SOLAR ZENITH ANGLE
!
! Output:
!    ratek() for photolysis reactions
!
! Called by:  either boxmain (as part of chem. solution, or for info)
!             or by quadchem (as part of chemistry solution)
!
! Calls to:   ZEN1 (sets zenith angle)
!             jvalpro (HV lookup for specified time and location)
!
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
!
    subroutine hvrates
!
      implicit none
!
!   hvrate(kk, 56)    Interpolated hv rates averaged over time interval.
!                                      (sec-1)
!   fhvint             Fraction assigned to each discrete time period
!                       used to interpolate for the full interval
!                       (=1/ihvint)
!   ihvint             Number of discrete times  used  to interpolate
!                       hv over the simulated time interval
!
!   jval(      56)       J-value output from jval2.f parameterization.
!                                      (sec-1)
!   zen(kk)           Zenith angle (degrees) (90=horizon)
!
      ! Interp hv over time interval
      real(dp) :: hvrate(c_kvec,56)
      ! j-value output from jval2.f
      real(dp) :: jval(56)
      ! zenith angle, degrees
      real(dp) :: zen(c_kvec)
      ! Partition btwen time intervals
      real(dp) :: fhvint
      ! Counter for time intervals
      integer :: ihour
      ! Number of time intervals
      integer :: ihvint
      ! Species index for hv array
      integer :: jc
      real(dp) :: xhouro
!
      kk = 1
      !
      ! HV RATES - HERE IS NEW INTERPOLATION FROM GRID (TUVGRID)
      !  PROGRAM WAS TESTED AS tuvtest2.f  and jvalmain1.f
      !   (see INFO in dhvmad02)
      ! INSERT REQUIRES DOBSON, XHM, IDATE, zen,
      !   HV RATES CALCULATED BASED ON AVERAGE OVER TIME PERIOD
      ! AVERAGE CALCULATED FROM HV RATES AT X TIME INTERVALS THROUGH THE PERIO
      ! X-TIME-INTERVAL ADDITION (1050 and 1600 AT START OF HV RATE CALC.).

      xhouro = c_hour
      ihvint = int(dabs(c_time)/1800.0D0)
      if ( ihvint < 1 ) ihvint = 1
      fhvint = d_one/(d_zero+ihvint)
      do ij = 1 , 56
        hvrate(kk,ij) = d_zero
      end do
      !
      ! NOTE OPTION:  hv FOR TIME INTERVAL ENDING AT THE  HOUR (preferred)
      !               OR hv FOR TIME INTERVAL CENTERED AT SPECIFIED HOUR.
      !    (CENTERED = allows exact comparison with measured OH at exact time
      !      But ENDING = dynamically correct calculation for time interval).
      !      Correct - ENDING with hv TIME INTERVAL CALCULATION = 1).
      !
      do ihour = 1 , ihvint
        c_hour = xhouro+(c_time/3600.0D0)*fhvint*(ihour-0.5D0-(1.0D0*ihvint))
        !
        ! END X-TIME-INTERVAL ADDITION
        !
        ! ESTABLISH COS(ZENITH) AS A FUNCTION OF TIME, DATE AND LATITUDE
        ! ZEN1 MODIFIED - USES XHOUR, DOES NOT ADJUST FOR TIME INTERPOLATION
        !      call ZEN1(zen)
        ! ENTER ZENITH ANGLE INTO PARAMETER MATRIX.
        !  NONVECTORIZED - NOTE ZENITH IS VECTOR.
        c_jparam(20) = c_idate
        !
        !  NOTE:  VALUES SHOULD BE ENTERED FROM BOX MODEL CALC (ALTITUDE)
        !  OR FROM STORED DATA:
        !  c_jparam(1) = zenith angle     2=altitude (kPA or km)
        !   3=ozone column,  4=SO2 column, 5=NO2 column
        !   6= aerosol optical depth,  7=albedo  (8-future ->ground in km)
        !   9=cloud-above optical depth, 10=cloud-below depth
        !   11=cloud-above altitude above current layer (km)
        !   12=cloud-below altitude above ground (km)
        !
        ! CALL JVALPRO FOR THESE CONDITIONS - sets jvalues
        !
        call jvalpro(c_nhv,c_hvmat,c_hvmatb,c_jarray,c_jparam,jval)
        !
        ! ADD TO ARRAY SUM FOR AVERAGE JVAL OVER TIME PERIOD WITH X INTERVALS.
        !
        do ij = 1 , 56
          hvrate(kk,ij) = hvrate(kk,ij)+fhvint*jval(ij)
          c_jval(kk,ij) = hvrate(kk,ij)
        end do
        !
        ! X-TIME-INTERVAL ADDITION.
        !
        !         END X-TIME-INTERVAL LOOP AFTER 400. RESET XHOUR.
      end do
      c_hour = xhouro
      !
      ! ENTER HV RATES FROM TUV TABLE INTO RATE CONSTANT ARRAY
      !   INCLUDING MULTIPLICATIVE FACTOR AND SPECIAL HV RATE CONSTANTS
      do j = 1 , c_nreac
        jc = c_nrk(j)
        if ( jc > 0 .and. jc < 56 ) then
          ratek(kk,j) = hvrate(kk,jc)*c_rk(1,j)
        end if
        !
        !  SPECIAL HV RATE CONSTANTS
        !  RATE FOR 9.O3+hv=2OH (IF PRODUCT=OH)
        !
        if ( jc == 102 ) then
          if ( c_nh2o > 0 ) then
            ratek(kk,j) = hvrate(kk,2)*c_rk(1,j) / &
              (d_one+0.13181D0*c_dens(kk)/c_h2ogas(kk))
          end if
        end if
      end do
    end subroutine hvrates
!
! -----------------------------------------------------------
!
!  this is used to set termolecular rates.
! Need TEMP,DENS from COMMON.  CUT; passed on in argument.
! MODIFIED TO FIT WITH NASA 1992 FORM
!
! Inputs:
!  Parameters (From DeMore et al., 1997 (JPL), p.8)
!    b = ko(300) = Low pressure limit at 300 K.
!    c = exponent for temperature adjustment, -n in NASA.
!        (ko(T)=ko(300)*(T/300)**c)
!    d = kinf(300) = High pressure limit at 300 K.
!    e = exponent for temperature admustment, -m in NASA.
!      (kinf(T)=kinf(300)*(T/300)**e)
!    u = Fc, base for exponent in log-T-P adjustment, always 0.6 in JPL.
!  Meteorology
!    tempx = temperature
!    denx  = density, molec/cm3
!
! Output:
!    bod = reaction rate cm3-sec
!
! Called by:  chemrates
! Calls to:   none
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
!
! Output reaction rate in cm3-sec
!
    real(dp) function bod(u,b,c,d,e,tempx,denx)
!
      implicit none
!
      ! ko(300) = Low pressure limit at 300 K.
      real(dp) , intent(in) :: b
      ! exponent for temperature, -n in JPL.
      real(dp) , intent(in) :: c
      ! kinf(300), High pressure limit at 300 K
      real(dp) , intent(in) :: d
      ! exponent for temperature, -m in JPL.
      real(dp) , intent(in) :: e
      ! Fc, base for log-T-P adjustment (0.6)
      real(dp) , intent(in) :: u
      ! temperature K
      real(dp) , intent(in) :: tempx
      ! density, molec/cm3
      real(dp) , intent(in) :: denx
!
      ! Interim parameter
      real(dp) :: f1
      ! Interim parameter
      real(dp) :: f2
      ! Interim parameter
      real(dp) :: ee
!
!     f1 = b*(tempx**c)*denx
!     f2 = f1 / (d*(tempx**e))
      !
      ! MODIFICATION TO TEMPX/300.
      !
      f1 = b*((tempx/300.0D0)**c)*denx
      f2 = f1 / (d*((tempx/300.0D0)**e))
      ee = d_one / ( d_one + (dlog10(f2))**d_two )
      bod = (f1/(d_one+f2) ) * u**ee
    end function bod
!
!  --------------------------------------------------------------
!
!  this sets the yield of alkyl nitrates from RO2+NO reactions.
!
!  Inputs:
!    c =  Number of carbon atoms in the RO2 (non-integer for lumped RO2)
!    tempx = temperature
!    denx  = density, molec/cm3
!
!  Output:
!    ytn = RNO3 yield (stoichiometry)
!
! Called by:  chemrates
! Calls to:   none
!
! ---------------------------------------------
! History:
!  12/06 Written by Sandy Sillman from boxchemv7.f
!
! -------------------------------------------------------------------
!
    real(dp) function ytn(c,tempx,denx)
!
      implicit none
!
      ! Number of carbon atoms in RO2
      real(dp) , intent(in) :: c
      ! temperature K
      real(dp) , intent(in) :: tempx
      ! density, molec/cm3
      real(dp) , intent(in) :: denx
!
      real(dp) :: x
      real(dp) :: y
      real(dp) :: z
      real(dp) , parameter :: par = 4.3D-25

      x = par*tempx*dexp(1.08D0*c)*denx*(300.0D0/tempx)**5.05D0
      y = 0.384D0*(300.0D0/tempx)**4.16D0
      z = d_one / ( d_one + (dlog10(x/y))**d_two )
      ytn = ( x / (d_one + x/y ))*(0.467D0**z)
    end function ytn
!
end module mod_cbmz_rates1
