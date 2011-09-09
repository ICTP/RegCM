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
!                                                                       
!                                                                       
      subroutine chemrates 
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
! ---------------------------------------------------------             
      IMPLICIT NONE 
                                                                        
      include 'chemmech.EXT' 
      INCLUDE 'chemvars.EXT' 
      include 'chemlocal.EXT' 
                                                                        
! LOCAL VARIABLES                                                       
                               ! temperature K                          
      double precision tempx 
                               ! density, molec/cm3                     
      double precision denx 
                               ! Function: 3-bodyreaction rate cm3-sec  
      double precision bod 
                               ! Function: RNO3 yield from RO2+NO       
      double precision ytn 
                                                                        
! KK = VECTORIZED LOOP VARIABLE, SET AT ONE                             
           kk=1 
!                                                                       
! pi (if not built-in special)                                          
! AND "RU" PARAMETER FOR  molecular speed (Barth, personal comm)        
!   (8.314e7 g-cm2/s2-mol-K)   (note: 8.314e0 kg-m2/s2-mol-K)           
!  Moved to CHEMVAR.EXT as PARAMETERS                                   
!     pii = 3.141592654                                                 
!      RUMOLEC = 8.314E7                                                
                                                                        
                                                                        
!  SET ALL RATE CONSTANTS (including initial index values for HV rates,)
                                                                        
! --------------------------                                            
! LOOP TO SET RATE CONSTANTS                                            
! --------------------------                                            
                 ! counter for parameterized ro2-ro2 reactions          
      nr1 = 0 
      do     j=1,c_nreac 
                                                                        
! Default                                                               
       ratek(kk, j) = c_rk(1,j) 
                                                                        
! Standard 2-body reaction format                                       
       if(   c_nrk(j).eq.0) then 
!        do  kk=1,c_kmax                                                
           ratek(kk ,j)= c_rk(1,j)*exp(-c_rk(2,j)/c_temp(kk)) 
!        enddo                                                          
       endif 
                                                                        
! Standard 2-body reaction format with input rate at 298 K.             
       if(   c_nrk(j).eq.-1) then 
!        do  kk=1,c_kmax                                                
           ratek(kk ,j)= c_rk(1,j)                                      &
     &       *exp(-c_rk(2,j)*(1./c_temp(kk)-0.0033557047))              
!        enddo                                                          
       endif 
                                                                        
! Standard format *T/300                                                
       if(   c_nrk(j).eq.-2) then 
!        do kk=1,c_kmax                                                 
           ratek(kk ,j)= c_rk(1,j)*exp(-c_rk(2,j)/c_temp(kk))           &
     &                 *(c_temp(kk)/300.)**c_rk(3,j)                    
!        enddo                                                          
       endif 
                                                                        
! Standard 3-body reaction format - bod function                        
       if(   c_nrk(j).eq.-3) then 
!        do  kk=1,c_kmax                                                
           tempx = c_temp(kk) 
           denx = c_dens(kk) 
           ratek(kk ,j)=bod( c_rk(1,j), c_rk(2,j), c_rk(3,j),           &
     &                       c_rk(4,j), c_rk(5,j), tempx, denx )        
!        enddo                                                          
       endif 
                                                                        
! Combined 2-body and  3-body reaction format                           
       if(   c_nrk(j).eq.-4) then 
!        do  kk=1,c_kmax                                                
           tempx = c_temp(kk) 
           denx = c_dens(kk) 
           ratek(kk ,j)=bod( c_rk(3,j), c_rk(4,j), c_rk(5,j),           &
     &                       c_rk(6,j), c_rk(7,j), tempx, denx )        
           ratek(kk ,j)= ratek(kk,j)                                    &
     &                    * c_rk(1,j)*exp(-c_rk(2,j)/tempx)             
                                                                        
!        enddo                                                          
       endif 
                                                                        
! 2-body reaction format with nitrate yield subtracted                  
       if(   c_nrk(j).eq.-5) then 
!        do  kk=1,c_kmax                                                
           ratek(kk ,j)= c_rk(1,j)*exp(-c_rk(2,j)/c_temp(kk)) 
           ratek(kk, j) = ratek(kk, j )                                 &
     &       * (1.- ytn(c_rk(3,j), c_temp(kk ), c_dens(kk )) )          
!        enddo                                                          
       endif 
                                                                        
! 2-body reaction format for  nitrate yield                             
       if(   c_nrk(j).eq.-6) then 
!        do  kk=1,c_kmax                                                
           ratek(kk ,j)= c_rk(1,j)*exp(-c_rk(2,j)/c_temp(kk)) 
           ratek(kk, j) = ratek(kk, j )                                 &
     &       * (    ytn(c_rk(3,j), c_temp(kk ), c_dens(kk )) )          
!        enddo                                                          
       endif 
                                                                        
! Density-dependent format A exp(-B/T)* (C + D*dens)/(E + F*dens)       
!    for CO, GLYX+OH                                                    
       if(   c_nrk(j).eq.-7) then 
!        do  kk=1,c_kmax                                                
           ratek(kk ,j)= c_rk(1,j)*exp(-c_rk(2,j)/c_temp(kk)) 
           ratek(kk,j) = ratek(kk,j)                                    &
     &      *(c_rk(3,j) + c_rk(4,j)*c_dens(kk))                         &
     &      /(c_rk(5,j) + c_rk(6,j)*c_dens(kk))                         
!        enddo                                                          
       endif 
                                                                        
                                                                        
! Combined format, interacting two-body reactions with density          
!   (for MO2+MO2 and DMS)                                               
! A exp(-B/T) / (1 + (C+E*dens)* exp(-D/T))                             
       if(   c_nrk(j).eq.-8) then 
!        do  kk=1,c_kmax                                                
           ratek(kk ,j)= c_rk(1,j)*exp(-c_rk(2,j)/c_temp(kk)) 
           ratek(kk,j) = ratek(kk,j) / ( 1. +                           &
     &              (c_rk(3,j) + c_rk(5,j)*c_dens(kk) )                 &
     &                   *exp(-c_rk(4,j)/c_temp(kk))    )               
!        enddo                                                          
       endif 
                                                                        
! Special HO2+HO2 rate: (Aexp(-B/t) + C*dens*exp(-D/t)) *(1+ E exp(-F/T)
!  from Chameides around 1996                                           
       if(   c_nrk(j).eq.-9) then 
!        do  kk=1,c_kmax                                                
           ratek(kk ,j)= c_rk(1,j)*exp(-c_rk(2,j)/c_temp(kk)) 
           ratek(kk,j) = ratek(kk,j) +                                  &
     &              c_rk(3,j)*exp(-c_rk(4,j)/c_temp(kk)) * c_dens(kk)   
           ratek(kk,j) = ratek(kk,j) *  (1. +                           &
     &        c_rk(3,j)*exp(-c_rk(4,j)/c_temp(kk)) * xc(kk, 14)  )      
                                                                        
!        enddo                                                          
       endif 
                                                                        
                                                                        
! Special HNO3+OH rate:                                                 
!   A exp(-B/t) + C exp(-D/t)*E exp(-F/t)/(Cexp(-D/T)+Eexp(-F/T))       
!  from Evans, Fiore 2003                                               
       if(   c_nrk(j).eq.-10) then 
!        do  kk=1,c_kmax                                                
           ratek(kk ,j)= c_rk(1,j)*exp(-c_rk(2,j)/c_temp(kk)) 
           beta(kk)    = c_rk(3,j)*exp(-c_rk(4,j)/c_temp(kk))           &
     &                              *c_dens(kk)                         
           gamma(kk)   = c_rk(5,j)*exp(-c_rk(6,j)/c_temp(kk)) 
                                                                        
           ratek(kk,j) = ratek(kk,j) +                                  &
     &          beta(kk)*gamma(kk)/(beta(kk)+gamma(kk))                 
                                                                        
!        enddo                                                          
       endif 
                                                                        
! Special C3H8 format (IUPAC 2002):  Aexp(-B/t)/(1+C*T/300^(D*exp(-E/t))
       if(   c_nrk(j).eq.-11) then 
!        do  kk=1,c_kmax                                                
           ratek(kk ,j)= c_rk(1,j)*exp(-c_rk(2,j)/c_temp(kk)) 
           beta(kk)    = c_rk(4,j)*exp(-c_rk(5,j)/c_temp(kk)) 
                                                                        
           ratek(kk,j) = ratek(kk,j) / (1. +                            &
     &         c_rk(3,j)* (c_temp(kk)/300.)**beta(kk)   )               
                                                                        
!        enddo                                                          
       endif 
                                                                        
! Sum of two reactions (for ACET+OH, JPL 2002):  A exp(-B/T) + C exp(-D/
       if(   c_nrk(j).eq.-12) then 
!        do  kk=1,c_kmax                                                
           ratek(kk ,j)= c_rk(1,j)*exp(-c_rk(2,j)/c_temp(kk))           &
     &                   + c_rk(3,j)*exp(-c_rk(4,j)/c_temp(kk))         
!        enddo                                                          
       endif 
                                                                        
! Parameterized RO2-RO2 reactions                                       
!      (CBMZ format: RO2->partial product from RO2-RO2 reaction,        
!        see setro2 in chemsolve.f for exact information.)              
       if(   c_nrk(j).eq.-13) then 
         nr1 = nr1+1 
!        do  kk=1,c_kmax                                                
           ratek(kk ,j)= c_rk(1,j)*exp(-c_rk(2,j)/c_temp(kk)) 
           ratero2(kk,nr1,1) = ratek(kk,j) 
           ratero2(kk ,nr1,2)=                                          &
     &          (c_rk(3,j)*exp(-c_rk(4,j)/c_temp(kk)))**0.5             
!        enddo                                                          
       endif 
                                                                        
! REACTION ON SULFATE AEROSOL SURFACES.  c_nrk(j) = -20                 
!   ratek = saersa*(gamma/4)* [8RT/(pi*MW)]^0.5 (Pandis and Seinfeld)   
!   for saersa, aerosol surface area in cm2/cm3                         
!   gamma = reaction efficiency (input)                                 
!   MW = molecular weight (input)                                       
!                                                                       
       if(   c_nrk(j).eq.-20) then 
!        do  kk=1,c_kmax                                                
           ratek(kk ,j)=  c_saersa(kk) * (0.25d0*c_rk(1,j))             &
     &       * (  (8.d0*RUMOLEC   * c_temp(kk))/(pii*c_rk(2,j))  )**0.5 
!        enddo                                                          
       endif 
                                                                        
                                                                        
                       ! END SET ALL RATE CONSTANTS LOOP                
      enddo 
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
                                                                        
!        (Replaced by parameter values in common)                       
!     atmos = 2.247E+19                                                 
!     avogadrl = 6.02E+20                                               
!     rtcon= 0.08206                                                    
                                                                        
!  FUTURE CHANGE: ADD ALT. SPECIAL FORMULAS FOR HENRYS LAW CONSTANTS    
!   AS NEEDED, CONTROLLED BY   c_nrkh()                                 
                                                                        
      do 200 j=1,c_nreach 
!      do 18051  kk=1,c_kmax                                            
         rateh(kk ,j)= c_rkh(1,j)                                       &
     &             *exp(0.-c_rkh(2,j)*(1./c_temp(kk)-1./298.) )         
         rateh(kk ,j) = rateh(kk, j)                                    &
     &    *rtcon*c_temp(kk)*c_h2oliq(kk)                                
!                                                                       
! (OLDER CONVERSION)                                                    
!        rateh(kk ,j) = rateh(kk, j)                                    
!    *       *c_h2oliq(kk) * avogadrl/atmos                             
!                                                                       
18051  continue 
  200 continue 
                                                                        
! H+ AQUEOUS EQUILIBRIUM COEFFICIENTS (MOLES/LITER)                     
!                                                                       
!  FUTURE CHANGE: ADD ALT. SPECIAL FORMULAS FOR EQUILIBRIUM COEFFICIENTS
!   AS NEEDED, CONTROLLED BY   c_nrkq()                                 
                                                                        
      do 250 j=1,c_nreacq 
!      do 18061  kk=1,c_kmax                                            
         rateq(kk ,j)= c_rkq(1,j)                                       &
     &          *exp(0.-c_rkq(2,j)*(1./c_temp(kk)-1./298.) )            
18061  continue 
  250 continue 
                                                                        
! (CONVERSION:  AQUEOUS K1 IS MULTIPLIED BY HENRY'S LAW CONSTANT        
!   AQUEOUS K2 IS MULTIPLIED BY K1KH, etc.                              
! DELETED )                                                             
                                                                        
!  SPECIAL AQUEOUS EQUILIBRIUM COEFFICIENTS (MOLES/LITER)               
!  FUTURE CHANGE: ADD ALT. SPECIAL FORMULAS FOR SPECIAL EQ  COEFFICIENTS
!   AS NEEDED, CONTROLLED BY   c_nrkqq()                                
                                                                        
                                                                        
      do 255 j=1,c_nreaqq 
!      do 18062  kk=1,c_kmax                                            
         rateqq(kk ,j)= c_rkqq(1,j)                                     &
     &           *exp(0.-c_rkqq(2,j)*(1./c_temp(kk)-1./298.) )          
18062  continue 
  255 continue 
                                                                        
! ---------------------------                                           
! END AQUEOUS RATE CONSTANTS                                            
! ---------------------------                                           
                                                                        
! END BRATES                                                            
 2000   return 
      END                                           
! ---------------------------------------------------------------       
! ----------------------------------------------------------            
                                                                        
                                                                        
                                                                        
      subroutine hvrates 
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
!                                                                       
! -------------------------------------------------------------         
!                                                                       
      IMPLICIT NONE 
                                                                        
      include 'chemmech.EXT' 
      INCLUDE 'chemvars.EXT' 
      include 'chemlocal.EXT' 
                                                                        
!  LOCAL VARIABLES FOR J-VALUE CALCULATION                              
                                                                        
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
      double precision hvrate(c_kvec ,56) 
                                        ! j-value output from jval2.f   
      double precision jval(       56) 
                                        ! zenith angle, degrees         
      double precision zen(c_kvec) 
                                                                        
                                   ! Partition btwen time intervals     
      double precision fhvint 
                                    ! Counter for time intervals        
      integer ihour 
                                    ! Number of time intervals          
      integer ihvint 
                                    ! Species index for hv array        
      integer jc 
      double precision xhouro 
                                                                        
! ---------------------------------------                               
             kk=1 
                                                                        
! HV RATES - HERE IS NEW INTERPOLATION FROM GRID (TUVGRID)              
!  PROGRAM WAS TESTED AS tuvtest2.f  and jvalmain1.f                    
!   (see INFO in dhvmad02)                                              
                                                                        
! INSERT REQUIRES DOBSON, XHM, IDATE, zen,                              
                                                                        
!   HV RATES CALCULATED BASED ON AVERAGE OVER TIME PERIOD               
! AVERAGE CALCULATED FROM HV RATES AT X TIME INTERVALS THROUGH THE PERIO
! X-TIME-INTERVAL ADDITION (1050 and 1600 AT START OF HV RATE CALC.).   
                                                                        
       xhouro = c_hour 
       ihvint = int(abs(c_time)/1800.) 
       if(ihvint.lt.1) ihvint = 1 
!      if(ihvint.lt.4) ihvint = 4                                       
       fhvint = 1./(0.+ihvint) 
!       if(c_hour.gt.12.and.c_hour.le.15) write(c_out,*) ihvint, fhvint 
        do 1050 ij=1,56 
!              do 18016 kk=1,c_kmax                                     
        hvrate( kk ,ij) = 0. 
18016          continue 
 1050    continue 
                                                                        
! NOTE OPTION:  hv FOR TIME INTERVAL ENDING AT THE  HOUR (preferred)    
!               OR hv FOR TIME INTERVAL CENTERED AT SPECIFIED HOUR.     
!    (CENTERED = allows exact comparison with measured OH at exact time 
!      But ENDING = dynamically correct calculation for time interval). 
!      Correct - ENDING with hv TIME INTERVAL CALCULATION = 1).         
                                                                        
       do 1600 ihour = 1,ihvint 
       c_hour = xhouro + (c_time/3600.)*fhvint*(ihour-0.5-(1.0*ihvint)) 
!      c_hour = xhouro + (c_time/3600.)*fhvint*(ihour-0.5-(0.5*ihvint)) 
                                                                        
!      write(c_out,*) ihour,xhouro,c_hour                               
                                                                        
! END X-TIME-INTERVAL ADDITION                                          
                                                                        
! ESTABLISH COS(ZENITH) AS A FUNCTION OF TIME, DATE AND LATITUDE        
!  ZEN1 MODIFIED - USES XHOUR, DOES NOT ADJUST FOR TIME INTERPOLATION (X
                                                                        
!      call ZEN1(zen) 
                                                                        
! ENTER ZENITH ANGLE INTO PARAMETER MATRIX.                             
!  NONVECTORIZED - NOTE ZENITH IS VECTOR.                               
                                                                        
!      c_jparam(1) = zen(1) 
      c_jparam(20) = c_IDATE 
                                                                        
!  NOTE:  VALUES SHOULD BE ENTERED FROM BOX MODEL CALC (ALTITUDE)       
!  OR FROM STORED DATA:                                                 
!  c_jparam(1) = zenith angle     2=altitude (kPA or km)                
!   3=ozone column,  4=SO2 column, 5=NO2 column                         
!   6= aerosol optical depth,  7=albedo  (8-future ->ground in km)      
!   9=cloud-above optical depth, 10=cloud-below depth                   
!   11=cloud-above altitude above current layer (km)                    
!   12=cloud-below altitude above ground (km)                           
                                                                        
! CALL JVALPRO FOR THESE CONDITIONS - sets jvalues                      
                                                                        
      call jvalpro(c_nhv,c_hvmat, c_hvmatb, c_jarray,c_jparam, jval) 
                                                                        
! ADD TO ARRAY SUM FOR AVERAGE JVAL OVER TIME PERIOD WITH X INTERVALS.  
                                                                        
      do ij=1,56 
!     if(zen( kk ).lt.94.)                                              
       hvrate(kk ,ij) = hvrate(kk  ,ij)                                 &
     &  +fhvint*jval(ij)                                                
       c_jval(kk,ij)  = hvrate(kk ,ij)
      enddo 
                                                                        
! TEST WRITES                                                           
      if(c_kkw.gt.0) then 
         write(c_out,1501) ihour, zen, fhvint,( hvrate(1,ij),ij=1,4) 
 1501    format(' IHOUR ZENITH JVAL(1-4)=',i5,f8.2, 5(1Pe10.3)) 
         write(c_out,1502) c_jparam 
 1502    format(8(1pe10.3)) 
      endif 
                                                                        
! X-TIME-INTERVAL ADDITION.                                             
!         END X-TIME-INTERVAL LOOP AFTER 400. RESET XHOUR.              
 1600  continue 
      c_hour = xhouro 
                                                                        
! ENTER HV RATES FROM TUV TABLE INTO RATE CONSTANT ARRAY                
!   INCLUDING MULTIPLICATIVE FACTOR AND SPECIAL HV RATE CONSTANTS       
                                                                        
      do     j=1,c_nreac 
        jc = c_nrk(j) 
        if(jc    .gt.0.and.jc    .lt.56) then 
! vec                                    do kk=1,c_kvec                 
         ratek( kk ,j) = hvrate( kk ,jc)                                &
     &     * c_rk(1,j)                                                  
! vec                                    enddo  !do kk=1,c_kvec         
                            !if(c_nrk(j).gt.0...                        
        endif 
                                                                        
!  SPECIAL HV RATE CONSTANTS                                            
                                                                        
!  RATE FOR 9.O3+hv=2OH (IF PRODUCT=OH)                                 
        if(jc.eq.102) then 
         if(c_nh2o.gt.0) then 
! vec                                    do kk=1,c_kvec                 
         ratek( kk ,j) = hvrate( kk , 2)                                &
     &     * c_rk(1,j)                                                  &
     &     /( 1.+0.13181*c_dens( kk )/c_h2ogas(kk)    )                 
!    *     /( 1.+c_rk(2,j)*c_dens( kk )/c_h2ogas(kk)    )               
!    *     /( 1.+0.13181*c_dens( kk )/xc( kk ,c_nh2o) )                 
                                                                        
! vec                                    enddo  !do kk=1,c_kvec         
                           !if(c_nh2o.gt.0) then                        
         endif 
                         !if(jc.eq.102) then                            
        endif 
                                                                        
! TEST WRITE                                                            
        if(c_kkw.gt.0) then 
          if(c_rk(1,j) .eq. 1.000)  then 
            if(jc.eq.4.or.jc.eq.2.or.jc.eq.102) then 
              write(c_out,802) c_hour,  c_treac(1,j), ratek(c_kkw,j) 
  802         format(/,' TEST XHOUR, SPEC+HV JVAL= ',                   &
     &                  f6.2,2x,a8,2x,1pe10.3 )                         
            if(jc.eq.102) write(43,*) 'hvrate, dens, h2o=',             &
     &          hvrate(c_kkw,2), c_dens(c_kkw), c_h2ogas(c_kkw)         
                                                                        
            endif 
          endif 
        endif 
                                                                        
! TEST WRITE                                                            
!      if(c_kkw.gt.0) write(c_out,*)'J JC SPEC JVAL=',                  
!    &         j,jc,c_treac(1,j), ratek(c_kkw,j)                        
                                                                        
                           !do     j=1,c_nreac                          
      enddo 
                                                                        
!  END HVRATES                                                          
 2000 return 
      END                                           
! -----------------------------------------------------------           
                                                                        
                                                                        
                                                                        
      FUNCTION bod(u,b,c,d,e, tempx, denx) 
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
!  --------------------------------------------------------------       
      IMPLICIT NONE 
                                                                        
!                                                                       
! LOCAL VARIABLES:  INPUT/OUTPUT                                        
                               ! ko(300) = Low pressure limit at 300 K. 
      double precision b 
                               ! exponent for temperature, -n in JPL.   
      double precision c 
                               ! kinf(300), High pressure limit at 300 K
      double precision d 
                               ! exponent for temperature, -m in JPL.   
      double precision e 
                               ! Fc, base for log-T-P adjustment (0.6)  
      double precision u 
                               ! temperature K                          
      double precision tempx 
                               ! density, molec/cm3                     
      double precision denx 
                               ! Output reaction rate cm3-sec           
      double precision bod 
!                                                                       
!                                                                       
! LOCAL VARIABLES:  INTERNAL                                            
                                     ! Interim parameter                
      double precision f1 
                                     ! Interim parameter                
      double precision f2 
                                     ! Interim parameter                
      double precision ee 
!                                                                       
!                                                                       
!     f1=b*(tempx**c)*denx                                              
!     f2=f1/( d * (tempx**e))                                           
! MODIFICATION TO TEMPX/300.                                            
      f1=b*((tempx/300.)**c)*denx 
      f2=f1/( d * ((tempx/300.)**e)) 
                                                                        
      ee= 1. / ( 1 + (log10(f2))**2. ) 
      bod = (f1/(1+f2) ) * u ** ee 
                                                                        
      return 
      END                                           
!  --------------------------------------------------------------       
                                                                        
      FUNCTION ytn(c, tempx, denx) 
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
!  --------------------------------------------------------------       
      IMPLICIT NONE 
                                                                        
!                                                                       
!                                                                       
! LOCAL VARIABLES:  INPUT/OUTPUT                                        
                               ! Number of carbon atoms in RO2          
      double precision c 
                               ! temperature K                          
      double precision tempx 
                               ! density, molec/cm3                     
      double precision denx 
                               ! Output: RNO3 yield from RO2+NO         
      double precision ytn 
!                                                                       
! LOCAL VARIABLES:  INTERNAL                                            
                                                                        
      double precision x 
      double precision y 
      double precision z 
      double precision par 
      PARAMETER (par=4.3D-25) 
                                                                        
      x=par*tempx*exp(1.08*c)*denx*(300./tempx)**5.05 
      y=0.384*(300./tempx)**4.16 
      z=1./ ( 1.+ (log10(x/y))**2. ) 
                                                                        
      ytn= ( x / (1.+ x/y ))*(0.467**z) 
                                                                        
      return 
      END                                           
!  --------------------------------------------------------------       
                                                                        
                                                                        
                                                                        
      subroutine ZEN1(ZENITH) 
                                                                        
! This calculates solar zenith and azimuth angles                       
!   for a particular time and location.                                 
                                                                        
!  Based on equations given by W.H.Smart,                               
!     "Textbook of Spherical Astronomy," 6th ed. Cambridge U. Press     
!                                                                       
!  Presented here in vectorized form                                    
!   (for vectorized latitude and longitude, scalar time)                
!   with common-block variables from the Sillman chemical solver.       
!                                                                       
! NOTE ALTERNATIVES: YYMMDD (standard)  or YYYYDDD                      
!                                                                       
! Inputs:                                                               
!       LAT - latitude in decimal degrees (LAT -> c_lat(kk))            
!       LONG - longitude in decimal degrees (LONG-> c_lon(kk))          
!       IDATE - Date at Greenwich:                                      
!           ALTERNATIVES:  YYMMDD (standard) or YYYYDDD                 
!                   specify year (19yy), month (mm), day (dd)           
!                   format is six-digit integer:  yymmdd                
!              ( OR enter as negative day number (1-366))               
!                                                                       
!       c_hour  Decimal time, EST.                                      
!          Note: exact time is set in hvrates based on time interval.   
!                                                                       
!       (c_time  time step, secs - not used here)                       
!       (original: GMT  - Greenwich mean time - decimal military        
!                         e.g.  22.75 = 45 min after ten pm gmt         
! OUTPUT                                                                
!       ZENITH(kk) = zenith angle, degrees (vectorized)                 
!       AZMUTH =  Azimuth                                               
!                                                                       
!                                                                       
! Called by:  hvrates                                                   
! Calls to:   none                                                      
!                                                                       
! ---------------------------------------------                         
! History:                                                              
! Written by Brian A. Ridley at York Univ., ca. 1978                    
!  based on equations given by W.H.Smart,                               
!    "Textbook of Spherical Astronomy,"	6th ed. Cambridge U. Press      
!                                                                       
!  1987 - given to Sandy Sillman from Michael Prather                   
!                                                                       
!  1994:  modified by Sandy Sillman for inclusion in chem. solver       
!  2005:  modified for vectorized form.                                 
!  12/06 Modified by Sandy Sillman from boxchemv7.f                     
!                                                                       
! -------------------------------------------------------------------   
!  --------------------------------------------------------------       
                                                                        
! -----------------------------------------------------                 
      IMPLICIT NONE 
                                                                        
                                                                        
                               ! Needed for output file number          
      include 'chemmech.EXT' 
                               ! Needed for c_kvec,                     
      INCLUDE 'chemvars.EXT' 
!                          also includes c_lat, c_lon, IDATE. c_hour    
                                                                        
! LOCAL VARIABLES:  INPUT   (all in chemvars.EXT)                       
!                                                                       
!   (LAT ->c_lat(c_kvec)  vectorized latitude)                          
!   (LON ->c_lon(c_kvec)  vectorized longitude)                         
!   (IDATE)                                                             
!   (c_hour = decimal hour EST)                                         
!                                                                       
! LOCAL VARIABLES:  OUTPUT                                              
                                         ! degrees, vectorized for CTM  
      double precision ZENITH(c_kvec) 
                                         ! azimuth angle, degrees       
      double precision AZMUTH 
!                                                                       
! LOCAL VARIABLES:  INTERNAL                                            
      DOUBLE PRECISION GMT, LBGMT,LZGMT 
                                  ! PI                                  
      DOUBLE PRECISION PI 
                                  ! PI/180, convert degrees to radians  
      DOUBLE PRECISION PIR 
      INTEGER NYEARS, LEAP, NOLEAP,  IN, I 
                                               ! Parse date             
      INTEGER IIYEAR,IYEAR,IMTH,IDAY, IIY 
      INTEGER   IMN(12) 
                          ! days in year from 1973 to present year      
      INTEGER YR 
                          ! DAY NUMBER                                  
      INTEGER IJD 
                          ! Number of days since JAN 1, 1973            
      INTEGER  JD 
                          ! Number of years since 1973                  
      INTEGER  IJ73 
                          ! Number of days with time fraction added.    
      DOUBLE PRECISION  D 
!                                                                       
                                     ! Geom. mean longitude             
      DOUBLE PRECISION ML 
                                     ! Geom mean longitude radians      
      DOUBLE PRECISION RML 
                                     ! mean long. of perigee,deg        
      DOUBLE PRECISION W 
                                     ! mean long. of peregee, rad.      
      DOUBLE PRECISION WR 
                                     ! cos(WR)                          
      DOUBLE PRECISION CW 
                                     ! sin(WR)                          
      DOUBLE PRECISION SW 
                                     ! sin(2WR)                         
      DOUBLE PRECISION SSW 
                                     ! Eccentricity                     
      DOUBLE PRECISION EC 
                                     ! Mean obliquity of ecliptic, deg  
      DOUBLE PRECISION EPSI 
                                     ! Mean obliquity of ecliptic, rad  
      DOUBLE PRECISION PEPSI 
                                     !  tangent calc                    
      DOUBLE PRECISION YT 
                                                                        
                                     ! 2.*EC*YT                         
      DOUBLE PRECISION EYT 
      DOUBLE PRECISION FEQT1 
      DOUBLE PRECISION FEQT2 
      DOUBLE PRECISION FEQT3 
      DOUBLE PRECISION FEQT4 
      DOUBLE PRECISION FEQT5 
      DOUBLE PRECISION FEQT6 
      DOUBLE PRECISION FEQT7 
      DOUBLE PRECISION FEQT 
                                 !equation of time, secs                
      DOUBLE PRECISION EQT 
                                 !equation of time in hours             
      DOUBLE PRECISION EQH 
                                 !equation of time , degrees            
      DOUBLE PRECISION REQT 
                                                                        
                                 ! right ascension, deg                 
      DOUBLE PRECISION RA 
                                 ! right ascension, rads                
      DOUBLE PRECISION RRA 
      DOUBLE PRECISION TAB 
                                 ! declination, rads                    
      DOUBLE PRECISION RDECL 
                                 ! declination, deg                     
      DOUBLE PRECISION DECL 
!                                                                       
                                 ! local hour angle, rad                
      DOUBLE PRECISION RLT 
                                 ! hour angle, rad                      
      DOUBLE PRECISION ZPT 
                                 ! cos(zen)                             
      DOUBLE PRECISION CSZ 
                                 ! zenith, radians                      
      DOUBLE PRECISION ZR 
                                                                        
                                         ! cos( azimuth)                
      DOUBLE PRECISION CAZ 
                                         ! azimuth (radians)            
      DOUBLE PRECISION RAZ 
                                                                        
! LOCAL VARIABLES ALSO IN chemmech, chemlocal                           
                                 ! vector index                         
      integer kk 
!                                                                       
      DATA IMN/31,28,31,30,31,30,31,31,30,31,30,31/ 
                                                                        
! KK = VECTORIZED LOOP VARIABLE, SET AT ONE                             
           kk=1 
                                                                        
! ENTER TIME INTERPOLATION PARAMETERS.  THESE ARE NOT VECTOR VBLS,      
!  HAVE BEEN MOVED AHEAD OF THE BIG VECTORIZED LOOP.                    
                                                                        
! *** NOTE: CURRENT VERSION ASSUMES XHOUR=EASTERN STANDARD TIME.        
!     THIS CAN BE CHANGED BY CHANGING GMT LINE.                         
                                                                        
!     GMT      = c_hour + 5. -c_time/3600.                              
      GMT      = c_hour + 5. 
      if(GMT     .gt.24) GMT      = GMT      - 24. 
                                                                        
! convert to radians                                                    
                                                                        
      PI = 3.1415926535590 
      PIR = PI/180. 
                                                                        
! parse date  YYMMDD                                                    
!  NOVEMBER 2005 CHANGE:  only if DATE>1000, else treat it as DAY NUMBER
!         if(c_IDATE.gt.0) then                                         
                                                                        
      IIYEAR = c_IDATE     /10000 
      IYEAR = 19*100 + IIYEAR 
      IMTH = (c_IDATE      - IIYEAR*10000)/100 
      IDAY = c_IDATE      - IIYEAR*10000 - IMTH*100 
                                                                        
! alternative algorithm for YYYYDDD                                     
!     IYEAR =  c_IDATE/1000                                             
!     IDAY = c_IDATE - IYEAR*1000                                       
! c     IDAY =   MOD ( IDATE  , 1000 )                                  
                                                                        
                                                                        
                                                                        
! identify and correct leap years                                       
                                                                        
      IIY = (IIYEAR/4)*4 
      IF(IIY.EQ.IIYEAR) IMN(2) = 29 
                                                                        
! count days from Dec.31,1973 to Jan 1                                  
                                                                        
      NYEARS = IYEAR - 1974 
      LEAP = (NYEARS+1)/4 
      IF(NYEARS.LE.-1) LEAP = (NYEARS-2)/4 
      NOLEAP = NYEARS - LEAP 
      YR = 365*NOLEAP + 366*LEAP 
                                                                        
                                                                        
! YYMMDD version continues here.                                        
! YYYYDDD     version would cut lines beginning here....                
      IJD = 0 
      IN = IMTH - 1 
      IF(IN.EQ.0) GO TO 40 
                                                                        
                                                                        
      DO 30 I=1,IN 
         IJD = IJD + IMN(I) 
   30 END DO 
                                                                        
      IJD = IJD + IDAY 
      GO TO 50 
   40 IJD = IDAY 
   50 IJ73 = IYEAR - 1973 
!                                                                       
! ... and ending here, and then insert  (for YYYYDDD version)           
!     IJD = IDAY                                                        
                                                                        
                                                                        
! julian days current "ijd"                                             
                                                                        
      JD = IJD + YR 
      D = JD + GMT     /24.0 
                                                                        
! NOTE:                                                                 
!  IJD = DAY NUMBER                                                     
!  JD = NUMBER OF DAYS SINCE JAN 1, 1973                                
!        = DAY NUMBER + YR (days in years from 1973 to present year)    
!  IJ = number of years since 1973                                      
!  D = number of days with GMT fraction added in.                       
                                                                        
! END  NOVEMBER 2005 CHANGE:  only if DATE>0,                           
!                                   else treat it as DAY NUMBER         
!         endif             !if(c_IDATE.gt.0) then                      
                                                                        
! TEST WRITE: IDATE D                                                   
!  980101   8.7671770833333339E+03                                      
!  Days from January 1, 1974 to January 1, 1998 = 24*365+6=8766.        
!  Calculation is all based on this number.                             
                                                                        
!      write(44,*) c_IDATE, D                                           
                                                                        
                                                                        
! calc geom mean longitude                                              
                                                                        
      ML = 279.2801988 + .9856473354*D + 2.267E-13*D*D 
      RML = ML*PIR 
                                                                        
                                                                        
! calc equation of time in sec                                          
!  w = mean long of perigee                                             
!  e = eccentricity                                                     
!  epsi = mean obliquity of ecliptic                                    
                                                                        
      W = 282.4932328 + 4.70684E-5*D + 3.39E-13*D*D 
      WR = W*PIR 
      EC = 1.6720041E-2 - 1.1444E-9*D - 9.4E-17*D*D 
      EPSI = 23.44266511 - 3.5626E-7*D - 1.23E-15*D*D 
      PEPSI = EPSI*PIR 
      YT = (TAN(PEPSI/2.0))**2 
      CW = COS(WR) 
      SW = SIN(WR) 
      SSW = SIN(2.0*WR) 
      EYT = 2.*EC*YT 
      FEQT1 = SIN(RML)*(-EYT*CW - 2.*EC*CW) 
      FEQT2 = COS(RML)*(2.*EC*SW - EYT*SW) 
      FEQT3 = SIN(2.*RML)*(YT - (5.*EC**2/4.)*(CW**2-SW**2)) 
      FEQT4 = COS(2.*RML)*(5.*EC**2*SSW/4.) 
      FEQT5 = SIN(3.*RML)*(EYT*CW) 
      FEQT6 = COS(3.*RML)*(-EYT*SW) 
      FEQT7 = -SIN(4.*RML)*(.5*YT**2) 
      FEQT = FEQT1 + FEQT2 + FEQT3 + FEQT4 + FEQT5 + FEQT6 + FEQT7 
      EQT = FEQT*13751.0 
                                                                        
! equation of time in hrs:                                              
                                                                        
      EQH = EQT/3600. 
                                                                        
! convert eq of time from sec to deg                                    
                                                                        
      REQT = EQT/240. 
                                                                        
! calc right ascension in rads                                          
                                                                        
      RA = ML - REQT 
      RRA = RA*PIR 
                                                                        
! calc declination in rads, deg                                         
                                                                        
      TAB = 0.43360*SIN(RRA) 
      RDECL = ATAN(TAB) 
      DECL = RDECL/PIR 
                                                                        
! CTM.  BEGIN VECTORIZED LOOP.  THE WHOLE LOOP    SHOULD VECTORIZE.     
                                                                        
!          do 18079  kk=1,c_kmax                                        
                                                                        
                                                                        
! calc local hour angle                                                 
                                                                        
      RLT =  c_lat( kk)*PIR 
      LBGMT = 12.0 - EQT/3600. -  c_lon( kk)*24./360. 
      LZGMT = 15.0*(GMT      - LBGMT) 
      ZPT = LZGMT*PIR 
      CSZ = SIN(RLT)*SIN(RDECL) + COS(RLT)*COS(RDECL)*COS(ZPT) 
      IF(CSZ .GT. 1.) CSZ = 1. 
      IF(CSZ .LT. -1.) CSZ = -1. 
                                                                        
! ZENITH AND AZIMUTH CALCULATIONS SKIPPED.  CSZ IS RETURNED VALUE.      
                                                                        
      ZR = ACOS(CSZ) 
      ZENITH( kk) = ZR/PIR 
                                                                        
       if(c_kkw.gt.0) then 
         write(c_out,11)  c_hour, GMT     ,                             &
     &    c_lat(c_kkw),  c_lon(c_kkw), ZENITH(c_kkw)                    
       endif 
   11  format(' IN ZENITH: XHOUR GMT LAT LONG,ZEN=',8f8.2) 
                                                                        
! calc local solar azimuth                                              
                                                                        
!     CAZ = (SIN(RDECL) - SIN(RLT)*COS(ZR))/(COS(RLT)*SIN(ZR))          
!     IF(CAZ .GT. 1.) CAZ = 1.                                          
!     IF(CAZ .LT. -1.) CAZ = -1.                                        
!     RAZ = ACOS(CAZ)                                                   
!     AZMUTH = RAZ/PIR                                                  
                                                                        
18079           continue 
                                                                        
      RETURN 
      END                                           
