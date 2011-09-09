! chemmain.f/chemmainbox.f    APRIL 2007                                
!                                                                       
! TEMPORARY xcemit = 0 for expo decay test                              
                                                                        
! This is the INTERFACE from the box model                              
!  to the 'radical balance' solver for photochemical production and loss
!  ( Sillman 1991 and  2007, JGR).                                      
!                                                                       
! It includes input and output of all variables                         
!  from the main model to and from the chemistry solver                 
!  and calls the solver.                                                
!                                                                       
! This program provides a 'shell' for inserting the chemistry solver    
!   into any any larger model.  (Instructions included).                
!                                                                       
! Input:  initial species concentrations,                               
!         meteorology and time information,                             
!         chemistry mechanism information                               
!                                                                       
! Output: Updated species concentrations and wet deposition.            
!                                                                       
! Called by:   boxmain.f (driver program for box model)                 
!                                                                       
! Calls to:   chemrates   (calculates reaction rate constants)          
!             hvrates     (photolysis rate constants)                   
!             quadchem    (radical balance chemistry solver)            
!                                                                       
! ---------------------------------------------                         
! History:                                                              
!  4/07  Written by Sandy Sillman                                       
! ---------------------------------------------                         
!                                                                       
! FUTURE ADD:                                                           
!   non-steady state gas/aqueous partitioning                           
!       with input flagged as gas or soluble aerosol.                   
!   Vectorize the remaining parts of the program:  hv rates.            
!   Clean up notes for CHEMSOLVE etc.                                   
!                                                                       
! ===================================================================   
!                                                                       
! Information for inserting the chemistry solver in other programs:     
!  (1) Read chemistry data:                                             
!      To read and process the chemistry (one time only)                
!         call the following subroutines:                               
!                                                                       
!         call chemread                                                 
!         call hvread                                                   
!         call cheminit                                                 
!                                                                       
!  (2) Data files:                                                      
!                                                                       
!        REACTION.DAT = chemical mechanism                              
!        TUVGRID2     = Photolysis rate data                            
!                                                                       
!  (3) Variable declarations:                                           
!       Include the following external file                             
!        to pass data into the chemistry solver:                        
!                                                                       
!     include 'chemmech.EXT'                                            
!                                                                       
!       This contains DECLARATIONS and COMMON BLOCKS                    
!        for  variables relating to the chemistry mechanism             
!        which are read at the start (by chemread, etc.).               
!       They must be passed through the larger program to the solver.   
!       This must be included in the larger program.                    
!                                                                       
!     INCLUDE 'chemvars.EXT'   - NOT necessary in the larger program.   
!                                                                       
!       This contains declarations and common blocks                    
!        for input/output variables for the chemistry solver.           
!        It is included in the interface program (chemmain).            
!        Values from the main program are entered into                  
!        variables which are declared in this file                      
!        and thus passed through to the solver.                         
!                                                                       
!    Note:  all common-block variables in the chemistry solver          
!      begin with the letters "c_".                                     
!                                                                       
! (4) Input/ouput and calling the chemistry solver:                     
!     Use the format of the program below to enter all needed input     
!       and call the solver.                                            
! ===================================================================   
! START OF SUBROUTINE                                                   
!                                                                       
      subroutine chemmain 
!                                                                       
      IMPLICIT NONE 
                               ! Varibles from the larger box model     
      include 'boxvars.EXT' 
                               ! Mechanism parameters in the solver.    
      include 'chemmech.EXT' 
                               ! I/O variables in the chem. solver      
      INCLUDE 'chemvars.EXT' 
                                                                        
!     include 'chemlocal.EXT'  ! Local variables in the solver          
!       (chemlocal has duplicated indices with boxvars -can't use)      
!                                                                       
! LOCAL VARIABLES                                                       
!     integer kk               ! vector loop  - in boxvars AND chemlocal
!     integer ic               ! chemistry counter - in boxvars AND chem
                               ! aqueous species counter                
      integer icq 
                               ! Aqueous link to gas species counter    
      integer neq 
                                                                        
                                      ! Avogadro Number * lit/cm3       
      double precision avogadrl 
      PARAMETER (avogadrl=6.0221367D20) 
                                                                        
                                                                        
! PRELIM:  set vector variable                                          
                                                                        
               ! Index for vector arrays, =1 unless vector loop entered 
      kk=1 
                                                                        
! ---------------------------------------                               
! ENTER INPUTS FOR CHEMISTRY SOLVER                                     
! ---------------------------------------                               
! NOTE:  ARRAY NUMBERS ( 1 ) can be used in VECTORIZED VERSION.         
!                                                                       
! -------------                                                         
! BASIC INPUTS:  initial species concentrations (c_xcout)               
! -------------                                                         
!                                                                       
!  ZERO ARRAYS                                                          
      do ic=1,c_nchem2 
       c_xcin( 1 ,ic) = 0. 
       c_xcemit( 1 , ic) = 0. 
       c_xcout( 1, ic) = 0. 
      enddo 
                                                                        
! ENTER INPUT SPECIES CONCENTRATIONS:                                   
!    c_xcin(ic)  = initial concentration + emissions, molec/cm3         
      do ic=1,64!nchem1 
       c_xcin( 1 ,ic) = xr( 1 ,ic) 
      enddo 
                                                                        
! OPTIONAL - ENTER EMISSIONS FOR TIME STEP. (Default = 0.)              
!    c_xcemit(kk,ic) = OPTIONAL emissions component only, molec/cm3     
!                                                                       
!    Note:  This does not affect species concentrations.                
!     It only affects variation vs time within the time interval        
!        in the exponential decay solution.                             
!                                                                       
!    Emissions component MUST be smaller than the input c_xcin.         
!                                                                       
!      do ic=1,nchem1 
!       c_xcemit( 1 ,ic) = xremit(ic) 
                                                                        
! TEMPORARY SET =0 FOR TEST                                             
!      c_xcemit( 1 ,ic) = 0.                                            
                                                                        
!                                                                       
! WARNING IF xcemit NOT LESS THAN xcin.                                 
!!        if(c_xcemit(1,ic).gt.0 .and.                                    &
!     &          c_xcemit(1,ic).gt.       c_xcin(1,ic)) then             
!    *          c_xcemit(1,ic).gt.0.9999*c_xcin(1,ic)) then             
!          write(6,101) iit, c_tchem(ic), c_xcemit(1,ic), c_xcin(1,ic) 
!  101     format(/,'WARNING:  EMISSIONS COMPONENT IN CHEMISTRY SOLVER ',&
!     &             ' MAY BE SET TOO HIGH.',/,                           &
!     &             'time species xcemis xcin =',                        &
!     &              i4,2x,a8,2(1pe10.3))                                
!                                                                        
!        c_xcemit( 1 ,ic) = 0.9999*c_xcin( 1 ,ic) 
!                                                                        
!        endif 
!                                                                       
!      enddo 
                                                                        
                                                                        
! FLAG FOR EXPONENTIAL DECAY SOLUTION                                   
!   Exponential decay solution should be used in remote locations       
!     where emissions are zero.                                         
!                                                                       
!   It can be used in locations if emissions are known.                 
!    Emissions must be "effective" emissions, including transport       
!     (for example, in the daytime mixed layer,                         
!      several model layers are effected by emissions                   
!      on scales shorter than the simulated time interval)              
                                                                        
!       do kk=1,c_kvec               ! kk vector loop                   
!     c_lexpo(kk) = .FALSE.                                             
      c_lexpo(kk) = .FALSE.! lexpo 
!       enddo                        ! kk vector loop                   
                                                                        
! TIME STEP (sec) FOR PHOTOCHEMICAL COMPONENT OF MODEL                  
                                                                        
      c_time = time 
                                                                        
! -------------                                                         
! AMBIENT CONDITIONS THAT AFFECT CHEMISTRY:                             
! -------------                                                         
!                                                                       
!  density (molec/cm3)                                                  
                                                                        
      c_dens( 1 ) = dens( 1 ) 
!                                                                       
! temperature (K)                                                       
                                                                        
      c_temp( 1 ) = temp( 1 ) 
!                                                                       
! absolute humitidy (H2O concentration, molec cm-3                      
                                                                        
      c_h2ogas( 1 ) = xh2o 
!                                                                       
! Liquid water concentration, grams cm-3                                
      c_h2oliq( 1 ) = acqua( 1 ) 
                                                                        
! Rainout fraction during time step for aqueous species                 
                                                                        
      c_rainfr( 1 ) = rainfr( 1 ) 
                                                                        
! PROTECTION AGAINST ERRONEOUS RAINOUT FRACTION (0<=rainfr<1)           
            if(c_rainfr( 1).lt.0)  c_rainfr( 1)=0. 
            if(c_rainfr( 1).gt.0.9999)c_rainfr( 1)=0.9999 
                                                                        
! Sulfate aerosol surface area, cm2 cm-3                                
!    (used for reaction rates on aerosol surfaces)                      
                                                                        
      c_saersa( 1 ) = saersa( 1 ) 
                                                                        
! Water droplet radius, cm  (default  .001 = 10 mcm)                    
                                                                        
      c_DROPLET( 1)  = DROPLET(1) 
      c_DROPLET( 1)  = 0.001 
                                                                        
! Gas-aqueous transfer rate                                             
!    enter -1 to cancel and use calculated value from Lelieveld, 1991   
! FUTURE OPTION, NOT CURRENTLY AVAILABLE.                               
                                                                        
!     c_rgasaq(1, 1 ) = RGASAQ                                          
!     c_rgasaq(1, 1 ) = -1                                              
                                                                        
                                                                        
!  FLAG FOR NON-STEADY-STATE AQUEOUS (T=steady state)                   
! FUTURE OPTION, NOT CURRENTLY AVAILABLE.                               
                                                                        
      c_lstsaq = lstsaq 
                                                                        
! NaCl OPTION:  Set NaCl (molec/cm3) to a fixed moles/liter value       
!   from input xnacl                                                    
                                                                        
!         if(xnacl.gt.0) then 
!           ic1=0 
!           ic2=0 
!           do ic=1,nchem2 
!            if(c_tchem(ic).eq.'    NAOH') ic1=c_npequil(ic) 
!            if(c_tchem(ic).eq.'     HCL') ic2=c_npequil(ic) 
!           enddo 
!           if(ic1.gt.0.and.ic2.gt.0) then 
!            c_xcin(1,  ic1)=xnacl*c_h2oliq(1)*avogadrl 
!            c_xcin(1,  ic2)=xnacl*c_h2oliq(1)*avogadrl 
!           endif 
!         endif 
                                                                        
!                                                                       
! -------------                                                         
! TIME AND PLACE (for photolysis rates)                                 
! -------------                                                         
!                                                                       
! LATITUDE, degrees                                                     
!      c_lat( 1 ) = xlat(1) 
!                                                                       
!                                                                       
! LONGITUDE, degrees                                                    
!      c_lon( 1 ) = xlon(1) 
!                                                                       
! HOUR AT THE END OF THE TIME INTERVAL, decimal hours(0-23 EST)         
      c_hour = xhour 
!                                                                       
! DATE (format yymmdd)                                                  
!  Note:  IDATE is used in ZEN1, as YYMMDD (YY=100 for 2000)            
!   but an alternative is commented out in ZEN1:  YYYYDDD.              
!  The date is also used in the photolysis array (c_jparam)             
!   and there it must be either YYMMDD or DDD.                          
!                                                                       
      c_IDATE = IDATE 
!                                                                       
! -------------                                                         
! OTHER AMBIENT CONDITIONS THAT AFFECT PHOTOLYSIS                       
! -------------                                                         
!                                                                       
!  The photolysis input array is c_jparam. Not yet vectorized.          
!                                                                       
!     c_jparam( 1) = zenith angle, calculated by the program.           
                                      ! Temperature, entered above.     
      c_jparam(1)  = zenith                                
      c_jparam(13) = c_temp( 1 ) 
!                                                                       
! DATE for the photolysis parameter.                                    
!   Format YYMMDD (YY=100 for 2000) or DDD.                             
!   Note, this format may be different from c_IDATE, used in zenith.    
!                                                                       
                                      ! Date, entered above.            
      c_jparam(20) = float(c_IDATE) 
!                                                                       
!  ALTITUDE AT THE MIDDLE OF THE LAYER (km). (future option, kPa)       
      c_jparam( 2) = altmid( 1 ) 
!                                                                       
! OZONE COLUMN (Dobson Units)                                           
      c_jparam( 3) = DOBSON( 1 ) 
!                                                                       
! SO2 COLUMN (DU).  0.1 = 1 ppb SO2 from 0 to 1 km.                     
      c_jparam( 4) = SO2COL 
!                                                                       
! NO2 COLUMN (DU).  0.1 = 1 ppb NO2 from 0 to 1 km.                     
      c_jparam( 5) = NO2COL 
!                                                                       
! Aerosol optical depth.  (0.36 = std from Elterman, 1968; 0.76=polluted
      c_jparam( 6) = AAERX( 1 ) 
!                                                                       
! Aerosol single scattering albedo (.75-.99, .75 Logan, .99 Madronich)  
      c_jparam( 7) = AERSSA 
!                                                                       
! Surface albedo  (typically 0.1)                                       
      c_jparam( 8) = ALBEDO 
!                                                                       
! Cloud-above optical depth:                                            
!   summed optical depth for clouds ABOVE the current altitude          
      c_jparam( 9) = DEPTHA 
!                                                                       
! Cloud-below optical depth:                                            
!   summed optical depth for clouds BELOW the current altitude          
      c_jparam(10) = DEPTHB 
!                                                                       
! Cloud-above altitude (km above sea level)                             
!   weighted average altitude of clouds ABOVE the current altitude      
!     weighted by optical depth                                         
      c_jparam(11) = ALTABOVE 
!                                                                       
! Cloud-below altitude (km above sea level)                             
!   weighted average altitude of clouds BELOW the current altitude      
!     weighted by optical depth                                         
      c_jparam(12) = ALTBELOW 

!                                                                       
! ---------------------------------------                               
!  CALLS TO RUN THE CHEMISTRY SOLVER                                    
! ---------------------------------------                               
!                                                                       
!  SET REACTION RATE CONSTANTS                                          
!                                                                       
      call chemrates 
!                                                                       
! SET PHOTOLYSIS RATE CONSTANTS                                         
      call hvrates 
!                                                                       
! CALL THE DRIVER SUBROUTINE FOR THE CHEMSTRY SOLUTION                  
      call quadchem 
!                                                                       
! ---------------------------------------                               
!  TRANSFER OUTPUT FROM THE CHEMISTRY SOLVER TO THE MAIN PROGRAM        
! ---------------------------------------                               
!                                                                       
! OUTPUT SPECIES CONCENTRATIONS                                         
                                                                        
!    c_xcout(ic)  = final concentration, molec/cm3                      
      do ic=1,64 !nchem1 
       xr( 1 ,ic) = c_xcout( 1 ,ic) 
      enddo 
!                                                                       
! OUTPUT WET DEPOSITION (entered into accumulator)                      
!                                                                       
!  Wet deposition is summed here for the GAS-PHASE SPECIES,             
!    and represents wet deposition for the aqueous species              
!    linked through Henry's law and aq. equilibria                      
!    to the gas phase species or pseudo-species).                       
!    (e.g. HNO3g for HNO3(aq) and NO3-).                                
!                                                                       
!  Wet deposition (M/cm2) is equal to                                   
!      the aqueous concentration (M/liter)                              
!   times the rainfall amount (cm) with unit conversion (.01 dm2/cm2)   
!                                                                       
!   Rainfall amount (cm) = rainfr (fraction) * c_h2oliq (g/cm3)         
!                   * alt  (layer thickness in cm)                      
!                                                                       
!  This loop also identifies AQUEOUS SPECIES CONCENTRATIONS             
!                                                                       
!      do icq=1,c_nchem2 
!                              ! Index for gas spec. linked to aqueous   
!        ic = c_npequil(ic) 
!                                                                        
!                               ! Identifies ic1 as aqueous species      
!        if(ic.ne.icq) then 
!           xrwdep( ic) = xrwdep(ic) + c_xcav( 1 ,icq)                   &
!     &      *c_rainfr( 1 )  * c_h2oliq( 1 )  * alt  * 0.01              
!                                                                        
!                    !if(ic.ne.icq) then                                 
!        endif 
!                  !do icq=1,c_nchem2                                    
!      enddo 
!                                                                        
!! (Alternative loop:                                                    
!!      do ic=1,c_nchem1                                                 
!!       do neq=1,nequil(ic)                                             
!        icq = ncequil(ic,neq)                                          
!        if(icq.gt.0) then                                              
!          xrwdep( ic) = xrwdep(ic) + c_xcav( 1 ,icq)                   
!    *      *c_rainfr( 1 )  * c_h2oliq( 1 )  * alt  * 0.01              
!        endif                                                          
!       enddo                                                           
!      enddo                                                            
                                                                        
!                                                                       
! OUTPUT H+ CONCENTRATION AND DEPOSITION                                
!  Output is special because H+ is not included                         
!    in the standard array of input/output species (ic=1,nchem1)        
!                                                                       
!  c_xcwdep( 1 ,c_nhplus) = H+ wet deposition in GAS UNITS (molec/cm3)  
!  c_xcout(c_nhplus) = H+ concentration in AQUEOUS UNITS (M/lit)        
!     (c_xcout and c_xcav should be identical  for H+)                  
                                                                        
!       if(c_nhplus.gt.0)  then 
!         xrwdep(c_nhplus) = xrwdep(c_nhplus) + c_xcav( 1, c_nhplus)     &
!     &      + c_xcav( 1 ,c_nhplus)                                      &
!     &       *c_rainfr( 1 ) * c_h2oliq( 1 )  * alt  * 0.01              
!       endif 
!                                                                       
! ---------------------------------------                               
! OPTION FOR WRITTEN CHEMISTRY OUTPUT                                   
! ---------------------------------------                               
!  Add controls to write output for specified times here.               
!   (The output uses common variables from the file chemvars.EXT)       
!                                                                       
!                                                                       
!  WRITE OPTIONS:                                                       
                                                                        
! OPTION TO SKIP WRITTEN OUTPUT                                         
!             if(iit.eq.0) then                                         
! OPTION TO WRITE FOR ALL TIME STEPS                                    
!              if(iit.gt.0) then 
! SAMPLE OPTION TO WRITE ONCE PER DAY                                   
!             if(c_hour.ge.12.and.c_hour.lt.(12.+c_time/3600.) then     
                                                                        
! TIME STEP AND NUMBER OF ITERATIONS                                    
!       write(lout,1139) iit, c_iter, c_ohtest, c_notest 
! 1139  format(/,' TIME STEP =', i5, '      ITERATION #',i3              &
!     &        ,  '     OH, NO tests= ',2(1pe10.3))                      
                                                                        
! WRITE SPECIES CONCENTRATIONS AND REACTION RATES                       
!       call chemwrit(1) 
                                                                        
! WRITE DETAILED ANALYSIS FOR SPECIFIED SPECIES.                        
!        call analyze('      O3', 1)                                    
                                                                        
! END WRITTEN OUTPUT OPTION                                             
!       endif 
                                                                        
                                                                        
 2000 return 
      END                                           
