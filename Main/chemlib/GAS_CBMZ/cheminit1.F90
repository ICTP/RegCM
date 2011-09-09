!  cheminit.f   April, 2007                                             
!                                                                       
!    4-2009: error with ifort, but not with -C compile option.          
!      indices c_noh, etc. are written incorrectly, possibly related to 
!       warning message about double precision in COMMON.               
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
!                                                                       
!                                                                       
          subroutine chemread 
                                                                        
                                                                        
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
      IMPLICIT NONE 
                                                                        
      include 'chemmech.EXT' 
!     INCLUDE 'chemvars.EXT'                                            
      include 'chemlocal.EXT' 
!                                                                       
! LOCAL VARIABLES                                                       
!   DUMMY VARIABLES FOR READ                                            
                                ! dummy input character variable        
      character*8 tdum(5) 
                                ! dummy vbl to identify READ            
      character*4 titl 
                                ! dummy input integer variable          
      integer ndum(5) 
                                ! dummy input real variable             
      real xdum(5) 
                                                                        
!   OTHER LOCAL VARIABLES                                               
                                                                        
!  caspair(ic,2)  Array to identify linked pair in cascade solver;      
!                    2nd species is linked in pair chain to 1st         
!                    used to establish pair chain array (nppair)        
                                                                        
                                      ! Array to ID paired species      
      integer caspair(c_cdim, 2    ) 
                                ! Index for species number in cascade   
      integer ncas 
                                ! Integer for spec. number in caspair   
      integer ncasp 
!                                                                       
!  lcastest: Flag set to TRUE when species is read in cascade.          
                                  ! Flag to check species in cascade    
      logical lcastest(c_cdim) 
!                                                                       
                            ! Counter for number of lumped species      
      integer nlump 
                            ! General counter                           
      integer jj 
                            ! Added reaction counter (nr in chemlocal)  
      integer nnr 
                                                                        
                            ! function to return chem species index     
      integer namechem 
                                                                        
! -------------------------------------------------------------------   
!                                                                       
                       ! Set vector variable for non-vectorized case    


                       write(*,*)'READDDDD',c_rin
      kk=1 
                                                                        
! GENERAL FORMATS                                                       
! 11    format(///////,a1)                                              
   11 format(a4) 
! 11    format(a8)                                                      
   12 format(i6) 
   13 format(1pe10.3) 
                                                                        
! PRELIMINARY ZERO FOR READ                                             
         do    nr=1,c_rdim 
         c_nnpro(nr) = 0 
          do    i=1,20 
           if (i.le.2) c_reactant(nr,i) = 0 
           if(i.le.6) c_treac(i,nr) = '        ' 
           if(i.le.2.and.nr.le.61) c_treach(i,nr) = '        ' 
           if(i.le.2.and.nr.le.61)  c_henry(nr,i) = 0 
           if(i.le.3.and.nr.le.61) c_treacq(i,nr) = '        ' 
           if(i.le.3.and.nr.le.61) c_aqueous(nr,i) = 0 
           c_product(nr,i) = 0 
           c_stoich(nr,i) = 0. 
          enddo 
          do ic=1,c_cdim 
            c_prodarr(nr,ic) = 0. 
          enddo 
         enddo 
                                                                        
         do  i=1,c_cdim 
          do  ii=1,c_cdim 
           c_cascade(i,ii) = 0 
          enddo 
          do ii=1,3 
            c_lump(i,ii)=0 
          enddo 
          do ii=1,2 
           caspair(i,ii) = 0 
          enddo 
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
          do ii=3,23 
            c_nppair(i,ii) = 0 
          enddo 
         enddo 
                                                                        
! ZERO SPECIAL SPECIES INDICES                                          
         c_nh2o  =0 
         c_nco2  =0 
         c_nco   = 0 
         c_nhplus   =0 
         c_nohmin  =0 
         c_noh   =0 
         c_nho2  =0 
         c_nh2o2 =0 
         c_no3   =0 
         c_nno2  =0 
         c_nno   =0 
         c_no3   =0 
         c_nn2o5 =0 
         c_nhno3 =0 
         c_nch4  =0 
         c_nco3  =0 
                                                                        
! ZERO SPECIAL INDICES FOR PARAMETERIZED RO2-RO2 REACTIONS (CBMZ)       
         c_nnrro2 = 0 
         do nr=1, c_rdim 
          c_nrro2 = 0 
         enddo 
                                                                        
                                                                        
! READ LOOP:  Read char*4 until you find 'START'.                       
!  Then BREAK and begin first read cycle.                               
!                                                                       
       do 9005 i=1,1111 
        read(c_rin,11) titl 
        if(titl.eq.'STAR') go to 9006 
 9005  continue 
 9006  continue 
!       write(c_out,11) titl                                            
                                                                        
                                                                        
! READ PRELIMINARY INDICES:  VECTOR KKW, KMAX; NUMITER AND CONVERGE     
       read(c_rin,12) c_kkw 
       read(c_rin,12) c_kmax 
       read(c_rin,12) c_numitr 
         if(c_kkw.eq.5) write(c_out,12) c_kkw, c_kmax,c_numitr 
       read(c_rin,13) c_converge 
         if(c_kkw.eq.5) write(c_out,13) c_converge 
                                                                        
! ------------------------                                              
! READ LOOP: NAMES AND CATEGORIES FOR TRANSPORTED (INPUT/OUTPUT) SPECIES
! ------------------------                                              
       ic=0 
                                                                        
       do 9010 i=1,1111 
        read(c_rin,11) titl 
        if(titl.eq.'STAR') go to 9011 
 9010  continue 
 9011  continue 
        if(c_kkw.eq.5) write(c_out,11) titl 
                                                                        
! OPTION: READ IN GROUPS OF 1 OR 5:  ii=1,1 vs ii=1,5 THROUGHOUT        
!   (main species read here and lumped species read below).             
!   (also can change in summary write)                                  
                                                                        
       do 30 i=1,c_cdim 
       read(c_rin,29) (tdum(ii),ndum(ii),ii=1,5) 
        if(c_kkw.eq.5) write(c_out,29) (tdum(ii),ndum(ii),ii=1,5) 
! 29     format(5(7x,a4,i3,2x))                                         
   29  format(5(3x,a8,i3,2x)) 
                                                                        
! WHILE NAME NE 'END', ENTER DUMMY READ INTO TCHEM ARRAY.               
! AT 'END' BREAK AND LEAVE LOOP                                         
        do 35 ii=1,5 
         if(tdum(ii).eq.'     END'.or.tdum(ii).eq.'    END ') go to 31 
         ic=ic+1 
            if(ic.eq.c_cdim) then 
              write(c_out,901) ic 
              write(6,901) ic 
  901         format(/,'WARNING: NCHEM MAY EXCEED DIMENSION'            &
     &        ,' IN COMMON BLOCK.  NCHEM=',i4)                          
            endif 
         c_tchem(ic)=tdum(ii) 
         c_icat(ic)=ndum(ii) 
! INSERT:  ASSIGN LSTS=T IF ICAT<0.  (Note, moved below after write)    
!        c_lsts(ic) = .FALSE.                                           
!        if(c_icat(ic).lt.0) c_lsts(ic) = .TRUE.                        
!        if(c_icat(ic).lt.0) c_icat(ic) = 0-c_icat(ic)                  
!          if(c_kkw.eq.5) write(c_out,*)                                
!    *          ic, c_tchem(ic),c_icat(ic), c_lsts(ic)                  
                                                                        
   35   continue 
   30  continue 
   31  continue 
!                                                                       
! READ INSERT:  RECORD NUMBER OF INPUT   SPECIES (NCHEM1)               
                                                                        
       c_nchem1=ic 
                                                                        
       if(c_kkw.gt.0) write(c_out,49) c_nchem1 
   49  format(' TOTAL NUMBER OF CHEMICALS READ =', i5) 
                                                                        
! SUMMARY WRITE:  INPUT/OUTPUT SPECIES                                  
       write(c_out,1301) 
 1301   format(/,'TRANSPORTED (INPUT) SPECIES:  ',                      &
     &     '  (negative=steady state)'   ,/)                            
                                                                        
       ncf = int(0.2*(float(c_nchem1)-0.01))+1 
       do      nc=1,ncf 
         nc1 = 5*(nc-1) + 1 
         nc2 = nc1 + 4 
         if(nc1.gt.c_nchem1)nc1=c_nchem1 
         if(nc2.gt.c_nchem1)nc2=c_nchem1 
         write(c_out, 1302) (nn, c_tchem(nn),c_icat(nn) , nn=nc1,nc2 ) 
! 1302      format(5(i3,4x,a4,i3,2x))                                   
 1302     format(5(i3,   a8,i3,2x)) 
! 1302      format(1(i3,   a8,i3,2x))                                   
       enddo 
                                                                        
! SET CATEGORY AND STEADY STATE INDEX (after WRITE)                     
!    ASSIGN LSTS=T IF ICAT<0, MAKE CATEGORY POSITIVE.                   
!    ALSO INITIALIZE LUMP FLAG=FALSE                                    
       do ic=1,c_nchem1 
                                 ! lump accumulator flag initially F    
         c_llump(ic) = .FALSE. 
         c_lsts(ic) = .FALSE. 
         if(c_icat(ic).lt.0) c_lsts(ic) = .TRUE. 
         if(c_icat(ic).lt.0) c_icat(ic) = 0-c_icat(ic) 
!          if(c_kkw.eq.5) write(c_out,*)                                
!    *            ic, c_tchem(ic),c_icat(ic), c_lsts(ic)                
                           !do ic=1,c_nchem1                            
       enddo 
                                                                        
! RECORD INTERIM NUMBER OF TOTAL SPECIES (INPUT + INTERNAL) (NCHEM2)    
!  (This will be increased as other species are added)                  
                                                                        
       c_nchem2=c_nchem1 
                                                                        
! ------------------------                                              
! END OF CHEMICAL NAME LOOP.  TOTAL NUMBER OF CHEMICALS=NCHEM2          
!  (to be moved below, after HENRY'S LAW and AQUEOUS EQULIBRIA)         
! ------------------------                                              
                                                                        
! -------------------------------------------------------               
! READ LUMPED SPECIES                                                   
! -------------------------------------------------------               
! LUMP(IC, IC1, IC2): SPECIES FROM IC1 TO IC2 ARE LUMPED                
! INTO ACCUMULATOR (IC), AND SHARE LSTS, ICAT WITH LUMPED SUM.          
! -------------------------------------------------------               
!                                                                       
!  11/22/06                                                             
! NEW READ: List of individual lumped species read here.                
                                                                        
! NOTE: Individual lumped species are automatically added               
!       to the species list.                                            
! The name of the lump group must be already included                   
!      in the list of 'primary' (i/o) species above.                    
                                                                        
       nlump=0 
                                                                        
       do 9020 i=1,1111 
        read(c_rin,11) titl 
        if(titl.eq.'STAR') go to 9021 
 9020  continue 
 9021  continue 
           if(c_kkw.eq.5) write(c_out,11) titl 
                                                                        
       do 210 i=1,20 
!       read(c_rin,209)(tdum(ii),ii=1,3)                                
        read(c_rin,209)(tdum(ii),ii=1,1) 
           if(c_kkw.eq.5) write(c_out,209)(tdum(ii),ii=1,1) 
! 209     format(3(4x,a4,3x))                                           
  209   format(3(   a8,3x)) 
        if(tdum(1).eq.'     END'.or.tdum(1).eq.'    END ')              &
     &                                              go to 211           
                                                                        
! ADVANCE COUNTER, SET SPECIES INDICES AT ZERO.                         
                                                                        
        nlump = nlump + 1 
                                                                        
        c_lump(nlump,1) = 0 
        c_lump(nlump,2) = 0 
        c_lump(nlump,3) = 0 
!                                                                       
!                                                                       
! ENTER LUMP SPECIES NUMBER INTO LUMP ARRAY.                            
! IF LUMPED SPECIES  NAME IS UNKNOWN, ERROR AND EXIT LOOP.              
         ic=namechem(tdum( 1)) 
         if(ic.eq.0) then 
          write(c_out,219) (tdum(j),j=1,1),ii 
          write(   6,219) (tdum(j),j=1,1),ii 
  219     format(/,'MAJOR ERROR: UNKNOWN SPECIES IN LUMPED SPECIES READ'&
     &      ,/, 'SPECIES LIST = ',1(a8,4x),'UNKNOWN #',i3)              
          c_lump(nlump,1) = 0 
          c_lump(nlump,2) = 0 
          c_lump(nlump,3) = 0 
          nlump=nlump-1 
          go to 210 
         endif 
                                                                        
        c_lump(nlump, 1) = ic 
        c_llump(ic) = .TRUE. 
                                                                        
! READ LIST OF SPECIES TO BE INCLUDED IN LUMP, ADD TO SPECIES LIST.     
                                                                        
        do 40 jj=1,111 
          read(c_rin,39) (tdum(ii),ii=1,5) 
          if(c_kkw.eq.5) write(c_out,39) (tdum(ii),ii=1,5) 
! 39        format(5(7x,a4,5x))                                         
   39    format(5(3x,a8,5x)) 
                                                                        
          do 45 ii=1,5 
           if(tdum(ii).eq.'       x'.or.tdum(ii).eq.'        '          &
     &      .or.tdum(ii).eq.'      x '                                  &
     &                                               ) go to 41         
                                                                        
!  ENTER SPECIES INTO SPECIES LIST AND INTO LUMP ARRAY.                 
!     (lump2 = 1st lumped species, lump3 = last lumped species)         
                                                                        
           ic=namechem(tdum(ii)) 
           if(ic.eq.0) then 
                                                                        
             c_nchem2 = c_nchem2+1 
             if(c_nchem2.eq.c_cdim) then 
               write(c_out,901) c_nchem2 
               write(6,901) c_nchem2 
             endif 
             c_tchem(c_nchem2)=tdum(ii) 
                                                                        
!  ADD TO LUMP ARRAY.                                                   
             if(c_lump(nlump,2).eq.0) c_lump(nlump, 2) = c_nchem2 
             c_lump(nlump, 3) = c_nchem2 
                                                                        
                    !  if(ic.eq.0)                                      
           else 
                                                                        
!  TEST IF SPECIES IS ALREADY INCLUDED IN SPECIES LIST -                
!    IF LUMPED SPECIES ALREADY INCLUDED, ERROR AND EXIT LUMP LOOP       
                                                                        
            write(c_out,218) (tdum(j),j=1,1) 
            write(   6,218) (tdum(j),j=1,1) 
  218     format(/,'MAJOR ERROR: PREVIOUSLY SET SPECIES IN LUMPED LIST '&
     &      ,/, 'SPECIES = ',1(a8   )               )                   
            c_lump(nlump,1) = 0 
            c_lump(nlump,2) = 0 
            c_lump(nlump,3) = 0 
            nlump=nlump-1 
            go to 210 
                                                                        
                  !  if(ic.eq.0)                                        
         endif 
                                                                        
   45   continue 
   40  continue 
   41  continue 
                                                                        
                                                                        
! ENTER CATEGORY AND STEADY STATE FOR LUMPED SPECIES.                   
! SET LUMP FLAG= F (T for accumulator only)                             
         do 225 ic=c_lump(nlump,2),c_lump(nlump,3) 
          c_icat(ic) = c_icat(c_lump(nlump,1)) 
          c_lsts(ic) = c_lsts(c_lump(nlump,1)) 
          c_llump(ic)=.FALSE. 
  225    continue 
                                                                        
! SUMMARY WRITE: LUMPED SPECIES  SUMMARY                                
       if(nlump.eq.1) write(c_out, 1310) 
 1310  format(/, 'LUMPED SPECIES: NAME AND SPECIES CONTENTS') 
       write(c_out,1311) c_tchem(c_lump(nlump,1)) 
 1311   format(/,   a8) 
                                                                        
       ncf = int(0.2*(float(c_lump(nlump,3)-c_lump(nlump,2)+1)-0.01))+1 
       do      nc=1,ncf 
         nc1 = 5*(nc-1) + c_lump(nlump,2) 
         nc2 = nc1 + 4 
         if(nc1.gt.c_lump(nlump,3)) nc1=c_lump(nlump,3) 
         if(nc2.gt.c_lump(nlump,3)) nc2=c_lump(nlump,3) 
         write(c_out, 1302) (nn, c_tchem(nn),c_icat(nn) , nn=nc1,nc2 ) 
       enddo 
                                                                        
                                                                        
  210  continue 
  211  continue 
! -------------------------------------------------------               
! END LUMPED SPECIES READ                                               
! -------------------------------------------------------               
                                                                        
! -------------------------------------------------------               
! END LUMPED SPECIES READ                                               
! -------------------------------------------------------               
!                                                                       
! ------------------------------------                                  
! READ LOOP:  HENRY'S LAW COEFFICIENTS                                  
! ------------------------------------                                  
        nrh=0 
                                                                        
       do 9030 i=1,1111 
        read(c_rin,11) titl 
        if(titl.eq.'STAR') go to 9031 
 9030  continue 
 9031  continue 
         if(c_kkw.eq.5) write(c_out,11) titl 
                                                                        
! LOOP FOR HENRY'S LAW READS: CONTINUE UNTIL 'END'.                     
      do  50 i=1,c_cdim 
                                                                        
! ADVANCE COUNTER.                                                      
       nrh = nrh+1 
! --------------------------------                                      
! READ NAMES OF HENRY'S LAW SPECIES.                                    
! --------------------------------                                      
            if(nrh.eq.161.or.nrh.eq.121.or.nrh.eq.c_cdim) then 
              write(c_out,903) ic 
              write(6,903) ic 
            endif 
  903       format(/,'WARNING: NHENRY MAY EXCEED DIMENSION'             &
     &      ,' IN COMMON BLOCK.  NHENRY=',i4)                           
                                                                        
        read(c_rin, 59) tdum(1),tdum(2) 
            if(c_kkw.eq.5) write(c_out, 59) tdum(1),tdum(2) 
! 59      format(3(4x,a4,3x))                                           
   59   format(3(   a8,3x)) 
                                                                        
!                                                                       
! ENTER NAMES INTO HENRY'S LAW MATRIX                                   
! 'END' BREAKS THE LOOP TO END THE HENRY'S LAW READ.                    
         do  65 iii=1,2 
          if(tdum(iii) .eq.'     END'.or.tdum(iii).eq.'    END ')       &
     &                                                  go to  51       
                                                                        
             c_treach(iii,nrh)=tdum(iii) 
                                                                        
          ic = namechem(tdum(iii)) 
!            write(c_out,*) nrh,iii, tdum(iii),ic                       
                                                                        
! ADD AQUEOUS HENRY'S LAW SPECIES INTO SPECIES LIST (if not already ther
          if(ic.eq.0.and.iii.eq.2) then 
           ic     = c_nchem2+1 
           c_nchem2 = ic 
           if(c_nchem2.eq.c_cdim) then 
              write(c_out,901) c_nchem2 
              write(6,901) c_nchem2 
            endif 
           c_tchem(c_nchem2)=tdum(iii) 
                             !if(ic.eq.0.and.iii.eq.2) then             
          endif 
                                                                        
! ENTER SPECIES INTO HENRY LAW ARRAY                                    
          if(ic.gt.0) then 
                  c_henry(nrh,iii)=ic 
          else 
! TEST FOR POSSIBLE ERROR IN CHEMICAL NAMES.                            
!  THE GAS-PHASE SPECIES SHOULD ALREADY BE IN THE SPECIES LIST.         
            if(                        tdum(iii).ne.'        '.and.     &
     &         tdum(iii).ne.'       x'.and.tdum(iii).ne.'       X'.and. &
     &         tdum(iii).ne.'      x '.and.tdum(iii).ne.'      X '.and. &
     &         tdum(iii).ne.'    MORE')  then                           
                 write(c_out, 69) tdum(iii), nrh,iii,(tdum(j),j=1,2) 
                 write(   6, 69) tdum(iii), nrh,iii,(tdum(j),j=1,2) 
            endif 
   69      format(/,'WARNING:  UNRECOGNIZED CHEMICAL NAME:',a8,/,       &
     &      ' NRH=',i4,'  NAME # ',i4,                                  &
     &      '  HENRY:  ',a8,'=',a8)                                     
          endif 
   65    continue 
                                                                        
! DEFAULT ACCOMODATION COEFFICIENT AND MOLECULAR WEIGHT                 
       c_accom(nrh) = 0.05 
       c_molwt(nrh) = 30. 
                                                                        
!  READ RATE PARAMETERS                                                 
!      write(c_out,*) nrh                                               
                                                                        
! CURRENT OPTION: ACCOMODATION COEFFICIENT AND MOLECULAR WT READ.       
! FUTURE MODIFICATION: ADD ALTERNATIVE RATE PARAMETER FORMATS IF NEEDED.
!                                                                       
!  READ LINES WITH THE FOLLOWING  FORMAT:                               
!    CO2 =     CO2L                                                     
!   1    0 3.400E-02     2400. 5.000E-02       44.                      
                                                                        
       read(c_rin, 33)   c_nrkh(nrh),  c_rkh(1,nrh),c_rkh(2,nrh),       &
     &    c_accom(nrh), c_molwt(nrh)                                    
          if(c_kkw.eq.5) write(c_out, 33)                               &
     &       c_nrkh(nrh), c_rkh(1,nrh),c_rkh(2,nrh),c_accom(nrh)        &
     &         , c_molwt(nrh)                                           
                                                                        
   33  format(5x,i5, 1pe10.3, 0pf10.0, 1pe10.3, 0pf10.0) 
                                                                        
   50 continue 
   51 continue 
! --------------------------------------------------------              
! END HENRY'S LAW LOOP.                                                 
! --------------------------------------------------------              
       c_nreach = nrh-1 
                                                                        
       if(c_kkw.gt.0) write(c_out,57) c_nreach 
   57  format(' TOTAL NUMBER OF HENRYS-LAW COEFFICIENTS READ =', i5) 
!                                                                       
!                                                                       
! ------------------------------------                                  
! READ LOOP:  AQUEUS EQUILIBRIUM CONSTANTS WITH H+                      
! ------------------------------------                                  
        nrq=0 
                                                                        
       do 9040 i=1,1111 
        read(c_rin,11) titl 
        if(titl.eq.'STAR') go to 9041 
 9040  continue 
 9041  continue 
         if(c_kkw.eq.5) write(c_out,11) titl 
                                                                        
! LOOP FOR EQUILIBRIUM READS: CONTINUE UNTIL 'END'.                     
      do  70 i=1,c_cdim 
                                                                        
! ADVANCE COUNTER.                                                      
       nrq = nrq+1 
                                                                        
            if(nrq.eq. 161.or.nrq.eq.121.or.nrq.eq. c_cdim)             &
     &      write(c_out,904) nrq                                        
  904       format(/,'WARNING: NAQUEOUS MAY EXCEED DIMENSION'           &
     &      ,' IN COMMON BLOCK.  NAQUEOUS=',i4)                         
                                                                        
! --------------------------------                                      
! READ NAMES OF H+ EQUILIBRIUM SPECIES:  A = B + C* (*C=H+or OH-)       
!    NOTE, 1st EQUILIBRIUM MUST BE 0<=>H+ + OH-                         
! --------------------------------                                      
                                                                        
!         write(c_out,*) nrq                                            
        read(c_rin, 59) tdum(1),tdum(2), tdum(3) 
            if(c_kkw.eq.5) write(c_out, 59) tdum(1),tdum(2), tdum(3) 
                                                                        
! SWITCH 2nd AND 3rd SPECIES if 2nd SPECIES = H+ or OH-                 
        if (nrq.gt.1) then 
          if (tdum(2).eq.c_treacq(2,1).or.tdum(2).eq.c_treacq(3,1)      &
     &        .or.tdum(2).eq.'      H+'.or.tdum(2).eq.'     OH-'        &
     &        .or.tdum(2).eq.'     H+ '.or.tdum(2).eq.'    OH- '        &
     &                                                     ) then       
            if(                        tdum( 3 ).ne.'        '.and.     &
     &         tdum( 3 ).ne.'       x'.and.tdum( 3 ).ne.'       X'.and. &
     &         tdum( 3 ).ne.'      x '.and.tdum( 3 ).ne.'      X '.and. &
     &         tdum( 3 ).ne.'    MORE')  then                           
                                                                        
              tdum(4) = tdum(2) 
              tdum(2) = tdum(3) 
              tdum(3) = tdum(4) 
              if(c_kkw.gt.0) write(c_out,88)                            &
     &                nrq, tdum(1),tdum(2),tdum(3)                      
   88         format(/,'AQUEOUS ORDER SWITCH: NUMBER:',i3,'  REACTION:',&
     &               a8,'=',a8,'&',a8)                                  
                                                                        
                      !if (tdum(3).ne.'        ')                       
            endif 
                    !if (tdum(2).eq. ...H+)                             
          endif 
                             !if (nrq.gt.1) then                        
        endif 
                                                                        
!                                                                       
! ENTER NAMES INTO AQUEOUS EQUILIBRIUM MATRIX                           
! 'END' BREAKS THE LOOP TO END THE EQUILIBRIUM READ.                    
         do  85 iii=1,3 
                                                                        
          if(tdum(iii) .eq.'     END'.or.tdum(iii).eq.'    END ')       &
     &                                                   go to  71      
                                                                        
          c_treacq(iii,nrq)=tdum(iii) 
                                                                        
          ic = namechem(tdum(iii)) 
!            write(c_out,*) nrq,iii, tdum(iii),ic                       
                                                                        
! ENTER NEW CATION/ANION INTO SPECIES LIST IF NOT ALREADY INCLUDED.     
!   (NOTE, FIRST EQUILIBRIUM SHOULD BE 0 <=> H+ + OH-;                  
!    SUBSEQUENT SHOULD BE A <=> B + H+ or OH-)                          
                                                                        
         if(ic.eq.0.and.(iii.eq.2.or.(iii.eq.3.and.nrq.eq.1))) then 
           ic     = c_nchem2+1 
           c_nchem2 = ic 
           if(c_nchem2.eq.c_cdim) then 
              write(c_out,901) c_nchem2 
              write(6,901) c_nchem2 
            endif 
           c_tchem(c_nchem2)=tdum(iii) 
                                                                        
                       !if(ic.eq.0.and.(iii.eq.2....))                  
         endif 
                                                                        
                                                                        
! ENTER SPECIES INTO AQUEOUS EQUILIBRIUM ARRAY                          
          if(ic.gt.0) then 
                c_aqueous(nrq,iii)=ic 
          else 
! TEST FOR POSSIBLE ERROR IN CHEMICAL NAMES                             
            if(                        tdum(iii).ne.'        '.and.     &
     &         tdum(iii).ne.'       x'.and.tdum(iii).ne.'       X'.and. &
     &         tdum(iii).ne.'      x '.and.tdum(iii).ne.'      X '.and. &
     &         tdum(iii).ne.'    MORE')  then                           
                 write(c_out, 89) tdum(iii), nrq,iii,(tdum(j),j=1,3) 
                 write(   6, 89) tdum(iii), nrq,iii,(tdum(j),j=1,3) 
            endif 
   89       format(/,'WARNING:  UNRECOGNIZED CHEMICAL NAME:',a8,/,      &
     &      ' NRQ=',i4,'  NAME # ',i4,                                  &
     &      ' AQUEOUS:  ',a8,'=',a8,' + ',a8)                           
          endif 
   85    continue 
                                                                        
!  END NAMEREAD LOOP.  READ RATE PARAMETERS                             
                                                                        
       read(c_rin, 33)   c_nrkq(nrq), c_rkq(1,nrq),c_rkq(2,nrq) 
          if(c_kkw.eq.5) write(c_out, 33)                               &
     &              c_nrkq(nrq), c_rkq(1,nrq),c_rkq(2,nrq)              
                                                                        
   70 continue 
   71 continue 
                                                                        
       c_nreacq = nrq-1 
                                                                        
       if(c_kkw.gt.0) write(c_out,79) c_nreacq 
   79  format(' TOTAL NUMBER OF AQUEOUS EQUILIBRIUM CONSTANTS =', i5) 
                                                                        
!                                                                       
! ------------------------------------                                  
! READ LOOP:  AQUEUS SPECIAL EQUILIBRIUM CONSTANTS WITHOUT H+/OH-       
! ------------------------------------                                  
! (note:  These are special.  A <->B+C; the combined species A are set  
!  in equilibria to be automatically solved when B is called.           
!  Species C MUST have much higher concentrations then either A or B,   
!  or else be zero.)                                                    
                                                                        
        nrqq=0 
                                                                        
       do 9045 i=1,1111 
        read(c_rin,11) titl 
        if(titl.eq.'STAR') go to 9046 
 9045  continue 
 9046  continue 
         if(c_kkw.eq.5) write(c_out,11) titl 
                                                                        
! LOOP FOR EQUILIBRIUM READS: CONTINUE UNTIL 'END'.                     
      do  170 i=1,c_cdim 
                                                                        
! ADVANCE COUNTER.                                                      
       nrqq = nrqq+1 
                                                                        
            if(nrqq.eq.121.or.nrqq.eq. 91.or.nrqq.eq.c_cdim)            &
     &      write(c_out,905) nrqq                                       
  905       format(/,'WARNING: NAQUEOUSQ MAY EXCEED DIMENSION'          &
     &      ,' IN COMMON BLOCK.  NAQUEOUSQ=',i4)                        
                                                                        
! --------------------------------                                      
! READ NAMES OF SPECIAL EQUILIBRIUM SPECIES:  A = B + C* (*C=common spec
! --------------------------------                                      
                                                                        
!         write(c_out,*) nrqq                                           
        read(c_rin, 59) tdum(1),tdum(2), tdum(3) 
            if(c_kkw.eq.5) write(c_out, 59) tdum(1),tdum(2), tdum(3) 
!                                                                       
! ENTER NAMES INTO AQUEOUS SPECIAL EQUILIBRIUM MATRIX                   
! 'END' BREAKS THE LOOP TO END THE EQUILIBRIUM READ.                    
         do  86 iii=1,3 
                                                                        
          if(tdum(iii) .eq.'     END'.or.tdum(iii).eq.'    END ')       &
     &                                                  go to  171      
                                                                        
          c_treaqq(iii,nrqq)=tdum(iii) 
                                                                        
          ic = namechem(tdum(iii)) 
!            write(c_out,*) nrqq,iii, tdum(iii),ic                      
          if(ic.gt.0) then 
                c_aqspec(nrqq,iii)=ic 
          else 
! TEST FOR POSSIBLE ERROR IN CHEMICAL NAMES                             
            if(                        tdum(iii).ne.'        '.and.     &
     &         tdum(iii).ne.'       x'.and.tdum(iii).ne.'       X'.and. &
     &         tdum(iii).ne.'      x '.and.tdum(iii).ne.'      X '.and. &
     &         tdum(iii).ne.'    MORE')  then                           
                 write(c_out, 87) tdum(iii), nrqq,iii,(tdum(j),j=1,3) 
                 write(   6, 87) tdum(iii), nrqq,iii,(tdum(j),j=1,3) 
            endif 
   87       format(/,'WARNING:  UNRECOGNIZED CHEMICAL NAME:',a8,/,      &
     &      ' NRQQ=',i4,'  NAME # ',i4,                                 &
     &      ' AQUEOUSQ:  ',a8,'=',a8,' + ',a8)                          
          endif 
   86    continue 
                                                                        
!  END NAMEREAD LOOP.  READ RATE PARAMETERS                             
!   NOTE FUTURE DO OPTION: Make rate=rate of return reaction            
!     with nn corresponding to nn for regular reactions (   nrk)        
                                                                        
       read(c_rin, 33)   c_nrkqq(nrqq),  c_rkqq(1,nrqq),c_rkqq(2,nrqq) 
          if(c_kkw.eq.5) write(c_out, 33)                               &
     &          c_nrkqq(nrqq), c_rkqq(1,nrqq),c_rkqq(2,nrqq)            
                                                                        
  170  continue 
  171  continue 
                                                                        
       c_nreaqq = nrqq-1 
                                                                        
       if(c_kkw.gt.0) write(c_out,99) c_nreaqq 
   99  format(' TOTAL NUMBER OF AQUEOUSQ EQUILIBRIUM CONSTANTS =', i5) 
                                                                        
                                                                        
! --------------------------------------                                
! READ LOOP:  REACTION RATES; REACTANTS, PRODUCTS AND STOICHIOMETRIES   
! --------------------------------------                                
       nr = 0 
                                                                        
       do 9050 i=1,1111 
        read(c_rin,11) titl 
        if(titl.eq.'STAR') go to 9051 
 9050  continue 
 9051  continue 
         if(c_kkw.eq.5) write(c_out,11) titl 
                                                                        
! LOOP FOR REACTION READS: CONTINUE UNTIL 'END'.                        
!     do 100 i=1,3211                                                   
      do 100 i=1,c_rdim 
                                                                        
! ADVANCE REACTION COUNTER                                              
       nr = nr+1 
!      write(c_out,*) nr                                                
                                                                        
            if(nr.eq.1711.or.nr.eq.1911.or.nr.eq.3211)                  &
     &      write(c_out,906) nr                                         
  906       format(/,'WARNING: NREACT   MAY EXCEED DIMENSION'           &
     &      ,' IN COMMON BLOCK.  NREACT  =',i4)                         
                                                                        
! --------------------------------------------------------              
! READ NAMES OF REACTANTS, PRODUCTS AND STOICHIOMETRIES.                
! --------------------------------------------------------              
       do 110 ii=1,10 
!         write(c_out,*) ii                                             
        read(c_rin,109) tdum(1),tdum(2),xdum(3),tdum(3),                &
     &   xdum(4),tdum(4),xdum(5),tdum(5)                                
          if(c_kkw.eq.5) write(c_out,109) tdum(1),tdum(2),xdum(3),      &
     &     tdum(3),xdum(4),tdum(4),xdum(5),tdum(5)                      
! 109     format(4x,a4,7x,a4,4x,2(f6.3,6x,a4,3x),f6.3,6x,a4)            
  109   format(   a8,3x,a8,4x,2(f6.3,2x,a8,3x),f6.3,2x,a8) 
!                                                                       
! ENTER NAMES INTO REACTANT OR PRODUCT MATRIX. (also TREAC)             
! 'NEXT' BREAKS THE LOOP TO ADVANCE TO NEXT REACTION                    
! 'END' BREAKS THE LOOP TO END THE REACTION READ.                       
         do 120 iii=1,5 
          if(ii.eq.1.and.tdum(iii).ne.'    NEXT'.and.                   &
     &                   tdum(iii).ne.'     END'.and.                   &
     &                   tdum(iii).ne.'   NEXT '.and.                   &
     &                   tdum(iii).ne.'    END '                        &
     &                                                     ) then       
             c_treac(iii,nr)=tdum(iii) 
          endif 
                                                                        
          if(tdum(iii).eq.'    NEXT'.or.tdum(iii).eq.'   NEXT '         &
     &                                                   ) go to 111    
          if(tdum(iii) .eq.'     END'.or.tdum(iii).eq.'    END ')       &
     &                                                     go to 101    
                                                                        
          ic = namechem(tdum(iii)) 
!            write(c_out,*) nr, tdum(iii),c_tchem(1),ic                 
          if(ic.gt.0) then 
            if(iii.le.2) then 
               c_reactant(nr,iii)=ic 
!              if(c_kkw.eq.5) write(c_out,*) nr,iii,c_reactant(nr,iii)  
            else 
              do 125 j=1,20 
               if (c_product(nr,j).eq.0) go to 126 
  125         continue 
  126         continue 
              c_product(nr,j) = ic 
              c_stoich(nr,j) = xdum(iii) 
              c_nnpro(nr) = j 
              c_prodarr(nr,ic) = xdum(iii) 
!             if(c_kkw.eq.5) write(c_out,*) nr,j,c_nnpro(nr),           
!    *          c_product(nr,j)                                         
                                                                        
! TEST WRITE PRODARR                                                    
!  (=stoichiometry for production of species ic from reaction nr)       
!                                                                       
            if(c_kkw.gt.0) write(c_out,127) nr, c_rk(1,nr), tdum 
            if(c_kkw.gt.0) write(c_out,128)                             &
     &         j, ic, c_stoich(nr,j), c_prodarr(nr,ic)                  
  127       format(' ID PRODUCTS:  NR, RK1=', i5,1pe10.3,/,             &
     &      '  TREACT= ',a8,'+',a8,'->',a8,'+',a8,'+',a8)               
  128       format(' J IC STOIC PRODARR = ', 2i5,2f10.3) 
                                                                        
            endif 
          else 
!                                                                       
! FLAG FOR HV REACTION:  Identified by 'hv' as Reactiant #2.            
!   Set REACTANT(2) =-1,                                                
!    and default hv rate index 4 (rate proportional to jNO2)            
!                                                                       
            if(iii.eq.2.and.                                            &
     &          (tdum(iii).eq.'      hv'.or.tdum(iii).eq.'     hv ')    &
     &                                         ) then                   
               c_reactant(nr,iii)=-1 
               if(c_nrk(nr).eq.0) c_nrk(nr) =  4 
                                                                        
!  TEST WRITE HV                                                        
!           write(c_out,138) nr, c_rk(1,nr), tdum                       
! 138         format(' ID HV RATES:  NR, RK1=', i5,1pe10.3,/,           
!    *      '  TREACT= ',a8,'+',a8,'->',a8,'+',a8,'+',a8)               
                                                                        
                        ! FLAG FOR HV                                   
            endif 
!                                                                       
! TEST FOR POSSIBLE ERROR IN CHEMICAL NAMES                             
            if(tdum(iii).ne.'      hv'.and.tdum(iii).ne.'        '.and. &
     &         tdum(iii).ne.'       x'.and.tdum(iii).ne.'       X'.and. &
     &         tdum(iii).ne.'     hv '.and.tdum(iii).ne.'      x '.and. &
     &         tdum(iii).ne.'    MORE')  then                           
                 write(c_out,139) tdum(iii), nr,iii,tdum 
                 write(   6,139) tdum(iii), nr,iii,tdum 
            endif 
  139       format(/,'WARNING:  UNRECOGNIZED CHEMICAL NAME:',a8,/,      &
     &      ' NR=',i4,'  NAME # ',i4,                                   &
     &      '  TREACT= ',a8,'+',a8,'->',a8,'+',a8,'+',a8)               
          endif 
  120    continue 
                                                                        
  110  continue 
  111  continue 
! -------------                                                         
! END READ REACTION NAMES                                               
! -------------                                                         
                                                                        
! --------------------------------------------------------              
!  READ RATE PARAMETERS (moved after NAMES)                             
! --------------------------------------------------------              
                                                                        
!      read(c_rin, 33) c_rk(1,nr),c_rk(2,nr), c_nrk(nr)                 
       read(c_rin,131) ii, c_nrk(nr) 
          if(c_kkw.eq.5) write(c_out,131) ii, c_nrk(nr) 
  131  format(2i5) 
                                                                        
! READ OPTIONS AND FLAG FOR READ                                        
       if(c_nrk(nr).eq.0.or.c_nrk(nr).eq.-1.or.c_nrk(nr).eq.-20) then 
          read(c_rin,132) c_rk(1,nr), c_rk(2,nr) 
            if(c_kkw.eq.5) write(c_out,132) c_rk(1,nr), c_rk(2,nr) 
       endif 
! 132     format(10x,1pe10.3,0pf8.0,0pf6.1)                             
  132   format(10x,1pe10.3,0pf8.0,0pf8.3) 
                                                                        
       if( c_nrk(nr).eq.-2.or.c_nrk(nr).eq.-5 .or.c_nrk(nr).eq.-6) then 
          read(c_rin,132) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr) 
            if(c_kkw.eq.5) write(c_out,132)                             &
     &                 c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)               
       endif 
                                                                        
       if(c_nrk(nr).eq.-3) then 
          read(c_rin,133) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)            &
     &                    ,c_rk(4,nr), c_rk(5,nr)                       
            if(c_kkw.eq.5) write(c_out,133)                             &
     &                 c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)               &
     &                    ,c_rk(4,nr), c_rk(5,nr)                       
       endif 
  133   format(10x,0pf6.2,1pe10.3,0pf6.1,1pe10.3,0pf6.1) 
                                                                        
       if(c_nrk(nr).eq.-4) then 
          read(c_rin,134) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)            &
     &                    ,c_rk(4,nr), c_rk(5,nr),c_rk(6,nr), c_rk(7,nr)
            if(c_kkw.eq.5) write(c_out,134)                             &
     &                 c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)               &
     &                    ,c_rk(4,nr), c_rk(5,nr),c_rk(6,nr), c_rk(7,nr)
       endif 
  134   format(10x,1pe10.3,0pf8.0,0pf6.2,1pe10.3,0pf6.1,1pe10.3,0pf6.1) 
                                                                        
       if(c_nrk(nr).eq.-7) then 
          read(c_rin,135) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)            &
     &                    ,c_rk(4,nr), c_rk(5,nr),c_rk(6,nr)            
            if(c_kkw.eq.5) write(c_out,135)                             &
     &                 c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)               &
     &                    ,c_rk(4,nr), c_rk(5,nr),c_rk(6,nr)            
       endif 
  135   format(10x,1pe10.3,0pf8.0,4(1pe10.3)) 
                                                                        
       if(c_nrk(nr).eq.-8) then 
          read(c_rin,136) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)            &
     &                    ,c_rk(4,nr), c_rk(5,nr)                       
            if(c_kkw.eq.5) write(c_out,136)                             &
     &                 c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)               &
     &                    ,c_rk(4,nr), c_rk(5,nr)                       
       endif 
                                               ! format used for others 
  136   format(10x,3(1pe10.3,0pf8.0) ) 
                                                                        
       if(c_nrk(nr).eq.-9.or.c_nrk(nr).eq.-10                 ) then 
          read(c_rin,136) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)            &
     &                    ,c_rk(4,nr), c_rk(5,nr),c_rk(6,nr)            
            if(c_kkw.eq.5) write(c_out,136)                             &
     &                 c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)               &
     &                    ,c_rk(4,nr), c_rk(5,nr),c_rk(6,nr)            
       endif 
                                                                        
       if(c_nrk(nr).eq.-11) then 
          read(c_rin,137) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)            &
     &                    ,c_rk(4,nr), c_rk(5,nr)                       
            if(c_kkw.eq.5) write(c_out,137)                             &
     &                 c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)               &
     &                    ,c_rk(4,nr), c_rk(5,nr)                       
       endif 
  137   format(10x,1pe10.3,0pf8.0,1pe10.3,0pf8.3,0pf8.0) 
                                                                        
       if(c_nrk(nr).eq.-12.or.c_nrk(nr).eq.-13) then 
          read(c_rin,136) c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)            &
     &                    ,c_rk(4,nr)                                   
            if(c_kkw.eq.5) write(c_out,136)                             &
     &                 c_rk(1,nr), c_rk(2,nr), c_rk(3,nr)               &
     &                    ,c_rk(4,nr), c_rk(5,nr),c_rk(6,nr)            
       endif 
                                                                        
! 138     format(10x,3(1pe10.3))                                        
                                                                        
!  Addition for secondary organic aerosols (Luis Olcese)                
       if(c_nrk(nr).eq.-21) then 
        read(c_rin,142) c_rk(1,nr),c_rk(2,nr),c_rk(3,nr),c_rk(4,nr) 
        if(c_kkw.eq.5) write(c_out,142) c_rk(1,nr),c_rk(2,nr),          &
     &   c_rk(3,nr),c_rk(4,nr)                                          
       endif 
                                                                        
       if(c_nrk(nr).eq.-22) then 
        read(c_rin,142) c_rk(1,nr),c_rk(2,nr),c_rk(3,nr),c_rk(4,nr) 
        if(c_kkw.eq.5) write(c_out,142) c_rk(1,nr),c_rk(2,nr),          &
     &   c_rk(3,nr),c_rk(4,nr)                                          
       endif 
                                                                        
  142  format(10x,1pe10.3,1x,1pe10.3,1x,1pe10.3,1x,1pe10.3) 
                                                                        
                                                                        
! READ FOR hv REACTIONS                                                 
       if(c_nrk(nr).gt.0.or.c_nrk(nr).le.-30  ) then 
          read(c_rin,132) c_rk(1,nr) 
            if(c_kkw.eq.5) write(c_out,132) c_rk(1,nr) 
       endif 
                                                                        
! COUNT NUMBER OF PARAMETERIZED RO2 REACTIONS AND RECORD REACTION NUMBER
      if(c_nrk(nr).eq.-13) then 
        c_nnrro2 = c_nnrro2 + 1 
        c_nrro2(c_nnrro2) = nr 
      endif 
                                                                        
  100 continue 
  101 continue 
! --------------------------------------------------------              
! END REACTION LOOP.                                                    
! --------------------------------------------------------              
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
       do 9070 i=1,1111 
        read(c_rin,11) titl 
        if(titl.eq.'STAR') go to 9071 
 9070  continue 
 9071  continue 
                                                                        
       ncas = 1 
       ncasp = 0 
       nc = 0 
                                                                        
                                       ! CASCADE READ LOOP              
       do  i=1,c_cdim 
        read(c_rin,309)(tdum(ii),ii=1,3) 
! 309     format(3(4x,a4,2x))                                           
  309   format(3(   a8,2x)) 
                                                                        
! LOOP TO READ THREE SPECIES AND CONVERT TO SPECIES NUMBERS             
!   NOTE:  This read loop repeats to read full group                    
!          until it reads a break mark (x)                              
                                                                        
                            ! LOOP TO READ GROUP OF THREE SPECIES (320) 
        do ii=1,3 
                                                                        
! IF NAME IS 'END', BREAK THE CASCADE READ AND EXIT                     
!     (ALSO ADJUST COUNTER)                                             
                                                                        
         if(tdum(1).eq.'     END'.or.tdum(1).eq.'    END ') then 
           if(nc .eq.0) ncas = ncas - 1 
           go to 311 
                        !if(tdum(1).eq.'     END'....                   
         endif 
                                                                        
! IF NAME IS  x, BREAK, PROCESS, AND GO TO THE NEXT CASCADE  GROUP.     
                                                                        
          if(tdum(ii).eq.'       x'.or.tdum(ii).eq.'        '           &
     &      .or. tdum(ii).eq.'     pre'.or.tdum(ii).eq.'    post'       &
     &                                               ) go to 310        
                                                                        
! IDENTIFY SPECIES NAME                                                 
         ic=namechem(tdum(ii)) 
                                                                        
! IF FIRST NAME IS 'pair' , SET FLAG TO ESTABLISH CASCADE PAIRS         
         if (ii.eq.1.and.                                               &
     &        (tdum(ii).eq.'    pair'                                   &
     &         .or.tdum(ii).eq.'    PAIR'                               &
     &         .or.tdum(ii).eq.'   pair '                               &
     &         .or.tdum(ii).eq.'   PAIR '                               &
     &                                 ) ) ic=-1                        
                                                                        
! WARNING FOR NO-NAME IN CASCADE. (no-name is regarded as break)        
         if(ic.eq.0) then 
          write(c_out,319) (tdum(j),j=1,3),ii 
          write(   6,319) (tdum(j),j=1,3),ii 
  319     format(/,'MAJOR ERROR: UNKNOWN SPECIES IN CASCADE SPEC  READ.'&
     &      ,/, 'SPECIES LIST = ',3(a8,4x),'UNKNOWN #',i3)              
          go to 310 
         endif 
                                                                        
! ENTER SPECIES IN CASCADE LIST (AND ADVANCE COUNTER)                   
!  (NOTE: LATER ADJUSTMENT FOR CASCADE PAIR FLAG)                       
        nc = nc + 1 
        c_cascade(ncas ,nc ) = ic 
                                                                        
                                                                        
                           !do ii=1,3  LOOP TO READ GROUP OF THREE SPECI
        enddo 
                                                                        
                           ! AFTER READ OF COMPLETE CASCADE GROUP.      
  310  continue 
                                                                        
! PROCESS CASCADE SPECIES GROUP  AFTER READ OF COMPLETE GROUP.          
                                                                        
! ESTABLISH CASCADE PAIRS BASED ON READ GROUP OF THREE SPECIES.         
                                                                        
! (1) IF CASCADE LIST CONSISTS OF JUST TWO SPECIES,                     
!     THEN IDENTIFY SECOND SPECIES AS PAIRED WITH THE FIRST             
!     AND REMOVE THE SECOND SPECIES  FROM THE CASCADE ARRAY             
!    (paired species will automatically be solved when primary is called
                                                                        
      if (c_cascade(ncas,1).gt.0.and.                                   &
     &     c_cascade(ncas,2).gt.0.and.c_cascade(ncas,3).eq.0) then      
                                                                        
         ncasp = ncasp + 1 
         caspair(ncasp ,1) = c_cascade(ncas,1) 
         caspair(ncasp ,2) = c_cascade(ncas,2) 
         c_cascade(ncas,2) = 0 
                                                                        
               !if (c_cascade(ncas,2).gt.0....                          
      endif 
                                                                        
! (2) IF CASCADE LIST HAS 'PAIR' FLAG (first species is -1)             
!     THEN IDENTIFY 2nd AND 3rd SPECIES AS PAIR                         
!     AND REMOVE BOTH FROM CASCADE ARRAY.                               
!      (2nd species should be in cascade or in cascade pair elsewhere)  
!    (with 1st species set to zero, this cascade number will be re-done)
                                                                        
      if(c_cascade(ncas,1).eq.-1                                        &
     &    .and.c_cascade(ncas,2).gt.0.and.c_cascade(ncas,3).gt.0) then  
                                                                        
         ncasp = ncasp + 1 
         caspair(ncasp ,1) = c_cascade(ncas,2) 
         caspair(ncasp ,2) = c_cascade(ncas,3) 
                                                                        
         c_cascade(ncas,1) = 0 
         c_cascade(ncas,2) = 0 
         c_cascade(ncas,3) = 0 
                                                                        
!   ADD ERROR WRITE if species is zero                                  
                                                                        
                !if(c_cascade(ncas,1).eq.-1) then                       
      endif 
                                                                        
!  ADVANCE THE SPECIES COUNTER UNLESS 1st SPECIES IS ZERO               
!    (note: do not advance if the line was just to identify a pair;     
!      or if it is a mistake.)                                          
                                                                        
              if(c_cascade(ncas,1).gt.0) ncas = ncas + 1 
              nc = 0 
                                                                        
                     !do  i=1,c_cdim  CASCADE READ LOOP                 
       enddo 
  311  continue 
                                                                        
! TEMPORARY TEST:  LIST CASCADE PAIRS AND LIST CASCADE                  
       if(c_kkw.eq.5) then 
         write(c_out,314) 
  314    format(/,' LIST OF CASCADE PAIRS AND CASCADE') 
         do i=1,ncasp 
           write(c_out, 315) i, c_tchem(caspair(i,1)),                  &
     &                      c_tchem(caspair(i,2))                       
  315      format(i4,   a8,   a8,   a8) 
         enddo 
         do i=1,  ncas 
           do ii=1,3 
             ic = c_cascade(i,ii) 
             if(ic.gt.0) then 
               write(c_out,315) i, c_tchem(ic) 
             endif 
           enddo 
         enddo 
                       !if(c_kkw.eq.5) then                             
       endif 
                                                                        
                                                                        
! -------------------------------------------------------               
! END CASCADE READ                                                      
! -------------------------------------------------------               
!                                                                       
! --------------------------                                            
! IDENTIFY SPECIAL SPECIES NUMBERS AND CATEGORIES                       
! --------------------------                                            
!  This replaces 'family'.                                              
                                                                        
! LOOP TO IDENTIFY SPECIAL SPECIES                                      
      do ic=1, c_nchem2 
                                                                        
! TEMPORARY DEBUG WRITE                                                 
!        if(c_kkw.gt.0) write(c_out,*) 'SPECIAL LOOP', ic,              
!    &       c_tchem(ic)                                                
                                                                        
         if(c_tchem(ic).eq.'     H2O'.or.c_tchem(ic).eq.'    H2O '      &
     &                                                 )    then        
            c_nh2o =ic 
            c_icat(ic) = 0 
         endif 
         if(c_tchem(ic).eq.'     CO2'.or.c_tchem(ic).eq.'    CO2 '      &
     &                                                 )      then      
            c_nco2 =ic 
            c_icat(ic) =20 
         endif 
         if(c_tchem(ic).eq.'      CO'.or.c_tchem(ic).eq.'     CO '      &
     &                                                 )      then      
            c_nco  =ic 
         endif 
         if(c_tchem(ic).eq.'     CH4'.or.c_tchem(ic).eq.'    CH4 '      &
     &                                                 )      then      
            c_nch4 =ic 
         endif 
         if(c_tchem(ic).eq.'     CO3'.or.c_tchem(ic).eq.'    CO3 '      &
     &                                                 )      then      
            c_nco3 =ic 
            c_icat(ic) = 8 
         endif 
         if(c_tchem(ic).eq.'      H+'.or.c_tchem(ic).eq.'     H+ '      &
     &                                                 )      then      
            c_nhplus  =ic 
            c_icat(ic) =17 
         endif 
         if(c_tchem(ic).eq.'     OH-'.or.c_tchem(ic).eq.'    OH- '      &
     &                                                 )      then      
            c_nohmin =ic 
            c_icat(ic) =18 
         endif 
         if(c_tchem(ic).eq.'      OH'.or.c_tchem(ic).eq.'     OH '      &
     &                                                 )      then      
            c_noh  =ic 
            c_icat(ic) = 9 
                                                                        
! TEMPORARY DEBUG WRITE                                                 
!            if(c_kkw.gt.0) write(c_out,*) 'c_noh =', c_noh             
                                                                        
         endif 
         if(c_tchem(ic).eq.'     HO2'.or.c_tchem(ic).eq.'    HO2 '      &
     &                                                 )      then      
            c_nho2 =ic 
            c_icat(ic) =10 
         endif 
         if(c_tchem(ic).eq.'    H2O2'.or.c_tchem(ic).eq.'   H2O2 '      &
     &                                                 )      then      
            c_nh2o2=ic 
            c_icat(ic) =19 
         endif 
         if(c_tchem(ic).eq.'      O3'.or.c_tchem(ic).eq.'     O3 '      &
     &                                                 )      then      
            c_no3  =ic 
            c_icat(ic) =11 
         endif 
         if(c_tchem(ic).eq.'     NO2'.or.c_tchem(ic).eq.'    NO2 '      &
     &                                                 )      then      
            c_nno2 =ic 
            c_icat(ic) =12 
         endif 
         if(c_tchem(ic).eq.'      NO'.or.c_tchem(ic).eq.'     NO '      &
     &                                                 )      then      
            c_nno  =ic 
            c_icat(ic) =13 
         endif 
         if(c_tchem(ic).eq.'     NO3'.or.c_tchem(ic).eq.'    NO3 '      &
     &                                                 )      then      
            c_nno3 =ic 
            c_icat(ic) =14 
         endif 
         if(c_tchem(ic).eq.'    N2O5'.or.c_tchem(ic).eq.'   N2O5 '      &
     &                                                 )      then      
            c_nn2o5=ic 
            c_icat(ic) =15 
         endif 
         if(c_tchem(ic).eq.'    HNO3'.or.c_tchem(ic).eq.'   HNO3 '      &
     &                                                 )      then      
            c_nhno3=ic 
            c_icat(ic) =16 
         endif 
                                                                        
                        !do ic=1, c_nchem2                              
      enddo 
! END  LOOP TO IDENTIFY SPECIAL SPECIES                                 
                                                                        
! WARN IF THESE SPECIES ARE MISSING                                     
                                                                        
      if(c_nh2o .eq.0) write(c_out, 1381) 
 1381 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING',         &
     &         ' OR MISNAMED:  H2O')                                    
      if(c_nco2 .eq.0) write(c_out, 1382) 
 1382 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING',         &
     &         ' OR MISNAMED:  CO2')                                    
      if(c_nhplus  .eq.0) write(c_out, 1383) 
 1383 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING',         &
     &         ' OR MISNAMED:   H+')                                    
      if(c_nohmin .eq.0) write(c_out, 1384) 
 1384 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING',         &
     &         ' OR MISNAMED:  OH-')                                    
      if(c_noh  .eq.0) write(c_out, 1385) 
 1385 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING',         &
     &         ' OR MISNAMED:   OH')                                    
      if(c_nho2 .eq.0) write(c_out, 1386) 
 1386 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING',         &
     &         ' OR MISNAMED:  HO2')                                    
      if(c_nh2o2.eq.0) write(c_out, 1387) 
 1387 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING',         &
     &         ' OR MISNAMED: H2O2')                                    
      if(c_no3  .eq.0) write(c_out, 1388) 
 1388 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING',         &
     &         ' OR MISNAMED:   O3')                                    
      if(c_nno2 .eq.0) write(c_out, 1389) 
 1389 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING',         &
     &         ' OR MISNAMED:  NO2')                                    
      if(c_nno  .eq.0) write(c_out, 1390) 
 1390 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING',         &
     &         ' OR MISNAMED:   NO')                                    
      if(c_nno3 .eq.0) write(c_out, 1391) 
 1391 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING',         &
     &         ' OR MISNAMED:  NO3')                                    
      if(c_nn2o5.eq.0) write(c_out, 1392) 
 1392 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING',         &
     &         ' OR MISNAMED: N2O5')                                    
      if(c_nhno3.eq.0) write(c_out, 1393) 
 1393 format(/,'WARNING. THE FOLLOWING SPECIES MAY BE MISSING',         &
     &         ' OR MISNAMED: HNO3')                                    
                                                                        
                                                                        
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
                                                                        
! PRELIMINARY ZERO                                                      
         do 610 ic=1,c_nchem2 
          c_nequil(ic) = 0 
          c_npequil(ic) = ic 
          c_ion(ic) = 0 
          do 620 i=1,6 
           c_ncequil(ic,i) = 0 
  620      continue 
  610     continue 
                                                                        
! SET ION INDEX FOR H+, OH-                                             
! H+                                                                    
       if(c_aqueous(1,2).gt.0) then 
         c_ion(c_aqueous(1,2))= 1 
       endif 
! OH-                                                                   
       if(c_aqueous(1,3).gt.0) then 
         c_ion(c_aqueous(1,3))=-1 
       endif 
                                                                        
! LOOP FOR HENRY'S LAW COEFFICIENTS                                     
                                                                        
       do 650 nrh=1,c_nreach 
        if(c_henry(nrh,1).eq.0) go to 651 
                                                                        
! FOR EACH HENRY'S LAW COEFFICIENT, ASSOCIATE AQUEOUS WITH GAS SPECIES. 
        ic = c_henry(nrh,1) 
        ic1 = c_henry(nrh,2) 
        c_npequil(ic1) = ic 
        c_icat(ic1)=c_icat(ic) 
        c_lsts(ic1)=c_lsts(ic) 
        c_nequil(ic) = c_nequil(ic) + 1 
        nne=c_nequil(ic) 
        c_ncequil(ic,nne) = ic1 
        c_nrequil(ic,nne) = nrh 
                                                                        
! CHANGE PRODUCT ARRAY TO GAS-MASTER                                    
!                     (add to avoid replacing already-set array)        
       do nr=1,c_nreac 
         c_prodarr(nr,ic) = c_prodarr(nr,ic) + c_prodarr(nr,ic1) 
         c_prodarr(nr,ic1) = 0. 
                                                                        
! TEST WRITE PRODARR                                                    
!        if(c_kkw.ne.0.and.c_prodarr(nr,ic).ne.0)                       
!    *    write(c_out,656) c_tchem(ic), c_tchem(ic1),                   
!    *      nr, ic,ic1, c_prodarr(nr,ic)                                
! 656      format('TEST PRODARR ADJUSTMENT: NR IC IC1 PRODARR= ',       
!    *       /,2(a8,2x),3i5, f10.3)                                     
       enddo 
                                                                        
! GO THROUGH AQUEOUS EQUILIBRIA.  ADD EQUILIBRIUM SPECIES               
! TO GAS-MASTER SPECIES LIST.                                           
! NOTE EQUILIBRIA MUST COME IN SEQUENCE:                                
!             GAS->AQUEOUS->FIRST ION->DOUBLE ION.                      
! ALSO SET ION INDEX.                                                   
                                                                        
        do 660 nrq=1,c_nreacq 
        if(c_aqueous(nrq,1).eq.0.and.c_aqueous(nrq,2).eq.0) go to 661 
                                                                        
!  H+ SPECIAL here: for H+, OH-  skip equilibrium index assignment      
        if(c_aqueous(nrq,1).eq.0) go to 660 
                                                                        
                                                                        
        if(c_npequil(c_aqueous(nrq,1)).eq.ic) then 
          ic1 = c_aqueous (nrq,2) 
          c_npequil(ic1) = ic 
          c_icat(ic1)=c_icat(ic) 
          c_lsts(ic1)=c_lsts(ic) 
          c_nequil(ic) = c_nequil(ic) + 1 
          nne=c_nequil(ic) 
          c_ncequil(ic,nne) = ic1 
          c_nrequil(ic,nne) = nrq 
                                                                        
          c_ion(ic1) = c_ion(c_aqueous(nrq,1)) 
          if(c_aqueous(nrq,3).gt.0) then 
            c_ion(ic1) = c_ion(ic1) - c_ion(c_aqueous(nrq,3)) 
            if(c_kkw.gt.0)  write(c_out,659) (c_treacq(j,nrq),j=1,3),   &
     &       c_ion(ic1), c_ion(c_aqueous(nrq,3))                        
  659        format(a8,'=',a8,'+',a8,2x, 'ION2=',i3, '  ION3=',i3) 
          endif 
                                                                        
! CHANGE PRODUCT ARRAY TO GAS-MASTER.                                   
!                (add to avoid replacing already-set array)             
          do nr=1,c_nreac 
            c_prodarr(nr,ic) = c_prodarr(nr,ic) + c_prodarr(nr,ic1) 
            c_prodarr(nr,ic1) = 0. 
! TEST WRITE PRODARR                                                    
!        if(c_kkw.ne.0.and.c_prodarr(nr,ic).ne.0)                       
!    *    write(c_out,656) c_tchem(ic), c_tchem(ic1),                   
!    *      nr, ic,ic1, c_prodarr(nr,ic)                                
                                                                        
          enddo 
                                                                        
        endif 
                                                                        
  660   continue 
  661   continue 
                                                                        
  650  continue 
  651  continue 
                                                                        
!      if(c_kkw.gt.0) write(c_out,669) (c_npequil(i),i=1,c_nchem2)      
! 669    format(/'TEST NPEQUIL:',/,(20i5))                              
                                                                        
! SUMMARY WRITE: HENRY'S LAW AND AQUEOUS EQUILIBRIUM SPECIES (ncf)      
       write(c_out,1321) 
 1321   format(/,'HENRYs LAW AND LINKED AQUEOUS EQUILIBRIUM SPECIES',/) 
                                                                        
        write(c_out,1322) c_tchem(c_aqueous(1,2)),                      &
     &                c_tchem(c_aqueous(1,3))                           
        do  nrh=1,c_nreach 
          ic = c_henry(nrh,1) 
          if(c_nequil(ic).gt.0)  then 
            write(c_out, 1322) c_tchem(ic),                             &
     &        ( c_tchem(c_ncequil(ic,i))                                &
     &                                      ,i=1,c_nequil(ic) )         
! 1322            format(5(4x,a4,        2x))                           
 1322           format(5(   a8,        2x)) 
                                   !if(c_nequil(ic).gt.0)  then         
          endif 
                                !do  nrh=1,c_nreach                     
        enddo 
                                                                        
! -----------------------------------------                             
!  END ASSIGN-AQUEOUS                                                   
! -----------------------------------------                             
!                                                                       
!                                                                       
! -----------------------------------------                             
!  CONVERT SPECIAL EQUILIBRIUM SPECIES TO BACK-FORTH REACTIONS          
! -----------------------------------------                             
!  Special equilibrium forward reaction 1e8 sec-1,                      
!   backward reaction based on equilibrium coefficient.                 
!                                                                       
!  (Special equilibrium reactions NOT entered into cascade pair.)       
                                                                        
! LOOP FOR SPECIAL AQUEOUS EQUILIBRIUM SPECIES                          
                                                                        
       nr = c_nreac 
                                                                        
       do     nrqq=1,c_nreaqq 
        if(c_aqspec(nrqq,1).ne.0.and.c_aqspec(nrqq,2).ne.0) then 
          icc1 = c_aqspec(nrqq,1) 
          icc2 = c_aqspec(nrqq,2) 
          icc3 = c_aqspec(nrqq,3) 
                                                                        
          ic1 = c_npequil(icc1) 
          ic2 = c_npequil(icc2) 
          ic3 = c_npequil(icc3) 
                                                                        
! ENTER INTO CASCADE PAIR - CUT.                                        
! Note:  1st SPECIAL EQUIL species is secondary paired species          
!        2nd species is primary paired species  (A <=> B + C)           
!        3rd species is totally independent                             
! If this is included, also change 'ncasp =' above.                     
                                                                        
!         caspair(nrqq,  1 ) = icc2                                     
!         caspair(nrqq,  2 ) = icc1                                     
                                                                        
! ENTER FORWARD REACTION                                                
          nr = nr+1 
          if(nr.eq.1711.or.nr.eq.1911.or.nr.eq.3211)                    &
     &      write(c_out,906) nr                                         
          c_treac(1,nr) = c_tchem(icc1) 
          c_treac(3,nr) = c_tchem(icc2) 
          c_reactant(nr,1) = icc1 
          c_product(nr,1) = icc2 
          c_stoich(nr,1) = 1. 
          c_nnpro(nr) = 1 
          c_prodarr(nr,ic2 ) = 1. 
          if(icc3.gt.0) then 
           c_treac(4,nr) = c_tchem(icc3) 
           c_product(nr,2) = icc3 
           c_stoich(nr,2) = 1. 
           c_nnpro(nr) = 2 
           c_prodarr(nr,ic3 ) = 1. 
          endif 
! FORWARD RATE CONSTANT:  1e8 s-1.                                      
! OPTION - maybe  this causes cncn-> 0.  alt 1e2 s-1                    
!         c_rk(1,nr) = 1.0D+08                                          
          c_rk(1,nr) = 1.0D+02 
                                                                        
! ENTER BACKWARD REACTION                                               
          nr = nr+1 
          if(nr.eq.711.or.nr.eq.911.or.nr.eq.1211)                      &
     &      write(c_out,906) nr                                         
          c_treac(1,nr) = c_tchem(icc2) 
          c_treac(3,nr) = c_tchem(icc1) 
          c_reactant(nr,1) = icc2 
          c_product(nr,1) = icc1 
          c_stoich(nr,1) = 1. 
          c_nnpro(nr) = 1 
          c_prodarr(nr,ic1 ) = 1. 
          if(icc3.gt.0) then 
           c_treac(2,nr) = c_tchem(icc3) 
           c_reactant(nr,2) = icc3 
          endif 
! BACKWARD RATE CONSTANT                                                
!  Equil  constant in  moles/liter A exp(-B * (1/temp - 1/298) )        
!  A<->B+C  KA = BC  K in moles/liter                                   
!  k1a = k2BC,  k1/k2=K,  k2=k1/K,  k2 as A2 exp(-B2/temp) also moles/li
!  k1=1e8 s-1.  A2=k1/(Aexp(+B/298),  B2=-B                             
          c_rk(1,nr) = 1.0D+08/                                         &
     &                  ( c_rkqq(1,nrqq)*exp(c_rkqq(2,nrqq)/298.) )     
          c_rk(2,nr) = 0. - c_rkqq(2, nrqq ) 
                                                                        
        endif 
       enddo 
! END LOOP FOR SPECIAL EQUILIBRIUM SPECIES                              
                                                                        
       c_nreac = nr 
! -----------------------------------------                             
!  END CONVERT SPECIAL EQUILIBRIUM SPECIES TO BACK-FORTH REACTIONS      
! -----------------------------------------                             
!                                                                       
! -----------------------------------------                             
!  ADD NO3-N2O5 (hard-wired) INTO CASCADE PAIR LIST                     
! -----------------------------------------                             
       ncasp = ncasp + 1 
       caspair(ncasp,  1 ) = c_nno3 
       caspair(ncasp,  2 ) = c_nn2o5 
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
      do i=1,ncasp 
        icc1 = caspair(i,1) 
        icc2 = caspair(i,2) 
        ic1 = 0 
        ic2 = 0 
        if(icc1.gt.0) then 
          ic1 = c_npequil(icc1) 
        endif 
        if(icc2.gt.0) then 
          ic2 = c_npequil(icc2) 
        endif 
        ic3 = 0 
        if(ic2.gt.0) then 
         ic3=c_nppair(ic2,1) 
        endif 
                                                                        
! ERROR CHECK:  IF 2ND SPECIES IS ALREADY PAIRED, MAJOR ERROR           
        if(ic3.ne.ic2.or.icc1.eq.0.or.icc2.eq.0.or.                     &
     &                       ic1.eq.0.or.ic2.eq.0)  then                
          write(c_out,667) i, icc1, icc2, ic1, ic2, c_nppair(ic2,1) 
          write( 6  ,667) i, icc1, icc2, ic1, ic2, c_nppair(ic2,1) 
  667     format(/,' MAJOR ERROR IN CASCADE PAIR SETUP:',/,             &
     &     '   EITHER ZERO OR SECOND CASCADE PAIR ALREADY PAIRED.',     &
     &    /,' PAIR NUMBER, icc1, icc2, ic1, ic2, c_nppair(ic2,1)=', 6i5)
          if(icc1.gt.0) write (c_out,309) c_tchem(icc1) 
          if(icc2.gt.0) write (c_out,309) c_tchem(icc2) 
          if(ic1.gt.0) write (c_out,309) c_tchem(ic1) 
          if(ic2.gt.0) write (c_out,309) c_tchem(ic2) 
          if(icc1.gt.0) write (c_out,309) c_tchem(icc1) 
          if(icc1.gt.0) write ( 6  ,309) c_tchem(icc1) 
          if(icc2.gt.0) write ( 6  ,309) c_tchem(icc2) 
          if(ic1.gt.0) write ( 6  ,309) c_tchem(ic1) 
          if(ic2.gt.0) write ( 6  ,309) c_tchem(ic2) 
          if(icc1.gt.0) write ( 6  ,309) c_tchem(icc1) 
        else 
                                                                        
! Assign all species and subspecies to the primary species              
!               of the first of the pair.                               
!  Note:  npair(ic2,3) = number of subspecies;                          
!         npair(ic2,4) = 1st subsp, etc.                                
!  Also assign direct pair species - just for this secondary species.   
                                                                        
          c_nppair(ic2,1) = ic1 
                                                                        
          ic3=c_nppair(ic1,2) 
          do ii = 1, (c_nppair(ic2,3)+1) 
            ic  = ic2 
            if(ii.gt.1) ic  = c_nppair(ic2,(ii+2)) 
            c_nppair(ic ,2) = ic3 
            c_nppair(ic3,3) = c_nppair(ic3,3) + 1 
            np = c_nppair(ic3,3) + 3 
            c_nppair(ic3,np) = ic 
                                                                        
! Assign CATEGORY and STEADY STATE INDEX to secondary species           
!  also assign to all its subspecies, including aqueous equil. species. 
                                                                        
!           do neq=1, (c_nequil(ic )+1)                                 
!             icc = ic                                                  
!             if(neq.gt.1) icc = c_ncequil(ic, (neq-1) )                
!             c_nppair(ic ,1) = ic3                                     
!             c_icat(icc) = c_icat(ic3)                                 
!             c_lsts(icc) = c_lsts(ic3)                                 
!           enddo                                                       
                                                                        
          enddo 
! End Assign loop                                                       
                                                                        
        endif 
! End error check if                                                    
                                                                        
      enddo 
! END LOOP THROUGH CASCADE PAIRS                                        
!                                                                       
! -----------------------------------------                             
! END GENERATE POINTERS FOR CASCADE PAIRS                               
! -----------------------------------------                             
                                                                        
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
      do ic=1,c_nchem2 
        lcastest(ic) = .FALSE. 
        if(ic.eq.c_aqueous(1,2).or.ic.eq.c_aqueous(1,3))                &
     &        lcastest(ic) = .TRUE.                                     
        if(c_tchem(c_npequil(ic )).eq.'     CO2') lcastest(ic) = .TRUE. 
        if(c_tchem(c_npequil(ic )).eq.'     H2O') lcastest(ic) = .TRUE. 
                       !do ic=1,c_nchem2                                
      enddo 
                                                                        
!   Set LCASTEST = true for LUMPED SPECIES - these are never included.  
      do     i=1,nlump 
         icc=c_lump(i,1) 
         lcastest(icc)=.TRUE. 
      enddo 
                                                                        
! SET COUNTERS (NSOL, NSOLV) AND MAXIMA (NCDIM, NSDIM)                  
      ncdim = 1 
      nsdim = 1 
                                                                        
! TEST WRITE CASCADE                                                    
      if(c_kkw.gt.0) write(c_out,347) 
  347 format(/,'CASCADE TEST: LIST OF SPECIES IN ORDER OF SOLUTION') 
                                                                        
!  RUN  CASCADE.  IDENTIFY INCLUDED SPECIES.                            
!   COUNT SPECIES IN CASCADE GROUP.   ALSO ASSIGN NMULTI                
!                                                                       
      do     i=1,c_nchem2 
        if(c_cascade(i,1).eq.0) go to 332 
        nsolv = 0 
        nsol = 0 
        do     ii=1,c_nchem2 
          ics =c_cascade(i,ii) 
          if(ics.eq.0) go to 331 
          ics =c_nppair(c_npequil(ics ) , 2) 
          if(ii.eq.1)  ic1=ics 
          nsol = nsol + 1 
          nsolv = nsolv + 1 + c_nppair(ics,3) 
                                                                        
          do n=1,(c_nppair(ics,3)+1) 
            icc = ics 
            if(n.gt.1) then 
              icc=c_nppair(ics,n+2) 
            endif 
                                                                        
! TEST WRITE CASCADE                                                    
               if(c_kkw.ge.1.and.icc.gt.0) write(c_out,348)n,           &
     &                                              c_tchem(icc)        
  348        format(i4,4x,a8) 
                                                                        
            do  neq=1,(c_nequil(icc)+1) 
              ic=icc 
              if(neq.gt.1) ic = c_ncequil(icc,(neq-1)) 
              lcastest(ic) = .TRUE. 
              c_npmulti(ic,1) = ic1 
                                   !do  neq=1,(c_nequil(icc)+1)         
            enddo 
                                 !do n=1,(c_nppair(ic,3)+1)             
          enddo 
                             !do     ii=1,c_nchem2                      
        enddo 
  331   continue 
        if(ncdim.lt.nsolv) ncdim=nsolv 
        if(nsdim.lt.nsol)  nsdim=nsol 
                           !do     i=1,c_nchem2                         
      enddo 
  332 continue 
                                                                        
!  INCLUDE FAMILY SPECIES, SAME AS FOR CASCADE.                         
!                                                                       
      do     i=1,7 
       nsolv = 0 
       nsol = 0 
       do ii=1,3 
                                        ! ADDED 2009                    
        ics = 0 
        if(i.eq.1.and.ii.eq.1) ics=c_nhno3 
        if(i.eq.2.and.ii.eq.1) ics=c_nno3 
!       if(i.eq.3.and.ii.eq.1) ics=c_nn2o5                              
        if(i.eq.4.and.ii.eq.1) ics=c_no3 
        if(i.eq.5.and.ii.eq.1) ics=c_nno2 
        if(i.eq.6.and.ii.eq.1) ics=c_nno 
        if(i.eq.7.and.ii.eq.1) ics=c_noh 
        if(i.eq.7.and.ii.eq.2) ics=c_nho2 
        if(i.eq.7.and.ii.eq.3) ics=c_nh2o2 
                                                                        
         if(ics.gt.0) then 
          ics =c_nppair(c_npequil(ics ) , 2) 
          nsol = nsol + 1 
          nsolv = nsolv + 1 + c_nppair(ics,3) 
                                                                        
          do n=1,(c_nppair(ics,3)+1) 
            icc = ics 
            if(n.gt.1) then 
              icc=c_nppair(ics,n+2) 
            endif 
                                                                        
! TEST WRITE CASCADE - SPECIAL ASSIGNED SPECIES                         
            if(c_kkw.ge.1) write(c_out, *   )                           &
     &       'AUTOMATIC CASCADE SPECIES:', n,c_tchem(icc)               
                                                                        
                                                                        
            do  neq=1,(c_nequil(icc)+1) 
              ic=icc 
              if(neq.gt.1) ic = c_ncequil(icc,(neq-1)) 
              lcastest(ic) = .TRUE. 
                                   !do  neq=1,(c_nequil(icc)+1)         
            enddo 
                                 !do n=1,(c_nppair(ic,3)+1)             
          enddo 
                                !if(ics.gt.0) then                      
         endif 
                        ! do ii=1,3                                     
        enddo 
                                                                        
        if(ncdim.lt.nsolv) ncdim=nsolv 
        if(nsdim.lt.nsol)  nsdim=nsol 
                           !do     i=1,9                                
      enddo 
!                                                                       
! TEST WRITE                                                            
      if(c_kkw.gt.0) write(c_out,336) ncdim, nsdim 
  336 format(/,' MAXIMUM CASCADE SPECIES, PAIR GROUPS = ', 2i5) 
                                                                        
       if(c_kkw.eq.1.or.c_kkw.eq.5) write(c_out,342) 
  342  format(/'TEST IC CHEM ICAT ISTS  NEQUIL PRIMARY  NPAIR PRIMARY') 
                                                                        
! CATEGORY TEST LOOP AMONG SPECIES.                                     
!   IDENTIFY SPECIES MISSING FROM CASCADE AND MISSING CATEGORIES.       
                                                                        
       do ic=1,c_nchem2 
!  Control for diagnostic write                                         
         if(c_kkw.eq.1.or.c_kkw.eq.5) then 
           write(c_out,343) ic, c_tchem(ic), c_icat(ic), c_lsts(ic),    &
     &       c_nequil(ic),c_npequil(ic),  c_nppair(ic,3)                &
     &       , c_nppair(ic,2)  ,c_nppair(ic,1)                          
  343       format(i4,2x,a8,2x,i5,l4,2x,8i5) 
           if(c_npequil(ic).gt.0) then 
             write(c_out,309) c_tchem(c_npequil(ic)) 
           endif 
           if(c_nppair(ic,2).gt.0) then 
             write(c_out,309) c_tchem(c_nppair(ic,2)) 
           endif 
           if(c_nppair(ic,1).gt.0) then 
             write(c_out,309) c_tchem(c_nppair(ic,1)) 
           endif 
           if(c_nequil(ic).gt.0) then 
             write(c_out,344) (c_ncequil(ic,i),i=1,c_nequil(ic)) 
  344        format(' c_ncequil(aqueous subspecies) = ',10i5) 
             do i=1, c_nequil(ic) 
               if(c_ncequil(ic,i).gt.0) then 
                 write(c_out,309) c_tchem(c_ncequil(ic,i)) 
               endif 
             enddo 
           endif 
           if(c_nppair(ic,3).gt.0) then 
             write(c_out,345) ( c_nppair(ic,i),i=4,(3+c_nppair(ic,3)) ) 
  345        format(' ncpair (chem. paired subspecies)  = ',10i5) 
             do i=4, (3+c_nppair(ic,3)) 
               if(c_nppair(ic,i).gt.0) then 
                 write(c_out,309) c_tchem(c_nppair(ic,i)) 
               endif 
             enddo 
           endif 
         endif 
! End control for diagnostic write.                                     
                                                                        
         if(c_icat(ic).le.0) then 
          write(c_out,341) ic,c_tchem(ic), c_icat(ic) 
  341     format(/,' WARNING:  MISSING SPECIES CATEGORY',/,             &
     &   'SPECIES =',i4,2x,a8,'  CATEGORY=',i4)                         
!         c_icat(ic) = 23                                               
         endif 
                                                                        
         if(.not.lcastest(ic))                                          &
     &    write(c_out,352) ic,c_tchem(ic), c_npequil(ic)                &
     &  ,c_tchem(c_npequil(ic))                                         
  352     format(/,' WARNING:  SPECIES OMITTED FROM CASCADE',/,         &
     &    'SPECIES =',i4,4x,a8,'  NPEQUIL =',i4,2x, a8)                 
                                                                        
                 !do ic=1,c_nchem2                                      
       enddo 
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
                                                                        
       do 350 i=1,c_nchem2 
        lcastest(i) = .FALSE. 
  350  continue 
                                                                        
! RUN CASCADE.  SET INDEX=1 AS EACH SPECIES IS SOLVED FOR.              
! ( DO WHILE NCASCADE(1).NE.0)                                          
                                                                        
! LOOP:  RUN CASCADE                                                    
      do     i=1,c_cdim 
       if(c_cascade(i,1).eq.0) go to 401 
                                                                        
! LOOP:  TEST REACTIONS CALLED BY THE CASCADE SPECIES                   
       do     ii=1,c_cdim 
        icc=c_cascade(i,ii) 
        if(icc.gt.0) then 
         icc=c_npequil(      (icc)) 
         ic=icc 
!                                                                       
! NOTE PROBLEM HERE (write 419): TEST OF REACTANTS AND PRODUCTS         
!   BELONGS IN QUADINIT                                                 
                                                                        
! FIRST TEST REACTANTS AND PRODUCTS FROM CASCADE REACTIONS.             
!  EXCEPT FOR 'POST' FLAG, PRODUCTS SHOULD BE ALL DOWN-CASCADE.         
!  CASCADE MODIFIED TO REFLECT ALL CHEMICAL PAIRS                       
!                                ->POINT TO PRIMARY SPEC.               
          if(c_cascade(i,2).ne.-2) then 
            nnr=c_nnrchem(ic) 
            if(nnr.gt.0) then 
              do 420 n=1,nnr 
               nr= c_nrchem(ic,n) 
               if(nr.gt.0) then 
                 icr1=c_reactant(nr,1) 
                 if(icr1.gt.0) then 
                   icr1=c_nppair( c_npequil(icr1), 2 ) 
                 endif 
                 icr2=c_reactant(nr,2) 
                 if(icr2.gt.0) then 
                   icr2=c_nppair( c_npequil(icr2), 2 ) 
                 endif 
                 if(c_nnpro(nr).gt.0) then 
                   do 425 nn=1,c_nnpro(nr) 
                     icp=c_product(nr,nn) 
                     if(icp.gt.0) then 
                       icp=c_nppair( c_npequil(icp), 2 ) 
                     endif 
                     if(icp.eq.icr1) then 
                       icr1=0 
                       icp=0 
                     endif 
                     if(icp.eq.icr2) then 
                       icr2=0 
                       icp=0 
                     endif 
! PRODUCT TEST                                                          
                      if(icp.gt.0) then 
                         if(lcastest(icp)) then 
                           write(c_out,419)                             &
     &                        c_tchem(icp),nr,c_tchem(ic)               &
     &                            ,(c_treac(jj,nr),jj=1,5)              
                           write(   6,419)                              &
     &                        c_tchem(icp),nr,c_tchem(ic)               &
     &                            ,(c_treac(jj,nr),jj=1,5)              
                        endif 
                      endif 
  425              continue 
                 endif 
! REACTANT TEST                                                         
                 if(icr1.gt.0) then 
                   if(lcastest(icr1)) then 
                      write(c_out,419)                                  &
     &                  c_tchem(icr1),nr,c_tchem(ic)                    &
     &                    ,(c_treac(jj,nr),jj=1,5)                      
                      write(   6,419)                                   &
     &                  c_tchem(icr1),nr,c_tchem(ic)                    &
     &                    ,(c_treac(jj,nr),jj=1,5)                      
                   endif 
                 endif 
                 if(icr2.gt.0) then 
                   if(lcastest(icr2))then 
                     write(c_out,419)                                   &
     &                  c_tchem(icr2),nr,c_tchem(ic)                    &
     &                    ,(c_treac(jj,nr),jj=1,5)                      
                     write(   6,419)                                    &
     &                  c_tchem(icr2),nr,c_tchem(ic)                    &
     &                    ,(c_treac(jj,nr),jj=1,5)                      
                   endif 
                 endif 
                                                                        
               endif 
  420         continue 
            endif 
          endif 
                                                                        
! LAST, ADD CASCADE REACTIONS TO REACTION LIST, UNLESS 'pre' FLAG       
          if(c_cascade(i,2).ne.-1) then 
             if(lcastest(ic)) then 
               write(c_out,418) c_tchem(ic),i, (c_cascade(i,jj),jj=1,3) 
               write(   6,418) c_tchem(ic),i, (c_cascade(i,jj),jj=1,3) 
             endif 
             lcastest(ic) = .TRUE. 
          endif 
                                                                        
        endif 
       enddo 
! END LOOP:  TEST REACTIONS CALLED BY THE CASCADE SPECIES               
                                                                        
      enddo 
! END LOOP:  RUN CASCADE                                                
                                                                        
  401 continue 
                                                                        
  419 format(/,'WARNING:  CASCADE OUT OF SEQUENCE FOR SPECIES NAME=  ', &
     &  a8,/,'PRODUCED AT REACTION #',i5,'  CALLED FOR SPECIES= ',      &
     &  a8,/,' REACTION: ',a8,'+',a8,'=>',a8,'+',a8,'+',a8)             
                                                                        
  418 format(/,'WARNING:  SPECIES CALLED TWICE IN CASCADE. SPECIES=',   &
     &  a8,/,' SECOND CALL AT CASCADE#',i3,' WITH INDICES=',3i3)        
                                                                        
                                                                        
! -------------------------------------------------------               
! END CASCADE TEST                                                      
! -------------------------------------------------------               
                                                                        
! OPTION:  call cheminit here to set initial chemistry parameters       
                                                                        
! ----------------------------------------------------------------      
!  END READ.                                                            
! ----------------------------------------------------------------      
                                                                        
! END CHEMREAD                                                          
 2000 return 
      END                                           
! ----------------------------------------------------------------      
                                                                        
                                                                        
                                                                        
       function namechem(titl) 
!      integer function namechem(titl)                                  
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
!                                                                       
      IMPLICIT NONE 
                                                                        
      include 'chemmech.EXT' 
!     INCLUDE 'chemvars.EXT'                                            
      include 'chemlocal.EXT' 
                                                                        
      character*8 titl 
      integer namechem 
                                                                        
       namechem = 0 
       if(titl.eq.'        ') return 
       do 100 i=1,c_nchem2 
        if(titl.eq.c_tchem(i)) then 
          namechem = i 
          go to 101 
        endif 
  100  continue 
  101  continue 
                                                                        
! END NAMECHEM                                                          
 2000   return 
      END                                           
! ----------------------------------------------------------------      
!                                                                       
                                                                        
      subroutine hvread 
                                                                        
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
      IMPLICIT NONE 
                                                                        
      include 'chemmech.EXT' 
!     INCLUDE 'chemvars.EXT'                                            
! c     include 'chemlocal.EXT'                                         
                                                                        
       call readhv(c_hvin,c_nhv,c_hvmat, c_hvmatb, c_jarray) 
                                                                        
! END HVREAD                                                            
 2000  return 
      END                                           
! --------------------------------------------------------------        
!                                                                       
!                                                                       
                                                                        
          subroutine cheminit 
!                   (quadinit)                                          
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
!                                                                       
!                                                                       
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
      IMPLICIT NONE 
                                                                        
      include 'chemmech.EXT' 
!     INCLUDE 'chemvars.EXT'                                            
      include 'chemlocal.EXT' 
!                                                                       
! LOCAL VARIABLES                                                       
! lloss        Local: Flag for identifying exchange loss reaction       
! lpro         Local: Flag for identifying exchange production reaction 
                                                                        
!  xoddhx      Counter for net change in odd hydrogen in reaction       
!                (intermediate for c_oddhx) - split into RO2; OH-HO2-CO3
!  xpronox     Counter for net production of odd nitrogen in reaction   
!                (intermediate for c_pronox)                            
                                                                        
                      ! Flag for identifying exchange loss reaction     
      logical lloss 
                      ! Flag for id exchange production reaction        
      logical  lpro 
                                                                        
                                   ! Counter for odd hydrogen RO2 only  
      double precision xoddhx 
                                   ! Counter for odd hydrogen w/o RO2   
      double precision xoddhx2 
                                   ! Counter for odd nitrogen           
      double precision xpronox 
                                   ! Counter for product reactions      
      integer nrp 
                                                                        
!                                                                       
! LOCAL DIMENSION                                                       
!         integer neqread(6)                                            
                                                                        
! GENERAL FORMATS                                                       
   11 format(///////,a1) 
! 12    format(i6)                                                      
! 13    format((1pe10.3))                                               
                                                                        
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
                                                                        
! THE IDENTIFICATION OF REACTIONS WITH SPECIES INSURES THAT             
! ALL REACTIONS WILL BE PROCESSED AT THE PROPER PLACE IN THE CASCADE.   
                                                                        
! (NEW) LOSS AND PRODUCTION REACTIONS ARE ID'D SEPARATELY.              
! -----------------------------------------                             
                                                                        
       kk=1 
                                                                        
                                                                        
!  ZERO THE IMPORTANT ARRAYS                                            
       do 160 ic=1,c_nchem2 
        c_nnrchem(ic) = 0 
        c_nnrchp(ic) = 0 
        do 157 i=1,25 
         c_nrchem(ic,i) = 0 
         c_nrchmp(ic,i) = 0 
  157   continue 
  160  continue 
                                                                        
       do 165 nr=1,c_nreac 
         c_stoiloss(nr) = 0. 
  165  continue 
                                                                        
! LOOP FOR EACH INDIVIDUAL REACTION                                     
       do 1000 nr=1,c_nreac 
                                                                        
! ZERO INDEX FOR IC ASSOCIATED WITH REACTION.                           
         ic = 0 
                                                                        
! ESTABLISH REACTANT CATEGORIES.                                        
!     NOTE: REACTION IS ASSIGNED ONLY FOR ICAT>0.  CHECK ERROR.         
!                                                                       
       icat1 = 0 
       icat2 = 0 
       if(c_reactant(nr,1).gt.0) then 
         icat1 = c_icat(c_reactant(nr,1)) 
       endif 
       if(c_reactant(nr,2).gt.0) then 
         icat2 = c_icat(c_reactant(nr,2)) 
       endif 
                                                                        
! FIRST ELIMINATE REDUNDANT REACTANTS:                                  
!   REACTANTS THAT ARE ALSO PRODUCTS,  AND ALSO H2O                     
!                                                                       
!  (This is for reactions RCO3+RO2=> products.                          
!    Entered as two reactions:                                          
!       RCO3+RO2=>RCO3+productB                                         
!       RCO3+RO2=>RO2 +productA                                         
!    Assign reaction with RO2 products to RO2, not RCO3.)               
                                                                        
      if(c_nnpro(nr).gt.0) then 
       do 305 n=1,c_nnpro(nr) 
        if(c_reactant(nr,1).eq.c_product(nr,n).and.c_stoich(nr,n).ge.1) &
     &    icat1 = 0 - icat1                                             
        if(c_reactant(nr,2).eq.c_product(nr,n).and.c_stoich(nr,n).ge.1) &
     &    icat2 = 0 - icat2                                             
  305  continue 
      endif 
                                                                        
! FIRST ASSOCIATE REACTION WITH REACTANTS                               
       if(icat1.lt.9.and.icat1.gt.0                ) then 
         ic = c_reactant(nr,1) 
       else 
         if(icat2.lt.9.and.icat2.gt.0                ) then 
           ic = c_reactant(nr,2) 
         endif 
       endif 
                                                                        
! RESTORE CATEGORY OF REDUNDANT REACTIONS (RCO3+RO2)                    
       icat1 = abs(icat1) 
       icat2 = abs(icat2) 
                                                                        
! (PRIOR:  Look for NO3-N2O5 reactants before non-special products.     
!  Now, Look for nonspecial products first, then NO3-N2O5 reactants.)   
                                                                        
                                                                        
! IF NEITHER REACTANT IS A NON-SPECIAL SPECIES, THEN CHECK PRODUCTS     
!  (e.g. HNO3, H2O2)                                                    
       if(ic.le.0) then 
         if(c_nnpro(nr).gt.0) then 
          do 310 n=1,c_nnpro(nr) 
            icatp = c_icat(c_npequil(c_product(nr,n))) 
                                                                        
                                                                        
            if(icatp.lt.9.and.icatp.gt.0                ) then 
              ic = c_product(nr,n) 
              go to 311 
            endif 
  310     continue 
  311     continue 
         endif 
       endif 
                                                                        
! IF NEITHER REACTANT IS A NON-SPECIAL SPECIES, LOOK FOR NO3-N2O5-HNO3. 
                                                                        
                                                                        
       if(ic.le.0) then 
         if(icat1.eq.14.or.icat1.eq.15.or.icat1.eq.16) then 
           ic = c_reactant(nr,1) 
         else 
           if(icat2.eq.14.or.icat2.eq.15.or.icat2.eq.16) then 
             ic = c_reactant(nr,2) 
           endif 
         endif 
       endif 
                                                                        
                                                                        
! IF NOTHING, LOOK FOR NO3/N2O5/HNO3 IN PRODUCTS                        
       if(ic.le.0) then 
            if(c_nnpro(nr).gt.0) then 
             do 320 n=1,c_nnpro(nr) 
              icatp = c_icat(c_npequil(c_product(nr,n))) 
              if(icatp.eq.14.or.icatp.eq.15.or.icatp.eq.16) then 
                ic = c_product(nr,n) 
                go to 321 
              endif 
  320        continue 
  321        continue 
            endif 
       endif 
                                                                        
! IF STILL NOTHING, LOOK FOR O3 or NOx IN REACTANTS OR PRODUCTS         
!    O3 (11), NO2 (12) or NO (13)                                       
                                                                        
! 2005 CHANGE   O3, NO, NO2 ORDER DOESN'T MATTER,                       
!   BUT DECLARE THROUGH REACTANTS NOT PRODUCTS - else brpro too early   
!   (HO2+NO=>NO2+OH must  be assigned to NO2)                           
!   MUST BE  11 or 12 or 13 here                                        
                                                                        
       if(ic.le.0) then 
!        if(icat1.eq.11.or.icat1.eq.12) then                            
         if(icat1.eq.11.or.icat1.eq.12.or.icat1.eq.13) then 
           ic = c_reactant(nr,1) 
         else 
!          if(icat2.eq.11.or.icat2.eq.12) then                          
           if(icat2.eq.11.or.icat2.eq.12.or.icat2.eq.13) then 
             ic = c_reactant(nr,2) 
           else 
            if(c_nnpro(nr).gt.0) then 
             do 330 n=1,c_nnpro(nr) 
              icatp = c_icat(c_npequil(c_product(nr,n))) 
!             if(icatp.eq.11.or.icatp.eq.12) then                       
              if(icatp.eq.11.or.icatp.eq.12.or.icatp.eq.13) then 
                ic = c_product(nr,n) 
                go to 331 
              endif 
  330        continue 
  331        continue 
            endif 
           endif 
         endif 
       endif 
!                                                                       
! H2O2 OPTION:  MAKE H2O2 A SPECIAL SPECIES FOR ASSIGNING REACTIONS:    
! LOOK FOR H2O2 (19 ) IN REACTANTS OR PRODUCTS.                         
       if(ic.le.0) then 
         if(               icat1.eq.19) then 
           ic = c_reactant(nr,1) 
         else 
           if(               icat2.eq.19) then 
             ic = c_reactant(nr,2) 
           else 
            if(c_nnpro(nr).gt.0) then 
             do 345 n=1,c_nnpro(nr) 
              icatp = c_icat(c_npequil(c_product(nr,n))) 
              if(               icatp.eq.19) then 
                ic = c_product(nr,n) 
                go to 346 
              endif 
  345        continue 
  346        continue 
            endif 
           endif 
         endif 
       endif 
                                                                        
! IF STILL NOTHING, THE REMAINDER HAD BETTER HAVE OH OR HO2!            
!  (MAYBE TRY TO USE OH RATHER THAN HO2?)                               
       if(ic.le.0) then 
         if(icat1.eq.9.or.icat1.eq.10) then 
           ic = c_reactant(nr,1) 
         else 
           if(icat2.eq.9.or.icat2.eq.10) then 
             ic = c_reactant(nr,2) 
           else 
            if(c_nnpro(nr).gt.0) then 
             do 350 n=1,c_nnpro(nr) 
              icatp = c_icat(c_npequil(c_product(nr,n))) 
              if(icatp.eq.9.or.icatp.eq.10) then 
                ic = c_product(nr,n) 
                go to 351 
              endif 
  350        continue 
  351        continue 
            endif 
           endif 
         endif 
       endif 
                                                                        
!                                                                       
! IF STILL IC=0, ERROR AND EXIT!! (unless reaction=0)                   
      if(ic.le.0.and.                                                   &
     &     (c_reactant(nr,1).gt.0.or.c_reactant(nr,2).gt.0)) then       
      write(c_out,355) nr, c_reactant(nr,1), icat1,                     &
     &                  c_reactant(nr,2), icat2,                        &
     &  (c_product(nr,n),n=1,c_nnpro(nr) )                              
      write(6,355) nr, c_reactant(nr,1), icat1,                         &
     &                  c_reactant(nr,2), icat2,                        &
     &  (c_product(nr,n),n=1,c_nnpro(nr) )                              
  355 format(//,'MAJOR ERROR:  REACTION NOT IDENTIFIED WITH SPECIES',   &
     & //,' REACTION NUMBER =',i4,/,                                    &
     &    ' FIRST REACTANT =', i4,'  CATEGORY =',i3,/,                  &
     &    ' SECOND REACTANT=', i4,'  CATEGORY =',i3,/,                  &
     &    ' PRODUCTS =', (10i5)  )                                      
 3000  stop 3000 
      endif 
!                                                                       
                                                                        
! FINALLY!  ENTER REACTION INTO PROPER SPECIES LIST                     
! NEW - SEPARATE LIST FOR SPECIES-LOSS REACTIONS AND PRODUCT REACTIONS, 
!       AND CALC OF STOILOSS FOR LOSS REACTIONS.                        
! 2000-REACTIONS ASSIGNED TO MAIN SPECIES FOR SPECIAL AQUEOUS LINKS.    
                                                                        
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
                                                                        
      if(ic.gt.0) then 
        icc = c_npequil(ic) 
        if(ic.eq.c_reactant(nr,1).or.ic.eq.c_reactant(nr,2)) then 
          c_nnrchem(icc) = c_nnrchem(icc) + 1 
          c_nrchem(icc,c_nnrchem(icc)) = nr 
          c_nicreac(nr) = icc 
          c_stoiloss(nr) = 1. 
          if(c_reactant(nr,1).eq.c_reactant(nr,2)) c_stoiloss(nr) = 2. 
          if(c_nnpro(nr).gt.0) then 
            do 360 n=1,c_nnpro(nr) 
             if(ic.eq.c_product(nr,n))                                  &
     &          c_stoiloss(nr) = c_stoiloss(nr)-c_stoich(nr,n)          
  360       continue 
          endif 
        else 
          c_nnrchp(icc) = c_nnrchp(icc) + 1 
          c_nrchmp(icc,c_nnrchp(icc)) = nr 
        endif 
      endif 
                                                                        
! END ANALYSIS OF SPECIES-LINK WITH REACTIONS.                          
                                                                        
                                                                        
!  CALCULATE CHANGE IN ODD HYDROGEN (ODDHX)                             
!  AND NET PRODUCTION OF ODD NITROGEN (PRONOX, >0).                     
!  (CURRENT OPTION - RECORD NET PRONOX, THEN IN PROGRAM                 
!   SUM NOx PRODUCTION-ONLY.)                                           
!  (ODDHX IS NET CHANGE, +/-.  BUT PRONOX IS PRODUCTION ONLY.           
!   NET CHANGE IN NOX CAN BE OBTAINED FROM RP, RL FOR NO AND NO2,       
!   BUT NOX ANALYSIS REQUIRES SEPARATION OF NOX SOURCES AND SINKS.      
!   NOTE - PAN->NOX AND HNO3->NOX OVERPRODUCTION IS FIXED               
!   BY PANCORR AND NO3CORR - DISCOUNTS BACK-FORTH NOX PRODUCTION.)      
                                                                        
      xoddhx = 0. 
      xoddhx2 = 0. 
      xpronox = 0. 
                                                                        
      if(icat1.eq.9.or.icat1.eq.10.or.icat1.eq.8)                       &
     &                                          xoddhx2 = xoddhx2 - 1.  
      if(icat2.eq.9.or.icat2.eq.10.or.icat1.eq.8)                       &
     &                                          xoddhx2 = xoddhx2 - 1.  
      if(icat1.eq.3)                                                    &
     &                                            xoddhx = xoddhx - 1.  
      if(icat2.eq.3)                                                    &
     &                                            xoddhx = xoddhx - 1.  
      if(icat1.eq.12.or.icat1.eq.13) xpronox = xpronox - 1. 
      if(icat2.eq.12.or.icat2.eq.13) xpronox = xpronox - 1. 
                                                                        
! TEMPORARY TEST                                                        
!     if(icat1.eq.3.and.icat2.eq.3)                                     
!    *     write(c_out,12001) nr, (c_treac(ir,nr),ir=1,5),              
!    *         xoddhx                                                   
                                                                        
          if(c_nnpro(nr).gt.0) then 
      do 410 n=1,c_nnpro(nr) 
       icatp = c_icat(c_npequil(c_product(nr,n))) 
                                                                        
       if(icatp.eq.9.or.icatp.eq.10.or.icatp.eq.8)                      &
     &     xoddhx2 = xoddhx2 + c_stoich(nr,n)                           
       if(icatp.eq.3)                                                   &
     &     xoddhx = xoddhx + c_stoich(nr,n)                             
       if(icatp.eq.12.or.icatp.eq.13)                                   &
     &     xpronox = xpronox + c_stoich(nr,n)                           
                                                                        
                                                                        
! TEMPORARY TEST                                                        
!     if(icat1.eq.3.and.icat2.eq.3) write(c_out,*) icatp                
!     if(icat1.eq.3.and.icat2.eq.3)                                     
!    *     write(c_out,12001) nr, (c_treac(ir,nr),ir=1,5),              
!    *         xoddhx                                                   
! 12001      format(i5,2x,a8,'+',a8,'=>',a8,'+',a8,'+',a8,/,3(1pe10.3)) 
                                                                        
                                                                        
  410 continue 
           endif 
                                                                        
! GOT HERE                                                              
       c_oddhx(nr,1) = xoddhx + xoddhx2 
       c_oddhx(nr,2) = xoddhx2 
       c_pronox(nr) = xpronox 
                                                                        
! DEBUG -  TEST WRITE                                                   
        if(c_kkw.gt.0) then 
          write(c_out,1001) nr,(c_treac(j,nr),j=1,5),ic,c_tchem(ic)     &
     &    , c_oddhx(nr,1),c_oddhx(nr,2), c_pronox(nr)                   
 1001   format('REACTION #',i4,': ',a8,'+',a8,'=>',a8,',',a8,',',a8,/,  &
     &     ' KEY SPECIES =',i4,2x,a8,'  ODD-H=',2f8.2,'  PRONOX=',f8.2) 
          write(c_out,1002) icc, c_nicreac(nr), c_nnrchem(icc)          &
     &                      , c_nnrchp(icc)                             
 1002     format('ICC,NICREAC, NNRCHEM NNRCHMP=',4i5) 
          if(c_nnrchem(icc).gt.0) write(c_out,1003)                     &
     &               c_nrchem(icc,c_nnrchem(icc)),c_stoiloss(nr)        
 1003     format('NRCHEM, STOILOSS=',i5,f8.2) 
          if(c_nnrchp(icc).gt.0)                                        &
     &       write(c_out,1004) c_nrchmp(icc,c_nnrchp(icc))              
 1004     format('NRCHEMP=',2i5) 
        endif 
                                                                        
! ADDITION FOR ODD HYDROGEN:  IDENTIFY THE HO2->OH DIRECT REACTIONS.    
!                                                                       
!        if( (c_npequil(c_reactant(nr,1)).eq.c_nho2     .or.            
!    *        c_npequil(c_reactant(nr,2)).eq.c_nho2     ).and.          
!    *        c_nnpro(nr).gt.0)    then                                 
!                                                                       
         if(c_nnpro(nr).gt.0)      then 
                if(c_npequil(c_reactant(nr,1)).ne.c_nho2     ) then 
                  if(c_reactant(nr,2).le.0) then 
                    go to 431 
                  else 
                    if(c_reactant(nr,2).ne.c_nho2     ) go to 431 
                  endif 
                endif 
          do 430 n=1,c_nnpro(nr) 
            if( c_npequil(c_product(nr,n)).eq.c_noh      ) then 
              go to 431 
            endif 
  430     continue 
  431     continue 
                                                                        
         endif 
                                                                        
 1000   continue 
! END REACTION LOOP - END MECHANISM ANALYSIS.                           
                                                                        
                                                                        
! ----------------------------------------------------------            
! IDENTIFY 'EXCHANGE' REACTIONS WITH THE FORM A+B=>C; C=>A+B.           
! AND PAIRFAC/MULTIFAC CONSERVATION INDICES.                            
!                                                                       
! THESE REACTION NUMBERS ARE RECORDED IN EXPRO(IC,I) AND EXLOSS(IC,I)   
! WHERE IC IS FOR THE PRODUCT C  (GAS-PHASE).                           
                                                                        
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
                                                                        
! ----------------------------------------------------------            
!  ZERO THE IMPORTANT ARRAYS                                            
       do  ic=1,c_nchem2 
        do i=1,20 
         c_exspec(ic,i) = 0 
         do ii=1,5 
           c_expro(ic,i,ii) = 0 
           c_exloss(ic,i,ii) = 0 
                             !do ii=1,5                                 
         enddo 
                          !do i=1,20                                    
        enddo 
                        !do  ic=1,c_nchem2                              
       enddo 
                                                                        
! INITIALIZE PAIRFAC (conservation parameter for pair species)          
!   = 1 unless adjusted by EXCHANGE REACTIONS                           
                                                                        
       do ic=1,c_nchem2 
         c_pairfac(ic) = 1. 
         c_multfac(ic) = 1. 
                         !do ic=1,c_nchem2                              
       enddo 
!                                                                       
! --------                                                              
! BEGIN REACTION LOOP FOR EXCHANGE REACTIONS (1200)                     
                                                                        
       do nrx=1,c_nreac 
                                                                        
! EXLOSS REACTION:                                                      
!  ICX (NICREAC) is key species for reaction (1st non-special reactant).
          icx = c_nicreac(nrx) 
                                                                        
!  Identify up to two REACTANTS ICY1, ICY2;                             
!               and up to two PRODUCTS ICX1,ICX2.                       
                                                                        
          icy1 = c_npequil(c_reactant(nrx,1)) 
          icy2 = 0 
          if (c_reactant(nrx,2).gt.0) then 
            icy2 = c_npequil(c_reactant(nrx,2)) 
          endif 
                                                                        
          icx1 = 0 
          icx2 = 0 
          if(c_product(nrx,1).gt.0) icx1=c_npequil(c_product(nrx,1)) 
          if(c_product(nrx,2).gt.0) icx2=c_npequil(c_product(nrx,2)) 
                                                                        
          icat1 = 0 
          icat2 = 0 
          icatp = 0 
          if(icy1.gt.0) then 
             icat1 = c_icat(icy1) 
          endif 
          if(icy2.gt.0) then 
             icat2 = c_icat(icy2) 
          endif 
          if(icx1.gt.0) then 
             icatp = c_icat(icx1) 
          endif 
                                                                        
                                                                        
! CONTROLS TO SKIP LOOP (1200):                                         
!                                                                       
! SKIP AND CONTINUE LOOP IF (a) TWO REACTANTS, ONE PRODUCT;             
!  OR IF (b) KEY SPECIES IS A PRODUCT SPECIES.  (IN WHICH CASE          
!  THIS REACTION WOULD BE THE EXCHANGE 'RETURN' - EXPRO, not EXLOSS)    
          if(icy2.gt.0.and.icx2.eq.0) then 
          else 
           if(icx.ne.icy1.and.icx.ne.icy2.and.icy2.gt.0) then 
           else 
                                                                        
! OPTIONAL CRITERIA:  EXCHANGE FOR FAST SPECIES ONLY.  FOR SLOW SPECIES,
! SKIP AND CONTINUE LOOP                                                
            if(c_icat(icy1).eq.1.or.c_icat(icy1).eq.4) then 
            else 
                                                                        
! SEARCH THROUGH REACTIONS TO FIND MATCHING EXPRO REACTION(s) (1300)    
              do nrp=1,c_nreac 
                                                                        
!  ICP (NICREAC) is key species for EXPRO reaction                      
               icp = c_nicreac(nrp) 
                                                                        
!  Identify EXPRO reactants IC1,  IC2 and products ICP1, ICP2.          
               ic1=0 
               ic2=0 
               if(c_reactant(nrp,1).gt.0) then 
                  ic1=c_npequil(c_reactant(nrp,1)) 
               endif 
               if(c_reactant(nrp,2).gt.0) then 
                  ic2=c_npequil(c_reactant(nrp,2)) 
               endif 
               icp1=0 
               icp2=0 
               if(c_product(nrp,1).gt.0) then 
                  icp1=c_npequil(c_product(nrp,1)) 
               endif 
               if(c_product(nrp,2).gt.0) then 
                  icp2=c_npequil(c_product(nrp,2)) 
               endif 
                                                                        
! CONTROLS TO IDENTIFY  EXCHANGE REACTION (EXLOSS-EXPRO):  A+B=C, C=A+B 
               if( (icy1.eq.icp1.and.icy2.eq.icp2).or.                  &
     &                  (icy1.eq.icp2.and.icy2.eq.icp1) )  then         
                 if( (icx1.eq.ic1.and.icx2.eq.ic2).or.                  &
     &                  (icx1.eq.ic2.and.icx2.eq.ic1) )  then           
                                                                        
! SET PAIRFAC (pair group conservation of mass)                         
!   FOR EXCHANGE REACTIONS A+B=C                                        
!    (note A+B=>C must be return reaction)                              
                                                                        
                 if(ic1.ne.0.and.ic2.ne.0.and.icp1.ne.0) then 
                   if(c_nppair(ic1,2).eq.c_nppair(ic2,2)) then 
                     if(c_nppair(icp1,2).eq.c_nppair(ic1,2)) then 
                       if(icp2.eq.0) then 
                          c_pairfac(icp1) =                             &
     &                       c_pairfac(ic1)+c_pairfac(ic2)              
! TEST WRITE - QUADINIT PAIRFAC (STANDARD)                              
                          if(c_kkw.gt.0) write(c_out,1209)              &
     &                     c_tchem(ic1), c_tchem(ic2), c_tchem(icp1),   &
     &                     c_pairfac(ic1), c_pairfac(ic2),              &
     &                     c_pairfac(icp1)                              
 1209                      format(' SET PAIRFAC:  species A + B <->C.', &
     &                            ' SPECIES AND PAIRFACS = ',/,         &
     &                            3(a8,2x), 3(1pe10.3)                ) 
                       else 
                          if(c_nppair(icp2,2).ne.c_nppair(ic1,2)) then 
                            c_pairfac(icp2) =                           &
     &                         c_pairfac(ic1)+c_pairfac(ic2)            
! TEST WRITE - QUADINIT PAIRFAC (STANDARD)                              
                            if(c_kkw.gt.0) write(c_out,1209)            &
     &                       c_tchem(ic1), c_tchem(ic2), c_tchem(icp1)  
                          endif 
                                           !if(icp2.eq.0) then          
                       endif 
                                  !if(c_nppair(icp1,2).eq....           
                     else 
                       if(icp2.gt.0) then 
                          if(c_nppair(icp2,2).eq.c_nppair(ic1,2)) then 
                            c_pairfac(icp2) =                           &
     &                         c_pairfac(ic1)+c_pairfac(ic2)            
! TEST WRITE - QUADINIT PAIRFAC (STANDARD)                              
                            if(c_kkw.gt.0) write(c_out,1209)            &
     &                       c_tchem(ic1), c_tchem(ic2), c_tchem(icp1)  
                          endif 
                                           !if(icp2.gt.0) then          
                       endif 
                                     !if(c_nppair(icp1,2).eq....        
                     endif 
                                  !if(c_nppair(ic1,2).eq....            
                   endif 
                              !if(ic1.ne.0.and.ic2.ne.0.and.icp1.ne.0)..
                 endif 
                                                                        
!                                                                       
! SET MULTIFAC (multi group conservation of mass)                       
!   FOR EXCHANGE REACTIONS A+B=C                                        
!    (note A+B=>C must be return reaction)                              
                                                                        
                 if(ic1.ne.0.and.ic2.ne.0.and.icp1.ne.0) then 
                   if(c_npmulti(ic1,1).eq.c_npmulti(ic2,1)) then 
                     if(c_npmulti(icp1,1).eq.c_npmulti(ic1,1)) then 
                       if(icp2.eq.0) then 
                          c_multfac(icp1) =                             &
     &                       c_multfac(ic1)+c_multfac(ic2)              
! TEST WRITE - QUADINIT MULTIFAC (STANDARD)                             
                          if(c_kkw.gt.0) write(c_out,1208)              &
     &                     c_tchem(ic1), c_tchem(ic2), c_tchem(icp1),   &
     &                     c_multfac(ic1), c_multfac(ic2),              &
     &                               c_multfac(icp1)                    
 1208                      format(' SET MULTIFAC:  species A + B <->C.',&
     &                            ' SPECIES AND MULTIFACS = ',/,        &
     &                            3(a8,2x), 3(1pe10.3)                ) 
                       else 
                          if(c_npmulti(icp2,1).ne.c_npmulti(ic1,1)) then 
                           c_multfac(icp2) =                            &
     &                          c_multfac(ic1)+c_multfac(ic2)           
! TEST WRITE - QUADINIT MULTIFAC (STANDARD)                             
                            if(c_kkw.gt.0) write(c_out,1208)            &
     &                       c_tchem(ic1), c_tchem(ic2), c_tchem(icp1)  
                          endif 
                                           !if(icp2.eq.0) then          
                       endif 
                                     !if(c_npmulti(icp1,1).eq....       
                     else 
                       if(icp2.gt.0) then 
                          if(c_npmulti(icp2,1).eq.c_npmulti(ic1,1)) then 
                           c_multfac(icp2) =                            &
     &                        c_multfac(ic1)+c_multfac(ic2)             
! TEST WRITE - QUADINIT MULTIFAC (STANDARD)                             
                            if(c_kkw.gt.0) write(c_out,1208)            &
     &                       c_tchem(ic1), c_tchem(ic2), c_tchem(icp1)  
                          endif 
                                           !if(icp2.gt.0) then          
                       endif 
                                     !if(c_npmulti(icp1,1).eq....       
                     endif 
                                  !if(c_npmulti(ic1,1).eq....           
                   endif 
                              !if(ic1.ne.0.and.ic2.ne.0.and.icp1.ne.0)  
                 endif 
                                                                        
! SET LOGIC FLAGS TO TEST IF REACTION IS ALREADY LISTED                 
                   lloss = .TRUE. 
                   lpro  = .TRUE. 
                                                                        
! LOOP TO TEST WHETHER EXCHANGE REACTION IS ALREADY LISTED IN REVERSE   
!   If listed, EXPRO reaction would be  EXLOSS for its key species      
                                                                        
                   if(icp.gt.0) then 
                    do i=1,20 
                     do ii=1,5 
! IF EXCHANGE RN IS LISTED IN REVERSE, USE LPRO, LLOSS to ID reaction as
                      if(c_exloss(icp,i,ii).eq.nrp) then 
                       lloss= .FALSE. 
                       lpro= .FALSE. 
                      endif 
                                    !do ii=1,5                          
                     enddo 
                                 !do i=1,20                             
                    enddo 
                               !if(icp.gt.0)                            
                   endif 
                                                                        
! LOOP TO TEST WHETHER EXLOSS, PRO REACTIONS ARE ALREADY IN FORWARD LIST
                                                                        
                   do i=1,20 
                    do ii=1,5 
! IF EXCHANGE RN IS LISTED IN REVERSE, USE LPRO, LLOSS to ID reaction as
                     if(c_exloss(icx,i,ii).eq.nrx) lloss = .FALSE. 
                     if( c_expro(icx,i,ii).eq.nrp)  lpro = .FALSE. 
                                   !do ii=1,5                           
                    enddo 
                                !do i=1,20                              
                   enddo 
                                                                        
! IDENTIFY EXCHANGED SPECIES (icp):                                     
!  Record species only if base and exchanged species                    
!                                         are in same pair group.       
!   and if exchanged  species is in reaction product list               
!  If no exchanged species, flag  species number as -1.                 
                   if(icp.gt.0) then 
                    if( icp.ne.icx1.and.icp.ne.icx2)   icp=0 
                    if(c_nppair(icx,2).ne.c_nppair(icp,2)) icp=0 
                                   !if(icp.gt.0)                        
                   endif 
                   if(icp.eq.0)  icp=-1 
                                                                        
! RECORD EXCHANGE SPECIES ONLY IF MAIN AND EXCHANGE SPECIES ARE PAIRED  
!  AND IF KEY SPECIES OF 2ND REACTION IS PRODUCT OF 1ST                 
!                                                                       
!  ( EXCORR must be called for single species also, e.g. HNO4.          
                                                                        
! LOOP TO ADD EXCHANGED SPECIES TO EXSPEC LIST                          
                                                                        
                  if(lloss.or.lpro) then 
                    do i=1,20 
                     is=i 
                     if(c_exspec(icx,i).eq.icp) go to 1205 
                     if(c_exspec(icx,i).eq.0) then 
                       c_exspec(icx,i)  = icp 
                       go to 1205 
                     endif 
                     if(i.eq.20) write(c_out,1206) 
                     if(i.eq.20) write(  6 ,1206) 
 1206                format(/,' MAJOR ERROR IN QUADINIT: ',             &
     &                'EXCHANGED SPECIES INDEX EXCEEDED LIMIT' ,        &
     &                ' (i=1,20, 1205)',/)                              
                                    !   do i=1,20                       
                    enddo 
 1205               continue 
                                 !if(lloss.or.lpro)                     
                  endif 
                                                                        
! LOOP TO ADD EXLOSS REACTION                                           
                  if(lloss) then 
                   do ii=1,5 
! ENTER THE EXCHANGE REACTION                                           
                    if(c_exloss(icx,is,ii).eq.0) then 
                      c_exloss(icx,is,ii)=nrx 
! TEST WRITE - QUADINIT EXCHANGE(STANDARD)                              
                      if(c_kkw.gt.0)                                    &
     &                   write(c_out,1212) icx,icp, is,i,nrx,           &
     &                   c_reactant(nrx,1),c_reactant(nrx,2),           &
     &                   c_product(nrx,1),c_product(nrx,2)              
 1212                    format (/,'EXLOSS(ICX,IS,I)=NRX. ',            &
     &                   'ICX,ICP,IS,I,NRX=', 5i5,/,                    &
     &                   ' REACTANTS=',2i5, '  PRODUCTS=',2i5   )       
                      if(c_kkw.gt.0.and.c_reactant(nrx,1).gt.0)         &
     &                 write(c_out,*) c_tchem(c_reactant(nrx,1))        
                      if(c_kkw.gt.0.and.c_reactant(nrx,2).gt.0)         &
     &                 write(c_out,*) c_tchem(c_reactant(nrx,2))        
                      if(c_kkw.gt.0.and.c_product(nrx,1).gt.0)          &
     &                 write(c_out,*) c_tchem( c_product(nrx,1))        
                      if(c_kkw.gt.0.and. c_product(nrx,2).gt.0)         &
     &                 write(c_out,*) c_tchem( c_product(nrx,2))        
                                                                        
! HAVING ENTERED THE REACTION, BREAK THE LOOP                           
                      go to 1211 
! ERROR WRITE                                                           
                      write(c_out,1213) 
                      write( 6  ,1213) 
 1213                format(/,' MAJOR ERROR IN QUADINIT: ',             &
     &                'EXLOSS INDEX EXCEEDED LIMIT' ,                   &
     &                ' (ii=1,5, 1213)',/)                              
                                                                        
                                    ! if(c_exloss(icx,is,ii).eq.0       
                    endif 
                              !do ii=1,5 (1210)                         
                   enddo 
                                !if(lloss)                              
                  endif 
! END LOOP TO ADD EXLOSS REACTION (1210)                                
!    CONTINUE AFTER BREAK LOOP                                          
 1211             continue 
!                                                                       
                                                                        
! LOOP TO ADD EXPRO  REACTION                                           
                  if(lpro ) then 
                   do ii=1,5 
! ENTER THE EXCHANGE REACTION                                           
                    if( c_expro(icx,is,ii).eq.0) then 
                       c_expro(icx,is,ii)=nrp 
! TEST WRITE - QUADINIT EXCHANGE(STANDARD)                              
                      if(c_kkw.gt.0)                                    &
     &                   write(c_out,1222) icx,icp, is,i,nrp,           &
     &                   c_reactant(nrp,1),c_reactant(nrp,2),           &
     &                   c_product(nrp,1),c_product(nrp,2)              
 1222                    format (/,' EXPRO(ICX,IS,I)=NRX.  ',           &
     &                   'ICX,ICP,IS,I,NRX=', 5i5,/,                    &
     &                   ' REACTANTS=',2i5,'  PRODUCTS=',2i5   )        
                      if(c_kkw.gt.0.and.c_reactant(nrp,1).gt.0)         &
     &                 write(c_out,*) c_tchem(c_reactant(nrp,1))        
                      if(c_kkw.gt.0.and.c_reactant(nrp,2).gt.0)         &
     &                 write(c_out,*) c_tchem(c_reactant(nrp,2))        
                      if(c_kkw.gt.0.and.c_product(nrp,1).gt.0)          &
     &                 write(c_out,*) c_tchem( c_product(nrp,1))        
                      if(c_kkw.gt.0.and. c_product(nrp,2).gt.0)         &
     &                 write(c_out,*) c_tchem( c_product(nrp,2))        
                                                                        
! ENTER WARNING FOR EXCHANGE REACTION AMONG UNPAIRED SPECIES            
                      if(icp.lt.0) then 
                        if ( (c_reactant(nrp,1).gt.11)                  &
     &                       .or. (c_reactant(nrp,2).gt.11)) then       
!                        write(6    ,1224) nrp,  (c_treac(j,nrp),j=1,5) 
                         write(c_out,1224) nrp,  (c_treac(j,nrp),j=1,5) 
 1224                    format(/,'WARNING: POSSIBLE EXCHANGE REACTION '&
     &                           ,'AMONG UNPAIRED SPECIES: nr =',       &
     &                          /, i5, a8,'+',a8,'=>',a8,'+',a8,'+', a8)
                        endif 
                      endif 
                                                                        
! HAVING ENTERED THE REACTION, BREAK THE LOOP                           
                      go to 1221 
! ERROR WRITE                                                           
                      write(c_out,1223) 
                      write( 6  ,1223) 
 1223                format(/,' MAJOR ERROR IN QUADINIT: ',             &
     &                'EXPRO INDEX EXCEEDED LIMIT' ,                    &
     &                ' (ii=1,5, 1223)',/)                              
                                                                        
                    endif 
                              !do ii=1,5 (1220)                         
                   enddo 
                               !if(lpro )                               
                  endif 
! END LOOP TO ADD EXPRO  REACTION (1220)                                
!    CONTINUE AFTER BREAK LOOP                                          
 1221             continue 
                                                                        
                                                                        
! END CONTROLS TO IDENTIFY EXCHANGE REACTIONS                           
                          !if( (icx1.eq.ic1.and.icx2.eq.ic2) etc.       
                 endif 
                        !if( (icy1.eq.icp1.and.icy2.eq.icp2) etc.       
               endif 
                                                                        
                        !do nrp=1,c_nreac                               
              enddo 
! END LOOP - SEARCH THROUGH REACTIONS TO FIND MATCHING EXPRO REACTION (1
                                                                        
                                                                        
! END - CONTROLS TO SKIP AND CONTINUE LOOP                              
                         !if(c_icat(icy1).eq.1.or.c_icat(icy1).eq.4)    
            endif 
                       !if(icx.ne.icy1.and.icx.ne.icy2.and.icy2.gt.0)   
           endif 
                     !if(icy2.gt.0.and.icx2.eq.0)                       
          endif 
                                                                        
                 !do nrx=1,c_nreac                                      
       enddo 
! END   REACTION LOOP FOR EXCHANGE REACTIONS (1200)                     
                                                                        
                                                                        
! END EXCHANGE ANALYSIS.                                                
! ----------------------                                                
!                                                                       
!                                                                       
!   IDENTIFY 'CROSS' REACTIONS FOR PAIRED SPECIES IN CASCADE            
! NNRCPRO, NRCPRO, STOICPRO WITH SAME FORMAT AS NRCHEM, ABOVE.          
!  THESE REPRESENT REACTIONS THAT PRODUCE THE SPECIES FROM 1/2 PAIR.    
                                                                        
! NNRCPRX, NRCPRX, STOICPRX REPRESENT PRODUCTION REACTIONS ASSOCIATED   
!  WITH THE THIRD PAIR IN A TRIPLET - USED ONLY FOR O3-NO-NO2.          
!                                                                       
! JANUARY 2005:  THIS HAS BEEN DELETED.   IT RELATES TO OLD TWOSOLVE    
!  AND NOXSOLVE - REPLACED BY NEW PAIR AND MULTISOLVE TREATMENTS.       
!                                                                       
!                                                                       
!                                                                       
! END CROSS-REACTION ANALYSIS.  END MECHANISM ANALYSIS.                 
! ----------------------                                                
!                                                                       
! ----------------------                                                
!  NOV 2007:  Identify net RP for NOx, Ox tracers                       
! IDENTIFY REACTION NUMBERS FOR 'SPECIAL' REACTIONS                     
!  (used in preliminary chem, oh solver and in NOx-Ox global tracers)   
! ----------------------                                                
                                                                        
! ZERO                                                                  
      do i=1,13 
       do nr=1,c_nreac 
        c_noxchem(i,nr) = 0. 
       enddo 
      enddo 
      c_nro3no = 0 
      c_nrno2x = 0 
      c_nro3hv=0 
      c_nrohno2 =0 
      c_nrohco =0 
      c_nrho2no=0 
      c_nrho22 =0 
      c_nro3no2 =0 
      c_nroho3 =0 
      c_nrohch4 =0 
                                                                        
! SPECIAL REACTIONS: LOOP FOR EACH INDIVIDUAL REACTION                  
!    moved from presolve                                                
                                                                        
      do nr=1,c_nreac 
!  c_nro3no                                                             
                                                                        
        if(                                                             &
     &    (c_reactant(nr,1).eq.c_nno.and.c_reactant(nr,2).eq.c_no3      &
     &     .and.c_product(nr,1).eq.c_nno2)                              &
     &    .or.                                                          &
     &    (c_reactant(nr,2).eq.c_nno.and.c_reactant(nr,1).eq.c_no3      &
     &     .and.c_product(nr,1).eq.c_nno2)                              &
     &      ) c_nro3no = nr                                             
                                                                        
! c_nrno2hv                                                             
        if(                                                             &
     &    (c_product (nr,1).eq.c_nno.and.c_product (nr,2).eq.c_no3      &
     &     .and.c_reactant(nr,1).eq.c_nno2)                             &
     &    .or.                                                          &
     &    (c_product (nr,2).eq.c_nno.and.c_product (nr,1).eq.c_no3      &
     &     .and.c_reactant(nr,1).eq.c_nno2)                             &
     &      )c_nrno2x = nr                                              
                                                                        
         if(c_reactant(nr,1).eq.c_no3.and.c_reactant(nr,2).eq.-1        &
     &       .and.c_product(nr,1).eq.c_noh)                             &
     &    c_nro3hv   =nr                                                
                                                                        
         if( (c_reactant(nr,1).eq.c_nno2.and.c_reactant(nr,2).eq.c_noh) &
     &   .or.(c_reactant(nr,1).eq.c_noh.and.c_reactant(nr,2).eq.c_nno2))&
     &    c_nrohno2   =nr                                               
                                                                        
         if( (c_reactant(nr,1).eq.c_nco .and.c_reactant(nr,2).eq.c_noh) &
     &   .or.(c_reactant(nr,1).eq.c_noh.and.c_reactant(nr,2).eq.c_nco ))&
     &    c_nrohco    =nr                                               
                                                                        
                                                                        
         if( (c_reactant(nr,1).eq.c_nno.and.c_reactant(nr,2).eq.c_nho2) &
     &   .or.(c_reactant(nr,1).eq.c_nho2.and.c_reactant(nr,2).eq.c_nno))&
     &    c_nrho2no   =nr                                               
                                                                        
         if(c_reactant(nr,1).eq.c_nho2.and.c_reactant(nr,2).eq.c_nho2)  &
     &    c_nrho22    =nr                                               
                                                                        
         if( (c_reactant(nr,1).eq.c_nno2.and.c_reactant(nr,2).eq.c_no3) &
     &   .or.(c_reactant(nr,1).eq.c_no3.and.c_reactant(nr,2).eq.c_nno2))&
     &    c_nro3no2    =nr                                              
                                                                        
         if( (c_reactant(nr,1).eq.c_no3 .and.c_reactant(nr,2).eq.c_noh) &
     &   .or.(c_reactant(nr,1).eq.c_noh.and.c_reactant(nr,2).eq.c_no3 ))&
     &    c_nroho3    =nr                                               
                                                                        
         if( (c_reactant(nr,1).eq.c_nch4.and.c_reactant(nr,2).eq.c_noh) &
     &   .or.(c_reactant(nr,1).eq.c_noh.and.c_reactant(nr,2).eq.c_nch4))&
     &    c_nrohch4   =nr                                               
                                                                        
                                                                        
                                                                        
                                                                        
                        !do nr=1,c_nreac                                
      enddo 
                                                                        
! TRACER ANALYSIS:                                                      
! LOOP FOR EACH INDIVIDUAL REACTION                                     
       do  nr=1,c_nreac 
                                                                        
                                                                        
!  noxchem(i,nr) = net rp O3 NOx PANs HNO3 RNO3;  NOx-PAN, PAN-NOX; HNO3
! REACTANT1                                                             
        if( c_icat(c_reactant(nr,1)).eq.11                              &
     &    .or.  c_icat(c_reactant(nr,1)).eq.12                          &
     &                   )  c_noxchem(1,nr) = c_noxchem(1,nr) - 1.      
                                                                        
        if( c_icat(c_reactant(nr,1)).eq.6                               &
     &    .or.  (c_icat(c_reactant(nr,1)).ge.12.and.                    &
     &           c_icat(c_reactant(nr,1)).le.15)                        &
     &                   )  c_noxchem(2,nr) = c_noxchem(2,nr) - 1.      
                                                                        
        if( c_icat(c_reactant(nr,1)).eq. 5.                             &
     &                   )  c_noxchem(3,nr) = c_noxchem(3,nr) - 1.      
                                                                        
        if( c_icat(c_reactant(nr,1)).eq.16.                             &
     &                   )  c_noxchem(4,nr) = c_noxchem(4,nr) - 1.      
                                                                        
        if( c_icat(c_reactant(nr,1)).eq. 7.                             &
     &                   )  c_noxchem(5,nr) = c_noxchem(5,nr) - 1.      
                                                                        
         if(c_tchem(c_reactant(nr,1)).eq.'    R4N1'                     &
     &                   )  c_noxchem(5,nr) = c_noxchem(5,nr) - 1.      
                                                                        
         if(c_tchem(c_reactant(nr,1)).eq.'    R3N1'                     &
     &                   )  c_noxchem(5,nr) = c_noxchem(5,nr) - 1.      
                                                                        
         if(c_tchem(c_reactant(nr,1)).eq.'    R7N1'                     &
     &                   )  c_noxchem(5,nr) = c_noxchem(5,nr) - 1.      
                                                                        
         if(c_tchem(c_reactant(nr,1)).eq.'    R6N1'                     &
     &                   )  c_noxchem(5,nr) = c_noxchem(5,nr) - 1.      
                                                                        
         if(c_tchem(c_reactant(nr,1)).eq.'    R5N1'                     &
     &                   )  c_noxchem(5,nr) = c_noxchem(5,nr) - 1.      
                                                                        
         if(c_tchem(c_reactant(nr,1)).eq.'    ISN1'                     &
     &                   )  c_noxchem(5,nr) = c_noxchem(5,nr) - 1.      
                                                                        
         if(c_tchem(c_reactant(nr,1)).eq.'    ISNR'                     &
     &                   )  c_noxchem(5,nr) = c_noxchem(5,nr) - 1.      
                                                                        
! REACTANT2                                                             
                                                                        
       if(c_reactant(nr,2).gt.0) then 
        if( c_icat(c_reactant(nr,2)).eq.11                              &
     &    .or.  c_icat(c_reactant(nr,2)).eq.12                          &
     &                   )  c_noxchem(1,nr) = c_noxchem(1,nr) - 1.      
                                                                        
        if( c_icat(c_reactant(nr,2)).eq.6                               &
     &    .or.  (c_icat(c_reactant(nr,2)).ge.12.and.                    &
     &           c_icat(c_reactant(nr,2)).le.15)                        &
     &                   )  c_noxchem(2,nr) = c_noxchem(2,nr) - 1.      
                                                                        
        if( c_icat(c_reactant(nr,2)).eq. 5.                             &
     &                   )  c_noxchem(3,nr) = c_noxchem(3,nr) - 1.      
                                                                        
        if( c_icat(c_reactant(nr,2)).eq.16.                             &
     &                   )  c_noxchem(4,nr) = c_noxchem(4,nr) - 1.      
                                                                        
        if( c_icat(c_reactant(nr,2)).eq. 7.                             &
     &                   )  c_noxchem(5,nr) = c_noxchem(5,nr) - 1.      
                                                                        
         if(c_tchem(c_reactant(nr,2)).eq.'    R4N1'                     &
     &                   )  c_noxchem(5,nr) = c_noxchem(5,nr) - 1.      
                                                                        
         if(c_tchem(c_reactant(nr,2)).eq.'    R3N1'                     &
     &                   )  c_noxchem(5,nr) = c_noxchem(5,nr) - 1.      
                                                                        
         if(c_tchem(c_reactant(nr,2)).eq.'    R7N1'                     &
     &                   )  c_noxchem(5,nr) = c_noxchem(5,nr) - 1.      
                                                                        
         if(c_tchem(c_reactant(nr,2)).eq.'    R6N1'                     &
     &                   )  c_noxchem(5,nr) = c_noxchem(5,nr) - 1.      
                                                                        
         if(c_tchem(c_reactant(nr,2)).eq.'    R5N1'                     &
     &                   )  c_noxchem(5,nr) = c_noxchem(5,nr) - 1.      
                                                                        
         if(c_tchem(c_reactant(nr,2)).eq.'    ISN1'                     &
     &                   )  c_noxchem(5,nr) = c_noxchem(5,nr) - 1.      
                                                                        
         if(c_tchem(c_reactant(nr,2)).eq.'    ISNR'                     &
     &                   )  c_noxchem(5,nr) = c_noxchem(5,nr) - 1.      
                                 !if(c_reactant(nr,2).gt.0) then        
        endif 
                                                                        
! PRODUCTS                                                              
                                                                        
       do  n=1,c_nnpro(nr) 
       if(c_product (nr,n).gt.0) then 
        if( c_icat(c_product (nr,n)).eq.11                              &
     &    .or.  c_icat(c_product (nr,n)).eq.12                          &
     &           )  c_noxchem(1,nr) = c_noxchem(1,nr) + c_stoich(nr,n)  
                                                                        
        if( c_icat(c_product (nr,n)).eq.6                               &
     &    .or.  (c_icat(c_product (nr,n)).ge.12.and.                    &
     &           c_icat(c_product (nr,n)).le.15)                        &
     &           )  c_noxchem(2,nr) = c_noxchem(2,nr) + c_stoich(nr,n)  
                                                                        
        if( c_icat(c_product (nr,n)).eq. 5.                             &
     &           )  c_noxchem(3,nr) = c_noxchem(3,nr) + c_stoich(nr,n)  
                                                                        
        if( c_icat(c_product (nr,n)).eq.16.                             &
     &           )  c_noxchem(4,nr) = c_noxchem(4,nr) + c_stoich(nr,n)  
                                                                        
        if( c_icat(c_product (nr,n)).eq. 7.                             &
     &           )  c_noxchem(5,nr) = c_noxchem(5,nr) + c_stoich(nr,n)  
                                                                        
         if(c_tchem(c_product (nr,n)).eq.'    R4N1'                     &
     &           )  c_noxchem(5,nr) = c_noxchem(5,nr) + c_stoich(nr,n)  
                                                                        
         if(c_tchem(c_product (nr,n)).eq.'    R3N1'                     &
     &           )  c_noxchem(5,nr) = c_noxchem(5,nr) + c_stoich(nr,n)  
                                                                        
         if(c_tchem(c_product (nr,n)).eq.'    R7N1'                     &
     &           )  c_noxchem(5,nr) = c_noxchem(5,nr) + c_stoich(nr,n)  
                                                                        
         if(c_tchem(c_product (nr,n)).eq.'    R6N1'                     &
     &           )  c_noxchem(5,nr) = c_noxchem(5,nr) + c_stoich(nr,n)  
                                                                        
         if(c_tchem(c_product (nr,n)).eq.'    R5N1'                     &
     &           )  c_noxchem(5,nr) = c_noxchem(5,nr) + c_stoich(nr,n)  
                                                                        
         if(c_tchem(c_product (nr,n)).eq.'    ISN1'                     &
     &           )  c_noxchem(5,nr) = c_noxchem(5,nr) + c_stoich(nr,n)  
                                                                        
         if(c_tchem(c_product (nr,n)).eq.'    ISNR'                     &
     &           )  c_noxchem(5,nr) = c_noxchem(5,nr) + c_stoich(nr,n)  
                                 !if(c_product (nr,n).gt.0) then        
        endif 
                                 !do  n=1,c_nnpro(nr)                   
        enddo 
                                                                        
! ENTER STOICHIOMETRIES FOR                                             
!    6 and 7:  NOx=>PAN, PAN=>NOx                                       
!    8 and 9:  NOx=>HNO3, HNO3=>NOX                                     
!   10 and 11:  NOx=>RNO3;  RNO3=>NOx                                   
!   12 and 13:  NOx=>H+RNO3; H+RNO3=>NOx                                
                                                                        
       if(c_noxchem(3,nr).gt.0) c_noxchem(6,nr) = c_noxchem(3,nr) 
       if(0.-c_noxchem(2,nr).lt.c_noxchem(6,nr) )                       &
     &       c_noxchem(6,nr) = 0.-c_noxchem(2,nr)                       
       if(c_noxchem(6,nr).lt.0.) c_noxchem(6,nr)=0. 
                                                                        
       if(c_noxchem(3,nr).lt.0) c_noxchem(7,nr) =0.- c_noxchem(3,nr) 
       if(   c_noxchem(2,nr).lt.c_noxchem(7,nr) )                       &
     &       c_noxchem(7,nr) =    c_noxchem(2,nr)                       
       if(c_noxchem(7,nr).lt.0.) c_noxchem(7,nr)=0. 
                                                                        
       if(c_noxchem(4,nr).gt.0) c_noxchem(8,nr) = c_noxchem(4,nr) 
       if(0.-c_noxchem(2,nr).lt.c_noxchem(8,nr) )                       &
     &       c_noxchem(8,nr) = 0.-c_noxchem(2,nr)                       
       if(c_noxchem(8,nr).lt.0.) c_noxchem(8,nr)=0. 
                                                                        
       if(c_noxchem(4,nr).lt.0) c_noxchem(9,nr) =0.- c_noxchem(4,nr) 
       if(   c_noxchem(2,nr).lt.c_noxchem(9,nr) )                       &
     &       c_noxchem(9,nr) =    c_noxchem(2,nr)                       
       if(c_noxchem(9,nr).lt.0.) c_noxchem(9,nr)=0. 
                                                                        
       if(c_noxchem(5,nr).gt.0) c_noxchem(10,nr) = c_noxchem(5,nr) 
       if(0.-c_noxchem(2,nr).lt.c_noxchem(10,nr) )                      &
     &       c_noxchem(10,nr) = 0.-c_noxchem(2,nr)                      
       if(c_noxchem(10,nr).lt.0.) c_noxchem(10,nr)=0. 
                                                                        
       if(c_noxchem(5,nr).lt.0) c_noxchem(11,nr) =0.- c_noxchem(5,nr) 
       if(   c_noxchem(2,nr).lt.c_noxchem(11,nr) )                      &
     &       c_noxchem(11,nr) =    c_noxchem(2,nr)                      
       if(c_noxchem(11,nr).lt.0.) c_noxchem(11,nr)=0. 
                                                                        
       if((c_noxchem(4,nr)+c_noxchem(5,nr).gt.0) )                      &
     &         c_noxchem(12,nr) = c_noxchem(4,nr) + c_noxchem(5,nr)     
       if(0.-c_noxchem(2,nr).lt.c_noxchem(12,nr) )                      &
     &       c_noxchem(12,nr) = 0.-c_noxchem(2,nr)                      
       if(c_noxchem(12,nr).lt.0.) c_noxchem(12,nr)=0. 
                                                                        
       if((c_noxchem(4,nr)+c_noxchem(5,nr).lt.0) )                      &
     &         c_noxchem(13,nr) =0.- c_noxchem(4,nr) - c_noxchem(5,nr)  
       if(   c_noxchem(2,nr).lt.c_noxchem(13,nr) )                      &
     &       c_noxchem(13,nr) =    c_noxchem(2,nr)                      
       if(c_noxchem(13,nr).lt.0.) c_noxchem(13,nr)=0. 
                                                                        
! TEST WRITE - DEBUG   treac                                            
!      if(c_kkw .gt.0) then                                             
!        if(nr.eq.1) write(c_out,1231)                                  
! 1231     format(/,'TEST REACTION INDICES FOR OX, NOX, PAN, HNO3, RNO3:
!        write(c_out, 1233) nr,  (c_treac(j,nr ),j=1,5)                 
! 1233     format( i5, a8,'+',a8,'=>',a8,'+',a8,'+', a8)                
!        write(c_out,1235) (c_noxchem(n,nr),n=1,13)                     
! 1235     format(13f6.2)                                               
!        if(c_nro3no.eq.nr) write(c_out,*) 'c_nro3no =', nr             
!        if(c_nrno2x.eq.nr) write(c_out,*) 'c_nrno2x =', nr             
!      endif          ! if(c_kkw .gt.0) then                            
                                                                        
                                                                        
                             !do  nr=1,c_nreac                          
       enddo 
                                                                        
!                                                                       
                                                                        
! ----------------------                                                
!  END QUADINIT                                                         
 2000     return 
      END                                           
                                                                        
                                                                        
                                                                        
                                                                        
       subroutine chemwrit( kw) 
                                                                        
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
                                                                        
                                                                        
                                                                        
! -------------------------------------------------------------------   
      IMPLICIT NONE 
                                                                        
      include 'chemmech.EXT' 
      INCLUDE 'chemvars.EXT' 
      include 'chemlocal.EXT' 
                                                                        
                                  ! 'SUM' name for output               
        character*8 tsum 
                                  ! Gas phase concentration             
        double precision xcgas 
                                  ! Aqueous   concentration             
        double precision acquacon 
                                                                        
        if( kw.le.0) return 
        kk=1 
                                                                        
!  SPECIES CONCENTRATIONS.                                              
       write(c_out,1141)  c_hour 
 1141  format(/,' SPECIES CONCENTRATIONS AT  XHOUR =',  f6.2) 
       ncf = int(0.2*(float(c_nchem2)-0.01))+1 
       do      nc=1,ncf 
         nc1 = 4*(nc-1) + 1 
          nc2 = nc1 + 3 
          if(nc1.gt.c_nchem2)nc1=c_nchem2 
          if(nc2.gt.c_nchem2)nc2=c_nchem2 
           write(c_out,1142)                                            &
     &        ( c_tchem(nn),c_xcout(kw ,(nn)) , nn=nc1,nc2 )            
 1142    format(4(a8,1pe10.3,2x )) 
       enddo 
                                                                        
!   AQUEOUS SPECIES CONCENTRATIONS                                      
                                                                        
! NOTE  c_xcout(kw,ic) for GAS SPECIES IS GAS+AQ SUM, molec/cm3 equiv.  
!       c_xcout(kw,icq) for AQUEOUS SPECIES is moles/liter              
!    .  GAS-ONLY IS CALCULATED HERE.                                    
                                                                        
! NOTE!  WRITE error for lumped species such as  MP (in rooh);          
!        c_xcout(kw,ic) is preserved as species fraction, not sum?      
                                                                        
          if(c_h2oliq(kw ).gt.0) then 
           acquacon = c_h2oliq( kw)*avogadrl 
           write(c_out,1301)   c_h2oliq(kw ) 
           write(c_out,1303) acquacon 
           if(c_aqueous(1,2).gt.0) then 
           write(c_out,1302) c_xcout(kw ,c_aqueous(1,2)) 
          endif 
                                                                        
 1301       format(/,'AQUEOUS CHEMISTRY ',/,'WATER (grams/cm3) =',      &
     &    1pe10.3)                                                      
 1302       format( ' [H+] (moles per liter) =',1pe10.3) 
 1303       format( ' GAS SPECIES and G-A SUM molec/cm3. ',             &
     &      ' AQUEOUS moles/liter. CONVERSION= ', 1pe10.3)              
                                                                        
! 1303        format( ' GAS SPECIES is G+A sum, molec/cm3. ',           
!    *      ' AQUEOUS moles/liter. CONVERSION= ', 1pe10.3)              
                                                                        
            tsum = '     SUM' 
!          do      ic=1,c_nchem2                                        
           do  nrh=1,c_nreach 
             ic = c_henry(nrh,1) 
             if(c_nequil(ic).gt.0)  then 
               xcgas = c_xcout(kw,ic) 
             do i=1,c_nequil(ic) 
               xcgas = xcgas - c_xcout(kw,c_ncequil(ic,i))*acquacon 
             enddo 
!             write(c_out,1305) c_tchem(ic), c_xcout(kw ,ic),           
!    *        ( c_tchem(c_ncequil(ic,i)),c_xcout(kw ,c_ncequil(ic,i))   
!    *                                      ,i=1,c_nequil(ic) )         
               write(c_out,1305) c_tchem(ic), xcgas,                    &
     &        ( c_tchem(c_ncequil(ic,i)),c_xcout(kw ,c_ncequil(ic,i))   &
     &                                      ,i=1,c_nequil(ic) )         &
     &        ,tsum, c_xcout(kw,ic)                                     
 1305          format(4(a8,1pe10.3,2x),/,                               &
     &          '         0.000E+00  ',3(a8,1pe10.3,2x))                
               if(rhdif(kw,nrh).gt.0)                                   &
     &         xcgas = (c_xcout(kw,c_ncequil(ic,1))*acquacon)           &
     &         /rhdif(kw,nrh)                                           
              write(c_out,1306)xcgas, rateh(kw ,nrh), rhdif(kw ,nrh) 
 1306        format('   (',1pe10.3,                                     &
     &       ')                      HENRY, Hw/DIFF=', 3(1pe10.3))      
!             write(c_out,1306) rateh(kw ,nrh), rhdif(kw ,nrh)          
! 1306         format('                                ',               
!    *       '     HENRY, Hw/DIFF=',2(1pe10.3))                         
            endif 
           enddo 
         endif 
                                                                        
                                                                        
! --------------                                                        
! CHEMWRIT  SUMMARY WRITE:                                              
!    CONCENTRATIONS, XCIN, RP, RL from QUADCHEM                         
! --------------                                                        
                                                                        
        write(c_out,1801) c_iter, kw 
 1801   format(//,' SUMMARY WRITE:  ITER =',i3, '    VECTOR KKW=',i3) 
       write(c_out,1804) c_IDATE, c_hour,  c_lat(1), c_lon(1),c_temp(1) 
 1804  format('IDATE xhour lat lon temp=',i8,4f8.2) 
       write(c_out,1805) (c_jparam(i),i=1,13),c_jparam(20) 
 1805  format('JPARAMS: zenith altmid dobson so2 no2 aaerx aerssa',     &
     &  ' albedo',/,'cld-below cld-above claltB claltA temp date',      &
     &   /,(8f10.3))                                                    
        write(c_out,1802) 
 1802   format(/,' #  IC     XCOUT     XCIN  XCF/AV       RL',          &
     &              '        RP        dXC      dR')                    
 1803    format(i4,a8,2(1pe10.3),0pf7.3,1x,(4(1pe10.3))) 
                                                                        
! WRITE GAS-AQUEOUS SUMS                                                
!   MAKE XC,  RP, RL=GAS+AQ SUM for gas species.                        
!    CONVERT AQUEOUS XC TO GAS UNITS                                    
!     (Note, sum aqueous into gas for xc, but not for xcout.            
!      xcout is made gas-aqueous sum in postlump)                       
!                                                                       
        do j=1,c_nchem2 
          if(c_npequil(j).eq.j) then 
            xrr( kw,1) = 0. 
            rpro( kw,1) = 0. 
            rloss( kw,1) = 0. 
            do  neq=1,(c_nequil(  j)+1) 
              ic=j 
              if(neq.gt.1) ic = c_ncequil(  j,(neq-1)) 
              alpha(  kw) = 1. 
              if(neq.gt.1) alpha( 1)= c_h2oliq( 1)*avogadrl 
              rloss( kw,1) = rloss( kw,1)+c_rl( kw,ic) 
              rpro( kw,1)  = rpro( kw,1) +c_rp( kw,ic) 
! GAS-AQUEOUS SUM + CONVERSION                                          
!             xrr( kw,1)   = xrr( kw,1)  +c_xcout( kw,ic)*alpha( kw)    
! AQUEOUS CONVERSION, NO GAS-AQ SUM                                     
              xrr( kw,1)   = c_xcout( kw,ic)*alpha( kw) 
                              !do  neq=1,(c_nequil(  j)+1)              
            enddo 
                                          ! CORRECTION- XR IS GAS-AQ SUM
            xrr( kw,1) = c_xcout( kw,j) 
            beta(1)  = xrr( kw,1) -   c_xcin( kw,j) 
            gamma(1) = rpro( kw,1) - rloss( kw,1) 
                                                                        
            write(c_out,1803) j, c_tchem(j),                            &
     &       xrr( kw,1),  c_xcin( kw,j), xcfinr(kw,j)                   &
     &       ,rloss( kw,1), rpro( kw,1)                                 &
     &       ,beta(1), gamma(1)                                         
                            !if(c_npequil(j).eq.j) then                 
          endif 
                         !do j=1,c_nchem2                               
        enddo 
                                                                        
                                                                        
! ------                                                                
! END SUMMARY WRITE OPTION                                              
! ------                                                                
                                                                        
                                                                        
! GAS-PHASE REACTION RATES                                              
         write(c_out,1143) 
 1143    format(/' GAS-PHASE REACTION RATES') 
         do      nr=1,c_nreac 
          write(c_out,1144) nr, (c_treac(j,nr),j=1,5), c_rr(kw ,nr),    &
     &    ratek(kw ,nr)                                                 
 1144     format(i4,2x,a8,'+',a8,'=>',a8,'+',a8,'+',a8,2(1pe10.3)) 
         enddo 
                                                                        
          write(c_out,1146) 
 1146      format(/,'HENRYS LAW CONSTANT + MODIFIED CONSTANT',          &
     &       ' WITH DIFFUSION ADJUSTMENT:')                             
          if(c_nreach.gt.0) then 
            do i=1,c_nreach 
               write(c_out,1147) (c_treach(j,i),j=1,2),                 &
     &         rateh(kw ,i),rhdif(kw ,i)                                
            enddo 
!           write(c_out,1147) ((c_treach(j,i),j=1,2),                   
!    *         rateh(kw ,i),rhdif(kw ,i),i=1,c_nreach)                  
 1147      format(a8,'=',a8,2x,2(1pe10.3)) 
          endif 
                                                                        
! AQUEOUS EQUILIBRIUM CONSTANTS                                         
         write(c_out,1153) 
 1153    format(/,' AQUEOUS EQUILIBRIUM CONSTANTS (H+)') 
         do  nrq = 1,c_nreacq 
          write(c_out,1154) nrq, (c_treacq(j,nrq),j=1,3),               &
     &    rateq(kw ,nrq)                                                
 1154     format(i4,2x,a8,'<=>',a8,'+',a8,2(1pe10.3)) 
         enddo 
                                                                        
! SPECIAL EQUILIBRIUM CONSTANTS                                         
         write(c_out,1156) 
 1156    format(/,' SPECIAL EQUILIBRIUM CONSTANTS ') 
         do  nrqq = 1,c_nreaqq 
          write(c_out,1154) nrqq, (c_treaqq(j,nrqq),j=1,3),             &
     &    rateqq(kw ,nrqq)                                              
         enddo 
                                                                        
                                                                        
 2000  return 
      END                                           
                                                                        
                                                                        
! --------------------------------------------------------------        
        subroutine analyze(titl, kw) 
                                                                        
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
!                                                                       
! -------------------------------------------------------------------   
      IMPLICIT NONE 
                                                                        
      include 'chemmech.EXT' 
      INCLUDE 'chemvars.EXT' 
      include 'chemlocal.EXT' 
                                                                        
                                   ! Name of specified chem species     
        character*8 titl 
                                   ! Gas phase concentration            
        double precision xcgas 
                                  ! Aqueous   concentration             
        double precision acquacon 
                                   ! Dimensionless Henry coefficient    
        double precision xcoeff 
                                   ! Droplet diffusion factor           
        double precision xcoeff2 
!                                                                       
                                   ! Chem. reaction rate molec/cm3/step 
        double precision tpro 
                                   ! Species production  molec/cm3/step 
        double precision xpro 
                                   ! Stoichiometry for spec. production 
        double precision stopro 
!                                                                       
                                   ! Chem. reaction loss molec/cm3/step 
        double precision tloss 
                                   ! Species loss rate   molec/cm3/step 
        double precision xloss 
                                   ! Stoichiometry for species loss     
        double precision stoloss 
!                                                                       
                                   ! Net production minus loss /cm3/step
        double precision tnetpro 
                                   ! Change in species conc molec/cm3   
        double precision tdelta 
                                                                        
                                   ! counter for  aqueous  spec         
        integer neq1 
                              ! function to return chem species index   
        integer namechem 
                                                                        
        if( kw.le.0) return 
        kk=1 
                                                                        
!        (Replaced by parameter values in common)                       
!       avogadrl = 6.02E20                                              
!       atmos = 2.247E+19                                               
                                                                        
        ics=namechem(titl) 
        if(ics.le.0) then 
           write(c_out,21) titl, ics 
   21     format(/,'WARNING:  UNKNOWN CHEMICAL NAME IN SUBROUTINE',     &
     &    'ANALYZE',/,'NAME = ',a8,'  IC=',i5)                          
           return 
        endif 
                                                                        
!                                                                       
! LOOP FOR ALL SPECIES LINKED THROUGH CHEMICAL PAIRS  - CUT             
        ic=ics 
                                                                        
        write(c_out,22) titl,ic 
   22   format(/,'CHEMICAL PRODUCTION AND LOSS ANALYSIS FOR:',/,        &
     &  'SPECIES = ',a8,'  IC=',i5)                                     
                                                                        
        if(              c_nequil(ic).gt.0)                             &
     &     write(c_out,23) c_h2oliq( kw)                                
   23  format(/,'LIQUID WATER (gr/cm3) =', 1pe10.3) 
                                                                        
!  IF ACQUA>0, WRITE ASSOCIATED HENRY'S LAW AND EQUILIBRIA.             
        acquacon = c_h2oliq( kw)*avogadrl 
      if(c_nequil(ic).gt.0.and.c_h2oliq( kw).gt.0)                      &
     &  write(c_out,31) acquacon                                        
   31   format(/,'AQUATIC CONVERSION FACTOR:',                          &
     &  ' MOLE/LITER per MOLEC/CM3 =', 1pe10.3)                         
                                                                        
                                                                        
!  IF ACQUA>0, LOOP FOR AQUEOUS PHASE SPECIES                           
        if(c_nequil(ic).gt.0.and.c_h2oliq( kw).gt.0) then 
                                                                        
! CONCENTRATION OF GAS-SPECIES ONLY                                     
         xcgas = c_xcout( kw,ic) 
         do 40 n=1,c_nequil(ic) 
          icc=c_ncequil(ic,n) 
          xcgas = xcgas - c_xcout( kw,icc)*acquacon 
                 write(c_out,10031) icc, c_tchem(icc),c_xcout(kw,icc),  &
     &              acquacon, xcgas                                     
10031            format(' TEST ICC TCHEM XR*AQUACON subtr fr XGSUM',/,  &
     &              i5,2x,a8,2x,8(1pe10.3))                             
   40    continue 
                                                                        
         write(c_out,32) c_xcout( kw,ic), xcgas 
   32    format(/,'  TOTAL SPECIES CONCENTRATION =',1pe10.3,/,          &
     &           '    GAS SPECIES CONCENTRATION =',1pe10.3)             
                                                                        
! AQUEOUS CONCENTRATIONS AND CORRESPONDING HENRY'S LAW COEFFICIENTS     
! OR EQUILIBRIUM CONSTANTS                                              
        do 50 n=1,c_nequil(ic) 
         icc = c_ncequil(ic,n) 
         nrq = c_nrequil(ic,n) 
         nrh = c_nrequil(ic,1) 
         if(n.eq.1) then 
           xcoeff = rateh( kw,nrh)*atmos/acquacon 
           xcoeff2 = 0. 
           if(rateh(kw,nrh).gt.0)                                       &
     &       xcoeff2 = rhdif(kw,nrh)/rateh(kw,nrh)                      
           write(c_out,51) c_tchem(icc),c_xcout( kw,icc) 
   51      format(/,'AQUEOUS SPECIES = ',a8,' MOLES/LITER=',1pe10.3) 
           write(c_out,52) (c_treach(i,nrh),i=1,2)                      &
     &        ,xcoeff, xcoeff2                                          
   52      format(a8,'=',a8,'  HENRYs LAW COEF=',1pe10.3,               &
     &       '  DROPLET DIFF FAC=', 1pe10.3)                            
         else 
!          if(xcoeff.gt.0) then                                         
!            xcoeff = rateq( kw,nrq)/xcoeff                             
!          else                                                         
!            xcoeff = rateq( kw,nrq)                                    
!          endif                                                        
             xcoeff = rateq( kw,nrq) 
           write(c_out,51) c_tchem(icc),c_xcout( kw,icc) 
           write(c_out,53) (c_treacq(i,nrq),i=1,3),xcoeff 
   53      format(a8,'=',a8,'+',a8,                                     &
     &      '  EQUILIBRIUM CONSTANT= ',1pe10.3)                         
         endif 
   50   continue 
! END IF FOR AQUEOUS SPECIES                                            
      endif 
                                                                        
                                                                        
!  SUM AND WRITE CHEMICAL PRODUCTION  (RATE MOL CM-3 PER TIME STEP)     
       tpro  = 0. 
                                                                        
       write(c_out,201) 
  201  format(/,'PHOTOCHEMICAL PRODUCTION:',/,                          &
     & '         REACTION                ',                             &
     &    ' PRODUCTION  RATE        RATEK     ')                        
                                                                        
! DO FOR ALL SPECIES LINKED THROUGH CHEMICAL PAIRS - CUT.               
        ic=ics 
                                                                        
! DO FOR GAS AND ALL ACQUEOUS EQUIVALENT SPECIES                        
       do 280 neq1=1,(c_nequil(ic)+1) 
        if(neq1.eq.1) then 
          acquacon = 1. 
          icc=ic 
        else 
          acquacon = c_h2oliq( kw)*avogadrl 
          neq=neq1-1 
          icc=c_ncequil(ic,neq) 
        endif 
        if(acquacon.eq.0) go to 280 
                                                                        
! LOOP THROUGH ALL REACTIONS TO FIND LOSSES FOR SPECIES                 
! NOTE THAT RR, REACTION RATE, IS ALREADY IN GAS UNITS FOR  AQ SPECIES. 
        do 270 nr=1,c_nreac 
         stopro  = 0. 
         if(c_reactant(nr,1).eq.icc) stopro  = stopro  - 1. 
         if(c_reactant(nr,2).eq.icc) stopro  = stopro  - 1. 
         do 265 i=1,20 
          if(c_product(nr,i).eq.icc) stopro  = stopro  + c_stoich(nr,i) 
  265     continue 
         if(stopro .gt.0.) then 
           xpro  = c_rr( kw,nr)*         stopro 
           tpro  = tpro  + xpro 
           write(c_out,102) nr, (c_treac(i,nr),i=1,5),                  &
     &       xpro, c_rr( kw,nr) , ratek( kw,nr)                         
         endif 
  270   continue 
  280  continue 
                                                                        
!  SUM AND WRITE CHEMICAL LOSSES                                        
       tloss = 0. 
                                                                        
       write(c_out,101) 
  101  format(/,'PHOTOCHEMICAL LOSSES:',/,                              &
     & '         REACTION                ',                             &
     &   ' LOSS        RATE        RATEK')                              
  102  format(i4,1x,a8,'+',a8,'=>',a8,'+',a8,'+',a8,3(1pe10.3   )) 
                                                                        
                                                                        
! DO FOR ALL SPECIES LINKED THROUGH CHEMICAL PAIRS - CUT                
        ic=ics 
! DO FOR GAS AND ALL ACQUEOUS EQUIVALENT SPECIES                        
       do 180 neq1=1,(c_nequil(ic)+1) 
        if(neq1.eq.1) then 
          acquacon = 1. 
          icc=ic 
        else 
          acquacon = c_h2oliq( kw)*avogadrl 
          neq=neq1-1 
          icc=c_ncequil(ic,neq) 
        endif 
        if(acquacon.eq.0) go to 180 
                                                                        
! LOOP THROUGH ALL REACTIONS TO FIND LOSSES FOR SPECIES                 
        do 170 nr=1,c_nreac 
         stoloss = 0. 
         if(c_reactant(nr,1).eq.icc) stoloss = stoloss + 1 
         if(c_reactant(nr,2).eq.icc) stoloss = stoloss + 1 
         do 165 i=1,20 
          if(c_product(nr,i).eq.icc) stoloss = stoloss - c_stoich(nr,i) 
  165     continue 
         if(stoloss.gt.0.) then 
           xloss = c_rr( kw,nr)*         stoloss 
           tloss = tloss + xloss 
           write(c_out,102) nr,(c_treac(i,nr),i=1,5),                   &
     &      xloss, c_rr( kw,nr) , ratek( kw,nr)                         
         endif 
  170   continue 
  180  continue 
                                                                        
! WRITE FINAL SUMMARY                                                   
       tnetpro = tpro-tloss 
                                                                        
       write(c_out,301) titl, tpro, tloss, tnetpro 
  301  format(/,'SPECIES = ',a8,/,                                      &
     &          ' SUMMED PHOTOCHEMICAL PRODUCTION = ',1pe10.3,/,        &
     &          ' SUMMED PHOTOCHEMICAL LOSS       = ',1pe10.3,/,        &
     &          ' NET    PHOTOCHEMICAL PRODUCTION = ',1pe10.3)          
                                                                        
!    FOR ALL SPECIES LINKED THROUGH CHEMICAL PAIRS - CUT                
        ic=ics 
!                                                                       
       tdelta = c_xcout(kw,ic) -   c_xcin(kw,ic) 
       write(c_out,302) c_tchem(ic),                                    &
     &    c_xcout( kw,ic),  c_xcin( kw,ic), tdelta,                     &
     &          c_rp( kw,ic),c_rl( kw,ic)                               
  302  format(/,                                                        &
     & ' SPEC GAS+AQ CONC. PRIOR CONC. NET CHANGE   ',                  &
     &     '(INTERNAL RPRO     RLOSS)',/                                &
     & , a8, 3(1pe12.4),2x,2(1pe12.4))                                  
                                                                        
 2000  return 
      END                                           
