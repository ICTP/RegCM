module aerosol_thermodynamics

  USE mo_kind, ONLY : dp
    !use chamber_Parameters
    !use chamber_Global, only: n_org_spc



   ! Prototype module for calculating water content and vapour pressures for 
   ! organic species condensing in a chamber study.


   ! Prototype module for efficient calculation of the partitioning of species between
   ! gas and aerosol. Built on a newly developed inorganic compound partitioning code
   ! developed by Dave Topping.
   !
   ! Written by Dave.Topping. Pure organic component properties predicted by Mark Barley
   ! based on VOCs predicted in MCM simulations performed by Mike Jenkin. Delivered by
   ! Gordon McFiggans as Deliverable D22 from WP1.4 in the EU FP6 EUCAARI Integrated Project.
   !
   ! Queries concerning the use of this code through Gordon McFiggans, g.mcfiggans@manchester.ac.uk,
   ! Ownership: D. Topping, Centre for Atmospheric Sciences, University of Manchester, 2007

 !implicit none   
 !!glc- delcare vapour pure here, so that it can be accessed from the rates module
 !real(kind=dp), dimension(n_org_spc) :: vapour_pure
 !!glc public
 !public :: vapour_pure
 contains


   subroutine inorganic_pdfite(rh,temp,ions,water_total,Press_HNO3,Press_HCl,Press_NH3,gamma_out,mols_out)

     !the following subroutine calculates the water content of a mixed inorganic/organic
     !particle as well as equilibrium vapour pressures above the solution (HNO3,HCl,NH3 and 
     !representative organic compounds)
     !Inorganic ions treated (H-NH4-Na-SO4-NO3-Cl) - bisulphate ion concentration calculated 
     !as described in a shortly to appear publication.
     !currently, one representative secondary organic compound is treated
     !- commented sections indicate how 16 components (to follow, after completion of D30 and 
     !D34). The one representative SOA compound currently carries the abundance weighted
     !average properties from >3000 potentially condensable organic compounds predicted in
     !the base case simulation using the Master Chemical Mechanism (MCM). This module accounts for
     !no interaction between the inorganic & organic components
     !
     !--CALCULATION PROCEDURE----
     !
     !the following code determines inorganic solutes and 'equivalent' mole fractions from the 
     !ionic inputs (ions), for use in calculating the water content (water_total) and activity 
     !coefficients. after activity coefficients are calculated these are then used for calculating
     !vapour pressures using appropriate equilibrium constants (Press_HNO3,Press_NH3,Press_HCl).
     !the one representative organic is treated as ideal for contributions to the total water and for
     !calculating its equilibrium vapour pressure using raoults law.
     !note: the organic pure component vapour pressures for >3000 potentially condensable organic 
     !compounds from the base case simulation of the Master Chemical Mechanism (MCM) are calculated 
     !from boiling point predictions. absorptive partitioning of all compounds at their total predicted
     !atmospheric burden to a similar organic mass is used to calculate the relative abundance of organic
     !species. the properties of the single representative organic component in this simulation are
     !calculated as the condensed phase abundance-weighted average of all condensing components
     !
     !the following text describe the sections of code individually which are labelled 1,2,3,4,5
     !
     !1) - COMPOSITION DEFINITIONS
     !
     ! a)Inorganic ion pairing
     !
     !  in order to calculate the water content, which is also used in calculating vapour 
     !  pressures, one needs to pair the anions and cations for use in the ZSR mixing rule
     !  The equation provided by Clegg et al (2001) is used for ion pairing.
     !  The solutes chosen comprise of 9 inorganic salts and acids which provide a pairing between 
     !  each anion and cation:
     !  (NH4)2SO4,NH4NO3,NH4Cl,Na2SO4,NaNO3,NaCl,H2SO4,HNO3,HCl
     !  the organic compound is treated as a seperate solute
     ! 
     ! b)Inorganic equivalent fractions
     !  These values are calculated so that activity coefficients can be expressed by
     !  a linear additive rule, thus allowing more efficient calculations and future expansion
     !  (see more detailed description below)
     !
     !2) - WATER CALCULATION
     !
     ! a)The water content is calculated using the ZSR rule with solute concentrations calculated 
     !  using 1a above
     !  Whilst the usual approximation of ZSR relies on binary data consisting of 5th or higher order 
     !  polynomials, in this code 4 different RH regimes are used, each housing cubic equations for 
     !  the water associated with each solute listed above.
     !  Binary water contents for inorganic components were calculated using AIM online (Clegg et al 1998)
     !  The water associated with the organic compound is calculated assuming ideality and that aw=RH
     !
     ! b)molality of each inorganic ion and organic solute (initial input) is calculated for use in 
     !  vapour pressure calculation
     !
     !3) - BISULPHATE DISSOCIATION CALCULATION
     !
     ! The dissociation of the bisulphate ion is calculated explicitly
     !  A solution to the equilibrium euqation between the bisulphate ion, hydrogen ion and sulphate ion
     !  is found using tabulated equilibrium constants (referenced)
     !  It is necessary to calculate the activity coefficients of HHSO4 and H2SO4 in a noniterative 
     !  manner. these are calulated using the same format as described in 4) below, where both activity
     !  coefficients were fit to the output from ADDEM (Topping et al 2005a,b) covering an extensive 
     !  composition space, providing the activity coefficients and bisulphate ion dissociation as a 
     !  function of equivalent mole fractions and relative humidity.
     !
     !4) ACTIVITY COEFFICIENTS -for vapour pressures of HNO3,HCl and NH3
     !
     ! This section evaluates activity coefficients and vapour pressures using the water content 
     !  calculated above) for each inorganic condensing species
     ! a - HNO3
     ! b - NH3
     ! c - HCl
     !  the following procedure is used:
     !
     !  Zaveri et al (2005) found that one could express the variation of activity coefficients
     !  linearly in log space if equivalent mole fractions were used. So, by a taylor series expansion
     !  log(activity coefficient)=log(binary activity coefficient at a given RH)+
     !   (equivalent mole fraction compound A)*('interaction' parameter between A and condensing specie)+
     !   (equivalent mole fraction compound B)*('interaction' parameter between B and condensing specie)+
     !  here, the interaction parameters have been fit to ADDEM by searching the whole compositon space 
     !  and fit usign the Levenberg-Marquardt nonlinear least squares algorithm.
     !
     !  They are given as a function of RH and vary with complexity ranging from linear to 5th order 
     !  polynomial expressions, the binary activity coefficients were calculated using AIM online.
     !  note: for NH3, no binary activity coefficient was used and the data were fit to the ratio of
     !  the activity coefficients for the ammonium and hydrogen ions
     !  Once the activity coefficients are obtained the vapour pressure can be easily calculated
     !  using tabulated equilibrium constants (referenced)
     !  This procedure differs from that of Zaveri et al (2005) in that it is not assumed one can carry
     !  behaviour from binary mixtures in multicomponent systems. To this end we have fit the 
     ! 'interaction' parameters explicitly to a general inorganic equilibrium model (ADDEM - Topping et 
     !  al 2005a,b). Such parameters take into account bisulphate ion dissociation and water content.
     !  This also allows us to consider one regime for all composition space, rather than defining 
     !  sulphate rich and sulphate poor regimes
     !
     !5) - The vapour pressure of the condensing organic species are calculated in this prototype by 
     !  assuming ideality and using Raoults law. In the final module, fits to organic composition space
     !  within UNIFAC will be used to estimate component activities (including water associated with the 
     !  organic compounds) and hence the organic vapour pressures.
     !
     !-----------------------------------------------------------------------------------

     !VARIABLE DECLARATIONS

     ! this array is for the output of the activity coefficients for calculating the 
     ! non-ideal dissociation constants
     ! 1: gamma_HNO3
     ! 2: gamma_HCl
     ! 3: gamma_NH4+/gamma_H+ ("gamma_NH3")
     ! 4: (gamma_HHSO4**2)/gamma_H2SO4 (for H2SO4 dissociation)
     ! 5: (gamma_H2SO4**3)/(gamma_HHSO4**2) (for HSO4 dissociation)

     ! Juha: Added the following
     ! 6: gamma_NH4HSO2
     ! 7: gamma_HHSO4

     real(kind=dp), dimension(7) :: gamma_out  ! Juha: Changed dimension 5->7

     real(kind=dp),dimension(7)::ions,ions_mol

     real(kind=dp),dimension(7) :: mols_out  ! Juha: put out ion molalities

     real(kind=dp)::charge_sum, nitric_acid, hydrochloric_acid, sulphuric_acid,&
          &ammonium_sulphate, ammonium_nitrate, ammonium_chloride, sodium_sulphate,&
          &sodium_nitrate, sodium_chloride, solutes, hydrochloric_acid_eq_frac,&
          &sulphuric_acid_eq_frac, ammonium_sulphate_eq_frac, ammonium_nitrate_eq_frac,&
          &ammonium_chloride_eq_frac, sodium_sulphate_eq_frac, sodium_nitrate_eq_frac,&
          &sodium_chloride_eq_frac, nitric_acid_eq_frac, RH, Temp

     real(kind=dp)::hno3_w, hcl_w, h2so4_w, nh42so4_w, nh4no3_w, nh4cl_w,&
          &na2so4_w, nano3_w, nacl_w, water_total, h_out, hso4_out, so4_out


     real(kind=dp)::nh43hso42_w, nh4hso4_w, na3hso42_w, nahso4_w

     real(kind=dp)::hno3_water, hcl_water, h2so4_water, nh42so4_water, nh4no3_water, nh4cl_water,&
          &na2so4_water, nano3_water, nacl_water, nh43hso42_water, nh4hso4_water, na3hso42_water, nahso4_water

     real(kind=dp)::nNa, nNH4, nSulf, Xt, fNa, fNH4


     integer::binary_case, full_complexity

     real(kind=dp)::binary_hno3

     real(kind=dp)::HCL_hno3, H2SO4_hno3, NH42SO4_hno3, NH4NO3_hno3, NH4Cl_hno3,&
          &Na2SO4_hno3, NaNO3_hno3, NaCl_hno3

     real(kind=dp)::Ln_HNO3_act, gamma_hno3, Press_HNO3, K_hno3

     real(kind=dp)::binary_hcl

     real(kind=dp)::HNO3_hcl, H2SO4_hcl, NH42SO4_hcl, NH4NO3_hcl, NH4Cl_hcl,&
          &Na2SO4_hcl, NaNO3_hcl, NaCl_hcl

     real(kind=dp)::Ln_HCl_act, gamma_hcl, Press_HCl

     real(kind=dp)::HNO3_nh3, HCl_nh3, H2SO4_nh3, NH42SO4_nh3, NH4NO3_nh3, NH4Cl_nh3,&
          &Na2SO4_nh3, NaNO3_nh3, NaCl_nh3

     real(kind=dp)::Ln_NH3_act, gamma_nh3, Press_NH3

     real(kind=dp)::K_hcl

     real(kind=dp)::Kh, Knh4, Kw, molality_ratio_nh3

     !----variables used in calculating bisulphate dissociation----

     real(kind=dp)::binary_hhso4

     real(kind=dp)::HNO3_hhso4, HCL_hhso4, NH42SO4_hhso4, NH4NO3_hhso4, NH4Cl_hhso4,&
          &Na2SO4_hhso4, NaNO3_hhso4, NaCl_hhso4

     real(kind=dp)::Ln_hhso4_act, gamma_hhso4

     real(kind=dp)::binary_h2so4

     real(kind=dp)::HNO3_h2so4, HCl_h2so4, NH42SO4_h2so4, NH4NO3_h2so4, NH4Cl_h2so4,&
          &Na2SO4_h2so4, NaNO3_h2so4, NaCl_h2so4

     real(kind=dp)::Ln_h2so4_act, gamma_h2so4

     real(kind=dp)::act_product, a, b, c, root1, root2

     real(kind=dp)::h_real, hso4_real, so4_real

     real(kind=dp)::Mol_frac_hno3,Mol_frac_hcl,Mol_frac_h2so4,Mol_frac_nh42so4,Mol_frac_nh4no3,&
          &Mol_frac_nh4cl,Mol_frac_na2so4,Mol_frac_nano3,Mol_frac_nacl


     real(kind=dp)::Mol_frac_nh43hso42, Mol_frac_nh4hso4, Mol_frac_na3hso42, Mol_frac_nahso4


     real(kind=dp):: henrys_temp_dep

     !---new nh3 variables---
     real(kind=dp)::binary_nh4hso4,HNO3_nh4hso4,HCL_nh4hso4,H2SO4_nh4hso4,NH42SO4_nh4hso4,&
          &NH4NO3_nh4hso4,NH4Cl_nh4hso4,Na2SO4_nh4hso4,NaNO3_nh4hso4,NaCl_nh4hso4
     real(kind=dp)::Ln_NH4HSO4_act,gamma_nh4hso4

     !----------------------------------------------------------------------------------
     !
     !VALUE INITIALISATION
     henrys_temp_dep = (1./temp - 1./298.d0)

     HCL_hno3=1.0d0;H2SO4_hno3=1.0d0;NH42SO4_hno3=1.0d0;NH4NO3_hno3=1.0d0;NH4Cl_hno3=1.0d0;
     Na2SO4_hno3=1.0d0;NaNO3_hno3=1.0d0;NaCl_hno3=1.0d0;
     HNO3_hcl=1.0d0;H2SO4_hcl=1.0d0;NH42SO4_hcl=1.0d0;NH4NO3_hcl=1.0d0;NH4Cl_hcl=1.0d0;
     Na2SO4_hcl=1.0d0;NaNO3_hcl=1.0d0;NaCl_hcl=1.0d0;
     HNO3_nh3=1.0d0;HCl_nh3=1.0d0;H2SO4_nh3=1.0d0;NH42SO4_nh3=1.0d0;NH4NO3_nh3=1.0d0;
     NH4Cl_nh3=1.0d0;Na2SO4_nh3=1.0d0;NaNO3_nh3=1.0d0;NaCl_nh3=1.0d0;
     HNO3_hhso4=1.0d0;HCL_hhso4=1.0d0;NH42SO4_hhso4=1.0d0;NH4NO3_hhso4=1.0d0;NH4Cl_hhso4=1.0d0;
     Na2SO4_hhso4=1.0d0;NaNO3_hhso4=1.0d0;NaCl_hhso4=1.0d0
     HNO3_h2so4=1.0d0;HCl_h2so4=1.0d0;NH42SO4_h2so4=1.0d0;NH4NO3_h2so4=1.0d0;NH4Cl_h2so4=1.0d0;
     Na2SO4_h2so4=1.0d0;NaNO3_h2so4=1.0d0;NaCl_h2so4=1.0d0;
     !--new nh3 variables--
     HNO3_nh4hso4=1.0d0;HCL_nh4hso4=1.0d0;H2SO4_nh4hso4=1.0d0;NH42SO4_nh4hso4=1.0d0;
     NH4NO3_nh4hso4=1.0d0;NH4Cl_nh4hso4=1.0d0;Na2SO4_nh4hso4=1.0d0;NaNO3_nh4hso4=1.0d0;NaCl_nh4hso4=1.0d0;


     ! Juha: added
     mols_out = 0.d0

     Press_HNO3=0.0d0
     Press_HCl=0.0d0
     Press_NH3=0.0d0       !Initialising vapour pressure over the multicomponent
     !particle

     gamma_out = 1d0   ! i.e. don't alter the ideal mixing ratios if there's nothing there.

     !------------------------------------------------------------------------
     !
     !1)-COMPOSITION DEFINITIONS
     !
     ! a)Inorganic ion pairing 
     !  pair cations and anions into solutes according to Clegg et al (2001)

     charge_sum=ions(1)+ions(2)+ions(3)+2.0d0*ions(4)+ions(5)+ions(6)+ions(7)
     nitric_acid=0.0d0;hydrochloric_acid=0.0d0;sulphuric_acid=0.0d0
     ammonium_sulphate=0.0d0;ammonium_nitrate=0.0d0;ammonium_chloride=0.0d0
     sodium_sulphate=0.0d0;sodium_nitrate=0.0d0;sodium_chloride=0.0d0
     nitric_acid=(2.0d0*ions(1)*ions(6)*((1.0d0/1.0d0)**0.5))/(charge_sum)
     hydrochloric_acid=(2.0d0*ions(1)*ions(7)*((1.0d0/1.0d0)**0.5))/(charge_sum)
     sulphuric_acid=(2.0d0*ions(1)*ions(4)*((2.0d0/2.0d0)**0.5))/(charge_sum)
     ammonium_sulphate=(2.0d0*ions(2)*ions(4)*((2.0d0/2.0d0)**0.5))/(charge_sum)  
     ammonium_nitrate=(2.0d0*ions(2)*ions(6)*((1.0d0/1.0d0)**0.5))/(charge_sum)   
     ammonium_chloride=(2.0d0*ions(2)*ions(7)*((1.0d0/1.0d0)**0.5))/(charge_sum) 
     sodium_sulphate=(2.0d0*ions(3)*ions(4)*((2.0d0/2.0d0)**0.5))/(charge_sum) 
     sodium_nitrate=(2.0d0*ions(3)*ions(6)*((1.0d0/1.0d0)**0.5))/(charge_sum)  
     sodium_chloride=(2.0d0*ions(3)*ions(7)*((1.0d0/1.0d0)**0.5))/(charge_sum)

     ! b) - Inorganic equivalent fractions

     solutes=0.0d0
     solutes=3.0d0*sulphuric_acid+2.0d0*hydrochloric_acid+2.0d0*nitric_acid+3.0d0*ammonium_sulphate+&
          &2.0d0*ammonium_nitrate+2.0d0*ammonium_chloride+3.0d0*sodium_sulphate+2.0d0*sodium_nitrate+&
          &2.0d0*sodium_chloride

     nitric_acid_eq_frac=2.0d0*nitric_acid/(solutes)
     hydrochloric_acid_eq_frac=2.0d0*hydrochloric_acid/(solutes)
     sulphuric_acid_eq_frac=3.0d0*sulphuric_acid/(solutes)
     ammonium_sulphate_eq_frac=3.0d0*ammonium_sulphate/(solutes)
     ammonium_nitrate_eq_frac=2.0d0*ammonium_nitrate/(solutes)
     ammonium_chloride_eq_frac=2.0d0*ammonium_chloride/(solutes)
     sodium_sulphate_eq_frac=3.0d0*sodium_sulphate/(solutes)
     sodium_nitrate_eq_frac=2.0d0*sodium_nitrate/(solutes)
     sodium_chloride_eq_frac=2.0d0*sodium_chloride/(solutes)


     !--inorganic ion molalities

     ions_mol(:)=0.0d0
     ions_mol(1)=ions(1)/(water_total*18.01528d-3)  !H
     ions_mol(2)=ions(2)/(water_total*18.01528d-3)  !NH4
     ions_mol(3)=ions(3)/(water_total*18.01528d-3)  !Na
     ions_mol(4)=ions(4)/(water_total*18.01528d-3)  !SO4
     ions_mol(5)=ions(5)/(water_total*18.01528d-3)  !HSO4
     ions_mol(6)=ions(6)/(water_total*18.01528d-3)  !NO3
     ions_mol(7)=ions(7)/(water_total*18.01528d-3)  !Cl

     !***
     !At this point we may need to introduce a method for prescribing H+ when
     !there is no 'real' value for H+..i.e. in the sulphate poor domain
     !this will give a value for
     !-solve quadratic proposed by Zaveri et al 2005




     !-----------------------------------------------------------------------------------------
     !
     !3) BISULPHATE ION DISSOCIATION
     !
     !Note: the flags "binary_case" and "full_complexity" are not used in this prototype
     !They are used for simplification of the fit expressions when using limited composition regions
     !
     !--this section of code calculates the bisulphate ion concentration

     if (ions(1) .GT. 0.0d0 .AND. ions(4) .GT. 0.0d0) then

        !----------HHSO4-------------
        binary_case=1;

        if (RH .GT. 0.1 .AND. RH .LT. 0.9) then
           binary_hhso4=- 4.9521*(RH**3)+9.2881*(RH**2)-10.777*RH+6.0534;
        elseif (RH .GE. 0.9 .AND. RH .LT. 0.955) then
           binary_hhso4=-6.3777*RH+5.962;
        elseif (RH .GE. 0.955 .AND. RH .LT. 0.99) then
           binary_hhso4=2367.2*(RH**3)-6849.7*(RH**2)+6600.9*RH-2118.7;
        elseif (RH .GE. 0.99 .AND. RH .LT. 0.9999) then
           binary_hhso4=3E-07*(RH**5)-2E-05*(RH**4)+0.0004*(RH**3)-0.0035*(RH**2)+0.0123*RH-0.3025;
        endif

        if (nitric_acid .gt. 0.0d0) then
           HNO3_hhso4=-4.2204*(RH**4)+12.193*(RH**3)-12.481*(RH**2)+6.459*RH-1.9004;
        endif

        if (hydrochloric_acid .gt. 0.0d0) then 
           HCL_hhso4=- 54.845*(RH**7)+209.54*(RH**6)-336.59*(RH**5)+294.21*(RH**4)-150.07*(RH**3)+43.767*(RH**2)-6.5495*RH+0.60048;
        endif

        if (ammonium_sulphate .gt. 0.0d0) then
           NH42SO4_hhso4=16.768*(RH**3)-28.75*(RH**2)+20.011*RH-8.3206;
        endif

        if (ammonium_nitrate .gt. 0.0d0) then
           NH4NO3_hhso4=-17.184*(RH**4)+56.834*(RH**3)-65.765*(RH**2)+35.321*RH-9.252;
        endif

        if (ammonium_chloride .gt. 0.0d0) then
           if (RH .LT. 0.2 .AND. RH .GE. 0.1) then
              NH4Cl_hhso4=3.2809*RH-2.0637;
           elseif (RH .GE. 0.2 .AND. RH .LT. 0.99) then
              NH4Cl_hhso4=-1.2981*(RH**3)+4.7461*(RH**2)-2.3269*RH-1.1259;
           endif
        endif

        if (sodium_sulphate .gt. 0.0d0) then
           Na2SO4_hhso4=118.87*(RH**6)-358.63*(RH**5)+435.85*(RH**4)-272.88*(RH**3)+94.411*(RH**2)-18.21*RH+0.45935;
        endif

        if (sodium_nitrate .gt. 0.0d0) then
           if (RH .LT. 0.2 .AND. RH .GE. 0.1) then
              NaNO3_hhso4=4.8456*RH-2.5773;    
           elseif (RH .GE. 0.2 .AND. RH .LT. 0.99) then
              NaNO3_hhso4=0.5964*(RH**3)-0.38967*(RH**2)+1.7918*RH-1.9691; 
           endif
        endif

        if (sodium_chloride .gt. 0.0d0) then
           if (RH .LT. 0.2) then
              NaCl_hhso4=0.51995*RH-1.3981;
           elseif (RH .GE. 0.2 .AND. RH .LT. 0.99) then
              NaCl_hhso4=1.6539*RH-1.6101;
           endif
        endif

        Ln_hhso4_act=binary_hhso4+nitric_acid_eq_frac*HNO3_hhso4+hydrochloric_acid_eq_frac*HCL_hhso4+&
             &ammonium_sulphate_eq_frac*NH42SO4_hhso4+ammonium_nitrate_eq_frac*NH4NO3_hhso4+ammonium_chloride_eq_frac*NH4Cl_hhso4+&
             &sodium_sulphate_eq_frac*Na2SO4_hhso4+sodium_nitrate_eq_frac*NaNO3_hhso4+sodium_chloride_eq_frac*NaCl_hhso4

        gamma_hhso4=exp(Ln_hhso4_act)


        !----------H2SO4--------------

        if (RH .GE. 0.1 .AND. RH .LT. 0.9) then
           binary_h2so4=2.4493*(RH**2)-6.2326*RH+2.1763;
        elseif (RH .GE. 0.9 .AND. RH .LT. 0.98) then
           binary_h2so4=914.68*(RH**3)-2502.3*(RH**2)+2281.9*RH-695.11;
        elseif (RH .GE. 0.98 .AND. RH .LT. 0.9999) then
           binary_h2so4=3E-08*(RH**4)-5E-06*(RH**3)+0.0003*(RH**2)-0.0022*RH-1.1305;
        endif

        if (nitric_acid .gt. 0.0d0) then
           HNO3_h2so4=-16.382*(RH**5)+46.677*(RH**4)-54.149*(RH**3)+34.36*(RH**2)-12.54*RH+2.1368;
        endif

        if (hydrochloric_acid .gt. 0.0d0) then 
           HCl_h2so4=-14.409*(RH**5)+42.804*(RH**4)-47.24*(RH**3)+24.668*(RH**2)-5.8015*RH+0.084627;
        endif

        if (ammonium_sulphate .gt. 0.0d0) then
           NH42SO4_h2so4=66.71*(RH**5)-187.5*(RH**4)+210.57*(RH**3)-121.04*(RH**2)+39.182*RH-8.0606;
        endif

        if (ammonium_nitrate .gt. 0.0d0) then
           NH4NO3_h2so4=-22.532*(RH**4)+66.615*(RH**3)-74.647*(RH**2)+37.638*RH-6.9711;
        endif

        if (ammonium_chloride .gt. 0.0d0) then
           if (RH .GE. 0.1 .AND. RH .LT. 0.2) then
              NH4Cl_h2so4=-0.32089*RH+0.57738;
           elseif (RH .GE. 0.2 .AND. RH .LT. 0.9) then
              NH4Cl_h2so4=18.089*(RH**5)-51.083*(RH**4)+50.32*(RH**3)-17.012*(RH**2)-0.93435*RH+1.0548;
           elseif (RH .GE. 0.9 .AND. RH .LT. 0.99) then
              NH4Cl_h2so4=-1.5749*RH+1.7002;
           endif
        endif

        if (sodium_sulphate .gt. 0.0d0) then
           Na2SO4_h2so4=29.843*(RH**4)-69.417*(RH**3)+61.507*(RH**2)-29.874*RH+7.7556;
        endif

        if (sodium_nitrate .gt. 0.0d0) then
           NaNO3_h2so4=-122.37*(RH**6)+427.43*(RH**5)-604.68*(RH**4)+443.08*(RH**3)-178.61*(RH**2)+37.242*RH-1.9564;
        endif

        if (sodium_chloride .gt. 0.0d0) then
           NaCl_h2so4=-40.288*(RH**5)+115.61*(RH**4)-129.99*(RH**3)+72.652*(RH**2)-22.124*RH+4.2676;
        endif

        Ln_h2so4_act=binary_h2so4+nitric_acid_eq_frac*HNO3_h2so4+hydrochloric_acid_eq_frac*HCl_h2so4+&
             &ammonium_sulphate_eq_frac*NH42SO4_h2so4+ammonium_nitrate_eq_frac*NH4NO3_h2so4+ammonium_chloride_eq_frac*NH4Cl_h2so4+&
             &sodium_sulphate_eq_frac*Na2SO4_h2so4+sodium_nitrate_eq_frac*NaNO3_h2so4+sodium_chloride_eq_frac*NaCl_h2so4;

        gamma_h2so4=exp(Ln_h2so4_act) 

        ! export activity coefficients
        if(gamma_H2SO4.gt.0d0) gamma_out(4) = (gamma_HHSO4**2d0)/gamma_H2SO4
        if(gamma_HHSO4.gt.0d0) gamma_out(5) = (gamma_H2SO4**3d0)/(gamma_HHSO4**2d0)


        !--solve the quadratic

        act_product=(gamma_h2so4**3d0)/(gamma_hhso4**2d0)

        a=1;
        b=-1*(ions(4)+ions(1)+((water_total*18.0d-3)/(99*act_product)));
        c=ions(4)*ions(1);
        root1=((-1*b)+(((b**2)-4*a*c)**(0.5)))/(2*a)
        root2=((-1*b)-(((b**2)-4*a*c)**(0.5)))/(2*a)


        if (root1 .GT. ions(1) .OR. root1 .LT. 0.0d0) then
           root1=0.0d0
        endif

        if (root2 .GT. ions(1) .OR. root2 .LT. 0.0d0) then
           root2=0.0d0
        endif


        !--calculate your new hydrogen ion, bisulphate ion and sulphate ion concentration

        hso4_real=0
        h_real=ions(1)
        so4_real=ions(4)
        if (root1 .EQ. 0.0d0) then
           hso4_real=root2
        elseif (root2 .EQ. 0.0d0) then
           hso4_real=root1
        endif
        h_real=ions(1)-hso4_real
        so4_real=ions(4)-hso4_real

        !--recalculate ion molalities
        ions_mol(1)=h_real/(water_total*18.01528d-3) !sulphuric_acid*2+hydrochloric_acid+nitric_acid  !H
        ions_mol(4)=so4_real/(water_total*18.01528d-3)  !sulphuric_acid+ammonium_sulphate+sodium_sulphate  !SO4
        ions_mol(5)=hso4_real/(water_total*18.01528d-3)  !0.0d0  !HSO4

        h_out=h_real
        hso4_out=hso4_real
        so4_out=so4_real

     elseif (ions(1) .EQ. 0.0d0 .OR. ions(4) .EQ. 0.0d0) then

        h_out=ions(1)
        hso4_out=0.0d0
        so4_out=ions(4)

     endif




     !-----------------------------------------------------------------------------------------
     !
     !4) ACTIVITY COEFFICIENTS -for vapour pressures of HNO3,HCl and NH3
     !
     !Note: the flags "binary_case" and "full_complexity" are not used in this prototype
     !They are used for simplification of the fit expressions when using limited composition regions
     !
     !a) - ACTIVITY COEFF/VAPOUR PRESSURE - HNO3
     !
     !-----------binary hno3 act coeff---------------
     if (ions(1) .gt. 0 .and. ions(6) .gt. 0) then

        binary_case=1;

        if (RH .gt. 0.1 .and. RH .lt. 0.98) then
           !10 to 98! RH full run
           if (binary_case.eq.1) then
              binary_hno3=1.8514*(RH**3)-4.6991*(RH**2)+1.5514*(RH)+0.90236
              !40 to 98 full
           elseif (binary_case.eq.2) then
              binary_hno3=-1.1751*(RH**2)-0.53794*(RH)+1.2808
           endif
        elseif (RH .ge. 0.98 .and. RH .lt. 0.9999) then
           !98 to 100!
           binary_hno3=1244.69635941351*(RH**3)-2613.93941099991*(RH**2)+1525.0684974546*(RH)-155.946764059316;
        endif

        !------Contributions from other solutes---------

        full_complexity=1

        if (hydrochloric_acid .gt. 0.0d0) then
           !---HCL---
           if (full_complexity.eq.1 .or. RH .le. 0.4) then
              HCL_hno3=16.051*(RH**4)-44.357*(RH**3)+45.141*(RH**2)-21.638*(RH)+ 4.8182
              !40! +
           elseif (full_complexity.eq.0 .and. RH .gt. 0.4) then
              HCL_hno3=-1.5833*(RH)+1.5569
              !y = 1.0392*x**{2} - 3.071*x + 2.0431
           endif
        endif

        if (sulphuric_acid .gt. 0.0d0) then
           !---H2SO4----
           if (full_complexity.eq.1 .or. RH .le. 0.4) then
              H2SO4_hno3=-3.0849*(RH**3)+5.9609*(RH**2)-4.468*(RH)+1.5658;
           elseif (full_complexity.eq.0 .and. RH .gt. 0.4) then
              !40! +
              H2SO4_hno3=-0.93473*RH+0.9363
              !y = - 0.33243*x**{2} - 0.45881*x + 0.78079
           endif
        endif

        if (ammonium_sulphate .gt. 0.0d0) then
           !---NH42SO4----
           NH42SO4_hno3=16.821*(RH**3)-28.391*(RH**2)+18.133*(RH)-6.7356
        endif

        if (ammonium_nitrate .gt. 0.0d0) then
           !---NH4NO3-----
           NH4NO3_hno3=11.01*(RH**3)-21.578*(RH**2)+14.808*RH-4.2593
        endif

        if (ammonium_chloride .gt. 0.0d0) then
           !---NH4Cl------
           if (full_complexity.eq.1 .or. RH .le. 0.4) then
              NH4Cl_hno3=-1.176*(RH**3)+5.0828*(RH**2)-3.8792*(RH)-0.05518
           elseif (full_complexity.eq.0 .and. RH .gt. 0.4) then
              !40! +
              NH4Cl_hno3=2.6219*(RH**2)-2.2609*RH-0.38436
           endif
        endif

        if (sodium_sulphate .gt. 0.0d0) then
           !---Na2SO4-----
           Na2SO4_hno3=35.504*(RH**4)-80.101*(RH**3)+67.326*(RH**2)-28.461*(RH)+5.6016
        endif

        if (sodium_nitrate .gt. 0.0d0) then
           !---NaNO3------
           if (full_complexity.eq.1 .or. RH .le. 0.4) then
              NaNO3_hno3=23.659*(RH**5)-66.917*(RH**4)+74.686*(RH**3)-40.795*(RH**2)+10.831*(RH)-1.4701
           elseif (full_complexity.eq.0 .and. RH .gt. 0.4) then
              !40!+
              NaNO3_hno3=14.749*(RH**4)-35.237*(RH**3)+31.196*(RH**2)-12.076*(RH)+1.3605
           endif
        endif

        if (sodium_chloride .gt. 0.0d0) then
           !---NaCl------
           if (full_complexity.eq.1 .or. RH .le. 0.4) then
              NaCl_hno3=13.682*(RH**4)-35.122*(RH**3)+33.397*(RH**2)-14.586*(RH)+2.6276
           elseif (full_complexity.eq.0 .and. RH .gt. 0.4) then
              !40!
              NaCl_hno3=1.1882*(RH**3)-1.1037*(RH**2)-0.7642*(RH)+0.6671
           endif
        endif

        Ln_HNO3_act=binary_hno3+hydrochloric_acid_eq_frac*HCL_hno3+sulphuric_acid_eq_frac*H2SO4_hno3+&
             &ammonium_sulphate_eq_frac*NH42SO4_hno3+ammonium_nitrate_eq_frac*NH4NO3_hno3+ammonium_chloride_eq_frac*NH4Cl_hno3+&
             &sodium_sulphate_eq_frac*Na2SO4_hno3+sodium_nitrate_eq_frac*NaNO3_hno3+sodium_chloride_eq_frac*NaCl_hno3

        gamma_hno3=exp(Ln_HNO3_act)
        gamma_out(1) = gamma_hno3


        !---PARTIAL PRESSURE CALCULATION
        !K_hno3=2.51*(10**6);  
        !K_hno3=2.628145923d6  !calculated using AIM online (Clegg et al 1998)
        K_hno3 = 2.6d6*exp(8700d0*henrys_temp_dep) !after Chameides (1984) (and NIST database)
        Press_HNO3=(ions_mol(1)*ions_mol(6)*(gamma_hno3**2))/K_hno3

     endif

     !b) - ACTIVITY COEFF/VAPOUR PRESSURE - NH3

     if (ions(2) .gt. 0 .and. ions(1) .gt. 0) then
        !follow the two solute approach of zaveri


        !A)-------------------------NH4HSO4-----------------------------
        binary_nh4hso4=56.907*(RH**6)-155.32*(RH**5)+142.94*(RH**4)-32.298*(RH**3)-27.936*(RH**2)+19.502*(RH)-4.2618
        !---------binary---------

        if (nitric_acid .gt. 0.0d0) then
           !---HNO3---
           HNO3_nh4hso4=104.8369*(RH**8)-288.8923*(RH**7)+129.3445*(RH**6)+373.0471*(RH**5)-571.0385*(RH**4)+ &
                326.3528*(RH**3)-74.169*(RH**2)-2.4999*(RH)+3.17
        endif

        if (hydrochloric_acid .gt. 0.0d0) then
           !---HCl---
           HCL_nh4hso4=-7.9133*(RH**8)+126.6648*(RH**7)-460.7425*(RH**6)+731.606*(RH**5)-582.7467*(RH**4)+ & 
                216.7197*(RH**3)-11.3934*(RH**2)-17.7728*(RH)+5.75
        endif

        if (sulphuric_acid .gt. 0.0d0) then
           !---H2SO4---
           H2SO4_nh4hso4=195.981*(RH**8)-779.2067*(RH**7)+1226.3647*(RH**6)-964.0261*(RH**5)+391.7911*(RH**4)- &
                84.1409*(RH**3)+20.0602*(RH**2)-10.2663*(RH)+3.5817
        endif

        if (ammonium_sulphate .gt. 0.0d0) then
           !---NH42SO4---
           NH42SO4_nh4hso4=617.777*(RH**8)-2547.427*(RH**7)+4361.6009*(RH**6)-4003.162*(RH**5)+ &
                2117.8281*(RH**4)-640.0678*(RH**3)+98.0902*(RH**2)-2.2615*(RH)-2.3811
        endif

        if (ammonium_nitrate .gt. 0.0d0) then
           !---NH4NO3---
           NH4NO3_nh4hso4=-104.4504*(RH**8)+539.5921*(RH**7)-1157.0498*(RH**6)+1322.4507*(RH**5)- &
                852.2475*(RH**4)+298.3734*(RH**3)-47.0309*(RH**2)+1.297*(RH)-0.8029
        endif

        if (ammonium_chloride .gt. 0.0d0) then
           !---NH4Cl---
           NH4Cl_nh4hso4=258.1792*(RH**8)-1019.3777*(RH**7)+1592.8918*(RH**6)-1221.0726*(RH**5)+ &
                442.2548*(RH**4)-43.6278*(RH**3)-7.5282*(RH**2)-3.8459*(RH)+2.2728
        endif

        if (sodium_sulphate .gt. 0.0d0) then
           !---Na2SO4---
           Na2SO4_nh4hso4=225.4238*(RH**8)-732.4113*(RH**7)+843.7291*(RH**6)-322.7328*(RH**5)- &
                88.6252*(RH**4)+72.4434*(RH**3)+22.9252*(RH**2)-25.3954*(RH)+4.6971
        endif

        if (sodium_nitrate .gt. 0.0d0) then
           !---NaNO3---
           NaNO3_nh4hso4=96.1348*(RH**8)-341.6738*(RH**7)+406.5314*(RH**6)-98.5777*(RH**5)- &
                172.8286*(RH**4)+149.3151*(RH**3)-38.9998*(RH**2)-0.2251*(RH)+0.4953
        endif

        if (sodium_chloride .gt. 0.0d0) then
           !---NaCl---
           NaCl_nh4hso4=91.7856*(RH**8)-316.6773*(RH**7)+358.2703*(RH**6)-68.9142*(RH**5)- &
                156.5031*(RH**4)+116.9592*(RH**3)-22.5271*(RH**2)-3.7716*(RH)+1.56
        endif

        Ln_NH4HSO4_act=binary_nh4hso4+nitric_acid_eq_frac*HNO3_nh4hso4+hydrochloric_acid_eq_frac*HCL_nh4hso4+ &
             sulphuric_acid_eq_frac*H2SO4_nh4hso4+ammonium_sulphate_eq_frac*NH42SO4_nh4hso4+&
             ammonium_nitrate_eq_frac*NH4NO3_nh4hso4+ammonium_chloride_eq_frac*NH4Cl_nh4hso4+&
             sodium_sulphate_eq_frac*Na2SO4_nh4hso4+sodium_nitrate_eq_frac*NaNO3_nh4hso4+&
             sodium_chloride_eq_frac*NaCl_nh4hso4

        gamma_nh4hso4=exp(Ln_NH4HSO4_act)
        ! Juha: Put this out as well
        gamma_out(6) = gamma_nh4hso4
        ! --

        !  -------------------------------------------------------------
        !B)-------------------------HHSO4-------------------------------
        binary_hhso4=-4.9521*(RH**3)+9.2881*(RH**2)-10.777*RH+6.0534;
        !---------binary---------

        if (nitric_acid .gt. 0.0d0) then
           !---HNO3---
           HNO3_hhso4=-4.2204*(RH**4)+12.193*(RH**3)-12.481*(RH**2)+6.459*RH-1.9004;
        endif

        if (hydrochloric_acid .gt. 0.0d0) then
           !---HCl---
           HCL_hhso4=- 54.845*(RH**7)+209.54*(RH**6)-336.59*(RH**5)+294.21*(RH**4)-150.07*(RH**3)+43.767*(RH**2)-6.5495*RH+0.60048;
        endif

        if (ammonium_sulphate .gt. 0.0d0) then
           !---NH42SO4---
           NH42SO4_hhso4=16.768*(RH**3)-28.75*(RH**2)+20.011*RH-8.3206;
        endif

        if (ammonium_nitrate .gt. 0.0d0) then
           !---NH4NO3---
           NH4NO3_hhso4=-17.184*(RH**4)+56.834*(RH**3)-65.765*(RH**2)+35.321*RH-9.252;
        endif

        if (ammonium_chloride .gt. 0.0d0) then
           !---NH4Cl---
           if (RH .lt. 0.2 .and. RH .ge. 0.1) then
              NH4Cl_hhso4=3.2809*RH-2.0637;
           elseif (RH .ge. 0.2 .and. RH .lt. 0.99) then 
              NH4Cl_hhso4=-1.2981*(RH**3)+4.7461*(RH**2)-2.3269*RH-1.1259;
           endif
        endif

        if (sodium_sulphate .gt. 0.0d0) then
           !---Na2SO4---
           Na2SO4_hhso4=118.87*(RH**6)-358.63*(RH**5)+435.85*(RH**4)-272.88*(RH**3)+94.411*(RH**2)-18.21*RH+0.45935;
        endif

        if (sodium_nitrate .gt. 0.0d0) then
           !---NaNO3---
           if (RH .lt. 0.2 .and. RH .ge. 0.1) then
              NaNO3_hhso4=4.8456*RH-2.5773;    
           elseif (RH .ge. 0.2 .and. RH .lt. 0.99) then 
              NaNO3_hhso4=0.5964*(RH**3)-0.38967*(RH**2)+1.7918*RH-1.9691; 
           endif
        endif

        if (sodium_chloride .gt. 0.0d0) then
           !---NaCl---
           if (RH .lt. 0.2) then
              NaCl_hhso4=0.51995*RH-1.3981;
           elseif (RH .ge. 0.2 .and. RH .lt. 0.99) then
              NaCl_hhso4=1.6539*RH-1.6101;
           endif
        endif

        Ln_HHSO4_act=binary_hhso4+nitric_acid_eq_frac*HNO3_hhso4+hydrochloric_acid_eq_frac*HCL_hhso4+&
             &ammonium_sulphate_eq_frac*NH42SO4_hhso4+ammonium_nitrate_eq_frac*NH4NO3_hhso4+ammonium_chloride_eq_frac*NH4Cl_hhso4+&
             &sodium_sulphate_eq_frac*Na2SO4_hhso4+sodium_nitrate_eq_frac*NaNO3_hhso4+sodium_chloride_eq_frac*NaCl_hhso4

        gamma_hhso4=exp(Ln_HHSO4_act)
        ! Juha: Put this out as well
        gamma_out(7) = gamma_hhso4
        ! --
        !  -------------------------------------------------------------

        gamma_nh3=(gamma_nh4hso4**2)/(gamma_hhso4**2)  !exp(Ln_NH3_act)
        gamma_out(3) = gamma_nh3


        !this actually represents the ratio of the ammonium to hydrogen ion activity coefficients 
        !(see Zaveri paper) - multiply this by the ratio of the ammonium to hydrogen ion molality 
        !and the ratio of appropriate equilibrium constants

        !equilibrium constants
        !Kh=57.64d0  !Zaveri et al (2005)
        Kh = 5.8d1*exp(4085d0*henrys_temp_dep) ! after Chameides (1984) (and NIST database)
        !Knh4=1.81d-5  !Zaveri et al (2005)
        Knh4 = 1.7d-5*exp(-4325*henrys_temp_dep) ! after Chameides (1984)
        !Kw=1.01d-14   !Zaveri et al (2005)
        Kw = 1d-14*exp(-6716*henrys_temp_dep) ! after Chameides (1984)

        molality_ratio_nh3=ions_mol(2)/ions_mol(1)
        Press_NH3=molality_ratio_nh3*gamma_nh3*(Kw/(Kh*Knh4))

     endif

     !c) - ACTIVITY COEFF/VAPOUR PRESSURE - HCl

     if (ions(1) .gt. 0 .and. ions(7) .gt. 0) then
        !-----------binary hcl act coeff---------------
        binary_case=1;
        if (RH .gt. 0.1 .and. RH .lt. 0.98) then
           !10 to 98! RH full run
           if (binary_case.eq.1) then
              !binary = 1.7526*(RH**2) - 6.9942*(RH) + 5.0046;
              binary_hcl=- 5.0179*(RH**3)+9.8816*(RH**2)-10.789*(RH)+5.4737
              !40 to 98 full
           elseif (binary_case.eq.2) then
              binary_hcl=- 4.6221*(RH)+4.2633
           endif
        elseif (RH .ge. 0.98 .and. RH .lt. 0.9999) then
           !98 to 100!
           binary_hcl=775.6111008626*(RH**3)-2146.01320888771*(RH**2)+1969.01979670259*(RH)-598.878230033926
        endif

        if (nitric_acid .gt. 0.0d0) then
           !---HNO3---
           if (full_complexity.eq.1 .or. RH .le. 0.4) then
              HNO3_hcl=9.6256*(RH**4)-26.507*(RH**3)+27.622*(RH**2)-12.958*(RH)+2.2193
           elseif (full_complexity.eq.0 .and. RH .gt. 0.4) then
              HNO3_hcl=1.3242*(RH**2)-1.8827*(RH)+0.55706
           endif
        endif

        if (sulphuric_acid .gt. 0.0d0) then
           !--h2so4--
           if (full_complexity.eq.1 .or. RH .le. 0.4) then
              H2SO4_hcl=1.4406*(RH**3)-2.7132*(RH**2)+1.014*RH+0.25226
           elseif (full_complexity.eq.0 .and. RH .gt. 0.4) then
              H2SO4_hcl=0.30993*(RH**2)-0.99171*RH+0.66913
           endif
        endif

        if (ammonium_sulphate .gt. 0.0d0) then
           !--NH42SO4---
           NH42SO4_hcl=22.071*(RH**3)-40.678*(RH**2)+27.893*RH-9.4338
        endif

        if (ammonium_nitrate .gt. 0.0d0) then
           !--nh4no3---
           NH4NO3_hcl=19.935*(RH**3)-42.335*(RH**2)+31.275*(RH)-8.8675
        endif

        if (ammonium_chloride .gt. 0.0d0) then
           !--nh4cl---
           if (full_complexity.eq.1 .or. RH .le. 0.4) then
              NH4Cl_hcl=2.8048*(RH**3)-4.3182*(RH**2)+3.1971*(RH)-1.6824
           elseif (full_complexity.eq.0 .and. RH .gt. 0.4) then
              NH4Cl_hcl=1.2304*(RH**2)-0.18262*(RH)-1.0643
           endif
        endif

        if (sodium_sulphate .gt. 0.0d0) then
           !--Na2SO4---
           Na2SO4_hcl=36.104*(RH**4)-78.658*(RH**3)+63.441*(RH**2)-26.727*RH+5.7007
        endif

        if (sodium_nitrate .gt. 0.0d0) then
           !--NaNO3---
           if (full_complexity.eq.1 .or. RH .le. 0.4) then
              NaNO3_hcl=54.471*(RH**5)-159.42*(RH**4)+180.25*(RH**3)-98.176*(RH**2)+25.309*RH-2.4275
           elseif (full_complexity.eq.0 .and. RH .gt. 0.4) then
              NaNO3_hcl=21.632*(RH**4)-53.088*(RH**3)+47.285*(RH**2)-18.519*RH+2.6846
           endif
        endif

        if (sodium_chloride .gt. 0.0d0) then
           !--NaCl---
           if (full_complexity.eq.1 .or. RH .le. 0.4) then
              NaCl_hcl=5.4138*(RH**4)-12.079*(RH**3)+9.627*(RH**2)-3.3164*RH+0.35224
           elseif (full_complexity.eq.0 .and. RH .gt. 0.4) then
              NaCl_hcl=2.432*(RH**3)-4.3453*(RH**2)+2.3834*RH-0.4762
           endif
        endif

        Ln_HCl_act=binary_hcl+nitric_acid_eq_frac*HNO3_hcl+sulphuric_acid_eq_frac*H2SO4_hcl+&
             &ammonium_sulphate_eq_frac*NH42SO4_hcl+ammonium_nitrate_eq_frac*NH4NO3_hcl+ammonium_chloride_eq_frac*NH4Cl_hcl+&
             &sodium_sulphate_eq_frac*Na2SO4_hcl+sodium_nitrate_eq_frac*NaNO3_hcl+sodium_chloride_eq_frac*NaCl_hcl


        gamma_hcl=exp(Ln_HCl_act)
        gamma_out(2) = gamma_hcl


        !equilibrium constant
        !K_hcl=1.97d6 !Zaveri et al (2005)
        K_hcl = 2d6*exp(9000d0*henrys_temp_dep) ! after Wagman et al (1982) (and NIST database)
        Press_HCl=(ions_mol(1)*ions_mol(7)*(gamma_hcl**2))/K_hcl

     endif


     ! Juha: for ion malility output
     mols_out = ions_mol


     !REFERENCES
     !Clegg et al (1998) A Thermodynamic Model of the System H+-NH4+-Na+-SO42- -NO3--Cl--H2O at 298.15 K, J. Phys. Chem., 102A, 2155-2171. 
     !Clegg et al (2001) Thermodynamic modelling of aqueous aerosols containing electrolytes and dissolved organic compounds. Journal of Aerosol Science 2001;32(6):713-738.
     !Topping et al (2005a) A curved multi-component aerosol hygroscopicity model framework: Part 1 - Inorganic compounds. Atmospheric Chemistry and Physics 2005;5:1205-1222.
     !Topping et al (2005b) A curved multi-component aerosol hygroscopicity model framework: Part 2 - Including organic compounds. Atmospheric Chemistry and Physics 2005;5:1223-1242.
     !Zaveri et al (2005). A new method for multicomponent activity coefficients of electrolytes in aqueous atmospheric aerosols, JGR, 110, D02201, 2005.
   end subroutine inorganic_pdfite



end module aerosol_thermodynamics















