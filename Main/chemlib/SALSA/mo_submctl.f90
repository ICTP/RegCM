!! mo_submctl contains parameters, switches and initialization routines for the m7 aerosol module.
!!
!=========================================================================================
! Modifications 05/2014 (J Tonttila):
! All arrays whose dimension depends on the number of size bins for aerosols/clouds are now 
! allocated and initialized in *salsa_initialize* in order to adjust the size distribution 
! properties and resolution by using a NAMELIST.
!=========================================================================================




MODULE mo_submctl

  USE mo_kind,             ONLY: dp

  IMPLICIT NONE

  PRIVATE

  !M7 and SALSA
  PUBLIC :: lspropupdate,lsdistupdate,nsnucl

  PUBLIC :: oldupdate, debug

  !SALSA:
  PUBLIC :: locgas, lsol2b, act_coeff,nj3, time

  PUBLIC :: in1a,in2a,in2b,fn1a,fn2a,fn2b,nbins
  PUBLIC :: nbin, nbin2, nbin3,reglim,nlim,prlim,nreg
  PUBLIC :: epsv,vhilim,vlolim,vratiohi,vratiolo,dpmid,sigma
  PUBLIC :: pi, pi6, rg, avog, boltz, cpa, mair, grav
  PUBLIC :: rhosu,rhooc, rhobc,rhoss, rhodu, rhowa, rhonh, rhono
  PUBLIC :: msu,mdu,mno,mnh,n3,massacc,d_sa,pstand,mss,mbc,moc,epsoc,mwa,slim,ions,&
            mvsu,mvoc,mvss,surfw0,mvwa,mvno,mvnh
  PUBLIC :: recalc
  PUBLIC :: csr_strat_wat,csr_strat_mix,csr_strat_ice,csr_conv,zbcr,cfracn,zfracn,zfracn_cv

  PUBLIC :: t_section,t_parallelbin
  PUBLIC :: ncldbin,ica,fca,icb,fcb,ira,fra,ncld,nprc,dmincld

  PUBLIC :: aerobins, cloudbins, precpbins

  PUBLIC :: nlcoag,                &
            nlcgaa,nlcgcc,nlcgpp,  &
            nlcgca,nlcgpa,nlcgpc,  &
            nlcnd,                 &
            nlcndgas,              &
            nlcndh2oae,nlcndh2ocl, &
            nlauto,                &
            nlactiv,               &
            nlactbase,nlactintst,  &
            
            lscoag,                &
            lscgaa,lscgcc,lscgpp,  &
            lscgca,lscgpa,lscgpc,  &
            lscnd,                 &
            lscndgas,              &
            lscndh2oae,lscndh2ocl, &
            lsauto,                &
            lsactiv,               &
            lsactbase,lsactintst

  PUBLIC :: nspec, listspec, maxspec, nmod, prcpfm
  PUBLIC :: sigmag, dpg, n, volDistA, volDistB, nf2a
  PUBLIC :: isdtyp

  TYPE t_section
     ! Used to store information for the sectional size distribution of cloud/drizzle
     ! droplets(and later aerosols?)
     REAL(dp) :: vhilim,     & ! bin volume at the high limit
                 vlolim,     & ! - '' - at the low limit
                 vratiohi,   & ! volume ratio between the center and high limit
                 vratiolo,   & ! - '' - and the low limit
                 dmid,       & ! bin middle diameter
                 !******************************************************
                 ! ^ Do NOT change the stuff above after initialization !
                 !******************************************************
                 dwet,       & ! Wet diameter or mean droplet diameter

                 volc(8),    & ! Volume concentrations of aerosol species + water
                               ! Since most of the stuff in SALSA is hard coded, these
                               ! *have to be* in the order: 1:SO4, 2:OC, 3:BC, 4:DU, 5:SS, 6:NO, 7:NH, 8:H2O

                 veqh2o,     & ! Equilibrium h2o concentration for each particle
                 numc,       & ! Number concentration of particles/droplets 
                 core          ! Volume of dry particle
  END TYPE t_section

  TYPE t_parallelbin
     ! Map bin indices between parallel size distributions
     INTEGER :: cur  ! Index for current distribution
     INTEGER :: par  ! Index for corresponding parallel distribution
  END TYPE t_parallelbin

  !--- 1) Define and pre-set switches for the processes of M7: -----------------------

  !--- Physical:
  !Switches for both SALSA aerosol microphysical processes
  LOGICAL :: oldupdate  = .FALSE.
  LOGICAL :: debug      = .FALSE.

  ! Process switches: nl* is read from the NAMELIST and NOT changed during runtime.
  !                   ls* is the switch actually used and will get the value of nl* 
  !                   except for special circumstances such as spinup period etc.
  LOGICAL :: nlcoag    = .FALSE., & ! Coagulation master switch
             lscoag
    LOGICAL :: nlcgaa  = .TRUE., & ! Coagulation between aerosols
               lscgaa
    LOGICAL :: nlcgcc  = .TRUE., & ! Collision-coalescence between cloud droplets
               lscgcc
    LOGICAL :: nlcgca  = .TRUE., & ! Cloud collection of aerosols
               lscgca
    LOGICAL :: nlcgpc  = .TRUE., & ! Collection of cloud droplets by rain
               lscgpc
    LOGICAL :: nlcgpa  = .TRUE., & ! Collection of aerosols by rain
               lscgpa
    LOGICAL :: nlcgpp  = .TRUE., & ! Collision between rain drops
               lscgpp

    LOGICAL :: nlcnd     = .TRUE.,   & ! Condensation 
               lscnd
    LOGICAL :: nlcndgas  = .TRUE.,   & ! Condensation of precursor gases
               lscndgas
    LOGICAL :: nlcndh2ocl  = .TRUE., & ! Condensation of water vapour on clouds (drizzle)
               lscndh2ocl
    LOGICAL :: nlcndh2oae = .TRUE.,  & ! Condensation of water vapour on aerosol particles (FALSE -> equilibrium calc.) 
               lscndh2oae = .TRUE.

    LOGICAL :: nlauto    = .FALSE.,   & ! Autoconversion of cloud droplets (needs activation)
               lsauto
    LOGICAL :: nlactiv   = .FALSE.,   & ! Master switch for cloud droplet activation
               lsactiv              
    LOGICAL :: nlactintst = .TRUE.,  & ! Switch for interstitial activation: Use particle wet size determined by
               lsactintst = .TRUE.     ! codensation equations and supersaturation directly from the host model
    LOGICAL :: nlactbase = .TRUE.,   & ! Switch for cloud base activation: Use the regular parameterized method
               lsactbase = .TRUE.      ! for maximum supersaturation and cloud activation.


  LOGICAL :: lsdistupdate = .TRUE.  ! Perform the size distribution update
  
  LOGICAL :: lspropupdate = .FALSE.  ! Update diagnostic particle properties between processes

  ! 1) Switches for aerosol microphysical processes ------------------------
  INTEGER, PARAMETER :: nmod = 7
  INTEGER :: nwater     = 1         ! Aerosol water uptake scheme:
                                    !
                                    ! nwater = 0 Jacobson et al., JGR 1996
                                    !        = 1 Kappa-Koehler theory based approach (Petters and Kreidenweis, ACP 2007)

  INTEGER :: nsnucl     = 0         ! Choice of the H2SO4/H2O nucleation scheme:
                                    ! SALSA:
                                    ! 0 = off   
                                    ! 1 = binary nucleation
                                    ! 2 = activation type nucleation
                                    ! 3 = kinetic nucleation
                                    ! 4 = ternary nucleation
                                    ! 5 = nucleation with ORGANICs
                                    ! 6 = activation type of nucleation with H2SO4+ORG
                                    ! 7 = heteromolecular nucleation with H2SO4*ORG
                                    ! 8 = homomolecular nucleation of  H2SO4 + 
                                    !           heteromolecular nucleation with H2SO4*ORG
                                    ! 9 = homomolecular nucleation of  H2SO4 and ORG + 
                                    !           heteromolecular nucleation with H2SO4*ORG

  LOGICAL :: prcpfm   = .TRUE.      ! Choice of precipitation formation, works with switched off
				    ! = .FALSE. no precipitation 
				    ! = .TRUE. coagulation based formation 

  LOGICAL :: lgcr       = .TRUE.    ! Calculate ionization due to galactic cosmic rays
  
  LOGICAL :: locgas = .FALSE.,&   ! emission of organic carbon in gas phase
             lsol2b = .FALSE.     ! repartitioning of insoluble material in 
                                  ! case of increase in solubility 

  LOGICAL :: recalc   = .FALSE.   ! recalculation of wet diameter between
                                  ! calculation of microphysical processes

  INTEGER ::                    & ! J3 parametrization
             nj3 = 1              ! 1 = condensational sink (Kerminen&Kulmala, 2002)
                                  ! 2 = coagulational sink (Lehtinen et al. 2007)
                                  ! 3 = coagS+self-coagulation (Anttila et al. 2010)
  REAL(dp) :: act_coeff=1.e-7_dp  ! activation coefficient

  ! Define which aerosol species used and initial size distributions
  INTEGER :: nspec = 1
  INTEGER, PARAMETER :: maxspec = 7
  CHARACTER(len=3) :: listspec(maxspec) = (/'SO4','   ','   ','   ','   ','   ','   '/)

  ! Volume fractions between aerosol species for A and B-bins
  REAL(dp) :: volDistA(maxspec) = (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
  REAL(dp) :: volDistB(maxspec) = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
  ! Number fraction allocated to a-bins in regime 2 (b-bins will get 1-nf2a)
  REAL(dp) :: nf2a = 1.0   

  REAL(dp) :: time 

  INTEGER :: isdtyp = 0  ! Type of input aerosol size distribution: 0 - Uniform 
                         !                                          1 - Read vertical profile of the mode 
                         !                                              number concentrations from an input file

  REAL(dp) :: sigmag(nmod) = (/2.0,2.0,2.0,2.0,2.0,2.0,2.0/),               & ! Stdev
              dpg(nmod) = (/0.03, 0.15, 0.2, 0.2, 0.2, 0.2, 0.2/),    & ! Mode diam in um
              n(nmod) = (/1600.,640.,0.,0.,0.,0.,0./)                   ! #/mg ~ #/cm3



  INTEGER, PARAMETER ::            &
   nreg = 2                          ! number of main size regimes


  REAL(dp) ::                       &
   reglim(nreg+2) =                            & ! low/high diameter limits
    (/ 3.e-9_dp, 5.e-8_dp, 7.e-7_dp, 1.e-5_dp /) ! of main size regimes [m]  

   INTEGER :: &
   nbin(nreg) = (/ 3, 7 /)   ! number of bins in each main regime 

  INTEGER ::      &
   nbin2, & != 4,                & ! number of bins in former 2-region
   nbin3 != nbin(2) - nbin2     ! number of bins in former 3-region

  INTEGER ::      & ! number/radius: start index
   in1a, & != 1,                 & ! regime 1a
   in2a, & != in1a + nbin(1),    & ! regime 2a
   in2b, & != in2a + nbin(2),    & ! regime 2b

                               ! number/radius: last index
   fn1a, & != in2a - 1,          & ! regime 1a
   fn2a, & != fn1a + nbin(2),    & ! regime 2a
   fn2b, & != fn2a + nbin(2),    & ! regime 2b

   nbins != fn2b                ! total number of size bins

   ! Juha: Cloud and rain bins:
  INTEGER :: ncldbin(2) = (/7,7/)        ! Number of bins for cloud bins in regime a and b 
                                         
  TYPE(t_parallelbin) ::   ica, & ! cloud droplets (first, regime a)
                           fca, & ! cloud droplets (last, regime a)
                           icb, & ! cloud droplets (first, regime b)
                           fcb    ! cloud droplets (last, regime b)
  INTEGER             ::   ira,fra! Rain/drizzle bin indices      
  INTEGER ::               ncld   ! Total number of cloud bins
  INTEGER ::               nprc   ! Totla number of precipitation bins
  REAL(dp) :: dmincld = 5.e-8_dp    !Minimum diameter for the cloud droplet regime in terms of the ccn dry radius. The first cloud droplet bin
                                    !is taken to coincide with the smallest full aerosol bin that conforms with this diameter **NAMELISTIIN        
  
  REAL(dp), ALLOCATABLE :: aerobins(:),  &  ! These are just to deliver information about the bin diameters if the 
                           cloudbins(:), &  ! host model needs it (lower limits).
                           precpbins(:)


  !!!! SOME OF THIS IS REPLACED BY THE T_SECTION DATATYPE
  REAL(dp), ALLOCATABLE :: epsv(:),           &
                           vhilim(:),         &
                           vlolim(:),         &
                           vratiohi(:),       &
                           vratiolo(:),       &
                           dpmid(:),          &
                           sigma(:),          &
                           csr_strat_wat(:),  &
                           csr_strat_mix(:),  &
                           csr_strat_ice(:),  &
                           csr_conv(:),       &
                           zbcr(:)


  REAL(dp), PARAMETER ::     &
   avog   = 6.0221e+23_dp,   & ! Avogadro number (#/mol)
   boltz  = 1.3807e-23_dp,   & ! Boltzmann constant (J/K)
   grav   = 9.81_dp,         & ! gravitational acceleration (m/s^2)
   pstand = 1.01325e+5_dp,   & ! standard pressure (Pa)
   rg     = 8.314_dp,        & ! molar gas constant (J/(mol K)) 
   pi     = 3.1415927_dp,    & ! self explanatory
   pi6    = 0.5235988_dp,    & ! pi/6
   cpa    = 1010._dp,        & ! specific heat of dry air, constant P (J/kg/K)
   mair   = 28.97e-3_dp,     & ! molar mass of air (mol/kg)
   deltav = 1.096e-7_dp,     & ! vapor jump length (m)
   deltaT = 2.16e-7_dp,      & ! thermal jump length (m)
   alphaT = 0.96_dp,         & ! thermal accomodation coefficient
   alphac = 1.0_dp             ! condensation coefficient

  REAL(dp), PARAMETER ::     & ! molar mass [kg/mol]
   msu = 98.08e-3_dp,        & ! sulphate   
   mno = 62.01e-3_dp,        & ! HNO3 
   mnh = 18.04e-3_dp,        & ! NH3 
   moc = 150.e-3_dp,         & ! organic carbon
   mbc = 12.e-3_dp,          & ! black carbon
   mss = 58.44e-3_dp,        & ! sea salt (NaCl)
   mdu = 100.e-3_dp,         & ! mineral dust
   mwa = 18.016e-3_dp,       & ! water
   mas = 132.14e-3_dp,       & ! ammoniums sulphate ((NH4)2SO4)
                               !
                               ! densities [kg/m3]
   rhosu = 1830._dp,         & ! sulphate
   rhono = 1479._dp,         & ! HNO3
   rhonh = 1530._dp,         & ! NH3 
   rhooc = 2000._dp,         & ! organic carbon
   rhobc = 2000._dp,         & ! black carbon
   rhoss = 2165._dp,         & ! sea salt (NaCl)
   rhodu = 2650._dp,         & ! mineral dust
   rhowa = 1000._dp,         & ! water 
                               !
                               ! volume of molecule [kg/#]
   mvsu = msu/avog/rhosu,    & ! sulphate
   mvno = mno/avog/rhono,    & ! HNO3 
   mvnh = mnh/avog/rhonh,    & ! NH3  
   mvoc = moc/avog/rhooc,    & ! organic carbon
   mvss = mss/avog/rhoss,    & ! sea salt
   mvwa = mwa/avog/rhowa,    &
                              !
   volratio =                & ! ratio of molecular volumes for
    (msu*rhoss)/(rhosu*mss), & ! sulphate and sea salt
                               !
   n3 = 158.79_dp               ! number of H2SO4 molecules in 3 nm cluster 
                               !  assuming d_sa = 5.54 ???     
  !-- 4.3) Properties of condensing vapours

  REAL(dp), PARAMETER ::                               & ! diameter of condensing molecule [m]
      d_sa   = 5.539376964394570e-10_dp,               &

      d_oc   = 6.195906936656752e-10_dp,               &

      d_h2o  = 3.851565216195334e-10_dp

  REAL(dp), PARAMETER :: &
       slim = 1.005_dp,  & ! water saturation used as limit
       ions = 3.0_dp,    & ! van't Hoff factor (ions produced upon dissociation)
       surfw0 = 0.073_dp,& ! surface tension of pure water @ ~ 293 K [J/m2]
       epsoc = 0.15_dp     ! water uptake of organic material

  !-- 7) Parameters for cloud activation
  ! OSAN NÄISTÄ VOIS POISTAA

  REAL(dp), PARAMETER :: crcut=0.035*1E-6_dp ! Assumed lower cut-off of the
                                             ! aerosol size distribution [m]

  !--- Ulrike: included for activation in convective clouds
  REAL(dp), PARAMETER :: crcut_cv=0.025*1E-6_dp ! Assumed lower cut-off of the
  


  REAL(dp), ALLOCATABLE :: cfracn(:)
  REAL(dp), ALLOCATABLE :: zfracn(:)
  REAL(dp), ALLOCATABLE :: zfracn_cv(:)
  REAL(dp), ALLOCATABLE :: massacc(:)



  REAL(dp), PARAMETER :: &
   nlim = 1.e-1_dp,         & ! number conc. limit below which bin empty  [#/m3] 
   prlim = 1.e-40_dp,     & ! The same for precipitation drops for which concentrations are normally much lower [#/m3]
   m3_2_um3 = 1.e+18_dp    ! conversion factor for volume from m3 to um3

  INTEGER, ALLOCATABLE, PUBLIC :: iso4b(:), inob(:), inhb(:),   &
                                  iocb(:),  ibcb(:),            &
                                  idub(:),  issb(:)
  

  !--- 12) Service routines for initialization and auxiliary computations ----------

END MODULE mo_submctl
