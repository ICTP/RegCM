!****************************************************************
!*                                                              *
!*   module MOD_AERO_INIT                                       *
!*                                                              *
!*   Contains subroutines and functions that are used           *
!*   to initialize the particle size grid and aerosol           *
!*   processes                                                  *
!*                                                              *
!****************************************************************
!----------------------------------------------------------------------------
!!$   Copyright 2014 Atmospheric Research Centre of Eastern Finland,
!!$         Finnish Meteorological Institute, Kuopio, Finland
!!$
!!$   Licensed under the Apache License, Version 2.0 (the "License");
!!$   you may not use this file except in compliance with the License.
!!$   You may obtain a copy of the License at
!!$
!!$       http://www.apache.org/licenses/LICENSE-2.0
!!$
!!$   Unless required by applicable law or agreed to in writing, software
!!$   distributed under the License is distributed on an "AS IS" BASIS,
!!$   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!$   See the License for the specific language governing permissions and
!!$   limitations under the License.
!----------------------------------------------------------------------------
MODULE mo_salsa_init

  IMPLICIT NONE

CONTAINS

  ! fxm: when dpmid is used for calculating coagulation coefficients
  ! (only?), would it make more sense to use approximate wet radii
  ! e.g. for sea salt particles?
  !********************************************************************
  !
  ! subroutine SET_SIZEBINS()
  !
  !********************************************************************
  !
  ! Purpose:
  ! --------
  ! Initializes particle size distribution grid by 
  ! calculating size bin limits and mid-size for
  ! *dry* particles in each bin
  !
  !
  ! Method:
  ! ------- 
  ! Size distribution described using
  !   1) moving center method (regimes 1 and 2)
  !   (Jacobson, Atmos. Env., 31, 131-144, 1997)
  !   2) fixed sectional method (regime 3)
  ! 
  ! Size bins in each regime are spaced logarithmically
  ! based on given regime size limits and bin number.
  !
  !
  ! Interface:
  ! ----------
  ! Called from model driver
  ! (only at the beginning of simulation)
  !
  !
  ! Coded by:
  ! ---------
  ! Hannele Korhonen (FMI) 2005 
  ! Harri Kokkola (FMI) 2006
  !
  ! Bug fixes for box model + updated for the new aerosol datatype:
  ! ----------
  ! Juha Tonttila (FMI) 2014
  !
  !---------------------------------------------------------------------

  SUBROUTINE set_sizebins()

    USE mo_kind, ONLY : dp

    ! POISTA MITÄ EI TARVII
    USE mo_submctl, ONLY : &
         pi6,         & ! pi/6
         reglim,      & ! diameter limits for size regimes [m]
         nbin,        & ! number of size bins in each (sub)regime
         in1a, fn1a,  & ! size regime bin indices: 1a
         in2a, fn2a,  & !     - " -       2a
         in2b, fn2b,  & !     - " -       2b
         nbin2, nbin3,& !
         nbins,  &
         sigma,  &
         aerobins

    USE mo_salsa_driver, ONLY : &
         kbdim,kproma,klev,  &
         aero


    IMPLICIT NONE

    !-- local variables ----------
    INTEGER :: ii, jj,cc,dd,vv ! loop indices

    REAL(dp) ::  &
         ratio,  & ! ratio of regime upper and lower diameter
         vlo,    & ! 
         vhi


    ! -------------------------
    ! Allocate aerosol tracers
    !--------------------------
    ALLOCATE(aero(kbdim,klev,nbins))

    ! Juha: Corrected some bugs with the section volume definitions
    !-------------------------------------------------------------------------------

    DO jj = 1,klev
       DO ii = 1,kproma

          !-- 1) size regime 1: --------------------------------------
          !  - minimum & maximum *dry* volumes [fxm] 
          !  - bin mid *dry* diameter [m]
          ratio = reglim(2)/reglim(1)   ! section spacing

          DO cc = in1a,fn1a
             aero(ii,jj,cc)%vlolim = pi6*(reglim(1)*ratio**(real(cc-1)/nbin(1)))**3
             aero(ii,jj,cc)%vhilim = pi6*(reglim(1)*ratio**(real(cc)/nbin(1)))**3
             aero(ii,jj,cc)%dmid = ( (aero(ii,jj,cc)%vhilim + aero(ii,jj,cc)%vlolim) /  &
                                     (2._dp*pi6) )**(1._dp/3._dp)
             aero(ii,jj,cc)%vratiohi = aero(ii,jj,cc)%vhilim/(pi6*aero(ii,jj,cc)%dmid**3)
             aero(ii,jj,cc)%vratiolo = aero(ii,jj,cc)%vlolim/(pi6*aero(ii,jj,cc)%dmid**3)
          END DO

          !-- 2) size regime 2: --------------------------------------
          !  - minimum & maximum *dry* volumes [fxm] 
          !  - bin mid *dry* diameter [m]

          !-- 2.1) first for subregime 2a
          ratio = reglim(3)/reglim(2)   ! section spacing

          DO dd = in2a,in2a+nbin2-1
             cc = dd - in2a
             aero(ii,jj,dd)%vlolim = pi6*(reglim(2)*ratio**(real(cc)/nbin2))**3
             aero(ii,jj,dd)%vhilim = pi6*(reglim(2)*ratio**(real(cc+1)/nbin2))**3
             aero(ii,jj,dd)%dmid = ( (aero(ii,jj,dd)%vhilim + aero(ii,jj,dd)%vlolim) /  &
                                     (2._dp*pi6) )**(1._dp/3._dp)
             aero(ii,jj,dd)%vratiohi = aero(ii,jj,dd)%vhilim/(pi6*aero(ii,jj,dd)%dmid**3)
             aero(ii,jj,dd)%vratiolo = aero(ii,jj,dd)%vlolim/(pi6*aero(ii,jj,dd)%dmid**3)
          END DO

          !-- 3) size regime 3: --------------------------------------
          !  - bin mid *dry* diameter [m]
          ratio = reglim(4)/reglim(3)   ! section spacing

          DO dd = in2a+nbin2,fn2a
             cc = dd - (fn2a-(nbin3-1))
 
             aero(ii,jj,dd)%vlolim = pi6*(reglim(3)*ratio**(real(cc)/nbin3))**3
             aero(ii,jj,dd)%vhilim = pi6*(reglim(3)*ratio**(real(cc+1)/nbin3))**3
             aero(ii,jj,dd)%dmid = ( (aero(ii,jj,dd)%vhilim + aero(ii,jj,dd)%vlolim) /  &
                                     (2._dp*pi6) )**(1._dp/3._dp)
             aero(ii,jj,dd)%vratiohi = aero(ii,jj,dd)%vhilim/(pi6*aero(ii,jj,dd)%dmid**3)
             aero(ii,jj,dd)%vratiolo = aero(ii,jj,dd)%vlolim/(pi6*aero(ii,jj,dd)%dmid**3)             
          END DO

          !-- 2.2) same values for subregime 2b
          aero(ii,jj,in2b:fn2b)%vlolim = aero(ii,jj,in2a:fn2a)%vlolim
          aero(ii,jj,in2b:fn2b)%vhilim = aero(ii,jj,in2a:fn2a)%vhilim
          aero(ii,jj,in2b:fn2b)%dmid = aero(ii,jj,in2a:fn2a)%dmid
          aero(ii,jj,in2b:fn2b)%vratiohi = aero(ii,jj,in2a:fn2a)%vratiohi
          aero(ii,jj,in2b:fn2b)%vratiolo = aero(ii,jj,in2a:fn2a)%vratiolo

          ! Initialize the wet diameter with the bin dry diameter to avoid numerical proplems later
          aero(ii,jj,:)%dwet = aero(ii,jj,:)%dmid

          ! Set volume and number concentrations to zero
          aero(ii,jj,:)%numc = 0._dp
          aero(ii,jj,:)%core = 0._dp
          DO vv=1,8
             aero(ii,jj,:)%volc(vv) = 0._dp
          END DO

       END DO !ii
    END DO !!jj

    ! Save bin limits to be delivered e.g. to host model if needed
    ALLOCATE(aerobins(nbins))
    DO cc = 1,nbins
       aerobins(cc) = (aero(1,1,cc)%vlolim/pi6)**(1./3.)
    END DO
    aerobins = 0.5_dp*aerobins ! to radius

  END SUBROUTINE set_sizebins

  !--------------------------------------------------------------------------
  !
  ! *******************************
  ! SUBROUTINE set_cloudbins
  ! *******************************
  !
  ! Setup of hydrometeor size bins similar to the subroutine *set_sizebins*
  !
  ! Juha Tonttila (FMI) 2014
  !
  !---------------------------------------------------------------------------

  SUBROUTINE set_cloudbins()
    USE mo_kind, ONLY : dp
    !POISTA MITÄ EI TARVII
    USE mo_submctl, ONLY : pi6,             &
                               ica,icb,         &
                               fca,fcb,         &
                               ira,fra,         &
                               ncld,nprc,        &
                               ncldbin,         & 
                               dmincld,         &
                               in1a,fn2a,       &
                               in2b,fn2b,       &
                               cloudbins,       &
                               precpbins
    USE mo_salsa_driver, ONLY : kproma, kbdim, klev, &
                                 cloud,precp,aero

                                      
    IMPLICIT NONE

    INTEGER :: ii,jj,cc,bb

    ! For cloud bins
    REAL(dp) :: zvoldiff(fn2a)
    INTEGER :: zindex(fn2a)
    INTEGER :: imin, nba,nbb
    INTEGER :: armin(1)
    LOGICAL :: l_min(fn2a)

    REAL(dp) :: tmplolim(7), tmphilim(7), tmpmid(7)

    ! Helper arrays to set up precipitation size bins
    tmplolim = (/50.,55.,65.,100.,200.,500.,1000./)*1.e-6_dp
    tmphilim = (/55.,65.,100.,200.,500.,1000.,2000./)*1.e-6_dp
    !tmpmid = (tmplolim + tmphilim)/2._dp ÄLÄ KÄYTÄ


    armin = 0
    
    ! For setting cloud droplet bins
    DO ii = in1a,fn2a
       zindex(ii) = ii
    END DO

    ! Determine the smallest full aerosol bin with which the smallest cloud droplet bin will coincide
    zvoldiff(in1a:fn2a) = ABS(aero(1,1,in1a:fn2a)%vlolim - pi6*dmincld**3)
    l_min = ( zvoldiff(in1a:fn2a) == MINVAL(zvoldiff(in1a:fn2a)) )
    IF (COUNT(l_min) /= 1) STOP 'Initialization error! Try changing the cloud droplet bin lower bound slightly'
    armin = PACK(zindex,l_min)
    imin = armin(1)
    IF (aero(1,1,imin)%vlolim < pi6*dmincld**3) imin=imin+1
    IF (imin < in1a .OR. imin > fn2a) STOP 'Invalid first cloud bin'
 
    ! Number of cloud bins in regime a (soluble nuclei)
    nba = fn2a-imin+1
    ! Number of cloud bins in regime b (insoluble nuclei)
    IF (aero(1,1,imin)%dmid > aero(1,1,in2b)%dmid) THEN
       ! All cloud bins are within regime 2
       nbb = nba
    ELSE
       ! The smallest regime a cloud bins go to regime 1
       nbb = fn2b-fn2a
    END IF

    ncldbin(1) = nba
    ncldbin(2) = nbb

    ! Reset cloud bin indices accordingly. The two components give the cloud regime index, 
    ! and the aerosol bin index with which they are parallel
    ica%cur = 1;                       ica%par = imin
    fca%cur = ica%cur + ncldbin(1)-1; fca%par = ica%par + ncldbin(1)-1
    icb%cur = fca%cur + 1;             icb%par = fn2b - ncldbin(2) + 1
    fcb%cur = icb%cur + ncldbin(2)-1; fcb%par = icb%par + ncldbin(2)-1

    ncld = fcb%cur

    ! Rain/drizzle bins
    ira = 1; fra = 7;
    nprc = fra

    ! ----------------------------------------
    ! Allocate cloud and precipitation arrays
    ! ----------------------------------------
    ALLOCATE(cloud(kbdim,klev,ncld), precp(kbdim,klev,nprc))

    DO jj = 1,klev
       DO ii = 1,kproma

          ! -------------------------------------------------
          ! Set cloud properties (parallel to aerosol bins)
          ! -------------------------------------------------
          cloud(ii,jj,ica%cur:fca%cur)%vhilim = aero(ii,jj,ica%par:fca%par)%vhilim
          cloud(ii,jj,icb%cur:fcb%cur)%vhilim = aero(ii,jj,icb%par:fcb%par)%vhilim

          cloud(ii,jj,ica%cur:fca%cur)%vlolim = aero(ii,jj,ica%par:fca%par)%vlolim
          cloud(ii,jj,icb%cur:fcb%cur)%vlolim = aero(ii,jj,icb%par:fcb%par)%vlolim

          cloud(ii,jj,ica%cur:fca%cur)%vratiohi = aero(ii,jj,ica%par:fca%par)%vratiohi
          cloud(ii,jj,icb%cur:fcb%cur)%vratiohi = aero(ii,jj,icb%par:fcb%par)%vratiohi

          cloud(ii,jj,ica%cur:fca%cur)%vratiolo = aero(ii,jj,ica%par:fca%par)%vratiolo
          cloud(ii,jj,icb%cur:fcb%cur)%vratiolo = aero(ii,jj,icb%par:fcb%par)%vratiolo

          cloud(ii,jj,ica%cur:fca%cur)%dmid = aero(ii,jj,ica%par:fca%par)%dmid
          cloud(ii,jj,icb%cur:fcb%cur)%dmid = aero(ii,jj,icb%par:fcb%par)%dmid

          ! Initialize the droplet diameter ("wet diameter") as the dry 
          ! mid diameter of the nucleus to avoid problems later.
          cloud(ii,jj,ica%cur:fcb%cur)%dwet = cloud(ii,jj,ica%cur:fcb%cur)%dmid

          ! Initialize the volume and number concentrations for clouds.
          ! First "real" values are only obtained upon the first calculation
          ! of the cloud droplet activation.
          DO cc = 1,8
             cloud(ii,jj,ica%cur:fcb%cur)%volc(cc) = 0._dp
          END DO
          cloud(ii,jj,ica%cur:fcb%cur)%numc = 0._dp
          cloud(ii,jj,ica%cur:fcb%cur)%core = 0._dp
          cloud(ii,jj,ica%cur:fcb%cur)%veqh2o = 0._dp

          ! ---------------------------------------------------------------------------------------
          ! Set the precipitation properties; unlike aerosol and cloud bins, the size distribution 
          ! goes according to the *wet* radius
          ! ---------------------------------------------------------------------------------------
          DO bb = ira,fra
             precp(ii,jj,bb)%vhilim = pi6*tmphilim(bb)**3
             precp(ii,jj,bb)%vlolim = pi6*tmplolim(bb)**3
             precp(ii,jj,bb)%dmid = ( (precp(ii,jj,bb)%vlolim + precp(ii,jj,bb)%vhilim) /  &
                                      (2._dp*pi6) )**(1._dp/3._dp)
             precp(ii,jj,bb)%vratiohi = precp(ii,jj,bb)%vhilim / ( pi6*precp(ii,jj,bb)%dmid**3 )
             precp(ii,jj,bb)%vratiolo = precp(ii,jj,bb)%vlolim / ( pi6*precp(ii,jj,bb)%dmid**3 )
             
             ! Initialize the wet diameter as the bin mid diameter
             precp(ii,jj,bb)%dwet = precp(ii,jj,bb)%dmid

             DO cc = 1,8
                precp(ii,jj,bb)%volc(cc) = 0._dp
             END DO
             precp(ii,jj,bb)%numc = 0._dp
             precp(ii,jj,bb)%core = 0._dp
             precp(ii,jj,bb)%veqh2o = 0._dp

          END DO
          
       END DO
    END DO

    ! Save bin limits to be delivered e.g. to host model if needed
    ALLOCATE(cloudbins(ncld))
    DO bb = 1,ncld
       cloudbins(bb) = (cloud(1,1,bb)%vlolim/pi6)**(1./3.)
    END DO
    cloudbins = 0.5_dp*cloudbins ! To radius
    ALLOCATE(precpbins(nprc))
    DO bb = 1,nprc
       precpbins(bb) = (precp(1,1,bb)%vlolim/pi6)**(1./3.)
    END DO
    precpbins = 0.5*precpbins ! To radius

  END SUBROUTINE set_cloudbins



  ! fxm: is it really the upper limit of the particles that we want to
  !  use as a criterion? how about supersat of 0.5%.
  !*********************************************************************
  !
  !  subroutine ACTCURVE()
  !
  !*********************************************************************
  !
  ! Purpose:
  ! --------
  ! Calculates the criterion based on which particles from
  !  insoluble sections are moved into soluble sections
  !
  !
  ! Method: 
  ! -------  
  ! Critical dry volume ratio of soluble matter calculated at
  !  upper volume boundary (in regime 3: fixed bin volume) of each bin
  !
  ! A fixed critical saturation ratio assumed for this calculation
  !
  ! NB! The critical soluble fraction is calculated according to
  !  sulphate properties
  !
  !
  ! Interface:
  ! ----------
  ! Called from model driver
  ! (only at the beginning of simulation)
  !
  !
  ! Coded by:
  ! ---------
  ! Hannele Korhonen (FMI) 2005 
  !
  ! Juha Tonttila (FMI) 2014: Updated for the new aerosol datatype
  !
  !----------------------------------------------------------------------

  SUBROUTINE actcurve()

    USE mo_kind, ONLY: dp

    USE mo_submctl, ONLY: &
         pi,        & 
         pi6,       & ! pi/6
         rg,        & ! universal gas constant [J/(mol K)]
         in1a,fn2b, & ! size bin indices
         mwa,       & ! molar mass of water [kg/mol]
         rhowa,     & ! density of water [kg/m3]
         surfw0,    & ! surface tension of water [J/m2]
         msu,       & ! molar mass of sulphate [kg/mol]
         rhosu,     & ! density of sulphate [kg/m3]
         ions,      & ! van t'Hoff factor
         slim,      & ! critical saturation value
         epsv         ! critical volume ratio of soluble material (*temp**3)

    USE  mo_salsa_driver, ONLY : &
         aero

    IMPLICIT NONE

    !-- local variables ------
    REAL(dp) :: aa, bb   ! constants in K?hler equation

    !--------------------------------------------------------------------

    aa = 4.*mwa*surfw0/(rg*rhowa)   ! curvature (Kelvin) effect
    bb = 6.*mwa*ions*rhosu/(pi*rhowa*msu) ! solute (Raoult) effect

    !-- Calculation of critical soluble material in each bin
    !   (must be divided with temp**3 when used) ------------------
    !
    !TB changed the calculation of  maximum soluble fraction from 
    !   high volume limit to mean diameter
    !

    !    epsv(in1a:fn2b) = 4._dp*aa**3/(27._dp*bb*(log(slim))**2*vhilim(in1a:fn2b))
    epsv(in1a:fn2b) = 4._dp*aa**3/(27._dp*bb*(log(slim))**2*(pi6*aero(1,1,in1a:fn2b)%dmid**3))

  END SUBROUTINE actcurve


  !----------------------------------------------------------------------
  !
  ! *************************
  ! SUBROUTINE define_salsa
  ! *************************
  !
  ! Reads logical switches and aerosol/hydrometeor size bin definitions
  ! from a namelist.
  !
  ! Juha Tonttila (FMI) 2014
  !
  !----------------------------------------------------------------------
  SUBROUTINE define_salsa()

    USE mo_kind, ONLY : dp
    USE mo_submctl, ONLY : nlcoag,                &
                               nlcgaa,nlcgcc,nlcgpp,  &
                               nlcgca,nlcgpa,nlcgpc,  &
                               nlcnd,                 &
                               nlcndgas,              &
                               nlcndh2oae,nlcndh2ocl, &
                               nlauto,nlactiv,        &
                               nlactintst,nlactbase,  &
                               dmincld,nbin,reglim,   &
                               nspec,listspec,        &
                               volDistA, volDistB,    &
                               nf2a, isdtyp,         &
                               sigmag,dpg,n
                               
    IMPLICIT NONE

    NAMELIST /salsa/    &
         nlcoag,        & ! Coagulation master switch
         nlcgaa,        & ! Coagulation between aerosols
         nlcgcc,        & ! Collision-coalescence between cloud droplets 
         nlcgpp,        & ! Collisions between rain drops
         nlcgca,        & ! Cloud collection of aerosols
         nlcgpa,        & ! Collection of aerosols by precip
         nlcgpc,        & ! Collection of cloud droplets by rain
         nlcnd,         & ! Switch for condensation subroutine
         nlcndgas,      & ! Condensation of precursor gases
         nlcndh2ocl,    & ! Condensation of water vapour on clouds (drizzle)
         nlcndh2oae,    & ! Condensation of water vapour on aerosols (FALSE -> equilibrium calc.)
         nlauto,        & ! Switch for autoconversion of cloud droplets to drizzle and rain
         nlactiv,       & ! Master witch for cloud droplet activation 
         nlactbase,     & ! Switch for parameterized cloud base activation
         nlactintst,    & ! Swithc for interstitial activation based on particle growth and host model S
         isdtyp,        & ! Type of initial size distribution: 0 - uniform; 1 - vertical profile, read from file
         reglim,        & ! Low/high diameter limits of the 2 aerosol size regimes (1d table with length 4)
         nbin,          & ! Number of bins used for each of the aerosol size regimes (1d table with length 2)
         dmincld,       & ! Minimum hydrometeor bin diameter: the lower limit of the smallest hydrometeor bin 
                          ! is set at or above. For now this is also the only way to control the number of
                          ! hydrometeor bins as they are assumed fully parallel with the aerosol bins.
                          ! Implementing different bin limits for hydromets than aerosols would require
                          ! somewhat extensive modification of the microphysical SALSA dynamics calculations
         nspec,         & ! Number of aerosol species used in the model
         listspec,      & ! List of strings specifying the names of the aerosol species that are active.
                          ! Must be an array of length 7, with empty strings for unused stuff.
         volDistA,      & ! Initial relative contribution [0-1] of each species to particle volume in a-bins. Must be
                          ! an array of length 7, with zero for unused species.
         volDistB,      & ! Same as above but for b-bins
         nf2a,          & ! Number fraction of particles allocated to a-bins in regime 2. b-bins will get 1-nf2a

         sigmag,        & ! Stdev for the 7 initial lognormal modes
         dpg,           & ! Mean diameter for the 7 initial lognormal modes
         n                ! Number concentration for the 7 initial lognormal modes

    
    OPEN(11,STATUS='old',FILE='namelist.salsa')
    READ(11,NML=salsa)
    CLOSE(11)
 
  END SUBROUTINE define_salsa



  !-------------------------------------------------------------------------------
  !
  ! *****************************
  ! SUBROUTINE salsa_initialize
  ! *****************************
  !
  ! SALSA initializations. Modified and rewritten for more dynamic control 
  ! and LES implementation.
  !
  ! Juha Tonttila (FMI) 2014
  !
  !-------------------------------------------------------------------------------
  SUBROUTINE salsa_initialize()
 
    !
    !-------------------------------------------------------------------------------

    USE mo_kind, ONLY : dp
    USE mo_submctl, ONLY : reglim, nbin, nbin2, nbin3,     &
                               in1a,fn1a,in2a,fn2a,in2b,fn2b,  &
                               nbins,                          &
                               
                               epsv,sigma,csr_strat_wat,csr_strat_mix, &
                               csr_strat_ice,csr_conv,zbcr,             &

                               cfracn,zfracn,zfracn_cv,massacc,          &
                               iso4b,inob,inhb,iocb,ibcb,idub,issb

    USE mo_salsa_driver, ONLY : kproma,kbdim

    IMPLICIT NONE

    ! Remember to call 'define_salsa' for namelist paramers before calling this subroutine!

    ! --1) Set derived indices
    nbin2 = 4 
    nbin3 = nbin(2) - nbin2  

    in1a = 1     
    in2a = in1a + nbin(1)   
    in2b = in2a + nbin(2)  

    fn1a = in2a - 1
    fn2a = fn1a + nbin(2)
    fn2b = fn2a + nbin(2)

    nbins = fn2b          

    ! --2) Allocate arrays
    ALLOCATE(epsv(nbins), &
             sigma(fn2b),csr_strat_wat(fn2b),csr_strat_mix(fn2b),     &
             csr_strat_ice(fn2b),csr_conv(fn2b),zbcr(fn2b)                        )

    ALLOCATE(cfracn(fn2b),zfracn(fn2b),zfracn_cv(fn2b),massacc(nbins))

    massacc = 1._dp

    ! Tarvitaanko näitä??
    ALLOCATE(iso4b(fn2b),inob(fn2b),inhb(fn2b),iocb(fn2b),ibcb(fn2b),idub(fn2b),issb(fn2b))
    

    ! -- Aerosol tracers are allocated in *set_sizebins*
    ! -- Hydrometeor tracer in *set_cloudbins*

    ! --3) Call other initialization routines
    CALL set_sizebins()

    CALL set_cloudbins()

    CALL actcurve()

  END SUBROUTINE salsa_initialize
    
  
END MODULE mo_salsa_init
