
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

module mod_che_salsa

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_runparams , only : dtsec, dtche 
  !use mod_runparams , only : rcmtimer
  use mod_che_common
  use mod_che_indices
  use mod_che_species
  USE mo_submctl
  USE mo_salsa, ONLY : salsa
  USE mo_salsa_properties, ONLY : equilibration    
  USE class_componentIndex, ONLY : ComponentIndex, GetIndex, GetNcomp, IsUsed, ComponentIndexConstructor

  implicit none

  private
    ! grid points for SALSA
  integer(ik4)   :: kproma 
  INTEGER(ik4)  :: kbdim 
  INTEGER(ik4)  :: klev 
  INTEGER(ik4)  :: krow
  integer(ik4), parameter  :: nsalsp = 4
  integer(ik4), parameter  :: nmode = 7

  !
  logical, save :: salsa_firstcall
  ! -- Local hydrometeor properties (set up in aero initialize)
  TYPE(t_section), ALLOCATABLE, SAVE :: cloud(:,:,:) ! cloud properties
  TYPE(t_section), ALLOCATABLE, SAVE :: aero(:,:,:)  ! Aerosol properties
  TYPE(t_section), ALLOCATABLE, SAVE :: precp(:,:,:) ! Precipitation properties
  TYPE(t_section), ALLOCATABLE       :: actd(:,:,:)  ! activated droplet
  TYPE(t_section), ALLOCATABLE       :: aero_old(:,:,:),aero_ten(:,:,:) 
  
  REAL(rkx), allocatable :: ini_aero(:,:,:), ini_rh(:,:)
  TYPE(ComponentIndex) :: prtcl ! Contains "getIndex" which gives the index for a
                                ! aerosol component name, i.e. 1:SO4, 2:OC, 3:BC, 4:DU,
                                ! 5:SS, 6:NO, 7:NH, 8:H2O

  real(rkx) ::ssigmag(nmode),sdpg(nmode),sng(nmode)   

  ! -- Local gas compound tracers [# m-3]
  REAL(rkx), allocatable :: zgso4(:,:),   &
              zghno3(:,:),  &
              zgnh3(:,:),   &
              zgocnv(:,:),  &
              zgocsv(:,:),  &
              in_t(:,:),    & 
              in_p(:,:),    &
              in_rv(:,:),   &
              in_rs(:,:),   &
              in_w(:,:)    

  public :: run_salsa, init_salsa 

  contains 

  subroutine init_salsa()
    implicit none
    integer(ik4) ::i,j,k,n,nc,vc,str,end ,offset
    


    kproma = (jci2-jci1+1)*(ici2-ici1+1)
    kbdim = kproma
    klev = kz 
    krow = 1 ! ca sert pas


    CALL ComponentIndexConstructor(prtcl, nsalsp, 4, complist)


    allocate(in_t(1:kbdim,1:klev))    
    allocate(in_p(1:kbdim,1:klev))    
    allocate(in_rv(1:kbdim,1:klev))   
    allocate(in_rs(1:kbdim,1:klev))   
    allocate(in_w(1:kbdim,1:klev)) 
    allocate(zgso4(1:kbdim,1:klev))   
    allocate(zghno3(1:kbdim,1:klev))  
    allocate(zgnh3(1:kbdim,1:klev))   
    allocate(zgocnv(1:kbdim,1:klev))  
    allocate(zgocsv(1:kbdim,1:klev))  
     
    call salsa_initialize()

  ! FAB :This part is to inititialize particle size distribution 
  ! temporary !
  ! Think about how to handle this in aerosol ICBC later 
  ! Stdev of modes
    ssigmag = (/ 2.0, 1.5, 1.59, 1.59, 2.0, 2.0, 2.0 /)*1._rkx
  ! Mean diameter of the modes (m)
    sdpg = (/0.15, 1.0, 0.2, 0.2, 6.0, 6.0, 6.0/)*1.e-6_rkx
  ! Number concentration of modes (#/cm3)
    sng = (/100.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0/)*1.e6_rkx

    allocate(ini_aero(kbdim,klev,fn2b))   
    allocate(ini_rh(kbdim,klev))
    allocate(aero_old(kbdim,klev,fn2b))
    allocate(aero_ten(kbdim,klev,fn2b))

    CALL size_distribution(kproma, kbdim,  klev,   &
       sng, sdpg, ssigmag, ini_aero)


    offset = 0
       
!initialise les masses maintenant hypothese just sulfate 
       if (IsUsed(prtcl,'SO4')) THEN
          nc = GetIndex(prtcl,'SO4')
          vc = 1
! aero%volc(vc) should be in m3/m3 , convertir les traceurs de kg/kg à m3/m3
!ini_aero en #/m3 
          aero(1:kproma,1:klev,1:nbins)%volc(vc) = ini_aero(1:kproma,1:klev,1:nbins) &
                                                      * pi6*aero(1:kproma,1:klev,1:nbins)%dmid**3
! 
          aero_old(1:kproma,1:klev,1:nbins)%volc(vc) = aero(1:kproma,1:klev,1:nbins)%volc(vc)
! 
! need to initialize also tracers : to be refined with CHBC
         
          str = offset + (nc-1)*nbins+1
          end = offset + nc*nbins
          do k = 1 , kz
          n = 1
            do i = ici1 , ici2
              do j = jci1 , jci2
              !back from m3/m3 to kg/kg
                chemt(j,i,k,str:end) = aero(n,k,1:nbins)%volc(vc)/crhob3d(j,i,k)*rhosu
                n = n + 1
           end do
           end do
           end do
     end if
     
 !intialise number 
       aero(1:kproma,1:klev,1:nbins)%numc  = ini_aero(1:kproma,1:klev,1:nbins)
       aero_old(1:kproma,1:klev,1:nbins)%numc = aero_old(1:kproma,1:klev,1:nbins)%numc
       str = offset + (nsalsp)*nbins+1
       end = offset + (nsalsp+1)*nbins
       do k = 1 , kz
          n = 1
            do i = ici1 , ici2
              do j = jci1 , jci2
              !back from #/m3 to #/kg
                chemt(j,i,k,str:end) = aero(n,k,1:nbins)%numc/crhob3d(j,i,k)
                n = n + 1
            end do
           end do
        end do


    ini_rh(:,:) = 0.3_rkx


    salsa_firstcall = .true.



  end subroutine init_salsa  

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

!    USE mo_salsa_driver, ONLY : &
!         kbdim,kproma,klev,  &
!         aero


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
             aero(ii,jj,cc)%dmid = ( (aero(ii,jj,cc)%vhilim +aero(ii,jj,cc)%vlolim) /  &
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

          ! Initialize the wet diameter with the bin dry diameter to avoid
          ! numerical proplems later
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
  ! 
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
 !FAB    USE mo_salsa_driver, ONLY : kproma, kbdim, klev, &
 !                                cloud,precp,aero

                                      
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

    ! Rain/drizzle /bins
    ira = 1; fra = 7;
    nprc = fra

    ! ----------------------------------------
    ! Allocate cloud and precipitation arrays
    ! ----------------------------------------
    ALLOCATE(cloud(kbdim,klev,ncld), precp(kbdim,klev,nprc),actd(kbdim,klev,ncld))

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

!    USE  mo_salsa_driver, ONLY : &
!         aero

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

    !    epsv(in1a:fn2b) =
    !    4._dp*aa**3/(27._dp*bb*(log(slim))**2*vhilim(in1a:fn2b))
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
         nlcndh2oae,    & ! Condensation of water vapour on aerosols (FALSE ->equilibrium calc.)
         nlauto,        & ! Switch for autoconversion of cloud droplets todrizzle and rain
         nlactiv,       & ! Master witch for cloud droplet activation 
         nlactbase,     & ! Switch for parameterized cloud base activation
         nlactintst,    & ! Swithc for interstitial activation based on particlegrowth and host model S
         isdtyp,        & ! Type of initial size distribution: 0 - uniform; 1 -vertical profile, read from file
         reglim,        & ! Low/high diameter limits of the 2 aerosol sizeregimes (1d table with length 4)
         nbin,          & ! Number of bins used for each of the aerosol sizeregimes (1d table with length 2)
         dmincld,       & ! Minimum hydrometeor bin diameter: the lower limit ofthe smallest hydrometeor bin 
                          ! is set at or above. For now this is also the only
                          ! way to control the number of
                          ! hydrometeor bins as they are assumed fully parallel
                          ! with the aerosol bins.
                          ! Implementing different bin limits for hydromets than
                          ! aerosols would require
                          ! somewhat extensive modification of the microphysical
                          ! SALSA dynamics calculations
         nspec,         & ! Number of aerosol species used in the model
         listspec,      & ! List of strings specifying the names of the aerosol species that are active.
                          ! Must be an array of length 7, with empty strings for
                          ! unused stuff.
         volDistA,      & ! Initial relative contribution [0-1] of each speciesto particle volume in a-bins. Must be
                          ! an array of length 7, with zero for unused species.
         volDistB,      & ! Same as above but for b-bins
         nf2a,          & ! Number fraction of particles allocated to a-bins inregime 2. b-bins will get 1-nf2a

         sigmag,        & ! Stdev for the 7 initial lognormal modes
         dpg,           & ! Mean diameter for the 7 initial lognormal modes
         n                ! Number concentration for the 7 initial lognormalmodes

    
    OPEN(11,STATUS='old',FILE='namelist.salsa')
    READ(11,NML=salsa)
    CLOSE(11)
 
  END SUBROUTINE define_salsa

    SUBROUTINE size_distribution(kproma, kbdim, klev, &
       n, dpg, sigmag, naero)

    USE mo_submctl, ONLY :      &
         nreg,                      &
         pi6,                       &
         pi,                        &
         in1a,                      &
         fn2a,                      &
         fn2b

!    USE mo_salsa_driver, ONLY : aero ! This is needed for size bins spacings

    USE mo_kind, ONLY : dp

    IMPLICIT NONE

    INTEGER, PARAMETER :: nmod = 7

    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid points 
         kbdim,                     & ! dimension for arrays 
         klev                         ! number of vertical levels 

    REAL(dp), INTENT(IN) ::         &
         n(nmod)                  , & ! total concentration of a mode
         dpg(nmod)                , & ! geometric-mean diameter of a mode
         sigmag(nmod)                 ! standard deviation of a mode

    REAL(dp), INTENT(OUT) ::        &
         naero(kbdim,klev,fn2b)      ! number concentration  [#/m3]

    !-- local variables
    REAL(dp) ::                     &
         deltadp                      ! bin width [m]

    REAL(dp) :: d1,d2,delta_d,dmid

    INTEGER :: ii, jj, kk,ib

    naero = 0._dp

    DO jj = 1,klev    ! vertical grid
       DO ii = 1,kproma ! horizontal grid

          DO kk = in1a, fn2a
             naero(ii,jj,kk) = 0.0

             d1 = (aero(ii,jj,kk)%vlolim/pi6)**(1._dp/3._dp)
             d2 = (aero(ii,jj,kk)%vhilim/pi6)**(1._dp/3._dp)
             delta_d=(d2-d1)/10
             DO ib = 1,10
                d1=(aero(ii,jj,kk)%vlolim/pi6)**(1._dp/3._dp)+(ib-1.)*delta_d
                d2=d1+delta_d
                dmid=(d1+d2)/2
                deltadp = log(d2/d1)
                
                !deltadp = (aero(ii,jj,kk)%vhilim**(1._dp/3._dp)-aero(ii,jj,kk)%vlolim**(1._dp/3._dp))/   &
                !     pi6**(1._dp/3._dp)
                
                !-- size distribution
                !   ntot = total number, total area, or total volume concentration
                !   dpg = geometric-mean number, area, or volume diameter
                !   n(kk) = number, area, or volume concentration in a bin
!                naero(ii,jj,kk) = naero(ii,jj,kk)+sum(n*deltadp/                        &
!                     (sqrt(2._dp*pi)*log(sigmag))*                 &
!                     exp(-log(aero(ii,jj,kk)%dmid/dpg)**2/(2._dp*log(sigmag)**2)))
                  naero(ii,jj,kk) = naero(ii,jj,kk)+sum(n*deltadp/                        &
                     (sqrt(2._dp*pi)*log(sigmag))*                 &
                     exp(-log(dmid/dpg)**2/(2._dp*log(sigmag)**2)))
                

 !                write(*,*) naero(ii,jj,kk)*1e-6,d1,d2,dmid
             END DO
 !            write(*,*) dpg(1:2),n(1:2)*1e-6
 !            write(*,*) kk,naero(ii,jj,kk)*1e-6,aero(ii,jj,kk)%dmid, (aero(ii,jj,kk)%vlolim/pi6)**(1._dp/3._dp), &
 !                 (aero(ii,jj,kk)%vhilim/pi6)**(1._dp/3._dp)
 !            write(*,*) sum(naero(ii,jj,:))
 !            write(*,*) 'here'
 !            pause
          END DO

       END DO

    END DO

  END SUBROUTINE size_distribution

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

!    USE mo_salsa_driver, ONLY : kproma,kbdim

    IMPLICIT NONE

    ! Remember to call 'define_salsa' for namelist paramers before calling this
    ! subroutine!

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
             csr_strat_ice(fn2b),csr_conv(fn2b),zbcr(fn2b))

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

 !---------------------------------------------------------------
  ! SET_SALSA_RUNTIME
  ! Set logical switches according to the host model state and
  ! user-specified NAMELIST options. 
  !
  ! Juha Tonttila, FMI, 2014
  !
  SUBROUTINE set_SALSA_runtime(prunmode)
    USE mo_submctl, ONLY : nlcoag,                 &
                               nlcgaa,nlcgcc,nlcgpp,   &
                               nlcgca,nlcgpa,nlcgpc,   &
                               nlcnd,                  &
                               nlcndgas,               &
                               nlcndh2oae, nlcndh2ocl, &
                               nlauto,nlactiv,         &
                               nlactbase,nlactintst,   &
                               
                               lscoag,                 &
                               lscgaa,lscgcc,lscgpp,   &
                               lscgca,lscgpa,lscgpc,   &
                               lscnd,                  &
                               lscndgas,               &
                               lscndh2oae, lscndh2ocl, &
                               lsauto,lsactiv,         &
                               lsactbase,lsactintst

    IMPLICIT NONE

    INTEGER, INTENT(in) :: prunmode

    SELECT CASE(prunmode)

       CASE(1) ! Initialization

          lscoag      = .FALSE.
          lscnd       = nlcnd
          lscndgas    = nlcndgas
          lscndh2oae  = nlcndh2oae
          lscndh2ocl  = nlcndh2ocl
          lsauto      = .FALSE.
          lsactiv     = nlactiv
          lsactbase   = .FALSE.
          lsactintst  = .TRUE.   

       CASE(2)  ! Spinup period
          
          lscoag      = ( .FALSE. .AND. nlcoag   ) 
          lscnd       = ( .TRUE.  .AND. nlcnd    )
          lscndgas    = ( .TRUE.  .AND. nlcndgas )
          lscndh2oae  = ( .TRUE.  .AND. nlcndh2oae )
          lscndh2ocl  = ( .TRUE.  .AND. nlcndh2ocl )
          lsauto      = ( .FALSE. .AND. nlauto   )
          lsactiv     = ( .TRUE. .AND. nlactiv  )
          lsactbase   = ( .TRUE. .AND. nlactbase )
          lsactintst  = ( .TRUE. .AND. nlactintst )

       CASE(3)  ! Run
          
          lscoag      = nlcoag
          lscgaa      = nlcgaa
          lscgcc      = nlcgcc
          lscgpp      = nlcgpp
          lscgca      = nlcgca
          lscgpa      = nlcgpa
          lscgpc      = nlcgpc
          lscnd       = nlcnd
          lscndgas    = nlcndgas
          lscndh2oae  = nlcndh2oae
          lscndh2ocl  = nlcndh2ocl
          lsauto      = nlauto
          lsactiv     = nlactiv
          lsactbase   = nlactbase
          lsactintst  = nlactintst

    END SELECT
       
  END SUBROUTINE set_SALSA_runtime


!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine run_salsa()
    USE mo_submctl, only: ncld 

    implicit none
    
    integer(ik4) :: i,j,k,n,nc,vc,offset
    logical, parameter :: dbg = .FALSE.
    integer(rk4) :: pruntime
! set input


   do k = 1 , kz
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        in_t(n,k) = ctb3d(j,i,k)
        in_p(n,k) = cpb3d(j,i,k)
        in_rv(n,k)= cqxb3d(j,i,k,iqv)
        in_rs(n,k) =  cqsb3d(j,i,k) 
        in_w(n,k) = cwpx3d(j,i,k) 
! FAB délicat !!
!
! for the gas convert to #/m3

        zgso4(n,k) = d_zero
        zghno3(n,k)= 8.4682e-15_rkx*avog*1.e6_rkx
        zgnh3(n,k) = d_zero
        zgocnv(n,k)= d_zero
        zgocsv(n,k)= d_zero 

        n = n + 1
      end do
    end do
   end do

!                str = (nc-1)*ncld+1
!                end = nc*ncld
!                cloud(1,1,1:ncld)%volc(vc) =
!pa_vcloudp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhosu
!                vcloud_old(kk,ii,jj,str:end) = cloud(1,1,1:ncld)%volc(vc)

!                str = (nc-1)*nprc+1
!                end = nc*nprc
!                precp(1,1,1:nprc)%volc(vc) =
!pa_vprecpp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhosu
!                vprecp_old(kk,ii,jj,str:end) = precp(1,1,1:nprc)%volc(vc)




! aero%numc lorsque traceurs actives  aero(n,k,1:nbins)%numc devra etre converty de kg/kg à 
! #/m3 
!
! 
! aero%volc(vc) should be in m3/m3 , convertir les traceurs de kg/kg à m3/m3
!ini_aero en #/m3 
        

   if (salsa_firstcall) then
     pruntime = 1
     salsa_firstcall = .false.
   else
     pruntime = 3
   end if

! aero%volc(vc) should be in m3/m3 , convertir les traceurs de kg/kg à m3/m3
!ini_aero en #/m3 
        
     if (pruntime==1) call equilibration(kproma,kbdim,klev, ini_rh,in_t,aero,.TRUE.)


! Set aerosol bins volume concentration in volume mixing ratio.
! Very important :  since the mechanisms are hard coded in salsa
! the value of vc is fixed in function of a determined aerosol mixture 
! respecting the indexing : 
! 1:SO4, 2:OC, 3:BC, 4:DU, 5:SS, 6:NO, 7:NH, 8:H2O
! You can define a subset of these components as tracer effectively used,
! and this is defined through  prtcl. The GetIndex will adjust accordingly
! for stacking in the tracer table chib3d.    
     offset = 0

     if (pruntime >1 ) then 


       if (IsUsed(prtcl,'SO4')) THEN
            nc = GetIndex(prtcl,'SO4')
            vc = 1
            call trac2salsa('mass',vc,nc,offset,rhosu)
       end if
       if (IsUsed(prtcl,'NO3')) THEN
            nc = GetIndex(prtcl,'NO3')
            vc = 6 ! cf what I was saying before
            call trac2salsa('mass',vc,nc,offset,rhono)
       end if
       if (IsUsed(prtcl,'NH4')) THEN
            nc = GetIndex(prtcl,'NO3')
            vc = 7 ! cf what I was saying before
            call trac2salsa('mass',vc,nc,offset,rhonh)
       end if
       ! H2O is always is the mixture
       nc = GetIndex(prtcl,'H2O')
       vc = 8 ! cf what I was saying before
       call trac2salsa('mass',vc,nc,offset,rhoh2o)

      ! set number distribution  
       nc =nsalsp +1
       vc = 0 ! not used 
       call trac2salsa('numb',vc,nc,offset,d_zero)

     end if ! pruntime  


     call set_salsa_runtime(pruntime) 

!
     if(1==0) then
     call salsa(kproma, kbdim,  klev,   krow,                  &
                        in_p,   in_rv,  in_rs,  in_t, dt,   &
                        zgso4,  zgocnv, zgocsv, zghno3,        &
                        zgnh3,  aero,   cloud,  precp,         &
                        actd,   in_w,   dbg,   prtcl          )
     end if
! ladies and gentlemen it is time to calculate the tendencies ..

       if (IsUsed(prtcl,'SO4')) THEN
            nc = GetIndex(prtcl,'SO4')
            vc = 1
            call salsa2tractend('mass',vc,nc,offset,rhosu)
       end if
       if (IsUsed(prtcl,'NO3')) THEN
            nc = GetIndex(prtcl,'NO3')
            vc = 6 ! cf what I was saying before
            call salsa2tractend('mass',vc,nc,offset,rhono)
       end if
       if (IsUsed(prtcl,'NH4')) THEN
            nc = GetIndex(prtcl,'NO3')
            vc = 7 ! cf what I was saying before
            call salsa2tractend('mass',vc,nc,offset,rhonh)
       end if
       ! H2O is always is the mixture
       nc = GetIndex(prtcl,'H2O')
       vc = 8 ! cf what I was saying before
       call salsa2tractend('mass',vc,nc,offset,rhoh2o)

      ! set number distribution   
       nc = nsalsp + 1
       vc = 0
       call salsa2tractend('numb',vc,nc,offset,d_zero)

end subroutine run_salsa

  subroutine trac2salsa(typ,vc,nc,offset,rhop)
     implicit none
  integer(ik4) , intent(in) :: vc,nc,offset
  character(len=4),intent(in) :: typ
  real(rkx) , intent(in) :: rhop
  integer(ik4) :: i,j,k,n,str,end

   if (typ  == "mass") then 
   str = offset + (nc-1)*nbins+1
   end = offset + nc*nbins       
   do k = 1 , kz
      n = 1
       do i = ici1 , ici2
         do j = jci1 , jci2   
           aero(n,k,1:nbins)%volc(vc) = chib3d(j,i,k,str:end) * crhob3d(j,i,k) /rhop
           aero_old(n,k,1:nbins)%volc(vc) = aero(n,k,1:nbins)%volc(vc)
           n = n + 1
      end do
    end do
   end do
   else if( typ == "numb") then
   str = offset + (nc-1)*nbins+1
   end = offset + nc*nbins
   do k = 1 , kz
      n = 1
       do i = ici1 , ici2
         do j = jci1 , jci2
           aero(n,k,1:nbins)%numc = chib3d(j,i,k,str:end) * crhob3d(j,i,k) 
           aero_old(n,k,1:nbins)%numc = aero(n,k,1:nbins)%numc
           n = n + 1
      end do
    end do
   end do
   end if
  end subroutine trac2salsa



  subroutine salsa2tractend(typ,vc,nc,offset,rhop)
     implicit none
  integer(ik4) , intent(in) :: vc,nc,offset
  character(len=4),intent(in) :: typ
  real(rkx) , intent(in) :: rhop
  integer(ik4) :: i,j,k,n,str,end


   if (typ  == "mass") then
   str = offset + (nc-1)*nbins+1
   end = offset + nc*nbins
   do k = 1 , kz
      n = 1
       do i = ici1 , ici2
         do j = jci1 , jci2
           aero_ten(n,k,1:nbins)%volc(vc) = (aero(n,k,1:nbins)%volc(vc) - &
                                             aero_old(n,k,1:nbins)%volc(vc)) / dt


           !back from m3/m3 to kg/kg
           chiten(j,i,k,str:end) = chiten(j,i,k,str:end) +  &
                                   aero_ten(n,k,1:nbins)%volc(vc)/crhob3d(j,i,k) *rhop

           n = n + 1
      end do
    end do
   end do
   else if( typ == "numb") then
   str = offset + (nc-1)*nbins+1
   end = offset + nc*nbins
   do k = 1 , kz
      n = 1
       do i = ici1 , ici2
         do j = jci1 , jci2
           aero_ten(n,k,1:nbins)%numc = (aero(n,k,1:nbins)%numc - &
                                         aero_old(n,k,1:nbins)%numc) / dt
           !back from #/m3 to #/kg
           chiten(j,i,k,str:end) = chiten(j,i,k,str:end) +  &
                                   aero_ten(n,k,1:nbins)%volc(vc) * crhob3d(j,i,k)

           n = n + 1
      end do
    end do
   end do
   end if
  end subroutine salsa2tractend




end module mod_che_salsa
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
