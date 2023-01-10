MODULE mo_salsa_driver
USE mo_kind, ONLY : dp
USE mo_submctl, ONLY : t_section
IMPLICIT NONE

!---------------------------------------------------------------
!
! MOD_SALSA_DRIVER:
! Contains the primary SALSA input/output variables as well as 
! subroutines used to call the main SALSA routine.
!
! Juha Tonttila, FMI, 2014
!
!---------------------------------------------------------------


  ! JT: Variables from SALSA
  ! --------------------------------------------
  ! grid points for SALSA
  INTEGER, PARAMETER :: kproma = 1
  INTEGER, PARAMETER :: kbdim = 1
  INTEGER, PARAMETER :: klev = 1
  INTEGER, PARAMETER :: krow = 1

  REAL(dp), PARAMETER :: init_rh(kbdim,klev) = 0.3_dp

  ! -- Local hydrometeor properties (set up in aero initialize)
  TYPE(t_section), ALLOCATABLE, SAVE :: cloud(:,:,:) ! cloud properties
  TYPE(t_section), ALLOCATABLE, SAVE :: aero(:,:,:)  ! Aerosol properties
  TYPE(t_section), ALLOCATABLE, SAVE :: precp(:,:,:) ! Precipitation properties

  ! -- Local gas compound tracers [# m-3]
  REAL(dp) :: zgso4(kbdim,klev),   &
              zghno3(kbdim,klev),  &
              zgnh3(kbdim,klev),   &
              zgocnv(kbdim,klev),  &
              zgocsv(kbdim,klev)

   ! --------------------------------------------


  CONTAINS

  !
  !----------------------------------------------------
  ! RUN_SALSA:
  ! Performs necessary unit and dimension conversion between
  ! the host model and SALSA module, and calls the main SALSA
  ! routine
  !
  ! Partially adobted form the original SALSA boxmodel version.
  !
  ! Now takes masses in as kg/kg from LES!! Converted to m3/m3 for SALSA
  !
  ! Juha Tonttila, FMI, 2014
  !
  SUBROUTINE run_SALSA(pnx, pny, pnz, n4, press, tk, rv, rt, rs, wp, pdn,   &
                       pa_naerop,  pa_naerot,  pa_vaerop,  pa_vaerot,   &
                       pa_ncloudp, pa_ncloudt, pa_vcloudp, pa_vcloudt,  &
                       pa_nprecpp, pa_nprecpt, pa_vprecpp, pa_vprecpt,  &
                       pa_nactd,   pa_vactd,   pa_gaerop,  pa_gaerot,   &
                       pa_Radry,   pa_Rcdry,   pa_Rpdry,   pa_Rawet,    &
                       pa_Rcwet,   pa_Rpwet,   pa_rhop,    prunmode,    &
                       prtcl,      tstep,      dbg2                     )

    USE mo_submctl, ONLY : nbins,ncld,nprc,pi6,          &
                               rhowa, rhosu, rhobc, rhooc,   &
                               rhono, rhonh, rhoss, rhodu, fn2a
    USE mo_salsa, ONLY : salsa
    USE mo_salsa_properties, ONLY : equilibration
    USE class_componentIndex, ONLY : ComponentIndex, GetIndex, GetNcomp, IsUsed
    IMPLICIT NONE

    INTEGER, INTENT(in) :: pnx,pny,pnz,n4
    REAL(dp), INTENT(in)    :: tstep                            ! Model timestep length

    REAL(dp), INTENT(in)    :: press(pnz,pnx,pny), &            ! Pressure (Pa)
                               tk(pnz,pnx,pny),    &            ! Temperature (K)
                               rv(pnz,pnx,pny),    &            ! Water vapor mixing ratio
                               rs(pnz,pnx,pny),    &            ! Water vapour saturation mixing ratio
                               wp(pnz,pnx,pny)                  ! Vertical velocity (m s-1)
    
    REAL(dp), INTENT(in)    :: pdn(pnz,pnx,pny)             ! Air density (for normalizing concentrations)

    REAL(dp), INTENT(in)    :: pa_naerop(pnz,pnx,pny,nbins),   & ! aerosol number concentration (# kg-1)
                               pa_vaerop(pnz,pnx,pny,n4*nbins), & ! aerosol volume concentration (kg kg-1)
                               pa_ncloudp(pnz,pnx,pny,ncld),   & ! Cloud droplet number concentration (# kg-1)
                               pa_vcloudp(pnz,pnx,pny,n4*ncld), & ! Cloud droplet volume concentration (kg kg-1)
                               pa_nprecpp(pnz,pnx,pny,nprc),   & ! Rain drop number concentration (# kg-1)
                               pa_vprecpp(pnz,pnx,pny,n4*nprc)    ! Rain drop volume concentration (kg kg-1)

    REAL(dp), INTENT(in)    :: pa_gaerop(pnz,pnx,pny,5)         ! Gaseous tracers [# kg]

    INTEGER, INTENT(in) :: prunmode                         ! 1: Initialization call
                                                            ! 2: Spinup period call
                                                            ! 3: Regular runtime call'

    LOGICAL, INTENT(in) :: dbg2

    TYPE(ComponentIndex), INTENT(in) :: prtcl ! Object containing the indices of different aerosol components for mass arrays

    REAL(dp), INTENT(inout)   :: pa_naerot(pnz,pnx,pny,nbins),   & ! Aerosol number tendency
                                 pa_vaerot(pnz,pnx,pny,n4*nbins), & ! Aerosol volume tendency
                                 pa_ncloudt(pnz,pnx,pny,ncld),   & ! Cloud droplet number tendency
                                 pa_vcloudt(pnz,pnx,pny,n4*ncld), & ! Cloud droplet volume tendency
                                 pa_nprecpt(pnz,pnx,pny,nprc),   & ! Rain drop number tendency
                                 pa_vprecpt(pnz,pnx,pny,n4*nprc)    ! Rain drop volume tendency

    REAL(dp), INTENT(inout)   :: pa_gaerot(pnz,pnx,pny,5)         ! Gaseous tracer tendency
    REAL(dp), INTENT(inout)   :: rt(pnz,pnx,pny)                  ! Water vapour tendency

    REAL(dp), INTENT(in)   :: pa_Radry(pnz,pnx,pny,nbins),   & ! Aerosol dry particle radius
                              pa_Rcdry(pnz,pnx,pny,ncld),    & ! Cloud dry radius
                              pa_Rpdry(pnz,pnx,pny,nprc),    & ! Rain dry radius
                              pa_Rawet(pnz,pnx,pny,nbins),   & ! Aerosol wet radius
                              pa_Rcwet(pnz,pnx,pny,ncld),    & ! Cloud wet radius
                              pa_Rpwet(pnz,pnx,pny,nprc),    & ! Rain drop wet radius
                              pa_rhop(pnz,pnx,pny,nbins)       ! Aerosol density (kg/m3)

    REAL(dp), INTENT(out)   :: pa_vactd(pnz,pnx,pny,n4*ncld) ! Volume concentrations of newly activated droplets for calculating the 
                                                         ! actual tendency due to new droplet formation.
    REAL(dp), INTENT(out)   :: pa_nactd(pnz,pnx,pny,ncld)   ! Same for number concentration

    TYPE(t_section) :: actd(kbdim,klev,ncld) ! Activated droplets - for interfacing with SALSA

    REAL(dp) :: vaero_old(pnz,pnx,pny,n4*nbins) ! Nääki vois muuttaa TYPE-muotoon
    REAL(dp) :: vcloud_old(pnz,pnx,pny,n4*ncld)
    REAL(dp) :: vprecp_old(pnz,pnx,pny,n4*nprc)
    REAL(dp) :: naero_old(pnz,pnx,pny,nbins)
    REAL(dp) :: ncloud_old(pnz,pnx,pny,ncld)
    REAL(dp) :: nprecp_old(pnz,pnx,pny,nprc)

    LOGICAL :: dbg3

    INTEGER :: jj,ii,kk,ss,str,end, nc,vc, cc
    REAL(dp) :: in_p(kbdim,klev), in_t(kbdim,klev), in_rv(kbdim,klev), in_rs(kbdim,klev), in_w(kbdim,klev) 


    vaero_old = 0._dp;  naero_old = 0._dp
    vcloud_old = 0._dp; ncloud_old = 0._dp
    vprecp_old = 0._dp; nprecp_old = 0._dp

    ! NÄiden "luokkien" alustamiseen tarttis jonku kätevämmän systeemin
    actd(1:kproma,:,:)%numc = 0._dp
    aero(1:kproma,:,:)%numc = 0._dp
    cloud(1:kproma,:,:)%numc = 0._dp
    precp(1:kproma,:,:)%numc = 0._dp
    DO ss = 1,8 !GetNcomp(prtcl)+1  !!!! FIXED, should be 1,8

       actd(1:kproma,:,:)%volc(ss) = 0._dp
       aero(1:kproma,:,:)%volc(ss) = 0._dp
       cloud(1:kproma,:,:)%volc(ss) = 0._dp
       precp(1:kproma,:,:)%volc(ss) = 0._dp

    END DO

    ! Set the SALSA runtime config (saisiko hoidettua tehokkaammin?)
    CALL set_salsa_runtime(prunmode)
    
    ! Convert input concentrations for SALSA into #/m3 or m3/m3 instead of kg/kg (multiplied by pdn/divided by subtance density)
    DO jj = 3,pny-2
       DO ii = 3,pnx-2 
          DO kk = pnz-1,2,-1
 
             IF ( ANY(pa_vaerot(kk,ii,jj,:) /= pa_vaerot(kk,ii,jj,:)) )THEN
                WRITE(*,*) 'NAN1hop'
                WRITE(*,*) kk,ii,jj,pa_vaerot(kk,ii,jj,:)
             END IF

             IF ( ANY(pa_vaerot(kk,ii,jj,:) /= pa_vaerot(kk,ii,jj,:)) ) STOP
             IF ( ANY(pa_naerop(kk,ii,jj,:) /= pa_naerop(kk,ii,jj,:)) ) WRITE(*,*) 'NAN2'
             IF ( ANY(pa_naerop(kk,ii,jj,:) /= pa_naerop(kk,ii,jj,:)) ) STOP
             IF ( ANY(pa_vcloudt(kk,ii,jj,:) /= pa_vcloudt(kk,ii,jj,:)) ) WRITE(*,*) 'NAN3'
             IF ( ANY(pa_vcloudt(kk,ii,jj,:) /= pa_vcloudt(kk,ii,jj,:)) ) STOP
             IF ( ANY(pa_vcloudp(kk,ii,jj,:) /= pa_vcloudp(kk,ii,jj,:)) ) WRITE(*,*) 'NAN4'
             IF ( ANY(pa_vcloudp(kk,ii,jj,:) /= pa_vcloudp(kk,ii,jj,:)) ) STOP
             IF ( ANY(pa_vprecpt(kk,ii,jj,:) /= pa_vprecpt(kk,ii,jj,:)) ) WRITE(*,*) 'NAN5'
             IF ( ANY(pa_vprecpt(kk,ii,jj,:) /= pa_vprecpt(kk,ii,jj,:)) ) STOP
             IF ( ANY(pa_vprecpp(kk,ii,jj,:) /= pa_vprecpp(kk,ii,jj,:)) ) WRITE(*,*) 'NAN6'
             IF ( ANY(pa_vprecpp(kk,ii,jj,:) /= pa_vprecpp(kk,ii,jj,:)) ) STOP

            ! Set inputs
             in_p(1,1) = press(kk,ii,jj)
             in_t(1,1) = tk(kk,ii,jj)
             in_rv(1,1) = rv(kk,ii,jj)
             in_rs(1,1) = rs(kk,ii,jj)
             in_w(1,1) = wp(kk,ii,jj)

             ! Set volume concentrations
             IF (IsUsed(prtcl,'SO4')) THEN
                nc = GetIndex(prtcl,'SO4')
                vc = 1
                str = (nc-1)*nbins+1
                end = nc*nbins
                aero(1,1,1:nbins)%volc(vc) = pa_vaerop(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhosu
                vaero_old(kk,ii,jj,str:end) = aero(1,1,1:nbins)%volc(vc)

                str = (nc-1)*ncld+1
                end = nc*ncld
                cloud(1,1,1:ncld)%volc(vc) = pa_vcloudp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhosu
                vcloud_old(kk,ii,jj,str:end) = cloud(1,1,1:ncld)%volc(vc)

                str = (nc-1)*nprc+1
                end = nc*nprc
                precp(1,1,1:nprc)%volc(vc) = pa_vprecpp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhosu
                vprecp_old(kk,ii,jj,str:end) = precp(1,1,1:nprc)%volc(vc)
             END IF
             
             IF (IsUsed(prtcl,'OC')) THEN
                nc = GetIndex(prtcl,'OC')
                vc = 2
                str = (nc-1)*nbins+1
                end = nc*nbins
                aero(1,1,1:nbins)%volc(vc) = pa_vaerop(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhooc
                vaero_old(kk,ii,jj,str:end) = aero(1,1,1:nbins)%volc(vc)

                str = (nc-1)*ncld+1
                end = nc*ncld
                cloud(1,1,1:ncld)%volc(vc) = pa_vcloudp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhooc
                vcloud_old(kk,ii,jj,str:end) = cloud(1,1,1:ncld)%volc(vc)

                str = (nc-1)*nprc+1
                end = nc*nprc
                precp(1,1,1:nprc)%volc(vc) = pa_vprecpp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhooc
                vprecp_old(kk,ii,jj,str:end) = precp(1,1,1:nprc)%volc(vc)
             END IF

             IF (IsUsed(prtcl,'BC')) THEN
                nc = GetIndex(prtcl,'BC')
                vc = 3
                str = (nc-1)*nbins+1
                end = nc*nbins
                aero(1,1,1:nbins)%volc(vc) = pa_vaerop(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhobc
                vaero_old(kk,ii,jj,str:end) = aero(1,1,1:nbins)%volc(vc)

                str = (nc-1)*ncld+1
                end = nc*ncld
                cloud(1,1,1:ncld)%volc(vc) = pa_vcloudp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhobc
                vcloud_old(kk,ii,jj,str:end) = cloud(1,1,1:ncld)%volc(vc)

                str = (nc-1)*nprc+1
                end = nc*nprc
                precp(1,1,1:nprc)%volc(vc) = pa_vprecpp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhobc
                vprecp_old(kk,ii,jj,str:end) = precp(1,1,1:nprc)%volc(vc)
             END IF

             IF (IsUsed(prtcl,'DU')) THEN
                nc = GetIndex(prtcl,'DU')
                vc = 4
                str = (nc-1)*nbins+1
                end = nc*nbins
                aero(1,1,1:nbins)%volc(vc) = pa_vaerop(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhodu
                vaero_old(kk,ii,jj,str:end) = aero(1,1,1:nbins)%volc(vc)

                str = (nc-1)*ncld+1
                end = nc*ncld
                cloud(1,1,1:ncld)%volc(vc) = pa_vcloudp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhodu
                vcloud_old(kk,ii,jj,str:end) = cloud(1,1,1:ncld)%volc(vc)

                str = (nc-1)*nprc+1
                end = nc*nprc
                precp(1,1,1:nprc)%volc(vc) = pa_vprecpp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhodu
                vprecp_old(kk,ii,jj,str:end) = precp(1,1,1:nprc)%volc(vc)
             END IF

             IF (IsUsed(prtcl,'SS')) THEN
                nc = GetIndex(prtcl,'SS')
                vc = 5
                str = (nc-1)*nbins+1
                end = nc*nbins
                aero(1,1,1:nbins)%volc(vc) = pa_vaerop(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhoss
                vaero_old(kk,ii,jj,str:end) = aero(1,1,1:nbins)%volc(vc)

                str = (nc-1)*ncld+1
                end = nc*ncld
                cloud(1,1,1:ncld)%volc(vc) = pa_vcloudp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhoss
                vcloud_old(kk,ii,jj,str:end) = cloud(1,1,1:ncld)%volc(vc)

                str = (nc-1)*nprc+1
                end = nc*nprc
                precp(1,1,1:nprc)%volc(vc) = pa_vprecpp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhoss
                vprecp_old(kk,ii,jj,str:end) = precp(1,1,1:nprc)%volc(vc)
             END IF

             IF (IsUsed(prtcl,'NO')) THEN
                nc = GetIndex(prtcl,'NO')
                vc = 6
                str = (nc-1)*nbins+1
                end = nc*nbins
                aero(1,1,1:nbins)%volc(vc) = pa_vaerop(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhono
                vaero_old(kk,ii,jj,str:end) = aero(1,1,1:nbins)%volc(vc)

                str = (nc-1)*ncld+1
                end = nc*ncld
                cloud(1,1,1:ncld)%volc(vc) = pa_vcloudp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhono
                vcloud_old(kk,ii,jj,str:end) = cloud(1,1,1:ncld)%volc(vc)

                str = (nc-1)*nprc+1
                end = nc*nprc
                precp(1,1,1:nprc)%volc(vc) = pa_vprecpp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhono
                vprecp_old(kk,ii,jj,str:end) = precp(1,1,1:nprc)%volc(vc)
             END IF

             IF (IsUsed(prtcl,'NH')) THEN
                nc = GetIndex(prtcl,'NH')
                vc = 7
                str = (nc-1)*nbins+1
                end = nc*nbins
                aero(1,1,1:nbins)%volc(vc) = pa_vaerop(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhonh
                vaero_old(kk,ii,jj,str:end) = aero(1,1,1:nbins)%volc(vc)

                str = (nc-1)*ncld+1
                end = nc*ncld
                cloud(1,1,1:ncld)%volc(vc) = pa_vcloudp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhonh
                vcloud_old(kk,ii,jj,str:end) = cloud(1,1,1:ncld)%volc(vc)

                str = (nc-1)*nprc+1
                end = nc*nprc
                precp(1,1,1:nprc)%volc(vc) = pa_vprecpp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhonh
                vprecp_old(kk,ii,jj,str:end) = precp(1,1,1:nprc)%volc(vc)
             END IF

             ! Water (always used)
             ! -----------------------------
             nc = GetIndex(prtcl,'H2O')
             vc = 8
             str = (nc-1)*nbins+1
             end = nc*nbins
             aero(1,1,1:nbins)%volc(vc) = pa_vaerop(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhowa
             vaero_old(kk,ii,jj,str:end) = aero(1,1,1:nbins)%volc(vc)

             str = (nc-1)*ncld+1
             end = nc*ncld
             cloud(1,1,1:ncld)%volc(vc) = pa_vcloudp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhowa
             vcloud_old(kk,ii,jj,str:end) = cloud(1,1,1:ncld)%volc(vc)

             str = (nc-1)*nprc+1
             end = nc*nprc
             precp(1,1,1:nprc)%volc(vc) = pa_vprecpp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rhowa
             vprecp_old(kk,ii,jj,str:end) = precp(1,1,1:nprc)%volc(vc)
             ! -----------------------------

             ! Number concentrations and particle sizes
             aero(1,1,1:nbins)%numc = pa_naerop(kk,ii,jj,1:nbins)*pdn(kk,ii,jj)
             naero_old(kk,ii,jj,1:nbins) = aero(1,1,1:nbins)%numc
             aero(1,1,1:nbins)%dwet = pa_Rawet(kk,ii,jj,1:nbins)*2._dp
             aero(1,1,1:nbins)%core = pi6*(pa_Radry(kk,ii,jj,1:nbins)*2._dp)**3._dp

             cloud(1,1,1:ncld)%numc = pa_ncloudp(kk,ii,jj,1:ncld)*pdn(kk,ii,jj)
             ncloud_old(kk,ii,jj,1:ncld) = cloud(1,1,1:ncld)%numc
             cloud(1,1,1:ncld)%dwet = pa_Rcwet(kk,ii,jj,1:ncld)*2._dp
             cloud(1,1,1:ncld)%core = pi6*(pa_Rcdry(kk,ii,jj,1:ncld)*2._dp)**3._dp

             precp(1,1,1:nprc)%numc = pa_nprecpp(kk,ii,jj,1:nprc)*pdn(kk,ii,jj)
             nprecp_old(kk,ii,jj,1:nprc) = precp(1,1,1:nprc)%numc
             precp(1,1,1:nprc)%dwet = pa_Rpwet(kk,ii,jj,1:nprc)*2._dp
             precp(1,1,1:nprc)%core = pi6*(pa_Rpdry(kk,ii,jj,1:nprc)*2._dp)**3._dp


             ! If this is an initialization call, calculate the equilibrium particle 
             ! size at 30 %. SIIRRÄ JOHONKIN FIKSUMPAAN PAIKKAAN
             If (prunmode == 1) CALL equilibration(kproma,kbdim,klev,   &
                                                    init_rh,in_t,aero,.TRUE.)

             ! Convert to #/m3
             zgso4(1,1) = pa_gaerop(kk,ii,jj,1)*pdn(kk,ii,jj)
             zghno3(1,1) = pa_gaerop(kk,ii,jj,2)*pdn(kk,ii,jj)
             zgnh3(1,1) = pa_gaerop(kk,ii,jj,3)*pdn(kk,ii,jj)
             zgocnv(1,1) = pa_gaerop(kk,ii,jj,4)*pdn(kk,ii,jj)
             zgocsv(1,1) = pa_gaerop(kk,ii,jj,5)*pdn(kk,ii,jj)

             ! ***************************************!
             !                Run SALSA               !
             ! ***************************************! 
             CALL salsa(kproma, kbdim,  klev,   krow,          &
                        in_p,   in_rv,  in_rs,  in_t, tstep,   &
                        zgso4,  zgocnv, zgocsv, zghno3,        &
                        zgnh3,  aero,   cloud,  precp,         &
                        actd,   in_w,   dbg3,   prtcl          )


             !IF ( SUM(aero(1,1,1:fn2a)%numc) < 1.e9_dp ) THEN
             !   WRITE(*,*) 'ahhhh'
             !   WRITE(*,*) aero(1,1,1:fn2a)
             !   WRITE(*,*) naero_old(kk,ii,jj,1:fn2a)
             !END IF

             ! Calculate tendencies (convert back to #/kg or kg/kg)
             pa_naerot(kk,ii,jj,1:nbins) = pa_naerot(kk,ii,jj,1:nbins) + &
                  ( aero(1,1,1:nbins)%numc - naero_old(kk,ii,jj,1:nbins) )/pdn(kk,ii,jj)/tstep
             pa_ncloudt(kk,ii,jj,1:ncld) = pa_ncloudt(kk,ii,jj,1:ncld) + &
                  ( cloud(1,1,1:ncld)%numc - ncloud_old(kk,ii,jj,1:ncld) )/pdn(kk,ii,jj)/tstep
             pa_nprecpt(kk,ii,jj,1:nprc) = pa_nprecpt(kk,ii,jj,1:nprc) + &
                  ( precp(1,1,1:nprc)%numc - nprecp_old(kk,ii,jj,1:nprc) )/pdn(kk,ii,jj)/tstep
             ! Activated droplets
             pa_nactd(kk,ii,jj,1:ncld) = actd(1,1,1:ncld)%numc/pdn(kk,ii,jj)

             IF (IsUsed(prtcl,'SO4')) THEN
                nc = GetIndex(prtcl,'SO4')
                vc = 1
                ! Aerosol bins
                str = (nc-1)*nbins+1
                end = nc*nbins
                pa_vaerot(kk,ii,jj,str:end) = pa_vaerot(kk,ii,jj,str:end) + &
                     ( aero(1,1,1:nbins)%volc(vc) - vaero_old(kk,ii,jj,str:end) )*rhosu/pdn(kk,ii,jj)/tstep
                ! Hydrometeor bins
                str = (nc-1)*ncld+1
                end = nc*ncld
                pa_vcloudt(kk,ii,jj,str:end) = pa_vcloudt(kk,ii,jj,str:end) + &
                     ( cloud(1,1,1:ncld)%volc(vc) - vcloud_old(kk,ii,jj,str:end) )*rhosu/pdn(kk,ii,jj)/tstep
                ! Rain drops
                str = (nc-1)*nprc+1
                end = nc*nprc
                pa_vprecpt(kk,ii,jj,str:end) = pa_vprecpt(kk,ii,jj,str:end) + &
                     ( precp(1,1,1:nprc)%volc(vc) - vprecp_old(kk,ii,jj,str:end) )*rhosu/pdn(kk,ii,jj)/tstep
                ! Activated droplets
                str = (nc-1)*ncld+1
                end = nc*ncld
                pa_vactd(kk,ii,jj,str:end) = actd(1,1,1:ncld)%volc(vc)*rhosu/pdn(kk,ii,jj)
             END IF

             IF (IsUsed(prtcl,'OC')) THEN
                nc = GetIndex(prtcl,'OC')
                vc = 2
                ! Aerosol bins
                str = (nc-1)*nbins+1
                end = nc*nbins
                pa_vaerot(kk,ii,jj,str:end) = pa_vaerot(kk,ii,jj,str:end) + &
                     ( aero(1,1,1:nbins)%volc(vc) - vaero_old(kk,ii,jj,str:end) )*rhooc/pdn(kk,ii,jj)/tstep
                ! Hydrometeor bins
                str = (nc-1)*ncld+1
                end = nc*ncld
                pa_vcloudt(kk,ii,jj,str:end) = pa_vcloudt(kk,ii,jj,str:end) + &
                     ( cloud(1,1,1:ncld)%volc(vc) - vcloud_old(kk,ii,jj,str:end) )*rhooc/pdn(kk,ii,jj)/tstep
                ! Rain drops
                str = (nc-1)*nprc+1
                end = nc*nprc
                pa_vprecpt(kk,ii,jj,str:end) = pa_vprecpt(kk,ii,jj,str:end) + &
                     ( precp(1,1,1:nprc)%volc(vc) - vprecp_old(kk,ii,jj,str:end) )*rhooc/pdn(kk,ii,jj)/tstep
                ! Activated droplets
                str = (nc-1)*ncld+1
                end = nc*ncld
                pa_vactd(kk,ii,jj,str:end) = actd(1,1,1:ncld)%volc(vc)*rhooc/pdn(kk,ii,jj)
             END IF

             IF (IsUsed(prtcl,'BC')) THEN
                nc = GetIndex(prtcl,'BC')
                vc = 3
                ! Aerosol bins
                str = (nc-1)*nbins+1
                end = nc*nbins
                pa_vaerot(kk,ii,jj,str:end) = pa_vaerot(kk,ii,jj,str:end) + &
                     ( aero(1,1,1:nbins)%volc(vc) - vaero_old(kk,ii,jj,str:end) )*rhobc/pdn(kk,ii,jj)/tstep
                ! Hydrometeor bins
                str = (nc-1)*ncld+1
                end = nc*ncld
                pa_vcloudt(kk,ii,jj,str:end) = pa_vcloudt(kk,ii,jj,str:end) + &
                     ( cloud(1,1,1:ncld)%volc(vc) - vcloud_old(kk,ii,jj,str:end) )*rhobc/pdn(kk,ii,jj)/tstep
                ! Rain drops
                str = (nc-1)*nprc+1
                end = nc*nprc
                pa_vprecpt(kk,ii,jj,str:end) = pa_vprecpt(kk,ii,jj,str:end) + &
                     ( precp(1,1,1:nprc)%volc(vc) - vprecp_old(kk,ii,jj,str:end) )*rhobc/pdn(kk,ii,jj)/tstep
                ! Activated droplets
                str = (nc-1)*ncld+1
                end = nc*ncld
                pa_vactd(kk,ii,jj,str:end) = actd(1,1,1:ncld)%volc(vc)*rhobc/pdn(kk,ii,jj)
             END IF

             IF (IsUsed(prtcl,'DU')) THEN
                nc = GetIndex(prtcl,'DU')
                vc = 4
                ! Aerosol bins
                str = (nc-1)*nbins+1
                end = nc*nbins
                pa_vaerot(kk,ii,jj,str:end) = pa_vaerot(kk,ii,jj,str:end) + &
                     ( aero(1,1,1:nbins)%volc(vc) - vaero_old(kk,ii,jj,str:end) )*rhodu/pdn(kk,ii,jj)/tstep
                ! Hydrometeor bins
                str = (nc-1)*ncld+1
                end = nc*ncld
                pa_vcloudt(kk,ii,jj,str:end) = pa_vcloudt(kk,ii,jj,str:end) + &
                     ( cloud(1,1,1:ncld)%volc(vc) - vcloud_old(kk,ii,jj,str:end) )*rhodu/pdn(kk,ii,jj)/tstep
                ! Rain drops
                str = (nc-1)*nprc+1
                end = nc*nprc
                pa_vprecpt(kk,ii,jj,str:end) = pa_vprecpt(kk,ii,jj,str:end) + &
                     ( precp(1,1,1:nprc)%volc(vc) - vprecp_old(kk,ii,jj,str:end) )*rhodu/pdn(kk,ii,jj)/tstep
                ! Activated droplets
                str = (nc-1)*ncld+1
                end = nc*ncld
                pa_vactd(kk,ii,jj,str:end) = actd(1,1,1:ncld)%volc(vc)*rhodu/pdn(kk,ii,jj)
             END IF

             IF (IsUsed(prtcl,'SS')) THEN
                nc = GetIndex(prtcl,'SS')
                vc = 5
                ! Aerosol bins
                str = (nc-1)*nbins+1
                end = nc*nbins
                pa_vaerot(kk,ii,jj,str:end) = pa_vaerot(kk,ii,jj,str:end) + &
                     ( aero(1,1,1:nbins)%volc(vc) - vaero_old(kk,ii,jj,str:end) )*rhoss/pdn(kk,ii,jj)/tstep
                ! Hydrometeor bins
                str = (nc-1)*ncld+1
                end = nc*ncld
                pa_vcloudt(kk,ii,jj,str:end) = pa_vcloudt(kk,ii,jj,str:end) + &
                     ( cloud(1,1,1:ncld)%volc(vc) - vcloud_old(kk,ii,jj,str:end) )*rhoss/pdn(kk,ii,jj)/tstep
                ! Rain drops
                str = (nc-1)*nprc+1
                end = nc*nprc
                pa_vprecpt(kk,ii,jj,str:end) = pa_vprecpt(kk,ii,jj,str:end) + &
                     ( precp(1,1,1:nprc)%volc(vc) - vprecp_old(kk,ii,jj,str:end) )*rhoss/pdn(kk,ii,jj)/tstep
                ! Activated droplets
                str = (nc-1)*ncld+1
                end = nc*ncld
                pa_vactd(kk,ii,jj,str:end) = actd(1,1,1:ncld)%volc(vc)*rhoss/pdn(kk,ii,jj)
             END IF

             IF (IsUsed(prtcl,'NO')) THEN
                nc = GetIndex(prtcl,'NO')
                vc = 6
                ! Aerosol bins
                str = (nc-1)*nbins+1
                end = nc*nbins
                pa_vaerot(kk,ii,jj,str:end) = pa_vaerot(kk,ii,jj,str:end) + &
                     ( aero(1,1,1:nbins)%volc(vc) - vaero_old(kk,ii,jj,str:end) )*rhono/pdn(kk,ii,jj)/tstep
                ! Hydrometeor bins
                str = (nc-1)*ncld+1
                end = nc*ncld
                pa_vcloudt(kk,ii,jj,str:end) = pa_vcloudt(kk,ii,jj,str:end) + &
                     ( cloud(1,1,1:ncld)%volc(vc) - vcloud_old(kk,ii,jj,str:end) )*rhono/pdn(kk,ii,jj)/tstep
                ! Rain drops
                str = (nc-1)*nprc+1
                end = nc*nprc
                pa_vprecpt(kk,ii,jj,str:end) = pa_vprecpt(kk,ii,jj,str:end) + &
                     ( precp(1,1,1:nprc)%volc(vc) - vprecp_old(kk,ii,jj,str:end) )*rhono/pdn(kk,ii,jj)/tstep
                ! Activated droplets
                str = (nc-1)*ncld+1
                end = nc*ncld
                pa_vactd(kk,ii,jj,str:end) = actd(1,1,1:ncld)%volc(vc)*rhono/pdn(kk,ii,jj)
             END IF

             IF (IsUsed(prtcl,'NH')) THEN
                nc = GetIndex(prtcl,'NH')
                vc = 7
                ! Aerosol bins
                str = (nc-1)*nbins+1
                end = nc*nbins
                pa_vaerot(kk,ii,jj,str:end) = pa_vaerot(kk,ii,jj,str:end) + &
                     ( aero(1,1,1:nbins)%volc(vc) - vaero_old(kk,ii,jj,str:end) )*rhonh/pdn(kk,ii,jj)/tstep
                ! Hydrometeor bins
                str = (nc-1)*ncld+1
                end = nc*ncld
                pa_vcloudt(kk,ii,jj,str:end) = pa_vcloudt(kk,ii,jj,str:end) + &
                     ( cloud(1,1,1:ncld)%volc(vc) - vcloud_old(kk,ii,jj,str:end) )*rhonh/pdn(kk,ii,jj)/tstep
                ! Rain drops
                str = (nc-1)*nprc+1
                end = nc*nprc
                pa_vprecpt(kk,ii,jj,str:end) = pa_vprecpt(kk,ii,jj,str:end) + &
                     ( precp(1,1,1:nprc)%volc(vc) - vprecp_old(kk,ii,jj,str:end) )*rhonh/pdn(kk,ii,jj)/tstep
                ! Activated droplets
                str = (nc-1)*ncld+1
                end = nc*ncld
                pa_vactd(kk,ii,jj,str:end) = actd(1,1,1:ncld)%volc(vc)*rhonh/pdn(kk,ii,jj)
             END IF

             ! Water (always used)
             ! ---------------------------------------
             nc = GetIndex(prtcl,'H2O')
             vc = 8
             ! Aerosol bins
             str = (nc-1)*nbins+1
             end = nc*nbins
             pa_vaerot(kk,ii,jj,str:end) = pa_vaerot(kk,ii,jj,str:end) + &
                  ( aero(1,1,1:nbins)%volc(vc) - vaero_old(kk,ii,jj,str:end) )*rhowa/pdn(kk,ii,jj)/tstep
             ! Hydrometeor bins
             str = (nc-1)*ncld+1
             end = nc*ncld
             pa_vcloudt(kk,ii,jj,str:end) = pa_vcloudt(kk,ii,jj,str:end) + &
                  ( cloud(1,1,1:ncld)%volc(vc) - vcloud_old(kk,ii,jj,str:end) )*rhowa/pdn(kk,ii,jj)/tstep
             ! Rain drops
             str = (nc-1)*nprc+1
             end = nc*nprc
             pa_vprecpt(kk,ii,jj,str:end) = pa_vprecpt(kk,ii,jj,str:end) + &
                  ( precp(1,1,1:nprc)%volc(vc) - vprecp_old(kk,ii,jj,str:end) )*rhowa/pdn(kk,ii,jj)/tstep
             ! Activated droplets
             str = (nc-1)*ncld+1
             end = nc*ncld
             pa_vactd(kk,ii,jj,str:end) = actd(1,1,1:ncld)%volc(vc)*rhowa/pdn(kk,ii,jj)
             ! ----------------------------------------
             
             !IF ( ANY(vcloud_old(kk,ii,jj,:) /= vcloud_old(kk,ii,jj,:)) ) WRITE(*,*) 'SALSASSA 1',kk,ii,jj
             !IF ( ANY(cloud(1,1,1:ncld)%volc(1) /= cloud(1,1,1:ncld)%volc(1)) ) WRITE(*,*) 'SALSASSA 2',kk,ii,jj
             IF ( ANY(precp(1,1,1:nprc)%volc(2) /= precp(1,1,1:nprc)%volc(2)) .OR.  &
                  ANY(precp(1,1,1:nprc)%volc(1) /= precp(1,1,1:nprc)%volc(1))) THEN
                WRITE(*,*) 'SALSASSA 3',kk,ii,jj
                WRITE(*,*) cloud(1,1,1:7)%numc
                WRITE(*,*) cloud(1,1,1:7)%volc(1)
                WRITE(*,*) cloud(1,1,1:ncld)%volc(2)
                WRITE(*,*) ncloud_old(kk,ii,jj,1:7)
                WRITE(*,*) vcloud_old(kk,ii,jj,1:ncld)
                WRITE(*,*) vcloud_old(kk,ii,jj,ncld+1:2*ncld)
                WRITE(*,*) '-----------------------------------'
                WRITE(*,*) precp(1,1,1:7)%numc
                WRITE(*,*) precp(1,1,1:nprc)%volc(1)
                WRITE(*,*) precp(1,1,1:nprc)%volc(2)
                WRITE(*,*) nprecp_old(kk,ii,jj,1:7)
                WRITE(*,*) vprecp_old(kk,ii,jj,1:nprc)
                WRITE(*,*) vprecp_old(kk,ii,jj,nprc+1:2*nprc)
                STOP
             END IF


             pa_gaerot(kk,ii,jj,1) = pa_gaerot(kk,ii,jj,1) + &
                  ( zgso4(1,1)/pdn(kk,ii,jj) - pa_gaerop(kk,ii,jj,1) )/tstep

             pa_gaerot(kk,ii,jj,2) = pa_gaerot(kk,ii,jj,2) + & 
                  ( zghno3(1,1)/pdn(kk,ii,jj) - pa_gaerop(kk,ii,jj,2) )/tstep

             pa_gaerot(kk,ii,jj,3) = pa_gaerot(kk,ii,jj,3) + &
                  ( zgnh3(1,1)/pdn(kk,ii,jj) - pa_gaerop(kk,ii,jj,3) )/tstep

             pa_gaerot(kk,ii,jj,4) = pa_gaerot(kk,ii,jj,4) + &
                  ( zgocnv(1,1)/pdn(kk,ii,jj) - pa_gaerop(kk,ii,jj,4) )/tstep

             pa_gaerot(kk,ii,jj,5) = pa_gaerot(kk,ii,jj,5) + &
                  ( zgocsv(1,1)/pdn(kk,ii,jj) - pa_gaerop(kk,ii,jj,5) )/tstep


             ! Tendency of water vapour mixing ratio is obtained from the change in RH during SALSA run.
             ! Assumes no temperature change during SALSA run.
             rt(kk,ii,jj) = rt(kk,ii,jj) + &
                  ( in_rv(1,1) - rv(kk,ii,jj) )/tstep

          END DO ! kk
       END DO ! ii
    END DO ! jj

  END SUBROUTINE run_SALSA

  !
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


END MODULE mo_salsa_driver
