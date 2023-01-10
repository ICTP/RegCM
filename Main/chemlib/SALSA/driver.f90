PROGRAM driver

  USE mo_submctl,  ONLY :  &
                   nmod,   &
                   pstand, &
                   pi6,    &
                   ncld,   &
                   nprc,   &
                   rhosu,  &
                   msu,    &
                   rhono,  &
                   rhonh,  &
                   rhowa,  &
                   mwa,    &
                   mair,   &
                   rg,     &
                   avog,   &
                   time
  USE mo_salsa_driver
  USE mo_salsa_init
  USE mo_salsa_sizedist
  USE mo_kind, ONLY: dp
  USE mo_salsa_cloud
  USE class_componentIndex
  USE mo_salsa_dynamics, ONLY: satvaph2o

  IMPLICIT NONE

  INTEGER :: ii, jj


  !**********************************************************
  !*                                                        *
  !* I) Tracers/subroutines to be coupled with host model   *
  !*                                                        *
  !*  NB: When coupling, check units carefully!!!!          *
  !*                                                        *
  !**********************************************************

  !-----------------------------------------------------------------
  !-- Tracers (provided by/to host model) --------------------------
  !-----------------------------------------------------------------

  !-- bin indices -----------------
   INTEGER, PARAMETER :: &
   nbin(2) = (/ 3, 7 /)   ! number of bins in each main regime

  INTEGER, PARAMETER ::    &
    in1a = 1,              &
    in2a = in1a + nbin(1), &
    in2b = in2a + nbin(2), &
    fn1a = in2a - 1,       &
    fn2a = fn1a + nbin(2), &
    fn2b = fn2a + nbin(2), &
    nbins = fn2b,          &
    pnx   = 5,             &
    pny   = 5,             &
    pnz   = 3,             &
    n4    = 4                ! number of compounds


  REAL(dp), ALLOCATABLE :: &
       pa_ncloudp(:,:,:,:), &           ! Cloud droplet number concentration (# kg-1)
       pa_vcloudp(:,:,:,:), &           ! Cloud droplet volume concentration (kg kg-1)
       pa_nprecpp(:,:,:,:), &           ! Rain drop number concentration (# kg-1)
       pa_vprecpp(:,:,:,:), &           ! Rain drop volume concentration (kg kg-1)
       pa_ncloudt(:,:,:,:), &           ! Cloud droplet number tendency
       pa_vcloudt(:,:,:,:), &           ! Cloud droplet volume tendency
       pa_nprecpt(:,:,:,:), &           ! Rain drop number tendency
       pa_vprecpt(:,:,:,:), &           ! Rain drop volume tendency
       pa_Rcdry(:,:,:,:),   &           ! Cloud dry radius
       pa_Rpdry(:,:,:,:),   &           ! Rain dry radius
       pa_Rcwet(:,:,:,:),   &           ! Cloud wet radius
       pa_Rpwet(:,:,:,:),   &           ! Rain drop wet radius
       pa_vactd(:,:,:,:),   &           ! actual tendency due to new droplet formation.
       pa_nactd(:,:,:,:)                ! Same for number concentration

  TYPE(ComponentIndex) :: prtcl ! Contains "getIndex" which gives the index for a given
                                ! aerosol component name, i.e. 1:SO4, 2:OC, 3:BC, 4:DU,
                                ! 5:SS, 6:NO, 7:NH, 8:H2O

  REAL(dp) :: &
       ng(nmod),sigmag(nmod),dpg(nmod) 

  !-- gas compound tracers ------------
  REAL(dp) :: &
       
       c_h2so4(kbdim,klev), & ! sulphuric acid concentration in gas phase
                              ! for each grid point (kbdim,klev) [#/m3]
       
       c_hno3(kbdim,klev),  & ! nitric acid concentration in gas phase
                              ! for each grid point (kbdim,klev) [#/m3]
       
       c_ocnv(kbdim,klev),  & ! non-volatile organic vapour concentration
                              ! for each grid point (kbdim,klev) [#/m3]
       
       c_ocsv(kbdim,klev)     ! semivolatile organic vapour concentration
                              ! for each grid point (kbdim,klev) [#/m3]


  !-- atmospheric conditions --------------
  REAL(dp) :: &
        dpmid(nbins),                   & ! mean diameter within a bin
        zcore(nbins),                   & !
        n_aero(kbdim,klev,nbins),       & ! number concentration of aerosol particles
        vols(kbdim,klev,nbins,n4),      &
        papp1(pnz,pnx,pny),             & ! atmospheric pressure at time t+1
                                          ! for each grid point (kbdim,klev) [Pa]
       
       prelhum(kbdim,klev),             & ! atmospheric relative humidity at time t+1
                                          ! for each grid point (kbdim,klev) [0-1]
       
       ptp1(pnz,pnx,pny),               & ! atmospheric temperature at time t+1
                                          ! for each grid point (kbdim,klev) [K]
       rv(pnz,pnx,pny),                 & ! Water vapor mixing ratio
       rs(pnz,pnx,pny),                 & ! Water vapour saturation mixing ratio
       wp(pnz,pnx,pny),                 & ! Vertical velocity (m s-1)
       pdn(pnz,pnx,pny),                & ! Air density (for normalizing concentrations)
       pa_naerop(pnz,pnx,pny,nbins),    & ! aerosol number concentration (# kg-1)
       pa_vaerop(pnz,pnx,pny,n4*nbins), & ! aerosol volume concentration (kg kg-1)
       pa_naerot(pnz,pnx,pny,nbins),    & ! Aerosol number tendency
       pa_vaerot(pnz,pnx,pny,n4*nbins), & ! Aerosol volume tendency
       pa_gaerot(pnz,pnx,pny,5),        & ! Gaseous tracer tendency
       rt(pnz,pnx,pny)                    ! Water vapour tendency


    REAL(dp)    :: pa_gaerop(pnz,pnx,pny,5),           &  ! Gaseous tracers [# kg]
                   tstep                                  ! time step

    INTEGER     :: str(n4), fnl(n4), nc(n4), iReason

    LOGICAL     :: dbg2

    REAL(dp)    :: pa_Radry(pnz,pnx,pny,nbins),   & ! Aerosol dry particle radius
                   pa_Rawet(pnz,pnx,pny,nbins),   & ! Aerosol wet radius
                   pa_rhop(pnz,pnx,pny,nbins)       ! Aerosol density (kg/m3)
                  

    CHARACTER(len=3) :: complist(n4)

    REAL(dp) :: temperature, pressure, relative_humidity, rempty, vempty, nempty

    REAL(dp) :: timein, timeout
    INTEGER :: i, icounter, &
         prunmode                     ! 1 = initialize, 2 = spinup, 3 = actual run

  OPEN(13,FILE='radii2.dat', STATUS='unknown')
  OPEN(14,FILE='output.dat',STATUS='unknown')
  OPEN(15,FILE='diagnose2.dat',STATUS='unknown')
  OPEN(16,FILE='diagnose3.dat',STATUS='unknown')

  OPEN(9, FILE='input.dat', STATUS='old')


  ! Set values to 'ambient tracers'. 

  ! NB! Later these variables provided by host model

  READ(9,*) temperature, pressure, relative_humidity ! read in from cloud parcel output
  
  ptp1    = temperature              ! ambient temperature
  papp1   = pressure*100._dp         ! ambient pressure

  ! initialization
  iReason = 1
  ! air density kg/m3
  pdn = mair*papp1(pnz-1,pnx-2,pny-2)/rg/ptp1(pnz-1,pnx-2,pny-2)
  ! Water vapor saturation mixing ratio (kg/m3)
  rs = mwa*satvaph2o(ptp1(pnz-1,pnx-2,pny-2))/(mair*papp1(pnz-1,pnx-2,pny-2))

  ! Water vapor mixing ratio (kg/m3)
  rv = relative_humidity*rs(pnz-1,pnx-2,pny-2)

  wp=.2_dp
  pa_rhop = 0._dp

  c_h2so4 = 0._dp!5.E14_dp
  c_hno3 = 8.4682e-15_dp*avog*1.e6_dp!5.E15_dp
  c_ocnv = 0._dp!5.E14_dp
  c_ocsv = 0._dp!1.E14_dp

  ! runmode for initialization
  prunmode = 1

  ! debug off
  dbg2 = .false.

  CALL salsa_initialize

  ! mean volume of the aerosol bins
  zcore(in1a:fn2b) = pi6 * aero(1,1,in1a:fn2b)%dmid**3

  ! Stdev of modes
  sigmag = (/ 2.0, 1.5, 1.59, 1.59, 2.0, 2.0, 2.0 /)
  ! Mean diameter of the modes (m)
  dpg = (/0.15, 1.0, 0.2, 0.2, 6.0, 6.0, 6.0/)*1.e-6_dp
  ! Number concentration of modes (#/cm3)
  ng = (/100.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0/)*1.e6_dp

  CALL size_distribution(kproma, kbdim,  klev,   &
       ng, dpg, sigmag, n_aero)

  ! number mixing ratio
  pa_naerop = 0._dp
  pa_naerop(pnz-1,pnx-2,pny-2,in1a:fn2b) = n_aero(1,1,in1a:fn2b)/pdn(pnz-1,pnx-2,pny-2)

  ! gas phase mixing ratios
  pa_gaerop = 0._dp
  pa_gaerop(pnz-1,pnx-2,pny-2,1) = c_h2so4(1,1)/pdn(pnz-1,pnx-2,pny-2)
  pa_gaerop(pnz-1,pnx-2,pny-2,2) = c_hno3(1,1)/pdn(pnz-1,pnx-2,pny-2)
  pa_gaerop(pnz-1,pnx-2,pny-2,3) = c_hno3(1,1)/pdn(pnz-1,pnx-2,pny-2)

  ! setup chemical components
  complist(1)='SO4'
  complist(2)='NO'
  complist(3)='NH'
  complist(4)='H2O'

  CALL ComponentIndexConstructor(prtcl, n4, 4, complist)


  ! sulfate aerosol mixing ratio
  pa_vaerop = 0._dp
  pa_vaerop(pnz-1,pnx-2,pny-2,in1a:fn2b) = &
       n_aero(1,1,in1a:fn2b)*zcore(in1a:fn2b)*rhosu/pdn(pnz-1,pnx-2,pny-2)

  str(3) = (nc(3)-1)*nbins
!  pa_vaerop(pnz-1,pnx-2,pny-2,str(3)+in1a:str(3)+fn2b)=&
!       2._dp/3._dp*n_aero(1,1,in1a:fn2b)*zcore(in1a:fn2b)*rhonh/pdn(pnz-1,pnx-2,pny-2)



  ! Allocate bins
  ALLOCATE(pa_ncloudp(pnz,pnx,pny,ncld))           ! Cloud droplet number mixing ratio
  ALLOCATE(pa_ncloudt(pnz,pnx,pny,ncld))           ! Cloud droplet number tendency
  ALLOCATE(pa_vcloudp(pnz,pnx,pny,n4*ncld))           ! Cloud droplet volume mixing ratio
  ALLOCATE(pa_vcloudt(pnz,pnx,pny,n4*ncld))           ! Cloud droplet volume tendency
  ALLOCATE(pa_nprecpp(pnz,pnx,pny,ncld))           ! Rain drop number tendency
  ALLOCATE(pa_nprecpt(pnz,pnx,pny,ncld))           ! Rain drop number tendency
  ALLOCATE(pa_vprecpp(pnz,pnx,pny,n4*nprc))           ! Rain drop volume tendency
  ALLOCATE(pa_vprecpt(pnz,pnx,pny,n4*nprc))           ! Rain drop volume tendency
  ALLOCATE(pa_Rcdry(pnz,pnx,pny,ncld))           ! Cloud dry radius
  ALLOCATE(pa_Rpdry(pnz,pnx,pny,nprc))           ! Rain dry radius
  ALLOCATE(pa_Rcwet(pnz,pnx,pny,ncld))           ! Cloud wet radius
  ALLOCATE(pa_Rpwet(pnz,pnx,pny,nprc))           ! Rain drop wet radius
  ALLOCATE(pa_vactd(pnz,pnx,pny,n4*ncld))           ! actual tendency due to new droplet formation.
  ALLOCATE(pa_nactd(pnz,pnx,pny,ncld))           ! Same for number concentration

  ! Initialize bins
  pa_ncloudp = 0._dp
  pa_vcloudp = 0._dp
  pa_nprecpp = 0._dp
  pa_vprecpp = 0._dp
  pa_Rcdry = 0._dp
  pa_Rpdry = 0._dp
  pa_Rcwet = 0._dp
  pa_Rpwet = 0._dp
  pa_vactd = 0._dp
  pa_nactd = 0._dp
  
  tstep = 1._dp
  
  wp=.2_dp
  pa_rhop = 0._dp
  rempty = 1.0000000000000000E-010
  nempty = 1._dp
  vempty = 0._dp

  CALL CPU_TIME(timein)
  
  icounter = 1

  DO ii = 1, 1000000

     time = real(ii,dp)

     IF(ii == icounter) THEN
        write(6,*) ii
        icounter=icounter + 100
     END IF

     ! Initialize tendencies
     pa_naerot  = 0._dp
     pa_vaerot  = 0._dp
     pa_ncloudt = 0._dp
     pa_vcloudt = 0._dp
     pa_nprecpt = 0._dp
     pa_vprecpt = 0._dp
     pa_gaerot  = 0._dp
     rt         = 0._dp

     ! loop over aerosol bins
     DO jj = in1a, fn2b 
        ! Calculate dry radius
        
        ! initialize
        pa_Radry(pnz-1,pnx-2,pny-2,jj) = aero(1,1,jj)%dmid/2._dp
        pa_Rawet(pnz-1,pnx-2,pny-2,jj) = aero(1,1,jj)%dmid/2._dp
        
        IF(pa_naerop(pnz-1,pnx-2,pny-2,jj) > nempty) THEN
           
           ! indices for aerosol bins
           nc (1) = GetIndex(prtcl,'SO4')
           str(1) = (nc(1)-1)*nbins
           fnl(1) = nc(1)*nbins

           nc(2) = GetIndex(prtcl,'NO')
           str(2) = (nc(2)-1)*nbins
           fnl(2) = nc(2)*nbins

           nc(3) = GetIndex(prtcl,'NH')
           str(3) = (nc(3)-1)*nbins
           fnl(3) = nc(3)*nbins

           nc(4) = GetIndex(prtcl,'H2O')
           str(4) = (nc(4)-1)*nbins
           fnl(4) = nc(4)*nbins

           ! calculate dry radius for non-empty bins
           pa_Radry(pnz-1,pnx-2,pny-2,jj) = ((pa_vaerop(pnz-1,pnx-2,pny-2,str(1)+jj)/rhosu + & ! sulfate
                                              pa_vaerop(pnz-1,pnx-2,pny-2,str(2)+jj)/rhono + & ! nitrate
                                              pa_vaerop(pnz-1,pnx-2,pny-2,str(3)+jj)/rhonh)/ & ! ammonia
                                              pa_naerop(pnz-1,pnx-2,pny-2,jj)/pi6)**(1._dp/3._dp)/2._dp
           
           ! calculate wet radius for non-empty bins
           ! ---------------------------------------
        
           pa_Rawet(pnz-1,pnx-2,pny-2,jj) = ((pa_vaerop(pnz-1,pnx-2,pny-2,str(1)+jj)/rhosu + & ! sulfate
                                              pa_vaerop(pnz-1,pnx-2,pny-2,str(2)+jj)/rhono + & ! nitrate
                                              pa_vaerop(pnz-1,pnx-2,pny-2,str(3)+jj)/rhonh + & ! ammonia
                                              pa_vaerop(pnz-1,pnx-2,pny-2,str(4)+jj)/rhowa)/ & ! water
                                              pa_naerop(pnz-1,pnx-2,pny-2,jj)/pi6)**(1._dp/3._dp)/2._dp
 ELSE

           pa_Radry(pnz-1,pnx-2,pny-2,jj)  = rempty
           pa_Rawet(pnz-1,pnx-2,pny-2,jj)  = rempty
           pa_naerop(pnz-1,pnx-2,pny-2,jj) = vempty
    pa_vaerop(pnz-1,pnx-2,pny-2,jj) = vempty

        END IF

     END DO

     ! loop over cloud bins

     DO jj = 1, ncld 

        ! indices for cloud bins
        nc (1) = GetIndex(prtcl,'SO4')
        str(1) = (nc(1)-1)*ncld
     
        nc(2) = GetIndex(prtcl,'NO')
        str(2) = (nc(2)-1)*ncld
        
        nc(3) = GetIndex(prtcl,'NH')
        str(3) = (nc(3)-1)*ncld

        nc(4) = GetIndex(prtcl,'H2O')
        str(4) = (nc(4)-1)*ncld

 IF(pa_ncloudp(pnz-1,pnx-2,pny-2,jj) > nempty) THEN
         ! calculate dry radius for cloud bins
         pa_Rcdry(pnz-1,pnx-2,pny-2,jj) = ((pa_vcloudp(pnz-1,pnx-2,pny-2,str(1)+jj)/rhosu + & ! sulfate
                                                   pa_vcloudp(pnz-1,pnx-2,pny-2,str(2)+jj)/rhono + & ! nitrate
                                                   pa_vcloudp(pnz-1,pnx-2,pny-2,str(3)+jj)/rhonh)/ & ! ammonium
                                                   pa_ncloudp(pnz-1,pnx-2,pny-2,jj)/pi6)**(1._dp/3._dp)/2._dp
        
         ! calculate wet radius for cloud bins
         pa_Rcwet(pnz-1,pnx-2,pny-2,jj) = ((pa_vcloudp(pnz-1,pnx-2,pny-2,str(1)+jj)/rhosu + & ! sulfate
                                                   pa_vcloudp(pnz-1,pnx-2,pny-2,str(2)+jj)/rhono + & ! water (!! hard coded index)
                                                   pa_vcloudp(pnz-1,pnx-2,pny-2,str(3)+jj)/rhonh + & ! water (!! hard coded index)
                                                   pa_vcloudp(pnz-1,pnx-2,pny-2,str(4)+jj)/rhowa) / & ! water (!! hard coded index)
               pa_ncloudp(pnz-1,pnx-2,pny-2,jj)/pi6)**(1._dp/3._dp)/2._dp
 ELSE
            pa_Rcdry(pnz-1,pnx-2,pny-2,jj) = rempty
            pa_Rcwet(pnz-1,pnx-2,pny-2,jj) = rempty
  pa_ncloudp(pnz-1,pnx-2,pny-2,jj) = vempty
  pa_vcloudp(pnz-1,pnx-2,pny-2,jj) = vempty
 END IF
     END DO

     ! loop over prec bins
     DO jj = 1, nprc

        ! indices for cloud bins
        nc (1) = GetIndex(prtcl,'SO4')
        str(1) = (nc(1)-1)*nprc
     
        nc(2) = GetIndex(prtcl,'NO')
        str(2) = (nc(2)-1)*nprc
        
        nc(3) = GetIndex(prtcl,'NH')
        str(3) = (nc(3)-1)*nprc

        nc(4) = GetIndex(prtcl,'H2O')
        str(4) = (nc(4)-1)*nprc


 IF(pa_nprecpp(pnz-1,pnx-2,pny-2,jj) > nempty) THEN
            ! calculate dry radius for precipitation bins
            pa_Rpdry(pnz-1,pnx-2,pny-2,jj) = ((pa_vprecpp(pnz-1,pnx-2,pny-2,str(1)+jj)/rhosu + & ! sulfate
                                                   pa_vprecpp(pnz-1,pnx-2,pny-2,str(2)+jj)/rhono + & ! nitrate
                                                   pa_vprecpp(pnz-1,pnx-2,pny-2,str(3)+jj)/rhonh)/ & ! nitrate
                                                   pa_nprecpp(pnz-1,pnx-2,pny-2,jj)/pi6)**(1._dp/3._dp)/2._dp
            ! calculate wet radius for precipitation bins
            pa_Rpwet(pnz-1,pnx-2,pny-2,jj) = ((pa_vprecpp(pnz-1,pnx-2,pny-2,str(1)+jj)/rhosu + & ! sulfate
                                                   pa_vprecpp(pnz-1,pnx-2,pny-2,str(2)+jj)/rhono + & ! nitrate
                                                   pa_vprecpp(pnz-1,pnx-2,pny-2,str(3)+jj)/rhonh + & ! nitrate
                                                   pa_vprecpp(pnz-1,pnx-2,pny-2,str(4)+jj)/rhowa)/ & ! water
                                                   pa_nprecpp(pnz-1,pnx-2,pny-2,jj)/pi6)**(1._dp/3._dp)/2._dp
 ELSE
            pa_Rpdry(pnz-1,pnx-2,pny-2,jj) = rempty
            pa_Rpwet(pnz-1,pnx-2,pny-2,jj) = rempty
  pa_nprecpp(pnz-1,pnx-2,pny-2,jj) = vempty
  pa_vprecpp(pnz-1,pnx-2,pny-2,jj) = vempty
END IF
     END DO

     rs = mwa*satvaph2o(ptp1(pnz-1,pnx-2,pny-2))/(mair*papp1(pnz-1,pnx-2,pny-2))
     str(2) = (nc(2)-1)
     str(4) = (nc(4)-1)
     write(13,*) pa_Radry(pnz-1,pnx-2,pny-2,1:nbins),pa_Rcwet(pnz-1,pnx-2,pny-2,1:ncld),pa_Rpwet(pnz-1,pnx-2,pny-2,1:nprc)
     write(14,*) time,pa_vaerop(pnz-1,pnx-2,pny-2,str(2)*nbins+1:str(2)*nbins+10), &
                 pa_vcloudp(pnz-1,pnx-2,pny-2,str(2)*ncld+1:str(2)*ncld+ncld),   &
                 pa_vcloudp(pnz-1,pnx-2,pny-2,str(2)*nprc+1:str(2)*nprc+nprc)
!     write(6,*) time, temperature, rv(pnz-1,pnx-2,pny-2)/rs(pnz-1,pnx-2,pny-2)

     CALL run_SALSA(pnx,        pny,        pnz,        n4,          &
                    papp1,      ptp1,       rv,         rt,          &
                    rs,         wp,         pdn,                     &
                    pa_naerop,  pa_naerot,  pa_vaerop,  pa_vaerot,   &
                    pa_ncloudp, pa_ncloudt, pa_vcloudp, pa_vcloudt,  &
                    pa_nprecpp, pa_nprecpt, pa_vprecpp, pa_vprecpt,  &
                    pa_nactd,   pa_vactd,   pa_gaerop,  pa_gaerot,   &
                    pa_Radry,   pa_Rcdry,   pa_Rpdry,   pa_Rawet,    &
                    pa_Rcwet,   pa_Rpwet,   pa_rhop,    prunmode,    &
                    prtcl,      tstep,      dbg2                     )

    READ(9,*,IOSTAT=iReason) temperature, pressure, relative_humidity

!     IF(time < 200._dp) THEN
!        temperature =temperature-0.02_dp
!     ELSE
!        temperature =temperature+0.02_dp
!     END IF

     IF(time == 400._dp) STOP 'Maximum time reached. Simulation ends normally.'
 

     IF(iReason < 0) STOP 'End of input.dat file reached. Simulation ends normally.'     

     ptp1    = temperature
     papp1   = pressure*100._dp

     pdn(pnz-1,pnx-2,pny-2) = mair*papp1(pnz-1,pnx-2,pny-2)/rg/ptp1(pnz-1,pnx-2,pny-2)
     ! Water vapor saturation mixing ratio (kg/m3)
     rs(pnz-1,pnx-2,pny-2) = 0.622_dp*exp(20.386_dp-5132/ptp1(pnz-1,pnx-2,pny-2))*133.32_dp/& 
          papp1(pnz-1,pnx-2,pny-2)*pdn(pnz-1,pnx-2,pny-2)

     pa_naerop  = pa_naerop  + pa_naerot  * tstep
     pa_vaerop  = pa_vaerop  + pa_vaerot  * tstep
     pa_ncloudp = pa_ncloudp + pa_ncloudt * tstep
     pa_vcloudp = pa_vcloudp + pa_vcloudt * tstep
     pa_nprecpp = pa_nprecpp + pa_nprecpt * tstep
     pa_vprecpp = pa_vprecpp + pa_vprecpt * tstep     
     pa_gaerop  = pa_gaerop  + pa_gaerot  * tstep
     rv         = rv         + rt         * tstep
     prunmode   = 3

  END DO
  
665 FORMAT(99(E11.4,1X))

  CALL CPU_TIME(timeout)

END PROGRAM driver

