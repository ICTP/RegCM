MODULE mo_salsa_cloud

  !*********************************************************
  !  MOD_AERO CLOUD
  !*********************************************************
  ! 
  ! Purpose: Calculates the number of activated cloud 
  ! droplets according to parameterizations by:
  !
  ! Abdul-Razzak et al: "A parameterization of aerosol activation - 
  !                      3. Sectional representation"
  !                      J. Geophys. Res. 107, 10.1029/2001JD000483, 2002. 
  !                      [Part 3]
  !
  ! Abdul Razzak et al: "A parameterization of aerosol activation - 
  !                      1. Single aerosol type"
  !                      J. Geophys. Res. 103, 6123-6130, 1998. 
  !                      [Part 1]
  !
  !
  ! Interface:
  ! ----------
  ! Called from main aerosol model
  ! 
  !
  ! Coded by:
  ! ---------
  ! T. Anttila (FMI)     2007
  ! H. Kokkola (FMI)     2007
  ! A.-I. Partanen (FMI) 2007
  !

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

  !*********************************************************

CONTAINS 

  SUBROUTINE cloud_activation(kproma, kbdim, klev,   &
                              temp,   pres,  rv,     &
                              rs,     w,     paero,  &
                              pcloud, pactd          )
    
    USE mo_kind,        ONLY : dp

    USE mo_constants,   ONLY : g

    USE mo_submctl, ONLY :                      &
         lsactbase,                                 & ! Activation at cloud base
         lsactintst,                                & ! Activation of interstitial particles
         rg,                                        & ! molar gas constant 
                                                      ! [J/(mol K)]
         slim,                                      &
         surfw0,                                    & ! surface tension 
                                                      ! of water [J/m2]
         nbin,                                      & ! number of size bins 
                                                      ! in subranges
         nlim,                                      & ! lowest possible particle conc. in a bin [#/m3]
         rhosu, msu,                                & ! properties of compounds
         rhooc, moc,                                &
         rhono, mno,                                &
         rhonh, mnh,                                &
         rhobc, mbc,                                &
         rhoss, mss,                                &
         rhowa, mwa,                                &
         rhodu, mdu,                                &
         pi,                                        &
         pi6,                                       &
         cpa,                                       &
         mair,                                      &
         in1a,in2a,in2b,fn1a,fn2a, fn2b,            & ! size regime bin indices
         t_section,                                 & ! Data type for cloud/rain drops
         t_parallelbin,                             & ! Data type for mapping indices between parallel size bins
         ica,fca,icb,fcb,ira,fra,                   & ! Cloud bin indices
         ncld                                         ! Total number of cloud bins

    IMPLICIT NONE

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) ::              &
             kproma,                    & ! number of horiz. grid points 
             kbdim,                     & ! dimension for arrays 
             klev                       ! number of vertical levels

    REAL(dp), INTENT(in) ::             &
             pres(kbdim,klev),          &
             temp(kbdim,klev),          &
             w(kbdim,klev)
    
    REAL(dp), INTENT(inout) :: rv(kbdim,klev) ! Water vapor mixing ratio
    REAL(dp), INTENT(in)    :: rs(kbdim,klev) ! Saturation vapor mixing ratio

    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld),  &
                                      paero(kbdim,klev,fn2b)

    ! Properties of newly activate particles
    TYPE(t_section), INTENT(out) :: pactd(kbdim,klev,ncld)

    ! Properties of newly activated interstitial particles (for diagnostics)
    TYPE(t_section) :: pactd_intst(kbdim,klev,ncld)


    !-- local variables --------------
    INTEGER :: ii, jj, kk, mm, cc             ! loop indices

    INTEGER :: hlp      ! helper variables for index calculations
    INTEGER :: bcrita(kbdim,klev),bcritb(kbdim,klev) ! Index of the critical aerosol bin for regimes a and b

    REAL(dp) :: sfty = 1.e-30_dp  ! For safety

    REAL(dp) ::                         &
             sil,                       & ! critical supersaturation 
                                          !     at the upper bound of the bin 
             siu,                       & !  "  at the lower bound of the bin
             scrit(fn2b),               & !  "  at the center of the bin
             aa(kbdim,klev),                        & ! curvature (Kelvin) effect [m]
             bb,                        & ! solute (Raoult) effect [m3]
             ns(kbdim,klev,fn2b),                  & ! number of moles of solute
             nshi,                      & !             " at the upper bound of the bin
             nslo,                      & !             " at the lower bound of the bin
             ratio,                     & ! volume ratio
             vmiddle(kbdim,klev,fn2b),  & ! volume in the middle of the bin [m3]
             s_max,                     & ! maximum supersaturation
             s_eff,                     & ! effective supersaturation
             x, x1, x2, x3, a1, sum1,   & ! technical variables
             ka1,                       & ! thermal conductivity
             dv1,                       & ! diffusion coefficient
             Gc,                        & ! growth coefficient
             alpha,                     & ! see Abdul-Razzak and Ghan, part 3
             gamma,                     & ! see Abdul-Razzak and Ghan, part 3
             L,                         & ! latent heat of evaporation
             ps,                        & ! saturation vapor pressure of water [Pa]
             khi,                       & ! see Abdul-Razzak and Ghan, part 3
             theta,                     & ! see Abdul-Razzak and Ghan, part 3
             frac(kbdim,klev,fn2b),                & ! fraction of activated droplets in a bin
             ntot,                      & ! total number conc of particles [#/m3]
             dinsol(fn2b),              & ! diameter of the insoluble fraction [m]
             dinsolhi,                  & !    "   at the upper bound of a bin [m]
             dinsollo,                  & !    "   at the lower bound of a bin [m]
             zdcrit(kbdim,klev,fn2b),   & ! critical diameter [m]
             zdcstar(kbdim,klev),                   & ! Critical diameter corresponding to Smax
             !zdcint,                    &  ! Critical diameter for calculating activation of interstitial particles
             !zvcrhi,zvcrlo,             & ! critical volume boundaries for partially activated bins
             zdcrhi(kbdim,klev,fn2b),   & ! Critical diameter at the high end of a bin
             zdcrlo(kbdim,klev,fn2b),   & ! Critical diameter at the low end of a bin
             V,                         & ! updraft velocity [m/s]
             rref,                      & ! reference radius [m]
             A, B, D_p0, dmx,           & !
             vlo, k,                    &
             zvact,               & ! Total volume concentration of the newly activated droplets
             zcoreact,            & ! Volume concentration of the newly activated dry CCN
             zvtot,               &
             zcore,               &
             zrho(8)

    REAL(dp) :: testi


    ! ------------------------------------------------------------------
    ! Initialization
    ! ------------------------------------------------------------------
    zrho(1:8) = (/rhosu,rhooc,rhono,rhonh,rhobc,rhodu,rhoss,rhowa/)

    bb = 6._dp*mwa/(pi*rhowa)             ! Raoult effect [m3]
                                          ! NOTE!
                                          ! bb must be multiplied
                                          ! by the number of moles of
                                          ! solute
    zdcrit(:,:,:) = 0._dp
    zdcrlo(:,:,:) = 0._dp
    zdcrhi(:,:,:) = 0._dp
    frac(:,:,:) = 0._dp
    bcrita(:,:) = fn2a
    bcritb(:,:) = fn2b

    DO jj = 1,klev    ! vertical grid
       DO ii = 1,kproma ! horizontal grid

          vmiddle(ii,jj,in1a:fn2b) = pi6*paero(ii,jj,in1a:fn2b)%dmid**3
          DO kk = 1,ncld
             pactd(ii,jj,kk)%volc(:) = 0._dp
             pactd(ii,jj,kk)%numc = 0._dp
          END DO
          aa(ii,jj) = 4._dp*mwa*surfw0/(rg*rhowa*temp(ii,jj)) ! Kelvin effect [m]

       END DO
    END DO
    ! Get moles of solute at the middle of the bin
    CALL getSolute(kproma,kbdim,klev,paero,ns)

    ! ----------------------------------------------------------------


    ! -------------------------------------
    ! Interstitial activation
    ! -------------------------------------
    IF ( lsactintst ) THEN

       !write(*,*) 'ACT'
       CALL actInterst(kproma,kbdim,klev,paero,pcloud,rv,rs,aa)
       
       ! Update the moles of solute after interstitial activation
       CALL getSolute(kproma,kbdim,klev,paero,ns)
    END IF 

    
    ! -----------------------------------
    ! Activation at cloud base
    ! -----------------------------------
    IF ( lsactbase ) THEN

       DO jj = 1,klev    ! vertical grid
          DO ii = 1,kproma ! horizontal grid
          
             A  = aa(ii,jj) * 1.e6_dp                            !     "         [um] 
             x  = 4._dp*aa(ii,jj)**3/(27._dp*bb)

             ! Get the critical supersaturation for aerosol bins for necessary summation terms (sum1 & ntot)
             ntot = 0._dp
             sum1 = 0._dp

             !-- sums in equation (8), part 3
             ntot = ntot + SUM(paero(ii,jj,in1a:fn2b)%numc)

             scrit(in1a:fn2b) = exp(sqrt(x/max(1.e-28,ns(ii,jj,in1a:fn2b)))) - 1._dp

             sum1 = sum1 + SUM(paero(ii,jj,in1a:fn2b)%numc/scrit(in1a:fn2b)**(2._dp/3._dp))

             IF(ntot < nlim) CYCLE
             IF (w(ii,jj) <= 0._dp) CYCLE
             V  = w(ii,jj)

             !-- latent heat of evaporation [J/kg]
             L     = 2.501e6_dp-2370._dp*(temp(ii,jj)-273.15_dp)

             !-- saturation vapor pressure of water [Pa]
             a1    = 1._dp-(373.15_dp/temp(ii,jj))
             ps    = 101325._dp*                                                 &
                  exp(13.3185_dp*a1-1.976_dp*a1**2-0.6445_dp*a1**3-0.1299_dp*a1**4)     

             !-- part 1, eq (11)
             alpha = g*mwa*L/(cpa*rg*temp(ii,jj)**2)-                            &
                  g*mair/(rg*temp(ii,jj))

             !-- part 1, eq (12)
             gamma = rg*temp(ii,jj)/(ps*mwa) &
                  + mwa*L**2/(cpa*pres(ii,jj)*mair*temp(ii,jj))

             !-- diffusivity [m2/s], Seinfeld and Pandis (15.65)
             x1 = pres(ii,jj) / 101325._dp
             dv1= 1.e-4_dp * (0.211_dp/x1) * ((temp(ii,jj)/273._dp)**1.94_dp)
             
             rref = 10.e-9_dp
             !-- corrected diffusivity, part 1, eq (17)
             ! dv = dv1 / (rref/(rref + deltaV) + (dv1/(rref * alphac)) *        &
             !     SQRT(2.*pi*mwa/(rg*temp(ii,jj))))
             
             !-- thermal conductivity [J/(m s K)], Seinfeld and Pandis (15.75)
             ka1= 1.e-3_dp * (4.39_dp + 0.071_dp * temp(ii,jj))

             !-- growth coefficient, part 1, eq (16)
             !-- (note: here uncorrected diffusivities and conductivities are used
             !    based on personal communication with H. Abdul-Razzak, 2007)
             Gc = 1._dp/(rhowa*rg*temp(ii,jj)/(ps*dv1*mwa) +                      &
                  L*rhowa/(ka1*temp(ii,jj)) * (L*mwa/(temp(ii,jj)*rg)-1._dp))

             !-- effective critical supersaturation: part 3, eq (8)
             s_eff = (ntot/sum1)**(3._dp/2._dp)

             !-- part 3, equation (5)
             
             theta = ((alpha*V/Gc)**(3._dp/2._dp))/(2._dp*pi*rhowa*gamma*ntot)
          
             !-- part 3, equation (6)
             khi = (2._dp/3._dp)*aa(ii,jj)*SQRT(alpha*V/Gc)

             !-- maximum supersaturation of the air parcel: part 3, equation (9)
             s_max = s_eff / SQRT(0.5_dp*(khi/theta)**(3._dp/2._dp)              &
                  + ((s_eff**2)/(theta+3._dp*khi))**(3._dp/4._dp))

             !-- Juha: Get the critical diameter corresponding to the maximum supersaturation
             ! ------ NÄISTÄ LASKUISTA PUUTTUU LIUKENEMATTOMAN COREN VAIKUTUS?
             zdcstar = 2._dp*aa(ii,jj)/(3._dp*s_max)

             !WRITE(*,*) s_max
             !WRITE(*,*) zdcstar

             DO kk = in1a, fn2b

                IF (paero(ii,jj,kk)%numc < nlim) CYCLE

                !-- moles of solute in particle at the upper bound of the bin
                nshi = ns(ii,jj,kk)*paero(ii,jj,kk)%vratiohi
                
                !-- critical supersaturation
                sil = exp(sqrt(x/nshi)) - 1._dp

                IF(s_max < sil) CYCLE 

                !-- moles of solute at the lower bound of the bin:
                nslo = ns(ii,jj,kk)*paero(ii,jj,kk)%vratiolo
                
                !-- critical supersaturation
                siu = exp(sqrt(x/nslo)) - 1._dp

                !-- fraction of activated in a bin, eq (13), part 3
                frac(ii,jj,kk) = min(1._dp,log(s_max/sil)/log(siu/sil))

                !-- Critical diameters for each bin and bin edges
                zdcrlo(ii,jj,kk) = 2._dp*sqrt(3._dp*nslo*bb/aa(ii,jj))
                zdcrhi(ii,jj,kk) = 2._dp*sqrt(3._dp*nshi*bb/aa(ii,jj))
                zdcrit(ii,jj,kk) = 2._dp*sqrt(3._dp*ns(ii,jj,kk)*bb/aa(ii,jj))
                
             END DO ! kk

             ! Find critical bin
             DO kk = in1a,fn2a
                IF (frac(ii,jj,kk) < 1._dp .AND. frac(ii,jj,kk) > 0._dp) THEN
                   bcrita(ii,jj) = kk
                   EXIT
                END IF
             END DO
             DO kk = in2b,fn2b
                IF (frac(ii,jj,kk) < 1._dp .AND. frac(ii,jj,kk) > 0._dp) THEN
                   bcritb(ii,jj) = kk
                   EXIT
                END IF
             END DO

          END DO ! ii
       
       END DO ! jj

       CALL activate3(kproma,kbdim,klev,paero,bcrita,bcritb,  &
                      zdcrit, zdcrlo, zdcrhi, zdcstar, pactd  )

    END IF ! lsactbase

    ! Regimes 1 and 2a
    !CALL activate2(kproma,  kbdim,  klev,     &
    !               ica,     fca,    paero,    &
    !               zdcrit,  zdcrlo, zdcrhi,   &
    !               zdcstar, pactd             )

    ! Regime 2b
    !CALL activate2(kproma, kbdim, klev, icb, fcb,   &
    !                  paero,  frac,  zdcrit, pactd         )

  END SUBROUTINE cloud_activation

! -----------------------------------------------------------------

  SUBROUTINE getSolute(kproma,kbdim,klev,paero,pns)
    USE mo_kind, ONLY : dp
    USE mo_submctl, ONLY : t_section,nlim,       &
                               in1a,fn1a,            &
                               in2a,fn2a,            &
                               in2b,fn2b,            &
                               rhosu, rhooc, rhobc,  &
                               rhonh, rhono, rhodu,  &
                               rhoss,                &
                               msu, moc, mbc,        &
                               mnh, mno, mdu,        &
                               mss
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kproma,kbdim,klev
    TYPE(t_section), INTENT(IN) :: paero(kbdim,klev,fn2b)
    REAL(dp), INTENT(OUT) :: pns(kbdim,klev,fn2b) 

    INTEGER :: ii,jj,kk

    pns = 0._dp

    DO jj = 1,klev

       DO ii = 1,kproma

          !-- subrange 1a

          !-- calculation of critical superaturation in the middle of the bin

          !-- volume in the middle of the bin

          DO kk = in1a, fn1a
             IF (paero(ii,jj,kk)%numc > nlim) THEN

                !-- number of moles of solute in one particle [mol]
                pns(ii,jj,kk) = (3._dp*paero(ii,jj,kk)%volc(1)*rhosu/msu  +   &
                     paero(ii,jj,kk)%volc(2)*rhooc/moc)/                      &
                     paero(ii,jj,kk)%numc

             END IF
          END DO

          !-- subrange 2a

          DO kk = in2a, fn2a
             IF (paero(ii,jj,kk)%numc > nlim) THEN

                pns(ii,jj,kk) = (3._dp*paero(ii,jj,kk)%volc(1)*rhosu/msu  +   &
                     paero(ii,jj,kk)%volc(2) * rhooc/moc +                    &
                     paero(ii,jj,kk)%volc(6) * rhono/mno +                    &
                     paero(ii,jj,kk)%volc(7) * rhonh/mnh +                    &
                     2._dp*paero(ii,jj,kk)%volc(5) * rhoss/mss)/              &
                     paero(ii,jj,kk)%numc  

             END IF
          END DO

          !-- subrange 2b

          DO kk = in2b, fn2b
             IF (paero(ii,jj,kk)%numc > nlim) THEN

                pns(ii,jj,kk) = (3._dp*paero(ii,jj,kk)%volc(1)*rhosu/msu  +   &   
                        paero(ii,jj,kk)%volc(2) * rhooc/moc +                 &
                        paero(ii,jj,kk)%volc(6) * rhono/mno +                 &
                        paero(ii,jj,kk)%volc(7) * rhonh/mnh +                 &
                        2._dp*paero(ii,jj,kk)%volc(5) * rhoss/mss)/           &
                        paero(ii,jj,kk)%numc 

             END IF
          END DO


          
       END DO

    END DO


  END SUBROUTINE getSolute

! -----------------------------------------------------------------

  SUBROUTINE actInterst(kproma,kbdim,klev,paero,pcloud,prv,prs,paa)
    !
    ! Activate interstitial aerosols if they've reached their critical size
    ! HUOM TEKEE VAAN A-BINIT NYT!!
    !
    !
    ! 1. Formulate the profiles of Dwet within bins
    !      - Get the slopes between adjacent bin mids (known)
    !      - Interpolate Dwet to bin edges
    ! 2. Based on Dcrit and Dwet slopes, estimate the Ddry where Dwet becomes > Dcrit
    ! 3. Formulate the slopes for number concentration
    ! 4. Use the Dry limits from (2) as the integration limits if they are defined

    USE mo_kind, ONLY : dp
    USE mo_submctl, ONLY : t_section,nlim,pi6,ica,fca,icb,fcb, &
                               in1a,fn1a,in2a,fn2a,                &
                               nbins,ncld
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kproma,kbdim,klev
    TYPE(t_section), INTENT(INOUT) :: paero(kbdim,klev,nbins),  &
                                      pcloud(kbdim,klev,ncld)
    REAL(dp), INTENT(IN) :: prv(kbdim,klev),prs(kbdim,klev),    &  ! Water vapour and saturation mixin ratios
                            paa(kbdim,klev)                        ! Coefficient for Kelvin effect

    REAL(dp) :: zdcstar,zvcstar   ! Critical diameter/volume corresponding to S_LES
    REAL(dp) :: zvcint            ! Integration limit volume based on above
    REAL(dp) :: zactvol           ! Total volume of the activated particles
    
    TYPE(t_section) :: zactd(ncld)   ! Temporary storage for activated particles
    REAL(dp) :: Nact, Vact(8)        ! Helper variables for transferring the activated particles

    REAL(dp) :: Nmid, Nim1, Nip1     ! Bin number concentrations in current and adjacent bins 
    REAL(dp) :: dNmid, dNim1, dNip1  ! Density function value of the number distribution for current and adjacent bins

    REAL(dp) :: Vmid, Vim1, Vip1     ! Dry particle volume in the middle of the bin
    REAL(dp) :: Vlo, Vhi             ! Dry particle volume scaled to bin edges
    REAL(dp) :: Vlom1,Vhim1          ! - '' - For adjacent bins
    REAL(dp) :: Vlop1,Vhip1          !

    REAL(dp) :: Vwmid, Vwim1, Vwip1  ! Wet particle volume in the middle of the bin
    REAL(dp) :: Vwlo,Vwhi            ! Wet particle volume at bin edges
    REAL(dp) :: Vwlom1, Vwhim1       ! Wet particle volumes for adjacent bins
    REAL(dp) :: Vwlop1, Vwhip1       ! - '' -

    REAL(dp) :: zs1,zs2           ! Slopes for number concetration distributions within bins

    REAL(dp) :: N01,N02           ! Origin values for number distribution slopes
    REAL(dp) :: V01,V02           ! Origin values for wet particle volume slopes
    REAL(dp) :: Nnorm, Vnorm      ! Normalization factors for number and volume integrals
    
    REAL(dp) :: vcut,vint1,vint2  ! cut volume, integration limit volumes
    LOGICAL(dp) :: intrange(4)    ! Logical table for integration ranges depending on the shape of the wet size profile:
                                  ! [Vlo -- vint1][vint1 -- Vmid][Vmid -- vint2][vint1 -- Vhi]
    INTEGER :: cb,ab, ii,jj,ss

    
    DO jj = 1,klev

       DO ii = 1,kproma

          IF ( prv(ii,jj)/prs(ii,jj) > 1.000_dp ) THEN 

             zactd(1:ncld)%numc = 0._dp
             DO ss = 1,8
                zactd(1:ncld)%volc(ss) = 0._dp
             END DO

             ! Determine Dstar == critical diameter corresponding to the host model S
             zdcstar = 2._dp*paa(ii,jj)/( 3._dp*( (prv(ii,jj)/prs(ii,jj))-1._dp ) )
             zvcstar = pi6*zdcstar**3

             ! Loop over cloud droplet (and aerosol) bins
             
             !IF (ANY(paero(ii,jj,:)%numc<nlim)) WRITE(*,*) 'aero', paero(ii,jj,:)%dmid
             !IF (ANY(pcloud(ii,jj,:)%numc<nlim))

             DO cb = ica%cur, fca%cur
                ab = ica%par + (cb-ica%cur)

                IF ( paero(ii,jj,ab)%numc < nlim) CYCLE

                intrange = .FALSE.

                ! Define some parameters
                Nmid = paero(ii,jj,ab)%numc     ! Number concentration at the current bin center
                Vwmid = SUM(paero(ii,jj,ab)%volc(:))/Nmid  ! Wet volume at the current bin center
                Vmid = SUM(paero(ii,jj,ab)%volc(1:7))/Nmid ! Dry volume at the current bin center
                Vlo = Vmid*paero(ii,jj,ab)%vratiolo        ! Dry vol at low limit
                Vhi = Vmid*paero(ii,jj,ab)%vratiohi        ! Dry vol at high limit
                
                ! Number concentrations and volumes at adjacent bins (check for sizedistribution boundaries)
                IF ( ab > in1a ) THEN
                   Nim1 = paero(ii,jj,ab-1)%numc
                   IF (Nim1 > nlim) THEN
                      Vim1 = SUM(paero(ii,jj,ab-1)%volc(1:7))/Nim1
                      Vwim1 = SUM(paero(ii,jj,ab-1)%volc(:))/Nim1
                   ELSE
                      Vim1 = pi6*paero(ii,jj,ab-1)%dmid**3
                      Vwim1 = pi6*paero(ii,jj,ab-1)%dmid**3
                   END IF
                   Vlom1 = Vim1*paero(ii,jj,ab-1)%vratiolo
                   Vhim1 = Vim1*paero(ii,jj,ab-1)%vratiohi
                ELSE ! ab == in1a
                   Nim1 = nlim
                   Vim1 = Vlo/2._dp
                   Vlom1 = 0._dp
                   Vhim1 = Vlo
                   Vwim1 = Vwmid/3._dp ! Tää ny o vähä tämmöne
                END IF
                IF ( ab < fn2a ) THEN
                   Nip1 = paero(ii,jj,ab+1)%numc
                   IF (Nip1 > nlim) THEN
                      Vip1 = SUM(paero(ii,jj,ab+1)%volc(1:7))/Nip1
                      Vwip1 = SUM(paero(ii,jj,ab+1)%volc(:))/Nip1
                   ELSE
                      Vip1 = pi6*paero(ii,jj,ab+1)%dmid**3
                      Vwip1 = pi6*paero(ii,jj,ab+1)%dmid**3
                   END IF
                   Vlop1 = Vip1*paero(ii,jj,ab+1)%vratiolo
                   Vhip1 = Vip1*paero(ii,jj,ab+1)%vratiohi
                ELSE ! ab == fn2a
                   Nip1 = nlim
                   Vip1 = Vhi + 0.5_dp*(Vhi-Vlo)
                   Vlop1 = Vhi
                   Vhip1 = Vhi + (Vhi-Vlo)
                   Vwip1 = Vhip1  ! ....
                END IF

                ! Keeping thins smooth...
                Vip1 = MAX(Vhi,Vip1)
                Vim1 = MIN(Vlo,Vim1)

                ! First, make profiles of particle wet radius in 
                ! order to determine the integration boundaries
                zs1 = (Vwmid - MAX(Vwim1,0._dp))/(Vmid - Vim1)
                zs2 = (MAX(Vwip1,0._dp) - Vwmid)/(Vip1 - Vmid)

                ! Get the origin values for slope equations
                V01 = Vwmid - zs1*Vmid
                V02 = Vwmid - zs2*Vmid

                ! Get the wet sizes at bins edges
                Vwlo = MAX(V01 + zs1*Vlo, 0._dp)
                Vwhi = MAX(V02 + zs2*Vhi, 0._dp)
                
                ! Find out dry vol integration boundaries based on *zvcstar*:
                IF ( zvcstar < Vwlo .AND. zvcstar < Vwmid .AND. zvcstar < Vwhi ) THEN
                   ! Whole bin activates
                   vint1 = Vlo
                   vint2 = Vhi
                   
                   intrange(1:4) = .TRUE.

                ELSE IF ( zvcstar > Vwlo .AND. zvcstar > Vwmid .AND. zvcstar > Vwhi) THEN
                   ! None activates
                   vint1 = 999.
                   vint2 = 999.

                   intrange(1:4) = .FALSE.

                ELSE
                   ! Partial activation:
                   ! Slope1
                   vcut = (zvcstar - V01)/zs1  ! Where the wet size profile intersects the critical size (slope 1)
                   IF (vcut < Vlo .OR. vcut > Vmid) THEN
                      ! intersection volume outside the current size range -> set as the lower limit
                      vint1 = Vlo
                   ELSE
                      vint1 = vcut
                   END IF

                   ! Slope2
                   vcut = (zvcstar - V02)/zs2  ! Where the wet size profile intersects the critical size (slope 2)
                   IF (vcut < Vmid .OR. vcut > Vhi) THEN
                      ! Intersection volume outside the current size range -> set as the lower limit
                      vint2 = Vmid
                   ELSE
                      vint2 = vcut
                   END IF

                   ! Determine which size ranges have wet volume larger than the critical
                   intrange(1) = ( Vwlo > zvcstar )
                   intrange(2) = ( Vwmid > zvcstar )
                   intrange(3) = ( Vwmid > zvcstar )
                   intrange(4) = ( Vwhi > zvcstar )
                   
                END IF
                ! DONE WITH INTEGRATION LIMITS
                ! -------------------------------------
             
                ! Number concentration profiles within bins and integration for number of activated:
                ! -----------------------------------------------------------------------------------
                ! get density distribution values for number concentration
                dNim1 = Nim1/(Vhim1-Vlom1)
                dNip1 = Nip1/(Vhip1-Vlop1)
                dNmid = Nmid/(Vhi-Vlo)

                ! Get slopes
                zs1 = ( dNmid - dNim1 )/( Vmid - Vim1 )
                zs2 = ( dNip1 - dNmid )/( Vip1 - Vmid )
                
                N01 = dNmid - zs1*Vmid  ! Origins
                N02 = dNmid - zs2*Vmid  !
                
                ! Define normalization factors
                Nnorm = intgN(zs1,N01,Vlo,Vmid) + intgN(zs2,N02,Vmid,Vhi)
                Vnorm = intgV(zs1,N01,Vlo,Vmid) + intgV(zs2,N02,Vmid,Vhi)
                
                ! Accumulated variables
                zactvol = 0._dp
                Nact = 0._dp
                Vact(:) = 0._dp                
                
                ! Integration over each size range within a bin
                IF ( intrange(1) ) THEN
                   Nact = Nact + (Nmid/Nnorm)*intgN(zs1,N01,Vlo,vint1)
                   zactvol = zactvol + (Nmid*Vmid/Vnorm)*intgV(zs1,N01,Vlo,vint1)
                END IF

                IF ( intrange(2) ) THEN
                   Nact = Nact + (Nmid/Nnorm)*intgN(zs1,N01,vint1,Vmid)
                   zactvol = zactvol + (Nmid*Vmid/Vnorm)*intgV(zs1,N01,vint1,Vmid)
                END IF

                IF ( intrange(3) ) THEN
                   Nact = Nact + (Nmid/Nnorm)*intgN(zs2,N02,Vmid,vint2)
                   zactvol = zactvol + (Nmid*Vmid/Vnorm)*intgV(zs2,N02,Vmid,vint2)
                END IF

                IF ( intrange(4) ) THEN
                   Nact = Nact + (Nmid/Nnorm)*intgN(zs2,N02,vint2,Vhi)
                   zactvol = zactvol + (Nmid*Vmid/Vnorm)*intgV(zs2,N02,vint2,Vhi)
                END IF

                DO ss = 1,8
                   Vact(ss) = zactvol*( paero(ii,jj,ab)%volc(ss)/(Vmid*Nmid) )
                END DO

                ! Store the number concentration and mass of activated particles for current bins
                zactd(cb)%numc = Nact
                zactd(cb)%volc(:) = Vact(:)

                IF (zactd(cb)%numc < 0._dp) THEN
                   WRITE(*,*) Nim1,Nmid,Nip1
                   WRITE(*,*) dNim1,dNmid,dNip1
                   WRITE(*,*) Vim1,Vlo,Vmid,Vhi,Vip1
                   WRITE(*,*) N01,N02,zs1,zs2
                   WRITE(*,*) vint1,vint2
                END IF

             END DO ! cb

             IF (ANY(zactd(:)%numc < 0._dp)) THEN
                WRITE(*,*) 'NEGATIVE ACTIVATED'
                WRITE(*,*) zactd(:)%numc
                WRITE(*,*) '---------'
                WRITE(*,*) 'Setting MAX(Nact,0)'
                zactd(:)%numc = MAX(zactd(:)%numc, 0._dp)
                DO ss = 1,8
                   WHERE(zactd(:)%numc == 0._dp) zactd(:)%volc(ss) = 0._dp
                END DO

             END IF
             
             ! Apply the number and mass activated to aerosol and cloud bins
             !WRITE(*,*) zactd(1:7)%numc/paero(ii,jj,4:10)%numc
             paero(ii,jj,ica%par:fca%par)%numc =   &
                  MAX(0._dp, paero(ii,jj,ica%par:fca%par)%numc - zactd(ica%cur:fca%cur)%numc)
             pcloud(ii,jj,ica%cur:fca%cur)%numc = pcloud(ii,jj,ica%cur:fca%cur)%numc + zactd(ica%cur:fca%cur)%numc
             DO ss = 1,8
                paero(ii,jj,ica%par:fca%par)%volc(ss) =    &
                     MAX(0._dp, paero(ii,jj,ica%par:fca%par)%volc(ss) - zactd(ica%cur:fca%cur)%volc(ss))
                pcloud(ii,jj,ica%cur:fca%cur)%volc(ss) = pcloud(ii,jj,ica%cur:fca%cur)%volc(ss) + zactd(ica%cur:fca%cur)%volc(ss)
             END DO
           
          END IF ! RH limit
             
       END DO ! ii
          
    END DO ! jj

  END SUBROUTINE actInterst


  ! ----------------------------------------------
  
  SUBROUTINE activate3(kproma,kbdim,klev,paero,pbcrita,pbcritb, &
                       pdcrit, pdcrlo, pdcrhi, pdcstar, pactd   )
    !
    ! Gets the number and mass activated in the critical aerosol size bin
    ! TEKEE NYT VAIN a-binit!!!!
    !
    USE mo_kind, ONLY : dp
    USE mo_submctl, ONLY : t_parallelbin, t_section, pi6, nlim, fn2b, ncld,  &
                               in1a,in2a,fn1a,fn2a, ica,fca,icb,fcb
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: kproma,kbdim,klev
    TYPE(t_section), INTENT(IN) :: paero(kbdim,klev,fn2b)
    INTEGER, INTENT(IN) :: pbcrita(kbdim,klev),         & ! Index of the critical aerosol bin in regime a 
                           pbcritb(kbdim,klev)            ! Index of the critical aerosol bin in regime b
    REAL(dp), INTENT(IN) :: pdcrit(kbdim,klev,fn2b),    & ! Bin middle critical diameter
                            pdcrlo(kbdim,klev,fn2b),    & ! Critical diameter at low limit
                            pdcrhi(kbdim,klev,fn2b)       ! Critical diameter at high limit
    REAL(dp), INTENT(IN) :: pdcstar(kbdim,klev)           ! Critical diameter corresponding to Smax
    TYPE(t_section), INTENT(OUT) :: pactd(kbdim,klev,ncld) ! Properties of the maximum amount of newly activated droplets
    

    REAL(dp) :: zvcstar, zvcint
    REAL(dp) :: zs1,zs2             ! Slopes
    REAL(dp) :: Nmid, Nim1,Nip1
    REAL(dp) :: dNmid, dNim1, dNip1
    REAL(dp) :: Vmid,Vlo,Vhi
    REAL(dp) :: Vim1,Vlom1,Vhim1
    REAL(dp) :: Vip1,Vlop1,Vhip1
    REAL(dp) :: Nnorm,Vnorm,N01,N02    
    REAL(dp) :: zactvol

    REAL(dp) :: zhlp
    REAL(dp) :: zratio

    INTEGER :: ii,jj,ss,cb,ab

    DO jj = 1,klev
       DO ii = 1,kproma
          
          ! Determine the dry volume of the critical sized particle
          IF ( pdcstar(ii,jj) >= pdcrit(ii,jj,pbcrita(ii,jj)) ) THEN 

             zhlp = ( pdcstar(ii,jj) - pdcrit(ii,jj,pbcrita(ii,jj)) ) /   &
                    ( pdcrhi(ii,jj,pbcrita(ii,jj)) - pdcrit(ii,jj,pbcrita(ii,jj)) )
             
             zratio = (1._dp-zhlp)*pdcrit(ii,jj,pbcrita(ii,jj)) /  &
                      (SUM(paero(ii,jj,pbcrita(ii,jj))%volc(1:7))/paero(ii,jj,pbcrita(ii,jj))%numc/pi6)**(1._dp/3._dp) + &
                      zhlp*pdcrhi(ii,jj,pbcrita(ii,jj))/(paero(ii,jj,pbcrita(ii,jj))%vhilim/pi6)**(1._dp/3._dp)

          ELSE IF ( pdcstar(ii,jj) < pdcrit(ii,jj,pbcrita(ii,jj)) ) THEN
             
             zhlp = ( pdcrit(ii,jj,pbcrita(ii,jj)) - pdcstar(ii,jj) ) /   & 
                    ( pdcrit(ii,jj,pbcrita(ii,jj)) - pdcrlo(ii,jj,pbcrita(ii,jj)) ) 

             zratio = (1._dp-zhlp)*pdcrit(ii,jj,pbcrita(ii,jj)) /  &
                      (SUM(paero(ii,jj,pbcrita(ii,jj))%volc(1:7))/paero(ii,jj,pbcrita(ii,jj))%numc/pi6)**(1._dp/3._dp) + &
                      zhlp*pdcrlo(ii,jj,pbcrita(ii,jj))/(paero(ii,jj,pbcrita(ii,jj))%vlolim/pi6)**(1._dp/3._dp)

          END IF
          ! zvcstar is the DRY particle volume corresponding to the critical diameter at S_max
          zvcstar = pi6*(pdcstar(ii,jj)/zratio)**3
          zvcstar = MAX( zvcstar, paero(ii,jj,pbcrita(ii,jj))%vlolim )
          zvcstar = MIN( zvcstar, paero(ii,jj,pbcrita(ii,jj))%vhilim )

          ! Loop over cloud droplet (and aerosol) bins
          DO cb = ica%cur,fca%cur
             ab = ica%par + (cb-ica%cur)
         
             IF ( paero(ii,jj,ab)%numc < nlim) CYCLE
                
             ! Formulate a slope for Wet particle size within bins and integrate over 
             ! the particles larger than zvcstar
                
             Nmid = MAX(paero(ii,jj,ab)%numc, nlim)
             Vmid = SUM(paero(ii,jj,ab)%volc(1:7))/Nmid ! Dry bin mid volume
             Vlo = paero(ii,jj,ab)%vlolim      ! Mid dry volume scaled to bin low limit (this is mostly an educated guess... )
             Vhi = paero(ii,jj,ab)%vhilim      ! Same for high limit
                
             IF ( ab > in1a ) THEN
                Nim1 = MAX(paero(ii,jj,ab-1)%numc, nlim)
                IF (Nim1 > nlim) THEN
                   Vim1 = SUM(paero(ii,jj,ab-1)%volc(1:7))/Nim1
                ELSE
                   Vim1 = pi6*paero(ii,jj,ab-1)%dmid**3
                END IF
                Vlom1 = paero(ii,jj,ab-1)%vlolim
                Vhim1 = paero(ii,jj,ab-1)%vhilim
             ELSE ! ab == in1a
                Nim1 = nlim
                Vim1 = Vlo/2._dp
                Vlom1 = 0._dp
                Vhim1 = Vlo
             END IF

             IF ( ab < fn2a ) THEN
                Nip1 = MAX(paero(ii,jj,ab+1)%numc, nlim)
                IF (Nip1 > nlim) THEN
                   Vip1 = SUM(paero(ii,jj,ab+1)%volc(1:7))/Nip1
                ELSE
                   Vip1 = pi6*paero(ii,jj,ab+1)%dmid**3
                END IF
                Vlop1 = paero(ii,jj,ab+1)%vlolim
                Vhip1 = paero(ii,jj,ab+1)%vhilim
             ELSE ! ab == fn2a
                Nip1 = nlim
                Vip1 = Vhi + 0.5_dp*(Vhi-Vlo)
                Vlop1 = Vhi
                Vhip1 = Vhi + (Vhi-Vlo)
             END IF
                
             ! get density distribution values for 
             dNim1 = Nim1/(Vhim1-Vlom1)
             dNip1 = Nip1/(Vhip1-Vlop1)
             dNmid = Nmid/(Vhi-Vlo)
          
             ! Get slopes
             zs1 = ( dNmid - dNim1 )/( Vmid - Vim1 )
             zs2 = ( dNip1 - dNmid )/( Vip1 - Vmid )

             N01 = dNim1 - zs1*Vim1  ! Origins
             N02 = dNmid - zs2*Vmid  !
          
             ! Define normalization factors
             Nnorm = intgN(zs1,N01,Vlo,Vmid) + intgN(zs2,N02,Vmid,Vhi)
             Vnorm = intgV(zs1,N01,Vlo,Vmid) + intgV(zs2,N02,Vmid,Vhi)
             
             IF (zvcstar < Vmid) THEN 
             
                ! Use actual critical volume only in the critical bin, otherwise current bin limits
                zvcint = MAX(zvcstar, Vlo)

                pactd(ii,jj,cb)%numc = (Nmid/Nnorm) * ( intgN(zs1,N01,zvcint,Vmid) + intgN(zs2,N02,Vmid,Vhi) )
                ! For different species, assume the mass distribution identical in particles within the bin
                zactvol = (Nmid*Vmid/Vnorm) * ( intgV(zs1,N01,zvcint,Vmid) + intgV(zs2,N02,Vmid,Vhi) )
                DO ss = 1,8
                   pactd(ii,jj,cb)%volc(ss) = zactvol*( paero(ii,jj,ab)%volc(ss)/(Vmid*Nmid) )
                END DO

             ELSE IF (zvcstar >= Vmid) THEN
                
                ! Use actual critical volume only in the critical bin, otherwise current bin limits
                zvcint = MIN(zvcstar,Vhi)

                pactd(ii,jj,cb)%numc = (Nmid/Nnorm) * ( intgN(zs2,N02,zvcint,Vhi) )
                zactvol = (Nmid*Vmid/Vnorm) * ( intgV(zs2,N02,zvcint,Vhi) )
                DO ss = 1,8
                   pactd(ii,jj,cb)%volc(ss) = zactvol*( paero(ii,jj,ab)%volc(ss)/(Vmid*Nmid) )
                END DO

             END IF

             IF ( pactd(ii,jj,cb)%numc < 0._dp .OR. ANY(pactd(ii,jj,cb)%volc(:) < 0._dp) ) WRITE(*,*) 'HEPHEPHPE'

             pactd(ii,jj,cb)%numc = MAX(0._dp, pactd(ii,jj,cb)%numc)
             DO ss = 1,8
                pactd(ii,jj,cb)%volc(ss) = MAX(0._dp, pactd(ii,jj,cb)%volc(ss))
             END DO

             ! "Artificially" adjust the wet size of newly activated a little bit to prevent them from being
             ! evaporated immediately

             pactd(ii,jj,cb)%volc(8) = pactd(ii,jj,cb)%numc*pi6*(pdcrit(ii,jj,ab)**3) *  &
!                                       MIN(2._dp,(3.e-6_dp/pdcrit(ii,jj,ab))**2)
                                       (3.e-6_dp/max(3.e-6_dp/sqrt(2._dp),pdcrit(ii,jj,ab)))**2
          END DO ! cb

       END DO ! ii
    END DO ! jj

  END SUBROUTINE activate3
  ! ------------------------------------------------
  REAL(dp) FUNCTION intgN(ikk,icc,ilow,ihigh)
    ! Gets the integral over a (linear) number concentration distribution
    ! 
    USE mo_kind, ONLY : dp
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: ikk,icc,ilow,ihigh
    intgN = 0.5_dp*ikk*MAX(ihigh**2 - ilow**2,0._dp) + icc*MAX(ihigh - ilow,0._dp)
  END FUNCTION intgN
  ! ------------------------------------------------
  REAL(dp) FUNCTION intgV(ikk,icc,ilow,ihigh)
    ! Gets the integral over a volume volume distribution based on a linear
    ! number concentration distribution
    USE mo_kind, ONLY : dp
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: ikk,icc,ilow,ihigh
    intgV = (1._dp/3._dp)*ikk*MAX(ihigh**3 - ilow**3,0._dp) + 0.5_dp*icc*MAX(ihigh**2 - ilow**2,0._dp)
  END FUNCTION intgV




  !----------------------------------------- 
  SUBROUTINE autoconv2(kproma,kbdim,klev,   &
                      pcloud,pprecp         )
  ! 
  ! Uses a more straightforward method for converting cloud droplets to drizzle.
  ! Assume a lognormal cloud droplet distribution for each bin. Sigma_g is an adjustable
  ! parameter and is set to 1.2 by default
  ! 
    USE mo_kind, ONLY : dp
    USE mo_submctl, ONLY : t_section,   &
                               ica,fca,     &
                               icb,fcb,     &
                               ira,fra,     &
                               ncld,        &
                               nprc,        &
                               rhowa,       &
                               pi6,         &
                               nlim
    USE mo_constants, ONLY : rd
    IMPLICIT NONE

    INTEGER, INTENT(in) :: kproma,kbdim,klev
    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld)
    TYPE(t_section), INTENT(inout) :: pprecp(kbdim,klev,nprc)

    REAL(dp) :: Vrem, Nrem, Vtot, Ntot
    REAL(dp) :: dvg,dg

    REAL(dp), PARAMETER :: zd0 = 50.e-6_dp
    REAL(dp), PARAMETER :: sigmag = 1.2_dp
 
    INTEGER :: ii,jj,cc,ss

    ! Find the cloud bins where the mean droplet diameter is above 50 um
    ! Do some fitting...
    DO jj = 1,klev
       DO ii = 1,kproma
          DO cc = 1,ncld

             Vrem = 0._dp
             Nrem = 0._dp
             Ntot = 0._dp
             Vtot = 0._dp
             dvg = 0._dp
             dg = 0._dp

             Ntot = pcloud(ii,jj,cc)%numc
             Vtot = SUM(pcloud(ii,jj,cc)%volc(:))

             IF ( Ntot > 0._dp .AND. Vtot > 0._dp ) THEN
                
                ! Volume geometric mean diameter
                dvg = ((Vtot/Ntot/pi6)**(1._dp/3._dp))*EXP( (3._dp*LOG(sigmag)**2)/2._dp )
                dg = dvg*EXP( -3._dp*LOG(sigmag)**2 )

                !testi = cumlognorm(2.19e-4_dp,sigmag,zd0)
                !WRITE(*,*) 'TESTAAN: ',testi

                Vrem = Vtot*( 1._dp - cumlognorm(dvg,sigmag,zd0) )
                Nrem = Ntot*( 1._dp - cumlognorm(dg,sigmag,zd0) )



                IF ( Vrem < 0._dp ) WRITE(*,*) 'ERROR Vrem < 0', Vrem, Vtot
                IF ( Vrem > Vtot ) WRITE(*,*) 'ERROR Vrem > Vtot', cumlognorm(dvg,sigmag,zd0),dvg,dg,Vtot
                IF ( Nrem < 0._dp ) WRITE(*,*) 'ERROR Nrem < 0', Nrem
                IF ( Nrem > Ntot ) WRITE(*,*) 'ERROR Nrem > Ntot', cumlognorm(dg,sigmag,zd0),dvg,dg,Ntot

                IF ( Vrem > 0._dp .AND. Nrem > 0._dp) THEN
                   
                   ! Put the mass and number to the first precipitation bin and remover from 
                   ! cloud droplets
                   DO ss = 1,7
                      pprecp(ii,jj,1)%volc(ss) = pprecp(ii,jj,1)%volc(ss) + pcloud(ii,jj,cc)%volc(ss)*(Nrem/Ntot)
                      pcloud(ii,jj,cc)%volc(ss) = pcloud(ii,jj,cc)%volc(ss)*(1._dp - (Nrem/Ntot))
                   END DO
                   pprecp(ii,jj,1)%volc(8) = pprecp(ii,jj,1)%volc(8) + pcloud(ii,jj,cc)%volc(8)*(Vrem/Vtot)
                   pcloud(ii,jj,cc)%volc(8) = pcloud(ii,jj,cc)%volc(8)*(1_dp - (Vrem/Vtot))

                   pprecp(ii,jj,1)%numc = pprecp(ii,jj,1)%numc + Nrem
                   pcloud(ii,jj,cc)%numc = pcloud(ii,jj,cc)%numc - Nrem
                
                END IF ! Nrem Vrem

             END IF ! Ntot Vtot

          END DO ! cc
       END DO ! ii
    END DO ! jj

  END SUBROUTINE autoconv2
  !
  ! -----------------------------------------------------------------
  !
  REAL(dp) FUNCTION cumlognorm(dg,sigmag,dpart)
    USE mo_kind, ONLY : dp
    USE mo_submctl, ONLY : pi
    IMPLICIT NONE
    ! Cumulative lognormal function
    REAL(dp), INTENT(in) :: dg
    REAL(dp), INTENT(in) :: sigmag
    REAL(dp), INTENT(in) :: dpart

    REAL(dp) :: hlp1,hlp2

    !cumlognorm = 0._dp

    hlp1 = ( LOG(dpart) - LOG(dg) )
    hlp2 = SQRT(2._dp)*LOG(sigmag)
    cumlognorm = 0.5_dp + 0.5_dp*ERF( hlp1/hlp2 )

  END FUNCTION cumlognorm
  !
  ! ----------------------------------------------------------------
  !
  REAL(dp) FUNCTION errf(x)
    USE mo_kind, ONLY : dp
    IMPLICIT NONE
    ! (Approximative) Error function.
    ! This is available as an intrinsic function as well but the implementation is somewhat compiler-specific
    ! so an approximation is given here explicitly...
    
    ! EI TOIMI OIKEIN ÄÄRIPÄISSÄ


    REAL(dp), INTENT(in) :: x
    REAL(dp) :: hlp

    hlp = 1._dp + 0.0705230784_dp*x + 0.0422820123_dp*(x**2) +  &
          0.0092705272_dp*(x**3) + 0.0001520143_dp*(x**4) +       &
          0.0002765672_dp*(x**5) + 0.0000430638_dp*(x**6)

    errf = 1._dp - 1._dp/(hlp**16)

    WRITE(*,*) 'TESTAAN 2: ', x, errf

  END FUNCTION errf




END MODULE mo_salsa_cloud
