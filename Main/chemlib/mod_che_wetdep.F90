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
!
module mod_che_wetdep

! This module contains routines for wet scavenging of gas and aerosols
 
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_che_common
  use mod_che_indices
  use mod_che_drydep

  public :: sethet , wetdepa

  contains
!
! Calculate the total wet deposition rate for each species that is able to be
! washed/rained out (1/s).  This module was originally from MOZART4 but has been
! adapted to work within the RegCM4.1-chem framework.
! note: the press array is in pascals and must be
! mutiplied by 10 to yield dynes/cm**2.
! set wet deposition for
!           1. h2o2         2. hno3
!           3. ch2o         4. ch3ooh
!           5. pooh         6. ch3coooh
!           7. ho2no2       8. onit
!           9. mvk         10. macr
!          11. c2h5ooh     12. c3h7ooh
!          13. rooh        14. ch3cocho
!          15. pb          16. macrooh
!          17. xooh        18. onitr
!          19. isopooh     20. ch3oh
!          21. c2h5oh      22. glyald
!          23. hyac        24. hydrald
!          25. ch3cho      26. isopno3
!
  subroutine sethet(j,zmid1,phis1,tfld1,cmfdqr1,nrain1, &
                    nevapr1,delt,xhnm1,qin1,ps2)
    implicit none
!
    integer, intent(in) :: j               ! longitude index
    ! surf geopotential (m2/s2)
    real(dp), dimension(:) , intent(in) :: phis1
    ! midpoint geopot convert (m) to (km)
    real(dp), dimension(:,:) , intent(in) :: zmid1
    ! temperature (K)
    real(dp), dimension(:,:) , intent(in) :: tfld1
    ! dq/dt for convection (kg/kg/s) convert to (1/s) below
    real(dp), dimension(:,:) , intent(in) :: cmfdqr1
    ! stratoform precip (1/s)
    real(dp), dimension(:,:) , intent(in) :: nrain1
    ! evaporation (kg/kg/s) convert to (1/s) below
    real(dp), dimension(:,:) , intent(in) :: nevapr1
    ! time step ( s )
    real(dp), intent(in) :: delt
    ! total atms density (kg/m3) convert to (#/cm^3) below
    real(dp), dimension(:,:) , intent(in) :: xhnm1
    ! exported species ( mmr )
    real(dp), dimension(:,:,:) , intent(in) :: qin1
    ! Pressure ( cb )
    real(dp), dimension(:) , intent(in) :: ps2
!
    real(dp) :: phis(ici1:ici2)        ! surf geopotential (m2/s2)
    real(dp) :: zmid(ici1:ici2,kz)     ! midpoint geopot convert (m) to (km)
    real(dp) :: tfld(ici1:ici2,kz)     ! temperature (K)
    real(dp) :: cmfdqr(ici1:ici2,kz)   ! dq/dt for convection (1/s) 
    real(dp) :: nrain(ici1:ici2,kz)    ! stratoform precip (1/s)
    real(dp) :: nevapr(ici1:ici2,kz)   ! evaporation (1/s)
    real(dp) :: xhnm(ici1:ici2,kz)     ! total atms density (#/cm^3)
    real(dp) :: qin(ici1:ici2,kz,ntr)  ! exported species ( mmr )
    real(dp) :: temp_dep(ici1:ici2)    ! temp var of wet dep rate
                                       ! (centibars_tracer/sec)
    real(dp) :: vmr_hno3(ici1:ici2,kz) ! volume mixing ratio
    real(dp) :: vmr_h2o2(ici1:ici2,kz) ! volume mixing ratio
    ! Effective henry's law constant (1/s)
    real(dp) :: het_rates(ici1:ici2,kz,ntr)

    ! mean diameter of rain drop (cm)
    real(dp) , parameter :: xrm = 0.189D0
    ! mean rain drop terminal velocity (cm/s)
    real(dp) , parameter :: xum = 748.D0
    ! kinetic viscosity (cm^2/s)
    real(dp) , parameter :: xvv = 6.18D-2
    ! mass transport coefficient (cm/s)
    real(dp) , parameter :: xdg = 0.112D0
    ! reference temperature (k)
    real(dp) , parameter :: t0 = 298.0D0
    real(dp) , parameter :: xph0 = 1.D-5 ! cloud [h+]
    ! saturation factor for hno3 in clouds 
    real(dp) , parameter :: satf_hno3 = 0.016D0
    ! saturation factor for hno3 in clouds 
    real(dp) , parameter :: satf_h2o2 = 0.016D0
    ! saturation factor for hno3 in clouds 
    real(dp) , parameter :: satf_ch2o = 0.1D0
    ! (atmospheres/deg k/cm^3)
    real(dp) , parameter :: const0 = boltzk * 1.D-6
    ! hno3 dissociation constant
    real(dp) , parameter :: hno3_diss = 15.4D0
    ! geometry factor (surf area/volume = geo_fac/diameter)
    real(dp) , parameter :: geo_fac = 6.0D0
    ! mass of background atmosphere (amu)
    real(dp) , parameter :: mass_air = 29.0D0
    ! mass of water vapor (amu)
    real(dp) , parameter :: mass_h2o = 18.0D0
    real(dp) , parameter :: h2o_mol = 1.D3/mass_h2o   ! (gm/mol water)
    real(dp) , parameter :: km2cm = 1.D5              ! convert km to cm
    real(dp) , parameter :: m2km = 1.D-3              ! convert m to km
    real(dp) , parameter :: cm3_2_m3 = 1.D-6          ! convert cm^3 to m^3
    real(dp) , parameter :: m3_2_cm3 = 1.D6           ! convert m^3 to cm^3
    real(dp) , parameter :: liter_per_gram = 1.D-3
    ! (liter/gm/mol*(m/cm)^3)
    real(dp), parameter :: avo2 = navgdr * liter_per_gram * cm3_2_m3

    integer :: ktop   ! index of top model layer that can have clouds
    integer :: i , k , kk , itr ! indicies
    real(dp) :: xkgm         ! mass flux on rain drop
    real(dp) :: all1 , all2  ! work variables
    real(dp) :: stay         ! fraction of layer traversed by falling drop
                             ! in timestep delt
    real(dp) :: xeqca1 , xeqca2 , xca1 , xca2 , xdtm
    real(dp) :: xxx1 , xxx2 , yhno3 , yh2o2
    real(dp) , dimension(ici1:ici2) :: xk0 , work1 , work2 , work3 , zsurf
    real(dp) , dimension(kz) :: xgas1, xgas2
    real(dp) , dimension(ici1:ici2) :: tmp0_rates , tmp1_rates
    ! layer depth about interfaces (cm)
    real(dp) , dimension(ici1:ici2,kz) :: delz
    ! hno3 concentration (molecules/cm^3)
    real(dp) , dimension(ici1:ici2,kz) :: xhno3
    ! h2o2 concentration (molecules/cm^3)
    real(dp) , dimension(ici1:ici2,kz) :: xh2o2
    ! liquid rain water content in a grid cell (gm/m^3)
    real(dp) , dimension(ici1:ici2,kz) :: xliq
    ! conversion rate of water vapor into rain water (molecules/cm^3/s)
    real(dp) , dimension(ici1:ici2,kz) :: rain
    real(dp) , dimension(ici1:ici2,kz) :: xhen_hno3 , xhen_h2o2 , xhen_ch2o , &
        xhen_ch3ooh , xhen_ch3co3h , xhen_ch3cocho , xhen_xooh , xhen_onitr , &
        xhen_ho2no2 , xhen_glyald , xhen_ch3cho , xhen_mvk , xhen_macr
    real(dp) , dimension(ici1:ici2,kz) :: xhen_nh3 , xhen_ch3cooh
    real(dp) , dimension(ici1:ici2,kz,2) :: tmp_hetrates
    real(dp) , dimension(ici1:ici2,kz) :: precip

    !-----------------------------------------------------------------
    !  Initialize rainout array
    !-----------------------------------------------------------------
    het_rates = d_zero
    temp_dep  = d_zero
    vmr_hno3  = d_zero
    vmr_h2o2  = d_zero

    !-----------------------------------------------------------------
    !  Set maximum cloud top to highest model level
    !  (e.g. clouds can exist anywhere)
    !-----------------------------------------------------------------
    ktop = 0

    !-----------------------------------------------------------------
    !  Perform a few conversion for xhnm, zmid, nevap, cmfdqr
    !-----------------------------------------------------------------
    ![kg_h2o/kg_air/s -> 1/s]
    cmfdqr = cmfdqr1(ici1:ici2,:)*(xhnm1(ici1:ici2,:)/rhoh2o)
    ![kg_h2o/kg_air/s -> 1/s]
    nevapr = nevapr1(ici1:ici2,:)*(xhnm1(ici1:ici2,:)/rhoh2o)
    ![m     -> km]
    zmid   = zmid1(ici1:ici2,:)*m2km
    ![kg/m3 -> #/cm3]
    xhnm   = xhnm1(ici1:ici2,:)*navgdr*cm3_2_m3/( amd*d_1000)

    phis   = phis1(ici1:ici2)
    qin    = qin1(ici1:ici2,:,:)
    tfld   = tfld1(ici1:ici2,:)
    nrain  = nrain1(ici1:ici2,:)

    !-----------------------------------------------------------------
    !	... the 2 and .6 multipliers are from a formula by frossling (1938)
    !-----------------------------------------------------------------
    xkgm = xdg/xrm * d_two + xdg/xrm * 0.6D0 * &
           dsqrt( xrm*xum/xvv ) * (xvv/xdg)**(d_one/d_three) 
    do k = ktop + 1 , kz
      precip(:,k) = nrain(:,k)        ! = cmfdqr(:,k) + nrain(:,k) - nevapr(:,k)
      rain(:,k)   = mass_air*precip(:,k)*xhnm(:,k) / mass_h2o
      xliq(:,k)   = precip(:,k) * delt * xhnm(:,k) / navgdr*mass_air * m3_2_cm3
      !---------------------------------------
      ! convert from cb_X to mass mixing ratio
      !---------------------------------------
      if ( ihno3 > 0 ) vmr_hno3(:,k) = (qin(:,k,ihno3)/ps2(:)) * (amd/63.012D0)
      if ( ih2o2 > 0 ) vmr_h2o2(:,k) = (qin(:,k,ih2o2)/ps2(:)) * (amd/34.0147D0)
    end do 
!
    if ( ihno3 > 0 ) then
      xhno3(:,:)  = vmr_hno3(:,:) * xhnm(:,:)
    else
      xhno3(:,:)  = d_zero
    end if
    if ( ih2o2 > 0 ) then
      xh2o2(:,:)  = vmr_h2o2(:,:) * xhnm(:,:)
    else
      xh2o2(:,:)  = d_zero
    end if

    zsurf(:) = m2km * phis(:) / egrav
    do k = ktop + 1 , kz - 1
      delz(:,k) = dabs( (zmid(:,k) - zmid(:,k+1))*km2cm ) 
    end do 
    delz(:,kz) = dabs( (zmid(:,kz) - zsurf(:) )*km2cm ) 

    !-----------------------------------------------------------------
    ! ... part 0b,  for temperature dependent of henrys
    !    xxhe1 = henry con for hno3
    !    xxhe2 = henry con for h2o2
    !lwh 10/00 -- take henry''s law constants from brasseur et al. [1999],
    !    appendix j. for hno3, also consider dissociation to
    !    get effective henry''s law constant; equilibrium
    !    constant for dissociation from brasseur et al. [1999],
    !    appendix k. assume ph=5 (set as xph0 above).
    !    heff = h*k/[h+] for hno3 (complete dissociation)
    !    heff = h for h2o2 (no dissociation)
    !    heff = h * (1 + k/[h+]) (in general)
    !-----------------------------------------------------------------
    do k = ktop + 1 , kz
      work1(:) = (t0 - tfld(:,k))/(t0*tfld(:,k))
    !-----------------------------------------------------------------
    ! 	... effective henry''s law constants:
    !	   hno3, h2o2, ch2o, ch3ooh, ch3coooh (brasseur et al., 1999)
    !    xooh, onitr, macrooh (j.-f. muller; brocheton, 1999)
    !    isopooh (equal to hno3, as for macrooh)
    !    ho2no2 (mozart-1)
    !    ch3cocho, hoch2cho (betterton and hoffman, environ. sci. technol., 1988)
    !    ch3cho (staudinger and roberts, crit. rev. sci. technol., 1996)
    !    mvk, macr (allen et al., environ. toxicol. chem., 1998)
    !-----------------------------------------------------------------
      xk0(:)             = 2.1D5 * dexp(8700.0D0*work1(:))
      xhen_hno3(:,k)     = xk0(:) * (d_one + hno3_diss / xph0)
      xhen_h2o2(:,k)     = 7.45D4 * dexp(6620.0D0*work1(:))
      xhen_ch2o(:,k)     = 6.3D3 * dexp(6460.0D0*work1(:))
      xhen_ch3ooh(:,k)   = 2.27D2 * dexp(5610.0D0*work1(:))
      xhen_ch3co3h(:,k)  = 4.73D2 * dexp(6170.0D0*work1(:))
      xhen_ch3cocho(:,k) = 3.70D3 * dexp(7275.0D0*work1(:))
      xhen_xooh(:,k)     = 90.5 * dexp(5607.0D0*work1(:))
      xhen_onitr(:,k)    = 7.51D3 * dexp(6485.0D0*work1(:))
      xhen_ho2no2(:,k)   = 2.0D4
      xhen_glyald(:,k)   = 4.1D4 * dexp(4600.0D0*work1(:))
      xhen_ch3cho(:,k)   = 1.4D1 * dexp(5600.0D0*work1(:))
      xhen_mvk(:,k)      = 21.0D0 * dexp(7800.0D0*work1(:))
      xhen_macr(:,k)     = 4.3D0 * dexp(5300.0D0*work1(:))
      xhen_ch3cooh(:,k)  = 4.1D3 * dexp(6300.0D0*work1(:))
      xhen_nh3 (:,k)     = 1.D6
      tmp_hetrates(:,k,:) = d_zero
    end do
   
    !-----------------------------------------------------------------
    !   part 1, solve for high henry constant ( hno3, h2o2)
    !-----------------------------------------------------------------
    long_loop : &
    do i = ici1 , ici2
      xgas1(:) = xhno3(i,:)  ! xgas will change during 
      xgas2(:) = xh2o2(i,:)  ! different levels wash 
    level_loop1 : &
      do kk = ktop + 1 , kz
        stay = d_one
        if ( dabs(rain(i,kk)) > d_zero ) then  ! finding rain cloud
          all1 = d_zero ! accumulation to justisfy saturation
          all2 = d_zero
          stay = ((zmid(i,kk) - zsurf(i))*km2cm)/(xum*delt)
          stay = dmin1(stay,d_one)
          !-----------------------------------------------------------
          !  calculate the saturation concentration eqca
          !-----------------------------------------------------------
          do k = kk , kz ! cal washout below cloud
            xeqca1 =  xgas1(k) / (xliq(i,kk)*avo2 + &
              d_one/(xhen_hno3(i,k)*const0*tfld(i,k))) * xliq(i,kk)*avo2
            xeqca2 =  xgas2(k) / (xliq(i,kk)*avo2 + &
              d_one/(xhen_h2o2(i,k)*const0*tfld(i,k))) * xliq(i,kk)*avo2
            !---------------------------------------------------------
            ! calculate ca; inside cloud concentration in #/cm3(air)
            !---------------------------------------------------------
            xca1 = geo_fac*xkgm*xgas1(k)/(xrm*xum)*delz(i,k)*xliq(i,kk)*cm3_2_m3
            xca2 = geo_fac*xkgm*xgas2(k)/(xrm*xum)*delz(i,k)*xliq(i,kk)*cm3_2_m3
            !---------------------------------------------------------
            ! if is not saturated
            !    hno3(gas)_new = hno3(gas)_old - hno3(h2o)
            ! otherwise
            !    hno3(gas)_new = hno3(gas)_old
            !---------------------------------------------------------
            all1 = all1 + xca1
            all2 = all2 + xca2
            if ( all1 < xeqca1 ) then
              xgas1(k) = dmax1(xgas1(k)-xca1,d_zero)
            end if
            if ( all2 < xeqca2 ) then
              xgas2(k) = dmax1(xgas2(k)-xca2,d_zero)
            end if
          end do
        end if
        !-----------------------------------------------------------------
        ! calculate the lifetime of washout (second)
        !  after all layers washout the concentration of hno3 is reduced 
        !  then the lifetime xtt is calculated by
        !
        !  xtt = (xhno3(ini) - xgas1(new))/(dt*xhno3(ini))
        !  where dt = passing time (s) in vertical path below the cloud
        !  dt = dz(cm)/um(cm/s)
        !-----------------------------------------------------------------
        xdtm = delz(i,kk) / xum ! the traveling time in each dz
        xxx1 = (xhno3(i,kk) - xgas1(kk))
        xxx2 = (xh2o2(i,kk) - xgas2(kk))
        if ( dabs(xxx1) > d_zero ) then  ! if no washout lifetime = 1.e29
          yhno3  = xhno3(i,kk)/xxx1 * xdtm
        else
          yhno3  = 1.D29
        end if
        if ( dabs(xxx2) > d_zero ) then  ! if no washout lifetime = 1.e29
          yh2o2  = xh2o2(i,kk)/xxx2 * xdtm
        else
          yh2o2  = 1.D29 
        end if 
        tmp_hetrates(i,kk,1) = dmax1(d_one/yh2o2,d_zero) * stay
        tmp_hetrates(i,kk,2) = dmax1(d_one/yhno3,d_zero) * stay
      end do level_loop1
    end do long_loop

    !-----------------------------------------------------------------
    ! part 2, in-cloud solve for low henry constant hno3 and h2o2 have
    ! both in and under cloud
    !   *** abt go directly to calculating the tracer tendency (chiten)
    !   *** MOZART code originally used het_rates and passed it out
    ! old way MOZART: het_rates(:,k,ihcho)  = work3(:)
    ! for RegCM4.1  : chiten(:,k,j,ihcho)  = chiten + work3(:) * chib
    !-----------------------------------------------------------------
    level_loop2 : &
    do k = ktop + 1 , kz
      work1(:) = avo2 * xliq(:,k)
      work2(:) = const0 * tfld(:,k)
      work3(:) = dmax1(rain(:,k)/(h2o_mol*(work1(:) + &
           d_one/(xhen_ch2o(:,k)*work2(:)))),d_zero) * satf_ch2o
      if ( ihcho > 0 ) then
        het_rates(:,k,ihcho)  = work3(:)
      end if
      if ( iisopno3 > 0 ) then
        het_rates(:,k,iisopno3) = work3(:)
      end if
      if ( ihyac > 0 ) then
        het_rates(:,k,ihyac) = work3(:)
      end if
      if ( ihydrald > 0 ) then
        het_rates(:,k,ihydrald) = work3(:)
      end if
      work3(:) = dmax1(rain(:,k)/(h2o_mol*(work1(:) + &
                 d_one/(xhen_ch3ooh(:,k)*work2(:)))),d_zero)
      if ( ich3ooh > 0 ) then
        het_rates(:,k,ich3ooh)  = work3(:)
      end if
      if ( ipooh > 0 ) then
        het_rates(:,k,ipooh)  = work3(:)
      end if
      if ( ic2h5ooh > 0 ) then
        het_rates(:,k,ic2h5ooh) = work3(:)
      end if
      if ( ic3h7ooh > 0 ) then
        het_rates(:,k,ic3h7ooh) = work3(:)
      end if
      if ( irooh > 0 ) then
        het_rates(:,k,irooh) = work3(:)
      end if
      if ( imoh > 0 ) then
        het_rates(:,k,imoh) = work3(:)
      end if
      if ( ieoh > 0 ) then
        het_rates(:,k,ieoh) = work3(:)
      end if

      if ( ich3coooh > 0 ) then
        het_rates(:,k,ich3coooh) = dmax1(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_ch3co3h(:,k)*work2(:)))),d_zero)
      end if
      if ( ihno4 > 0 ) then
        het_rates(:,k,ihno4) = dmax1(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_ho2no2(:,k)*work2(:)))),d_zero)
      end if
      if ( ich3cocho > 0 ) then
        het_rates(:,k,ich3cocho) = dmax1(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_ch3cocho(:,k)*work2(:)))),d_zero)
      end if
      if ( ixooh > 0 ) then
        het_rates(:,k,ixooh) = dmax1(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_xooh(:,k)*work2(:)))),d_zero)
      end if
      if ( ionitr > 0 ) then
        het_rates(:,k,ionitr) = dmax1(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_onitr(:,k)*work2(:)))),d_zero)
      end if
      if ( iglyald > 0 ) then
        het_rates(:,k,iglyald) = dmax1(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_glyald(:,k)*work2(:)))),d_zero)
      end if
      if ( iald2 > 0 ) then
        het_rates(:,k,iald2) = dmax1(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_ch3cho(:,k)*work2(:)))),d_zero)
      end if
      if ( imvk > 0 ) then
        het_rates(:,k,imvk) = dmax1(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_mvk(:,k)*work2(:)))),d_zero)
      end if
      if ( imacr > 0 ) then
        het_rates(:,k,imacr) = dmax1(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_macr(:,k)*work2(:)))),d_zero)
      end if
      if ( ih2o2 > 0 ) then
        work3(:) = satf_h2o2 * dmax1(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_h2o2(:,k)*work2(:)))),d_zero)
        het_rates(:,k,ih2o2) =  work3(:) + tmp_hetrates(:,k,1)
      end if

      work3(:) = tmp_hetrates(:,k,2) + satf_hno3 * &
           dmax1(rain(:,k)/(h2o_mol*(work1(:) +    &
           d_one/(xhen_hno3(:,k)*work2(:)))),d_zero)
      tmp0_rates(:) = work3(:)
      tmp1_rates(:) = 0.2D0*work3(:)
      if ( ihno3 > 0 ) then
        het_rates(:,k,ihno3) = work3(:)
      end if
      ! ONIT is organic nitrate
      if ( ionit > 0 ) then
        het_rates(:,k,ionit) = work3(:)
      end if
       !**** PB is lead
      if ( ipb > 0 ) then
        het_rates(:,k,ipb) = work3(:)
      end if
      !  MARCROOH is oxidation product of methacrolein
      if ( imacrooh > 0 ) then
        het_rates(:,k,imacrooh) = work3(:)
      end if
      ! ISOPOOH is hydroxyhydroperoxide (unsaturated)
      if ( iisopooh > 0 ) then
        het_rates(:,k,iisopooh) = work3(:)
      end if
      if ( iso2 > 0 ) then
        het_rates(:,k,iso2) = het_rates(:,k,ih2o2)
      end if
      ! SOA is Secondary Organic aerosols
      if ( isoa > 0 ) then
         het_rates(:,k,isoa) = tmp1_rates(:)
      end if
      ! OC2 is organic carbon, hydrophylic
!     if ( ioc2 > 0 ) then
!	      het_rates(:,k,ioc2) = tmp1_rates(:)
!	    end if
      ! CB2 is black carbon, hydrophylic
!     if ( icb2 > 0 ) then
!	      het_rates(:,k,icb2) = tmp1_rates(:)
!	    end if
      if ( iso4 > 0 ) then
        het_rates(:,k,iso4) = tmp1_rates(:)
      end if
      ! SA1 is sea salt 0.1-0.5 microns
!     if ( sa1_ndx > 0 ) then
!	      het_rates(:,k,sa1_ndx) = tmp1_rates(:)
!	    end if
      ! SA2 is sea salt 0.5-1.5 microns
!     if ( sa2_ndx > 0 ) then
!	      het_rates(:,k,sa2_ndx) = tmp1_rates(:)
!	    end if
      ! SA3 is sea salt 1.5-5 microns
!     if ( sa3_ndx > 0 ) then
!	      het_rates(:,k,sa3_ndx) = tmp1_rates(:)
!	    end if
      ! SA4 is sea salt 5-10 microns
!     if ( sa4_ndx > 0 ) then
!	      het_rates(:,k,sa4) = tmp1_rates(:)
!	    end if
      if ( inh3 > 0 ) then
         het_rates(:,k,inh3) = dmax1(rain(:,k)/(h2o_mol*(work1(:) + &
               d_one/(xhen_nh3(:,k)*work2(:)))),d_zero) * satf_hno3
      end if
      if ( inh4 > 0 ) then
        het_rates(:,k,inh4) = tmp1_rates(:)
      end if
      if ( inh4no3 > 0 ) then
        het_rates(:,k,inh4no3 ) = tmp1_rates(:)
      end if
      ! ALKOOH is C5H11OOH
      if ( ialkooh  > 0 ) then
        het_rates(:,k,ialkooh)  = het_rates(:,k,ich3ooh)
      end if
      ! MEKOOH is CH3COCH(OOH)CH3
      if ( imekooh  > 0 ) then
        het_rates(:,k,imekooh)  = het_rates(:,k,ich3ooh)
      end if
      ! TOLOOH is C7H10O6
      if ( itolooh  > 0 ) then
        het_rates(:,k,itolooh)  = het_rates(:,k,ich3ooh)
      end if
      ! TERPOOH is C10H16(OH)(OOH)
      if ( iterpooh > 0 ) then
        het_rates(:,k,iterpooh) = het_rates(:,k,ich3ooh)
      end if
      ! CH3COOH is acetic acid
      if ( ich3cooh > 0 ) then
        het_rates(:,k,ich3cooh) = dmax1(rain(:,k)/(h2o_mol*(work1(:) + &
                  d_one/(xhen_ch3cooh(:,k)*work2(:)))),d_zero)
      end if
      !-----------------------------------------------------------------
      ! *** abt added to work with RegCM4.1
      ! Directly calculate the tracer tendency (chiten)
      ! Loop over all non-aerosol tracers
      ! old way MOZART: het_rates(:,k,ihcho)  = work3(:)
      ! for RegCM4.1  : chiten(:,k,j,ihcho)  = chiten + work3(:) * chib
      !-----------------------------------------------------------------
      tracer_loop1 : &
      do itr = 1 , ntr
        if (chtrname(itr).ne.'DUST'  .or. chtrname(itr).ne.'BC_HB' .or. &
            chtrname(itr).ne.'SSALT' .or. chtrname(itr).ne.'BC_HL' .or. &
            chtrname(itr).ne.'DMS'   .or. chtrname(itr).ne.'OC_HB' .or. &
            chtrname(itr).ne.'OC_HL' ) then
          temp_dep(:) = (d_one-dexp(-het_rates(:,k,itr)*delt))*qin(:,k,itr)
          chiten(j,ici1:ici2,k,itr) = chiten(j,ici1:ici2,k,itr) - &
                                      temp_dep(:)/delt 
          remcvc(ici1:ici2,k,j,itr) = remcvc(ici1:ici2,k,j,itr) + &
                                      temp_dep(:)/d_two
          remlsc(ici1:ici2,k,j,itr) = remlsc(ici1:ici2,k,j,itr) + &
                                      temp_dep(:)/d_two
        end if
      end do tracer_loop1
    end do level_loop2
  end subroutine sethet

  subroutine wetdepa(j,mbin,indp,beffdiam,rhoaer,t,wl,fracloud,fracum, &
                     pressg,shj,rho,totppt,pdepv)
    implicit none
!
    integer , intent(in) :: j , mbin
    integer , dimension(mbin) , intent(in) :: indp
    real(dp) , dimension(mbin) , intent(in) :: beffdiam
    real(dp) , intent(in) :: rhoaer ! specific aerosol density
    real(dp) , dimension(iy,kz) , intent(in) :: wl , fracloud , fracum , &
                                                totppt , t , rho
    real(dp) , dimension(iy) ,intent(in) :: pressg
    real(dp) , dimension(kz) ,intent(in) :: shj
    real(dp) , dimension(iy,kz,ntr) , intent (in) :: pdepv
    ! size of the aerosol bin
    ! index of the correponding aerosol in the chi table

    real(dp) , dimension(ici1:ici2,kz,mbin) :: colef , wetdep , rhsize , rhop
    real(dp) , dimension(ntr) :: wetrem , wetrem_cvc 
    real(dp) :: wtend
    integer :: n , k , i

    ! rain out parametrisation 
    ! clmin = non-precipitating cloud
    ! conversion threshold, clmin=0.01g/m3
    real(dp) , parameter :: clmin = 0.01D0
    ! remcum= removal rate for cumulus
    ! cloud scavenging (s-1)
    real(dp) , parameter :: remcum = 1.0D-3

    do n = 1 , mbin
      ! wet deposition term
      ! Wet removal at resolvable scale (fcc)
      ! add non-precipitating cloud conversion (threshold clmin=0.01g/m3)
      ! the same as that in subroutine exmois clmin = 0.01
      ! FAB: 
      ! Care, remrat and rembc as calculated in precip are already grid
      ! level average removal rates
      ! here remrat is divided by fracloud (large scale cloud fraction)
      ! to get the correct in cloud removal rate 
      if ( ichremlsc == 1 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            if ( wl(i,k) > clmin ) then
              wetrem(indp(n)) = d_zero
              if ( cremrat(j,i,k) > d_zero ) then
                wetrem(indp(n)) = fracloud(i,k)*chtrsol(indp(n)) * &
                   chib(j,i,k,indp(n))* &
                   (dexp(-cremrat(j,i,k)/fracloud(i,k)*dtche)-d_one)
                chiten(j,i,k,indp(n)) = chiten(j,i,k,indp(n)) + &
                   wetrem(indp(n))/dtche
                remlsc(i,k,j,indp(n)) = remlsc(i,k,j,indp(n)) - &
                   wetrem(indp(n))/d_two
              end if
            end if
          end do
        end do
      end if
   
      if ( ichremcvc == 1 ) then
        ! sub-scale wet removal, cumulus cloud (fracum)
        ! remcum = in cloud removal rate for cumulus cloud scavenging (s-1)
        ! remcum = 1.e-3
        do i = ici1 , ici2
          if ( kcumtop(j,i) >  0 ) then
            do k = kcumtop(j,i) , kz
              wetrem_cvc(indp(n)) = fracum(i,k)*chtrsol(indp(n)) * &
                   chib(j,i,k,indp(n))*(dexp(-remcum*dtche)-d_one)
              chiten(j,i,k,indp(n)) = chiten(j,i,k,indp(n)) + &
                   wetrem_cvc(indp(n))/dtche
              remcvc(i,k,j,indp(n)) = remcvc(i,k,j,indp(n)) - &
                   wetrem_cvc(indp(n))/d_two
            end do
          end if
        end do
      end if
    end do ! end tracer loop

    ! below cloud scavenging for aerosol:
    ! calculate now the below cloud (washout scavenging rate) for
    ! aerosol  WETDEP (in s-1)
    ! calculate the effective diameter from beffdiam, correction for
    ! humidity not done yet.
    ! take particule effective diameter as the average bin
    ! kept as dry diameter . multiply by 0.5 since the radius is
    ! used in collection efficiency calculation.
    do i = ici1 , ici2
      do k = 1 , kz
        rhsize(i,k,1:mbin) = d_half * 1.D-06 * beffdiam(1:mbin)
      end do
    end do
    ! dry density for now
    rhop(:,:,:) = rhoaer
   
    call blcld(mbin,indp,rhsize,t,pressg,shj,rho,totppt,pdepv,rhop,wetdep,colef)

    ! calculate the tendency due to wahsout dep here.

    do n = 1 , mbin
      do k = 1 , kz
        do i = ici1 , ici2
          wtend = chib(j,i,k,indp(n))*(d_one-dexp(-wetdep(i,k,n)*dtche))/dtche
          chiten(j,i,k,indp(n)) = chiten(j,i,k,indp(n)) - wtend 
          ! wet deposition diagnostic, adding to rainout contribution ! 
          ! nod differenciation between conv and large scale yet
          remcvc(i,k,j,indp(n)) = remcvc(i,k,j,indp(n)) + wtend*dtche/d_two
          remlsc(i,k,j,indp(n)) = remlsc(i,k,j,indp(n)) + wtend*dtche/d_two
        end do
      end do
    end do

  end subroutine wetdepa

  subroutine blcld(mbin,indp,rhsize,t,pressg,shj,rho,totppt,pdepv, &
                   rhop,wetdep,colef)
    implicit none 
 
    integer , intent(in) :: mbin
    real(dp) , dimension(ici1:ici2,kz,mbin) , intent(in) :: rhsize , rhop
    real(dp) , dimension(iy,kz) , intent(in) :: t , rho , totppt
    real(dp) , dimension(iy,kz,ntr) , intent(in) :: pdepv ! care ,ntr dimension
    real(dp) , dimension(iy) , intent(in) :: pressg
    real(dp) , dimension(kz) , intent(in) :: shj
    real(dp) , dimension(ici1:ici2,kz,mbin) , intent(out) :: colef , wetdep 

    ! index of the correponding aerosol in the chi table
    integer , dimension(mbin) , intent(in) :: indp

    real(dp) :: dm , tl , rrm
    integer :: i , n , k

    real(dp) , parameter :: bcrain = 0.5D0
    real(dp) , parameter :: bcsnow = 0.8D0

    !----------------------------------------------------------------------c
    ! call to compute collection efficiency coefficients                   c
    !----------------------------------------------------------------------c 
 
    call cas(t,colef,indp,rhop,rho,rhsize,pressg,mbin,totppt,pdepv,shj) 

    wetdep(:,:,:) = d_zero
    do n = 1 , mbin 
      do k = 1 , kz
        do i = ici1 , ici2
          if ( totppt(i,k) > 1.D-20 ) then 
            tl = t(i,k)-tzero
            !----------------------------------------------------------c
            ! rain scavenging rate : wetdep in s-1   ?                 c
            ! units of totppt (mm.s-1) converted to m.s-1 (by 1.0e-3)  c
            !----------------------------------------------------------c
            if ( tl > d_zero ) then 
              rrm = 0.35D0*(totppt(i,k)*3600.0D0)**0.25D0*1.0D-3 
              ! formulation for monoddipesre rain of radius rrm
              ! (cf seinfeld), gives a scavenging rate wetdep in
              ! s-1 ?? bcrain is 0.5 when it is 1.5 in seinfeld ?? 
              wetdep(i,k,n) = bcrain*totppt(i,k)*1.0D-3*colef(i,k,n)/rrm 
            end if 
            !---------------------------------------------------------c
            ! snow scavenging rate                                    c
            ! * density of snow is set as 1/10 of liquid water.       c
            ! * factor of 1.0e-2 in wetdep calculation takes this and c
            ! * unit                                                  c
            !      change (to m.s-1) into account                     c
            !---------------------------------------------------------c 
            ! . . . . snow scavenging
            if ( tl <= d_zero .and. tl >= -8.0D0 ) then
              dm = 3.8D-5 ! characteristic length scale [m] 
              wetdep(i,k,n) = bcsnow*totppt(i,k)*1.0D-3*colef(i,k,n)/dm 
            end if
            ! . . . . steller snow scavenging 
            if ( tl < -8. .and. tl >= -25.0D0 ) then 
              dm = 2.7D-5
              wetdep(i,k,n) = bcsnow*totppt(i,k)*1.0D-3*colef(i,k,n)/dm 
            end if
            ! . . . . graupel scavenging 
            if ( tl > -25.0D0 ) then 
              dm = 1.4D-4 
              wetdep(i,k,n) = bcsnow*totppt(i,k)*1.0D-3*colef(i,k,n)/dm 
            endif
          end if 
        end do 
      end do 
    end do
  end subroutine blcld ! below-cloud scavenging subroutine 
!
!*********************************************************************** 
!****  a e r o s o l   c o l l e c t i o n   e f f i c i e n c y    ****
!***********************************************************************
! . . . calculation of the collection efficiency of the aerosols of  . .
! . . . type n by collector droplets of radius rcol  . . . . . . . . . . 
! . . . . . . . . . . . . . . . . . . . . . . . . . . . .  . . . . . . . 
! . . . dec 19/96 - s.l.gong - vectorised the whole program & added  . .
! . . . . . . . . . . . . . .  some working space  . . . . . . . . . . . 
! . . . jul 8/94  - s.l.gong - first version . . . . . . . . . . . . . . 
!***********************************************************************
!
  subroutine cas(t,colef,indp,rhop,rho,rhsize,pressg,mbin,totppt,pdepv,shj)
  !------------------------------------------------------------------------
  !  purpose:
  !  --------
  !  this is a module of calculating the collection efficiency of
  !  the aerosols of type n by collector droplets of radius rcol.
  ! 
  implicit none

  integer , intent(in) :: mbin
  integer , dimension(mbin) , intent(in) :: indp 
  real(dp) , dimension(iy,kz) , intent(in) :: t , rho , totppt
  real(dp) , dimension(ici1:ici2,kz,mbin) , intent(in) :: rhop , rhsize
  real(dp) , dimension(iy) , intent(in) :: pressg
  real(dp) , dimension(iy,kz,ntr) , intent(in) :: pdepv
  real(dp) , dimension(kz) , intent(in) :: shj
  ! collection efficiency
  real(dp) , dimension(ici1:ici2,kz,mbin) , intent(out) :: colef
!
  real(dp) :: pdiff
  real(dp) :: amu , anu !dynamic viscosity of air
  real(dp) :: amfp   ! mean molecular free path
  real(dp) :: schm   ! schmidt number
  real(dp) :: prii
  real(dp) :: priiv 
  real(dp) :: cfac
  real(dp) :: re  ! reynolds number
  real(dp) :: rr  ! ratio of collected particle radius:collector particle radius
  real(dp) :: st  ! stokes number of collected particles
  real(dp) :: vr
  real(dp) :: sstar 
  real(dp) :: amob 
  real(dp) :: colimp , pre
  real(dp) :: vpr ! average settling velocity
  real(dp) :: rrm ! mass mean raindrop radius (for rain)
                  ! characteristic capture length (for snow)
  real(dp) :: alpha
  integer :: n , k , i

  ! Cunningham slip correction factor parameters for rain

  real(dp) , parameter :: aa1r = 1.249D0
  real(dp) , parameter :: aa2r = 0.42D0
  real(dp) , parameter :: aa3r = 0.87D0
  real(dp) , parameter :: rhorain = 1000.0D0
  real(dp) , parameter :: amuw = 1.002D-3 ! at 20*c [kg/m/sec]
 
  colef(:,:,:) = d_zero

  do n = 1 , mbin 
    do k = 1 , kz
      do i = ici1 , ici2
        ! . . . for precipitation: 
        vpr = d_zero
        if ( totppt(i,k) > 1.D-15 ) then
          !********************************************************
          !* air's dynamic viscosity                           ****
          !********************************************************
          amu = a1*1.D-8*t(i,k)**a2/(t(i,k)+a3)
          !. . . . mid layer pressure in [pascal].
          pre = pressg(i)*shj(k)
          !********************************************************
          !* mean molecular free path.                         ****
          !*     k.v. beard [1976], j atm. sci., 33            ****
          !********************************************************
          amfp = c1*(amu/c2)*(c3/pre)*(t(i,k)/c4)**d_half
          !********************************************************
          !* cunningham slip correction factor and             ****
          !* relaxation time = vg/grav.                        ****
          !********************************************************
          cfac = d_one + amfp/rhsize(i,k,n) * &
                    (aa1+aa2*dexp(-aa3*rhsize(i,k,n)/amfp))
          !*****************************************************
          !* the schmidt number is the ratio of the         ****
          !* kinematic viscosity of air to the particle     ****
          !* brownian diffusivity ===> sc=v/d               ****
          !*****************************************************
          anu = amu/rho(i,k)
          amob = 6.0D0*mathpi*amu*rhsize(i,k,n)/cfac
          pdiff = boltzk*t(i,k)/amob
          schm = anu/pdiff 
          !----------------------------------------------------c
          ! collector drop size                                c
          !----------------------------------------------------c
          !----------------------------------------------------c
          ! rain scavenging rate                               c
          ! units of totppt (mm.s-1) converted to m.s-1        c
          ! rrm = mass mean raindrop radius                    c
          !----------------------------------------------------c
          if ( t(i,k) > tzero ) then
            rrm = 0.35D0*(totppt(i,k)*3600.0D0)**0.25D0*1.D-3 
            ! what empirical parametization is that ?
            ! needs apparently rain rate in mm/hr 
            prii = 2.0D0/9.0D0*egrav/amu
            priiv = prii*(rhorain-rho(i,k)) 
            ! cunningham settling slip-correction factor
            ! settling velocity
            cfac = d_one+amfp/rrm*(aa1r+aa2r*dexp(-aa3r*rrm/amfp))
            vpr = priiv*rrm**d_two*cfac
          end if 
          !----------------------------------------------------c
          ! snow scavenging                                    c
          ! data from slinn (1984) in atmospheric science &    c
          ! power production, ed. darryl randerson             c
          ! density of snow is set as 1/10 of liquid water.    c
          ! factor of 1.0e-2 in wetdep calculation takes this  c
          ! and unit change (to m.s-1) into account            c
          ! rrm = characteristic capture length                c
          !----------------------------------------------------c
          ! needle-snow scavenging:
          if ( t(i,k) <= tzero .and. t(i,k)>= tzero-8.0D0 ) then
            vpr = 50.0D-2
            rrm = 10.D-6
            alpha = 1.0D0
          end if
          ! . . . steller-snow scavenging:
          if ( t(i,k) < tzero-8.0D0 .and. t(i,k) >= tzero-25.0D0) then
            vpr = 57.0D-2
            rrm = 100.D-6
            alpha = 0.5D0
          end if
          ! . . . graupel scavenging:
          if ( t(i,k) < tzero-25.0D0 ) then
            vpr = 180.0D-2
            rrm = 1000.D-6
            alpha = 2.0D0/3.0D0
          endif

          ! reynolds number: 
          re = rrm*vpr*rho(i,k)/amu
          ! stokes number of collected particles: 
          st = d_two*pdepv(i,k,indp(n))/egrav*(vpr - &
               pdepv(i,k,indp(n)))/(d_two*rrm)
          ! ratio of radius of collected particle to colector drop: 
          rr = rhsize(i,k,n)/rrm
          vr = amuw/amu
          sstar = (1.2D0+(d_one/12.0D0)* &
                  dlog(d_one+re))/(d_one+dlog(d_one+re))
          if ( st > sstar ) then
            colimp = ((st-sstar)/(st-sstar+d_two/d_three))**(d_three/d_two) * &
                    dsqrt(d_1000/rhop(i,k,n))
          else
            colimp = d_zero
          end if
          ! rain scavenging efficiency:
          if ( t(i,k) > tzero ) then
            colef(i,k,n) = d_four/(re*schm)*(d_one+0.4D0*dsqrt(re) * &
                 schm**(d_one/d_three)+0.16D0*dsqrt(re*schm)) +  &
                 d_four*rr*(d_one/vr + (d_one+d_two*dsqrt(re))*rr) + colimp
          else
          ! snow scavenging efficiency: 
          colef(i,k,n) = (d_one/schm)**alpha +  &
                 (d_one-dexp(-(d_one+dsqrt(re))*rr**d_two)) + colimp
          endif
          ! setting the upper-bound for collection efficiency: 
          colef(i,k,n) = dmax1(d_zero,dmin1(d_one,colef(i,k,n)))
        end if
      end do
    end do 
  end do 
  end subroutine cas
!
end module mod_che_wetdep
