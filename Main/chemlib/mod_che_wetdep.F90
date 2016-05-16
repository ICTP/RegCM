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

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_che_common
  use mod_che_indices
  use mod_che_drydep

  implicit none

  private

  public :: sethet , wetdepa

  contains
  !
  ! Calculate the total wet deposition rate for each species that is able to be
  ! washed/rained out (1/s).
  !  This module was originally from MOZART4 but has been adapted to work
  ! within the RegCM4.1-chem framework.
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
  subroutine sethet(j,zmid1,phis,tfld,strappt,convppt, &
                    nevapr1,delt,xhnm1,qin,ps2)
    implicit none
    integer, intent(in) :: j               ! longitude index
    ! surf geopotential (m2/s2)
    real(rkx), dimension(ici1:ici2) , intent(in) :: phis
    ! midpoint geopot convert (m) to (km)
    real(rkx), dimension(ici1:ici2,kz) , intent(in) :: zmid1
    ! temperature (K)
    real(rkx), dimension(ici1:ici2,kz) , intent(in) :: tfld
    ! dq/dt for convection (kg/kg/s) converted to (1/s) below
    ! real(rkx), dimension(ici1:ici2,kz) , intent(in) :: cmfdqr1
    ! rainwater formation tendency kg/kg/s1 converted to (1/s) :
    !          PASSED by cremrat already
    !real(rkx), dimension(ici1:ici2,kz) , intent(in) :: nrain1
    ! precip rate mm/s (=kg/m2/s), stratif and convec / grid level already,
    ! converted next to 1/s
    real(rkx), dimension(ici1:ici2,kz) , intent(in) :: strappt, convppt
    ! evaporation (kg/kg/s) convert to (1/s) below
    real(rkx), dimension(ici1:ici2,kz) , intent(in) :: nevapr1
    ! time step ( s )
    real(rkx), intent(in) :: delt
    ! total atms density (kg/m3) convert to (#/cm^3) below
    real(rkx), dimension(ici1:ici2,kz) , intent(in) :: xhnm1
    ! exported species ( mmr )
    real(rkx), dimension(ici1:ici2,kz,ntr) , intent(in) :: qin
    ! Pressure ( cb )
    real(rkx), dimension(ici1:ici2) , intent(in) :: ps2
    real(rkx) :: zmid(ici1:ici2,kz)     ! midpoint geopot convert (m) to (km)
    real(rkx) :: cmfdqr(ici1:ici2,kz)   ! dq/dt for convection (1/s)
    real(rkx) :: nrain(ici1:ici2,kz)    ! stratoform precip (1/s)
    real(rkx) :: nevapr(ici1:ici2,kz)   ! evaporation (1/s)
    real(rkx) :: xhnm(ici1:ici2,kz)     ! total atms density (#/cm^3)
    real(rkx) :: temp_dep(ici1:ici2)    ! temp var of wet dep rate
                                        ! (centibars_tracer/sec)
    ! temp var of wet dep rate
    real(rkx) :: temp_rain(ici1:ici2) , temp_wash(ici1:ici2)
                                        ! (centibars_tracer/sec)
    real(rkx) :: vmr_hno3(ici1:ici2,kz) ! volume mixing ratio
    real(rkx) :: vmr_h2o2(ici1:ici2,kz) ! volume mixing ratio
    ! Effective henry's law constant (1/s)
    real(rkx) :: het_rates(ici1:ici2,kz,ntr)

    ! mean diameter of rain drop (cm)
    real(rkx) , parameter :: xrm = 0.189_rkx
    ! mean rain drop terminal velocity (cm/s)
    real(rkx) , parameter :: xum = 748._rkx
    ! kinetic viscosity (cm^2/s)
    real(rkx) , parameter :: xvv = 6.18e-2_rkx
    ! mass transport coefficient (cm/s)
    real(rkx) , parameter :: xdg = 0.112_rkx
    ! reference temperature (k)
    real(rkx) , parameter :: t0 = 298.0_rkx
    real(rkx) , parameter :: xph0 = 1.e-5_rkx ! cloud [h+]
    ! saturation factor for hno3 in clouds
    real(rkx) , parameter :: satf_hno3 = 0.016_rkx
    ! saturation factor for hno3 in clouds
    real(rkx) , parameter :: satf_h2o2 = 0.016_rkx
    ! saturation factor for hno3 in clouds
    real(rkx) , parameter :: satf_ch2o = 0.1_rkx
    ! (atmospheres/deg k/cm^3)
    real(rkx) , parameter :: const0 = boltzk * 1.e-6_rkx
    ! hno3 dissociation constant
    real(rkx) , parameter :: hno3_diss = 15.4_rkx
    ! geometry factor (surf area/volume = geo_fac/diameter)
    real(rkx) , parameter :: geo_fac = 6.0_rkx
    ! mass of background atmosphere (amu)
    real(rkx) , parameter :: mass_air = 29.0_rkx
    ! mass of water vapor (amu)
    real(rkx) , parameter :: mass_h2o = 18.0_rkx
    real(rkx) , parameter :: h2o_mol = 1.e3_rkx/mass_h2o   ! (gm/mol water)
    real(rkx) , parameter :: km2cm = 1.e5_rkx              ! convert km to cm
    real(rkx) , parameter :: m2km = 1.e-3_rkx              ! convert m to km
    real(rkx) , parameter :: cm3_2_m3 = 1.e-6_rkx          ! convert cm^3 to m^3
    real(rkx) , parameter :: m3_2_cm3 = 1.e6_rkx           ! convert m^3 to cm^3
    real(rkx) , parameter :: liter_per_gram = 1.e-3_rkx
    ! (liter/gm/mol*(m/cm)^3)
    real(rkx), parameter :: avo2 = navgdr * liter_per_gram * cm3_2_m3

    integer(ik4) :: ktop   ! index of top model layer that can have clouds
    integer(ik4) :: i , k , kk , itr ! indicies
    real(rkx) :: xkgm         ! mass flux on rain drop
    real(rkx) :: all1 , all2  ! work variables
    real(rkx) :: stay         ! fraction of layer traversed by falling drop
                             ! in timestep delt
    real(rkx) :: xeqca1 , xeqca2 , xca1 , xca2 , xdtm
    real(rkx) :: xxx1 , xxx2 , yhno3 , yh2o2
    real(rkx) , dimension(ici1:ici2) :: xk0 , work1 , work2 , work3 , zsurf
    real(rkx) , dimension(kz) :: xgas1, xgas2
    real(rkx) , dimension(ici1:ici2) :: tmp0_rates , tmp1_rates
    ! layer depth about interfaces (cm)
    real(rkx) , dimension(ici1:ici2,kz) :: delz
    ! hno3 concentration (molecules/cm^3)
    real(rkx) , dimension(ici1:ici2,kz) :: xhno3
    ! h2o2 concentration (molecules/cm^3)
    real(rkx) , dimension(ici1:ici2,kz) :: xh2o2
    ! liquid rain water content in a grid cell (gm/m^3)
    real(rkx) , dimension(ici1:ici2,kz) :: xliq
    ! conversion rate of water vapor into rain water (molecules/cm^3/s)
    real(rkx) , dimension(ici1:ici2,kz) :: rain
    real(rkx) , dimension(ici1:ici2,kz) :: xhen_hno3 , xhen_h2o2 , xhen_ch2o , &
        xhen_ch3ooh , xhen_ch3co3h , xhen_ch3cocho , xhen_xooh , xhen_onitr ,  &
        xhen_ho2no2 , xhen_glyald , xhen_ch3cho , xhen_mvk , xhen_macr
    real(rkx) , dimension(ici1:ici2,kz) :: xhen_nh3 , xhen_ch3cooh
    real(rkx) , dimension(ici1:ici2,kz,2) :: tmp_hetrates
    real(rkx) , dimension(ici1:ici2,kz) :: precip

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

    ! nrain(ici1:ici2,:)  = nrain1(ici1:ici2,:) * (xhnm1(ici1:ici2,:)/rhoh2o)
    !
    !FOR NOW LARGE SCALE ONLY
    ! already in s-1 / only large scale precip tendency
    nrain(ici1:ici2,:) = cremrat(j,ici1:ici2,:)

    cmfdqr(:,:) = d_zero
    nevapr(:,:) = d_zero
    !FAB:  improve that by calculating or passing a grid level conversion
    ! rate for convective processes
    ![kg_h2o/kg_air/s -> 1/s]
    !cmfdqr(ici1:ici2,:) = cmfdqr1(ici1:ici2,:)*(xhnm1(ici1:ici2,:)/rhoh2o)
    ![kg_h2o/kg_air/s -> 1/s]
    ! nevapr(ici1:ici2,:) = nevapr1(ici1:ici2,:)*(xhnm1(ici1:ici2,:)/rhoh2o)

    ![m     -> km]
    zmid(ici1:ici2,:) = zmid1(ici1:ici2,:)*m2km
    ![kg/m3 -> #/cm3]
    xhnm(ici1:ici2,:) = xhnm1(ici1:ici2,:)*navgdr*cm3_2_m3/( amd*d_1000)


    !-----------------------------------------------------------------
    !	... the 2 and .6 multipliers are from a formula by frossling (1938)
    !-----------------------------------------------------------------
    xkgm = xdg/xrm * d_two + xdg/xrm * 0.6_rkx * &
           sqrt( xrm*xum/xvv ) * (xvv/xdg)**(d_one/d_three)

    do k = ktop + 1 , kz

      precip(:,k) = cmfdqr(:,k) + nrain(:,k) - nevapr(:,k)

      rain(:,k)   = mass_air*precip(:,k)*xhnm(:,k) / mass_h2o

      ! Note : xliq calculated from cloud to rain convesrion rate, shouldn't it
      ! be also from precip rate  ?
      ! In fact xliq is the rainwater generated at each level but is also be
      ! used to scavenge gas in level below in washout process after iteration
      ! on all cloudy levels, it is equivalent to have an approach with
      ! precipitation rate since: a given layer is sucessively washd out
      ! by rainwater converted in the different cloudy layer abvove.

      ! follows the original mozart code here
      ! xliq is ig g (rainwater fromed during dt) / m3
      xliq(:,k)   = precip(:,k) * delt * xhnm(:,k) / navgdr*mass_air * m3_2_cm3

      !---------------------------------------
      ! convert from cb_X to mass mixing ratio
      !---------------------------------------
      if ( ihno3 > 0 ) vmr_hno3(:,k) = (qin(:,k,ihno3)/ps2(:)) * (amd/63.012_rkx)
      if ( ih2o2 > 0 ) vmr_h2o2(:,k) = (qin(:,k,ih2o2)/ps2(:)) * (amd/34.0147_rkx)
    end do

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

    zsurf(:) = m2km * phis(:) * regrav
    do k = ktop + 1 , kz - 1
      delz(:,k) = abs( (zmid(:,k) - zmid(:,k+1))*km2cm )
    end do
    delz(:,kz) = abs( (zmid(:,kz) - zsurf(:) )*km2cm )

    !-----------------------------------------------------------------
    ! ... part 0b,  for temperature dependent of henrys
    ! xxhe1 = henry con for hno3
    ! xxhe2 = henry con for h2o2
    ! lwh 10/00 -- take henry''s law constants from brasseur et al. [1999],
    ! appendix j. for hno3, also consider dissociation to
    ! get effective henry''s law constant; equilibrium
    ! constant for dissociation from brasseur et al. [1999],
    ! appendix k. assume ph=5 (set as xph0 above).
    ! heff = h*k/[h+] for hno3 (complete dissociation)
    ! heff = h for h2o2 (no dissociation)
    ! heff = h * (1 + k/[h+]) (in general)
    !-----------------------------------------------------------------
    do k = ktop + 1 , kz
      work1(:) = (t0 - tfld(:,k))/(t0*tfld(:,k))
      !-----------------------------------------------------------------
      ! 	... effective henry''s law constants:
      !	hno3, h2o2, ch2o, ch3ooh, ch3coooh (brasseur et al., 1999)
      ! xooh, onitr, macrooh (j.-f. muller; brocheton, 1999)
      ! isopooh (equal to hno3, as for macrooh)
      ! ho2no2 (mozart-1)
      ! ch3cocho, hoch2cho (betterton and hoffman, environ. sci. technol., 1988)
      ! ch3cho (staudinger and roberts, crit. rev. sci. technol., 1996)
      ! mvk, macr (allen et al., environ. toxicol. chem., 1998)
      !-----------------------------------------------------------------
      xk0(:)             = 2.1e5_rkx * exp(8700.0_rkx*work1(:))
      xhen_hno3(:,k)     = xk0(:) * (d_one + hno3_diss / xph0)
      xhen_h2o2(:,k)     = 7.45e4_rkx * exp(6620.0_rkx*work1(:))
      xhen_ch2o(:,k)     = 6.3e3_rkx * exp(6460.0_rkx*work1(:))
      xhen_ch3ooh(:,k)   = 2.27e2_rkx * exp(5610.0_rkx*work1(:))
      xhen_ch3co3h(:,k)  = 4.73e2_rkx * exp(6170.0_rkx*work1(:))
      xhen_ch3cocho(:,k) = 3.70e3_rkx * exp(7275.0_rkx*work1(:))
      xhen_xooh(:,k)     = 90.5 * exp(5607.0_rkx*work1(:))
      xhen_onitr(:,k)    = 7.51e3_rkx * exp(6485.0_rkx*work1(:))
      xhen_ho2no2(:,k)   = 2.0e4_rkx
      xhen_glyald(:,k)   = 4.1e4_rkx * exp(4600.0_rkx*work1(:))
      xhen_ch3cho(:,k)   = 1.4e1_rkx * exp(5600.0_rkx*work1(:))
      xhen_mvk(:,k)      = 21.0_rkx * exp(7800.0_rkx*work1(:))
      xhen_macr(:,k)     = 4.3_rkx * exp(5300.0_rkx*work1(:))
      xhen_ch3cooh(:,k)  = 4.1e3_rkx * exp(6300.0_rkx*work1(:))
      xhen_nh3 (:,k)     = 1.e6_rkx
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
        if ( abs(rain(i,kk)) > d_zero ) then  ! finding rain cloud
          all1 = d_zero ! accumulation to justisfy saturation
          all2 = d_zero
          stay = ((zmid(i,kk) - zsurf(i))*km2cm)/(xum*delt)
          stay = min(stay,d_one)
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
              xgas1(k) = max(xgas1(k)-xca1,d_zero)
            end if
            if ( all2 < xeqca2 ) then
              xgas2(k) = max(xgas2(k)-xca2,d_zero)
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
        if ( abs(xxx1) > d_zero ) then  ! if no washout lifetime = 1.e29
          yhno3  = xhno3(i,kk)/xxx1 * xdtm
        else
          yhno3  = 1.e29_rkx
        end if
        if ( abs(xxx2) > d_zero ) then  ! if no washout lifetime = 1.e29
          yh2o2  = xh2o2(i,kk)/xxx2 * xdtm
        else
          yh2o2  = 1.e29_rkx
        end if
        tmp_hetrates(i,kk,1) = max(d_one/yh2o2,d_zero) * stay
        tmp_hetrates(i,kk,2) = max(d_one/yhno3,d_zero) * stay
      end do level_loop1
    end do long_loop

    !-----------------------------------------------------------------
    ! part 2, in-cloud solve for low henry constant hno3 and h2o2 have
    ! both in and under cloud.
    ! FAB : tmp_hetrates is so the washout of h2O2 and hno3: used for diagnostic
    !
    !   *** abt go directly to calculating the tracer tendency (chiten)
    !   *** MOZART code originally used het_rates and passed it out
    ! old way MOZART: het_rates(:,k,ihcho)  = work3(:)
    ! for RegCM4.1  : chiten(:,k,j,ihcho)  = chiten + work3(:) * chib
    !-----------------------------------------------------------------
    level_loop2 : &
    do k = ktop + 1 , kz
      work1(:) = avo2 * xliq(:,k)
      work2(:) = const0 * tfld(:,k)
      work3(:) = max(rain(:,k)/(h2o_mol*(work1(:) + &
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
      work3(:) = max(rain(:,k)/(h2o_mol*(work1(:) + &
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
      if ( ich3oh > 0 ) then
        het_rates(:,k,ich3oh) = work3(:)
      end if
      if ( ieoh > 0 ) then
        het_rates(:,k,ieoh) = work3(:)
      end if

      if ( ich3coooh > 0 ) then
        het_rates(:,k,ich3coooh) = max(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_ch3co3h(:,k)*work2(:)))),d_zero)
      end if
      if ( ihno4 > 0 ) then
        het_rates(:,k,ihno4) = max(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_ho2no2(:,k)*work2(:)))),d_zero)
      end if
      if ( ich3cocho > 0 ) then
        het_rates(:,k,ich3cocho) = max(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_ch3cocho(:,k)*work2(:)))),d_zero)
      end if
      if ( ixooh > 0 ) then
        het_rates(:,k,ixooh) = max(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_xooh(:,k)*work2(:)))),d_zero)
      end if
      if ( ionitr > 0 ) then
        het_rates(:,k,ionitr) = max(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_onitr(:,k)*work2(:)))),d_zero)
      end if
      if ( iglyald > 0 ) then
        het_rates(:,k,iglyald) = max(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_glyald(:,k)*work2(:)))),d_zero)
      end if
      if ( iald2 > 0 ) then
        het_rates(:,k,iald2) = max(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_ch3cho(:,k)*work2(:)))),d_zero)
      end if
      if ( imvk > 0 ) then
        het_rates(:,k,imvk) = max(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_mvk(:,k)*work2(:)))),d_zero)
      end if
      if ( imacr > 0 ) then
        het_rates(:,k,imacr) = max(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_macr(:,k)*work2(:)))),d_zero)
      end if
      if ( ih2o2 > 0 ) then
        work3(:) = satf_h2o2 * max(rain(:,k)/(h2o_mol*(work1(:) + &
                   d_one/(xhen_h2o2(:,k)*work2(:)))),d_zero)
        het_rates(:,k,ih2o2) =  work3(:) + tmp_hetrates(:,k,1)
      end if

      work3(:) = tmp_hetrates(:,k,2) + satf_hno3 * &
           max(rain(:,k)/(h2o_mol*(work1(:) +    &
           d_one/(xhen_hno3(:,k)*work2(:)))),d_zero)
      tmp0_rates(:) = work3(:)
      tmp1_rates(:) = 0.2_rkx*work3(:)
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
!      if ( iso4 > 0 ) then
!        het_rates(:,k,iso4) = tmp1_rates(:)
!      end if
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
         het_rates(:,k,inh3) = max(rain(:,k)/(h2o_mol*(work1(:) + &
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
        het_rates(:,k,ich3cooh) = max(rain(:,k)/(h2o_mol*(work1(:) + &
                  d_one/(xhen_ch3cooh(:,k)*work2(:)))),d_zero)
      end if

    end do level_loop2

    !-----------------------------------------------------------------
    ! *** abt added to work with RegCM4.1
    ! Directly calculate the tracer tendency (chiten)
    ! Loop over all non-aerosol tracers
    ! old way MOZART: het_rates(:,k,ihcho)  = work3(:)
    ! for RegCM4.1  : chiten(:,k,j,ihcho)  = chiten + work3(:) * chib
    ! FAB: add diagnostic
    !-----------------------------------------------------------------

    do k= 1 , kz
      do itr = 1 , ntr
        do i = ici1 , ici2
          if ( het_rates(i,k,itr) <= 2.e-29_rkx ) cycle

          ! tendency due to rainout and washout
          temp_dep(i) = (d_one-exp(-het_rates(i,k,itr)*delt))*qin(i,k,itr)
          chiten(j,i,k,itr) = chiten(j,i,k,itr) - temp_dep(i)/delt

          ! Disagnostic
          ! rainout: rainout
          ! washout washout

          if ( itr == ih2o2 ) then
            temp_rain(i) = (d_one- &
              exp(-(het_rates(i,k,itr)-tmp_hetrates(i,k,1))*delt))*qin(i,k,itr)
            temp_wash(i) = (d_one- &
              exp(-tmp_hetrates(i,k,1)*delt))*qin(i,k,itr)
          else if ( itr == ihno3 ) then
            temp_rain(i) = (d_one- &
              exp(-(het_rates(i,k,itr)-tmp_hetrates(i,k,2))*delt))*qin(i,k,itr)
            temp_wash(i) = (d_one- &
              exp(-tmp_hetrates(i,k,2)*delt))*qin(i,k,itr)
          else
            temp_rain(i) = temp_dep(i)
            temp_wash(i) = d_zero
          end if
          rainout(j,i,k,itr) = rainout(j,i,k,itr) - temp_rain(i)*cfdout
          washout(j,i,k,itr) = washout(j,i,k,itr) - temp_wash(i)*cfdout
        end do
      end do
    end do

    ! diagnostic for durface fluxes
    do itr = 1 , ntr
      do i = ici1 , ici2
        wdrout(j,i,itr) = d_zero
        wdwout(j,i,itr) = d_zero
        do k = 1 , kz
          ! sum on the vertical to get total surface flux diag fo rain out
          ! and washout (already weighted for time average cfdout !),
          ! also change sign convention normalise by psb to get the right
          ! flux unit
          wdrout(j,i,itr) = wdrout(j,i,itr) - &
            rainout(j,i,k,itr)*cdzq(j,i,k) *crhob3d(j,i,k) /cpsb(j,i)
          wdwout(j,i,itr) = wdwout(j,i,itr) - &
            washout(j,i,k,itr)*cdzq(j,i,k) *crhob3d(j,i,k)/cpsb(j,i)
        end do
      end do
    end do

  end subroutine sethet

  subroutine wetdepa(j,mbin,indp,beffdiam,rhoaer,t,wl,fracloud,fracum, &
                     pressg,shj,rho,strappt,convppt,pdepv)
    implicit none
    integer(ik4) , intent(in) :: j , mbin
    integer(ik4) , dimension(mbin) , intent(in) :: indp
    real(rkx) , dimension(mbin) , intent(in) :: beffdiam
    real(rkx) , intent(in) :: rhoaer ! specific aerosol density
    real(rkx) , dimension(ici1:ici2,kz) , intent(in) :: wl , t , rho , &
                                                       strappt
     real(rkx) , dimension(ici1:ici2,kz) , intent(in)  :: convppt
    real(rkx) , dimension(ici1:ici2,kz) , intent(in) :: fracloud , fracum
    real(rkx) , dimension(ici1:ici2) ,intent(in) :: pressg
    real(rkx) , dimension(kz) ,intent(in) :: shj
    real(rkx) , dimension(ici1:ici2,kz,ntr) , intent (in) :: pdepv
    ! size of the aerosol bin
    ! index of the correponding aerosol in the chi table
    real(rkx) , dimension(ici1:ici2,kz) :: totppt
    real(rkx) , dimension(ici1:ici2,kz,mbin) :: colef , wetdep , rhsize , rhop
    real(rkx) , dimension(ntr) :: wetrem , wetrem_cvc
    real(rkx) :: wtend
    integer(ik4) :: n , k , i, nk,nkh

    ! rain out parametrisation
    ! clmin = non-precipitating cloud
    ! conversion threshold, clmin=0.01g/m3
    real(rkx) , parameter :: clmin = 0.01_rkx
    ! remcum= removal rate for cumulus
    ! cloud scavenging (s-1)
    real(rkx) , parameter :: remcum = 1.0e-3_rkx

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

      ! Notes :
      ! For the CLM interface,
      ! We calculate now  an instantaneous surface wet depoistion fluxcwet dep
      ! representing the sum of convective and large scale washout and rainout.
      !
      ! For the diag, convective and large scale rainout ( washout) are
      ! added

      ! Important: convective, large scale, washout and rainout fluxes
      ! of every vertical level are summed for surface flux calculation.
      ! cwetdepflx(:,:,indp(n)) is also a time accumulated array which is
      ! set to 0 when the surface (CLM45) scheme is called. The average between
      ! surface time steps is calculated in the atm to surface interface
      if ( ichremlsc == 1 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            if ( wl(i,k) > clmin ) then
              wetrem(indp(n)) = d_zero
              if ( cremrat(j,i,k) > d_zero .and. fracloud(i,k) > d_zero ) then
                wetrem(indp(n)) = fracloud(i,k)*chtrsol(indp(n)) * &
                   chib(j,i,k,indp(n))* &
                   (exp(-cremrat(j,i,k)/fracloud(i,k)*dt)-d_one)
                chiten(j,i,k,indp(n)) = chiten(j,i,k,indp(n)) + &
                   wetrem(indp(n))/dt

                ! sum up the flux on the vertical to get instantaneous surface
                ! deposition flux
                ! in (Kg/m2/s) passed to CLM surface scheme ( change sign
                ! convention to get a positive deposition flux)
                cwetdepflx(j,i,indp(n))  = cwetdepflx(j,i,indp(n)) - &
                   wetrem(indp(n))/ dt *cdzq(j,i,k)*crhob3d(j,i,k) /cpsb(j,i)
                ! save the tendency as a diag (accumulated for average at output
                ! time step)
                rainout(j,i,k,indp(n)) = rainout(j,i,k,indp(n)) + &
                   wetrem(indp(n))/dt  *cfdout
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
                   chib(j,i,k,indp(n))*(exp(-remcum*dt)-d_one)
              chiten(j,i,k,indp(n)) = chiten(j,i,k,indp(n)) + &
                   wetrem_cvc(indp(n))/dt
              ! add to large scale and sum up the flux on the vertical to
              ! get instantaneous surface flux in (Kg/m2/s) passed to CLM
              ! surface scheme ( change sign convention to get a positive flux)
              cwetdepflx(j,i,indp(n))  = cwetdepflx(j,i,indp(n)) - &
                wetrem_cvc(indp(n))/ dt *cdzq(j,i,k)*crhob3d(j,i,k) /cpsb(j,i)
              ! add the concvetive rainout to large scale save the
              ! tendency as a diag
              rainout(j,i,k,indp(n)) = rainout(j,i,k,indp(n)) + &
                   wetrem_cvc(indp(n))/ dt *cfdout
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
        rhsize(i,k,1:mbin) = d_half * 1.e-6_rkx * beffdiam(1:mbin)
      end do
    end do
    ! dry density for now
    rhop(:,:,:) = rhoaer

    ! convppt is a pseudo-3d array of con prec rate conatining constant rate
    ! (equals to surface cumul precipi) from surface to cumtop.
    ! inside the cloud, need itnterpolate linearly between the prec rate and
    ! o ( top of the cloud). Consider that half of the convective
    ! column is cloud .. improve this in the future !
    totppt(:,:) = d_zero
    do i =  ici1 , ici2
      if ( kcumtop(j,i) > 0 ) then
        nk = kcumbot(j,i) - kcumtop(j,i) + 1
        n = 1
        do k =  kcumtop(j,i) , kz
          nkh = nk / 2
          if ( n <= nkh ) then
            totppt(i,k) = (real(n,rkx) / real(nkh,rkx)) * convppt(i,kz)
          else
            totppt(i,k)  = convppt(i,kz)
          end if
          n = n + 1
        end do
      end if
    end do
    ! add the contribution of strat prec
    totppt(:,:) = strappt(:,:) + totppt(:,:)
    call blcld(mbin,indp,rhsize,t,pressg,shj,rho,totppt,pdepv,rhop,wetdep,colef)

    ! calculate the tendency due to wahsout dep here.

    do n = 1 , mbin
      do k = 1 , kz
        do i = ici1 , ici2
          wtend = chib(j,i,k,indp(n))*(d_one-exp(-wetdep(i,k,n)*dt))/dt
          chiten(j,i,k,indp(n)) = chiten(j,i,k,indp(n)) - wtend

          ! add to rainout flux  and sum up on the vertical to get
          ! instantaneous surface flux (in Kg/m2/s) passed to CLM
          ! surface scheme.
          cwetdepflx(j,i,indp(n))  = cwetdepflx(j,i,indp(n)) + &
               wtend * cdzq(j,i,k) * crhob3d(j,i,k) / cpsb(j,i)
          ! wet deposition washout diagnostic ( both from conv and large scale)
          washout(j,i,k,indp(n)) = washout(j,i,k,indp(n)) - wtend * cfdout
        end do
      end do
    end do

    do n = 1, mbin
      do i = ici1 , ici2
        wdrout(j,i,indp(n)) = d_zero
        wdwout(j,i,indp(n)) = d_zero
        do k = 1 , kz
          ! sum on the vertical to get total surface flux diag for rain out
          ! and washout (already weighted for time average cfdout !),
          ! also change sign convention
          ! normalise by psb to get the right flux unit
          wdrout(j,i,indp(n)) = wdrout(j,i,indp(n)) - &
            rainout(j,i,k,indp(n))*cdzq(j,i,k) *crhob3d(j,i,k) /cpsb(j,i)
          wdwout(j,i,indp(n)) = wdwout(j,i,indp(n)) - &
            washout(j,i,k,indp(n))*cdzq(j,i,k) *crhob3d(j,i,k)/cpsb(j,i)
        end do
      end do
    end do

  end subroutine wetdepa

  subroutine blcld(mbin,indp,rhsize,t,pressg,shj,rho,totppt,pdepv, &
                   rhop,wetdep,colef)
    implicit none

    integer(ik4) , intent(in) :: mbin
    real(rkx) , dimension(ici1:ici2,kz,mbin) , intent(in) :: rhsize , rhop
    real(rkx) , dimension(ici1:ici2,kz) , intent(in) :: t , rho , totppt
    ! care ,ntr dimension
    real(rkx) , dimension(ici1:ici2,kz,ntr) , intent(in) :: pdepv
    real(rkx) , dimension(ici1:ici2) , intent(in) :: pressg
    real(rkx) , dimension(kz) , intent(in) :: shj
    real(rkx) , dimension(ici1:ici2,kz,mbin) , intent(out) :: colef , wetdep

    ! index of the correponding aerosol in the chi table
    integer(ik4) , dimension(mbin) , intent(in) :: indp

    real(rkx) :: dm , tl , rrm
    integer(ik4) :: i , n , k

    real(rkx) , parameter :: bcrain = 0.5_rkx
    real(rkx) , parameter :: bcsnow = 0.8_rkx

    !----------------------------------------------------------------------c
    ! call to compute collection efficiency coefficients                   c
    !----------------------------------------------------------------------c

    call cas(t,colef,indp,rhop,rho,rhsize,pressg,mbin,totppt,pdepv,shj)

    wetdep(:,:,:) = d_zero
    do n = 1 , mbin
      do k = 1 , kz
        do i = ici1 , ici2
          if ( totppt(i,k) > 1.e-20_rkx ) then
            tl = t(i,k)-tzero
            !----------------------------------------------------------c
            ! rain scavenging rate : wetdep in s-1   ?                 c
            ! units of totppt (mm.s-1) converted to m.s-1 (by 1.0e-3)  c
            !----------------------------------------------------------c
            if ( tl > d_zero ) then
              rrm = 0.35_rkx*(totppt(i,k)*3600.0_rkx)**0.25_rkx*1.0e-3_rkx
              ! formulation for monoddipesre rain of radius rrm
              ! (cf seinfeld), gives a scavenging rate wetdep in
              ! s-1 ?? bcrain is 0.5 when it is 1.5 in seinfeld ??
              wetdep(i,k,n) = bcrain*totppt(i,k)*1.0e-3_rkx*colef(i,k,n)/rrm
            end if
            !---------------------------------------------------------c
            ! snow scavenging rate                                    c
            ! * density of snow is set as 1/10 of liquid water.       c
            ! * factor of 1.0e-2 in wetdep calculation takes this and c
            ! * unit                                                  c
            !      change (to m.s-1) into account                     c
            !---------------------------------------------------------c
            ! . . . . snow scavenging
            if ( tl <= d_zero .and. tl >= -8.0_rkx ) then
              dm = 3.8e-5_rkx ! characteristic length scale [m]
              wetdep(i,k,n) = bcsnow*totppt(i,k)*1.0e-3_rkx*colef(i,k,n)/dm
            end if
            ! . . . . steller snow scavenging
            if ( tl < -8. .and. tl >= -25.0_rkx ) then
              dm = 2.7e-5_rkx
              wetdep(i,k,n) = bcsnow*totppt(i,k)*1.0e-3_rkx*colef(i,k,n)/dm
            end if
            ! . . . . graupel scavenging
            if ( tl > -25.0_rkx ) then
              dm = 1.4e-4_rkx
              wetdep(i,k,n) = bcsnow*totppt(i,k)*1.0e-3_rkx*colef(i,k,n)/dm
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
    !---------------------------------------------------------------------
    !  purpose:
    !  --------
    !  this is a module of calculating the collection efficiency of
    !  the aerosols of type n by collector droplets of radius rcol.
    !
    implicit none

    integer(ik4) , intent(in) :: mbin
    integer(ik4) , dimension(mbin) , intent(in) :: indp
    real(rkx) , dimension(ici1:ici2,kz) , intent(in) :: t , rho , totppt
    real(rkx) , dimension(ici1:ici2,kz,mbin) , intent(in) :: rhop , rhsize
    real(rkx) , dimension(ici1:ici2) , intent(in) :: pressg
    real(rkx) , dimension(ici1:ici2,kz,ntr) , intent(in) :: pdepv
    real(rkx) , dimension(kz) , intent(in) :: shj
    ! collection efficiency
    real(rkx) , dimension(ici1:ici2,kz,mbin) , intent(out) :: colef

    real(rkx) :: pdiff
    real(rkx) :: amu , anu !dynamic viscosity of air
    real(rkx) :: amfp   ! mean molecular free path
    real(rkx) :: schm   ! schmidt number
    real(rkx) :: prii
    real(rkx) :: priiv
    ! real(rkx) :: cfac
    real(rkx) :: cfaca
    real(rkx) :: re ! reynolds number
    real(rkx) :: rr ! (collected particle radius)/(collector particle radius)
    real(rkx) :: st ! stokes number of collected particles
    real(rkx) :: vr
    real(rkx) :: sstar
    real(rkx) :: amob
    real(rkx) :: colimp , pre
    real(rkx) :: vpr ! average settling velocity
    real(rkx) :: rrm ! mass mean raindrop radius (for rain)
                    ! characteristic capture length (for snow)
    real(rkx) :: alpha
    integer(ik4) :: n , k , i

    ! Cunningham slip correction factor parameters for rain

    !real(rkx) , parameter :: aa1r = 1.249_rkx
    !real(rkx) , parameter :: aa2r = 0.42_rkx
    !real(rkx) , parameter :: aa3r = 0.87_rkx
    real(rkx) , parameter :: rhorain = 1000.0_rkx
    real(rkx) , parameter :: amuw = 1.002e-3_rkx ! at 20*c [kg/m/sec]

    colef(:,:,:) = d_zero

    do n = 1 , mbin
      do k = 1 , kz
        do i = ici1 , ici2
          ! . . . for precipitation:
          vpr = d_zero
          if ( totppt(i,k) > 1.e-15_rkx ) then
            !********************************************************
            !* air's dynamic viscosity                           ****
            !********************************************************
            amu = a1*1.e-8_rkx*t(i,k)**a2/(t(i,k)+a3)
            !. . . . mid layer pressure in [pascal].
            pre = pressg(i)*shj(k)
            !********************************************************
            !* mean molecular free path.                         ****
            !*     k.v. beard [1976], j atm. sci., 33            ****
            !********************************************************
            amfp = c1*(amu/c2)*(c3/pre)*sqrt(t(i,k)/c4)
            !********************************************************
            !* cunningham slip correction factor and             ****
            !* relaxation time = vg/grav.                        ****
            !********************************************************
            cfaca = d_one + amfp/rhsize(i,k,n) * &
                      (aa1+aa2*exp(-aa3*rhsize(i,k,n)/amfp))
            !*****************************************************
            !* the schmidt number is the ratio of the         ****
            !* kinematic viscosity of air to the particle     ****
            !* brownian diffusivity ===> sc=v/d               ****
            !*****************************************************
            anu = amu/rho(i,k)
            amob = 6.0_rkx*mathpi*amu*rhsize(i,k,n)/cfaca
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
              rrm = 0.35_rkx*(totppt(i,k)*3600.0_rkx)**0.25_rkx*1.e-3_rkx
              ! what empirical parametization is that ?
              ! needs apparently rain rate in mm/hr
              prii = 2.0_rkx/9.0_rkx*egrav/amu
              priiv = prii*(rhorain-rho(i,k))
              ! cunningham settling slip-correction factor
              ! settling velocity
              ! FAB : wrong this formulation apply only for small
              ! particles settling but not for big rain drops (high reynolds)!!
              ! cf Seinfeld / Vpr would be overestimated
              !  cfac = d_one+amfp/rrm*(aa1r+aa2r*exp(-aa3r*rrm/amfp))
              ! vpr = priiv*rrm**2*cfac
              ! try with a mean rainfall velcoity of 3 m/s
              vpr = 3._rkx
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
            if ( t(i,k) <= tzero .and. t(i,k)>= tzero-8.0_rkx ) then
              vpr = 50.0e-2_rkx
              rrm = 10.e-6_rkx
              alpha = 1.0_rkx
            end if
            ! . . . steller-snow scavenging:
            if ( t(i,k) < tzero-8.0_rkx .and. t(i,k) >= tzero-25.0_rkx) then
              vpr = 57.0e-2_rkx
              rrm = 100.e-6_rkx
              alpha = 0.5_rkx
            end if
            ! . . . graupel scavenging:
            if ( t(i,k) < tzero-25.0_rkx ) then
              vpr = 180.0e-2_rkx
              rrm = 1000.e-6_rkx
              alpha = 2.0_rkx/3.0_rkx
            endif

            ! reynolds number:
            re = rrm*vpr*rho(i,k)/amu
            ! stokes number of collected particles:
            st = d_two*pdepv(i,k,indp(n))*regrav*(vpr - &
                 pdepv(i,k,indp(n)))/(d_two*rrm)
            ! ratio of radius of collected particle to colector drop:
            rr = rhsize(i,k,n)/rrm
            vr = amuw/amu
            sstar = (1.2_rkx+(d_one/12.0_rkx)* &
                    log(d_one+re))/(d_one+log(d_one+re))
            if ( st > sstar ) then
              colimp = ((st-sstar)/ &
                        (st-sstar+d_two/d_three))**(d_three/d_two) * &
                         sqrt(d_1000/rhop(i,k,n))
            else
              colimp = d_zero
            end if
            ! rain scavenging efficiency:
            if ( t(i,k) > tzero ) then
              colef(i,k,n) = d_four/(re*schm)*(d_one+0.4_rkx*sqrt(re) * &
                   schm**(d_one/d_three)+0.16_rkx*sqrt(re*schm)) +  &
                   d_four*rr*(d_one/vr + (d_one+d_two*sqrt(re))*rr) + colimp
            else
            ! snow scavenging efficiency:
            colef(i,k,n) = (d_one/schm)**alpha +  &
                   (d_one-exp(-(d_one+sqrt(re))*rr**2)) + colimp
            endif
            ! setting the upper-bound for collection efficiency:
            colef(i,k,n) = max(d_zero,min(d_one,colef(i,k,n)))
          end if
        end do
      end do
    end do
  end subroutine cas

end module mod_che_wetdep
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
