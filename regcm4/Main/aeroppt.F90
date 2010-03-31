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
 
      subroutine aeroppt(rh,pint,tauxar_mix,tauasc_mix,gtota_mix,       &
                       & ftota_mix,tauxar_mix_cs,tauasc_mix_cs,         &
                       & gtota_mix_cs,ftota_mix_cs)
 
      use mod_regcm_param
      use mod_param2
      use mod_aerosol , only : nspi , ksdust , wsdust , gsdust ,        &
                       & kscoef , wscoef , gscoef , ksoc_hl , wsoc_hl , &
                       & gsoc_hl , ksbc_hl , wsbc_hl , gsbc_hl ,        &
                       & ksoc_hb , wsoc_hb, gsoc_hb, ksbc_hb , wsbc_hb ,&
                       & gsbc_hb , dextmix , dssamix , dgmix , ksbase , &
                       & wsbase , gsbase , aermmr
      use mod_trachem
      use mod_message
      use mod_constants , only : rhoso4 , rhobc, rhooc, rhodust , gtigts
      implicit none
!
! Dummy arguments
!
!     Aerosol extinction optical depth
      real(8) , dimension(iym1,0:kz,nspi) :: ftota_mix , gtota_mix ,    &
           & tauasc_mix , tauxar_mix
      real(8) , dimension(iym1,nspi) :: ftota_mix_cs , gtota_mix_cs ,   &
           & tauasc_mix_cs , tauxar_mix_cs
!     Interface pressure, relative humidity
      real(8) , dimension(iym1,kzp1) :: pint
      real(8) , dimension(iym1,kz) :: rh
      intent (in) pint , rh
      intent (inout) ftota_mix , ftota_mix_cs , gtota_mix ,             &
                   & gtota_mix_cs , tauasc_mix , tauasc_mix_cs ,        &
                   & tauxar_mix , tauxar_mix_cs
!
! Local variables
!
      real(8) , dimension(iym1,kz) :: aermtot , aervtot
      real(8) , dimension(iym1,0:kz,ntr) :: fa , ga , tauxar , uaer ,   &
           & wa
      real(8) , dimension(iym1,ntr) :: faer , gaer , tauaer , utaer ,   &
                                      & waer
      real(8) , dimension(4) :: frac , prop
      integer :: i , i1 , i2 , i3 , i4 , ibin , itr , k , ns
      real(8) :: path
!
! uaer, tauxar  - Aerosol radiative properties (local arrays)
! wa            - Aerosol single scattering albedo
! ga            - Aerosol asimmetry parameter
! fa            - Aerosol forward scattered fraction
! utaer, tauaer - Total column aerosol extinction
! waer          - Aerosol single scattering albedo
! gaer          - Aerosol asymmetry parameter
! faer          - Aerosol forward scattered fraction
!
! Visible band
!     real(kind=8)  aertau(iym1,ntr)
!     real(kind=8)  aerprf(iym1,0:kz,ntr)
!     real(kind=8)  tauprf

!trapuv
      tauxar = 0.0
      wa = 0.0
      ga = 0.0
      fa = 0.0
!trapuv_
 
!     Spectral loop
      do ns = 1 , nspi
!---------------
        do itr = 1 , ntr
          do i = 1 , iym1
!           set above top value
            uaer(i,0,itr) = 0.0
            tauxar(i,0,itr) = 0.0
            wa(i,0,itr) = 0.0
            ga(i,0,itr) = 0.0
            fa(i,0,itr) = 0.0
            utaer(i,itr) = 0.0
            tauaer(i,itr) = 0.0
            waer(i,itr) = 0.0
            gaer(i,itr) = 0.0
            faer(i,itr) = 0.0
          end do
        end do
 
        do i = 1 , iym1
          tauxar_mix_cs(i,ns) = 0.0
          tauasc_mix_cs(i,ns) = 0.0
          gtota_mix_cs(i,ns) = 0.0
          ftota_mix_cs(i,ns) = 0.0
        end do
 
        do k = 0 , kz
          do i = 1 , iym1
            tauxar_mix(i,k,ns) = 0.0
            tauasc_mix(i,k,ns) = 0.0
            gtota_mix(i,k,ns) = 0.0
            ftota_mix(i,k,ns) = 0.0
          end do
        end do
 
        if ( idirect.ge.1 ) then
!
!         calculate optical properties of each aerosol component
!
 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!         Option 1  melange externe
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          if ( mixtype.eq.1 ) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            do k = 1 , kz
              do i = 1 , iym1
                path = (pint(i,k+1)-pint(i,k))/gtigts
                ibin = 0
                do itr = 1 , ntr
                  uaer(i,k,itr) = 0.
                  if ( rh(i,k).lt.0.0 .or. rh(i,k).gt.1.0 ) print * ,   &
                     & i , k , rh(i,k) , '  RH WARNING !!!!!'
 
                  if ( chtrname(itr).eq.'DUST' ) then
                    uaer(i,k,itr) = aermmr(i,k,itr)*path
!                   gaffe au facteur !!
!                   gaffe au ntr/bins
                    ibin = ibin + 1
                    if ( ibin.gt.4 ) print * , 'DUST OP PBLEME !!!!'
 
                    tauxar(i,k,itr) = 1.D5*uaer(i,k,itr)*ksdust(ns,ibin)
                    wa(i,k,itr) = wsdust(ns,ibin)
                    ga(i,k,itr) = gsdust(ns,ibin)
                    fa(i,k,itr) = gsdust(ns,ibin)*gsdust(ns,ibin)
 
                  else if ( chtrname(itr).eq.'SO4' ) then
                    uaer(i,k,itr) = aermmr(i,k,itr)*path
                    tauxar(i,k,itr) = 1.D5*uaer(i,k,itr)*ksbase(ns)     &
                                    & *exp(kscoef(ns,1)+kscoef(ns,2)    &
                                    & /(rh(i,k)+kscoef(ns,3))           &
                                    & +kscoef(ns,4)                     &
                                    & /(rh(i,k)+kscoef(ns,5)))
!
                    wa(i,k,itr) = 1.0 - wsbase(ns)                      &
                                & *exp(wscoef(ns,1)+wscoef(ns,2)        &
                                & /(rh(i,k)+wscoef(ns,3))+wscoef(ns,4)  &
                                & /(rh(i,k)+wscoef(ns,5)))
!
                    ga(i,k,itr) = gsbase(ns)                            &
                                & *exp(gscoef(ns,1)+gscoef(ns,2)        &
                                & /(rh(i,k)+gscoef(ns,3))+gscoef(ns,4)  &
                                & /(rh(i,k)+gscoef(ns,5)))
!
                    fa(i,k,itr) = ga(i,k,itr)*ga(i,k,itr)
 
                  else if ( chtrname(itr).eq.'OC_HL' ) then
                    uaer(i,k,itr) = aermmr(i,k,itr)*path
!                   Humidity effect !
                    tauxar(i,k,itr) = 1.D5*uaer(i,k,itr)*ksoc_hl(ns)    &
                                    & *(1-rh(i,k))**(-0.2)
                    wa(i,k,itr) = wsoc_hl(ns)
                    ga(i,k,itr) = gsoc_hl(ns)
                    fa(i,k,itr) = ga(i,k,itr)*ga(i,k,itr)
                  else if ( chtrname(itr).eq.'BC_HL' ) then
                    uaer(i,k,itr) = aermmr(i,k,itr)*path
!                   Humidity effect !
                    tauxar(i,k,itr) = 1.D5*uaer(i,k,itr)*ksbc_hl(ns)    &
                                    & *(1-rh(i,k))**(-0.25)
                    wa(i,k,itr) = wsbc_hl(ns)
                    ga(i,k,itr) = gsbc_hl(ns)
                    fa(i,k,itr) = ga(i,k,itr)*ga(i,k,itr)
                  else if ( chtrname(itr).eq.'OC_HB' ) then
                    uaer(i,k,itr) = aermmr(i,k,itr)*path
                    tauxar(i,k,itr) = 1.D5*uaer(i,k,itr)*ksoc_hb(ns)
                    wa(i,k,itr) = wsoc_hb(ns)
                    ga(i,k,itr) = gsoc_hb(ns)
                    fa(i,k,itr) = gsoc_hb(ns)*gsoc_hb(ns)
                  else if ( chtrname(itr).eq.'BC_HB' ) then
                    uaer(i,k,itr) = aermmr(i,k,itr)*path
!                   Absorbing aerosols (soot type)
                    tauxar(i,k,itr) = 1.D5*uaer(i,k,itr)*ksbc_hb(ns)
                    wa(i,k,itr) = wsbc_hb(ns)
                    ga(i,k,itr) = gsbc_hb(ns)
                    fa(i,k,itr) = gsbc_hb(ns)*gsbc_hb(ns)
                  else
                  end if
                end do  ! end tracer loop
              end do
            end do
!           optical properties for the clear sky diagnostic
            do i = 1 , iym1
              do itr = 1 , ntr
                do k = 1 , kz
                  utaer(i,itr) = utaer(i,itr) + uaer(i,k,itr)
                  tauaer(i,itr) = tauaer(i,itr) + tauxar(i,k,itr)
                  waer(i,itr) = waer(i,itr) + wa(i,k,itr)*uaer(i,k,itr)
                  gaer(i,itr) = gaer(i,itr) + ga(i,k,itr)*uaer(i,k,itr)
                  faer(i,itr) = faer(i,itr) + fa(i,k,itr)*uaer(i,k,itr)
                end do
                if ( utaer(i,itr).le.1.D-10 ) utaer(i,itr) = 1.D-10
                waer(i,itr) = waer(i,itr)/utaer(i,itr)
                gaer(i,itr) = gaer(i,itr)/utaer(i,itr)
                faer(i,itr) = faer(i,itr)/utaer(i,itr)
              end do
            end do
!
!           Calculate the EXTERNAL Mixing of aerosols
!           melange externe
!
            do i = 1 , iym1
              do itr = 1 , ntr
!               only for climatic feedback allowed
                do k = 0 , kz
                  tauxar_mix(i,k,ns) = tauxar_mix(i,k,ns)               &
                                     & + tauxar(i,k,itr)
                  tauasc_mix(i,k,ns) = tauasc_mix(i,k,ns)               &
                                     & + tauxar(i,k,itr)*wa(i,k,itr)
                  gtota_mix(i,k,ns) = gtota_mix(i,k,ns) + ga(i,k,itr)   &
                                    & *tauxar(i,k,itr)*wa(i,k,itr)
                  ftota_mix(i,k,ns) = ftota_mix(i,k,ns) + fa(i,k,itr)   &
                                    & *tauxar(i,k,itr)*wa(i,k,itr)
                end do
!
!               Clear sky (always calcuated if idirect >=1 for
!               diagnostic radiative forcing)
!
                tauxar_mix_cs(i,ns) = tauxar_mix_cs(i,ns)               &
                                    & + tauaer(i,itr)
                tauasc_mix_cs(i,ns) = tauasc_mix_cs(i,ns)               &
                                    & + tauaer(i,itr)*waer(i,itr)
                gtota_mix_cs(i,ns) = gtota_mix_cs(i,ns) + gaer(i,itr)   &
                                   & *tauaer(i,itr)*waer(i,itr)
                ftota_mix_cs(i,ns) = ftota_mix_cs(i,ns) + faer(i,itr)   &
                                   & *tauaer(i,itr)*waer(i,itr)
              end do
            end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!           Option 2  melange interne
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          else if ( mixtype.eq.2 ) then
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            do i = 1 , iym1
              do k = 1 , kz
                path = (pint(i,k+1)-pint(i,k))/gtigts
                if ( rh(i,k).lt.0.0 .or. rh(i,k).gt.1.0 ) print * , i , &
                   & k , rh(i,k) , '  RH WARNING !!!!!'
 
!               sum of hydrophilic aerosols
                aervtot(i,k) = 0.
                aermtot(i,k) = 0.
 
                if ( iso4.ne.0 ) then
                  aervtot(i,k) = aervtot(i,k) + aermmr(i,k,iso4)/rhoso4
                  aermtot(i,k) = aermtot(i,k) + aermmr(i,k,iso4)
                end if
 
                if ( ibchl.ne.0 ) then
                  aervtot(i,k) = aervtot(i,k) + aermmr(i,k,ibchl)/rhobc
                  aermtot(i,k) = aermtot(i,k) + aermmr(i,k,ibchl)
                end if
 
                if ( iochl.ne.0 ) then
                  aervtot(i,k) = aervtot(i,k) + aermmr(i,k,iochl)/rhooc
                  aermtot(i,k) = aermtot(i,k) + aermmr(i,k,iochl)
                end if
 
                if ( idust(1).ne.0 ) then
                  aervtot(i,k) = aervtot(i,k) + aermmr(i,k,idust(1))    &
                               & /rhodust
                  aermtot(i,k) = aermtot(i,k) + aermmr(i,k,idust(1))
                end if
 
!               minimum quantity of total aerosol
 
                if ( aermtot(i,k).gt.1.D-14 ) then
!                 indexes in the internal mixing table
                  prop(1) = (aermmr(i,k,iso4)/rhoso4)/aervtot(i,k)
                  prop(2) = (aermmr(i,k,ibchl)/rhobc)/aervtot(i,k)
                  prop(3) = (aermmr(i,k,iochl)/rhooc)/aervtot(i,k)
                  prop(4) = (aermmr(i,k,idust(1))/rhodust)/aervtot(i,k)
                  frac(1) = fraction(prop(1))
                  frac(2) = fraction(prop(2))
                  frac(3) = fraction(prop(3))
                  frac(4) = fraction(prop(4))
!                 FIND THE GREATEST FRACTIONAL PART
 
                  if ( iso4.ne.0 ) then
                    i1 = nint(10*prop(1)) + 1
                  else
                    i1 = 0 + 1
                  end if
 
                  if ( ibchl.ne.0 ) then
                    i2 = nint(10*prop(2)) + 1
                  else
                    i2 = 0 + 1
                  end if
 
                  if ( iochl.ne.0 ) then
                    i3 = nint(10*prop(3)) + 1
                  else
                    i3 = 0 + 1
                  end if
 
                  if ( idust(1).ne.0 ) then
                    i4 = nint(10*prop(4)) + 1
                  else
                    i4 = 0 + 1
                  end if
!
!                 final optical parameters
!
                  if ( i1+i2+i3+i4.eq.13 ) i4 = i4 + 1
                  if ( i1+i2+i3+i4.eq.15 ) then
                    if ( i4.ne.1 ) i4 = i4 - 1
                  end if
 
                  if ( i1+i2+i3+i4.eq.15 ) then
                    if ( i1.ne.1 ) i1 = i1 - 1
                  end if
 
                  if ( i1+i2+i3+i4.eq.15 ) then
                    if ( i3.ne.1 ) i3 = i3 - 1
                  end if
 
                  if ( i1+i2+i3+i4.eq.15 ) call fatal(__FILE__,__LINE__,&
                      &'WRONG COMBINATION. SHOULD NEVER HAPPEN')
 
                  if ( i1+i2+i3+i4.ne.14 ) then
                    print * , i1 , i2 , i3 , i4 , i1 + i2 + i3 + i4
                    print * , idust(1) , iochl , ibchl , iso4
                    print * , 'OC HL' , aermmr(i,k,iochl)/rhooc
                    print * , 'BC HL' , aermmr(i,k,ibchl)/rhobc
                    print * , 'SO4' , aermmr(i,k,iso4)/rhoso4
                    print * , 'DUST' , aermmr(i,k,idust(1))/rhodust
                    print * , 'VOL TOT' , aervtot(i,k)
                    print * , 'OC HL%' , 10*(aermmr(i,k,iochl)/rhooc)   &
                        & /aervtot(i,k)
                    print * , 'BC HL%' , 10*(aermmr(i,k,ibchl)/rhobc)   &
                        & /aervtot(i,k)
                    print * , 'SO4 %' , 10*(aermmr(i,k,iso4)/rhoso4)    &
                        & /aervtot(i,k)
                    print * , 'SO4 %' ,                                 &
                        & nint(10*(aermmr(i,k,iso4)/rhoso4)/aervtot(i,k)&
                        & )
                    print * , 'DUST %' ,                                &
                        & 10*(aermmr(i,k,idust(1))/rhodust)/aervtot(i,k)
                    print * , 'DUST %' ,                                &
                        & nint(10*(aermmr(i,k,idust(1))/rhodust)        &
                        & /aervtot(i,k))
                    call fatal(__FILE__,__LINE__,                       &
                              &'SOMETHING WRONG ON SPECIES ABUNDANCE')
                  end if
                  tauxar_mix(i,k,ns) = dextmix(1,ns,i4,i2,i3,i1)        &
                                     & *aermtot(i,k)*path*1D5
                  tauasc_mix(i,k,ns) = dssamix(1,ns,i4,i2,i3,i1)        &
                                     & *tauxar_mix(i,k,ns)
                  gtota_mix(i,k,ns) = dgmix(1,ns,i4,i2,i3,i1)           &
                                    & *tauasc_mix(i,k,ns)               &
                                    & *tauxar_mix(i,k,ns)
                  ftota_mix(i,k,ns) = dgmix(1,ns,i4,i2,i3,i1)           &
                                    & *dgmix(1,ns,i4,i2,i3,i1)          &
                                    & *tauasc_mix(i,k,ns)               &
                                    & *tauxar_mix(i,k,ns)
 
!                 clear sky dignostic
                  utaer(i,1) = utaer(i,1) + aermtot(i,k)*path
 
                  tauaer(i,1) = tauaer(i,1) + dextmix(1,ns,i4,i2,i3,i1) &
                              & *aermtot(i,k)*path*1D5
                  waer(i,1) = waer(i,1) + dssamix(1,ns,i4,i2,i3,i1)     &
                            & *aermtot(i,k)*path
                  gaer(i,1) = gaer(i,1) + dgmix(1,ns,i4,i2,i3,i1)       &
                            & *aermtot(i,k)*path
                  faer(i,1) = gaer(i,1) + dgmix(1,ns,i4,i2,i3,i1)       &
                            & *dgmix(1,ns,i4,i2,i3,i1)*aermtot(i,k)*path
 
                end if ! end minimum concentration conditions
              end do ! end k loop
 
              if ( utaer(i,1).gt.1.D-12 ) then
                waer(i,1) = waer(i,1)/utaer(i,1)
                gaer(i,1) = gaer(i,1)/utaer(i,1)
                faer(i,1) = faer(i,1)/utaer(i,1)
              end if
!             clear sky final effective optical properties
 
              tauxar_mix_cs(i,ns) = tauaer(i,1)
              tauasc_mix_cs(i,ns) = waer(i,1)*tauaer(i,1)
              gtota_mix_cs(i,ns) = gaer(i,1)*waer(i,1)*tauaer(i,1)
              ftota_mix_cs(i,ns) = faer(i,1)*waer(i,1)*tauaer(i,1)
            end do ! end i loop
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          else
          end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        end if
 
!       end spectral loop
      end do
 
      end subroutine aeroppt
