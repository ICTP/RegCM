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

  module mod_che_tend
!
! Tendency and budget for tracer transport and chemicals
!
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams , only : iqv , iqc , syncro_che
  !use mod_runparams , only : rcmtimer , syncro_srf
  use mod_mppparam
  use mod_che_common
  use mod_che_indices
  use mod_che_param
  use mod_che_sox
  use mod_che_drydep
  use mod_che_wetdep
  use mod_che_emission
  use mod_che_dust
  use mod_che_seasalt
  use mod_che_carbonaer
  use mod_che_mppio
  use mod_che_chemistry
  use mod_che_isorropia
  use mod_che_pollen
  use mod_che_bionit
  use mod_che_ccn
  use mod_che_linox
  implicit none

  private

  public :: tractend2

  contains

!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! This subroutine computes the tendencies for tracer transport and chemistry
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine tractend2(lmonth,lday,declin)
      implicit none
      integer(ik4) , intent(in) :: lmonth , lday
      real(rk8) , intent(in) :: declin

#ifndef CLM45
      real(rkx) :: facb , facs , facv , pres10 , qsat10 , shu10
#endif
      real(rkx) :: fact , u10 , v10
      real(rkx) , dimension(jci1:jci2,kz,ici1:ici2) :: rho , ttb,  wl , prec , &
                                                      convprec
      real(rkx) , dimension(jci1:jci2,kz,ici1:ici2) :: hgt , ph
      real(rkx) , dimension(jci1:jci2,kz,ici1:ici2) :: fracloud, fracum
      integer(ik4) , dimension(jci1:jci2,ici1:ici2) :: ivegcov
      real(rkx) , dimension(jci1:jci2,kz,ntr,ici1:ici2) :: pdepv
      real(rkx) , dimension(jci1:jci2,ntr,ici1:ici2) :: ddepa ! , ddepg
      real(rkx) , dimension(jci1:jci2,ici1:ici2) :: psurf , rh10 , soilw , &
       srad , temp10 , tsurf , vegfrac , snowfrac , wid10 , zeff , hsurf 
      real(rkx) , dimension(jci1:jci2,ici1:ici2,kz) :: ncpc
      real(rkx) , dimension(jci1:jci2,kz,ntr,ici1:ici2) :: bchi
      real(rkx) , dimension(jci1:jci2,ici1:ici2) :: ustar
      real(rkx) , dimension(jci1:jci2,ici1:ici2) :: xra
      real(rkx) , dimension(ntr) :: xrho
      ! evap of l-s precip (see mod_precip.f90; [kg_h2o/kg_air/s)
      ! cum h2o vapor tendency for cum precip (kg_h2o/kg_air/s)
      real(rkx) , dimension(jci1:jci2,kz,ici1:ici2) :: chevap
!     real(rkx) , dimension(jci1:jci2,kz,ici1:ici2) :: checum
      real(rkx) , dimension (1) :: polrftab
      integer(ik4) , dimension (1) :: poltab
      integer(ik4) :: i , j , k , n , ibin
      !
      !*********************************************************************
      ! A : PRELIMINARY CALCULATIONS
      !     For historical reasons and to avoid increase memory usage,
      !
      !*********************************************************************
      !
      rho         = d_zero
      wl          = d_zero
      ttb         = d_zero
      prec        = d_zero
      convprec    = d_zero
      ustar       = d_zero
      fracloud    = d_zero
      fracum      = d_zero
      psurf       = d_zero
      ivegcov     = 0
      !
      ! the unit: rho - kg/m3, wl - g/m3
      !
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! rho(j,i,k) = (sfs%psb(j,i)*a(k)+r8pt)* &
            ! what the hell   1000./287./atm2%t(j,k,i)*sfs%psb(j,i)
            hgt(j,k,i)  = cza(j,i,k)
            ph(j,k,i)   = cpb3d(j,i,k)
            rho(j,k,i)  = crhob3d(j,i,k)
            wl(j,k,i)   = cqxb3d(j,i,k,iqc)*crhob3d(j,i,k)*d_1000
            ttb(j,k,i)  = ctb3d(j,i,k)
            ! precipiation rate is a rquired variable for deposition routines.
            ! It is directly taken as rembc (saved in precip routine) in mm/hr !
            ncpc(j,i,k) =  crembc(j,i,k) / 3600._rkx !passed in mm/s
            prec(j,k,i) = ncpc(j,i,k)
            !and the quivalent for convective prec
            convprec(j,k,i) = cconvpr(j,i,k) ! already in mm/s
          end do
        end do
      end do

      if ( count(icarb > 0) > 0 .and. carb_aging_control ) then
        call carb_prepare( )
      end if
      !
      ! cloud fractionnal cover for wet deposition
      ! large scale : fracloud, calculated from fcc coming from
      ! pcp.f = large scale fraction
      ! cumulus scale : fracum, now directly saved (== cloudfra right after
      ! call to convection scheme)
      ! N.B : the difference with 'radiative cldfra' is that cldfra has an
      ! upper threshold of 0.8 and accound
      ! both large scale and cumulus cloud)
      ! here cfcc + convcldfra can be > 1  !! Overestimation of removal ?
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2

            fracloud(j,k,i)  =  cfcc(j,i,k)
            fracum(j,k,i)    =  convcldfra (j,i,k)

!            if ( kcumtop(j,i) > 0 ) then
!              do kk = kcumtop(j,i) , kz
!                fracum(j,kk,i) = ccldfra(j,i,kk) - cfcc(j,i,kk)
!                 fracum(j,kk,i) = convcldfra (j,i,kk)
!              end do
!            end if
          end do
        end do
      end do
      !
      ! variables used for natural fluxes and deposition velocities
      !
      do i = ici1 , ici2
        do j = jci1 , jci2
          !
          ! care ocean-lake in veg2d is now back type 14-15 !!
          ! ivegcov mask is also defined when CLM45 is used
          if ( cveg2d(j,i) == 14 .or. cveg2d(j,i) == 15 ) then
            ivegcov(j,i) = 0
          else
            ivegcov(j,i) = cveg2d(j,i)
          end if
          psurf(j,i) = cps2d(j,i)
          ! incoming solar radiation (for stb criteria used to calculate
          ! aerodynamic resistance)
          srad(j,i) = csol2d(j,i)
          hsurf(j,i) = cht(j,i)
          if ( idynamic == 3 ) then
            bchi(j,:,:,i) = chemt(j,i,:,:)
          else
            bchi(j,:,:,i) = chib(j,i,:,:)
          end if
          ! fraction of vegetation
#ifdef CLM45
          vegfrac(j,i) = d_one - csfracb2d(j,i)
          ! csfracb2d points on clm dust emitting ground fraction
          ! it already accouts for snow and lake mask , and is calculated
          ! as 1 - LAI/0.3 ( cf clm doc ).
#else
          vegfrac(j,i) = cvegfrac(j,i)
          snowfrac(j,i) = csfracs2d(j,i)
#endif
        end do
      end do
      ! Roughness lenght, 10 m wind           !

      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( ivegcov(j,i) /= 0 ) then
            zeff(j,i) = czo(j,i) 
          else
            zeff(j,i) = zoce
          end if
          ustar(j,i) = custar(j,i)
          xra(j,i) = cra(j,i)
          wid10(j,i) = cw10m(j,i)
          tsurf(j,i) = ctg(j,i)
          ! use simplifcations to approximate t and hr at 10 m from
          ! first model level
          temp10(j,i) =  ctb3d(j,i,kz) + lrate * (cza(j,i,kz) - 10._rkx)
          rh10(j,i) = crhb3d(j,i,kz)
          ! when CLM45 we use directly the volumetic soil water of the
          ! first CLM layer
          ! converted to gravimetric water conrtent using a sandy soil
          ! bulk density ratio to water density of 0.45
          ! consistent with bionox data
          !soilw(j,i) = csw_vol(j,i,1)*0.45_rkx
          soilw(j,i) = cssw2da(j,i)/cdepuv(nint(clndcat(j,i))) / &
                       1.e-3_rkx/(2650.0_rkx * &
                       (d_one-cxmopor(ciexsol(nint(clndcat(j,i))))))
          !
        end do
      end do
      !
      ! END of preliminary calculations
      !
      !*****************************************************************
      ! CALCULATION OF TRACER TENDENCY
      ! (except advection and convection transport)
      !*****************************************************************
      !
      ! SOX CHEMSITRY ( from offline oxidant)
      !
      if ( iso2 > 0 .and. (iso4 > 0 .or. ih2so4 > 0) ) then
        do i = ici1 , ici2
          call chemsox(i,wl(:,:,i),fracloud(:,:,i), &
                       fracum(:,:,i),rho(:,:,i),ttb(:,:,i))
        end do
      end if
      !
      ! aging of carboneaceous aerosols
      !
      if ( (ibchb > 0 .and. nbchl > 0) .or. &
           (iochb > 0 .and. nochl > 0) .or. &
           (ism1 > 0  .and. ism2 > 0) ) then
        call aging_carb
      end if
      !
      do i = ici1 , ici2
        do j = jci1 , jci2
        ! if ( ivegcov(j,i) == 8) print*, 'HE desert',custar(j,i),ustar(j,i), wid10(j,i), cw10m(j,i) 
     !FAB TEST
        end do
      end do
      !
      ! Before emission and deposition routine set the surfecae netflux
      ! used by BL schems to zero
      !
      chifxuw = d_zero
      !
      ! NATURAL EMISSIONS FLUX and tendencies  (dust -sea salt)
      !
      if ( nbin > 0 .and. ichsursrc == 1 ) then
        if ( ichdustemd /= 3 ) then
                call sfflux(lmonth,ivegcov,vegfrac,snowfrac,ustar,zeff, &
                            soilw,wid10,crho2d,dustbsiz)
        else
          ! OPTION for using CLM45 dust emission scheme
          ! if flux calculated by clm45 / update the tendency if ichdustemd == 3
          call clm_dust_tend
        end if
      end if
      !sea salt
      if ( isslt(1) > 0 .and. ichsursrc == 1 ) then
        call sea_salt(wid10,ivegcov)
      end if
      !
      ! pollen emission
      !
      if ( ipollen > 0 ) then
        call pollen_emission(ustar,wid10,rh10,ncpc(:,:,kz),convprec(:,:,kz))
      end if
      !
      ! biogenic nox emission
      if ( ichsursrc == 1 .and. ino > 0 .and. ichbion == 1 ) then
        call soilnitro_emissions(ivegcov,wid10)
      end if
      !
      ! linox emissions
      if ( ichlinox == 1 .and. ino > 0 ) then
        call linox_em(ivegcov)
      end if

      ! update emission tendencies from external inventories
      ! handle biogenic emission fluxes coming from CLM45
      if ( ichsursrc == 1 ) then
        do i = ici1 , ici2
          call emis_tend(i,declin)
        end do
      end if
      !
      ! aerosol settling and drydep
      ! include calculation of dry dep/settling velocities and
      ! updating tendencies
      !
      pdepv(:,:,:,:) = d_zero
      ddepa(:,:,:)   = d_zero
      ! ddepg(:,:,:)   = d_zero
      if ( nbin > 0 .and. ichdrdepo > 0 ) then
        do i = ici1 , ici2
          call drydep_aero(i,nbin,idust,rhodust,ivegcov(:,i),       &
                           ttb(:,:,i),rho(:,:,i),ph(:,:,i),         &
                           temp10(:,i),tsurf(:,i),srad(:,i),        &
                           rh10(:,i),wid10(:,i),zeff(:,i),dustbed,  &
                           pdepv(:,:,:,i),ddepa(:,:,i),ustar(:,i),&
                           xra(:,i))
        end do
        ! mineralogical tracers
        if ( nmine > 0 ) then
          do n = 1 , nmine
            do i = ici1 , ici2
              call drydep_aero(i,nbin,imine(:,n),rhodust,ivegcov(:,i), &
                           ttb(:,:,i),rho(:,:,i),ph(:,:,i),            &
                           temp10(:,i),tsurf(:,i),srad(:,i),           &
                           rh10(:,i),wid10(:,i),zeff(:,i),dustbed,     &
                           pdepv(:,:,:,i),ddepa(:,:,i),ustar(:,i),   &
                           xra(:,i))

            end do
          end do
        end if
      end if
      if ( isslt(1) > 0 .and. ichdrdepo > 0 ) then
        do i = ici1 , ici2
          call drydep_aero(i,sbin,isslt,rhosslt,ivegcov(:,i),       &
                           ttb(:,:,i),rho(:,:,i),ph(:,:,i),         &
                           temp10(:,i),tsurf(:,i),srad(:,i),        &
                           rh10(:,i),wid10(:,i),zeff(:,i),ssltbed,  &
                           pdepv(:,:,:,i),ddepa(:,:,i),ustar(:,i),&
                           xra(:,i))
        end do
      end if
      if ( icarb(1) > 0 .and. ichdrdepo > 0 ) then
        ibin = count( icarb > 0 )
        do i = ici1 , ici2
          call drydep_aero(i,ibin,icarb(1:ibin),rhooc,ivegcov(:,i), &
                           ttb(:,:,i),rho(:,:,i),ph(:,:,i),         &
                           temp10(:,i),tsurf(:,i),srad(:,i),        &
                           rh10(:,i),wid10(:,i),zeff(:,i),          &
                           carbed(1:ibin),pdepv(:,:,:,i),           &
                           ddepa(:,:,i),ustar(:,i),xra(:,i))
        end do
      end if
      if ( ipollen > 0 .and. ichdrdepo > 0 ) then
        ibin = 1
        poltab(1) = ipollen
        polrftab(1) = reffpollen
        do i = ici1 , ici2
          call drydep_aero(i,ibin,poltab,rhopollen,ivegcov(:,i), &
                           ttb(:,:,i),rho(:,:,i),ph(:,:,i),      &
                           temp10(:,i),tsurf(:,i),srad(:,i),     &
                           rh10(:,i),wid10(:,i),zeff(:,i),       &
                           polrftab,pdepv(:,:,:,i),ddepa(:,:,i), &
                           ustar(:,i),xra(:,i))
        end do
      end if
      !
      ! GAS phase dry deposition velocity + tendencies
      ! option compatible with BATS and CLM
      ! dry deposition for SO2  is calculated also in non gaschem simulations
      !
      if ( (iso2 > 0 .or. igaschem == 1) .and. ichdrdepo > 0 ) then
        do i = ici1 , ici2
          call drydep_gas(i,lmonth,lday,ivegcov(:,i),rh10(:,i), &
                          srad(:,i),tsurf(:,i),prec(:,kz,i),    &
                          temp10(:,i),cxlai2d(:,i),     &
                          ustar(:,i),xra(:,i))
        end do
      end if
      !
      ! WET deposition (rainout and washout) for aerosol
      !
      if ( nbin > 0 .and. ichremlsc == 1 ) then
        xrho(1:nbin) = rhodust
        do i = ici1 , ici2
          call wetdepa(i,nbin,idust,dustbed,xrho(1:nbin),ttb(:,:,i), &
                       wl(:,:,i),fracloud(:,:,i),fracum(:,:,i),      &
                       psurf(:,i),hsigma,rho(:,:,i),prec(:,:,i),     &
                       convprec(:,:,i), pdepv(:,:,:,i))
        end do
        ! mineralogical tracers
        if ( nmine > 0 ) then
          xrho(1:nbin) = rhodust
          do n = 1 , nmine
            do i = ici1 , ici2
              call wetdepa(i,nbin,imine(:,n),dustbed,xrho(1:nbin),ttb(:,:,i), &
                           wl(:,:,i),fracloud(:,:,i),fracum(:,:,i),           &
                           psurf(:,i),hsigma,rho(:,:,i),prec(:,:,i),          &
                           convprec(:,:,i), pdepv(:,:,:,i))
            end do
          end do
        end if
      end if
      if ( isslt(1) > 0 .and. ichremlsc == 1 )  then
        xrho(1:sbin) = rhosslt
        do i = ici1 , ici2
          call wetdepa(i,sbin,isslt,ssltbed,xrho(1:sbin),ttb(:,:,i),  &
                       wl(:,:,i),fracloud(:,:,i),fracum(:,:,i),       &
                       psurf(:,i),hsigma,rho(:,:,i),prec(:,:,i),      &
                       convprec(:,:,i), pdepv(:,:,:,i))
        end do
      end if
      if ( icarb(1) > 0 .and. ichremlsc == 1 )  then
        ibin = count( icarb > 0 )
        ! Yep ?
        xrho(1:ibin) = rhobchl(1)
        do i = ici1 , ici2
          call wetdepa(i,ibin,icarb(1:ibin),carbed(1:ibin),xrho(1:ibin),   &
                       ttb(:,:,i),wl(:,:,i),fracloud(:,:,i),fracum(:,:,i), &
                       psurf(:,i),hsigma,rho(:,:,i),prec(:,:,i),           &
                       convprec(:,:,i),pdepv(:,:,:,i))
        end do
      end if
      !
      if ( ipollen > 0 .and. ichremlsc == 1 )  then
        ibin = 1
        poltab(1) = ipollen
        polrftab(1) = reffpollen
        xrho(1) = rhopollen
        do i = ici1 , ici2
          call wetdepa(i,ibin,poltab,polrftab,xrho(1:1),ttb(:,:,i), &
                       wl(:,:,i),fracloud(:,:,i),fracum(:,:,i),     &
                       psurf(:,i),hsigma,rho(:,:,i),prec(:,:,i),    &
                       convprec(:,:,i),pdepv(:,:,:,i))
        end do
      end if
      !
      ! Wet Deposition for gasphase species
      !
      if ( igaschem == 1 .and. ichremlsc == 1 ) then
        do i = ici1 , ici2
          call sethet(i,bchi(:,:,:,i), dt,psurf(:,i))
        end do
      end if
      !
      ! Gas phase solver
      ! note : solver is called every dtchsolv (900s)- chemten
      ! chemistry raecation tendendy is calculated but chemical tracer
      ! tendency is still updated every dt ( =dt) time step
      ! ( insure smoothness)
      !
      if ( igaschem == 1 .and. ichsolver > 0 ) then
        if ( syncro_che%act( ) ) then
          chemten(:,:,:,:) = d_zero
          call chemistry
          if ( myid == italk .and. syncro_rep%will_act( ) ) then
            write(stdout,'(a,2g12.5)') ' $$$ Jvalue min/max NO2surf : ', &
              minval(jphoto(:,:,kz,jvNO2 )),  maxval(jphoto(:,:,kz,jvNO2 ))
          end if
          ! secondary inorganic aerosol solver ( modify also chemten !)
          if ( iisoropia == 1 ) then
            call aerodriver
          endif
        end if ! end chem timestep

        ! add tendency due to chemistry reaction + thermo equilibrium (every dt)
        chiten(:,:,:,:) = chiten(:,:,:,:) + chemten(:,:,:,:)
        if ( ichdiag > 0 ) then
          chemdiag(:,:,:,:) = chemdiag(:,:,:,:) + chemten(:,:,:,:) * cfdout
        end if
      end if
      !
      ! FAB try to put a criteria to avoid grid point storm and explosion
      ! in the chemistry.
      ! normally grid point storm produces a crazy precipitation in the column
      ! set tracer total tracer tendency to zero if so ..
      ! threshold 0.005 mm/s = 432 mm/d
      !
      !do i = ici1 , ici2
      !  do j = jci1 , jci2
      !    if ( maxval(convprec(j,:,i) + prec(j,:,i)) > 0.001_rkx ) then
      !      chiten (j,i,:,:) = d_zero
      !    end if
      !  end do
      !end do
      !
      ! diagnostics
      ! tracer instantaneous burden for diag
      !
      dtrace(:,:,:) = d_zero
      do n = 1 , ntr
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              dtrace(j,i,n) = dtrace(j,i,n) +  &
                  chib3d(j,i,k,n)*cdzq(j,i,k)*crhob3d(j,i,k)
            end do
          end do
        end do
      end do
      ! calculate ccn number for use in precip autoconversion
      ! calculation ( 2nd indirect effect)
      if ( iindirect > 0 .and. iaerosol == 1 ) then
        call ccn
      end if

      contains

#include <pfesat.inc>
#include <pfwsat.inc>

    end subroutine tractend2

end module mod_che_tend
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
