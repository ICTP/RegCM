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
  use mod_runparams , only : iqv , iqc
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
    subroutine tractend2(ktau,lyear,lmonth,lday,calday,declin)
      implicit none
      integer(ik4) , intent(in) :: lmonth , lday , lyear
      real(rk8) , intent(in) :: calday , declin
      integer(ik8) , intent(in) :: ktau

      real(rk8) :: facb , facs , fact , facv , pres10 , qsat10 , &
                  shu10 , u10 , v10
      real(rk8) , dimension(ici1:ici2,kz,jci1:jci2) :: rho , ttb,  wl , prec , &
                                                      convprec
      real(rk8) , dimension(ici1:ici2,kz,jci1:jci2) :: hgt
      real(rk8) , dimension(ici1:ici2,kz,jci1:jci2) :: fracloud, fracum
      integer(ik4) , dimension(ici1:ici2,jci1:jci2) :: ivegcov
      real(rk8) , dimension(ici1:ici2,kz,ntr,jci1:jci2) :: pdepv
      real(rk8) , dimension(ici1:ici2,ntr,jci1:jci2) :: ddepa ! , ddepg
      real(rk8) , dimension(ici1:ici2,jci1:jci2) :: psurf , rh10 , soilw , &
                 srad , temp10 , tsurf , vegfrac , wid10 , zeff , ustar , &
                 hsurf
      real(rk8) , dimension(ici1:ici2,kz,ntr,jci1:jci2) :: bchi
      real(rk8) , dimension(ici1:ici2,1) :: xra
      real(rk8) , dimension(ici1:ici2,nbin,jci1:jci2) :: dust_flx
      real(rk8) , dimension(ici1:ici2,sbin,jci1:jci2) :: seasalt_flx
      ! evap of l-s precip (see mod_precip.f90; [kg_h2o/kg_air/s)
      ! cum h2o vapor tendency for cum precip (kg_h2o/kg_air/s)
      real(rk8) , dimension(ici1:ici2,kz,jci1:jci2) :: chevap
!     real(rk8) , dimension(ici1:ici2,kz,jci1:jci2) :: checum
      real(rk8) , dimension (1) :: polrftab
      integer(ik4) , dimension (1) :: poltab
      integer(ik4) :: i , j , ibin , k
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
      dust_flx    = d_zero
      seasalt_flx = d_zero
      ivegcov     = 0
      !
      ! the unit: rho - kg/m3, wl - g/m3
      !
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! rho(j,i,k) = (sfs%psb(i,j)*a(k)+r8pt)* &
            ! what the hell   1000./287./atm2%t(i,k,j)*sfs%psb(i,j)
            hgt(i,k,j)  = cza(j,i,k)
            rho(i,k,j)  = crhob3d(j,i,k)
            wl(i,k,j)   = cqxb3d(j,i,k,iqc)*crhob3d(j,i,k)*d_1000
            ttb(i,k,j)  = ctb3d(j,i,k)
            ! precipiation rate is a rquired variable for deposition routines.
            ! It is directly taken as rembc (saved in precip routine) in mm/hr !
            prec(i,k,j) = crembc(j,i,k) / 3600.D0 !passed in mm/s
            !and the quivalent for convective prec
            convprec(i,k,j) = cconvpr(j,i,k) ! already in mm/s
          end do
        end do
      end do
      !
      ! cloud fractionnal cover for wet deposition
      ! large scale : fracloud, calculated from fcc coming from pcp.f = large scale fraction
      ! cumulus scale : fracum, now directly saved (== cloudfra right after call to convection scheme)
      ! N.B : the difference with 'radiative cldfra' is that cldfra has an upper threshold of 0.8 and accound
      ! both large scale and cumulus cloud)
      ! here cfcc + convcldfra can be > 1  !! Overestimation of removal ?
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2

            fracloud(i,k,j)  =  cfcc(j,i,k)
            fracum(i,k,j)    =  convcldfra (j,i,k)

!            if ( kcumtop(j,i) > 0 ) then
!              do kk = kcumtop(j,i) , kz
!                fracum(i,kk,j) = ccldfra(j,i,kk) - cfcc(j,i,kk)
!                 fracum(i,kk,j) = convcldfra (j,i,kk)
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
          !
          if ( cveg2d(j,i) == 14 .or. cveg2d(j,i) == 15 ) then
            ivegcov(i,j) = 0
          else
            ivegcov(i,j) = cveg2d(j,i)
          end if
          psurf(i,j) = cps2d(j,i)
          !
          ! method based on bats diagnostic in routine interf.
          !
          if ( ivegcov(i,j) /= 0 ) then
            facv = dlog(cza(j,i,kz)/d_10) / &
                   dlog(cza(j,i,kz)/crough(ivegcov(i,j)))
            facb = dlog(cza(j,i,kz)/d_10)/dlog(cza(j,i,kz)/zlnd)
            facs = dlog(cza(j,i,kz)/d_10)/dlog(cza(j,i,kz)/zsno)
            fact = csfracv2d(j,i)*facv + csfracb2d(j,i)*facb + &
                   csfracs2d(j,i)*facs
            !
            ! grid level effective roughness lenght (linear averaging for now)
            zeff(i,j) = crough(ivegcov(i,j))*csfracv2d(j,i) + &
                       zlnd * csfracb2d(j,i) + zsno * csfracs2d(j,i)
          else
            ! water surface
            fact = dlog(cza(j,i,kz)/d_10)/dlog(cza(j,i,kz)/zoce)
            zeff(i,j) = zoce
          end if
          ! 10 m wind
          u10 = (cubx3d(j,i,kz))*(d_one-fact)
          v10 = (cvbx3d(j,i,kz))*(d_one-fact)
          wid10(i,j) = sqrt(u10**2+v10**2)
          ! 10 m air temperature
          temp10(i,j) = ctb3d(j,i,kz) - csdeltk2d(j,i)*fact
          ! specific  humidity at 10m
          shu10 = cqxb3d(j,i,kz,iqv)/(d_one+cqxb3d(j,i,kz,iqv)) - &
                  csdelqk2d(j,i)*fact
          ! back to mixing ratio
          shu10 = shu10/(d_one-shu10)
          ! saturation mixing ratio at 10m
          pres10 = psurf(i,j) - 98.0D0
          qsat10 = pfqsat(temp10(i,j),pres10)
          ! relative humidity at 10m
          rh10(i,j) = d_zero
          if ( qsat10 > d_zero ) rh10(i,j) = shu10/qsat10
          !
          ! soil wetness
          soilw(i,j) = cssw2da(j,i)/cdepuv(idnint(clndcat(j,i)))/(2650.0D0 * &
                       (d_one-cxmopor(ciexsol(idnint(clndcat(j,i))))))
          ! fraction of vegetation
          vegfrac(i,j) = cvegfrac(j,i)
          ! surface temperature
          ! over land recalculated from the BATS as deltk air/ surface
          ! temperature account for a composite temperature between
          ! bare ground and vegetation
          if ( ivegcov(i,j) /= 0 ) then
            tsurf(i,j) = ctb3d(j,i,kz) - csdeltk2d(j,i)
          else
            ! ocean temperature in this case
            tsurf(i,j) = ctg(j,i)
          end if
          !
          ! aborbed solar radiation (for stb criteria used to calculate
          ! aerodynamic resistance)
          !
          srad(i,j) = csol2d(j,i)
          hsurf(i,j) = cht(j,i)
          bchi(i,:,:,j) = chib(j,i,:,:)
        end do
      end do
      !
      ! END of preliminary calculations)
      !
      !*****************************************************************
      ! CALCULATION OF TRACER TENDENCY
      !       (except full gas phase chemistry solver)
      !*****************************************************************
      !
      ! SOX CHEMSITRY ( from offline oxidant)
      !
      if ( igaschem == 0 ) then
        if ( iso2 > 0 .and. iso4 > 0 ) then
          do j = jci1 , jci2
            call chemsox(j,wl(:,:,j),fracloud(:,:,j), &
                         fracum(:,:,j),rho(:,:,j),ttb(:,:,j))
          end do
        end if
      end if
      !
      ! aging of carboneaceous aerosols
      !
      if ( (ibchb > 0 .and. ibchl > 0 ) .or. &
           (iochb > 0 .and. iochl > 0) ) then
        do j = jci1 , jci2
          call aging_carb(j)
        end do
      end if
      !
      ! Before emission and deposition routine set the surfecae netflux
      ! used by BL schems to zero
      !
      chifxuw = d_zero
      !
      ! NATURAL EMISSIONS FLUX and tendencies  (dust -sea salt)
      !
      if ( idust(1) > 0 .and. ichdustemd < 3 .and.  ichsursrc == 1 ) then
        do j = jci1 , jci2
          where (ivegcov(:,j) == 11)
            zeff(:,j) = 0.01D0 ! value set to desert type for semi-arid)
          end where
          call aerodyresis(zeff(:,j),wid10(:,j),temp10(:,j),tsurf(:,j), &
            rh10(:,j),srad(:,j),ivegcov(:,j),ustar(:,j),xra(:,1))
          call sfflux(j,ivegcov(:,j),vegfrac(:,j),ustar(:,j),      &
                      zeff(:,j),soilw(:,j),wid10(:,j),rho(:,kz,j), &
                      dustbsiz,dust_flx(:,:,j))
        end do
      end if
      ! if flux calculated by clm45 / update the tendency if ichdustemd == 3
#if defined CLM45
      if (idust(1) > 0 .and.  ichdustemd == 3 .and. ichsursrc == 1 ) then
        if ( ktau == 0 .or. mod(ktau+1,ntsrf) == 0 ) call clm_dust_tend
      end if
#endif
      !sea salt
      if ( isslt(1) > 0 .and. ichsursrc == 1 ) then
        do j = jci1 , jci2
          call sea_salt(j,wid10(:,j),ivegcov(:,j),seasalt_flx(:,:,j))
        end do
      end if
      !
      ! pollen emission
      !
      if ( ipollen > 0 ) then
        do j = jci1 , jci2
         call aerodyresis(zeff(:,j),wid10(:,j),temp10(:,j),tsurf(:,j), &
           rh10(:,j),srad(:,j),ivegcov(:,j),ustar(:,j),xra(:,1))
         call pollen_emission(j,ustar(:,j),wid10(:,j),rh10(:,j), &
           prec(:,kz,j), convprec(:,kz,j))
        end do
      end if
      !
      ! biogenic nox emission
      if ( ichsursrc == 1 .and. ino > 0 .and. ichbion == 1 ) then
       do j = jci1 , jci2
          call soilnitro_emissions(j,ivegcov(:,j),wid10(:,j))
       end do
      end if
      !
      ! update emission tendencies from inventories
      !
      if ( ichsursrc == 1 ) then
        do j = jci1 , jci2
          call emis_tend(ktau,j,lmonth,declin)
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
      if ( idust(1) > 0 .and. ichdrdepo > 0 ) then
        do j = jci1 , jci2
          call drydep_aero(j,nbin,idust,rhodust,ivegcov(:,j),      &
                           ttb(:,:,j),rho(:,:,j),hsigma,psurf(:,j),  &
                           temp10(:,j),tsurf(:,j),srad(:,j),       &
                           rh10(:,j),wid10(:,j),zeff(:,j),dustbed, &
                           pdepv(:,:,:,j),ddepa(:,:,j))
        end do
      end if
      if ( isslt(1) > 0  .and. ichdrdepo > 0 ) then
        do j = jci1 , jci2
          call drydep_aero(j,sbin,isslt,rhosslt,ivegcov(:,j),      &
                           ttb(:,:,j),rho(:,:,j),hsigma,psurf(:,j),  &
                           temp10(:,j),tsurf(:,j),srad(:,j),       &
                           rh10(:,j),wid10(:,j),zeff(:,j),ssltbed, &
                           pdepv(:,:,:,j),ddepa(:,:,j))
        end do
      end if
      if ( icarb(1) > 0  .and. ichdrdepo > 0 ) then
        ibin = count( icarb > 0 )
        do j = jci1 , jci2
          call drydep_aero(j,ibin,icarb(1:ibin),rhooc,ivegcov(:,j), &
                           ttb(:,:,j),rho(:,:,j),hsigma,psurf(:,j),   &
                           temp10(:,j),tsurf(:,j),srad(:,j),        &
                           rh10(:,j),wid10(:,j),zeff(:,j),          &
                           carbed(1:ibin),pdepv(:,:,:,j),ddepa(:,:,j))
        end do
      end if
      if ( ipollen > 0 .and. ichdrdepo > 0 ) then
        ibin = 1
        poltab(1) = ipollen
        polrftab(1) = reffpollen
        do j = jci1 , jci2
          call drydep_aero(j,ibin,poltab,rhopollen,ivegcov(:,j), &
                           ttb(:,:,j),rho(:,:,j),hsigma,psurf(:,j),   &
                           temp10(:,j),tsurf(:,j),srad(:,j),        &
                           rh10(:,j),wid10(:,j),zeff(:,j),          &
                           polrftab,pdepv(:,:,:,j),ddepa(:,:,j))
        end do
      end if
      !
      ! GAS phase dry deposition velocity + tendencies
      ! option compatible with BATS and CLM
      ! dry deposition for SO2  is calculated also in non gaschem simulations
      !
      if ( (iso2 > 0 .or. igaschem == 1) .and. ichdrdepo > 0 ) then
        do j = jci1 , jci2
          call drydep_gas(j,lmonth,lday,ivegcov(:,j),rh10(:,j), &
                          srad(:,j),tsurf(:,j),prec(:,kz,j),    &
                          temp10(:,j),wid10(:,j),zeff(:,j))
        end do
      end if
      !
      ! WET deposition (rainout and washout) for aerosol
      !
      if ( idust(1) > 0 .and. ichremlsc == 1 ) then
        do j = jci1 , jci2
          call wetdepa(j,nbin,idust,dustbed,rhodust,ttb(:,:,j),  &
                       wl(:,:,j),fracloud(:,:,j),fracum(:,:,j),  &
                       psurf(:,j),hsigma,rho(:,:,j),prec(:,:,j), &
                       convprec(:,:,j), pdepv(:,:,:,j))
        end do
      end if
      if ( isslt(1) > 0 .and.   ichremlsc == 1 )  then
        do j = jci1 , jci2
          call wetdepa(j,sbin,isslt,ssltbed,rhosslt,ttb(:,:,j),  &
                       wl(:,:,j),fracloud(:,:,j),fracum(:,:,j),  &
                       psurf(:,j),hsigma,rho(:,:,j),prec(:,:,j), &
                       convprec(:,:,j), pdepv(:,:,:,j))
        end do
      end if
      if ( icarb(1) > 0 .and.  ichremlsc == 1 )  then
        ibin = count( icarb > 0 )
        do j = jci1 , jci2
          call wetdepa(j,ibin,icarb(1:ibin),carbed(1:ibin),rhobchl,        &
                       ttb(:,:,j),wl(:,:,j),fracloud(:,:,j),fracum(:,:,j), &
                       psurf(:,j),hsigma,rho(:,:,j),prec(:,:,j),           &
                       convprec(:,:,j),pdepv(:,:,:,j))
        end do
      end if
      !
      if ( ipollen > 0 .and.  ichremlsc == 1 )  then
        ibin = 1
        poltab(1) = ipollen
        polrftab(1) = reffpollen
        do j = jci1 , jci2
          call wetdepa(j,ibin,poltab,polrftab,rhopollen,        &
                       ttb(:,:,j),wl(:,:,j),fracloud(:,:,j),fracum(:,:,j), &
                       psurf(:,j),hsigma,rho(:,:,j),prec(:,:,j),           &
                       convprec(:,:,j),pdepv(:,:,:,j))
        end do
      end if
      !
      ! Wet Deposition for gasphase species
      !
      if ( igaschem == 1 .and. ichremlsc == 1 ) then
        ! fix the interface for this variable
        ! no effect of cumulus scavenging now
        ! checum = d_zero
        chevap = d_zero
        do j = jci1 , jci2
          call sethet(j,hgt(:,:,j),hsurf(:,j),ttb(:,:,j),        &
                      prec(:,:,j),convprec(:,:,j),chevap(:,:,j), &
                      dt,rho(:,:,j),bchi(:,:,:,j),psurf(:,j))
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
        if ( mod(ktau+1,kchsolv) == 0 ) then
          chemten(:,:,:,:) = d_zero
          do j = jci1 , jci2
            call chemistry(j)
          end do
          if ( myid == italk .and. mod(ktau+1,krep) == 0 ) then
            write(stdout,'(a,2g12.5)') ' $$$ Jvalue min/max NO2surf : ', &
              minval(jphoto(:,:,kz,jvNO2 )),  maxval(jphoto(:,:,kz,jvNO2 ))
          end if
          ! secondary inorganic aerosol solver ( modify also chemten !)
          if ( iisoropia == 1 ) then
            call aerodriver
          endif
        end if ! end chem timestep

        ! add tendency due to chemistry reaction + thermo equilibrium (every dt)
        chiten(jci1:jci2,:,:,:) = chiten(jci1:jci2,:,:,:) + &
                                  chemten(jci1:jci2,:,:,:)
        if ( ichdiag > 0 ) then
          chemdiag(jci1:jci2,:,:,:) = chemdiag(jci1:jci2,:,:,:) + &
              chemten(jci1:jci2,:,:,:) * cfdout
        end if
      end if
      !
      ! Finally save tarcer instantaneous burden for diag
      !
      dtrace(:,:,:) = d_zero
      do k=1,kz
        do i = ici1 , ici2
          do j = jci1 , jci2
           dtrace(j,i,:) = dtrace(j,i,:) +  &
             chib3d(j,i,k,:)*cdzq(j,i,k)*crhob3d(j,i,k)
          end do
        end do
      end do

    end subroutine tractend2
!
end module mod_che_tend
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
