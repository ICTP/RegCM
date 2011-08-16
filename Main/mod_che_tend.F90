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
  use mod_runparams
  use mod_atm_interface
  use mod_che_interface
  use mod_che_trac
  use mod_pbldim
  use mod_bats
  use mod_rad
  use mod_pmoist
  use mod_che_dust
  use mod_message
  use mod_che_semdde
  use mod_diffusion
  use mod_advection
  use mod_diagnosis
  use mod_slice
  use mod_tcm_interface
  use mod_che_indices
  use mod_mppio
  private

  public :: tractend2 , tracbud

  real(8) , parameter :: d10e6 = 1.0D+06

  contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! This subroutine computes the tendencies for tracer transport and
! chemistry
!
! j:             index of j slice in current computation
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine tractend2(j,xkc)
    implicit none
!
    integer , intent(in) :: j
    real(8) , dimension(iy,kz,jxp) :: xkc
!
    real(8) :: agct , ak00t , ak0tm , akval , chimol , cldno , clmin ,&
               facb , facs , fact , facv , oh1 , pres10 , qsat10 ,    &
               remcum , rxs1 , rxs11 , rxs2 , rxs21 , satvp , shu10 , &
               u10 , v10
#ifdef CHEMTEST
    real(8) :: h2o2mol
#endif
    real(8) , dimension(ntr) :: agingtend , wetrem , wetrem_cvc
    real(8) , dimension(iy,kz) :: concmin , cutend_dwd , cutend_up ,  &
                                  fracloud , fracum , rho , settend , &
                                  ttb , wk , wl
    logical :: gfcall , gfcall2 , gfcall3 , gfcall4
    integer :: i , ibin , itr , k , kb , kdwd
    integer , dimension(iy) :: ivegcov
    real(8) , dimension(iy,kz,nbin) :: pdepv
    real(8) , dimension(iy) :: psurf , rh10 , soilw , srad , temp10 , &
                               tsurf , vegfrac , wid10 , zeff , ustar
    real(8) , dimension(iy,nbin) :: rsfrow
!
!   real(kind=8) :: ustar(iy)
!   real(kind=8) :: zza(iy,kz)
!
!   clmin = non-precipitating cloud conversion threshold,
!   clmin=0.01g/m3
    clmin = 0.01D0
!   remcum= removal rate for cumulus cloud scavenging (s-1)
    remcum = 1.0D-3
!
!   Preliminary calculations independant of tracer nature
!   the unit: rho - kg/m3, wl - g/m3
    do k = 1 , kz
      do i = 2 , iym2
        rho(i,k) = (sps2%ps(i,j)*a(k)+r8pt)*d_1000 / &
                    287.0D0/atm2%t(i,k,j)*sps2%ps(i,j)
        wl(i,k) = atm2%qc(i,k,j)/sps2%ps(i,j)*d_1000*rho(i,k)
      end do
    end do
   
    ivegcov = 0
   
!   cloud fractionnal cover for wet deposition
!   large scale : fracloud, calculated from fcc coming from pcp.f
   
!   cumulus scale : fracum, calculated from the total cloud fraction
   
!   (as defined for the radiation scheme in cldfrac.f routine)
    do i = 2 , iym2
      do k = 1 , kz
        fracloud(i,k) = dmin1(fcc(i,k,j),fcmax)
        fracum(i,k) = d_zero
      end do
      if ( icumtop(i,j) > 0 ) then
        do k = icumtop(i,j) , kz
          fracum(i,k) = cldfra(i,k) - fracloud(i,k)
        end do
      end if
    end do
   
   
!   TRANSPORT OF TRACERS
!   initialize tracer tendencies, and scratch arrays
    do itr = 1 , ntr
      do k = 1 , kz
        do i = 1 , iym1
          chiten(i,k,j,itr) = d_zero
        end do
      end do
    end do
   
!   horizontal and vertical advection
   
    do itr = 1 , ntr
!
      call hadv_x(chiten(:,:,j,itr),chi(:,:,:,itr),dx,j,2)
      call vadv(chiten(:,:,j,itr),qdot,chia(:,:,j,itr),j,5,kpbl(:,j))
!     horizontal diffusion: initialize scratch vars to 0.
!     need to compute tracer tendencies due to diffusion
      call diffu_x(chiten(:,:,j,itr),atms%chib3d(:,:,:,itr), &
                   sps2%ps,xkc(:,:,j),j,kz)
   
    end do ! end tracer loop
   
!   subgrid vertical transport by convective mass flux : a modifier !
   
    if ( ichcumtra == 2 ) then
      do k = 2 , kz
        do i = 2 , iym2
          wk(i,k) = (d_one/sps1%ps(i,j))*(twt(k,1)*chib(i,k,j,itr)+ &
                      twt(k,2)*chib(i,k-1,j,itr))
          cutend_up(i,k) = d_zero
          cutend_dwd(i,k) = d_zero
        end do
      end do
   
      do i = 2 , iym2
   
        if ( icumtop(i,j) > 0 ) then
   
          kt = max0(icumtop(i,j),3)
          kb = icumbot(i,j)
          kdwd = icumdwd(i,j)
   
!         cutend(i,kt) =  mflx(i) * g * 1.e-3*
!                         (wk(i,kb)-wk(i,kt))/(sigma(kb)-sigma(kt))
   
!         transport linked to updraft
!         betwwen kt et kdwd , the tendancy is averaged (mixing)
   
          if ( kdwd < kt ) then
            write (aline, *) 'Problem in tractend2 !'
            call say(myid)
          end if
          do k = kt , kdwd
            cutend_up(i,k) = mflx(i,1)*egrav*d_r1000*wk(i,kb)/(sigma(kdwd)-sigma(kt))
          end do
          cutend_up(i,kb) = -mflx(i,1)*egrav*d_r1000*wk(i,kb)/(dsigma(kb))

!         transport linked to downdraft
   
          cutend_dwd(i,kdwd) = -mflx(i,2)*egrav*d_r1000*wk(i,kdwd)/(dsigma(kdwd))
          cutend_dwd(i,kz) = +mflx(i,2)*egrav*d_r1000*wk(i,kdwd)/(dsigma(kz))
          do k = kt , kz
            chiten(i,k,j,itr) = chiten(i,k,j,itr) + cutend_up(i,k)+cutend_dwd(i,k)
          end do
        end if
      end do
    end if
   
!
!   SOURCE AND SINKS TERMS ( dependant on the nature of tracers)
    gfcall = .true.  ! logical call for SOx
    gfcall2 = .true. ! logical for DUST
    gfcall3 = .true. ! logical for BC aging
    gfcall4 = .true. ! logical for OC aging
   
    ibin = 0
!*******************************************
!   begining of tracer loop
!*******************************************
    do itr = 1 , ntr
!****************************
!     print*,'gfcall', itr, gfcall
   
!
!---------------------------------------------
!     SOX CHEMSITRY / IN-CLOUD PROCESS
!--------------------------------------------
      if ( gfcall .and. &
          (chtrname(itr) == 'SO2' .and. iso4 > 0) ) then
        gfcall = .false.
   
!
!       GAZEOUS CONVERSION
   
!       calculate the k value for the gaseous conversion by OH
   
!       ohconc
   
        do k = 1 , kz
          do i = 2 , iym2
            cldno = d_one    ! no cloud fraction
   
!           if (coszrs(i,j) < 0.001) ohconc(i,j,k)=ohconc(i,j,k)*0.01
!           oh1=ohconc(i,j,k)*rho(i,k)*2.084e13             !
!           molecules/cm3 test j eprends directement une valeur de oh1
#ifdef CHEMTEST 
            oh1 = oh(i,k,j)                            ! molecules/cm3
#else
            oh1 = 15.0D5
#endif
            if ( coszrs(i) < 0.001D0 ) oh1 = oh1*0.01D0 ! diurnal evolution
   
            ak0tm = 3.0D-31*(atm1%t(i,k,j)/sps1%ps(i,j) / &
                    300.0D0)**(-3.3D0)*rho(i,k)*2.084D19 ! k0(T)*[M]
            ak00t = 1.5D-12                            ! K00(T)
            akval = ak0tm/(d_one+ak0tm/ak00t) * &
                    0.6D0**(d_one/(d_one+(dlog10(ak0tm/ak00t))**d_two))
   
!           tendencies
!           here p1 unit: s^-1  and the ratio of molar mass of SO4 to
   
!           SO2 is 96/64 = 1.5
            chiten(i,k,j,iso2) = chiten(i,k,j,iso2) - &
                                 chib(i,k,j,iso2)*akval*oh1*cldno
            chiten(i,k,j,iso4) = chiten(i,k,j,iso4) + &
                                 chib(i,k,j,iso2)*akval*oh1*cldno*1.5D0
   
!           gazeous conversion diagnostic
   
            rxsg(i,k,j,iso2) = rxsg(i,k,j,iso2) + chib(i,k,j,iso2) * &
                               akval*oh1*cldno*dt*d_half
            rxsg(i,k,j,iso4) = rxsg(i,k,j,iso4) + chib(i,k,j,iso2) * &
                               akval*oh1*cldno*1.5D0*dt*d_half
   
          end do
        end do
   
!       AQUEOUS CONVERSION IN CLOUDS AND WET REMOVAL
   
!       Aqueous conversion from so2 to so4 ;control by h2o2
   
        do k = 1 , kz
          do i = 2 , iym2
#ifdef CHEMTEST
            h2o2mol = h2o2(i,k,j)
#endif
            chimol = 28.9D0/64.0D0 * chib(i,k,j,iso2)/sps2%ps(i,j) ! kg/kg to mole
!           concmin(i,k)=dmin1(h2o2mol,chimol)*64./28.9*sps2%ps(i,j)  !
!           cb*kg/kg do tests, suppose h2o2 always enough
            concmin(i,k) = chimol*64.0D0/28.9D0*sps2%ps(i,j)     ! cb*kg/kg
          end do
        end do
   
!       Large scale clouds
   
        do k = 1 , kz
          do i = 2 , iym2
            rxs1 = d_zero
            rxs11 = d_zero ! fraction of conversion, not removed, as SO4 src
            wetrem(iso2) = d_zero ! scavenging for SO2, below lsc
            wetrem(iso4) = d_zero
   
            if ( wl(i,k) > clmin ) then
!             conversion from so2 to so4
              rxs1 = fracloud(i,k)*chtrsol(iso2)*concmin(i,k) * &
                     (dexp(-wl(i,k)/360.0D0*dt)-d_one)
   
              rxs11 = rxs1*1.5D0
              ! SO4 src term and the ratio of molar
              ! mass of SO4 to SO2 is 96/64 = 1.5
   
!             if removal occurs, a fraction of SO4 src term is also
!             removed and accounted for in the term  wetrem(iso4)
   
              if ( remrat(i,k) > d_zero ) then
                wetrem(iso4) = (fracloud(i,k)*chtrsol(iso4)* &
                     chib(i,k,j,iso4)-rxs11)*(dexp(-remrat(i,k)*dt)-d_one)
              end if
            end if
   
!           Below cloud scavenging only for SO2
   
            if ( rembc(i,k) > d_zero ) then
              wetrem(iso2) = fracloud(i,k) * &
                    chtrsol(iso2)*concmin(i,k)*(dexp(-rembc(i,k)*dt)-d_one)
            end if
   
!           Tendancies large scale cloud
            chiten(i,k,j,iso2) = chiten(i,k,j,iso2) + rxs1/dt + wetrem(iso2)/dt
            chiten(i,k,j,iso4) = chiten(i,k,j,iso4) - rxs11/dt + wetrem(iso4)/dt
   
!           and wetdep diagnostics
            remlsc(i,k,j,iso2) = remlsc(i,k,j,iso2)-wetrem(iso2)*d_half
            remlsc(i,k,j,iso4) = remlsc(i,k,j,iso4)-wetrem(iso4)*d_half
   
!           chemical aqueous conversion diagnostic
            rxsaq1(i,k,j,iso2) = rxsaq1(i,k,j,iso2) - rxs1*d_half
            rxsaq1(i,k,j,iso4) = rxsaq1(i,k,j,iso4) - rxs11*d_half
   
          end do
        end do
   
!       cumulus clouds
!       wet removal by cumulus clouds (over the fraction of grid box
!       fracum) assume the cloud water content = 2 g/m3  (ref.
!       Kasibhatla )
        do i = 2 , iym2
          if ( icumtop(i,j) > 0 ) then
            do k = icumtop(i,j) , kz
              rxs2 = d_zero
              rxs21 = d_zero  ! fraction of conversion, not removed, as SO4 src
              wetrem_cvc(iso2) = d_zero   ! scavenging for SO2, below lsc
              wetrem_cvc(iso4) = d_zero
   
!             conversion from so2 to so4
              rxs2 = fracum(i,k)*chtrsol(iso2)*concmin(i,k) * &
                     (dexp(-d_two/360.0D0*dt)-d_one)
              rxs21 = rxs2*1.5D0
   
!             removal (including theremoval on the rxs21 term)
              wetrem_cvc(iso4) = (fracum(i,k)*chtrsol(iso4)* &
                      chib(i,k,j,iso4)-rxs21)*(dexp(-remcum*dt)-d_one)
   
!             tendancies due to convective processes
              chiten(i,k,j,iso2) = chiten(i,k,j,iso2) + rxs2/dt
              chiten(i,k,j,iso4) = chiten(i,k,j,iso4) + wetrem_cvc(iso4)/dt - rxs21/dt
   
!             diagnostic of wet deposition
!             remcvc(i,k,j,1) = remcvc(i,k,j,1)-wetrem_cvc(iso2)/2.0D0
              remcvc(i,k,j,iso4) = remcvc(i,k,j,iso4)  - wetrem_cvc(iso4)*d_half
!             chemical aquesous conversion diagnostic
              rxsaq2(i,k,j,iso2) = rxsaq2(i,k,j,iso2) - rxs2*d_half
              rxsaq2(i,k,j,iso4) = rxsaq2(i,k,j,iso4) - rxs21*d_half
   
            end do
          end if
        end do
   
      end if ! end of SOX chemistry
   
!---------------------------------------
!     Other than sulfate CARBON AEROSOL, DUST
!----------------------------------------
      if ( chtrname(itr) == 'BC_HB' .or. &
           chtrname(itr) == 'BC_HL' .or. &
           chtrname(itr) == 'OC_HB' .or. &
           chtrname(itr) == 'OC_HL' .or. &
           chtrname(itr) == 'DUST' ) then
   
!       wet deposition term
   
        if ( ichremlsc == 1 ) then
!         Wet removal at resolvable scale (fcc)
!         add non-precipitating cloud conversion (threshold
!         clmin=0.01g/m3) the same as that in subroutine exmois
!         clmin = 0.01
   
          do k = 1 , kz
            do i = 2 , iym2
              if ( wl(i,k) > clmin ) then
                wetrem(itr) = d_zero
                if ( remrat(i,k) > d_zero ) then
                  wetrem(itr) = fracloud(i,k)*chtrsol(itr) * &
                                chib(i,k,j,itr) * (dexp(-remrat(i,k)*dt)-d_one)
                  chiten(i,k,j,itr) = chiten(i,k,j,itr) + wetrem(itr)/dt
                  remlsc(i,k,j,itr) = remlsc(i,k,j,itr) - wetrem(itr)*d_half
                end if
              end if
            end do
          end do
        end if
   
        if ( ichremcvc == 1 ) then
!         sub-scale wet removal, cumulus cloud (fracum)
!         remcum = removal rate for cumulus cloud scavenging (s-1)
!         remcum = 1.e-3
          do i = 2 , iym2
            if ( icumtop(i,j) > 0 ) then
              do k = icumtop(i,j) , kz
                wetrem_cvc(itr) = fracum(i,k)*chtrsol(itr) * &
                                  chib(i,k,j,itr) * (dexp(-remcum*dt)-d_one)
                chiten(i,k,j,itr) = chiten(i,k,j,itr) + wetrem_cvc(itr)/dt
                remcvc(i,k,j,itr) = remcvc(i,k,j,itr) - wetrem_cvc(itr)*d_half
              end do
            end if
          end do
        end if
   
      end if ! end wet removal DUST, CARBON
   
!     Conversion from hydrophobic to hydrophilic: Carbonaceopus
!     species time constant ( 1.15 day cooke et al.,1999)
   
      if ( gfcall3 .and. chtrname(itr) == 'BC_HB' .and. ibchl > 0 )  then
        gfcall3 = .false.
        agct = 1.15D0*secpd
!       agct = 2.30D0*secpd
!bxq    do itr = 1 , ntr
!bxq      agingtend(itr) = d_zero
!bxq    end do
   
        do k = 1 , kz
          do i = 2 , iym2
            agingtend(ibchb) = -chib(i,k,j,ibchb)*(d_one-dexp(-dt/agct))/dt
            agingtend(ibchl) = -agingtend(ibchb)
            chiten(i,k,j,ibchb) = chiten(i,k,j,ibchb) + agingtend(ibchb)
            chiten(i,k,j,ibchl) = chiten(i,k,j,ibchl) + agingtend(ibchl)
          end do
        end do
      end if
   
   
      if ( gfcall4 .and. chtrname(itr) == 'OC_HB' .and. iochl > 0 ) then
   
        gfcall4 = .false.
        agct = 1.15D0*secpd
!       agct = 2.30D0*secpd
!bxq    do itr = 1 , ntr
!bxq      agingtend(itr) = d_zero
!bxq    end do
   
        do k = 1 , kz
          do i = 2 , iym2
            agingtend(iochb) = -chib(i,k,j,iochb)*(d_one-dexp(-dt/agct))/dt
            agingtend(iochl) = -agingtend(iochb)
            chiten(i,k,j,iochb) = chiten(i,k,j,iochb) + agingtend(iochb)
            chiten(i,k,j,iochl) = chiten(i,k,j,iochl) + agingtend(iochl)
          end do
        end do
   
      end if ! end aging
   
!*********************************
!     SURFACE source terms
!********************************
   
!     care chiten must be consistent with chia,b,c (= chi * pstar i.e)
!     1.e3 comes from Psurf/Pstar
   
!     en chantier DUST
!
!     print*,'before SFLUX',j,maxval(SFLT)
!     calculation of 10m wind ggaffe au facteur sur le vent
!     gaffe au kl ??
   
!     1 wind at 10 m
   
!     define the bin size
   
      if ( chtrname(itr) == 'DUST' .and. gfcall2 ) then
   
        do i = 2 , iym2
          ivegcov(i) = veg2d(i,j)
          if ( ivegcov(i) == 14 .or. ivegcov(i) == 15 ) then
            ivegcov(i) = 0
          end if
          psurf(i) = sps2%ps(i,j)*d_1000 + r8pt
   
          do k = 1 , kz
            ttb(i,k) = atm2%t(i,k,j)/sps2%ps(i,j)
!           zza(i,k) = za(i,k,j)
          end do
   
!         calculate 10 M input for wind erosion and dry deposition
!         method based on bats diagnostic in routine interf.
   
          if ( (ivegcov(i) /= 0) ) then
            facv = dlog(za(i,kz,j)*d_r10)/dlog(za(i,kz,j)/rough(ivegcov(i)))
            facb = dlog(za(i,kz,j)*d_r10)/dlog(za(i,kz,j)/zlnd)
            facs = dlog(za(i,kz,j)*d_r10)/dlog(za(i,kz,j)/zsno)
   
            fact = sfracv2d(i,j)*facv+sfracb2d(i,j)*facb + sfracs2d(i,j)*facs
   
!           grid level effective roughness lenght
!           (linear averaging for now)

            zeff(i) = rough(ivegcov(i))*sfracv2d(i,j) + &
                      zlnd * sfracb2d(i,j) + zsno * sfracs2d(i,j)

          else
!           water surface
            fact = dlog(za(i,kz,j)*d_r10)/dlog(za(i,kz,j)/zoce)
            zeff(i) = zoce
          end if
   
!         10 m wind
          u10 = (atm2%u(i,kz,j)/sps2%ps(i,j))*(d_one-fact)
          v10 = (atm2%v(i,kz,j)/sps2%ps(i,j))*(d_one-fact)
          wid10(i) = dsqrt(u10**d_two+v10**d_two)
!         wid10(5) = 15
!         10 m air temperature
   
          temp10(i) = ttb(i,kz) - sdeltk2d(i,j)*fact
   
!         specific  humidity at 10m
          shu10 = (atm2%qv(i,kz,j)/sps2%ps(i,j))/ &
            (d_one+atm2%qv(i,kz,j)/sps2%ps(i,j))-sdelqk2d(i,j)*fact
   
!         retransform in mixing ratio
   
          shu10 = shu10/(d_one-shu10)
   
!         saturation mixing ratio at 10m
          if ( temp10(i) > tzero ) then
            satvp = svp1*d_1000*dexp(svp2*(temp10(i)-tzero)/(temp10(i)-svp3))
          else
            satvp = svp4*d_1000*dexp(svp5-svp6/temp10(i))
          end if
          pres10 = psurf(i) - 98.0D0
          qsat10 = ep2*satvp/(pres10-satvp)
   
!         relative humidity at 10m
          rh10(i) = d_zero
          if ( qsat10 > d_zero ) rh10(i) = shu10/qsat10
!
!         friction velocity ( not used fo the moment)
!
          ustar(i) = sqrt ( sfsta%uvdrag(i,j)             *           &
                     sqrt ( (atm2%u(i,kz,j)/sps2%ps(i,j) )**d_two   + &
                            (atm2%v(i,kz,j)/sps2%ps(i,j) )**d_two ) / &
                             rho(i,kz) )
   
!         soil wetness
!         soilw(i) = ssw2da(i,j)                                       &
!                    /(xmopor(iexsol(idnint(mddom%lndcat(i,j))))*depuv &
!                    (idnint(mddom%lndcat(i,j))))
   
          soilw(i) = ssw2da(i,j)/depuv(idnint(mddom%lndcat(i,j)))/(2650.0D0 * &
                (d_one-xmopor(iexsol(idnint(mddom%lndcat(i,j))))))
   
!         soilw(i) = ssw2da(i,j) /(xmopor(iexsol(ivegcov(i)) )
!                                   * depuv(ivegcov(i))      )
   
!         fraction of vegetation
          vegfrac(i) = svegfrac2d(i,j)
   
!         surface temperature
!         over land recalculated from the BATS  deltk air/ surface
!         temperature account for a composite temperature between
   
!         bare ground and vegetation
          if ( ivegcov(i) /= 0 ) then
            tsurf(i) = ttb(i,kz) - sdeltk2d(i,j)
          else
!           ocean temperature in this case
            tsurf(i) = sts2%tg(i,j)
          end if
   
!         absorbed solar radiation ( for stb criteria)
   
          srad(i) = sol2d(i,j)
   
        end do

        call sfflux(iy,2,iym2,j,ivegcov,vegfrac,ustar, &
                    zeff,soilw,wid10,rho(:,kz),dustbsiz,rsfrow)
   
        call chdrydep(iy,2,iym2,kz,1,nbin,ivegcov,ttb,rho,a,psurf,    &
                      temp10,tsurf,srad,rh10,wid10,zeff,dustbsiz,pdepv)
   
   
        gfcall2 = .false.
      end if ! end dust flux and deposition calculated for bins
!     just calculated for the first case of itr = DUST (save time man)
   
      if ( chtrname(itr) == 'DUST' ) then
!       define the corresponding index between itr and DUST bins
!       (can be different if we use DUST + other particles)
   
        ibin = ibin + 1
   
!       calculate the source tendancy
        do i = 2 , iym2
          chemsrc(i,j,idatex%month,itr) = rsfrow(i,ibin)
          chiten(i,kz,j,itr) = chiten(i,kz,j,itr) + &
                               rsfrow(i,ibin)*egrav/(dsigma(kz)*d_1000)
!         diagnostique source
          cemtr(i,j,itr) = cemtr(i,j,itr)+chemsrc(i,j,idatex%month,itr)*dt*d_half
        end do
   
!       calculate the tendancy du to gravitationnal settling and dry deposition
        do k = 2 , kz
          do i = 2 , iym2
            wk(i,k) = (d_one/sps2%ps(i,j)) * &
                      (twt(k,1)*chib(i,k,j,itr)+twt(k,2)*chib(i,k-1,j,itr))
          end do
        end do
   
!       remember PDEPV is defined for ibin which is not necessarly itr
   
        do i = 2 , iym2
          do k = 2 , kz - 1
            ! do not apply to the first level
            settend(i,k) = (wk(i,k+1)*pdepv(i,k+1,ibin) - &
                            wk(i,k)*pdepv(i,k,ibin))*egrav*d_r1000/dsigma(k)
            chiten(i,k,j,itr) = chiten(i,k,j,itr) - settend(i,k)
          end do
!
          settend(i,kz) = -(wk(i,kz)*pdepv(i,kz,ibin)*egrav*d_r1000)/dsigma(kz)
          chiten(i,kz,j,itr) = chiten(i,kz,j,itr) + settend(i,kz)
   
!         dignoctic for dry deposition
          remdrd(i,j,itr) = remdrd(i,j,itr) - settend(i,kz)*dt*d_half
        end do
     
      end if !( end calculation of dust tendancies)
   
!     Source tendencies
   
      do i = 2 , iym2
        if ( chtrname(itr) /= 'DUST' ) then
          chiten(i,kz,j,itr) = chiten(i,kz,j,itr) + &
                     chemsrc(i,j,idatex%month,itr)*egrav*0.7D0/(dsigma(kz)*d_1000)
          chiten(i,kz-1,j,itr) = chiten(i,kz-1,j,itr) + &
                     chemsrc(i,j,idatex%month,itr)*egrav*0.15D0/(dsigma(kz-1)*d_1000)
          chiten(i,kzm2,j,itr) = chiten(i,kzm2,j,itr) + &
                     chemsrc(i,j,idatex%month,itr)*egrav*0.15D0/(dsigma(kzm2)*d_1000)
!         diagnostic for source, cumul
          cemtr(i,j,itr) = cemtr(i,j,itr) + chemsrc(i,j,idatex%month,itr)*dt*d_half
        end if
      end do
!   
!     end loop on tracers
!
    end do
!
  end subroutine tractend2
!
  subroutine tracbud
!
#ifndef IBM
    use mpi
#else 
    include 'mpif.h'
#endif 
    implicit none
!
    integer :: i , itr , j , k
!
    do itr = 1 , ntr
      do j = 1 , jendx
        do i = 1 , iym1
          dtrace(i,j,itr) = d_zero
          wdlsc(i,j,itr)  = d_zero
          wdcvc(i,j,itr)  = d_zero
          wxsg(i,j,itr)   = d_zero
          wxaq(i,j,itr)   = d_zero
          ddsfc(i,j,itr)  = d_zero
        end do
      end do
    end do
   
!   tracers (unit = kg):
    do itr = 1 , ntr
      do j = 1 , jendx
        do i = 1 , iym1
          do k = 1 , kz
            dtrace(i,j,itr) = dtrace(i,j,itr) + chia(i,k,j,itr)*dsigma(k)
            wdlsc(i,j,itr) = wdlsc(i,j,itr) + remlsc(i,k,j,itr)*dsigma(k)
            wdcvc(i,j,itr) = wdcvc(i,j,itr) + remcvc(i,k,j,itr)*dsigma(k)
            wxsg(i,j,itr) = wxsg(i,j,itr) + rxsg(i,k,j,itr)*dsigma(k)
!           sum ls and conv contribution
            wxaq(i,j,itr) = wxaq(i,j,itr) + &
                            (rxsaq1(i,k,j,itr)+rxsaq2(i,k,j,itr))*dsigma(k)
          end do
          ddsfc(i,j,itr) = ddsfc(i,j,itr) + remdrd(i,j,itr)*dsigma(kz)
!         Source cumulated diag(care the unit are alredy .m-2)
          cemtrac(i,j,itr) = cemtr(i,j,itr)
        end do
      end do
    end do

#ifndef BAND
    if (debug_level > 2) call contrac
#endif

    do itr = 1 , ntr
      do j = 1 , jendx
        do i = 1 , iym1
          dtrace(i,j,itr) = d10e6*dtrace(i,j,itr)*d_1000*regrav ! unit: mg/m2
          wdlsc(i,j,itr) = d10e6*wdlsc(i,j,itr)*d_1000*regrav
          wdcvc(i,j,itr) = d10e6*wdcvc(i,j,itr)*d_1000*regrav
          ddsfc(i,j,itr) = d10e6*ddsfc(i,j,itr)*d_1000*regrav
          wxsg(i,j,itr)  = d10e6*wxsg(i,j,itr)*d_1000*regrav
          wxaq(i,j,itr)  = d10e6*wxaq(i,j,itr)*d_1000*regrav
!         emtrac isbuilt from chsurfem so just need the 1e6*dt/2
!         factor to to pass im mg/m2
          cemtrac(i,j,itr) = d10e6*cemtrac(i,j,itr)
        end do
      end do
    end do

  end subroutine tracbud
!
end module mod_che_tend
