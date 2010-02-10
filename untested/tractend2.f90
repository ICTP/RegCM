!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine tractend2(j)
!
!     This subroutine computes the tendencies for tracer transport and
!     chemistry
!
!     ntr:           dimension of tracer arrays in species index
!     j:             index of j slice in current computation
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use regcm_param
      use param1
      use param3
      use main
      use mainchem
      use cvaria
      use trachem
      use pbldim
      use mod_bats
      use rad
      use pmoist
      use dust
      use date
      use message
      implicit none
!
! Dummy arguments
!
      integer :: j
!
! Local variables
!
      real(8) :: agct , ak00t , ak0tm , akval , chimol , cldno , clmin ,&
               & facb , facs , fact , facv , oh1 , pres10 , qsat10 ,    &
               & remcum , rxs1 , rxs11 , rxs2 , rxs21 , satvp , shu10 , &
               & u10 , v10
      real(8) , dimension(ntr) :: agingtend , wetrem , wetrem_cvc
      real(8) , dimension(ix,kx) :: concmin , cutend_dwd , cutend_up ,  &
                                  & fracloud , fracum , rho , settend , &
                                  & ttb , wk , wl
      logical :: gfcall , gfcall2 , gfcall3 , gfcall4
      integer :: i , ibin , itr , k , kb , kdwd
      integer , dimension(ix) :: ivegcov , soilt
      real(8) , dimension(ix,kx,nbin) :: pdepv
      real(8) , dimension(ix) :: psurf , rh10 , soilw , srad ,  &
                               & temp10 , tsurf , vegfrac , wid10 , zeff
      real(8) , dimension(ix,nbin) :: rsfrow
!
!bxq  real(kind=8)  h2o2mol
!     real(kind=8)  ustar(ix)
!     real(kind=8)  zza(ix,kx)
!
!     clmin = non-precipitating cloud conversion threshold,
!     clmin=0.01g/m3
      clmin = 0.01
!     remcum= removal rate for cumulus cloud scavenging (s-1)
      remcum = 1.E-3
!
!     Preliminary calculations independant of tracer nature
 
!     the unit: rho - kg/m3, wl - g/m3
      do k = 1 , kx
        do i = 2 , ilxm
          rho(i,k) = (psb(i,j)*a(k)+ptop)*1000./287./tb(i,k,j)*psb(i,j)
          wl(i,k) = qcb(i,k,j)/psb(i,j)*1000.*rho(i,k)
        end do
      end do
 
      ivegcov = 0
 
!     cloud fractionnal cover for wet deposition
!     large scale : fracloud, calculated from fcc coming from pcp.f
 
!     cumulus scale : fracum, calculated from the total cloud fraction
 
!     (as defined for the radiation scheme in cldfrac.f routine)
      do i = 2 , ilxm
        do k = 1 , kx
          fracloud(i,k) = dmin1(fcc(i,k,j),fcmax)
          fracum(i,k) = 0.
        end do
        if ( icumtop(i,j).ne.0 ) then
          do k = icumtop(i,j) , kx
            fracum(i,k) = cldfra(i,k) - fracloud(i,k)
          end do
        end if
      end do
 
 
!     TRANSPORT OF TRACERS
!----initialize tracer tendencies, and scratch arrays
      do itr = 1 , ntr
        do k = 1 , kx
          do i = 1 , ilx
            chiten(i,k,j,itr) = 0.
          end do
        end do
      end do
 
!-----horizontal and vertical advection
 
      do itr = 1 , ntr
!
        call hadvch(chiten(1,1,j,itr),dx,itr,j,2)
 
        call vadv(chiten(1,1,j,itr),chia(1,1,j,itr),j,5)
 
!       call vadv(chiten(1,1,j,itr),chia(1,1,j,itr),j,3)
 
!----horizontal diffusion: initialize scratch vars to 0.
!       need to compute tracer tendencies due to diffusion
        call diffutch(chiten(1,1,j,itr),xkc(1,1,j),c203,itr,j)
 
      end do ! end tracer loop
 
!     subgrid vertical transport by convective mass flux : a modifier !
 
      if ( ichcumtra.eq.2 ) then
        do k = 2 , kx
          do i = 2 , ilxm
            wk(i,k) = (1./psa(i,j))                                     &
                    & *(twt(k,1)*chib(i,k,j,itr)+twt(k,2)*chib(i,k-1,j, &
                    & itr))
 
            cutend_up(i,k) = 0.
            cutend_dwd(i,k) = 0.
          end do
        end do
 
        do i = 2 , ilxm
 
          if ( icumtop(i,j).ne.0 ) then
 
            kt = max0(icumtop(i,j),3)
            kb = icumbot(i,j)
            kdwd = icumdwd(i,j)
 
!           cutend(i,kt) =  mflx(i) * g * 1.e-3*
!           &               (wk(i,kb)-wk(i,kt))/(sigma(kb)-sigma(kt))
 
!           transport linked to updraft
!           betwwen kt et kdwd , the tendancy is averaged (mixing)
 
            if ( kdwd.lt.kt ) then
              write (aline, *) 'Problem in tractend2 !'
              call say
            end if
            do k = kt , kdwd
              cutend_up(i,k) = mflx(i,1)*g*1.E-3*wk(i,kb)               &
                             & /(sigma(kdwd)-sigma(kt))
            end do
 
            cutend_up(i,kb) = -mflx(i,1)*g*1.E-3*wk(i,kb)/(dsigma(kb))
!           transport linked to downdraft
 
            cutend_dwd(i,kdwd) = -mflx(i,2)*g*1.E-3*wk(i,kdwd)          &
                               & /(dsigma(kdwd))
 
            cutend_dwd(i,kx) = +mflx(i,2)*g*1.E-3*wk(i,kdwd)            &
                             & /(dsigma(kx))
 
            do k = kt , kx
              chiten(i,k,j,itr) = chiten(i,k,j,itr) + cutend_up(i,k)    &
                                & + cutend_dwd(i,k)
            end do
          end if
        end do
      end if
 
!
!     SOURCE AND SINKS TERMS ( dependant on the nature of tracers)
      gfcall = .true.  ! logical call for SOx
      gfcall2 = .true. ! logical for DUST
      gfcall3 = .true. ! logical for BC aging
      gfcall4 = .true. ! logical for OC aging
 
      ibin = 0
!*******************************************
!     begining of tracer loop
!*******************************************
      do itr = 1 , ntr
!****************************
!       print*,'gfcall', itr, gfcall
 
!
!---------------------------------------------
!       SOX CHEMSITRY / IN-CLOUD PROCESS
!--------------------------------------------
        if ( gfcall .and. (chtrname(itr).eq.'SO2' .and. iso4.gt.0) )    &
           & then
          gfcall = .false.
 
!
!         GAZEOUS CONVERSION
 
!         calculate the k value for the gaseous conversion by OH
 
!         ohconc
 
          do k = 1 , kx
            do i = 2 , ilxm
              cldno = 1.    ! no cloud fraction
 
!             if(coszrs(i,j).lt.0.001) ohconc(i,j,k)=ohconc(i,j,k)*0.01
!             oh1=ohconc(i,j,k)*rho(i,k)*2.084e13             !
!             molecules/cm3 test j eprends directement une valeur de oh1
 
              oh1 = 15.E5                                ! molecules/cm3
              if ( coszrs(i).lt.0.001 ) oh1 = oh1*0.01   ! diurnal evolution
 
              ak0tm = 3.E-31*(ta(i,k,j)/psa(i,j)/300.)**(-3.3)*rho(i,k) &
                    & *2.084E19                          ! k0(T)*[M]
              ak00t = 1.5E-12                            ! K00(T)
              akval = ak0tm/(1.+ak0tm/ak00t)                            &
                    & *0.6**(1./(1.+(dlog10(ak0tm/ak00t))**2.))
 
!             tendencies
!             here p1 unit: s^-1  and the ratio of molar mass of SO4 to
 
!             SO2 is 96/64 = 1.5
              chiten(i,k,j,iso2) = chiten(i,k,j,iso2) - chib(i,k,j,iso2)&
                                 & *akval*oh1*cldno
              chiten(i,k,j,iso4) = chiten(i,k,j,iso4) + chib(i,k,j,iso2)&
                                 & *akval*oh1*cldno*1.5
 
!             gazeous conversion diagnostic
 
              rxsg(i,k,j,iso2) = rxsg(i,k,j,iso2) + chib(i,k,j,iso2)    &
                               & *akval*oh1*cldno*dt/2.
              rxsg(i,k,j,iso4) = rxsg(i,k,j,iso4) + chib(i,k,j,iso2)    &
                               & *akval*oh1*cldno*1.5*dt/2.
 
            end do
          end do
 
!         AQUEOUS CONVERSION IN CLOUDS AND WET REMOVAL
 
!         Aqueous conversion from so2 to so4 ;control by h2o2
 
          do k = 1 , kx
            do i = 2 , ilxm
!bxq          h2o2mol = 1.e-6 * h2o2conc(i,j,k)
              chimol = 28.9/64.*chib(i,k,j,iso2)/psb(i,j)      ! kg/kg to mole
!             concmin(i,k)=dmin1(h2o2mol,chimol)*64./28.9*psb(i,j)  !
!             cb*kg/kg do tests, suppose h2o2 always enough
              concmin(i,k) = chimol*64./28.9*psb(i,j)       ! cb*kg/kg
            end do
          end do
 
!         Large scale clouds
 
          do k = 1 , kx
            do i = 2 , ilxm
              rxs1 = 0.0
              rxs11 = 0.0      ! fraction of conversion, not removed, as SO4 src
              wetrem(iso2) = 0.
                               ! scavenging for SO2, below lsc
              wetrem(iso4) = 0.
 
              if ( wl(i,k).gt.clmin ) then
!               conversion from so2 to so4
                rxs1 = fracloud(i,k)*chtrsol(iso2)*concmin(i,k)         &
                     & *(dexp(-wl(i,k)/360.*dt)-1.)
 
                rxs11 = rxs1*1.5
                                ! SO4 src term and the ratio of molar
                                ! mass of SO4 to SO2 is 96/64 = 1.5
 
!               if removal occurs, a fraction of SO4 src term is also
!               removed and accounted for in the term  wetrem(iso4)
 
                if ( remrat(i,k).gt.0. ) wetrem(iso4)                   &
                   & = (fracloud(i,k)*chtrsol(iso4)*chib(i,k,j,iso4)    &
                   & -rxs11)*(dexp(-remrat(i,k)*dt)-1.)
              end if
 
!             Below cloud scavenging only for SO2
 
              if ( rembc(i,k).gt.0. ) wetrem(iso2) = fracloud(i,k)      &
                 & *chtrsol(iso2)*concmin(i,k)*(dexp(-rembc(i,k)*dt)-1.)
 
!             Tendancies large scale cloud
              chiten(i,k,j,iso2) = chiten(i,k,j,iso2) + rxs1/dt +       &
                                 & wetrem(iso2)/dt
              chiten(i,k,j,iso4) = chiten(i,k,j,iso4) - rxs11/dt +      &
                                 & wetrem(iso4)/dt
 
!             and wetdep diagnostics
              remlsc(i,k,j,iso2) = remlsc(i,k,j,iso2) - wetrem(iso2)/2.
              remlsc(i,k,j,iso4) = remlsc(i,k,j,iso4) - wetrem(iso4)/2.
 
!             chemical aqueous conversion diagnostic
              rxsaq1(i,k,j,iso2) = rxsaq1(i,k,j,iso2) - rxs1/2.
              rxsaq1(i,k,j,iso4) = rxsaq1(i,k,j,iso4) - rxs11/2.
 
            end do
          end do
 
!         cumulus clouds
!         wet removal by cumulus clouds (over the fraction of grid box
!         fracum) assume the cloud water content = 2 g/m3  (ref.
 
 
!         Kasibhatla )
          do i = 2 , ilxm
            if ( icumtop(i,j).ne.0 ) then
              do k = icumtop(i,j) , kx
                rxs2 = 0.0
                rxs21 = 0.0    ! fraction of conversion, not removed, as SO4 src
                wetrem_cvc(iso2) = 0.   ! scavenging for SO2, below lsc
                wetrem_cvc(iso4) = 0.
 
!               conversion from so2 to so4
                rxs2 = fracum(i,k)*chtrsol(iso2)*concmin(i,k)           &
                     & *(dexp(-2./360.*dt)-1.)
                rxs21 = rxs2*1.5
 
!               removal (including theremoval on the rxs21 term)
                wetrem_cvc(iso4) = (fracum(i,k)*chtrsol(iso4)*chib(i,k,j&
                                 & ,iso4)-rxs21)*(dexp(-remcum*dt)-1.)
 
!               tendancies due to convective processes
                chiten(i,k,j,iso2) = chiten(i,k,j,iso2) + rxs2/dt
                chiten(i,k,j,iso4) = chiten(i,k,j,iso4)                 &
                                   & + wetrem_cvc(iso4)/dt - rxs21/dt
 
!               diagnostic of wet deposition
!               remcvc(i,k,j,1) = remcvc(i,k,j,1) - wetrem_cvc(iso2)/2.
                remcvc(i,k,j,iso4) = remcvc(i,k,j,iso4)                 &
                                   & - wetrem_cvc(iso4)/2.
!               chemical aquesous conversion diagnostic
                rxsaq2(i,k,j,iso2) = rxsaq2(i,k,j,iso2) - rxs2/2.
                rxsaq2(i,k,j,iso4) = rxsaq2(i,k,j,iso4) - rxs21/2.
 
              end do
            end if
          end do
 
        end if ! end of SOX chemistry
 
!---------------------------------------
!       Other than sulfate CARBON AEROSOL, DUST
!----------------------------------------
        if ( chtrname(itr).eq.'BC_HB' .or. chtrname(itr).eq.'BC_HL' .or.&
           & chtrname(itr).eq.'OC_HB' .or. chtrname(itr).eq.'OC_HL' .or.&
           & chtrname(itr).eq.'DUST' ) then
 
!         wet deposition term
 
          if ( ichremlsc.eq.1 ) then
!           Wet removal at resolvable scale (fcc)
!           add non-precipitating cloud conversion (threshold
!           clmin=0.01g/m3) the same as that in subroutine exmois
!           clmin = 0.01
 
            do k = 1 , kx
              do i = 2 , ilxm
                if ( wl(i,k).gt.clmin ) then
                  wetrem(itr) = 0.
                  if ( remrat(i,k).gt.0. ) then
                    wetrem(itr) = fracloud(i,k)*chtrsol(itr)            &
                                & *chib(i,k,j,itr)                      &
                                & *(dexp(-remrat(i,k)*dt)-1.)
                    chiten(i,k,j,itr) = chiten(i,k,j,itr) + wetrem(itr) &
                                      & /dt
                    remlsc(i,k,j,itr) = remlsc(i,k,j,itr) - wetrem(itr) &
                                      & /2.
                  end if
                end if
              end do
            end do
          end if
 
          if ( ichremcvc.eq.1 ) then
!           sub-scale wet removal, cumulus cloud (fracum)
!           remcum = removal rate for cumulus cloud scavenging (s-1)
!           remcum = 1.e-3
            do i = 2 , ilxm
              if ( icumtop(i,j).ne.0 ) then
                do k = icumtop(i,j) , kx
                  wetrem_cvc(itr) = fracum(i,k)*chtrsol(itr)            &
                                  & *chib(i,k,j,itr)                    &
                                  & *(dexp(-remcum*dt)-1.)
                  chiten(i,k,j,itr) = chiten(i,k,j,itr)                 &
                                    & + wetrem_cvc(itr)/dt
                  remcvc(i,k,j,itr) = remcvc(i,k,j,itr)                 &
                                    & - wetrem_cvc(itr)/2.
                end do
              end if
            end do
          end if
 
        end if ! end wet removal DUST, CARBON
 
!       Conversion from hydrophobic to hydrophilic: Carbonaceopus
!       species time constant ( 1.15 day cooke et al.,1999)
 
        if ( gfcall3 .and. chtrname(itr).eq.'BC_HB' .and. ibchl.gt.0 )  &
           & then
          gfcall3 = .false.
          agct = 1.15*86400
!         agct = 2.30 *86400
!bxq      do itr=1,ntr
!bxq      agingtend(itr) = 0.
!bxq      end do
 
          do k = 1 , kx
            do i = 2 , ilxm
              agingtend(ibchb) = -chib(i,k,j,ibchb)*(1.-dexp(-dt/agct)) &
                               & /dt
              agingtend(ibchl) = -agingtend(ibchb)
 
              chiten(i,k,j,ibchb) = chiten(i,k,j,ibchb)                 &
                                  & + agingtend(ibchb)
              chiten(i,k,j,ibchl) = chiten(i,k,j,ibchl)                 &
                                  & + agingtend(ibchl)
            end do
          end do
        end if
 
 
        if ( gfcall4 .and. chtrname(itr).eq.'OC_HB' .and. iochl.gt.0 )  &
           & then
 
          gfcall4 = .false.
          agct = 1.15*86400
!         agct = 2.30 *86400
!bxq      do itr=1,ntr
!bxq      agingtend(itr) = 0.
!bxq      end do
 
          do k = 1 , kx
            do i = 2 , ilxm
              agingtend(iochb) = -chib(i,k,j,iochb)*(1-dexp(-dt/agct))  &
                               & /dt
              agingtend(iochl) = -agingtend(iochb)
 
              chiten(i,k,j,iochb) = chiten(i,k,j,iochb)                 &
                                  & + agingtend(iochb)
              chiten(i,k,j,iochl) = chiten(i,k,j,iochl)                 &
                                  & + agingtend(iochl)
            end do
          end do
 
        end if ! end aging
 
!*********************************
!       SURFACE source terms
!********************************
 
!       care chiten must be consistent with chia,b,c (= chi * pstar i.e)
!       1.e3 comes from Psurf/Pstar
 
!       en chantier DUST
!
!       print*,'before SFLUX',j,maxval(SFLT)
!       calculation of 10m wind ggaffe au facteur sur le vent
!       gaffe au kl ??
 
!       1 wind at 10 m
 
!       define the bin size
 
        if ( chtrname(itr).eq.'DUST' .and. gfcall2 ) then
 
          do i = 2 , ilxm
            ivegcov(i) = nint(veg2d(i,j))
            psurf(i) = psb(i,j)*1000. + ptop
 
            do k = 1 , kx
              ttb(i,k) = tb(i,k,j)/psb(i,j)
!             zza(i,k) = za(i,k,j)
            end do
 
!           calculate 10 M input for wind erosion and dry deposition
!           method based on bats diagnostic in routine interf.
 
            if ( (ivegcov(i).ne.0) ) then
              facv = dlog(za(i,kx,j)/10.)                               &
                   & /dlog(za(i,kx,j)/rough(ivegcov(i)))
              facb = dlog(za(i,kx,j)/10.)/dlog(za(i,kx,j)/zlnd)
              facs = dlog(za(i,kx,j)/10.)/dlog(za(i,kx,j)/zsno)
 
              fact = sfracv2d(i,j)*facv + sfracb2d(i,j)                 &
                   & *facb + sfracs2d(i,j)*facs
 
!             grid level effective roughness lenght ( a faire !)
              zeff(i) = rough(ivegcov(i)) ! ajouter contrib sol nu et snow
            else
!             water surface
              fact = dlog(za(i,kx,j)/10.)/dlog(za(i,kx,j)/zoce)
 
              zeff(i) = zoce
            end if
 
!           10 m wind
            u10 = (ub(i,kx,j)/psb(i,j))*(1-fact)
            v10 = (vb(i,kx,j)/psb(i,j))*(1-fact)
            wid10(i) = sqrt(u10**2+v10**2)
!           wid10(5) = 15
!           10 m air temperature
 
            temp10(i) = ttb(i,kx) - sdeltk2d(i,j)*fact
 
!           specific  humidity at 10m
            shu10 = (qvb(i,kx,j)/psb(i,j))/(1.+qvb(i,kx,j)/psb(i,j))    &
                  & - sdelqk2d(i,j)*fact
 
!           retransform in mixing ratio
 
            shu10 = shu10/(1-shu10)
 
!           saturation mixing ratio at 10m
            if ( temp10(i).gt.273.15 ) then
              satvp = svp1*1.E3*dexp(svp2*(temp10(i)-273.15)            &
                    & /(temp10(i)-svp3))
            else
              satvp = .611*1.E3*dexp(22.514-6.15E3/temp10(i))
            end if
            pres10 = psurf(i) - 98
            qsat10 = ep2*satvp/(pres10-satvp)
 
!           relative humidity at 10m
            rh10(i) = 0.
            if ( qsat10.gt.0. ) rh10(i) = shu10/qsat10
!
!           friction velocity ( not used fo the moment)
!
!           ustar(i) = sqrt ( uvdrag(i,j)                 *
!           &              sqrt ( (ub(i,kx,j)/psb(i,j) )**2 +
!           &                     (vb(i,kx,j)/psb(i,j) )**2 ) /
!           &                rho(i,kx)                          )
 
!           soil wetness
            soilw(i) = ssw2da(i,j)                                      &
                     & /(xmopor(iexsol(nint(satbrt(i,j))))*depuv        &
                     & (nint(satbrt(i,j))))
 
!           soilw(i) = ssw2da(i,j) /(xmopor(iexsol(ivegcov(i)) )
!           &                            * depuv(ivegcov(i))      )
 
!           fraction of vegetation
            vegfrac(i) = svegfrac2d(i,j)
 
!           soil texture ( from external preoproc)
            soilt(i) = nint(dustsotex(i,j))
 
!           surface temperature
!           over land recalculated from the BATS  deltk air/ surface
!           temperature account for a composite temperature between
 
!           bare ground and vegetation
            if ( ivegcov(i).ne.0 ) then
              tsurf(i) = ttb(i,kx) - sdeltk2d(i,j)
            else
!             ocean temperature in this case
              tsurf(i) = tgb(i,j)
            end if
 
!           aborbed solar radiation ( for stb criteria)
 
            srad(i) = sol2d(i,j)
 
          end do
 
          call sfflux(ix,2,ilxm,j,20,ivegcov(1),vegfrac(1),soilt(1),    &
                    & zeff,soilw(1),wid10(1),rho(1,kx),dustbsiz,rsfrow)
 
          call chdrydep(ix,2,ilxm,kx,1,nbin,ivegcov,ttb,rho,a,psurf,    &
                      & temp10,tsurf,srad,rh10,wid10,zeff,dustbsiz,     &
                      & pdepv)
 
 
          gfcall2 = .false.
        end if ! end dust flux and deposition calculated for bins
!       just calculated for the first case of itr = DUST (save time man)
 
        if ( chtrname(itr).eq.'DUST' ) then
!         define the corresponding index between itr and DUST bins
!         (can be different if we use DUST + other particles)
 
          ibin = ibin + 1
 
!         calculate the source tendancy
          do i = 2 , ilxm
            chemsrc(i,j,lmonth,itr) = rsfrow(i,ibin)
            chiten(i,kx,j,itr) = chiten(i,kx,j,itr) + rsfrow(i,ibin)    &
                               & *g/(dsigma(kx)*1.E3)
!           diagnostique source
            cemtr(i,j,itr) = cemtr(i,j,itr) + chemsrc(i,j,lmonth,itr)   &
                           & *dt/2.
          end do
 
!         calculate the tendancy du to gravitationnal settling and dry
 
!         deposition
          do k = 2 , kx
            do i = 2 , ilxm
              wk(i,k) = (1./psb(i,j))                                   &
                      & *(twt(k,1)*chib(i,k,j,itr)+twt(k,2)*chib(i,k-1, &
                      & j,itr))
            end do
          end do
 
!         remember PDEPV is defined for ibin which is not necessarly itr
 
          do i = 2 , ilxm
            do k = 2 , kx - 1
                        ! do not apply to the first level
              settend(i,k) = (wk(i,k+1)*pdepv(i,k+1,ibin)-wk(i,k)*pdepv(&
                           & i,k,ibin))*g*1.E-3/dsigma(k)
              chiten(i,k,j,itr) = chiten(i,k,j,itr) - settend(i,k)
            end do
!
            settend(i,kx) = -(wk(i,kx)*pdepv(i,kx,ibin)*g*1.E-3)        &
                          & /dsigma(kx)
            chiten(i,kx,j,itr) = chiten(i,kx,j,itr) + settend(i,kx)
 
!           dignoctic for dry deposition
            remdrd(i,j,itr) = remdrd(i,j,itr) - settend(i,kx)*dt/2.
          end do
 
        end if !( end calculation of dust tendancies)
 
!CCCC   Source tendenciesCCCC
 
        do i = 2 , ilxm
          if ( chtrname(itr).ne.'DUST' ) then
            chiten(i,kx,j,itr) = chiten(i,kx,j,itr)                     &
                               & + chemsrc(i,j,lmonth,itr)              &
                               & *g*0.7/(dsigma(kx)*1.E3)
            chiten(i,kx-1,j,itr) = chiten(i,kx-1,j,itr)                 &
                                 & + chemsrc(i,j,lmonth,itr)            &
                                 & *g*0.15/(dsigma(kx-1)*1.E3)
            chiten(i,kx-2,j,itr) = chiten(i,kx-2,j,itr)                 &
                                 & + chemsrc(i,j,lmonth,itr)            &
                                 & *g*0.15/(dsigma(kx-2)*1.E3)
!           diagnostic for source, cumul
            cemtr(i,j,itr) = cemtr(i,j,itr) + chemsrc(i,j,lmonth,itr)   &
                           & *dt/2.
          end if
        end do
 
!       end loop on tracers
!******************
      end do
!*******************
      end subroutine tractend2
