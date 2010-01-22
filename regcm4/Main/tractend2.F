C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tractend2(j)
c
c     This subroutine computes the tendencies for tracer transport and
c     chemistry
c     
c     ntr:           dimension of tracer arrays in species index
c     j:             index of j slice in current computation
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer j
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'mpiregcm.h'
#endif
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"

#     include "include/iunits.cb"
#     include "include/main.cb"
#     include "include/mainchem.cb"
#     include "include/bdycod.cb"
#     include "include/cvaria.cb"
#     include "include/trachem.cb"
#     include "include/pbldim.cb"
#     include "include/bats.cb2"
#     include "include/rad.cb"
#     include "include/pmoist.cb"
#     include "include/slice.cb"
#     include "include/dust.cb"
#     include "include/bats.cb"
#     include "include/date.cb"
c
      real(kind=8)  wetrem(ntr), wetrem_cvc(ntr)
      real(kind=8)  agingtend(ntr)
      real(kind=8)  agct
      real(kind=8)  fracum(ix,kx),fracloud(ix,kx)

      real(kind=8)  wl(ix,kx),rho(ix,kx)
      real(kind=8)  concmin(ix,kx)
      real(kind=8)  remcum,clmin, chimol
cbxq  real(kind=8)  h2o2mol
      real(kind=8)  rxs1, rxs11 ,rxs2, rxs21
      real(kind=8)  cldno,oh1,ak0tm,ak00t,akval
      real(kind=8)  cutend_up(ix,kx), cutend_dwd(ix,kx), wk(ix,kx)
      real(kind=8)  soilw(ix), wid10(ix),temp10(ix),rh10(ix),tsurf(ix)
      real(kind=8)  psurf(ix),srad(ix),vegfrac(ix),zeff(ix)
      real(kind=8)  ustar(ix)
      real(kind=8)  facv,facb,facs,fact, u10, v10
      real(kind=8)  satvp,shu10,pres10,qsat10 
      real(kind=8)  PDEPV(ix,kx,nbin), ttb(ix,kx),zza(ix,kx)
      real(kind=8)  settend(ix,kx)
      integer ivegcov(ix), soilt(ix)
      integer i,k,itr,kb,kdwd, ibin
      logical gfcall, gfcall2,gfcall3,gfcall4 
      real(kind=8)  RSFROW(ix,nbin),smask(ix)     
c
c clmin = non-precipitating cloud conversion threshold, clmin=0.01g/m3
      clmin = 0.01
c remcum= removal rate for cumulus cloud scavenging (s-1)
      remcum = 1.e-3
C
c Preliminary calculations independant of tracer nature
     
c   the unit: rho - kg/m3, wl - g/m3
      do 12 k=1,kx
      do 12 i=2,ilxm
        rho(i,k)=(psb(i,j)*a(k)+ptop)*1000./287./tb(i,k,j)*psb(i,j)
        wl(i,k)=qcb(i,k,j)/psb(i,j)*1000.*rho(i,k)
12    continue

      ivegcov = 0

c cloud fractionnal cover for wet deposition
c large scale : fracloud, calculated from fcc coming from pcp.f

c cumulus scale : fracum, calculated from the total cloud fraction (as defined for the radiation scheme in cldfrac.f routine)

c
      do  i=2,ilxm
       do  k=1,kx
        fracloud(i,k)= dmin1( fcc(i,k,j), fcmax)
        fracum(i,k) = 0.  
       end do
       if (icumtop(i,j).ne.0) then
        do k=icumtop(i,j),kx
          fracum(i,k)= cldfra(i,k) - fracloud(i,k)         
        end do
       endif
      end do


c TRANSPORT OF TRACERS
c----initialize tracer tendencies, and scratch arrays
      do 10 itr=1,ntr
      do 10 k=1,kx
      do 10 i=1,ilx
        chiten(i,k,j,itr)= 0.
10    continue

c-----horizontal and vertical advection

      do itr=1,ntr     
c   
       call hadvCH(chiten(1,1,j,itr),dx,itr,j,2)

       call vadv(chiten(1,1,j,itr),chia(1,1,j,itr),j,5)

c      call vadv(chiten(1,1,j,itr),chia(1,1,j,itr),j,3)

c----horizontal diffusion: initialize scratch vars to 0.
c  need to compute tracer tendencies due to diffusion
       call diffutCH(chiten(1,1,j,itr),xkc(1,1,j),c203,itr,j)
 
      end do ! end tracer loop

c subgrid vertical transport by convective mass flux : a modifier !
     
      if(ichcumtra .eq. 2) then
        do 210  k=2,kx
        do 210  i=2,ilxm
          wk(i,k)= (1./psa(i,j))*(twt(k,1)*chib(i,k,j,itr)
     &           + twt(k,2)*chib(i,k-1,j,itr))

          cutend_up(i,k)=0.
          cutend_dwd(i,k)=0.
210     continue        

        do  i=2,ilxm

          if (icumtop(i,j).ne.0) then

            kt=max0(icumtop(i,j),3)           
            kb=icumbot(i,j) 
            kdwd = icumdwd(i,j) 

c          cutend(i,kt) =  mflx(i) * g * 1.e-3*
c     &               (wk(i,kb)-wk(i,kt))/(sigma(kb)-sigma(kt))

c transport linked to updraft 
c betwwen kt et kdwd , the tendancy is averaged (mixing)

            if (kdwd .lt. kt) print*,'PBLEME !'
            do k=kt, kdwd
              cutend_up(i,k) =  mflx(i,1) * g * 1.e-3*
     &               wk(i,kb)/(sigma(kdwd) - sigma(kt))
            end do

            cutend_up(i,kb) =  -  mflx(i,1) * g * 1.e-3*
     &               wk(i,kb)/(dsigma(kb))
c transport linked to downdraft 
             
            cutend_dwd(i,kdwd) =  -  mflx(i,2) * g * 1.e-3*
     &               wk(i,kdwd)/(dsigma(kdwd))

            cutend_dwd(i,kx) =  +  mflx(i,2) * g * 1.e-3*
     &               wk(i,kdwd)/(dsigma(kx))

            do k=kt,kx           
              chiten(i,k,j,itr) = chiten(i,k,j,itr)+cutend_up(i,k) 
     &                        + cutend_dwd(i,k)
            end do
          end if
        end do
      end if

c
c SOURCE AND SINKS TERMS ( dependant on the nature of tracers) 
      gfcall =  .true. ! logical call for SOx 
      gfcall2 = .true. ! logical for DUST  
      gfcall3 = .true. ! logical for BC aging
      gfcall4 = .true. ! logical for OC aging

      ibin = 0
c*******************************************
c begining of tracer loop 
c*******************************************
      DO 99 itr = 1,ntr
c****************************
c       print*,'gfcall', itr, gfcall

c
c---------------------------------------------
c SOX CHEMSITRY / IN-CLOUD PROCESS 
c--------------------------------------------
        IF(gfcall .and. ( chtrname(itr).eq.'SO2' .and. iso4.gt.0))then
          gfcall = .false. 

c
c  GAZEOUS CONVERSION   

c calculate the k value for the gaseous conversion by OH

c ohconc

          do 50 k = 1,kx
          do 50 i = 2,ilxm
            cldno=1.        ! no cloud fraction

c        if(coszrs(i,j).lt.0.001) ohconc(i,j,k)=ohconc(i,j,k)*0.01
c        oh1=ohconc(i,j,k)*rho(i,k)*2.084e13             ! molecules/cm3
c test j eprends directement une valeur de oh1

            oh1= 15.e5                                   ! molecules/cm3 
            if(coszrs(i).lt.0.001) oh1=oh1*0.01          ! diurnal evolution

            ak0tm=3.e-31*(ta(i,k,j)/psa(i,j)/300.)**(-3.3)
     &           *rho(i,k)*2.084e19                      ! k0(T)*[M]
            ak00t=1.5e-12                                ! K00(T)
            akval=ak0tm/(1.+ak0tm/ak00t)
     &           *0.6**(1./(1.+(dlog10(ak0tm/ak00t))**2.))

c tendencies 
c here p1 unit: s^-1  and the ratio of molar mass of SO4 to SO2 is 96/64 = 1.5

            chiten(i,k,j,iso2) =chiten(i,k,j,iso2) 
     &                       -chib(i,k,j,iso2)*akval*oh1*cldno
            chiten(i,k,j,iso4) =chiten(i,k,j,iso4) 
     &                       +chib(i,k,j,iso2)*akval*oh1*cldno*1.5

c gazeous conversion diagnostic 

            rxsg(i,k,j,iso2) =rxsg(i,k,j,iso2)
     &                       +chib(i,k,j,iso2)*akval*oh1*cldno*dt/2.
            rxsg(i,k,j,iso4) =rxsg(i,k,j,iso4)
     &                       +chib(i,k,j,iso2)*akval*oh1*cldno*1.5*dt/2.

 50       continue

c AQUEOUS CONVERSION IN CLOUDS AND WET REMOVAL

c Aqueous conversion from so2 to so4 ;control by h2o2  

          do k = 1,kx
          do i = 2,ilxm
cbxq        h2o2mol = 1.e-6 * h2o2conc(i,j,k)
            chimol  = 28.9/64.*chib(i,k,j,iso2)/psb(i,j)       ! kg/kg to mole
c concmin(i,k)=dmin1(h2o2mol,chimol)*64./28.9*psb(i,j)  ! cb*kg/kg
c do tests, suppose h2o2 always enough
            concmin(i,k) = chimol * 64./28.9*psb(i,j)       ! cb*kg/kg
          enddo
          enddo
 
c Large scale clouds     

          do 250 k=1,kx
          do 250 i=2,ilxm
            rxs1 = 0.0
            rxs11= 0.0         ! fraction of conversion, not removed, as SO4 src
            wetrem(iso2)=0.    ! scavenging for SO2, below lsc
            wetrem(iso4)=0.        

            if(wl(i,k).gt.clmin) then
c conversion from so2 to so4
              rxs1 = fracloud(i,k)*chtrsol(iso2)*concmin(i,k)
     &             *(dexp(-wl(i,k)/360.*dt) - 1.)         

              rxs11 = rxs1 *1.5 ! SO4 src term and the ratio of molar
                                ! mass of SO4 to SO2 is 96/64 = 1.5

c if removal occurs, a fraction of SO4 src term is also removed
c and accounted for in the term  wetrem(iso4)
              if(remrat(i,k).gt.0.) then

                wetrem(iso4) =
     &             (fracloud(i,k)*chtrsol(iso4)*chib(i,k,j,iso4)-rxs11)
     &                       * (dexp(-remrat(i,k)*dt) - 1.)
              endif
            endif

c Below cloud scavenging only for SO2

            if(rembc(i,k).gt.0.) then
              wetrem(iso2) = fracloud(i,k)*chtrsol(iso2)*concmin(i,k)
     &                     *(dexp(-rembc(i,k)*dt) - 1.)
            endif

c Tendancies large scale cloud 
            chiten(i,k,j,iso2)=chiten(i,k,j,iso2)
     &                        +rxs1/dt +wetrem(iso2)/dt
            chiten(i,k,j,iso4)=chiten(i,k,j,iso4)
     &                        -rxs11/dt+wetrem(iso4)/dt

c and wetdep diagnostics
            remlsc(i,k,j,iso2) = remlsc(i,k,j,iso2) - wetrem(iso2)/2.
            remlsc(i,k,j,iso4) = remlsc(i,k,j,iso4) - wetrem(iso4)/2.

c chemical aqueous conversion diagnostic
            rxsaq1(i,k,j,iso2) = rxsaq1(i,k,j,iso2) - rxs1/2.
            rxsaq1(i,k,j,iso4) = rxsaq1(i,k,j,iso4) - rxs11/2.

 250      continue   

c cumulus clouds 
c  wet removal by cumulus clouds (over the fraction of grid box fracum)
c assume the cloud water content = 2 g/m3  (ref. Kasibhatla )


          do 300 i=2,ilxm
            if (icumtop(i,j).ne.0) then
              do 290 k=icumtop(i,j),kx
                rxs2 = 0.0
                rxs21= 0.0     ! fraction of conversion, not removed, as SO4 src
                wetrem_cvc(iso2)=0.     ! scavenging for SO2, below lsc
                wetrem_cvc(iso4)=0. 

c conversion from so2 to so4
                rxs2 = fracum(i,k)*chtrsol(iso2) * concmin(i,k) *
     &                  (dexp(-2./360.*dt) - 1.)
                rxs21 = rxs2 *1.5

c removal (including theremoval on the rxs21 term) 
                wetrem_cvc(iso4)=
     &             (fracum(i,k)*chtrsol(iso4)*chib(i,k,j,iso4)-rxs21)
     &                          * (dexp(-remcum*dt) - 1.)

c tendancies due to convective processes
                chiten(i,k,j,iso2) = chiten(i,k,j,iso2) +rxs2/dt
                chiten(i,k,j,iso4) = chiten(i,k,j,iso4)
     &                             +wetrem_cvc(iso4)/dt -rxs21/dt

c diagnostic of wet deposition
c              remcvc(i,k,j,1) = remcvc(i,k,j,1) - wetrem_cvc(iso2)/2.
               remcvc(i,k,j,iso4)=remcvc(i,k,j,iso4)-wetrem_cvc(iso4)/2.
c chemical aquesous conversion diagnostic
                rxsaq2(i,k,j,iso2) = rxsaq2(i,k,j,iso2) - rxs2/2.
                rxsaq2(i,k,j,iso4) = rxsaq2(i,k,j,iso4) - rxs21/2.

  290         continue
            endif
  300     continue

        END IF ! end of SOX chemistry

c---------------------------------------
c Other than sulfate CARBON AEROSOL, DUST  
c----------------------------------------
        IF (chtrname(itr) .eq.'BC_HB' .or.  chtrname(itr) .eq. 'BC_HL'
     & .or. chtrname(itr) .eq.'OC_HB' .or.  chtrname(itr) .eq. 'OC_HL'
     & .or. chtrname(itr) .eq.'DUST') then

c wet deposition term 

        if (ichremlsc .eq. 1) then 
c Wet removal at resolvable scale (fcc)
c add non-precipitating cloud conversion (threshold clmin=0.01g/m3)
c the same as that in subroutine exmois
c         clmin = 0.01

          do k=1,kx
          do i=2,ilxm
            if(wl(i,k).gt.clmin) then
              wetrem(itr)=0. 
              if(remrat(i,k).gt.0.) then 
                wetrem(itr)= fracloud(i,k)*chtrsol(itr)*chib(i,k,j,itr)
     &                     *(dexp(-remrat(i,k)*dt) - 1.) 
                chiten(i,k,j,itr) = chiten(i,k,j,itr) + wetrem(itr)/dt
                remlsc(i,k,j,itr) = remlsc(i,k,j,itr) - wetrem(itr)/2. 
              end if         
            end if
          end do
          end do
        end if 

        if ( ichremcvc .eq. 1) then  
c sub-scale wet removal, cumulus cloud (fracum)
c remcum = removal rate for cumulus cloud scavenging (s-1)
c         remcum = 1.e-3
          do i=2,ilxm
            if (icumtop(i,j).ne.0) then
              do k=icumtop(i,j),kx
                wetrem_cvc(itr)=fracum(i,k)*chtrsol(itr)*chib(i,k,j,itr)
     &                         *(dexp(-remcum*dt)- 1.)    
                chiten(i,k,j,itr)= chiten(i,k,j,itr)+ wetrem_cvc(itr)/dt
                remcvc(i,k,j,itr)= remcvc(i,k,j,itr)- wetrem_cvc(itr)/2.
              end do
            end if
          end do
        end if 

        END IF ! end wet removal DUST, CARBON

c Conversion from hydrophobic to hydrophilic: Carbonaceopus species
c time constant ( 1.15 day cooke et al.,1999)       

        if(gfcall3 .and. chtrname(itr).eq.'BC_HB' .and. ibchl.gt.0 )then
          gfcall3 = .false.
          agct =    1.15 *86400
c         agct = 2.30 *86400
Cbxq      do itr=1,ntr
Cbxq        agingtend(itr) = 0.
Cbxq      enddo

          do k=1,kx  
          do i=2,ilxm
            agingtend(ibchb) = -chib(i,k,j,ibchb)*(1.-dexp(-dt/agct))/dt
            agingtend(ibchl) = -agingtend(ibchb)
                
            chiten(i,k,j,ibchb) = chiten(i,k,j,ibchb) + agingtend(ibchb)
            chiten(i,k,j,ibchl) = chiten(i,k,j,ibchl) + agingtend(ibchl)
          end do
          end do
        end if 
 

        if(gfcall4 .and. chtrname(itr).eq.'OC_HB' .and. iochl.gt.0 )then

          gfcall4 = .false.
          agct =    1.15 *86400
c         agct = 2.30 *86400
Cbxq      do itr=1,ntr
Cbxq        agingtend(itr) = 0.
Cbxq      enddo
 
          do k=1,kx  
          do i=2,ilxm
            agingtend(iochb) = -chib(i,k,j,iochb)*(1-dexp(-dt/agct))/dt
            agingtend(iochl) = -agingtend(iochb)

            chiten(i,k,j,iochb) = chiten(i,k,j,iochb) + agingtend(iochb)
            chiten(i,k,j,iochl) = chiten(i,k,j,iochl) + agingtend(iochl)
          end do   
          end do

        end if ! end aging

c*********************************
c SURFACE source terms
c********************************

c care chiten must be consistent with chia,b,c (= chi * pstar i.e)
c 1.e3 comes from Psurf/Pstar 

c en chantier DUST
c 
c      print*,'before SFLUX',j,maxval(SFLT)
c calculation of 10m wind ggaffe au facteur sur le vent        
c gaffe au kl ?? 

c 1 wind at 10 m 

c define the bin size   

        IF ( chtrname(itr) .eq.'DUST' .and. gfcall2 ) then
                        
          do i= 2,ilxm
            ivegcov(i) = nint(veg2d(i,j))
            psurf(i)  =  psb(i,j)*1000.+ptop
  
            do k=1,kx
              ttb(i,k) = tb(i,k,j)/ psb(i,j)
              zza(i,k) = za(i,k,j) 
            end do 

c        calculate 10 M input for wind erosion and dry deposition
c        method based on bats diagnostic in routine interf.

            if ( (ivegcov(i) .ne. 0)) then 
            facv=dlog(za(i,kx,j)/10.)/dlog(za(i,kx,j)/rough(ivegcov(i)))
            facb=dlog(za(i,kx,j)/10.)/dlog(za(i,kx,j)/zlnd)
            facs=dlog(za(i,kx,j)/10.)/dlog(za(i,kx,j)/zsno)

              fact = sfracv2d(i,j)*facv +sfracb2d(i,j)*facb
     &              +sfracs2d(i,j)*facs

c grid level effective roughness lenght ( a faire !) 
              zeff(i) = rough(ivegcov(i)) ! ajouter contrib sol nu et snow
            else
c    water surface 
              fact= dlog(za(i,kx,j)/10.)/dlog(za(i,kx,j)/zoce)

              zeff(i)  = zoce
            end if

c 10 m wind 
            u10 = ( ub(i,kx,j)/psb(i,j)) * (1 -fact) 
            v10 = ( vb(i,kx,j)/psb(i,j)) * (1-fact)
            wid10(i) = sqrt ( u10**2 + v10**2)
c           wid10(5) = 15
c 10 m air temperature        

            temp10(i) = ttb(i,kx) - sdeltk2d(i,j) * fact

c specific  humidity at 10m
            shu10  = (qvb(i,kx,j)/psb(i,j))/ (1.+qvb(i,kx,j)/psb(i,j)) 
     &             -   sdelqk2d(i,j) * fact

c retransform in mixing ratio

            shu10  = shu10 / (1-shu10)  

c saturation mixing ratio at 10m
            if(temp10(i).gt.273.15)then
              satvp=svp1*1.e3*dexp(svp2*(temp10(i)-273.15)
     &                         /(temp10(i)-svp3))
            else
              satvp=.611*1.e3*dexp(22.514-6.15e3/temp10(i))
            endif
            pres10 = psurf(i) - 98
            qsat10  =ep2*satvp/(pres10-satvp)

c relative humidity at 10m
            rh10(i)=0.
            if (qsat10.gt.0.) rh10 (i) = shu10/qsat10
c
c friction velocity ( not used fo the moment)

            ustar(i) = sqrt ( uvdrag(i,j)                 * 
     &              sqrt ( (ub(i,kx,j)/psb(i,j) )**2 + 
     &                     (vb(i,kx,j)/psb(i,j) )**2 ) /
     &                rho(i,kx)                          )

c soil wetness
            soilw(i) = ssw2da(i,j) /(xmopor(iexsol(nint(satbrt(i,j))))
     &                            * depuv(nint(satbrt(i,j))))

c           soilw(i) = ssw2da(i,j) /(xmopor(iexsol(ivegcov(i)) ) 
c     &                            * depuv(ivegcov(i))      )

c fraction of vegetation
            vegfrac(i) = svegfrac2d(i,j) 

c soil texture ( from external preoproc) 
            soilt(i) = nint( dustsotex(i,j)) 

c surface temperature
c over land recalculated from the BATS  deltk air/ surface temperature
c account for a composite temperature between bare ground and vegetation

            if ( ivegcov(i) .ne. 0) then 
              tsurf(i) = ttb(i,kx) -  sdeltk2d(i,j)
            else
c ocean temperature in this case
              tsurf(i) = tgb(i,j) 
            end if

c aborbed solar radiation ( for stb criteria) 

            srad(i) = sol2d(i,j)

          end do

          call SFFLUX( ix, 2, ilxm,j,20,ivegcov(1),vegfrac(1),  
     &              soilt(1),smask(1),zeff,soilw(1), wid10(1),rho(1,kx),
     &              dustbsiz,RSFROW )

          call chdrydep( ix, 2,ilxm,kx,1,nbin,
     &                  ivegcov, ttb, rho, a, psurf,
     &                temp10, tsurf,srad,rh10, wid10,zeff,dustbsiz,
     &                PDEPV)  

         
          gfcall2 = .false.
        END IF ! end dust flux and deposition calculated for bins
c just calculated for the first case of itr = DUST (save time man)

        IF ( chtrname (itr) .eq. 'DUST') then 
c define the corresponding index between itr and DUST bins 
c (can be different if we use DUST + other particles) 
 
          ibin = ibin+1

c calculate the source tendancy 
          do i=2,ilxm
            chemsrc(i,j,lmonth,itr) =  RSFROW(i,ibin) 
            chiten(i,kx,j,itr)=chiten(i,kx,j,itr) + RSFROW(i,ibin) 
     &                          * g / (dsigma(kx)*1.e3)
c diagnostique source
            cemtr(i,j,itr)=cemtr(i,j,itr)+chemsrc(i,j,lmonth,itr)*dt/2.
          end do

c  calculate the tendancy du to gravitationnal settling and dry deposition 

          do k=2,kx
          do i=2,ilxm
            wk(i,k)= (1./psb(i,j))*(twt(k,1)*chib(i,k,j,itr)
     &            + twt(k,2)*chib(i,k-1,j,itr))
          end do
          end do       

c  remember PDEPV is defined for ibin which is not necessarly itr

          do i=2,ilxm 
            do k=2,kx-1 ! do not apply to the first level 
              settend(i,k) = ( wk(i,k+1)*PDEPV(i,k+1,ibin) 
     &                        -wk(i,k)  *PDEPV(i,k,ibin)  )
     &                   * g * 1.e-3/dsigma(k)       
              chiten(i,k,j,itr)=chiten(i,k,j,itr) - settend(i,k)
            end do
c    
            settend(i,kx) = - (wk(i,kx) *PDEPV(i,kx,ibin)*g* 1.e-3)/
     &                   dsigma(kx)       
            chiten(i,kx,j,itr)=chiten(i,kx,j,itr) + settend(i,kx)
        
c dignoctic for dry deposition 
            remdrd(i,j,itr)=remdrd(i,j,itr) - settend(i,kx)* dt/2.
          end do 

        END IF !( end calculation of dust tendancies)      

CCCCC Source tendenciesCCCC
     
        do i=2,ilxm
          if (chtrname(itr) .ne. 'DUST'   ) then  
            chiten(i,kx,j,itr)=chiten(i,kx,j,itr)
     &         + chemsrc(i,j,lmonth,itr) * g *0.7/ (dsigma(kx)*1.e3)
            chiten(i,kx-1,j,itr)=chiten(i,kx-1,j,itr)
     &         + chemsrc(i,j,lmonth,itr) * g *0.15/ (dsigma(kx-1)*1.e3)
            chiten(i,kx-2,j,itr)=chiten(i,kx-2,j,itr)
     &         + chemsrc(i,j,lmonth,itr) * g *0.15/ (dsigma(kx-2)*1.e3)
c      diagnostic for source, cumul
            cemtr(i,j,itr)=cemtr(i,j,itr)+chemsrc(i,j,lmonth,itr)*dt/2.
          end if 
        end do

c end loop on tracers
c******************
99    continue    
c*******************
      return
      end
