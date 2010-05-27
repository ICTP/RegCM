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
 
      subroutine radclw(jslc,ts,tnm,qnm,o3vmr,pmid,pint,pmln,piln,plco2,&
                      & plh2o,n2o,ch4,cfc11,cfc12,cld,tclrsf,qrl,flns,  &
                      & flnt,flnsc,flntc,flwds,emiss1d,aerlwfo,aerlwfos)

!-----------------------------------------------------------------------
!
! Compute longwave radiation heating rates and boundary fluxes
!
! Uses broad band absorptivity/emissivity method to compute clear sky;
! assumes randomly overlapped clouds with variable cloud emissivity to
! include effects of clouds.
!
! Computes clear sky absorptivity/emissivity at lower frequency (in
! general) than the model radiation frequency; uses previously computed
! and stored values for efficiency
!
! Note: This subroutine contains vertical indexing which proceeds
!       from bottom to top rather than the top to bottom indexing
!       used in the rest of the model.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
      use mod_dynparam
      use mod_param1 , only : ifrabe
      use mod_param2 , only : iemiss , idirect
      use mod_radbuf
      use mod_aerosol, only : aerlwtr
      use mod_date , only : jyear , jyear0 , ktau
      use mod_constants , only : gocp , stebol
      use mod_when , only : whenne , whenflt
      implicit none
!
!     Input arguments
!
! ts      - Ground (skin) temperature
! emiss1d - Emissivity of surface
!
!     Input arguments which are only passed to other routines
!
! tnm     - Level temperature
! qnm     - Level moisture field
! o3vmr   - ozone volume mixing ratio
! pmid    - Level pressure
! pint    - Model interface pressure
! pmln    - Ln(pmid)
! piln    - Ln(pint)
! plco2   - Path length co2
! plh2o   - Path length h2o
! n2o     - nitrous oxide mass mixing ratio
! ch4     - methane mass mixing ratio
! cfc11   - cfc11 mass mixing ratio
! cfc12   - cfc12 mass mixing ratio
!
!     Input/Output arguments
!
! cld     - Cloud cover
! tclrsf  - Clear sky fraction
!
!     Output arguments
!
! qrl     - Longwave heating rate
! flns    - Surface cooling flux
! flnt    - Net outgoing flux
! flnsc   - Clear sky surface cooing
! flntc   - Net clear sky outgoing flux
! flwds   - Down longwave flux at surface
!
! Dummy arguments
!
      integer :: jslc
      real(8) , dimension(iym1,kz) :: cfc11 , cfc12 , ch4 , n2o ,    &
           & o3vmr , pmid , pmln , qnm , qrl , tnm
      real(8) , dimension(iym1,kzp1) :: cld , piln , pint , plco2 ,   &
           & plh2o , tclrsf
      real(8) , dimension(iym1) :: emiss1d , flns , flnsc , flnt ,     &
                                  & flntc , flwds , ts
      real(8), dimension(iym1) :: aerlwfo , aerlwfos


      intent (in) cld , emiss1d
      intent (out) flns , flnsc , flnt , flntc , flwds , qrl , aerlwfo, &
                 & aerlwfos

      intent (inout) tclrsf
!
!---------------------------Local variables-----------------------------
!
! i       - Longitude index
! k       - Level index
! k1      - Level index
! k2      - Level index
! k3      - Level index
! km      - Level index
! km1     - Level index
! km2     - Level index
! km3     - Level index
! km4     - Level index
! tmp     - Temporary
! tmp1    - Temporary 1
! absbt   - Downward emission at model top
! plol    - O3 pressure wghted path length
! plos    - O3 path length
! co2em   - Layer co2 normalized planck funct. derivative
! co2eml  - Interface co2 normalized planck funct. deriv.
! delt    - Diff t**4 mid layer to top interface
! delt1   - Diff t**4 lower intrfc to mid layer
! bk1     - Absrptvty for vertical quadrature
! bk2     - Absrptvty for vertical quadrature
! ful     - Total upwards longwave flux
! fsul    - Clear sky upwards longwave flux
! fdl     - Total downwards longwave flux
! fsdl    - Clear sky downwards longwv flux
! fclb4   - Sig t**4 for cld bottom interfc
! fclt4   - Sig t**4 for cloud top interfc
! s       - Flx integral sum
! tplnka  - Planck fnctn temperature
! s2c     - H2o cont amount
! s2t     - H2o cont temperature
! w       - H2o path
! tplnke  - Planck fnctn temperature
! h2otr   - H2o trnmsn for o3 overlap
! co2t    - Prs wghted temperature path
! tint    - Interface temperature
! tint4   - Interface temperature**4
! tlayr   - Level temperature
! tlayr4  - Level temperature**4
! rtclrsf - 1./tclrsf(i,k)
! klov    - Cloud lowest level index
! khiv    - Cloud highest level index
! khivm   - khiv(i) - 1
!
!     Trace gas variables
!
! ucfc11  - CFC11 path length
! ucfc12  - CFC12 path length
! un2o0   - N2O path length
! un2o1   - N2O path length (hot band)
! uch4    - CH4 path length
! uco211  - CO2 9.4 micron band path length
! uco212  - CO2 9.4 micron band path length
! uco213  - CO2 9.4 micron band path length
! uco221  - CO2 10.4 micron band path length
! uco222  - CO2 10.4 micron band path length
! uco223  - CO2 10.4 micron band path length
! bn2o0   - pressure factor for n2o
! bn2o1   - pressure factor for n2o
! bch4    - pressure factor for ch4
! uptype  - p-type continuum path length
! abplnk1 - non-nearest layer Plack factor
! abplnk2 - nearest layer factor
!
! Local variables
!
      real(8) , dimension(14,iym1,kzp1) :: abplnk1 , abplnk2
      real(8) , dimension(iym1) :: absbt , bk1 , bk2 , delt , delt1 ,  &
                                  & tmp , tplnke
      real(8) , dimension(iym1,kzp1) :: bch4 , bn2o0 , bn2o1 , co2em ,&
           & co2t , fdl , fsdl , fsul , ful , h2otr , plol , plos ,     &
           & rtclrsf , s2c , s2t , tint , tint4 , tlayr , tlayr4 ,      &
           & tplnka , ucfc11 , ucfc12 , uch4 , uco211 , uco212 ,        &
           & uco213 , uco221 , uco222 , uco223 , un2o0 , un2o1 ,        &
           & uptype , w , fsul0 , fsdl0
      real(8) , dimension(iym1,kz) :: co2eml , fclb4 , fclt4
      logical , dimension(iym1) :: done , start
      real(8) :: tmp1
      integer :: i , ii , k , k1 , k2 , k3 , khighest , km , km1 , km2 ,&
               & km3 , km4 , iym1c , rad , n
      integer , dimension(iym1) :: indx , khiv , khivm , klov
      real(8) , dimension(iym1,kzp1,kzp1) :: s
      real(8) , dimension(iym1,kzp1,kzp1) :: tone
!
      integer , external :: intmax

      tone(:,:,:)=1. 
!
      do i = 1 , iym1
        rtclrsf(i,1) = 1.0/tclrsf(i,1)
      end do
!
      do k = 1 , kz
        do i = 1 , iym1
          fclb4(i,k) = 0.
          fclt4(i,k) = 0.
          tclrsf(i,k+1) = tclrsf(i,k)*(1.-cld(i,k+1))
          rtclrsf(i,k+1) = 1./tclrsf(i,k+1)
        end do
      end do
!
!     Calculate some temperatures needed to derive absorptivity and
!     emissivity, as well as some h2o path lengths
!
      call radtpl(tnm,ts,qnm,pint,plh2o,tplnka,s2c,s2t,w,tplnke,tint,   &
                & tint4,tlayr,tlayr4,pmln,piln)
 
!     do emissivity and absorptivity calculations
!     only if abs/ems computation
!
      if ( (jyear.eq.jyear0 .and. ktau.eq.0) .or.                       &
         & (mod(ktau+1,ifrabe).eq.0) ) then
 
!
!       Compute ozone path lengths at frequency of a/e calculation.
!
        call radoz2(o3vmr,pint,plol,plos)
!
!       Compute trace gas path lengths
!
        call trcpth(tnm,pint,cfc11,cfc12,n2o,ch4,qnm,ucfc11,ucfc12,     &
                  & un2o0,un2o1,uch4,uco211,uco212,uco213,uco221,uco222,&
                  & uco223,bn2o0,bn2o1,bch4,uptype)
!
!
!       Compute total emissivity:
!
        call radems(s2c,s2t,w,tplnke,plh2o,pint,plco2,tint,tint4,tlayr, &
                  & tlayr4,plol,plos,ucfc11,ucfc12,un2o0,un2o1,uch4,    &
                  & uco211,uco212,uco213,uco221,uco222,uco223,uptype,   &
                  & bn2o0,bn2o1,bch4,co2em,co2eml,co2t,h2otr,abplnk1,   &
                  & abplnk2,jslc)
!
!       Compute total absorptivity:
!
        call radabs(pmid,pint,co2em,co2eml,tplnka,s2c,s2t,w,h2otr,plco2,&
                  & plh2o,co2t,tint,tlayr,plol,plos,pmln,piln,ucfc11,   &
                  & ucfc12,un2o0,un2o1,uch4,uco211,uco212,uco213,uco221,&
                  & uco222,uco223,uptype,bn2o0,bn2o1,bch4,abplnk1,      &
                  & abplnk2,jslc)

         abstot0(:,:,:,jslc) = abstot(:,:,:,jslc)
         emstot0(:,:,jslc) = emstot(:,:,jslc)
         absnxt0(:,:,:,jslc) = absnxt(:,:,:,jslc)   
 
      end if
!
!     Find the lowest and highest level cloud for each grid point
!     Note: Vertical indexing here proceeds from bottom to top
!
      do i = 1 , iym1
        klov(i) = 0
        done(i) = .false.
      end do
      do k = 1 , kz
        do i = 1 , iym1
          if ( .not.done(i) .and. cld(i,kzp2-k).gt.0.0 ) then
            done(i) = .true.
            klov(i) = k
          end if
        end do
      end do
      call whenne(iym1,klov,1,0,indx,iym1c)
      do i = 1 , iym1
        khiv(i) = klov(i)
        done(i) = .false.
      end do
      do k = kz , 1 , -1
        do ii = 1 , iym1c
          i = indx(ii)
          if ( .not.done(i) .and. cld(i,kzp2-k).gt.0.0 ) then
            done(i) = .true.
            khiv(i) = k
          end if
        end do
      end do
      do i = 1 , iym1
        khivm(i) = khiv(i) - 1
      end do
!
!     Note: Vertical indexing here proceeds from bottom to top
!
      do ii = 1 , iym1c
        i = indx(ii)
        do k = klov(i) , khiv(i)
          fclt4(i,kzp1-k) = stebol*tint4(i,kzp2-k)
          fclb4(i,kzp1-k) = stebol*tint4(i,kzp3-k)
        end do
      end do

!
! option to calculate LW aerosol radiative forcing

!     CCFAB LW radiative forcing ( rad=1 : avec dust)
      fsul0(:,:) = 0.
      fsdl0(:,:) = 0.
      abstot(:,:,:,jslc) = abstot0(:,:,:,jslc)
      emstot(:,:,jslc) = emstot0(:,:,jslc)
      absnxt(:,:,:,jslc) = absnxt0(:,:,:,jslc)

      do rad = 1 , 2
        if ( rad==2 ) then

          abstot(:,:,:,jslc) = 1-(1-abstot0(:,:,:,jslc))*aerlwtr(:,:,:)

          emstot(:,:,jslc) = 1-(1-emstot0(:,:,jslc))*aerlwtr(:,:,1)
       
          do k = 1 , kz  ! aerlwtr defined on plev levels
            do n = 1 , 4
              absnxt(:,k,n,jslc) = 1-(1-absnxt0(:,k,n,jslc))*           &
                                & (aerlwtr(:,k,k+1)**xuinpl(:,k,n,jslc))
!             print* , 'z' , jslc , xuinpl(:,k,n,jslc)
            end do
          end do
      
!         print * , jslc, minval(absnxt(:,:,:,jslc)),                   &
!        &   maxval(absnxt(:,:,:,jslc)), minval(absnxt0(:,:,:,jslc)),   &
!        &   maxval(absnxt0(:,:,:,jslc))

        end if
!
!     Compute sums used in integrals (all longitude points)
!
!     Definition of bk1 & bk2 depends on finite differencing.  for
!     trapezoidal rule bk1=bk2. trapezoidal rule applied for nonadjacent
!     layers only.
!
!     delt=t**4 in layer above current sigma level km.
!     delt1=t**4 in layer below current sigma level km.
!
        do i = 1 , iym1
          delt(i) = tint4(i,kz) - tlayr4(i,kzp1)
          delt1(i) = tlayr4(i,kzp1) - tint4(i,kzp1)
          s(i,kzp1,kzp1) = stebol*(delt1(i)*absnxt(i,kz,1,jslc)         &
                         & +delt(i)*absnxt(i,kz,4,jslc))
          s(i,kz,kzp1) = stebol*(delt(i)*absnxt(i,kz,2,jslc)+delt1(i)   &
                        & *absnxt(i,kz,3,jslc))
        end do
        do k = 1 , kz - 1
          do i = 1 , iym1
            bk2(i) = (abstot(i,k,kz,jslc)+abstot(i,k,kzp1,jslc))*0.5
            bk1(i) = bk2(i)
            s(i,k,kzp1) = stebol*(bk2(i)*delt(i)+bk1(i)*delt1(i))
          end do
        end do
!
!       All k, km>1
!
        do km = kz , 2 , -1
          do i = 1 , iym1
            delt(i) = tint4(i,km-1) - tlayr4(i,km)
            delt1(i) = tlayr4(i,km) - tint4(i,km)
          end do
          do k = kzp1 , 1 , -1
            if ( k.eq.km ) then
              do i = 1 , iym1
                bk2(i) = absnxt(i,km-1,4,jslc)
                bk1(i) = absnxt(i,km-1,1,jslc)
              end do
            else if ( k.eq.km-1 ) then
              do i = 1 , iym1
                bk2(i) = absnxt(i,km-1,2,jslc)
                bk1(i) = absnxt(i,km-1,3,jslc)
              end do
            else
              do i = 1 , iym1
                bk2(i) = (abstot(i,k,km-1,jslc)+abstot(i,k,km,jslc))*0.5
                bk1(i) = bk2(i)
              end do
            end if
            do i = 1 , iym1
              s(i,k,km) = s(i,k,km+1)                                   &
                      & + stebol*(bk2(i)*delt(i)+bk1(i)*delt1(i))
            end do
          end do
        end do
!
!       Computation of clear sky fluxes always set first level of fsul
!
        do i = 1 , iym1
          if ( iemiss.eq.1 ) then
            fsul(i,kzp1) = emiss1d(i)*(stebol*(ts(i)**4))
          else
            fsul(i,kzp1) = stebol*(ts(i)**4)
          end if
        end do
!
!       Downward clear sky fluxes store intermediate quantities in down
!       flux Initialize fluxes to clear sky values.
!
        do i = 1 , iym1
          tmp(i) = fsul(i,kzp1) - stebol*tint4(i,kzp1)
          fsul(i,1) = fsul(i,kzp1) - abstot(i,1,kzp1,jslc)*tmp(i)       &
                  & + s(i,1,2)
          fsdl(i,1) = stebol*(tplnke(i)**4)*emstot(i,1,jslc)
          ful(i,1) = fsul(i,1)
          fdl(i,1) = fsdl(i,1)
        end do
!
!       fsdl(i,kzp1) assumes isothermal layer
!
        do k = 2 , kz
          do i = 1 , iym1
            fsul(i,k) = fsul(i,kzp1) - abstot(i,k,kzp1,jslc)*tmp(i)     &
                    & + s(i,k,k+1)
            ful(i,k) = fsul(i,k)
            fsdl(i,k) = stebol*(tplnke(i)**4)*emstot(i,k,jslc)          &
                    & - (s(i,k,2)-s(i,k,k+1))
            fdl(i,k) = fsdl(i,k)
          end do
        end do
!
!       Store the downward emission from level 1 = total gas emission *
!       sigma t**4.  fsdl does not yet include all terms
!
        do i = 1 , iym1
          ful(i,kzp1) = fsul(i,kzp1)
          absbt(i) = stebol*(tplnke(i)**4)*emstot(i,kzp1,jslc)
          fsdl(i,kzp1) = absbt(i) - s(i,kzp1,2)
          fdl(i,kzp1) = fsdl(i,kzp1)
        end do

!     FAB radiative forcing sur fsul

        if ( rad==1 ) then
          fsul0(:,:) = fsul(:,:)! save fsul0 = no dust
          fsdl0(:,:) = fsdl(:,:)!
        end if

      end do ! end rad loop

!     FAB after this DO loop fsul account for dust LW effect
!     which is OK in case of idirect=2

      aerlwfo(:) = fsul0(:,1) - fsul(:,1)

!     surface lw net ! fsul(i,plevp) - fsdl(i,plevp)
!     aerlwfos(:)= fsdl0(:,kz)-fsdl(:,kz)
      aerlwfos(:) = (fsul0(:,kzp1)-fsdl0(:,kzp1))-                      &
               &    (fsul(:,kzp1) - fsdl(:,kzp1))
         
!     return to no aerosol LW effect  situation if idirect ==1
      if ( idirect==1 ) then
        fsul(:,:) = fsul0(:,:)
        fsdl(:,:) = fsdl0(:,:)
      end if 

!
!     Modifications for clouds
!
!     Further qualify longitude subset for computations.  Select only
!     those locations where there are clouds (total cloud fraction <=
!     1.e-3 treated as clear)
!
      call whenflt(iym1,tclrsf(1,kzp1),1,0.999D0,indx,iym1c)
!
!     Compute downflux at level 1 for cloudy sky
!
      do ii = 1 , iym1c
        i = indx(ii)
!
!       First clear sky flux plus flux from cloud at level 1
!
        fdl(i,kzp1) = fsdl(i,kzp1)*tclrsf(i,kz)                     &
                     & *rtclrsf(i,kzp1-khiv(i)) + fclb4(i,kz-1)      &
                     & *cld(i,kz)
      end do
!
!     Flux emitted by other layers
!     Note: Vertical indexing here proceeds from bottom to top
!
      khighest = khiv(intmax(iym1,khiv,1))
      do km = 3 , khighest
        km1 = kzp1 - km
        km2 = kzp2 - km
        km4 = kzp4 - km
        do ii = 1 , iym1c
          i = indx(ii)
          if ( km.le.khiv(i) ) then
            tmp1 = cld(i,km2)*tclrsf(i,kz)*rtclrsf(i,km2)
            fdl(i,kzp1) = fdl(i,kzp1) + (fclb4(i,km1)-s(i,kzp1,km4)) &
                         & *tmp1
          end if
        end do
      end do
!
!     Note: Vertical indexing here proceeds from bottom to top
!
      do k = 1 , khighest - 1
        k1 = kzp1 - k
        k2 = kzp2 - k
        k3 = kzp3 - k
        do ii = 1 , iym1c
          i = indx(ii)
          if ( k.ge.klov(i) .and. k.le.khivm(i) ) ful(i,k2) = fsul(i,k2)&
             & *(tclrsf(i,kzp1)*rtclrsf(i,k1))
        end do
        do km = 1 , k
          km1 = kzp1 - km
          km2 = kzp2 - km
          km3 = kzp3 - km
          do ii = 1 , iym1c
            i = indx(ii)
!
            if ( k.le.khivm(i) .and. km.ge.klov(i) .and. km.le.khivm(i) &
               & ) ful(i,k2) = ful(i,k2)                                &
                             & + (fclt4(i,km1)+s(i,k2,k3)-s(i,k2,km3))  &
                             & *cld(i,km2)*(tclrsf(i,km1)*rtclrsf(i,k1))
          end do
        end do             ! km=1,k
      end do               ! k=1,khighest-1
!
      do k = 1 , kzp1
        k2 = kzp2 - k
        k3 = kzp3 - k
        do i = 1 , iym1
          start(i) = .false.
        end do
        do ii = 1 , iym1c
          i = indx(ii)
          if ( k.ge.khiv(i) ) then
            start(i) = .true.
            ful(i,k2) = fsul(i,k2)*tclrsf(i,kzp1)                      &
                      & *rtclrsf(i,kzp1-khiv(i))
          end if
        end do
        do km = 1 , khighest
          km1 = kzp1 - km
          km2 = kzp2 - km
          km3 = kzp3 - km
          do ii = 1 , iym1c
            i = indx(ii)
            if ( start(i) .and. km.ge.klov(i) .and. km.le.khiv(i) )     &
               & ful(i,k2) = ful(i,k2)                                  &
                           & + (cld(i,km2)*tclrsf(i,km1)*rtclrsf(i,     &
                           & kzp1-khiv(i)))                            &
                           & *(fclt4(i,km1)+s(i,k2,k3)-s(i,k2,km3))
          end do
        end do          ! km=1,khighest
      end do            ! k=1,kzp1
!
!     Computation of the downward fluxes
!
      do k = 2 , khighest - 1
        k1 = kzp1 - k
        k2 = kzp2 - k
        k3 = kzp3 - k
        do ii = 1 , iym1c
          i = indx(ii)
          if ( k.le.khivm(i) ) fdl(i,k2) = 0.
        end do
        do km = k + 1 , khighest
          km1 = kzp1 - km
          km2 = kzp2 - km
          km4 = kzp4 - km
          do ii = 1 , iym1c
            i = indx(ii)
!
            if ( k.le.khiv(i) .and. km.ge.max0(k+1,klov(i)) .and.       &
               & km.le.khiv(i) ) fdl(i,k2) = fdl(i,k2)                  &
               & + (cld(i,km2)*tclrsf(i,k1)*rtclrsf(i,km2))             &
               & *(fclb4(i,km1)-s(i,k2,km4)+s(i,k2,k3))
          end do
        end do             ! km=k+1,khighest
        do ii = 1 , iym1c
          i = indx(ii)
          if ( k.le.khivm(i) ) fdl(i,k2) = fdl(i,k2) + fsdl(i,k2)       &
             & *(tclrsf(i,k1)*rtclrsf(i,kzp1-khiv(i)))
        end do
      end do               ! k=1,khighest-1
!
!     End cloud modification loops
!
!     All longitudes: store history tape quantities
!
      do i = 1 , iym1
!
!       Downward longwave flux
!
        flwds(i) = fdl(i,kzp1)
!
!       Net flux
!
        flns(i) = ful(i,kzp1) - fdl(i,kzp1)
!
!       Clear sky flux at top of atmosphere
!
        flntc(i) = fsul(i,1)
        flnsc(i) = fsul(i,kzp1) - fsdl(i,kzp1)
!
!       Outgoing ir
!
        flnt(i) = ful(i,1) - fdl(i,1)
      end do
!
!     Computation of longwave heating (k per sec)
!
      do k = 1 , kz
        do i = 1 , iym1
          qrl(i,k) = (ful(i,k)-fdl(i,k)-ful(i,k+1)+fdl(i,k+1))          &
                   & *gocp/((pint(i,k)-pint(i,k+1)))
        end do
      end do
!
      end subroutine radclw
