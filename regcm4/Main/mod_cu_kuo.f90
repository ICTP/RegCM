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
 
      module mod_cu_kuo

      use mod_constants
      use mod_dynparam
      use mod_runparams
      use mod_main
      use mod_cvaria
      use mod_pmoist
      use mod_rad
      use mod_bats
      use mod_trachem
      use mod_date

      private

      public :: cupara

!     qdcrit is the precipitation threshold for moisture convergence.
      real(8) :: qdcrit
      data qdcrit /3.0E-7/

      contains

      subroutine cupara(j)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine performs cumulus parameterization scheme.       c
!     the basic method follows anthes and keyser (1979) and           c
!     kuo (1983).                                                     c
!                                                                     c
!     all the other arguments are passed from subroutine "tend" and   c
!     explained in "tend".                                            c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
!
! Dummy arguments
!
      integer :: j
!
! Local variables
!
      real(8) :: akclth , apcnt , aprdiv , arh , c301 , cdscld , dalr , &
               & deqt , dlnp , dlt , dplr , dsc , e1 , eddyf , emax ,   &
               & eqt , eqtm , es , perq , pert , plcl , pmax , prainx , &
               & psg , psx , pux , q , qmax , qs , rh , rsht , rswt ,   &
               & sca , siglcl , suma , sumb , t1 , tdmax , tlcl , tmax ,&
               & tmean , ttconv , ttp , ttsum , xsav , zlcl
      integer :: i , k , kbase , kbaseb , kclth , kk , ktop
      real(8) , dimension(kz) :: seqt
      real(8) , dimension(iy,kz) :: tmp3
!
!
!----------------------------------------------------------------------
!
!     pert   : perturbation temperature
!     perq   : perturbation mixing ratio
!     dlt    : temperature difference used to allow over shooting.
!     cdscld : critical cloud depth in delta sigma.
!
      data pert , perq/1. , 1.E-3/
      data dlt , cdscld/3.0 , 0.3/
!
!
      pmax = 0.0
      qmax = 0.0
      tmax = 0.0
      do k = 1 , kz
        do i = 1 , iym1
          cldlwc(i,k) = 0.
          cldfra(i,k) = 0.
        end do
      end do
      do k = 1 , kz
        do i = 1 , iym1
          qvten(i,k,j) = 0.
        end do
      end do
!
!-----compute the horizontal advection terms:
!
      call hadvqv(qvten(1,1,j),dx4,j,1)
!---------------
!chem2
      if ( ichem.eq.1 ) then
!
!       icumtop = top level of cumulus clouds
!       icumbot = bottom level of cumulus clouds
!       (calculated in cupara and stored for tractend)
!       before do 100 put
        do i = 2 , iym2
          icumtop(i,j) = 0
          icumbot(i,j) = 0
        end do
      end if
!chem2__
!
!-----compute the moisture convergence in a column:
!     at this stage, qvten(i,k,j) only includes horizontal advection.
!     sca: is the amount of total moisture convergence
!
      do i = 2 , iym2
!
        sca = 0.0
        do k = 1 , kz
          sca = sca + qvten(i,k,j)*dsigma(k)
        end do
!
!-----determine if moist convection exists:
!
        if ( sca.ge.qdcrit ) then
!
!-----check for stability
!
!--1--compute eqt (equivalent potential temperature)
!         between surface and 700 mb, with perturbation temperature
!         and moisture added. the maximum eqt will be regarded
!         as the origin of air parcel that produce cloud.
!
          eqtm = 0.0
          do k = k700 , kz
            ttp = ta(i,k,j)/psa(i,j) + pert
            q = qva(i,k,j)/psa(i,j) + perq
            psg = psa(i,j)*a(k) + r8pt
            t1 = ttp*(100./psg)**rovcp
            eqt = t1*dexp(wlhvocp*q/ttp)
            if ( eqt.gt.eqtm ) then
              eqtm = eqt
              tmax = ttp
              qmax = q
              pmax = psg
            end if
          end do
!
!--2--compute lcl, get the sigma and p of lcl
!
          emax = qmax*pmax/(ep2+qmax)
          tdmax = 5418.12/(19.84659-dlog(emax/.611))
          dalr = gti*rcpd
          dplr = (gti*tdmax*tdmax)/(ep2*wlhv*tmax)
          zlcl = (tmax-tdmax)/(dalr-dplr)
          tlcl = tmax - dalr*zlcl
          tmean = 0.5*(tmax+tlcl)
          dlnp = (gti*zlcl)/(rgas*tmean)
          plcl = pmax*dexp(-dlnp)
          siglcl = (plcl-r8pt)/psa(i,j)
!
!--3--compute seqt (saturation equivalent potential temperature)
!         of all the levels that are above the lcl
!
          do k = 1 , kz
            if ( a(k).ge.siglcl ) exit
          end do
          kbase = k
          if ( kbase.gt.kz ) kbase = kz
!
!.....kbase is the layer where lcl is located.
!
          do k = 1 , kbase
            ttp = ta(i,k,j)/psa(i,j)
            psg = psa(i,j)*a(k) + r8pt
            es = .611*dexp(19.84659-5418.12/ttp)
            qs = ep2*es/(psg-es)
            t1 = ttp*(100./psg)**rovcp
            seqt(k) = t1*dexp(wlhvocp*qs/ttp)
          end do
!
!--4--when seqt = eqt + dt, cloud top is reached.
!         eqt is the eqt of cloud (same as lcl eqt).
!
          do kk = 1 , kbase
            k = kbase + 1 - kk
            deqt = seqt(k) - eqtm
            if ( deqt.gt.dlt ) exit
          end do
!
!.....cloud top has been reached
!
          ktop = k
!
!--5--check cloud depth
!         if cloud depth is less than critical depth (cdscld = 0.3),
!         the convection is killed
!
          dsc = (siglcl-a(ktop))
          if ( dsc.ge.cdscld ) then
!
!--6--check negative area
!           if negative area is larger than the positive area
!           convection is killed.
!
            ttsum = 0.
            do k = ktop , kbase
              ttsum = (eqtm-seqt(k))*dsigma(k) + ttsum
            end do
            if ( ttsum.ge.0. ) then
!
!.....you     are here if stability was found.
!
!.....if      values dont already exist in array twght,vqflx for this
!             kbase/ktop, then flag it, and set kbase/ktop to standard
!
              if ( (kbase.lt.5) .or. (ktop.gt.kbase-3) ) then
                print 99001 , ktau , jyear , i , j , kbase , ktop
                if ( kbase.lt.5 ) kbase = 5
                if ( ktop.gt.kbase-3 ) ktop = kbase - 3
              end if
!
!.....convection exist, compute convective flux of water vapor and
!             latent heating
!             icon   : is a counter which keep track the total points
!             where deep convection occurs.
!             c301   : is the 'b' factor in kuo's scheme.
!
              icon(j) = icon(j) + 1
              suma = 0.
              sumb = 0.
              arh = 0.
              psx = psa(i,j)
              do k = 1 , kz
                qwght(k) = 0.0
              end do
              do k = ktop , kz
                pux = psx*a(k) + r8pt
                e1 = .611*dexp(19.84659-5418.12/(ta(i,k,j)/psx))
                qs = ep2*e1/(pux-e1)
                rh = qva(i,k,j)/(qs*psx)
                rh = dmin1(rh,1.D0)
                xsav = (1.0-rh)*qs
                qwght(k) = xsav
                sumb = sumb + qs*dsigma(k)
                arh = arh + rh*qs*dsigma(k)
                suma = suma + xsav*dsigma(k)
              end do
              arh = arh/sumb
              c301 = 2.0*(1.0-arh)
              if ( c301.lt.0.0 ) c301 = 0.0
              if ( c301.gt.1.0 ) c301 = 1.0
              if ( suma.le.0.0 ) then
                c301 = 0.0
                suma = 1.0
              end if
              do k = ktop , kz
                qwght(k) = qwght(k)/suma
              end do
              do k = 1 , kz
                ttconv = wlhvocp*(1.0-c301)*twght(k,kbase,ktop)*sca
                rsheat(i,k,j) = rsheat(i,k,j) + ttconv*dt/2.
!x              if (ttconv*2. .gt. 0.01) write(18,1234) i,j,k,ttconv*2.
!1234           format(1x,'cupara, i=',i4,' j=',i4,' k=',i4,'
!               qteva=',e12.4)
                apcnt = (1.0-c301)*sca/4.3E-3
                eddyf = apcnt*vqflx(k,kbase,ktop)
                qvten(i,k,j) = eddyf
                rswat(i,k,j) = rswat(i,k,j) + c301*qwght(k)*sca*dt/2.
              end do
!
!             find cloud fractional cover and liquid water content
!
              kbaseb = min0(kbase,kzm2)
              if ( ktop.le.kbaseb ) then
                kclth = kbaseb - ktop + 1
                akclth = 1./dble(kclth)
                do k = ktop , kbaseb
                  cldlwc(i,k) = cllwcv
                  cldfra(i,k) = 1. - (1.-clfrcv)**akclth
                end do
              end if
!.....the     unit for rainfall is mm.
              prainx = (1.-c301)*sca*dtmin*60000.*rgti
              rainc(i,j) = rainc(i,j) + prainx
!             instantaneous precipitation rate for use in bats (mm/s)
              aprdiv = dble(nbatst)
              if ( jyear.eq.jyear0 .and. ktau.eq.0 ) aprdiv = 1.
              pptc(i,j) = pptc(i,j) + prainx/(dtmin*60.)/aprdiv
!
!chem2
              if ( ichem.eq.1 ) then
!               before go to 100 put
                icumtop(i,j) = ktop
                icumbot(i,j) = kbaseb
              end if
!chem2_
 
              cycle
            end if
          end if
        end if
!
!.....convection not exist, compute the vertical advection term:
!
        tmp3(i,1) = 0.
        do k = 2 , kz
          if ( qva(i,k,j).lt.1.E-15 ) then
            tmp3(i,k) = 0.0
          else
            tmp3(i,k) = qva(i,k,j)*(qva(i,k-1,j)/qva(i,k,j))**qcon(k)
          end if
        end do
        qvten(i,1,j) = qvten(i,1,j) - qdot(i,2,j)*tmp3(i,2)/dsigma(1)
        do k = 2 , kzm1
          qvten(i,k,j) = qvten(i,k,j)                                   &
                       & - (qdot(i,k+1,j)*tmp3(i,k+1)-qdot(i,k,j)       &
                       & *tmp3(i,k))/dsigma(k)
        end do
        qvten(i,kz,j) = qvten(i,kz,j) + qdot(i,kz,j)*tmp3(i,kz)         &
                      & /dsigma(kz)
!
      end do           !end i=2,iym2 loop
!
      do k = 1 , kz
        do i = 2 , iym2
          rsheat(i,k,j) = dmax1(rsheat(i,k,j),0.D0)
          rswat(i,k,j) = dmax1(rswat(i,k,j),0.D0)
          rsht = rsheat(i,k,j)/tauht
          rswt = rswat(i,k,j)/tauht
          tten(i,k,j) = tten(i,k,j) + rsht
          qvten(i,k,j) = qvten(i,k,j) + rswt
          rsheat(i,k,j) = rsheat(i,k,j)*(1.-dt/(2.*tauht))
          rswat(i,k,j) = rswat(i,k,j)*(1.-dt/(2.*tauht))
!bxq if(rsht/psb(i,j).gt..0002)write(18,1222)ktau,jyear,i,j,k,rsht/psb(i
!1222     format (1x,'ktau= ',i7,' jyear= ',i5,' i= ',i5,' j= ',i5,
!         1        ' k= ',i5,' ttconv =',e15.7)
        end do
      end do
99001 format (/,' >>in **cupara**: at ktau=',i8,' in year=',i5,         &
             &' & (i,j)=(',i2,',',i2,'),   ',                           &
             &' kbase/ktop are non-standard:',2I3,                      &
             &'  & will be set to closest standard.')
!
      end subroutine cupara
!
!
!
      end module mod_cu_kuo
