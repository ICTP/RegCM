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
!
!     This module implements Kuo cumulus parameterization scheme.
!     The basic method follows Anthes and Keyser (1979) and Kuo (1983).
!
      use mod_runparams
      use mod_main
      use mod_cvaria
      use mod_pmoist
      use mod_rad
      use mod_bats
      use mod_trachem
      use mod_date
      use mod_advection

      private

      public :: cupara , htdiff
!
!     qdcrit : the precipitation threshold for moisture convergence.
!     pert   : perturbation temperature
!     perq   : perturbation mixing ratio
!     dlt    : temperature difference used to allow over shooting.
!     cdscld : critical cloud depth in delta sigma.
!
      real(8) , parameter :: qdcrit = 3.0D-7
      real(8) , parameter :: pert   = 1.0D0
      real(8) , parameter :: perq   = 1.0D-3
      real(8) , parameter :: dlt    = 3.0D0
      real(8) , parameter :: cdscld = 0.3D0
!
      contains
!
      subroutine cupara(j)
!
!     All the other arguments are passed from subroutine "tend" and
!     explained in "tend".
!
      implicit none
!
      integer :: j
!
      real(8) :: akclth , apcnt , aprdiv , arh , c301 , dalr ,    &
               & deqt , dlnp , dplr , dsc , e1 , eddyf , emax ,   &
               & eqt , eqtm , es , plcl , pmax , prainx , psg ,   &
               & psx , pux , q , qmax , qs , rh , rsht , rswt ,   &
               & sca , siglcl , suma , sumb , t1 , tdmax , tlcl , &
               & tmax , tmean , ttconv , ttp , ttsum , xsav , zlcl
      integer :: i , k , kbase , kbaseb , kclth , kk , ktop
      real(8) , dimension(kz) :: seqt
      real(8) , dimension(iy,kz) :: tmp3
!
!----------------------------------------------------------------------
!
!
      pmax = d_zero
      qmax = d_zero
      tmax = d_zero
      do k = 1 , kz
        do i = 1 , iym1
          cldlwc(i,k) = d_zero
          cldfra(i,k) = d_zero
        end do
      end do
      do k = 1 , kz
        do i = 1 , iym1
          aten%qv(i,k,j) = d_zero
        end do
      end do
!
!-----compute the horizontal advection terms:
!
      call hadv_x(aten%qv(:,:,j),atmx%qv,dx4,j,1)
!---------------
!chem2
      if ( ichem == 1 ) then
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
!     at this stage, aten%qv(i,k,j) only includes horizontal advection.
!     sca: is the amount of total moisture convergence
!
      do i = 2 , iym2
!
        sca = d_zero
        do k = 1 , kz
          sca = sca + aten%qv(i,k,j)*dsigma(k)
        end do
!
!-----determine if moist convection exists:
!
        if ( sca >= qdcrit ) then
!
!-----check for stability
!
!--1--compute eqt (equivalent potential temperature)
!         between surface and 700 mb, with perturbation temperature
!         and moisture added. the maximum eqt will be regarded
!         as the origin of air parcel that produce cloud.
!
          eqtm = d_zero
          do k = k700 , kz
            ttp = atm1%t(i,k,j)/sps1%ps(i,j) + pert
            q = atm1%qv(i,k,j)/sps1%ps(i,j) + perq
            psg = sps1%ps(i,j)*a(k) + r8pt
            t1 = ttp*(d_100/psg)**rovcp
            eqt = t1*dexp(wlhvocp*q/ttp)
            if ( eqt > eqtm ) then
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
          tdmax = 5418.12D0/(19.84659D0-dlog(emax/0.611D0))
          dalr = egrav*rcpd
          dplr = (egrav*tdmax*tdmax)/(ep2*wlhv*tmax)
          zlcl = (tmax-tdmax)/(dalr-dplr)
          tlcl = tmax - dalr*zlcl
          tmean = (tmax+tlcl)*d_half
          dlnp = (egrav*zlcl)/(rgas*tmean)
          plcl = pmax*dexp(-dlnp)
          siglcl = (plcl-r8pt)/sps1%ps(i,j)
!
!--3--compute seqt (saturation equivalent potential temperature)
!         of all the levels that are above the lcl
!
          do k = 1 , kz
            if ( a(k) >= siglcl ) exit
          end do
          kbase = k
          if ( kbase > kz ) kbase = kz
!
!.....kbase is the layer where lcl is located.
!
          do k = 1 , kbase
            ttp = atm1%t(i,k,j)/sps1%ps(i,j)
            psg = sps1%ps(i,j)*a(k) + r8pt
            es = 0.611D0*dexp(19.84659D0-5418.12D0/ttp)
            qs = ep2*es/(psg-es)
            t1 = ttp*(d_100/psg)**rovcp
            seqt(k) = t1*dexp(wlhvocp*qs/ttp)
          end do
!
!--4--when seqt = eqt + dt, cloud top is reached.
!         eqt is the eqt of cloud (same as lcl eqt).
!
          do kk = 1 , kbase
            k = kbase + 1 - kk
            deqt = seqt(k) - eqtm
            if ( deqt > dlt ) exit
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
          if ( dsc >= cdscld ) then
!
!--6--check negative area
!           if negative area is larger than the positive area
!           convection is killed.
!
            ttsum = d_zero
            do k = ktop , kbase
              ttsum = (eqtm-seqt(k))*dsigma(k) + ttsum
            end do
            if ( ttsum >= d_zero ) then
!
!.....you     are here if stability was found.
!
!.....if      values dont already exist in array twght,vqflx for this
!             kbase/ktop, then flag it, and set kbase/ktop to standard
!
              if ( (kbase < 5) .or. (ktop > kbase-3) ) then
                print 99001 , ktau , jyear , i , j , kbase , ktop
                if ( kbase < 5 ) kbase = 5
                if ( ktop > kbase-3 ) ktop = kbase - 3
              end if
!
!.....convection exist, compute convective flux of water vapor and
!             latent heating
!             icon   : is a counter which keep track the total points
!             where deep convection occurs.
!             c301   : is the 'b' factor in kuo's scheme.
!
              icon(j) = icon(j) + 1
              suma = d_zero
              sumb = d_zero
              arh = d_zero
              psx = sps1%ps(i,j)
              do k = 1 , kz
                qwght(k) = d_zero
              end do
              do k = ktop , kz
                pux = psx*a(k) + r8pt
                e1 = 0.611D0*dexp(19.84659D0-5418.12D0/ &
                                  (atm1%t(i,k,j)/psx))
                qs = ep2*e1/(pux-e1)
                rh = atm1%qv(i,k,j)/(qs*psx)
                rh = dmin1(rh,d_one)
                xsav = (d_one-rh)*qs
                qwght(k) = xsav
                sumb = sumb + qs*dsigma(k)
                arh = arh + rh*qs*dsigma(k)
                suma = suma + xsav*dsigma(k)
              end do
              arh = arh/sumb
              c301 = d_two*(d_one-arh)
              if ( c301 < d_zero ) c301 = d_zero
              if ( c301 > d_one ) c301 = d_one
              if ( suma <= d_zero ) then
                c301 = d_zero
                suma = d_one
              end if
              do k = ktop , kz
                qwght(k) = qwght(k)/suma
              end do
              do k = 1 , kz
                ttconv = wlhvocp*(d_one-c301)*twght(k,kbase,ktop)*sca
                rsheat(i,k,j) = rsheat(i,k,j) + ttconv*dt*d_half
!x              if (ttconv*2. > 0.01) write(18,1234) i,j,k,ttconv*2.
!1234           format(1x,'cupara, i=',i4,' j=',i4,' k=',i4,'
!               qteva=',e12.4)
                apcnt = (d_one-c301)*sca/4.3D-3
                eddyf = apcnt*vqflx(k,kbase,ktop)
                aten%qv(i,k,j) = eddyf
                rswat(i,k,j) = rswat(i,k,j) + &
                               c301*qwght(k)*sca*dt*d_half
              end do
!
!             find cloud fractional cover and liquid water content
!
              kbaseb = min0(kbase,kzm2)
              if ( ktop <= kbaseb ) then
                kclth = kbaseb - ktop + 1
                akclth = d_one/dble(kclth)
                do k = ktop , kbaseb
                  cldlwc(i,k) = cllwcv
                  cldfra(i,k) = d_one - (d_one-clfrcv)**akclth
                end do
              end if
!.....the     unit for rainfall is mm.
              prainx = (d_one-c301)*sca*dtsec*1000.0D0*regrav
              if ( prainx > dlowval ) then
                sfsta%rainc(i,j) = sfsta%rainc(i,j) + prainx
!               instantaneous precipitation rate for use in bats (mm/s)
                aprdiv = dble(nbatst)
                if ( jyear == jyear0 .and. ktau == 0 ) aprdiv = d_one
                pptc(i,j) = pptc(i,j) + prainx/dtsec/aprdiv
              end if
!
!chem2
              if ( ichem == 1 ) then
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
        tmp3(i,1) = d_zero
        do k = 2 , kz
          if ( atm1%qv(i,k,j) < 1.0D-15 ) then
            tmp3(i,k) = d_zero
          else
            tmp3(i,k) = atm1%qv(i,k,j)*(atm1%qv(i,k-1,j)/ &
                        atm1%qv(i,k,j))**qcon(k)
          end if
        end do
        aten%qv(i,1,j) = aten%qv(i,1,j)-qdot(i,2,j)*tmp3(i,2)/dsigma(1)
        do k = 2 , kzm1
          aten%qv(i,k,j) = aten%qv(i,k,j)                               &
                       & - (qdot(i,k+1,j)*tmp3(i,k+1)-qdot(i,k,j)       &
                       & *tmp3(i,k))/dsigma(k)
        end do
        aten%qv(i,kz,j) = aten%qv(i,kz,j) + qdot(i,kz,j)*tmp3(i,kz)     &
                      & /dsigma(kz)
!
      end do           !end i=2,iym2 loop
!
      do k = 1 , kz
        do i = 2 , iym2
          rsheat(i,k,j) = dmax1(rsheat(i,k,j),d_zero)
          rswat(i,k,j) = dmax1(rswat(i,k,j),d_zero)
          rsht = rsheat(i,k,j)/tauht
          rswt = rswat(i,k,j)/tauht
          aten%t(i,k,j) = aten%t(i,k,j) + rsht
          aten%qv(i,k,j) = aten%qv(i,k,j) + rswt
          rsheat(i,k,j) = rsheat(i,k,j)*(d_one-dt/(d_two*tauht))
          rswat(i,k,j) = rswat(i,k,j)*(d_one-dt/(d_two*tauht))
!bxq if (rsht/psb(i,j) > .0002)write(18,1222)ktau,jyear,i,j,k,rsht/psb(i
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
      subroutine htdiff(dxsq,akht1)

#ifdef MPP1
#ifndef IBM
      use mpi
#else 
      include 'mpif.h'
#endif 
#endif
      implicit none
!
      real(8) :: akht1 , dxsq
      intent (in) akht1 , dxsq
!
      integer :: i , im1 , ip1 , j , jm1 , jp1 , k
#ifdef MPP1
      integer :: ierr
      real(8) , dimension(iy,0:jxp+1) :: wr
#else
      real(8) , dimension(iy,jx) :: wr
#endif
!
#ifdef MPP1
      do k = 1 , kz
        do j = 1 , jendl
          do i = 1 , iy
            wr(i,j) = rsheat(i,k,j)
          end do
        end do
        call mpi_sendrecv(wr(1,jxp),iy,mpi_real8,ieast,1,               &
                        & wr(1,0),iy,mpi_real8,iwest,1,                 &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(wr(1,1),iy,mpi_real8,iwest,2,                 &
                        & wr(1,jxp+1),iy,mpi_real8,ieast,2,             &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        do j = jbegin , jendm
#ifdef BAND
          jm1 = j - 1
          jp1 = j + 1
#else
          if ( myid == 0 ) then
            jm1 = max0(j-1,2)
          else
            jm1 = j - 1
          end if
          if ( myid == nproc-1 ) then
            jp1 = min0(j+1,jxp-2)
          else
            jp1 = j + 1
          end if
#endif
          do i = 2 , iym2
            im1 = max0(i-1,2)
            ip1 = min0(i+1,iym2)
            rsheat(i,k,j) = rsheat(i,k,j)                               &
                          & + akht1*dt*d_half/dxsq*(wr(im1,j)+wr(ip1,j) &
                          & +wr(i,jm1)+wr(i,jp1)-d_four*wr(i,j))
          end do
        end do
      end do
#else
      do k = 1 , kz
        do j = 1 , jx
          do i = 1 , iy
            wr(i,j) = rsheat(i,k,j)
          end do
        end do
#ifdef BAND
        do j = 1 , jx
          jm1 = j - 1
          jp1 = j + 1
          if (jm1 == 0) jm1 = jx
          if (jp1 == jx+1) jp1 = 1
#else
        do j = 2 , jxm2
          jm1 = max0(j-1,2)
          jp1 = min0(j+1,jxm2)
#endif
          do i = 2 , iym2
            im1 = max0(i-1,2)
            ip1 = min0(i+1,iym2)
            rsheat(i,k,j) = rsheat(i,k,j)                               &
                          & + akht1*dt*d_half/dxsq*(wr(im1,j)+wr(ip1,j) &
                          & +wr(i,jm1)+wr(i,jp1)-d_four*wr(i,j))
          end do
        end do
      end do
#endif
!
#ifdef MPP1
      do k = 1 , kz
        do j = 1 , jendl
          do i = 1 , iy
            wr(i,j) = rswat(i,k,j)
          end do
        end do
        call mpi_sendrecv(wr(1,jxp),iy,mpi_real8,ieast,1,               &
                        & wr(1,0),iy,mpi_real8,iwest,1,                 &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(wr(1,1),iy,mpi_real8,iwest,2,                 &
                        & wr(1,jxp+1),iy,mpi_real8,ieast,2,             &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        do j = jbegin , jendm
#ifdef BAND
          jm1 = j - 1
          jp1 = j + 1
#else
          if ( myid == 0 ) then
            jm1 = max0(j-1,2)
          else
            jm1 = j - 1
          end if
          if ( myid == nproc-1 ) then
            jp1 = min0(j+1,jxp-2)
          else
            jp1 = j + 1
          end if
#endif
          do i = 2 , iym2
            im1 = max0(i-1,2)
            ip1 = min0(i+1,iym2)
            rswat(i,k,j) = rswat(i,k,j)                                 &
                         & + akht1*dt*d_half/dxsq*(wr(im1,j)+wr(ip1,j)  &
                         & +wr(i,jm1)+wr(i,jp1)-d_four*wr(i,j))
          end do
        end do
      end do
#else
      do k = 1 , kz
        do j = 1 , jx
          do i = 1 , iy
            wr(i,j) = rswat(i,k,j)
          end do
        end do
#ifdef BAND
        do j = 1 , jx
          jm1 = j - 1
          jp1 = j + 1
          if (jm1 == 0) jm1 = jx
          if (jp1 == jx+1) jp1 = 1
#else
        do j = 2 , jxm2
          jm1 = max0(j-1,2)
          jp1 = min0(j+1,jxm2)
#endif
          do i = 2 , iym2
            im1 = max0(i-1,2)
            ip1 = min0(i+1,iym2)
            rswat(i,k,j) = rswat(i,k,j)                                 &
                         & + akht1*dt*d_half/dxsq*(wr(im1,j)+wr(ip1,j)  &
                         & +wr(i,jm1)+wr(i,jp1)-d_four*wr(i,j))
          end do
        end do
      end do
#endif
      end subroutine htdiff
!
      end module mod_cu_kuo
