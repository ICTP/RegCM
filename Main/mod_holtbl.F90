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
 
      module mod_holtbl
!
! Holtslag planetary boundary layer scheme
! Reference : Holtslag, De Bruijn and Pan - MWR - 8/90
!
      use mod_runparams
      use mod_main
      use mod_mainchem
      use mod_pbldim
      use mod_cvaria
      use mod_pmoist
      use mod_bats
      use mod_slice
      use mod_trachem
      use mod_diagnosis
#ifdef MPP1
      use mod_mppio
#endif
!
      private
!
      public :: allocate_mod_holtbl , holtbl
!
      real(8) , allocatable , dimension(:,:,:) :: cgh , kvc , kvh ,   &
                                                  kvm , kvq ! , cgq
      real(8) , allocatable, dimension(:,:) :: hfxv , obklen , th10 , &
                                               ustr , xhfx , xqfx
!
      real(8) , allocatable , dimension(:,:) :: alphak , betak , chix , &
                               coef1 , coef2 , coef3 , coefe , coeff1 , &
                               coeff2 , tpred1 , tpred2
      real(8) , allocatable , dimension(:,:) :: kzm , rc , ttnp
      real(8) , allocatable , dimension(:,:) :: vdep
      real(8) , allocatable , dimension(:) :: govrth , hydf
!
      real(8) , allocatable , dimension(:,:,:) :: auxx , avxx , dza , qcx
      real(8) , allocatable , dimension(:,:,:) :: akzz1 , akzz2
#ifdef MPP1
      real(8) , allocatable , dimension(:) :: wkrecv , wksend
#endif
      real(8) , allocatable , dimension(:,:,:) :: rhohf
!
      real(8) , allocatable , dimension(:,:) :: ri
      real(8) , allocatable , dimension(:) :: therm
!
!     minimum eddy diffusivity
      real(8) , parameter :: kzo = d_one
      real(8) , parameter :: szkm = 1600.0D0
!     coef. of proportionality and lower % of bl in sfc layer
      real(8) , parameter :: fak = 8.5D0
      real(8) , parameter :: sffrac = 0.1D0
!     beta coefs. for momentum, stable conditions and heat
      real(8) , parameter :: betam = 15.0D0
      real(8) , parameter :: betas = 5.0D0
      real(8) , parameter :: betah = 15.0D0
      real(8) , parameter :: mult = 0.61D0
      real(8) , parameter :: ccon = fak*sffrac*vonkar
      real(8) , parameter :: gvk = egrav*vonkar
      real(8) , parameter :: gpcf = egrav/d_1000
      real(8) , parameter :: binm = betam*sffrac
      real(8) , parameter :: binh = betah*sffrac
!     power in formula for k and critical ri for judging stability
      real(8) , parameter :: pink = d_two
      real(8) , parameter :: ricr = d_rfour
! 
      contains
!
      subroutine allocate_mod_holtbl
      implicit none
      allocate(alphak(iy,kz))
      allocate(betak(iy,kz))
      allocate(coef1(iy,kz))
      allocate(coef2(iy,kz))
      allocate(coef3(iy,kz))
      allocate(coefe(iy,kz))
      allocate(coeff1(iy,kz))
      allocate(coeff2(iy,kz))
      allocate(tpred1(iy,kz))
      allocate(tpred2(iy,kz))
      allocate(kzm(iym1,kz))
      allocate(rc(iym1,kz))
      allocate(ttnp(iym1,kz))
      allocate(govrth(iym1))
      allocate(hydf(kz))
      allocate(ri(iy,kz))
      allocate(therm(iy))
#ifdef MPP1
      allocate(cgh(iy,kz,jxp))
!     allocate(cgq(iy,kz,jxp))
      allocate(kvh(iy,kz,jxp))        
      allocate(kvm(iy,kz,jxp))        
      allocate(kvq(iy,kz,jxp))
      allocate(hfxv(iy,jxp))        
      allocate(obklen(iy,jxp))        
      allocate(th10(iy,jxp))        
      allocate(ustr(iy,jxp))        
      allocate(xhfx(iy,jxp))        
      allocate(xqfx(iy,jxp))        
      allocate(auxx(iym1,kz,jxp))
      allocate(avxx(iym1,kz,jxp))
      allocate(dza(iym1,kz,jxp))
      allocate(qcx(iym1,kz,jxp))
      allocate(akzz1(iym1,kz,0:jxp+1))
      allocate(akzz2(iym1,kz,0:jxp+1))
      allocate(wkrecv(2*iym2*kz))
      allocate(wksend(2*iym2*kz))
      allocate(rhohf(iy,kz,jxp))
#else 
#ifdef BAND
      allocate(cgh(iy,kz,jx))        
      allocate(kvh(iy,kz,jx))        
      allocate(kvm(iy,kz,jx))        
      allocate(kvq(iy,kz,jx))
      allocate(auxx(iym1,kz,jx))
      allocate(avxx(iym1,kz,jx))
      allocate(dza(iym1,kz,jx))
      allocate(qcx(iym1,kz,jx))
      allocate(akzz1(iym1,kz,0:jx))
      allocate(akzz2(iym1,kz,0:jx))
      allocate(rhohf(iy,kz,jx))
#else
      allocate(cgh(iy,kz,jxm1))        
      allocate(kvh(iy,kz,jxm1))        
      allocate(kvm(iy,kz,jxm1))        
      allocate(kvq(iy,kz,jxm1))
      allocate(auxx(iym1,kz,jxm1))
      allocate(avxx(iym1,kz,jxm1))
      allocate(dza(iym1,kz,jxm1))
      allocate(qcx(iym1,kz,jxm1))
      allocate(akzz1(iym1,kz,0:jxm1))
      allocate(akzz2(iym1,kz,0:jxm1))
      allocate(rhohf(iy,kz,jxm1))
#endif
      allocate(hfxv(iy,jx))        
      allocate(obklen(iy,jx))        
      allocate(th10(iy,jx))        
      allocate(ustr(iy,jx))        
      allocate(xhfx(iy,jx))        
      allocate(xqfx(iy,jx))        
#endif 
      if ( ichem == 1 ) then
        if ( ichdrdepo == 1 ) then
          allocate(chix(iy,kz))
          allocate(vdep(iym1,ntr))
          chix = d_zero
          vdep = d_zero
        end if
#ifdef MPP1
        allocate(kvc(iy,kz,jxp))        
#else
#ifdef BAND
        allocate(kvc(iy,kz,jx))        
#else
        allocate(kvc(iy,kz,jxm1))        
#endif
#endif
        kvc = d_zero
      end if
      cgh = d_zero
      kvh = d_zero
      kvm = d_zero
      kvq = d_zero
      hfxv = d_zero
      obklen = d_zero
      th10 = d_zero
      ustr = d_zero
      xhfx = d_zero
      xqfx = d_zero
      end  subroutine allocate_mod_holtbl
!
      subroutine holtbl
!
#ifdef MPP1
#ifndef IBM
      use mpi
#else 
      include 'mpif.h'
#endif 
#endif
!
      implicit none
!
      real(8) :: drgdot , dumr , kzmax , oblen , xps , ps2 , ri , &
                 sf , sh10 , ss , uflxsf , uflxsfx , vflxsf ,     &
                 vflxsfx , rdtx
      integer :: jdx , jm1
#ifndef BAND
      integer :: jdxm1
#endif
      integer :: i , idx , idxm1 , itr , j , k
#ifdef MPP1
      integer :: ierr , ii
#endif
!
#ifndef BAND
!
! *********************************************************************
!     diagnostic on total evaporation
! *********************************************************************
!
      if (debug_level > 2) call conqeva
!
#endif
!
!----------------------------------------------------------------------
!-----some of the storage spaces for high-resolution pbl
!     are use mod_to store the variables in this subroutine.
!     difft(i,k,j)   : temperature tendency (tten)
!     diffq(i,k,j)   : water vapor tendency (qvten)
!
!-----decouple flux-form variables to give u,v,t,theta,theta-vir.,
!     t-vir., qv, and qc at cross points and at ktau-1.
!
!     *** note ***
!     the boundary winds may not be adequately affected by friction,
!     so use only interior values of ubx3d and vbx3d to calculate
!     tendencies.
!
#ifdef MPP1
      call mpi_sendrecv(sps2%ps(1,jxp),iy,mpi_real8,ieast,1,       &
                        sps2%ps(1,0),  iy,mpi_real8,iwest,1,       &
                        mpi_comm_world,mpi_status_ignore,ierr)
      call mpi_sendrecv(sfsta%uvdrag(1,jxp),iy,mpi_real8,ieast,1,  &
                        sfsta%uvdrag(1,0),  iy,mpi_real8,iwest,1,  &
                        mpi_comm_world,mpi_status_ignore,ierr)
#endif 
      rdtx = d_one/dt
      do k = 1 , kz
        hydf(k) = gpcf/dsigma(k)
      end do
#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif 
        jm1 = j-1
#if defined(BAND) && (!defined(MPP1))
        if (jm1 == 0) jm1 = jx
#endif 
        do k = 1 , kz
          do i = 2 , iym1
            dumr = d_four/(sps2%ps(i,j)  +sps2%ps(i,jm1) +    &
                           sps2%ps(i-1,j)+sps2%ps(i-1,jm1))
            auxx(i,k,j) = atm2%u(i,k,j)*dumr
            avxx(i,k,j) = atm2%v(i,k,j)*dumr
          end do
        end do
      end do
!
#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif 
        do k = 1 , kz
          do i = 2 , iym1
            thvx(i,k,j) = thx3d(i,k,j)*(d_one+ep1*qvb3d(i,k,j))
            qcx(i,k,j)  = atm2%qc(i,k,j)/sps2%ps(i,j)
          end do
        end do
!
!.....density at surface is stored in rhox2d(i,j), at half levels in
!       rhohf(i,k,j).
!
        do k = 1 , kzm1
          do i = 2 , iym1
            dza(i,k,j) = za(i,k,j) - za(i,k+1,j)
            xps = (a(k)*sps2%ps(i,j)+r8pt)*d_1000
            ps2 = (a(k+1)*sps2%ps(i,j)+r8pt)*d_1000
            rhohf(i,k,j) = (ps2-xps)/(egrav*dza(i,k,j))
          end do
        end do
!
        do i = 2 , iym1
          govrth(i) = egrav/thx3d(i,kz,j)
        end do
!
! *********************************************************************
!
!     compute the vertical diffusion term:
!
        do k = 2 , kz
          do i = 2 , iym1
            rc(i,k) = 0.257D0*dzq(i,k,j)**0.175D0
          end do
        end do
!
!-----compute the diffusion coefficient:
!
!       blackadar scheme above boundary layer top
!
        do k = 2 , kz
          do i = 2 , iym1
            kzmax = 0.8D0*dza(i,k-1,j)*dzq(i,k,j)*rdtx
            ss = ((ubx3d(i,k-1,j)-ubx3d(i,k,j))*   &
                  (ubx3d(i,k-1,j)-ubx3d(i,k,j))+   &
                  (vbx3d(i,k-1,j)-vbx3d(i,k,j))*   &
                  (vbx3d(i,k-1,j)-vbx3d(i,k,j)))/  &
                  (dza(i,k-1,j)*dza(i,k-1,j)) + 1.0D-9
            ri = govrth(i)*(thvx(i,k-1,j)-thvx(i,k,j))/(ss*dza(i,k-1,j))
            if ( (ri-rc(i,k)) >= d_zero ) then
              kzm(i,k) = kzo
            else
              kzm(i,k) = kzo + dsqrt(ss)*(rc(i,k)-ri)*szkm/rc(i,k)
            end if
            kzm(i,k) = dmin1(kzm(i,k),kzmax)
          end do
        end do
!
! *********************************************************************
!
!       holtslag pbl
!
!       initialize bl diffusion coefficients and counter-gradient terms
!       with free atmosphere values and make specific humidity
!
        do k = 2 , kz
          do i = 2 , iym1
!           eddy diffusivities for momentum, heat and moisture
            kvm(i,k,j) = kzm(i,k)
            kvh(i,k,j) = kzm(i,k)
            kvq(i,k,j) = kzm(i,k)
!chem
            if ( ichem == 1 ) then
              kvc(i,k,j) = kzm(i,k)
            end if
!chem_
!           counter gradient terms for heat and moisture
            cgh(i,k,j) = d_zero
          end do
        end do
 
        do i = 2 , iym1
!         compute friction velocity
          idx = i
          jdx = j
#ifndef BAND
          idx = min0(idx,iym1)
          idxm1 = i - 1
          idxm1 = max0(idxm1,2)
          jdxm1 = j - 1
#ifdef MPP1
          if ( myid == nproc-1 ) jdx = min0(jdx,jendx)
          if ( myid == 0 ) jdxm1 = max0(jdxm1,2)
#else
          jdx = min0(jdx,jxm1)
          jdxm1 = max0(jdxm1,2)
#endif
#endif
          uflxsfx = sfsta%uvdrag(idx,jdx)*ubx3d(i,kz,j)
          vflxsfx = sfsta%uvdrag(idx,jdx)*vbx3d(i,kz,j)

          ustr(i,j) = dsqrt(dsqrt(uflxsfx*uflxsfx+vflxsfx*vflxsfx) / &
                             rhox2d(i,j))
 
!         convert surface fluxes to kinematic units
          xhfx(i,j) = sfsta%hfx(i,j)/(cpd*rhox2d(i,j))
          xqfx(i,j) = sfsta%qfx(i,j)/rhox2d(i,j)
!         compute virtual heat flux at surface
          hfxv(i,j) = xhfx(i,j) + mult*thx3d(i,kz,j)*xqfx(i,j)
        end do
!
!       estimate potential temperature at 10m via log temperature
!       profile in the surface layer (brutsaert, p. 63).
!       calculate mixing ratio at 10m by assuming a constant
!       value from the surface to the lowest model level.
!
 
        do i = 2 , iym1
          sh10 = qvb3d(i,kz,j)/(qvb3d(i,kz,j)+d_one)
!         th10(i,j) = ((thx3d(i,kz,j)+sts2%tg(i,j))*d_half)*(d_one+mult*sh10)
!         th10(i,j) = thvx(i,kz,j) + hfxv(i,j)/(vonkar*ustr(i,j)* &
!                     dlog(za(i,kz,j)*d_r10)
 
!         "virtual" potential temperature
          if ( hfxv(i,j) >= d_zero ) then
            th10(i,j) = thvx(i,kz,j)
          else
!           th10(i,j) =
!----       (0.25*thx3d(i,kz,j)+0.75*sts2%tg(i,j))*(d_one+mult*sh10) first
!           approximation for obhukov length
            oblen = -d_half*(thx3d(i,kz,j)+sts2%tg(i,j)) *  &
                    (d_one+mult*sh10)*ustr(i,j)**d_three /  &
                    (gvk*(hfxv(i,j)+dsign(1.0D-10,hfxv(i,j))))
            if ( oblen >= za(i,kz,j) ) then
              th10(i,j) = thvx(i,kz,j) + hfxv(i,j)/(vonkar*ustr(i,j))*  &
                 (dlog(za(i,kz,j)*d_r10)+d_five/oblen*(za(i,kz,j)-d_10))
            else if ( oblen < za(i,kz,j) .and. oblen > d_10 ) then
              th10(i,j) = thvx(i,kz,j) + hfxv(i,j)/(vonkar*ustr(i,j))*  &
                  (dlog(oblen*d_r10)+d_five/oblen*(oblen-d_10)+         &
                  6.0D0*dlog(za(i,kz,j)/oblen))
            else if ( oblen <= d_10 ) then
              th10(i,j) = thvx(i,kz,j) + hfxv(i,j)/(vonkar*ustr(i,j))   &
                          *6.0D0*dlog(za(i,kz,j)*d_r10)
            end if
            th10(i,j) = dmax1(th10(i,j),sts2%tg(i,j))
          end if
!gtb      th10(i,j) = dmin1(th10(i,j),sts2%tg(i,j))  ! gtb add to minimize
 
!         obklen compute obukhov length
          obklen(i,j) = -th10(i,j)*ustr(i,j)**d_three / &
                  (gvk*(hfxv(i,j)+dsign(1.0D-10,hfxv(i,j))))
        end do
!
!       compute diffusivities and counter gradient terms
!
      end do

      call blhnew

#ifdef MPP1
      do j = jbegin , jendx
#ifndef BAND
      if ( (myid /= nproc-1) .or. (myid == nproc-1 .and. j < jendx)) then
#endif
        do k = 2 , kz
          do i = 2 , iym1
            akzz1(i,k,j) = rhohf(i,k-1,j)*kvm(i,k,j)/dza(i,k-1,j)
          end do
        end do
        do k = 1 , kz
          do i = 2 , iym1
            akzz2(i,k,j) = hydf(k)/sps2%ps(i,j)
          end do
        end do
#ifndef BAND
      end if
#endif
      end do
      ii = 0
      do k = 1 , kz
        do i = 2 , iym1
          ii = ii + 1
          wksend(ii) = akzz1(i,k,jxp)
        end do
      end do
      do k = 1 , kz
        do i = 2 , iym1
          ii = ii + 1
          wksend(ii) = akzz2(i,k,jxp)
        end do
      end do
      call mpi_sendrecv(wksend,iym2*kz*2,mpi_real8,ieast,1, &
                        wkrecv,iym2*kz*2,mpi_real8,iwest,1, &
                        mpi_comm_world,mpi_status_ignore,ierr)
      ii = 0
      do k = 1 , kz
        do i = 2 , iym1
          ii = ii + 1
          akzz1(i,k,0) = wkrecv(ii)
        end do
      end do
      do k = 1 , kz
        do i = 2 , iym1
          ii = ii + 1
          akzz2(i,k,0) = wkrecv(ii)
        end do
      end do
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm2
#endif
        do k = 1 , kz
          do i = 2 , iym1
            if ( k > 1 ) then
              akzz1(i,k,j) = rhohf(i,k-1,j)*kvm(i,k,j)/dza(i,k-1,j)
            end if
            akzz2(i,k,j) = egrav/(sps2%ps(i,j)*d_1000)/dsigma(k)
          end do
        end do
      end do
#endif

#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif
        jm1 = j-1
#if defined(BAND) && (!defined(MPP1))
        if (jm1 == 0) jm1 = jx
#endif
 
!       calculate coefficients at dot points for u and v wind
 
#ifndef BAND
#ifdef MPP1
        if ( myid == 0 .and. j == 2 ) then
#else
        if ( j == 2 ) then
#endif
          do k = 2 , kz
            do i = 2 , iym1
              idx = i
              idx = min0(idx,iym2)
              idxm1 = i - 1
              idxm1 = max0(idxm1,2)
              betak(i,k) = d_half*(akzz1(idx,k,j)+akzz1(idxm1,k,j))
            end do
          end do
          do k = 1 , kz
            do i = 2 , iym1
              idx = i
              idx = min0(idx,iym2)
              idxm1 = i - 1
              idxm1 = max0(idxm1,2)
              alphak(i,k) = d_half*(akzz2(idx,k,j)+akzz2(idxm1,k,j))
            end do
          end do
#ifdef MPP1
        else if ( myid == nproc-1 .and. j == jendx ) then
#else
        else if ( j == jxm1 ) then
#endif
          do k = 2 , kz
            do i = 2 , iym1
              idx = i
              idx = min0(idx,iym2)
              idxm1 = i - 1
              idxm1 = max0(idxm1,2)
              betak(i,k) = d_half*(akzz1(idx,k,jm1)+akzz1(idxm1,k,jm1))
            end do
          end do
          do k = 1 , kz
            do i = 2 , iym1
              idx = i
              idx = min0(idx,iym2)
              idxm1 = i - 1
              idxm1 = max0(idxm1,2)
              alphak(i,k) = d_half*(akzz2(idx,k,jm1)+akzz2(idxm1,k,jm1))
            end do
          end do
        else
#endif
         do k = 2 , kz
           do i = 2 , iym1
             idx = i
             idx = min0(idx,iym2)
             idxm1 = i - 1
             idxm1 = max0(idxm1,2)
             betak(i,k) = (akzz1(idx,k,jm1)+akzz1(idxm1,k,jm1)+ &
                           akzz1(idx,k,j)+akzz1(idxm1,k,j))*d_rfour
           end do
         end do
         do k = 1 , kz
           do i = 2 , iym1
             idx = i
             idx = min0(idx,iym2)
             idxm1 = i - 1
             idxm1 = max0(idxm1,2)
             alphak(i,k) = (akzz2(idx,k,jm1)+akzz2(idxm1,k,jm1)+ &
                            akzz2(idx,k,j)+akzz2(idxm1,k,j))*d_rfour
           end do
         end do
#ifndef BAND
        end if
#endif
!
! **********************************************************************
!
!       start now procedure for implicit diffusion calculations
!       performed separately for wind (dot points)
!       and temperature and water vapor (cross points)
!       countergradient term is not included in the implicit diffusion
!       scheme its effect is included as in the old explicit scheme
!       calculations assume fluxes positive upward, so the sign in front
!       of uflxsf and vflxsf has been changed in the various terms
!
!       wind components
!
!       first compute coefficients of the tridiagonal matrix
!
! **********************************************************************
!
        do i = 2 , iym1
          coef1(i,1) = dt*alphak(i,1)*betak(i,2)
          coef2(i,1) = d_one + dt*alphak(i,1)*betak(i,2)
          coef3(i,1) = d_zero
          coefe(i,1) = coef1(i,1)/coef2(i,1)
          coeff1(i,1) = auxx(i,1,j)/coef2(i,1)
          coeff2(i,1) = avxx(i,1,j)/coef2(i,1)
        end do
 
        do k = 2 , kz - 1
          do i = 2 , iym1
            coef1(i,k) = dt*alphak(i,k)*betak(i,k+1)
            coef2(i,k) = d_one+dt*alphak(i,k)*(betak(i,k+1)+betak(i,k))
            coef3(i,k) = dt*alphak(i,k)*betak(i,k)
            coefe(i,k) = coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
            coeff1(i,k) = (auxx(i,k,j)+coef3(i,k)*coeff1(i,k-1))/       &
                          (coef2(i,k)-coef3(i,k)*coefe(i,k-1))
            coeff2(i,k) = (avxx(i,k,j)+coef3(i,k)*coeff2(i,k-1))/       &
                          (coef2(i,k)-coef3(i,k)*coefe(i,k-1))
          end do
        end do
 
        do i = 2 , iym1
          idx = i
          idx = min0(idx,iym1)
          idxm1 = i - 1
          idxm1 = max0(idxm1,2)

#ifdef BAND
          drgdot = (sfsta%uvdrag(idxm1,jm1)+sfsta%uvdrag(idxm1,j) + &
                    sfsta%uvdrag(idx,jm1)  +sfsta%uvdrag(idx,j))*d_rfour
#else
          jdx = j
          jdxm1 = j - 1
#ifdef MPP1
          if ( myid == nproc-1 ) jdx = min0(jdx,jendx)
          if ( myid == 0 ) jdxm1 = max0(jdxm1,2)
#else
          jdx = min0(jdx,jxm1)
          jdxm1 = max0(jdxm1,2)
#endif
          drgdot = (sfsta%uvdrag(idxm1,jdxm1)+sfsta%uvdrag(idxm1,jdx)+  &
                  sfsta%uvdrag(idx,jdxm1)+sfsta%uvdrag(idx,jdx))*d_rfour
#endif
          uflxsf = drgdot*auxx(i,kz,j)
          vflxsf = drgdot*avxx(i,kz,j)
          coef1(i,kz) = d_zero
          coef2(i,kz) = d_one + dt*alphak(i,kz)*betak(i,kz)
          coef3(i,kz) = dt*alphak(i,kz)*betak(i,kz)
          coefe(i,kz) = d_zero
          coeff1(i,kz) = (auxx(i,kz,j)-dt*alphak(i,kz)*uflxsf+          &
                          coef3(i,kz)*coeff1(i,kz-1))/                  &
                         (coef2(i,kz)-coef3(i,kz)*coefe(i,kz-1))
          coeff2(i,kz) = (avxx(i,kz,j)-dt*alphak(i,kz)*vflxsf+          &
                          coef3(i,kz)*coeff2(i,kz-1))/                  &
                         (coef2(i,kz)-coef3(i,kz)*coefe(i,kz-1))
 
        end do
!
!       all coefficients have been computed, predict field and put it in
!       temporary work space tpred
!
        do i = 2 , iym1
          tpred1(i,kz) = coeff1(i,kz)
          tpred2(i,kz) = coeff2(i,kz)
        end do
 
        do k = kz - 1 , 1 , -1
          do i = 2 , iym1
            tpred1(i,k) = coefe(i,k)*tpred1(i,k+1) + coeff1(i,k)
            tpred2(i,k) = coefe(i,k)*tpred2(i,k+1) + coeff2(i,k)
          end do
        end do
 
!
!       calculate tendency due to vertical diffusion using temporary
!       predicted field
!
        do k = 1 , kz
          do i = 2 , iym1
            dumr = (sps2%ps(i,j)  +sps2%ps(i,jm1) +  &
                    sps2%ps(i-1,j)+sps2%ps(i-1,jm1))*d_rfour
            aten%u(i,k,j) = aten%u(i,k,j) + &
                            (tpred1(i,k)-auxx(i,k,j))*rdtx*dumr
            aten%v(i,k,j) = aten%v(i,k,j) + &
                            (tpred2(i,k)-avxx(i,k,j))*rdtx*dumr
          end do
        end do
 
!
!       Common coefficients
!
        do k = 1 , kz
          do i = 2 , iym1
            alphak(i,k) = hydf(k)/sps2%ps(i,j)
          end do
        end do
!
!       temperature
!       calculate coefficients at cross points for temperature
!
        do k = 2 , kz
          do i = 2 , iym1
            betak(i,k) = rhohf(i,k-1,j)*kvh(i,k,j)/dza(i,k-1,j)
          end do
        end do
 
        do i = 2 , iym1
          coef1(i,1) = dt*alphak(i,1)*betak(i,2)
          coef2(i,1) = d_one + dt*alphak(i,1)*betak(i,2)
          coef3(i,1) = d_zero
          coefe(i,1) = coef1(i,1)/coef2(i,1)
          coeff1(i,1) = thx3d(i,1,j)/coef2(i,1)
        end do
 
        do k = 2 , kz - 1
          do i = 2 , iym1
            coef1(i,k) = dt*alphak(i,k)*betak(i,k+1)
            coef2(i,k) = d_one+dt*alphak(i,k)*(betak(i,k+1)+betak(i,k))
            coef3(i,k) = dt*alphak(i,k)*betak(i,k)
            coefe(i,k) = coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
            coeff1(i,k) = (thx3d(i,k,j)+coef3(i,k)*coeff1(i,k-1)) / &
                          (coef2(i,k)-coef3(i,k)*coefe(i,k-1))
          end do
        end do
 
        do i = 2 , iym1
          coef1(i,kz) = d_zero
          coef2(i,kz) = d_one + dt*alphak(i,kz)*betak(i,kz)
          coef3(i,kz) = dt*alphak(i,kz)*betak(i,kz)
          coefe(i,kz) = d_zero
          coeff1(i,kz) = (thx3d(i,kz,j) + &
                 dt*alphak(i,kz)*sfsta%hfx(i,j)*rcpd + &
                   coef3(i,kz)*coeff1(i,kz-1)) /       &
                   (coef2(i,kz)-coef3(i,kz)*coefe(i,kz-1))
        end do
 
!
!       all coefficients have been computed, predict field and put it in
!       temporary work space tpred
!
 
        do i = 2 , iym1
          tpred1(i,kz) = coeff1(i,kz)
        end do
 
        do k = kz - 1 , 1 , -1
          do i = 2 , iym1
            tpred1(i,k) = coefe(i,k)*tpred1(i,k+1) + coeff1(i,k)
          end do
        end do
 
!
!       calculate tendency due to vertical diffusion using temporary
!       predicted field
!
        do k = 1 , kz
          do i = 2 , iym1
            sf = atm2%t(i,k,j)/thx3d(i,k,j)
            difft(i,k,j) = difft(i,k,j) + &
                           (tpred1(i,k)-thx3d(i,k,j))*rdtx*sf
          end do
        end do
!
!       water vapor
!
 
!       calculate coefficients at cross points for water vapor
 
        do k = 2 , kz
          do i = 2 , iym1
            betak(i,k) = rhohf(i,k-1,j)*kvq(i,k,j)/dza(i,k-1,j)
          end do
        end do
 
        do i = 2 , iym1
          coef1(i,1) = dt*alphak(i,1)*betak(i,2)
          coef2(i,1) = d_one + dt*alphak(i,1)*betak(i,2)
          coef3(i,1) = d_zero
          coefe(i,1) = coef1(i,1)/coef2(i,1)
          coeff1(i,1) = qvb3d(i,1,j)/coef2(i,1)
        end do
 
        do k = 2 , kz - 1
          do i = 2 , iym1
            coef1(i,k) = dt*alphak(i,k)*betak(i,k+1)
            coef2(i,k) = d_one+dt*alphak(i,k)*(betak(i,k+1)+betak(i,k))
            coef3(i,k) = dt*alphak(i,k)*betak(i,k)
            coefe(i,k) = coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
            coeff1(i,k) = (qvb3d(i,k,j)+coef3(i,k)*coeff1(i,k-1)) / &
                           (coef2(i,k)-coef3(i,k)*coefe(i,k-1))
          end do
        end do
 
        do i = 2 , iym1
          coef1(i,kz) = d_zero
          coef2(i,kz) = d_one + dt*alphak(i,kz)*betak(i,kz)
          coef3(i,kz) = dt*alphak(i,kz)*betak(i,kz)
          coefe(i,kz) = d_zero
          coeff1(i,kz) = (qvb3d(i,kz,j) + &
                   dt*alphak(i,kz)*sfsta%qfx(i,j) + &
                   coef3(i,kz)*coeff1(i,kz-1)) /    &
                   (coef2(i,kz)-coef3(i,kz)*coefe(i,kz-1))
        end do
 
!
!       all coefficients have been computed, predict field and put it in
!       temporary work space tpred
 
        do i = 2 , iym1
          tpred1(i,kz) = coeff1(i,kz)
        end do
 
        do k = kz - 1 , 1 , -1
          do i = 2 , iym1
            tpred1(i,k) = coefe(i,k)*tpred1(i,k+1) + coeff1(i,k)
          end do
        end do
 
!
!       calculate tendency due to vertical diffusion using temporary
!       predicted field
!
        do k = 1 , kz
          do i = 2 , iym1
            diffq(i,k,j) = diffq(i,k,j) + &
                           (tpred1(i,k)-atm2%qv(i,k,j) / &
                           sps2%ps(i,j))*rdtx*sps2%ps(i,j)
          end do
        end do
 
!       calculate coefficients at cross points for cloud vater
 
        do k = 2 , kz
          do i = 2 , iym1
            betak(i,k) = rhohf(i,k-1,j)*kvq(i,k,j)/dza(i,k-1,j)
          end do
        end do
 
        do i = 2 , iym1
          coef1(i,1) = dt*alphak(i,1)*betak(i,2)
          coef2(i,1) = d_one + dt*alphak(i,1)*betak(i,2)
          coef3(i,1) = d_zero
          coefe(i,1) = coef1(i,1)/coef2(i,1)
          coeff1(i,1) = qcx(i,1,j)/coef2(i,1)
        end do
 
        do k = 2 , kz - 1
          do i = 2 , iym1
            coef1(i,k) = dt*alphak(i,k)*betak(i,k+1)
            coef2(i,k) = d_one+dt*alphak(i,k)*(betak(i,k+1)+betak(i,k))
            coef3(i,k) = dt*alphak(i,k)*betak(i,k)
            coefe(i,k) = coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
            coeff1(i,k) = (qcx(i,k,j)+coef3(i,k)*coeff1(i,k-1)) / &
                          (coef2(i,k)-coef3(i,k)*coefe(i,k-1))
          end do
        end do
 
        do i = 2 , iym1
          coef1(i,kz) = d_zero
          coef2(i,kz) = d_one + dt*alphak(i,kz)*betak(i,kz)
          coef3(i,kz) = dt*alphak(i,kz)*betak(i,kz)
          coefe(i,kz) = d_zero
          coeff1(i,kz) = (qcx(i,kz,j)+coef3(i,kz)*coeff1(i,kz-1)) / &
                         (coef2(i,kz)-coef3(i,kz)*coefe(i,kz-1))
        end do
 
!
!       all coefficients have been computed, predict field and put it in
!       temporary work space tpred
!
 
        do i = 2 , iym1
          tpred1(i,kz) = coeff1(i,kz)
        end do
 
        do k = kz - 1 , 1 , -1
          do i = 2 , iym1
            tpred1(i,k) = coefe(i,k)*tpred1(i,k+1) + coeff1(i,k)
          end do
        end do
 
!
!       calculate tendency due to vertical diffusion using temporary
!       predicted field
!
        do k = 1 , kz
          do i = 2 , iym1
            aten%qc(i,k,j) = aten%qc(i,k,j) +             &
                            (tpred1(i,k)-atm2%qc(i,k,j) / &
                            sps2%ps(i,j))*rdtx*sps2%ps(i,j)
          end do
        end do
 
!
! **********************************************************************
!
!       now add countergradient term to temperature and water vapor
!       equation
!trapuv
        do i = 2 , iym1
          ttnp(i,1) = d_zero
        end do
!trapuv_
        do k = 2 , kz
          do i = 2 , iym1
            sf = atm2%t(i,k,j)/(sps2%ps(i,j)*thx3d(i,k,j))
            ttnp(i,k) = sf*cpd*rhohf(i,k-1,j)*kvh(i,k,j)*cgh(i,k,j)
          end do
        end do
!
!-----compute the tendencies:
!
        do i = 2 , iym1
          difft(i,kz,j) = difft(i,kz,j) - hydf(kz)*ttnp(i,kz)/cpd
        end do
!
        do k = 1 , kzm1
          do i = 2 , iym1
            difft(i,k,j) = difft(i,k,j) + &
                           hydf(k)*(ttnp(i,k+1)-ttnp(i,k))/cpd
          end do
        end do
 
!
!chem2
        if ( ichem == 1 .and. ichdrdepo == 1 ) then
!
!         coef1, coef2, coef3 and coefe are the same as for water vapor
!         and cloud water so they do not need to be recalculated
 
!
!         recalculation of coef1,2,3  with tracer diffusivity kvc
 
          do k = 2 , kz
            do i = 2 , iym1
              betak(i,k) = rhohf(i,k-1,j)*kvc(i,k,j)/dza(i,k-1,j)
            end do
          end do
 
          do k = 2 , kz - 1
            do i = 2 , iym1
              coef1(i,k) = dt*alphak(i,k)*betak(i,k+1)
              coef2(i,k) = d_one+dt*alphak(i,k)*(betak(i,k+1)+betak(i,k))
              coef3(i,k) = dt*alphak(i,k)*betak(i,k)
            end do
          end do
 
          do i = 2 , iym1
            coef1(i,1) = dt*alphak(i,1)*betak(i,2)
            coef2(i,1) = d_one + dt*alphak(i,1)*betak(i,2)
            coef3(i,1) = d_zero
            coef1(i,kz) = d_zero
            coef2(i,kz) = d_one + dt*alphak(i,kz)*betak(i,kz)
            coef3(i,kz) = dt*alphak(i,kz)*betak(i,kz)
          end do
!
!         set the Vd for in case of prescribed deposition velocities
!
          do itr = 1 , ntr
            do i = 2 , iym1
#ifdef CLM
              if ( ocld2d(1,i,j) == 0 ) then
#else
              if ( veg2d(i,j) == 14 .or. veg2d(i,j) == 15 ) then
#endif
                vdep(i,itr) = chtrdpv(itr,2)
              else
                vdep(i,itr) = chtrdpv(itr,1)
              end if
!             provisoire test de la routine chdrydep pour les dust
 
              if ( chtrname(itr) == 'DUST' ) vdep(i,itr) = d_zero
            end do
          end do
!
          do itr = 1 , ntr
!
            do k = 1 , kz
              do i = 2 , iym1
                chix(i,k) = chib(i,k,j,itr)/sps2%ps(i,j)
              end do
            end do
!
            do i = 2 , iym1
              coefe(i,1) = coef1(i,1)/coef2(i,1)
              coeff1(i,1) = chix(i,1)/coef2(i,1)
            end do
!
            do k = 2 , kz - 1
              do i = 2 , iym1
                coefe(i,k) = coef1(i,k)                                 &
                             /(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
                coeff1(i,k) = (chix(i,k)+coef3(i,k)*coeff1(i,k-1))      &
                              /(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
              end do
            end do
 
            do i = 2 , iym1
              coefe(i,kz) = d_zero
 
!             add dry deposition option1
              coeff1(i,kz) = (chix(i,kz)-dt*alphak(i,kz)*chix(i,kz)     &
                             *vdep(i,itr)*rhox2d(i,j)+coef3(i,kz)       &
                             *coeff1(i,kz-1))                           &
                             /(coef2(i,kz)-coef3(i,kz)*coefe(i,kz-1))
            end do
!
!           all coefficients have been computed, predict field and put
!           it in temporary work space tpred1
!
            do i = 2 , iym1
              tpred1(i,kz) = coeff1(i,kz)
            end do
!
            do k = kz - 1 , 1 , -1
              do i = 2 , iym1
                tpred1(i,k) = coefe(i,k)*tpred1(i,k+1) + coeff1(i,k)
              end do
            end do
!
!
!           calculate tendency due to vertical diffusion using temporary
!           predicted field
!           Dry deposition option 1 is included
 
!
            do k = 1 , kz
              do i = 2 , iym1
!qian           chiten(i,k,j,itr)=chiten(i,k,j,itr)
!CGAFFE         TEST diffusion/10
                chiten(i,k,j,itr) = chiten(i,k,j,itr)                   &
                                    + (tpred1(i,k)-chix(i,k))           &
                                    *rdtx*sps2%ps(i,j)
!               chiten(i,k,j,itr)=chiten(i,k,j,itr)+0.1 *(tpred1(i,k)-
!               1  chix(i,k))*rdtx *sps2%ps(i,j)
 
              end do
            end do
            do i = 2 , iym1
 
              if ( chtrname(itr) /= 'DUST' ) &
                remdrd(i,j,itr) = remdrd(i,j,itr) + chix(i,kz)* &
                    vdep(i,itr)*sps2%ps(i,j)*dt*d_half*rhox2d(i,j)* &
                    hydf(kz)/sps2%ps(i,j)
 
            end do
          end do
        end if
!chem2_
 
      end do
      end subroutine holtbl
!
! ------------------------------------------------------------
! this routine computes the boundary layer eddy diffusivities
! for momentum, heat and moisture and the counter-gradient
! terms for heat and moisture.
!
! reference : holtslag, de bruijn and pan - mwr - 8/90
!
! input arguments :  j       longitudinal position index
!                    ubx3d   u wind component
!                    vbx3d   v wind component
!                    thx3d   potential temperature
!                    thvx    virtual potential temperature
!                    za      height of half sigma levels
!                    f       coriolis parameter
!                    shum    specific humidity
!                    xhfx    sensible heat flux
!                    xqfx    sfc kinematic moisture flux
!                    th10    virt. pot. temp. at 10m
!                    hfxv    surface virtual heat flux
!                    obklen  monin obukov length
!                    ustr    friction velocity
!
! input/output
! arguments :        therm   thermal temperature excess
!
! output arguments : cgh     counter-gradient term for heat
!                    cgq     counter-gradient term for moisture
!                    kvm     eddy diffusivity for momentum
!                    kvh     eddy diffusivity for heat
!                    kvq     eddy diffusivity for moisture
!                    zpbl     boundary layer height
! ------------------------------------------------------------
!
      subroutine blhnew
!
      implicit none
!
      real(8) :: fak1 , fak2 , fht , xfmt , pblk , pblk1 , pblk2 ,  &
                 pfcor , phpblm , pr , therm2 , tkv , tlv , vv ,    &
                 vvl , wsc , z , zh , zl , zm , zp , zzh , zzhnew , &
                 zzhnew2
      integer :: i , j , k , k2
!
      pblk2 = d_zero
      zzhnew2 = d_zero
!
#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif
!       ****note: kt, max no. of pbl levels, calculated in param
!       ******   compute richardson number
!
        therm(:) = d_zero
 
        do k = kz , kt , -1
          do i = 2 , iym1
            vv = ubx3d(i,k,j)*ubx3d(i,k,j) + vbx3d(i,k,j)*vbx3d(i,k,j)
            ri(i,k) = egrav*(thvx(i,k,j)-th10(i,j))*za(i,k,j)/ &
                            (th10(i,j)*vv)
          end do
        end do
!       ******   first, set bl height to height of lowest model level
        do i = 2 , iym1
          sfsta%zpbl(i,j) = za(i,kz,j)
        end do
!       ******   looking for bl top
        do k = kz , kt + 1 , -1
          k2 = k - 1
          do i = 2 , iym1
!     ******   bl height lies between this level and the last
!     ******   use linear interp. of rich. no. to height of ri=ricr
            if ( (ri(i,k) < ricr) .and. (ri(i,k2) >= ricr) )       &
              sfsta%zpbl(i,j) = za(i,k,j) + (za(i,k2,j)-za(i,k,j)) &
                  *((ricr-ri(i,k))/(ri(i,k2)-ri(i,k)))
          end do
        end do
 
        do i = 2 , iym1
!     ******   set bl top to highest allowable model layer
          if ( ri(i,kt) < ricr ) sfsta%zpbl(i,j) = za(i,kt,j)
        end do
 
!       ******   recompute richardson no. at lowest model level
        do i = 2 , iym1
          if ( hfxv(i,j) > d_zero ) then
!           ******   estimate of convective velocity scale
            xfmt = (d_one-(binm*sfsta%zpbl(i,j)/obklen(i,j)))**onet
            wsc = ustr(i,j)*xfmt
!           ******   thermal temperature excess
            therm(i) = (xhfx(i,j)+mult*thx3d(i,kz,j)*xqfx(i,j))*fak/wsc
            vvl = ubx3d(i,kz,j)*ubx3d(i,kz,j) + &
                  vbx3d(i,kz,j)*vbx3d(i,kz,j)
            ri(i,kz) = -egrav*therm(i)*za(i,kz,j)/(th10(i,j)*vvl)
          end if
        end do
 
!       ******   recompute richardson no. at other model levels
        do k = kz - 1 , kt , -1
          do i = 2 , iym1
            if ( hfxv(i,j) > d_zero ) then
              tlv = th10(i,j) + therm(i)
              tkv = thx3d(i,k,j) * &
                    (d_one+mult*(qvb3d(i,k,j)/(qvb3d(i,k,j)+d_one)))
              vv = ubx3d(i,k,j)*ubx3d(i,k,j)+vbx3d(i,k,j)*vbx3d(i,k,j)
              ri(i,k) = egrav*(tkv-tlv)*za(i,k,j)/(th10(i,j)*vv)
            end if
          end do
        end do
 
!       ******   improve estimate of bl height under convective
!       conditions ******   using convective temperature excess (therm)
        do k = kz , kt + 1 , -1
          k2 = k - 1
          do i = 2 , iym1
            if ( hfxv(i,j) > d_zero ) then
!     ******   bl height lies between this level and the last
!     ******   use linear interp. of rich. no. to height of ri=ricr
              if ( (ri(i,k) < ricr) .and. (ri(i,k2) >= ricr) ) then
                sfsta%zpbl(i,j) = za(i,k,j) + (za(i,k2,j)-za(i,k,j)) &
                               *((ricr-ri(i,k))/(ri(i,k2)-ri(i,k)))
              end if
            end if
          end do
        end do
 
        do i = 2 , iym1
          if ( hfxv(i,j) > d_zero ) then
!     ******   set bl top to highest allowable model layer
            if ( ri(i,kt) < ricr ) then
              sfsta%zpbl(i,j) = za(i,kt,j)
            end if
          end if
        end do
 
!       ******   limit bl height to be at least mech. mixing depth
        do i = 2 , iym1
!         ******   limit coriolis parameter to value at 10 deg. latitude
          pfcor = dmax1(dabs(mddom%f(i,j)),2.546D-5)
!         ******   compute mechanical mixing depth,
!         ******   set to lowest model level if lower
          phpblm = 0.07D0*ustr(i,j)/pfcor
          phpblm = dmax1(phpblm,za(i,kz,j))
          sfsta%zpbl(i,j) = dmax1(sfsta%zpbl(i,j),phpblm)
        end do
 
        do k = kz , kt + 1 , -1
          k2 = k - 1
          do i = 2 , iym1
            pblk = d_zero
            zm = za(i,k,j)
            zp = za(i,k2,j)
            if ( zm < sfsta%zpbl(i,j) ) then
              zp = dmin1(zp,sfsta%zpbl(i,j))
              z = (zm+zp)*d_half
              zh = z/sfsta%zpbl(i,j)
              zl = z/obklen(i,j)
              if ( zh <= d_one ) then
                zzh = d_one - zh
                zzh = zzh**pink
!xexp4          zzhnew = sfsta%zpbl(i,j)*(d_one-zh)*zh**1.5
!xexp5          zzhnew = 0.5*sfsta%zpbl(i,j)*(d_one-zh)*zh**1.5
!xexp6          zzhnew = d_one - zh
!xexp7          zzhnew =0.5* (d_one - zh)
!Sara
!               zzhnew =0.25* (d_one - zh)
!               zzhnew =0.75* (d_one - zh)
!Sara_
                zzhnew = (d_one-zh)*d_rfour
!xexp10         zzhnew =zh * (d_one - zh)**2
!chem
                if ( ichem == 1 ) then
                  zzhnew2 = (d_one-zh)**d_two
                end if
!chem_
              else
                zzh = d_zero
                zzhnew = d_zero
!chem
                zzhnew2 = d_zero
!chem_
              end if
              fak1 = ustr(i,j)*sfsta%zpbl(i,j)*vonkar
              if ( hfxv(i,j) <= d_zero ) then
!**             stable and neutral conditions
!**             igroup = 1
 
!**             prevent pblk from becoming too small in very stable
!               conditions
                if ( zl <= d_one ) then
                  pblk = fak1*zh*zzh/(d_one+betas*zl)
!xexp5            pblk1 = vonkar * ustr(i,j) / (d_one+betas*zl) * zzhnew
                  pblk1 = fak1*zh*zzhnew/(d_one+betas*zl)
!chem
                  if ( ichem == 1 ) then
                    pblk2 = fak1*zh*zzhnew2/(d_one+betas*zl)
                  end if
!chem_
                else
                  pblk = fak1*zh*zzh/(betas+zl)
!xexp5            pblk1 = vonkar * ustr(i,j) / (betas+zl) * zzhnew
                  pblk1 = fak1*zh*zzhnew/(betas+zl)
!chem
                  if ( ichem == 1 ) then
                    pblk2 = fak1*zh*zzhnew2/(betas+zl)
                  end if
!chem_
                end if
!**             compute eddy diffusivities
                kvm(i,k,j) = dmax1(pblk,kzo)
                kvh(i,k,j) = kvm(i,k,j)
                kvq(i,k,j) = dmax1(pblk1,kzo)
! Erika put k=0 in very stable conditions
                if ( zl <= 0.1D0 ) then
                  kvm(i,k,j) = d_zero
                  kvh(i,k,j) = kvm(i,k,j)*d_zero
                  kvq(i,k,j) = d_zero
                end if
! Erika put k=0 in very stable conditions

!chem
                if ( ichem == 1 ) then
                  kvc(i,k,j) = dmax1(pblk2,kzo)
                end if
!chem_
!**             compute counter-gradient term
                cgh(i,k,j) = d_zero
!               cgq(i,k) = d_zero
              else
!**             unstable conditions
 
!**             compute counter gradient term
                if ( zh >= sffrac ) then
!**               igroup = 2
                  xfmt = (d_one-binm*sfsta%zpbl(i,j)/obklen(i,j))**onet
                  fht = dsqrt(d_one-binh*sfsta%zpbl(i,j)/obklen(i,j))
                  wsc = ustr(i,j)*xfmt
                  pr = (xfmt/fht) + ccon
                  fak2 = wsc*sfsta%zpbl(i,j)*vonkar
                  pblk = fak2*zh*zzh
!xexp5            pblk1 = vonkar * wsc * zzhnew
                  pblk1 = fak2*zh*zzhnew
!chem
                  if ( ichem == 1 ) then
                    pblk2 = fak2*zh*zzhnew2
                  end if
!chem_
                  therm2 = fak/(sfsta%zpbl(i,j)*wsc)
                  cgh(i,k,j) = hfxv(i,j)*therm2
!                 cgq(i,k) = xqfx(i,j)*therm2
                else
!**               igroup = 3
                  pblk = fak1*zh*zzh*(d_one-betam*zl)**onet
!xexp5            pblk1 = vonkar * ustr(i,j) * zzhnew *
!                 (d_one-betam*zl)**onet
                  pblk1 = fak1*zh*zzhnew*(d_one-betam*zl)**onet
!chem
                  if ( ichem == 1 ) then
                    pblk2 = fak1*zh*zzhnew2*(d_one-betam*zl)**onet
                  end if
!chem_
                  pr = ((d_one-betam*zl)**onet)/dsqrt(d_one-betah*zl)
                  cgh(i,k,j) = d_zero
!                 cgq(i,k) = d_zero
                end if
 
!**             compute eddy diffusivities
                kvm(i,k,j) = dmax1(pblk,kzo)
                kvh(i,k,j) = dmax1((pblk/pr),kzo)
!               kvq(i,k,j) = kvh(i,k,j)
                kvq(i,k,j) = dmax1(pblk1,kzo)
!chem
                if ( ichem == 1 ) then
                  kvc(i,k,j) = dmax1(pblk2,kzo)
                end if
!chem_
              end if
            end if
          end do
        end do
      end do
 
      end subroutine blhnew
!
      end module mod_holtbl
