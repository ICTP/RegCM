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
 
      subroutine holtbl

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      use mod_regcm_param
      use mod_param1
      use mod_param2
      use mod_param3 , only : dsigma , ptop , a
      use mod_main
      use mod_mainchem
      use mod_pbldim
      use mod_cvaria
      use mod_pmoist
      use mod_bats , only : veg2d
      use mod_slice
      use mod_trachem
      use mod_constants , only : gti , vonkar , cpd , rcpd , ep1
#ifdef MPP1
      use mod_mppio
      use mpi
#endif
#ifdef DIAG
      use mod_diagnosis
#endif
      use mod_blh_tmp
      implicit none
!
! Local variables
!
      real(8) , dimension(ix,kx) :: alphak , betak , chix , coef1 ,     &
                                  & coef2 , coef3 , coefe , coeff1 ,    &
                                  & coeff2 , tpred1 , tpred2
      real(8) :: drgdot , dumr , kzmax , kzo , oblen , xps , ps2 , ri , &
               & sf , sh10 , ss , szkm , tvcon , uflxsf , uflxsfx ,     &
               & vflxsf , vflxsfx
      real(8) , dimension(ixm1) :: govrth
      integer :: i , idx , idxm1 , itr , j , jdx , jdxm1 , k
      real(8) , dimension(ixm1,kx) :: kzm , rc , ttnp
      real(8) , dimension(ixm1,ntr) :: vdep
#ifdef MPP1
      integer :: ierr , ii
      real(8) , dimension(ixm1,kx,jxp) :: auxx , avxx , dza , qcx
      real(8) , dimension(ixm1,kx,0:jxp+1) :: akxx1 , akxx2
      real(8) , dimension(2*(ixm2)*kx) :: wkrecv , wksend
      real(8) , dimension(ix,kx,jxp) :: rhohf
#else
      real(8) , dimension(ixm1,kx,jxm1) :: auxx , avxx , dza , qcx
      real(8) , dimension(ixm1,kx,0:jxm1) :: akxx1 , akxx2
      real(8) , dimension(ix,kx,jxm1) :: rhohf
#endif
!
      data kzo/1./
      data szkm/1600./
!
! *********************************************************************
!
!     diagnostic on total evaporation
!
#ifdef DIAG
#ifdef MPP1
      call mpi_gather(qfx(1,1),ix*jxp,mpi_real8,qfx_io(1,1),            &
                    & ix*jxp,mpi_real8,0,mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do j = 2 , jxm2
          do i = 2 , ixm2
            tqeva = tqeva + qfx_io(i,j)*dx*dx*dtmin*60.
          end do
        end do
      end if
      call mpi_bcast(tqeva,1,mpi_real8,0,mpi_comm_world,ierr)
#else
      do j = 2 , jxm2
        do i = 2 , ixm2
          tqeva = tqeva + qfx(i,j)*dx*dx*dtmin*60.
        end do
      end do
#endif
#endif

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
      call mpi_sendrecv(psb(1,jxp),ix,mpi_real8,ieast,1,                &
                      & psb(1,0),ix,mpi_real8,iwest,1,                  &
                      & mpi_comm_world,mpi_status_ignore,ierr)
      call mpi_sendrecv(uvdrag(1,jxp),ix,mpi_real8,ieast,1,             &
                      & uvdrag(1,0),ix,mpi_real8,iwest,1,               &
                      & mpi_comm_world,mpi_status_ignore,ierr)
      do j = jbegin , jendx
        do k = 1 , kx
          do i = 2 , ixm1
            dumr = 4./(psb(i,j)+psb(i,j-1)+psb(i-1,j)+psb(i-1,j-1))
            auxx(i,k,j) = ub(i,k,j)*dumr
            avxx(i,k,j) = vb(i,k,j)*dumr
          end do
        end do
#else
      do j = 2 , jxm1
        do k = 1 , kx
          do i = 2 , ixm1
            dumr = 4./(psb(i,j)+psb(i,j-1)+psb(i-1,j)+psb(i-1,j-1))
            auxx(i,k,j) = ub(i,k,j)*dumr
            avxx(i,k,j) = vb(i,k,j)*dumr
          end do
        end do
#endif 
!
        do k = 1 , kx
          do i = 2 , ixm1
            tvcon = (1.+ep1*qvb3d(i,k,j))
            thvx(i,k,j) = thx3d(i,k,j)*tvcon
          end do
        end do
!
        do k = 1 , kx
          do i = 2 , ixm1
            qcx(i,k,j) = qcb(i,k,j)/psb(i,j)
          end do
        end do
!
!.....density at surface is stored in rhox2d(i,j), at half levels in
!       rhohf(i,k,j).
!
        do k = 1 , kxm1
          do i = 2 , ixm1
            dza(i,k,j) = za(i,k,j) - za(i,k+1,j)
            xps = (a(k)*psb(i,j)+ptop)*1000.
            ps2 = (a(k+1)*psb(i,j)+ptop)*1000.
            rhohf(i,k,j) = (ps2-xps)/(gti*dza(i,k,j))
          end do
        end do
!
        do i = 2 , ixm1
          govrth(i) = gti/thx3d(i,kx,j)
        end do
!
! *********************************************************************
!
!-----compute the vertical diffusion term:
!
        do k = 2 , kx
          do i = 2 , ixm1
            rc(i,k) = 0.257*dzq(i,k,j)**0.175
          end do
        end do
!
!-----compute the diffusion coefficient:
!
!       blackadar scheme above boundary layer top
!
        do k = 2 , kx
          do i = 2 , ixm1
            kzmax = 0.8*dza(i,k-1,j)*dzq(i,k,j)/dt
            ss = ((ubx3d(i,k-1,j)-ubx3d(i,k,j))                         &
               & *(ubx3d(i,k-1,j)-ubx3d(i,k,j))                         &
               & +(vbx3d(i,k-1,j)-vbx3d(i,k,j))                         &
               & *(vbx3d(i,k-1,j)-vbx3d(i,k,j)))                        &
               & /(dza(i,k-1,j)*dza(i,k-1,j)) + 1.E-9
            ri = govrth(i)*(thvx(i,k-1,j)-thvx(i,k,j))/(ss*dza(i,k-1,j))
            if ( (ri-rc(i,k)).ge.0. ) then
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
        do k = 2 , kx
          do i = 2 , ixm1
!           eddy diffusivities for momentum, heat and moisture
            kvm(i,k,j) = kzm(i,k)
            kvh(i,k,j) = kzm(i,k)
            kvq(i,k,j) = kzm(i,k)
!chem
            if ( ichem.eq.1 ) kvc(i,k,j) = kzm(i,k)
!chem_
!           counter gradient terms for heat and moisture
            cgh(i,k,j) = 0.0
          end do
        end do
 
        do i = 2 , ixm1
!         compute friction velocity
          idx = i
          idx = min0(idx,ixm1)
          idxm1 = i - 1
          idxm1 = max0(idxm1,2)
          jdx = j
          jdxm1 = j - 1
#ifdef MPP1
          if ( myid.eq.nproc-1 ) jdx = min0(jdx,jendx)
          if ( myid.eq.0 ) jdxm1 = max0(jdxm1,2)
#else
          jdx = min0(jdx,jxm1)
          jdxm1 = max0(jdxm1,2)
#endif
          uflxsfx = uvdrag(idx,jdx)*ubx3d(i,kx,j)
          vflxsfx = uvdrag(idx,jdx)*vbx3d(i,kx,j)
          ustr(i,j) = dsqrt(dsqrt(uflxsfx*uflxsfx+vflxsfx*vflxsfx)      &
                    & /rhox2d(i,j))
 
!         convert surface fluxes to kinematic units
          xhfx(i,j) = hfx(i,j)/(cpd*rhox2d(i,j))
          xqfx(i,j) = qfx(i,j)/rhox2d(i,j)
!         compute virtual heat flux at surface
          hfxv(i,j) = xhfx(i,j) + 0.61*thx3d(i,kx,j)*xqfx(i,j)
        end do
!
!       estimate potential temperature at 10m via log temperature
!       profile in the surface layer (brutsaert, p. 63).
!       calculate mixing ratio at 10m by assuming a constant
!       value from the surface to the lowest model level.
!
 
        do i = 2 , ixm1
          sh10 = qvb3d(i,kx,j)/(qvb3d(i,kx,j)+1)
!         th10(i,j) = ((thx3d(i,kx,j)+tgb(i,j))/2.0)*(1.0+0.61*sh10)
!         th10(i,j) = thvx(i,kx,j) + hfxv(i,j)/(vonkar*ustr(i,j))
!         1            *dlog(za(i,kx,j)/10.)
 
!         "virtual" potential temperature
          if ( hfxv(i,j).ge.0. ) then
            th10(i,j) = thvx(i,kx,j)
          else
!           th10(i,j) =
!----       (0.25*thx3d(i,kx,j)+0.75*tgb(i,j))*(1.0+0.61*sh10) first
!           approximation for obhukov length
            oblen = -0.5*(thx3d(i,kx,j)+tgb(i,j))*(1.0+0.61*sh10)       &
                  & *ustr(i,j)                                          &
                  & **3/(gti*vonkar*(hfxv(i,j)+dsign(1.D-10,hfxv(i,j))))
            if ( oblen.ge.za(i,kx,j) ) then
              th10(i,j) = thvx(i,kx,j) + hfxv(i,j)/(vonkar*ustr(i,j))   &
                        & *(dlog(za(i,kx,j)/10.)                        &
                        & +5./oblen*(za(i,kx,j)-10.))
            else if ( oblen.lt.za(i,kx,j) .and. oblen.gt.10. ) then
              th10(i,j) = thvx(i,kx,j) + hfxv(i,j)/(vonkar*ustr(i,j))   &
                        & *(dlog(oblen/10.)+5./oblen*(oblen-10.)        &
                        & +6*dlog(za(i,kx,j)/oblen))
            else if ( oblen.le.10. ) then
              th10(i,j) = thvx(i,kx,j) + hfxv(i,j)/(vonkar*ustr(i,j))   &
                        & *6*dlog(za(i,kx,j)/10.)
            else
            end if
            th10(i,j) = dmax1(th10(i,j),tgb(i,j))
          end if
!gtb      th10(i,j) = dmin1(th10(i,j),tgb(i,j))  ! gtb add to minimize
 
!         obklen compute obukhov length
          obklen(i,j) = -th10(i,j)*ustr(i,j)                            &
                      & **3/(gti*vonkar*(hfxv(i,j)+dsign(1.D-10,        &
                      & hfxv(i,j))))
        end do
!
!       compute diffusivities and counter gradient terms
!
      end do
      call blhnew
#ifdef MPP1
      do j = jbegin , jendx
        if ( (myid.ne.nproc-1) .or. (myid.eq.nproc-1 .and. j.lt.jendx) )&
           & then
          do k = 1 , kx
            do i = 2 , ixm1
              if ( k.gt.1 ) akxx1(i,k,j) = rhohf(i,k-1,j)*kvm(i,k,j)    &
                 & /dza(i,k-1,j)
              akxx2(i,k,j) = gti/(psb(i,j)*1000.)/dsigma(k)
            end do
          end do
        end if
      end do
      ii = 0
      do k = 1 , kx
        do i = 2 , ixm1
          ii = ii + 1
          wksend(ii) = akxx1(i,k,jxp)
        end do
      end do
      do k = 1 , kx
        do i = 2 , ixm1
          ii = ii + 1
          wksend(ii) = akxx2(i,k,jxp)
        end do
      end do
      call mpi_sendrecv(wksend(1),(ixm2)*kx*2,mpi_real8,                &
                      & ieast,1,wkrecv(1),(ixm2)*kx*2,                  &
                      & mpi_real8,iwest,1,mpi_comm_world,               &
                      & mpi_status_ignore,ierr)
      ii = 0
      do k = 1 , kx
        do i = 2 , ixm1
          ii = ii + 1
          akxx1(i,k,0) = wkrecv(ii)
        end do
      end do
      do k = 1 , kx
        do i = 2 , ixm1
          ii = ii + 1
          akxx2(i,k,0) = wkrecv(ii)
        end do
      end do
#else
      do j = 2 , jxm2
        do k = 1 , kx
          do i = 2 , ixm1
            if ( k.gt.1 ) akxx1(i,k,j) = rhohf(i,k-1,j)*kvm(i,k,j)      &
               & /dza(i,k-1,j)
            akxx2(i,k,j) = gti/(psb(i,j)*1000.)/dsigma(k)
          end do
        end do
      end do
#endif

#ifdef MPP1
      do j = jbegin , jendx
#else
      do j = 2 , jxm1
#endif
 
!       calculate coefficients at dot points for u and v wind
 
#ifdef MPP1
        if ( myid.eq.0 .and. j.eq.2 ) then
#else
        if ( j.eq.2 ) then
#endif
          do k = 1 , kx
            do i = 2 , ixm1
              idx = i
              idx = min0(idx,ixm2)
              idxm1 = i - 1
              idxm1 = max0(idxm1,2)
              if ( k.gt.1 ) betak(i,k)                                  &
                 & = 0.5*(akxx1(idx,k,j)+akxx1(idxm1,k,j))
              alphak(i,k) = 0.5*(akxx2(idx,k,j)+akxx2(idxm1,k,j))
            end do
          end do
#ifdef MPP1
        else if ( myid.eq.nproc-1 .and. j.eq.jendx ) then
#else
        else if ( j.eq.jxm1 ) then
#endif
          do k = 1 , kx
            do i = 2 , ixm1
              idx = i
              idx = min0(idx,ixm2)
              idxm1 = i - 1
              idxm1 = max0(idxm1,2)
              if ( k.gt.1 ) betak(i,k)                                  &
                 & = 0.5*(akxx1(idx,k,j-1)+akxx1(idxm1,k,j-1))
              alphak(i,k) = 0.5*(akxx2(idx,k,j-1)+akxx2(idxm1,k,j-1))
            end do
          end do
        else
          do k = 1 , kx
            do i = 2 , ixm1
              idx = i
              idx = min0(idx,ixm2)
              idxm1 = i - 1
              idxm1 = max0(idxm1,2)
              if ( k.gt.1 ) betak(i,k)                                  &
                 & = 0.25*(akxx1(idx,k,j-1)+akxx1(idxm1,k,j-1)          &
                 & +akxx1(idx,k,j)+akxx1(idxm1,k,j))
              alphak(i,k) = 0.25*(akxx2(idx,k,j-1)+akxx2(idxm1,k,j-1)   &
                          & +akxx2(idx,k,j)+akxx2(idxm1,k,j))
            end do
          end do
        end if
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
 
!
!       wind components
!
 
!
!       first compute coefficients of the tridiagonal matrix
!
 
        do k = 2 , kx - 1
          do i = 2 , ixm1
            coef1(i,k) = dt*alphak(i,k)*betak(i,k+1)
            coef2(i,k) = 1. + dt*alphak(i,k)*(betak(i,k+1)+betak(i,k))
            coef3(i,k) = dt*alphak(i,k)*betak(i,k)
          end do
        end do
 
        do i = 2 , ixm1
          coef1(i,1) = dt*alphak(i,1)*betak(i,2)
          coef2(i,1) = 1. + dt*alphak(i,1)*betak(i,2)
          coef3(i,1) = 0.
          coef1(i,kx) = 0.
          coef2(i,kx) = 1. + dt*alphak(i,kx)*betak(i,kx)
          coef3(i,kx) = dt*alphak(i,kx)*betak(i,kx)
        end do
 
        do i = 2 , ixm1
          coefe(i,1) = coef1(i,1)/coef2(i,1)
          coeff1(i,1) = auxx(i,1,j)/coef2(i,1)
          coeff2(i,1) = avxx(i,1,j)/coef2(i,1)
        end do
 
        do k = 2 , kx - 1
          do i = 2 , ixm1
            coefe(i,k) = coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
            coeff1(i,k) = (auxx(i,k,j)+coef3(i,k)*coeff1(i,k-1))        &
                        & /(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
            coeff2(i,k) = (avxx(i,k,j)+coef3(i,k)*coeff2(i,k-1))        &
                        & /(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
          end do
        end do
 
        do i = 2 , ixm1
          idx = i
          idx = min0(idx,ixm1)
          idxm1 = i - 1
          idxm1 = max0(idxm1,2)
          jdx = j
          jdxm1 = j - 1
#ifdef MPP1
          if ( myid.eq.nproc-1 ) jdx = min0(jdx,jendx)
          if ( myid.eq.0 ) jdxm1 = max0(jdxm1,2)
#else
          jdx = min0(jdx,jxm1)
          jdxm1 = max0(jdxm1,2)
#endif
          drgdot = 0.25*(uvdrag(idxm1,jdxm1)+uvdrag(idxm1,jdx)          &
                 & +uvdrag(idx,jdxm1)+uvdrag(idx,jdx))
          uflxsf = drgdot*auxx(i,kx,j)
          vflxsf = drgdot*avxx(i,kx,j)
 
          coefe(i,kx) = 0.
          coeff1(i,kx) = (auxx(i,kx,j)-dt*alphak(i,kx)*uflxsf+coef3(i,kx&
                       & )*coeff1(i,kx-1))                              &
                       & /(coef2(i,kx)-coef3(i,kx)*coefe(i,kx-1))
          coeff2(i,kx) = (avxx(i,kx,j)-dt*alphak(i,kx)*vflxsf+coef3(i,kx&
                       & )*coeff2(i,kx-1))                              &
                       & /(coef2(i,kx)-coef3(i,kx)*coefe(i,kx-1))
 
        end do
!
!       all coefficients have been computed, predict field and put it in
!       temporary work space tpred
!
        do i = 2 , ixm1
          tpred1(i,kx) = coeff1(i,kx)
          tpred2(i,kx) = coeff2(i,kx)
        end do
 
        do k = kx - 1 , 1 , -1
          do i = 2 , ixm1
            tpred1(i,k) = coefe(i,k)*tpred1(i,k+1) + coeff1(i,k)
            tpred2(i,k) = coefe(i,k)*tpred2(i,k+1) + coeff2(i,k)
          end do
        end do
 
!
!       calculate tendency due to vertical diffusion using temporary
!       predicted field
!
        do k = 1 , kx
          do i = 2 , ixm1
            dumr = 0.25*(psb(i,j)+psb(i,j-1)+psb(i-1,j)+psb(i-1,j-1))
            uten(i,k,j) = uten(i,k,j) + (tpred1(i,k)-auxx(i,k,j))       &
                        & /dt*dumr
            vten(i,k,j) = vten(i,k,j) + (tpred2(i,k)-avxx(i,k,j))       &
                        & /dt*dumr
          end do
        end do
 
!       temperature
!
 
!       calculate coefficients at cross points for temperature
 
        do k = 1 , kx
          do i = 2 , ixm1
            if ( k.gt.1 ) betak(i,k) = rhohf(i,k-1,j)*kvh(i,k,j)        &
                                     & /dza(i,k-1,j)
            alphak(i,k) = gti/(psb(i,j)*1000.)/dsigma(k)
          end do
        end do
 
        do k = 2 , kx - 1
          do i = 2 , ixm1
            coef1(i,k) = dt*alphak(i,k)*betak(i,k+1)
            coef2(i,k) = 1. + dt*alphak(i,k)*(betak(i,k+1)+betak(i,k))
            coef3(i,k) = dt*alphak(i,k)*betak(i,k)
          end do
        end do
 
        do i = 2 , ixm1
          coef1(i,1) = dt*alphak(i,1)*betak(i,2)
          coef2(i,1) = 1. + dt*alphak(i,1)*betak(i,2)
          coef3(i,1) = 0.
          coef1(i,kx) = 0.
          coef2(i,kx) = 1. + dt*alphak(i,kx)*betak(i,kx)
          coef3(i,kx) = dt*alphak(i,kx)*betak(i,kx)
        end do
 
        do i = 2 , ixm1
          coefe(i,1) = coef1(i,1)/coef2(i,1)
          coeff1(i,1) = thx3d(i,1,j)/coef2(i,1)
        end do
 
        do k = 2 , kx - 1
          do i = 2 , ixm1
            coefe(i,k) = coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
            coeff1(i,k) = (thx3d(i,k,j)+coef3(i,k)*coeff1(i,k-1))       &
                        & /(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
          end do
        end do
 
        do i = 2 , ixm1
          coefe(i,kx) = 0.
          coeff1(i,kx) = (thx3d(i,kx,j)+dt*alphak(i,kx)*hfx(i,j)        &
                       & *rcpd+coef3(i,kx)*coeff1(i,kx-1))              &
                       & /(coef2(i,kx)-coef3(i,kx)*coefe(i,kx-1))
        end do
 
!
!       all coefficients have been computed, predict field and put it in
!       temporary work space tpred
!
 
        do i = 2 , ixm1
          tpred1(i,kx) = coeff1(i,kx)
        end do
 
        do k = kx - 1 , 1 , -1
          do i = 2 , ixm1
            tpred1(i,k) = coefe(i,k)*tpred1(i,k+1) + coeff1(i,k)
          end do
        end do
 
!
!       calculate tendency due to vertical diffusion using temporary
!       predicted field
!
        do k = 1 , kx
          do i = 2 , ixm1
            sf = tb(i,k,j)/thx3d(i,k,j)
            difft(i,k,j) = difft(i,k,j) + (tpred1(i,k)-thx3d(i,k,j))    &
                         & /dt*sf
          end do
        end do
!
!       water vapor
!
 
!       calculate coefficients at cross points for water vapor
 
        do k = 1 , kx
          do i = 2 , ixm1
            if ( k.gt.1 ) betak(i,k) = rhohf(i,k-1,j)*kvq(i,k,j)        &
                                     & /dza(i,k-1,j)
            alphak(i,k) = gti/(psb(i,j)*1000.)/dsigma(k)
          end do
        end do
 
        do k = 2 , kx - 1
          do i = 2 , ixm1
            coef1(i,k) = dt*alphak(i,k)*betak(i,k+1)
            coef2(i,k) = 1. + dt*alphak(i,k)*(betak(i,k+1)+betak(i,k))
            coef3(i,k) = dt*alphak(i,k)*betak(i,k)
          end do
        end do
 
        do i = 2 , ixm1
          coef1(i,1) = dt*alphak(i,1)*betak(i,2)
          coef2(i,1) = 1. + dt*alphak(i,1)*betak(i,2)
          coef3(i,1) = 0.
          coef1(i,kx) = 0.
          coef2(i,kx) = 1. + dt*alphak(i,kx)*betak(i,kx)
          coef3(i,kx) = dt*alphak(i,kx)*betak(i,kx)
        end do
 
        do i = 2 , ixm1
          coefe(i,1) = coef1(i,1)/coef2(i,1)
          coeff1(i,1) = qvb3d(i,1,j)/coef2(i,1)
        end do
 
        do k = 2 , kx - 1
          do i = 2 , ixm1
            coefe(i,k) = coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
            coeff1(i,k) = (qvb3d(i,k,j)+coef3(i,k)*coeff1(i,k-1))       &
                        & /(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
          end do
        end do
 
        do i = 2 , ixm1
          coefe(i,kx) = 0.
          coeff1(i,kx) = (qvb3d(i,kx,j)+dt*alphak(i,kx)*qfx(i,j)        &
                       & +coef3(i,kx)*coeff1(i,kx-1))                   &
                       & /(coef2(i,kx)-coef3(i,kx)*coefe(i,kx-1))
        end do
 
!
!       all coefficients have been computed, predict field and put it in
!       temporary work space tpred
 
        do i = 2 , ixm1
          tpred1(i,kx) = coeff1(i,kx)
        end do
 
        do k = kx - 1 , 1 , -1
          do i = 2 , ixm1
            tpred1(i,k) = coefe(i,k)*tpred1(i,k+1) + coeff1(i,k)
          end do
        end do
 
!
!       calculate tendency due to vertical diffusion using temporary
!       predicted field
!
        do k = 1 , kx
          do i = 2 , ixm1
            diffq(i,k,j) = diffq(i,k,j)                                 &
                         & + (tpred1(i,k)-qvb(i,k,j)/psb(i,j))          &
                         & /dt*psb(i,j)
          end do
        end do
 
!       calculate coefficients at cross points for cloud vater
 
        do k = 1 , kx
          do i = 2 , ixm1
            if ( k.gt.1 ) betak(i,k) = rhohf(i,k-1,j)*kvq(i,k,j)        &
                                     & /dza(i,k-1,j)
            alphak(i,k) = gti/(psb(i,j)*1000.)/dsigma(k)
          end do
        end do
 
        do k = 2 , kx - 1
          do i = 2 , ixm1
            coef1(i,k) = dt*alphak(i,k)*betak(i,k+1)
            coef2(i,k) = 1. + dt*alphak(i,k)*(betak(i,k+1)+betak(i,k))
            coef3(i,k) = dt*alphak(i,k)*betak(i,k)
          end do
        end do
 
        do i = 2 , ixm1
          coef1(i,1) = dt*alphak(i,1)*betak(i,2)
          coef2(i,1) = 1. + dt*alphak(i,1)*betak(i,2)
          coef3(i,1) = 0.
          coef1(i,kx) = 0.
          coef2(i,kx) = 1. + dt*alphak(i,kx)*betak(i,kx)
          coef3(i,kx) = dt*alphak(i,kx)*betak(i,kx)
        end do
 
        do i = 2 , ixm1
          coefe(i,1) = coef1(i,1)/coef2(i,1)
          coeff1(i,1) = qcx(i,1,j)/coef2(i,1)
        end do
 
        do k = 2 , kx - 1
          do i = 2 , ixm1
            coefe(i,k) = coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
            coeff1(i,k) = (qcx(i,k,j)+coef3(i,k)*coeff1(i,k-1))         &
                        & /(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
          end do
        end do
 
        do i = 2 , ixm1
          coefe(i,kx) = 0.
          coeff1(i,kx) = (qcx(i,kx,j)+coef3(i,kx)*coeff1(i,kx-1))       &
                       & /(coef2(i,kx)-coef3(i,kx)*coefe(i,kx-1))
        end do
 
!
!       all coefficients have been computed, predict field and put it in
!       temporary work space tpred
!
 
        do i = 2 , ixm1
          tpred1(i,kx) = coeff1(i,kx)
        end do
 
        do k = kx - 1 , 1 , -1
          do i = 2 , ixm1
            tpred1(i,k) = coefe(i,k)*tpred1(i,k+1) + coeff1(i,k)
          end do
        end do
 
!
!       calculate tendency due to vertical diffusion using temporary
!       predicted field
!
        do k = 1 , kx
          do i = 2 , ixm1
            qcten(i,k,j) = qcten(i,k,j)                                 &
                         & + (tpred1(i,k)-qcb(i,k,j)/psb(i,j))          &
                         & /dt*psb(i,j)
          end do
        end do
 
!
! **********************************************************************
!
!       now add countergradient term to temperature and water vapor
!       equation
!trapuv
        do i = 2 , ixm1
          ttnp(i,1) = 0.0D0
        end do
!trapuv_
        do k = 2 , kx
          do i = 2 , ixm1
            sf = tb(i,k,j)/(psb(i,j)*thx3d(i,k,j))
            ttnp(i,k) = sf*cpd*rhohf(i,k-1,j)*kvh(i,k,j)*cgh(i,k,j)
          end do
        end do
!
!-----compute the tendencies:
!
        do i = 2 , ixm1
          difft(i,kx,j) = difft(i,kx,j) - gti*ttnp(i,kx)                &
                        & /(1000.*cpd*dsigma(kx))
        end do
!
        do k = 1 , kxm1
          do i = 2 , ixm1
            difft(i,k,j) = difft(i,k,j) + gti*(ttnp(i,k+1)-ttnp(i,k))   &
                         & /(1000.*cpd*dsigma(k))
          end do
        end do
 
!
!chem2
        if ( ichem.eq.1 .and. ichdrdepo.eq.1 ) then
!
!         coef1, coef2, coef3 and coefe are the same as for water vapor
!         and cloud water so they do not need to be recalculated
 
!
!         recalculation of coef1,2,3  with tracer diffusivity kvc
 
          do k = 1 , kx
            do i = 2 , ixm1
              if ( k.gt.1 ) betak(i,k) = rhohf(i,k-1,j)*kvc(i,k,j)      &
                 & /dza(i,k-1,j)
              alphak(i,k) = gti/(psb(i,j)*1000.)/dsigma(k)
            end do
          end do
 
          do k = 2 , kx - 1
            do i = 2 , ixm1
              coef1(i,k) = dt*alphak(i,k)*betak(i,k+1)
              coef2(i,k) = 1. + dt*alphak(i,k)*(betak(i,k+1)+betak(i,k))
              coef3(i,k) = dt*alphak(i,k)*betak(i,k)
            end do
          end do
 
          do i = 2 , ixm1
            coef1(i,1) = dt*alphak(i,1)*betak(i,2)
            coef2(i,1) = 1. + dt*alphak(i,1)*betak(i,2)
            coef3(i,1) = 0.
            coef1(i,kx) = 0.
            coef2(i,kx) = 1. + dt*alphak(i,kx)*betak(i,kx)
            coef3(i,kx) = dt*alphak(i,kx)*betak(i,kx)
          end do
!
!         set the Vd for in case of prescribed deposition velocities
!
 
          do itr = 1 , ntr
            do i = 2 , ixm1
              if ( veg2d(i,j).le.0.00001 ) then
                vdep(i,itr) = chtrdpv(itr,2)
              else
                vdep(i,itr) = chtrdpv(itr,1)
              end if
!             provisoire test de la routine chdrydep pour les dust
 
              if ( chtrname(itr).eq.'DUST' ) vdep(i,itr) = 0.
            end do
          end do
!
          do itr = 1 , ntr
!
            do k = 1 , kx
              do i = 2 , ixm1
                chix(i,k) = chib(i,k,j,itr)/psb(i,j)
              end do
            end do
!
            do i = 2 , ixm1
              coefe(i,1) = coef1(i,1)/coef2(i,1)
              coeff1(i,1) = chix(i,1)/coef2(i,1)
            end do
!
            do k = 2 , kx - 1
              do i = 2 , ixm1
                coefe(i,k) = coef1(i,k)                                 &
                           & /(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
                coeff1(i,k) = (chix(i,k)+coef3(i,k)*coeff1(i,k-1))      &
                            & /(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
              end do
            end do
 
            do i = 2 , ixm1
              coefe(i,kx) = 0.
 
!             add dry deposition option1
              coeff1(i,kx) = (chix(i,kx)-dt*alphak(i,kx)*chix(i,kx)     &
                           & *vdep(i,itr)*rhox2d(i,j)+coef3(i,kx)       &
                           & *coeff1(i,kx-1))                           &
                           & /(coef2(i,kx)-coef3(i,kx)*coefe(i,kx-1))
            end do
!
!           all coefficients have been computed, predict field and put
!           it in temporary work space tpred1
!
            do i = 2 , ixm1
              tpred1(i,kx) = coeff1(i,kx)
            end do
!
            do k = kx - 1 , 1 , -1
              do i = 2 , ixm1
                tpred1(i,k) = coefe(i,k)*tpred1(i,k+1) + coeff1(i,k)
              end do
            end do
!
!
!           calculate tendency due to vertical diffusion using temporary
!           predicted field
!           Dry deposition option 1 is included
 
!
            do k = 1 , kx
              do i = 2 , ixm1
!qian           chiten(i,k,j,itr)=chiten(i,k,j,itr)
!CGAFFE         TEST diffusion/10
                chiten(i,k,j,itr) = chiten(i,k,j,itr)                   &
                                  & + (tpred1(i,k)-chix(i,k))           &
                                  & /dt*psb(i,j)
!               chiten(i,k,j,itr)=chiten(i,k,j,itr)+0.1 *(tpred1(i,k)-
!               1  chix(i,k))/dt *psb(i,j)
 
              end do
            end do
            do i = 2 , ixm1
 
              if ( chtrname(itr).ne.'DUST' ) remdrd(i,j,itr)            &
                 & = remdrd(i,j,itr) + chix(i,kx)*vdep(i,itr)*psb(i,j)  &
                 & *dt/2.*rhox2d(i,j)*gti/(psb(i,j)*1000.*dsigma(kx))
 
            end do
          end do
        end if
!chem2_
 
      end do
!
      end subroutine holtbl
