!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! Historical, not compiled: the setup of the Kuo scheme (mod_cu_kuo) once
! in param of mod_params (icup = 1), with vqrang = 5.0e-4 of
! mod_runparams: the heating profile twght and the weights of the
! convective moisture flux vqflx by cloud base and top, and the level
! k700 of 700 hPa. In the hydrostatic core the vertical advection of the
! water vapour (vadv of qxdyn in mod_tendency) was skipped with Kuo, which
! used the vertical moisture convergence (aten%qx, pc_dynamic) as its
! closure.
!
    if ( any(icup == 1) ) then
      !
      ! specify heating profile (twght) and weighted function
      ! for moisture fluxes due to convection (vqflx)
      ! assume base of cloud varies as  < kbase = 5,kz >
      ! top  of cloud varies as  < ktop  = 1,kbase-3 >
      ! exceptions to this are treated explicitly in subroutine
      ! "cupara".
      !
      do kbase = 5, kz
        do ktop = 1, kbase - 3
          do k = 1, kz
            twght(k,kbase,ktop) = d_zero
            vqflx(k,kbase,ktop) = d_zero
          end do
          !
          ! get twght from 1/2 level sigma values
          !
          bb = log(hsigma(ktop)) + log(hsigma(kbase))
          cc = log(hsigma(ktop))*log(hsigma(kbase))
          ssum = d_zero
          do k = ktop, kbase
            xx = log(hsigma(k))
            twght(k,kbase,ktop) = (xx*xx) - (bb*xx) + cc
            ssum = ssum + twght(k,kbase,ktop)*dsigma(k)
          end do
          do k = ktop, kbase
            twght(k,kbase,ktop) = twght(k,kbase,ktop)/ssum
          end do
          !
          ! get vqflx from  d(w*q) / dsigma on full levels
          ! do computations in p to avoid sigma=0. discontinuity
          !
          xtop = log((d_100-ptop)*sigma(ktop)+ptop)
          xbot = log((d_100-ptop)*sigma(kbase+1)+ptop)
          bb = xtop + xbot
          cc = xtop*xbot
          vqmax = d_zero
          ssum = d_zero
          xx = xtop
          yy = xbot
          wk = (xx*xx) - (bb*xx) + cc
          qk = -((yy*yy)-(bb*yy)+cc)
          do k = ktop, kbase
            xx = log((d_100-ptop)*sigma(k+1)+ptop)
            yy = log((d_100-ptop) * &
                 (sigma(ktop)+sigma(kbase+1)-sigma(k+1))+ptop)
            wkp1 = (xx*xx) - (bb*xx) + cc
            qkp1 = -((yy*yy)-(bb*yy)+cc)
            vqflx(k,kbase,ktop) = -((wkp1*qkp1)-(wk*qk))/dsigma(k)
            ssum = ssum + vqflx(k,kbase,ktop)
            if ( abs(vqflx(k,kbase,ktop)) > vqmax ) &
                 vqmax = abs(vqflx(k,kbase,ktop))
            wk = wkp1
            qk = qkp1
          end do
          do k = ktop, kbase
            vqflx(k,kbase,ktop) = vqflx(k,kbase,ktop)*vqrang/vqmax
          end do
        end do
      end do
      if ( myid == italk ) then
        write(stdout, *) 'Anthes-Kuo Convection Scheme used.'
      end if
    end if
    if ( any(icup == 1) ) then
      sig700 = (70.0_rkx-ptop)/(d_100-ptop)
      do k = 1, kz
        k700 = k
        if ( sig700 <= sigma(k+1) .and. sig700 > sigma(k) ) exit
      end do
    end if
