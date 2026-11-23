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
! Historical, not compiled: the condensation of the SUBEX microphysics
! (ipptls = 1, mod_micro_subex), once in mod_micro_interface, called by the
! tendency of the hydrostatic and non hydrostatic cores (mod_tendency,
! after the sum of the tendencies) and of MOLOCH (mod_moloch). The
! threshold rh0adj was computed in cldfrac from the large scale cloud
! fraction fcc:
!
!   rh0adj = 1 - (1-rh)/max(1-fcc,1e-5)**2
!
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine condtq
    use mod_atm_interface, only : mo_atm, atm0, atm2, sfs, aten
    implicit none
    integer(ik4) :: i, j, k
    real(rkx) :: qccs, qvcs, tmp1, tmp2, tmp3
    real(rkx) :: dqv, exces, pres, qvc_cld, qvs, fccc, &
                 r1, rhc, rlv, cpm

    !---------------------------------------------------------------------
    !     1.  Compute t, qv, and qc at tau+1 without condensational term
    !---------------------------------------------------------------------
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      if ( idynamic == 3 ) then
        tmp3 = mo_atm%t(j,i,k) + dt*mo_atm%tten(j,i,k)
        qvcs = max(mo_atm%qx(j,i,k,iqv) + dt*mo_atm%qxten(j,i,k,iqv),minqq)
        qccs = max(mo_atm%qx(j,i,k,iqc) + dt*mo_atm%qxten(j,i,k,iqc),d_zero)
        pres = mo_atm%p(j,i,k)
        qvc_cld = max((mo2mc%qs(j,i,k) + &
                  dt * mc2mo%qxten(j,i,k,iqv)),minqq)
      else
        tmp3 = (atm2%t(j,i,k)+dt*aten%t(j,i,k,pc_total))/sfs%psc(j,i)
        qvcs = atm2%qx(j,i,k,iqv) + dt*aten%qx(j,i,k,iqv,pc_total)
        qccs = atm2%qx(j,i,k,iqc) + dt*aten%qx(j,i,k,iqc,pc_total)
        qvc_cld = max((mo2mc%qs(j,i,k) + &
                   dt * mc2mo%qxten(j,i,k,iqv)/sfs%psc(j,i)),minqq)
        if ( idynamic == 1 ) then
          pres = (hsigma(k)*sfs%psc(j,i)+ptop)*d_1000
        else
          pres = atm0%pr(j,i,k) + &
             (atm2%pp(j,i,k)+dt*aten%pp(j,i,k,pc_total))/sfs%psc(j,i)
        end if
        if ( qvcs < minqq * sfs%psc(j,i) ) then
          qvcs = minqq * sfs%psc(j,i)
        end if
        if ( qccs < dlowval * sfs%psc(j,i) ) then
          qccs = d_zero
        end if
        qvcs = qvcs /sfs%psc(j,i)
        qccs = qccs /sfs%psc(j,i)
      end if
#ifndef OPENACC
#ifdef DEBUG
      if ( tmp3 < d_zero ) then
        write(stderr,*) 'Time step = ', rcmtimer%lcount
        write(stderr,*) 'Consistency TEMPERATURE ERROR in condtq (T < 0K)'
        write(stderr,*) 'At global J : ',j
        write(stderr,*) 'At global I : ',i
        write(stderr,*) 'At global K : ',k
      end if
#endif
#endif
      !
      ! 2.  Compute the cloud condensation/evaporation term.
      !
      ! 2a. Calculate the saturation mixing ratio and relative humidity
      qvs = pfwsat(tmp3,pres)
      rlv = wlh(tmp3)
      cpm = cpd*(d_one-qvcs) + cpv*qvcs
      r1 = d_one/(d_one+rlv*rlv*qvs/(rwat*cpm*tmp3*tmp3))
      rhc = min(max(qvcs/qvs,rhmin),rhmax)
      ! 2b. Compute the relative humidity threshold at ktau+1
      if ( rhc < rh0adj(j,i,k) ) then  ! Low cloud cover
        dqv = conf * (qvcs - qvs)
      else if ( rhc > 0.99999_rkx ) then
        dqv = conf * (qvcs - qvs)      ! High cloud cover
      else
        fccc = d_one-sqrt((rhmax-rhc)/(rhmax-rh0adj(j,i,k)))
        fccc = min(max(fccc,d_zero),d_one)
        ! qv diff between predicted qv_c
        dqv = conf * fccc * (qvc_cld - qvs)
      end if

      ! 2c. Compute the water vapor in excess of saturation
      tmp1 = r1*dqv               ! grid cell average

      ! 2d. Compute the new cloud water + old cloud water
      exces = qccs + tmp1
      if ( exces >= d_zero ) then ! Some cloud is left
        tmp2 = tmp1/dt
      else                        ! The cloud evaporates
        tmp2 = -qccs/dt
      end if
      !
      ! 3. Compute the tendencies.
      !
      if ( abs(tmp2) > dlowval ) then
        if ( idynamic == 3 ) then
          mo_atm%qxten(j,i,k,iqv) = mo_atm%qxten(j,i,k,iqv) - tmp2
          mo_atm%qxten(j,i,k,iqc) = mo_atm%qxten(j,i,k,iqc) + tmp2
          mo_atm%tten(j,i,k) = mo_atm%tten(j,i,k) + tmp2*rlv/cpm
        else
          aten%qx(j,i,k,iqv,pc_physic) = &
              aten%qx(j,i,k,iqv,pc_physic) - sfs%psc(j,i)*tmp2
          aten%qx(j,i,k,iqc,pc_physic) = &
              aten%qx(j,i,k,iqc,pc_physic) + sfs%psc(j,i)*tmp2
          aten%t(j,i,k,pc_physic) = &
              aten%t(j,i,k,pc_physic) + sfs%psc(j,i)*tmp2*rlv/cpm
        end if
      end if
    end do

    contains

#include <saturation.inc>
#include <wlh.inc>

  end subroutine condtq
