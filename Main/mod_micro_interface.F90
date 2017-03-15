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

module mod_micro_interface

  use mod_realkinds
  use mod_service
  use mod_constants
  use mod_dynparam
  use mod_memutil
  use mod_mppparam
  use mod_regcm_types
  use mod_runparams
  use mod_micro_nogtom
  use mod_micro_subex
  use mod_micro_wsm5
  use mod_cloud_subex
  use mod_cloud_xuran
  use mod_cloud_thomp

  implicit none

  private

  public :: allocate_micro , init_micro , microscheme , cldfrac , condtq

  type(nogtom_stats) , public :: ngs

  type(mod_2_micro) :: mo2mc
  type(micro_2_mod) :: mc2mo

  real(rkx) , pointer , dimension(:,:,:) :: totc

  real(rkx) , parameter :: qvmin = 1.0e-8_rkx
  real(rkx) , parameter :: alphaice = d_one
  real(rkx) , parameter :: eps = 1.0e-7_rkx

  integer(ik4) , parameter :: nchi = 256
  real(rkx) , dimension(0:nchi-1) :: chis

  public :: qck1 , cgul , rh0 , cevap , xcevap , caccr

  contains

#include <pfesat.inc>
#include <pfwsat.inc>
#include <clwfromt.inc>

  subroutine allocate_micro
    implicit none
    integer(ik4) :: i
    real(rkx) :: cf
    call allocate_subex
    if ( ipptls == 2 ) then
      call allocate_mod_nogtom
#ifdef DEBUG
      if ( stats ) then
        call getmem3d(ngs%statssupw,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statssupw')
        call getmem3d(ngs%statssupc,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statssupc')
        call getmem3d(ngs%statsdetrw,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsdetrw')
        call getmem3d(ngs%statsdetrc,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsdetrc')
        call getmem3d(ngs%statserosw,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statserosw')
        call getmem3d(ngs%statserosc,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statserosc')
        call getmem3d(ngs%statsevapw,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsevapw')
        call getmem3d(ngs%statsevapc,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsevapc')
        call getmem3d(ngs%statscond1w,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statscond1w')
        call getmem3d(ngs%statscond1c,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statscond1c')
        call getmem3d(ngs%statscond2w,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statscond2w')
        call getmem3d(ngs%statscond2c,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statscond2c')
        call getmem3d(ngs%statsdepos,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsdepos')
        call getmem3d(ngs%statsmelt,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsmelt')
        call getmem3d(ngs%statsfrz,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsfrz')
        call getmem3d(ngs%statsrainev,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsrainev')
        call getmem3d(ngs%statssnowev,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statssnowev')
        call getmem3d(ngs%statsautocvw,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsautocvw')
        call getmem3d(ngs%statsautocvc,jci1,jci2, &
                                    ici1,ici2,1,kz,'micro:statsautocvc')
      end if
#endif
    else if ( ipptls == 3 ) then
      call allocate_mod_wsm5
    end if
    call getmem3d(totc,jci1,jci2,ici1,ici2,1,kz,'subex:totc')
    do i = 1 , nchi
      cf = real(i-1)/real(nchi-1)
      chis(i-1) = 0.97_rkx*exp(-((cf-0.098_rkx)**2)/0.0365_rkx)+0.255_rkx
    end do
  end subroutine allocate_micro

  subroutine init_micro
    use mod_atm_interface
    use mod_che_interface
    implicit none

    call assignpnt(mddom%xlat,mo2mc%xlat)
    call assignpnt(sfs%psb,mo2mc%psb)
    call assignpnt(atms%pb3d,mo2mc%phs)
    call assignpnt(atms%pf3d,mo2mc%pfs)
    call assignpnt(atms%tb3d,mo2mc%t)
    call assignpnt(atms%dzq,mo2mc%delz)
    call assignpnt(atms%wpx3d,mo2mc%pverv)
    call assignpnt(atms%wb3d,mo2mc%verv)
    call assignpnt(atms%qxb3d,mo2mc%qxx)
    call assignpnt(atms%rhob3d,mo2mc%rho)
    call assignpnt(atms%rhb3d,mo2mc%rh)
    call assignpnt(atms%qsb3d,mo2mc%qs)
    call assignpnt(heatrt,mo2mc%heatrt)
    call assignpnt(q_detr,mo2mc%qdetr)
    call assignpnt(fcc,mo2mc%fcc)

    call assignpnt(atms%qxb3d,mo2mc%qvn,iqv)
    call assignpnt(atms%qxb3d,mo2mc%qcn,iqc)
    if ( ipptls > 1 ) then
      call assignpnt(atms%qxb3d,mo2mc%qin,iqi)
      call assignpnt(atms%qxb3d,mo2mc%qsn,iqs)
      call assignpnt(atms%qxb3d,mo2mc%qrn,iqr)
    end if

    if ( ichem == 1 ) then
      if ( iaerosol == 1 .and. iindirect == 2 ) then
        call assignpnt(ccn,mo2mc%ccn)
      end if
      if ( idiag == 1 ) then
        call assignpnt(qdiag%qcr,mc2mo%dia_qcr)
        call assignpnt(qdiag%qcl,mc2mo%dia_qcl)
        call assignpnt(qdiag%acr,mc2mo%dia_acr)
      end if
    end if

    call assignpnt(aten%qx,mc2mo%qxten,pc_physic)
    call assignpnt(aten%t,mc2mo%tten,pc_physic)
    call assignpnt(sfs%rainnc,mc2mo%rainnc)
    call assignpnt(sfs%snownc,mc2mo%snownc)
    call assignpnt(pptnc,mc2mo%lsmrnc)
    call assignpnt(rain_ls,mc2mo%rainls)
    call assignpnt(remrat,mc2mo%remrat)
    call assignpnt(rembc,mc2mo%rembc)

    select case ( ipptls )
      case (1)
        call init_subex(mddom%xlat)
      case (2)
        call init_nogtom(mddom%ldmsk)
      case(3)
        call init_wsm5
      case default
        return
    end select
  end subroutine init_micro

  subroutine microscheme
    use mod_atm_interface
    implicit none
    select case ( ipptls )
      case (1)
        call subex(mo2mc,mc2mo)
      case (2)
        call nogtom(mo2mc,ngs,mc2mo)
      case (3)
        call wsm5(mo2mc,mc2mo)
      case default
        return
    end select
  end subroutine microscheme
  !
  ! This subroutine computes the fractional cloud coverage and
  ! liquid water content (in cloud value).  Both are use in
  ! radiation.
  !
  subroutine cldfrac
    use mod_atm_interface , only : mddom , atms , cldlwc , cldfra
    implicit none
    real(rkx) :: exlwc
    integer(ik4) :: i , j , k , ichi

    select case ( icldfrac )
      case (1)
        call xuran_cldfrac(mo2mc%phs,mo2mc%qcn,mo2mc%qvn, &
                           mo2mc%qs,mo2mc%rh,mo2mc%fcc)
      case (2)
        call thomp_cldfrac(mo2mc%phs,mo2mc%t,mo2mc%rho,mo2mc%qvn,     &
                           mo2mc%qcn,mo2mc%qsn,mo2mc%qin,mddom%ldmsk, &
                           ds,mo2mc%fcc)
      case default
        call subex_cldfrac(mo2mc%t,mo2mc%phs,mo2mc%qvn, &
                           mo2mc%rh,tc0,rh0,mo2mc%fcc)
    end select

    if ( ipptls > 1 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            totc(j,i,k) = (mo2mc%qcn(j,i,k) + alphaice*mo2mc%qin(j,i,k))
          end do
        end do
      end do
    else
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            totc(j,i,k) = mo2mc%qcn(j,i,k)
          end do
        end do
      end do
    end if

    !------------------------------------------
    ! 1a. Determine Marine stratocumulus clouds
    !------------------------------------------

    if ( icldmstrat == 1 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( mddom%ldmsk(j,i) == 0 ) then
              if ( mo2mc%phs(j,i,k) >= 70000.0_rkx ) then
                ! Klein, S. A., and D. L. Hartmann,
                ! The seasonal cycle of low stratiform clouds,
                ! J. Climate, 6, 1587-1606, 1993
                mo2mc%fcc(j,i,k) = max(mo2mc%fcc(j,i,k), &
                      (atms%th700(j,i)-atms%th3d(j,i,k)) * &
                           0.057_rkx - 0.5573_rkx)
              end if
            end if
          end do
        end do
      end do
    end if

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          mo2mc%fcc(j,i,k) = max(min(mo2mc%fcc(j,i,k),hicld),lowcld)
        end do
      end do
    end do

    !-----------------------------------------------------------------
    ! 2.  Combine large-scale and convective fraction and liquid water
    !     to be passed into radiation.
    !-----------------------------------------------------------------

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          exlwc = dlowval
          ! Cloud Water Volume
          ! kg gq / kg dry air * kg dry air / m3 * 1000 = g qc / m3
          if ( iconvlwp == 1 ) then
            ! Apply the parameterisation based on temperature to the
            ! the large scale clouds.
            ! Scaling for CF
            ! Implements CF scaling as in Liang GRL 32, 2005
            ! doi: 10.1029/2004GL022301
            if ( mo2mc%fcc(j,i,k) > lowcld+eps ) then
              ichi = int(mo2mc%fcc(j,i,k)*real(nchi-1,rkx))
              exlwc = clwfromt(mo2mc%t(j,i,k)) * chis(ichi)
            end if
          else
            ! NOTE : IN CLOUD HERE IS NEEDED !!!
            if ( mo2mc%fcc(j,i,k) > lowcld+eps ) then
              ! In g / m^3
              exlwc = ((totc(j,i,k)*d_1000)/mo2mc%fcc(j,i,k))*mo2mc%rho(j,i,k)
            end if
          end if
          if ( cldfra(j,i,k) > lowcld+eps .or. &
               mo2mc%fcc(j,i,k) > lowcld+eps ) then
            ! get maximum cloud fraction between cumulus and large scale
            cldlwc(j,i,k) = (exlwc * mo2mc%fcc(j,i,k) + &
                            cldlwc(j,i,k) * cldfra(j,i,k)) / &
                            (cldfra(j,i,k) + mo2mc%fcc(j,i,k))
            cldfra(j,i,k) = max(cldfra(j,i,k),mo2mc%fcc(j,i,k))
          else
            cldfra(j,i,k) = lowcld
            cldlwc(j,i,k) = dlowval
          end if
          cldfra(j,i,k) = min(max(cldfra(j,i,k),lowcld),cftotmax)
        end do
      end do
    end do
  end subroutine cldfrac
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !                                                                 c
  ! This subroutine computes the condensational or evaporational    c
  ! heating term from the explicit moisture scheme.                 c
  !                                                                 c
  ! ---the condensational or evaporational term are one step        c
  !    adjustment based on asai (1965, j. meteo. soc. japan).       c
  !                                                                 c
  ! ---modified to include the effects of partial cloud cover       c
  !    (see Pal et al 2000).  When partial clouds exist, the qxten  c
  !    in/out of the clear and cloudy portions of the grid cell is  c
  !    assumed to be at the same rate (i.e., if there is 80% cloud  c
  !    cover, .2 of qxten goes to raising qv in the clear region    c
  !    and .8 goes to condensation or evaporation of qc in the      c
  !    cloudy portion).                                             c
  !                                                                 c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine condtq
    use mod_atm_interface , only : atm0 , atm2 , sfs , aten
    implicit none
    !
    ! rhc    - Relative humidity at ktau+1
    ! rh0adj - Adjusted relative humidity threshold at ktau+1
    ! fccc   - Cloud fraction at ktau+1
    !
    real(rkx) :: qccs , qvcs , tmp1 , tmp2 , tmp3
    real(rkx) :: dqv , exces , fccc , pres , qvc_cld , qvs , &
               r1 , rh0adj , rhc
    integer(ik4) :: i , j , k

    !---------------------------------------------------------------------
    !     1.  Compute t, qv, and qc at tau+1 without condensational term
    !---------------------------------------------------------------------
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          tmp3 = (atm2%t(j,i,k)+dt*aten%t(j,i,k,pc_total))/sfs%psc(j,i)
#ifdef DEBUG
          if ( tmp3 < d_zero ) then
            write(stderr,*) 'Time ktau = ', ktau
            write(stderr,*) 'Consistency TEMPERATURE ERROR in condtq (T < 0K)'
            write(stderr,*) 'At global J : ',j
            write(stderr,*) 'At global I : ',i
            write(stderr,*) 'At global K : ',k
          end if
#endif
          qvcs = max((atm2%qx(j,i,k,iqv) + &
                   dt*aten%qx(j,i,k,iqv,pc_total)),qvmin)/sfs%psc(j,i)
          qccs = max((atm2%qx(j,i,k,iqc) + &
                   dt*aten%qx(j,i,k,iqc,pc_total)),d_zero)/sfs%psc(j,i)
          !
          ! 2.  Compute the cloud condensation/evaporation term.
          !
          ! 2a. Calculate the saturation mixing ratio and relative humidity
          if ( idynamic == 1 ) then
            pres = (hsigma(k)*sfs%psc(j,i)+ptop)*d_1000
          else
            pres = atm0%pr(j,i,k) + &
               (atm2%pp(j,i,k)+dt*aten%pp(j,i,k,pc_total))/sfs%psc(j,i)
          end if
          qvs = pfwsat(tmp3,pres)
          rhc = min(max(qvcs/qvs,rhmin),rhmax)

          r1 = d_one/(d_one+wlhv*wlhv*qvs/(rwat*cpd*tmp3*tmp3))

          ! 2b. Compute the relative humidity threshold at ktau+1
          if ( tmp3 > tc0 ) then
            rh0adj = rh0(j,i)
          else ! high cloud (less subgrid variability)
            rh0adj = rhmax - (rhmax-rh0(j,i))/(d_one+0.15_rkx*(tc0-tmp3))
          end if
          rh0adj = max(rhmin,min(rh0adj,rhmax))
          if ( rhc < rh0adj ) then      ! Low cloud cover
            dqv = conf * (qvcs - qvs)
          else if ( rhc >= rhmax ) then ! Full cloud cover
            dqv = conf * (qvcs - qvs)
          else
            fccc = d_one-sqrt(d_one-(rhc-rh0adj)/(rhmax-rh0adj))
            if ( pres >= 75000.0_rkx .and. qvcs <= 0.003_rkx ) then
              fccc = fccc * max(0.15_rkx, min(d_one,qvcs/0.003_rkx))
            end if
            fccc = min(max(fccc,d_zero),d_one)
            qvc_cld = max((mo2mc%qs(j,i,k) + &
                     dt * mc2mo%qxten(j,i,k,iqv)/sfs%psc(j,i)),d_zero)
            dqv = conf * fccc * (qvc_cld - qvs) ! qv diff between predicted qv_c
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
            aten%qx(j,i,k,iqv,pc_physic) = &
                aten%qx(j,i,k,iqv,pc_physic) - sfs%psc(j,i)*tmp2
            aten%qx(j,i,k,iqc,pc_physic) = &
                aten%qx(j,i,k,iqc,pc_physic) + sfs%psc(j,i)*tmp2
            aten%t(j,i,k,pc_physic) = &
                aten%t(j,i,k,pc_physic) + sfs%psc(j,i)*tmp2*wlhvocp
          end if
        end do
      end do
    end do
  end subroutine condtq

end module mod_micro_interface
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
