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

module mod_micro_interface

  use mod_realkinds
  use mod_service
  use mod_stdio
  use mod_constants
  use mod_dynparam
  use mod_memutil
  use mod_mppparam
  use mod_regcm_types
  use mod_runparams
  use mod_micro_nogtom
  use mod_micro_subex
  use mod_micro_wsm5
  use mod_micro_wsm7
  use mod_micro_wdm7
  use mod_cloud_subex
  use mod_cloud_xuran
  use mod_cloud_thomp
  use mod_cloud_guli2007
  use mod_cloud_texeira
  use mod_cloud_tompkins
  use mod_cloud_echam5

  implicit none

  private

  public :: allocate_micro, init_micro, microscheme, cldfrac, condtq

  type(nogtom_stats), public :: ngs

  type(mod_2_micro) :: mo2mc
  type(micro_2_mod) :: mc2mo

  real(rkx), pointer, contiguous, dimension(:,:) :: rh0 => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: qtcrit => null( )
  ! rh0adj - Adjusted relative humidity threshold
  real(rkx), pointer, contiguous, dimension(:,:,:) :: totc => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: rh0adj => null( )

  integer(ik4), parameter :: nchi = 256
  real(rkx), dimension(0:nchi-1) :: chis

  logical, parameter :: do_cfscaling = .true.
  real(rkx), parameter :: qccrit_lnd = 1.0e-9_rkx
  real(rkx), parameter :: qccrit_oce = 5.0e-9_rkx

  public :: qck1, cgul, rh0, cevap, xcevap, caccr

  contains

  subroutine allocate_micro
    implicit none
    integer(ik4) :: i
    real(rkx) :: cf
    if ( ipptls == 1 ) then
      call allocate_subex
      call getmem(rh0adj,jci1,jci2,ici1,ici2,1,kz,'micro:rh0adj')
    else if ( ipptls == 2 ) then
      call allocate_mod_nogtom
#ifdef DEBUG
      if ( stats ) then
        call getmem(ngs%statssupw,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statssupw')
        call getmem(ngs%statssupc,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statssupc')
        call getmem(ngs%statsdetrw,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statsdetrw')
        call getmem(ngs%statsdetrc,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statsdetrc')
        call getmem(ngs%statserosw,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statserosw')
        call getmem(ngs%statserosc,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statserosc')
        call getmem(ngs%statsevapw,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statsevapw')
        call getmem(ngs%statsevapc,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statsevapc')
        call getmem(ngs%statscond1w,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statscond1w')
        call getmem(ngs%statscond1c,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statscond1c')
        call getmem(ngs%statsdepos,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statsdepos')
        call getmem(ngs%statsmelt,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statsmelt')
        call getmem(ngs%statsfrz,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statsfrz')
        call getmem(ngs%statsrainev,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statsrainev')
        call getmem(ngs%statssnowev,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statssnowev')
        call getmem(ngs%statsautocvw,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statsautocvw')
        call getmem(ngs%statsautocvc,1,kz,jci1,jci2, &
                                    ici1,ici2,'micro:statsautocvc')
      end if
#endif
    else if ( ipptls == 3 ) then
      call allocate_mod_wsm5
    else if ( ipptls == 4 ) then
      call allocate_mod_wsm7
    else if ( ipptls == 5 ) then
      call allocate_mod_wdm7
    end if
    call getmem(rh0,jci1,jci2,ici1,ici2,'micro:rh0')
    call getmem(qtcrit,jci1,jci2,ici1,ici2,'micro:qtcrit')
    call getmem(totc,jci1,jci2,ici1,ici2,1,kz,'micro:totc')
    do i = 1, nchi
      cf = real(i-1,rkx)/real(nchi-1,rkx)
      chis(i-1) = 0.97_rkx*exp(-((cf-0.098_rkx)**2)/0.0365_rkx)+0.255_rkx
    end do
  end subroutine allocate_micro

  subroutine init_micro
    use mod_atm_interface
    use mod_che_interface
    implicit none

    call assignpnt(mddom%ldmsk,mo2mc%ldmsk)
    call assignpnt(mddom%iveg,mo2mc%iveg)
    call assignpnt(mddom%xlat,mo2mc%xlat)
    call assignpnt(sfs%psb,mo2mc%psb)
    call assignpnt(atms%pb3d,mo2mc%phs)
    call assignpnt(atms%pf3d,mo2mc%pfs)
    call assignpnt(atms%tb3d,mo2mc%t)
    call assignpnt(atms%za,mo2mc%z)
    call assignpnt(atms%dzq,mo2mc%delz)
    call assignpnt(atms%wpx3d,mo2mc%pverv)
    call assignpnt(atms%wb3d,mo2mc%verv)
    call assignpnt(atms%qxb3d,mo2mc%qxx)
    call assignpnt(atms%rhob3d,mo2mc%rho)
    call assignpnt(atms%rhb3d,mo2mc%rh)
    call assignpnt(atms%qsb3d,mo2mc%qs)
    call assignpnt(atms%ps2d,mo2mc%ps2)
    call assignpnt(heatrt,mo2mc%heatrt)
    call assignpnt(q_detr,mo2mc%qdetr)
    call assignpnt(fcc,mo2mc%cldf)

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
      if ( idiag > 0 ) then
        call assignpnt(qdiag%qcr,mc2mo%dia_qcr)
        call assignpnt(qdiag%qcl,mc2mo%dia_qcl)
        call assignpnt(qdiag%acr,mc2mo%dia_acr)
      end if
    end if

    call assignpnt(fcc,mc2mo%fcc)
    if ( idynamic == 3 ) then
      call assignpnt(mo_atm%tten,mc2mo%tten)
      call assignpnt(mo_atm%qxten,mc2mo%qxten)
    else
      call assignpnt(aten%t,mc2mo%tten,pc_physic)
      call assignpnt(aten%qx,mc2mo%qxten,pc_physic)
    end if
    call assignpnt(sfs%rainnc,mc2mo%rainnc)
    call assignpnt(sfs%snownc,mc2mo%snownc)
    call assignpnt(sfs%grplnc,mc2mo%grplnc)
    call assignpnt(sfs%hailnc,mc2mo%hailnc)
    call assignpnt(pptnc,mc2mo%lsmrnc)
    call assignpnt(rain_ls,mc2mo%rainls)
    call assignpnt(remrat,mc2mo%remrat)
    call assignpnt(rembc,mc2mo%rembc)
    call assignpnt(ncrrate,mc2mo%trrate)

    select case ( ipptls )
      case (1)
        call init_subex(mddom%xlat)
      case (2)
        call init_nogtom(mddom%ldmsk)
      case(3)
        call init_wsm5
      case(4)
        call init_wsm7
      case(5)
        call init_wdm7
      case default
        return
    end select
  end subroutine init_micro

  subroutine microscheme
    implicit none
    select case ( ipptls )
      case (1)
        call subex(mo2mc,mc2mo)
      case (2)
#ifdef DEBUG
        call nogtom(mo2mc,mc2mo,ngs)
#else
        call nogtom(mo2mc,mc2mo)
#endif
      case (3)
        call wsm5(mo2mc,mc2mo)
      case (4)
        call wsm7(mo2mc,mc2mo)
      case (5)
        call wdm7(mo2mc,mc2mo)
      case default
        return
    end select
  end subroutine microscheme
  !
  ! This subroutine computes the fractional cloud coverage and
  ! liquid water content (in cloud value).  Both are use in
  ! radiation.
  !
  subroutine cldfrac(cldlwc,cldfra)
    use mod_atm_interface, only : atms
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: cldlwc, cldfra
    integer(ik4) :: i, j, k
    real(rkx) :: ls_exlwc, conv_exlwc, exlwc
    integer(ik4) :: ichi

    if ( ipptls > 1 ) then
      if ( icldfrac == 3 ) then
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          totc(j,i,k) = mo2mc%qcn(j,i,k) + mo2mc%qin(j,i,k) + &
                        mo2mc%qrn(j,i,k) + mo2mc%qsn(j,i,k)
        end do
      else
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          totc(j,i,k) = mo2mc%qcn(j,i,k) + mo2mc%qin(j,i,k)
        end do
      end if
    else
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        totc(j,i,k) = mo2mc%qcn(j,i,k)
      end do
    end if
    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      totc(j,i,k) = totc(j,i,k)/(d_one+totc(j,i,k))
    end do

    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      qtcrit(j,i) = mo2mc%ldmsk(j,i)*qccrit_lnd - &
                   (mo2mc%ldmsk(j,i)-1)*qccrit_oce
    end do

    select case ( icldfrac )
      case (0)
        call subex_cldfrac(mo2mc%t,mo2mc%phs,mo2mc%qvn, &
                           totc,mo2mc%rh,tc0,rh0,qtcrit,mc2mo%fcc)
      case (1)
        call xuran_cldfrac(mo2mc%phs,totc,mo2mc%qs,mo2mc%rh, &
                           qtcrit,mc2mo%fcc)
      case (2)
        call thomp_cldfrac(mo2mc%phs,mo2mc%t,mo2mc%rho,mo2mc%qvn, &
                           totc,mo2mc%qsn,mo2mc%qin,mo2mc%ldmsk,  &
                           ds,mc2mo%fcc)
      case (3)
        call gulisa_cldfrac(totc,mo2mc%z,mc2mo%fcc)
      case (4)
        call texeira_cldfrac(totc,mo2mc%qs,mo2mc%rh,rh0, &
                             qtcrit,mc2mo%fcc)
      case (5)
        call tompkins_cldfrac(mo2mc%rh,totc,mo2mc%phs,mo2mc%ps2, &
                              qtcrit,mc2mo%fcc)
      case (6)
        call echam5_cldfrac(totc,mo2mc%rh,mo2mc%phs,mo2mc%ps2, &
                            qtcrit,mc2mo%fcc)
      case (7)
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          if ( totc(j,i,k) * mo2mc%delz(j,i,k) > 1.0e-2_rkx ) then
            mc2mo%fcc(j,i,k) = 1.0_rkx
          else
            mc2mo%fcc(j,i,k) = 0.0_rkx
          end if
        end do
      case default
        mc2mo%fcc = d_zero
        if ( myid == italk ) then
          write(stdout,*) ' No cloud fraction algorithm used.'
        end if
    end select

    if ( ipptls == 1 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        rh0adj(j,i,k) = d_one - (d_one-mo2mc%rh(j,i,k)) / &
                        max((d_one-mc2mo%fcc(j,i,k)),1.0e-5_rkx)**2
      end do
    end if

    !------------------------------------------
    ! 1a. Determine Marine stratocumulus clouds
    !------------------------------------------

    if ( icldmstrat == 1 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        if ( mo2mc%iveg(j,i) == 15 ) then
          if ( mo2mc%phs(j,i,k) >= 70000.0_rkx ) then
            ! Klein, S. A., and D. L. Hartmann,
            ! The seasonal cycle of low stratiform clouds,
            ! J. Climate, 6, 1587-1606, 1993
            mc2mo%fcc(j,i,k) = max(mc2mo%fcc(j,i,k), &
                      min(((atms%th700(j,i)-atms%th3d(j,i,k)) * &
                           0.057_rkx) - 0.5573_rkx,1.0_rkx))
          end if
        end if
      end do
    end if

    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      mc2mo%fcc(j,i,k) = min(max(mc2mo%fcc(j,i,k),d_zero),d_one)
    end do

    !-----------------------------------------------------------------
    ! 2.  Combine large-scale and convective fraction and liquid water
    !     to be passed into radiation.
    !-----------------------------------------------------------------

    if ( iconvlwp == 1 ) then
#ifdef STDPAR_FIXED
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
#else
      !$acc parallel loop collapse(3)
      do k = 1, kz
      do i = ici1, ici2
      do j = jci1, jci2
#endif
        ! get overlap of clouds
        cldfra(j,i,k) = rov2(cldfra(j,i,k),mc2mo%fcc(j,i,k))
        if ( cldfra(j,i,k) > lowcld ) then
          ! Cloud Water Volume
          ! Apply the parameterisation based on temperature to the
          ! the large scale clouds. This is an in-cloud here.
          cldlwc(j,i,k) = clwfromt(mo2mc%t(j,i,k))
        else
          cldfra(j,i,k) = d_zero
          cldlwc(j,i,k) = d_zero
        end if
#ifndef STDPAR_FIXED
      end do
      end do
#endif
      end do
    else
      if ( any(icup > 1) ) then
#ifdef STDPAR_FIXED
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
#else
        !$acc parallel loop collapse(3)
        do k = 1, kz
        do i = ici1, ici2
        do j = jci1, jci2
#endif
          ! Cloud Water Volume
          ! kg gq / kg dry air * kg dry air / m3 * 1000 = g qc / m3
          ls_exlwc = (totc(j,i,k)*d_1000)*mo2mc%rho(j,i,k)
          conv_exlwc = clwfromt(mo2mc%t(j,i,k))
          exlwc = conv_exlwc*cldfra(j,i,k) + ls_exlwc
          ! get overlap of clouds
          cldfra(j,i,k) = rov2(cldfra(j,i,k),mc2mo%fcc(j,i,k))
          if ( cldfra(j,i,k) > lowcld ) then
            ! NOTE : IN CLOUD LWC IS NEEDED IN THE RADIATION !!!
            exlwc = exlwc/cldfra(j,i,k)
            ! Scaling for CF
            ! Implements CF scaling as in Liang GRL 32, 2005
            ! doi: 10.1029/2004GL022301
            if ( do_cfscaling ) then
              ichi = int(cldfra(j,i,k)*real(nchi-1,rkx))
              exlwc = exlwc * chis(ichi)
            end if
            cldlwc(j,i,k) = exlwc
          else
            cldfra(j,i,k) = d_zero
            cldlwc(j,i,k) = d_zero
          end if
#ifndef STDPAR_FIXED
        end do
        end do
#endif
        end do
      else
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          ! Cloud Water Volume
          ! kg gq / kg dry air * kg dry air / m3 * 1000 = g qc / m3
          cldfra(j,i,k) = min(max(mc2mo%fcc(j,i,k),d_zero),d_one)
          if ( cldfra(j,i,k) > lowcld ) then
            exlwc = (totc(j,i,k)*d_1000)*mo2mc%rho(j,i,k)
            ! Scaling for CF
            ! Implements CF scaling as in Liang GRL 32, 2005
            ! doi: 10.1029/2004GL022301
            if ( do_cfscaling ) then
              ichi = int(cldfra(j,i,k)*real(nchi-1,rkx))
              exlwc = exlwc * chis(ichi)
            end if
            ! NOTE : IN CLOUD HERE IS NEEDED !!!
            cldlwc(j,i,k) = exlwc/cldfra(j,i,k)
          else
            cldfra(j,i,k) = d_zero
            cldlwc(j,i,k) = d_zero
          end if
        end do
      end if
    end if

    contains

#include <clwfromt.inc>

  end subroutine cldfrac

  pure real(rkx) function max_overlap(a)
    implicit none
    real(rkx), dimension(:), intent(in) :: a
    max_overlap = maxval(a)
  end function max_overlap

  pure real(rkx) function min_overlap(a)
    implicit none
    real(rkx), dimension(:), intent(in) :: a
    min_overlap = minval(a)
  end function min_overlap

  pure real(rkx) function mean_overlap(a)
    implicit none
    real(rkx), dimension(:), intent(in) :: a
    mean_overlap = sum(a)/size(a)
  end function mean_overlap

  pure real(rkx) function random_overlap(a)
    implicit none
    real(rkx), dimension(:), intent(in) :: a
    real(rkx) :: fac
    integer :: i
    fac = d_one
    do i = 1, size(a)
      fac = fac * (d_one-max(min(a(i),d_one),d_zero))
    end do
    random_overlap = d_one - fac
  end function random_overlap

  pure real(rkx) function max_random_overlap(a)
    !$acc routine seq
    implicit none
    real(rkx), dimension(:), intent(in) :: a
    real(rkx) :: fac
    integer :: i
    fac = d_one - a(1)
    do i = 2, size(a)
      fac = fac * (d_one-max(a(i-1),a(i)))/(d_one-a(i-1))
    end do
    max_random_overlap = d_one - fac
  end function max_random_overlap

  elemental real(rkx) function rov2(a,b)
    !$acc routine seq
    implicit none
    real(rkx), intent(in) :: a, b
    rov2 = max_random_overlap([a,b])
  end function rov2
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
        fccc = d_one-sqrt((d_one-rhc)/(d_one-rh0adj(j,i,k)))
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

#include <pfwsat.inc>
#include <wlh.inc>

  end subroutine condtq

end module mod_micro_interface
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
