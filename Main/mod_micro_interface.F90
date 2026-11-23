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

  use mod_realkinds, only : rk8, rkx
  use mod_stdio, only : stderr, stdout
  use mod_constants, only : cpd, cpv, d_1000, d_one, d_zero, dlowval, ep2, &
    lowcld, minqq, rtice, rtwat, rtwat_rtice_r, rwat, tzero, wlhs, wlhv
  use mod_dynparam, only : ds, ici1, ici2, idynamic, jci1, jci2, kz, myid, &
    ptop
  use mod_memutil, only : assignpnt, getmem
  use mod_mppparam, only : italk
  use mod_regcm_types, only : ik4, micro_2_mod, mod_2_micro, nogtom_stats
  use mod_runparams, only : dt, hsigma, iaerosol, ichem, icldfrac, &
    icldmstrat, icup, idiag, iindirect, ipptls, iqc, iqi, iqr, iqs, &
    iqv, pc_physic, pc_total, rcmtimer, rhmax, rhmin, stats
  use mod_micro_nogtom, only : allocate_mod_nogtom, init_nogtom, nogtom
  use mod_micro_wsm5, only : allocate_mod_wsm5, init_wsm5, wsm5
  use mod_micro_wsm7, only : allocate_mod_wsm7, init_wsm7, wsm7
  use mod_micro_wdm7, only : allocate_mod_wdm7, init_wdm7, wdm7
  use mod_cloud_sundqvist, only : sundqvist_cldfrac
  use mod_cloud_xuran, only : xuran_cldfrac
  use mod_cloud_thomp, only : thomp_cldfrac
  use mod_cloud_guli2007, only : gulisa_cldfrac
  use mod_cloud_texeira, only : texeira_cldfrac
  use mod_cloud_tompkins, only : tompkins_cldfrac
  use mod_cloud_echam5, only : echam5_cldfrac

  implicit none

  private

  public :: allocate_micro, init_micro, microscheme, cldfrac

  type(nogtom_stats), public :: ngs

  type(mod_2_micro) :: mo2mc
  type(micro_2_mod) :: mc2mo

  ! Threshold relative humidity of the Sundqvist cloud fraction
  real(rkx), pointer, contiguous, dimension(:,:) :: rh0 => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: rhcrit => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: totc => null( )

  integer(ik4), parameter :: nchi = 256
  ! Allocate through getmem so GPU stdpar can access the table in managed mode.
  real(rkx), pointer, contiguous, dimension(:) :: chis => null()
  real(rkx), parameter :: xchi = nchi-1
  real(rkx), parameter :: rchi = 1.0_rkx/xchi

  logical, parameter :: do_cfscaling = .true.

  public :: rh0

  contains

  subroutine allocate_micro
    implicit none
    integer(ik4) :: i
    real(rkx) :: cf
    if ( ipptls == 2 ) then
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
    if ( icldfrac == 0 ) then
      call getmem(rhcrit,jci1,jci2,ici1,ici2,1,kz,'micro:rhcrit')
    end if
    call getmem(totc,jci1,jci2,ici1,ici2,1,kz,'micro:totc')
    call getmem(chis,0,nchi-1,'micro:chis')
    do i = 1, nchi
      cf = real(i-1,rkx)*rchi
      chis(i-1) = 0.97_rkx*exp(-((cf-0.098_rkx)**2)/0.0365_rkx)+0.255_rkx
    end do
  end subroutine allocate_micro

  subroutine init_micro
    use mod_atm_interface, only : aten, atms, ccn, fcc, heatrt, mddom, &
      mo_atm, ncrrate, pptnc, q_detr, qdiag, rain_ls, rembc, remrat, sfs
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
    call assignpnt(atms%qxb3d,mo2mc%qin,iqi)
    call assignpnt(atms%qxb3d,mo2mc%qsn,iqs)
    call assignpnt(atms%qxb3d,mo2mc%qrn,iqr)

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
    real(rkx), pointer, contiguous, &
      dimension(:,:,:), intent(inout) :: cldlwc, cldfra
    real(rkx) :: w1, w2, liq, ice, fcc_ls, cwc
    integer(ik4) :: i, j, k
    integer(ik4) :: ichi

    if ( icldfrac == 3 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        liq = mo2mc%qcn(j,i,k)
        ice = mo2mc%qin(j,i,k)
        if ( liq < 2.0e-7_rkx ) liq = 0.0_rkx
        if ( ice < 2.0e-7_rkx ) ice = 0.0_rkx
        totc(j,i,k) = liq + ice + mo2mc%qrn(j,i,k) + mo2mc%qsn(j,i,k)
      end do
    else
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        liq = mo2mc%qcn(j,i,k)
        ice = mo2mc%qin(j,i,k)
        if ( liq < 2.0e-7_rkx ) liq = 0.0_rkx
        if ( ice < 2.0e-7_rkx ) ice = 0.0_rkx
        totc(j,i,k) = liq + ice
      end do
    end if

    select case ( icldfrac )
      case (0)
        call sundqvist_cldfrac(mo2mc%t,mo2mc%phs,mo2mc%qvn, &
                           mo2mc%rh,rh0,rhcrit,mc2mo%fcc)
      case (1)
        call xuran_cldfrac(totc,mo2mc%qs,mo2mc%rh,mc2mo%fcc)
      case (2)
        call thomp_cldfrac(mo2mc%phs,mo2mc%t,mo2mc%qvn,totc, &
                           mo2mc%qsn,mo2mc%qin,mo2mc%ldmsk,  &
                           ds,mc2mo%fcc)
      case (3)
        call gulisa_cldfrac(totc,mo2mc%z,mc2mo%fcc)
      case (4)
        call texeira_cldfrac(totc,mo2mc%qs,mo2mc%rh,rh0,mc2mo%fcc)
      case (5)
        call tompkins_cldfrac(mo2mc%rh,totc,mo2mc%phs,mo2mc%ps2,mc2mo%fcc)
      case (6)
        call echam5_cldfrac(totc,mo2mc%rh,mo2mc%phs,mo2mc%ps2,mc2mo%fcc)
      case (7)
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          if ( totc(j,i,k)*mo2mc%rho(j,i,k) * &
               mo2mc%delz(j,i,k) > 1.0e-2_rkx ) then
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

    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      ! Transform in g/m^3
      ! kg gq / kg dry air * kg dry air / m3 * 1000 = g qc / m3
      totc(j,i,k) = totc(j,i,k) * 1000.0_rkx * mo2mc%rho(j,i,k)
    end do

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
        ! In cloud convective condensate: from the cumulus scheme if it
        ! provides it (cucloud stores it in cldlwc), else a function of
        ! the temperature.
        cwc = cldlwc(j,i,k)
        if ( cwc <= d_zero ) cwc = clwfromt(mo2mc%t(j,i,k))
        ! get overlap of clouds
        if ( cldfra(j,i,k) > lowcld ) then
          if ( mc2mo%fcc(j,i,k) > lowcld ) then
            fcc_ls = mc2mo%fcc(j,i,k)
            cldfra(j,i,k) = max(fcc_ls, cldfra(j,i,k))
            w2 = fcc_ls / cldfra(j,i,k)
            w1 = 1.0_rkx - w2
            cldlwc(j,i,k) = w2 * totc(j,i,k)/fcc_ls + w1 * cwc
          else
            cldlwc(j,i,k) = cwc
          end if
        else
          if ( mc2mo%fcc(j,i,k) > lowcld ) then
            cldfra(j,i,k) = mc2mo%fcc(j,i,k)
            cldlwc(j,i,k) = totc(j,i,k)/mc2mo%fcc(j,i,k)
          else
            cldfra(j,i,k) = 0.0_rkx
            cldlwc(j,i,k) = 0.0_rkx
          end if
        end if
        ! Scaling for CF
        ! Implements CF scaling as in Liang GRL 32, 2005
        ! doi: 10.1029/2004GL022301
        if ( cldlwc(j,i,k) > 0.0_rkx .and. do_cfscaling ) then
          ichi = int(cldfra(j,i,k)*xchi)
          cldlwc(j,i,k) = cldlwc(j,i,k) * chis(ichi)
        end if
#ifndef STDPAR_FIXED
      end do
      end do
#endif
      end do
    else
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        ! Cloud Water Volume
        cldfra(j,i,k) = mc2mo%fcc(j,i,k)
        ! kg gq / kg dry air * kg dry air / m3 * 1000 = g qc / m3
        if ( cldfra(j,i,k) > lowcld ) then
          ! NOTE : IN CLOUD HERE IS NEEDED !!!
          cldlwc(j,i,k) = totc(j,i,k)/cldfra(j,i,k)
          ! Scaling for CF
          ! Implements CF scaling as in Liang GRL 32, 2005
          ! doi: 10.1029/2004GL022301
          if ( do_cfscaling ) then
            ichi = int(cldfra(j,i,k)*xchi)
            cldlwc(j,i,k) = cldlwc(j,i,k) * chis(ichi)
          end if
        else
          cldfra(j,i,k) = d_zero
          cldlwc(j,i,k) = d_zero
        end if
      end do
    end if

    contains

#include <clwfromt.inc>

  end subroutine cldfrac

end module mod_micro_interface
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
