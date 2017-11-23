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

!
! Link atmospheric model and cumulus schemes
!
module mod_cu_interface
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_runparams
  use mod_memutil
  use mod_regcm_types
  use mod_mppparam , only : exchange , uvcross2dot , italk

  use mod_cu_common , only : cuscheme , total_precip_points , cevapu ,     &
      model_cumulus_cloud , init_mod_cumulus
  use mod_cu_common , only : avg_tten , avg_uten , avg_vten , avg_qten , &
      avg_chiten
  use mod_cu_common , only : cu_uten , cu_vten , cu_tten , cu_qten , &
      cu_prate , cu_ktop , cu_kbot , cu_cldfrc , cu_qdetr , cu_raincc , &
      cu_convpr , cu_chiten , avg_ww
  use mod_cu_tiedtke , only : allocate_mod_cu_tiedtke , tiedtkedrv
  use mod_cu_tables , only : init_convect_tables
  use mod_cu_bm , only : allocate_mod_cu_bm , bmpara , lutbl , cldefi ,    &
      tbase
  use mod_cu_em , only : allocate_mod_cu_em , cupemandrv , cbmf2d ,        &
      elcrit2d , epmax2d
  use mod_cu_kuo , only : allocate_mod_cu_kuo , cupara , twght , vqflx , k700
  use mod_cu_grell , only : allocate_mod_cu_grell , cuparan , mincld2d ,   &
      shrmax2d , shrmin2d , edtmax2d , edtmin2d ,  edtmaxo2d , edtmaxx2d , &
      edtmino2d , edtminx2d , htmax2d , htmin2d , pbcmax2d , dtauc2d ,     &
      pbcmax2d , kbmax2d
  use mod_cu_kf
  use mod_cu_shallow

  implicit none

  private

  public :: allocate_cumulus
  public :: init_cumulus
  public :: cumulus , shallow_convection
  public :: cucloud

  public :: lutbl

  public :: cuscheme
  public :: cbmf2d
  public :: cldefi
  public :: tbase
  public :: avg_ww
  public :: avg_tten
  public :: avg_uten
  public :: avg_vten
  public :: avg_qten
  public :: avg_chiten
  public :: twght
  public :: vqflx
  public :: shrmax2d
  public :: shrmin2d
  public :: edtmax2d
  public :: edtmin2d
  public :: edtmaxo2d
  public :: edtmino2d
  public :: edtmaxx2d
  public :: edtminx2d
  public :: pbcmax2d
  public :: mincld2d
  public :: kbmax2d
  public :: htmax2d
  public :: htmin2d
  public :: dtauc2d
  public :: elcrit2d
  public :: epmax2d
  public :: cevapu
  public :: k700
  public :: total_precip_points

  type(mod_2_cum) :: m2c
  type(cum_2_mod) :: c2m

  contains

  subroutine allocate_cumulus
    use mod_atm_interface
    implicit none
    integer(ik4) :: i , j
    call getmem2d(cuscheme,jci1,jci2,ici1,ici2,'cumulus:cuscheme')
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( mddom%lndcat(j,i) > 14.5_rkx .and. &
             mddom%lndcat(j,i) < 15.5_rkx ) then
          cuscheme(j,i) = icup_ocn
        else
          cuscheme(j,i) = icup_lnd
        end if
      end do
    end do
    if ( any(icup == 1) ) then
      call allocate_mod_cu_kuo
    end if
    if ( any(icup == 2) ) then
      call allocate_mod_cu_grell
    end if
    if ( any(icup == 3) ) then
      call allocate_mod_cu_bm
    end if
    if ( any(icup == 4) ) then
      call allocate_mod_cu_em
    end if
    if ( any(icup == 5) ) then
      if ( iconv /= 4 ) call init_convect_tables
      call allocate_mod_cu_tiedtke
    end if
    if ( any(icup == 6) ) then
      call allocate_mod_cu_kf
      call kf_lutab
    end if
  end subroutine allocate_cumulus

  subroutine init_cumulus
    use mod_atm_interface
    use mod_che_interface
    implicit none
    ! INPUT
    call assignpnt(mddom%ht,m2c%ht)
    call assignpnt(mddom%ldmsk,m2c%ldmsk)
    call assignpnt(sfs%psb,m2c%psb)
    call assignpnt(sfs%psa,m2c%psa)
    call assignpnt(sfs%psdotb,m2c%psdotb)
    call assignpnt(atms%za,m2c%zas)
    call assignpnt(atms%zq,m2c%zfs)
    call assignpnt(atms%dzq,m2c%dzq)
    call assignpnt(atms%ps2d,m2c%psf)
    call assignpnt(atms%pb3d,m2c%pas)
    call assignpnt(atms%pf3d,m2c%pasf)
    call assignpnt(atms%tb3d,m2c%tas)
    call assignpnt(atms%ubx3d,m2c%uas)
    call assignpnt(atms%vbx3d,m2c%vas)
    call assignpnt(atms%wpx3d,m2c%wpas)
    call assignpnt(atms%wb3d,m2c%was)
    call assignpnt(atm2%tke,m2c%tkeas) ! Not coupled here.
    call assignpnt(atms%qsb3d,m2c%qsas)
    call assignpnt(atms%qxb3d,m2c%qxas)
    call assignpnt(atms%rhob3d,m2c%rhoas)
    call assignpnt(atms%chib3d,m2c%chias)
    call assignpnt(qdot,m2c%qdot)
    call assignpnt(atm1%qx,m2c%qq1,iqv)
    call assignpnt(sfs%qfx,m2c%qfx)
    call assignpnt(sfs%hfx,m2c%hfx)
    call assignpnt(ktrop,m2c%ktrop)
    call assignpnt(ccn,m2c%ccn)
    call assignpnt(aten%t,m2c%tten,pc_physic)
    call assignpnt(aten%qx,m2c%qxten,pc_physic)
    call assignpnt(aten%u,m2c%uten,pc_physic)
    call assignpnt(aten%v,m2c%vten,pc_physic)
    call assignpnt(aten%chi,m2c%chiten,pc_physic)
    call assignpnt(heatrt,m2c%heatrt)
    ! OUTPUT
    call assignpnt(aten%t,c2m%tten,pc_physic)
    call assignpnt(aten%u,c2m%uten,pc_physic)
    call assignpnt(aten%v,c2m%vten,pc_physic)
    call assignpnt(aten%qx,c2m%qxten,pc_physic)
    call assignpnt(aten%chi,c2m%chiten,pc_physic)
    call assignpnt(sfs%rainc,c2m%rainc)
    call assignpnt(pptc,c2m%pcratec)
    call assignpnt(cldfra,c2m%cldfrc)
    call assignpnt(cldlwc,c2m%cldlwc)
    call assignpnt(icumtop,c2m%kcumtop)
    call assignpnt(icumbot,c2m%kcumbot)
    call assignpnt(convpr,c2m%convpr)
    call assignpnt(q_detr,c2m%q_detr)
    call assignpnt(rain_cc,c2m%rain_cc)
    call init_mod_cumulus
  end subroutine init_cumulus

  subroutine cucloud
    implicit none
    integer(ik4) :: i , j , k
    if ( any(icup == 1) .or. any(icup == 3) ) then
      call model_cumulus_cloud(m2c)
    end if
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          c2m%cldfrc(j,i,k) = max(cu_cldfrc(j,i,k),0.0_rkx)
          if ( cu_cldfrc(j,i,k) > 0.001_rkx ) then
            c2m%cldlwc(j,i,k) = clwfromt(m2c%tas(j,i,k))
          else
            c2m%cldlwc(j,i,k) = d_zero
          end if
        end do
      end do
    end do

    contains

#include <clwfromt.inc>

  end subroutine cucloud

  subroutine cumulus
    implicit none
    integer(ik4) :: i , j , k , n
    real(rkx) :: w1

    if ( any(icup == 6) ) then
      w1 = d_one/real(max(int(900.0_rkx/dtsec),1),rkx)
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            avg_ww(j,i,k) = (d_one-w1) * avg_ww(j,i,k) + &
                          w1 * d_half * (m2c%was(j,i,k)+m2c%was(j,i,k+1))
          end do
        end do
      end do
    end if

    w1 = d_one/real(max(int(dtcum/dtsec),1),rkx)

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          avg_tten(j,i,k) = (d_one-w1) * avg_tten(j,i,k) + &
                        w1 * m2c%tten(j,i,k)/m2c%psb(j,i)+m2c%heatrt(j,i,k)
        end do
      end do
    end do

    if ( any(icup == 5) ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            avg_ww(j,i,k) = (d_one-w1)*avg_ww(j,i,k) + w1*m2c%wpas(j,i,k)
          end do
        end do
      end do
      call exchange(m2c%uten,1,jdi1,jdi2,idi1,idi2,1,kz)
      call exchange(m2c%vten,1,jdi1,jdi2,idi1,idi2,1,kz)
      do k = 1 , kz
        do i = icii1 , icii2
          do j = jcii1 , jcii2
            avg_uten(j,i,k) = (d_one-w1) * avg_uten(j,i,k) + &
                w1 * d_rfour * &
                    (m2c%uten(j,i,k)   + m2c%uten(j+1,i,k) + &
                     m2c%uten(j,i+1,k) + m2c%uten(j+1,i+1,k)) / m2c%psb(j,i)
            avg_vten(j,i,k) = (d_one-w1) * avg_vten(j,i,k) + &
                w1 * d_rfour * &
                    (m2c%vten(j,i,k)   + m2c%vten(j+1,i,k) + &
                     m2c%vten(j,i+1,k) + m2c%vten(j+1,i+1,k)) / m2c%psb(j,i)
          end do
        end do
      end do
    end if
    do n = 1 , nqx
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            avg_qten(j,i,k,n) = (d_one-w1) * avg_qten(j,i,k,n) + &
                          w1 * m2c%qxten(j,i,k,n)/m2c%psb(j,i)
          end do
        end do
      end do
    end do
    if ( ichem == 1 ) then
      do n = 1 , ntr
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              avg_chiten(j,i,k,n) = (d_one-w1) * avg_chiten(j,i,k,n) + &
                             w1 * m2c%chiten(j,i,k,n)/m2c%psb(j,i)
            end do
          end do
        end do
      end do
    end if

    if ( rcmtimer%integrating( ) ) then

      if ( syncro_cum%act( ) ) then

        if ( debug_level > 3 .and. myid == italk ) then
          write(stdout,*) 'Calling cumulus scheme at ',trim(rcmtimer%str())
        end if
        ! Update input cumulus tendencies

        cu_prate(:,:) = d_zero
        cu_ktop(:,:) = 0
        cu_kbot(:,:) = 0
        cu_tten(:,:,:) = d_zero
        cu_uten(:,:,:) = d_zero
        cu_vten(:,:,:) = d_zero
        cu_qten(:,:,:,:) = d_zero
        cu_cldfrc(:,:,:) = d_zero
        if ( ichem == 1 ) then
          cu_chiten(:,:,:,:) = d_zero
          cu_convpr(:,:,:) = d_zero
        end if
        if ( any(icup == 5) ) then
          cu_qdetr(:,:,:) = d_zero
          cu_raincc(:,:,:) = d_zero
        end if

        total_precip_points = 0
        if ( icup_lnd == icup_ocn ) then
          select case ( icup_lnd )
            case (1)
              call cupara(m2c)
            case (2)
              call cuparan(m2c)
            case (3)
              call bmpara(m2c)
            case (4)
              call cupemandrv(m2c)
            case (5)
              call tiedtkedrv(m2c)
            case (6)
              call kfdrv(m2c)
          end select
        else
          select case ( icup_lnd )
            case (2)
              call cuparan(m2c)
            case (4)
              call cupemandrv(m2c)
            case (5)
              call tiedtkedrv(m2c)
            case (6)
              call kfdrv(m2c)
          end select
          select case ( icup_ocn )
            case (2)
              call cuparan(m2c)
            case (4)
              call cupemandrv(m2c)
            case (5)
              call tiedtkedrv(m2c)
            case (6)
              call kfdrv(m2c)
          end select
        end if

        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              cu_uten(j,i,k) = cu_uten(j,i,k) * m2c%psb(j,i)
              cu_vten(j,i,k) = cu_vten(j,i,k) * m2c%psb(j,i)
            end do
          end do
        end do

      end if

      ! Sum cumulus tendencies

      do i = ici1 , ici2
        do j = jci1 , jci2
          c2m%pcratec(j,i) = c2m%pcratec(j,i) + cu_prate(j,i)
          c2m%rainc(j,i) = c2m%rainc(j,i) + cu_prate(j,i) * dtsec
          c2m%kcumtop(j,i) = cu_ktop(j,i)
          c2m%kcumbot(j,i) = cu_kbot(j,i)
        end do
      end do

      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            c2m%tten(j,i,k) = c2m%tten(j,i,k) + cu_tten(j,i,k) * m2c%psb(j,i)
          end do
        end do
      end do

      call uvcross2dot(cu_uten,cu_vten,c2m%uten,c2m%vten)

      do n = 1 , nqx
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              c2m%qxten(j,i,k,n) = c2m%qxten(j,i,k,n) + &
                       cu_qten(j,i,k,n) * m2c%psb(j,i)
            end do
          end do
        end do
      end do

      if ( ichem == 1 ) then
        do n = 1 , ntr
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                c2m%chiten(j,i,k,n) = c2m%chiten(j,i,k,n) + &
                                 cu_chiten(j,i,k,n) * m2c%psb(j,i)
              end do
            end do
          end do
        end do
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              c2m%convpr(j,i,k) = cu_convpr(j,i,k)
            end do
          end do
        end do
      end if

      if ( any(icup == 5) ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              c2m%q_detr(j,i,k) = cu_qdetr(j,i,k)
              c2m%rain_cc(j,i,k) = cu_raincc(j,i,k)
            end do
          end do
        end do
        if ( ipptls /= 2 ) then
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                c2m%qxten(j,i,k,iqv) = c2m%qxten(j,i,k,iqv) + &
                                       c2m%q_detr(j,i,k) * m2c%psb(j,i) / dt
              end do
            end do
          end do
        end if
      end if

    end if

  end subroutine cumulus

  subroutine shallow_convection
    use mod_atm_interface , only : aten
    implicit none
    integer(ik4) :: i , j , k

    if ( rcmtimer%integrating( ) ) then

      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            avg_tten(j,i,k) = aten%t(j,i,k,pc_total)/m2c%psb(j,i)
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            avg_qten(j,i,k,iqv) = aten%qx(j,i,k,iqv,pc_total) / m2c%psb(j,i)
          end do
        end do
      end do

      call shallcu(m2c)

      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            aten%t(j,i,k,pc_total) = aten%t(j,i,k,pc_total) + &
                        cu_tten(j,i,k) * m2c%psb(j,i)
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            aten%qx(j,i,k,iqv,pc_total) = aten%qx(j,i,k,iqv,pc_total) + &
                        cu_qten(j,i,k,iqv) * m2c%psb(j,i)
          end do
        end do
      end do

    end if
  end subroutine shallow_convection

end module mod_cu_interface

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
