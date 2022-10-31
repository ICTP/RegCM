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
  use mod_constants
  use mod_stdio
  use mod_dynparam
  use mod_runparams
  use mod_memutil
  use mod_regcm_types
  use mod_mppparam , only : exchange_lrbt , uvcross2dot , uvdot2cross , italk
  use mod_mppparam , only : tenxtouvten , uvtentotenx
  use mod_mppparam , only : meanall , minall , maxall

  use mod_cu_common , only : cuscheme , total_precip_points , &
      model_cumulus_cloud , init_mod_cumulus
  use mod_cu_common , only : cu_uten , cu_vten , cu_tten , cu_qten , &
      cu_prate , cu_ktop , cu_kbot , cu_cldfrc , cu_qdetr , cu_raincc , &
      cu_convpr , cu_chiten , avg_ww
  use mod_cu_tiedtke , only : allocate_mod_cu_tiedtke , tiedtkedrv , &
      pmean , nmctop
  use mod_cu_tables , only : init_convect_tables
  use mod_cu_bm , only : allocate_mod_cu_bm , bmpara , cldefi
  use mod_cu_em , only : allocate_mod_cu_em , cupemandrv , cbmf2d
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

  public :: cuscheme
  public :: cbmf2d
  public :: avg_ww
  public :: cldefi
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
  public :: k700
  public :: total_precip_points

  type(mod_2_cum) :: m2c
  type(cum_2_mod) :: c2m

  real(rkx) , pointer , dimension(:,:,:) :: utenx , vtenx
  real(rkx) , pointer , dimension(:,:,:) :: utend , vtend

  ! Midlevel convection top pressure for Tiedtke iconv = 1
  real(rkx) , parameter :: cmcptop = 30000.0_rkx

  contains

  subroutine allocate_cumulus
    use mod_atm_interface
    implicit none
    integer(ik4) :: i , j
    call getmem2d(cuscheme,jci1,jci2,ici1,ici2,'cumulus:cuscheme')
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( isocean(mddom%lndcat(j,i)) ) then
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
    if ( any(icup == 4) .or. any(icup == 5) ) then
      if ( idynamic == 3 ) then
        call getmem3d(utend,jdi1gb,jdi2gb,ici1,ici2,1,kz,'pbl_common:utend')
        call getmem3d(vtend,jci1,jci2,idi1gb,idi2gb,1,kz,'pbl_common:vtend')
      else
        call getmem3d(utend,jdi1ga,jdi2ga,idi1ga,idi2ga,1,kz,'pbl_common:utend')
        call getmem3d(vtend,jdi1ga,jdi2ga,idi1ga,idi2ga,1,kz,'pbl_common:vtend')
      end if
    end if
    if ( any(icup == 4) ) then
      call allocate_mod_cu_em
    end if
    if ( any(icup == 5) ) then
      if ( iconv /= 4 ) call init_convect_tables
      call allocate_mod_cu_tiedtke
      call getmem3d(utenx,jci1,jci2,ici1,ici2,1,kz,'pbl_common:utenx')
      call getmem3d(vtenx,jci1,jci2,ici1,ici2,1,kz,'pbl_common:vtenx')
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
    if ( idynamic == 3 ) then
      call assignpnt(mo_atm%tke,m2c%tkeas)
      call assignpnt(mo_atm%qx,m2c%qq1,iqv)
    else
      call assignpnt(atm2%tke,m2c%tkeas) ! Not coupled here.
      call assignpnt(atm1%qx,m2c%qq1,iqv)
    end if
    call assignpnt(atms%qsb3d,m2c%qsas)
    call assignpnt(atms%qxb3d,m2c%qxas)
    call assignpnt(atms%rhob3d,m2c%rhoas)
    call assignpnt(atms%chib3d,m2c%chias)
    call assignpnt(qdot,m2c%qdot)
    call assignpnt(sfs%qfx,m2c%qfx)
    call assignpnt(sfs%hfx,m2c%hfx)
    call assignpnt(ktrop,m2c%ktrop)
    call assignpnt(ccn,m2c%ccn)
    call assignpnt(heatrt,m2c%heatrt)
    ! Pass dynamic tendencies as input.
    if ( idynamic == 3 ) then
      call assignpnt(mo_atm%tten,m2c%tten)
      call assignpnt(mo_atm%qxten,m2c%qxten)
      call assignpnt(mo_atm%uten,m2c%uten)
      call assignpnt(mo_atm%vten,m2c%vten)
      if ( ichem == 1 ) call assignpnt(mo_atm%chiten,m2c%chiten)
      call assignpnt(mo_atm%tten,c2m%tten)
      call assignpnt(mo_atm%uten,c2m%uten)
      call assignpnt(mo_atm%vten,c2m%vten)
      call assignpnt(mo_atm%qxten,c2m%qxten)
      if ( ichem == 1 ) call assignpnt(mo_atm%chiten,c2m%chiten)
    else
      call assignpnt(aten%t,m2c%tten,pc_physic)
      call assignpnt(aten%qx,m2c%qxten,pc_physic)
      call assignpnt(aten%qx,m2c%dynqx,pc_dynamic)
      call assignpnt(aten%u,m2c%uten,pc_physic)
      call assignpnt(aten%v,m2c%vten,pc_physic)
      if ( ichem == 1 ) call assignpnt(aten%chi,m2c%chiten,pc_physic)
      ! OUTPUT
      call assignpnt(aten%t,c2m%tten,pc_physic)
      call assignpnt(aten%u,c2m%uten,pc_physic)
      call assignpnt(aten%v,c2m%vten,pc_physic)
      call assignpnt(aten%qx,c2m%qxten,pc_physic)
      if ( ichem == 1 ) call assignpnt(aten%chi,c2m%chiten,pc_physic)
    end if
    call assignpnt(sfs%rainc,c2m%rainc)
    call assignpnt(pptc,c2m%pcratec)
    call assignpnt(cldfra,c2m%cldfrc)
    call assignpnt(cldlwc,c2m%cldlwc)
    call assignpnt(icumtop,c2m%kcumtop)
    call assignpnt(icumbot,c2m%kcumbot)
    if ( ichem == 1 ) call assignpnt(convpr,c2m%convpr)
    call assignpnt(q_detr,c2m%q_detr)
    call assignpnt(rain_cc,c2m%rain_cc)
    call assignpnt(crrate,c2m%trrate)
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
    use mod_atm_interface
    implicit none
    integer(ik4) :: i , j , k , n
    integer(ik4) :: iplmlc , mintop , maxtop
    real(rkx) :: mymean , w1

    if ( any(icup == 6) ) then
      w1 = d_one/real(max(int(max(dtcum,900.0_rkx)/dtsec),1),rkx)
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( cuscheme(j,i) == 6 ) then
              avg_ww(j,i,k) = (d_one - w1) * avg_ww(j,i,k) + &
                            w1 * d_half * (m2c%was(j,i,k)+m2c%was(j,i,k+1))
            end if
          end do
        end do
      end do
    end if

    if ( any(icup == 5) ) then
      if ( iconv == 1 ) then
        ! Calculate average elevation of cmcptop level
        nmctop = 0
        do i = ici1 , ici2
          do j = jci1 , jci2
            iplmlc = 1
            do k = 1 , kzp1
              iplmlc = k
              if ( m2c%pasf(j,i,k) >= cmcptop ) exit
            end do
            nmctop = nmctop + iplmlc
          end do
        end do
        iplmlc = nmctop / ((jci2-jci1)*(ici2-ici1))
        call minall(iplmlc,mintop)
        call maxall(iplmlc,maxtop)
        nmctop = (mintop+maxtop)/2
      else if ( iconv == 4 ) then
        do k = 1 , kz
          mymean = sum(m2c%pas(:,:,k))/real(((jci2-jci1)*(ici2-ici1)),rkx)
          call meanall(mymean,pmean(k))
        end do
      end if

      w1 = d_one/real(max(int(max(dtcum,300.0_rkx)/dtsec),1),rkx)
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( cuscheme(j,i) == 5 ) then
              avg_ww(j,i,k) = (d_one-w1)*avg_ww(j,i,k) + w1*m2c%wpas(j,i,k)
            end if
          end do
        end do
      end do
    end if

    ! Skip first timestep

    if ( rcmtimer%integrating( ) ) then

      if ( syncro_cum%act( ) ) then

        if ( debug_level > 3 .and. myid == italk ) then
          write(stdout,*) 'Calling cumulus scheme at ',trim(rcmtimer%str())
        end if

        cu_prate(:,:) = d_zero
        cu_ktop(:,:) = 0
        cu_kbot(:,:) = 0
        cu_tten(:,:,:) = d_zero
        if ( any(icup == 4) .or. any(icup == 5) ) then
          cu_uten(:,:,:) = d_zero
          cu_vten(:,:,:) = d_zero
          if ( any(icup == 5) ) then
            if ( idynamic == 3 ) then
              utend(jdi1:jdi2,ici1:ici2,:) = m2c%uten
              vtend(jci1:jci2,idi1:idi2,:) = m2c%vten
              call uvtentotenx(utend,vtend,utenx,vtenx)
            else
              utend(jdi1:jdi2,idi1:idi2,:) = m2c%uten
              vtend(jdi1:jdi2,idi1:idi2,:) = m2c%vten
              call uvdot2cross(utend,vtend,utenx,vtenx)
            end if
          end if
          utend = d_zero
          vtend = d_zero
        end if
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
              call tiedtkedrv(m2c,utenx,vtenx)
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
              call tiedtkedrv(m2c,utenx,vtenx)
            case (6)
              call kfdrv(m2c)
          end select
          select case ( icup_ocn )
            case (2)
              call cuparan(m2c)
            case (4)
              call cupemandrv(m2c)
            case (5)
              call tiedtkedrv(m2c,utenx,vtenx)
            case (6)
              call kfdrv(m2c)
          end select
        end if

        ! Update output wind cumulus tendencies (cross to dot points)

        if ( any(icup == 5) .or. any(icup == 4) ) then
          if ( idynamic == 3 ) then
            call tenxtouvten(cu_uten,cu_vten,utend,vtend)
          else
            call uvcross2dot(cu_uten,cu_vten,utend,vtend)
          end if
        end if

        if ( ipptls == 2 ) then
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( cuscheme(j,i) == 5 ) then
                  cu_qten(j,i,k,iqc) = d_zero
                  cu_qten(j,i,k,iqi) = d_zero
                end if
              end do
            end do
          end do
        end if

      end if

      ! Sum cumulus tendencies

      do i = ici1 , ici2
        do j = jci1 , jci2
          c2m%pcratec(j,i) = c2m%pcratec(j,i) + cu_prate(j,i)
          c2m%trrate(j,i) = cu_prate(j,i)
          c2m%rainc(j,i) = c2m%rainc(j,i) + cu_prate(j,i) * dtsec
          c2m%kcumtop(j,i) = cu_ktop(j,i)
          c2m%kcumbot(j,i) = cu_kbot(j,i)
        end do
      end do

      if ( idynamic == 3 ) then

        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              c2m%tten(j,i,k) = c2m%tten(j,i,k) + cu_tten(j,i,k)
            end do
          end do
        end do

        if ( any(icup == 5) .or. any(icup == 4) ) then
          do k = 1 , kz
            do i = idi1 , idi2
              do j = jci1 , jci2
                c2m%vten(j,i,k) = c2m%vten(j,i,k) + vtend(j,i,k)
              end do
            end do
          end do
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jdi1 , jdi2
                c2m%uten(j,i,k) = c2m%uten(j,i,k) + utend(j,i,k)
              end do
            end do
          end do
        end if

        do n = 1 , nqx
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                c2m%qxten(j,i,k,n) = c2m%qxten(j,i,k,n) + cu_qten(j,i,k,n)
              end do
            end do
          end do
        end do

        if ( ichem == 1 ) then
          do n = 1 , ntr
            do k = 1 , kz
              do i = ici1 , ici2
                do j = jci1 , jci2
                  c2m%chiten(j,i,k,n) = &
                      c2m%chiten(j,i,k,n) + cu_chiten(j,i,k,n)
                end do
              end do
            end do
          end do
        end if

      else

        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              c2m%tten(j,i,k) = c2m%tten(j,i,k) + cu_tten(j,i,k) * m2c%psb(j,i)
            end do
          end do
        end do

        if ( any(icup == 5) .or. any(icup == 4) ) then
          do k = 1 , kz
            do i = idi1 , idi2
              do j = jdi1 , jdi2
                c2m%uten(j,i,k) = c2m%uten(j,i,k) + utend(j,i,k)*m2c%psdotb(j,i)
                c2m%vten(j,i,k) = c2m%vten(j,i,k) + vtend(j,i,k)*m2c%psdotb(j,i)
              end do
            end do
          end do
        end if

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
        end if
      end if

      if ( ichem == 1 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              c2m%convpr(j,i,k) = cu_convpr(j,i,k)
            end do
          end do
        end do
      end if

      if ( any(icup == 5) ) then
        if ( ipptls == 2 ) then
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                c2m%q_detr(j,i,k) = cu_qdetr(j,i,k) * dt
              end do
            end do
          end do
        end if
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              c2m%rain_cc(j,i,k) = cu_raincc(j,i,k)
            end do
          end do
        end do
      end if
    end if
  end subroutine cumulus

  subroutine shallow_convection
    use mod_atm_interface , only : aten , mo_atm
    implicit none
    integer(ik4) :: i , j , k

    if ( rcmtimer%integrating( ) ) then

      call shallcu(m2c)

      if ( idynamic == 3 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              mo_atm%tten(j,i,k) = cu_tten(j,i,k)
            end do
          end do
        end do
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              mo_atm%qxten(j,i,k,iqv) = cu_qten(j,i,k,iqv)
            end do
          end do
        end do
      else
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

    end if
  end subroutine shallow_convection

end module mod_cu_interface

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
