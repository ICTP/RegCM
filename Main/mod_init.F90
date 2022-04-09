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

module mod_init
  !
  ! RegCM Init module
  !
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_stdio
  use mod_runparams
  use mod_mppparam
  use mod_lm_interface
  use mod_atm_interface
  use mod_che_interface
  use mod_cu_interface
  use mod_pbl_interface
  use mod_rad_interface
  use mod_slabocean
  use rrtmg_sw_init
  use rrtmg_lw_init
  use mod_pbl_interface
  use mod_diffusion , only : initialize_diffusion
  use mod_micro_interface
  use mod_bdycod
  use mod_mpmessage
  use mod_sun
  use mod_ncio
  use mod_savefile
  use mod_slice
  use mod_constants
  use mod_outvars
  use mod_service
  use mod_massck
  use mod_zita
  use mod_sound , only : init_sound
  use mod_moloch , only : init_moloch

  implicit none

  private

  public :: init

  real(rkx) , parameter :: tlp = 50.0_rkx
  real(rkx) , parameter :: ts00 = 288.0_rkx

  contains

#include <pfesat.inc>
#include <pfwsat.inc>

  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !                                                                   c
  !  This subroutine reads in the initial and boundary conditions.    c
  !                                                                   c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine init
    implicit none
    integer(ik4) :: i , j , k , n
    real(rkx) :: rdnnsg
    real(rkx) :: zzi , zfilt
    real(rkx) , dimension(kzp1) :: ozprnt
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'init'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! For an initial run -- not a restart
    !
    if ( .not. ifrest ) then
      if ( irceideal == 1 ) then
        call initideal( )
      end if
      !
      ! Initialize model atmospheric status variables
      ! Data are from the ICBC input at first timestep.
      !
      if ( idynamic < 3 ) then
        do k = 1 , kz
          do i = ide1 , ide2
            do j = jde1 , jde2
              atm1%u(j,i,k) = xub%b0(j,i,k)
              atm1%v(j,i,k) = xvb%b0(j,i,k)
              atm2%u(j,i,k) = xub%b0(j,i,k)
              atm2%v(j,i,k) = xvb%b0(j,i,k)
            end do
          end do
        end do
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              atm1%t(j,i,k) = xtb%b0(j,i,k)
              atm2%t(j,i,k) = xtb%b0(j,i,k)
              atm1%qx(j,i,k,iqv) = xqb%b0(j,i,k)
              atm2%qx(j,i,k,iqv) = xqb%b0(j,i,k)
            end do
          end do
        end do
        if ( ipptls > 0 ) then
          if ( is_present_qc( ) ) then
            do k = 1 , kz
              do i = ice1 , ice2
                do j = jce1 , jce2
                  atm1%qx(j,i,k,iqc) = xlb%b0(j,i,k)
                  atm2%qx(j,i,k,iqc) = xlb%b0(j,i,k)
                end do
              end do
            end do
          end if
          if ( ipptls > 1 ) then
            if ( is_present_qi( ) ) then
              do k = 1 , kz
                do i = ice1 , ice2
                  do j = jce1 , jce2
                    atm1%qx(j,i,k,iqi) = xib%b0(j,i,k)
                    atm2%qx(j,i,k,iqi) = xib%b0(j,i,k)
                  end do
                end do
              end do
            end if
          end if
        end if
        if ( idynamic == 1 ) then
          do i = ice1 , ice2
            do j = jce1 , jce2
              sfs%psa(j,i) = xpsb%b0(j,i)
              sfs%psb(j,i) = xpsb%b0(j,i)
            end do
          end do
        else
          do i = ice1 , ice2
            do j = jce1 , jce2
              sfs%psa(j,i) = atm0%ps(j,i) * d_r1000
              sfs%psb(j,i) = sfs%psa(j,i)
              sfs%psc(j,i) = sfs%psa(j,i)
            end do
          end do
          do k = 1 , kz
            do i = ice1 , ice2
              do j = jce1 , jce2
                atm1%pp(j,i,k) = xppb%b0(j,i,k)
                atm2%pp(j,i,k) = xppb%b0(j,i,k)
              end do
            end do
          end do
          do k = 1 , kzp1
            do i = ice1 , ice2
              do j = jce1 , jce2
                atm1%w(j,i,k) = xwwb%b0(j,i,k)
                atm2%w(j,i,k) = xwwb%b0(j,i,k)
              end do
            end do
          end do
        end if
        call exchange(sfs%psa,1,jce1,jce2,ice1,ice2)
        call psc2psd(sfs%psa,sfs%psdota)
        call exchange(sfs%psdota,1,jde1,jde2,ide1,ide2)
        call exchange(sfs%psb,idif,jce1,jce2,ice1,ice2)
        call psc2psd(sfs%psb,sfs%psdotb)
        call exchange(sfs%psdotb,idif,jde1,jde2,ide1,ide2)
      else
        do k = 1 , kz
          do i = ide1 , ide2
            do j = jce1 , jce2
              mo_atm%v(j,i,k) = xvb%b0(j,i,k)
            end do
          end do
        end do
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jde1 , jde2
              mo_atm%u(j,i,k) = xub%b0(j,i,k)
            end do
          end do
        end do
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              mo_atm%t(j,i,k) = xtb%b0(j,i,k)
              mo_atm%qx(j,i,k,iqv) = xqb%b0(j,i,k)
            end do
          end do
        end do
        if ( ipptls > 1 ) then
          if ( is_present_qc( ) ) then
            do k = 1 , kz
              do i = ice1 , ice2
                do j = jce1 , jce2
                  mo_atm%qx(j,i,k,iqc) = xlb%b0(j,i,k)
                end do
              end do
            end do
          else
            do k = 1 , kz
              do i = ice1 , ice2
                do j = jce1 , jce2
                  if ( mo_atm%t(j,i,k) > 253.15 ) then
                    mo_atm%qx(j,i,k,iqc) = minqc
                  end if
                end do
              end do
            end do
          end if
          if ( is_present_qi( ) ) then
            do k = 1 , kz
              do i = ice1 , ice2
                do j = jce1 , jce2
                  mo_atm%qx(j,i,k,iqi) = xib%b0(j,i,k)
                end do
              end do
            end do
          else
            do k = 1 , kz
              do i = ice1 , ice2
                do j = jce1 , jce2
                  if ( mo_atm%t(j,i,k) < 253.15 ) then
                    mo_atm%qx(j,i,k,iqi) = minqc
                  end if
                end do
              end do
            end do
          end if
        else if ( ipptls == 1 ) then
          if ( is_present_qc( ) ) then
            do k = 1 , kz
              do i = ice1 , ice2
                do j = jce1 , jce2
                  mo_atm%qx(j,i,k,iqc) = xlb%b0(j,i,k)
                end do
              end do
            end do
          else
            do k = 1 , kz
              do i = ice1 , ice2
                do j = jce1 , jce2
                  mo_atm%qx(j,i,k,iqc) = minqc
                end do
              end do
            end do
          end if
        end if
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              mo_atm%tvirt(j,i,k) = mo_atm%t(j,i,k) * &
                           (d_one + ep1*mo_atm%qx(j,i,k,iqv))
              mo_atm%pai(j,i,k) = xpaib%b0(j,i,k)
            end do
          end do
        end do
        ! Compute pressure and density
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              mo_atm%p(j,i,k) = (mo_atm%pai(j,i,k)**cpovr) * p00
              mo_atm%rho(j,i,k) = mo_atm%p(j,i,k)/(rgas* mo_atm%t(j,i,k))
            end do
          end do
        end do
        do i = ice1 , ice2
          do j = jce1 , jce2
            sfs%psa(j,i) = xpsb%b0(j,i)
          end do
        end do
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              mo_atm%qs(j,i,k) = pfwsat(mo_atm%t(j,i,k),mo_atm%p(j,i,k))
              ! Remove excessive supersaturation
              mo_atm%qx(j,i,k,iqv) = &
                  min(mo_atm%qx(j,i,k,iqv),mo_atm%qs(j,i,k))
            end do
          end do
        end do
        do i = ice1 , ice2
          do j = jce1 , jce2
            mo_atm%pf(j,i,kzp1) = sfs%psa(j,i)
          end do
        end do
        do k = kz , 2 , -1
          do i = ice1 , ice2
            do j = jce1 , jce2
              mo_atm%pf(j,i,k) = p00 * &
                (d_half*(mo_atm%pai(j,i,k)+mo_atm%pai(j,i,k-1)))**cpovr
            end do
          end do
        end do
        do i = ice1 , ice2
          do j = jce1 , jce2
            mo_atm%pf(j,i,1) = mo_atm%p(j,i,1) - egrav * mo_atm%rho(j,i,1) * &
              (mo_atm%zetaf(j,i,1)-mo_atm%zeta(j,i,1))
          end do
        end do

        call exchange_lr(mo_atm%u,1,jde1,jde2,ice1,ice2,1,kz)
        call exchange_bt(mo_atm%v,1,jce1,jce2,ide1,ide2,1,kz)
        mo_atm%w(:,:,1) = d_zero
        mo_atm%w(:,:,kzp1) = d_zero
        do k = 2 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              mo_atm%w(j,i,k) = mo_atm%w(j,i,k-1) + &
                 (mo_atm%zetaf(j,i,k-1)-mo_atm%zetaf(j,i,k)) * rdx * &
                 ((mo_atm%u(j+1,i,k)-mo_atm%u(j,i,k)) + &
                  (mo_atm%v(j,i+1,k)-mo_atm%v(j,i,k)))
            end do
          end do
        end do
      end if

      do i = ici1 , ici2
        do j = jci1 , jci2
          sfs%tg(j,i) = xtsb%b0(j,i)
        end do
      end do
      !
      ! Initialize PBL Hgt
      !
      zpbl(:,:) = 500.0_rkx
      !
      ! Inizialize the surface atmospheric temperature
      !
      do i = ici1 , ici2
        do j = jci1 , jci2
          sfs%tgbb(j,i) = sfs%tg(j,i)
        end do
      end do
      !
      ! Initialize surface parameters for aerosol scheme
      !
      if ( ichem == 1 ) then
        sfracv2d(:,:)  = d_half
        sfracb2d(:,:)  = d_half
      end if
      !
      ! Set the TKE variables for UW PBL to a default value
      !
      if ( idynamic == 3 ) then
        if ( ibltyp == 2 ) then
          mo_atm%tke(:,:,:) = tkemin
        else if ( ibltyp == 4 ) then
          atms%tkepbl(:,:,:) = tkemin
          sfs%uz0 = d_zero
          sfs%vz0 = d_zero
          sfs%thz0 = d_zero
          sfs%qz0 = d_zero
        end if
      else
        if ( ibltyp == 2 ) then
          atm1%tke(:,:,:) = tkemin
          atm2%tke(:,:,:) = tkemin
        else if ( ibltyp == 4 ) then
          atms%tkepbl = tkemin
          sfs%uz0 = d_zero
          sfs%vz0 = d_zero
          sfs%thz0 = d_zero
          sfs%qz0 = d_zero
        end if
      end if
      !
      ! Init the diurnal cycle SST scheme
      !
      if ( idcsst == 1 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            do n = 1 , nnsg
              lms%sst(n,j,i) = xtsb%b0(j,i)
            end do
          end do
        end do
        lms%tskin(:,:,:) = lms%sst
        lms%deltas(:,:,:) = 0.001_rkx
        lms%tdeltas(:,:,:) = lms%sst(:,:,:) - lms%deltas(:,:,:)
      end if
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nnsg
            lms%um10(n,j,i) = 1.0_rkx
          end do
        end do
      end do
      !
      ! Some output fields rely on atms structure.
      !
      call mkslice
      !
      ! End of initial run case
      !
    else
      if ( ichem == 1 .and. ichecold == 1 ) then
        if ( myid == italk ) then
          write(stdout,*) 'Starting tracer run in a previous run with ichem==0'
        end if
        if ( ichem == 1 ) then
          sfracv2d(:,:)  = d_half
          sfracb2d(:,:)  = d_half
        end if
      end if
      !
      ! When restarting, read in the data saved from previous run
      !
      call read_savefile(rcmtimer%idate)
      !
      ! Comunicate the data to other processors
      !
      if ( do_parallel_save ) then
        if ( idynamic == 3 ) then
          mo_atm%u(jde1:jde2,ice1:ice2,:) = atm_u_io(jde1:jde2,ice1:ice2,:)
          mo_atm%v(jce1:jce2,ide1:ide2,:) = atm_v_io(jce1:jce2,ide1:ide2,:)
          mo_atm%w(jce1:jce2,ice1:ice2,:) = atm_w_io(jce1:jce2,ice1:ice2,:)
          mo_atm%t(jce1:jce2,ice1:ice2,:) = atm_t_io(jce1:jce2,ice1:ice2,:)
          mo_atm%pai(jce1:jce2,ice1:ice2,:) = atm_pai_io(jce1:jce2,ice1:ice2,:)
          mo_atm%qx(jce1:jce2,ice1:ice2,:,:) = &
                             atm_qx_io(jce1:jce2,ice1:ice2,:,:)
          if ( ibltyp == 2 ) then
            mo_atm%tke(jce1:jce2,ice1:ice2,:) = &
                     atm_tke_io(jce1:jce2,ice1:ice2,:)
          end if
          if ( ichem == 1 .and. ichecold == 0 ) then
            mo_atm%trac(jce1:jce2,ice1:ice2,:,:) = &
                       trac_io(jce1:jce2,ice1:ice2,:,:)
          end if
          sfs%psa(jce1:jce2,ice1:ice2) = ps_io(jce1:jce2,ice1:ice2)
        else
          atm1%u(jde1:jde2,ide1:ide2,:) = atm1_u_io(jde1:jde2,ide1:ide2,:)
          atm1%v(jde1:jde2,ide1:ide2,:) = atm1_v_io(jde1:jde2,ide1:ide2,:)
          atm1%t(jce1:jce2,ice1:ice2,:) = atm1_t_io(jce1:jce2,ice1:ice2,:)
          atm1%qx(jce1:jce2,ice1:ice2,:,:) = &
                             atm1_qx_io(jce1:jce2,ice1:ice2,:,:)
          atm2%u(jde1:jde2,ide1:ide2,:) = atm2_u_io(jde1:jde2,ide1:ide2,:)
          atm2%v(jde1:jde2,ide1:ide2,:) = atm2_v_io(jde1:jde2,ide1:ide2,:)
          atm2%t(jce1:jce2,ice1:ice2,:) = atm2_t_io(jce1:jce2,ice1:ice2,:)
          atm2%qx(jce1:jce2,ice1:ice2,:,:) = &
                             atm2_qx_io(jce1:jce2,ice1:ice2,:,:)
          if ( ibltyp == 2 ) then
            atm1%tke(jce1:jce2,ice1:ice2,:) = &
                     atm1_tke_io(jce1:jce2,ice1:ice2,:)
            atm2%tke(jce1:jce2,ice1:ice2,:) = &
                     atm2_tke_io(jce1:jce2,ice1:ice2,:)
          end if
          if ( ichem == 1 .and. ichecold == 0 ) then
            atm1%chi(jce1:jce2,ice1:ice2,:,:) = &
                     chia_io(jce1:jce2,ice1:ice2,:,:)
            atm2%chi(jce1:jce2,ice1:ice2,:,:) = &
                     chib_io(jce1:jce2,ice1:ice2,:,:)
          end if
          sfs%psa(jce1:jce2,ice1:ice2) = psa_io(jce1:jce2,ice1:ice2)
          sfs%psb(jce1:jce2,ice1:ice2) = psb_io(jce1:jce2,ice1:ice2)
        end if
        if ( ibltyp == 2 ) then
          kpbl = kpbl_io
        end if
        if ( ibltyp == 4 ) then
          atms%tkepbl = tke_pbl_io
          kpbl = kpbl_io
          sfs%uz0 = myjsf_uz0_io
          sfs%vz0 = myjsf_vz0_io
          sfs%thz0 = myjsf_thz0_io
          sfs%qz0 = myjsf_qz0_io
        end if
        if ( idynamic == 2 ) then
          atm1%pp(jce1:jce2,ice1:ice2,:) = atm1_pp_io(jce1:jce2,ice1:ice2,:)
          atm2%pp(jce1:jce2,ice1:ice2,:) = atm2_pp_io(jce1:jce2,ice1:ice2,:)
          atm1%w(jce1:jce2,ice1:ice2,:) = atm1_w_io(jce1:jce2,ice1:ice2,:)
          atm2%w(jce1:jce2,ice1:ice2,:) = atm2_w_io(jce1:jce2,ice1:ice2,:)
        end if
        sfs%hfx = hfx_io
        sfs%qfx = qfx_io
        sfs%tgbb = tgbb_io
        sfs%zo = zo_io
        sfs%uvdrag = uvdrag_io
        sfs%ram1 = ram_io
        sfs%rah1 = rah_io
        sfs%br = br_io
        sfs%q2m = q2m_io
        sfs%u10m = u10m_io
        sfs%v10m = v10m_io
        sfs%w10m = w10m_io
        sfs%ustar = ustar_io
        if ( ipptls > 0 ) then
          fcc = fcc_io
        end if
        heatrt = heatrt_io
        o3prof = o3prof_io
        if ( iocnflx == 2 ) then
          zpbl = zpbl_io
        end if
        if ( any(icup == 3) ) then
          cldefi = cldefi_io
        end if
        if ( any(icup == 4) ) then
          cbmf2d = cbmf2d_io
        end if
        if ( any(icup == 6) .or. any(icup == 5) ) then
          avg_ww = cu_avg_ww_io
        end if
        if ( irrtm == 0 ) then
          gasabsnxt = gasabsnxt_io
          gasabstot = gasabstot_io
          gasemstot = gasemstot_io
        end if
        lms%sw = sw_io
#ifdef CLM45
        if ( ichem == 1 .and. ichecold == 0 ) then
          tsoi = tsoi_io
          sw_vol = swvol_io
        end if
#else
        lms%gwet = gwet_io
        lms%ldew = ldew_io
        lms%taf = taf_io
#endif
        lms%tgrd = tgrd_io
        lms%tgbrd = tgbrd_io
        lms%tlef = tlef_io
        lms%sncv = sncv_io
        lms%sfice = sfice_io
        lms%snag = snag_io
        lms%emisv = emisv_io
        lms%um10 = um10_io
        lms%swalb = swalb_io
        lms%lwalb = lwalb_io
        lms%swdiralb = swdiralb_io
        lms%swdifalb = swdifalb_io
        lms%lwdiralb = lwdiralb_io
        lms%lwdifalb = lwdifalb_io
        mdsub%ldmsk(:,jci1:jci2,ici1:ici2) = ldmsk1_io(:,jci1:jci2,ici1:ici2)
        solis = solis_io
        solvs = solvs_io
        solvsd = solvsd_io
        solvl = solvl_io
        solvld = solvld_io
        sabveg = sabveg_io
        flw = flw_io
        flwd = flwd_io
        fsw = fsw_io
        sinc = sinc_io
        mddom%ldmsk = ldmsk_io
#ifndef CLM
        if ( lakemod == 1 ) then
          lms%eta = eta_io
          lms%hi = hi_io
          lms%tlake = tlak_io
        end if
#else
        if ( imask == 2 ) then
          mddom%lndcat(jce1:jce2,ice1:ice2) = lndcat_io(jce1:jce2,ice1:ice2)
        end if
#endif
        if ( idcsst == 1 ) then
          lms%sst = sst_io
          lms%tskin = tskin_io
          lms%deltas = deltas_io
          lms%tdeltas = tdeltas_io
        end if
        if ( idynamic == 1 ) then
          dstor(jde1:jde2,ide1:ide2,:) = dstor_io(jde1:jde2,ide1:ide2,:)
          hstor(jde1:jde2,ide1:ide2,:) = hstor_io(jde1:jde2,ide1:ide2,:)
        end if
        if ( ichem == 1 .and. ichecold == 0 ) then
          convpr = convpr_io
          rainout = rainout_io
          washout = washout_io
          remdrd = remdrd_io
          if ( igaschem == 1 .and. ichsolver > 0 ) then
            chemall = chemall_io
            taucldsp = taucldsp_io
          end if
          ssw2da = ssw2da_io
#ifdef CLM45
          dustflx_clm = duflux_io
          voc_em_clm = voflux_io
#else
          sdelt = sdelt_io
          sdelq = sdelq_io
          svegfrac2d = svegfrac2d_io
#endif
          sfracv2d = sfracv2d_io
          sfracb2d = sfracb2d_io
          sfracs2d = sfracs2d_io
        end if
        if ( islab_ocean == 1 .and. do_restore_sst ) then
          qflux_restore_sst = qflux_restore_sst_io
        end if
      else
        if ( idynamic == 3 ) then
          call grid_distribute(atm_u_io,mo_atm%u,jde1,jde2,ice1,ice2,1,kz)
          call grid_distribute(atm_v_io,mo_atm%v,jce1,jce2,ide1,ide2,1,kz)
          call grid_distribute(atm_w_io,mo_atm%w,jce1,jce2,ice1,ice2,1,kzp1)
          call grid_distribute(atm_t_io,mo_atm%t,jce1,jce2,ice1,ice2,1,kz)
          call grid_distribute(atm_pai_io,mo_atm%pai,jce1,jce2,ice1,ice2,1,kz)
          call grid_distribute(atm_qx_io,mo_atm%qx, &
                               jce1,jce2,ice1,ice2,1,kz,1,nqx)
          if ( ibltyp == 2 ) then
            call grid_distribute(atm_tke_io,mo_atm%tke, &
                                 jce1,jce2,ice1,ice2,1,kzp1)
          end if
          if ( ichem == 1 ) then
            call grid_distribute(trac_io,mo_atm%trac, &
                                 jce1,jce2,ice1,ice2,1,kz,1,ntr)
          end if
          call grid_distribute(ps_io,sfs%psa,jce1,jce2,ice1,ice2)
        else
          call grid_distribute(atm1_u_io,atm1%u,jde1,jde2,ide1,ide2,1,kz)
          call grid_distribute(atm1_v_io,atm1%v,jde1,jde2,ide1,ide2,1,kz)
          call grid_distribute(atm1_t_io,atm1%t,jce1,jce2,ice1,ice2,1,kz)
          call grid_distribute(atm1_qx_io,atm1%qx, &
                               jce1,jce2,ice1,ice2,1,kz,1,nqx)
          call grid_distribute(atm2_u_io,atm2%u,jde1,jde2,ide1,ide2,1,kz)
          call grid_distribute(atm2_v_io,atm2%v,jde1,jde2,ide1,ide2,1,kz)
          call grid_distribute(atm2_t_io,atm2%t,jce1,jce2,ice1,ice2,1,kz)
          call grid_distribute(atm2_qx_io,atm2%qx, &
                               jce1,jce2,ice1,ice2,1,kz,1,nqx)
          if ( ibltyp == 2 ) then
            call grid_distribute(atm1_tke_io,atm1%tke, &
                                 jce1,jce2,ice1,ice2,1,kzp1)
            call grid_distribute(atm2_tke_io,atm2%tke, &
                                 jce1,jce2,ice1,ice2,1,kzp1)
          end if
          if ( ichem == 1 ) then
            call grid_distribute(chia_io,atm1%chi, &
                                 jce1,jce2,ice1,ice2,1,kz,1,ntr)
            call grid_distribute(chib_io,atm2%chi, &
                                 jce1,jce2,ice1,ice2,1,kz,1,ntr)
          end if
          call grid_distribute(psa_io,sfs%psa,jce1,jce2,ice1,ice2)
          call grid_distribute(psb_io,sfs%psb,jce1,jce2,ice1,ice2)
        end if
        if ( ibltyp == 2 ) then
          call grid_distribute(kpbl_io,kpbl,jci1,jci2,ici1,ici2)
        end if
        if ( ibltyp == 4 ) then
          call grid_distribute(tke_pbl_io,atms%tkepbl, &
                               jci1,jci2,ici1,ici2,1,kz)
          call grid_distribute(kpbl_io,kpbl,jci1,jci2,ici1,ici2)
          call grid_distribute(myjsf_uz0_io,sfs%uz0,jci1,jci2,ici1,ici2)
          call grid_distribute(myjsf_vz0_io,sfs%vz0,jci1,jci2,ici1,ici2)
          call grid_distribute(myjsf_thz0_io,sfs%thz0,jci1,jci2,ici1,ici2)
          call grid_distribute(myjsf_qz0_io,sfs%qz0,jci1,jci2,ici1,ici2)
        end if
        if ( idynamic == 2 ) then
          call grid_distribute(atm1_w_io,atm1%w,jce1,jce2,ice1,ice2,1,kzp1)
          call grid_distribute(atm2_w_io,atm2%w,jce1,jce2,ice1,ice2,1,kzp1)
          call grid_distribute(atm1_pp_io,atm1%pp,jce1,jce2,ice1,ice2,1,kz)
          call grid_distribute(atm2_pp_io,atm2%pp,jce1,jce2,ice1,ice2,1,kz)
        end if
        call grid_distribute(hfx_io,sfs%hfx,jci1,jci2,ici1,ici2)
        call grid_distribute(qfx_io,sfs%qfx,jci1,jci2,ici1,ici2)
        call grid_distribute(tgbb_io,sfs%tgbb,jci1,jci2,ici1,ici2)
        call grid_distribute(zo_io,sfs%zo,jci1,jci2,ici1,ici2)
        call grid_distribute(uvdrag_io,sfs%uvdrag,jci1,jci2,ici1,ici2)
        call grid_distribute(ram_io,sfs%ram1,jci1,jci2,ici1,ici2)
        call grid_distribute(rah_io,sfs%rah1,jci1,jci2,ici1,ici2)
        call grid_distribute(br_io,sfs%br,jci1,jci2,ici1,ici2)
        call grid_distribute(q2m_io,sfs%q2m,jci1,jci2,ici1,ici2)
        call grid_distribute(u10m_io,sfs%u10m,jci1,jci2,ici1,ici2)
        call grid_distribute(v10m_io,sfs%v10m,jci1,jci2,ici1,ici2)
        call grid_distribute(w10m_io,sfs%w10m,jci1,jci2,ici1,ici2)
        call grid_distribute(ustar_io,sfs%ustar,jci1,jci2,ici1,ici2)
        if ( ipptls > 0 ) then
          call grid_distribute(fcc_io,fcc,jci1,jci2,ici1,ici2,1,kz)
        end if
        call grid_distribute(heatrt_io,heatrt,jci1,jci2,ici1,ici2,1,kz)
        call grid_distribute(o3prof_io,o3prof,jci1,jci2,ici1,ici2,1,kzp1)
        if ( iocnflx == 2 .or. ibltyp == 3 ) then
          call grid_distribute(zpbl_io,zpbl,jci1,jci2,ici1,ici2)
        end if
        if ( any(icup == 3) ) then
          call grid_distribute(cldefi_io,cldefi,jci1,jci2,ici1,ici2)
        end if
        if ( any(icup == 4) ) then
          call grid_distribute(cbmf2d_io,cbmf2d,jci1,jci2,ici1,ici2)
        end if
        if ( any(icup == 6) .or. any(icup == 5 ) ) then
          call grid_distribute(cu_avg_ww_io,avg_ww,jci1,jci2,ici1,ici2,1,kz)
        end if
        if ( irrtm == 0 ) then
          call grid_distribute(gasabsnxt_io,gasabsnxt, &
                               jci1,jci2,ici1,ici2,1,kz,1,4)
          call grid_distribute(gasabstot_io,gasabstot, &
                               jci1,jci2,ici1,ici2,1,kzp1,1,kzp1)
          call grid_distribute(gasemstot_io,gasemstot, &
                               jci1,jci2,ici1,ici2,1,kzp1)
        end if
        call subgrid_distribute(sw_io,lms%sw,jci1,jci2, &
                                ici1,ici2,1,num_soil_layers)
#ifdef CLM45
        if ( ichem == 1 ) then
          call grid_distribute(tsoi_io,tsoi,jci1,jci2, &
                               ici1,ici2,1,num_soil_layers)
          call grid_distribute(swvol_io,sw_vol,jci1,jci2, &
                               ici1,ici2,1,num_soil_layers)
        end if
#else
        call subgrid_distribute(gwet_io,lms%gwet,jci1,jci2,ici1,ici2)
        call subgrid_distribute(ldew_io,lms%ldew,jci1,jci2,ici1,ici2)
        call subgrid_distribute(taf_io,lms%taf,jci1,jci2,ici1,ici2)
#endif
        call subgrid_distribute(tgrd_io,lms%tgrd,jci1,jci2,ici1,ici2)
        call subgrid_distribute(tgbrd_io,lms%tgbrd,jci1,jci2,ici1,ici2)
        call subgrid_distribute(tlef_io,lms%tlef,jci1,jci2,ici1,ici2)
        call subgrid_distribute(sncv_io,lms%sncv,jci1,jci2,ici1,ici2)
        call subgrid_distribute(sfice_io,lms%sfice,jci1,jci2,ici1,ici2)
        call subgrid_distribute(snag_io,lms%snag,jci1,jci2,ici1,ici2)
        call subgrid_distribute(emisv_io,lms%emisv,jci1,jci2,ici1,ici2)
        call subgrid_distribute(um10_io,lms%um10,jci1,jci2,ici1,ici2)
        call subgrid_distribute(swalb_io,lms%swalb,jci1,jci2,ici1,ici2)
        call subgrid_distribute(lwalb_io,lms%lwalb,jci1,jci2,ici1,ici2)
        call subgrid_distribute(swdiralb_io,lms%swdiralb,jci1,jci2,ici1,ici2)
        call subgrid_distribute(swdifalb_io,lms%swdifalb,jci1,jci2,ici1,ici2)
        call subgrid_distribute(lwdiralb_io,lms%lwdiralb,jci1,jci2,ici1,ici2)
        call subgrid_distribute(lwdifalb_io,lms%lwdifalb,jci1,jci2,ici1,ici2)
        call subgrid_distribute(ldmsk1_io,mdsub%ldmsk,jci1,jci2,ici1,ici2)
        call grid_distribute(solis_io,solis,jci1,jci2,ici1,ici2)
        call grid_distribute(solvs_io,solvs,jci1,jci2,ici1,ici2)
        call grid_distribute(solvsd_io,solvsd,jci1,jci2,ici1,ici2)
        call grid_distribute(solvl_io,solvl,jci1,jci2,ici1,ici2)
        call grid_distribute(solvld_io,solvld,jci1,jci2,ici1,ici2)
        call grid_distribute(sabveg_io,sabveg,jci1,jci2,ici1,ici2)
        call grid_distribute(flw_io,flw,jci1,jci2,ici1,ici2)
        call grid_distribute(fsw_io,fsw,jci1,jci2,ici1,ici2)
        call grid_distribute(flwd_io,flwd,jci1,jci2,ici1,ici2)
        call grid_distribute(sinc_io,sinc,jci1,jci2,ici1,ici2)
        call grid_distribute(ldmsk_io,mddom%ldmsk,jci1,jci2,ici1,ici2)
#ifndef CLM
        if ( lakemod == 1 ) then
          call subgrid_distribute(eta_io,lms%eta,jci1,jci2,ici1,ici2)
          call subgrid_distribute(hi_io,lms%hi,jci1,jci2,ici1,ici2)
          call subgrid_distribute(tlak_io,lms%tlake, &
                                  jci1,jci2,ici1,ici2,1,ndpmax)
        endif
#else
        if ( imask == 2 ) then
          call grid_distribute(lndcat_io,mddom%lndcat,jci1,jci2,ici1,ici2)
        end if
#endif
        if ( idcsst == 1 ) then
          call subgrid_distribute(sst_io,lms%sst,jci1,jci2,ici1,ici2)
          call subgrid_distribute(tskin_io,lms%tskin,jci1,jci2,ici1,ici2)
          call subgrid_distribute(deltas_io,lms%deltas,jci1,jci2,ici1,ici2)
          call subgrid_distribute(tdeltas_io,lms%tdeltas,jci1,jci2,ici1,ici2)
        end if
        if ( idynamic == 1 ) then
          call grid_distribute(dstor_io,dstor,jde1,jde2,ide1,ide2,1,nsplit)
          call grid_distribute(hstor_io,hstor,jde1,jde2,ide1,ide2,1,nsplit)
        end if
        if ( ichem == 1 ) then
          call grid_distribute(convpr_io,convpr,jci1,jci2,ici1,ici2,1,kz)
          call grid_distribute(rainout_io,rainout, &
                               jci1,jci2,ici1,ici2,1,kz,1,ntr)
          call grid_distribute(washout_io,washout, &
                               jci1,jci2,ici1,ici2,1,kz,1,ntr)
          call grid_distribute(remdrd_io,remdrd,jci1,jci2,ici1,ici2,1,ntr)
          if ( igaschem == 1 .and. ichsolver > 0 ) then
            call grid_distribute(chemall_io,chemall, &
                                  jci1,jci2,ici1,ici2,1,kz,1,totsp)
            call grid_distribute(taucldsp_io,taucldsp, &
                                 jci1,jci2,ici1,ici2,0,kz,1,nspi)
          end if
          call grid_distribute(ssw2da_io,ssw2da,jci1,jci2,ici1,ici2)
#ifdef CLM45
          call grid_distribute(duflux_io,dustflx_clm,jci1,jci2,ici1,ici2,1,4)
          call grid_distribute(voflux_io,voc_em_clm,jci1,jci2,ici1,ici2,1,ntr)
#else
          call grid_distribute(sdelt_io,sdelt,jci1,jci2,ici1,ici2)
          call grid_distribute(sdelq_io,sdelq,jci1,jci2,ici1,ici2)
          call grid_distribute(svegfrac2d_io,svegfrac2d,jci1,jci2,ici1,ici2)
#endif
          call grid_distribute(sfracv2d_io,sfracv2d,jci1,jci2,ici1,ici2)
          call grid_distribute(sfracb2d_io,sfracb2d,jci1,jci2,ici1,ici2)
          call grid_distribute(sfracs2d_io,sfracs2d,jci1,jci2,ici1,ici2)
        end if
      end if

      if ( myid == italk ) then
        ozprnt = o3prof(3,3,:)
        call vprntv(ozprnt,kzp1,'Ozone profiles restart')
      end if

      if ( idynamic /= 3 ) then
        call exchange(sfs%psb,idif,jce1,jce2,ice1,ice2)
        call exchange(sfs%psa,1,jce1,jce2,ice1,ice2)
        call psc2psd(sfs%psa,sfs%psdota)
        call psc2psd(sfs%psb,sfs%psdotb)
        call exchange(sfs%psdota,1,jde1,jde2,ide1,ide2)
        call exchange(sfs%psdotb,idif,jde1,jde2,ide1,ide2)
      end if

      if ( idynamic == 2 ) then
        do i = ice1 , ice2
          do j = jce1 , jce2
            sfs%psc(j,i) = sfs%psa(j,i)
          end do
        end do
      end if
#ifdef CLM
      !
      ! CLM modifies landuse table. Get the modified one from restart file
      !
      if ( imask == 2 ) then
        do n = 1 , nnsg
          mdsub%lndcat(n,jci1:jci2,ici1:ici2) = &
                              mddom%lndcat(jci1:jci2,ici1:ici2)
        end do
      end if
#endif
      rdnnsg = d_one/real(nnsg,rkx)
      aldirs = sum(lms%swdiralb,1)*rdnnsg
      aldirl = sum(lms%lwdiralb,1)*rdnnsg
      aldifs = sum(lms%swdifalb,1)*rdnnsg
      aldifl = sum(lms%lwdifalb,1)*rdnnsg
      albvs = sum(lms%swalb,1)*rdnnsg
      albvl = sum(lms%lwalb,1)*rdnnsg
      sfs%tg = sum(lms%tgrd,1)*rdnnsg

      call bcast(declin)
      call bcast(solcon)

      if ( islab_ocean == 1 .and. do_restore_sst ) then
        call grid_distribute(qflux_restore_sst_io,qflux_restore_sst, &
          jci1,jci2,ici1,ici2,1,12)
        call bcast(stepcount)
      end if
      if ( idynamic == 2 .and. ifupr == 1 ) then
        call bcast(tmask)
      end if

      if ( debug_level > 0 ) then
        call bcast(dryini)
        call bcast(watini)
        call bcast(dryerror)
        call bcast(waterror)
      end if
      !
      ! Init boundary
      !
      call init_bdy
      !
      ! Report success
      !
      if ( myid == italk ) then
        write(stdout,*) 'Successfully read restart file at time = ', &
                rcmtimer%str( )
      end if
      !
      ! Setup all timeseps for a restart
      !
      dtbat = dtsrf
      if ( idynamic /= 3 ) then
        dt = dt2
        rdt = d_one/dt
        dtsq = dt*dt
        dtcb = dt*dt*dt
      end if
      !
      ! End of restart phase
      !
    end if

    if ( idynamic == 3 ) then
      call init_moloch
    end if

    if ( idynamic == 1 ) then
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            atm1%pr(j,i,k) = (hsigma(k)*sfs%psa(j,i) + ptop)*d_1000
            atm2%pr(j,i,k) = (hsigma(k)*sfs%psb(j,i) + ptop)*d_1000
          end do
        end do
      end do
    else if ( idynamic == 2 ) then
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            atm1%pr(j,i,k) = atm0%pr(j,i,k) + atm1%pp(j,i,k)/sfs%psa(j,i)
            atm1%rho(j,i,k) = atm1%pr(j,i,k) /        &
                (rgas*atm1%t(j,i,k)/sfs%psa(j,i) *    &
                (d_one+ep1*atm1%qx(j,i,k,iqv)/sfs%psa(j,i)))
            atm2%pr(j,i,k) = atm0%pr(j,i,k) + atm2%pp(j,i,k)/sfs%psb(j,i)
          end do
        end do
      end do
    else if ( idynamic == 3 ) then
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            mo_atm%p(j,i,k) = (mo_atm%pai(j,i,k)**cpovr) * p00
            mo_atm%rho(j,i,k) = mo_atm%p(j,i,k)/(rgas*mo_atm%t(j,i,k))
            mo_atm%qs(j,i,k) = pfwsat(mo_atm%t(j,i,k),mo_atm%p(j,i,k))
          end do
        end do
      end do
      if ( ipptls > 0 ) then
        if ( ipptls > 1 ) then
          do k = 1 , kz
            do i = ice1 , ice2
              do j = jce1 , jce2
                mo_atm%tvirt(j,i,k) = mo_atm%t(j,i,k) * &
                       (d_one + ep1*mo_atm%qx(j,i,k,iqv) - &
                        mo_atm%qx(j,i,k,iqc) - mo_atm%qx(j,i,k,iqi) - &
                        mo_atm%qx(j,i,k,iqr) - mo_atm%qx(j,i,k,iqs))
              end do
            end do
          end do
        else
          do k = 1 , kz
            do i = ice1 , ice2
              do j = jce1 , jce2
                mo_atm%tvirt(j,i,k) = mo_atm%t(j,i,k) * &
                               (d_one + ep1*mo_atm%qx(j,i,k,iqv) - &
                                mo_atm%qx(j,i,k,iqc))
              end do
            end do
          end do
        end if
      else
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              mo_atm%tvirt(j,i,k) = mo_atm%t(j,i,k) * &
                             (d_one + ep1*mo_atm%qx(j,i,k,iqv))
            end do
          end do
        end do
      end if
      do i = ice1 , ice2
        do j = jce1 , jce2
          mo_atm%pf(j,i,kzp1) = sfs%psa(j,i)
        end do
      end do
      do k = kz , 2 , -1
        do i = ice1 , ice2
          do j = jce1 , jce2
            mo_atm%pf(j,i,k) = p00 * &
              (d_half*(mo_atm%pai(j,i,k)+mo_atm%pai(j,i,k-1)))**cpovr
          end do
        end do
      end do
      ! Top pressure
      do i = ice1 , ice2
        do j = jce1 , jce2
          mo_atm%pf(j,i,1) = mo_atm%p(j,i,1) - egrav * mo_atm%rho(j,i,1) * &
            (mo_atm%zetaf(j,i,1)-mo_atm%zeta(j,i,1))
        end do
      end do
    end if
    !
    ! pressure of tropopause
    !
    if ( irceideal == 1 ) then
      ptrop = 10000_rkx
    else
      do i = ici1 , ici2
        do j = jci1 , jci2
          ptrop(j,i) = 17.0e3_rkx - 6.0e3_rkx*cos(mddom%xlat(j,i)*degrad)**2
        end do
      end do
    end if

    if ( .not. ifrest ) then
      if ( any(icup == 6)  .or. any(icup == 5) ) then
        if ( idynamic == 2 ) then
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( cuscheme(j,i) == 6 ) then
                  avg_ww(j,i,k) = 0.5_rkx * &
                    (atm1%w(j,i,k) + atm1%w(j,i,k+1)) / sfs%psa(j,i)
                end if
                if ( cuscheme(j,i) == 5 ) then
                  avg_ww(j,i,k) = -0.5_rkx * egrav * atm1%rho(j,i,k) * &
                    (atm1%w(j,i,k) + atm1%w(j,i,k+1)) / sfs%psa(j,i)
                end if
              end do
            end do
          end do
        end if
      end if
    end if
    !
    ! The following allows to change landuse on restart.
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        mddom%iveg(j,i) = nint(mddom%lndcat(j,i))
        mddom%itex(j,i) = nint(mddom%lndtex(j,i))
        do n = 1 , nnsg
          mdsub%iveg(n,j,i) = nint(mdsub%lndcat(n,j,i))
          mdsub%itex(n,j,i) = nint(mdsub%lndtex(n,j,i))
        end do
      end do
    end do
    !
    ! Initialize solar elevation (zenith angle)
    !
    call zenitm(mddom%xlat,mddom%xlon,coszrs)
    !
    ! Initialize the Surface Model
    !
    if ( idynamic == 3 ) then
      if ( mo_nzfilt > 0 ) then
        ! Sponge layer at the top of the atmosphere
        zfilt = (kzp1-mo_nzfilt)*mo_dzita
        do k = 1 , kz
          if ( k > mo_nzfilt ) then
            ffilt(k) = d_zero
          else
            zzi = (mo_dzita*(kzp1-k)-zfilt)/(mo_ztop-zfilt)
            ffilt(k) = 0.8_rkx*sin(d_half*mathpi*zzi)**2
          end if
        end do
      else
        ffilt(:) = d_zero
      end if
    end if
    call initialize_surface_model
    if ( idynamic /= 3 ) then
      call initialize_diffusion
      if ( idynamic == 2 ) then
        call init_sound
      end if
    end if
    !
    ! RRTM_SW gas / abs constant initialisation
    !
    if ( irrtm == 1 ) then
      call rrtmg_sw_ini(cpd)
      call rrtmg_lw_ini(cpd)
    end if
    !
    ! chemistry initialisation
    !
    if ( ichem == 1 ) then
      call start_chem
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

  end subroutine init

end module mod_init

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
