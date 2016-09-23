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

module mod_output

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_runparams
  use mod_header
  use mod_mpmessage
  use mod_mppparam
  use mod_service
  use mod_atm_interface
  use mod_che_interface
  use mod_che_output
  use mod_lm_interface
  use mod_rad_interface
  use mod_cu_interface
  use mod_pbl_interface
  use mod_ncout
  use mod_bdycod
  use mod_precip
  use mod_split
  use mod_savefile
  use mod_slabocean
  use mod_cloud_s1

  implicit none

  private

  public :: output

  contains

  subroutine output
    implicit none
    logical :: ldoatm , ldosrf , ldorad , ldoche
    logical :: ldosav , ldolak , ldosub
    logical :: ldoslab
    logical :: lstartup
    integer(ik4) :: i , j , k , itr
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'output'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    lstartup = .false.
    if ( ktau == 0 .or. doing_restart ) then
      !
      ! Set up static variables (first time in)
      !
      if ( associated(xlon_out) ) then
        xlon_out = mddom%xlon(jci1:jci2,ici1:ici2)
        xlat_out = mddom%xlat(jci1:jci2,ici1:ici2)
        mask_out = mddom%mask(jci1:jci2,ici1:ici2)
        topo_out = mddom%ht(jci1:jci2,ici1:ici2)
        topo_out = topo_out*regrav
      end if
      if ( associated(sub_xlon_out) ) then
        call reorder_subgrid(mdsub%xlon,sub_xlon_out)
        call reorder_subgrid(mdsub%xlat,sub_xlat_out)
        call reorder_subgrid(mdsub%mask,sub_mask_out)
        call reorder_subgrid(mdsub%ht,sub_topo_out)
        sub_topo_out = sub_topo_out*regrav
      end if
      if ( idynamic == 2 ) then
        if ( associated(p0_out) ) then
          p0_out = atm0%ps(jci1:jci2,ici1:ici2) + ptop*d_1000
        end if
      end if
      !
      ! Reset the accumulation arrays
      !
      if ( associated(sts_tgmax_out) )  sts_tgmax_out  = -1.e30_rkx
      if ( associated(sts_tgmin_out) )  sts_tgmin_out  =  1.e30_rkx
      if ( associated(sts_t2max_out) )  sts_t2max_out  = -1.e30_rkx
      if ( associated(sts_t2min_out) )  sts_t2min_out  =  1.e30_rkx
      if ( associated(sts_w10max_out) ) sts_w10max_out = -1.e30_rkx
      if ( associated(sts_psmin_out) )  sts_psmin_out  =  1.e30_rkx
      if ( associated(sts_pcpmax_out) ) sts_pcpmax_out = -1.e30_rkx
      call newoutfiles(idatex)
      lstartup = .true.
      if ( doing_restart ) then
        doing_restart = .false.
#ifdef DEBUG
        call time_end(subroutine_name,idindx)
#endif
        return
      end if
    end if

    ldoatm = .false.
    ldosrf = .false.
    ldolak = .false.
    ldosub = .false.
    ldorad = .false.
    ldoche = .false.
    ldosav = .false.
    ldoslab = .false.

    if ( ktau > 0 ) then
      if ( ksav > 0 ) then
        if ( ktau == mtau .or. mod(ktau,ksav) == 0 ) then
          ldosav = .true.
        end if
      else
        if ( ksav == 0 ) then
          if ( ( ktau == mtau ) .or. &
               (lfdomonth(idatex) .and. lmidnight(idatex)) ) then
            ldosav = .true.
          end if
        else if ( ksav < 0 ) then
          if ( ktau == mtau .or. mod(ktau,-ksav) == 0 .or. &
               (lfdomonth(idatex) .and. lmidnight(idatex)) ) then
            ldosav = .true.
          end if
        end if
      end if
      if ( mod(ktau,katm) == 0 ) then
        ldoatm = .true.
      end if
      if ( mod(ktau,ksrf) == 0 ) then
        ldosrf = .true.
      end if
      if ( mod(ktau,klak) == 0 ) then
        ldolak= .true.
      end if
      if ( mod(ktau,ksub) == 0 ) then
        ldosub= .true.
      end if
      if ( mod(ktau,krad) == 0 ) then
        ldorad = .true.
      end if
      if ( mod(ktau,kche) == 0 ) then
        ldoche = .true.
      end if
      if ( idatex == idate2 ) then
        ldoslab = .true.
      end if
    end if

    if ( ktau == 0 ) then
      if ( debug_level > 2 ) then
        ldoatm = .true.
      end if
    end if

    if ( atm_stream > 0 ) then
      if ( ldoatm ) then
        ps_out = sfs%psa(jci1:jci2,ici1:ici2)
        if ( associated(atm_t_out) ) then
          do k = 1 , kz
            atm_t_out(:,:,k) = atm1%t(jci1:jci2,ici1:ici2,k)/ps_out
          end do
        end if
        if ( associated(atm_u_out) .and. associated(atm_v_out) ) then
          call exchange(atm1%u,1,jde1,jde2,ide1,ide2,1,kz)
          call exchange(atm1%v,1,jde1,jde2,ide1,ide2,1,kz)
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                atm_u_out(j,i,k) = d_rfour*(atm1%u(j,i,k)+atm1%u(j+1,i,k) + &
                                 atm1%u(j,i+1,k)+atm1%u(j+1,i+1,k))/ps_out(j,i)
                atm_v_out(j,i,k) = d_rfour*(atm1%v(j,i,k)+atm1%v(j+1,i,k) + &
                                 atm1%v(j,i+1,k)+atm1%v(j+1,i+1,k))/ps_out(j,i)
              end do
            end do
          end do
        end if
        if ( associated(atm_omega_out) ) &
          atm_omega_out = omega(jci1:jci2,ici1:ici2,:)*d_10
        if ( associated(atm_w_out) ) then
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                atm_w_out(j,i,k) = d_half * (atm1%w(j,i,k+1) + &
                                             atm1%w(j,i,k)) / ps_out(j,i)
              end do
            end do
          end do
        end if
        if ( associated(atm_pp_out) ) then
          do k = 1 , kz
            atm_pp_out(:,:,k) = atm1%pp(jci1:jci2,ici1:ici2,k)/ps_out
          end do
        end if
        if ( associated(atm_qv_out) ) then
          do k = 1 , kz
            atm_qv_out(:,:,k) = atm1%qx(jci1:jci2,ici1:ici2,k,iqv)/ps_out
          end do
          ! Specific humidity in the output, not mixing ratio
          atm_qv_out = atm_qv_out/(d_one+atm_qv_out)
        end if
        if ( associated(atm_qc_out) ) then
          do k = 1 , kz
            atm_qc_out(:,:,k) = atm1%qx(jci1:jci2,ici1:ici2,k,iqc)/ps_out
          end do
        end if
        if ( associated(atm_qr_out) ) then
          do k = 1 , kz
            atm_qr_out(:,:,k) = atm1%qx(jci1:jci2,ici1:ici2,k,iqr)/ps_out
          end do
        end if
        if ( associated(atm_qi_out) ) then
          do k = 1 , kz
            atm_qi_out(:,:,k) = atm1%qx(jci1:jci2,ici1:ici2,k,iqi)/ps_out
          end do
        end if
        if ( associated(atm_qs_out) ) then
           do k = 1 , kz
            atm_qs_out(:,:,k) = atm1%qx(jci1:jci2,ici1:ici2,k,iqs)/ps_out
          end do
        end if
        if ( associated(atm_rh_out) ) then
          do k = 1 , kz
            atm_rh_out(:,:,k) = atms%rhb3d(jci1:jci2,ici1:ici2,k)*d_100
          end do
        end if
        if ( associated(atm_zf_out) ) then
          do k = 1 , kz
            atm_zf_out(:,:,k) = atms%zq(jci1:jci2,ici1:ici2,k)
          end do
        end if
        if ( associated(atm_zh_out) ) then
          do k = 1 , kz
            atm_zh_out(:,:,k) = atms%za(jci1:jci2,ici1:ici2,k)
          end do
        end if
        if ( associated(atm_pf_out) ) then
          do k = 1 , kz
            atm_pf_out(:,:,k) = atms%pb3d(jci1:jci2,ici1:ici2,k)
          end do
        end if
        if ( associated(atm_ph_out) ) then
          do k = 1 , kz
            atm_ph_out(:,:,k) = atms%pf3d(jci1:jci2,ici1:ici2,k)
          end do
        end if
        if ( ipptls == 2 ) then
          if ( associated(atm_rainls_out) ) then
            do k = 1 , kz
              atm_rainls_out(:,:,k) = rain_ls(jci1:jci2,ici1:ici2,k)
            end do
          end if
        end if
        if ( associated(atm_raincc_out) ) then
          do k = 1 , kz
            atm_raincc_out(:,:,k) = rain_cc(jci1:jci2,ici1:ici2,k)
          end do
        end if
        if ( associated(atm_q_detr_out) ) then
          do k = 1 , kz
            atm_q_detr_out(:,:,k) = q_detr(jci1:jci2,ici1:ici2,k)
          end do
        end if

#ifdef DEBUG
        if ( ipptls == 2 .and. stats ) then
          if ( associated(atm_stats_supw_out) ) then
            do k = 1 , kz
              atm_stats_supw_out(:,:,k) = statssupw(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_supc_out) ) then
            do k = 1 , kz
              atm_stats_supc_out(:,:,k) = statssupc(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_detw_out) ) then
            do k = 1 , kz
              atm_stats_detw_out(:,:,k) = statserosw(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_detc_out) ) then
            do k = 1 , kz
              atm_stats_detc_out(:,:,k) = statserosc(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_erow_out) ) then
            do k = 1 , kz
              atm_stats_erow_out(:,:,k) = statsdetrw(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_eroc_out) ) then
            do k = 1 , kz
              atm_stats_eroc_out(:,:,k) = statsdetrc(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_evw_out) ) then
            do k = 1 , kz
              atm_stats_evw_out(:,:,k) = statsevapw(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_evc_out) ) then
            do k = 1 , kz
              atm_stats_evc_out(:,:,k) = statsevapc(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_con1w_out) ) then
            do k = 1 , kz
              atm_stats_con1w_out(:,:,k) = statscond1w(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_con1c_out) ) then
            do k = 1 , kz
              atm_stats_con1c_out(:,:,k) = statscond1c(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_con2w_out) ) then
            do k = 1 , kz
              atm_stats_con2w_out(:,:,k) = statscond2w(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_con2c_out) ) then
            do k = 1 , kz
              atm_stats_con2c_out(:,:,k) = statscond2c(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_dep_out) ) then
            do k = 1 , kz
              atm_stats_dep_out(:,:,k) = statsdepos(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_melt_out) ) then
            do k = 1 , kz
              atm_stats_melt_out(:,:,k) = statsmelt(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_frz_out) ) then
            do k = 1 , kz
              atm_stats_frz_out(:,:,k) = statsfrz(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_rainev_out) ) then
            do k = 1 , kz
              atm_stats_rainev_out(:,:,k) = statsrainev(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_snowev_out) ) then
            do k = 1 , kz
              atm_stats_snowev_out(:,:,k) = statssnowev(jci1:jci2,ici1:ici2,k)
            end do
          end if
        end if
#endif

        if ( ibltyp == 2 ) then
          if ( associated(atm_tke_out) ) &
            atm_tke_out = atm1%tke(jci1:jci2,ici1:ici2,1:kz)
          if ( associated(atm_kth_out) ) &
            atm_kth_out = uwstateb%kth(jci1:jci2,ici1:ici2,1:kz)
          if ( associated(atm_kzm_out) ) &
            atm_kzm_out = uwstateb%kzm(jci1:jci2,ici1:ici2,1:kz)
        end if

        if ( ichem == 1 .and. iaerosol == 1 .and. iindirect == 2 ) then
          if ( associated(atm_ccnnum_out) ) then
            do k = 1 , kz
              ! convert to 1/cm3
              atm_ccnnum_out(:,:,k) = ccn(jci1:jci2,ici1:ici2,k) * 1.e-6_rkx
            end do
          end if
          if ( idiag == 1 ) then
            if ( associated(atm_qcrit_out) ) then
              do k = 1 , kz
                atm_qcrit_out(:,:,k) = qdiag%qcr(jce1:jce2,ice1:ice2,k)
              end do
            end if
            if ( associated(atm_qincl_out) ) then
              do k = 1 , kz
                atm_qincl_out(:,:,k) = qdiag%qcl(jce1:jce2,ice1:ice2,k)
              end do
            end if
            if ( associated(atm_autoconvr_out) ) then
              do k = 1 , kz
                atm_autoconvr_out(:,:,k) = qdiag%acr(jce1:jce2,ice1:ice2,k)
              end do
            end if
          end if
        end if

        if ( associated(atm_tpr_out) ) then
          atm_tpr_out = (sfs%rainc+sfs%rainnc)/(atmfrq*secph)
        end if
        if ( associated(atm_tsn_out) ) &
          atm_tsn_out = sfs%snownc/(atmfrq*secph)

        if ( associated(atm_tgb_out) ) &
          atm_tgb_out = atm_tgb_out * rsrf_in_atm
        if ( associated(atm_tsw_out) ) then
          where ( mddom%ldmsk > 0 )
            atm_tsw_out = atm_tsw_out * rsrf_in_atm
          elsewhere
            atm_tsw_out = dmissval
          end where
        end if

        ! FAB add tendency diagnostic here
        if ( idiag == 1 ) then
          if ( associated(atm_tten_adh_out) ) then
            do k = 1 , kz
              atm_tten_adh_out(:,:,k) = tdiag%adh(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            tdiag%adh = d_zero
          end if
          if ( associated(atm_tten_adv_out) ) then
            do k = 1 , kz
              atm_tten_adv_out(:,:,k) = tdiag%adv(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            tdiag%adv = d_zero
          end if
          if ( associated(atm_tten_tbl_out) ) then
            do k = 1 , kz
              atm_tten_tbl_out(:,:,k) = tdiag%tbl(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            tdiag%tbl = d_zero
          end if
          if ( associated(atm_tten_dif_out) ) then
            do k = 1 , kz
              atm_tten_dif_out(:,:,k) = tdiag%dif(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            tdiag%dif = d_zero
          end if
          if ( associated(atm_tten_bdy_out) ) then
            do k = 1 , kz
              atm_tten_bdy_out(:,:,k) = tdiag%bdy(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            tdiag%bdy = d_zero
          end if
          if ( associated(atm_tten_con_out) ) then
            do k = 1 , kz
              atm_tten_con_out(:,:,k) = tdiag%con(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            tdiag%con = d_zero
          end if
          if ( associated(atm_tten_adi_out) ) then
            do k = 1 , kz
              atm_tten_adi_out(:,:,k) = tdiag%adi(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            tdiag%adi = d_zero
          end if
          if ( associated(atm_tten_rad_out) ) then
            do k = 1 , kz
              atm_tten_rad_out(:,:,k) = tdiag%rad(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            tdiag%rad = d_zero
          end if
          if ( associated(atm_tten_lsc_out) ) then
            do k = 1 , kz
              atm_tten_lsc_out(:,:,k) = tdiag%lsc(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            tdiag%lsc = d_zero
          end if
          if ( associated(atm_qten_adh_out) ) then
            do k = 1 , kz
              atm_qten_adh_out(:,:,k) = qdiag%adh(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            qdiag%adh = d_zero
          end if
          if ( associated(atm_qten_adv_out) ) then
            do k = 1 , kz
              atm_qten_adv_out(:,:,k) = qdiag%adv(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            qdiag%adv = d_zero
          end if
          if ( associated(atm_qten_tbl_out) ) then
            do k = 1 , kz
              atm_qten_tbl_out(:,:,k) = qdiag%tbl(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            qdiag%tbl = d_zero
          end if
          if ( associated(atm_qten_dif_out) ) then
            do k = 1 , kz
              atm_qten_dif_out(:,:,k) = qdiag%dif(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            qdiag%dif = d_zero
          end if
          if ( associated(atm_qten_bdy_out) ) then
            do k = 1 , kz
              atm_qten_bdy_out(:,:,k) = qdiag%bdy(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            qdiag%bdy = d_zero
          end if
          if ( associated(atm_qten_con_out) ) then
            do k = 1 , kz
              atm_qten_con_out(:,:,k) = qdiag%con(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            qdiag%con = d_zero
          end if
          if ( associated(atm_qten_adi_out) ) then
            do k = 1 , kz
              atm_qten_adi_out(:,:,k) = qdiag%adi(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            qdiag%adi = d_zero
          end if
          if ( associated(atm_qten_rad_out) ) then
            do k = 1 , kz
              atm_qten_rad_out(:,:,k) = qdiag%rad(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            qdiag%rad = d_zero
          end if
          if ( associated(atm_qten_lsc_out) ) then
            do k = 1 , kz
              atm_qten_lsc_out(:,:,k) = qdiag%lsc(jci1:jci2,ici1:ici2,k)/ps_out
            end do
            qdiag%lsc = d_zero
          end if
        end if

        if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + atms%ppb3d(j,i,kz)
            end do
          end do
        else
          ps_out = d_1000*(ps_out+ptop)
        end if

        call write_record_output_stream(atm_stream,idatex)
        if ( myid == italk ) &
          write(stdout,*) 'ATM variables written at ' , tochar(idatex)

        if ( associated(atm_tgb_out) ) atm_tgb_out = d_zero
        if ( associated(atm_tsw_out) ) atm_tsw_out = d_zero
      end if
    end if

    if ( srf_stream > 0 ) then
      if ( ldosrf ) then

        if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + atms%ppb3d(j,i,kz)
            end do
          end do
        else
          ps_out = d_1000*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        end if
        ! Averaged values
        if ( associated(srf_tpr_out) ) &
          srf_tpr_out = srf_tpr_out*rnsrf_for_srffrq
        if ( associated(srf_prcv_out) ) &
          srf_prcv_out = srf_prcv_out*rnsrf_for_srffrq
        if ( associated(srf_zpbl_out) ) &
          srf_zpbl_out = srf_zpbl_out*rnsrf_for_srffrq
        if ( associated(srf_dew_out) .and. associated(srf_evp_out) ) then
          srf_dew_out = -(srf_evp_out*rnsrf_for_srffrq)
          srf_dew_out = max(srf_dew_out, d_zero)
        end if
        if ( associated(srf_evp_out) ) then
          srf_evp_out = srf_evp_out*rnsrf_for_srffrq
          srf_evp_out = max(srf_evp_out, d_zero)
        end if
        if ( associated(srf_scv_out) ) then
          where ( mddom%ldmsk > 0 )
            srf_scv_out = srf_scv_out*rnsrf_for_srffrq
          elsewhere
            srf_scv_out = dmissval
          end where
        end if
        if ( associated(srf_srunoff_out) ) then
          where ( srf_srunoff_out < dlowval )
            srf_srunoff_out = d_zero
          end where
          where ( mddom%ldmsk > 0 )
            srf_srunoff_out = srf_srunoff_out*rnsrf_for_srffrq
          elsewhere
            srf_srunoff_out = dmissval
          end where
        end if
        if ( associated(srf_trunoff_out) ) then
          where ( mddom%ldmsk > 0 )
            where ( srf_trunoff_out < dlowval )
              srf_trunoff_out = d_zero
            end where
            srf_trunoff_out = srf_trunoff_out*rnsrf_for_srffrq
          elsewhere
            srf_trunoff_out = dmissval
          end where
        end if
        if ( associated(srf_sena_out) ) &
          srf_sena_out = srf_sena_out*rnsrf_for_srffrq
        if ( associated(srf_flw_out) ) &
          srf_flw_out = srf_flw_out*rnsrf_for_srffrq
        if ( associated(srf_fsw_out) ) &
          srf_fsw_out = srf_fsw_out*rnsrf_for_srffrq
        if ( associated(srf_fld_out) ) &
          srf_fld_out = srf_fld_out*rnsrf_for_srffrq
        if ( associated(srf_sina_out) ) &
          srf_sina_out = srf_sina_out*rnsrf_for_srffrq
        if ( associated(srf_snowmelt_out) ) &
          srf_snowmelt_out = srf_snowmelt_out*rnsrf_for_srffrq

        call write_record_output_stream(srf_stream,idatex)
        if ( myid == italk ) &
          write(stdout,*) 'SRF variables written at ' , tochar(idatex)

        if ( associated(srf_tpr_out) ) srf_tpr_out = d_zero
        if ( associated(srf_prcv_out) ) srf_prcv_out = d_zero
        if ( associated(srf_zpbl_out) ) srf_zpbl_out = d_zero
        if ( associated(srf_evp_out) ) srf_evp_out = d_zero
        if ( associated(srf_scv_out) ) srf_scv_out = d_zero
        if ( associated(srf_srunoff_out) ) srf_srunoff_out = d_zero
        if ( associated(srf_trunoff_out) ) srf_trunoff_out = d_zero
        if ( associated(srf_sena_out) ) srf_sena_out = d_zero
        if ( associated(srf_flw_out) ) srf_flw_out = d_zero
        if ( associated(srf_fsw_out) ) srf_fsw_out = d_zero
        if ( associated(srf_fld_out) ) srf_fld_out = d_zero
        if ( associated(srf_sina_out) ) srf_sina_out = d_zero
        if ( associated(srf_sund_out) ) srf_sund_out = d_zero
        if ( associated(srf_snowmelt_out) ) srf_snowmelt_out = d_zero
      end if
    end if

    if ( sub_stream > 0 ) then
      if ( ldosub ) then

        sub_ps_out = sub_ps_out*rnsrf_for_subfrq

        if ( associated(sub_evp_out) ) then
          sub_evp_out = sub_evp_out*rnsrf_for_subfrq
          sub_evp_out = max(sub_evp_out, d_zero)
        end if
        if ( associated(sub_scv_out) ) then
          where ( sub_scv_out < dmissval )
            sub_scv_out = sub_scv_out*rnsrf_for_subfrq
          end where
        end if
        if ( associated(sub_sena_out) ) &
          sub_sena_out = sub_sena_out*rnsrf_for_subfrq
        if ( associated(sub_srunoff_out) ) then
          where ( sub_srunoff_out < dmissval )
            sub_srunoff_out = sub_srunoff_out*rnsrf_for_subfrq
          end where
        end if
        if ( associated(sub_trunoff_out) ) then
          where ( sub_trunoff_out < dmissval )
            sub_trunoff_out = sub_trunoff_out*rnsrf_for_subfrq
          end where
        end if

        call write_record_output_stream(sub_stream,idatex)
        if ( myid == italk ) &
          write(stdout,*) 'SUB variables written at ' , tochar(idatex)

        if ( associated(sub_evp_out) ) sub_evp_out = d_zero
        if ( associated(sub_scv_out) ) sub_scv_out = d_zero
        if ( associated(sub_sena_out) ) sub_sena_out = d_zero
        if ( associated(sub_srunoff_out) ) sub_srunoff_out = d_zero
        if ( associated(sub_trunoff_out) ) sub_trunoff_out = d_zero
      end if
    end if

    if ( lak_stream > 0 ) then
      if ( ldolak ) then

        if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + atms%ppb3d(j,i,kz)
            end do
          end do
        else
          ps_out = d_1000*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        end if
        if ( associated(lak_tpr_out) ) &
          lak_tpr_out = lak_tpr_out*rnsrf_for_lakfrq
        if ( associated(lak_scv_out) ) then
          where ( mddom%ldmsk > 0 )
            lak_scv_out = lak_scv_out*rnsrf_for_lakfrq
          elsewhere
            lak_scv_out = dmissval
          end where
        end if
        if ( associated(lak_sena_out) ) &
          lak_sena_out = lak_sena_out*rnsrf_for_lakfrq
        if ( associated(lak_flw_out) ) &
          lak_flw_out = lak_flw_out*rnsrf_for_lakfrq
        if ( associated(lak_fsw_out) ) &
          lak_fsw_out = lak_fsw_out*rnsrf_for_lakfrq
        if ( associated(lak_fld_out) ) &
          lak_fld_out = lak_fld_out*rnsrf_for_lakfrq
        if ( associated(lak_sina_out) ) &
          lak_sina_out = lak_sina_out*rnsrf_for_lakfrq
        if ( associated(lak_evp_out) ) then
          lak_evp_out = lak_evp_out*rnsrf_for_lakfrq
          lak_evp_out = max(lak_evp_out, d_zero)
        end if

        call write_record_output_stream(lak_stream,idatex)
        if ( myid == italk ) &
          write(stdout,*) 'LAK variables written at ' , tochar(idatex)

        if ( associated(lak_tpr_out) )    lak_tpr_out = d_zero
        if ( associated(lak_scv_out) )    lak_scv_out = d_zero
        if ( associated(lak_sena_out) )   lak_sena_out = d_zero
        if ( associated(lak_flw_out) )    lak_flw_out = d_zero
        if ( associated(lak_fsw_out) )    lak_fsw_out = d_zero
        if ( associated(lak_fld_out) )    lak_fld_out = d_zero
        if ( associated(lak_sina_out) )   lak_sina_out = d_zero
        if ( associated(lak_evp_out) )    lak_evp_out = d_zero
      end if
    end if

    if ( opt_stream > 0 ) then
      if ( ldoche ) then
        if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + atms%ppb3d(j,i,kz)
            end do
          end do
        else
          ps_out = d_1000*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        end if
        if ( associated(opt_pp_out) ) then
          do k = 1 , kz
            opt_pp_out(:,:,k) = atm1%pp(jci1:jci2,ici1:ici2,k) / &
                                  sfs%psa(jci1:jci2,ici1:ici2)
          end do
        end if
        if ( associated(opt_acstoarf_out) ) &
          opt_acstoarf_out = opt_acstoarf_out * rnrad_for_chem
        if ( associated(opt_acstsrrf_out) ) &
          opt_acstsrrf_out = opt_acstsrrf_out * rnrad_for_chem
        if ( associated(opt_acstalrf_out) ) &
          opt_acstalrf_out = opt_acstalrf_out * rnrad_for_chem
        if ( associated(opt_acssrlrf_out) ) &
          opt_acssrlrf_out = opt_acssrlrf_out * rnrad_for_chem
        if ( associated(opt_aastoarf_out) ) &
          opt_aastoarf_out = opt_aastoarf_out * rnrad_for_chem
        if ( associated(opt_aastsrrf_out) ) &
          opt_aastsrrf_out = opt_aastsrrf_out * rnrad_for_chem
        if ( associated(opt_aastalrf_out) ) &
          opt_aastalrf_out = opt_aastalrf_out * rnrad_for_chem
        if ( associated(opt_aassrlrf_out) ) &
          opt_aassrlrf_out = opt_aassrlrf_out * rnrad_for_chem

        call write_record_output_stream(opt_stream,idatex)
        if ( myid == italk ) &
          write(stdout,*) 'OPT variables written at ' , tochar(idatex)
        if ( associated(opt_acstoarf_out) ) opt_acstoarf_out = d_zero
        if ( associated(opt_acstsrrf_out) ) opt_acstsrrf_out = d_zero
        if ( associated(opt_acstalrf_out) ) opt_acstalrf_out = d_zero
        if ( associated(opt_acssrlrf_out) ) opt_acssrlrf_out = d_zero
        if ( associated(opt_aastoarf_out) ) opt_acstoarf_out = d_zero
        if ( associated(opt_aastsrrf_out) ) opt_acstsrrf_out = d_zero
        if ( associated(opt_aastalrf_out) ) opt_aastalrf_out = d_zero
        if ( associated(opt_aassrlrf_out) ) opt_aassrlrf_out = d_zero
      end if
    end if

    if ( che_stream > 0 ) then
      if ( ldoche ) then
        if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + atms%ppb3d(j,i,kz)
            end do
          end do
        else
          ps_out = d_1000*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        end if
        if ( associated(che_pp_out) ) then
          do k = 1 , kz
            che_pp_out(:,:,k) = atm1%pp(jci1:jci2,ici1:ici2,k) / &
                                 sfs%psa(jci1:jci2,ici1:ici2)
          end do
        end if
        do itr = 1 , ntr
          call fill_chem_outvars(itr)
          call write_record_output_stream(che_stream,idatex,itr)
        end do
        if ( myid == italk ) &
          write(stdout,*) 'CHE variables written at ' , tochar(idatex)
      end if
    end if

    if ( sts_stream > 0 ) then
      if ( mod(ktau+kstsoff,ksts) == 0 .and. ktau > kstsoff+2 ) then

        if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + atms%ppb3d(j,i,kz)
            end do
          end do
        else
          ps_out = d_1000*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        end if
        if ( associated(sts_pcpavg_out) ) &
          sts_pcpavg_out = sts_pcpavg_out*rnsrf_for_day
        if ( associated(sts_t2avg_out) ) &
          sts_t2avg_out = sts_t2avg_out*rnsrf_for_day
        if ( associated(sts_psavg_out) ) &
          sts_psavg_out = sts_psavg_out*rnsrf_for_day
        if ( associated(sts_srunoff_out) ) then
          where ( mddom%ldmsk > 0 )
            sts_srunoff_out = sts_srunoff_out*rnsrf_for_day
          elsewhere
            sts_srunoff_out = dmissval
          end where
        end if
        if ( associated(sts_trunoff_out) ) then
          where ( mddom%ldmsk > 0 )
            sts_trunoff_out = sts_trunoff_out*rnsrf_for_day
          elsewhere
            sts_trunoff_out = dmissval
          end where
        end if

        call write_record_output_stream(sts_stream,idatex)
        if ( myid == italk ) &
          write(stdout,*) 'STS variables written at ' , tochar(idatex)

        if ( associated(sts_pcpavg_out) )  sts_pcpavg_out  = d_zero
        if ( associated(sts_t2avg_out) )   sts_t2avg_out   = d_zero
        if ( associated(sts_psavg_out) )   sts_psavg_out   = d_zero
        if ( associated(sts_tgmax_out) )   sts_tgmax_out   = -1.e30_rkx
        if ( associated(sts_tgmin_out) )   sts_tgmin_out   =  1.e30_rkx
        if ( associated(sts_t2max_out) )   sts_t2max_out   = -1.e30_rkx
        if ( associated(sts_t2min_out) )   sts_t2min_out   =  1.e30_rkx
        if ( associated(sts_w10max_out) )  sts_w10max_out  = -1.e30_rkx
        if ( associated(sts_psmin_out) )   sts_psmin_out   =  1.e30_rkx
        if ( associated(sts_pcpmax_out) )  sts_pcpmax_out  = -1.e30_rkx
        if ( associated(sts_sund_out) )    sts_sund_out    = d_zero
        if ( associated(sts_srunoff_out) ) sts_srunoff_out = d_zero
        if ( associated(sts_trunoff_out) ) sts_trunoff_out = d_zero

      end if
    end if

    if ( rad_stream > 0 ) then
      if ( ldorad ) then
        if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + atms%ppb3d(j,i,kz)
            end do
          end do
        else
          ps_out = d_1000*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        end if
        if ( associated(rad_pp_out) ) then
          do k = 1 , kz
            rad_pp_out(:,:,k) = atm1%pp(jci1:jci2,ici1:ici2,k)/ &
                             sfs%psa(jci1:jci2,ici1:ici2)
          end do
        end if
        call write_record_output_stream(rad_stream,idatex)
        if ( myid == italk ) &
          write(stdout,*) 'RAD variables written at ' , tochar(idatex)
      end if
    end if

    if ( slaboc_stream > 0 ) then
      if ( ldoslab ) then
        if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + atms%ppb3d(j,i,kz)
            end do
          end do
        else
          ps_out = d_1000*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        end if
        call fill_slaboc_outvars
        call writevar_output_stream(slaboc_stream,v3dvar_slaboc(slab_qflx))
        if ( myid == italk ) then
          write(stdout,*) 'SOM variables written at ' , tochar(idatex)
        end if
      end if
    end if

    if ( ifsave ) then
      if ( ldosav ) then
        call grid_collect(atm1%u,atm1_u_io,jde1,jde2,ide1,ide2,1,kz)
        call grid_collect(atm1%v,atm1_v_io,jde1,jde2,ide1,ide2,1,kz)
        call grid_collect(atm1%t,atm1_t_io,jce1,jce2,ice1,ice2,1,kz)
        call grid_collect(atm1%qx,atm1_qx_io,jce1,jce2,ice1,ice2,1,kz,1,nqx)

        call grid_collect(atm2%u,atm2_u_io,jde1,jde2,ide1,ide2,1,kz)
        call grid_collect(atm2%v,atm2_v_io,jde1,jde2,ide1,ide2,1,kz)
        call grid_collect(atm2%t,atm2_t_io,jce1,jce2,ice1,ice2,1,kz)
        call grid_collect(atm2%qx,atm2_qx_io,jce1,jce2,ice1,ice2,1,kz,1,nqx)

        if ( ibltyp == 2 ) then
          call grid_collect(atm1%tke,atm1_tke_io,jce1,jce2,ice1,ice2,1,kzp1)
          call grid_collect(atm2%tke,atm2_tke_io,jce1,jce2,ice1,ice2,1,kzp1)
          call grid_collect(kpbl,kpbl_io,jci1,jci2,ici1,ici2)
        end if

        if ( idynamic == 2 ) then
          call grid_collect(atm1%pp,atm1_pp_io,jce1,jce2,ice1,ice2,1,kz)
          call grid_collect(atm2%pp,atm2_pp_io,jce1,jce2,ice1,ice2,1,kz)
          call grid_collect(atm1%w,atm1_w_io,jce1,jce2,ice1,ice2,1,kzp1)
          call grid_collect(atm2%w,atm2_w_io,jce1,jce2,ice1,ice2,1,kzp1)
        end if

        call grid_collect(sfs%psa,psa_io,jce1,jce2,ice1,ice2)
        call grid_collect(sfs%psb,psb_io,jce1,jce2,ice1,ice2)
        call grid_collect(sfs%tga,tga_io,jci1,jci2,ici1,ici2)
        call grid_collect(sfs%tgb,tgb_io,jci1,jci2,ici1,ici2)

        call grid_collect(sfs%hfx,hfx_io,jci1,jci2,ici1,ici2)
        call grid_collect(sfs%qfx,qfx_io,jci1,jci2,ici1,ici2)
        call grid_collect(sfs%tgbb,tgbb_io,jci1,jci2,ici1,ici2)
        call grid_collect(sfs%uvdrag,uvdrag_io,jci1,jci2,ici1,ici2)

        if ( ipptls > 0 ) then
          call grid_collect(fcc,fcc_io,jci1,jci2,ici1,ici2,1,kz)
          if ( ipptls == 2 ) then
            call grid_collect(sfs%snownc,snownc_io,jci1,jci2,ici1,ici2)
          end if
        end if
        call grid_collect(heatrt,heatrt_io,jci1,jci2,ici1,ici2,1,kz)
        call grid_collect(o3prof,o3prof_io,jci1,jci2,ici1,ici2,1,kzp1)

        if ( iocnflx == 2 ) then
          call grid_collect(zpbl,zpbl_io,jci1,jci2,ici1,ici2)
        end if
        if ( any(icup == 3) ) then
          call grid_collect(tbase,tbase_io,jci1,jci2,ici1,ici2,1,kz)
          call grid_collect(cldefi,cldefi_io,jci1,jci2,ici1,ici2)
        end if
        if ( any(icup == 4) ) then
          call grid_collect(cbmf2d,cbmf2d_io,jci1,jci2,ici1,ici2)
        end if
        if ( any(icup == 6) ) then
          call grid_collect(kfwavg,kfwavg_io,jci1,jci2,ici1,ici2,1,kz)
        end if

        if ( irrtm == 0 ) then
          call grid_collect(gasabsnxt,gasabsnxt_io,jci1,jci2,ici1,ici2,1,kz,1,4)
          call grid_collect(gasabstot,gasabstot_io, &
                            jci1,jci2,ici1,ici2,1,kzp1,1,kzp1)
          call grid_collect(gasemstot,gasemstot_io,jci1,jci2,ici1,ici2,1,kzp1)
        end if

        call subgrid_collect(lms%sw,sw_io,jci1,jci2,ici1,ici2,1,num_soil_layers)
        call subgrid_collect(lms%gwet,gwet_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(lms%ldew,ldew_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(lms%tgrd,tgrd_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(lms%tgbrd,tgbrd_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(lms%taf,taf_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(lms%tlef,tlef_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(lms%sncv,sncv_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(lms%snag,snag_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(lms%sfice,sfice_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(lms%emisv,emisv_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(lms%scvk,scvk_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(lms%um10,um10_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(lms%swdiralb,swdiralb_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(lms%swdifalb,swdifalb_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(lms%lwdiralb,lwdiralb_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(lms%lwdifalb,lwdifalb_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(mdsub%ldmsk,ldmsk1_io,jci1,jci2,ici1,ici2)

        call grid_collect(solis,solis_io,jci1,jci2,ici1,ici2)
        call grid_collect(solvs,solvs_io,jci1,jci2,ici1,ici2)
        call grid_collect(solvsd,solvsd_io,jci1,jci2,ici1,ici2)
        call grid_collect(solvl,solvl_io,jci1,jci2,ici1,ici2)
        call grid_collect(solvld,solvld_io,jci1,jci2,ici1,ici2)
        call grid_collect(sabveg,sabveg_io,jci1,jci2,ici1,ici2)
        call grid_collect(flw,flw_io,jci1,jci2,ici1,ici2)
        call grid_collect(flwd,flwd_io,jci1,jci2,ici1,ici2)
        call grid_collect(fsw,fsw_io,jci1,jci2,ici1,ici2)
        call grid_collect(sinc,sinc_io,jci1,jci2,ici1,ici2)
        call grid_collect(mddom%ldmsk,ldmsk_io,jci1,jci2,ici1,ici2)

#ifndef CLM
        if ( lakemod == 1 ) then
          call subgrid_collect(lms%eta,eta_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%hi,hi_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%tlake,tlak_io,jci1,jci2,ici1,ici2,1,ndpmax)
        end if
#else
        if ( imask == 2 ) then
          call grid_collect(mddom%lndcat,lndcat_io,jci1,jci2,ici1,ici2)
        end if
#endif
        if ( idcsst == 1 ) then
          call subgrid_collect(lms%sst,sst_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%tskin,tskin_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%deltas,deltas_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%tdeltas,tdeltas_io,jci1,jci2,ici1,ici2)
        end if

        if ( idynamic == 1 ) then
          call grid_collect(dstor,dstor_io,jde1,jde2,ide1,ide2,1,nsplit)
          call grid_collect(hstor,hstor_io,jde1,jde2,ide1,ide2,1,nsplit)
        end if

        if ( ichem == 1 ) then
          call grid_collect(chia,chia_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
          call grid_collect(chib,chib_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
          call grid_collect(rainout,rainout_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
          call grid_collect(washout,washout_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
          call grid_collect(remdrd,remdrd_io,jce1,jce2,ice1,ice2,1,ntr)
          if ( igaschem == 1 .and. ichsolver > 0 ) then
            call grid_collect(chemall,chemall_io,jci1,jci2,ici1,ici2, &
                              1,kz,1,totsp)
            call grid_collect(taucldsp,taucldsp_io,jci1,jci2,ici1,ici2, &
                              0,kz,1,nspi)
          end if

          call grid_collect(ssw2da,ssw2da_io,jci1,jci2,ici1,ici2)
          call grid_collect(sdelt,sdelt_io,jci1,jci2,ici1,ici2)
          call grid_collect(sdelq,sdelq_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfracv2d,sfracv2d_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfracb2d,sfracb2d_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfracs2d,sfracs2d_io,jci1,jci2,ici1,ici2)
          call grid_collect(svegfrac2d,svegfrac2d_io,jci1,jci2,ici1,ici2)
        end if

        if ( islab_ocean == 1 .and. do_restore_sst ) then
          call grid_collect(qflux_restore_sst,qflux_restore_sst_io, &
            jci1,jci2,ici1,ici2,1,12)
        end if
        call write_savefile(idatex)
      end if
    end if

    if ( lfdomonth(idatex) .and. lmidnight(idatex) ) then
      if ( .not. lstartup .and. idatex /= idate2 ) then
        call newoutfiles(idatex)

        ! This must be removed
        ! if ( ifchem .and. myid == iocpu ) then
        !   call prepare_chem_out(idatex,ifrest)
        ! end if

        call checktime(myid)
      end if
    end if

    if ( ldoatm ) then
      sfs%rainc   = d_zero
      sfs%rainnc  = d_zero
      if ( ipptls == 2 ) sfs%snownc  = d_zero
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine output

end module mod_output

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
