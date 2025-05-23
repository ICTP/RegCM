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

module mod_outvars

  use mod_realkinds

  implicit none

  public

  real(rkx), dimension(:,:), pointer, contiguous :: xlon_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: xlat_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: topo_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: mask_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: area_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: ps_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: p0_out => null()

  real(rkx), dimension(:,:), pointer, contiguous :: sub_xlon_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sub_xlat_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sub_topo_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sub_mask_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sub_area_out => null()

  real(rkx), dimension(:,:), pointer, contiguous :: sub_ps_out => null()

  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_u_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_v_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_w_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_t_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_pp_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_pai_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_omega_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qv_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qc_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qr_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qi_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qs_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qg_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qh_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_nn_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_nc_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_nr_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_rh_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_zf_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_zh_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_pf_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_ph_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_q_detr_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_rainls_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_raincc_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_tke_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_kth_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_kzm_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_tten_adh_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_tten_adv_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_tten_tbl_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_tten_dif_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_tten_bdy_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_tten_con_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_tten_adi_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_tten_rad_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_tten_lsc_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qten_adh_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qten_adv_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qten_tbl_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qten_dif_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qten_bdy_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qten_con_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qten_adi_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qten_rad_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qten_lsc_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qcrit_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_ccnnum_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_qincl_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_autoconvr_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_smw_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_tsoil_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: atm_tgb_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: atm_tpr_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: atm_mrso_out => null()
  ! stats
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_supw_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_supc_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_detw_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_detc_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_erow_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_eroc_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_evw_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_evc_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_con1w_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_con1c_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_dep_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_melt_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_frz_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_rainev_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_snowev_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_autocw_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: atm_stats_autocc_out => null()

  real(rkx), dimension(:,:), pointer, contiguous :: shf_pcpmax_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: shf_twetb_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: shf_pcpavg_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: shf_pcprcv_out => null()

  real(rkx), dimension(:,:), pointer, contiguous :: srf_uvdrag_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_taux_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_tauy_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_ustar_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_zo_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_rhoa_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_tg_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_tlef_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_evp_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_dew_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_scv_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_sena_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_lena_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_flw_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_fsw_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_uflw_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_ufsw_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_fld_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_sina_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_tpr_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_prcv_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_snow_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_hail_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_grau_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_zpbl_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_aldirs_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_aldifs_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_sund_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_snowmelt_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_seaice_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_srunoff_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_trunoff_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_totcf_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_wspd_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_mslp_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_evpot_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_pcpmax_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_twetb_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_tprw_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_cape_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_cin_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_li_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_mrsos_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_htindx_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: srf_hfso_out => null()

  real(rkx), dimension(:,:,:), pointer, contiguous :: srf_u10m_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: srf_v10m_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: srf_t2m_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: srf_q2m_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: srf_rh2m_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: srf_smw_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: srf_tsoil_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: srf_ua50_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: srf_va50_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: srf_ta50_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: srf_hus50_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: srf_ua100_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: srf_va100_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: srf_ua150_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: srf_va150_out => null()

  real(rkx), dimension(:,:), pointer, contiguous :: sts_tgmax_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sts_tgmin_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sts_pcpmax_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sts_pcpavg_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sts_sund_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sts_psavg_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sts_psmin_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sts_srunoff_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sts_trunoff_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sts_wsgsmax_out => null()

  real(rkx), dimension(:,:,:), pointer, contiguous :: sts_t2max_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: sts_t2min_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: sts_t2avg_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: sts_w10max_out => null()

  real(rkx), dimension(:,:), pointer, contiguous :: sub_uvdrag_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sub_tg_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sub_tlef_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sub_evp_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sub_scv_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sub_sena_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sub_srunoff_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: sub_trunoff_out => null()

  real(rkx), dimension(:,:,:), pointer, contiguous :: sub_u10m_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: sub_v10m_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: sub_t2m_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: sub_q2m_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: sub_smw_out => null()

  real(rkx), dimension(:,:), pointer, contiguous :: rad_frsa_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: rad_frla_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: rad_clrst_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: rad_clrss_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: rad_clrls_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: rad_clrlt_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: rad_solin_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: rad_solout_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: rad_totwv_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: rad_totcl_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: rad_totci_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: rad_lwout_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: rad_higcl_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: rad_midcl_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: rad_lowcl_out => null()

  real(rkx), dimension(:,:,:), pointer, contiguous :: rad_pp_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: rad_pai_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: rad_cld_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: rad_clwp_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: rad_qrs_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: rad_qrl_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: rad_o3_out => null()
  real(rkx), dimension(:,:,:,:), pointer, contiguous :: rad_taucl_out => null()
  real(rkx), dimension(:,:,:,:), pointer, contiguous :: rad_tauci_out => null()

  real(rkx), dimension(:,:), pointer, contiguous :: lak_tg_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: lak_tpr_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: lak_scv_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: lak_sena_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: lak_sina_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: lak_fsw_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: lak_flw_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: lak_fld_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: lak_evp_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: lak_ice_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: lak_aldirs_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: lak_aldifs_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: lak_tlake_out => null()

  real(rkx), dimension(:,:), pointer, contiguous :: opt_acstoarf_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: opt_acstsrrf_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: opt_aastoarf_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: opt_aastsrrf_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: opt_acstalrf_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: opt_acssrlrf_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: opt_aastalrf_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: opt_aassrlrf_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: opt_aod_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: opt_pp_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: opt_pai_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: opt_aext8_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: opt_assa8_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: opt_agfu8_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: opt_deltaz_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: opt_ncon_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: opt_surf_out => null()

  real(rkx), dimension(:,:), pointer, contiguous :: che_wdrflx_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: che_wdcflx_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: che_ddflx_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: che_emflx_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: che_ddvel_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: che_burden_out => null()
  real(rkx), dimension(:,:), pointer, contiguous :: che_pblten_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_pp_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_pai_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_mixrat_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_cheten_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_advhten_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_advvten_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_difhten_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_cuten_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_tuten_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_raiten_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_wasten_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_bdyten_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_sedten_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_emten_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_chgact_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_ncon_out => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: che_airden_out => null()

  real(rkx), dimension(:,:,:), pointer, contiguous :: slab_qflx_out => null()

end module mod_outvars
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
