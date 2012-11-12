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
module mod_outvars

  use mod_realkinds

  public

  real(rk8) , dimension(:,:) , pointer :: xlon_out => null()
  real(rk8) , dimension(:,:) , pointer :: xlat_out => null()
  real(rk8) , dimension(:,:) , pointer :: topo_out => null()
  real(rk8) , dimension(:,:) , pointer :: mask_out => null()
  real(rk8) , dimension(:,:) , pointer :: ps_out => null()

  real(rk8) , dimension(:,:) , pointer :: xlon_sub_out => null()
  real(rk8) , dimension(:,:) , pointer :: xlat_sub_out => null()
  real(rk8) , dimension(:,:) , pointer :: topo_sub_out => null()
  real(rk8) , dimension(:,:) , pointer :: mask_sub_out => null()

  real(rk8) , dimension(:,:) , pointer :: ps_sub_out => null()

  real(rk8) , dimension(:,:,:) , pointer :: atm_u_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: atm_v_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: atm_t_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: atm_omega_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: atm_qv_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: atm_qc_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: atm_tke_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: atm_kth_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: atm_kzm_out => null()
  real(rk8) , dimension(:,:) , pointer :: atm_tgb_out => null()
  real(rk8) , dimension(:,:) , pointer :: atm_tpr_out => null()
  real(rk8) , dimension(:,:) , pointer :: atm_tsw_out => null()

  real(rk8) , dimension(:,:) , pointer :: uvdrag_out => null()
  real(rk8) , dimension(:,:) , pointer :: tg_out => null()
  real(rk8) , dimension(:,:) , pointer :: tlef_out => null()
  real(rk8) , dimension(:,:) , pointer :: evp_out => null()
  real(rk8) , dimension(:,:) , pointer :: scv_out => null()
  real(rk8) , dimension(:,:) , pointer :: sena_out => null()
  real(rk8) , dimension(:,:) , pointer :: flw_out => null()
  real(rk8) , dimension(:,:) , pointer :: fsw_out => null()
  real(rk8) , dimension(:,:) , pointer :: fld_out => null()
  real(rk8) , dimension(:,:) , pointer :: sina_out => null()
  real(rk8) , dimension(:,:) , pointer :: tpr_out => null()
  real(rk8) , dimension(:,:) , pointer :: prcv_out => null()
  real(rk8) , dimension(:,:) , pointer :: zpbl_out => null()
  real(rk8) , dimension(:,:) , pointer :: aldirs_out => null()
  real(rk8) , dimension(:,:) , pointer :: aldifs_out => null()
  real(rk8) , dimension(:,:) , pointer :: sund_out => null()
  real(rk8) , dimension(:,:) , pointer :: seaice_out => null()

  real(rk8) , dimension(:,:,:) , pointer :: u10m_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: v10m_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: t2m_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: q2m_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: smw_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: runoff_out => null()

  real(rk8) , dimension(:,:) , pointer :: tgmax_out => null()
  real(rk8) , dimension(:,:) , pointer :: tgmin_out => null()
  real(rk8) , dimension(:,:) , pointer :: pcpmax_out => null()
  real(rk8) , dimension(:,:) , pointer :: pcpavg_out => null()
  real(rk8) , dimension(:,:) , pointer :: sundd_out => null()
  real(rk8) , dimension(:,:) , pointer :: psmin_out => null()

  real(rk8) , dimension(:,:,:) , pointer :: t2max_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: t2min_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: t2avg_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: w10max_out => null()

  real(rk8) , dimension(:,:) , pointer :: uvdrag_sub_out => null()
  real(rk8) , dimension(:,:) , pointer :: tg_sub_out => null()
  real(rk8) , dimension(:,:) , pointer :: tlef_sub_out => null()
  real(rk8) , dimension(:,:) , pointer :: evp_sub_out => null()
  real(rk8) , dimension(:,:) , pointer :: scv_sub_out => null()
  real(rk8) , dimension(:,:) , pointer :: sena_sub_out => null()
  real(rk8) , dimension(:,:) , pointer :: tlake_sub_out => null()

  real(rk8) , dimension(:,:,:) , pointer :: u10m_sub_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: v10m_sub_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: t2m_sub_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: q2m_sub_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: smw_sub_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: runoff_sub_out => null()

  real(rk8) , dimension(:,:) , pointer :: frsa_out => null()
  real(rk8) , dimension(:,:) , pointer :: frla_out => null()
  real(rk8) , dimension(:,:) , pointer :: clrst_out => null()
  real(rk8) , dimension(:,:) , pointer :: clrss_out => null()
  real(rk8) , dimension(:,:) , pointer :: clrls_out => null()
  real(rk8) , dimension(:,:) , pointer :: clrlt_out => null()
  real(rk8) , dimension(:,:) , pointer :: solin_out => null()
  real(rk8) , dimension(:,:) , pointer :: sabtp_out => null()
  real(rk8) , dimension(:,:) , pointer :: totcf_out => null()
  real(rk8) , dimension(:,:) , pointer :: totcl_out => null()
  real(rk8) , dimension(:,:) , pointer :: totci_out => null()
  real(rk8) , dimension(:,:) , pointer :: firtp_out => null()

  real(rk8) , dimension(:,:,:) , pointer :: cld_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: clwp_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: qrs_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: qrl_out => null()

  real(rk8) , dimension(:,:) , pointer :: aveice_out => null()
  real(rk8) , dimension(:,:) , pointer :: hsnow_out => null()
  real(rk8) , dimension(:,:,:) , pointer :: tlake_out => null()

end module mod_outvars
