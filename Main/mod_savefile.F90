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

module mod_savefile

  use mod_nchelper , only : regcm_vartype
  use mod_intkinds
  use mod_realkinds
  use mod_date
  use mod_stdio
  use mod_dynparam
  use mod_runparams
  use mod_mpmessage
  use mod_mppparam
  use mod_memutil
  use mod_atm_interface , only : tmask
  use mod_lm_interface
  use mod_che_interface
  use mod_che_mppio
  use mod_massck
  use netcdf

  implicit none

  private

#ifndef QUAD_PRECISION
  integer , parameter :: wrkp = rk8
#else
  integer , parameter :: wrkp = rk16
#endif

  public :: allocate_mod_savefile
  public :: read_savefile
  public :: write_savefile

  integer(ik4) , parameter :: maxdims = 18
  integer(ik4) , parameter :: maxvars = 128
  integer(ik4) :: ncstatus

  integer(ik4) , parameter :: idjcross = 1
  integer(ik4) , parameter :: idicross = 2
  integer(ik4) , parameter :: idjdot = 3
  integer(ik4) , parameter :: ididot = 4
  integer(ik4) , parameter :: idkh = 5
  integer(ik4) , parameter :: idkf = 6
  integer(ik4) , parameter :: idnsplit = 7
  integer(ik4) , parameter :: idnnsg = 8
  integer(ik4) , parameter :: idmonth = 9
  integer(ik4) , parameter :: idnqx = 10
  integer(ik4) , parameter :: idspi = 11
  integer(ik4) , parameter :: idspw = 12
  integer(ik4) , parameter :: idntr = 13
  integer(ik4) , parameter :: idtotsp = 14
  integer(ik4) , parameter :: iddpt = 15
  integer(ik4) , parameter :: ikern = 16
  integer(ik4) , parameter :: indep = 17
  integer(ik4) , parameter :: idndu = 18

  integer(ik4) , public , pointer , dimension(:,:,:) :: ldmsk1_io
  integer(ik4) , public , pointer , dimension(:,:) :: ldmsk_io

  real(rkx) , public , pointer , dimension(:,:,:) :: atm_u_io
  real(rkx) , public , pointer , dimension(:,:,:) :: atm_v_io
  real(rkx) , public , pointer , dimension(:,:,:) :: atm_w_io
  real(rkx) , public , pointer , dimension(:,:,:) :: atm_pai_io
  real(rkx) , public , pointer , dimension(:,:,:) :: atm_t_io
  real(rkx) , public , pointer , dimension(:,:,:) :: atm_tke_io
  real(rkx) , public , pointer , dimension(:,:,:,:) :: atm_qx_io
  real(rkx) , public , pointer , dimension(:,:) :: ps_io

  real(rkx) , public , pointer , dimension(:,:,:) :: atm1_u_io
  real(rkx) , public , pointer , dimension(:,:,:) :: atm2_u_io
  real(rkx) , public , pointer , dimension(:,:,:) :: atm1_v_io
  real(rkx) , public , pointer , dimension(:,:,:) :: atm2_v_io
  real(rkx) , public , pointer , dimension(:,:,:) :: atm1_t_io
  real(rkx) , public , pointer , dimension(:,:,:) :: atm2_t_io
  real(rkx) , public , pointer , dimension(:,:,:) :: atm1_w_io
  real(rkx) , public , pointer , dimension(:,:,:) :: atm2_w_io
  real(rkx) , public , pointer , dimension(:,:,:) :: atm1_pp_io
  real(rkx) , public , pointer , dimension(:,:,:) :: atm2_pp_io
  real(rkx) , public , pointer , dimension(:,:,:,:) :: atm1_qx_io
  real(rkx) , public , pointer , dimension(:,:,:,:) :: atm2_qx_io
  real(rkx) , public , pointer , dimension(:,:,:) :: atm1_tke_io
  real(rkx) , public , pointer , dimension(:,:,:) :: atm2_tke_io

  real(rkx) , public , pointer , dimension(:,:,:) :: tke_pbl_io

  real(rkx) , public , pointer , dimension(:,:) :: psa_io
  real(rkx) , public , pointer , dimension(:,:) :: psb_io

  real(rkx) , public , pointer , dimension(:,:) :: hfx_io
  real(rkx) , public , pointer , dimension(:,:) :: qfx_io
  real(rkx) , public , pointer , dimension(:,:) :: tgbb_io
  real(rkx) , public , pointer , dimension(:,:) :: uvdrag_io
  real(rkx) , public , pointer , dimension(:,:) :: q2m_io
  real(rkx) , public , pointer , dimension(:,:) :: u10m_io
  real(rkx) , public , pointer , dimension(:,:) :: v10m_io
  real(rkx) , public , pointer , dimension(:,:) :: w10m_io
  real(rkx) , public , pointer , dimension(:,:) :: br_io
  real(rkx) , public , pointer , dimension(:,:) :: ram_io
  real(rkx) , public , pointer , dimension(:,:) :: rah_io
  real(rkx) , public , pointer , dimension(:,:) :: ustar_io
  real(rkx) , public , pointer , dimension(:,:) :: zo_io
  real(rkx) , public , pointer , dimension(:,:) :: dsrnof_io
  real(rkx) , public , pointer , dimension(:,:) :: dtrnof_io

  real(rkx) , public , pointer , dimension(:,:,:) :: ldew_io
  real(rkx) , public , pointer , dimension(:,:,:) :: snag_io
  real(rkx) , public , pointer , dimension(:,:,:) :: sncv_io
  real(rkx) , public , pointer , dimension(:,:,:) :: sfice_io
  real(rkx) , public , pointer , dimension(:,:,:) :: gwet_io
  real(rkx) , public , pointer , dimension(:,:,:,:) :: sw_io
  real(rkx) , public , pointer , dimension(:,:,:) :: tsoi_io
  real(rkx) , public , pointer , dimension(:,:,:) :: swvol_io
  real(rkx) , public , pointer , dimension(:,:,:) :: taf_io
  real(rkx) , public , pointer , dimension(:,:,:) :: tgrd_io
  real(rkx) , public , pointer , dimension(:,:,:) :: tgbrd_io
  real(rkx) , public , pointer , dimension(:,:,:) :: tlef_io
  real(rkx) , public , pointer , dimension(:,:,:) :: emisv_io
  real(rkx) , public , pointer , dimension(:,:,:) :: um10_io
  real(rkx) , public , pointer , dimension(:,:,:) :: eta_io
  real(rkx) , public , pointer , dimension(:,:,:) :: hi_io
  real(rkx) , public , pointer , dimension(:,:,:,:) :: tlak_io

  real(rkx) , public , pointer , dimension(:,:) :: flw_io
  real(rkx) , public , pointer , dimension(:,:) :: flwd_io
  real(rkx) , public , pointer , dimension(:,:) :: fsw_io
  real(rkx) , public , pointer , dimension(:,:) :: sabveg_io
  real(rkx) , public , pointer , dimension(:,:) :: totcf_io
  real(rkx) , public , pointer , dimension(:,:) :: sinc_io
  real(rkx) , public , pointer , dimension(:,:) :: solis_io
  real(rkx) , public , pointer , dimension(:,:) :: solvs_io
  real(rkx) , public , pointer , dimension(:,:) :: solvsd_io
  real(rkx) , public , pointer , dimension(:,:) :: solvl_io
  real(rkx) , public , pointer , dimension(:,:) :: solvld_io

  real(rkx) , public , pointer , dimension(:,:,:) :: sst_io
  real(rkx) , public , pointer , dimension(:,:,:) :: tskin_io
  real(rkx) , public , pointer , dimension(:,:,:) :: tdeltas_io
  real(rkx) , public , pointer , dimension(:,:,:) :: deltas_io

  integer(ik4) , public , pointer , dimension(:,:) :: kpbl_io
  real(rkx) , public , pointer , dimension(:,:) :: zpbl_io
  real(rkx) , public , pointer , dimension(:,:) :: myjsf_uz0_io
  real(rkx) , public , pointer , dimension(:,:) :: myjsf_vz0_io
  real(rkx) , public , pointer , dimension(:,:) :: myjsf_thz0_io
  real(rkx) , public , pointer , dimension(:,:) :: myjsf_qz0_io

  real(rkx) , public , pointer , dimension(:,:) :: cbmf2d_io

  real(rkx) , public , pointer , dimension(:,:,:) :: fcc_io

  real(rkx) , public , pointer , dimension(:,:,:,:) :: gasabsnxt_io
  real(rkx) , public , pointer , dimension(:,:,:,:) :: gasabstot_io
  real(rkx) , public , pointer , dimension(:,:,:) :: gasemstot_io

  real(rkx) , public , pointer , dimension(:,:,:) :: heatrt_io
  real(rkx) , public , pointer , dimension(:,:,:) :: o3prof_io

  real(rkx) , public , pointer , dimension(:,:,:) :: dstor_io
  real(rkx) , public , pointer , dimension(:,:,:) :: hstor_io

  real(rkx) , public , pointer , dimension(:,:) :: cldefi_io

  real(rkx) , public , pointer , dimension(:,:,:) :: cu_avg_ww_io

  real(rkx) , public , pointer , dimension(:,:,:) :: qflux_restore_sst_io

#ifdef CLM
  real(rkx) , public , pointer , dimension(:,:) :: lndcat_io
#endif

  real(rkx) , public , pointer , dimension(:,:,:) :: swalb_io
  real(rkx) , public , pointer , dimension(:,:,:) :: lwalb_io
  real(rkx) , public , pointer , dimension(:,:,:) :: swdiralb_io
  real(rkx) , public , pointer , dimension(:,:,:) :: swdifalb_io
  real(rkx) , public , pointer , dimension(:,:,:) :: lwdiralb_io
  real(rkx) , public , pointer , dimension(:,:,:) :: lwdifalb_io

  real(rkx) , public , pointer , dimension(:,:,:,:) :: tmp_io

  interface myputvar
    module procedure myputvar2dd
    module procedure myputvar2ddf
    module procedure myputvar3dd
    module procedure myputvar4dd
    module procedure myputvar1dif
    module procedure myputvar2di
    module procedure myputvar3di
  end interface myputvar

  interface mygetvar
    module procedure mygetvar2dd
    module procedure mygetvar2ddf
    module procedure mygetvar3dd
    module procedure mygetvar4dd
    module procedure mygetvar1dif
    module procedure mygetvar2di
    module procedure mygetvar3di
  end interface mygetvar

  contains

  subroutine allocate_mod_savefile
    implicit none

    if ( do_parallel_save ) then
      if ( idynamic == 3 ) then
        call getmem3d(atm_u_io,jde1,jde2,ice1,ice2,1,kz,'atm_u_io')
        call getmem3d(atm_v_io,jce1,jce2,ide1,ide2,1,kz,'atm_v_io')
        call getmem3d(atm_t_io,jce1,jce2,ice1,ice2,1,kz,'atm_t_io')
        call getmem3d(atm_pai_io,jce1,jce2,ice1,ice2,1,kz,'atm_pai_io')
        call getmem3d(atm_w_io,jce1,jce2,ice1,ice2,1,kzp1,'atm_w_io')
        call getmem4d(atm_qx_io,jce1,jce2,ice1,ice2,1,kz,1,nqx,'atm_qx_io')
        if ( ibltyp == 2 ) then
          call getmem3d(atm_tke_io,jce1,jce2,ice1,ice2,1,kzp1,'atm_tke_io')
        end if
        call getmem2d(ps_io,jce1,jce2,ice1,ice2,'ps_io')
      else
        call getmem3d(atm1_u_io,jde1,jde2,ide1,ide2,1,kz,'atm1_u_io')
        call getmem3d(atm1_v_io,jde1,jde2,ide1,ide2,1,kz,'atm1_v_io')
        call getmem3d(atm1_t_io,jce1,jce2,ice1,ice2,1,kz,'atm1_t_io')
        call getmem4d(atm1_qx_io,jce1,jce2,ice1,ice2,1,kz,1,nqx,'atm1_qx_io')
        call getmem3d(atm2_u_io,jde1,jde2,ide1,ide2,1,kz,'atm2_u_io')
        call getmem3d(atm2_v_io,jde1,jde2,ide1,ide2,1,kz,'atm2_v_io')
        call getmem3d(atm2_t_io,jce1,jce2,ice1,ice2,1,kz,'atm2_t_io')
        call getmem4d(atm2_qx_io,jce1,jce2,ice1,ice2,1,kz,1,nqx,'atm2_qx_io')
        if ( ibltyp == 2 ) then
          call getmem3d(atm1_tke_io,jce1,jce2,ice1,ice2,1,kzp1,'atm1_tke_io')
          call getmem3d(atm2_tke_io,jce1,jce2,ice1,ice2,1,kzp1,'atm2_tke_io')
        end if
        call getmem2d(psa_io,jce1,jce2,ice1,ice2,'psa_io')
        call getmem2d(psb_io,jce1,jce2,ice1,ice2,'psb_io')
      end if
      if ( ibltyp == 4 ) then
        call getmem3d(tke_pbl_io,jci1,jci2,ici1,ici2,1,kz,'tke_pbl_io')
      end if
      if ( idynamic == 2 ) then
        call getmem3d(atm1_w_io,jce1,jce2,ice1,ice2,1,kzp1,'atm1_w_io')
        call getmem3d(atm2_w_io,jce1,jce2,ice1,ice2,1,kzp1,'atm2_w_io')
        call getmem3d(atm1_pp_io,jce1,jce2,ice1,ice2,1,kz,'atm1_pp_io')
        call getmem3d(atm2_pp_io,jce1,jce2,ice1,ice2,1,kz,'atm2_pp_io')
      end if

      call getmem2d(hfx_io,jci1,jci2,ici1,ici2,'hfx_io')
      call getmem2d(qfx_io,jci1,jci2,ici1,ici2,'qfx_io')
      call getmem2d(tgbb_io,jci1,jci2,ici1,ici2,'tgbb_io')
      call getmem2d(uvdrag_io,jci1,jci2,ici1,ici2,'uvdrag_io')
      call getmem2d(q2m_io,jci1,jci2,ici1,ici2,'q2m_io')
      call getmem2d(u10m_io,jci1,jci2,ici1,ici2,'u10m_io')
      call getmem2d(v10m_io,jci1,jci2,ici1,ici2,'v10m_io')
      call getmem2d(w10m_io,jci1,jci2,ici1,ici2,'w10m_io')
      call getmem2d(br_io,jci1,jci2,ici1,ici2,'br_io')
      call getmem2d(ram_io,jci1,jci2,ici1,ici2,'ram_io')
      call getmem2d(rah_io,jci1,jci2,ici1,ici2,'rah_io')
      call getmem2d(ustar_io,jci1,jci2,ici1,ici2,'ustar_io')
      call getmem2d(zo_io,jci1,jci2,ici1,ici2,'zo_io')
      if ( iocncpl == 1 .or. iwavcpl == 1 ) then
        call getmem2d(dsrnof_io,jci1,jci2,ici1,ici2,'dsrnof_io')
        call getmem2d(dtrnof_io,jci1,jci2,ici1,ici2,'dtrnof_io')
      end if
#ifdef CLM45
      if ( ichem == 1 ) then
        call getmem3d(tsoi_io,jci1,jci2,ici1,ici2, &
                            1,num_soil_layers,'tsoi_io')
        call getmem3d(swvol_io,jci1,jci2,ici1,ici2, &
                            1,num_soil_layers,'swvol_io')
      end if
#else
      call getmem3d(gwet_io,1,nnsg,jci1,jci2,ici1,ici2,'gwet_io')
      call getmem3d(ldew_io,1,nnsg,jci1,jci2,ici1,ici2,'ldew_io')
      call getmem3d(taf_io,1,nnsg,jci1,jci2,ici1,ici2,'taf_io')
#endif
      call getmem3d(sncv_io,1,nnsg,jci1,jci2,ici1,ici2,'sncv_io')
      call getmem3d(sfice_io,1,nnsg,jci1,jci2,ici1,ici2,'sfice_io')
      call getmem3d(snag_io,1,nnsg,jci1,jci2,ici1,ici2,'snag_io')
      call getmem4d(sw_io,1,nnsg,jci1,jci2,ici1,ici2, &
                          1,num_soil_layers,'sw_io')
      call getmem3d(tgrd_io,1,nnsg,jci1,jci2,ici1,ici2,'tgrd_io')
      call getmem3d(tgbrd_io,1,nnsg,jci1,jci2,ici1,ici2,'tgbrd_io')
      call getmem3d(tlef_io,1,nnsg,jci1,jci2,ici1,ici2,'tlef_io')
      call getmem3d(emisv_io,1,nnsg,jci1,jci2,ici1,ici2,'emisv_io')
      call getmem3d(um10_io,1,nnsg,jci1,jci2,ici1,ici2,'um10_io')
      call getmem3d(ldmsk1_io,1,nnsg,jci1,jci2,ici1,ici2,'ldmsk1')
      call getmem2d(ldmsk_io,jci1,jci2,ici1,ici2,'ldmsk_io')
      call getmem2d(flw_io,jci1,jci2,ici1,ici2,'flw_io')
      call getmem2d(flwd_io,jci1,jci2,ici1,ici2,'flwd_io')
      call getmem2d(fsw_io,jci1,jci2,ici1,ici2,'fsw_io')
      call getmem2d(sabveg_io,jci1,jci2,ici1,ici2,'sabveg_io')
      call getmem2d(totcf_io,jci1,jci2,ici1,ici2,'totcf_io')
      call getmem2d(sinc_io,jci1,jci2,ici1,ici2,'sinc_io')
      call getmem2d(solis_io,jci1,jci2,ici1,ici2,'solis_io')
      call getmem2d(solvs_io,jci1,jci2,ici1,ici2,'solvs_io')
      call getmem2d(solvsd_io,jci1,jci2,ici1,ici2,'solvsd_io')
      call getmem2d(solvl_io,jci1,jci2,ici1,ici2,'solvl_io')
      call getmem2d(solvld_io,jci1,jci2,ici1,ici2,'solvld_io')
      if ( ibltyp == 2 ) then
        call getmem2d(kpbl_io,jci1,jci2,ici1,ici2,'kpbl_io')
      else if ( ibltyp == 4 ) then
        call getmem2d(myjsf_uz0_io,jci1,jci2,ici1,ici2,'myjsf_uz0_io')
        call getmem2d(myjsf_vz0_io,jci1,jci2,ici1,ici2,'myjsf_vz0_io')
        call getmem2d(myjsf_thz0_io,jci1,jci2,ici1,ici2,'myjsf_thz0_io')
        call getmem2d(myjsf_qz0_io,jci1,jci2,ici1,ici2,'myjsf_qz0_io')
        call getmem2d(kpbl_io,jci1,jci2,ici1,ici2,'kpbl_io')
      end if
      if ( idcsst == 1 ) then
        call getmem3d(sst_io,1,nnsg,jci1,jci2,ici1,ici2,'sst_io')
        call getmem3d(tskin_io,1,nnsg,jci1,jci2,ici1,ici2,'tskin_io')
        call getmem3d(deltas_io,1,nnsg,jci1,jci2,ici1,ici2,'deltas_io')
        call getmem3d(tdeltas_io,1,nnsg,jci1,jci2,ici1,ici2,'tdeltas_io')
      end if
      if ( any(icup == 4) ) then
        call getmem2d(cbmf2d_io,jci1,jci2,ici1,ici2,'cbmf2d_io')
      end if
      if ( ipptls > 0 ) then
        call getmem3d(fcc_io,jci1,jci2,ici1,ici2,1,kz,'fcc_io')
      end if

      if ( idynamic == 1 ) then
        call getmem3d(dstor_io,jde1,jde2,ide1,ide2,1,nsplit,'dstor_io')
        call getmem3d(hstor_io,jde1,jde2,ide1,ide2,1,nsplit,'hstor_io')
      end if

      if ( irrtm == 0 ) then
        call getmem4d(gasabsnxt_io,jci1,jci2,ici1,ici2,1,kz,1,4,'gasabsnxt_io')
        call getmem4d(gasabstot_io,jci1,jci2,ici1,ici2, &
                      1,kzp1,1,kzp1,'gasabstot_io')
        call getmem3d(gasemstot_io,jci1,jci2,ici1,ici2,1,kzp1,'gasemstot_io')
      end if

      call getmem3d(heatrt_io,jci1,jci2,ici1,ici2,1,kz,'heatrt_io')
      call getmem3d(o3prof_io,jci1,jci2,ici1,ici2,1,kzp1,'o3prof_io')

      if ( islab_ocean == 1 .and. do_restore_sst ) then
        call getmem3d(qflux_restore_sst_io,jci1,jci2,ici1,ici2, &
                      1,12,'qflux_restore_sst_io')
      end if

      if ( iocnflx == 2 .or. ibltyp == 3 ) then
        call getmem2d(zpbl_io,jci1,jci2,ici1,ici2,'zpbl_io')
      end if
      if ( any(icup == 3) ) then
        call getmem2d(cldefi_io,jci1,jci2,ici1,ici2,'cldefi_io')
      end if
      if ( any(icup == 6) .or. any(icup == 5) ) then
        call getmem3d(cu_avg_ww_io,jci1,jci2,ici1,ici2,1,kz,'cu_avg_ww_io')
      end if
#ifdef CLM
      if ( imask == 2 ) then
        call getmem2d(lndcat_io,jce1,jce2,ice1,ice2,'lndcat_io')
      end if
#else
      if ( lakemod == 1 ) then
        call getmem3d(eta_io,1,nnsg,jci1,jci2,ici1,ici2,'eta_io')
        call getmem3d(hi_io,1,nnsg,jci1,jci2,ici1,ici2,'hi_io')
        call getmem4d(tlak_io,1,nnsg,jci1,jci2,ici1,ici2,1,ndpmax,'tlak_io')
      end if
#endif
      call getmem3d(swalb_io,1,nnsg,jci1,jci2,ici1,ici2,'swalb')
      call getmem3d(lwalb_io,1,nnsg,jci1,jci2,ici1,ici2,'lwalb')
      call getmem3d(swdiralb_io,1,nnsg,jci1,jci2,ici1,ici2,'swdiralb')
      call getmem3d(swdifalb_io,1,nnsg,jci1,jci2,ici1,ici2,'swdifalb')
      call getmem3d(lwdiralb_io,1,nnsg,jci1,jci2,ici1,ici2,'lwdiralb')
      call getmem3d(lwdifalb_io,1,nnsg,jci1,jci2,ici1,ici2,'lwdifalb')
      return
    end if

    if ( myid == iocpu ) then
      if ( idynamic == 3 ) then
        call getmem3d(atm_u_io,jdot1,jdot2,icross1,icross2,1,kz,'atm_u_io')
        call getmem3d(atm_v_io,jcross1,jcross2,idot1,idot2,1,kz,'atm_v_io')
        call getmem3d(atm_t_io,jcross1,jcross2,icross1,icross2,1,kz,'atm_t_io')
        call getmem3d(atm_pai_io,jcross1,jcross2, &
                                 icross1,icross2,1,kz,'atm_pai_io')
        call getmem3d(atm_w_io,jcross1,jcross2, &
                               icross1,icross2,1,kzp1,'atm_w_io')
        call getmem4d(atm_qx_io,jcross1,jcross2, &
                      icross1,icross2,1,kz,1,nqx,'atm_qx_io')
        if ( ibltyp == 2 ) then
          call getmem3d(atm_tke_io,jcross1,jcross2, &
                        icross1,icross2,1,kzp1,'atm_tke_io')
        end if
        call getmem2d(ps_io,jcross1,jcross2,icross1,icross2,'ps_io')
      else
        call getmem3d(atm1_u_io,jdot1,jdot2,idot1,idot2,1,kz,'atm1_u_io')
        call getmem3d(atm1_v_io,jdot1,jdot2,idot1,idot2,1,kz,'atm1_v_io')
        call getmem3d(atm1_t_io,jcross1,jcross2, &
                      icross1,icross2,1,kz,'atm1_t_io')
        call getmem4d(atm1_qx_io,jcross1,jcross2, &
                      icross1,icross2,1,kz,1,nqx,'atm1_qx_io')
        call getmem3d(atm2_u_io,jdot1,jdot2,idot1,idot2,1,kz,'atm2_u_io')
        call getmem3d(atm2_v_io,jdot1,jdot2,idot1,idot2,1,kz,'atm2_v_io')
        call getmem3d(atm2_t_io,jcross1,jcross2, &
                      icross1,icross2,1,kz,'atm2_t_io')
        call getmem4d(atm2_qx_io,jcross1,jcross2, &
                      icross1,icross2,1,kz,1,nqx,'atm2_qx_io')
        if ( ibltyp == 2 ) then
          call getmem3d(atm1_tke_io,jcross1,jcross2, &
                        icross1,icross2,1,kzp1,'atm1_tke_io')
          call getmem3d(atm2_tke_io,jcross1,jcross2, &
                        icross1,icross2,1,kzp1,'atm2_tke_io')
        end if
        call getmem2d(psa_io,jcross1,jcross2,icross1,icross2,'psa_io')
        call getmem2d(psb_io,jcross1,jcross2,icross1,icross2,'psb_io')
      end if
      if ( ibltyp == 4 ) then
        call getmem3d(tke_pbl_io,jcross1,jcross2, &
                      icross1,icross2,1,kz,'tke_pbl_io')
      end if
      if ( idynamic == 2 ) then
        call getmem3d(atm1_w_io,jcross1,jcross2, &
                      icross1,icross2,1,kzp1,'atm1_w_io')
        call getmem3d(atm2_w_io,jcross1,jcross2, &
                      icross1,icross2,1,kzp1,'atm2_w_io')
        call getmem3d(atm1_pp_io,jcross1,jcross2, &
                      icross1,icross2,1,kz,'atm1_pp_io')
        call getmem3d(atm2_pp_io,jcross1,jcross2, &
                      icross1,icross2,1,kz,'atm2_pp_io')
      end if
      call getmem2d(hfx_io,jcross1,jcross2,icross1,icross2,'hfx_io')
      call getmem2d(qfx_io,jcross1,jcross2,icross1,icross2,'qfx_io')
      call getmem2d(tgbb_io,jcross1,jcross2,icross1,icross2,'tgbb_io')
      call getmem2d(uvdrag_io,jcross1,jcross2,icross1,icross2,'uvdrag_io')
      call getmem2d(q2m_io,jcross1,jcross2,icross1,icross2,'q2m_io')
      call getmem2d(u10m_io,jcross1,jcross2,icross1,icross2,'u10m_io')
      call getmem2d(v10m_io,jcross1,jcross2,icross1,icross2,'v10m_io')
      call getmem2d(w10m_io,jcross1,jcross2,icross1,icross2,'w10m_io')
      call getmem2d(br_io,jcross1,jcross2,icross1,icross2,'br_io')
      call getmem2d(ram_io,jcross1,jcross2,icross1,icross2,'ram_io')
      call getmem2d(rah_io,jcross1,jcross2,icross1,icross2,'rah_io')
      call getmem2d(ustar_io,jcross1,jcross2,icross1,icross2,'ustar_io')
      call getmem2d(zo_io,jcross1,jcross2,icross1,icross2,'zo_io')
      if ( iocncpl == 1 .or. iwavcpl == 1 ) then
        call getmem2d(dsrnof_io,jcross1,jcross2,icross1,icross2,'dsrnof_io')
        call getmem2d(dtrnof_io,jcross1,jcross2,icross1,icross2,'dtrnof_io')
      end if

#ifdef CLM45
      if ( ichem == 1 ) then
        call getmem3d(tsoi_io,jcross1,jcross2,icross1,icross2, &
                            1,num_soil_layers,'tsoi_io')
        call getmem3d(swvol_io,jcross1,jcross2,icross1,icross2, &
                            1,num_soil_layers,'swvol_io')
      end if
#else
      call getmem3d(gwet_io,1,nnsg,jcross1,jcross2,icross1,icross2,'gwet_io')
      call getmem3d(ldew_io,1,nnsg,jcross1,jcross2,icross1,icross2,'ldew_io')
      call getmem3d(taf_io,1,nnsg,jcross1,jcross2,icross1,icross2,'taf_io')
#endif
      call getmem3d(sncv_io,1,nnsg,jcross1,jcross2,icross1,icross2,'sncv_io')
      call getmem3d(sfice_io,1,nnsg,jcross1,jcross2,icross1,icross2,'sfice_io')
      call getmem3d(snag_io,1,nnsg,jcross1,jcross2,icross1,icross2,'snag_io')
      call getmem4d(sw_io,1,nnsg,jcross1,jcross2,icross1,icross2, &
                          1,num_soil_layers,'sw_io')
      call getmem3d(tgrd_io,1,nnsg,jcross1,jcross2,icross1,icross2,'tgrd_io')
      call getmem3d(tgbrd_io,1,nnsg,jcross1,jcross2,icross1,icross2,'tgbrd_io')
      call getmem3d(tlef_io,1,nnsg,jcross1,jcross2,icross1,icross2,'tlef_io')
      call getmem3d(emisv_io,1,nnsg,jcross1,jcross2,icross1,icross2,'emisv_io')
      call getmem3d(um10_io,1,nnsg,jcross1,jcross2,icross1,icross2,'um10_io')
      call getmem3d(ldmsk1_io,1,nnsg,jcross1,jcross2,icross1,icross2,'ldmsk1')
      call getmem2d(ldmsk_io,jcross1,jcross2,icross1,icross2,'ldmsk_io')
      call getmem2d(flw_io,jcross1,jcross2,icross1,icross2,'flw_io')
      call getmem2d(flwd_io,jcross1,jcross2,icross1,icross2,'flwd_io')
      call getmem2d(fsw_io,jcross1,jcross2,icross1,icross2,'fsw_io')
      call getmem2d(sabveg_io,jcross1,jcross2,icross1,icross2,'sabveg_io')
      call getmem2d(totcf_io,jcross1,jcross2,icross1,icross2,'totcf_io')
      call getmem2d(sinc_io,jcross1,jcross2,icross1,icross2,'sinc_io')
      call getmem2d(solis_io,jcross1,jcross2,icross1,icross2,'solis_io')
      call getmem2d(solvs_io,jcross1,jcross2,icross1,icross2,'solvs_io')
      call getmem2d(solvsd_io,jcross1,jcross2,icross1,icross2,'solvsd_io')
      call getmem2d(solvl_io,jcross1,jcross2,icross1,icross2,'solvl_io')
      call getmem2d(solvld_io,jcross1,jcross2,icross1,icross2,'solvld_io')
      if ( ibltyp == 2 ) then
        call getmem2d(kpbl_io,jcross1,jcross2,icross1,icross2,'kpbl_io')
      else if ( ibltyp == 4 ) then
        call getmem2d(myjsf_uz0_io,jcross1,jcross2, &
                                   icross1,icross2,'myjsf_uz0_io')
        call getmem2d(myjsf_vz0_io,jcross1,jcross2, &
                                   icross1,icross2,'myjsf_vz0_io')
        call getmem2d(myjsf_thz0_io,jcross1,jcross2, &
                                    icross1,icross2,'myjsf_thz0_io')
        call getmem2d(myjsf_qz0_io,jcross1,jcross2, &
                                   icross1,icross2,'myjsf_qz0_io')
        call getmem2d(kpbl_io,jcross1,jcross2,icross1,icross2,'kpbl_io')
      end if
      if ( idcsst == 1 ) then
        call getmem3d(sst_io,1,nnsg,jcross1,jcross2, &
                                    icross1,icross2,'sst_io')
        call getmem3d(tskin_io,1,nnsg,jcross1,jcross2, &
                                      icross1,icross2,'tskin_io')
        call getmem3d(deltas_io,1,nnsg,jcross1,jcross2, &
                                       icross1,icross2,'deltas_io')
        call getmem3d(tdeltas_io,1,nnsg,jcross1,jcross2, &
                                        icross1,icross2,'tdeltas_io')
      end if

      if ( any(icup == 4) ) then
        call getmem2d(cbmf2d_io,jcross1,jcross2,icross1,icross2,'cbmf2d_io')
      end if
      if ( ipptls > 0 ) then
        call getmem3d(fcc_io,jcross1,jcross2,icross1,icross2,1,kz,'fcc_io')
      end if

      if ( idynamic == 1 ) then
        call getmem3d(dstor_io,jdot1,jdot2,idot1,idot2,1,nsplit,'dstor_io')
        call getmem3d(hstor_io,jdot1,jdot2,idot1,idot2,1,nsplit,'hstor_io')
      end if

      if ( irrtm == 0 ) then
        call getmem4d(gasabsnxt_io,jcross1,jcross2, &
                      icross1,icross2,1,kz,1,4,'gasabsnxt_io')
        call getmem4d(gasabstot_io,jcross1,jcross2, &
                      icross1,icross2,1,kzp1,1,kzp1,'gasabstot_io')
        call getmem3d(gasemstot_io,jcross1,jcross2, &
                      icross1,icross2,1,kzp1,'gasemstot_io')
      end if

      call getmem3d(heatrt_io,jcross1,jcross2, &
                    icross1,icross2,1,kz,'heatrt_io')
      call getmem3d(o3prof_io,jcross1,jcross2, &
                    icross1,icross2,1,kzp1,'o3prof_io')

      if ( islab_ocean == 1 .and. do_restore_sst ) then
        call getmem3d(qflux_restore_sst_io,jcross1,jcross2, &
                      icross1,icross2,1,12,'qflux_restore_sst_io')
      end if

      if ( iocnflx == 2 .or. ibltyp == 3 ) then
        call getmem2d(zpbl_io,jcross1,jcross2,icross1,icross2,'zpbl_io')
      end if
      if ( any(icup == 3) ) then
        call getmem2d(cldefi_io,jcross1,jcross2,icross1,icross2,'cldefi_io')
      end if
      if ( any(icup == 6) .or. any(icup == 5) ) then
        call getmem3d(cu_avg_ww_io,jcross1,jcross2, &
                                   icross1,icross2,1,kz,'cu_avg_ww_io')
      end if
#ifdef CLM
      if ( imask == 2 ) then
        call getmem2d(lndcat_io,jcross1,jcross2,icross1,icross2,'lndcat_io')
      end if
#else
      if ( lakemod == 1 ) then
        call getmem3d(eta_io,1,nnsg,jcross1,jcross2,icross1,icross2,'eta_io')
        call getmem3d(hi_io,1,nnsg,jcross1,jcross2,icross1,icross2,'hi_io')
        call getmem4d(tlak_io,1,nnsg,jcross1,jcross2, &
                              icross1,icross2,1,ndpmax,'tlak_io')
      end if
#endif
      call getmem3d(swalb_io,1,nnsg,jcross1,jcross2,icross1,icross2,'swalb')
      call getmem3d(lwalb_io,1,nnsg,jcross1,jcross2,icross1,icross2,'lwalb')
      call getmem3d(swdiralb_io,1,nnsg,jcross1,jcross2, &
                                       icross1,icross2,'swdiralb')
      call getmem3d(swdifalb_io,1,nnsg,jcross1,jcross2, &
                                       icross1,icross2,'swdifalb')
      call getmem3d(lwdiralb_io,1,nnsg,jcross1,jcross2, &
                                       icross1,icross2,'lwdiralb')
      call getmem3d(lwdifalb_io,1,nnsg,jcross1,jcross2, &
                                       icross1,icross2,'lwdifalb')
    endif
  end subroutine allocate_mod_savefile

  subroutine read_savefile(idate)
#ifdef PNETCDF
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    type (rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: ncid
    character(256) :: ffin
    character(32) :: fbname

    if ( .not. do_parallel_save ) then
      if ( myid /= iocpu ) then
        return
      end if
    end if

    write (fbname, '(a,a)') 'SAV.', trim(tochar10(idate))
    ffin = trim(dirout)//pthsep//trim(prestr)//trim(domname)// &
           '_'//trim(fbname)//'.nc'

    call saveopen(ffin,ncid)

    if ( idynamic == 3 ) then
      call mygetvar(ncid,'atm_u',atm_u_io)
      call mygetvar(ncid,'atm_v',atm_v_io)
      call mygetvar(ncid,'atm_w',atm_w_io)
      call mygetvar(ncid,'atm_t',atm_t_io)
      call mygetvar(ncid,'atm_pai',atm_pai_io)
      call mygetvar(ncid,'atm_qx',atm_qx_io)
      if ( ibltyp == 2 ) then
        call mygetvar(ncid,'atm_tke',atm_tke_io)
      end if
      call mygetvar(ncid,'ps',ps_io)
    else
      call mygetvar(ncid,'atm1_u',atm1_u_io)
      call mygetvar(ncid,'atm2_u',atm2_u_io)
      call mygetvar(ncid,'atm1_v',atm1_v_io)
      call mygetvar(ncid,'atm2_v',atm2_v_io)
      call mygetvar(ncid,'atm1_t',atm1_t_io)
      call mygetvar(ncid,'atm2_t',atm2_t_io)
      call mygetvar(ncid,'atm1_qx',atm1_qx_io)
      call mygetvar(ncid,'atm2_qx',atm2_qx_io)
      if ( ibltyp == 2 ) then
        call mygetvar(ncid,'atm1_tke',atm1_tke_io)
        call mygetvar(ncid,'atm2_tke',atm2_tke_io)
      end if
      call mygetvar(ncid,'psa',psa_io)
      call mygetvar(ncid,'psb',psb_io)
    end if
    if ( ibltyp == 2 ) then
      call mygetvar(ncid,'kpbl',kpbl_io)
    else if ( ibltyp == 4 ) then
      call mygetvar(ncid,'tke_pbl',tke_pbl_io)
      call mygetvar(ncid,'kpbl',kpbl_io)
      call mygetvar(ncid,'myjsf_uz0',myjsf_uz0_io)
      call mygetvar(ncid,'myjsf_vz0',myjsf_vz0_io)
      call mygetvar(ncid,'myjsf_thz0',myjsf_thz0_io)
      call mygetvar(ncid,'myjsf_qz0',myjsf_qz0_io)
    end if
    if ( idynamic == 2 ) then
      call mygetvar(ncid,'atm1_w',atm1_w_io)
      call mygetvar(ncid,'atm2_w',atm2_w_io)
      call mygetvar(ncid,'atm1_pp',atm1_pp_io)
      call mygetvar(ncid,'atm2_pp',atm2_pp_io)
    end if
    call mygetvar(ncid,'hfx',hfx_io)
    call mygetvar(ncid,'qfx',qfx_io)
    call mygetvar(ncid,'tgbb',tgbb_io)
    call mygetvar(ncid,'uvdrag',uvdrag_io)
    call mygetvar(ncid,'q2m',q2m_io)
    call mygetvar(ncid,'u10m',u10m_io)
    call mygetvar(ncid,'v10m',v10m_io)
    call mygetvar(ncid,'w10m',w10m_io)
    call mygetvar(ncid,'br',br_io)
    call mygetvar(ncid,'ram',ram_io)
    call mygetvar(ncid,'rah',rah_io)
    call mygetvar(ncid,'ustar',ustar_io)
    call mygetvar(ncid,'zo',zo_io)
    if ( iocncpl == 1 .or. iwavcpl == 1 ) then
      call mygetvar(ncid,'dsrnof',dsrnof_io)
      call mygetvar(ncid,'dtrnof',dtrnof_io)
    end if
    if ( any(icup == 3) ) then
      call mygetvar(ncid,'cldefi',cldefi_io)
    end if
    if ( any(icup == 4) ) then
      call mygetvar(ncid,'cbmf2d',cbmf2d_io)
    end if
    if ( any(icup == 6) .or. any(icup == 5) ) then
      call mygetvar(ncid,'cu_avg_ww',cu_avg_ww_io)
    end if
    if ( idcsst == 1 ) then
      call mygetvar(ncid,'sst',sst_io)
      call mygetvar(ncid,'tskin',tskin_io)
      call mygetvar(ncid,'deltas',deltas_io)
      call mygetvar(ncid,'tdeltas',tdeltas_io)
    end if
    if ( irrtm == 0 ) then
      call mygetvar(ncid,'gasabsnxt',gasabsnxt_io)
      call mygetvar(ncid,'gasabstot',gasabstot_io)
      call mygetvar(ncid,'gasemstot',gasemstot_io)
    end if
    if ( ipptls > 0 ) then
      call mygetvar(ncid,'fcc',fcc_io)
    end if
    call mygetvar(ncid,'solis',solis_io)
    call mygetvar(ncid,'solvs',solvs_io)
    call mygetvar(ncid,'solvsd',solvsd_io)
    call mygetvar(ncid,'solvl',solvl_io)
    call mygetvar(ncid,'solvld',solvld_io)
    call mygetvar(ncid,'sabveg',sabveg_io)
    call mygetvar(ncid,'totcf',totcf_io,.true.)
    call mygetvar(ncid,'sw',sw_io)
    call mygetvar(ncid,'tlef',tlef_io)
    call mygetvar(ncid,'tgrd',tgrd_io)
    call mygetvar(ncid,'tgbrd',tgbrd_io)
    call mygetvar(ncid,'sncv',sncv_io)
#ifdef CLM45
    if ( ichem == 1 .and. ichecold == 0 ) then
      call mygetvar(ncid,'tsoi',tsoi_io)
      call mygetvar(ncid,'swvol',swvol_io)
    end if
#else
    call mygetvar(ncid,'gwet',gwet_io)
    call mygetvar(ncid,'ldew',ldew_io)
    call mygetvar(ncid,'taf',taf_io)
#endif
    call mygetvar(ncid,'sfice',sfice_io)
    call mygetvar(ncid,'snag',snag_io)
    call mygetvar(ncid,'ldmsk1',ldmsk1_io)
    call mygetvar(ncid,'emiss',emisv_io)
    call mygetvar(ncid,'um10',um10_io)
#ifndef CLM
    if ( lakemod == 1 ) then
      call mygetvar(ncid,'eta',eta_io)
      call mygetvar(ncid,'hi',hi_io)
      call mygetvar(ncid,'tlak',tlak_io)
    end if
#endif
    call mygetvar(ncid,'heatrt',heatrt_io)
    call mygetvar(ncid,'o3prof',o3prof_io)
    call mygetvar(ncid,'flw',flw_io)
    call mygetvar(ncid,'flwd',flwd_io)
    call mygetvar(ncid,'fsw',fsw_io)
    call mygetvar(ncid,'sinc',sinc_io)
    call mygetvar(ncid,'ldmsk',ldmsk_io)
    if ( iocnflx == 2 .or. ibltyp == 3 ) then
      call mygetvar(ncid,'zpbl',zpbl_io)
    end if
    if ( ichem == 1 .and. ichecold == 0 ) then
      if ( idynamic == 3 ) then
        call mygetvar(ncid,'trac',trac_io)
      else
        call mygetvar(ncid,'chia',chia_io)
        call mygetvar(ncid,'chib',chib_io)
      end if
      if ( igaschem == 1 .and. ichsolver > 0 ) then
        call mygetvar(ncid,'chemall',chemall_io)
        call remappnt4(taucldsp_io,tmp_io)
        call mygetvar(ncid,'taucldsp',tmp_io)
      end if
      call mygetvar(ncid,'convpr',convpr_io)
      call mygetvar(ncid,'rainout',rainout_io)
      call mygetvar(ncid,'washout',washout_io)
      call mygetvar(ncid,'remdrd',remdrd_io)
      call mygetvar(ncid,'ssw2da',ssw2da_io)
#ifdef CLM45
      call mygetvar(ncid,'duflux',duflux_io)
      call mygetvar(ncid,'voflux',voflux_io)
#else
      call mygetvar(ncid,'sdelq',sdelq_io)
      call mygetvar(ncid,'sdelt',sdelt_io)
      call mygetvar(ncid,'svegfrac2d',svegfrac2d_io)
#endif
      call mygetvar(ncid,'sfracb2d',sfracb2d_io)
      call mygetvar(ncid,'sfracs2d',sfracs2d_io)
      call mygetvar(ncid,'sfracv2d',sfracv2d_io)
    end if
    if ( idynamic == 1 ) then
      call mygetvar(ncid,'dstor',dstor_io)
      call mygetvar(ncid,'hstor',hstor_io)
    end if
    if ( islab_ocean == 1 .and. do_restore_sst ) then
      call mygetvar(ncid,'qflux_restore_sst',qflux_restore_sst_io)
      call mygetvar(ncid,'stepcount',stepcount,get_varid(ncid,'stepcount'))
    end if
#ifdef CLM
    if ( imask == 2 ) then
      call mygetvar(ncid,'lndcat',lndcat_io)
    end if
#endif
    call mygetvar(ncid,'swalb',swalb_io)
    call mygetvar(ncid,'lwalb',lwalb_io)
    call mygetvar(ncid,'swdiralb',swdiralb_io)
    call mygetvar(ncid,'swdifalb',swdifalb_io)
    call mygetvar(ncid,'lwdiralb',lwdiralb_io)
    call mygetvar(ncid,'lwdifalb',lwdifalb_io)
    if ( idynamic == 2 .and. ifupr == 1 ) then
      call mygetvar(ncid,'tmask',tmask,get_varid(ncid,'tmask'))
    end if
    call saveclose(ffin,ncid)
  end subroutine read_savefile

  subroutine write_savefile(idate)
    use netcdf
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: ncid , ivcc
    integer(ik4) , dimension(maxdims) :: dimids
    integer(ik4) , dimension(maxdims) :: wrkdim
    integer(ik4) , dimension(maxvars) :: varids
    character(256) :: ffout
    character(32) :: fbname
#ifdef CLM
    integer(ik4) :: ioff
#endif

    if ( .not. do_parallel_save ) then
      if ( myid /= iocpu ) then
#ifdef CLM
        ioff = dtsrf-dtsec
        filer_rest = restFile_filename(type='netcdf',offset=ioff)
        call restFile_write(filer_rest)
        filer_rest = restFile_filename(type='binary',offset=ioff)
        call restFile_write_binary(filer_rest)
#endif
        return
      end if
    end if

    write (fbname, '(a,a)') 'SAV.', trim(tochar10(idate))
    ffout = trim(dirout)//pthsep//trim(prestr)//trim(domname)// &
            '_'//trim(fbname)//'.nc'
    call savecreate(ffout,ncid)
    dimids(idjcross) = savedefdim(ncid,'jcross',jcross2-jcross1+1)
    dimids(idicross) = savedefdim(ncid,'icross',icross2-icross1+1)
    dimids(idjdot) = savedefdim(ncid,'jdot',jdot2-jdot1+1)
    dimids(ididot) = savedefdim(ncid,'idot',idot2-idot1+1)
    dimids(idkh) = savedefdim(ncid,'khalf',kz)
    dimids(idkf) = savedefdim(ncid,'kfull',kz+1)
    if ( idynamic == 1 ) then
      dimids(idnsplit) = savedefdim(ncid,'nsplit',nsplit)
    end if
    dimids(idnnsg) = savedefdim(ncid,'nnsg',nnsg)
    dimids(idmonth) = savedefdim(ncid,'month',12)
    dimids(idnqx) = savedefdim(ncid,'nqx',nqx)
    if ( irrtm == 0 ) then
      dimids(idspw) = savedefdim(ncid,'spw',4)
    end if
    if ( ichem == 1 ) then
      dimids(idntr) = savedefdim(ncid,'ntr',ntr)
      if ( igaschem == 1 .and. ichsolver > 0 ) then
        dimids(idspi) = savedefdim(ncid,'nspi',nspi)
        dimids(idtotsp) = savedefdim(ncid,'totsp',totsp)
      end if
#ifdef CLM45
      dimids(idndu) = savedefdim(ncid,'ndust',4)
#endif
    end if
    if ( lakemod == 1 ) then
      dimids(iddpt) = savedefdim(ncid,'dpt',ndpmax)
    end if
    if ( idynamic == 2 .and. ifupr == 1 ) then
      dimids(ikern) = savedefdim(ncid,'ikern',13)
    end if
    dimids(indep) = savedefdim(ncid,'soil_layers',num_soil_layers)

    ivcc = 0
    if ( idynamic == 3 ) then
      wrkdim(1) = dimids(idjdot)
      wrkdim(2) = dimids(idicross)
      wrkdim(3) = dimids(idkh)
      call savedefvar(ncid,'atm_u',regcm_vartype,wrkdim,1,3,varids,ivcc)
      wrkdim(1) = dimids(idjcross)
      wrkdim(2) = dimids(ididot)
      call savedefvar(ncid,'atm_v',regcm_vartype,wrkdim,1,3,varids,ivcc)
      wrkdim(1) = dimids(idjcross)
      wrkdim(2) = dimids(idicross)
      wrkdim(3) = dimids(idkf)
      call savedefvar(ncid,'atm_w',regcm_vartype,wrkdim,1,3,varids,ivcc)
      wrkdim(3) = dimids(idkh)
      call savedefvar(ncid,'atm_t',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call savedefvar(ncid,'atm_pai',regcm_vartype,wrkdim,1,3,varids,ivcc)
      wrkdim(4) = dimids(idnqx)
      call savedefvar(ncid,'atm_qx',regcm_vartype,wrkdim,1,4,varids,ivcc)
      if ( ibltyp == 2 ) then
        wrkdim(3) = dimids(idkf)
        call savedefvar(ncid,'atm_tke',regcm_vartype,wrkdim,1,3,varids,ivcc)
      end if
      call savedefvar(ncid,'ps',regcm_vartype,wrkdim,1,2,varids,ivcc)
    else
      wrkdim(1) = dimids(idjdot)
      wrkdim(2) = dimids(ididot)
      wrkdim(3) = dimids(idkh)
      call savedefvar(ncid,'atm1_u',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call savedefvar(ncid,'atm2_u',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call savedefvar(ncid,'atm1_v',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call savedefvar(ncid,'atm2_v',regcm_vartype,wrkdim,1,3,varids,ivcc)
      wrkdim(1) = dimids(idjcross)
      wrkdim(2) = dimids(idicross)
      wrkdim(3) = dimids(idkh)
      call savedefvar(ncid,'atm1_t',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call savedefvar(ncid,'atm2_t',regcm_vartype,wrkdim,1,3,varids,ivcc)
      wrkdim(4) = dimids(idnqx)
      call savedefvar(ncid,'atm1_qx',regcm_vartype,wrkdim,1,4,varids,ivcc)
      call savedefvar(ncid,'atm2_qx',regcm_vartype,wrkdim,1,4,varids,ivcc)
      if ( ibltyp == 2 ) then
        wrkdim(3) = dimids(idkf)
        call savedefvar(ncid,'atm1_tke',regcm_vartype,wrkdim,1,3,varids,ivcc)
        call savedefvar(ncid,'atm2_tke',regcm_vartype,wrkdim,1,3,varids,ivcc)
      end if
      call savedefvar(ncid,'psa',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call savedefvar(ncid,'psb',regcm_vartype,wrkdim,1,2,varids,ivcc)
    end if
    if ( ibltyp == 2 ) then
      call savedefvar(ncid,'kpbl',nf90_int,wrkdim,1,2,varids,ivcc)
    else if ( ibltyp == 4 ) then
      wrkdim(3) = dimids(idkh)
      call savedefvar(ncid,'tke_pbl',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call savedefvar(ncid,'kpbl',nf90_int,wrkdim,1,2,varids,ivcc)
      call savedefvar(ncid,'myjsf_uz0',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call savedefvar(ncid,'myjsf_vz0',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call savedefvar(ncid,'myjsf_thz0',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call savedefvar(ncid,'myjsf_qz0',regcm_vartype,wrkdim,1,2,varids,ivcc)
    end if
    if ( idynamic == 2 ) then
      wrkdim(3) = dimids(idkf)
      call savedefvar(ncid,'atm1_w',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call savedefvar(ncid,'atm2_w',regcm_vartype,wrkdim,1,3,varids,ivcc)
      wrkdim(3) = dimids(idkh)
      call savedefvar(ncid,'atm1_pp',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call savedefvar(ncid,'atm2_pp',regcm_vartype,wrkdim,1,3,varids,ivcc)
    end if
    call savedefvar(ncid,'hfx',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'qfx',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'tgbb',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'uvdrag',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'q2m',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'u10m',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'v10m',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'w10m',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'br',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'ram',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'rah',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'ustar',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'zo',regcm_vartype,wrkdim,1,2,varids,ivcc)
    wrkdim(3) = dimids(idkh)
    if ( any(icup == 3) ) then
      call savedefvar(ncid,'cldefi',regcm_vartype,wrkdim,1,2,varids,ivcc)
    end if
    if ( any(icup == 6) .or. any(icup == 5) ) then
      call savedefvar(ncid,'cu_avg_ww',regcm_vartype,wrkdim,1,3,varids,ivcc)
    end if
    if ( any(icup == 4) ) then
      call savedefvar(ncid,'cbmf2d',regcm_vartype,wrkdim,1,2,varids,ivcc)
    end if
    if ( irrtm == 0 ) then
      wrkdim(3) = dimids(idkh)
      wrkdim(4) = dimids(idspw)
      call savedefvar(ncid,'gasabsnxt',regcm_vartype,wrkdim,1,4,varids,ivcc)
      wrkdim(3) = dimids(idkf)
      wrkdim(4) = dimids(idkf)
      call savedefvar(ncid,'gasabstot',regcm_vartype,wrkdim,1,4,varids,ivcc)
      call savedefvar(ncid,'gasemstot',regcm_vartype,wrkdim,1,3,varids,ivcc)
    end if
    if ( ipptls > 0 ) then
      wrkdim(3) = dimids(idkh)
      call savedefvar(ncid,'fcc',regcm_vartype,wrkdim,1,3,varids,ivcc)
    end if
    call savedefvar(ncid,'solis',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'solvs',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'solvsd',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'solvl',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'solvld',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'sabveg',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'totcf',regcm_vartype,wrkdim,1,2,varids,ivcc)
    wrkdim(1) = dimids(idnnsg)
    wrkdim(2) = dimids(idjcross)
    wrkdim(3) = dimids(idicross)
    wrkdim(4) = dimids(indep)
    call savedefvar(ncid,'sw',regcm_vartype,wrkdim,1,4,varids,ivcc)
    call savedefvar(ncid,'tlef',regcm_vartype,wrkdim,1,3,varids,ivcc)
    call savedefvar(ncid,'tgrd',regcm_vartype,wrkdim,1,3,varids,ivcc)
    call savedefvar(ncid,'tgbrd',regcm_vartype,wrkdim,1,3,varids,ivcc)
    call savedefvar(ncid,'sncv',regcm_vartype,wrkdim,1,3,varids,ivcc)
#ifdef CLM45
    wrkdim(1) = dimids(idjcross)
    wrkdim(2) = dimids(idicross)
    wrkdim(3) = dimids(indep)
    if ( ichem == 1 ) then
      call savedefvar(ncid,'tsoi',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call savedefvar(ncid,'swvol',regcm_vartype,wrkdim,1,3,varids,ivcc)
    end if
#else
    call savedefvar(ncid,'gwet',regcm_vartype,wrkdim,1,3,varids,ivcc)
    call savedefvar(ncid,'ldew',regcm_vartype,wrkdim,1,3,varids,ivcc)
    call savedefvar(ncid,'taf',regcm_vartype,wrkdim,1,3,varids,ivcc)
#endif
    wrkdim(1) = dimids(idnnsg)
    wrkdim(2) = dimids(idjcross)
    wrkdim(3) = dimids(idicross)
    wrkdim(4) = dimids(indep)
    call savedefvar(ncid,'sfice',regcm_vartype,wrkdim,1,3,varids,ivcc)
    call savedefvar(ncid,'snag',regcm_vartype,wrkdim,1,3,varids,ivcc)
    call savedefvar(ncid,'ldmsk1',nf90_int,wrkdim,1,3,varids,ivcc)
    call savedefvar(ncid,'emiss',regcm_vartype,wrkdim,1,3,varids,ivcc)
    call savedefvar(ncid,'um10',regcm_vartype,wrkdim,1,3,varids,ivcc)
    if ( idcsst == 1 ) then
      call savedefvar(ncid,'sst',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call savedefvar(ncid,'tskin',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call savedefvar(ncid,'deltas',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call savedefvar(ncid,'tdeltas',regcm_vartype,wrkdim,1,3,varids,ivcc)
    end if
#ifndef CLM
    if ( lakemod == 1 ) then
      call savedefvar(ncid,'eta',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call savedefvar(ncid,'hi',regcm_vartype,wrkdim,1,3,varids,ivcc)
      wrkdim(4) = dimids(iddpt)
      call savedefvar(ncid,'tlak',regcm_vartype,wrkdim,1,4,varids,ivcc)
    end if
#endif
    wrkdim(1) = dimids(idjcross)
    wrkdim(2) = dimids(idicross)
    wrkdim(3) = dimids(idkh)
    call savedefvar(ncid,'heatrt',regcm_vartype,wrkdim,1,3,varids,ivcc)
    wrkdim(3) = dimids(idkf)
    call savedefvar(ncid,'o3prof',regcm_vartype,wrkdim,1,3,varids,ivcc)
    call savedefvar(ncid,'flw',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'flwd',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'fsw',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'sinc',regcm_vartype,wrkdim,1,2,varids,ivcc)
    call savedefvar(ncid,'ldmsk',nf90_int,wrkdim,1,2,varids,ivcc)
    if ( iocnflx == 2 .or. ibltyp == 3 ) then
      call savedefvar(ncid,'zpbl',regcm_vartype,wrkdim,1,2,varids,ivcc)
    end if
    if ( ichem == 1 ) then
      wrkdim(3) = dimids(idkh)
      wrkdim(4) = dimids(idntr)
      if ( idynamic == 3 ) then
        call savedefvar(ncid,'trac',regcm_vartype,wrkdim,1,4,varids,ivcc)
      else
        call savedefvar(ncid,'chia',regcm_vartype,wrkdim,1,4,varids,ivcc)
        call savedefvar(ncid,'chib',regcm_vartype,wrkdim,1,4,varids,ivcc)
      end if
      if ( igaschem == 1 .and. ichsolver > 0 ) then
        wrkdim(3) = dimids(idkh)
        wrkdim(4) = dimids(idtotsp)
        call savedefvar(ncid,'chemall',regcm_vartype,wrkdim,1,4,varids,ivcc)
        wrkdim(3) = dimids(idkf)
        wrkdim(4) = dimids(idspi)
        call savedefvar(ncid,'taucldsp',regcm_vartype,wrkdim,1,4,varids,ivcc)
      end if
      wrkdim(3) = dimids(idkh)
      wrkdim(4) = dimids(idntr)
      call savedefvar(ncid,'convpr',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call savedefvar(ncid,'rainout',regcm_vartype,wrkdim,1,4,varids,ivcc)
      call savedefvar(ncid,'washout',regcm_vartype,wrkdim,1,4,varids,ivcc)
      wrkdim(3) = dimids(idntr)
      call savedefvar(ncid,'remdrd',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call savedefvar(ncid,'ssw2da',regcm_vartype,wrkdim,1,2,varids,ivcc)
#ifdef CLM45
      wrkdim(3) = dimids(idndu)
      call savedefvar(ncid,'duflux',regcm_vartype,wrkdim,1,3,varids,ivcc)
      wrkdim(3) = dimids(idntr)
      call savedefvar(ncid,'voflux',regcm_vartype,wrkdim,1,3,varids,ivcc)
#else
      call savedefvar(ncid,'sdelq',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call savedefvar(ncid,'sdelt',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call savedefvar(ncid,'svegfrac2d',regcm_vartype,wrkdim,1,2,varids,ivcc)
#endif
      call savedefvar(ncid,'sfracb2d',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call savedefvar(ncid,'sfracs2d',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call savedefvar(ncid,'sfracv2d',regcm_vartype,wrkdim,1,2,varids,ivcc)
    end if
    if ( idynamic == 1 ) then
      wrkdim(1) = dimids(idjdot)
      wrkdim(2) = dimids(ididot)
      wrkdim(3) = dimids(idnsplit)
      call savedefvar(ncid,'dstor',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call savedefvar(ncid,'hstor',regcm_vartype,wrkdim,1,3,varids,ivcc)
    end if
    if ( islab_ocean == 1 .and. do_restore_sst ) then
      wrkdim(1) = dimids(idjcross)
      wrkdim(2) = dimids(idicross)
      wrkdim(3) = dimids(idmonth)
      call savedefvar(ncid,'qflux_restore_sst',regcm_vartype, &
                    wrkdim,1,3,varids,ivcc)
      wrkdim(1) = dimids(idmonth)
      call savedefvar(ncid,'stepcount',nf90_int,wrkdim,1,1,varids,ivcc)
    end if
#ifdef CLM
    wrkdim(1) = dimids(idjcross)
    wrkdim(2) = dimids(idicross)
    if ( imask == 2 ) then
      call savedefvar(ncid,'lndcat',regcm_vartype,wrkdim,1,2,varids,ivcc)
    end if
#endif
    wrkdim(1) = dimids(idnnsg)
    wrkdim(2) = dimids(idjcross)
    wrkdim(3) = dimids(idicross)
    call savedefvar(ncid,'swalb',regcm_vartype,wrkdim,1,3,varids,ivcc)
    call savedefvar(ncid,'lwalb',regcm_vartype,wrkdim,1,3,varids,ivcc)
    call savedefvar(ncid,'swdiralb',regcm_vartype,wrkdim,1,3,varids,ivcc)
    call savedefvar(ncid,'swdifalb',regcm_vartype,wrkdim,1,3,varids,ivcc)
    call savedefvar(ncid,'lwdiralb',regcm_vartype,wrkdim,1,3,varids,ivcc)
    call savedefvar(ncid,'lwdifalb',regcm_vartype,wrkdim,1,3,varids,ivcc)
    if ( idynamic == 2 .and. ifupr == 1 ) then
      wrkdim(1) = dimids(ikern)
      wrkdim(2) = dimids(ikern)
      call savedefvar(ncid,'tmask',regcm_vartype,wrkdim,1,2,varids,ivcc)
    end if

    call saveready(ffout,idate,ncid)

    ivcc = 0

    if ( idynamic == 3 ) then
      call myputvar(ncid,'atm_u',atm_u_io,varids,ivcc)
      call myputvar(ncid,'atm_v',atm_v_io,varids,ivcc)
      call myputvar(ncid,'atm_w',atm_w_io,varids,ivcc)
      call myputvar(ncid,'atm_t',atm_t_io,varids,ivcc)
      call myputvar(ncid,'atm_pai',atm_pai_io,varids,ivcc)
      call myputvar(ncid,'atm_qx',atm_qx_io,varids,ivcc)
      if ( ibltyp == 2 ) then
        call myputvar(ncid,'atm_tke',atm_tke_io,varids,ivcc)
      end if
      call myputvar(ncid,'ps',ps_io,varids,ivcc)
    else
      call myputvar(ncid,'atm1_u',atm1_u_io,varids,ivcc)
      call myputvar(ncid,'atm2_u',atm2_u_io,varids,ivcc)
      call myputvar(ncid,'atm1_v',atm1_v_io,varids,ivcc)
      call myputvar(ncid,'atm2_v',atm2_v_io,varids,ivcc)
      call myputvar(ncid,'atm1_t',atm1_t_io,varids,ivcc)
      call myputvar(ncid,'atm2_t',atm2_t_io,varids,ivcc)
      call myputvar(ncid,'atm1_qx',atm1_qx_io,varids,ivcc)
      call myputvar(ncid,'atm2_qx',atm2_qx_io,varids,ivcc)
      if ( ibltyp == 2 ) then
        call myputvar(ncid,'atm1_tke',atm1_tke_io,varids,ivcc)
        call myputvar(ncid,'atm2_tke',atm2_tke_io,varids,ivcc)
      end if
      call myputvar(ncid,'psa',psa_io,varids,ivcc)
      call myputvar(ncid,'psb',psb_io,varids,ivcc)
    end if
    if ( ibltyp == 2 ) then
      call myputvar(ncid,'kpbl',kpbl_io,varids,ivcc)
    else if ( ibltyp == 4 ) then
      call myputvar(ncid,'tke_pbl',tke_pbl_io,varids,ivcc)
      call myputvar(ncid,'kpbl',kpbl_io,varids,ivcc)
      call myputvar(ncid,'myjsf_uz0',myjsf_uz0_io,varids,ivcc)
      call myputvar(ncid,'myjsf_vz0',myjsf_vz0_io,varids,ivcc)
      call myputvar(ncid,'myjsf_thz0',myjsf_thz0_io,varids,ivcc)
      call myputvar(ncid,'myjsf_qz0',myjsf_qz0_io,varids,ivcc)
    end if
    if ( idynamic == 2 ) then
      call myputvar(ncid,'atm1_w',atm1_w_io,varids,ivcc)
      call myputvar(ncid,'atm2_w',atm2_w_io,varids,ivcc)
      call myputvar(ncid,'atm1_pp',atm1_pp_io,varids,ivcc)
      call myputvar(ncid,'atm2_pp',atm2_pp_io,varids,ivcc)
    end if
    call myputvar(ncid,'hfx',hfx_io,varids,ivcc)
    call myputvar(ncid,'qfx',qfx_io,varids,ivcc)
    call myputvar(ncid,'tgbb',tgbb_io,varids,ivcc)
    call myputvar(ncid,'uvdrag',uvdrag_io,varids,ivcc)
    call myputvar(ncid,'q2m',q2m_io,varids,ivcc)
    call myputvar(ncid,'u10m',u10m_io,varids,ivcc)
    call myputvar(ncid,'v10m',v10m_io,varids,ivcc)
    call myputvar(ncid,'w10m',w10m_io,varids,ivcc)
    call myputvar(ncid,'br',br_io,varids,ivcc)
    call myputvar(ncid,'ram',ram_io,varids,ivcc)
    call myputvar(ncid,'rah',rah_io,varids,ivcc)
    call myputvar(ncid,'ustar',ustar_io,varids,ivcc)
    call myputvar(ncid,'zo',zo_io,varids,ivcc)
    if ( iocncpl == 1 .or. iwavcpl == 1 ) then
      call myputvar(ncid,'dsrnof',dsrnof_io,varids,ivcc)
      call myputvar(ncid,'dtrnof',dtrnof_io,varids,ivcc)
    end if
    if ( any(icup == 3) ) then
      call myputvar(ncid,'cldefi',cldefi_io,varids,ivcc)
    end if
    if ( any(icup == 6) .or. any(icup == 5) ) then
      call myputvar(ncid,'cu_avg_ww',cu_avg_ww_io,varids,ivcc)
    end if
    if ( any(icup == 4) ) then
      call myputvar(ncid,'cbmf2d',cbmf2d_io,varids,ivcc)
    end if
    if ( irrtm == 0 ) then
      call myputvar(ncid,'gasabsnxt',gasabsnxt_io,varids,ivcc)
      call myputvar(ncid,'gasabstot',gasabstot_io,varids,ivcc)
      call myputvar(ncid,'gasemstot',gasemstot_io,varids,ivcc)
    end if
    if ( ipptls > 0 ) then
      call myputvar(ncid,'fcc',fcc_io,varids,ivcc)
    end if
    call myputvar(ncid,'solis',solis_io,varids,ivcc)
    call myputvar(ncid,'solvs',solvs_io,varids,ivcc)
    call myputvar(ncid,'solvsd',solvsd_io,varids,ivcc)
    call myputvar(ncid,'solvl',solvl_io,varids,ivcc)
    call myputvar(ncid,'solvld',solvld_io,varids,ivcc)
    call myputvar(ncid,'sabveg',sabveg_io,varids,ivcc)
    call myputvar(ncid,'totcf',totcf_io,varids,ivcc)
    call myputvar(ncid,'sw',sw_io,varids,ivcc)
    call myputvar(ncid,'tlef',tlef_io,varids,ivcc)
    call myputvar(ncid,'tgrd',tgrd_io,varids,ivcc)
    call myputvar(ncid,'tgbrd',tgbrd_io,varids,ivcc)
    call myputvar(ncid,'sncv',sncv_io,varids,ivcc)
#ifdef CLM45
    if ( ichem == 1 ) then
      call myputvar(ncid,'tsoi',tsoi_io,varids,ivcc)
      call myputvar(ncid,'swvol',swvol_io,varids,ivcc)
    end if
#else
    call myputvar(ncid,'gwet',gwet_io,varids,ivcc)
    call myputvar(ncid,'ldew',ldew_io,varids,ivcc)
    call myputvar(ncid,'taf',taf_io,varids,ivcc)
#endif
    call myputvar(ncid,'sfice',sfice_io,varids,ivcc)
    call myputvar(ncid,'snag',snag_io,varids,ivcc)
    call myputvar(ncid,'ldmsk1',ldmsk1_io,varids,ivcc)
    call myputvar(ncid,'emiss',emisv_io,varids,ivcc)
    call myputvar(ncid,'um10',um10_io,varids,ivcc)
    if ( idcsst == 1 ) then
      call myputvar(ncid,'sst',sst_io,varids,ivcc)
      call myputvar(ncid,'tskin',tskin_io,varids,ivcc)
      call myputvar(ncid,'deltas',deltas_io,varids,ivcc)
      call myputvar(ncid,'tdeltas',tdeltas_io,varids,ivcc)
    end if
#ifndef CLM
    if ( lakemod == 1 ) then
      call myputvar(ncid,'eta',eta_io,varids,ivcc)
      call myputvar(ncid,'hi',hi_io,varids,ivcc)
      call myputvar(ncid,'tlak',tlak_io,varids,ivcc)
    end if
#endif
    call myputvar(ncid,'heatrt',heatrt_io,varids,ivcc)
    call myputvar(ncid,'o3prof',o3prof_io,varids,ivcc)
    call myputvar(ncid,'flw',flw_io,varids,ivcc)
    call myputvar(ncid,'flwd',flwd_io,varids,ivcc)
    call myputvar(ncid,'fsw',fsw_io,varids,ivcc)
    call myputvar(ncid,'sinc',sinc_io,varids,ivcc)
    call myputvar(ncid,'ldmsk',ldmsk_io,varids,ivcc)
    if ( iocnflx == 2 .or. ibltyp == 3 ) then
      call myputvar(ncid,'zpbl',zpbl_io,varids,ivcc)
    end if
    if ( ichem == 1 ) then
      if ( idynamic == 3 ) then
        call myputvar(ncid,'trac',trac_io,varids,ivcc)
      else
        call myputvar(ncid,'chia',chia_io,varids,ivcc)
        call myputvar(ncid,'chib',chib_io,varids,ivcc)
      end if
      if ( igaschem == 1 .and. ichsolver > 0 ) then
        call myputvar(ncid,'chemall',chemall_io,varids,ivcc)
        call remappnt4(taucldsp_io,tmp_io)
        call myputvar(ncid,'taucldsp',tmp_io,varids,ivcc)
      end if
      call myputvar(ncid,'convpr',convpr_io,varids,ivcc)
      call myputvar(ncid,'rainout',rainout_io,varids,ivcc)
      call myputvar(ncid,'washout',washout_io,varids,ivcc)
      call myputvar(ncid,'remdrd',remdrd_io,varids,ivcc)
      call myputvar(ncid,'ssw2da',ssw2da_io,varids,ivcc)
#ifdef CLM45
      call myputvar(ncid,'duflux',duflux_io,varids,ivcc)
      call myputvar(ncid,'voflux',voflux_io,varids,ivcc)
#else
      call myputvar(ncid,'sdelq',sdelq_io,varids,ivcc)
      call myputvar(ncid,'sdelt',sdelt_io,varids,ivcc)
      call myputvar(ncid,'svegfrac2d',svegfrac2d_io,varids,ivcc)
#endif
      call myputvar(ncid,'sfracb2d',sfracb2d_io,varids,ivcc)
      call myputvar(ncid,'sfracs2d',sfracs2d_io,varids,ivcc)
      call myputvar(ncid,'sfracv2d',sfracv2d_io,varids,ivcc)
    end if
    if ( idynamic == 1 ) then
      call myputvar(ncid,'dstor',dstor_io,varids,ivcc)
      call myputvar(ncid,'hstor',hstor_io,varids,ivcc)
    end if
    if ( islab_ocean == 1 .and. do_restore_sst ) then
      call myputvar(ncid,'qflux_restore_sst',qflux_restore_sst_io,varids,ivcc)
      call myputvar(ncid,'stepcount',stepcount,size(stepcount),varids,ivcc)
    end if
#ifdef CLM
    if ( imask == 2 ) then
      call myputvar(ncid,'lndcat',lndcat_io,varids,ivcc)
    end if
#endif
    call myputvar(ncid,'swalb',swalb_io,varids,ivcc)
    call myputvar(ncid,'lwalb',lwalb_io,varids,ivcc)
    call myputvar(ncid,'swdiralb',swdiralb_io,varids,ivcc)
    call myputvar(ncid,'swdifalb',swdifalb_io,varids,ivcc)
    call myputvar(ncid,'lwdiralb',lwdiralb_io,varids,ivcc)
    call myputvar(ncid,'lwdifalb',lwdifalb_io,varids,ivcc)
    if ( idynamic == 2 .and. ifupr == 1 ) then
      call myputvar(ncid,'tmask',tmask, &
                    size(tmask,1),size(tmask,2),varids,ivcc)
    end if

    call saveclose(ffout,ncid)

#ifdef CLM
    ioff = dtsrf-dtsec
    filer_rest = restFile_filename(type='netcdf',offset=ioff)
    call restFile_write(filer_rest)
    filer_rest = restFile_filename(type='binary',offset=ioff)
    call restFile_write_binary(filer_rest)
#endif

    if ( myid == iocpu ) then
      write(stdout,*) 'SAV variables written at ', rcmtimer%str( )
    end if
  end subroutine write_savefile

  subroutine check_ok(f,l,m1)
#ifdef PNETCDF
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    character(*) , intent(in) :: f, m1
    integer(ik4) , intent(in) :: l
    if ( ncstatus /= nf90_noerr ) then
      write (stderr,*) trim(m1)
#ifdef PNETCDF
      write (stderr,*) nf90mpi_strerror(ncstatus)
#else
      write (stderr,*) nf90_strerror(ncstatus)
#endif
      call fatal(f,l,'SAVEFILE')
    end if
  end subroutine check_ok

  integer(ik4) function get_varid(ncid,vname) result(varid)
#ifdef PNETCDF
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) :: vname
#ifdef PNETCDF
    ncstatus = nf90mpi_inq_varid(ncid,vname,varid)
#else
    ncstatus = nf90_inq_varid(ncid,vname,varid)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot find variable '//trim(vname))
  end function get_varid

  subroutine savedefvar(ncid,str,ityp,idims,i1,i2,ivar,iivar)
#ifdef PNETCDF
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid , ityp , i1 , i2
    integer(ik4) , intent(inout) :: iivar
    integer(ik4) , dimension(:) , intent(in) :: idims
    integer(ik4) , dimension(:) , intent(inout) :: ivar
    character(len=*) , intent(in) :: str

    iivar = iivar + 1
#ifdef PNETCDF
    ncstatus = nf90mpi_def_var(ncid,str,ityp,idims(i1:i2),ivar(iivar))
#else
    ncstatus = nf90_def_var(ncid,str,ityp,idims(i1:i2),ivar(iivar))
#endif
    call check_ok(__FILE__,__LINE__,'Cannot create var '//trim(str))
#if defined (NETCDF4_HDF5)
#if defined (NETCDF4_COMPRESS)
    ncstatus = nf90_def_var_deflate(ncid,ivar(iivar),1,1,deflate_level)
    call check_ok(__FILE__,__LINE__, &
      'Cannot set compression level to variable '//trim(str))
#endif
#endif
  end subroutine savedefvar

  subroutine myputvar2dd(ncid,str,var,ivar,iivar)
#ifdef PNETCDF
    use mpi , only : mpi_offset_kind
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    real(rkx) , pointer , dimension(:,:) , intent(in) :: var
    integer(ik4) , dimension(:) , intent(in) :: ivar
    integer(ik4) , intent(inout) :: iivar
#ifdef PNETCDF
    integer(kind=mpi_offset_kind) , dimension(2) :: istart , icount
#else
    integer(ik4) , dimension(2) :: istart , icount
#endif
    iivar = iivar + 1
    istart(1) = lbound(var,1)
    istart(2) = lbound(var,2)
    icount(1) = ubound(var,1)-istart(1)+1
    icount(2) = ubound(var,2)-istart(2)+1
#ifdef PNETCDF
    ncstatus = nf90mpi_put_var_all(ncid,ivar(iivar),var,istart,icount)
#else
    ncstatus = nf90_put_var(ncid,ivar(iivar),var,istart,icount)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot write var '//trim(str))
  end subroutine myputvar2dd

  subroutine mygetvar2dd(ncid,str,var,skippable)
#ifdef PNETCDF
    use mpi , only : mpi_offset_kind
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: var
    logical , optional :: skippable
#ifdef PNETCDF
    integer(kind=mpi_offset_kind) , dimension(2) :: istart , icount
#else
    integer(ik4) , dimension(2) :: istart , icount
#endif
    istart(1) = lbound(var,1)
    istart(2) = lbound(var,2)
    icount(1) = ubound(var,1)-istart(1)+1
    icount(2) = ubound(var,2)-istart(2)+1
#ifdef PNETCDF
    ncstatus = nf90mpi_get_var_all(ncid,get_varid(ncid,str),var,istart,icount)
#else
    ncstatus = nf90_get_var(ncid,get_varid(ncid,str),var,istart,icount)
#endif
    if ( present(skippable) ) then
      if ( skippable ) then
        var(:,:) = 0.0_rk8
        return
      else
        call check_ok(__FILE__,__LINE__,'Cannot read var '//trim(str))
      end if
    else
      call check_ok(__FILE__,__LINE__,'Cannot read var '//trim(str))
    end if
  end subroutine mygetvar2dd

  subroutine myputvar2ddf(ncid,str,var,nx,ny,ivar,iivar)
#ifdef PNETCDF
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid , nx , ny
    character(len=*) , intent(in) :: str
    real(rkx) , dimension(nx,ny) , intent(inout) :: var
    integer(ik4) , dimension(:) , intent(in) :: ivar
    integer(ik4) , intent(inout) :: iivar
    iivar = iivar + 1
#ifdef PNETCDF
    ncstatus = nf90mpi_put_var_all(ncid,ivar(iivar),var)
#else
    ncstatus = nf90_put_var(ncid,ivar(iivar),var)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot write var '//trim(str))
  end subroutine myputvar2ddf

  subroutine mygetvar2ddf(ncid,str,var,ivar)
#ifdef PNETCDF
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    real(rkx) , dimension(:,:) , intent(inout) :: var
    integer(ik4) , intent(in) :: ivar
#ifdef PNETCDF
    ncstatus = nf90mpi_get_var_all(ncid,ivar,var)
#else
    ncstatus = nf90_get_var(ncid,ivar,var)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot read var '//trim(str))
  end subroutine mygetvar2ddf

  subroutine myputvar3dd(ncid,str,var,ivar,iivar)
#ifdef PNETCDF
    use mpi , only : mpi_offset_kind
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: var
    integer(ik4) , dimension(:) , intent(in) :: ivar
    integer(ik4) , intent(inout) :: iivar
#ifdef PNETCDF
    integer(kind=mpi_offset_kind) , dimension(3) :: istart , icount
#else
    integer(ik4) , dimension(3) :: istart , icount
#endif
    iivar = iivar + 1
    istart(1) = lbound(var,1)
    istart(2) = lbound(var,2)
    istart(3) = lbound(var,3)
    icount(1) = ubound(var,1)-istart(1)+1
    icount(2) = ubound(var,2)-istart(2)+1
    icount(3) = ubound(var,3)-istart(3)+1
#ifdef PNETCDF
    ncstatus = nf90mpi_put_var_all(ncid,ivar(iivar),var,istart,icount)
#else
    ncstatus = nf90_put_var(ncid,ivar(iivar),var,istart,icount)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot write var '//trim(str))
  end subroutine myputvar3dd

  subroutine mygetvar3dd(ncid,str,var)
#ifdef PNETCDF
    use mpi , only : mpi_offset_kind
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: var
#ifdef PNETCDF
    integer(kind=mpi_offset_kind) , dimension(3) :: istart , icount
#else
    integer(ik4) , dimension(3) :: istart , icount
#endif
    istart(1) = lbound(var,1)
    istart(2) = lbound(var,2)
    istart(3) = lbound(var,3)
    icount(1) = ubound(var,1)-istart(1)+1
    icount(2) = ubound(var,2)-istart(2)+1
    icount(3) = ubound(var,3)-istart(3)+1
#ifdef PNETCDF
    ncstatus = nf90mpi_get_var_all(ncid,get_varid(ncid,str),var,istart,icount)
#else
    ncstatus = nf90_get_var(ncid,get_varid(ncid,str),var,istart,icount)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot read var '//trim(str))
  end subroutine mygetvar3dd

  subroutine myputvar4dd(ncid,str,var,ivar,iivar)
#ifdef PNETCDF
    use mpi , only : mpi_offset_kind
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    real(rkx) , pointer , dimension(:,:,:,:) , intent(in) :: var
    integer(ik4) , dimension(:) , intent(in) :: ivar
    integer(ik4) , intent(inout) :: iivar
#ifdef PNETCDF
    integer(kind=mpi_offset_kind) , dimension(4) :: istart , icount
#else
    integer(ik4) , dimension(4) :: istart , icount
#endif
    iivar = iivar + 1
    istart(1) = lbound(var,1)
    istart(2) = lbound(var,2)
    istart(3) = lbound(var,3)
    istart(4) = lbound(var,4)
    icount(1) = ubound(var,1)-istart(1)+1
    icount(2) = ubound(var,2)-istart(2)+1
    icount(3) = ubound(var,3)-istart(3)+1
    icount(4) = ubound(var,4)-istart(4)+1
#ifdef PNETCDF
    ncstatus = nf90mpi_put_var_all(ncid,ivar(iivar),var,istart,icount)
#else
    ncstatus = nf90_put_var(ncid,ivar(iivar),var,istart,icount)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot write var '//trim(str))
  end subroutine myputvar4dd

  subroutine mygetvar4dd(ncid,str,var)
#ifdef PNETCDF
    use mpi , only : mpi_offset_kind
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    real(rkx) , pointer , dimension(:,:,:,:) , intent(inout) :: var
#ifdef PNETCDF
    integer(kind=mpi_offset_kind) , dimension(4) :: istart , icount
#else
    integer(ik4) , dimension(4) :: istart , icount
#endif
    istart(1) = lbound(var,1)
    istart(2) = lbound(var,2)
    istart(3) = lbound(var,3)
    istart(4) = lbound(var,4)
    icount(1) = ubound(var,1)-istart(1)+1
    icount(2) = ubound(var,2)-istart(2)+1
    icount(3) = ubound(var,3)-istart(3)+1
    icount(4) = ubound(var,4)-istart(4)+1
#ifdef PNETCDF
    ncstatus = nf90mpi_get_var_all(ncid,get_varid(ncid,str),var,istart,icount)
#else
    ncstatus = nf90_get_var(ncid,get_varid(ncid,str),var,istart,icount)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot read var '//trim(str))
  end subroutine mygetvar4dd

  subroutine myputvar1dif(ncid,str,var,nx,ivar,iivar)
#ifdef PNETCDF
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid , nx
    character(len=*) , intent(in) :: str
    integer(ik4) , dimension(nx) , intent(inout) :: var
    integer(ik4) , dimension(:) , intent(in) :: ivar
    integer(ik4) , intent(inout) :: iivar
    iivar = iivar + 1
#ifdef PNETCDF
    ncstatus = nf90mpi_put_var_all(ncid,ivar(iivar),var)
#else
    ncstatus = nf90_put_var(ncid,ivar(iivar),var)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot write var '//trim(str))
  end subroutine myputvar1dif

  subroutine mygetvar1dif(ncid,str,var,ivar)
#ifdef PNETCDF
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    integer(ik4) , dimension(:) , intent(inout) :: var
    integer(ik4) , intent(in) :: ivar
#ifdef PNETCDF
    ncstatus = nf90mpi_get_var_all(ncid,ivar,var)
#else
    ncstatus = nf90_get_var(ncid,ivar,var)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot read var '//trim(str))
  end subroutine mygetvar1dif

  subroutine myputvar2di(ncid,str,var,ivar,iivar)
#ifdef PNETCDF
    use mpi , only : mpi_offset_kind
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    integer(ik4) , pointer , dimension(:,:) , intent(in) :: var
    integer(ik4) , dimension(:) , intent(in) :: ivar
    integer(ik4) , intent(inout) :: iivar
#ifdef PNETCDF
    integer(kind=mpi_offset_kind) , dimension(2) :: istart , icount
#else
    integer(ik4) , dimension(2) :: istart , icount
#endif
    iivar = iivar + 1
    istart(1) = lbound(var,1)
    istart(2) = lbound(var,2)
    icount(1) = ubound(var,1)-istart(1)+1
    icount(2) = ubound(var,2)-istart(2)+1
#ifdef PNETCDF
    ncstatus = nf90mpi_put_var_all(ncid,ivar(iivar),var,istart,icount)
#else
    ncstatus = nf90_put_var(ncid,ivar(iivar),var,istart,icount)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot write var '//trim(str))
  end subroutine myputvar2di

  subroutine mygetvar2di(ncid,str,var)
#ifdef PNETCDF
    use mpi , only : mpi_offset_kind
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    integer(ik4) , pointer , dimension(:,:) , intent(inout) :: var
#ifdef PNETCDF
    integer(kind=mpi_offset_kind) , dimension(2) :: istart , icount
#else
    integer(ik4) , dimension(2) :: istart , icount
#endif
    istart(1) = lbound(var,1)
    istart(2) = lbound(var,2)
    icount(1) = ubound(var,1)-istart(1)+1
    icount(2) = ubound(var,2)-istart(2)+1
#ifdef PNETCDF
    ncstatus = nf90mpi_get_var_all(ncid,get_varid(ncid,str),var,istart,icount)
#else
    ncstatus = nf90_get_var(ncid,get_varid(ncid,str),var,istart,icount)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot read var '//trim(str))
  end subroutine mygetvar2di

  subroutine myputvar3di(ncid,str,var,ivar,iivar)
#ifdef PNETCDF
    use mpi , only : mpi_offset_kind
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    integer(ik4) , pointer , dimension(:,:,:) , intent(in) :: var
    integer(ik4) , dimension(:) , intent(in) :: ivar
    integer(ik4) , intent(inout) :: iivar
#ifdef PNETCDF
    integer(kind=mpi_offset_kind) , dimension(3) :: istart , icount
#else
    integer(ik4) , dimension(3) :: istart , icount
#endif
    iivar = iivar + 1
    istart(1) = lbound(var,1)
    istart(2) = lbound(var,2)
    istart(3) = lbound(var,3)
    icount(1) = ubound(var,1)-istart(1)+1
    icount(2) = ubound(var,2)-istart(2)+1
    icount(3) = ubound(var,3)-istart(3)+1
#ifdef PNETCDF
    ncstatus = nf90mpi_put_var_all(ncid,ivar(iivar),var,istart,icount)
#else
    ncstatus = nf90_put_var(ncid,ivar(iivar),var,istart,icount)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot write var '//trim(str))
  end subroutine myputvar3di

  subroutine mygetvar3di(ncid,str,var)
#ifdef PNETCDF
    use mpi , only : mpi_offset_kind
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    integer(ik4) , pointer , dimension(:,:,:) , intent(inout) :: var
#ifdef PNETCDF
    integer(kind=mpi_offset_kind) , dimension(3) :: istart , icount
#else
    integer(ik4) , dimension(3) :: istart , icount
#endif
    istart(1) = lbound(var,1)
    istart(2) = lbound(var,2)
    istart(3) = lbound(var,3)
    icount(1) = ubound(var,1)-istart(1)+1
    icount(2) = ubound(var,2)-istart(2)+1
    icount(3) = ubound(var,3)-istart(3)+1
#ifdef PNETCDF
    ncstatus = nf90mpi_get_var_all(ncid,get_varid(ncid,str),var,istart,icount)
#else
    ncstatus = nf90_get_var(ncid,get_varid(ncid,str),var,istart,icount)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot read var '//trim(str))
  end subroutine mygetvar3di

  subroutine saveopen(sname,ncid)
#ifdef PNETCDF
    use pnetcdf
    use mpi , only : mpi_comm_self
#else
    use netcdf
#endif
    implicit none
    character(len=*) , intent(in) :: sname
    integer(ik4) , intent(out) :: ncid
    integer(ik4) :: imode
    integer(ik4) :: int10d , ical
    type (rcm_time_and_date) :: idatex
    real(rk8) :: rtmp
    imode = nf90_nowrite
    if ( do_parallel_save ) then
#ifdef PNETCDF
      ncstatus = nf90mpi_open(get_cartcomm( ),sname,imode, &
                              ncout_mpi_info,ncid)
#else
#ifdef NETCDF4_HDF5
      imode = ior(imode,nf90_share)
      ncstatus = nf90_open(sname,imode,ncid,comm=get_cartcomm( ), &
                           info=ncout_mpi_info)
#endif
#ifdef PNETCDF_IN_NETCDF
      ncstatus = nf90_open(sname,imode,ncid=ncid,comm=get_cartcomm( ), &
                           info=ncout_mpi_info)
#endif
#endif
    else
#ifdef PNETCDF
      ncstatus = nf90mpi_open(mpi_comm_self,sname,imode, &
                              ncout_mpi_info,ncid)
#else
      ncstatus = nf90_open(sname, imode, ncid)
#endif
    end if
    call check_ok(__FILE__,__LINE__,'Cannot open savefile '//trim(sname))

#ifdef PNETCDF
    ncstatus = nf90mpi_get_att(ncid,nf90_global,'idatex',int10d)
#else
    ncstatus = nf90_get_att(ncid,nf90_global,'idatex',int10d)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot get attribute idatex')
#ifdef PNETCDF
    ncstatus = nf90mpi_get_att(ncid,nf90_global,'calendar',ical)
#else
    ncstatus = nf90_get_att(ncid,nf90_global,'calendar',ical)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot get attribute calendar')
#ifdef PNETCDF
    ncstatus = nf90mpi_get_att(ncid,nf90_global,'declin',declin)
#else
    ncstatus = nf90_get_att(ncid,nf90_global,'declin',declin)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot get attribute declin')
#ifdef PNETCDF
    ncstatus = nf90mpi_get_att(ncid,nf90_global,'solcon',solcon)
#else
    ncstatus = nf90_get_att(ncid,nf90_global,'solcon',solcon)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot get attribute solcon')
    if ( debug_level > 0 ) then
      ! Do no give any error. User may have increased debug just now.
#ifdef PNETCDF
      ncstatus = nf90mpi_get_att(ncid,nf90_global,'dryini',rtmp)
#else
      ncstatus = nf90_get_att(ncid,nf90_global,'dryini',rtmp)
#endif
      if ( ncstatus == nf90_noerr ) then
        dryini = real(rtmp,wrkp)
      else
        dryini = 0.0_wrkp
      end if
#ifdef PNETCDF
      ncstatus = nf90mpi_get_att(ncid,nf90_global,'watini',rtmp)
#else
      ncstatus = nf90_get_att(ncid,nf90_global,'watini',rtmp)
#endif
      if ( ncstatus == nf90_noerr ) then
        watini = real(rtmp,wrkp)
      else
        watini = 0.0_wrkp
      end if
#ifdef PNETCDF
      ncstatus = nf90mpi_get_att(ncid,nf90_global,'dryerror',dryerror)
#else
      ncstatus = nf90_get_att(ncid,nf90_global,'dryerror',dryerror)
#endif
      if ( ncstatus /= nf90_noerr ) dryerror = 0.0_rk8
#ifdef PNETCDF
      ncstatus = nf90mpi_get_att(ncid,nf90_global,'waterror',waterror)
#else
      ncstatus = nf90_get_att(ncid,nf90_global,'waterror',waterror)
#endif
      if ( ncstatus /= nf90_noerr ) waterror = 0.0_rk8
    end if
    idatex = i4wcal(int10d,ical)
    if ( idatex /= rcmtimer%idate ) then
      write(stderr,*) 'Mismatch in dates namelist vs SAV file'
      write(stderr,*) 'idate1 in namelist is ', rcmtimer%str( )
      write(stderr,*) 'idate  in SAV file is ', tochar(idatex)
      call fatal(__FILE__,__LINE__,'SAV FILE DO NOT MATCH IDATE1')
    end if
  end subroutine saveopen

  subroutine savecreate(sname,ncid)
#ifdef PNETCDF
    use pnetcdf
    use mpi , only : mpi_comm_self
#else
    use netcdf
#endif
    implicit none
    character(len=*) , intent(in) :: sname
    integer(ik4) , intent(out) :: ncid
    integer(ik4) :: imode
#ifndef PNETCDF
#ifdef NETCDF4_HDF5
    imode = ior(nf90_clobber, nf90_netcdf4)
#else
#ifdef NETCDF_CDF5
    imode = ior(nf90_clobber, nf90_cdf5)
#else
    imode = nf90_clobber
#endif
#endif
#else
    imode = iomode
#endif
    if ( do_parallel_save ) then
#ifdef PNETCDF
      ncstatus = nf90mpi_create(get_cartcomm( ),sname,imode, &
                                ncout_mpi_info,ncid)
#else
#ifdef NETCDF4_HDF5
      imode = ior(imode,nf90_mpiio)
      ncstatus = nf90_create(sname,imode,ncid,comm=get_cartcomm( ), &
                             info=ncout_mpi_info)
#endif
#ifdef PNETCDF_IN_NETCDF
      imode = ior(imode,nf90_mpiio)
      ncstatus = nf90_create_par(sname,imode, &
                             get_cartcomm( ),ncout_mpi_info,ncid)
#endif
#endif
    else
#ifdef PNETCDF
      ncstatus = nf90mpi_create(mpi_comm_self,sname,imode, &
                                ncout_mpi_info,ncid)
#else
      ncstatus = nf90_create(sname, imode, ncid)
#endif
    end if
    call check_ok(__FILE__,__LINE__,'Cannot create savefile '//trim(sname))
  end subroutine savecreate

  subroutine saveready(sname,idate,ncid)
#ifdef PNETCDF
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    character(len=*) , intent(in) :: sname
    type (rcm_time_and_date) , intent(in) :: idate
    integer(ik4) , intent(in) :: ncid
#ifndef NETCDF_CDF5
    integer(ik4) :: itemp
    itemp = int(toint10(idate),ik4)
#else
    integer(ik8) :: itemp
    itemp = toint10(idate)
#endif
#ifdef PNETCDF
    ncstatus = nf90mpi_put_att(ncid,nf90_global,'idatex',itemp)
#else
    ncstatus = nf90_put_att(ncid,nf90_global,'idatex',itemp)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot save idatex')
#ifdef PNETCDF
    ncstatus = nf90mpi_put_att(ncid,nf90_global,'calendar',idate%calendar)
#else
    ncstatus = nf90_put_att(ncid,nf90_global,'calendar',idate%calendar)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot save calendar')
#ifdef PNETCDF
    ncstatus = nf90mpi_put_att(ncid,nf90_global,'declin',declin)
#else
    ncstatus = nf90_put_att(ncid,nf90_global,'declin',declin)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot save declin')
#ifdef PNETCDF
    ncstatus = nf90mpi_put_att(ncid,nf90_global,'solcon',solcon)
#else
    ncstatus = nf90_put_att(ncid,nf90_global,'solcon',solcon)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot save solcon')
    if ( debug_level > 0 ) then
#ifdef PNETCDF
      ncstatus = nf90mpi_put_att(ncid,nf90_global,'dryini',real(dryini,rk8))
#else
      ncstatus = nf90_put_att(ncid,nf90_global,'dryini',real(dryini,rk8))
#endif
      call check_ok(__FILE__,__LINE__,'Cannot save dryini')
#ifdef PNETCDF
      ncstatus = nf90mpi_put_att(ncid,nf90_global,'watini',real(watini,rk8))
#else
      ncstatus = nf90_put_att(ncid,nf90_global,'watini',real(watini,rk8))
#endif
      call check_ok(__FILE__,__LINE__,'Cannot save watini')
#ifdef PNETCDF
      ncstatus = nf90mpi_put_att(ncid,nf90_global,'dryerror',dryerror)
#else
      ncstatus = nf90_put_att(ncid,nf90_global,'dryerror',dryerror)
#endif
      call check_ok(__FILE__,__LINE__,'Cannot save dryerror')
#ifdef PNETCDF
      ncstatus = nf90mpi_put_att(ncid,nf90_global,'waterror',waterror)
#else
      ncstatus = nf90_put_att(ncid,nf90_global,'waterror',waterror)
#endif
      call check_ok(__FILE__,__LINE__,'Cannot save waterror')
    end if

#ifdef PNETCDF
    ncstatus = nf90mpi_enddef(ncid)
    call check_ok(__FILE__,__LINE__,'Cannot enable savefile '//trim(sname))
#else
    ncstatus = nf90_enddef(ncid)
    call check_ok(__FILE__,__LINE__,'Cannot enable savefile '//trim(sname))
#endif
  end subroutine saveready

  integer(ik4) function savedefdim(ncid,dname,dlen) result(dimid)
#ifdef PNETCDF
    use pnetcdf
    use mpi , only : mpi_offset_kind
#else
    use netcdf
#endif
    implicit none
    integer(ik4) , intent(in) :: ncid , dlen
    character(len=*) , intent(in) :: dname
#ifdef PNETCDF
    integer(kind=mpi_offset_kind) :: xlen
    xlen = dlen
    ncstatus = nf90mpi_def_dim(ncid,dname,xlen,dimid)
#else
    ncstatus = nf90_def_dim(ncid,dname,dlen,dimid)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot create dimension '//trim(dname))
  end function savedefdim

  subroutine saveclose(sname,ncid)
#ifdef PNETCDF
    use pnetcdf
#else
    use netcdf
#endif
    implicit none
    character(len=*) , intent(in) :: sname
    integer(ik4) , intent(inout) :: ncid
#ifdef PNETCDF
    ncstatus = nf90mpi_close(ncid)
#else
    ncstatus = nf90_close(ncid)
#endif
    call check_ok(__FILE__,__LINE__,'Cannot close savefile '//trim(sname))
  end subroutine saveclose

end module mod_savefile
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
