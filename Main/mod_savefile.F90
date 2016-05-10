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
  use mod_dynparam
  use mod_runparams
  use mod_mpmessage
  use mod_mppparam
  use mod_memutil
  use mod_atm_interface , only : tmask
  use mod_lm_interface
  use mod_che_interface
  use mod_che_mppio

  implicit none

  private

  public :: allocate_mod_savefile
  public :: read_savefile
  public :: write_savefile

  integer(ik4) , parameter :: maxdims = 16
  integer(ik4) , parameter :: maxvars = 108
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

  integer(ik4) , public , pointer , dimension(:,:,:) :: ldmsk1_io
  integer(ik4) , public , pointer , dimension(:,:) :: ldmsk_io

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

  real(rkx) , public , pointer , dimension(:,:) :: psa_io
  real(rkx) , public , pointer , dimension(:,:) :: psb_io
  real(rkx) , public , pointer , dimension(:,:) :: tga_io
  real(rkx) , public , pointer , dimension(:,:) :: tgb_io
  real(rkx) , public , pointer , dimension(:,:) :: hfx_io
  real(rkx) , public , pointer , dimension(:,:) :: qfx_io
  real(rkx) , public , pointer , dimension(:,:) :: rainc_io
  real(rkx) , public , pointer , dimension(:,:) :: rainnc_io
  real(rkx) , public , pointer , dimension(:,:) :: snownc_io
  real(rkx) , public , pointer , dimension(:,:) :: tgbb_io
  real(rkx) , public , pointer , dimension(:,:) :: uvdrag_io

  real(rkx) , public , pointer , dimension(:,:,:) :: ldew_io
  real(rkx) , public , pointer , dimension(:,:,:) :: snag_io
  real(rkx) , public , pointer , dimension(:,:,:) :: sncv_io
  real(rkx) , public , pointer , dimension(:,:,:) :: sfice_io
  real(rkx) , public , pointer , dimension(:,:,:) :: gwet_io
  real(rkx) , public , pointer , dimension(:,:,:) :: rsw_io
  real(rkx) , public , pointer , dimension(:,:,:) :: ssw_io
  real(rkx) , public , pointer , dimension(:,:,:) :: tsw_io
  real(rkx) , public , pointer , dimension(:,:,:) :: taf_io
  real(rkx) , public , pointer , dimension(:,:,:) :: tgrd_io
  real(rkx) , public , pointer , dimension(:,:,:) :: tgbrd_io
  real(rkx) , public , pointer , dimension(:,:,:) :: tlef_io
  real(rkx) , public , pointer , dimension(:,:,:) :: emisv_io
  real(rkx) , public , pointer , dimension(:,:,:) :: scvk_io
  real(rkx) , public , pointer , dimension(:,:,:) :: um10_io
  real(rkx) , public , pointer , dimension(:,:,:) :: eta_io
  real(rkx) , public , pointer , dimension(:,:,:) :: hi_io
  real(rkx) , public , pointer , dimension(:,:,:,:) :: tlak_io

  real(rkx) , public , pointer , dimension(:,:) :: flw_io
  real(rkx) , public , pointer , dimension(:,:) :: flwd_io
  real(rkx) , public , pointer , dimension(:,:) :: fsw_io
  real(rkx) , public , pointer , dimension(:,:) :: sabveg_io
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
  real(rkx) , public , pointer , dimension(:,:,:) :: tbase_io

  real(rkx) , public , pointer , dimension(:,:,:) :: kfwavg_io

  real(rkx) , public , pointer , dimension(:,:,:) :: qflux_restore_sst_io

#ifdef CLM
  real(rkx) , public , pointer , dimension(:,:) :: lndcat_io
#endif

#ifdef CLM45
  real(rkx) , public , pointer , dimension(:,:,:) :: swdiralb_io
  real(rkx) , public , pointer , dimension(:,:,:) :: swdifalb_io
  real(rkx) , public , pointer , dimension(:,:,:) :: lwdiralb_io
  real(rkx) , public , pointer , dimension(:,:,:) :: lwdifalb_io
#endif

  interface myputvar
    module procedure myputvar2dd
    module procedure myputvar2ddf
    module procedure myputvar3dd
    module procedure myputvar4dd
    module procedure myputvar1dif
    module procedure myputvar2di
    module procedure myputvar3di
  end interface myputvar

  contains

  subroutine allocate_mod_savefile
    implicit none

    if ( myid == iocpu ) then
      call getmem3d(atm1_u_io,jdot1,jdot2,idot1,idot2,1,kz,'atm1_u_io')
      call getmem3d(atm1_v_io,jdot1,jdot2,idot1,idot2,1,kz,'atm1_v_io')
      call getmem3d(atm1_t_io,jcross1,jcross2,icross1,icross2,1,kz,'atm1_t_io')
      call getmem4d(atm1_qx_io,jcross1,jcross2, &
                    icross1,icross2,1,kz,1,nqx,'atm1_qx_io')
      call getmem3d(atm2_u_io,jdot1,jdot2,idot1,idot2,1,kz,'atm2_u_io')
      call getmem3d(atm2_v_io,jdot1,jdot2,idot1,idot2,1,kz,'atm2_v_io')
      call getmem3d(atm2_t_io,jcross1,jcross2,icross1,icross2,1,kz,'atm2_t_io')
      call getmem4d(atm2_qx_io,jcross1,jcross2, &
                    icross1,icross2,1,kz,1,nqx,'atm2_qx_io')
      if ( ibltyp == 2 ) then
        call getmem3d(atm1_tke_io,jcross1,jcross2, &
                      icross1,icross2,1,kzp1,'atm1_tke_io')
        call getmem3d(atm2_tke_io,jcross1,jcross2, &
                      icross1,icross2,1,kzp1,'atm2_tke_io')
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
      call getmem2d(psa_io,jcross1,jcross2,icross1,icross2,'psa_io')
      call getmem2d(psb_io,jcross1,jcross2,icross1,icross2,'psb_io')
      call getmem2d(tga_io,jcross1,jcross2,icross1,icross2,'tga_io')
      call getmem2d(tgb_io,jcross1,jcross2,icross1,icross2,'tgb_io')
      call getmem2d(hfx_io,jcross1,jcross2,icross1,icross2,'hfx_io')
      call getmem2d(qfx_io,jcross1,jcross2,icross1,icross2,'qfx_io')
      call getmem2d(rainc_io,jcross1,jcross2,icross1,icross2,'rainc_io')
      call getmem2d(rainnc_io,jcross1,jcross2,icross1,icross2,'rainnc_io')
      if ( ipptls == 2 ) then
        call getmem2d(snownc_io,jcross1,jcross2,icross1,icross2,'snownc_io')
      end if
      call getmem2d(tgbb_io,jcross1,jcross2,icross1,icross2,'tgbb_io')
      call getmem2d(uvdrag_io,jcross1,jcross2,icross1,icross2,'uvdrag_io')

      call getmem3d(ldew_io,1,nnsg,jcross1,jcross2,icross1,icross2,'ldew_io')
      call getmem3d(gwet_io,1,nnsg,jcross1,jcross2,icross1,icross2,'gwet_io')
      call getmem3d(snag_io,1,nnsg,jcross1,jcross2,icross1,icross2,'snag_io')
      call getmem3d(sncv_io,1,nnsg,jcross1,jcross2,icross1,icross2,'sncv_io')
      call getmem3d(sfice_io,1,nnsg,jcross1,jcross2,icross1,icross2,'sfice_io')
      call getmem3d(rsw_io,1,nnsg,jcross1,jcross2,icross1,icross2,'rsw_io')
      call getmem3d(ssw_io,1,nnsg,jcross1,jcross2,icross1,icross2,'ssw_io')
      call getmem3d(tsw_io,1,nnsg,jcross1,jcross2,icross1,icross2,'tsw_io')
      call getmem3d(tgrd_io,1,nnsg,jcross1,jcross2,icross1,icross2,'tgrd_io')
      call getmem3d(tgbrd_io,1,nnsg,jcross1,jcross2,icross1,icross2,'tgbrd_io')
      call getmem3d(taf_io,1,nnsg,jcross1,jcross2,icross1,icross2,'taf_io')
      call getmem3d(tlef_io,1,nnsg,jcross1,jcross2,icross1,icross2,'tlef_io')
      call getmem3d(emisv_io,1,nnsg,jcross1,jcross2,icross1,icross2,'emisv_io')
      call getmem3d(scvk_io,1,nnsg,jcross1,jcross2,icross1,icross2,'scvk_io')
      call getmem3d(um10_io,1,nnsg,jcross1,jcross2,icross1,icross2,'um10_io')
      call getmem3d(ldmsk1_io,1,nnsg,jcross1,jcross2,icross1,icross2,'ldmsk1')
      call getmem2d(ldmsk_io,jcross1,jcross2,icross1,icross2,'ldmsk_io')
      call getmem2d(flw_io,jcross1,jcross2,icross1,icross2,'flw_io')
      call getmem2d(flwd_io,jcross1,jcross2,icross1,icross2,'flwd_io')
      call getmem2d(fsw_io,jcross1,jcross2,icross1,icross2,'fsw_io')
      call getmem2d(sabveg_io,jcross1,jcross2,icross1,icross2,'sabveg_io')
      call getmem2d(sinc_io,jcross1,jcross2,icross1,icross2,'sinc_io')
      call getmem2d(solis_io,jcross1,jcross2,icross1,icross2,'solis_io')
      call getmem2d(solvs_io,jcross1,jcross2,icross1,icross2,'solvs_io')
      call getmem2d(solvsd_io,jcross1,jcross2,icross1,icross2,'solvsd_io')
      call getmem2d(solvl_io,jcross1,jcross2,icross1,icross2,'solvl_io')
      call getmem2d(solvld_io,jcross1,jcross2,icross1,icross2,'solvld_io')
      if ( ibltyp == 2 ) then
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

      if ( iocnflx == 2 ) then
        call getmem2d(zpbl_io,jcross1,jcross2,icross1,icross2,'zpbl_io')
      end if
      if ( any(icup == 3) ) then
        call getmem2d(cldefi_io,jcross1,jcross2,icross1,icross2,'cldefi_io')
        call getmem3d(tbase_io,jcross1,jcross2,icross1,icross2,1,kz,'tbase_io')
      end if
      if ( any(icup == 6) ) then
        call getmem3d(kfwavg_io,jcross1,jcross2, &
                                icross1,icross2,1,kz,'kfwavg_io')
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
#ifdef CLM45
      call getmem3d(swdiralb_io,1,nnsg,jcross1,jcross2, &
                                       icross1,icross2,'swdiralb')
      call getmem3d(swdifalb_io,1,nnsg,jcross1,jcross2, &
                                       icross1,icross2,'swdifalb')
      call getmem3d(lwdiralb_io,1,nnsg,jcross1,jcross2, &
                                       icross1,icross2,'lwdiralb')
      call getmem3d(lwdifalb_io,1,nnsg,jcross1,jcross2, &
                                       icross1,icross2,'lwdifalb')
#endif
    endif
  end subroutine allocate_mod_savefile

  subroutine read_savefile(idate)
    use netcdf
    implicit none
    type (rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: ncid
    integer(ik4) :: int10d , ical
    integer(ik8) :: idt1 , idt2
    real(rkx) :: odtsec
    character(256) :: ffin
    character(32) :: fbname

    if ( myid == iocpu ) then
      write (fbname, '(a,i10)') 'SAV.', toint10(idate)
      ffin = trim(dirout)//pthsep//trim(domname)//'_'//trim(fbname)//'.nc'
      ncstatus = nf90_open(ffin,nf90_nowrite,ncid)
      call check_ok(__FILE__,__LINE__,'Cannot open savefile '//trim(ffin))

      ncstatus = nf90_get_att(ncid,nf90_global,'ktau',ktau)
      call check_ok(__FILE__,__LINE__,'Cannot get attribute ktau')
      ncstatus = nf90_get_att(ncid,nf90_global,'dtsec',odtsec)
      call check_ok(__FILE__,__LINE__,'Cannot get attribute dtsec')
      ncstatus = nf90_get_att(ncid,nf90_global,'idatex',int10d)
      call check_ok(__FILE__,__LINE__,'Cannot get attribute idatex')
      ncstatus = nf90_get_att(ncid,nf90_global,'calendar',ical)
      call check_ok(__FILE__,__LINE__,'Cannot get attribute calendar')
      idt1 = nint(odtsec)
      idt2 = nint(dtsec)
      if ( idt1 /= idt2 ) then
        write (stdout,*) 'Recalculating ktau for the new dt'
        write (stdout,*) 'Restart file ktau is       = ', ktau
        ktau = (ktau * idt1) / idt2
        write (stdout,*) 'Actual ktau with new dt is = ', ktau
        write (stdout,*) 'Done Recalculating ktau for the new dt'
      end if
      idatex = int10d
      call setcal(idatex,ical)

      ncstatus = nf90_get_var(ncid,get_varid(ncid,'atm1_u'),atm1_u_io)
      call check_ok(__FILE__,__LINE__,'Cannot read atm1_u')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'atm2_u'),atm2_u_io)
      call check_ok(__FILE__,__LINE__,'Cannot read atm2_u')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'atm1_v'),atm1_v_io)
      call check_ok(__FILE__,__LINE__,'Cannot read atm1_v')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'atm2_v'),atm2_v_io)
      call check_ok(__FILE__,__LINE__,'Cannot read atm2_v')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'atm1_t'),atm1_t_io)
      call check_ok(__FILE__,__LINE__,'Cannot read atm1_t')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'atm2_t'),atm2_t_io)
      call check_ok(__FILE__,__LINE__,'Cannot read atm2_t')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'atm1_qx'),atm1_qx_io)
      call check_ok(__FILE__,__LINE__,'Cannot read atm1_qx')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'atm2_qx'),atm2_qx_io)
      call check_ok(__FILE__,__LINE__,'Cannot read atm2_qx')
      if ( ibltyp == 2 ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'atm1_tke'),atm1_tke_io)
        call check_ok(__FILE__,__LINE__,'Cannot read atm1_tke')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'atm2_tke'),atm2_tke_io)
        call check_ok(__FILE__,__LINE__,'Cannot read atm2_tke')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'kpbl'),kpbl_io)
        call check_ok(__FILE__,__LINE__,'Cannot read kpbl')
      end if
      if ( idynamic == 2 ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'atm1_w'),atm1_w_io)
        call check_ok(__FILE__,__LINE__,'Cannot read atm1_w')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'atm2_w'),atm2_w_io)
        call check_ok(__FILE__,__LINE__,'Cannot read atm2_w')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'atm1_pp'),atm1_pp_io)
        call check_ok(__FILE__,__LINE__,'Cannot read atm1_pp')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'atm2_pp'),atm2_pp_io)
        call check_ok(__FILE__,__LINE__,'Cannot read atm2_pp')
      end if
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'psa'),psa_io)
      call check_ok(__FILE__,__LINE__,'Cannot read psa')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'psb'),psb_io)
      call check_ok(__FILE__,__LINE__,'Cannot read psb')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'tga'),tga_io)
      call check_ok(__FILE__,__LINE__,'Cannot read tga')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'tgb'),tgb_io)
      call check_ok(__FILE__,__LINE__,'Cannot read tgb')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'hfx'),hfx_io)
      call check_ok(__FILE__,__LINE__,'Cannot read hfx')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'qfx'),qfx_io)
      call check_ok(__FILE__,__LINE__,'Cannot read qfx')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'rainc'),rainc_io)
      call check_ok(__FILE__,__LINE__,'Cannot read rainc')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'rainnc'),rainnc_io)
      call check_ok(__FILE__,__LINE__,'Cannot read rainnc')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'tgbb'),tgbb_io)
      call check_ok(__FILE__,__LINE__,'Cannot read tgbb')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'uvdrag'),uvdrag_io)
      call check_ok(__FILE__,__LINE__,'Cannot read uvdrag')
      if ( any(icup == 3) ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'cldefi'),cldefi_io)
        call check_ok(__FILE__,__LINE__,'Cannot read cldefi')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'tbase'),tbase_io)
        call check_ok(__FILE__,__LINE__,'Cannot read tbase')
      end if
      if ( any(icup == 4) ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'cbmf2d'),cbmf2d_io)
        call check_ok(__FILE__,__LINE__,'Cannot read cbmf2d')
      end if
      if ( any(icup == 6) ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'kfwavg'),kfwavg_io)
        call check_ok(__FILE__,__LINE__,'Cannot read kfwavg')
      end if
      if ( idcsst == 1 ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'sst'),sst_io)
        call check_ok(__FILE__,__LINE__,'Cannot read sst')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'tskin'),tskin_io)
        call check_ok(__FILE__,__LINE__,'Cannot read tskin')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'deltas'),deltas_io)
        call check_ok(__FILE__,__LINE__,'Cannot read deltas')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'tdeltas'),tdeltas_io)
        call check_ok(__FILE__,__LINE__,'Cannot read tdeltas')
      end if
      if ( irrtm == 0 ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'gasabsnxt'),gasabsnxt_io)
        call check_ok(__FILE__,__LINE__,'Cannot read gasabsnxt')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'gasabstot'),gasabstot_io)
        call check_ok(__FILE__,__LINE__,'Cannot read gasabstot')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'gasemstot'),gasemstot_io)
        call check_ok(__FILE__,__LINE__,'Cannot read gasemstot')
      end if
      if ( ipptls > 0 ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'fcc'),fcc_io)
        call check_ok(__FILE__,__LINE__,'Cannot read fcc')
        if ( ipptls == 2 ) then
          ncstatus = nf90_get_var(ncid,get_varid(ncid,'snownc'),snownc_io)
          call check_ok(__FILE__,__LINE__,'Cannot read snownc')
        end if
      end if
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'solis'),solis_io)
      call check_ok(__FILE__,__LINE__,'Cannot read solis')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'solvs'),solvs_io)
      call check_ok(__FILE__,__LINE__,'Cannot read solvs')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'solvsd'),solvsd_io)
      call check_ok(__FILE__,__LINE__,'Cannot read solvsd')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'solvl'),solvl_io)
      call check_ok(__FILE__,__LINE__,'Cannot read solvl')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'solvld'),solvld_io)
      call check_ok(__FILE__,__LINE__,'Cannot read solvld')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'sabveg'),sabveg_io)
      call check_ok(__FILE__,__LINE__,'Cannot read sabveg')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'tlef'),tlef_io)
      call check_ok(__FILE__,__LINE__,'Cannot read tlef')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'ssw'),ssw_io)
      call check_ok(__FILE__,__LINE__,'Cannot read ssw')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'rsw'),rsw_io)
      call check_ok(__FILE__,__LINE__,'Cannot read rsw')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'tsw'),tsw_io)
      call check_ok(__FILE__,__LINE__,'Cannot read tsw')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'tgrd'),tgrd_io)
      call check_ok(__FILE__,__LINE__,'Cannot read tgrd')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'tgbrd'),tgbrd_io)
      call check_ok(__FILE__,__LINE__,'Cannot read tgbrd')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'sncv'),sncv_io)
      call check_ok(__FILE__,__LINE__,'Cannot read sncv')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'gwet'),gwet_io)
      call check_ok(__FILE__,__LINE__,'Cannot read gwet')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'snag'),snag_io)
      call check_ok(__FILE__,__LINE__,'Cannot read snag')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'sfice'),sfice_io)
      call check_ok(__FILE__,__LINE__,'Cannot read sfice')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'ldew'),ldew_io)
      call check_ok(__FILE__,__LINE__,'Cannot read ldew')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'ldmsk1'),ldmsk1_io)
      call check_ok(__FILE__,__LINE__,'Cannot read ldmsk1')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'taf'),taf_io)
      call check_ok(__FILE__,__LINE__,'Cannot read taf')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'emiss'),emisv_io)
      call check_ok(__FILE__,__LINE__,'Cannot read emiss')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'scvk'),scvk_io)
      call check_ok(__FILE__,__LINE__,'Cannot read scvk')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'um10'),um10_io)
      call check_ok(__FILE__,__LINE__,'Cannot read um10')
#ifndef CLM
      if ( lakemod == 1 ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'eta'),eta_io)
        call check_ok(__FILE__,__LINE__,'Cannot read eta')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'hi'),hi_io)
        call check_ok(__FILE__,__LINE__,'Cannot read hi')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'tlak'),tlak_io)
        call check_ok(__FILE__,__LINE__,'Cannot read tlak')
      end if
#endif
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'heatrt'),heatrt_io)
      call check_ok(__FILE__,__LINE__,'Cannot read heatrt')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'o3prof'),o3prof_io)
      call check_ok(__FILE__,__LINE__,'Cannot read o3prof')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'flw'),flw_io)
      call check_ok(__FILE__,__LINE__,'Cannot read flw')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'flwd'),flwd_io)
      call check_ok(__FILE__,__LINE__,'Cannot read flwd')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'fsw'),fsw_io)
      call check_ok(__FILE__,__LINE__,'Cannot read fsw')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'sinc'),sinc_io)
      call check_ok(__FILE__,__LINE__,'Cannot read sinc')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'ldmsk'),ldmsk_io)
      call check_ok(__FILE__,__LINE__,'Cannot read ldmsk')
      if ( iocnflx == 2 ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'zpbl'),zpbl_io)
        call check_ok(__FILE__,__LINE__,'Cannot read zpbl')
      end if
      if ( ichem == 1 ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'chia'),chia_io)
        call check_ok(__FILE__,__LINE__,'Cannot read chia')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'chib'),chib_io)
        call check_ok(__FILE__,__LINE__,'Cannot read chib')
        if ( igaschem == 1 .and. ichsolver > 0 ) then
          ncstatus = nf90_get_var(ncid,get_varid(ncid,'chemall'),chemall_io)
          call check_ok(__FILE__,__LINE__,'Cannot read chemall')
          ncstatus = nf90_get_var(ncid,get_varid(ncid,'taucldsp'),taucldsp_io)
          call check_ok(__FILE__,__LINE__,'Cannot read taucldsp')
        end if
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'rainout'),rainout_io)
        call check_ok(__FILE__,__LINE__,'Cannot read rainout')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'washout'),washout_io)
        call check_ok(__FILE__,__LINE__,'Cannot read washout')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'remdrd'),remdrd_io)
        call check_ok(__FILE__,__LINE__,'Cannot read remdrd')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'ssw2da'),ssw2da_io)
        call check_ok(__FILE__,__LINE__,'Cannot read ssw2da')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'sdelq'),sdelq_io)
        call check_ok(__FILE__,__LINE__,'Cannot read sdelq')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'sdelt'),sdelt_io)
        call check_ok(__FILE__,__LINE__,'Cannot read sdelt')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'sfracb2d'),sfracb2d_io)
        call check_ok(__FILE__,__LINE__,'Cannot read sfracb2d')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'sfracs2d'),sfracs2d_io)
        call check_ok(__FILE__,__LINE__,'Cannot read sfracs2d')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'sfracv2d'),sfracv2d_io)
        call check_ok(__FILE__,__LINE__,'Cannot read sfracv2d')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'svegfrac2d'),svegfrac2d_io)
        call check_ok(__FILE__,__LINE__,'Cannot read svegfrac2d')
      end if
      if ( idynamic == 1 ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'dstor'),dstor_io)
        call check_ok(__FILE__,__LINE__,'Cannot read dstor')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'hstor'),hstor_io)
        call check_ok(__FILE__,__LINE__,'Cannot read hstor')
      end if
      if ( islab_ocean == 1 .and. do_restore_sst ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'qflux_restore_sst'), &
                                qflux_restore_sst_io)
        call check_ok(__FILE__,__LINE__,'Cannot read qflux_restore_sst')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'stepcount'),stepcount)
        call check_ok(__FILE__,__LINE__,'Cannot read stepcount')
      end if
#ifdef CLM
      if ( imask == 2 ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'lndcat'),lndcat_io)
        call check_ok(__FILE__,__LINE__,'Cannot read lndcat')
      end if
#endif
#ifdef CLM45
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'swdiralb'),swdiralb_io)
      call check_ok(__FILE__,__LINE__,'Cannot read swdiralb')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'swdifalb'),swdifalb_io)
      call check_ok(__FILE__,__LINE__,'Cannot read swdifalb')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'lwdiralb'),lwdiralb_io)
      call check_ok(__FILE__,__LINE__,'Cannot read lwdiralb')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'lwdifalb'),lwdifalb_io)
      call check_ok(__FILE__,__LINE__,'Cannot read lwdifalb')
#endif
      if ( idynamic == 2 .and. ifupr == 1 ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'tmask'),tmask)
        call check_ok(__FILE__,__LINE__,'Cannot read tmask')
      end if
      ncstatus = nf90_close(ncid)
      call check_ok(__FILE__,__LINE__,'Cannot close savefile '//trim(ffin))
    end if
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

    if ( myid == iocpu ) then
      write (fbname, '(a,i10)') 'SAV.', toint10(idate)
      ffout = trim(dirout)//pthsep//trim(domname)//'_'//trim(fbname)//'.nc'

      ! Use 64-bit offset format file, instead of a netCDF classic format file.
      ! Sometimes SAV files can be really big...
      ncstatus = nf90_create(ffout, iomode, ncid)
      call check_ok(__FILE__,__LINE__,'Cannot create savefile '//trim(ffout))

      ncstatus = nf90_def_dim(ncid,'jcross',jcross2-jcross1+1,dimids(idjcross))
      call check_ok(__FILE__,__LINE__,'Cannot create dimension jcross')
      ncstatus = nf90_def_dim(ncid,'icross',icross2-icross1+1,dimids(idicross))
      call check_ok(__FILE__,__LINE__,'Cannot create dimension icross')
      ncstatus = nf90_def_dim(ncid,'jdot',jdot2-jdot1+1,dimids(idjdot))
      call check_ok(__FILE__,__LINE__,'Cannot create dimension jdot')
      ncstatus = nf90_def_dim(ncid,'idot',idot2-idot1+1,dimids(ididot))
      call check_ok(__FILE__,__LINE__,'Cannot create dimension idot')
      ncstatus = nf90_def_dim(ncid,'khalf',kz,dimids(idkh))
      call check_ok(__FILE__,__LINE__,'Cannot create dimension khalf')
      ncstatus = nf90_def_dim(ncid,'kfull',kz+1,dimids(idkf))
      call check_ok(__FILE__,__LINE__,'Cannot create dimension kfull')
      if ( idynamic == 1 ) then
        ncstatus = nf90_def_dim(ncid,'nsplit',nsplit,dimids(idnsplit))
        call check_ok(__FILE__,__LINE__,'Cannot create dimension nsplit')
      end if
      ncstatus = nf90_def_dim(ncid,'nnsg',nnsg,dimids(idnnsg))
      call check_ok(__FILE__,__LINE__,'Cannot create dimension nnsg')
      ncstatus = nf90_def_dim(ncid,'month',12,dimids(idmonth))
      call check_ok(__FILE__,__LINE__,'Cannot create dimension month')
      ncstatus = nf90_def_dim(ncid,'nqx',nqx,dimids(idnqx))
      call check_ok(__FILE__,__LINE__,'Cannot create dimension nqx')
      if ( irrtm == 0 ) then
        ncstatus = nf90_def_dim(ncid,'spw',4,dimids(idspw))
        call check_ok(__FILE__,__LINE__,'Cannot create dimension spw')
      end if
      if ( ichem == 1 ) then
        ncstatus = nf90_def_dim(ncid,'ntr',ntr,dimids(idntr))
        call check_ok(__FILE__,__LINE__,'Cannot create dimension ntr')
        if ( igaschem == 1 .and. ichsolver > 0 ) then
          ncstatus = nf90_def_dim(ncid,'nspi',nspi,dimids(idspi))
          call check_ok(__FILE__,__LINE__,'Cannot create dimension nspi')
          ncstatus = nf90_def_dim(ncid,'totsp',totsp,dimids(idtotsp))
          call check_ok(__FILE__,__LINE__,'Cannot create dimension totsp')
        end if
      end if
      if ( lakemod == 1 ) then
        ncstatus = nf90_def_dim(ncid,'dpt',ndpmax,dimids(iddpt))
        call check_ok(__FILE__,__LINE__,'Cannot create dimension dpt')
      end if
      if ( idynamic == 2 .and. ifupr == 1 ) then
        ncstatus = nf90_def_dim(ncid,'ikern',13,dimids(ikern))
        call check_ok(__FILE__,__LINE__,'Cannot create dimension ikern')
      end if

      ivcc = 0
      wrkdim(1) = dimids(idjdot)
      wrkdim(2) = dimids(ididot)
      wrkdim(3) = dimids(idkh)
      call mydefvar(ncid,'atm1_u',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'atm2_u',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'atm1_v',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'atm2_v',regcm_vartype,wrkdim,1,3,varids,ivcc)
      wrkdim(1) = dimids(idjcross)
      wrkdim(2) = dimids(idicross)
      wrkdim(3) = dimids(idkh)
      call mydefvar(ncid,'atm1_t',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'atm2_t',regcm_vartype,wrkdim,1,3,varids,ivcc)
      wrkdim(4) = dimids(idnqx)
      call mydefvar(ncid,'atm1_qx',regcm_vartype,wrkdim,1,4,varids,ivcc)
      call mydefvar(ncid,'atm2_qx',regcm_vartype,wrkdim,1,4,varids,ivcc)
      if ( ibltyp == 2 ) then
        wrkdim(3) = dimids(idkf)
        call mydefvar(ncid,'atm1_tke',regcm_vartype,wrkdim,1,3,varids,ivcc)
        call mydefvar(ncid,'atm2_tke',regcm_vartype,wrkdim,1,3,varids,ivcc)
        call mydefvar(ncid,'kpbl',nf90_int,wrkdim,1,2,varids,ivcc)
      end if
      if ( idynamic == 2 ) then
        wrkdim(3) = dimids(idkf)
        call mydefvar(ncid,'atm1_w',regcm_vartype,wrkdim,1,3,varids,ivcc)
        call mydefvar(ncid,'atm2_w',regcm_vartype,wrkdim,1,3,varids,ivcc)
        wrkdim(3) = dimids(idkh)
        call mydefvar(ncid,'atm1_pp',regcm_vartype,wrkdim,1,3,varids,ivcc)
        call mydefvar(ncid,'atm2_pp',regcm_vartype,wrkdim,1,3,varids,ivcc)
      end if
      call mydefvar(ncid,'psa',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'psb',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'tga',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'tgb',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'hfx',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'qfx',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'rainc',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'rainnc',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'tgbb',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'uvdrag',regcm_vartype,wrkdim,1,2,varids,ivcc)
      wrkdim(3) = dimids(idkh)
      if ( any(icup == 3) ) then
        call mydefvar(ncid,'cldefi',regcm_vartype,wrkdim,1,2,varids,ivcc)
        call mydefvar(ncid,'tbase',regcm_vartype,wrkdim,1,3,varids,ivcc)
      end if
      if ( any(icup == 6) ) then
        call mydefvar(ncid,'kfwavg',regcm_vartype,wrkdim,1,3,varids,ivcc)
      end if
      if ( any(icup == 4) ) then
        call mydefvar(ncid,'cbmf2d',regcm_vartype,wrkdim,1,2,varids,ivcc)
      end if
      if ( irrtm == 0 ) then
        wrkdim(3) = dimids(idkh)
        wrkdim(4) = dimids(idspw)
        call mydefvar(ncid,'gasabsnxt',regcm_vartype,wrkdim,1,4,varids,ivcc)
        wrkdim(3) = dimids(idkf)
        wrkdim(4) = dimids(idkf)
        call mydefvar(ncid,'gasabstot',regcm_vartype,wrkdim,1,4,varids,ivcc)
        call mydefvar(ncid,'gasemstot',regcm_vartype,wrkdim,1,3,varids,ivcc)
      end if
      if ( ipptls > 0 ) then
        wrkdim(3) = dimids(idkh)
        call mydefvar(ncid,'fcc',regcm_vartype,wrkdim,1,3,varids,ivcc)
        if ( ipptls == 2 ) then
          call mydefvar(ncid,'snownc',regcm_vartype,wrkdim,1,2,varids,ivcc)
        end if
      end if
      call mydefvar(ncid,'solis',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'solvs',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'solvsd',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'solvl',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'solvld',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'sabveg',regcm_vartype,wrkdim,1,2,varids,ivcc)
      wrkdim(1) = dimids(idnnsg)
      wrkdim(2) = dimids(idjcross)
      wrkdim(3) = dimids(idicross)
      call mydefvar(ncid,'tlef',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'ssw',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'rsw',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'tsw',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'tgrd',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'tgbrd',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'sncv',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'gwet',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'snag',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'sfice',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'ldew',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'ldmsk1',nf90_int,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'taf',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'emiss',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'scvk',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'um10',regcm_vartype,wrkdim,1,3,varids,ivcc)
      if ( idcsst == 1 ) then
        call mydefvar(ncid,'sst',regcm_vartype,wrkdim,1,3,varids,ivcc)
        call mydefvar(ncid,'tskin',regcm_vartype,wrkdim,1,3,varids,ivcc)
        call mydefvar(ncid,'deltas',regcm_vartype,wrkdim,1,3,varids,ivcc)
        call mydefvar(ncid,'tdeltas',regcm_vartype,wrkdim,1,3,varids,ivcc)
      end if
#ifndef CLM
      if ( lakemod == 1 ) then
        call mydefvar(ncid,'eta',regcm_vartype,wrkdim,1,3,varids,ivcc)
        call mydefvar(ncid,'hi',regcm_vartype,wrkdim,1,3,varids,ivcc)
        wrkdim(4) = dimids(iddpt)
        call mydefvar(ncid,'tlak',regcm_vartype,wrkdim,1,4,varids,ivcc)
      end if
#endif
      wrkdim(1) = dimids(idjcross)
      wrkdim(2) = dimids(idicross)
      wrkdim(3) = dimids(idkh)
      call mydefvar(ncid,'heatrt',regcm_vartype,wrkdim,1,3,varids,ivcc)
      wrkdim(3) = dimids(idkf)
      call mydefvar(ncid,'o3prof',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'flw',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'flwd',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'fsw',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'sinc',regcm_vartype,wrkdim,1,2,varids,ivcc)
      call mydefvar(ncid,'ldmsk',nf90_int,wrkdim,1,2,varids,ivcc)
      if ( iocnflx == 2 ) then
        call mydefvar(ncid,'zpbl',regcm_vartype,wrkdim,1,2,varids,ivcc)
      end if
      if ( ichem == 1 ) then
        wrkdim(3) = dimids(idkh)
        wrkdim(4) = dimids(idntr)
        call mydefvar(ncid,'chia',regcm_vartype,wrkdim,1,4,varids,ivcc)
        call mydefvar(ncid,'chib',regcm_vartype,wrkdim,1,4,varids,ivcc)
        if ( igaschem == 1 .and. ichsolver > 0 ) then
          wrkdim(3) = dimids(idkh)
          wrkdim(4) = dimids(idtotsp)
          call mydefvar(ncid,'chemall',regcm_vartype,wrkdim,1,4,varids,ivcc)
          wrkdim(3) = dimids(idkf)
          wrkdim(4) = dimids(idspi)
          call mydefvar(ncid,'taucldsp',regcm_vartype,wrkdim,1,4,varids,ivcc)
        end if
        wrkdim(3) = dimids(idkh)
        wrkdim(4) = dimids(idntr)
        call mydefvar(ncid,'rainout',regcm_vartype,wrkdim,1,4,varids,ivcc)
        call mydefvar(ncid,'washout',regcm_vartype,wrkdim,1,4,varids,ivcc)
        wrkdim(3) = dimids(idntr)
        call mydefvar(ncid,'remdrd',regcm_vartype,wrkdim,1,3,varids,ivcc)
        call mydefvar(ncid,'ssw2da',regcm_vartype,wrkdim,1,2,varids,ivcc)
        call mydefvar(ncid,'sdelq',regcm_vartype,wrkdim,1,2,varids,ivcc)
        call mydefvar(ncid,'sdelt',regcm_vartype,wrkdim,1,2,varids,ivcc)
        call mydefvar(ncid,'sfracb2d',regcm_vartype,wrkdim,1,2,varids,ivcc)
        call mydefvar(ncid,'sfracs2d',regcm_vartype,wrkdim,1,2,varids,ivcc)
        call mydefvar(ncid,'sfracv2d',regcm_vartype,wrkdim,1,2,varids,ivcc)
        call mydefvar(ncid,'svegfrac2d',regcm_vartype,wrkdim,1,2,varids,ivcc)
      end if
      if ( idynamic == 1 ) then
        wrkdim(1) = dimids(idjdot)
        wrkdim(2) = dimids(ididot)
        wrkdim(3) = dimids(idnsplit)
        call mydefvar(ncid,'dstor',regcm_vartype,wrkdim,1,3,varids,ivcc)
        call mydefvar(ncid,'hstor',regcm_vartype,wrkdim,1,3,varids,ivcc)
      end if
      if ( islab_ocean == 1 .and. do_restore_sst ) then
        wrkdim(1) = dimids(idjcross)
        wrkdim(2) = dimids(idicross)
        call mydefvar(ncid,'qflux_restore_sst',regcm_vartype, &
                      wrkdim,1,2,varids,ivcc)
        wrkdim(1) = dimids(idmonth)
        call mydefvar(ncid,'stepcount',nf90_int,wrkdim,1,1,varids,ivcc)
      end if
#ifdef CLM
      wrkdim(1) = dimids(idjcross)
      wrkdim(2) = dimids(idicross)
      if ( imask == 2 ) then
        call mydefvar(ncid,'lndcat',regcm_vartype,wrkdim,1,2,varids,ivcc)
      end if
#endif
#ifdef CLM45
      wrkdim(1) = dimids(idnnsg)
      wrkdim(2) = dimids(idjcross)
      wrkdim(3) = dimids(idicross)
      call mydefvar(ncid,'swdiralb',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'swdifalb',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'lwdiralb',regcm_vartype,wrkdim,1,3,varids,ivcc)
      call mydefvar(ncid,'lwdifalb',regcm_vartype,wrkdim,1,3,varids,ivcc)
#endif
      if ( idynamic == 2 .and. ifupr == 1 ) then
        wrkdim(1) = dimids(ikern)
        wrkdim(2) = dimids(ikern)
        call mydefvar(ncid,'tmask',regcm_vartype,wrkdim,1,2,varids,ivcc)
      end if

      ncstatus = nf90_put_att(ncid,nf90_global,'ktau',ktau)
      call check_ok(__FILE__,__LINE__,'Cannot save ktau')
      ncstatus = nf90_put_att(ncid,nf90_global,'dtsec',dtsec)
      call check_ok(__FILE__,__LINE__,'Cannot save dtsec')
      ncstatus = nf90_put_att(ncid,nf90_global,'idatex',toint10(idatex))
      call check_ok(__FILE__,__LINE__,'Cannot save idatex')
      ncstatus = nf90_put_att(ncid,nf90_global,'calendar',idatex%calendar)
      call check_ok(__FILE__,__LINE__,'Cannot save calendar')

      ncstatus = nf90_enddef(ncid)
      call check_ok(__FILE__,__LINE__,'Cannot setup savefile '//trim(ffout))

      ivcc = 0
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
        call myputvar(ncid,'atm2_tke',atm1_tke_io,varids,ivcc)
        call myputvar(ncid,'kpbl',kpbl_io,varids,ivcc)
      end if
      if ( idynamic == 2 ) then
        call myputvar(ncid,'atm1_w',atm1_w_io,varids,ivcc)
        call myputvar(ncid,'atm2_w',atm2_w_io,varids,ivcc)
        call myputvar(ncid,'atm1_pp',atm1_pp_io,varids,ivcc)
        call myputvar(ncid,'atm2_pp',atm2_pp_io,varids,ivcc)
      end if
      call myputvar(ncid,'psa',psa_io,varids,ivcc)
      call myputvar(ncid,'psb',psb_io,varids,ivcc)
      call myputvar(ncid,'tga',tga_io,varids,ivcc)
      call myputvar(ncid,'tgb',tgb_io,varids,ivcc)
      call myputvar(ncid,'hfx',hfx_io,varids,ivcc)
      call myputvar(ncid,'qfx',qfx_io,varids,ivcc)
      call myputvar(ncid,'rainc',rainc_io,varids,ivcc)
      call myputvar(ncid,'rainnc',rainnc_io,varids,ivcc)
      call myputvar(ncid,'tgbb',tgbb_io,varids,ivcc)
      call myputvar(ncid,'uvdrag',uvdrag_io,varids,ivcc)
      if ( any(icup == 3) ) then
        call myputvar(ncid,'cldefi',cldefi_io,varids,ivcc)
        call myputvar(ncid,'tbase',tbase_io,varids,ivcc)
      end if
      if ( any(icup == 6) ) then
        call myputvar(ncid,'kfwavg',kfwavg_io,varids,ivcc)
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
        if ( ipptls == 2 ) then
          call myputvar(ncid,'snownc',snownc_io,varids,ivcc)
        end if
      end if
      call myputvar(ncid,'solis',solis_io,varids,ivcc)
      call myputvar(ncid,'solvs',solvs_io,varids,ivcc)
      call myputvar(ncid,'solvsd',solvsd_io,varids,ivcc)
      call myputvar(ncid,'solvl',solvl_io,varids,ivcc)
      call myputvar(ncid,'solvld',solvld_io,varids,ivcc)
      call myputvar(ncid,'sabveg',sabveg_io,varids,ivcc)
      call myputvar(ncid,'tlef',tlef_io,varids,ivcc)
      call myputvar(ncid,'ssw',ssw_io,varids,ivcc)
      call myputvar(ncid,'rsw',rsw_io,varids,ivcc)
      call myputvar(ncid,'tsw',tsw_io,varids,ivcc)
      call myputvar(ncid,'tgrd',tgrd_io,varids,ivcc)
      call myputvar(ncid,'tgbrd',tgbrd_io,varids,ivcc)
      call myputvar(ncid,'sncv',sncv_io,varids,ivcc)
      call myputvar(ncid,'gwet',gwet_io,varids,ivcc)
      call myputvar(ncid,'snag',snag_io,varids,ivcc)
      call myputvar(ncid,'sfice',sfice_io,varids,ivcc)
      call myputvar(ncid,'ldew',ldew_io,varids,ivcc)
      call myputvar(ncid,'ldmsk1',ldmsk1_io,varids,ivcc)
      call myputvar(ncid,'taf',taf_io,varids,ivcc)
      call myputvar(ncid,'emiss',emisv_io,varids,ivcc)
      call myputvar(ncid,'scvk',scvk_io,varids,ivcc)
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
      if ( iocnflx == 2 ) then
        call myputvar(ncid,'zpbl',zpbl_io,varids,ivcc)
      end if
      if ( ichem == 1 ) then
        call myputvar(ncid,'chia',chia_io,varids,ivcc)
        call myputvar(ncid,'chib',chib_io,varids,ivcc)
        if ( igaschem == 1 .and. ichsolver > 0 ) then
          call myputvar(ncid,'chemall',chemall_io,varids,ivcc)
          call myputvar(ncid,'taucldsp',taucldsp_io,varids,ivcc)
        end if
        call myputvar(ncid,'rainout',rainout_io,varids,ivcc)
        call myputvar(ncid,'washout',washout_io,varids,ivcc)
        call myputvar(ncid,'remdrd',remdrd_io,varids,ivcc)
        call myputvar(ncid,'ssw2da',ssw2da_io,varids,ivcc)
        call myputvar(ncid,'sdelq',sdelq_io,varids,ivcc)
        call myputvar(ncid,'sdelt',sdelt_io,varids,ivcc)
        call myputvar(ncid,'sfracb2d',sfracb2d_io,varids,ivcc)
        call myputvar(ncid,'sfracs2d',sfracs2d_io,varids,ivcc)
        call myputvar(ncid,'sfracv2d',sfracv2d_io,varids,ivcc)
        call myputvar(ncid,'svegfrac2d',svegfrac2d_io,varids,ivcc)
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
#ifdef CLM45
      call myputvar(ncid,'swdiralb',swdiralb_io,varids,ivcc)
      call myputvar(ncid,'swdifalb',swdifalb_io,varids,ivcc)
      call myputvar(ncid,'lwdiralb',lwdiralb_io,varids,ivcc)
      call myputvar(ncid,'lwdifalb',lwdifalb_io,varids,ivcc)
#endif
      if ( idynamic == 2 .and. ifupr == 1 ) then
        call myputvar(ncid,'tmask',tmask, &
                      size(tmask,1),size(tmask,2),varids,ivcc)
      end if

      ncstatus = nf90_close(ncid)
      call check_ok(__FILE__,__LINE__,'Cannot close savefile '//trim(ffout))
    end if

#ifdef CLM
    ioff = int(dtsec*dble(ntsrf-1))
    filer_rest = restFile_filename(type='netcdf',offset=ioff)
    call restFile_write(filer_rest)
    filer_rest = restFile_filename(type='binary',offset=ioff)
    call restFile_write_binary(filer_rest)
#endif

    if ( myid == iocpu ) then
      write(stdout,*) 'SAV variables written at ', tochar(idate)
    end if
  end subroutine write_savefile

  subroutine check_ok(f,l,m1)
    use netcdf
    implicit none
    character(*) , intent(in) :: f, m1
    integer(ik4) , intent(in) :: l
    if ( ncstatus /= nf90_noerr ) then
      write (stderr,*) trim(m1)
      write (stderr,*) nf90_strerror(ncstatus)
      call fatal(f,l,'SAVEFILE')
    end if
  end subroutine check_ok

  integer(ik4) function get_varid(ncid,vname) result(varid)
    use netcdf
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) :: vname
    ncstatus = nf90_inq_varid(ncid,vname,varid)
    call check_ok(__FILE__,__LINE__,'Cannot find variable '//trim(vname))
  end function get_varid

  subroutine mydefvar(ncid,str,ityp,idims,i1,i2,ivar,iivar)
    use netcdf
    implicit none
    integer(ik4) , intent(in) :: ncid , ityp , i1 , i2
    integer(ik4) , intent(inout) :: iivar
    integer(ik4) , dimension(:) , intent(in) :: idims
    integer(ik4) , dimension(:) , intent(inout) :: ivar
    character(len=*) , intent(in) :: str
    iivar = iivar + 1
    ncstatus = nf90_def_var(ncid,str,ityp,idims(i1:i2),ivar(iivar))
    call check_ok(__FILE__,__LINE__,'Cannot create var '//trim(str))
  end subroutine mydefvar

  subroutine myputvar2dd(ncid,str,var,ivar,iivar)
    use netcdf
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    real(rkx) , pointer , dimension(:,:) , intent(in) :: var
    integer(ik4) , dimension(:) , intent(in) :: ivar
    integer(ik4) , intent(inout) :: iivar
    iivar = iivar + 1
    ncstatus = nf90_put_var(ncid,ivar(iivar),var)
    call check_ok(__FILE__,__LINE__,'Cannot write var '//trim(str))
  end subroutine myputvar2dd

  subroutine myputvar2ddf(ncid,str,var,nx,ny,ivar,iivar)
    use netcdf
    implicit none
    integer(ik4) , intent(in) :: ncid , nx , ny
    character(len=*) , intent(in) :: str
    real(rkx) , dimension(nx,ny) , intent(in) :: var
    integer(ik4) , dimension(:) , intent(in) :: ivar
    integer(ik4) , intent(inout) :: iivar
    iivar = iivar + 1
    ncstatus = nf90_put_var(ncid,ivar(iivar),var)
    call check_ok(__FILE__,__LINE__,'Cannot write var '//trim(str))
  end subroutine myputvar2ddf

  subroutine myputvar3dd(ncid,str,var,ivar,iivar)
    use netcdf
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: var
    integer(ik4) , dimension(:) , intent(in) :: ivar
    integer(ik4) , intent(inout) :: iivar
    iivar = iivar + 1
    ncstatus = nf90_put_var(ncid,ivar(iivar),var)
    call check_ok(__FILE__,__LINE__,'Cannot write var '//trim(str))
  end subroutine myputvar3dd

  subroutine myputvar4dd(ncid,str,var,ivar,iivar)
    use netcdf
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    real(rkx) , pointer , dimension(:,:,:,:) , intent(in) :: var
    integer(ik4) , dimension(:) , intent(in) :: ivar
    integer(ik4) , intent(inout) :: iivar
    iivar = iivar + 1
    ncstatus = nf90_put_var(ncid,ivar(iivar),var)
    call check_ok(__FILE__,__LINE__,'Cannot write var '//trim(str))
  end subroutine myputvar4dd

  subroutine myputvar1dif(ncid,str,var,nx,ivar,iivar)
    use netcdf
    implicit none
    integer(ik4) , intent(in) :: ncid , nx
    character(len=*) , intent(in) :: str
    integer(ik4) , dimension(nx) , intent(in) :: var
    integer(ik4) , dimension(:) , intent(in) :: ivar
    integer(ik4) , intent(inout) :: iivar
    iivar = iivar + 1
    ncstatus = nf90_put_var(ncid,ivar(iivar),var)
    call check_ok(__FILE__,__LINE__,'Cannot write var '//trim(str))
  end subroutine myputvar1dif

  subroutine myputvar2di(ncid,str,var,ivar,iivar)
    use netcdf
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    integer(ik4) , pointer , dimension(:,:) , intent(in) :: var
    integer(ik4) , dimension(:) , intent(in) :: ivar
    integer(ik4) , intent(inout) :: iivar
    iivar = iivar + 1
    ncstatus = nf90_put_var(ncid,ivar(iivar),var)
    call check_ok(__FILE__,__LINE__,'Cannot write var '//trim(str))
  end subroutine myputvar2di

  subroutine myputvar3di(ncid,str,var,ivar,iivar)
    use netcdf
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: str
    integer(ik4) , pointer , dimension(:,:,:) , intent(in) :: var
    integer(ik4) , dimension(:) , intent(in) :: ivar
    integer(ik4) , intent(inout) :: iivar
    iivar = iivar + 1
    ncstatus = nf90_put_var(ncid,ivar(iivar),var)
    call check_ok(__FILE__,__LINE__,'Cannot write var '//trim(str))
  end subroutine myputvar3di

end module mod_savefile
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
