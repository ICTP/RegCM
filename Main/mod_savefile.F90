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

  use mod_runparams
  use mod_mpmessage
  use mod_mppparam
  use mod_memutil
  use mod_lm_interface
  use mod_che_mppio

  private

  public :: allocate_mod_savefile
  public :: read_savefile
  public :: write_savefile

  integer(ik4) , parameter :: maxdims = 16
  integer(ik4) , parameter :: maxvars = 96
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

  integer(ik4) :: isavlast

#ifdef CLM
  character(len=256) :: thisclmrest
  character(len=256) :: lastclmrest
#endif

  data isavlast /-1/

  integer(ik4) , public , pointer , dimension(:,:,:) :: ldmsk1_io
  integer(ik4) , public , pointer , dimension(:,:) :: ldmsk_io

  real(rk8) , public , pointer , dimension(:,:,:) :: atm1_u_io
  real(rk8) , public , pointer , dimension(:,:,:) :: atm2_u_io
  real(rk8) , public , pointer , dimension(:,:,:) :: atm1_v_io
  real(rk8) , public , pointer , dimension(:,:,:) :: atm2_v_io
  real(rk8) , public , pointer , dimension(:,:,:) :: atm1_t_io
  real(rk8) , public , pointer , dimension(:,:,:) :: atm2_t_io
  real(rk8) , public , pointer , dimension(:,:,:,:) :: atm1_qx_io
  real(rk8) , public , pointer , dimension(:,:,:,:) :: atm2_qx_io
  real(rk8) , public , pointer , dimension(:,:,:) :: atm1_tke_io
  real(rk8) , public , pointer , dimension(:,:,:) :: atm2_tke_io

  real(rk8) , public , pointer , dimension(:,:) :: psa_io
  real(rk8) , public , pointer , dimension(:,:) :: psb_io
  real(rk8) , public , pointer , dimension(:,:) :: tga_io
  real(rk8) , public , pointer , dimension(:,:) :: tgb_io
  real(rk8) , public , pointer , dimension(:,:) :: hfx_io
  real(rk8) , public , pointer , dimension(:,:) :: qfx_io
  real(rk8) , public , pointer , dimension(:,:) :: rainc_io
  real(rk8) , public , pointer , dimension(:,:) :: rainnc_io
  real(rk8) , public , pointer , dimension(:,:) :: snownc_io
  real(rk8) , public , pointer , dimension(:,:) :: tgbb_io
  real(rk8) , public , pointer , dimension(:,:) :: uvdrag_io

  real(rk8) , public , pointer , dimension(:,:,:) :: ldew_io
  real(rk8) , public , pointer , dimension(:,:,:) :: snag_io
  real(rk8) , public , pointer , dimension(:,:,:) :: sncv_io
  real(rk8) , public , pointer , dimension(:,:,:) :: sfice_io
  real(rk8) , public , pointer , dimension(:,:,:) :: gwet_io
  real(rk8) , public , pointer , dimension(:,:,:) :: rsw_io
  real(rk8) , public , pointer , dimension(:,:,:) :: ssw_io
  real(rk8) , public , pointer , dimension(:,:,:) :: tsw_io
  real(rk8) , public , pointer , dimension(:,:,:) :: taf_io
  real(rk8) , public , pointer , dimension(:,:,:) :: tgrd_io
  real(rk8) , public , pointer , dimension(:,:,:) :: tgbrd_io
  real(rk8) , public , pointer , dimension(:,:,:) :: tlef_io
  real(rk8) , public , pointer , dimension(:,:,:) :: emisv_io
  real(rk8) , public , pointer , dimension(:,:,:) :: scvk_io
  real(rk8) , public , pointer , dimension(:,:,:) :: eta_io
  real(rk8) , public , pointer , dimension(:,:,:) :: hi_io
  real(rk8) , public , pointer , dimension(:,:,:,:) :: tlak_io

  real(rk8) , public , pointer , dimension(:,:) :: flw_io
  real(rk8) , public , pointer , dimension(:,:) :: flwd_io
  real(rk8) , public , pointer , dimension(:,:) :: fsw_io
  real(rk8) , public , pointer , dimension(:,:) :: sabveg_io
  real(rk8) , public , pointer , dimension(:,:) :: sinc_io
  real(rk8) , public , pointer , dimension(:,:) :: solis_io
  real(rk8) , public , pointer , dimension(:,:) :: solvs_io
  real(rk8) , public , pointer , dimension(:,:) :: solvsd_io
  real(rk8) , public , pointer , dimension(:,:) :: solvl_io
  real(rk8) , public , pointer , dimension(:,:) :: solvld_io

  real(rk8) , public , pointer , dimension(:,:,:) :: sst_io
  real(rk8) , public , pointer , dimension(:,:,:) :: tskin_io
  real(rk8) , public , pointer , dimension(:,:,:) :: tdeltas_io
  real(rk8) , public , pointer , dimension(:,:,:) :: deltas_io

  integer(ik4) , public , pointer , dimension(:,:) :: kpbl_io
  real(rk8) , public , pointer , dimension(:,:) :: zpbl_io

  real(rk8) , public , pointer , dimension(:,:) :: cbmf2d_io

  real(rk8) , public , pointer , dimension(:,:,:) :: fcc_io

  real(rk8) , public , pointer , dimension(:,:,:) :: rsheat_io
  real(rk8) , public , pointer , dimension(:,:,:) :: rswat_io

  real(rk8) , public , pointer , dimension(:,:,:,:) :: gasabsnxt_io
  real(rk8) , public , pointer , dimension(:,:,:,:) :: gasabstot_io
  real(rk8) , public , pointer , dimension(:,:,:) :: gasemstot_io

  real(rk8) , public , pointer , dimension(:,:,:) :: heatrt_io
  real(rk8) , public , pointer , dimension(:,:,:) :: o3prof_io

  real(rk8) , public , pointer , dimension(:,:,:) :: dstor_io
  real(rk8) , public , pointer , dimension(:,:,:) :: hstor_io

  real(rk8) , public , pointer , dimension(:,:) :: cldefi_io
  real(rk8) , public , pointer , dimension(:,:,:) :: tbase_io

  real(rk8) , public , pointer , dimension(:,:,:) :: qflux_restore_sst_io

#ifdef CLM
  real(rk8) , public , pointer , dimension(:,:) :: lndcat_io
#endif

  contains

  subroutine allocate_mod_savefile
    implicit none

    if ( myid == iocpu ) then
      call getmem3d(atm1_u_io,jdot1,jdot2,idot1,idot2,1,kz,'atm1_u_io')
      call getmem3d(atm1_v_io,jdot1,jdot2,idot1,idot2,1,kz,'atm1_v_io')
      call getmem3d(atm1_t_io,jcross1,jcross2,icross1,icross2,1,kz,'atm1_t_io')
      call getmem4d(atm1_qx_io,jcross1,jcross2, &
                    icross1,icross2,1,kz,1,nqx,'atm1_qx_io')
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call getmem3d(atm1_tke_io,jcross1,jcross2, &
                      icross1,icross2,1,kzp1,'atm1_tke_io')
      end if
      call getmem3d(atm2_u_io,jdot1,jdot2,idot1,idot2,1,kz,'atm2_u_io')
      call getmem3d(atm2_v_io,jdot1,jdot2,idot1,idot2,1,kz,'atm2_v_io')
      call getmem3d(atm2_t_io,jcross1,jcross2,icross1,icross2,1,kz,'atm2_t_io')
      call getmem4d(atm2_qx_io,jcross1,jcross2, &
                    icross1,icross2,1,kz,1,nqx,'atm2_qx_io')
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call getmem3d(atm2_tke_io,jcross1,jcross2, &
                      icross1,icross2,1,kzp1,'atm2_tke_io')
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
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
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
      if ( any(icup == 1) ) then
        call getmem3d(rsheat_io,jcross1,jcross2, &
                                icross1,icross2,1,kz,'rsheat_io')
        call getmem3d(rswat_io,jcross1,jcross2,icross1,icross2,1,kz,'rswat_io')
      end if

      call getmem3d(dstor_io,jdot1,jdot2,idot1,idot2,1,nsplit,'dstor_io')
      call getmem3d(hstor_io,jdot1,jdot2,idot1,idot2,1,nsplit,'hstor_io')

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
    endif
  end subroutine allocate_mod_savefile

  subroutine read_savefile(idate)
    use netcdf
    implicit none
    type (rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: ncid
    integer(ik4) :: int10d , ical
    integer(ik8) :: idt1 , idt2
    real(rk8) :: odtsec
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
      idt1 = idnint(odtsec)
      idt2 = idnint(dtsec)
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
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'atm1_tke'),atm1_tke_io)
        call check_ok(__FILE__,__LINE__,'Cannot read atm1_tke')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'atm2_tke'),atm2_tke_io)
        call check_ok(__FILE__,__LINE__,'Cannot read atm2_tke')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'kpbl'),kpbl_io)
        call check_ok(__FILE__,__LINE__,'Cannot read kpbl')
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
      if ( any(icup == 1) ) then
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'rsheat'),rsheat_io)
        call check_ok(__FILE__,__LINE__,'Cannot read rsheat')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'rswat'),rswat_io)
        call check_ok(__FILE__,__LINE__,'Cannot read rswat')
      end if
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
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'remlsc'),remlsc_io)
        call check_ok(__FILE__,__LINE__,'Cannot read remlsc')
        ncstatus = nf90_get_var(ncid,get_varid(ncid,'remcvc'),remcvc_io)
        call check_ok(__FILE__,__LINE__,'Cannot read remcvc')
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
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'dstor'),dstor_io)
      call check_ok(__FILE__,__LINE__,'Cannot read dstor')
      ncstatus = nf90_get_var(ncid,get_varid(ncid,'hstor'),hstor_io)
      call check_ok(__FILE__,__LINE__,'Cannot read hstor')
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
      ncstatus = nf90_close(ncid)
      call check_ok(__FILE__,__LINE__,'Cannot close savefile '//trim(ffin))
    end if
  end subroutine read_savefile

  subroutine write_savefile(idate,ltmp)
    use netcdf
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    logical , intent(in) :: ltmp
    integer(ik4) :: ncid
    integer(ik4) , dimension(maxdims) :: dimids
    integer(ik4) , dimension(maxdims) :: wrkdim
    integer(ik4) , dimension(maxvars) :: varids
    character(256) :: ffout
    character(32) :: fbname
#ifdef CLM
    integer(ik4) :: ioff
#endif

    if ( myid == iocpu ) then
      if (ltmp) then
        write (fbname, '(a,i10)') 'TMPSAV.', toint10(idate)
      else
        write (fbname, '(a,i10)') 'SAV.', toint10(idate)
      end if
      ffout = trim(dirout)//pthsep//trim(domname)//'_'//trim(fbname)//'.nc'
      ncstatus = nf90_create(ffout,nf90_clobber,ncid)
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
      ncstatus = nf90_def_dim(ncid,'nsplit',nsplit,dimids(idnsplit))
      call check_ok(__FILE__,__LINE__,'Cannot create dimension nsplit')
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

      wrkdim(1) = dimids(idjdot)
      wrkdim(2) = dimids(ididot)
      wrkdim(3) = dimids(idkh)
      ncstatus = nf90_def_var(ncid,'atm1_u',nf90_double,wrkdim(1:3),varids(1))
      call check_ok(__FILE__,__LINE__,'Cannot create var atm1_u')
      ncstatus = nf90_def_var(ncid,'atm2_u',nf90_double,wrkdim(1:3),varids(2))
      call check_ok(__FILE__,__LINE__,'Cannot create var atm2_u')
      ncstatus = nf90_def_var(ncid,'atm1_v',nf90_double,wrkdim(1:3),varids(3))
      call check_ok(__FILE__,__LINE__,'Cannot create var atm1_v')
      ncstatus = nf90_def_var(ncid,'atm2_v',nf90_double,wrkdim(1:3),varids(4))
      call check_ok(__FILE__,__LINE__,'Cannot create var atm2_v')
      wrkdim(1) = dimids(idjcross)
      wrkdim(2) = dimids(idicross)
      wrkdim(3) = dimids(idkh)
      ncstatus = nf90_def_var(ncid,'atm1_t',nf90_double,wrkdim(1:3),varids(5))
      call check_ok(__FILE__,__LINE__,'Cannot create var atm1_t')
      ncstatus = nf90_def_var(ncid,'atm2_t',nf90_double,wrkdim(1:3),varids(6))
      call check_ok(__FILE__,__LINE__,'Cannot create var atm2_t')
      wrkdim(4) = dimids(idnqx)
      ncstatus = nf90_def_var(ncid,'atm1_qx',nf90_double,wrkdim(1:4),varids(7))
      call check_ok(__FILE__,__LINE__,'Cannot create var atm1_qx')
      ncstatus = nf90_def_var(ncid,'atm2_qx',nf90_double,wrkdim(1:4),varids(8))
      call check_ok(__FILE__,__LINE__,'Cannot create var atm2_qx')
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        wrkdim(3) = dimids(idkf)
        ncstatus = nf90_def_var(ncid,'atm1_tke',nf90_double, &
                                wrkdim(1:3),varids(9))
        call check_ok(__FILE__,__LINE__,'Cannot create var atm1_tke')
        ncstatus = nf90_def_var(ncid,'atm2_tke',nf90_double, &
                                wrkdim(1:3),varids(10))
        call check_ok(__FILE__,__LINE__,'Cannot create var atm2_tke')
        ncstatus = nf90_def_var(ncid,'kpbl',nf90_int,wrkdim(1:2),varids(11))
        call check_ok(__FILE__,__LINE__,'Cannot create var kpbl')
      end if
      ncstatus = nf90_def_var(ncid,'psa',nf90_double,wrkdim(1:2),varids(12))
      call check_ok(__FILE__,__LINE__,'Cannot create var psa')
      ncstatus = nf90_def_var(ncid,'psb',nf90_double,wrkdim(1:2),varids(13))
      call check_ok(__FILE__,__LINE__,'Cannot create var psb')
      ncstatus = nf90_def_var(ncid,'tga',nf90_double,wrkdim(1:2),varids(14))
      call check_ok(__FILE__,__LINE__,'Cannot create var tga')
      ncstatus = nf90_def_var(ncid,'tgb',nf90_double,wrkdim(1:2),varids(15))
      call check_ok(__FILE__,__LINE__,'Cannot create var tgb')
      ncstatus = nf90_def_var(ncid,'hfx',nf90_double,wrkdim(1:2),varids(16))
      call check_ok(__FILE__,__LINE__,'Cannot create var hfx')
      ncstatus = nf90_def_var(ncid,'qfx',nf90_double,wrkdim(1:2),varids(17))
      call check_ok(__FILE__,__LINE__,'Cannot create var qfx')
      ncstatus = nf90_def_var(ncid,'rainc',nf90_double,wrkdim(1:2),varids(18))
      call check_ok(__FILE__,__LINE__,'Cannot create var rainc')
      ncstatus = nf90_def_var(ncid,'rainnc',nf90_double,wrkdim(1:2),varids(19))
      call check_ok(__FILE__,__LINE__,'Cannot create var rainnc')
      ncstatus = nf90_def_var(ncid,'tgbb',nf90_double,wrkdim(1:2),varids(20))
      call check_ok(__FILE__,__LINE__,'Cannot create var tgbb')
      ncstatus = nf90_def_var(ncid,'uvdrag',nf90_double,wrkdim(1:2),varids(21))
      call check_ok(__FILE__,__LINE__,'Cannot create var uvdrag')
      wrkdim(3) = dimids(idkh)
      if ( any(icup == 1) ) then
        ncstatus = nf90_def_var(ncid,'rsheat',nf90_double, &
                                wrkdim(1:3),varids(22))
        call check_ok(__FILE__,__LINE__,'Cannot create var rsheat')
        ncstatus = nf90_def_var(ncid,'rswat',nf90_double, &
                                wrkdim(1:3),varids(23))
        call check_ok(__FILE__,__LINE__,'Cannot create var rswat')
      end if
      if ( any(icup == 3) ) then
        ncstatus = nf90_def_var(ncid,'cldefi',nf90_double, &
                                wrkdim(1:2),varids(24))
        call check_ok(__FILE__,__LINE__,'Cannot create var cldefi')
        ncstatus = nf90_def_var(ncid,'tbase',nf90_double, &
                                wrkdim(1:3),varids(25))
        call check_ok(__FILE__,__LINE__,'Cannot create var tbase')
      end if
      if ( any(icup == 4) ) then
        ncstatus = nf90_def_var(ncid,'cbmf2d',nf90_double, &
                                wrkdim(1:2),varids(26))
        call check_ok(__FILE__,__LINE__,'Cannot create var cbmf2d')
      end if
      if ( irrtm == 0 ) then
        wrkdim(3) = dimids(idkh)
        wrkdim(4) = dimids(idspw)
        ncstatus = nf90_def_var(ncid,'gasabsnxt',nf90_double, &
                                wrkdim(1:4),varids(27))
        call check_ok(__FILE__,__LINE__,'Cannot create var gasabsnxt')
        wrkdim(3) = dimids(idkf)
        wrkdim(4) = dimids(idkf)
        ncstatus = nf90_def_var(ncid,'gasabstot',nf90_double, &
                                wrkdim(1:4),varids(28))
        call check_ok(__FILE__,__LINE__,'Cannot create var gasabstot')
        ncstatus = nf90_def_var(ncid,'gasemstot',nf90_double, &
                                wrkdim(1:3),varids(29))
        call check_ok(__FILE__,__LINE__,'Cannot create var gasemstot')
      end if
      if ( ipptls > 0 ) then
        wrkdim(3) = dimids(idkh)
        ncstatus = nf90_def_var(ncid,'fcc',nf90_double,wrkdim(1:3),varids(30))
        call check_ok(__FILE__,__LINE__,'Cannot create var fcc')
        if ( ipptls == 2 ) then
          ncstatus = nf90_def_var(ncid,'snownc',nf90_double, &
                                  wrkdim(1:2),varids(31))
          call check_ok(__FILE__,__LINE__,'Cannot create var snownc')
        end if
      end if
      ncstatus = nf90_def_var(ncid,'solis',nf90_double,wrkdim(1:2),varids(32))
      call check_ok(__FILE__,__LINE__,'Cannot create var solis')
      ncstatus = nf90_def_var(ncid,'solvs',nf90_double,wrkdim(1:2),varids(33))
      call check_ok(__FILE__,__LINE__,'Cannot create var solvs')
      ncstatus = nf90_def_var(ncid,'solvsd',nf90_double,wrkdim(1:2),varids(34))
      call check_ok(__FILE__,__LINE__,'Cannot create var solvsd')
      ncstatus = nf90_def_var(ncid,'solvl',nf90_double,wrkdim(1:2),varids(35))
      call check_ok(__FILE__,__LINE__,'Cannot create var solvl')
      ncstatus = nf90_def_var(ncid,'solvld',nf90_double,wrkdim(1:2),varids(36))
      call check_ok(__FILE__,__LINE__,'Cannot create var solvld')
      ncstatus = nf90_def_var(ncid,'sabveg',nf90_double,wrkdim(1:2),varids(37))
      call check_ok(__FILE__,__LINE__,'Cannot create var sabveg')
      wrkdim(1) = dimids(idnnsg)
      wrkdim(2) = dimids(idjcross)
      wrkdim(3) = dimids(idicross)
      ncstatus = nf90_def_var(ncid,'tlef',nf90_double,wrkdim(1:3),varids(38))
      call check_ok(__FILE__,__LINE__,'Cannot create var tlef')
      ncstatus = nf90_def_var(ncid,'ssw',nf90_double,wrkdim(1:3),varids(39))
      call check_ok(__FILE__,__LINE__,'Cannot create var ssw')
      ncstatus = nf90_def_var(ncid,'rsw',nf90_double,wrkdim(1:3),varids(40))
      call check_ok(__FILE__,__LINE__,'Cannot create var rsw')
      ncstatus = nf90_def_var(ncid,'tsw',nf90_double,wrkdim(1:3),varids(41))
      call check_ok(__FILE__,__LINE__,'Cannot create var tsw')
      ncstatus = nf90_def_var(ncid,'tgrd',nf90_double,wrkdim(1:3),varids(42))
      call check_ok(__FILE__,__LINE__,'Cannot create var tgrd')
      ncstatus = nf90_def_var(ncid,'tgbrd',nf90_double,wrkdim(1:3),varids(43))
      call check_ok(__FILE__,__LINE__,'Cannot create var tgbrd')
      ncstatus = nf90_def_var(ncid,'sncv',nf90_double,wrkdim(1:3),varids(44))
      call check_ok(__FILE__,__LINE__,'Cannot create var sncv')
      ncstatus = nf90_def_var(ncid,'gwet',nf90_double,wrkdim(1:3),varids(45))
      call check_ok(__FILE__,__LINE__,'Cannot create var gwet')
      ncstatus = nf90_def_var(ncid,'snag',nf90_double,wrkdim(1:3),varids(46))
      call check_ok(__FILE__,__LINE__,'Cannot create var snag')
      ncstatus = nf90_def_var(ncid,'sfice',nf90_double,wrkdim(1:3),varids(47))
      call check_ok(__FILE__,__LINE__,'Cannot create var sfice')
      ncstatus = nf90_def_var(ncid,'ldew',nf90_double,wrkdim(1:3),varids(48))
      call check_ok(__FILE__,__LINE__,'Cannot create var ldew')
      ncstatus = nf90_def_var(ncid,'ldmsk1',nf90_int,wrkdim(1:3),varids(49))
      call check_ok(__FILE__,__LINE__,'Cannot create var ldmsk1')
      ncstatus = nf90_def_var(ncid,'taf',nf90_double,wrkdim(1:3),varids(50))
      call check_ok(__FILE__,__LINE__,'Cannot create var taf')
      ncstatus = nf90_def_var(ncid,'emiss',nf90_double,wrkdim(1:3),varids(51))
      call check_ok(__FILE__,__LINE__,'Cannot create var emiss')
      ncstatus = nf90_def_var(ncid,'scvk',nf90_double,wrkdim(1:3),varids(52))
      call check_ok(__FILE__,__LINE__,'Cannot create var scvk')
      if ( idcsst == 1 ) then
        ncstatus = nf90_def_var(ncid,'sst',nf90_double, &
                                wrkdim(1:3),varids(53))
        ncstatus = nf90_def_var(ncid,'tskin',nf90_double, &
                                wrkdim(1:3),varids(54))
        call check_ok(__FILE__,__LINE__,'Cannot create var tskin')
        ncstatus = nf90_def_var(ncid,'deltas',nf90_double, &
                                wrkdim(1:3),varids(55))
        call check_ok(__FILE__,__LINE__,'Cannot create var deltas')
        ncstatus = nf90_def_var(ncid,'tdeltas',nf90_double, &
                                wrkdim(1:3),varids(56))
        call check_ok(__FILE__,__LINE__,'Cannot create var tdeltas')
      end if
#ifndef CLM
      if ( lakemod == 1 ) then
        ncstatus = nf90_def_var(ncid,'eta',nf90_double,wrkdim(1:3),varids(57))
        call check_ok(__FILE__,__LINE__,'Cannot create var eta')
        ncstatus = nf90_def_var(ncid,'hi',nf90_double,wrkdim(1:3),varids(58))
        call check_ok(__FILE__,__LINE__,'Cannot create var hi')
        wrkdim(4) = dimids(iddpt)
        ncstatus = nf90_def_var(ncid,'tlak',nf90_double,wrkdim(1:4),varids(59))
        call check_ok(__FILE__,__LINE__,'Cannot create var tlak')
      end if
#endif
      wrkdim(1) = dimids(idjcross)
      wrkdim(2) = dimids(idicross)
      wrkdim(3) = dimids(idkh)
      ncstatus = nf90_def_var(ncid,'heatrt',nf90_double,wrkdim(1:3),varids(60))
      call check_ok(__FILE__,__LINE__,'Cannot create var heatrt')
      wrkdim(3) = dimids(idkf)
      ncstatus = nf90_def_var(ncid,'o3prof',nf90_double,wrkdim(1:3),varids(61))
      call check_ok(__FILE__,__LINE__,'Cannot create var o3prof')
      ncstatus = nf90_def_var(ncid,'flw',nf90_double,wrkdim(1:2),varids(62))
      call check_ok(__FILE__,__LINE__,'Cannot create var flw')
      ncstatus = nf90_def_var(ncid,'flwd',nf90_double,wrkdim(1:2),varids(63))
      call check_ok(__FILE__,__LINE__,'Cannot create var flwd')
      ncstatus = nf90_def_var(ncid,'fsw',nf90_double,wrkdim(1:2),varids(64))
      call check_ok(__FILE__,__LINE__,'Cannot create var fsw')
      ncstatus = nf90_def_var(ncid,'sinc',nf90_double,wrkdim(1:2),varids(65))
      call check_ok(__FILE__,__LINE__,'Cannot create var sinc')
      ncstatus = nf90_def_var(ncid,'ldmsk',nf90_int,wrkdim(1:2),varids(66))
      call check_ok(__FILE__,__LINE__,'Cannot create var ldmsk')
      if ( iocnflx == 2 ) then
        ncstatus = nf90_def_var(ncid,'zpbl',nf90_double,wrkdim(1:2),varids(67))
        call check_ok(__FILE__,__LINE__,'Cannot create var zpbl')
      end if
      if ( ichem == 1 ) then
        wrkdim(3) = dimids(idkh)
        wrkdim(4) = dimids(idntr)
        ncstatus = nf90_def_var(ncid,'chia',nf90_double,wrkdim(1:4),varids(68))
        call check_ok(__FILE__,__LINE__,'Cannot create var chia')
        ncstatus = nf90_def_var(ncid,'chib',nf90_double,wrkdim(1:4),varids(69))
        call check_ok(__FILE__,__LINE__,'Cannot create var chib')
        if ( igaschem == 1 .and. ichsolver > 0 ) then
          wrkdim(3) = dimids(idkh)
          wrkdim(4) = dimids(idtotsp)
          ncstatus = nf90_def_var(ncid,'chemall',nf90_double, &
                                  wrkdim(1:4),varids(70))
          call check_ok(__FILE__,__LINE__,'Cannot create var chemall')
          wrkdim(3) = dimids(idkf)
          wrkdim(4) = dimids(idspi)
          ncstatus = nf90_def_var(ncid,'taucldsp',nf90_double, &
                                  wrkdim(1:4),varids(71))
          call check_ok(__FILE__,__LINE__,'Cannot create var taucldsp')
        end if
        wrkdim(3) = dimids(idkh)
        wrkdim(4) = dimids(idntr)
        ncstatus = nf90_def_var(ncid,'remlsc',nf90_double, &
                                wrkdim(1:4),varids(72))
        call check_ok(__FILE__,__LINE__,'Cannot create var remlsc')
        ncstatus = nf90_def_var(ncid,'remcvc',nf90_double, &
                                wrkdim(1:4),varids(73))
        call check_ok(__FILE__,__LINE__,'Cannot create var remcvc')
        wrkdim(3) = dimids(idntr)
        ncstatus = nf90_def_var(ncid,'remdrd',nf90_double, &
                                wrkdim(1:3),varids(74))
        call check_ok(__FILE__,__LINE__,'Cannot create var remdrd')
        ncstatus = nf90_def_var(ncid,'ssw2da',nf90_double, &
                                wrkdim(1:2),varids(75))
        call check_ok(__FILE__,__LINE__,'Cannot create var ssw2da')
        ncstatus = nf90_def_var(ncid,'sdelq',nf90_double, &
                                wrkdim(1:2),varids(76))
        call check_ok(__FILE__,__LINE__,'Cannot create var sdelq')
        ncstatus = nf90_def_var(ncid,'sdelt',nf90_double, &
                                wrkdim(1:2),varids(77))
        call check_ok(__FILE__,__LINE__,'Cannot create var sdelt')
        ncstatus = nf90_def_var(ncid,'sfracb2d',nf90_double, &
                                wrkdim(1:2),varids(78))
        call check_ok(__FILE__,__LINE__,'Cannot create var sfracb2d')
        ncstatus = nf90_def_var(ncid,'sfracs2d',nf90_double, &
                                wrkdim(1:2),varids(79))
        call check_ok(__FILE__,__LINE__,'Cannot create var sfracs2d')
        ncstatus = nf90_def_var(ncid,'sfracv2d',nf90_double, &
                                wrkdim(1:2),varids(80))
        call check_ok(__FILE__,__LINE__,'Cannot create var sfracv2d')
        ncstatus = nf90_def_var(ncid,'svegfrac2d',nf90_double, &
                                wrkdim(1:2),varids(81))
        call check_ok(__FILE__,__LINE__,'Cannot create var svegfrac2d')
      end if
      wrkdim(1) = dimids(idjdot)
      wrkdim(2) = dimids(ididot)
      wrkdim(3) = dimids(idnsplit)
      ncstatus = nf90_def_var(ncid,'dstor',nf90_double,wrkdim(1:3),varids(82))
      call check_ok(__FILE__,__LINE__,'Cannot create var dstor')
      ncstatus = nf90_def_var(ncid,'hstor',nf90_double,wrkdim(1:3),varids(83))
      call check_ok(__FILE__,__LINE__,'Cannot create var hstor')
      if ( islab_ocean == 1 .and. do_restore_sst ) then
        wrkdim(1) = dimids(idjcross)
        wrkdim(2) = dimids(idicross)
        ncstatus = nf90_def_var(ncid,'qflux_restore_sst',nf90_double, &
                                wrkdim(1:2),varids(84))
        call check_ok(__FILE__,__LINE__,'Cannot create var qflux_restore_sst')
        wrkdim(1) = dimids(idmonth)
        ncstatus = nf90_def_var(ncid,'stepcount',nf90_double, &
                                wrkdim(1:1),varids(85))
        call check_ok(__FILE__,__LINE__,'Cannot create var stepcount')
      end if
#ifdef CLM
      wrkdim(1) = dimids(idjcross)
      wrkdim(2) = dimids(idicross)
      if ( imask == 2 ) then
        ncstatus = nf90_def_var(ncid,'lndcat',nf90_double,wrkdim(1:2), &
                                varids(86))
        call check_ok(__FILE__,__LINE__,'Cannot create var lndcat')
      end if
#endif

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

      ncstatus = nf90_put_var(ncid,varids(1),atm1_u_io)
      call check_ok(__FILE__,__LINE__,'Cannot write atm1_u')
      ncstatus = nf90_put_var(ncid,varids(2),atm2_u_io)
      call check_ok(__FILE__,__LINE__,'Cannot write atm2_u')
      ncstatus = nf90_put_var(ncid,varids(3),atm1_v_io)
      call check_ok(__FILE__,__LINE__,'Cannot write atm1_v')
      ncstatus = nf90_put_var(ncid,varids(4),atm2_v_io)
      call check_ok(__FILE__,__LINE__,'Cannot write atm2_v')
      ncstatus = nf90_put_var(ncid,varids(5),atm1_t_io)
      call check_ok(__FILE__,__LINE__,'Cannot write atm1_t')
      ncstatus = nf90_put_var(ncid,varids(6),atm2_t_io)
      call check_ok(__FILE__,__LINE__,'Cannot write atm2_t')
      ncstatus = nf90_put_var(ncid,varids(7),atm1_qx_io)
      call check_ok(__FILE__,__LINE__,'Cannot write atm1_qx')
      ncstatus = nf90_put_var(ncid,varids(8),atm2_qx_io)
      call check_ok(__FILE__,__LINE__,'Cannot write atm2_qx')
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        ncstatus = nf90_put_var(ncid,varids(9),atm1_tke_io)
        call check_ok(__FILE__,__LINE__,'Cannot write atm1_tke')
        ncstatus = nf90_put_var(ncid,varids(10),atm2_tke_io)
        call check_ok(__FILE__,__LINE__,'Cannot write atm2_tke')
        ncstatus = nf90_put_var(ncid,varids(11),kpbl_io)
        call check_ok(__FILE__,__LINE__,'Cannot write kpbl')
      end if
      ncstatus = nf90_put_var(ncid,varids(12),psa_io)
      call check_ok(__FILE__,__LINE__,'Cannot write psa')
      ncstatus = nf90_put_var(ncid,varids(13),psb_io)
      call check_ok(__FILE__,__LINE__,'Cannot write psb')
      ncstatus = nf90_put_var(ncid,varids(14),tga_io)
      call check_ok(__FILE__,__LINE__,'Cannot write tga')
      ncstatus = nf90_put_var(ncid,varids(15),tgb_io)
      call check_ok(__FILE__,__LINE__,'Cannot write tgb')
      ncstatus = nf90_put_var(ncid,varids(16),hfx_io)
      call check_ok(__FILE__,__LINE__,'Cannot write hfx')
      ncstatus = nf90_put_var(ncid,varids(17),qfx_io)
      call check_ok(__FILE__,__LINE__,'Cannot write qfx')
      ncstatus = nf90_put_var(ncid,varids(18),rainc_io)
      call check_ok(__FILE__,__LINE__,'Cannot write rainc')
      ncstatus = nf90_put_var(ncid,varids(19),rainnc_io)
      call check_ok(__FILE__,__LINE__,'Cannot write rainnc')
      ncstatus = nf90_put_var(ncid,varids(20),tgbb_io)
      call check_ok(__FILE__,__LINE__,'Cannot write tgbb')
      ncstatus = nf90_put_var(ncid,varids(21),uvdrag_io)
      call check_ok(__FILE__,__LINE__,'Cannot write uvdrag')
      if ( any(icup == 1) ) then
        ncstatus = nf90_put_var(ncid,varids(22),rsheat_io)
        call check_ok(__FILE__,__LINE__,'Cannot write rsheat')
        ncstatus = nf90_put_var(ncid,varids(23),rswat_io)
        call check_ok(__FILE__,__LINE__,'Cannot write rswat')
      end if
      if ( any(icup == 3) ) then
        ncstatus = nf90_put_var(ncid,varids(24),cldefi_io)
        call check_ok(__FILE__,__LINE__,'Cannot write cldefi')
        ncstatus = nf90_put_var(ncid,varids(25),tbase_io)
        call check_ok(__FILE__,__LINE__,'Cannot write tbase')
      end if
      if ( any(icup == 4) ) then
        ncstatus = nf90_put_var(ncid,varids(26),cbmf2d_io)
        call check_ok(__FILE__,__LINE__,'Cannot write cbmf2d')
      end if
      if ( irrtm == 0 ) then
        ncstatus = nf90_put_var(ncid,varids(27),gasabsnxt_io)
        call check_ok(__FILE__,__LINE__,'Cannot write gasabsnxt')
        ncstatus = nf90_put_var(ncid,varids(28),gasabstot_io)
        call check_ok(__FILE__,__LINE__,'Cannot write gasabstot')
        ncstatus = nf90_put_var(ncid,varids(29),gasemstot_io)
        call check_ok(__FILE__,__LINE__,'Cannot write gasemstot')
      end if
      if ( ipptls > 0 ) then
        ncstatus = nf90_put_var(ncid,varids(30),fcc_io)
        call check_ok(__FILE__,__LINE__,'Cannot write fcc')
        if ( ipptls == 2 ) then
          ncstatus = nf90_put_var(ncid,varids(31),snownc_io)
          call check_ok(__FILE__,__LINE__,'Cannot write snownc')
        end if
      end if
      ncstatus = nf90_put_var(ncid,varids(32),solis_io)
      call check_ok(__FILE__,__LINE__,'Cannot write solis')
      ncstatus = nf90_put_var(ncid,varids(33),solvs_io)
      call check_ok(__FILE__,__LINE__,'Cannot write solvs')
      ncstatus = nf90_put_var(ncid,varids(34),solvsd_io)
      call check_ok(__FILE__,__LINE__,'Cannot write solvsd')
      ncstatus = nf90_put_var(ncid,varids(35),solvl_io)
      call check_ok(__FILE__,__LINE__,'Cannot write solvl')
      ncstatus = nf90_put_var(ncid,varids(36),solvld_io)
      call check_ok(__FILE__,__LINE__,'Cannot write solvld')
      ncstatus = nf90_put_var(ncid,varids(37),sabveg_io)
      call check_ok(__FILE__,__LINE__,'Cannot write sabveg')
      ncstatus = nf90_put_var(ncid,varids(38),tlef_io)
      call check_ok(__FILE__,__LINE__,'Cannot write tlef')
      ncstatus = nf90_put_var(ncid,varids(39),ssw_io)
      call check_ok(__FILE__,__LINE__,'Cannot write ssw')
      ncstatus = nf90_put_var(ncid,varids(40),rsw_io)
      call check_ok(__FILE__,__LINE__,'Cannot write rsw')
      ncstatus = nf90_put_var(ncid,varids(41),tsw_io)
      call check_ok(__FILE__,__LINE__,'Cannot write tsw')
      ncstatus = nf90_put_var(ncid,varids(42),tgrd_io)
      call check_ok(__FILE__,__LINE__,'Cannot write tgrd')
      ncstatus = nf90_put_var(ncid,varids(43),tgbrd_io)
      call check_ok(__FILE__,__LINE__,'Cannot write tgbrd')
      ncstatus = nf90_put_var(ncid,varids(44),sncv_io)
      call check_ok(__FILE__,__LINE__,'Cannot write sncv')
      ncstatus = nf90_put_var(ncid,varids(45),gwet_io)
      call check_ok(__FILE__,__LINE__,'Cannot write gwet')
      ncstatus = nf90_put_var(ncid,varids(46),snag_io)
      call check_ok(__FILE__,__LINE__,'Cannot write snag')
      ncstatus = nf90_put_var(ncid,varids(47),sfice_io)
      call check_ok(__FILE__,__LINE__,'Cannot write sfice')
      ncstatus = nf90_put_var(ncid,varids(48),ldew_io)
      call check_ok(__FILE__,__LINE__,'Cannot write ldew')
      ncstatus = nf90_put_var(ncid,varids(49),ldmsk1_io)
      call check_ok(__FILE__,__LINE__,'Cannot write ldmsk1')
      ncstatus = nf90_put_var(ncid,varids(50),taf_io)
      call check_ok(__FILE__,__LINE__,'Cannot write taf')
      ncstatus = nf90_put_var(ncid,varids(51),emisv_io)
      call check_ok(__FILE__,__LINE__,'Cannot write emiss')
      ncstatus = nf90_put_var(ncid,varids(52),scvk_io)
      call check_ok(__FILE__,__LINE__,'Cannot write scvk')
      if ( idcsst == 1 ) then
        ncstatus = nf90_put_var(ncid,varids(53),sst_io)
        call check_ok(__FILE__,__LINE__,'Cannot write sst')
        ncstatus = nf90_put_var(ncid,varids(54),tskin_io)
        call check_ok(__FILE__,__LINE__,'Cannot write tskin')
        ncstatus = nf90_put_var(ncid,varids(55),deltas_io)
        call check_ok(__FILE__,__LINE__,'Cannot write deltas')
        ncstatus = nf90_put_var(ncid,varids(56),tdeltas_io)
        call check_ok(__FILE__,__LINE__,'Cannot write tdeltas')
      end if
#ifndef CLM
      if ( lakemod == 1 ) then
        ncstatus = nf90_put_var(ncid,varids(57),eta_io)
        call check_ok(__FILE__,__LINE__,'Cannot write eta')
        ncstatus = nf90_put_var(ncid,varids(58),hi_io)
        call check_ok(__FILE__,__LINE__,'Cannot write hi')
        ncstatus = nf90_put_var(ncid,varids(59),tlak_io)
        call check_ok(__FILE__,__LINE__,'Cannot write tlak')
      end if
#endif
      ncstatus = nf90_put_var(ncid,varids(60),heatrt_io)
      call check_ok(__FILE__,__LINE__,'Cannot write heatrt')
      ncstatus = nf90_put_var(ncid,varids(61),o3prof_io)
      call check_ok(__FILE__,__LINE__,'Cannot write o3prof')
      ncstatus = nf90_put_var(ncid,varids(62),flw_io)
      call check_ok(__FILE__,__LINE__,'Cannot write flw')
      ncstatus = nf90_put_var(ncid,varids(63),flwd_io)
      call check_ok(__FILE__,__LINE__,'Cannot write flwd')
      ncstatus = nf90_put_var(ncid,varids(64),fsw_io)
      call check_ok(__FILE__,__LINE__,'Cannot write fsw')
      ncstatus = nf90_put_var(ncid,varids(65),sinc_io)
      call check_ok(__FILE__,__LINE__,'Cannot write sinc')
      ncstatus = nf90_put_var(ncid,varids(66),ldmsk_io)
      call check_ok(__FILE__,__LINE__,'Cannot write ldmsk')
      if ( iocnflx == 2 ) then
        ncstatus = nf90_put_var(ncid,varids(67),zpbl_io)
        call check_ok(__FILE__,__LINE__,'Cannot write zpbl')
      end if
      if ( ichem == 1 ) then
        ncstatus = nf90_put_var(ncid,varids(68),chia_io)
        call check_ok(__FILE__,__LINE__,'Cannot write chia')
        ncstatus = nf90_put_var(ncid,varids(69),chib_io)
        call check_ok(__FILE__,__LINE__,'Cannot write chib')
        if ( igaschem == 1 .and. ichsolver > 0 ) then
          ncstatus = nf90_put_var(ncid,varids(70),chemall_io)
          call check_ok(__FILE__,__LINE__,'Cannot write chemall')
          ncstatus = nf90_put_var(ncid,varids(71),taucldsp_io)
          call check_ok(__FILE__,__LINE__,'Cannot write taucldsp')
        end if
        ncstatus = nf90_put_var(ncid,varids(72),remlsc_io)
        call check_ok(__FILE__,__LINE__,'Cannot write remlsc')
        ncstatus = nf90_put_var(ncid,varids(73),remcvc_io)
        call check_ok(__FILE__,__LINE__,'Cannot write remcvc')
        ncstatus = nf90_put_var(ncid,varids(74),remdrd_io)
        call check_ok(__FILE__,__LINE__,'Cannot write remdrd')
        ncstatus = nf90_put_var(ncid,varids(75),ssw2da_io)
        call check_ok(__FILE__,__LINE__,'Cannot write ssw2da')
        ncstatus = nf90_put_var(ncid,varids(76),sdelq_io)
        call check_ok(__FILE__,__LINE__,'Cannot write sdelq')
        ncstatus = nf90_put_var(ncid,varids(77),sdelt_io)
        call check_ok(__FILE__,__LINE__,'Cannot write sdelt')
        ncstatus = nf90_put_var(ncid,varids(78),sfracb2d_io)
        call check_ok(__FILE__,__LINE__,'Cannot write sfracb2d')
        ncstatus = nf90_put_var(ncid,varids(79),sfracs2d_io)
        call check_ok(__FILE__,__LINE__,'Cannot write sfracs2d')
        ncstatus = nf90_put_var(ncid,varids(80),sfracv2d_io)
        call check_ok(__FILE__,__LINE__,'Cannot write sfracv2d')
        ncstatus = nf90_put_var(ncid,varids(81),svegfrac2d_io)
        call check_ok(__FILE__,__LINE__,'Cannot write svegfrac2d')
      end if
      ncstatus = nf90_put_var(ncid,varids(82),dstor_io)
      call check_ok(__FILE__,__LINE__,'Cannot write dstor')
      ncstatus = nf90_put_var(ncid,varids(83),hstor_io)
      call check_ok(__FILE__,__LINE__,'Cannot write hstor')
      if ( islab_ocean == 1 .and. do_restore_sst ) then
        ncstatus = nf90_put_var(ncid,varids(84),qflux_restore_sst_io)
        call check_ok(__FILE__,__LINE__,'Cannot write qflux_restore_sst')
        ncstatus = nf90_put_var(ncid,varids(85),stepcount)
        call check_ok(__FILE__,__LINE__,'Cannot write stepcount')
      end if
#ifdef CLM
      if ( imask == 2 ) then
        ncstatus = nf90_put_var(ncid,varids(86),lndcat_io)
        call check_ok(__FILE__,__LINE__,'Cannot write lndcat')
      end if
#endif
      ncstatus = nf90_close(ncid)
      call check_ok(__FILE__,__LINE__,'Cannot close savefile '//trim(ffout))
    end if

#ifdef CLM
    ioff = int(dtsec*dble(ntsrf-1))
    filer_rest = restFile_filename(type='netcdf',offset=ioff)
    call restFile_write(filer_rest)
    filer_rest = restFile_filename(type='binary',offset=ioff)
    call restFile_write_binary(filer_rest)
    thisclmrest = filer_rest(1:256)
#endif

    if ( myid == iocpu ) then
      write(stdout,*) 'SAV variables written at ', tochar(idate)
      if ( isavlast > 0 ) then
        write (fbname, '(a,i10)') 'TMPSAV.', isavlast
        ffout = trim(dirout)//pthsep//trim(domname)//'_'//trim(fbname)//'.nc'
        call unlink(ffout)
#ifdef CLM
        call unlink(trim(lastclmrest))
        call unlink((trim(lastclmrest)//'.nc'))
#endif            
      end if
      if (ltmp) then
        isavlast = toint10(idate)
      else
        isavlast = 0
      end if
#ifdef CLM
      lastclmrest = thisclmrest
#endif
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
    
end module mod_savefile
