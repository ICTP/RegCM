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

module mod_mppio

  use mod_runparams
  use mod_atm_interface
  use mod_lm_interface
  use mod_cu_interface
  use mod_rad_interface
  use mod_pbl_interface , only : tcm_state , allocate_tcm_state
  use mod_mppparam
  use mod_memutil
  use mod_mpmessage
!
  integer(ik4) , pointer , dimension(:,:,:) :: ldmsk1_io
  integer(ik4) , pointer , dimension(:,:,:) :: iveg1_io
  integer(ik4) , pointer , dimension(:,:) :: iveg_io
  integer(ik4) , pointer , dimension(:,:) :: ldmsk_io

  real(rk8) , pointer , dimension(:,:,:) :: ldew_io
  real(rk8) , pointer , dimension(:,:,:) :: gwet_io
  real(rk8) , pointer , dimension(:,:,:) :: snag_io
  real(rk8) , pointer , dimension(:,:,:) :: sncv_io
  real(rk8) , pointer , dimension(:,:,:) :: sfice_io
  real(rk8) , pointer , dimension(:,:,:) :: rsw_io
  real(rk8) , pointer , dimension(:,:,:) :: ssw_io
  real(rk8) , pointer , dimension(:,:,:) :: tsw_io
  real(rk8) , pointer , dimension(:,:,:) :: taf_io
  real(rk8) , pointer , dimension(:,:,:) :: tgrd_io
  real(rk8) , pointer , dimension(:,:,:) :: tgbrd_io
  real(rk8) , pointer , dimension(:,:,:) :: tlef_io
  real(rk8) , pointer , dimension(:,:,:) :: emiss_io

  real(rk8) , pointer , dimension(:,:,:) :: ht1_io
  real(rk8) , pointer , dimension(:,:,:) :: lndcat1_io
  real(rk8) , pointer , dimension(:,:,:) :: mask1_io
  real(rk8) , pointer , dimension(:,:,:) :: xlat1_io
  real(rk8) , pointer , dimension(:,:,:) :: xlon1_io
!
  real(rk8) , pointer , dimension(:,:) :: flw_io
  real(rk8) , pointer , dimension(:,:) :: flwd_io
  real(rk8) , pointer , dimension(:,:) :: fsw_io
  real(rk8) , pointer , dimension(:,:) :: sabveg_io
  real(rk8) , pointer , dimension(:,:) :: sinc_io
  real(rk8) , pointer , dimension(:,:) :: solis_io
  real(rk8) , pointer , dimension(:,:) :: solvd_io
  real(rk8) , pointer , dimension(:,:) :: solvs_io

  real(rk8) , pointer , dimension(:,:) :: dtskin_io
  real(rk8) , pointer , dimension(:,:) :: tdeltas_io
  real(rk8) , pointer , dimension(:,:) :: deltas_io

  integer(ik4) , pointer , dimension(:,:) :: kpbl_io
  real(rk8) , pointer , dimension(:,:) :: zpbl_io
!
  real(rk4) , pointer , dimension(:,:,:) :: fbat_io
  real(rk4) , pointer , dimension(:,:,:,:) :: fsub_io
  real(rk4) , pointer , dimension(:,:,:) :: frad2d_io
  real(rk4) , pointer , dimension(:,:,:,:) :: frad3d_io

  real(rk8) , pointer , dimension(:,:) :: cbmf2d_io
  real(rk8) , pointer , dimension(:,:,:) :: fcc_io
  real(rk8) , pointer , dimension(:,:,:) :: rsheat_io
  real(rk8) , pointer , dimension(:,:,:) :: rswat_io

  real(rk8) , pointer , dimension(:,:,:,:) :: gasabsnxt_io
  real(rk8) , pointer , dimension(:,:,:,:) :: gasabstot_io
  real(rk8) , pointer , dimension(:,:,:) :: gasemstot_io

  real(rk8) , pointer , dimension(:,:,:) :: heatrt_io
  real(rk8) , pointer , dimension(:,:,:) :: o3prof_io

  real(rk8) , pointer , dimension(:,:,:) :: dstor_io
  real(rk8) , pointer , dimension(:,:,:) :: hstor_io

  real(rk8) , pointer , dimension(:,:) :: ps0_io
  real(rk8) , pointer , dimension(:,:) :: ps1_io
  real(rk8) , pointer , dimension(:,:) :: ts0_io
  real(rk8) , pointer , dimension(:,:) :: ts1_io

  real(rk8) , pointer , dimension(:,:,:) :: qb0_io
  real(rk8) , pointer , dimension(:,:,:) :: qb1_io
  real(rk8) , pointer , dimension(:,:,:) :: tb0_io
  real(rk8) , pointer , dimension(:,:,:) :: tb1_io
  real(rk8) , pointer , dimension(:,:,:) :: ub0_io
  real(rk8) , pointer , dimension(:,:,:) :: ub1_io
  real(rk8) , pointer , dimension(:,:,:) :: vb0_io
  real(rk8) , pointer , dimension(:,:,:) :: vb1_io

  real(rk8) , pointer , dimension(:,:) :: cldefi_io

  real(rk8) , pointer , dimension(:,:,:) :: omega_io
  real(rk8) , pointer , dimension(:,:,:) :: tbase_io

  type(atmstate) :: atm1_io
  type(atmstate) :: atm2_io
  type(surfstate) :: sfs_io
  type(tcm_state) :: tcmstate_io

#ifdef CLM
  real(rk8) , pointer , dimension(:,:) :: sols2d_io
  real(rk8) , pointer , dimension(:,:) :: soll2d_io
  real(rk8) , pointer , dimension(:,:) :: solsd2d_io
  real(rk8) , pointer , dimension(:,:) :: solld2d_io
  real(rk8) , pointer , dimension(:,:) :: aldifl2d_io
  real(rk8) , pointer , dimension(:,:) :: aldirs2d_io
  real(rk8) , pointer , dimension(:,:) :: aldirl2d_io
  real(rk8) , pointer , dimension(:,:) :: aldifs2d_io
  real(rk8) , pointer , dimension(:,:) :: lndcat2d_io
#endif
!
!---------- DATA init section--------------------------------------------
!
  contains 
!
!     This routines allocate all the arrays contained in the module
!
  subroutine allocate_mod_mppio
    implicit none

    if ( myid == iocpu ) then
      call allocate_atmstate(atm1_io,ibltyp,.false.,zero_exchange_point)
      call allocate_atmstate(atm2_io,ibltyp,.false.,zero_exchange_point)
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call allocate_tcm_state(tcmstate_io,.false.)
      end if
      call allocate_surfstate(sfs_io,.false.)

      call getmem3d(ldew_io,1,nnsg,jcross1,jcross2,icross1,icross2,'ldew_io')
      call getmem3d(gwet_io,1,nnsg,jcross1,jcross2,icross1,icross2,'gwet_io')
      call getmem3d(snag_io,1,nnsg,jcross1,jcross2,icross1,icross2,'snag_io')
      call getmem3d(sncv_io,1,nnsg,jcross1,jcross2,icross1,icross2,'sncv_io')
      call getmem3d(sfice_io,1,nnsg,jcross1,jcross2,icross1,icross2,'sfice_io')
      call getmem3d(rsw_io,1,nnsg,jcross1,jcross2,icross1,icross2,'rsw_io')
      call getmem3d(ssw_io,1,nnsg,jcross1,jcross2,icross1,icross2,'ssw_io')
      call getmem3d(tsw_io,1,nnsg,jcross1,jcross2,icross1,icross2,'tsw_io')
      call getmem3d(taf_io,1,nnsg,jcross1,jcross2,icross1,icross2,'taf_io')
      call getmem3d(tgrd_io,1,nnsg,jcross1,jcross2,icross1,icross2,'tgrd_io')
      call getmem3d(tgbrd_io,1,nnsg,jcross1,jcross2,icross1,icross2,'tgbrd_io')
      call getmem3d(tlef_io,1,nnsg,jcross1,jcross2,icross1,icross2,'tlef_io')
      call getmem3d(emiss_io,1,nnsg,jcross1,jcross2,icross1,icross2,'emiss_io')
      call getmem3d(ldmsk1_io,1,nnsg,jcross1,jcross2,icross1,icross2,'ldmsk1')
      call getmem3d(iveg1_io,1,nnsg,jcross1,jcross2,icross1,icross2,'iveg1_io')
      call getmem2d(iveg_io,jcross1,jcross2,icross1,icross2,'iveg_io')
      call getmem2d(ldmsk_io,jcross1,jcross2,icross1,icross2,'ldmsk_io')
      call getmem2d(flw_io,jcross1,jcross2,icross1,icross2,'flw_io')
      call getmem2d(flwd_io,jcross1,jcross2,icross1,icross2,'flwd_io')
      call getmem2d(fsw_io,jcross1,jcross2,icross1,icross2,'fsw_io')
      call getmem2d(sabveg_io,jcross1,jcross2,icross1,icross2,'sabveg_io')
      call getmem2d(sinc_io,jcross1,jcross2,icross1,icross2,'sinc_io')
      call getmem2d(solis_io,jcross1,jcross2,icross1,icross2,'solis_io')
      call getmem2d(solvd_io,jcross1,jcross2,icross1,icross2,'solvd_io')
      call getmem2d(solvs_io,jcross1,jcross2,icross1,icross2,'solvs_io')
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call getmem2d(kpbl_io,jcross1,jcross2,icross1,icross2,'kpbl_io')
      end if
      if ( idcsst == 1 ) then
        call getmem2d(dtskin_io,jcross1,jcross2,icross1,icross2,'dtskin_io')
        call getmem2d(deltas_io,jcross1,jcross2,icross1,icross2,'deltas_io')
        call getmem2d(tdeltas_io,jcross1,jcross2,icross1,icross2,'tdeltas_io')
      end if

      call getmem3d(ht1_io,1,nnsg,jdot1,jdot2,idot1,idot2,'ht1_io')
      call getmem3d(lndcat1_io,1,nnsg,jdot1,jdot2,idot1,idot2,'lndcat1_io')
      call getmem3d(mask1_io,1,nnsg,jdot1,jdot2,idot1,idot2,'mask1_io')
      call getmem3d(xlat1_io,1,nnsg,jdot1,jdot2,idot1,idot2,'xlat1_io')
      call getmem3d(xlon1_io,1,nnsg,jdot1,jdot2,idot1,idot2,'xlon1_io')

      call getmem2d(cbmf2d_io,jcross1,jcross2,icross1,icross2,'cbmf2d_io')
      call getmem3d(fcc_io,jcross1,jcross2,icross1,icross2,1,kz,'fcc_io')
      call getmem3d(rsheat_io,jcross1,jcross2,icross1,icross2,1,kz,'rsheat_io')
      call getmem3d(rswat_io,jcross1,jcross2,icross1,icross2,1,kz,'rswat_io')

      call getmem3d(dstor_io,jdot1,jdot2,idot1,idot2,1,nsplit,'dstor_io')
      call getmem3d(hstor_io,jdot1,jdot2,idot1,idot2,1,nsplit,'hstor_io')

      call getmem4d(gasabsnxt_io,jcross1,jcross2, &
                    icross1,icross2,1,kz,1,4,'gasabsnxt_io')
      call getmem4d(gasabstot_io,jcross1,jcross2, &
                    icross1,icross2,1,kzp1,1,kzp1,'gasabstot_io')
      call getmem3d(gasemstot_io,jcross1,jcross2, &
                    icross1,icross2,1,kzp1,'gasemstot_io')

      call getmem3d(heatrt_io,jcross1,jcross2, &
                    icross1,icross2,1,kz,'heatrt_io')
      call getmem3d(o3prof_io,jcross1,jcross2, &
                    icross1,icross2,1,kzp1,'o3prof_io')

      call getmem2d(ps0_io,jdot1,jdot2,idot1,idot2,'ps0_io')
      call getmem2d(ps1_io,jdot1,jdot2,idot1,idot2,'ps1_io')
      call getmem2d(ts0_io,jdot1,jdot2,idot1,idot2,'ts0_io')
      call getmem2d(ts1_io,jdot1,jdot2,idot1,idot2,'ts1_io')
      call getmem3d(qb0_io,jdot1,jdot2,idot1,idot2,1,kz,'qb0_io')
      call getmem3d(qb1_io,jdot1,jdot2,idot1,idot2,1,kz,'qb1_io')
      call getmem3d(tb0_io,jdot1,jdot2,idot1,idot2,1,kz,'tb0_io')
      call getmem3d(tb1_io,jdot1,jdot2,idot1,idot2,1,kz,'tb1_io')
      call getmem3d(ub0_io,jdot1,jdot2,idot1,idot2,1,kz,'ub0_io')
      call getmem3d(ub1_io,jdot1,jdot2,idot1,idot2,1,kz,'ub1_io')
      call getmem3d(vb0_io,jdot1,jdot2,idot1,idot2,1,kz,'vb0_io')
      call getmem3d(vb1_io,jdot1,jdot2,idot1,idot2,1,kz,'vb1_io')

      if ( iocnflx == 2 ) then
        call getmem2d(zpbl_io,jcross1,jcross2,icross1,icross2,'zpbl_io')
      end if
      call getmem3d(omega_io,jcross1,jcross2, &
                    icross1,icross2,1,kz,'omega_io')
      if (icup == 3) then
        call getmem2d(cldefi_io,jcross1,jcross2,icross1,icross2,'cldefi_io')
        call getmem3d(tbase_io,jcross1,jcross2,icross1,icross2,1,kz,'tbase_io')
      end if
#ifdef CLM
      call getmem2d(sols2d_io,jcross1,jcross2,icross1,icross2,'sols2d_io')
      call getmem2d(soll2d_io,jcross1,jcross2,icross1,icross2,'soll2d_io')
      call getmem2d(solsd2d_io,jcross1,jcross2,icross1,icross2,'solsd2d_io')
      call getmem2d(solld2d_io,jcross1,jcross2,icross1,icross2,'solld2d_io')
      call getmem2d(aldifl2d_io,jcross1,jcross2,icross1,icross2,'aldifl2d_io')
      call getmem2d(aldifs2d_io,jcross1,jcross2,icross1,icross2,'aldifs2d_io')
      call getmem2d(aldirl2d_io,jcross1,jcross2,icross1,icross2,'aldirl2d_io')
      call getmem2d(aldirs2d_io,jcross1,jcross2,icross1,icross2,'aldirs2d_io')
      call getmem2d(lndcat2d_io,jcross1,jcross2,icross1,icross2,'lndcat2d_io')
#endif
      !
      ! Output array for SRF , RAD , SUB
      !
      call getmem3d(fbat_io,jout1,jout2,iout1,iout2,1,numbat,'fbat_io')
      call getmem4d(fsub_io,1,nnsg,jout1,jout2,iout1,iout2,1,numsub,'fsub_io')
      call getmem3d(frad2d_io,jout1,jout2,iout1,iout2,1,nrad2d,'frad2d_io')
      call getmem4d(frad3d_io,jout1,jout2,iout1,iout2,1,kz,1,nrad3d,'frad3d_io')
    endif

  end subroutine allocate_mod_mppio
!
end module mod_mppio
