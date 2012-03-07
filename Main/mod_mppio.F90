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
  use mod_che_interface
  use mod_rad_interface
  use mod_pbl_interface , only : ibltyp , tcm_state , allocate_tcm_state
  use mod_memutil
  use mod_mpmessage
!
  integer , pointer , dimension(:,:,:) :: ldmsk1_io
  integer , pointer , dimension(:,:,:) :: iveg1_io
  integer , pointer , dimension(:,:) :: iveg_io
  integer , pointer , dimension(:,:) :: ldmsk_io

  real(dp) , pointer , dimension(:,:,:) :: ldew_io
  real(dp) , pointer , dimension(:,:,:) :: gwet_io
  real(dp) , pointer , dimension(:,:,:) :: snag_io
  real(dp) , pointer , dimension(:,:,:) :: sncv_io
  real(dp) , pointer , dimension(:,:,:) :: sfice_io
  real(dp) , pointer , dimension(:,:,:) :: rsw_io
  real(dp) , pointer , dimension(:,:,:) :: ssw_io
  real(dp) , pointer , dimension(:,:,:) :: tsw_io
  real(dp) , pointer , dimension(:,:,:) :: taf_io
  real(dp) , pointer , dimension(:,:,:) :: tgrd_io
  real(dp) , pointer , dimension(:,:,:) :: tgbrd_io
  real(dp) , pointer , dimension(:,:,:) :: tlef_io
  real(dp) , pointer , dimension(:,:,:) :: emiss_io

  real(dp) , pointer , dimension(:,:,:) :: ht1_io
  real(dp) , pointer , dimension(:,:,:) :: lndcat1_io
  real(dp) , pointer , dimension(:,:,:) :: xlat1_io
  real(dp) , pointer , dimension(:,:,:) :: xlon1_io
!
  real(dp) , pointer , dimension(:,:) :: flw_io
  real(dp) , pointer , dimension(:,:) :: flwd_io
  real(dp) , pointer , dimension(:,:) :: fsw_io
  real(dp) , pointer , dimension(:,:) :: sabveg_io
  real(dp) , pointer , dimension(:,:) :: sinc_io
  real(dp) , pointer , dimension(:,:) :: solis_io
  real(dp) , pointer , dimension(:,:) :: solvd_io
  real(dp) , pointer , dimension(:,:) :: solvs_io

  integer , pointer , dimension(:,:) :: kpbl_io
  real(dp) , pointer , dimension(:,:) :: zpbl_io
!
  real(sp) , pointer , dimension(:,:,:) :: fbat_io
  real(sp) , pointer , dimension(:,:,:,:) :: fsub_io
  real(sp) , pointer , dimension(:,:,:) :: frad2d_io
  real(sp) , pointer , dimension(:,:,:,:) :: frad3d_io
  real(sp) , pointer , dimension(:,:) :: radpsa_io

  real(dp) , pointer , dimension(:,:) :: cbmf2d_io
  real(dp) , pointer , dimension(:,:,:) :: fcc_io
  real(dp) , pointer , dimension(:,:,:) :: rsheat_io
  real(dp) , pointer , dimension(:,:,:) :: rswat_io

  real(dp) , pointer , dimension(:,:,:,:) :: gasabsnxt_io
  real(dp) , pointer , dimension(:,:,:,:) :: gasabstot_io
  real(dp) , pointer , dimension(:,:,:) :: gasemstot_io

  real(dp) , pointer , dimension(:,:,:) :: heatrt_io
  real(dp) , pointer , dimension(:,:,:) :: o3prof_io

  real(dp) , pointer , dimension(:,:,:) :: dstor_io
  real(dp) , pointer , dimension(:,:,:) :: hstor_io

  real(dp) , pointer , dimension(:,:) :: ps0_io
  real(dp) , pointer , dimension(:,:) :: ps1_io
  real(dp) , pointer , dimension(:,:) :: ts0_io
  real(dp) , pointer , dimension(:,:) :: ts1_io

  real(dp) , pointer , dimension(:,:,:) :: qb0_io
  real(dp) , pointer , dimension(:,:,:) :: qb1_io
  real(dp) , pointer , dimension(:,:,:) :: tb0_io
  real(dp) , pointer , dimension(:,:,:) :: tb1_io
  real(dp) , pointer , dimension(:,:,:) :: ub0_io
  real(dp) , pointer , dimension(:,:,:) :: ub1_io
  real(dp) , pointer , dimension(:,:,:) :: vb0_io
  real(dp) , pointer , dimension(:,:,:) :: vb1_io

  real(dp) , pointer , dimension(:,:) :: cldefi_io

  real(dp) , pointer , dimension(:,:,:) :: omega_io
  real(dp) , pointer , dimension(:,:,:) :: tbase_io

  type(domain) :: mddom_io
  type(atmstate) :: atm1_io
  type(atmstate) :: atm2_io
  type(surfstate) :: sfs_io
  type(tcm_state) :: tcmstate_io

#ifdef CLM
  real(dp) , pointer , dimension(:,:) :: sols2d_io
  real(dp) , pointer , dimension(:,:) :: soll2d_io
  real(dp) , pointer , dimension(:,:) :: solsd2d_io
  real(dp) , pointer , dimension(:,:) :: solld2d_io
  real(dp) , pointer , dimension(:,:) :: aldifl2d_io
  real(dp) , pointer , dimension(:,:) :: aldirs2d_io
  real(dp) , pointer , dimension(:,:) :: aldirl2d_io
  real(dp) , pointer , dimension(:,:) :: aldifs2d_io
  real(dp) , pointer , dimension(:,:) :: lndcat2d_io
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

    if (myid == 0) then
      call allocate_domain(mddom_io,.false.)
      call allocate_atmstate(atm1_io,ibltyp,.false.,0,0)
      call allocate_atmstate(atm2_io,ibltyp,.false.,0,0)
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
      call getmem3d(ldmsk1_io,1,nnsg,jcross1,jcross2,icross1,icross2,'ldmsk1_io')
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

      call getmem3d(ht1_io,1,nnsg,jdot1,jdot2,idot1,idot2,'ht1_io')
      call getmem3d(lndcat1_io,1,nnsg,jdot1,jdot2,idot1,idot2,'lndcat1_io')
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

      call getmem2d(zpbl_io,jcross1,jcross2,icross1,icross2,'zpbl_io')
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
      call getmem2d(radpsa_io,jout1,jout2,iout1,iout2,'radpsa_io')
    endif

  end subroutine allocate_mod_mppio
!
  subroutine free_mpp_initspace
    implicit none
    if (ichem == 1) then
      call relmem4d(src0)
      call relmem3d(src1)
    end if
    if (myid == 0) then
      if (ichem == 1) then
        call relmem4d(src_0)
        call relmem3d(src_1)
      end if
    end if
  end subroutine free_mpp_initspace
!
end module mod_mppio
