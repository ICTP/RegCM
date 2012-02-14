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
  real(dp) , pointer , dimension(:) :: r8vector1
  real(dp) , pointer , dimension(:) :: r8vector2
  real(sp) , pointer , dimension(:) :: r4vector1
  real(sp) , pointer , dimension(:) :: r4vector2
  integer , pointer , dimension(:) :: i4vector1
  integer , pointer , dimension(:) :: i4vector2

  integer , pointer , dimension(:,:,:) :: ocld_io
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
  real(dp) , pointer , dimension(:,:,:) :: text2d_io
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

  real(dp) , pointer , dimension(:,:) :: sue_io
  real(dp) , pointer , dimension(:,:) :: sui_io
  real(dp) , pointer , dimension(:,:) :: nui_io
  real(dp) , pointer , dimension(:,:) :: nue_io
  real(dp) , pointer , dimension(:,:) :: sve_io
  real(dp) , pointer , dimension(:,:) :: svi_io
  real(dp) , pointer , dimension(:,:) :: nve_io
  real(dp) , pointer , dimension(:,:) :: nvi_io

  real(dp) , pointer , dimension(:,:) :: pptc_io
  real(dp) , pointer , dimension(:,:) :: pptnc_io

  real(dp) , pointer , dimension(:,:) :: cldefi_io
  real(dp) , pointer , dimension(:,:) :: hfx_io
                  
  real(dp) , pointer , dimension(:,:) :: zpbl_io

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
  interface deco1_scatter
    module procedure deco1_1d_real8_scatter ,   &
                     deco1_2d_real8_scatter ,   &
                     deco1_3d_real8_scatter ,   &
                     deco1_4d_real8_scatter ,   &
                     deco1_2d_real4_scatter ,   &
                     deco1_3d_real4_scatter ,   &
                     deco1_4d_real4_scatter ,   &
                     deco1_2d_integer_scatter , &
                     deco1_3d_integer_scatter , &
                     deco1_4d_integer_scatter
  end interface deco1_scatter

  interface deco1_gather
    module procedure deco1_1d_real8_gather ,   &
                     deco1_2d_real8_gather ,   &
                     deco1_3d_real8_gather ,   &
                     deco1_4d_real8_gather ,   &
                     deco1_2d_real4_gather ,   &
                     deco1_3d_real4_gather ,   &
                     deco1_4d_real4_gather ,   &
                     deco1_2d_integer_gather , &
                     deco1_3d_integer_gather , &
                     deco1_4d_integer_gather
  end interface deco1_gather

  interface subgrid_deco1_scatter
    module procedure subgrid_deco1_2d_real8_scatter ,   &
                     subgrid_deco1_3d_real8_scatter ,   &
                     subgrid_deco1_2d_real4_scatter ,   &
                     subgrid_deco1_3d_real4_scatter ,   &
                     subgrid_deco1_2d_integer_scatter , &
                     subgrid_deco1_3d_integer_scatter
  end interface subgrid_deco1_scatter

  interface subgrid_deco1_gather
    module procedure subgrid_deco1_2d_real8_gather ,   &
                     subgrid_deco1_3d_real8_gather ,   &
                     subgrid_deco1_2d_real4_gather ,   &
                     subgrid_deco1_3d_real4_gather ,   &
                     subgrid_deco1_2d_integer_gather , &
                     subgrid_deco1_3d_integer_gather
  end interface subgrid_deco1_gather

  interface deco1_exchange_left
    module procedure deco1_2d_real8_exchange_left ,  &
                     deco1_3d_real8_exchange_left,   &
                     deco1_4d_real8_exchange_left
  end interface deco1_exchange_left

  interface deco1_exchange_right
    module procedure deco1_2d_real8_exchange_right , &
                     deco1_3d_real8_exchange_right,  &
                     deco1_4d_real8_exchange_right
  end interface deco1_exchange_right

  public :: deco1_scatter , deco1_gather
  public :: subgrid_deco1_scatter , subgrid_deco1_gather
  public :: deco1_exchange_left , deco1_exchange_right
  public :: uvcross2dot , psc2psd

!---------- DATA init section--------------------------------------------

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
      call getmem3d(ocld_io,1,nnsg,jcross1,jcross2,icross1,icross2,'ocld_io')
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
      call getmem2d(kpbl_io,jcross1,jcross2,icross1,icross2,'kpbl_io')

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

      call getmem2d(sui_io,jdot1,jdot2,1,kz,'sue_io')
      call getmem2d(sue_io,jdot1,jdot2,1,kz,'sui_io')
      call getmem2d(nui_io,jdot1,jdot2,1,kz,'nui_io')
      call getmem2d(nue_io,jdot1,jdot2,1,kz,'nue_io')
      call getmem2d(sve_io,jdot1,jdot2,1,kz,'sve_io')
      call getmem2d(svi_io,jdot1,jdot2,1,kz,'svi_io')
      call getmem2d(nve_io,jdot1,jdot2,1,kz,'nve_io')
      call getmem2d(nvi_io,jdot1,jdot2,1,kz,'nvi_io')

      call getmem2d(pptc_io,jcross1,jcross2,icross1,icross2,'pptc_io')
      call getmem2d(pptnc_io,jcross1,jcross2,icross1,icross2,'pptnc_io')
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
  subroutine deco1_1d_real8_scatter(mg,ml)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:) , intent(out) :: ml ! model local
    integer :: ierr
    call mpi_scatter(mg,jxp,mpi_real8,ml,jxp,mpi_real8,0,mycomm,ierr)
  end subroutine deco1_1d_real8_scatter
!
  subroutine deco1_1d_real8_gather(ml,mg)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:) , intent(out) :: mg ! model global
    integer :: ierr
    call mpi_gather(ml,jxp,mpi_real8,mg,jxp,mpi_real8,0,mycomm,ierr)
  end subroutine deco1_1d_real8_gather
!
  subroutine deco1_2d_real8_scatter(mg,ml,j1,j2,i1,i2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ierr , ib , i , j , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx
    lsize = isize*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'deco1_2d_real8_scatter')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'deco1_2d_real8_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'deco1_2d_real8_scatter')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'deco1_2d_real8_scatter')
      end if
      ib = (j1-1)*isize+1
      do j = j1 , j2
        do i = i1 , i2
          r8vector1(ib) = mg(j,i)
          ib = ib + 1
        end do
      end do
    end if
    call mpi_scatter(r8vector1,lsize,mpi_real8, &
                     r8vector2,lsize,mpi_real8, &
                     0,mycomm,ierr)
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        ml(j,i) = r8vector2(ib)
        ib = ib + 1
      end do
    end do
  end subroutine deco1_2d_real8_scatter
!
  subroutine deco1_2d_real8_gather(ml,mg,j1,j2,i1,i2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ierr , ib , i , j , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx
    lsize = isize*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'deco1_2d_real8_gather')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'deco1_2d_real8_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        r8vector2(ib) = ml(j,i)
        ib = ib + 1
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'deco1_2d_real8_gather')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'deco1_2d_real8_gather')
      end if
    end if
    call mpi_gather(r8vector2,lsize,mpi_real8, &
                    r8vector1,lsize,mpi_real8, &
                    0,mycomm,ierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize+1
      do j = j1 , j2
        do i = i1 , i2
          mg(j,i) = r8vector1(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine deco1_2d_real8_gather
!
  subroutine deco1_3d_real8_scatter(mg,ml,j1,j2,i1,i2,k1,k2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ierr , ib , i , j , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx
    lsize = isize*ksize*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'deco1_3d_real8_scatter')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'deco1_3d_real8_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'deco1_3d_real8_scatter')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'deco1_3d_real8_scatter')
      end if
      ib = (j1-1)*isize*ksize+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            r8vector1(ib) = mg(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatter(r8vector1,lsize,mpi_real8, &
                     r8vector2,lsize,mpi_real8, &
                     0,mycomm,ierr)
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          ml(j,i,k) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end do
  end subroutine deco1_3d_real8_scatter
!
  subroutine deco1_3d_real8_gather(ml,mg,j1,j2,i1,i2,k1,k2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ierr , ib , i , j , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx
    lsize = isize*ksize*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'deco1_3d_real8_gather')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'deco1_3d_real8_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          r8vector2(ib) = ml(j,i,k)
          ib = ib + 1
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'deco1_3d_real8_gather')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'deco1_3d_real8_gather')
      end if
    end if
    call mpi_gather(r8vector2,lsize,mpi_real8, &
                    r8vector1,lsize,mpi_real8, &
                    0,mycomm,ierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            mg(j,i,k) = r8vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine deco1_3d_real8_gather
!
  subroutine deco1_4d_real8_scatter(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ierr , ib , i , j , k , n , js , je
    integer :: isize , ksize , nsize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    gsize = isize*ksize*nsize*jx
    lsize = isize*ksize*nsize*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'deco1_4d_real8_scatter')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'deco1_4d_real8_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'deco1_4d_real8_scatter')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'deco1_4d_real8_scatter')
      end if
      ib = (j1-1)*isize*ksize*nsize+1
      do j = j1 , j2
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              r8vector1(ib) = mg(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_scatter(r8vector1,lsize,mpi_real8, &
                     r8vector2,lsize,mpi_real8, &
                     0,mycomm,ierr)
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize*nsize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            ml(j,i,k,n) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end do
  end subroutine deco1_4d_real8_scatter
!
  subroutine deco1_4d_real8_gather(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ierr , ib , i , j , k , n , js , je
    integer :: isize , ksize , nsize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    gsize = isize*ksize*nsize*jx
    lsize = isize*ksize*nsize*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'deco1_4d_real8_gather')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'deco1_4d_real8_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize*nsize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            r8vector2(ib) = ml(j,i,k,n)
            ib = ib + 1
          end do
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'deco1_4d_real8_gather')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'deco1_4d_real8_gather')
      end if
    end if
    call mpi_gather(r8vector2,lsize,mpi_real8, &
                    r8vector1,lsize,mpi_real8, &
                    0,mycomm,ierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize*nsize+1
      do j = j1 , j2
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              mg(j,i,k,n) = r8vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine deco1_4d_real8_gather
!
  subroutine deco1_2d_real4_scatter(mg,ml,j1,j2,i1,i2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(sp) , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    real(sp) , pointer , dimension(:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ierr , ib , i , j , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx
    lsize = isize*jxp
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'deco1_2d_real4_scatter')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'deco1_2d_real4_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'deco1_2d_real4_scatter')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'deco1_2d_real4_scatter')
      end if
      ib = (j1-1)*isize+1
      do j = j1 , j2
        do i = i1 , i2
          r4vector1(ib) = mg(j,i)
          ib = ib + 1
        end do
      end do
    end if
    call mpi_scatter(r4vector1,lsize,mpi_real4, &
                     r4vector2,lsize,mpi_real4, &
                     0,mycomm,ierr)
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        ml(j,i) = r4vector2(ib)
        ib = ib + 1
      end do
    end do
  end subroutine deco1_2d_real4_scatter
!
  subroutine deco1_2d_real4_gather(ml,mg,j1,j2,i1,i2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(sp) , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    real(sp) , pointer , dimension(:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ierr , ib , i , j , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx
    lsize = isize*jxp
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'deco1_2d_real4_gather')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'deco1_2d_real4_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        r4vector2(ib) = ml(j,i)
        ib = ib + 1
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'deco1_2d_real4_gather')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'deco1_2d_real4_gather')
      end if
    end if
    call mpi_gather(r4vector2,lsize,mpi_real4, &
                    r4vector1,lsize,mpi_real4, &
                    0,mycomm,ierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize+1
      do j = j1 , j2
        do i = i1 , i2
          mg(j,i) = r4vector1(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine deco1_2d_real4_gather
!
  subroutine deco1_3d_real4_scatter(mg,ml,j1,j2,i1,i2,k1,k2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(sp) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ierr , ib , i , j , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx
    lsize = isize*ksize*jxp
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'deco1_3d_real4_scatter')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'deco1_3d_real4_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'deco1_3d_real4_scatter')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'deco1_3d_real4_scatter')
      end if
      ib = (j1-1)*isize*ksize+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            r4vector1(ib) = mg(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatter(r4vector1,lsize,mpi_real4, &
                     r4vector2,lsize,mpi_real4, &
                     0,mycomm,ierr)
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          ml(j,i,k) = r4vector2(ib)
          ib = ib + 1
        end do
      end do
    end do
  end subroutine deco1_3d_real4_scatter
!
  subroutine deco1_3d_real4_gather(ml,mg,j1,j2,i1,i2,k1,k2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(sp) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ierr , ib , i , j , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx
    lsize = isize*ksize*jxp
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'deco1_3d_real4_gather')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'deco1_3d_real4_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          r4vector2(ib) = ml(j,i,k)
          ib = ib + 1
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'deco1_3d_real4_gather')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'deco1_3d_real4_gather')
      end if
    end if
    call mpi_gather(r4vector2,lsize,mpi_real4, &
                    r4vector1,lsize,mpi_real4, &
                    0,mycomm,ierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            mg(j,i,k) = r4vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine deco1_3d_real4_gather
!
  subroutine deco1_4d_real4_scatter(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(sp) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ierr , ib , i , j , k , n , js , je
    integer :: isize , ksize , nsize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    gsize = isize*ksize*nsize*jx
    lsize = isize*ksize*nsize*jxp
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'deco1_4d_real4_scatter')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'deco1_4d_real4_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'deco1_4d_real4_scatter')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'deco1_4d_real4_scatter')
      end if
      ib = (j1-1)*isize*ksize*nsize+1
      do j = j1 , j2
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              r4vector1(ib) = mg(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_scatter(r4vector1,lsize,mpi_real4, &
                     r4vector2,lsize,mpi_real4, &
                     0,mycomm,ierr)
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize*nsize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            ml(j,i,k,n) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end do
  end subroutine deco1_4d_real4_scatter
!
  subroutine deco1_4d_real4_gather(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(sp) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ierr , ib , i , j , k , n , js , je
    integer :: isize , ksize , nsize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    gsize = isize*ksize*nsize*jx
    lsize = isize*ksize*nsize*jxp
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'deco1_4d_real4_gather')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'deco1_4d_real4_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize*nsize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            r4vector2(ib) = ml(j,i,k,n)
            ib = ib + 1
          end do
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'deco1_4d_real4_gather')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'deco1_4d_real4_gather')
      end if
    end if
    call mpi_gather(r4vector2,lsize,mpi_real4, &
                    r4vector1,lsize,mpi_real4, &
                    0,mycomm,ierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize*nsize+1
      do j = j1 , j2
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              mg(j,i,k,n) = r4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine deco1_4d_real4_gather
!
  subroutine deco1_2d_integer_scatter(mg,ml,j1,j2,i1,i2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    integer , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    integer , pointer , dimension(:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ierr , ib , i , j , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx
    lsize = isize*jxp
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'deco1_2d_integer_scatter')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'deco1_2d_integer_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'deco1_2d_integer_scatter')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'deco1_2d_integer_scatter')
      end if
      ib = (j1-1)*isize+1
      do j = j1 , j2
        do i = i1 , i2
          i4vector1(ib) = mg(j,i)
          ib = ib + 1
        end do
      end do
    end if
    call mpi_scatter(i4vector1,lsize,mpi_integer, &
                     i4vector2,lsize,mpi_integer, &
                     0,mycomm,ierr)
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        ml(j,i) = i4vector2(ib)
        ib = ib + 1
      end do
    end do
  end subroutine deco1_2d_integer_scatter
!
  subroutine deco1_2d_integer_gather(ml,mg,j1,j2,i1,i2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    integer , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    integer , pointer , dimension(:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ierr , ib , i , j , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx
    lsize = isize*jxp
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'deco1_2d_integer_gather')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'deco1_2d_integer_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        i4vector2(ib) = ml(j,i)
        ib = ib + 1
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'deco1_2d_integer_gather')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'deco1_2d_integer_gather')
      end if
    end if
    call mpi_gather(i4vector2,lsize,mpi_integer, &
                    i4vector1,lsize,mpi_integer, &
                    0,mycomm,ierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize+1
      do j = j1 , j2
        do i = i1 , i2
          mg(j,i) = i4vector1(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine deco1_2d_integer_gather
!
  subroutine deco1_3d_integer_scatter(mg,ml,j1,j2,i1,i2,k1,k2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    integer , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    integer , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ierr , ib , i , j , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx
    lsize = isize*ksize*jxp
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'deco1_3d_integer_scatter')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'deco1_3d_integer_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'deco1_3d_integer_scatter')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'deco1_3d_integer_scatter')
      end if
      ib = (j1-1)*isize*ksize+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            i4vector1(ib) = mg(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatter(i4vector1,lsize,mpi_integer, &
                     i4vector2,lsize,mpi_integer, &
                     0,mycomm,ierr)
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          ml(j,i,k) = i4vector2(ib)
          ib = ib + 1
        end do
      end do
    end do
  end subroutine deco1_3d_integer_scatter
!
  subroutine deco1_3d_integer_gather(ml,mg,j1,j2,i1,i2,k1,k2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    integer , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    integer , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ierr , ib , i , j , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx
    lsize = isize*ksize*jxp
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'deco1_3d_integer_gather')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'deco1_3d_integer_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          i4vector2(ib) = ml(j,i,k)
          ib = ib + 1
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'deco1_3d_integer_gather')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'deco1_3d_integer_gather')
      end if
    end if
    call mpi_gather(i4vector2,lsize,mpi_integer, &
                    i4vector1,lsize,mpi_integer, &
                    0,mycomm,ierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            mg(j,i,k) = i4vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine deco1_3d_integer_gather
!
  subroutine deco1_4d_integer_scatter(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    integer , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    integer , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ierr , ib , i , j , k , n , js , je
    integer :: isize , ksize , nsize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    gsize = isize*ksize*nsize*jx
    lsize = isize*ksize*nsize*jxp
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'deco1_4d_integer_scatter')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'deco1_4d_integer_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'deco1_4d_integer_scatter')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'deco1_4d_integer_scatter')
      end if
      ib = (j1-1)*isize*ksize*nsize+1
      do j = j1 , j2
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              i4vector1(ib) = mg(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_scatter(i4vector1,lsize,mpi_integer, &
                     i4vector2,lsize,mpi_integer, &
                     0,mycomm,ierr)
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize*nsize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            ml(j,i,k,n) = i4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end do
  end subroutine deco1_4d_integer_scatter
!
  subroutine deco1_4d_integer_gather(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    integer , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    integer , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ierr , ib , i , j , k , n , js , je
    integer :: isize , ksize , nsize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    gsize = isize*ksize*nsize*jx
    lsize = isize*ksize*nsize*jxp
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'deco1_4d_integer_gather')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'deco1_4d_integer_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize*nsize+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            i4vector2(ib) = ml(j,i,k,n)
            ib = ib + 1
          end do
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'deco1_4d_integer_gather')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'deco1_4d_integer_gather')
      end if
    end if
    call mpi_gather(i4vector2,lsize,mpi_integer, &
                    i4vector1,lsize,mpi_integer, &
                    0,mycomm,ierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize*nsize+1
      do j = j1 , j2
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              mg(j,i,k,n) = i4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine deco1_4d_integer_gather
!
  subroutine subgrid_deco1_2d_real8_scatter(mg,ml,j1,j2,i1,i2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ierr , ib , i , j , nn , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*nnsg*jx
    lsize = isize*nnsg*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'subgrid_deco1_2d_real8_scatter')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'subgrid_deco1_2d_real8_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'subgrid_deco1_2d_real8_scatter')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'subgrid_deco1_2d_real8_scatter')
      end if
      ib = (j1-1)*isize*nnsg+1
      do j = j1 , j2
        do i = i1 , i2
          do nn = 1 , nnsg
            r8vector1(ib) = mg(nn,j,i)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatter(r8vector1,lsize,mpi_real8, &
                     r8vector2,lsize,mpi_real8, &
                     0,mycomm,ierr)
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        do nn = 1 , nnsg
          ml(nn,j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end do
  end subroutine subgrid_deco1_2d_real8_scatter
!
  subroutine subgrid_deco1_2d_real8_gather(ml,mg,j1,j2,i1,i2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ierr , ib , i , j , nn , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*nnsg*jx
    lsize = isize*nnsg*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'subgrid_deco1_2d_real8_gather')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'subgrid_deco1_2d_real8_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        do nn = 1 , nnsg
          r8vector2(ib) = ml(nn,j,i)
          ib = ib + 1
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'subgrid_deco1_2d_real8_gather')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'subgrid_deco1_2d_real8_gather')
      end if
    end if
    call mpi_gather(r8vector2,lsize,mpi_real8, &
                    r8vector1,lsize,mpi_real8, &
                    0,mycomm,ierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*nnsg+1
      do j = j1 , j2
        do i = i1 , i2
          do nn = 1 , nnsg
            mg(nn,j,i) = r8vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine subgrid_deco1_2d_real8_gather
!
  subroutine subgrid_deco1_3d_real8_scatter(mg,ml,j1,j2,i1,i2,k1,k2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ierr , ib , i , j , nn , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*nnsg*jx
    lsize = isize*ksize*nnsg*jxp
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'subgrid_deco1_3d_real8_scatter')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'subgrid_deco1_3d_real8_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'subgrid_deco1_3d_real8_scatter')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'subgrid_deco1_3d_real8_scatter')
      end if
      ib = (j1-1)*isize*ksize*nnsg+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            do nn = 1 , nnsg
              r8vector1(ib) = mg(nn,j,i,k)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_scatter(r8vector1,lsize,mpi_real8, &
                     r8vector2,lsize,mpi_real8, &
                     0,mycomm,ierr)
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          do nn = 1 , nnsg
            ml(nn,j,i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end do
  end subroutine subgrid_deco1_3d_real8_scatter
!
  subroutine subgrid_deco1_3d_real8_gather(ml,mg,j1,j2,i1,i2,k1,k2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ierr , ib , i , j , k , nn , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx*nnsg
    lsize = isize*ksize*jxp*nnsg
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,lsize,'subgrid_deco1_3d_real8_gather')
    else if ( size(r8vector2) < lsize ) then
      call getmem1d(r8vector2,1,lsize,'subgrid_deco1_3d_real8_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          do nn = 1 , nnsg
            r8vector2(ib) = ml(nn,j,i,k)
            ib = ib + 1
          end do
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r8vector1) ) then
        call getmem1d(r8vector1,1,gsize,'subgrid_deco1_3d_real8_gather')
      else if ( size(r8vector1) < gsize ) then
        call getmem1d(r8vector1,1,gsize,'subgrid_deco1_3d_real8_gather')
      end if
    end if
    call mpi_gather(r8vector2,lsize,mpi_real8, &
                    r8vector1,lsize,mpi_real8, &
                    0,mycomm,ierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize*nnsg+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            do nn = 1 , nnsg
              mg(nn,j,i,k) = r8vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine subgrid_deco1_3d_real8_gather
!
  subroutine subgrid_deco1_2d_real4_scatter(mg,ml,j1,j2,i1,i2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(sp) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ierr , ib , i , j , nn , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx*nnsg
    lsize = isize*jxp*nnsg
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'subgrid_deco1_2d_real4_scatter')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'subgrid_deco1_2d_real4_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'subgrid_deco1_2d_real4_scatter')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'subgrid_deco1_2d_real4_scatter')
      end if
      ib = (j1-1)*isize*nnsg+1
      do j = j1 , j2
        do i = i1 , i2
          do nn = 1 , nnsg
            r4vector1(ib) = mg(nn,j,i)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatter(r4vector1,lsize,mpi_real4, &
                     r4vector2,lsize,mpi_real4, &
                     0,mycomm,ierr)
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        do nn = 1 , nnsg
          ml(nn,j,i) = r4vector2(ib)
          ib = ib + 1
        end do
      end do
    end do
  end subroutine subgrid_deco1_2d_real4_scatter
!
  subroutine subgrid_deco1_2d_real4_gather(ml,mg,j1,j2,i1,i2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(sp) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(sp) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ierr , ib , i , j , nn , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx*nnsg
    lsize = isize*jxp*nnsg
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'subgrid_deco1_2d_real4_gather')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'subgrid_deco1_2d_real4_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        do nn = 1 , nnsg
          r4vector2(ib) = ml(nn,j,i)
          ib = ib + 1
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'subgrid_deco1_2d_real4_gather')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'subgrid_deco1_2d_real4_gather')
      end if
    end if
    call mpi_gather(r4vector2,lsize,mpi_real4, &
                    r4vector1,lsize,mpi_real4, &
                    0,mycomm,ierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*nnsg+1
      do j = j1 , j2
        do i = i1 , i2
          do nn = 1 , nnsg
            mg(nn,j,i) = r4vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine subgrid_deco1_2d_real4_gather
!
  subroutine subgrid_deco1_3d_real4_scatter(mg,ml,j1,j2,i1,i2,k1,k2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(sp) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ierr , ib , i , j , nn , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx*nnsg
    lsize = isize*ksize*jxp*nnsg
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'subgrid_deco1_3d_real4_scatter')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'subgrid_deco1_3d_real4_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'subgrid_deco1_3d_real4_scatter')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'subgrid_deco1_3d_real4_scatter')
      end if
      ib = (j1-1)*isize*ksize*nnsg+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            do nn = 1 , nnsg
              r4vector1(ib) = mg(nn,j,i,k)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_scatter(r4vector1,lsize,mpi_real4, &
                     r4vector2,lsize,mpi_real4, &
                     0,mycomm,ierr)
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          do nn = 1 , nnsg
            ml(nn,j,i,k) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end do
  end subroutine subgrid_deco1_3d_real4_scatter
!
  subroutine subgrid_deco1_3d_real4_gather(ml,mg,j1,j2,i1,i2,k1,k2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(sp) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(sp) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ierr , ib , i , j , k , nn , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx*nnsg
    lsize = isize*ksize*jxp*nnsg
    if ( .not. associated(r4vector2) ) then
      call getmem1d(r4vector2,1,lsize,'subgrid_deco1_3d_real4_gather')
    else if ( size(r4vector2) < lsize ) then
      call getmem1d(r4vector2,1,lsize,'subgrid_deco1_3d_real4_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          do nn = 1 , nnsg
            r4vector2(ib) = ml(nn,j,i,k)
            ib = ib + 1
          end do
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(r4vector1) ) then
        call getmem1d(r4vector1,1,gsize,'subgrid_deco1_3d_real4_gather')
      else if ( size(r4vector1) < gsize ) then
        call getmem1d(r4vector1,1,gsize,'subgrid_deco1_3d_real4_gather')
      end if
    end if
    call mpi_gather(r4vector2,lsize,mpi_real4, &
                    r4vector1,lsize,mpi_real4, &
                    0,mycomm,ierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize*nnsg+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            do nn = 1 , nnsg
              mg(nn,j,i,k) = r4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine subgrid_deco1_3d_real4_gather
!
  subroutine subgrid_deco1_2d_integer_scatter(mg,ml,j1,j2,i1,i2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    integer , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    integer , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ierr , ib , i , j , nn , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx*nnsg
    lsize = isize*jxp*nnsg
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'subgrid_deco1_2d_integer_scatter')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'subgrid_deco1_2d_integer_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'subgrid_deco1_2d_integer_scatter')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'subgrid_deco1_2d_integer_scatter')
      end if
      ib = (j1-1)*isize*nnsg+1
      do j = j1 , j2
        do i = i1 , i2
          do nn = 1 , nnsg
            i4vector1(ib) = mg(nn,j,i)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatter(i4vector1,lsize,mpi_integer, &
                     i4vector2,lsize,mpi_integer, &
                     0,mycomm,ierr)
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        do nn = 1 , nnsg
          ml(nn,j,i) = i4vector2(ib)
          ib = ib + 1
        end do
      end do
    end do
  end subroutine subgrid_deco1_2d_integer_scatter
!
  subroutine subgrid_deco1_2d_integer_gather(ml,mg,j1,j2,i1,i2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    integer , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    integer , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2
    integer :: ierr , ib , i , j , nn , isize , gsize , lsize , js , je
    isize = i2-i1+1
    gsize = isize*jx*nnsg
    lsize = isize*jxp*nnsg
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'subgrid_deco1_2d_integer_gather')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'subgrid_deco1_2d_integer_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do i = i1 , i2
        do nn = 1 , nnsg
          i4vector2(ib) = ml(nn,j,i)
          ib = ib + 1
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'subgrid_deco1_2d_integer_gather')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'subgrid_deco1_2d_integer_gather')
      end if
    end if
    call mpi_gather(i4vector2,lsize,mpi_integer, &
                    i4vector1,lsize,mpi_integer, &
                    0,mycomm,ierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*nnsg+1
      do j = j1 , j2
        do i = i1 , i2
          do nn = 1 , nnsg
            mg(nn,j,i) = i4vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine subgrid_deco1_2d_integer_gather
!
  subroutine subgrid_deco1_3d_integer_scatter(mg,ml,j1,j2,i1,i2,k1,k2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    integer , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    integer , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ierr , ib , i , j , nn , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx*nnsg
    lsize = isize*ksize*jxp*nnsg
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'subgrid_deco1_3d_integer_scatter')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'subgrid_deco1_3d_integer_scatter')
    end if
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'subgrid_deco1_3d_integer_scatter')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'subgrid_deco1_3d_integer_scatter')
      end if
      ib = (j1-1)*isize*ksize*nnsg+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            do nn = 1 , nnsg
              i4vector1(ib) = mg(nn,j,i,k)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_scatter(i4vector1,lsize,mpi_integer, &
                     i4vector2,lsize,mpi_integer, &
                     0,mycomm,ierr)
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          do nn = 1 , nnsg
            ml(nn,j,i,k) = i4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end do
  end subroutine subgrid_deco1_3d_integer_scatter
!
  subroutine subgrid_deco1_3d_integer_gather(ml,mg,j1,j2,i1,i2,k1,k2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    integer , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    integer , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer :: ierr , ib , i , j , nn , k , js , je
    integer :: isize , ksize , gsize , lsize
    isize = i2-i1+1
    ksize = k2-k1+1
    gsize = isize*ksize*jx*nnsg
    lsize = isize*ksize*jxp*nnsg
    if ( .not. associated(i4vector2) ) then
      call getmem1d(i4vector2,1,lsize,'subgrid_deco1_3d_integer_gather')
    else if ( size(i4vector2) < lsize ) then
      call getmem1d(i4vector2,1,lsize,'subgrid_deco1_3d_integer_gather')
    end if
    if ( myid == 0 ) then
      js = j1
      je = jxp
      ib = (j1-1)*isize*ksize*nnsg+1
    else if ( myid == nproc-1 ) then
      js = 1
      je = jxp-(jx-j2)
      ib = 1
    else
      js = 1
      je = jxp
      ib = 1
    end if
    do j = js , je
      do k = k1 , k2
        do i = i1 , i2
          do nn = 1 , nnsg
            i4vector2(ib) = ml(nn,j,i,k)
            ib = ib + 1
          end do
        end do
      end do
    end do
    if ( myid == 0 ) then
      if ( .not. associated(i4vector1) ) then
        call getmem1d(i4vector1,1,gsize,'subgrid_deco1_3d_integer_gather')
      else if ( size(i4vector1) < gsize ) then
        call getmem1d(i4vector1,1,gsize,'subgrid_deco1_3d_integer_gather')
      end if
    end if
    call mpi_gather(i4vector2,lsize,mpi_integer, &
                    i4vector1,lsize,mpi_integer, &
                    0,mycomm,ierr)
    if ( myid == 0 ) then
      ib = (j1-1)*isize*ksize*nnsg+1
      do j = j1 , j2
        do k = k1 , k2
          do i = i1 , i2
            do nn = 1 , nnsg
              mg(nn,j,i,k) = i4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine subgrid_deco1_3d_integer_gather
!
  subroutine deco1_2d_real8_exchange_right(ml,nex,i1,i2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:,:) , intent(inout) :: ml
    integer , intent(in) :: nex , i1 , i2
    integer :: isize , ssize , i , j , ib
    integer :: ierr
    isize = i2-i1+1
    ssize = nex*isize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'deco1_2d_real8_exchange_right')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'deco1_2d_real8_exchange_right')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'deco1_2d_real8_exchange_right')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'deco1_2d_real8_exchange_right')
    end if
    if ( iwest /= mpi_proc_null ) then
      ib = 1
      do j = 1 , nex
        do i = i1 , i2
          r8vector1(ib) = ml(j,i)
          ib = ib + 1
        end do
      end do
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,iwest,2, &
                      r8vector2,ssize,mpi_real8,ieast,2, &
                      mycomm,mpi_status_ignore,ierr)
    if ( ieast /= mpi_proc_null ) then
      ib = 1
      do j = 1 , nex
        do i = i1 , i2
          ml(jxp+j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine deco1_2d_real8_exchange_right
!
  subroutine deco1_2d_real8_exchange_left(ml,nex,i1,i2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:,:) , intent(inout) :: ml
    integer , intent(in) :: nex , i1 , i2
    integer :: isize , ssize , j , i , ib
    integer :: ierr
    isize = i2-i1+1
    ssize = nex*isize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'deco1_2d_real8_exchange_left')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'deco1_2d_real8_exchange_left')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'deco1_2d_real8_exchange_left')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'deco1_2d_real8_exchange_left')
    end if
    if ( ieast /= mpi_proc_null ) then
      ib = 1 
      do j = 1 , nex
        do i = i1 , i2
          r8vector1(ib) = ml((jxp-j)+1,i)
          ib = ib + 1
        end do
      end do
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,ieast,1, &
                      r8vector2,ssize,mpi_real8,iwest,1, &
                      mycomm,mpi_status_ignore,ierr)
    if ( iwest /= mpi_proc_null ) then
      ib = 1
      do j = 1 , nex
        do i = i1 , i2
          ml((1-j),i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine deco1_2d_real8_exchange_left
!
  subroutine deco1_3d_real8_exchange_right(ml,nex,i1,i2,k1,k2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer , intent(in) :: nex , i1 , i2 , k1 , k2
    integer :: isize , ksize , ssize , hsize , i , j , k , ib
    integer :: ierr
    isize = i2-i1+1
    ksize = k2-k1+1
    hsize = isize*ksize
    ssize = nex*hsize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'deco1_3d_real8_exchange_right')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'deco1_3d_real8_exchange_right')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'deco1_3d_real8_exchange_right')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'deco1_3d_real8_exchange_right')
    end if
    if ( iwest /= mpi_proc_null ) then
      ib = 1
      do j = 1 , nex
        do k = k1 , k2
          do i = i1 , i2
            r8vector1(ib) = ml(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,iwest,2, &
                      r8vector2,ssize,mpi_real8,ieast,2, &
                      mycomm,mpi_status_ignore,ierr)
    if ( ieast /= mpi_proc_null ) then
      ib = 1
      do j = 1 , nex
        do k = k1 , k2
          do i = i1 , i2
            ml(jxp+j,i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine deco1_3d_real8_exchange_right
!
  subroutine deco1_3d_real8_exchange_left(ml,nex,i1,i2,k1,k2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer , intent(in) :: nex , i1 , i2 , k1 , k2
    integer :: isize , ksize , ssize , hsize , i , j , k , ib
    integer :: ierr
    isize = i2-i1+1
    ksize = k2-k1+1
    hsize = isize*ksize
    ssize = nex*hsize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'deco1_3d_real8_exchange_left')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'deco1_3d_real8_exchange_left')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'deco1_3d_real8_exchange_left')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'deco1_3d_real8_exchange_left')
    end if
    if ( ieast /= mpi_proc_null ) then
      ib = 1
      do j = 1 , nex
        do k = k1 , k2
          do i = i1 , i2
            r8vector1(ib) = ml((jxp-j)+1,i,k)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,ieast,1, &
                      r8vector2,ssize,mpi_real8,iwest,1, &
                      mycomm,mpi_status_ignore,ierr)
    if ( iwest /= mpi_proc_null ) then
      ib = 1
      do j = 1 , nex
        do k = k1 , k2
          do i = i1 , i2
            ml((1-j),i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine deco1_3d_real8_exchange_left
!
  subroutine deco1_4d_real8_exchange_right(ml,nex,i1,i2,k1,k2,n1,n2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer , intent(in) :: nex , i1 , i2 , k1 , k2 , n1 , n2
    integer :: isize , ssize , ksize , nsize , vsize , hsize , ib
    integer :: i , j , k , n
    integer :: ierr
    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    vsize = isize*ksize
    hsize = vsize*nsize
    ssize = nex*hsize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'deco1_4d_real8_exchange_right')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'deco1_4d_real8_exchange_right')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'deco1_4d_real8_exchange_right')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'deco1_4d_real8_exchange_right')
    end if
    if ( iwest /= mpi_proc_null ) then
      ib = 1
      do j = 1 , nex
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              r8vector1(ib) = ml(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,iwest,2, &
                      r8vector2,ssize,mpi_real8,ieast,2, &
                      mycomm,mpi_status_ignore,ierr)
    if ( ieast /= mpi_proc_null ) then
      ib = 1
      do j = 1 , nex
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              ml(jxp+j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine deco1_4d_real8_exchange_right
!
  subroutine deco1_4d_real8_exchange_left(ml,nex,i1,i2,k1,k2,n1,n2)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer , intent(in) :: nex , i1 , i2 , k1 , k2 , n1 , n2
    integer :: isize , ssize , ksize , nsize , vsize , hsize , ib
    integer :: i , j , k , n
    integer :: ierr
    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    vsize = isize*ksize
    hsize = vsize*nsize
    ssize = nex*hsize
    if ( .not. associated(r8vector1) ) then
      call getmem1d(r8vector1,1,ssize,'deco1_4d_real8_exchange_left')
    else if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'deco1_4d_real8_exchange_left')
    end if
    if ( .not. associated(r8vector2) ) then
      call getmem1d(r8vector2,1,ssize,'deco1_4d_real8_exchange_left')
    else if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'deco1_4d_real8_exchange_left')
    end if
    if ( ieast /= mpi_proc_null ) then
      ib = 1
      do j = 1 , nex
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              r8vector1(ib) = ml((jxp-j)+1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_sendrecv(r8vector1,ssize,mpi_real8,ieast,1, &
                      r8vector2,ssize,mpi_real8,iwest,1, &
                      mycomm,mpi_status_ignore,ierr)
    if ( iwest /= mpi_proc_null ) then
      ib = 1
      do j = 1 , nex
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              ml((1-j),i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine deco1_4d_real8_exchange_left
!
! Takes u and v on the cross grid (the same grid as t, qv, qc, etc.)
! and interpolates the u and v to the dot grid.
! This routine sheilds the user of the function from the need to worry
! about the details of the domain decomposition.  
!
! Written by Travis A. O'Brien 01/04/11.
!
  subroutine uvcross2dot(ux,vx,ud,vd)
    implicit none
    real(dp) , pointer , dimension(:,:,:) , intent(inout) :: ux , vx
    real(dp) , pointer , dimension(:,:,:) , intent(inout) :: ud , vd
    integer :: ib , ie , jb , je , i , j

    ! TODO:  It might make sense to encapsulate the following code
    ! in to a standard routine, since this boundary sending code is
    ! ubiquitous throughout the RegCM code and it is domain
    ! decomposition-dependent.

    ! Send the right-edge of the u/v tendencies to the left
    ! edge of the next process's u/v tendencies (so that
    ! invar%u(i,k,0) holds invar%u(i,k,jxp) of the parallel
    ! chunk next door)

    call deco1_exchange_left(ux,1,1,iy,1,kz)
    call deco1_exchange_left(vx,1,1,iy,1,kz)

    ! Set j-loop boundaries
    jb = jbegin
    je = jendx
    ! Set i-loop boundaries
    ib = 2
    ie = iym1

    !
    !     x     x     x     x     x     x
    !
    !        o     o     o     o     o 
    !         (i-1,j-1)     (i,j-1)            
    !     x     x     x-----x     x     x
    !                 |(i,j)|
    !        o     o  |  o  |  o     o 
    !                 |     |
    !     x     x     x-----x     x     x
    !           (i-1,j)     (i,j)
    !
    !        o     o     o     o     o 
    !
    !     x     x     x     x     x     x
    !

    ! Perform the bilinear interpolation necessary
    ! to put the u and v variables on the dot grid.

    do i = ib , ie
      do j = jb , je
        ud(j,i,:) =  ud(j,i,:) +             &
          d_rfour*(ux(j,i,:) + ux(j-1,i,:) +   &
                   ux(j,i-1,:) + ux(j-1,i-1,:))
        vd(j,i,:) =  vd(j,i,:) +             &
          d_rfour*(vx(j,i,:) + vx(j-1,i,:) +   &
                   vx(j,i-1,:) + vx(j-1,i-1,:))
      end do
    end do
  end subroutine uvcross2dot
!
  subroutine psc2psd(pc,pd)
    implicit none
    real(dp) , pointer , dimension(:,:) , intent(in)  :: pc
    real(dp) , pointer , dimension(:,:) , intent(out) :: pd
    integer :: i , j

    do i = idi1 , idi2
      do j = jdi1 , jdi2
        pd(j,i) = (pc(j,i)+pc(j,i-1)+pc(j-1,i)+pc(j-1,i-1))*d_rfour
      end do
    end do
    if ( .not. ma%bandflag ) then
      if ( ma%hasleft ) then
        do i = idi1 , idi2
          pd(jde1,i) = (pc(jce1,i)+pc(jce1,i-1))*d_half
        end do
        pd(jde1,ide1) = pc(jce1,ice1)
        pd(jde1,ide2) = pc(jce1,ice2)
      end if
      if ( ma%hasright ) then
        do i = idi1 , idi2
          pd(jde2,i) = (pc(jce2,i)+pc(jce2,i-1))*d_half
        end do
        pd(jde2,ide1) = pc(jce2,ice1)
        pd(jde2,ide2) = pc(jce2,ice2)
      end if
      if ( ma%hasbottom ) then
        do j = jdi1 , jdi2
          pd(j,ide1)  = (pc(j,ice1)+pc(j-1,ice1))*d_half
        end do
      end if
      if ( ma%hastop ) then
        do j = jdi1 , jdi2
          pd(j,ide2) = (pc(j,ice2)+pc(j-1,ice2))*d_half
        end do
      end if
    else
      if ( ma%hasbottom ) then
        do j = jde1 , jde2
          pd(j,ide1)  = (pc(j,ice1)+pc(j-1,ice1))*d_half
        end do
      end if
      if ( ma%hastop ) then
        do j = jde1 , jde2
          pd(j,ide2) = (pc(j,ice2)+pc(j-1,ice2))*d_half
        end do
      end if
    end if
  end subroutine psc2psd

end module mod_mppio
