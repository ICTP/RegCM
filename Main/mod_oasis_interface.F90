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

#ifdef OASIS

module mod_oasis_interface

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_stdio
  use mod_message
  use mod_service
  use mod_mppparam
  use mod_runparams , only : isocean , dtsec !, rcmtimer
  use mod_bats_common , only : rdnnsg
  use mod_atm_interface , only : atms , sfs , flwd , flw , fsw , sinc , mddom
  use mod_lm_interface , only : lms

  use mod_oasis
  use mod_oasis_params
  use mod_oasis_signature
  use mod_oasis_generic

  implicit none

  private

  public :: comp_name , comp_id ! -> mod_oasis_params
  public :: oasis_lag           ! -> mod_oasis_params

  !--------------------------------------------------------------------
  ! for grid and partition activation
  ! (dot/cross , external/internal)
  logical :: l_make_grdde , l_make_grddi , &
             l_make_grdce , l_make_grdci
  type(infogrd) , target , allocatable :: grdde , grddi , &
                                          grdce , grdci

  !--------------------------------------------------------------------
  ! below the namelist parameters (namelist oasisparam)
  public :: l_write_grids , write_restart_option ! -> mod_oasis_params
  public :: oasis_sync_lag                       ! -> mod_oasis_params
  logical , public :: l_cpl_im_sst , &
!                      l_cpl_im_sit , & ! not coded
                      l_cpl_im_wz0 , & ! not tested
                      l_cpl_im_wust , & ! not tested
                      l_cpl_ex_u10m , & ! not tested
                      l_cpl_ex_v10m , & ! not tested
                      l_cpl_ex_wspd , & ! not tested
                      l_cpl_ex_wdir , & ! not tested
                      l_cpl_ex_t2m , & ! not tested
!                      l_cpl_ex_t10m , & ! not coded
                      l_cpl_ex_q2m , & ! not tested
!                      l_cpl_ex_q10m , & ! not coded
                      l_cpl_ex_slp , &
                      l_cpl_ex_taux , &
                      l_cpl_ex_tauy , &
                      l_cpl_ex_z0 , & ! not tested
                      l_cpl_ex_ustr , & ! not tested
                      l_cpl_ex_evap , &
                      l_cpl_ex_prec , &
                      l_cpl_ex_nuwa , &
                      l_cpl_ex_ulhf , &
                      l_cpl_ex_ushf , &
                      l_cpl_ex_uwlw , &
                      l_cpl_ex_dwlw , &
                      l_cpl_ex_nulw , &
                      l_cpl_ex_uwsw , &
                      l_cpl_ex_dwsw , &
                      l_cpl_ex_ndsw , &
                      l_cpl_ex_rhoa ! not tested
  ! OASIS field +++

  !--------------------------------------------------------------------
  ! below the variables for the fields' information
  type(infofld) , allocatable :: im_sst , &
!                                 im_sit , &
                                 im_wz0 , &
                                 im_wust , &
                                 ex_u10m , &
                                 ex_v10m , &
                                 ex_wspd , &
                                 ex_wdir , &
                                 ex_t2m , &
!                                 ex_t10m , &
                                 ex_q2m , &
!                                 ex_q10m , &
                                 ex_slp , &
                                 ex_taux , &
                                 ex_tauy , &
                                 ex_z0 , &
                                 ex_ustr , &
                                 ex_evap , &
                                 ex_prec , &
                                 ex_nuwa , &
                                 ex_ulhf , &
                                 ex_ushf , &
                                 ex_uwlw , &
                                 ex_dwlw , &
                                 ex_nulw , &
                                 ex_uwsw , &
                                 ex_dwsw , &
                                 ex_ndsw , &
                                 ex_rhoa
  ! OASIS field +++

  ! OASIS needs double precision arrays. These variables are optional.
  ! Required for the fields to import; recommended for the fields to
  ! export which need a bit of work from the model variables.
  real(rkx) , dimension(:,:) , allocatable :: cpl_sst , &
!                                              cpl_sit , &
                                              cpl_wz0 , &
                                              cpl_wust , &
                                              cpl_wdir
  ! OASIS field +++

  ! OASIS needs double precision arrays. These variables are optional.
  ! They are used when calling the subroutine oasisxregcm_setup_field_array(),
  ! when an array is to be initialized. If not given, the array will be
  ! initialized with zeros.
  real(rkx) , parameter :: init_sst = tzero + 25.0
  ! OASIS field +++

  !--------------------------------------------------------------------
  ! accessibility from oasislib:
  public :: oasisxregcm_init , oasisxregcm_finalize ! -> mod_oasis_generic
  public :: oasisxregcm_header                      ! -> mod_oasis_signature

  !--------------------------------------------------------------------
  ! from this module:
  public :: oasisxregcm_params , oasisxregcm_release
  public :: oasisxregcm_def
  public :: oasisxregcm_rcv_all , oasisxregcm_snd_all
  public :: oasisxregcm_sync_wait

  contains

  subroutine oasisxregcm_params
    implicit none
    !--------------------------------------------------------------------------
    oasis_lag = 0
    ! state which grids have to be defined,
    ! depending on what field is activated.
    l_make_grdde   = .false.
    ! OASIS field +++
    l_make_grddi   = .false.
    ! OASIS field +++
    l_make_grdce   = l_cpl_ex_slp
    ! OASIS field +++
    l_make_grdci   = l_cpl_im_sst  .or. &
!                     l_cpl_im_sit  .or. &
                     l_cpl_im_wz0  .or. &
                     l_cpl_im_wust .or. &
                     l_cpl_ex_u10m .or. &
                     l_cpl_ex_v10m .or. &
                     l_cpl_ex_wspd .or. &
                     l_cpl_ex_wdir .or. &
                     l_cpl_ex_t2m  .or. &
!                     l_cpl_ex_t10m .or. &
                     l_cpl_ex_q2m  .or. &
!                     l_cpl_ex_q10m .or. &
                     l_cpl_ex_taux .or. &
                     l_cpl_ex_tauy .or. &
                     l_cpl_ex_z0   .or. &
                     l_cpl_ex_ustr .or. &
                     l_cpl_ex_evap .or. &
                     l_cpl_ex_prec .or. &
                     l_cpl_ex_nuwa .or. &
                     l_cpl_ex_ulhf .or. &
                     l_cpl_ex_ushf .or. &
                     l_cpl_ex_uwlw .or. &
                     l_cpl_ex_dwlw .or. &
                     l_cpl_ex_nulw .or. &
                     l_cpl_ex_uwsw .or. &
                     l_cpl_ex_dwsw .or. &
                     l_cpl_ex_ndsw .or. &
                     l_cpl_ex_rhoa
    ! OASIS field +++
    !
    ! initialize grids: grid variable, name with no mask, name with mask,
    !                   local j1, j2, i1, i2, global j1, j2, i1, i2, number of
    !                   corners per cell
    if ( l_make_grdde )   call oasisxregcm_setup_grid(grdde, 'rden', 'rdem', &
                               jde1, jde2, ide1, ide2, 1, jx, 1, iy, 4)
    if ( l_make_grddi )   call oasisxregcm_setup_grid(grddi, 'rdin', 'rdim', &
                               jdi1, jdi2, idi1, idi2, 2, jx-1, 2, iy-1, 4)
    if ( l_make_grdce )   call oasisxregcm_setup_grid(grdce, 'rcen', 'rcem', &
                               jce1, jce2, ice1, ice2, 1, jx-1, 1, iy-1, 4)
    if ( l_make_grdci )   call oasisxregcm_setup_grid(grdci, 'rcin', 'rcim', &
                               jci1, jci2, ici1, ici2, 2, jx-2, 2, iy-2, 4)
    !
    ! initialize fields: field variable, name, grid, field array (optional), initialization value
    !                                                                        (optional; 0 otherwise)
    if ( l_cpl_im_sst )  call oasisxregcm_setup_field(im_sst,  'RCM_SST',  grdci, cpl_sst, init_sst)
!    if ( l_cpl_im_sit )  call oasisxregcm_setup_field(im_sit,  'RCM_SIT',  grdci, cpl_sit)
    if ( l_cpl_im_wz0 )  call oasisxregcm_setup_field(im_wz0,  'RCM_WZ0',  grdci, cpl_wz0)
    if ( l_cpl_im_wust ) call oasisxregcm_setup_field(im_wust, 'RCM_WUST', grdci, cpl_wust)
    if ( l_cpl_ex_u10m ) call oasisxregcm_setup_field(ex_u10m, 'RCM_U10M', grdci)
    if ( l_cpl_ex_v10m ) call oasisxregcm_setup_field(ex_v10m, 'RCM_V10M', grdci)
    if ( l_cpl_ex_wspd ) call oasisxregcm_setup_field(ex_wspd, 'RCM_WSPD', grdci)
    if ( l_cpl_ex_wdir ) call oasisxregcm_setup_field(ex_wdir, 'RCM_WDIR', grdci, cpl_wdir)
    if ( l_cpl_ex_t2m )  call oasisxregcm_setup_field(ex_t2m,  'RCM_T2M',  grdci)
!    if ( l_cpl_ex_t10m ) call oasisxregcm_setup_field(ex_t10m, 'RCM_T10M', grdci)
    if ( l_cpl_ex_q2m )  call oasisxregcm_setup_field(ex_q2m,  'RCM_Q2M',  grdci)
!    if ( l_cpl_ex_q10m ) call oasisxregcm_setup_field(ex_q10m, 'RCM_Q10M', grdci)
    if ( l_cpl_ex_slp )  call oasisxregcm_setup_field(ex_slp,  'RCM_SLP',  grdce)
    if ( l_cpl_ex_taux ) call oasisxregcm_setup_field(ex_taux, 'RCM_TAUX', grdci)
    if ( l_cpl_ex_tauy ) call oasisxregcm_setup_field(ex_tauy, 'RCM_TAUY', grdci)
    if ( l_cpl_ex_z0 )   call oasisxregcm_setup_field(ex_z0,   'RCM_Z0',   grdci)
    if ( l_cpl_ex_ustr ) call oasisxregcm_setup_field(ex_ustr, 'RCM_USTR', grdci)
    if ( l_cpl_ex_evap ) call oasisxregcm_setup_field(ex_evap, 'RCM_EVAP', grdci)
    if ( l_cpl_ex_prec ) call oasisxregcm_setup_field(ex_prec, 'RCM_PREC', grdci)
    if ( l_cpl_ex_nuwa ) call oasisxregcm_setup_field(ex_nuwa, 'RCM_NUWA', grdci)
    if ( l_cpl_ex_ulhf ) call oasisxregcm_setup_field(ex_ulhf, 'RCM_ULHF', grdci)
    if ( l_cpl_ex_ushf ) call oasisxregcm_setup_field(ex_ushf, 'RCM_USHF', grdci)
    if ( l_cpl_ex_uwlw ) call oasisxregcm_setup_field(ex_uwlw, 'RCM_UWLW', grdci)
    if ( l_cpl_ex_dwlw ) call oasisxregcm_setup_field(ex_dwlw, 'RCM_DWLW', grdci)
    if ( l_cpl_ex_nulw ) call oasisxregcm_setup_field(ex_nulw, 'RCM_NULW', grdci)
    if ( l_cpl_ex_uwsw ) call oasisxregcm_setup_field(ex_uwsw, 'RCM_UWSW', grdci)
    if ( l_cpl_ex_dwsw ) call oasisxregcm_setup_field(ex_dwsw, 'RCM_DWSW', grdci)
    if ( l_cpl_ex_ndsw ) call oasisxregcm_setup_field(ex_ndsw, 'RCM_NDSW', grdci)
    if ( l_cpl_ex_rhoa ) call oasisxregcm_setup_field(ex_rhoa, 'RCM_RHOA', grdci)
    ! OASIS field +++
  end subroutine oasisxregcm_params

  ! call all definition subroutines for setting up OASIS
  subroutine oasisxregcm_def
    implicit none
    !--------------------------------------------------------------------------
    ! partition definition
#ifdef DEBUG
    write(ndebug,*) oasis_prefix, 'definition phase: partitions'
#endif
    if ( l_make_grdde ) call oasisxregcm_def_partition(grdde)
    if ( l_make_grddi ) call oasisxregcm_def_partition(grddi)
    if ( l_make_grdce ) call oasisxregcm_def_partition(grdce)
    if ( l_make_grdci ) call oasisxregcm_def_partition(grdci)
    ! grid definition
#ifdef DEBUG
    write(ndebug,*) oasis_prefix, 'definition phase: grids'
#endif
    if ( l_write_grids ) then
      call oasisxregcm_def_grid
    else
#ifdef DEBUG
      write(ndebug,*) oasis_prefix, '               >> skipped!'
#endif
      if ( myid == italk ) write(stdout,*) 'Note: OASIS grid files', &
                                   ' have not been written online.'
    end if
    ! variable definitions
#ifdef DEBUG
    write(ndebug,*) oasis_prefix, 'definition phase: fields'
#endif
    ! field variable, OASIS integer for the direction:
    !                 importing(OASIS_In)/exporting(OASIS_Out)
    if ( l_cpl_im_sst )  call oasisxregcm_def_field(im_sst,  OASIS_In)
!    if ( l_cpl_im_sit )  call oasisxregcm_def_field(im_sit,  OASIS_In)
    if ( l_cpl_im_wz0 )  call oasisxregcm_def_field(im_wz0,  OASIS_In)
    if ( l_cpl_im_wust ) call oasisxregcm_def_field(im_wust, OASIS_In)
    if ( l_cpl_ex_u10m ) call oasisxregcm_def_field(ex_u10m, OASIS_Out)
    if ( l_cpl_ex_v10m ) call oasisxregcm_def_field(ex_v10m, OASIS_Out)
    if ( l_cpl_ex_wspd ) call oasisxregcm_def_field(ex_wspd, OASIS_Out)
    if ( l_cpl_ex_wdir ) call oasisxregcm_def_field(ex_wdir, OASIS_Out)
    if ( l_cpl_ex_t2m )  call oasisxregcm_def_field(ex_t2m,  OASIS_Out)
!    if ( l_cpl_ex_t10m ) call oasisxregcm_def_field(ex_t10m, OASIS_Out)
    if ( l_cpl_ex_q2m )  call oasisxregcm_def_field(ex_q2m,  OASIS_Out)
!    if ( l_cpl_ex_q10m ) call oasisxregcm_def_field(ex_q10m, OASIS_Out)
    if ( l_cpl_ex_slp )  call oasisxregcm_def_field(ex_slp,  OASIS_Out)
    if ( l_cpl_ex_taux ) call oasisxregcm_def_field(ex_taux, OASIS_Out)
    if ( l_cpl_ex_tauy ) call oasisxregcm_def_field(ex_tauy, OASIS_Out)
    if ( l_cpl_ex_z0 )   call oasisxregcm_def_field(ex_z0,   OASIS_Out)
    if ( l_cpl_ex_ustr ) call oasisxregcm_def_field(ex_ustr, OASIS_Out)
    if ( l_cpl_ex_evap ) call oasisxregcm_def_field(ex_evap, OASIS_Out)
    if ( l_cpl_ex_prec ) call oasisxregcm_def_field(ex_prec, OASIS_Out)
    if ( l_cpl_ex_nuwa ) call oasisxregcm_def_field(ex_nuwa, OASIS_Out)
    if ( l_cpl_ex_ulhf ) call oasisxregcm_def_field(ex_ulhf, OASIS_Out)
    if ( l_cpl_ex_ushf ) call oasisxregcm_def_field(ex_ushf, OASIS_Out)
    if ( l_cpl_ex_uwlw ) call oasisxregcm_def_field(ex_uwlw, OASIS_Out)
    if ( l_cpl_ex_dwlw ) call oasisxregcm_def_field(ex_dwlw, OASIS_Out)
    if ( l_cpl_ex_nulw ) call oasisxregcm_def_field(ex_nulw, OASIS_Out)
    if ( l_cpl_ex_uwsw ) call oasisxregcm_def_field(ex_uwsw, OASIS_Out)
    if ( l_cpl_ex_dwsw ) call oasisxregcm_def_field(ex_dwsw, OASIS_Out)
    if ( l_cpl_ex_ndsw ) call oasisxregcm_def_field(ex_ndsw, OASIS_Out)
    if ( l_cpl_ex_rhoa ) call oasisxregcm_def_field(ex_rhoa, OASIS_Out)
    ! OASIS field +++
    ! termination of definition phase
#ifdef DEBUG
    write(ndebug,*) oasis_prefix, 'definition phase: end'
#endif
    call oasisxregcm_end_def

  contains

  ! define OASIS grids
  subroutine oasisxregcm_def_grid
    implicit none
    real(rkx) , pointer , dimension(:,:) :: dlon , dlat ! dot degree coordinates
    real(rkx) , pointer , dimension(:,:) :: xlon , xlat ! cross degree coordinates
    real(rkx) , pointer , dimension(:,:) :: lndcat ! land category (15 for ocean)
    real(rkx) , pointer , dimension(:,:) :: oasisgrid_lon , oasisgrid_lat ! lon, lat
    real(rkx) , pointer , dimension(:,:,:) :: oasisgrid_clon , & ! lon &
                                              oasisgrid_clat     ! lat of the corners
    ! note: 4 corners (always the case?) in the third dimension, counterclockwisely.
    real(rkx) , pointer , dimension(:,:) :: oasisgrid_srf ! surface of the grid meshes m2
    integer(ik4) , pointer , dimension(:,:) :: oasisgrid_mask ! mask, 0 = valid, 1 = mask
    !                                                          (OASIS convention)
    integer(ik4) :: il_flag ! flag for grid writing by proc 0
    integer(ik4) :: ierror
#ifdef DEBUG
    !--------------------------------------------------------------------------
    write(ndebug,*) oasis_prefix, 'start collecting grid data'
    write(ndebug,*) oasis_prefix, 'grids writing is done by the in/out cpu only'
#endif
    ! collect the domain definition information in RegCM
    allocate(dlon(jx,iy),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating dlon'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(dlat(jx,iy),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating dlat'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(xlon(jx,iy),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating xlon'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(xlat(jx,iy),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating xlat'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(lndcat(jx,iy),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating lndcat'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    call grid_collect(mddom%dlon,dlon,jde1,jde2,ide1,ide2)
    call grid_collect(mddom%dlat,dlat,jde1,jde2,ide1,ide2)
    call grid_collect(mddom%xlon,xlon,jde1,jde2,ide1,ide2)
    call grid_collect(mddom%xlat,xlat,jde1,jde2,ide1,ide2)
    call grid_collect(mddom%lndcat,lndcat,jde1,jde2,ide1,ide2)
    if ( myid == iocpu ) then
      ! start the grid writing process
#ifdef DEBUG
      write(ndebug,*) oasis_prefix, 'start grids writing'
#endif
      call oasis_start_grids_writing(il_flag)
      !
      if ( l_make_grdde .or. l_make_grddi) then ! dots
        !
        call oasisxregcm_make_oasisgrids_d( &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask, &
             dlon,dlat,xlon,xlat,lndcat) ! subroutine writen below
        !
        if ( l_make_grdde ) call oasisxregcm_write_oasisgrids(grdde, &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask)
        !
        if ( l_make_grddi ) call oasisxregcm_write_oasisgrids(grddi, &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask)
        !
        call oasisxregcm_deallocate_oasisgrids( &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask)
        !
      end if
      !
      if ( l_make_grdce .or. l_make_grdci) then ! crosses
        !
        call oasisxregcm_make_oasisgrids_c( &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask, &
             dlon,dlat,xlon,xlat,lndcat) ! subroutine writen below
        !
        if ( l_make_grdce ) call oasisxregcm_write_oasisgrids(grdce, &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask)
        !
        if ( l_make_grdci ) call oasisxregcm_write_oasisgrids(grdci, &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask)
        !
        call oasisxregcm_deallocate_oasisgrids( &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask)
        !
      end if
      ! terminate the grid writing process
#ifdef DEBUG
      write(ndebug,*) oasis_prefix, 'terminate grids writing'
#endif
      call oasis_terminate_grids_writing()
    end if
    if ( associated(dlon) ) nullify(dlon)
    if ( associated(dlat) ) nullify(dlat)
    if ( associated(xlon) ) nullify(xlon)
    if ( associated(xlat) ) nullify(xlat)
    if ( associated(lndcat) ) nullify(lndcat)
#ifdef DEBUG
    write(ndebug,"(' ',A,A,I3,A)") oasis_prefix, 'please refer to cpu ', iocpu, ' for grid writing'&
                           //' debug statements'
#endif
  end subroutine oasisxregcm_def_grid

  subroutine oasisxregcm_make_oasisgrids_d(lon,lat,clon,clat,srf,mask, &
                                           dlon,dlat,xlon,xlat,lndcat)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(in) :: dlon , dlat , &
                                                         xlon , xlat , lndcat
    real(rkx) , pointer , dimension(:,:) , intent(out) :: lon , lat , srf
    real(rkx) , pointer , dimension(:,:,:) , intent(out) :: clon , clat
    integer(ik4) , pointer , dimension(:,:) , intent(out) :: mask
    integer(ik4) :: i , j
    !--------------------------------------------------------------------------
    !
    !    ^    corner 2        corner 1
    !    |    x(j-1,i)         x(j,i)
    !  i |
    !    |             center
    !  a |             d(j,i)
    !  x |
    !  i |    corner 3        corner 4
    !  s |   x(j-1,i-1)       x(j,i-1)
    !    |
    !    +-------------------------------->
    !                 j axis
    !
    ! note: concretely xlat and xlon are defined on (1:jx,1:iy) as for dlat and dlon
    !       though the cross grid is read on (1:jx-1,1:iy-1) only
    ! consequence: xlon and xlat exist at (jx,jy) -> no exception needed
    !              however, exceptions are needed for extrem west and south values
    !                                                         0         0
    call oasisxregcm_allocate_oasisgrids(lon,lat,clon,clat,srf,mask, &
                                         jx,iy,4)
    do i = 1 , iy
      do j = 1 , jx
        ! center
        lon(j,i) = dlon(j,i)
        lat(j,i) = dlat(j,i)
        ! upper right corner: no problem
        clon(j,i,1) = xlon(j,i)
        clat(j,i,1) = xlat(j,i)
        ! left corners:
        !   upper left corner: use of j-1 -> exception for j=1
        !                      lon -> use of the last lon step to the West
        !                      lat -> use of the same value as for j
        !                             (expected to be the same)
        !   lower left corner: same + use of i-1 -> exception for i=1
        !                      lon -> use of the same value as for i
        !                             (expected to be the same)
        !                      lat -> use of the last lat step to the South
        if ( j == 1 ) then
          clon(j,i,2) = 2 * xlon(j,i) - xlon(j+1,i)
          clat(j,i,2) = xlat(j,i)
          if ( i == 1 ) then
            clon(j,i,3) = 2 * xlon(j,i) - xlon(j+1,i)
            clat(j,i,3) = 2 * xlat(j,i) - xlat(j,i+1)
          else
            clon(j,i,3) = 2 * xlon(j,i-1) - xlon(j+1,i-1)
            clat(j,i,3) = xlat(j,i-1)
          end if
        else
          clon(j,i,2) = xlon(j-1,i)
          clat(j,i,2) = xlat(j-1,i)
          if ( i == 1 ) then
            clon(j,i,3) = xlon(j-1,i)
            clat(j,i,3) = 2 * xlat(j-1,i) - xlat(j-1,i+1)
          else
            clon(j,i,3) = xlon(j-1,i-1)
            clat(j,i,3) = xlat(j-1,i-1)
          end if
        end if
        ! lower right corner: use of i-1 -> exception for i=1
        if ( i == 1 ) then
          clon(j,i,4) = xlon(j,i)
          clat(j,i,4) = 2 * xlat(j,i) - xlat(j,i+1)
        else
          clon(j,i,4) = xlon(j,i-1)
          clat(j,i,4) = xlat(j,i-1)
        end if
        ! surface
        srf(j,i) = srf_sqm(clon(j,i,:),clat(j,i,:))
        ! mask
        if ( isocean(lndcat(j,i)) ) then
          mask(j,i) = 0
        else
          mask(j,i) = 1
        end if
      end do
    end do
  end subroutine oasisxregcm_make_oasisgrids_d

  subroutine oasisxregcm_make_oasisgrids_c(lon,lat,clon,clat,srf,mask, &
                                           dlon,dlat,xlon,xlat,lndcat)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(in) :: dlon , dlat , &
                                                         xlon , xlat , lndcat
    real(rkx) , pointer , dimension(:,:) , intent(out) :: lon , lat , srf
    real(rkx) , pointer , dimension(:,:,:) , intent(out) :: clon , clat
    integer(ik4) , pointer , dimension(:,:) , intent(out) :: mask
    integer(ik4) :: i , j
    !--------------------------------------------------------------------------
    !
    !    ^    corner 2        corner 1
    !    |    d(j,i+1)       d(j+1,i+1)
    !  i |
    !    |             center
    !  a |             x(j,i)
    !  x |
    !  i |    corner 3        corner 4
    !  s |     d(j,i)         d(j+1,i)
    !    |
    !    +-------------------------------->
    !                 j axis
    !
    ! note: the cross case is easier as the cross grid is contained
    !       in the dot one thus all corners are already defined in dlon
    !
    call oasisxregcm_allocate_oasisgrids(lon,lat,clon,clat,srf,mask, &
                                         jx-1,iy-1,4)
    do i = 1 , iy-1
      do j = 1 , jx-1
        ! center
        lon(j,i) = xlon(j,i)
        lat(j,i) = xlat(j,i)
        ! upper right corner
        clon(j,i,1) = dlon(j+1,i+1)
        clat(j,i,1) = dlat(j+1,i+1)
        ! upper left corner
        clon(j,i,2) = dlon(j,i+1)
        clat(j,i,2) = dlat(j,i+1)
        ! lower left corner
        clon(j,i,3) = dlon(j,i)
        clat(j,i,3) = dlat(j,i)
        ! lower right corner
        clon(j,i,4) = dlon(j+1,i)
        clat(j,i,4) = dlat(j+1,i)
        ! surface
        srf(j,i) = srf_sqm(clon(j,i,:),clat(j,i,:))
        ! mask
        if ( isocean(lndcat(j,i)) ) then
          mask(j,i) = 0
        else
          mask(j,i) = 1
        end if
      end do
    end do
  end subroutine oasisxregcm_make_oasisgrids_c

  end subroutine oasisxregcm_def

  ! call all subroutines consisting of receiving OASIS fields
  ! and optionally reworking them
  subroutine oasisxregcm_rcv_all(time)
    implicit none
    integer(ik4) , intent(in) :: time ! execution time
    logical :: l_act
    !--------------------------------------------------------------------------
#ifdef DEBUG
    if ( time == 0 ) then
      write(ndebug,*) oasis_prefix, 'Initialization:'
    else if ( debug_level > 0 ) then
      write(ndebug,*) oasis_prefix, 'OASIS time: ', time
    end if
#endif
    !
    if ( l_cpl_im_sst ) then ! sea surface temperature [K]
      call oasisxregcm_rcv(cpl_sst,im_sst,time,l_act)
      !if ( l_act ) then
        ! fill the ocean parts: out array (2d or 3d with nsg dimension)
        !                       in array, mask, field grid
        call fill_ocean(sfs%tg,cpl_sst,mddom%lndcat,im_sst%grd)
        call fill_ocean(sfs%tgbb,cpl_sst,mddom%lndcat,im_sst%grd)
      !end if
    end if
    ! NOT IMPLEMENTED YET, BELOW IS A COPY OF THE PROCESS IN
    ! MOD_LM_INTERFACE.F90; IMPORT_DATA_INTO_SURFACE.
!    if ( l_cpl_im_sit ) then ! sea ice thickness [m]
!      call oasisxregcm_rcv(cpl_sit,im_sit,time,l_act)
!!      if ( l_act ) then
!        grd => im_sit%grd
!        ! dans le mask
!        ! tol = 0.5 E 20
!          if ( impfie%sit(j,i) < tol .and. lm%ldmsk(j,i) /= 1 ) then
!            if ( impfie%sit(j,i) > iceminh ) then
!              flag = .false.
!              if ( lm%ldmsk(j,i) == 0 ) flag = .true.
!              ! set land-sea mask
!              lm%ldmsk(j,i) = 2
!              do n = 1, nnsg
!                lm%ldmsk1(n,j,i) = 2
!                ! set sea ice thikness (in meter)
!                lms%sfice(n,j,i) = impfie%sit(j,i)
!              end do
!              ! write debug info
!              if ( flag ) then
!                write(*,f99002) j, i, 'water', 'ice  ', &
!                   lm%ldmsk(j,i), lms%sfice(1,j,i)
!              end if
!            else
!              if ( ldmskb(j,i) == 0 .and. lm%ldmsk(j,i) == 2 ) then
!                do n = 1, nnsg
!                  ! reduce to one tenth surface ice: it should melt away
!                  lms%sfice(n,j,i) = lms%sfice(n,j,i)*d_r10
!                  ! check that sea ice is melted or not
!                  if ( lms%sfice(n,j,i) <= iceminh ) then
!                    if ( ldmskb(j,i) /= lm%ldmsk(j,i) ) flag = .true.
!                    ! set land-sea mask to its original value
!                    lm%ldmsk(j,i) = ldmskb(j,i)
!                    lm%ldmsk1(n,j,i) = ldmskb(j,i)
!                    ! set land-use type to its original value
!                    ! set sea ice thikness (in meter)
!                    lms%sfice(n,j,i) = d_zero
!                  else
!                    flag = .false.
!                  end if
!                end do
!                ! write debug info
!                if ( flag ) then
!                  write(*,f99003) j, i, 'ice  ', 'water',  &
!                    lm%ldmsk(j,i), lms%sfice(1,j,i)
!                end if
!              end if
!            end if
!          end if
!        ???
!        nullify(grd)
!!      end if
!    end if
    !
    if ( l_cpl_im_wz0 ) then ! surface roughness length [m]
      call oasisxregcm_rcv(cpl_wz0,im_wz0,time,l_act)
!      if ( l_act ) then
        call fill_ocean(sfs%zo,cpl_wz0,mddom%lndcat,im_wz0%grd)
!      end if
    end if
    !
    if ( l_cpl_im_wust ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_wust,im_wust,time,l_act)
!      if ( l_act ) then
        call fill_ocean(sfs%ustar ,cpl_wust,mddom%lndcat,im_wust%grd)
!      end if
    end if
    ! OASIS field +++
    !
  end subroutine oasisxregcm_rcv_all

  ! call all subroutines consisting of sending OASIS fields
  ! with optional prior reworking
  subroutine oasisxregcm_snd_all(time,l_last_time)
    implicit none
    integer(ik4) , intent(in) :: time ! execution time
    logical , intent(in) :: l_last_time
    type(infogrd) , pointer :: grd
    integer(ik4) :: i , j , ishift , jshift
    logical :: l_write_restart
    !--------------------------------------------------------------------------
    l_write_restart = (write_restart_option == 1 .and. time == 0) .or. &
                      (write_restart_option == 2 .and. l_last_time) .or. &
                      (write_restart_option == 3 .and. (time == 0 .or. l_last_time))
    !
    if ( l_cpl_ex_u10m ) then ! eatward near-surface wind (10m-height) [m.s-1]
      grd => ex_u10m%grd
      call oasisxregcm_snd( &
           sfs%u10m(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_u10m, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_v10m ) then ! northward near-surface wind (10m-height) [m.s-1]
      grd => ex_v10m%grd
      call oasisxregcm_snd( &
           sfs%v10m(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_v10m, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_wspd ) then ! near-surface wind speed (10m-height) [m.s-1]
      grd => ex_wspd%grd
      call oasisxregcm_snd( &
           sqrt(sfs%u10m(grd%j1:grd%j2 , grd%i1:grd%i2)**2   &
              + sfs%v10m(grd%j1:grd%j2 , grd%i1:grd%i2)**2), &
           ex_wspd, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_wdir ) then ! near-surface wind direction (10m-height) [degree]
      grd => ex_wdir%grd
      do i = grd%i1 , grd%i2
        ishift = i - grd%i1 + 1
        do j = grd%j1 , grd%j2
          jshift = j - grd%j1 + 1
          cpl_wdir(jshift,ishift) = atan2( sfs%u10m(j,i) , sfs%v10m(j,i) ) * raddeg
          if ( cpl_wdir(jshift,ishift) < d_zero ) &
               cpl_wdir(jshift,ishift) = cpl_wdir(jshift,ishift) + 360_rkx
        end do
      end do
      call oasisxregcm_snd(cpl_wdir, ex_wdir, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_t2m ) then ! near-surface air temperature (2m-height) [K]
      grd => ex_t2m%grd
      call oasisxregcm_snd( &
           sum( lms%t2m(: , grd%j1:grd%j2 , grd%i1:grd%i2) , 1 ) * rdnnsg, &
           ex_t2m, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
!    if ( l_cpl_ex_t10m ) then ! near-surface air temperature (10m-height) [K]
!      grd => ex_t10m%grd
!      call oasisxregcm_snd( &
!           ! ???
!           ex_t10m, time, .false. .or. l_write_restart)
!      nullify(grd)
!    end if
    !
    if ( l_cpl_ex_q2m ) then ! near-surface air specific humidity (2m-height) [1]
      grd => ex_q2m%grd
      call oasisxregcm_snd( &
           sfs%q2m(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_q2m, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
!    if ( l_cpl_ex_q10m ) then ! near-surface air specific humidity (10m-height) [1]
!      grd => ex_q10m%grd
!      call oasisxregcm_snd( &
!           ! ???
!           ex_q10m, time, .false. .or. l_write_restart)
!      nullify(grd)
!    end if
    !
    if ( l_cpl_ex_slp ) then ! sea level pressure [Pa]
      grd => ex_slp%grd
      call oasisxregcm_snd( &
           atms%ps2d(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_slp, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_taux) then ! surface eastward wind stress [Pa]
      grd => ex_taux%grd
      call oasisxregcm_snd( &
           sum( lms%taux(: , grd%j1:grd%j2 , grd%i1:grd%i2) , 1 ) * rdnnsg, &
           ex_taux, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_tauy) then ! surface northward wind stress [Pa]
      grd => ex_tauy%grd
      call oasisxregcm_snd( &
           sum( lms%tauy(: , grd%j1:grd%j2 , grd%i1:grd%i2) , 1 ) * rdnnsg, &
           ex_tauy, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_z0 ) then ! surface roughness length [m]
      grd => ex_z0%grd
      call oasisxregcm_snd( &
           sfs%zo(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_z0, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_ustr ) then ! surface friction velocity [s-1]
      grd => ex_ustr%grd
      call oasisxregcm_snd( &
           sfs%ustar(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_ustr, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_evap ) then ! water evaporation flux [kg.m-2.s-1]
      grd => ex_evap%grd
      call oasisxregcm_snd( &
           sfs%qfx(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_evap, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_prec ) then ! precipitation flux [kg.m-2.s-1]
      grd => ex_prec%grd
      call oasisxregcm_snd( &
           sum( lms%prcp(: , grd%j1:grd%j2 , grd%i1:grd%i2) , 1 ) * rdnnsg, &
           ex_prec, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_nuwa ) then ! net upward water flux [kg.m-2.s-1]
      grd => ex_nuwa%grd
      call oasisxregcm_snd( &
           sum( lms%evpr(: , grd%j1:grd%j2 , grd%i1:grd%i2) - &
                lms%prcp(: , grd%j1:grd%j2 , grd%i1:grd%i2) , 1 ) * rdnnsg, &
           ex_nuwa, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_ulhf ) then ! surface upward latent heat flux [W.m-2]
      grd => ex_ulhf%grd
      call oasisxregcm_snd( &
           sum( lms%evpr(: , grd%j1:grd%j2 , grd%i1:grd%i2) &
           *wlh(lms%prcp(: , grd%j1:grd%j2 , grd%i1:grd%i2)) , 1 ) * rdnnsg, &
           ex_ulhf, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_ushf ) then ! surface upward sensible heat flux [W.m-2]
      grd => ex_ushf%grd
      call oasisxregcm_snd( &
           sfs%hfx(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_ushf, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_uwlw ) then ! surface upwelling long-wave radiation flux [W.m-2]
      grd => ex_uwlw%grd
      call oasisxregcm_snd( &
           flwd(grd%j1:grd%j2 , grd%i1:grd%i2) &
           + flw(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_uwlw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_dwlw ) then ! surface downwelling long-wave radiation flux [W.m-2]
      grd => ex_dwlw%grd
      call oasisxregcm_snd( &
           flwd(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_dwlw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_nulw ) then ! surface net upward long-wave radiation flux [W.m-2]
      grd => ex_nulw%grd
      call oasisxregcm_snd( &
           flw(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_nulw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_uwsw ) then ! surface upwelling short-wave radiation flux [W.m-2]
      grd => ex_uwsw%grd
      call oasisxregcm_snd( &
           sinc(grd%j1:grd%j2 , grd%i1:grd%i2) &
           - fsw(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_uwsw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_dwsw ) then ! surface downwelling short-wave radiation flux [W.m-2]
      grd => ex_dwsw%grd
      call oasisxregcm_snd( &
           sinc(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_dwsw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_ndsw ) then ! surface net downward short-wave radiation flux [W.m-2]
      grd => ex_ndsw%grd
      call oasisxregcm_snd( &
           fsw(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_ndsw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_rhoa ) then ! surface air density [kg.m-3]
      grd => ex_rhoa%grd
      call oasisxregcm_snd( &
           sum( lms%rhoa(:, grd%j1:grd%j2 , grd%i1:grd%i2) , 1 ) * rdnnsg, &
           ex_rhoa, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    ! OASIS field +++
    !
    if ( myid == italk ) then
      if ((write_restart_option == 1 .or. write_restart_option == 3) &
          .and. time == 0) then
        write(stdout,*) 'Note: OASIS restart files written at the' &
                      //' first time.'
 else if ((write_restart_option == 2 .or. write_restart_option == 3) &
          .and. l_last_time) then
        write(stdout,*) 'Note: OASIS restart files written at the' &
                      //' last time.'
      end if
    end if

    contains

#include <wlh.inc>

  end subroutine oasisxregcm_snd_all

  ! wait during the time indicated through oasis_sync_lag
  ! to make RegCM model time not synchronised with
  ! other coupled components
  ! update oasis_lag such that
  !   OASIS time = RegCM time + oasis_lag
  subroutine oasisxregcm_sync_wait(time)
    implicit none
    integer(ik4) , intent(in) :: time ! execution time
    !--------------------------------------------------------------------------
    if ( oasis_lag /= 0 ) then
      write(stderr,*) 'error initiating the lag with OASIS'
      call fatal(__FILE__,__LINE__,'SYNC WAIT')
    end if
    ! case number 1: oasis_sync_lag > 0
    ! RegCM actually starts oasis_sync_lag seconds after OASIS
    ! >> at the beginning of the run, OASIS starts while RegCM runs
    ! dummy loops with oasis_lag increasing like 0 -> oasis_sync_lag
    ! at the end of the run, OASIS and RegCM stop synchronously
    if ( oasis_sync_lag > 0 ) then
      do while ( oasis_lag < oasis_sync_lag )
        call oasisxregcm_rcv_all(time + oasis_lag)
        ! dummy loop
#ifdef DEBUG
        if ( oasis_lag == 0 ) then
          write(ndebug,*) oasis_prefix, 'RegCM is creating the lag ', &
                                        'with the other components'
          write(ndebug,*) oasis_prefix, ' oasis_lag: ', oasis_lag
        end if
#endif
        call oasisxregcm_snd_all(time + oasis_lag, .false.)
        oasis_lag = oasis_lag + int(dtsec,ik4)
#ifdef DEBUG
        write(ndebug,*) oasis_prefix, ' oasis_lag: ', oasis_lag
#endif
      end do
    ! case number 2: oasis_sync_lag < 0
    ! RegCM starts synchronously with OASIS
    ! however, some other components start with a lag
    ! at the end of the run, OASIS needs oasis_sync_lag seconds to run
    ! with the other components while RegCM has already finished
    ! >> at the end, OASIS keeps going while RegCM runs dummy loops
    ! with oasis_lag increasing like end_time + dt -> -oasis_sync_lag
    else if ( oasis_sync_lag < 0 ) then
#ifdef DEBUG
      write(ndebug,*) oasis_prefix, 'RegCM is catching up ', &
                                    'with the other components'
      write(ndebug,*) oasis_prefix, ' oasis_lag: ', oasis_lag
#endif
      do while ( oasis_lag < -oasis_sync_lag )
        call oasisxregcm_rcv_all(time + int(dtsec,ik4) + oasis_lag)
        ! dummy loop
        call oasisxregcm_snd_all(time + int(dtsec,ik4) + oasis_lag, .false.)
        oasis_lag = oasis_lag + int(dtsec,ik4)
#ifdef DEBUG
        write(ndebug,*) oasis_prefix, ' oasis_lag: ', oasis_lag
#endif
      end do
    end if
  end subroutine oasisxregcm_sync_wait

  ! call all subroutines linked with some deallocation of variables used in
  ! the OASIS coupling
  subroutine oasisxregcm_release
    implicit none
    !--------------------------------------------------------------------------
    call oasisxregcm_deallocate_field(im_sst, cpl_sst)
!    call oasisxregcm_deallocate_field(im_sit, cpl_sit)
    call oasisxregcm_deallocate_field(im_wz0, cpl_wz0)
    call oasisxregcm_deallocate_field(im_wust, cpl_wust)
    call oasisxregcm_deallocate_field(ex_u10m)
    call oasisxregcm_deallocate_field(ex_v10m)
    call oasisxregcm_deallocate_field(ex_wspd)
    call oasisxregcm_deallocate_field(ex_wdir, cpl_wdir)
    call oasisxregcm_deallocate_field(ex_t2m)
!    call oasisxregcm_deallocate_field(ex_t10m)
    call oasisxregcm_deallocate_field(ex_q2m)
!    call oasisxregcm_deallocate_field(ex_q10m)
    call oasisxregcm_deallocate_field(ex_slp)
    call oasisxregcm_deallocate_field(ex_taux)
    call oasisxregcm_deallocate_field(ex_tauy)
    call oasisxregcm_deallocate_field(ex_z0)
    call oasisxregcm_deallocate_field(ex_ustr)
    call oasisxregcm_deallocate_field(ex_evap)
    call oasisxregcm_deallocate_field(ex_prec)
    call oasisxregcm_deallocate_field(ex_nuwa)
    call oasisxregcm_deallocate_field(ex_ulhf)
    call oasisxregcm_deallocate_field(ex_ushf)
    call oasisxregcm_deallocate_field(ex_uwlw)
    call oasisxregcm_deallocate_field(ex_dwlw)
    call oasisxregcm_deallocate_field(ex_nulw)
    call oasisxregcm_deallocate_field(ex_uwsw)
    call oasisxregcm_deallocate_field(ex_dwsw)
    call oasisxregcm_deallocate_field(ex_ndsw)
    call oasisxregcm_deallocate_field(ex_rhoa)
    ! OASIS field +++
    if ( allocated(grdde) ) deallocate(grdde)
    if ( allocated(grddi) ) deallocate(grddi)
    if ( allocated(grdce) ) deallocate(grdce)
    if ( allocated(grdci) ) deallocate(grdci)
  end subroutine oasisxregcm_release

end module mod_oasis_interface
!
#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
