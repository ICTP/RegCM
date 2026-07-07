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

module mod_ncio

  use, intrinsic :: iso_c_binding
  use mod_intkinds
  use mod_realkinds
  use mod_date
  use mod_stdio
  use mod_constants
  use netcdf
  use mod_runparams
  use mod_dynparam
  use mod_ensemble
  use mod_mpmessage
  use mod_mppparam
  use mod_memutil
  use mod_nchelper
  use mod_ncstream, only : outstream_async_flush, outstream_netcdf_lock, &
                           outstream_netcdf_unlock
#ifdef ASYNC_NETCDF
  use mod_async_netcdf, only : async_netcdf_buffer, async_netcdf_counter, &
                               async_netcdf_configured_cap_bytes, &
                               async_netcdf_counter_fail, &
                               async_netcdf_counter_poll, &
                               async_netcdf_counter_reset, &
                               async_netcdf_counter_wait, &
                               async_netcdf_get_var_rkx, &
                               async_netcdf_initialize, &
                               async_netcdf_release_buffer, &
                               async_netcdf_try_acquire_buffer
#endif
  use mod_domain
  implicit none

  private

  public :: read_domain_info, read_subdomain_info
  public :: open_icbc, icbc_search, read_icbc, close_icbc
#ifdef ASYNC_NETCDF
  public :: warmup_icbc_prefetch, prefetch_icbc, consume_icbc_prefetch
#endif
  public :: open_som, som_search, read_som, close_som
  public :: open_clmbc, clmbc_search, close_clmbc, read_clmbc
  public :: fixqcqi, we_have_qc, we_have_qi
  public :: read_ccn, tccn_data

  integer(ik4) :: ibcin, somin, clmbcin
  integer(ik4) :: istatus
  integer(ik4) :: ibcrec, ibcnrec, clmbcrec, clmbcnrec
  integer(ik4) :: somrec
  type(rcm_time_and_date), dimension(:), allocatable :: icbc_idate
  type(rcm_time_and_date), dimension(:), allocatable :: clmbc_idate
  integer(ik4), dimension(16) :: icbc_ivar
  integer(ik4), dimension(4) :: clmbc_ivar
  integer(ik4), dimension(1) :: som_ivar
  logical :: has_qc = .false.
  logical :: has_qi = .false.

  data ibcin   /-1/
  data clmbcin /-1/
  data ibcrec  / 1/
  data ibcnrec / 0/
  data clmbcrec  / 1/
  data clmbcnrec / 0/
  data somin   /-1/
  data somrec  / 1/

  real(rkx), dimension(:,:), pointer, contiguous :: rspacesom => null()
  real(rkx), dimension(:,:), pointer, contiguous :: rspace2 => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: rspace3 => null()
#ifdef OPENACC
  real(rkx), dimension(:,:), pointer, contiguous :: rspace2_raw => null()
  real(rkx), dimension(:,:,:), pointer, contiguous :: rspace3_raw => null()
  type(c_ptr) :: cptr
  integer(c_int) :: istat
  integer(c_size_t) :: size_bytes
#endif
  real(rkx), dimension(:,:,:), pointer, contiguous :: tempw => null()
  real(rkx), dimension(:,:), pointer, contiguous :: tempwtop => null()

#ifdef ASYNC_NETCDF
  integer(c_int64_t), parameter :: icbc_prefetch_alignment = 128_c_int64_t

  type icbc_prefetch_slot
    logical :: active = .false.
    type(rcm_time_and_date) :: target_date
    integer(ik4) :: ncid = -1
    integer(ik4) :: record = -1
    integer(ik4) :: nitems = 0
    type(async_netcdf_buffer) :: buffer
    type(async_netcdf_counter) :: counter
    real(rkx), pointer, contiguous, dimension(:) :: ps_flat => null()
    real(rkx), pointer, contiguous, dimension(:) :: ts_flat => null()
    real(rkx), pointer, contiguous, dimension(:) :: u_flat => null()
    real(rkx), pointer, contiguous, dimension(:) :: v_flat => null()
    real(rkx), pointer, contiguous, dimension(:) :: t_flat => null()
    real(rkx), pointer, contiguous, dimension(:) :: qv_flat => null()
    real(rkx), pointer, contiguous, dimension(:) :: qc_flat => null()
    real(rkx), pointer, contiguous, dimension(:) :: qi_flat => null()
    real(rkx), pointer, contiguous, dimension(:) :: pp_flat => null()
    real(rkx), pointer, contiguous, dimension(:) :: ww_flat => null()
    real(rkx), pointer, contiguous, dimension(:) :: wtop_flat => null()
    real(rkx), pointer, contiguous, dimension(:,:) :: ps => null()
    real(rkx), pointer, contiguous, dimension(:,:) :: ts => null()
    real(rkx), pointer, contiguous, dimension(:,:,:) :: u => null()
    real(rkx), pointer, contiguous, dimension(:,:,:) :: v => null()
    real(rkx), pointer, contiguous, dimension(:,:,:) :: t => null()
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qv => null()
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qc => null()
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qi => null()
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pp => null()
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ww => null()
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pai => null()
    real(rkx), pointer, contiguous, dimension(:,:) :: wtop => null()
  end type icbc_prefetch_slot

  type(icbc_prefetch_slot), save, target :: icbc_prefetch
#endif

  type :: tccn_data
    integer(ik4) :: nlon, nlat, nlev
    real(rkx), pointer, dimension(:), contiguous :: lon
    real(rkx), pointer, dimension(:), contiguous :: lat
    real(rkx), pointer, dimension(:), contiguous :: altitude
    real(rkx), pointer, dimension(:,:,:,:), contiguous :: ccn
  end type

  contains

#include <pfesat.inc>

  subroutine fixqcqi( )
    implicit none
    call bcast(has_qc)
    call bcast(has_qi)
  end subroutine fixqcqi

  logical function we_have_qc( )
    implicit none
    we_have_qc = has_qc
  end function we_have_qc
  logical function we_have_qi( )
    implicit none
    we_have_qi = has_qi
  end function we_have_qi

  subroutine read_domain_info(ht,lnd,tex,mask,area,xlat,xlon,dlat,dlon, &
                              ulat,ulon,vlat,vlon,msfx,msfd,msfu,msfv,  &
                              coriol,snowam,smoist,rmoist,rts,hlake,ts0,&
                              htu,htv,rlat,rlon)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: ht
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: lnd
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: tex
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: mask
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: area
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: xlat
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: xlon
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: dlat
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: dlon
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: ulat
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: ulon
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: vlat
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: vlon
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: msfx
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: msfd
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: msfu
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: msfv
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: coriol
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: snowam
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: smoist
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: rmoist
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: rts
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: hlake
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: htu
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: htv
    real(rkx), pointer, contiguous, dimension(:), intent(inout) :: rlat
    real(rkx), pointer, contiguous, dimension(:), intent(inout) :: rlon
    real(rkx), intent(out) :: ts0
    real(rkx), dimension(:,:), pointer, contiguous :: tempmoist
    character(len=256) :: dname
    character(len=8) :: csmoist, cstemp
    integer(ik4) :: idmin, ilev
    integer(ik4), dimension(2) :: istart, icount
    integer(ik4), dimension(3) :: istart3, icount3
    real(rkx), dimension(:,:), pointer, contiguous :: rspace
    logical :: has_snow
    logical :: has_dhlake
    logical :: lerror

    call outstream_async_flush()

    has_snow = .true.
    has_dhlake = .true.
    dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
    if ( myid == italk ) then
      write(stdout,*) 'Reading Domain file : ', trim(dname)
    end if

    if ( do_parallel_netcdf_in ) then
      istart(1) = jde1
      istart(2) = ide1
      icount(1) = jde2-jde1+1
      icount(2) = ide2-ide1+1
      allocate(rspace(icount(1),icount(2)))
      call openfile_withname(dname,idmin)
      call check_domain(idmin)
      call get_attribute(idmin,'initialized_soil_moisture',csmoist)
      if ( csmoist == 'Yes' ) then
        replacemoist = .true.
      else
        replacemoist = .false.
      end if
      call get_attribute(idmin,'initialized_soil_temperature',cstemp)
      if ( cstemp == 'Yes' ) then
        replacetemp = .true.
      else
        replacetemp = .false.
      end if
      lerror = .true.
      call read_var1d_static(idmin,'kz',sigma,lerror)
      if ( .not. lerror ) then
        call read_var1d_static(idmin,'sigma',sigma)
      end if
      if ( idynamic == 3 ) then
        call read_var1d_static(idmin,'rlat',rlat,lerror)
        call read_var1d_static(idmin,'rlon',rlon,lerror)
      end if
      call read_var2d_static(idmin,'xlat',rspace,istart=istart,icount=icount)
      xlat(jde1:jde2,ide1:ide2) = rspace
      call read_var2d_static(idmin,'xlon',rspace,istart=istart,icount=icount)
      xlon(jde1:jde2,ide1:ide2) = rspace
      call read_var2d_static(idmin,'dlat',rspace,istart=istart,icount=icount)
      dlat(jde1:jde2,ide1:ide2) = rspace
      call read_var2d_static(idmin,'dlon',rspace,istart=istart,icount=icount)
      dlon(jde1:jde2,ide1:ide2) = rspace
      if ( idynamic == 3 ) then
        call read_var2d_static(idmin,'ulat',rspace,istart=istart,icount=icount)
        ulat(jde1:jde2,ide1:ide2) = rspace
        call read_var2d_static(idmin,'ulon',rspace,istart=istart,icount=icount)
        ulon(jde1:jde2,ide1:ide2) = rspace
        call read_var2d_static(idmin,'vlat',rspace,istart=istart,icount=icount)
        vlat(jde1:jde2,ide1:ide2) = rspace
        call read_var2d_static(idmin,'vlon',rspace,istart=istart,icount=icount)
        vlon(jde1:jde2,ide1:ide2) = rspace
        call read_var2d_static(idmin,'umap',rspace,istart=istart,icount=icount)
        msfu(jde1:jde2,ide1:ide2) = rspace
        call read_var2d_static(idmin,'vmap',rspace,istart=istart,icount=icount)
        msfv(jde1:jde2,ide1:ide2) = rspace
        call read_var2d_static(idmin,'xmap',rspace,istart=istart,icount=icount)
        msfx(jde1:jde2,ide1:ide2) = rspace
        call read_var2d_static(idmin,'topou',rspace,istart=istart,icount=icount)
        htu(jde1:jde2,ide1:ide2) = rspace
        call read_var2d_static(idmin,'topov',rspace,istart=istart,icount=icount)
        htv(jde1:jde2,ide1:ide2) = rspace
      else
        call read_var2d_static(idmin,'xmap',rspace,istart=istart,icount=icount)
        msfx(jde1:jde2,ide1:ide2) = rspace
        call read_var2d_static(idmin,'dmap',rspace,istart=istart,icount=icount)
        msfd(jde1:jde2,ide1:ide2) = rspace
      end if
      call read_var2d_static(idmin,'topo',rspace,istart=istart,icount=icount)
      if ( ensemble_run ) then
        if ( myid == italk ) then
          write(stdout,*) 'Applying perturbation to input dataset:'
        end if
        if ( lperturb_topo ) then
          if ( idynamic == 3 ) then
            write(stderr, *) 'Cannot apply perturbation to topography ', &
                'in model for dynamic == 3.'
            write(stderr, *) 'Apply it to the netcdf DOMAIN file.'
            call fatal(__FILE__,__LINE__,'DOMAIN READ')
          end if
          if ( myid == italk ) then
            write(stdout,'(a,f7.2,a)') 'Topo with value ', &
              perturb_frac_topo*d_100,'%'
          end if
          call randify(rspace,perturb_frac_topo,icount(1),icount(2))
        end if
      end if
      ht(jde1:jde2,ide1:ide2) = rspace
      call read_var2d_static(idmin,'areacella',rspace, &
                             istart=istart,icount=icount)
      area(jde1:jde2,ide1:ide2) = rspace
      call read_var2d_static(idmin,'mask',rspace,istart=istart,icount=icount)
      mask(jde1:jde2,ide1:ide2) = rspace
      call read_var2d_static(idmin,'landuse',rspace,istart=istart,icount=icount)
      lnd(jde1:jde2,ide1:ide2) = rspace
      call read_var2d_static(idmin,'texture',rspace,istart=istart,icount=icount)
      tex(jde1:jde2,ide1:ide2) = rspace
      call read_var2d_static(idmin,'coriol',rspace,istart=istart,icount=icount)
      coriol(jde1:jde2,ide1:ide2) = rspace
      call read_var2d_static(idmin,'smoist',rspace,istart=istart,icount=icount)
      smoist(jde1:jde2,ide1:ide2) = rspace
      if ( replacemoist ) then
        istart3(1:2) = istart
        icount3(1:2) = icount
        icount3(3) = 1
        do ilev = 1, num_soil_layers
          istart3(3) = ilev
          call read_var2d_static(idmin,'rmoist',rspace, &
                                 istart=istart3,icount=icount3)
          rmoist(jde1:jde2,ide1:ide2,ilev) = rspace
        end do
      end if
      if ( replacetemp ) then
        istart3(1:2) = istart
        icount3(1:2) = icount
        icount3(3) = 1
        do ilev = 1, num_soil_layers
          istart3(3) = ilev
          call read_var2d_static(idmin,'rts',rspace, &
                                 istart=istart3,icount=icount3)
          rts(jde1:jde2,ide1:ide2,ilev) = rspace
        end do
      end if
      rspace = d_zero
      call read_var2d_static(idmin,'snowam',rspace,has_snow, &
           istart=istart,icount=icount)
      if ( has_snow ) snowam(jde1:jde2,ide1:ide2) = rspace
      if ( lakemod == 1 ) then
        call read_var2d_static(idmin,'dhlake',rspace,has_dhlake, &
             istart=istart,icount=icount)
        if (has_dhlake) hlake(jde1:jde2,ide1:ide2) = rspace
      end if
      if ( idynamic == 2 ) then
        call read_reference_surface_temp(idmin,ts0)
      end if
      call closefile(idmin)
      deallocate(rspace)
    else
      if ( myid == iocpu ) then
        istart(1) = 1
        istart(2) = 1
        icount(1) = jx
        icount(2) = iy
        allocate(rspace(icount(1),icount(2)))
        call openfile_withname(dname,idmin)
        call check_domain(idmin)
        if ( idynamic == 2 ) then
          call read_reference_surface_temp(idmin,ts0)
        end if
        call bcast(ts0)
        call get_attribute(idmin,'initialized_soil_moisture',csmoist)
        if ( csmoist == 'Yes' ) then
          replacemoist = .true.
        else
          replacemoist = .false.
        end if
        call bcast(replacemoist)
        call get_attribute(idmin,'initialized_soil_moisture',cstemp)
        if ( cstemp == 'Yes' ) then
          replacetemp = .true.
        else
          replacetemp = .false.
        end if
        call bcast(replacetemp)
        lerror = .true.
        call read_var1d_static(idmin,'kz',sigma,lerror)
        if ( .not. lerror ) then
          call read_var1d_static(idmin,'sigma',sigma)
        end if
        call bcast(sigma)
        if ( idynamic == 3 ) then
          call read_var1d_static(idmin,'rlat',rlat,lerror)
          call read_var1d_static(idmin,'rlon',rlon,lerror)
          call bcast(rlat)
          call bcast(rlon)
        end if
        call read_var2d_static(idmin,'xlat',rspace,istart=istart,icount=icount)
        call grid_distribute(rspace,xlat,jde1,jde2,ide1,ide2)
        call read_var2d_static(idmin,'xlon',rspace,istart=istart,icount=icount)
        call grid_distribute(rspace,xlon,jde1,jde2,ide1,ide2)
        call read_var2d_static(idmin,'dlat',rspace,istart=istart,icount=icount)
        call grid_distribute(rspace,dlat,jde1,jde2,ide1,ide2)
        call read_var2d_static(idmin,'dlon',rspace,istart=istart,icount=icount)
          call grid_distribute(rspace,dlon,jde1,jde2,ide1,ide2)
        if ( idynamic == 3 ) then
          call read_var2d_static(idmin,'ulat',rspace, &
                                 istart=istart,icount=icount)
          call grid_distribute(rspace,ulat,jde1,jde2,ide1,ide2)
          call read_var2d_static(idmin,'ulon',rspace, &
                                 istart=istart,icount=icount)
          call grid_distribute(rspace,ulon,jde1,jde2,ide1,ide2)
          call read_var2d_static(idmin,'vlat',rspace, &
                                 istart=istart,icount=icount)
          call grid_distribute(rspace,vlat,jde1,jde2,ide1,ide2)
          call read_var2d_static(idmin,'vlon',rspace, &
                                 istart=istart,icount=icount)
          call grid_distribute(rspace,vlon,jde1,jde2,ide1,ide2)
          call read_var2d_static(idmin,'umap',rspace, &
                                 istart=istart,icount=icount)
          call grid_distribute(rspace,msfu,jde1,jde2,ide1,ide2)
          call read_var2d_static(idmin,'vmap',rspace, &
                                 istart=istart,icount=icount)
          call grid_distribute(rspace,msfv,jde1,jde2,ide1,ide2)
          call read_var2d_static(idmin,'xmap',rspace, &
                                 istart=istart,icount=icount)
          call grid_distribute(rspace,msfx,jde1,jde2,ide1,ide2)
          call read_var2d_static(idmin,'topou',rspace, &
                                 istart=istart,icount=icount)
          call grid_distribute(rspace,htu,jde1,jde2,ide1,ide2)
          call read_var2d_static(idmin,'topov',rspace, &
                                 istart=istart,icount=icount)
          call grid_distribute(rspace,htv,jde1,jde2,ide1,ide2)
        else
          call read_var2d_static(idmin,'xmap',rspace, &
                                 istart=istart,icount=icount)
          call grid_distribute(rspace,msfx,jde1,jde2,ide1,ide2)
          call read_var2d_static(idmin,'dmap',rspace, &
                                 istart=istart,icount=icount)
          call grid_distribute(rspace,msfd,jde1,jde2,ide1,ide2)
        end if
        call read_var2d_static(idmin,'topo',rspace,istart=istart,icount=icount)
        if ( ensemble_run ) then
          if ( idynamic == 3 ) then
            write(stderr, *) 'Cannot apply perturbation to topography ', &
                'in model for dynamic == 3.'
            write(stderr, *) 'Apply it to the netcdf DOMAIN file.'
            call fatal(__FILE__,__LINE__,'DOMAIN READ')
          end if
          if ( myid == italk ) then
            write(stdout,*) 'Applying perturbation to input dataset:'
          end if
          if ( lperturb_topo ) then
            if ( myid == italk ) then
              write(stdout,'(a,f7.2,a)') 'Topo with value ', &
                perturb_frac_topo*d_100,'%'
            end if
            call randify(rspace,perturb_frac_topo,icount(1),icount(2))
          end if
        end if
        call grid_distribute(rspace,ht,jde1,jde2,ide1,ide2)
        call read_var2d_static(idmin,'areacella',rspace, &
                               istart=istart,icount=icount)
        call grid_distribute(rspace,area,jde1,jde2,ide1,ide2)
        call read_var2d_static(idmin,'mask',rspace,istart=istart,icount=icount)
        call grid_distribute(rspace,mask,jde1,jde2,ide1,ide2)
        call read_var2d_static(idmin,'landuse',rspace, &
                istart=istart,icount=icount)
        call grid_distribute(rspace,lnd,jde1,jde2,ide1,ide2)
        call read_var2d_static(idmin,'texture',rspace, &
                istart=istart,icount=icount)
        call grid_distribute(rspace,tex,jde1,jde2,ide1,ide2)
        call read_var2d_static(idmin,'coriol',rspace, &
                istart=istart,icount=icount)
        call grid_distribute(rspace,coriol,jde1,jde2,ide1,ide2)
        call read_var2d_static(idmin,'smoist',rspace, &
                istart=istart,icount=icount)
        call grid_distribute(rspace,smoist,jde1,jde2,ide1,ide2)
        if ( replacemoist ) then
          allocate(tempmoist(jde1:jde2,ide1:ide2))
          istart3(1:2) = istart
          icount3(1:2) = icount
          icount3(3) = 1
          do ilev = 1, num_soil_layers
            istart3(3) = ilev
            call read_var2d_static(idmin,'rmoist',rspace, &
                    istart=istart3,icount=icount3)
            call grid_distribute(rspace,tempmoist,jde1,jde2,ide1,ide2)
            rmoist(jde1:jde2,ide1:ide2,ilev) = tempmoist
          end do
          deallocate(tempmoist)
        end if
        if ( replacetemp ) then
          allocate(tempmoist(jde1:jde2,ide1:ide2))
          istart3(1:2) = istart
          icount3(1:2) = icount
          icount3(3) = 1
          do ilev = 1, num_soil_layers
            istart3(3) = ilev
            call read_var2d_static(idmin,'rts',rspace, &
                    istart=istart3,icount=icount3)
            call grid_distribute(rspace,tempmoist,jde1,jde2,ide1,ide2)
            rts(jde1:jde2,ide1:ide2,ilev) = tempmoist
          end do
          deallocate(tempmoist)
        end if
        rspace = d_zero
        call read_var2d_static(idmin,'snowam',rspace,has_snow, &
                istart=istart,icount=icount)
        call grid_distribute(rspace,snowam,jde1,jde2,ide1,ide2)
        if ( lakemod == 1 ) then
          call read_var2d_static(idmin,'dhlake',rspace,has_dhlake, &
                   istart=istart,icount=icount)
          call grid_distribute(rspace,hlake,jde1,jde2,ide1,ide2)
        end if
        call closefile(idmin)
        deallocate(rspace)
      else
        call bcast(ts0)
        call bcast(replacemoist)
        call bcast(replacetemp)
        call bcast(sigma)
        if ( idynamic == 3 ) then
          call bcast(rlat)
          call bcast(rlon)
        end if
        call grid_distribute(rspace,xlat,jde1,jde2,ide1,ide2)
        call grid_distribute(rspace,xlon,jde1,jde2,ide1,ide2)
        call grid_distribute(rspace,dlat,jde1,jde2,ide1,ide2)
        call grid_distribute(rspace,dlon,jde1,jde2,ide1,ide2)
        if ( idynamic == 3 ) then
          call grid_distribute(rspace,ulat,jde1,jde2,ide1,ide2)
          call grid_distribute(rspace,ulon,jde1,jde2,ide1,ide2)
          call grid_distribute(rspace,vlat,jde1,jde2,ide1,ide2)
          call grid_distribute(rspace,vlon,jde1,jde2,ide1,ide2)
          call grid_distribute(rspace,msfu,jde1,jde2,ide1,ide2)
          call grid_distribute(rspace,msfv,jde1,jde2,ide1,ide2)
          call grid_distribute(rspace,msfx,jde1,jde2,ide1,ide2)
          call grid_distribute(rspace,htu,jde1,jde2,ide1,ide2)
          call grid_distribute(rspace,htv,jde1,jde2,ide1,ide2)
        else
          call grid_distribute(rspace,msfx,jde1,jde2,ide1,ide2)
          call grid_distribute(rspace,msfd,jde1,jde2,ide1,ide2)
        end if
        call grid_distribute(rspace,ht,jde1,jde2,ide1,ide2)
        call grid_distribute(rspace,area,jde1,jde2,ide1,ide2)
        call grid_distribute(rspace,mask,jde1,jde2,ide1,ide2)
        call grid_distribute(rspace,lnd,jde1,jde2,ide1,ide2)
        call grid_distribute(rspace,tex,jde1,jde2,ide1,ide2)
        call grid_distribute(rspace,coriol,jde1,jde2,ide1,ide2)
        call grid_distribute(rspace,smoist,jde1,jde2,ide1,ide2)
        if ( replacemoist ) then
          allocate(tempmoist(jde1:jde2,ide1:ide2))
          do ilev = 1, num_soil_layers
            call grid_distribute(rspace,tempmoist,jde1,jde2,ide1,ide2)
            rmoist(jde1:jde2,ide1:ide2,ilev) = tempmoist
          end do
          deallocate(tempmoist)
        end if
        if ( replacetemp ) then
          allocate(tempmoist(jde1:jde2,ide1:ide2))
          do ilev = 1, num_soil_layers
            call grid_distribute(rspace,tempmoist,jde1,jde2,ide1,ide2)
            rts(jde1:jde2,ide1:ide2,ilev) = tempmoist
          end do
          deallocate(tempmoist)
        end if
        call grid_distribute(rspace,snowam,jde1,jde2,ide1,ide2)
        if ( lakemod == 1 ) then
          call grid_distribute(rspace,hlake,jde1,jde2,ide1,ide2)
        end if
      end if
    end if
    if ( idynamic == 2 ) then
      if ( .not. do_parallel_netcdf_in ) then
        call getmem(tempw,jce1,jce2,ice1,ice2,1,kz,'read_domain:tempw')
        call getmem(tempwtop,jce1,jce2,ice1,ice2,'read_domain:tempwtop')
      end if
    end if
  end subroutine read_domain_info

  subroutine read_subdomain_info(ht,lnd,tex,mask,area,xlat,xlon,hlake)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: ht
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: lnd
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: tex
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: mask
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: area
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: xlat
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: xlon
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: hlake
    character(len=256) :: dname
    integer(ik4) :: idmin
    integer(ik4), dimension(2) :: istart, icount
    real(rkx), dimension(:,:), pointer, contiguous :: rspace
    real(rkx), dimension(:,:,:), pointer, contiguous :: rspace0
    logical :: has_dhlake
    character(len=3) :: sbstring

    call outstream_async_flush()

    has_dhlake = .true.
    write (sbstring,'(i0.3)') nsg
    dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN'//sbstring//'.nc'
    if ( myid == italk ) then
      write(stdout,*) 'Reading Sub-Domain file : ', trim(dname)
    end if

    if ( do_parallel_netcdf_in ) then
      istart(1) = jde1sg
      istart(2) = ide1sg
      icount(1) = jde2sg-jde1sg+1
      icount(2) = ide2sg-ide1sg+1
      allocate(rspace(icount(1),icount(2)))
      call openfile_withname(dname,idmin)
      call read_var2d_static(idmin,'xlat',rspace,istart=istart,icount=icount)
      call input_reorder(rspace,xlat,jde1,jde2,ide1,ide2)
      call read_var2d_static(idmin,'xlon',rspace,istart=istart,icount=icount)
      call input_reorder(rspace,xlon,jde1,jde2,ide1,ide2)
      call read_var2d_static(idmin,'topo',rspace,istart=istart,icount=icount)
      if ( ensemble_run ) then
        if ( myid == italk ) then
          write(stdout,*) 'Applying perturbation to input dataset:'
        end if
        if ( lperturb_topo ) then
          if ( myid == italk ) then
            write(stdout,'(a,f7.2,a)') 'Topo with value ', &
              perturb_frac_topo*d_100,'%'
          end if
          call randify(rspace,perturb_frac_topo,icount(1),icount(2))
        end if
      end if
      call input_reorder(rspace,ht,jde1,jde2,ide1,ide2)
      call read_var2d_static(idmin,'areacella',rspace, &
                             istart=istart,icount=icount)
      call input_reorder(rspace,area,jde1,jde2,ide1,ide2)
      call read_var2d_static(idmin,'mask',rspace,istart=istart,icount=icount)
      call input_reorder(rspace,mask,jde1,jde2,ide1,ide2)
      call read_var2d_static(idmin,'landuse',rspace,istart=istart,icount=icount)
      call input_reorder(rspace,lnd,jde1,jde2,ide1,ide2)
      call read_var2d_static(idmin,'texture',rspace,istart=istart,icount=icount)
      call input_reorder(rspace,tex,jde1,jde2,ide1,ide2)
      if ( lakemod == 1 ) then
        call read_var2d_static(idmin,'dhlake',rspace,has_dhlake, &
                               istart=istart,icount=icount)
        if ( has_dhlake ) call input_reorder(rspace,hlake,jde1,jde2,ide1,ide2)
      end if
      call closefile(idmin)
      deallocate(rspace)
    else
      if ( myid == iocpu ) then
        istart(1) = 1
        istart(2) = 1
        icount(1) = jxsg
        icount(2) = iysg
        allocate(rspace(icount(1),icount(2)))
        allocate(rspace0(nnsg,jx,iy))
        call openfile_withname(dname,idmin)
        call read_var2d_static(idmin,'xlat',rspace,istart=istart,icount=icount)
        call input_reorder(rspace,rspace0,1,jx,1,iy)
        call subgrid_distribute(rspace0,xlat,jde1,jde2,ide1,ide2)
        call read_var2d_static(idmin,'xlon',rspace,istart=istart,icount=icount)
        call input_reorder(rspace,rspace0,1,jx,1,iy)
        call subgrid_distribute(rspace0,xlon,jde1,jde2,ide1,ide2)
        call read_var2d_static(idmin,'topo',rspace,istart=istart,icount=icount)
        if ( ensemble_run ) then
          if ( myid == italk ) then
            write(stdout,*) 'Applying perturbation to input dataset:'
          end if
          if ( lperturb_topo ) then
            if ( myid == italk ) then
              write(stdout,'(a,f7.2,a)') 'Topo with value ', &
                perturb_frac_topo*d_100,'%'
            end if
            call randify(rspace,perturb_frac_topo,icount(1),icount(2))
          end if
        end if
        call input_reorder(rspace,rspace0,1,jx,1,iy)
        call subgrid_distribute(rspace0,ht,jde1,jde2,ide1,ide2)
        call read_var2d_static(idmin,'areacella',rspace, &
                               istart=istart,icount=icount)
        call input_reorder(rspace,rspace0,1,jx,1,iy)
        call subgrid_distribute(rspace0,area,jde1,jde2,ide1,ide2)
        call read_var2d_static(idmin,'mask',rspace,istart=istart,icount=icount)
        call input_reorder(rspace,rspace0,1,jx,1,iy)
        call subgrid_distribute(rspace0,mask,jde1,jde2,ide1,ide2)
        call read_var2d_static(idmin,'landuse', &
                rspace,istart=istart,icount=icount)
        call input_reorder(rspace,rspace0,1,jx,1,iy)
        call subgrid_distribute(rspace0,lnd,jde1,jde2,ide1,ide2)
        call read_var2d_static(idmin,'texture', &
                rspace,istart=istart,icount=icount)
        call input_reorder(rspace,rspace0,1,jx,1,iy)
        call subgrid_distribute(rspace0,tex,jde1,jde2,ide1,ide2)
        if ( lakemod == 1 ) then
          call read_var2d_static(idmin,'dhlake',rspace,has_dhlake, &
                                 istart=istart,icount=icount)
          if ( has_dhlake ) call input_reorder(rspace,rspace0,1,jx,1,iy)
          call subgrid_distribute(rspace0,hlake,jci1,jci2,ici1,ici2)
        end if
        call closefile(idmin)
        deallocate(rspace)
        deallocate(rspace0)
      else
        call subgrid_distribute(rspace0,xlat,jde1,jde2,ide1,ide2)
        call subgrid_distribute(rspace0,xlon,jde1,jde2,ide1,ide2)
        call subgrid_distribute(rspace0,ht,jde1,jde2,ide1,ide2)
        call subgrid_distribute(rspace0,area,jde1,jde2,ide1,ide2)
        call subgrid_distribute(rspace0,mask,jde1,jde2,ide1,ide2)
        call subgrid_distribute(rspace0,lnd,jde1,jde2,ide1,ide2)
        call subgrid_distribute(rspace0,tex,jde1,jde2,ide1,ide2)
        if ( lakemod == 1 ) then
          call subgrid_distribute(rspace0,hlake,jci1,jci2,ici1,ici2)
        end if
      end if
    end if
  end subroutine read_subdomain_info

  integer(ik4) function clmbc_search(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    type(rcm_time_interval) :: tdif
    character(len=32) :: appdat1, appdat2
    if ( .not. do_parallel_netcdf_in ) then
      if ( myid /= iocpu ) then
        clmbc_search = 1
        return
      end if
    end if
    if (idate > clmbc_idate(ibcnrec) .or. idate < clmbc_idate(1)) then
      clmbc_search = -1
    else
      tdif = idate-clmbc_idate(1)
      clmbcrec = nint(tohours(tdif))+1
      if ( clmbcrec < 1 .or. clmbcrec > clmbcnrec ) then
        appdat1 = tochar(idate)
        write (stderr,*) 'Record is not found in SFBC file for ',appdat1
        appdat1 = tochar(clmbc_idate(1))
        appdat2 = tochar(clmbc_idate(clmbcnrec))
        write (stderr,*) 'Range is : ', appdat1, '-', appdat2
        call fatal(__FILE__,__LINE__,'SFBC READ')
      end if
      clmbc_search = clmbcrec
    end if
  end function clmbc_search

  integer(ik4) function icbc_search(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    type(rcm_time_interval) :: tdif
    character(len=32) :: appdat1, appdat2
    if ( .not. do_parallel_netcdf_in ) then
      if ( myid /= iocpu ) then
        icbc_search = 1
        return
      end if
    end if
    if (idate > icbc_idate(ibcnrec) .or. idate < icbc_idate(1)) then
      icbc_search = -1
    else
      tdif = idate-icbc_idate(1)
      ibcrec = (nint(tohours(tdif))/ibdyfrq)+1
      if ( ibcrec < 1 .or. ibcrec > ibcnrec ) then
        appdat1 = tochar(idate)
        write (stderr,*) 'Record is not found in ICBC file for ',appdat1
        appdat1 = tochar(icbc_idate(1))
        appdat2 = tochar(icbc_idate(ibcnrec))
        write (stderr,*) 'Range is : ', appdat1, '-', appdat2
        call fatal(__FILE__,__LINE__,'ICBC READ')
      end if
      icbc_search = ibcrec
    end if
  end function icbc_search

#ifdef ASYNC_NETCDF
  subroutine warmup_icbc_prefetch()
    implicit none
    integer(c_int64_t) :: required_bytes, configured_bytes
    logical :: async_ready

#ifdef PNETCDF
    return
#endif
    if ( .not. do_parallel_netcdf_in ) return
    if ( ibcin < 0 .or. .not. allocated(icbc_idate) ) return

    required_bytes = icbc_prefetch_required_bytes()
    if ( required_bytes <= 0_c_int64_t ) return

    if ( myid == iocpu ) then
      async_ready = async_netcdf_initialize()
    else
      configured_bytes = async_netcdf_configured_cap_bytes()
      if ( configured_bytes < required_bytes ) return
      async_ready = async_netcdf_initialize(required_bytes)
    end if
    if ( .not. async_ready ) return
  end subroutine warmup_icbc_prefetch

  subroutine prefetch_icbc(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    integer(ik4), dimension(4) :: istart, icount
    integer(ik4) :: rec, stat, pending, pstatus
    integer(ik4) :: n2, n3
    integer(c_int64_t) :: required_bytes, configured_bytes
    logical :: failed
    logical :: async_ready

#ifdef PNETCDF
    return
#endif
    if ( .not. do_parallel_netcdf_in ) return
    if ( ibcin < 0 .or. .not. allocated(icbc_idate) ) return

    if ( icbc_prefetch%active ) then
      if ( icbc_prefetch%target_date == idate ) return
      call async_netcdf_counter_poll(icbc_prefetch%counter,pending,pstatus)
      if ( pending > 0 ) return
      call release_icbc_prefetch(.false.)
    end if

    rec = icbc_record_for_date(idate)
    if ( rec < 1 ) return

    required_bytes = icbc_prefetch_required_bytes()
    if ( required_bytes <= 0_c_int64_t ) return
    configured_bytes = async_netcdf_configured_cap_bytes()
    if ( configured_bytes < required_bytes ) return

    if ( myid == iocpu ) then
      async_ready = async_netcdf_initialize()
    else
      async_ready = async_netcdf_initialize(required_bytes)
    end if
    if ( .not. async_ready ) return

    n2 = (jde2-jde1+1)*(ide2-ide1+1)
    n3 = n2*kz

    call async_netcdf_try_acquire_buffer(icbc_prefetch%buffer, &
      required_bytes,stat)
    if ( stat /= nf90_noerr ) return

    call async_netcdf_counter_reset(icbc_prefetch%counter)
    icbc_prefetch%active = .true.
    icbc_prefetch%target_date = idate
    icbc_prefetch%ncid = ibcin
    icbc_prefetch%record = rec
    icbc_prefetch%nitems = 0
    call map_icbc_prefetch_buffers(n2,n3)

    failed = .false.
    istart(1) = jde1
    istart(2) = ide1
    istart(3) = rec
    icount(1) = jde2-jde1+1
    icount(2) = ide2-ide1+1
    icount(3) = 1
    call enqueue_icbc_prefetch_read(icbc_ivar(1),icbc_prefetch%ps_flat, &
      istart(1:3),icount(1:3),failed)
    call enqueue_icbc_prefetch_read(icbc_ivar(2),icbc_prefetch%ts_flat, &
      istart(1:3),icount(1:3),failed)

    istart(1) = jde1
    istart(2) = ide1
    istart(3) = 1
    istart(4) = rec
    icount(1) = jde2-jde1+1
    icount(2) = ide2-ide1+1
    icount(3) = kz
    icount(4) = 1
    call enqueue_icbc_prefetch_read(icbc_ivar(3),icbc_prefetch%u_flat, &
      istart,icount,failed)
    call enqueue_icbc_prefetch_read(icbc_ivar(4),icbc_prefetch%v_flat, &
      istart,icount,failed)
    call enqueue_icbc_prefetch_read(icbc_ivar(5),icbc_prefetch%t_flat, &
      istart,icount,failed)
    call enqueue_icbc_prefetch_read(icbc_ivar(6),icbc_prefetch%qv_flat, &
      istart,icount,failed)
    if ( has_qc ) then
      call enqueue_icbc_prefetch_read(icbc_ivar(10),icbc_prefetch%qc_flat, &
        istart,icount,failed)
    end if
    if ( has_qi ) then
      call enqueue_icbc_prefetch_read(icbc_ivar(11),icbc_prefetch%qi_flat, &
        istart,icount,failed)
    end if
    if ( idynamic == 2 ) then
      call enqueue_icbc_prefetch_read(icbc_ivar(7),icbc_prefetch%pp_flat, &
        istart,icount,failed)
      call enqueue_icbc_prefetch_read(icbc_ivar(8),icbc_prefetch%ww_flat, &
        istart,icount,failed)
      istart(1) = jde1
      istart(2) = ide1
      istart(3) = rec
      icount(1) = jde2-jde1+1
      icount(2) = ide2-ide1+1
      icount(3) = 1
      call enqueue_icbc_prefetch_read(icbc_ivar(9), &
        icbc_prefetch%wtop_flat,istart(1:3),icount(1:3),failed)
    else if ( idynamic == 3 ) then
      call enqueue_icbc_prefetch_read(icbc_ivar(7),icbc_prefetch%pai_flat, &
        istart,icount,failed)
    end if
    if ( failed .and. icbc_prefetch%nitems == 0 ) then
      call release_icbc_prefetch(.false.)
    end if
  end subroutine prefetch_icbc

  subroutine consume_icbc_prefetch(idate,ps,ts,ilnd,u,v,t,qv,qc,qi,pp,ww, &
      pai,consumed)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: ps
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: ts
    integer(ik4), pointer, contiguous, dimension(:,:), intent(in) :: ilnd
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: u
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: v
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: t
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: qv
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: qc
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: qi
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: pp
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: ww
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: pai
    logical, intent(out) :: consumed
    integer(ik4) :: pending, status

    consumed = .false.
    if ( .not. icbc_prefetch%active ) return
    if ( icbc_prefetch%target_date /= idate ) then
      call async_netcdf_counter_poll(icbc_prefetch%counter,pending,status)
      if ( pending == 0 ) call release_icbc_prefetch(.false.)
      return
    end if

    call async_netcdf_counter_wait(icbc_prefetch%counter,status)
    if ( status == nf90_noerr ) then
      ibcrec = icbc_prefetch%record
      call copy_prefetched_icbc(ps,ts,ilnd,u,v,t,qv,qc,qi,pp,ww,pai)
      consumed = .true.
    end if
    call release_icbc_prefetch(.false.)
  end subroutine consume_icbc_prefetch

  integer(ik4) function icbc_record_for_date(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    type(rcm_time_interval) :: tdif

    icbc_record_for_date = -1
    if ( .not. allocated(icbc_idate) ) return
    if ( ibcnrec < 1 ) return
    if ( idate > icbc_idate(ibcnrec) .or. idate < icbc_idate(1) ) return
    tdif = idate-icbc_idate(1)
    icbc_record_for_date = (nint(tohours(tdif))/ibdyfrq)+1
    if ( icbc_record_for_date < 1 .or. &
         icbc_record_for_date > ibcnrec ) icbc_record_for_date = -1
  end function icbc_record_for_date

  subroutine release_icbc_prefetch(wait)
    implicit none
    logical, intent(in) :: wait
    integer(ik4) :: status

    if ( .not. icbc_prefetch%active ) return
    if ( wait ) call async_netcdf_counter_wait(icbc_prefetch%counter,status)
    call async_netcdf_release_buffer(icbc_prefetch%buffer)
    call clear_icbc_prefetch_pointers()
    icbc_prefetch%active = .false.
    icbc_prefetch%ncid = -1
    icbc_prefetch%record = -1
    icbc_prefetch%nitems = 0
    call async_netcdf_counter_reset(icbc_prefetch%counter)
  end subroutine release_icbc_prefetch

  subroutine enqueue_icbc_prefetch_read(varid,values,start,count,failed)
    implicit none
    integer(ik4), intent(in) :: varid
    real(rkx), pointer, contiguous, dimension(:), intent(inout) :: values
    integer(ik4), dimension(:), intent(in) :: start, count
    logical, intent(inout) :: failed
    integer(ik4) :: stat

    if ( failed ) return
    stat = async_netcdf_get_var_rkx(icbc_prefetch%ncid,varid,values,start, &
      count,icbc_prefetch%counter)
    if ( stat == nf90_noerr ) then
      icbc_prefetch%nitems = icbc_prefetch%nitems + 1
    else
      failed = .true.
      call async_netcdf_counter_fail(icbc_prefetch%counter,stat)
    end if
  end subroutine enqueue_icbc_prefetch_read

  subroutine map_icbc_prefetch_buffers(n2,n3)
    implicit none
    integer(ik4), intent(in) :: n2, n3
    integer(ik4) :: offset

    call clear_icbc_prefetch_pointers()
    offset = 1
    call map_prefetch_2d(icbc_prefetch%ps_flat,icbc_prefetch%ps,offset,n2)
    call map_prefetch_2d(icbc_prefetch%ts_flat,icbc_prefetch%ts,offset,n2)
    call map_prefetch_3d(icbc_prefetch%u_flat,icbc_prefetch%u,offset,n3)
    call map_prefetch_3d(icbc_prefetch%v_flat,icbc_prefetch%v,offset,n3)
    call map_prefetch_3d(icbc_prefetch%t_flat,icbc_prefetch%t,offset,n3)
    call map_prefetch_3d(icbc_prefetch%qv_flat,icbc_prefetch%qv,offset,n3)
    if ( has_qc ) then
      call map_prefetch_3d(icbc_prefetch%qc_flat,icbc_prefetch%qc,offset,n3)
    end if
    if ( has_qi ) then
      call map_prefetch_3d(icbc_prefetch%qi_flat,icbc_prefetch%qi,offset,n3)
    end if
    if ( idynamic == 2 ) then
      call map_prefetch_3d(icbc_prefetch%pp_flat,icbc_prefetch%pp,offset,n3)
      call map_prefetch_3d(icbc_prefetch%ww_flat,icbc_prefetch%ww,offset,n3)
      call map_prefetch_2d(icbc_prefetch%wtop_flat,icbc_prefetch%wtop, &
        offset,n2)
    else if ( idynamic == 3 ) then
      call map_prefetch_3d(icbc_prefetch%pai_flat,icbc_prefetch%pai,offset,n3)
    end if
  end subroutine map_icbc_prefetch_buffers

  subroutine map_prefetch_2d(flat,field,offset,nvals)
    implicit none
    real(rkx), pointer, contiguous, dimension(:), intent(inout) :: flat
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: field
    integer(ik4), intent(inout) :: offset
    integer(ik4), intent(in) :: nvals

    flat(1:nvals) => icbc_prefetch%buffer%rkx(offset:offset+nvals-1)
    field(jde1:jde2,ide1:ide2) => &
      icbc_prefetch%buffer%rkx(offset:offset+nvals-1)
    offset = offset + aligned_icbc_elements(nvals)
  end subroutine map_prefetch_2d

  subroutine map_prefetch_3d(flat,field,offset,nvals)
    implicit none
    real(rkx), pointer, contiguous, dimension(:), intent(inout) :: flat
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: field
    integer(ik4), intent(inout) :: offset
    integer(ik4), intent(in) :: nvals

    flat(1:nvals) => icbc_prefetch%buffer%rkx(offset:offset+nvals-1)
    field(jde1:jde2,ide1:ide2,1:kz) => &
      icbc_prefetch%buffer%rkx(offset:offset+nvals-1)
    offset = offset + aligned_icbc_elements(nvals)
  end subroutine map_prefetch_3d

  integer(c_int64_t) function icbc_prefetch_required_bytes()
    implicit none
    integer(ik4) :: n2, n3

    n2 = (jde2-jde1+1)*(ide2-ide1+1)
    n3 = n2*kz
    icbc_prefetch_required_bytes = 2_c_int64_t*aligned_icbc_bytes(n2) + &
                                   4_c_int64_t*aligned_icbc_bytes(n3)
    if ( has_qc ) icbc_prefetch_required_bytes = &
      icbc_prefetch_required_bytes + aligned_icbc_bytes(n3)
    if ( has_qi ) icbc_prefetch_required_bytes = &
      icbc_prefetch_required_bytes + aligned_icbc_bytes(n3)
    if ( idynamic == 2 ) then
      icbc_prefetch_required_bytes = icbc_prefetch_required_bytes + &
        2_c_int64_t*aligned_icbc_bytes(n3) + aligned_icbc_bytes(n2)
    end if
  end function icbc_prefetch_required_bytes

  integer(c_int64_t) function aligned_icbc_bytes(nvals)
    implicit none
    integer(ik4), intent(in) :: nvals
    integer(c_int64_t) :: bytes

    bytes = int(nvals,c_int64_t)*icbc_rkx_bytes()
    aligned_icbc_bytes = ((bytes+icbc_prefetch_alignment-1_c_int64_t) / &
      icbc_prefetch_alignment) * icbc_prefetch_alignment
  end function aligned_icbc_bytes

  integer(ik4) function aligned_icbc_elements(nvals)
    implicit none
    integer(ik4), intent(in) :: nvals

    aligned_icbc_elements = int(aligned_icbc_bytes(nvals) / &
      icbc_rkx_bytes(),ik4)
  end function aligned_icbc_elements

  integer(c_int64_t) function icbc_rkx_bytes()
    implicit none

#ifdef SINGLE_PRECISION_REAL
    icbc_rkx_bytes = 4_c_int64_t
#else
    icbc_rkx_bytes = 8_c_int64_t
#endif
  end function icbc_rkx_bytes

  subroutine clear_icbc_prefetch_pointers()
    implicit none

    nullify(icbc_prefetch%ps_flat,icbc_prefetch%ts_flat)
    nullify(icbc_prefetch%u_flat,icbc_prefetch%v_flat)
    nullify(icbc_prefetch%t_flat,icbc_prefetch%qv_flat)
    nullify(icbc_prefetch%qc_flat,icbc_prefetch%qi_flat)
    nullify(icbc_prefetch%pp_flat,icbc_prefetch%ww_flat)
    nullify(icbc_prefetch%wtop_flat)
    nullify(icbc_prefetch%pai_flat)
    nullify(icbc_prefetch%ps,icbc_prefetch%ts)
    nullify(icbc_prefetch%u,icbc_prefetch%v)
    nullify(icbc_prefetch%t,icbc_prefetch%qv)
    nullify(icbc_prefetch%qc,icbc_prefetch%qi)
    nullify(icbc_prefetch%pp,icbc_prefetch%ww)
    nullify(icbc_prefetch%wtop)
    nullify(icbc_prefetch%pai)
  end subroutine clear_icbc_prefetch_pointers

  subroutine copy_prefetched_icbc(ps,ts,ilnd,u,v,t,qv,qc,qi,pp,ww,pai)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: ps
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: ts
    integer(ik4), pointer, contiguous, dimension(:,:), intent(in) :: ilnd
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: u
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: v
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: t
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: qv
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: qc
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: qi
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: pp
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: ww
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: pai
    real(rkx), pointer, contiguous, dimension(:,:) :: psbuf, tsbuf, wtopbuf
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ubuf, vbuf, tbuf
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qvbuf, qcbuf, qibuf
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ppbuf, wwbuf, paibuf
    integer(ik4) :: i, j, k
    real(rkx) :: told, pold, rhold, satvp, tnew, pnew

    psbuf => icbc_prefetch%ps
    tsbuf => icbc_prefetch%ts
    ubuf => icbc_prefetch%u
    vbuf => icbc_prefetch%v
    tbuf => icbc_prefetch%t
    qvbuf => icbc_prefetch%qv
    qcbuf => icbc_prefetch%qc
    qibuf => icbc_prefetch%qi
    ppbuf => icbc_prefetch%pp
    wwbuf => icbc_prefetch%ww
    wtopbuf => icbc_prefetch%wtop
    paibuf => icbc_prefetch%pai

    if ( ensemble_run ) then
      if ( lperturb_ps ) then
        if ( myid == italk ) then
          write(stdout,'(a,f7.2,a)') &
                  'PS with value ',perturb_frac_ps*d_100,'%'
        end if
        call randify(psbuf,perturb_frac_ps, &
          jde2-jde1+1,ide2-ide1+1)
      end if
    end if
    !$acc kernels deviceptr(psbuf)
    ps(jce1:jce2,ice1:ice2) = psbuf(jce1:jce2,ice1:ice2)
    !$acc end kernels

    if ( ensemble_run ) then
      if ( myid == italk ) then
        write(stdout,*) 'Applying perturbation to input dataset:'
      end if
      if ( lperturb_ts ) then
        if ( myid == italk ) then
          write(stdout,'(a,f7.2,a)') &
                  'TS with value ',perturb_frac_ts*d_100,'%'
        end if
        call randify(tsbuf,perturb_frac_ts, &
          jde2-jde1+1,ide2-ide1+1)
      end if
    end if
    !$acc kernels deviceptr(tsbuf)
    ts(jce1:jce2,ice1:ice2) = tsbuf(jce1:jce2,ice1:ice2)
    !$acc end kernels

    if ( ensemble_run ) then
      if ( lperturb_u ) then
        if ( myid == italk ) then
          write(stdout,'(a,f7.2,a)') 'U with value ',perturb_frac_u*d_100,'%'
        end if
        call randify(ubuf,perturb_frac_u, &
          jde2-jde1+1,ide2-ide1+1,kz)
      end if
    end if
    !$acc kernels deviceptr(ubuf)
    u(jde1:jde2,ide1:ide2,1:kz) = ubuf
    !$acc end kernels

    if ( ensemble_run ) then
      if ( lperturb_v ) then
        if ( myid == italk ) then
          write(stdout,'(a,f7.2,a)') 'V with value ',perturb_frac_v*d_100,'%'
        end if
        call randify(vbuf,perturb_frac_v, &
          jde2-jde1+1,ide2-ide1+1,kz)
      end if
    end if
    !$acc kernels deviceptr(vbuf)
    v(jde1:jde2,ide1:ide2,1:kz) = vbuf
    !$acc end kernels

    if ( ensemble_run ) then
      if ( lperturb_t ) then
        if ( myid == italk ) then
          write(stdout,'(a,f7.2,a)') 'T with value ',perturb_frac_t*d_100,'%'
        end if
        call randify(tbuf,perturb_frac_t, &
          jde2-jde1+1,ide2-ide1+1,kz)
      end if
    end if
    !$acc kernels deviceptr(tbuf)
    t(jce1:jce2,ice1:ice2,1:kz) = &
      tbuf(jce1:jce2,ice1:ice2,1:kz)
    !$acc end kernels

    if ( ensemble_run ) then
      if ( lperturb_q ) then
        if ( myid == italk ) then
          write(stdout,'(a,f7.2,a)') 'Q with value ',perturb_frac_q*d_100,'%'
        end if
        call randify(qvbuf,perturb_frac_q, &
          jde2-jde1+1,ide2-ide1+1,kz)
        !$acc kernels deviceptr(qvbuf)
        where ( qvbuf < 1.0e-8_rkx )
          qvbuf = 1.0e-8_rkx
        end where
        !$acc end kernels
      end if
    end if
    !$acc kernels deviceptr(qvbuf)
    qv(jce1:jce2,ice1:ice2,1:kz) = &
      qvbuf(jce1:jce2,ice1:ice2,1:kz)
    !$acc end kernels

    if ( has_qc ) then
      !$acc kernels deviceptr(qcbuf)
      qc(jce1:jce2,ice1:ice2,1:kz) = &
        qcbuf(jce1:jce2,ice1:ice2,1:kz)
      !$acc end kernels
    end if
    if ( has_qi ) then
      !$acc kernels deviceptr(qibuf)
      qi(jce1:jce2,ice1:ice2,1:kz) = &
        qibuf(jce1:jce2,ice1:ice2,1:kz)
      !$acc end kernels
    end if
    if ( idynamic == 2 ) then
      !$acc kernels deviceptr(ppbuf)
      pp(jce1:jce2,ice1:ice2,1:kz) = &
        ppbuf(jce1:jce2,ice1:ice2,1:kz)
      !$acc end kernels
      !$acc kernels deviceptr(wwbuf)
      ww(jce1:jce2,ice1:ice2,2:kzp1) = &
        wwbuf(jce1:jce2,ice1:ice2,1:kz)
      !$acc end kernels
      !$acc kernels deviceptr(wtopbuf)
      ww(jce1:jce2,ice1:ice2,1) = &
        wtopbuf(jce1:jce2,ice1:ice2)
      !$acc end kernels
    else if ( idynamic == 3 ) then
      !$acc kernels deviceptr(paibuf)
      pai(jce1:jce2,ice1:ice2,1:kz) = &
        paibuf(jce1:jce2,ice1:ice2,1:kz)
      !$acc end kernels
    end if

    if ( itweak == 1 ) then
      if ( itweak_sst == 1 ) then
        !$acc kernels
        do concurrent ( j = jce1:jce2, i = ice1:ice2 )
          if ( ilnd(j,i) == 0 ) then
            ts(j,i) = ts(j,i) + sst_tweak
          end if
        end do
        !$acc end kernels
      end if
      if ( itweak_temperature == 1 ) then
        !$acc kernels
        do k = 1, kz
          do i = ice1, ice2
            do j = jce1, jce2
              told = t(j,i,k)
              pold = sigma(k)*ps(j,i)*d_100
              satvp = pfesat(told,ps(j,i))
              rhold = max((qv(j,i,k)/(ep2*satvp/(pold-satvp))),d_zero)
              tnew = t(j,i,k) + temperature_tweak
              pnew = pold*(tnew/told)
              satvp = pfesat(tnew,ps(j,i))
              qv(j,i,k) = max(rhold*ep2*satvp/(pnew-satvp),d_zero)
              t(j,i,k) = tnew
            end do
          end do
        end do
        !$acc end kernels
      end if
    end if

  end subroutine copy_prefetched_icbc
#endif

  integer(ik4) function som_search(imon)
    implicit none
    integer(ik4), intent(in) :: imon
    if ( .not. do_parallel_netcdf_in ) then
      if ( myid /= iocpu ) then
        som_search = 1
        return
      end if
    end if
    if ( imon < 0 .or. imon > 13 ) then
      somrec = -1
    else
      if ( imon > 0 .and. imon < 13 ) then
        somrec = imon
      else if ( imon == 0 ) then
        somrec = 12
      else
        somrec = 1
      end if
    end if
    som_search = somrec
  end function som_search

  subroutine open_icbc(idate)
#ifdef OPENACC
    use cudafor
#endif
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    character(len=11) :: ctime
    integer(ik4) :: idimid, itvar, i, chkdiff
    real(rkx), dimension(:), allocatable :: icbc_nctime
    character(len=64) :: icbc_timeunits, icbc_timecal
    character(len=256) :: icbcname
    if ( .not. do_parallel_netcdf_in .and. myid /= iocpu ) then
      return
    end if
    call close_icbc
    write (ctime, '(a)') tochar10(idate)
    icbcname = trim(dirglob)//pthsep//trim(domname)// &
               '_ICBC.'//trim(ctime)//'.nc'
    call outstream_netcdf_lock()
    call openfile_withname(icbcname,ibcin)
    ibcrec = 1
    ibcnrec = 0
    call check_domain(ibcin,.true.)
    istatus = nf90_inq_dimid(ibcin, 'time', idimid)
    call check_ok(__FILE__,__LINE__,'Dimension time miss', 'ICBC FILE')
    istatus = nf90_inquire_dimension(ibcin, idimid, len=ibcnrec)
    call check_ok(__FILE__,__LINE__,'Dimension time read error', 'ICBC FILE')
    if ( ibcnrec < 1 ) then
      write (stderr,*) 'Time var in ICBC has zero dim.'
      call fatal(__FILE__,__LINE__,'ICBC READ')
    end if
    istatus = nf90_inq_varid(ibcin, 'time', itvar)
    call check_ok(__FILE__,__LINE__,'variable time miss', 'ICBC FILE')
    istatus = nf90_get_att(ibcin, itvar, 'units', icbc_timeunits)
    call check_ok(__FILE__,__LINE__,'variable time units miss','ICBC FILE')
    istatus = nf90_get_att(ibcin, itvar, 'calendar', icbc_timecal)
    call check_ok(__FILE__,__LINE__,'variable time calendar miss','ICBC FILE')
    allocate(icbc_nctime(ibcnrec), stat=istatus)
    if ( istatus /= 0 ) then
      write(stderr,*) 'Memory allocation error in ICBC for time real values'
      call fatal(__FILE__,__LINE__,'ICBC READ')
    end if
    allocate(icbc_idate(ibcnrec), stat=istatus)
    if ( istatus /= 0 ) then
      write(stderr,*) 'Memory allocation error in ICBC for time array'
      call fatal(__FILE__,__LINE__,'ICBC READ')
    end if
    istatus = nf90_get_var(ibcin, itvar, icbc_nctime)
    call check_ok(__FILE__,__LINE__,'variable time read error', 'ICBC FILE')
    do i = 1, ibcnrec
      icbc_idate(i) = timeval2date(icbc_nctime(i),icbc_timeunits,icbc_timecal)
    end do
    if ( ibcnrec > 1 ) then
      chkdiff = nint(icbc_nctime(2) - icbc_nctime(1))
      if (chkdiff /= ibdyfrq) then
        write (stderr,*) 'Time var in ICBC inconsistency.'
        write (stderr,*) 'Expecting ibdyfrq = ', ibdyfrq
        write (stderr,*) 'Found     ibdyfrq = ', chkdiff
        call fatal(__FILE__,__LINE__,'ICBC READ')
      end if
    end if
    deallocate(icbc_nctime)
    istatus = nf90_inq_varid(ibcin, 'ps', icbc_ivar(1))
    call check_ok(__FILE__,__LINE__,'variable ps miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'ts', icbc_ivar(2))
    call check_ok(__FILE__,__LINE__,'variable ts miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'u', icbc_ivar(3))
    call check_ok(__FILE__,__LINE__,'variable u miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'v', icbc_ivar(4))
    call check_ok(__FILE__,__LINE__,'variable v miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 't', icbc_ivar(5))
    call check_ok(__FILE__,__LINE__,'variable t miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'qv', icbc_ivar(6))
    call check_ok(__FILE__,__LINE__,'variable qv miss', 'ICBC FILE')
    if ( idynamic == 2 ) then
      istatus = nf90_inq_varid(ibcin, 'pp', icbc_ivar(7))
      call check_ok(__FILE__,__LINE__,'variable pp miss', 'ICBC FILE')
      istatus = nf90_inq_varid(ibcin, 'w', icbc_ivar(8))
      call check_ok(__FILE__,__LINE__,'variable w miss', 'ICBC FILE')
      istatus = nf90_inq_varid(ibcin, 'wtop', icbc_ivar(9))
      call check_ok(__FILE__,__LINE__,'variable wtop miss', 'ICBC FILE')
    else if ( idynamic == 3 ) then
      istatus = nf90_inq_varid(ibcin, 'pai', icbc_ivar(7))
      call check_ok(__FILE__,__LINE__,'variable pai miss', 'ICBC FILE')
    end if
    istatus = nf90_inq_varid(ibcin, 'qc', icbc_ivar(10))
    if ( istatus == nf90_noerr ) then
      has_qc = .true.
    end if
    istatus = nf90_inq_varid(ibcin, 'qi', icbc_ivar(11))
    if ( istatus == nf90_noerr ) then
      if ( ipptls > 1 ) then
        has_qi = .true.
      end if
    end if
    call outstream_netcdf_unlock()
    if ( do_parallel_netcdf_in ) then
#ifdef OPENACC
      size_bytes = (jde2-jde1+1)*(ide2-ide1+1)*kz*8
      if ( size_bytes > 0 ) then
        istat = cudaMallocHost(cptr, size_bytes)
        if ( istat /= cudaSuccess ) then
          write(*,*) 'cudaMallocHost failed with status: ',istat
          stop
        end if
        call c_f_pointer(cptr, rspace2_raw, [jde2-jde1+1,ide2-ide1+1])
        call c_f_pointer(cptr, rspace3_raw, [jde2-jde1+1,ide2-ide1+1,kz])
        rspace2(jde1:jde2,ide1:ide2) => rspace2_raw(:,:)
        rspace3(jde1:jde2,ide1:ide2,1:kz) => rspace3_raw(:,:,:)
      end if
#else
      allocate(rspace2(jde1:jde2,ide1:ide2))
      allocate(rspace3(jde1:jde2,ide1:ide2,kz))
#endif
    else
      allocate(rspace2(jx,iy))
      allocate(rspace3(jx,iy,kz))
    end if
  end subroutine open_icbc

  subroutine open_clmbc(idate)
    type(rcm_time_and_date), intent(in) :: idate
    character(len=11) :: ctime
    integer(ik4) :: idimid, itvar, i
    real(rkx), dimension(:), allocatable :: clmbc_nctime
    character(len=64) :: clmbc_timeunits, clmbc_timecal
    character(len=256) :: clmbcname
    if ( .not. do_parallel_netcdf_in .and. myid /= iocpu ) then
      return
    end if
    call close_clmbc
    write (ctime, '(a)') tochar10(idate)
    clmbcname = trim(dirglob)//pthsep//trim(domname)// &
               '_SFBC.'//trim(ctime)//'.nc'
    if ( myid == italk ) then
      write(stdout,*) 'Open ',trim(clmbcname)
    end if
    call outstream_netcdf_lock()
    call openfile_withname(clmbcname,clmbcin)
    clmbcrec = 1
    clmbcnrec = 0
    call check_domain(clmbcin,.true.)
    istatus = nf90_inq_dimid(clmbcin, 'time', idimid)
    call check_ok(__FILE__,__LINE__,'Dimension time miss', 'SFBC FILE')
    istatus = nf90_inquire_dimension(clmbcin, idimid, len=clmbcnrec)
    call check_ok(__FILE__,__LINE__,'Dimension time read error', 'SFBC FILE')
    if ( clmbcnrec < 1 ) then
      write (stderr,*) 'Time var in SFBC has zero dim.'
      call fatal(__FILE__,__LINE__,'SFBC READ')
    end if
    istatus = nf90_inq_varid(clmbcin, 'time', itvar)
    call check_ok(__FILE__,__LINE__,'variable time miss', 'SFBC FILE')
    istatus = nf90_get_att(clmbcin, itvar, 'units', clmbc_timeunits)
    call check_ok(__FILE__,__LINE__,'variable time units miss','SFBC FILE')
    istatus = nf90_get_att(clmbcin, itvar, 'calendar', clmbc_timecal)
    call check_ok(__FILE__,__LINE__,'variable time calendar miss','SFBC FILE')
    allocate(clmbc_nctime(clmbcnrec), stat=istatus)
    if ( istatus /= 0 ) then
      write(stderr,*) 'Memory allocation error in SFBC for time real values'
      call fatal(__FILE__,__LINE__,'SFBC READ')
    end if
    allocate(clmbc_idate(clmbcnrec), stat=istatus)
    if ( istatus /= 0 ) then
      write(stderr,*) 'Memory allocation error in SFBC for time array'
      call fatal(__FILE__,__LINE__,'SFBC READ')
    end if
    istatus = nf90_get_var(clmbcin, itvar, clmbc_nctime)
    call check_ok(__FILE__,__LINE__,'variable time read error', 'SFBC FILE')
    do i = 1, clmbcnrec
      clmbc_idate(i) = timeval2date(clmbc_nctime(i), &
                                    clmbc_timeunits,clmbc_timecal)
    end do
    deallocate(clmbc_nctime)
    istatus = nf90_inq_varid(clmbcin, 'pr', clmbc_ivar(1))
    call check_ok(__FILE__,__LINE__,'variable pr miss', 'SFBC FILE')
    istatus = nf90_inq_varid(clmbcin, 'ssr', clmbc_ivar(2))
    call check_ok(__FILE__,__LINE__,'variable ssr miss', 'SFBC FILE')
    istatus = nf90_inq_varid(clmbcin, 'strd', clmbc_ivar(3))
    call check_ok(__FILE__,__LINE__,'variable strd miss', 'SFBC FILE')
    istatus = nf90_inq_varid(clmbcin, 'clt', clmbc_ivar(4))
    call check_ok(__FILE__,__LINE__,'variable clt miss', 'SFBC FILE')
    call outstream_netcdf_unlock()
  end subroutine open_clmbc

  subroutine open_som( )
    character(len=10) :: ctime
    character(len=256) :: somname
    if ( .not. do_parallel_netcdf_in .and. myid /= iocpu ) then
      return
    end if
    call close_som
    ctime = 'YYYYMMDDHH'
    somname = trim(dirglob)//pthsep//trim(domname)//'_SOM.'//ctime//'.nc'
    call outstream_netcdf_lock()
    call openfile_withname(somname,somin)
    call check_domain(somin,.true.,.true.)
    istatus = nf90_inq_varid(somin, 'qflx', som_ivar(1))
    call check_ok(__FILE__,__LINE__,'variable qflx miss', 'SOM FILE')
    call outstream_netcdf_unlock()
    if ( do_parallel_netcdf_in ) then
      allocate(rspacesom(jci1:jci2,ici1:ici2))
    else
      allocate(rspacesom(2:jx-2,2:iy-2))
    end if
  end subroutine open_som

  subroutine read_clmbc(pr,ssr,strd,clt)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: pr
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: ssr
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: strd
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: clt

    integer(ik4), dimension(3) :: istart, icount

    if ( clmbcrec > clmbcnrec ) then
      call open_clmbc(rcmtimer%idate)
    end if
    if ( myid == italk ) then
      write(stdout,*) 'Reading SF values for ',tochar(clmbc_idate(clmbcrec))
    end if
    if ( do_parallel_netcdf_in .or. myid == iocpu ) then
      call outstream_netcdf_lock()
    end if
    if ( do_parallel_netcdf_in ) then
      istart(1) = jde1
      istart(2) = ide1
      istart(3) = clmbcrec
      icount(1) = jde2-jde1+1
      icount(2) = ide2-ide1+1
      icount(3) = 1
      istatus = nf90_get_var(clmbcin,clmbc_ivar(1),rspace2,istart,icount)
      call check_ok(__FILE__,__LINE__,'variable pr read error', 'SFBC FILE')
      pr(jci1:jci2,ici1:ici2) = rspace2(jci1:jci2,ici1:ici2)
      istatus = nf90_get_var(clmbcin,clmbc_ivar(2),rspace2,istart,icount)
      call check_ok(__FILE__,__LINE__,'variable ssr read error', 'SFBC FILE')
      ssr(jci1:jci2,ici1:ici2) = rspace2(jci1:jci2,ici1:ici2)
      istatus = nf90_get_var(clmbcin,clmbc_ivar(3),rspace2,istart,icount)
      call check_ok(__FILE__,__LINE__,'variable ssr read error', 'SFBC FILE')
      strd(jci1:jci2,ici1:ici2) = rspace2(jci1:jci2,ici1:ici2)
      istatus = nf90_get_var(clmbcin,clmbc_ivar(4),rspace2,istart,icount)
      call check_ok(__FILE__,__LINE__,'variable ssr read error', 'SFBC FILE')
      clt(jci1:jci2,ici1:ici2) = rspace2(jci1:jci2,ici1:ici2)
    else
      if ( myid == iocpu ) then
        istart(1) = 1
        istart(2) = 1
        istart(3) = clmbcrec
        icount(1) = jx
        icount(2) = iy
        icount(3) = 1
        istatus = nf90_get_var(clmbcin,clmbc_ivar(1),rspace2,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable pr read error', 'SFBC FILE')
        call grid_distribute(rspace2,pr,jci1,jci2,ici1,ici2)
        istatus = nf90_get_var(clmbcin,clmbc_ivar(2),rspace2,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable ssr read error', 'SFBC FILE')
        call grid_distribute(rspace2,ssr,jci1,jci2,ici1,ici2)
        istatus = nf90_get_var(clmbcin,clmbc_ivar(3),rspace2,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable strd read error', 'SFBC FILE')
        call grid_distribute(rspace2,strd,jci1,jci2,ici1,ici2)
        istatus = nf90_get_var(clmbcin,clmbc_ivar(4),rspace2,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable clt read error', 'SFBC FILE')
        call grid_distribute(rspace2,clt,jci1,jci2,ici1,ici2)
      else
        call grid_distribute(rspace2,pr,jci1,jci2,ici1,ici2)
        call grid_distribute(rspace2,ssr,jci1,jci2,ici1,ici2)
        call grid_distribute(rspace2,strd,jci1,jci2,ici1,ici2)
        call grid_distribute(rspace2,clt,jci1,jci2,ici1,ici2)
      end if
    end if
    if ( do_parallel_netcdf_in .or. myid == iocpu ) then
      call outstream_netcdf_unlock()
    end if
    clmbcrec = clmbcrec + 1
  end subroutine read_clmbc

  subroutine read_icbc(ps,ts,ilnd,u,v,t,qv,qc,qi,pp,ww,pai)
    !@acc use nvtx
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: ps
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: ts
    integer(ik4), pointer, contiguous, dimension(:,:), intent(in) :: ilnd
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: u
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: v
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: t
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: qv
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: qc
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: qi
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: pp
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: ww
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: pai

    integer(ik4), dimension(4) :: istart, icount
    integer(ik4) :: i, j, k
    real(rkx) :: told, pold, rhold, satvp, tnew, pnew
    !@acc call nvtxStartRange("read_icbc")
    if ( do_parallel_netcdf_in .or. myid == iocpu ) then
      call outstream_netcdf_lock()
    end if
    if ( do_parallel_netcdf_in ) then
      !@acc call nvtxStartRange("do_parallel_netcdf_in")
      istart(1) = jde1
      istart(2) = ide1
      istart(3) = ibcrec
      icount(1) = jde2-jde1+1
      icount(2) = ide2-ide1+1
      icount(3) = 1
      istatus = nf90_get_var(ibcin,icbc_ivar(1),rspace2,istart(1:3),icount(1:3))
      call check_ok(__FILE__,__LINE__,'variable ps read error', 'ICBC FILE')
      if ( ensemble_run ) then
        if ( lperturb_ps ) then
          if ( myid == italk ) then
            write(stdout,'(a,f7.2,a)') &
                    'PS with value ',perturb_frac_ps*d_100,'%'
          end if
          call randify(rspace2,perturb_frac_ps,icount(1),icount(2))
        end if
      end if
      !$acc kernels deviceptr(rspace2)
      ps(jce1:jce2,ice1:ice2) = rspace2(jce1:jce2,ice1:ice2)
      !$acc end kernels
      istatus = nf90_get_var(ibcin,icbc_ivar(2),rspace2,istart(1:3),icount(1:3))
      call check_ok(__FILE__,__LINE__,'variable ts read error', 'ICBC FILE')
      if ( ensemble_run ) then
        if ( myid == italk ) then
          write(stdout,*) 'Applying perturbation to input dataset:'
        end if
        if ( lperturb_ts ) then
          if ( myid == italk ) then
            write(stdout,'(a,f7.2,a)') &
                    'TS with value ',perturb_frac_ts*d_100,'%'
          end if
          call randify(rspace2,perturb_frac_ts,icount(1),icount(2))
        end if
      end if
      !$acc kernels deviceptr(rspace2)
      ts(jce1:jce2,ice1:ice2) = rspace2(jce1:jce2,ice1:ice2)
      !$acc end kernels
      istart(1) = jde1
      istart(2) = ide1
      istart(3) = 1
      istart(4) = ibcrec
      icount(1) = jde2-jde1+1
      icount(2) = ide2-ide1+1
      icount(3) = kz
      icount(4) = 1
      istatus = nf90_get_var(ibcin,icbc_ivar(3),rspace3,istart,icount)
      call check_ok(__FILE__,__LINE__,'variable u read error', 'ICBC FILE')
      if ( ensemble_run ) then
        if ( lperturb_u ) then
          if ( myid == italk ) then
            write(stdout,'(a,f7.2,a)') 'U with value ',perturb_frac_u*d_100,'%'
          end if
          call randify(rspace3,perturb_frac_u,icount(1),icount(2),icount(3))
        end if
      end if
      !$acc kernels deviceptr(rspace3)
      u(jde1:jde2,ide1:ide2,1:kz) = rspace3
      !$acc end kernels
      istatus = nf90_get_var(ibcin,icbc_ivar(4),rspace3,istart,icount)
      call check_ok(__FILE__,__LINE__,'variable v read error', 'ICBC FILE')
      if ( ensemble_run ) then
        if ( lperturb_v ) then
          if ( myid == italk ) then
            write(stdout,'(a,f7.2,a)') 'V with value ',perturb_frac_v*d_100,'%'
          end if
          call randify(rspace3,perturb_frac_v,icount(1),icount(2),icount(3))
        end if
      end if
      !$acc kernels deviceptr(rspace3)
      v(jde1:jde2,ide1:ide2,1:kz) = rspace3
      !$acc end kernels
      istatus = nf90_get_var(ibcin,icbc_ivar(5),rspace3,istart,icount)
      call check_ok(__FILE__,__LINE__,'variable t read error', 'ICBC FILE')
      if ( ensemble_run ) then
        if ( lperturb_t ) then
          if ( myid == italk ) then
            write(stdout,'(a,f7.2,a)') 'T with value ',perturb_frac_t*d_100,'%'
          end if
          call randify(rspace3,perturb_frac_t,icount(1),icount(2),icount(3))
        end if
      end if
      !$acc kernels deviceptr(rspace3)
      t(jce1:jce2,ice1:ice2,1:kz) = rspace3(jce1:jce2,ice1:ice2,1:kz)
      !$acc end kernels
      istatus = nf90_get_var(ibcin,icbc_ivar(6),rspace3,istart,icount)
      call check_ok(__FILE__,__LINE__,'variable qx read error', 'ICBC FILE')
      if ( ensemble_run ) then
        if ( lperturb_q ) then
          if ( myid == italk ) then
            write(stdout,'(a,f7.2,a)') 'Q with value ',perturb_frac_q*d_100,'%'
          end if
          call randify(rspace3,perturb_frac_q,icount(1),icount(2),icount(3))
          where ( rspace3 < 1.0e-8_rkx ) rspace3 = 1.0e-8_rkx
        end if
      end if
      !$acc kernels deviceptr(rspace3)
      qv(jce1:jce2,ice1:ice2,1:kz) = rspace3(jce1:jce2,ice1:ice2,1:kz)
      !$acc end kernels
      if ( has_qc ) then
        istatus = nf90_get_var(ibcin,icbc_ivar(10),rspace3,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable qc read error', 'ICBC FILE')
        qc(jce1:jce2,ice1:ice2,1:kz) = rspace3(jce1:jce2,ice1:ice2,1:kz)
      end if
      if ( has_qi ) then
        istatus = nf90_get_var(ibcin,icbc_ivar(11),rspace3,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable qi read error', 'ICBC FILE')
        qi(jce1:jce2,ice1:ice2,1:kz) = rspace3(jce1:jce2,ice1:ice2,1:kz)
      end if
      if ( idynamic == 2 ) then
        istatus = nf90_get_var(ibcin,icbc_ivar(7),rspace3,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable pp read error', 'ICBC FILE')
        !$acc kernels deviceptr(rspace3)
        pp(jce1:jce2,ice1:ice2,1:kz) = rspace3(jce1:jce2,ice1:ice2,1:kz)
        !$acc end kernels
        istatus = nf90_get_var(ibcin,icbc_ivar(8),rspace3,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable w read error', 'ICBC FILE')
        !$acc kernels deviceptr(rspace3)
        ww(jce1:jce2,ice1:ice2,2:kzp1) = rspace3(jce1:jce2,ice1:ice2,1:kz)
        !$acc end kernels
        istart(1) = jde1
        istart(2) = ide1
        istart(3) = ibcrec
        icount(1) = jde2-jde1+1
        icount(2) = ide2-ide1+1
        icount(3) = 1
        istatus = nf90_get_var(ibcin,icbc_ivar(9),rspace2, &
                               istart(1:3),icount(1:3))
        call check_ok(__FILE__,__LINE__, &
                      'variable wtop read error', 'ICBC FILE')
        !$acc kernels deviceptr(rspace2)
        ww(jce1:jce2,ice1:ice2,1) = rspace2(jce1:jce2,ice1:ice2)
        !$acc end kernels
      else if ( idynamic == 3 ) then
        istatus = nf90_get_var(ibcin,icbc_ivar(7),rspace3,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable pai read error', 'ICBC FILE')
        !$acc kernels deviceptr(rspace3)
        pai(jce1:jce2,ice1:ice2,1:kz) = rspace3(jce1:jce2,ice1:ice2,1:kz)
        !$acc end kernels
      end if
      !@acc call nvtxEndRange
    else
      if ( myid == iocpu ) then
        istart(1) = 1
        istart(2) = 1
        istart(3) = ibcrec
        icount(1) = jx
        icount(2) = iy
        icount(3) = 1
        istatus = nf90_get_var(ibcin,icbc_ivar(1), &
                rspace2,istart(1:3),icount(1:3))
        call check_ok(__FILE__,__LINE__,'variable ps read error', 'ICBC FILE')
        if ( ensemble_run ) then
          if ( lperturb_ps ) then
            if ( myid == italk ) then
              write(stdout,'(a,f7.2,a)') &
                    'PS with value ',perturb_frac_ps*d_100,'%'
            end if
            call randify(rspace2,perturb_frac_ps,icount(1),icount(2))
          end if
        end if
        call grid_distribute(rspace2,ps,jce1,jce2,ice1,ice2)
        istatus = nf90_get_var(ibcin,icbc_ivar(2), &
                rspace2,istart(1:3),icount(1:3))
        call check_ok(__FILE__,__LINE__,'variable ts read error', 'ICBC FILE')
        if ( ensemble_run ) then
          if ( myid == italk ) then
            write(stdout,*) 'Applying perturbation to input dataset:'
          end if
          if ( lperturb_ts ) then
            if ( myid == italk ) then
              write(stdout,'(a,f7.2,a)') &
                    'TS with value ',perturb_frac_ts*d_100,'%'
            end if
            call randify(rspace2,perturb_frac_ts,icount(1),icount(2))
          end if
        end if
        call grid_distribute(rspace2,ts,jce1,jce2,ice1,ice2)
        istart(1) = 1
        istart(2) = 1
        istart(3) = 1
        istart(4) = ibcrec
        icount(1) = jx
        icount(2) = iy
        icount(3) = kz
        icount(4) = 1
        istatus = nf90_get_var(ibcin,icbc_ivar(3),rspace3,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable u read error', 'ICBC FILE')
        if ( ensemble_run ) then
          if ( lperturb_u ) then
            if ( myid == italk ) then
              write(stdout,'(a,f7.2,a)') &
                'U with value ',perturb_frac_u*d_100,'%'
            end if
            call randify(rspace3,perturb_frac_u,icount(1),icount(2),icount(3))
          end if
        end if
        call grid_distribute(rspace3,u,jde1,jde2,ide1,ide2,1,kz)
        istatus = nf90_get_var(ibcin,icbc_ivar(4),rspace3,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable v read error', 'ICBC FILE')
        if ( ensemble_run ) then
          if ( lperturb_v ) then
            if ( myid == italk ) then
              write(stdout,'(a,f7.2,a)') &
                'V with value ',perturb_frac_v*d_100,'%'
            end if
            call randify(rspace3,perturb_frac_v,icount(1),icount(2),icount(3))
          end if
        end if
        call grid_distribute(rspace3,v,jde1,jde2,ide1,ide2,1,kz)
        istatus = nf90_get_var(ibcin,icbc_ivar(5),rspace3,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable t read error', 'ICBC FILE')
        if ( ensemble_run ) then
          if ( lperturb_t ) then
            if ( myid == italk ) then
              write(stdout,'(a,f7.2,a)') &
                'T with value ',perturb_frac_t*d_100,'%'
            end if
            call randify(rspace3,perturb_frac_t,icount(1),icount(2),icount(3))
          end if
        end if
        call grid_distribute(rspace3,t,jce1,jce2,ice1,ice2,1,kz)
        istatus = nf90_get_var(ibcin,icbc_ivar(6),rspace3,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable qx read error', 'ICBC FILE')
        if ( ensemble_run ) then
          if ( lperturb_q ) then
            if ( myid == italk ) then
              write(stdout,'(a,f7.2,a)') &
                'Q with value ',perturb_frac_q*d_100,'%'
            end if
            call randify(rspace3,perturb_frac_q,icount(1),icount(2),icount(3))
            where ( rspace3 < 1.0e-8_rkx ) rspace3 = 1.0e-8_rkx
          end if
        end if
        call grid_distribute(rspace3,qv,jce1,jce2,ice1,ice2,1,kz)
        if ( has_qc ) then
          istatus = nf90_get_var(ibcin,icbc_ivar(10),rspace3,istart,icount)
          call check_ok(__FILE__,__LINE__,'variable qc read error', 'ICBC FILE')
          call grid_distribute(rspace3,qc,jce1,jce2,ice1,ice2,1,kz)
        end if
        if ( has_qi ) then
          istatus = nf90_get_var(ibcin,icbc_ivar(11),rspace3,istart,icount)
          call check_ok(__FILE__,__LINE__,'variable qi read error', 'ICBC FILE')
          call grid_distribute(rspace3,qi,jce1,jce2,ice1,ice2,1,kz)
        end if
        if ( idynamic == 2 ) then
          istatus = nf90_get_var(ibcin,icbc_ivar(7),rspace3,istart,icount)
          call check_ok(__FILE__,__LINE__,'variable pp read error', 'ICBC FILE')
          call grid_distribute(rspace3,pp,jce1,jce2,ice1,ice2,1,kz)
          istatus = nf90_get_var(ibcin,icbc_ivar(8),rspace3,istart,icount)
          call check_ok(__FILE__,__LINE__,'variable w read error', 'ICBC FILE')
          call grid_distribute(rspace3,tempw,jce1,jce2,ice1,ice2,1,kz)
          ww(jce1:jce2,ice1:ice2,2:kzp1) = tempw(jce1:jce2,ice1:ice2,1:kz)
          istart(1) = 1
          istart(2) = 1
          istart(3) = ibcrec
          icount(1) = jx
          icount(2) = iy
          icount(3) = 1
          istatus = nf90_get_var(ibcin,icbc_ivar(9),rspace2, &
                                 istart(1:3),icount(1:3))
          call check_ok(__FILE__,__LINE__, &
                        'variable wtop read error', 'ICBC FILE')
          call grid_distribute(rspace2,tempwtop,jce1,jce2,ice1,ice2)
          ww(jce1:jce2,ice1:ice2,1) = tempwtop(jce1:jce2,ice1:ice2)
        else if ( idynamic == 3 ) then
          istatus = nf90_get_var(ibcin,icbc_ivar(7),rspace3,istart,icount)
          call check_ok(__FILE__,__LINE__, &
                        'variable pai read error', 'ICBC FILE')
          call grid_distribute(rspace3,pai,jce1,jce2,ice1,ice2,1,kz)
        end if
      else
        call grid_distribute(rspace2,ps,jce1,jce2,ice1,ice2)
        call grid_distribute(rspace2,ts,jce1,jce2,ice1,ice2)
        call grid_distribute(rspace3,u,jde1,jde2,ide1,ide2,1,kz)
        call grid_distribute(rspace3,v,jde1,jde2,ide1,ide2,1,kz)
        call grid_distribute(rspace3,t,jce1,jce2,ice1,ice2,1,kz)
        call grid_distribute(rspace3,qv,jce1,jce2,ice1,ice2,1,kz)
        if ( has_qc ) then
          call grid_distribute(rspace3,qc,jce1,jce2,ice1,ice2,1,kz)
        end if
        if ( has_qi ) then
          call grid_distribute(rspace3,qi,jce1,jce2,ice1,ice2,1,kz)
        end if
        if ( idynamic == 2 ) then
          call grid_distribute(rspace3,pp,jce1,jce2,ice1,ice2,1,kz)
          call grid_distribute(rspace3,tempw,jce1,jce2,ice1,ice2,1,kz)
          ww(jce1:jce2,ice1:ice2,2:kzp1) = tempw(jce1:jce2,ice1:ice2,1:kz)
          call grid_distribute(rspace2,tempwtop,jce1,jce2,ice1,ice2)
          ww(jce1:jce2,ice1:ice2,1) = tempwtop(jce1:jce2,ice1:ice2)
        else if ( idynamic == 3 ) then
          call grid_distribute(rspace3,pai,jce1,jce2,ice1,ice2,1,kz)
        end if
      end if
    end if
    if ( do_parallel_netcdf_in .or. myid == iocpu ) then
      call outstream_netcdf_unlock()
    end if
    if ( itweak == 1 ) then
      if ( itweak_sst == 1 ) then
        do concurrent ( j = jce1:jce2, i = ice1:ice2 )
          if ( ilnd(j,i) == 0 ) then
            ts(j,i) = ts(j,i) + sst_tweak
          end if
        end do
      end if
      if ( itweak_temperature == 1 ) then
        do k = 1, kz
          do i = ice1, ice2
            do j = jce1, jce2
              told = t(j,i,k)
              ! The below is not correct for non-hydro !
              pold = sigma(k)*ps(j,i)*d_100
              satvp = pfesat(told,ps(j,i))
              rhold = max((qv(j,i,k)/(ep2*satvp/(pold-satvp))),d_zero)
              tnew = t(j,i,k) + temperature_tweak
              pnew = pold*(tnew/told)
              satvp = pfesat(tnew,ps(j,i))
              qv(j,i,k) = max(rhold*ep2*satvp/(pnew-satvp),d_zero)
              t(j,i,k) = tnew
            end do
          end do
        end do
      end if
    end if
    !@acc call nvtxEndRange
  end subroutine read_icbc

  subroutine read_som(qflx)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: qflx
    integer(ik4), dimension(4) :: istart, icount
    character(len=3), dimension(12), parameter :: cmon = &
      ['jan','feb','mar','apr','may','jun', &
       'jul','aug','sep','oct','nov','dec']

    if ( myid == italk ) then
      write(stdout,*) 'Reading SOM data for ',cmon(somrec)
    end if
    if ( do_parallel_netcdf_in .or. myid == iocpu ) then
      call outstream_netcdf_lock()
    end if
    if ( do_parallel_netcdf_in ) then
      istart(1) = jci1-1
      istart(2) = ici1-1
      istart(3) = somrec
      icount(1) = jci2-jci1
      icount(2) = ici2-ici1
      icount(3) = 1
      istatus = nf90_get_var(somin,som_ivar(1),rspacesom, &
              istart(1:3),icount(1:3))
      call check_ok(__FILE__,__LINE__,'variable qflx read error', 'SOM FILE')
      qflx(jci1:jci2,ici1:ici2) = rspacesom
      else
        if ( myid == iocpu ) then
          istart(1) = 1
          istart(2) = 1
          istart(3) = somrec
          icount(1) = jx-3
          icount(2) = iy-3
          icount(3) = 1
          istatus = nf90_get_var(somin,som_ivar(1), &
                  rspacesom,istart(1:3),icount(1:3))
          call check_ok(__FILE__,__LINE__, &
                  'variable qflx read error', 'SOM FILE')
          call grid_distribute(rspacesom,qflx,jci1,jci2,ici1,ici2)
        else
          call grid_distribute(rspacesom,qflx,jci1,jci2,ici1,ici2)
        end if
      end if
    if ( do_parallel_netcdf_in .or. myid == iocpu ) then
      call outstream_netcdf_unlock()
    end if
  end subroutine read_som

  subroutine close_icbc
#ifdef OPENACC
    use cudafor
#endif
    implicit none
#ifdef ASYNC_NETCDF
    call release_icbc_prefetch(.true.)
#endif
    if (ibcin >= 0) then
      call outstream_netcdf_lock()
      istatus = nf90_close(ibcin)
      call check_ok(__FILE__,__LINE__,'Error Close ICBC file','ICBC FILE')
      call outstream_netcdf_unlock()
      if ( allocated(icbc_idate) ) deallocate(icbc_idate)
      ibcin = -1
    end if
#ifdef OPENACC
    if ( associated(rspace2) .or. associated(rspace3) ) then
      istat = cudaFreeHost(cptr)
      if ( istat /= cudaSuccess ) then
        write(*,*) 'cudaFreeHost failed with status: ',istat
        stop
      end if
      rspace2 => null()
      rspace3 => null()
    end if
#else
    if ( associated(rspace2) ) deallocate(rspace2)
    if ( associated(rspace3) ) deallocate(rspace3)
#endif
  end subroutine close_icbc

  subroutine close_clmbc
    implicit none
    if (clmbcin >= 0) then
      call outstream_netcdf_lock()
      istatus = nf90_close(clmbcin)
      call check_ok(__FILE__,__LINE__,'Error Close SFBC file','SFBC FILE')
      call outstream_netcdf_unlock()
      if ( allocated(clmbc_idate) ) deallocate(clmbc_idate)
      clmbcin = -1
    end if
  end subroutine close_clmbc

  subroutine close_som
    implicit none
    if (somin >= 0) then
      call outstream_netcdf_lock()
      istatus = nf90_close(somin)
      call check_ok(__FILE__,__LINE__,'Error Close SOM file','SOM FILE')
      call outstream_netcdf_unlock()
      somin = -1
    end if
    if ( associated(rspacesom) ) deallocate(rspacesom)
  end subroutine close_som

  subroutine check_ok(f,l,m1,mf)
    implicit none
    character(*), intent(in) :: f, m1, mf
    integer(ik4), intent(in) :: l
    if (istatus /= nf90_noerr) then
      write (stderr,*) trim(m1)
      write (stderr,*) nf90_strerror(istatus)
      call fatal(f,l,trim(mf))
    end if
  end subroutine check_ok

  subroutine read_ccn(ccn)
    implicit none
    type(tccn_data), intent(inout) :: ccn
    character(len=256) :: ccnfile
    integer(ik4) :: ncid, idim, ivar, i, j, k, n
    integer(ik4), parameter :: nseas = 4
    character(len=*), parameter :: omcam = 'CCN_climatology_cloudfree_8km.nc'

    ccnfile = trim(inpglob)//pthsep//'CCN'//pthsep//omcam
    call outstream_netcdf_lock()
    call openfile_withname(ccnfile,ncid)
    istatus = nf90_inq_dimid(ncid, 'lon', idim)
    call check_ok(__FILE__,__LINE__,'Dimension lon miss', 'CCN FILE')
    istatus = nf90_inquire_dimension(ncid, idim, len=ccn%nlon)
    call check_ok(__FILE__,__LINE__,'Dimension lon read error', 'CCN FILE')
    istatus = nf90_inq_dimid(ncid, 'lat', idim)
    call check_ok(__FILE__,__LINE__,'Dimension lat miss', 'CCN FILE')
    istatus = nf90_inquire_dimension(ncid, idim, len=ccn%nlat)
    call check_ok(__FILE__,__LINE__,'Dimension lat read error', 'CCN FILE')
    istatus = nf90_inq_dimid(ncid, 'altitude', idim)
    call check_ok(__FILE__,__LINE__,'Dimension altitude miss', 'CCN FILE')
    istatus = nf90_inquire_dimension(ncid, idim, len=ccn%nlev)
    call check_ok(__FILE__,__LINE__,'Dimension altitude read error', 'CCN FILE')
    allocate(ccn%lon(ccn%nlon), ccn%lat(ccn%nlat), ccn%altitude(ccn%nlev))
    allocate(ccn%ccn(ccn%nlon,ccn%nlat,ccn%nlev,nseas))
    istatus = nf90_inq_varid(ncid, 'lon', ivar)
    call check_ok(__FILE__,__LINE__,'variable lon miss', 'CCN FILE')
    istatus = nf90_get_var(ncid, ivar, ccn%lon)
    call check_ok(__FILE__,__LINE__,'variable lon read error', 'CCN FILE')
    istatus = nf90_inq_varid(ncid, 'lat', ivar)
    call check_ok(__FILE__,__LINE__,'variable lat miss', 'CCN FILE')
    istatus = nf90_get_var(ncid, ivar, ccn%lat)
    call check_ok(__FILE__,__LINE__,'variable lat read error', 'CCN FILE')
    istatus = nf90_inq_varid(ncid, 'altitude', ivar)
    call check_ok(__FILE__,__LINE__,'variable altitude miss', 'CCN FILE')
    istatus = nf90_get_var(ncid, ivar, ccn%altitude)
    call check_ok(__FILE__,__LINE__,'variable altitude read error', 'CCN FILE')
    istatus = nf90_inq_varid(ncid, 'CCN_cl_sn', ivar)
    call check_ok(__FILE__,__LINE__,'variable CCN_cl_sn miss', 'CCN FILE')
    istatus = nf90_get_var(ncid, ivar, ccn%ccn)
    call check_ok(__FILE__,__LINE__,'variable CCN_cl_sn read error', 'CCN FILE')
    istatus = nf90_close(ncid)
    call check_ok(__FILE__,__LINE__,'Error Close CCN file','CCN FILE')
    call outstream_netcdf_unlock()
    ! Unit of measure: set altitude in meters, CCN in #/m3
    ccn%altitude(:) = ccn%altitude(:) * 1000.0_rkx  ! In file km amsl
    do concurrent ( j = 1:ccn%nlon, i = 1:ccn%nlat, &
                    k = 1:ccn%nlev, n = 1:nseas )
      if ( ccn%ccn(j,i,k,n) < 1.0E20_rkx ) then
        ccn%ccn(j,i,k,n) = ccn%ccn(j,i,k,n) * 1.0e6_rkx ! In file #/cm3
      else
        ccn%ccn(j,i,k,n) = -9999.0_rkx
      end if
    end do
  end subroutine read_ccn

end module mod_ncio
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
