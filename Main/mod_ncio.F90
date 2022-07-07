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
module mod_ncio
!
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
  use mod_domain

  implicit none

  private

  public :: read_domain_info , read_subdomain_info
  public :: open_icbc , icbc_search , read_icbc , close_icbc
  public :: open_som , som_search , read_som , close_som
  public :: open_clmbc , clmbc_search , close_clmbc , read_clmbc
  public :: fixqcqi , we_have_qc , we_have_qi

  integer(ik4) :: ibcin , somin , clmbcin
  integer(ik4) :: istatus
  integer(ik4) :: ibcrec , ibcnrec , clmbcrec , clmbcnrec
  integer(ik4) :: somrec
  type(rcm_time_and_date) , dimension(:) , allocatable :: icbc_idate
  type(rcm_time_and_date) , dimension(:) , allocatable :: clmbc_idate
  integer(ik4) , dimension(16) :: icbc_ivar
  integer(ik4) , dimension(4) :: clmbc_ivar
  integer(ik4) , dimension(1) :: som_ivar
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

  real(rkx) , dimension(:,:) , pointer :: rspacesom => null()
  real(rkx) , dimension(:,:) , pointer :: rspace2 => null()
  real(rkx) , dimension(:,:,:) , pointer :: rspace3 => null()

  real(rkx) , dimension(:,:,:) , pointer :: tempw => null()
  real(rkx) , dimension(:,:) , pointer :: tempwtop => null()

  contains

  subroutine fixqcqi( )
    implicit none
    call bcast(has_qc)
    call bcast(has_qi)
  end subroutine

  logical function we_have_qc( )
    implicit none
    we_have_qc = has_qc
  end function we_have_qc
  logical function we_have_qi( )
    implicit none
    we_have_qi = has_qi
  end function we_have_qi

  subroutine read_domain_info(ht,lnd,tex,mask,area,xlat,xlon,dlat,dlon, &
                              ulat,ulon,vlat,vlon,msfx,msfd,msfu,  &
                              msfv,coriol,snowam,smoist,rmoist,hlake,ts0)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: ht
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: lnd
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: tex
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: mask
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: area
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: xlat
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: xlon
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: dlat
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: dlon
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: ulat
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: ulon
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: vlat
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: vlon
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: msfx
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: msfd
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: msfu
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: msfv
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: coriol
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: snowam
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: smoist
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: rmoist
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: hlake
    real(rkx) , intent(out) :: ts0
    real(rkx) , dimension(:,:) , pointer :: tempmoist
    character(len=256) :: dname
    character(len=8) :: csmoist
    integer(ik4) :: idmin , ilev
    integer(ik4) , dimension(2) :: istart , icount
    integer(ik4) , dimension(3) :: istart3 , icount3
    real(rkx) , dimension(:,:) , pointer :: rspace
    logical :: has_snow = .true.
    logical :: has_dhlake = .true.
    logical :: lerror

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
      lerror = .true.
      call read_var1d_static(idmin,'kz',sigma,lerror)
      if ( .not. lerror ) then
        call read_var1d_static(idmin,'sigma',sigma)
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
        do ilev = 1 , num_soil_layers
          istart3(3) = ilev
          call read_var2d_static(idmin,'rmoist',rspace, &
                                 istart=istart3,icount=icount3)
          rmoist(jde1:jde2,ide1:ide2,ilev) = rspace
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
        lerror = .true.
        call read_var1d_static(idmin,'kz',sigma,lerror)
        if ( .not. lerror ) then
          call read_var1d_static(idmin,'sigma',sigma)
        end if
        call bcast(sigma)
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
          write(stdout,*) 'Applying perturbation to input dataset:'
          if ( lperturb_topo ) then
            write(stdout,'(a,f7.2,a)') 'Topo with value ', &
              perturb_frac_topo*d_100,'%'
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
          do ilev = 1 , num_soil_layers
            istart3(3) = ilev
            call read_var2d_static(idmin,'rmoist',rspace, &
                    istart=istart3,icount=icount3)
            call grid_distribute(rspace,tempmoist,jde1,jde2,ide1,ide2)
            rmoist(jde1:jde2,ide1:ide2,ilev) = tempmoist
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
        call bcast(sigma)
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
          do ilev = 1 , num_soil_layers
            call grid_distribute(rspace,tempmoist,jde1,jde2,ide1,ide2)
            rmoist(jde1:jde2,ide1:ide2,ilev) = tempmoist
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
        call getmem3d(tempw,jce1,jce2,ice1,ice2,1,kz,'read_domain:tempw')
        call getmem2d(tempwtop,jce1,jce2,ice1,ice2,'read_domain:tempwtop')
      end if
    end if
  end subroutine read_domain_info

  subroutine read_subdomain_info(ht,lnd,tex,mask,area,xlat,xlon,hlake)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: ht
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: lnd
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: tex
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: mask
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: area
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: xlat
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: xlon
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: hlake
    character(len=256) :: dname
    integer(ik4) :: idmin
    integer(ik4) , dimension(2) :: istart , icount
    real(rkx) , dimension(:,:) , pointer :: rspace
    real(rkx) , dimension(:,:,:) , pointer :: rspace0
    logical :: has_dhlake = .true.
    character(len=3) :: sbstring

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
          write(stdout,*) 'Applying perturbation to input dataset:'
          if ( lperturb_topo ) then
            write(stdout,'(a,f7.2,a)') 'Topo with value ', &
              perturb_frac_topo*d_100,'%'
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

  integer function clmbc_search(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    type(rcm_time_interval) :: tdif
    character(len=32) :: appdat1 , appdat2
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

  integer function icbc_search(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    type(rcm_time_interval) :: tdif
    character(len=32) :: appdat1 , appdat2
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

  integer function som_search(imon)
    implicit none
    integer(ik4) , intent(in) :: imon
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
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=11) :: ctime
    integer(ik4) :: idimid , itvar , i , chkdiff
    real(rkx) , dimension(:) , allocatable :: icbc_nctime
    character(len=64) :: icbc_timeunits , icbc_timecal
    character(len=256) :: icbcname
    if ( .not. do_parallel_netcdf_in .and. myid /= iocpu ) then
      return
    end if
    call close_icbc
    write (ctime, '(a)') tochar10(idate)
    icbcname = trim(dirglob)//pthsep//trim(domname)// &
               '_ICBC.'//trim(ctime)//'.nc'
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
    do i = 1 , ibcnrec
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
    else
      istatus = nf90_inq_varid(ibcin, 'pp', icbc_ivar(7))
      if ( istatus == nf90_noerr ) then
        write(stderr,*) 'ERROR: Hydrostatic core and non-hydrostatic ICBC !'
        call fatal(__FILE__,__LINE__,'ICBC READ')
      end if
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
    if ( do_parallel_netcdf_in ) then
      allocate(rspace2(jde1:jde2,ide1:ide2))
      allocate(rspace3(jde1:jde2,ide1:ide2,kz))
    else
      allocate(rspace2(jx,iy))
      allocate(rspace3(jx,iy,kz))
    end if
  end subroutine open_icbc

  subroutine open_clmbc(idate)
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=11) :: ctime
    integer(ik4) :: idimid , itvar , i
    real(rkx) , dimension(:) , allocatable :: clmbc_nctime
    character(len=64) :: clmbc_timeunits , clmbc_timecal
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
    do i = 1 , clmbcnrec
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
    call openfile_withname(somname,somin)
    call check_domain(somin,.true.,.true.)
    istatus = nf90_inq_varid(somin, 'qflx', som_ivar(1))
    call check_ok(__FILE__,__LINE__,'variable qflx miss', 'SOM FILE')
    if ( do_parallel_netcdf_in ) then
      allocate(rspacesom(jci1:jci2,ici1:ici2))
    else
      allocate(rspacesom(2:jx-2,2:iy-2))
    end if
  end subroutine open_som

  subroutine read_clmbc(pr,ssr,strd,clt)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: pr
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: ssr
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: strd
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: clt

    integer(ik4) , dimension(3) :: istart , icount

    if ( clmbcrec > clmbcnrec ) then
      call open_clmbc(rcmtimer%idate)
    end if
    if ( myid == italk ) then
      write(stdout,*) 'Reading SF values for ',tochar(clmbc_idate(clmbcrec))
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
    clmbcrec = clmbcrec + 1
  end subroutine read_clmbc

  subroutine read_icbc(ps,ts,ilnd,u,v,t,qv,qc,qi,pp,ww)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: ps
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: ts
    integer(ik4) , pointer , dimension(:,:) , intent(in) :: ilnd
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: v
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: t
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: qv
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: qc
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: qi
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: pp
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: ww

    integer(ik4) , dimension(4) :: istart , icount
    integer(ik4) :: i , j , k
    real(rkx) :: told , pold , rhold , satvp , tnew , pnew

    if ( do_parallel_netcdf_in ) then
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
      ps(jce1:jce2,ice1:ice2) = rspace2(jce1:jce2,ice1:ice2)
      istatus = nf90_get_var(ibcin,icbc_ivar(2),rspace2,istart(1:3),icount(1:3))
      call check_ok(__FILE__,__LINE__,'variable ts read error', 'ICBC FILE')
      if ( ensemble_run ) then
        write(stdout,*) 'Applying perturbation to input dataset:'
        if ( lperturb_ts ) then
          if ( myid == italk ) then
            write(stdout,'(a,f7.2,a)') &
                    'TS with value ',perturb_frac_ts*d_100,'%'
          end if
          call randify(rspace2,perturb_frac_ts,icount(1),icount(2))
        end if
      end if
      ts(jce1:jce2,ice1:ice2) = rspace2(jce1:jce2,ice1:ice2)
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
      u(jde1:jde2,ide1:ide2,1:kz) = rspace3
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
      v(jde1:jde2,ide1:ide2,1:kz) = rspace3
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
      t(jce1:jce2,ice1:ice2,1:kz) = rspace3(jce1:jce2,ice1:ice2,1:kz)
      istatus = nf90_get_var(ibcin,icbc_ivar(6),rspace3,istart,icount)
      call check_ok(__FILE__,__LINE__,'variable qx read error', 'ICBC FILE')
      if ( ensemble_run ) then
        if ( lperturb_q ) then
          if ( myid == italk ) then
            write(stdout,'(a,f7.2,a)') 'Q with value ',perturb_frac_q*d_100,'%'
          end if
          call randify(rspace3,perturb_frac_q,icount(1),icount(2),icount(3))
        end if
      end if
      qv(jce1:jce2,ice1:ice2,1:kz) = rspace3(jce1:jce2,ice1:ice2,1:kz)
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
        pp(jce1:jce2,ice1:ice2,1:kz) = rspace3(jce1:jce2,ice1:ice2,1:kz)
        istatus = nf90_get_var(ibcin,icbc_ivar(8),rspace3,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable w read error', 'ICBC FILE')
        ww(jce1:jce2,ice1:ice2,2:kzp1) = rspace3(jce1:jce2,ice1:ice2,1:kz)
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
        ww(jce1:jce2,ice1:ice2,1) = rspace2(jce1:jce2,ice1:ice2)
      end if
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
            write(stdout,'(a,f7.2,a)') &
                    'PS with value ',perturb_frac_ps*d_100,'%'
            call randify(rspace2,perturb_frac_ps,icount(1),icount(2))
          end if
        end if
        call grid_distribute(rspace2,ps,jce1,jce2,ice1,ice2)
        istatus = nf90_get_var(ibcin,icbc_ivar(2), &
                rspace2,istart(1:3),icount(1:3))
        call check_ok(__FILE__,__LINE__,'variable ts read error', 'ICBC FILE')
        if ( ensemble_run ) then
          write(stdout,*) 'Applying perturbation to input dataset:'
          if ( lperturb_ts ) then
            write(stdout,'(a,f7.2,a)') &
                    'TS with value ',perturb_frac_ts*d_100,'%'
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
            write(stdout,'(a,f7.2,a)') 'U with value ',perturb_frac_u*d_100,'%'
            call randify(rspace3,perturb_frac_u,icount(1),icount(2),icount(3))
          end if
        end if
        call grid_distribute(rspace3,u,jde1,jde2,ide1,ide2,1,kz)
        istatus = nf90_get_var(ibcin,icbc_ivar(4),rspace3,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable v read error', 'ICBC FILE')
        if ( ensemble_run ) then
          if ( lperturb_v ) then
            write(stdout,'(a,f7.2,a)') 'V with value ',perturb_frac_v*d_100,'%'
            call randify(rspace3,perturb_frac_v,icount(1),icount(2),icount(3))
          end if
        end if
        call grid_distribute(rspace3,v,jde1,jde2,ide1,ide2,1,kz)
        istatus = nf90_get_var(ibcin,icbc_ivar(5),rspace3,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable t read error', 'ICBC FILE')
        if ( ensemble_run ) then
          if ( lperturb_t ) then
            write(stdout,'(a,f7.2,a)') 'T with value ',perturb_frac_t*d_100,'%'
            call randify(rspace3,perturb_frac_t,icount(1),icount(2),icount(3))
          end if
        end if
        call grid_distribute(rspace3,t,jce1,jce2,ice1,ice2,1,kz)
        istatus = nf90_get_var(ibcin,icbc_ivar(6),rspace3,istart,icount)
        call check_ok(__FILE__,__LINE__,'variable qx read error', 'ICBC FILE')
        if ( ensemble_run ) then
          if ( lperturb_q ) then
            write(stdout,'(a,f7.2,a)') 'Q with value ',perturb_frac_q*d_100,'%'
            call randify(rspace3,perturb_frac_q,icount(1),icount(2),icount(3))
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
        end if
      end if
    end if
    if ( itweak == 1 ) then
      if ( itweak_sst == 1 ) then
        do i = ice1 , ice2
          do j = jce1 , jce2
            if ( ilnd(j,i) == 0 ) then
              ts(j,i) = ts(j,i) + sst_tweak
            end if
          end do
        end do
      end if
      if ( itweak_temperature == 1 ) then
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              told = t(j,i,k)
              ! The below is not correct for non-hydro !
              pold = sigma(k)*ps(j,i)*d_100
              satvp = pfesat(told)
              rhold = max((qv(j,i,k)/(ep2*satvp/(pold-satvp))),d_zero)
              tnew = t(j,i,k) + temperature_tweak
              pnew = pold*(tnew/told)
              satvp = pfesat(tnew)
              qv(j,i,k) = max(rhold*ep2*satvp/(pnew-satvp),d_zero)
              t(j,i,k) = tnew
            end do
          end do
        end do
      end if
    end if

    contains

#include <pfesat.inc>

  end subroutine read_icbc

  subroutine read_som(qflx)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: qflx
    integer(ik4) , dimension(4) :: istart , icount
    character(len=3) , dimension(12) :: cmon = &
      ['jan','feb','mar','apr','may','jun', &
        'jul','aug','sep','oct','nov','dec']

    if ( myid == italk ) then
      write(stdout,*) 'Reading SOM data for ',cmon(somrec)
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
  end subroutine read_som

  subroutine close_icbc
    implicit none
    if (ibcin >= 0) then
      istatus = nf90_close(ibcin)
      call check_ok(__FILE__,__LINE__,'Error Close ICBC file','ICBC FILE')
      if ( allocated(icbc_idate) ) deallocate(icbc_idate)
      ibcin = -1
    end if
    if ( associated(rspace2) ) deallocate(rspace2)
    if ( associated(rspace3) ) deallocate(rspace3)
  end subroutine close_icbc

  subroutine close_clmbc
    implicit none
    if (clmbcin >= 0) then
      istatus = nf90_close(clmbcin)
      call check_ok(__FILE__,__LINE__,'Error Close SFBC file','SFBC FILE')
      if ( allocated(clmbc_idate) ) deallocate(clmbc_idate)
      clmbcin = -1
    end if
  end subroutine close_clmbc

  subroutine close_som
    implicit none
    if (somin >= 0) then
      istatus = nf90_close(somin)
      call check_ok(__FILE__,__LINE__,'Error Close SOM file','SOM FILE')
      somin = -1
    end if
    if ( associated(rspacesom) ) deallocate(rspacesom)
  end subroutine close_som

  subroutine check_ok(f,l,m1,mf)
    implicit none
    character(*) , intent(in) :: f, m1 , mf
    integer(ik4) , intent(in) :: l
    if (istatus /= nf90_noerr) then
      write (stderr,*) trim(m1)
      write (stderr,*) nf90_strerror(istatus)
      call fatal(f,l,trim(mf))
    end if
  end subroutine check_ok

end module mod_ncio
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
