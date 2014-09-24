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

subroutine myabort
  implicit none
  call abort
end subroutine myabort

module mod_remap
  use mod_realkinds
  use mod_intkinds

  implicit none

  private

  public :: remap

  interface remap
    module procedure remap_real4
    module procedure remap_real8
    module procedure remap_int4
  end interface remap

  contains

  subroutine remap_int4(nsg,iv,mask,var)
    implicit none
    integer(ik4) , intent(in) :: nsg
    integer(ik4) , intent(inout) , dimension(:) :: iv
    integer(ik4) , intent(in) , dimension(:,:) :: mask
    integer(ik4) , intent(out) , dimension(:,:) :: var
    integer(ik4) :: ib , n1 , n2 , i , j , ni , nj , ii , jj
    var(:,:) = bigint
    if ( nsg == 1 ) then
      n1 = size(mask,1)
      n2 = size(mask,2)
      ib = 1
      do i = 1 , n2
        do j = 1 , n1
          if ( mask(j,i) > 0 ) then
            var(j,i) = iv(ib)
            ib = ib + 1
          end if
        end do
      end do
    else
      n1 = size(mask,1)/nsg
      n2 = size(mask,2)/nsg
      ib = 1
      do i = 1 , n2
        do j = 1 , n1
          do ni = 1 , nsg
            ii = (i-1)*nsg + ni
            do nj = 1 , nsg
              jj = (j-1)*nsg + nj
              if ( mask(jj,ii) > 0 ) then
                var(jj,ii) = iv(ib)
                ib = ib + 1
              end if
            end do
          end do
        end do
      end do
    end if
  end subroutine remap_int4

  subroutine remap_real4(nsg,iv,mask,var)
    implicit none
    integer(ik4) , intent(in) :: nsg
    real(rk4) , intent(inout) , dimension(:) :: iv
    integer(ik4) , intent(in) , dimension(:,:) :: mask
    real(rk4) , dimension(:,:) , intent(out) :: var
    integer(ik4) :: ib , n1 , n2 , i , j , ni , nj , ii , jj
    var(:,:) = 1E+36
    if ( nsg == 1 ) then
      n1 = size(mask,1)
      n2 = size(mask,2)
      ib = 1
      do i = 1 , n2
        do j = 1 , n1
          if ( mask(j,i) > 0 ) then
            var(j,i) = iv(ib)
            ib = ib + 1
          end if
        end do
      end do
    else
      n1 = size(mask,1)/nsg
      n2 = size(mask,2)/nsg
      ib = 1
      do i = 1 , n2
        do j = 1 , n1
          do ni = 1 , nsg
            ii = (i-1)*nsg + ni
            do nj = 1 , nsg
              jj = (j-1)*nsg + nj
              if ( mask(jj,ii) > 0 ) then
                var(jj,ii) = iv(ib)
                ib = ib + 1
              end if
            end do
          end do
        end do
      end do
    end if
  end subroutine remap_real4

  subroutine remap_real8(nsg,iv,mask,var)
    implicit none
    integer(ik4) , intent(in) :: nsg
    real(rk8) , intent(inout) , dimension(:) :: iv
    integer(ik4) , intent(in) , dimension(:,:) :: mask
    real(rk8) , dimension(:,:) , intent(out) :: var
    integer(ik4) :: ib , n1 , n2 , i , j , ni , nj , ii , jj
    var(:,:) = 1D+36
    if ( nsg == 1 ) then
      n1 = size(mask,1)
      n2 = size(mask,2)
      ib = 1
      do i = 1 , n2
        do j = 1 , n1
          if ( mask(j,i) > 0 ) then
            var(j,i) = iv(ib)
            ib = ib + 1
          end if
        end do
      end do
    else
      n1 = size(mask,1)/nsg
      n2 = size(mask,2)/nsg
      ib = 1
      do i = 1 , n2
        do j = 1 , n1
          do ni = 1 , nsg
            ii = (i-1)*nsg + ni
            do nj = 1 , nsg
              jj = (j-1)*nsg + nj
              if ( mask(jj,ii) > 0 ) then
                var(jj,ii) = iv(ib)
                ib = ib + 1
              end if
            end do
          end do
        end do
      end do
    end if
  end subroutine remap_real8

end module mod_remap

program clm45_1dto2d
  use mod_realkinds
  use mod_intkinds
  use mod_stdio
  use mod_nchelper
  use mod_dynparam , only : iomode
  use mod_remap
  use netcdf

  implicit none

  character(len=256) :: prgname , ncfile , ncoutfile
  integer(ik4) :: numarg , istatus
  integer(ik4) :: ncid , ndims , nvars , natts , udimid
  integer(ik4) :: ncoutid
  integer(ik4) :: varid , idtime
  integer(ik4) , allocatable , dimension(:) :: dimids , dsize
  integer(ik4) , allocatable , dimension(:) :: varids
  integer(ik4) , allocatable , dimension(:) :: outdimids
  integer(ik4) , allocatable , dimension(:) :: mapids
  integer(ik4) , allocatable , dimension(:) :: vtype , vndims
  integer(ik4) , allocatable , dimension(:,:) :: vshape
  integer(ik4) :: iy , jx
  integer(ik4) :: id , iv , ia , iid1 , iid2
  integer(ik4) :: n , m
  integer(ik4) :: lndgrid , iydim , jxdim
  integer(ik4) , allocatable , dimension(:,:) :: mask
  real(rk4) , allocatable , dimension(:, :) :: var2d_single
  real(rk8) , allocatable , dimension(:, :) :: var2d_double
  integer(ik4) , allocatable , dimension(:, :) :: var2d_int
  logical , allocatable , dimension(:) :: lexpand
  character(len=32) :: vname , dname
  character(len=64) :: aname
  real(rk8) , allocatable , dimension(:) :: double_var_1d
  real(rk4) , allocatable , dimension(:) :: single_var_1d
  integer(ik4) , allocatable , dimension(:) :: int_var_1d
  real(rk8) , allocatable , dimension(:,:) :: double_var_2d
  real(rk4) , allocatable , dimension(:,:) :: single_var_2d
  integer(ik4) , allocatable , dimension(:,:) :: int_var_2d
  real(rk8) , allocatable , dimension(:,:,:) :: double_var_3d
  real(rk4) , allocatable , dimension(:,:,:) :: single_var_3d
  integer(ik4) , allocatable , dimension(:,:,:) :: int_var_3d
  integer(ik4) :: nsg

  call get_command_argument(0,value=prgname)
  numarg = command_argument_count()
  if (numarg < 1) then
    write (stderr,*) 'Not enough arguments.'
    write (stderr,*) ' '
    write (stderr,*) 'Usage : ', trim(prgname), ' Rcmfile.clm.h.nc'
    stop
  end if

  call get_command_argument(1,value=ncfile)

  iid1 = scan(ncfile, '/', .true.)+1
  iid2 = scan(ncfile, '.', .true.)-1
  ncoutfile = ncfile(iid1:iid2)//'_2d.nc'

  istatus = nf90_open(ncfile, nf90_nowrite, ncid)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error Open file '//trim(ncfile))

  istatus = nf90_create(ncoutfile, iomode, ncoutid)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error Create file '//trim(ncoutfile))

  istatus = nf90_get_att(ncid,nf90_global,'regcm_subgrid',nsg)
  if ( istatus /= nf90_noerr ) then
    write(stdout,*) 'Assuming nsg = 1'
    nsg = 1
  end if

  istatus = nf90_inquire(ncid,ndims,nvars,natts,udimid)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error inquire file '//trim(ncfile))
  allocate(dimids(ndims))
  allocate(dsize(ndims))
  allocate(outdimids(ndims))
  allocate(mapids(ndims))
  allocate(lexpand(nvars))
  allocate(varids(nvars))
  allocate(vndims(nvars))
  allocate(vtype(nvars))
  allocate(vshape(nvars,ndims))

  do ia = 1 , natts
    istatus = nf90_inq_attname(ncid, nf90_global, ia, aname)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire attribute')
    istatus = nf90_copy_att(ncid, nf90_global, aname, ncoutid, nf90_global)
    call checkncerr(istatus,__FILE__,__LINE__,'Error copy attribute')
  end do

  istatus = nf90_inq_dimid(ncid, 'jx', jxdim)
  call checkncerr(istatus,__FILE__,__LINE__,'Dimension x missing')
  istatus = nf90_inquire_dimension(ncid, jxdim, len=jx)
  call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dimension x')
  istatus = nf90_inq_dimid(ncid, 'iy', iydim)
  call checkncerr(istatus,__FILE__,__LINE__,'Dimension y missing')
  istatus = nf90_inquire_dimension(ncid, iydim, len=iy)
  call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dimension y')

  allocate(mask(jx,iy))
  allocate(var2d_single(jx,iy))
  allocate(var2d_double(jx,iy))
  allocate(var2d_int(jx,iy))
  istatus = nf90_inq_varid(ncid, "regcm_mask", varid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error find variable regcm_mask')
  istatus = nf90_get_var(ncid, varid, mask)
  call checkncerr(istatus,__FILE__,__LINE__,'Error read variable regcm_mask')

  do id = 1 , ndims
    istatus = nf90_inquire_dimension(ncid, id, dname, dsize(id))
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dimension')

    if ( dname == 'gridcell' .or. dname == 'landunit' .or. &
         dname == 'column' .or. dname == 'pft' .or. &
         dname == 'string_length') cycle

    if ( dname /= 'lndgrid' ) then
      if ( id == udimid ) then
        istatus = nf90_def_dim(ncoutid, dname, nf90_unlimited, outdimids(id))
      else
        istatus = nf90_def_dim(ncoutid, dname, dsize(id), outdimids(id))
      end if
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error define dimension '//trim(dname))
    else
      lndgrid = id
    end if
  end do

  vshape(:,:) = -1
  do iv = 1 , nvars
    lexpand(iv) = .false.
    istatus = nf90_inquire_variable(ncid,iv,name=vname,xtype=vtype(iv), &
                                    ndims=vndims(iv),dimids=dimids,natts=natts)
    ndims = vndims(iv)
    vshape(iv,1:ndims) = dsize(dimids(1:vndims(iv)))
    if ( vname == 'time' ) then
      idtime = iv
    end if
    do id = 1 , vndims(iv)
      if ( dimids(id) == lndgrid ) then
        lexpand(iv) = .true.
        mapids(id) = outdimids(jxdim)
        mapids(id+1) = outdimids(iydim)
        if ( vndims(iv) > id ) then
          mapids(id+2:vndims(iv)+1) = outdimids(dimids(id+1:ndims))
        end if
        ndims = ndims + 1
        exit
      end if
      mapids(id) = outdimids(dimids(id))
    end do

    istatus = nf90_def_var(ncoutid, vname, vtype(iv), &
                           mapids(1:ndims), varids(iv))
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error define variable '//trim(vname))
    do ia = 1 , natts
      istatus = nf90_inq_attname(ncid, iv, ia, aname)
      call checkncerr(istatus,__FILE__,__LINE__,'Error inquire attribute')
      istatus = nf90_copy_att(ncid, iv, aname, ncoutid, varids(iv))
      call checkncerr(istatus,__FILE__,__LINE__,'Error copy attribute')
    end do
    if ( vname == 'topo' .or. vname == 'lon' .or. &
         vname == 'lat' .or. vname == 'landmask' .or. &
         vname == 'landfrac' .or. vname == 'pftmask' .or. &
         vname == 'area' ) then
      istatus = nf90_put_att(ncoutid, varids(iv), '_FillValue', 1.D+36)
      call checkncerr(istatus,__FILE__,__LINE__,'Error set attribute')
    end if
  end do

  istatus = nf90_enddef(ncoutid)

  ! Put times
  allocate(double_var_1d(vshape(idtime,1)))
  istatus = nf90_get_var(ncid,idtime,double_var_1d)
  call checkncerr(istatus,__FILE__,__LINE__,'Error read time')
  istatus = nf90_put_var(ncoutid,idtime,double_var_1d)
  call checkncerr(istatus,__FILE__,__LINE__,'Error write time')
  deallocate(double_var_1d)

  writeloop: &
  do iv = 1 , nvars
    if ( iv == idtime ) cycle writeloop
    if ( .not. lexpand(iv) ) then
      select case (vtype(iv))
        case (nf90_real)
          select case (vndims(iv))
            case(1)
              allocate(single_var_1d(vshape(iv,1)))
              istatus = nf90_get_var(ncid,iv,single_var_1d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              istatus = nf90_put_var(ncoutid,iv,single_var_1d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
              deallocate(single_var_1d)
            case(2)
              allocate(single_var_2d(vshape(iv,1),vshape(iv,2)))
              istatus = nf90_get_var(ncid,iv,single_var_2d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              istatus = nf90_put_var(ncoutid,iv,single_var_2d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
              deallocate(single_var_2d)
            case(3)
              allocate(single_var_3d(vshape(iv,1),vshape(iv,2),vshape(iv,3)))
              istatus = nf90_get_var(ncid,iv,single_var_3d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              istatus = nf90_put_var(ncoutid,iv,single_var_3d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
              deallocate(single_var_3d)
          end select
        case (nf90_double)
          select case (vndims(iv))
            case(1)
              allocate(double_var_1d(vshape(iv,1)))
              istatus = nf90_get_var(ncid,iv,double_var_1d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              istatus = nf90_put_var(ncoutid,iv,double_var_1d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
              deallocate(double_var_1d)
            case(2)
              allocate(double_var_2d(vshape(iv,1),vshape(iv,2)))
              istatus = nf90_get_var(ncid,iv,double_var_2d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              istatus = nf90_put_var(ncoutid,iv,double_var_2d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
              deallocate(double_var_2d)
            case(3)
              allocate(double_var_3d(vshape(iv,1),vshape(iv,2),vshape(iv,3)))
              istatus = nf90_get_var(ncid,iv,double_var_3d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              istatus = nf90_put_var(ncoutid,iv,double_var_3d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
              deallocate(double_var_3d)
          end select
        case (nf90_int)
          select case (vndims(iv))
            case(1)
              allocate(int_var_1d(vshape(iv,1)))
              istatus = nf90_get_var(ncid,iv,int_var_1d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              istatus = nf90_put_var(ncoutid,iv,int_var_1d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
              deallocate(int_var_1d)
            case(2)
              allocate(int_var_2d(vshape(iv,1),vshape(iv,2)))
              istatus = nf90_get_var(ncid,iv,int_var_2d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              istatus = nf90_put_var(ncoutid,iv,int_var_2d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
              deallocate(int_var_2d)
            case(3)
              allocate(int_var_3d(vshape(iv,1),vshape(iv,2),vshape(iv,3)))
              istatus = nf90_get_var(ncid,iv,int_var_3d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              istatus = nf90_put_var(ncoutid,iv,int_var_3d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
              deallocate(int_var_3d)
          end select
        case default
          cycle writeloop
      end select
    else
      select case (vtype(iv))
        case (nf90_real)
          select case (vndims(iv))
            case(1)
              allocate(single_var_1d(vshape(iv,1)))
              istatus = nf90_get_var(ncid,iv,single_var_1d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              call remap(nsg,single_var_1d,mask,var2d_single)
              istatus = nf90_put_var(ncoutid,iv,var2d_single)
              call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
              deallocate(single_var_1d)
            case(2)
              allocate(single_var_2d(vshape(iv,1),vshape(iv,2)))
              istatus = nf90_get_var(ncid,iv,single_var_2d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              do n = 1 , vshape(iv,2)
                call remap(nsg,single_var_2d(:,n),mask,var2d_single)
                istatus = nf90_put_var(ncoutid,iv,var2d_single, &
                        start=(/1,1,n/), count=(/jx,iy,1/))
                call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
              end do
              deallocate(single_var_2d)
            case(3)
              allocate(single_var_3d(vshape(iv,1),vshape(iv,2),vshape(iv,3)))
              istatus = nf90_get_var(ncid,iv,single_var_3d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              do n = 1 , vshape(iv,3)
                do m = 1 , vshape(iv,2)
                  call remap(nsg,single_var_3d(:,m,n),mask,var2d_single)
                  istatus = nf90_put_var(ncoutid,iv,var2d_single, &
                          start=(/1,1,m,n/), count=(/jx,iy,1,1/))
                  call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
                end do
              end do
              deallocate(single_var_3d)
          end select
        case (nf90_double)
          select case (vndims(iv))
            case(1)
              allocate(double_var_1d(vshape(iv,1)))
              istatus = nf90_get_var(ncid,iv,double_var_1d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              call remap(nsg,double_var_1d,mask,var2d_double)
              istatus = nf90_put_var(ncoutid,iv,var2d_double)
              call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
              deallocate(double_var_1d)
            case(2)
              allocate(double_var_2d(vshape(iv,1),vshape(iv,2)))
              istatus = nf90_get_var(ncid,iv,double_var_2d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              do n = 1 , vshape(iv,2)
                call remap(nsg,double_var_2d(:,n),mask,var2d_double)
                istatus = nf90_put_var(ncoutid,iv,var2d_double, &
                        start=(/1,1,n/), count=(/jx,iy,1/))
                call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
              end do
              deallocate(double_var_2d)
            case(3)
              allocate(double_var_3d(vshape(iv,1),vshape(iv,2),vshape(iv,3)))
              istatus = nf90_get_var(ncid,iv,double_var_3d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              do n = 1 , vshape(iv,3)
                do m = 1 , vshape(iv,2)
                  call remap(nsg,double_var_3d(:,m,n),mask,var2d_double)
                  istatus = nf90_put_var(ncoutid,iv,var2d_double, &
                          start=(/1,1,m,n/), count=(/jx,iy,1,1/))
                  call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
                end do
              end do
              deallocate(double_var_3d)
          end select
        case (nf90_int)
          select case (vndims(iv))
            case(1)
              allocate(int_var_1d(vshape(iv,1)))
              istatus = nf90_get_var(ncid,iv,int_var_1d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              call remap(nsg,int_var_1d,mask,var2d_int)
              istatus = nf90_put_var(ncoutid,iv,var2d_int)
              call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
              deallocate(int_var_1d)
            case(2)
              allocate(int_var_2d(vshape(iv,1),vshape(iv,2)))
              istatus = nf90_get_var(ncid,iv,int_var_2d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              do n = 1 , vshape(iv,2)
                call remap(nsg,int_var_2d(:,n),mask,var2d_int)
                istatus = nf90_put_var(ncoutid,iv,var2d_int, &
                        start=(/1,1,n/), count=(/jx,iy,1/))
                call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
              end do
              deallocate(int_var_2d)
            case(3)
              allocate(int_var_3d(vshape(iv,1),vshape(iv,2),vshape(iv,3)))
              istatus = nf90_get_var(ncid,iv,int_var_3d)
              call checkncerr(istatus,__FILE__,__LINE__,'Error read var')
              do n = 1 , vshape(iv,3)
                do m = 1 , vshape(iv,2)
                  call remap(nsg,int_var_3d(:,m,n),mask,var2d_int)
                  istatus = nf90_put_var(ncoutid,iv,var2d_int, &
                          start=(/1,1,m,n/), count=(/jx,iy,1,1/))
                  call checkncerr(istatus,__FILE__,__LINE__,'Error write var')
                end do
              end do
              deallocate(int_var_3d)
          end select
        case default
          cycle writeloop
      end select
    end if
  end do writeloop

  istatus = nf90_close(ncid)
  istatus = nf90_close(ncoutid)

end program clm45_1dto2d

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
