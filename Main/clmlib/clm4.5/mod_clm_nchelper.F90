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

module mod_clm_nchelper

  use netcdf
  use mod_realkinds
  use mod_intkinds
  use mod_stdio
  use mod_memutil
  use mod_message
  use mod_mppparam
  use mod_dynparam
  use mod_regcm_types
  use mod_clm_decomp
  use mpi

  implicit none

  private

  integer(ik4) , parameter :: clm_maxdims = 64
  integer(ik4) , parameter :: clm_maxvars = 512

  type clm_filetype
    integer(ik4) :: ncid = -1
    character(len=256) :: fname
    integer(ik4) :: idimlast = 1
    integer(ik4) :: ivarlast = 1
    integer(ik4) , dimension(clm_maxdims) :: dimids
    integer(ik4) , dimension(clm_maxdims) :: dimhash
    character(len=32) , dimension(clm_maxdims) :: dimname
    integer(ik4) , dimension(clm_maxvars) :: varids
    integer(ik4) , dimension(clm_maxdims) :: varhash
    character(len=32) , dimension(clm_maxvars) :: varname
    integer(ik4) , dimension(:,:) , pointer :: i4buf => null()
    real(rk4) , dimension(:,:) , pointer :: r4buf => null()
    real(rk8) , dimension(:,:) , pointer :: r8buf => null()
  end type clm_filetype

  character , public , parameter :: clmvar_text       = 'c'
  logical , public , parameter :: clmvar_logical      = .false.
  integer(ik4) , public , parameter :: clmvar_integer = 1
  real(rk4) , public , parameter :: clmvar_real       = 1.0
  real(rk8) , public , parameter :: clmvar_double     = 1.0D0
  integer(ik4) , public , parameter :: clmvar_unlim   = -1

  interface assignment(=)
    module procedure copy_filetype
  end interface assignment(=)

  public :: clm_filetype

  public :: clm_createfile
  public :: clm_openfile
  public :: clm_closefile
  public :: clm_enddef
  public :: clm_checkncerr

  public :: clm_inqdim
  public :: clm_check_dim
  public :: clm_check_dims
  public :: clm_check_var
  public :: clm_check_dimlen

  public :: clm_adddim
  public :: clm_addatt
  public :: clm_addvar

  public :: clm_readvar
  public :: clm_writevar
  public :: clm_writevar_par

  public test_clmhelper

  integer(ik4) :: incstat

  integer(ik4) , dimension(clm_maxdims) :: usedims
  integer(ik4) , dimension(4) :: istart
  integer(ik4) , dimension(4) :: icount

  interface clm_addvar
    module procedure clm_addvar_char
    module procedure clm_addvar_int
    module procedure clm_addvar_real4
    module procedure clm_addvar_real8
    module procedure clm_addvar_logical
  end interface clm_addvar

  interface clm_addatt
    module procedure clm_addatt_text
    module procedure clm_addatt_integer
    module procedure clm_addatt_single
    module procedure clm_addatt_double
  end interface clm_addatt

  interface clm_readvar
    module procedure clm_readvar_text_0d
    module procedure clm_readvar_text_1d
    module procedure clm_readvar_logical_0d
    module procedure clm_readvar_logical_1d
    module procedure clm_readvar_logical_2d
    module procedure clm_readvar_logical_3d
    module procedure clm_readvar_logical_4d
    module procedure clm_readvar_integer_0d
    module procedure clm_readvar_integer_1d
    module procedure clm_readvar_integer_2d
    module procedure clm_readvar_integer_3d
    module procedure clm_readvar_integer_4d
    module procedure clm_readvar_real4_0d
    module procedure clm_readvar_real4_1d
    module procedure clm_readvar_real4_2d
    module procedure clm_readvar_real4_3d
    module procedure clm_readvar_real4_4d
    module procedure clm_readvar_real8_0d
    module procedure clm_readvar_real8_1d
    module procedure clm_readvar_real8_2d
    module procedure clm_readvar_real8_3d
    module procedure clm_readvar_real8_4d
    module procedure clm_readrec_logical_0d
    module procedure clm_readrec_logical_1d
    module procedure clm_readrec_logical_2d
    module procedure clm_readrec_logical_3d
    module procedure clm_readrec_integer_0d
    module procedure clm_readrec_integer_1d
    module procedure clm_readrec_integer_2d
    module procedure clm_readrec_integer_3d
    module procedure clm_readrec_real4_0d
    module procedure clm_readrec_real4_1d
    module procedure clm_readrec_real4_2d
    module procedure clm_readrec_real4_3d
    module procedure clm_readrec_real8_0d
    module procedure clm_readrec_real8_1d
    module procedure clm_readrec_real8_2d
    module procedure clm_readrec_real8_3d
  end interface clm_readvar

  interface clm_writevar
    module procedure clm_writevar_text_0d
    module procedure clm_writevar_text_1d
    module procedure clm_writevar_logical_0d
    module procedure clm_writevar_logical_1d
    module procedure clm_writevar_logical_2d
    module procedure clm_writevar_logical_3d
    module procedure clm_writevar_logical_4d
    module procedure clm_writevar_integer_0d
    module procedure clm_writevar_integer_1d
    module procedure clm_writevar_integer_2d
    module procedure clm_writevar_integer_3d
    module procedure clm_writevar_integer_4d
    module procedure clm_writevar_real4_0d
    module procedure clm_writevar_real4_1d
    module procedure clm_writevar_real4_2d
    module procedure clm_writevar_real4_3d
    module procedure clm_writevar_real4_4d
    module procedure clm_writevar_real8_0d
    module procedure clm_writevar_real8_1d
    module procedure clm_writevar_real8_2d
    module procedure clm_writevar_real8_3d
    module procedure clm_writevar_real8_4d
    module procedure clm_writerec_logical_0d
    module procedure clm_writerec_logical_1d
    module procedure clm_writerec_logical_2d
    module procedure clm_writerec_logical_3d
    module procedure clm_writerec_integer_0d
    module procedure clm_writerec_integer_1d
    module procedure clm_writerec_integer_2d
    module procedure clm_writerec_integer_3d
    module procedure clm_writerec_real4_0d
    module procedure clm_writerec_real4_1d
    module procedure clm_writerec_real4_2d
    module procedure clm_writerec_real4_3d
    module procedure clm_writerec_real8_0d
    module procedure clm_writerec_real8_1d
    module procedure clm_writerec_real8_2d
    module procedure clm_writerec_real8_3d
  end interface clm_writevar

  interface clm_writevar_par
    module procedure clm_writevar_logical_0d_par_sg
    module procedure clm_writevar_logical_1d_par_sg
    module procedure clm_writevar_logical_2d_par_sg
    module procedure clm_writevar_logical_3d_par_sg
    module procedure clm_writevar_logical_4d_par_sg
    module procedure clm_writevar_integer_0d_par_sg
    module procedure clm_writevar_integer_1d_par_sg
    module procedure clm_writevar_integer_2d_par_sg
    module procedure clm_writevar_integer_3d_par_sg
    module procedure clm_writevar_integer_4d_par_sg
    module procedure clm_writevar_real4_0d_par_sg
    module procedure clm_writevar_real4_1d_par_sg
    module procedure clm_writevar_real4_2d_par_sg
    module procedure clm_writevar_real4_3d_par_sg
    module procedure clm_writevar_real4_4d_par_sg
    module procedure clm_writevar_real8_0d_par_sg
    module procedure clm_writevar_real8_1d_par_sg
    module procedure clm_writevar_real8_2d_par_sg
    module procedure clm_writevar_real8_3d_par_sg
    module procedure clm_writevar_real8_4d_par_sg
    module procedure clm_writerec_logical_0d_par_sg
    module procedure clm_writerec_logical_1d_par_sg
    module procedure clm_writerec_logical_2d_par_sg
    module procedure clm_writerec_logical_3d_par_sg
    module procedure clm_writerec_integer_0d_par_sg
    module procedure clm_writerec_integer_1d_par_sg
    module procedure clm_writerec_integer_2d_par_sg
    module procedure clm_writerec_integer_3d_par_sg
    module procedure clm_writerec_real4_0d_par_sg
    module procedure clm_writerec_real4_1d_par_sg
    module procedure clm_writerec_real4_2d_par_sg
    module procedure clm_writerec_real4_3d_par_sg
    module procedure clm_writerec_real8_0d_par_sg
    module procedure clm_writerec_real8_1d_par_sg
    module procedure clm_writerec_real8_2d_par_sg
    module procedure clm_writerec_real8_3d_par_sg
    module procedure clm_writevar_logical_2d_par_gg
    module procedure clm_writevar_logical_3d_par_gg
    module procedure clm_writevar_logical_4d_par_gg
    module procedure clm_writevar_integer_2d_par_gg
    module procedure clm_writevar_integer_3d_par_gg
    module procedure clm_writevar_integer_4d_par_gg
    module procedure clm_writevar_real4_2d_par_gg
    module procedure clm_writevar_real4_3d_par_gg
    module procedure clm_writevar_real4_4d_par_gg
    module procedure clm_writevar_real8_2d_par_gg
    module procedure clm_writevar_real8_3d_par_gg
    module procedure clm_writevar_real8_4d_par_gg
    module procedure clm_writerec_logical_2d_par_gg
    module procedure clm_writerec_logical_3d_par_gg
    module procedure clm_writerec_integer_2d_par_gg
    module procedure clm_writerec_integer_3d_par_gg
    module procedure clm_writerec_real4_2d_par_gg
    module procedure clm_writerec_real4_3d_par_gg
    module procedure clm_writerec_real8_2d_par_gg
    module procedure clm_writerec_real8_3d_par_gg
  end interface clm_writevar_par

  contains
!
  subroutine clm_createfile(fname,ncid)
    implicit none
    character(len=*) , intent(in) :: fname
    type(clm_filetype) , intent(out) :: ncid

    if ( myid /= iocpu ) return
#ifdef NETCDF4_HDF5
    incstat = nf90_create(fname, &
             ior(ior(nf90_clobber,nf90_hdf5),nf90_classic_model),ncid%ncid)
#else
    incstat = nf90_create(fname, nf90_clobber, ncid%ncid)
#endif
    call clm_checkncerr(__FILE__,__LINE__, &
                    'Error creating NetCDF output '//trim(fname))
    ncid%fname = fname
    call getmem2d(ncid%i4buf,jout1,jout2,iout1,iout2,'clm_createfile')
    call getmem2d(ncid%r4buf,jout1,jout2,iout1,iout2,'clm_createfile')
    call getmem2d(ncid%r8buf,jout1,jout2,iout1,iout2,'clm_createfile')
    ncid%dimhash = huge(1)
    ncid%varhash = huge(1)
  end subroutine clm_createfile

  subroutine clm_openfile(fname,ncid)
    implicit none
    character(len=*) , intent(in) :: fname
    type(clm_filetype) , intent(out) :: ncid
    incstat = nf90_open(fname, nf90_nowrite, ncid%ncid)
    call clm_checkncerr(__FILE__,__LINE__, &
                    'Error open NetCDF input '//trim(fname))
    ncid%fname = fname
  end subroutine clm_openfile

  subroutine clm_enddef(ncid)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    if ( myid /= iocpu ) return
    incstat =  nf90_enddef(ncid%ncid)
    call clm_checkncerr(__FILE__,__LINE__,'Error enddef NetCDF output')
  end subroutine clm_enddef

  subroutine clm_inqdim(ncid,dname,dlen,lerror,lexist)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: dname
    integer(ik4) , intent(out) , optional :: dlen
    logical , intent(in) , optional :: lerror
    logical , intent(out) , optional :: lexist
    integer(ik4) :: idimid
    incstat = nf90_inq_dimid(ncid%ncid, dname, idimid)
    if ( incstat /= nf90_noerr ) then
      if ( present(lexist) ) then
        lexist = .false.
        return
      end if
      if ( present(lerror) ) then
        if ( lerror ) then
          if ( present(dlen) ) then
            dlen = -1
          end if
          return
        end if
      end if
    else
      if ( present(lexist) ) then
        lexist = .true.
        return
      end if
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search dimension '//dname//' in file '//trim(ncid%fname))
    if ( present(dlen) ) then
      incstat = nf90_inquire_dimension(ncid%ncid, idimid, len=dlen)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error read dimension '//dname//' in file '//trim(ncid%fname))
    end if
  end subroutine clm_inqdim

  subroutine clm_addatt_text(ncid,aname,aval,ivar,cvar)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: aname
    character(len=*) , intent(in) :: aval
    integer(ik4) , intent(in) , optional :: ivar
    character(len=*) , intent(in) , optional :: cvar
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    if ( present(ivar) ) then
      incstat = nf90_put_att(ncid%ncid,ivar,aname,aval)
    else if ( present(cvar) ) then
      ivarid = searchvar(ncid,cvar)
      if ( ivarid > 0 ) then
        incstat = nf90_put_att(ncid%ncid,ivarid,aname,aval)
      else
        incstat = nf90_enotvar
      end if
    else
      incstat = nf90_put_att(ncid%ncid,nf90_global,aname,aval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
       'Error adding attribute '//aname//' in file '//trim(ncid%fname))
  end subroutine clm_addatt_text

  subroutine clm_addatt_integer(ncid,aname,aval,ivar,cvar)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: aname
    integer(ik4) , intent(in) :: aval
    integer(ik4) , intent(in) , optional :: ivar
    character(len=*) , intent(in) , optional :: cvar
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    if ( present(ivar) ) then
      incstat = nf90_put_att(ncid%ncid,ivar,aname,aval)
    else if ( present(cvar) ) then
      ivarid = searchvar(ncid,cvar)
      if ( ivarid > 0 ) then
        incstat = nf90_put_att(ncid%ncid,ivarid,aname,aval)
      else
        incstat = nf90_enotvar
      end if
    else
      incstat = nf90_put_att(ncid%ncid,nf90_global,aname,aval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
       'Error adding attribute '//aname//' in file '//trim(ncid%fname))
  end subroutine clm_addatt_integer

  subroutine clm_addatt_single(ncid,aname,aval,ivar,cvar)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: aname
    real(rk4) , intent(in) :: aval
    integer(rk4) , intent(in) , optional :: ivar
    character(len=*) , intent(in) , optional :: cvar
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    if ( present(ivar) ) then
      incstat = nf90_put_att(ncid%ncid,ivar,aname,aval)
    else if ( present(cvar) ) then
      ivarid = searchvar(ncid,cvar)
      if ( ivarid > 0 ) then
        incstat = nf90_put_att(ncid%ncid,ivarid,aname,aval)
      else
        incstat = nf90_enotvar
      end if
    else
      incstat = nf90_put_att(ncid%ncid,nf90_global,aname,aval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
       'Error adding attribute '//aname//' in file '//trim(ncid%fname))
  end subroutine clm_addatt_single

  subroutine clm_addatt_double(ncid,aname,aval,ivar,cvar)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: aname
    integer(rk8) , intent(in) :: aval
    integer(ik4) , intent(in) , optional :: ivar
    character(len=*) , intent(in) , optional :: cvar
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    if ( present(ivar) ) then
      incstat = nf90_put_att(ncid%ncid,ivar,aname,aval)
    else if ( present(cvar) ) then
      ivarid = searchvar(ncid,cvar)
      if ( ivarid > 0 ) then
        incstat = nf90_put_att(ncid%ncid,ivarid,aname,aval)
      else
        incstat = nf90_enotvar
      end if
    else
      incstat = nf90_put_att(ncid%ncid,nf90_global,aname,aval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
       'Error adding attribute '//aname//' in file '//trim(ncid%fname))
  end subroutine clm_addatt_double

  logical function clm_check_dimlen(ncid,dname,ival)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: dname
    integer(ik4) , intent(in) :: ival
    integer(ik4) :: idimid , dlen
    clm_check_dimlen = .false.
    incstat = nf90_inq_dimid(ncid%ncid, dname, idimid)
    if ( incstat /= nf90_noerr ) return
    incstat = nf90_inquire_dimension(ncid%ncid, idimid, len=dlen)
    if ( incstat /= nf90_noerr ) return
    if ( dlen /= ival ) return
    clm_check_dimlen = .true.
  end function clm_check_dimlen

  subroutine clm_check_dims(ncid,ni,nj)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    integer(ik4) , intent(out) :: ni , nj
    call clm_inqdim(ncid,'lon',dlen=ni,lerror=.true.)
    if ( ni < 0 ) then
      call clm_inqdim(ncid,'lsmlon',dlen=ni,lerror=.true.)
      if ( ni < 0 ) then
        call clm_inqdim(ncid,'ni',dlen=ni,lerror=.true.)
        if ( ni < 0 ) then
          call clm_inqdim(ncid,'gridcell',dlen=ni)
          nj = 1
          return
        else
          call clm_inqdim(ncid,'nj',dlen=nj)
        end if
      else
        call clm_inqdim(ncid,'lsmlat',dlen=nj)
      end if
    else
      call clm_inqdim(ncid,'lat',dlen=nj)
    end if
  end subroutine clm_check_dims

  subroutine clm_adddim(ncid,dnam,nd)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: dnam
    integer(ik4) , intent(in) :: nd
    integer(ik4) :: nval
    if ( myid /= iocpu ) return
    nval = nd
    if ( nd == clmvar_unlim ) nval = nf90_unlimited
    if ( searchdim(ncid,dnam) > 0 ) return ! Already here
    call add_dimhash(ncid,dnam)
    incstat = nf90_def_dim(ncid%ncid, dnam, nval, ncid%dimids(ncid%idimlast))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error adding dimension '//dnam//' to file '//trim(ncid%fname))
    ncid%idimlast = ncid%idimlast + 1
  end subroutine clm_adddim

  subroutine clm_addvar_char(ctype,ncid,varname,cdims,long_name,units)
    implicit none
    character , intent(in) :: ctype
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: varname
    character(len=*) , intent(in) , optional , dimension(:) :: cdims
    character(len=*) , intent(in) , optional :: long_name
    character(len=*) , intent(in) , optional :: units
    integer(ik4) :: nd , i , varid
    if ( myid /= iocpu ) return
    call add_varhash(ncid,varname)
    if ( present(cdims) ) then
      nd = size(cdims)
      do i = 1 , nd
        usedims(i) = searchdim(ncid,cdims(i))
      end do
      incstat = nf90_def_var(ncid%ncid,varname,nf90_char,usedims(1:nd),varid)
    else
      incstat = nf90_def_var(ncid%ncid,varname,nf90_char,varid)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error adding variable '//varname//' to file '//trim(ncid%fname))
    if ( present(long_name) ) then
      incstat = nf90_put_att(ncid%ncid, varid, 'long_name',long_name)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error long_name to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(units) ) then
      incstat = nf90_put_att(ncid%ncid, varid, 'units', units)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error units to '//varname//' to file '//trim(ncid%fname))
    end if
    ncid%varids(ncid%ivarlast) = varid
    ncid%ivarlast = ncid%ivarlast + 1
  end subroutine clm_addvar_char

  subroutine clm_addvar_int(itype,ncid,varname,cdims,long_name,units, &
                            cell_method,comment,flag_meanings,missing_value,  &
                            fill_value,flag_values,valid_range)
    implicit none
    integer(ik4) , intent(in) :: itype
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: varname
    character(len=*) , intent(in) , optional , dimension(:) :: cdims
    character(len=*) , intent(in) , optional :: long_name
    character(len=*) , intent(in) , optional :: units
    character(len=*) , intent(in) , optional :: cell_method
    character(len=*) , intent(in) , optional :: comment
    character(len=*) , intent(in) , dimension(:) , optional :: flag_meanings
    integer(ik4) , intent(in) , optional :: missing_value
    integer(ik4) , intent(in) , optional :: fill_value
    integer(ik4) , intent(in) , dimension(:) , optional :: flag_values
    integer(ik4) , intent(in) , dimension(2) , optional :: valid_range
    integer(ik4) :: nd , i , varid
    character(len=256) :: str
    if ( myid /= iocpu ) return
    call add_varhash(ncid,varname)
    if ( present(cdims) ) then
      nd = size(cdims)
      do i = 1 , nd
        usedims(i) = searchdim(ncid,cdims(i))
      end do
      incstat = nf90_def_var(ncid%ncid,varname,nf90_int4,usedims(1:nd),varid)
    else
      incstat = nf90_def_var(ncid%ncid,varname,nf90_int4,varid)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error adding variable '//varname//' to file '//trim(ncid%fname))
    incstat = nf90_put_att(ncid%ncid, varid, 'long_name',long_name)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error long_name to '//varname//' to file '//trim(ncid%fname))
    incstat = nf90_put_att(ncid%ncid, varid, 'units', units)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error units to '//varname//' to file '//trim(ncid%fname))
    if ( present(fill_value) ) then
      incstat = nf90_put_att(ncid%ncid, varid, '_FillValue', fill_value)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error fill_value to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(missing_value) ) then
      incstat = nf90_put_att(ncid%ncid, varid, 'missing_value', missing_value)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error missing_value to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(comment) ) then
      incstat = nf90_put_att(ncid%ncid, varid, 'comment', comment)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error comment to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(cell_method) ) then
      str = 'time: ' // trim(cell_method)
      incstat = nf90_put_att(ncid%ncid, varid, 'cell_method', trim(str))
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error cell_method to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(valid_range) ) then
      incstat = nf90_put_att(ncid%ncid, varid, 'valid_range', valid_range)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error valid_range to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(flag_values) ) then
      incstat = nf90_put_att(ncid%ncid, varid, 'flag_values', flag_values)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error flag_values to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(flag_meanings) ) then
      str = ' '
      do i = 1 , size(flag_meanings)
        str = str//flag_meanings(i)
      end do
      incstat = nf90_put_att(ncid%ncid, varid, 'flag_meanings', trim(str))
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error flag_meanings to '//varname//' to file '//trim(ncid%fname))
    end if
    ncid%varids(ncid%ivarlast) = varid
    ncid%ivarlast = ncid%ivarlast + 1
  end subroutine clm_addvar_int

  subroutine clm_addvar_logical(ltype,ncid,varname,cdims,long_name,units, &
                            cell_method,comment,missing_value,fill_value)
    implicit none
    logical , intent(in) :: ltype
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: varname
    character(len=*) , intent(in) , optional , dimension(:) :: cdims
    character(len=*) , intent(in) , optional :: long_name
    character(len=*) , intent(in) , optional :: units
    character(len=*) , intent(in) , optional :: cell_method
    character(len=*) , intent(in) , optional :: comment
    integer(ik4) , intent(in) , optional :: missing_value
    integer(ik4) , intent(in) , optional :: fill_value
    integer(ik4) :: nd , i , varid
    character(len=256) :: str
    if ( myid /= iocpu ) return
    call add_varhash(ncid,varname)
    if ( present(cdims) ) then
      nd = size(cdims)
      do i = 1 , nd
        usedims(i) = searchdim(ncid,cdims(i))
      end do
      incstat = nf90_def_var(ncid%ncid,varname,nf90_int,usedims(1:nd),varid)
    else
      incstat = nf90_def_var(ncid%ncid,varname,nf90_int,varid)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error adding variable '//varname//' to file '//trim(ncid%fname))
    incstat = nf90_put_att(ncid%ncid, varid, 'long_name',long_name)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error long_name to '//varname//' to file '//trim(ncid%fname))
    incstat = nf90_put_att(ncid%ncid, varid, 'units', units)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error units to '//varname//' to file '//trim(ncid%fname))
    if ( present(fill_value) ) then
      incstat = nf90_put_att(ncid%ncid, varid, '_FillValue', fill_value)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error fill_value to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(missing_value) ) then
      incstat = nf90_put_att(ncid%ncid, varid, 'missing_value', missing_value)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error missing_value to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(comment) ) then
      incstat = nf90_put_att(ncid%ncid, varid, 'comment', comment)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error comment to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(cell_method) ) then
      str = 'time: ' // trim(cell_method)
      incstat = nf90_put_att(ncid%ncid, varid, 'cell_method', trim(str))
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error cell_method to '//varname//' to file '//trim(ncid%fname))
    end if
    incstat = nf90_put_att(ncid%ncid, varid, 'valid_range', (/0, 1/))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error valid_range to '//varname//' to file '//trim(ncid%fname))
    incstat = nf90_put_att(ncid%ncid, varid, 'flag_values', (/0, 1/))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error flag_values to '//varname//' to file '//trim(ncid%fname))
    incstat = nf90_put_att(ncid%ncid, varid, 'flag_meanings', 'FALSE TRUE')
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error flag_meanings to '//varname//' to file '//trim(ncid%fname))
    ncid%varids(ncid%ivarlast) = varid
    ncid%ivarlast = ncid%ivarlast + 1
  end subroutine clm_addvar_logical

  subroutine clm_addvar_real4(rtype,ncid,varname,cdims,long_name,units, &
                            cell_method,comment,flag_meanings,missing_value,  &
                            fill_value,flag_values,valid_range)
    implicit none
    real(rk4) , intent(in) :: rtype
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: varname
    character(len=*) , intent(in) , optional , dimension(:) :: cdims
    character(len=*) , intent(in) , optional :: long_name
    character(len=*) , intent(in) , optional :: units
    character(len=*) , intent(in) , optional :: cell_method
    character(len=*) , intent(in) , optional :: comment
    character(len=*) , intent(in) , dimension(:) , optional :: flag_meanings
    real(rk4) , intent(in) , optional :: missing_value
    real(rk4) , intent(in) , optional :: fill_value
    real(rk4) , intent(in) , dimension(:) , optional :: flag_values
    real(rk4) , intent(in) , dimension(2) , optional :: valid_range
    integer(ik4) :: nd , i , varid
    character(len=256) :: str
    if ( myid /= iocpu ) return
    call add_varhash(ncid,varname)
    if ( present(cdims) ) then
      nd = size(cdims)
      do i = 1 , nd
        usedims(i) = searchdim(ncid,cdims(i))
      end do
      incstat = nf90_def_var(ncid%ncid,varname,nf90_real4,usedims(1:nd),varid)
    else
      incstat = nf90_def_var(ncid%ncid,varname,nf90_real4,varid)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error adding variable '//varname//' to file '//trim(ncid%fname))
    incstat = nf90_put_att(ncid%ncid, varid, 'long_name',long_name)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error long_name to '//varname//' to file '//trim(ncid%fname))
    incstat = nf90_put_att(ncid%ncid, varid, 'units', units)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error units to '//varname//' to file '//trim(ncid%fname))
    if ( present(fill_value) ) then
      incstat = nf90_put_att(ncid%ncid, varid, '_FillValue', fill_value)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error fill_value to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(missing_value) ) then
      incstat = nf90_put_att(ncid%ncid, varid, 'missing_value', missing_value)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error missing_value to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(comment) ) then
      incstat = nf90_put_att(ncid%ncid, varid, 'comment', comment)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error comment to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(cell_method) ) then
      str = 'time: ' // trim(cell_method)
      incstat = nf90_put_att(ncid%ncid, varid, 'cell_method', trim(str))
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error cell_method to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(valid_range) ) then
      incstat = nf90_put_att(ncid%ncid, varid, 'valid_range', valid_range)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error valid_range to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(flag_values) ) then
      incstat = nf90_put_att(ncid%ncid, varid, 'flag_values', flag_values)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error flag_values to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(flag_meanings) ) then
      str = ' '
      do i = 1 , size(flag_meanings)
        str = str//flag_meanings(i)
      end do
      incstat = nf90_put_att(ncid%ncid, varid, 'flag_meanings', trim(str))
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error flag_meanings to '//varname//' to file '//trim(ncid%fname))
    end if
    ncid%varids(ncid%ivarlast) = varid
    ncid%ivarlast = ncid%ivarlast + 1
  end subroutine clm_addvar_real4

  subroutine clm_addvar_real8(rtype,ncid,varname,cdims,long_name,units, &
                            cell_method,comment,flag_meanings,missing_value,  &
                            fill_value,flag_values,valid_range)
    implicit none
    real(rk8) , intent(in) :: rtype
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: varname
    character(len=*) , intent(in) , optional , dimension(:) :: cdims
    character(len=*) , intent(in) , optional :: long_name
    character(len=*) , intent(in) , optional :: units
    character(len=*) , intent(in) , optional :: cell_method
    character(len=*) , intent(in) , optional :: comment
    character(len=*) , intent(in) , dimension(:) , optional :: flag_meanings
    real(rk8) , intent(in) , optional :: missing_value
    real(rk8) , intent(in) , optional :: fill_value
    real(rk8) , intent(in) , dimension(:) , optional :: flag_values
    real(rk8) , intent(in) , dimension(2) , optional :: valid_range
    integer(ik4) :: nd , i , varid
    character(len=256) :: str
    if ( myid /= iocpu ) return
    call add_varhash(ncid,varname)
    if ( present(cdims) ) then
      nd = size(cdims)
      do i = 1 , nd
        usedims(i) = searchdim(ncid,cdims(i))
      end do
      incstat = nf90_def_var(ncid%ncid,varname,nf90_real8,usedims(1:nd),varid)
    else
      incstat = nf90_def_var(ncid%ncid,varname,nf90_real8,varid)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error adding variable '//varname//' to file '//trim(ncid%fname))
    incstat = nf90_put_att(ncid%ncid, varid, 'long_name',long_name)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error long_name to '//varname//' to file '//trim(ncid%fname))
    incstat = nf90_put_att(ncid%ncid, varid, 'units', units)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error units to '//varname//' to file '//trim(ncid%fname))
    if ( present(fill_value) ) then
      incstat = nf90_put_att(ncid%ncid, varid, '_FillValue', fill_value)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error fill_value to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(missing_value) ) then
      incstat = nf90_put_att(ncid%ncid, varid, 'missing_value', missing_value)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error missing_value to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(comment) ) then
      incstat = nf90_put_att(ncid%ncid, varid, 'comment', comment)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error comment to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(cell_method) ) then
      str = 'time: ' // trim(cell_method)
      incstat = nf90_put_att(ncid%ncid, varid, 'cell_method', trim(str))
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error cell_method to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(valid_range) ) then
      incstat = nf90_put_att(ncid%ncid, varid, 'valid_range', valid_range)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error valid_range to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(flag_values) ) then
      incstat = nf90_put_att(ncid%ncid, varid, 'flag_values', flag_values)
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error flag_values to '//varname//' to file '//trim(ncid%fname))
    end if
    if ( present(flag_meanings) ) then
      str = ' '
      do i = 1 , size(flag_meanings)
        str = str//flag_meanings(i)
      end do
      incstat = nf90_put_att(ncid%ncid, varid, 'flag_meanings', trim(str))
      call clm_checkncerr(__FILE__,__LINE__, &
        'Error flag_meanings to '//varname//' to file '//trim(ncid%fname))
    end if
    ncid%varids(ncid%ivarlast) = varid
    ncid%ivarlast = ncid%ivarlast + 1
  end subroutine clm_addvar_real8

  logical function clm_check_dim(ncid,dname)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: dname
    integer(ik4) :: idimid
    incstat = nf90_inq_dimid(ncid%ncid,dname,idimid)
    clm_check_dim = ( incstat == nf90_noerr )
  end function clm_check_dim

  logical function clm_check_var(ncid,vname)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    clm_check_var = ( incstat == nf90_noerr )
  end function clm_check_var

  subroutine clm_closefile(ncid)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    if ( ncid%ncid < 0 ) return
    incstat = nf90_close(ncid%ncid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error closing file '//trim(ncid%fname))
    ncid%ncid = -1
    ncid%fname = ' '
    ncid%idimlast = -1
    ncid%ivarlast = -1
    if ( associated(ncid%i4buf) ) deallocate(ncid%i4buf)
    if ( associated(ncid%r4buf) ) deallocate(ncid%r4buf)
    if ( associated(ncid%r8buf) ) deallocate(ncid%r8buf)
  end subroutine clm_closefile

  subroutine clm_checkncerr(filename,line,arg)
    implicit none
    integer(ik4) , intent(in) :: line
    character(len=8) :: cline
    character(*) , intent(in) :: filename , arg
    if ( incstat /= nf90_noerr ) then
      write (cline,'(i8)') line
      write (stderr,*) nf90_strerror(incstat)
      call die(filename,trim(cline)//':'//arg,incstat)
    end if
  end subroutine clm_checkncerr

  integer(ik4) function hash(text) result(hashed) 
    implicit none 
    character(len=*) , intent(in) :: text
    integer(ik4) , parameter :: magic_numb = z'5d7a9f43'
    integer(ik4) :: i, j 
    hashed = 0
    do i = 1, len_trim(text) 
      j = mod(i-1, 4) * 8 
      hashed = ieor( hashed, ishft( ichar( text(i:i) ), j ) ) 
    end do 
    hashed = abs( ieor( hashed, magic_numb ) ) 
    return 
  end function hash

  subroutine add_dimhash(ncid,dname)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: dname
    ncid%dimhash(ncid%idimlast) = hash(dname)
    ncid%dimname(ncid%idimlast) = dname
  end subroutine add_dimhash

  integer(ik4) function searchdim(ncid,dname)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: dname
    integer(ik4) :: i , hashed
    hashed = hash(dname)
    searchdim = -1
    do i = 1 , ncid%idimlast
      if ( ncid%dimhash(i) == hashed ) then
        searchdim = i
        return
      end if
    end do
  end function searchdim

  subroutine add_varhash(ncid,vname)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    ncid%varhash(ncid%ivarlast) = hash(vname)
    ncid%varname(ncid%ivarlast) = vname
  end subroutine add_varhash

  integer(ik4) function searchvar(ncid,vname)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) :: i , hashed
    hashed = hash(vname)
    searchvar = -1
    do i = 1 , ncid%ivarlast
      if ( ncid%varhash(i) == hashed ) then
        searchvar = i
        return
      end if
    end do
  end function searchvar

  subroutine clm_readvar_text_0d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    character(len=*) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_text_0d

  subroutine clm_readvar_text_1d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    character(len=*) , dimension(:) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_text_1d

  subroutine clm_readvar_logical_0d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) :: rval
    logical , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
    xval = (rval > 0)
  end subroutine clm_readvar_logical_0d

  subroutine clm_readvar_logical_1d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:) , intent(out) :: xval
    integer(ik4) , dimension(:) , allocatable :: rval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    allocate(rval(size(xval)))
    incstat = nf90_get_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
    xval = (rval > 0)
    deallocate(rval)
  end subroutine clm_readvar_logical_1d

  subroutine clm_readvar_logical_2d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:) , intent(out) :: xval
    integer(ik4) , dimension(:,:) , allocatable :: rval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    allocate(rval(size(xval,1),size(xval,2)))
    incstat = nf90_get_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
    xval = (rval > 0)
    deallocate(rval)
  end subroutine clm_readvar_logical_2d

  subroutine clm_readvar_logical_3d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:,:) , intent(out) :: xval
    integer(ik4) , dimension(:,:,:) , allocatable :: rval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    allocate(rval(size(xval,1),size(xval,2),size(xval,3)))
    incstat = nf90_get_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
    xval = (rval > 0)
    deallocate(rval)
  end subroutine clm_readvar_logical_3d

  subroutine clm_readvar_logical_4d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:,:,:) , intent(out) :: xval
    integer(ik4) , dimension(:,:,:,:) , allocatable :: rval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    allocate(rval(size(xval,1),size(xval,2),size(xval,3),size(xval,4)))
    incstat = nf90_get_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
    xval = (rval > 0)
    deallocate(rval)
  end subroutine clm_readvar_logical_4d

  subroutine clm_readvar_integer_0d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_integer_0d

  subroutine clm_readvar_integer_1d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_integer_1d

  subroutine clm_readvar_integer_2d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_integer_2d

  subroutine clm_readvar_integer_3d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:,:) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_integer_3d

  subroutine clm_readvar_integer_4d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:,:,:) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_integer_4d

  subroutine clm_readvar_real4_0d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_real4_0d

  subroutine clm_readvar_real4_1d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_real4_1d

  subroutine clm_readvar_real4_2d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_real4_2d

  subroutine clm_readvar_real4_3d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:,:) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_real4_3d

  subroutine clm_readvar_real4_4d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:,:,:) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_real4_4d

  subroutine clm_readvar_real8_0d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_real8_0d

  subroutine clm_readvar_real8_1d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_real8_1d

  subroutine clm_readvar_real8_2d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_real8_2d

  subroutine clm_readvar_real8_3d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:,:) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_real8_3d

  subroutine clm_readvar_real8_4d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:,:,:) , intent(out) :: xval
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    incstat = nf90_get_var(ncid%ncid,ivarid,xval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readvar_real8_4d

  subroutine clm_readrec_logical_0d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , intent(out) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid
    integer(ik4) , dimension(1) :: rval
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    istart(1) = nt
    icount(1) = 1
    incstat = nf90_get_var(ncid%ncid,ivarid,rval,istart(1:1),icount(1:1))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
    xval = (rval(1) > 0)
  end subroutine clm_readrec_logical_0d

  subroutine clm_readrec_logical_1d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:) , intent(out) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1
    integer(ik4) , dimension(:) , allocatable :: rval
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    nv1 = size(xval,1)
    allocate(rval(nv1))
    istart(2) = nt
    istart(1) = 1
    icount(2) = 1
    icount(1) = nv1
    incstat = nf90_get_var(ncid%ncid,ivarid,rval,istart(1:2),icount(1:2))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
    xval = (rval > 0)
    deallocate(rval)
  end subroutine clm_readrec_logical_1d

  subroutine clm_readrec_logical_2d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:) , intent(out) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) , dimension(:,:) , allocatable :: rval
    integer(ik4) :: ivarid , nv1 , nv2
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    nv1 = size(xval,1)
    nv2 = size(xval,2)
    allocate(rval(nv1,nv2))
    istart(3) = nt
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = nv2
    icount(1) = nv1
    incstat = nf90_get_var(ncid%ncid,ivarid,rval,istart(1:3),icount(1:3))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
    xval = (rval > 0)
    deallocate(rval)
  end subroutine clm_readrec_logical_2d

  subroutine clm_readrec_logical_3d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:,:) , intent(out) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) , dimension(:,:,:) , allocatable :: rval
    integer(ik4) :: ivarid , nv1 , nv2 , nv3
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    nv1 = size(xval,1)
    nv2 = size(xval,2)
    nv3 = size(xval,3)
    allocate(rval(nv1,nv2,nv3))
    istart(4) = nt
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = nv3
    icount(2) = nv2
    icount(1) = nv1
    incstat = nf90_get_var(ncid%ncid,ivarid,rval,istart(1:4),icount(1:4))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
    xval = (rval > 0)
    deallocate(rval)
  end subroutine clm_readrec_logical_3d

  subroutine clm_readrec_integer_0d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , intent(out) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid
    integer(ik4) , dimension(1) :: rval
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    istart(1) = nt
    icount(1) = 1
    incstat = nf90_get_var(ncid%ncid,ivarid,rval,istart(1:1),icount(1:1))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
    xval = rval(1)
  end subroutine clm_readrec_integer_0d

  subroutine clm_readrec_integer_1d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:) , intent(out) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    nv1 = size(xval,1)
    istart(2) = nt
    istart(1) = 1
    icount(2) = 1
    icount(1) = nv1
    incstat = nf90_get_var(ncid%ncid,ivarid,xval,istart(1:2),icount(1:2))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readrec_integer_1d

  subroutine clm_readrec_integer_2d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:) , intent(out) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1 , nv2
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    nv1 = size(xval,1)
    nv2 = size(xval,2)
    istart(3) = nt
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = nv2
    icount(1) = nv1
    incstat = nf90_get_var(ncid%ncid,ivarid,xval,istart(1:3),icount(1:3))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readrec_integer_2d

  subroutine clm_readrec_integer_3d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:,:) , intent(out) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1 , nv2 , nv3
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    nv1 = size(xval,1)
    nv2 = size(xval,2)
    nv3 = size(xval,3)
    istart(4) = nt
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = nv3
    icount(2) = nv2
    icount(1) = nv1
    incstat = nf90_get_var(ncid%ncid,ivarid,xval,istart(1:4),icount(1:4))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readrec_integer_3d

  subroutine clm_readrec_real4_0d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , intent(out) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid
    real(rk4) , dimension(1) :: rval
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    istart(1) = nt
    icount(1) = 1
    incstat = nf90_get_var(ncid%ncid,ivarid,rval,istart(1:1),icount(1:1))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
    xval = rval(1)
  end subroutine clm_readrec_real4_0d

  subroutine clm_readrec_real4_1d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:) , intent(out) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    nv1 = size(xval,1)
    istart(2) = nt
    istart(1) = 1
    icount(2) = 1
    icount(1) = nv1
    incstat = nf90_get_var(ncid%ncid,ivarid,xval,istart(1:2),icount(1:2))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readrec_real4_1d

  subroutine clm_readrec_real4_2d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:) , intent(out) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1 , nv2
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    nv1 = size(xval,1)
    nv2 = size(xval,2)
    istart(3) = nt
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = nv2
    icount(1) = nv1
    incstat = nf90_get_var(ncid%ncid,ivarid,xval,istart(1:3),icount(1:3))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readrec_real4_2d

  subroutine clm_readrec_real4_3d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:,:) , intent(out) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1 , nv2 , nv3
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    nv1 = size(xval,1)
    nv2 = size(xval,2)
    nv3 = size(xval,3)
    istart(4) = nt
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = nv3
    icount(2) = nv2
    icount(1) = nv1
    incstat = nf90_get_var(ncid%ncid,ivarid,xval,istart(1:4),icount(1:4))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readrec_real4_3d

  subroutine clm_readrec_real8_0d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , intent(out) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid
    real(rk8) , dimension(1) :: rval
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    istart(1) = nt
    icount(1) = 1
    incstat = nf90_get_var(ncid%ncid,ivarid,rval,istart(1:1),icount(1:1))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
    xval = rval(1)
  end subroutine clm_readrec_real8_0d

  subroutine clm_readrec_real8_1d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:) , intent(out) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    nv1 = size(xval,1)
    istart(2) = nt
    istart(1) = 1
    icount(2) = 1
    icount(1) = nv1
    incstat = nf90_get_var(ncid%ncid,ivarid,xval,istart(1:2),icount(1:2))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readrec_real8_1d

  subroutine clm_readrec_real8_2d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:) , intent(out) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1 , nv2
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    nv1 = size(xval,1)
    nv2 = size(xval,2)
    istart(3) = nt
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = nv2
    icount(1) = nv1
    incstat = nf90_get_var(ncid%ncid,ivarid,xval,istart(1:3),icount(1:3))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readrec_real8_2d

  subroutine clm_readrec_real8_3d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:,:) , intent(out) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1 , nv2 , nv3
    incstat = nf90_inq_varid(ncid%ncid,vname,ivarid)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error search '//vname//' to file '//trim(ncid%fname))
    nv1 = size(xval,1)
    nv2 = size(xval,2)
    nv3 = size(xval,3)
    istart(4) = nt
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = nv3
    icount(2) = nv2
    icount(1) = nv1
    incstat = nf90_get_var(ncid%ncid,ivarid,xval,istart(1:4),icount(1:4))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_readrec_real8_3d

  subroutine clm_writevar_text_0d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    character(len=*) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_text_0d

  subroutine clm_writevar_text_1d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    character(len=*) , dimension(:) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_text_1d

  subroutine clm_writevar_logical_0d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , intent(in) :: xval
    integer(ik4) , dimension(1) :: rval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      rval = 0
      if ( xval ) rval = 1
      incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_logical_0d

  subroutine clm_writevar_logical_1d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:) , intent(in) :: xval
    integer(ik4) :: ivarid
    integer(ik4) , dimension(:) , allocatable :: rval
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      allocate(rval(size(xval)))
      rval = 0
      where (xval)
        rval = 1
      end where
      incstat = nf90_put_var(ncid%ncid,ivarid,rval)
      deallocate(rval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_logical_1d

  subroutine clm_writevar_logical_2d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:) , intent(in) :: xval
    integer(ik4) , dimension(:,:) , allocatable :: rval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      allocate(rval(size(xval,1),size(xval,2)))
      rval = 0
      where (xval)
        rval = 1
      end where
      incstat = nf90_put_var(ncid%ncid,ivarid,rval)
      deallocate(rval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_logical_2d

  subroutine clm_writevar_logical_3d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:,:) , intent(in) :: xval
    integer(ik4) , dimension(:,:,:) , allocatable :: rval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      allocate(rval(size(xval,1),size(xval,2),size(xval,3)))
      rval = 0
      where (xval)
        rval = 1
      end where
      incstat = nf90_put_var(ncid%ncid,ivarid,rval)
      deallocate(rval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_logical_3d

  subroutine clm_writevar_logical_4d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:,:,:) , intent(in) :: xval
    integer(ik4) , dimension(:,:,:,:) , allocatable :: rval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      allocate(rval(size(xval,1),size(xval,2),size(xval,3),size(xval,4)))
      rval = 0
      where (xval)
        rval = 1
      end where
      incstat = nf90_put_var(ncid%ncid,ivarid,rval)
      deallocate(rval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_logical_4d

  subroutine clm_writevar_integer_0d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_integer_0d

  subroutine clm_writevar_integer_1d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_integer_1d

  subroutine clm_writevar_integer_2d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_integer_2d

  subroutine clm_writevar_integer_3d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:,:) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_integer_3d

  subroutine clm_writevar_integer_4d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:,:,:) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_integer_4d

  subroutine clm_writevar_real4_0d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_real4_0d

  subroutine clm_writevar_real4_1d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_real4_1d

  subroutine clm_writevar_real4_2d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_real4_2d

  subroutine clm_writevar_real4_3d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:,:) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_real4_3d

  subroutine clm_writevar_real4_4d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:,:,:) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_real4_4d

  subroutine clm_writevar_real8_0d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_real8_0d

  subroutine clm_writevar_real8_1d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_real8_1d

  subroutine clm_writevar_real8_2d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_real8_2d

  subroutine clm_writevar_real8_3d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:,:) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_real8_3d

  subroutine clm_writevar_real8_4d(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:,:,:) , intent(in) :: xval
    integer(ik4) :: ivarid
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_real8_4d

  subroutine clm_writerec_logical_0d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid
    integer(ik4) , dimension(1) :: rval
    if ( myid /= iocpu ) return
    istart(1) = nt
    icount(1) = 1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      rval = 0
      if ( xval ) rval = 1
      incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:1),icount(1:1))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_logical_0d

  subroutine clm_writerec_logical_1d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) , dimension(:) , allocatable :: rval
    integer(ik4) :: ivarid , nv1
    if ( myid /= iocpu ) return
    nv1 = size(xval,1)
    istart(2) = nt
    istart(1) = 1
    icount(2) = 1
    icount(1) = nv1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      allocate(rval(nv1))
      rval = 0
      where ( xval )
        rval = 1
      end where
      incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:2),icount(1:2))
      deallocate(rval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_logical_1d

  subroutine clm_writerec_logical_2d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) , dimension(:,:) , allocatable :: rval
    integer(ik4) :: ivarid , nv1 , nv2
    if ( myid /= iocpu ) return
    nv1 = size(xval,1)
    nv2 = size(xval,2)
    istart(3) = nt
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = nv2
    icount(1) = nv1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      allocate(rval(nv1,nv2))
      rval = 0
      where ( xval )
        rval = 1
      end where
      incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:3),icount(1:3))
      deallocate(rval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_logical_2d

  subroutine clm_writerec_logical_3d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:,:) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) , dimension(:,:,:) , allocatable :: rval
    integer(ik4) :: ivarid , nv1 , nv2 , nv3
    if ( myid /= iocpu ) return
    nv1 = size(xval,1)
    nv2 = size(xval,2)
    nv3 = size(xval,3)
    istart(4) = nt
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = nv3
    icount(2) = nv2
    icount(1) = nv1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      allocate(rval(nv1,nv2,nv3))
      rval = 0
      where ( xval )
        rval = 1
      end where
      incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:4),icount(1:4))
      deallocate(rval)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_logical_3d

  subroutine clm_writerec_integer_0d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid
    integer(ik4) , dimension(1) :: rval
    if ( myid /= iocpu ) return
    istart(1) = nt
    icount(1) = 1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      rval(:) = xval
      incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:1),icount(1:1))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_integer_0d

  subroutine clm_writerec_integer_1d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1
    if ( myid /= iocpu ) return
    nv1 = size(xval,1)
    istart(2) = nt
    istart(1) = 1
    icount(2) = 1
    icount(1) = nv1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval,istart(1:2),icount(1:2))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_integer_1d

  subroutine clm_writerec_integer_2d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1 , nv2
    if ( myid /= iocpu ) return
    nv1 = size(xval,1)
    nv2 = size(xval,2)
    istart(3) = nt
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = nv2
    icount(1) = nv1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval,istart(1:3),icount(1:3))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_integer_2d

  subroutine clm_writerec_integer_3d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:,:) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1 , nv2 , nv3
    if ( myid /= iocpu ) return
    nv1 = size(xval,1)
    nv2 = size(xval,2)
    nv3 = size(xval,3)
    istart(4) = nt
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = nv3
    icount(2) = nv2
    icount(1) = nv1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval,istart(1:4),icount(1:4))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_integer_3d

  subroutine clm_writerec_real4_0d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid
    real(rk4) , dimension(1) :: rval
    if ( myid /= iocpu ) return
    istart(1) = nt
    icount(1) = 1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      rval(:) = xval
      incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:1),icount(1:1))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real4_0d

  subroutine clm_writerec_real4_1d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1
    if ( myid /= iocpu ) return
    nv1 = size(xval,1)
    istart(2) = nt
    istart(1) = 1
    icount(2) = 1
    icount(1) = nv1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval,istart(1:2),icount(1:2))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real4_1d

  subroutine clm_writerec_real4_2d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1 , nv2
    if ( myid /= iocpu ) return
    nv1 = size(xval,1)
    nv2 = size(xval,2)
    istart(3) = nt
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = nv2
    icount(1) = nv1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval,istart(1:3),icount(1:3))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real4_2d

  subroutine clm_writerec_real4_3d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:,:) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1 , nv2 , nv3
    if ( myid /= iocpu ) return
    nv1 = size(xval,1)
    nv2 = size(xval,2)
    nv3 = size(xval,3)
    istart(4) = nt
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = nv3
    icount(2) = nv2
    icount(1) = nv1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval,istart(1:4),icount(1:4))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real4_3d

  subroutine clm_writerec_real8_0d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid
    real(rk8) , dimension(1) :: rval
    if ( myid /= iocpu ) return
    istart(1) = nt
    icount(1) = 1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      rval(:) = xval
      incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:1),icount(1:1))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real8_0d

  subroutine clm_writerec_real8_1d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1
    if ( myid /= iocpu ) return
    nv1 = size(xval,1)
    istart(2) = nt
    istart(1) = 1
    icount(2) = 1
    icount(1) = nv1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval,istart(1:2),icount(1:2))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real8_1d

  subroutine clm_writerec_real8_2d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1 , nv2
    if ( myid /= iocpu ) return
    nv1 = size(xval,1)
    nv2 = size(xval,2)
    istart(3) = nt
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = nv2
    icount(1) = nv1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval,istart(1:3),icount(1:3))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real8_2d

  subroutine clm_writerec_real8_3d(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:,:) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , nv1 , nv2 , nv3
    if ( myid /= iocpu ) return
    nv1 = size(xval,1)
    nv2 = size(xval,2)
    nv3 = size(xval,3)
    istart(4) = nt
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = nv3
    icount(2) = nv2
    icount(1) = nv1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      incstat = nf90_put_var(ncid%ncid,ivarid,xval,istart(1:4),icount(1:4))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real8_3d

  subroutine copy_filetype(ncid2,ncid1)
    implicit none
    type(clm_filetype) , intent(in) :: ncid1
    type(clm_filetype) , intent(out) :: ncid2
    ncid2%ncid = ncid1%ncid
    ncid2%fname = ncid1%fname
    ncid2%idimlast = ncid1%idimlast
    ncid2%ivarlast = ncid1%ivarlast
    ncid2%dimids = ncid1%dimids
    ncid2%dimhash = ncid1%dimhash
    ncid2%dimname = ncid1%dimname
    ncid2%varids = ncid1%varids
    ncid2%varhash = ncid1%varhash
    ncid2%varname = ncid1%varname
  end subroutine copy_filetype

  subroutine clm_writevar_logical_0d_par_sg(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , intent(in) :: xval
    integer(ik4) :: ivarid
    integer(ik4) , dimension(1) :: ivar
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      if ( xval ) then
        ivar(1) = 1
      else
        ivar(1) = 0
      end if
      incstat = nf90_put_var(ncid%ncid,ivarid,ivar)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_logical_0d_par_sg

  subroutine clm_writevar_logical_1d_par_sg(ncid,vname,xval,sg)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) :: ivarid , mpierr
    integer(ik4) , dimension(:) , allocatable :: rval
    logical , dimension(:) , allocatable :: lval
    if ( myid == iocpu ) then
      allocate(lval(sg%ns))
      allocate(rval(sg%ns))
      ivarid = searchvar(ncid,vname)
    end if
    call mpi_gatherv(xval,sg%ic(myid+1),mpi_logical, &
                     lval,sg%ic,sg%id,mpi_logical,iocpu,sg%icomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
    if ( myid /= iocpu ) return
    where (lval)
      rval = 1
    else where
      rval = 0
    end where
    deallocate(lval)
    incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
    deallocate(rval)
  end subroutine clm_writevar_logical_1d_par_sg

  subroutine clm_writevar_logical_2d_par_sg(ncid,vname,xval,sg)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) :: ivarid , mpierr , nk , k
    integer(ik4) , dimension(:,:) , allocatable :: rval
    logical , dimension(:) , allocatable :: lval
    nk = size(xval,2)
    if ( myid == iocpu ) then
      allocate(lval(sg%ns))
      allocate(rval(sg%ns,nk))
      ivarid = searchvar(ncid,vname)
    end if
    do k = 1 , nk
      call mpi_gatherv(xval(:,k),sg%ic(myid+1),mpi_logical, &
                       lval,sg%ic,sg%id,mpi_logical,iocpu,sg%icomm,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
      where (lval)
        rval(:,k) = 1
      else where
        rval(:,k) = 0
      end where
    end do
    if ( myid /= iocpu ) return
    deallocate(lval)
    incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
    deallocate(rval)
  end subroutine clm_writevar_logical_2d_par_sg

  subroutine clm_writevar_logical_3d_par_sg(ncid,vname,xval,sg)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) :: ivarid , mpierr , nk , nn , n , k
    integer(ik4) , dimension(:,:,:) , allocatable :: rval
    logical , dimension(:) , allocatable :: lval
    nk = size(xval,2)
    nn = size(xval,3)
    if ( myid == iocpu ) then
      allocate(lval(sg%ns))
      allocate(rval(sg%ns,nk,nn))
      ivarid = searchvar(ncid,vname)
    end if
    do n = 1 , nn
      do k = 1 , nk
        call mpi_gatherv(xval(:,k,n),sg%ic(myid+1),mpi_logical, &
                         lval,sg%ic,sg%id,mpi_logical,iocpu,sg%icomm,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
        end if
        where (lval)
          rval(:,k,n) = 1
        else where
          rval(:,k,n) = 0
        end where
      end do
    end do
    if ( myid /= iocpu ) return
    deallocate(lval)
    incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
    deallocate(rval)
  end subroutine clm_writevar_logical_3d_par_sg

  subroutine clm_writevar_logical_4d_par_sg(ncid,vname,xval,sg)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:,:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) :: ivarid , mpierr , nk , nn , nl , n , k , l
    integer(ik4) , dimension(:,:,:,:) , allocatable :: rval
    logical , dimension(:) , allocatable :: lval
    nk = size(xval,2)
    nn = size(xval,3)
    nl = size(xval,4)
    if ( myid == iocpu ) then
      allocate(lval(sg%ns))
      allocate(rval(sg%ns,nk,nn,nl))
      ivarid = searchvar(ncid,vname)
    end if
    do l = 1 , nl
      do n = 1 , nn
        do k = 1 , nk
          call mpi_gatherv(xval(:,k,n,l),sg%ic(myid+1),mpi_logical, &
                           lval,sg%ic,sg%id,mpi_logical,iocpu,sg%icomm,mpierr)
          if ( mpierr /= mpi_success ) then
            call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
          end if
          where (lval)
            rval(:,k,n,l) = 1
          else where
            rval(:,k,n,l) = 0
          end where
        end do
      end do
    end do
    if ( myid /= iocpu ) return
    deallocate(lval)
    incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
    deallocate(rval)
  end subroutine clm_writevar_logical_4d_par_sg

  subroutine clm_writevar_integer_0d_par_sg(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , intent(in) :: xval
    integer(ik4) :: ivarid
    integer(ik4) , dimension(1) :: ivar
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      ivar(1) = xval
      incstat = nf90_put_var(ncid%ncid,ivarid,ivar)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_integer_0d_par_sg

  subroutine clm_writevar_integer_1d_par_sg(ncid,vname,xval,sg)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) :: ivarid , mpierr
    integer(ik4) , dimension(:) , allocatable :: rval
    if ( myid == iocpu ) then
      allocate(rval(sg%ns))
      ivarid = searchvar(ncid,vname)
    end if
    call mpi_gatherv(xval,sg%ic(myid+1),mpi_integer4, &
                     rval,sg%ic,sg%id,mpi_integer4,   &
                     iocpu,sg%icomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
    if ( myid /= iocpu ) return
    incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
    deallocate(rval)
  end subroutine clm_writevar_integer_1d_par_sg

  subroutine clm_writevar_integer_2d_par_sg(ncid,vname,xval,sg)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) :: ivarid , mpierr , nk , k
    integer(ik4) , dimension(:,:) , allocatable :: rval
    nk = size(xval,2)
    if ( myid == iocpu ) then
      allocate(rval(sg%ns,nk))
      ivarid = searchvar(ncid,vname)
    end if
    do k = 1 , nk
      call mpi_gatherv(xval(:,k),sg%ic(myid+1),mpi_integer4, &
                       rval(:,k),sg%ic,sg%id,mpi_integer4,   &
                       iocpu,sg%icomm,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
    end do
    if ( myid /= iocpu ) return
    incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
    deallocate(rval)
  end subroutine clm_writevar_integer_2d_par_sg

  subroutine clm_writevar_integer_3d_par_sg(ncid,vname,xval,sg)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) :: ivarid , mpierr , nk , nn , n , k
    integer(ik4) , dimension(:,:,:) , allocatable :: rval
    nk = size(xval,2)
    nn = size(xval,3)
    if ( myid == iocpu ) then
      allocate(rval(sg%ns,nk,nn))
      ivarid = searchvar(ncid,vname)
    end if
    do n = 1 , nn
      do k = 1 , nk
        call mpi_gatherv(xval(:,k,n),sg%ic(myid+1),mpi_integer4, &
                         rval(:,k,n),sg%ic,sg%id,mpi_integer4,   &
                         iocpu,sg%icomm,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
        end if
      end do
    end do
    if ( myid /= iocpu ) return
    incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
    deallocate(rval)
  end subroutine clm_writevar_integer_3d_par_sg

  subroutine clm_writevar_integer_4d_par_sg(ncid,vname,xval,sg)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:,:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) :: ivarid , mpierr , nk , nn , nl , n , k , l
    integer(ik4) , dimension(:,:,:,:) , allocatable :: rval
    nk = size(xval,2)
    nn = size(xval,3)
    nl = size(xval,4)
    if ( myid == iocpu ) then
      allocate(rval(sg%ns,nk,nn,nl))
      ivarid = searchvar(ncid,vname)
    end if
    do l = 1 , nl
      do n = 1 , nn
        do k = 1 , nk
          call mpi_gatherv(xval(:,k,n,l),sg%ic(myid+1),mpi_integer4, &
                           rval(:,k,n,l),sg%ic,sg%id,mpi_integer4,   &
                           iocpu,sg%icomm,mpierr)
          if ( mpierr /= mpi_success ) then
            call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
          end if
        end do
      end do
    end do
    if ( myid /= iocpu ) return
    incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
    deallocate(rval)
  end subroutine clm_writevar_integer_4d_par_sg

  subroutine clm_writevar_real4_0d_par_sg(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , intent(in) :: xval
    integer(ik4) :: ivarid
    real(rk4) , dimension(1) :: ivar
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      ivar(1) = xval
      incstat = nf90_put_var(ncid%ncid,ivarid,ivar)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_real4_0d_par_sg

  subroutine clm_writevar_real4_1d_par_sg(ncid,vname,xval,sg)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) :: ivarid , mpierr
    real(rk4) , dimension(:) , allocatable :: rval
    if ( myid == iocpu ) then
      allocate(rval(sg%ns))
      ivarid = searchvar(ncid,vname)
    end if
    call mpi_gatherv(xval,sg%ic(myid+1),mpi_real4, &
                     rval,sg%ic,sg%id,mpi_real4,   &
                     iocpu,sg%icomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
    if ( myid /= iocpu ) return
    incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
    deallocate(rval)
  end subroutine clm_writevar_real4_1d_par_sg

  subroutine clm_writevar_real4_2d_par_sg(ncid,vname,xval,sg)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) :: ivarid , mpierr , nk , k
    real(rk4) , dimension(:,:) , allocatable :: rval
    nk = size(xval,2)
    if ( myid == iocpu ) then
      allocate(rval(sg%ns,nk))
      ivarid = searchvar(ncid,vname)
    end if
    do k = 1 , nk
      call mpi_gatherv(xval(:,k),sg%ic(myid+1),mpi_real4, &
                       rval(:,k),sg%ic,sg%id,mpi_real4,   &
                       iocpu,sg%icomm,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
    end do
    if ( myid /= iocpu ) return
    incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
    deallocate(rval)
  end subroutine clm_writevar_real4_2d_par_sg

  subroutine clm_writevar_real4_3d_par_sg(ncid,vname,xval,sg)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) :: ivarid , mpierr , nk , nn , n , k
    real(rk4) , dimension(:,:,:) , allocatable :: rval
    nk = size(xval,2)
    nn = size(xval,3)
    if ( myid == iocpu ) then
      allocate(rval(sg%ns,nk,nn))
      ivarid = searchvar(ncid,vname)
    end if
    do n = 1 , nn
      do k = 1 , nk
        call mpi_gatherv(xval(:,k,n),sg%ic(myid+1),mpi_real4, &
                         rval(:,k,n),sg%ic,sg%id,mpi_real4,   &
                         iocpu,sg%icomm,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
        end if
      end do
    end do
    if ( myid /= iocpu ) return
    incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
    deallocate(rval)
  end subroutine clm_writevar_real4_3d_par_sg

  subroutine clm_writevar_real4_4d_par_sg(ncid,vname,xval,sg)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:,:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) :: ivarid , mpierr , nk , nn , nl , n , k , l
    real(rk4) , dimension(:,:,:,:) , allocatable :: rval
    nk = size(xval,2)
    nn = size(xval,3)
    nl = size(xval,4)
    if ( myid == iocpu ) then
      allocate(rval(sg%ns,nk,nn,nl))
      ivarid = searchvar(ncid,vname)
    end if
    do l = 1 , nl
      do n = 1 , nn
        do k = 1 , nk
          call mpi_gatherv(xval(:,k,n,l),sg%ic(myid+1),mpi_real4, &
                           rval(:,k,n,l),sg%ic,sg%id,mpi_real4,   &
                           iocpu,sg%icomm,mpierr)
          if ( mpierr /= mpi_success ) then
            call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
          end if
        end do
      end do
    end do
    if ( myid /= iocpu ) return
    incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
    deallocate(rval)
  end subroutine clm_writevar_real4_4d_par_sg

  subroutine clm_writevar_real8_0d_par_sg(ncid,vname,xval)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , intent(in) :: xval
    integer(ik4) :: ivarid
    real(rk8) , dimension(1) :: ivar
    if ( myid /= iocpu ) return
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      ivar(1) = xval
      incstat = nf90_put_var(ncid%ncid,ivarid,ivar)
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_real8_0d_par_sg

  subroutine clm_writevar_real8_1d_par_sg(ncid,vname,xval,sg)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) :: ivarid , mpierr
    real(rk8) , dimension(:) , allocatable :: rval
    if ( myid == iocpu ) then
      allocate(rval(sg%ns))
      ivarid = searchvar(ncid,vname)
    end if
    call mpi_gatherv(xval,sg%ic(myid+1),mpi_real8, &
                     rval,sg%ic,sg%id,mpi_real8,   &
                     iocpu,sg%icomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
    if ( myid /= iocpu ) return
    incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
    deallocate(rval)
  end subroutine clm_writevar_real8_1d_par_sg

  subroutine clm_writevar_real8_2d_par_sg(ncid,vname,xval,sg)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) :: ivarid , mpierr , nk , k
    real(rk8) , dimension(:,:) , allocatable :: rval
    nk = size(xval,2)
    if ( myid == iocpu ) then
      allocate(rval(sg%ns,nk))
      ivarid = searchvar(ncid,vname)
    end if
    do k = 1 , nk
      call mpi_gatherv(xval(:,k),sg%ic(myid+1),mpi_real8, &
                       rval(:,k),sg%ic,sg%id,mpi_real8,   &
                       iocpu,sg%icomm,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
    end do
    if ( myid /= iocpu ) return
    incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
    deallocate(rval)
  end subroutine clm_writevar_real8_2d_par_sg

  subroutine clm_writevar_real8_3d_par_sg(ncid,vname,xval,sg)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) :: ivarid , mpierr , nk , nn , n , k
    real(rk8) , dimension(:,:,:) , allocatable :: rval
    nk = size(xval,2)
    nn = size(xval,3)
    if ( myid == iocpu ) then
      allocate(rval(sg%ns,nk,nn))
      ivarid = searchvar(ncid,vname)
    end if
    do n = 1 , nn
      do k = 1 , nk
        call mpi_gatherv(xval(:,k,n),sg%ic(myid+1),mpi_real8, &
                         rval(:,k,n),sg%ic,sg%id,mpi_real8,   &
                         iocpu,sg%icomm,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
        end if
      end do
    end do
    if ( myid /= iocpu ) return
    incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
    deallocate(rval)
  end subroutine clm_writevar_real8_3d_par_sg

  subroutine clm_writevar_real8_4d_par_sg(ncid,vname,xval,sg)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:,:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) :: ivarid , mpierr , nk , nn , nl , n , k , l
    real(rk8) , dimension(:,:,:,:) , allocatable :: rval
    nk = size(xval,2)
    nn = size(xval,3)
    nl = size(xval,4)
    if ( myid == iocpu ) then
      allocate(rval(sg%ns,nk,nn,nl))
      ivarid = searchvar(ncid,vname)
    end if
    do l = 1 , nl
      do n = 1 , nn
        do k = 1 , nk
          call mpi_gatherv(xval(:,k,n,l),sg%ic(myid+1),mpi_real8, &
                           rval(:,k,n,l),sg%ic,sg%id,mpi_real8,   &
                           iocpu,sg%icomm,mpierr)
          if ( mpierr /= mpi_success ) then
            call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
          end if
        end do
      end do
    end do
    if ( myid /= iocpu ) return
    incstat = nf90_put_var(ncid%ncid,ivarid,rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error write '//vname//' to file '//trim(ncid%fname))
    deallocate(rval)
  end subroutine clm_writevar_real8_4d_par_sg

  subroutine clm_writerec_logical_0d_par_sg(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid , mpierr
    integer(ik4) , dimension(1) :: rval
    if ( myid /= iocpu ) return
    istart(1) = nt
    icount(1) = 1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      rval = 0
      if ( xval ) rval = 1
      incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:1),icount(1:1))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_logical_0d_par_sg

  subroutine clm_writerec_logical_1d_par_sg(ncid,vname,xval,sg,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) , intent(in) :: nt
    logical , dimension(:) , allocatable :: lval
    integer(ik4) , dimension(:) , allocatable :: rval
    integer(ik4) :: ivarid , mpierr
    if ( myid == iocpu ) then
      allocate(lval(sg%ns))
      allocate(rval(sg%ns))
      ivarid = searchvar(ncid,vname)
    end if
    call mpi_gatherv(xval,sg%ic(myid+1),mpi_logical, &
                     lval,sg%ic,sg%id,mpi_logical,   &
                     iocpu,sg%icomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
    if ( myid /= iocpu ) return
    where ( lval )
      rval = 1
    else where
      rval = 0
    end where
    deallocate(lval)
    istart(2) = nt
    istart(1) = 1
    icount(2) = 1
    icount(1) = sg%ns
    incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:2),icount(1:2))
    deallocate(rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_logical_1d_par_sg

  subroutine clm_writerec_logical_2d_par_sg(ncid,vname,xval,sg,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) , intent(in) :: nt
    logical , dimension(:) , allocatable :: lval
    integer(ik4) , dimension(:,:) , allocatable :: rval
    integer(ik4) :: ivarid , mpierr , nk , k
    nk = size(xval,2)
    if ( myid == iocpu ) then
      allocate(lval(sg%ns))
      allocate(rval(sg%ns,nk))
      ivarid = searchvar(ncid,vname)
    end if
    do k = 1 , nk
      call mpi_gatherv(xval(:,k),sg%ic(myid+1),mpi_logical, &
                       lval,sg%ic,sg%id,mpi_logical,   &
                       iocpu,sg%icomm,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
      where ( lval )
        rval(:,k) = 1
      else where
        rval(:,k) = 0
      end where
    end do
    if ( myid /= iocpu ) return
    deallocate(lval)
    istart(3) = nt
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = nk
    icount(1) = sg%ns
    incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:3),icount(1:3))
    deallocate(rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_logical_2d_par_sg

  subroutine clm_writerec_logical_3d_par_sg(ncid,vname,xval,sg,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) , intent(in) :: nt
    logical , dimension(:) , allocatable :: lval
    integer(ik4) , dimension(:,:,:) , allocatable :: rval
    integer(ik4) :: ivarid , mpierr , nk , nn , n , k
    nk = size(xval,2)
    nn = size(xval,3)
    if ( myid == iocpu ) then
      allocate(lval(sg%ns))
      allocate(rval(sg%ns,nk,nn))
      ivarid = searchvar(ncid,vname)
    end if
    do n = 1 , nn
      do k = 1 , nk
        call mpi_gatherv(xval(:,k,n),sg%ic(myid+1),mpi_logical, &
                         lval,sg%ic,sg%id,mpi_logical,   &
                         iocpu,sg%icomm,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
        end if
        where ( lval )
          rval(:,k,n) = 1
        else where
          rval(:,k,n) = 0
        end where
      end do
    end do
    if ( myid /= iocpu ) return
    deallocate(lval)
    istart(4) = nt
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = nn
    icount(2) = nk
    icount(1) = sg%ns
    incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:4),icount(1:4))
    deallocate(rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_logical_3d_par_sg

  subroutine clm_writerec_integer_0d_par_sg(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid
    integer(ik4) , dimension(1) :: rval
    if ( myid /= iocpu ) return
    istart(1) = nt
    icount(1) = 1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      rval(1) = xval
      incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:1),icount(1:1))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_integer_0d_par_sg

  subroutine clm_writerec_integer_1d_par_sg(ncid,vname,xval,sg,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) , intent(in) :: nt
    integer(ik4) , dimension(:) , allocatable :: rval
    integer(ik4) :: ivarid , mpierr
    if ( myid == iocpu ) then
      allocate(rval(sg%ns))
      ivarid = searchvar(ncid,vname)
    end if
    call mpi_gatherv(xval,sg%ic(myid+1),mpi_integer4, &
                     rval,sg%ic,sg%id,mpi_integer4,   &
                     iocpu,sg%icomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
    if ( myid /= iocpu ) return
    istart(2) = nt
    istart(1) = 1
    icount(2) = 1
    icount(1) = sg%ns
    incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:2),icount(1:2))
    deallocate(rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_integer_1d_par_sg

  subroutine clm_writerec_integer_2d_par_sg(ncid,vname,xval,sg,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) , intent(in) :: nt
    integer(ik4) , dimension(:,:) , allocatable :: rval
    integer(ik4) :: ivarid , mpierr , nk , k
    nk = size(xval,2)
    if ( myid == iocpu ) then
      allocate(rval(sg%ns,nk))
      ivarid = searchvar(ncid,vname)
    end if
    do k = 1 , nk
      call mpi_gatherv(xval(:,k),sg%ic(myid+1),mpi_integer4, &
                       rval(:,k),sg%ic,sg%id,mpi_integer4,   &
                       iocpu,sg%icomm,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
    end do
    if ( myid /= iocpu ) return
    istart(3) = nt
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = nk
    icount(1) = sg%ns
    incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:3),icount(1:3))
    deallocate(rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_integer_2d_par_sg

  subroutine clm_writerec_integer_3d_par_sg(ncid,vname,xval,sg,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) , intent(in) :: nt
    integer(ik4) , dimension(:,:,:) , allocatable :: rval
    integer(ik4) :: ivarid , mpierr , nk , nn , n , k
    nk = size(xval,2)
    nn = size(xval,3)
    if ( myid == iocpu ) then
      allocate(rval(sg%ns,nk,nn))
      ivarid = searchvar(ncid,vname)
    end if
    do n = 1 , nn
      do k = 1 , nk
        call mpi_gatherv(xval(:,k,n),sg%ic(myid+1),mpi_integer4, &
                         rval(:,k,n),sg%ic,sg%id,mpi_integer4,   &
                         iocpu,sg%icomm,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
        end if
      end do
    end do
    if ( myid /= iocpu ) return
    istart(4) = nt
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = nn
    icount(2) = nk
    icount(1) = sg%ns
    incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:4),icount(1:4))
    deallocate(rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_integer_3d_par_sg

  subroutine clm_writerec_real4_0d_par_sg(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid
    real(rk4) , dimension(1) :: rval
    if ( myid /= iocpu ) return
    istart(1) = nt
    icount(1) = 1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      rval(1) = xval
      incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:1),icount(1:1))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real4_0d_par_sg

  subroutine clm_writerec_real4_1d_par_sg(ncid,vname,xval,sg,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) , intent(in) :: nt
    real(rk4) , dimension(:) , allocatable :: rval
    integer(ik4) :: ivarid , mpierr
    if ( myid == iocpu ) then
      allocate(rval(sg%ns))
      ivarid = searchvar(ncid,vname)
    end if
    call mpi_gatherv(xval,sg%ic(myid+1),mpi_real4, &
                     rval,sg%ic,sg%id,mpi_real4,   &
                     iocpu,sg%icomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
    if ( myid /= iocpu ) return
    istart(2) = nt
    istart(1) = 1
    icount(2) = 1
    icount(1) = sg%ns
    incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:2),icount(1:2))
    deallocate(rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real4_1d_par_sg

  subroutine clm_writerec_real4_2d_par_sg(ncid,vname,xval,sg,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) , intent(in) :: nt
    real(rk4) , dimension(:,:) , allocatable :: rval
    integer(ik4) :: ivarid , mpierr , nk , k
    nk = size(xval,2)
    if ( myid == iocpu ) then
      allocate(rval(sg%ns,nk))
      ivarid = searchvar(ncid,vname)
    end if
    do k = 1 , nk
      call mpi_gatherv(xval(:,k),sg%ic(myid+1),mpi_real4, &
                       rval(:,k),sg%ic,sg%id,mpi_real4,   &
                       iocpu,sg%icomm,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
    end do
    if ( myid /= iocpu ) return
    istart(3) = nt
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = nk
    icount(1) = sg%ns
    incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:3),icount(1:3))
    deallocate(rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real4_2d_par_sg

  subroutine clm_writerec_real4_3d_par_sg(ncid,vname,xval,sg,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) , intent(in) :: nt
    real(rk4) , dimension(:,:,:) , allocatable :: rval
    integer(ik4) :: ivarid , mpierr , nk , nn , n , k
    nk = size(xval,2)
    nn = size(xval,3)
    if ( myid == iocpu ) then
      allocate(rval(sg%ns,nk,nn))
      ivarid = searchvar(ncid,vname)
    end if
    do n = 1 , nn
      do k = 1 , nk
        call mpi_gatherv(xval(:,k,n),sg%ic(myid+1),mpi_real4, &
                         rval(:,k,n),sg%ic,sg%id,mpi_real4,   &
                         iocpu,sg%icomm,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
        end if
      end do
    end do
    if ( myid /= iocpu ) return
    istart(4) = nt
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = nn
    icount(2) = nk
    icount(1) = sg%ns
    incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:4),icount(1:4))
    deallocate(rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real4_3d_par_sg

  subroutine clm_writerec_real8_0d_par_sg(ncid,vname,xval,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , intent(in) :: xval
    integer(ik4) , intent(in) :: nt
    integer(ik4) :: ivarid
    real(rk8) , dimension(1) :: rval
    if ( myid /= iocpu ) return
    istart(1) = nt
    icount(1) = 1
    ivarid = searchvar(ncid,vname)
    if ( ivarid < 0 ) then
      incstat = nf90_enotvar
    else
      rval(1) = xval
      incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:1),icount(1:1))
    end if
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real8_0d_par_sg

  subroutine clm_writerec_real8_1d_par_sg(ncid,vname,xval,sg,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) , intent(in) :: nt
    real(rk8) , dimension(:) , allocatable :: rval
    integer(ik4) :: ivarid , mpierr
    if ( myid == iocpu ) then
      allocate(rval(sg%ns))
      ivarid = searchvar(ncid,vname)
    end if
    call mpi_gatherv(xval,sg%ic(myid+1),mpi_real8, &
                     rval,sg%ic,sg%id,mpi_real8,   &
                     iocpu,sg%icomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
    if ( myid /= iocpu ) return
    istart(2) = nt
    istart(1) = 1
    icount(2) = 1
    icount(1) = sg%ns
    incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:2),icount(1:2))
    deallocate(rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real8_1d_par_sg

  subroutine clm_writerec_real8_2d_par_sg(ncid,vname,xval,sg,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) , intent(in) :: nt
    real(rk8) , dimension(:,:) , allocatable :: rval
    integer(ik4) :: ivarid , mpierr , nk , k
    nk = size(xval,2)
    if ( myid == iocpu ) then
      allocate(rval(sg%ns,nk))
      ivarid = searchvar(ncid,vname)
    end if
    do k = 1 , nk
      call mpi_gatherv(xval(:,k),sg%ic(myid+1),mpi_real8, &
                       rval(:,k),sg%ic,sg%id,mpi_real8,   &
                       iocpu,sg%icomm,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
    end do
    if ( myid /= iocpu ) return
    istart(3) = nt
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = nk
    icount(1) = sg%ns
    incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:3),icount(1:3))
    deallocate(rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real8_2d_par_sg

  subroutine clm_writerec_real8_3d_par_sg(ncid,vname,xval,sg,nt)
    implicit none
    type(clm_filetype) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:,:) , intent(in) :: xval
    type(subgrid_type) , intent(in) :: sg
    integer(ik4) , intent(in) :: nt
    real(rk8) , dimension(:,:,:) , allocatable :: rval
    integer(ik4) :: ivarid , mpierr , nk , nn , n , k
    nk = size(xval,2)
    nn = size(xval,3)
    if ( myid == iocpu ) then
      allocate(rval(sg%ns,nk,nn))
      ivarid = searchvar(ncid,vname)
    end if
    do n = 1 , nn
      do k = 1 , nk
        call mpi_gatherv(xval(:,k,n),sg%ic(myid+1),mpi_real8, &
                         rval(:,k,n),sg%ic,sg%id,mpi_real8,   &
                         iocpu,sg%icomm,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
        end if
      end do
    end do
    if ( myid /= iocpu ) return
    istart(4) = nt
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = nn
    icount(2) = nk
    icount(1) = sg%ns
    incstat = nf90_put_var(ncid%ncid,ivarid,rval,istart(1:4),icount(1:4))
    deallocate(rval)
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real8_3d_par_sg

  subroutine clm_writevar_logical_2d_par_gg(ncid,vname,xval,gg)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    logical , dimension(:) , allocatable :: lval
    integer(ik4) :: i , j , ib , ivarid , mpierr
    if ( myid == iocpu ) then
      allocate(lval(numg))
      ivarid = searchvar(ncid,vname)
    end if
    call mpi_gatherv(xval,gg%gc(myid+1),mpi_logical, &
                     lval,gg%gc,gg%gd,mpi_logical,   &
                     iocpu,gg%icomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
    if ( myid /= iocpu ) return
    ib = 1
    do i = iout1 , iout2
      do j = jout1 , jout2
        if ( gg%gcmask(j,i) ) then
          if ( lval(ib) ) then
            ncid%i4buf(j,i) = 1
          else
            ncid%i4buf(j,i) = 0
          end if
          ib = ib + 1
        end if
      end do
    end do
    deallocate(lval)
    istart(2) = 1
    istart(1) = 1
    icount(2) = niout
    icount(1) = njout
    incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf,istart(1:2),icount(1:2))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_logical_2d_par_gg

  subroutine clm_writevar_logical_3d_par_gg(ncid,vname,xval,gg)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    logical , dimension(:) , allocatable :: lval
    integer(ik4) :: i , j , k , nk , ib , ivarid , mpierr
    nk = size(xval,2)
    if ( myid == iocpu ) then
      allocate(lval(numg))
      ivarid = searchvar(ncid,vname)
    end if
    do k = 1 , nk
      call mpi_gatherv(xval(:,k),gg%gc(myid+1),mpi_logical, &
                       lval,gg%gc,gg%gd,mpi_logical,   &
                       iocpu,gg%icomm,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
      if ( myid == iocpu ) then
        ib = 1
        do i = iout1 , iout2
          do j = jout1 , jout2
            if ( gg%gcmask(j,i) ) then
              if ( lval(ib) ) then
                ncid%i4buf(j,i) = 1
              else
                ncid%i4buf(j,i) = 0
              end if
              ib = ib + 1
            end if
          end do
        end do
        istart(3) = k
        istart(2) = 1
        istart(1) = 1
        icount(3) = 1
        icount(2) = niout
        icount(1) = njout
        incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf, &
                               istart(1:3),icount(1:3))
        call clm_checkncerr(__FILE__,__LINE__, &
          'Error read '//vname//' to file '//trim(ncid%fname))
      end if
    end do
    if ( myid /= iocpu ) return
    deallocate(lval)
  end subroutine clm_writevar_logical_3d_par_gg

  subroutine clm_writevar_logical_4d_par_gg(ncid,vname,xval,gg)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:,:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    logical , dimension(:) , allocatable :: lval
    integer(ik4) :: i , j , k , l , nk , nl , ib , ivarid , mpierr
    nk = size(xval,2)
    nl = size(xval,3)
    if ( myid == iocpu ) then
      allocate(lval(numg))
      ivarid = searchvar(ncid,vname)
    end if
    do l = 1 , nl
      do k = 1 , nk
        call mpi_gatherv(xval(:,k,l),gg%gc(myid+1),mpi_logical, &
                         lval,gg%gc,gg%gd,mpi_logical,   &
                         iocpu,gg%icomm,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
        end if
        if ( myid == iocpu ) then
          ib = 1
          do i = iout1 , iout2
            do j = jout1 , jout2
              if ( gg%gcmask(j,i) ) then
                if ( lval(ib) ) then
                  ncid%i4buf(j,i) = 1
                else
                  ncid%i4buf(j,i) = 0
                end if
                ib = ib + 1
              end if
            end do
          end do
          istart(4) = l
          istart(3) = k
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = 1
          icount(2) = niout
          icount(1) = njout
          incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf, &
                                 istart(1:4),icount(1:4))
          call clm_checkncerr(__FILE__,__LINE__, &
            'Error read '//vname//' to file '//trim(ncid%fname))
        end if
      end do
    end do
    if ( myid /= iocpu ) return
    deallocate(lval)
  end subroutine clm_writevar_logical_4d_par_gg

  subroutine clm_writevar_integer_2d_par_gg(ncid,vname,xval,gg)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    integer(ik4) , dimension(:) , allocatable :: ival
    integer(ik4) :: i , j , ib , ivarid , mpierr
    if ( myid == iocpu ) then
      allocate(ival(numg))
      ivarid = searchvar(ncid,vname)
    end if
    call mpi_gatherv(xval,gg%gc(myid+1),mpi_integer4, &
                     ival,gg%gc,gg%gd,mpi_integer4,   &
                     iocpu,gg%icomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
    if ( myid /= iocpu ) return
    ib = 1
    do i = iout1 , iout2
      do j = jout1 , jout2
        if ( gg%gcmask(j,i) ) then
          ncid%i4buf(j,i) = ival(ib)
          ib = ib + 1
        end if
      end do
    end do
    deallocate(ival)
    istart(2) = 1
    istart(1) = 1
    icount(2) = niout
    icount(1) = njout
    incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf,istart(1:2),icount(1:2))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_integer_2d_par_gg

  subroutine clm_writevar_integer_3d_par_gg(ncid,vname,xval,gg)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    integer(ik4) , dimension(:) , allocatable :: ival
    integer(ik4) :: i , j , k , nk , ib , ivarid , mpierr
    nk = size(xval,2)
    if ( myid == iocpu ) then
      allocate(ival(numg))
      ivarid = searchvar(ncid,vname)
    end if
    do k = 1 , nk
      call mpi_gatherv(xval(:,k),gg%gc(myid+1),mpi_integer4, &
                       ival,gg%gc,gg%gd,mpi_integer4,   &
                       iocpu,gg%icomm,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
      if ( myid == iocpu ) then
        ib = 1
        do i = iout1 , iout2
          do j = jout1 , jout2
            if ( gg%gcmask(j,i) ) then
              ncid%i4buf(j,i) = ival(ib)
              ib = ib + 1
            end if
          end do
        end do
        istart(3) = k
        istart(2) = 1
        istart(1) = 1
        icount(3) = 1
        icount(2) = niout
        icount(1) = njout
        incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf, &
                               istart(1:3),icount(1:3))
        call clm_checkncerr(__FILE__,__LINE__, &
          'Error read '//vname//' to file '//trim(ncid%fname))
      end if
    end do
    if ( myid /= iocpu ) return
    deallocate(ival)
  end subroutine clm_writevar_integer_3d_par_gg

  subroutine clm_writevar_integer_4d_par_gg(ncid,vname,xval,gg)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:,:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    integer(ik4) , dimension(:) , allocatable :: ival
    integer(ik4) :: i , j , k , l , nk , nl , ib , ivarid , mpierr
    nk = size(xval,2)
    nl = size(xval,3)
    if ( myid == iocpu ) then
      allocate(ival(numg))
      ivarid = searchvar(ncid,vname)
    end if
    do l = 1 , nl
      do k = 1 , nk
        call mpi_gatherv(xval(:,k,l),gg%gc(myid+1),mpi_integer4, &
                         ival,gg%gc,gg%gd,mpi_integer4,          &
                         iocpu,gg%icomm,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
        end if
        if ( myid == iocpu ) then
          ib = 1
          do i = iout1 , iout2
            do j = jout1 , jout2
              if ( gg%gcmask(j,i) ) then
                ncid%i4buf(j,i) = ival(ib)
                ib = ib + 1
              end if
            end do
          end do
          istart(4) = l
          istart(3) = k
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = 1
          icount(2) = niout
          icount(1) = njout
          incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf, &
                                 istart(1:4),icount(1:4))
          call clm_checkncerr(__FILE__,__LINE__, &
            'Error read '//vname//' to file '//trim(ncid%fname))
        end if
      end do
    end do
    if ( myid /= iocpu ) return
    deallocate(ival)
  end subroutine clm_writevar_integer_4d_par_gg

  subroutine clm_writevar_real4_2d_par_gg(ncid,vname,xval,gg)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    real(rk4) , dimension(:) , allocatable :: ival
    integer(ik4) :: i , j , ib , ivarid , mpierr
    if ( myid == iocpu ) then
      allocate(ival(numg))
      ivarid = searchvar(ncid,vname)
    end if
    call mpi_gatherv(xval,gg%gc(myid+1),mpi_real4, &
                     ival,gg%gc,gg%gd,mpi_real4,   &
                     iocpu,gg%icomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
    if ( myid /= iocpu ) return
    ib = 1
    do i = iout1 , iout2
      do j = jout1 , jout2
        if ( gg%gcmask(j,i) ) then
          ncid%i4buf(j,i) = ival(ib)
          ib = ib + 1
        end if
      end do
    end do
    deallocate(ival)
    istart(2) = 1
    istart(1) = 1
    icount(2) = niout
    icount(1) = njout
    incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf,istart(1:2),icount(1:2))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_real4_2d_par_gg

  subroutine clm_writevar_real4_3d_par_gg(ncid,vname,xval,gg)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    real(rk4) , dimension(:) , allocatable :: ival
    integer(ik4) :: i , j , k , nk , ib , ivarid , mpierr
    nk = size(xval,2)
    if ( myid == iocpu ) then
      allocate(ival(numg))
      ivarid = searchvar(ncid,vname)
    end if
    do k = 1 , nk
      call mpi_gatherv(xval(:,k),gg%gc(myid+1),mpi_real4, &
                       ival,gg%gc,gg%gd,mpi_real4,   &
                       iocpu,gg%icomm,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
      if ( myid == iocpu ) then
        ib = 1
        do i = iout1 , iout2
          do j = jout1 , jout2
            if ( gg%gcmask(j,i) ) then
              ncid%i4buf(j,i) = ival(ib)
              ib = ib + 1
            end if
          end do
        end do
        istart(3) = k
        istart(2) = 1
        istart(1) = 1
        icount(3) = 1
        icount(2) = niout
        icount(1) = njout
        incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf, &
                               istart(1:3),icount(1:3))
        call clm_checkncerr(__FILE__,__LINE__, &
          'Error read '//vname//' to file '//trim(ncid%fname))
      end if
    end do
    if ( myid /= iocpu ) return
    deallocate(ival)
  end subroutine clm_writevar_real4_3d_par_gg

  subroutine clm_writevar_real4_4d_par_gg(ncid,vname,xval,gg)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:,:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    real(rk4) , dimension(:) , allocatable :: ival
    integer(ik4) :: i , j , k , l , nk , nl , ib , ivarid , mpierr
    nk = size(xval,2)
    nl = size(xval,3)
    if ( myid == iocpu ) then
      allocate(ival(numg))
      ivarid = searchvar(ncid,vname)
    end if
    do l = 1 , nl
      do k = 1 , nk
        call mpi_gatherv(xval(:,k,l),gg%gc(myid+1),mpi_real4, &
                         ival,gg%gc,gg%gd,mpi_real4,          &
                         iocpu,gg%icomm,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
        end if
        if ( myid == iocpu ) then
          ib = 1
          do i = iout1 , iout2
            do j = jout1 , jout2
              if ( gg%gcmask(j,i) ) then
                ncid%i4buf(j,i) = ival(ib)
                ib = ib + 1
              end if
            end do
          end do
          istart(4) = l
          istart(3) = k
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = 1
          icount(2) = niout
          icount(1) = njout
          incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf, &
                                 istart(1:4),icount(1:4))
          call clm_checkncerr(__FILE__,__LINE__, &
            'Error read '//vname//' to file '//trim(ncid%fname))
        end if
      end do
    end do
    if ( myid /= iocpu ) return
    deallocate(ival)
  end subroutine clm_writevar_real4_4d_par_gg

  subroutine clm_writevar_real8_2d_par_gg(ncid,vname,xval,gg)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    real(rk8) , dimension(:) , allocatable :: ival
    integer(ik4) :: i , j , ib , ivarid , mpierr
    if ( myid == iocpu ) then
      allocate(ival(numg))
      ivarid = searchvar(ncid,vname)
    end if
    call mpi_gatherv(xval,gg%gc(myid+1),mpi_real8, &
                     ival,gg%gc,gg%gd,mpi_real8,   &
                     iocpu,gg%icomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
    if ( myid /= iocpu ) return
    ib = 1
    do i = iout1 , iout2
      do j = jout1 , jout2
        if ( gg%gcmask(j,i) ) then
          ncid%i4buf(j,i) = ival(ib)
          ib = ib + 1
        end if
      end do
    end do
    deallocate(ival)
    istart(2) = 1
    istart(1) = 1
    icount(2) = niout
    icount(1) = njout
    incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf,istart(1:2),icount(1:2))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writevar_real8_2d_par_gg

  subroutine clm_writevar_real8_3d_par_gg(ncid,vname,xval,gg)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    real(rk8) , dimension(:) , allocatable :: ival
    integer(ik4) :: i , j , k , nk , ib , ivarid , mpierr
    nk = size(xval,2)
    if ( myid == iocpu ) then
      allocate(ival(numg))
      ivarid = searchvar(ncid,vname)
    end if
    do k = 1 , nk
      call mpi_gatherv(xval(:,k),gg%gc(myid+1),mpi_real8, &
                       ival,gg%gc,gg%gd,mpi_real8,   &
                       iocpu,gg%icomm,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
      if ( myid == iocpu ) then
        ib = 1
        do i = iout1 , iout2
          do j = jout1 , jout2
            if ( gg%gcmask(j,i) ) then
              ncid%i4buf(j,i) = ival(ib)
              ib = ib + 1
            end if
          end do
        end do
        istart(3) = k
        istart(2) = 1
        istart(1) = 1
        icount(3) = 1
        icount(2) = niout
        icount(1) = njout
        incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf, &
                               istart(1:3),icount(1:3))
        call clm_checkncerr(__FILE__,__LINE__, &
          'Error read '//vname//' to file '//trim(ncid%fname))
      end if
    end do
    if ( myid /= iocpu ) return
    deallocate(ival)
  end subroutine clm_writevar_real8_3d_par_gg

  subroutine clm_writevar_real8_4d_par_gg(ncid,vname,xval,gg)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:,:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    real(rk8) , dimension(:) , allocatable :: ival
    integer(ik4) :: i , j , k , l , nk , nl , ib , ivarid , mpierr
    nk = size(xval,2)
    nl = size(xval,3)
    if ( myid == iocpu ) then
      allocate(ival(numg))
      ivarid = searchvar(ncid,vname)
    end if
    do l = 1 , nl
      do k = 1 , nk
        call mpi_gatherv(xval(:,k,l),gg%gc(myid+1),mpi_real8, &
                         ival,gg%gc,gg%gd,mpi_real8,          &
                         iocpu,gg%icomm,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
        end if
        if ( myid == iocpu ) then
          ib = 1
          do i = iout1 , iout2
            do j = jout1 , jout2
              if ( gg%gcmask(j,i) ) then
                ncid%i4buf(j,i) = ival(ib)
                ib = ib + 1
              end if
            end do
          end do
          istart(4) = l
          istart(3) = k
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = 1
          icount(2) = niout
          icount(1) = njout
          incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf, &
                                 istart(1:4),icount(1:4))
          call clm_checkncerr(__FILE__,__LINE__, &
            'Error read '//vname//' to file '//trim(ncid%fname))
        end if
      end do
    end do
    if ( myid /= iocpu ) return
    deallocate(ival)
  end subroutine clm_writevar_real8_4d_par_gg

  subroutine clm_writerec_logical_2d_par_gg(ncid,vname,xval,gg,nt)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    integer(ik4) , intent(in) :: nt
    logical , dimension(:) , allocatable :: lval
    integer(ik4) :: i , j , ib , ivarid , mpierr
    if ( myid == iocpu ) then
      allocate(lval(numg))
      ivarid = searchvar(ncid,vname)
    end if
    call mpi_gatherv(xval,gg%gc(myid+1),mpi_logical, &
                     lval,gg%gc,gg%gd,mpi_logical,   &
                     iocpu,gg%icomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
    if ( myid /= iocpu ) return
    ib = 1
    do i = iout1 , iout2
      do j = jout1 , jout2
        if ( gg%gcmask(j,i) ) then
          if ( lval(ib) ) then
            ncid%i4buf(j,i) = 1
          else
            ncid%i4buf(j,i) = 0
          end if
          ib = ib + 1
        end if
      end do
    end do
    deallocate(lval)
    istart(3) = nt
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = niout
    icount(1) = njout
    incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf,istart(1:3),icount(1:3))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_logical_2d_par_gg

  subroutine clm_writerec_integer_2d_par_gg(ncid,vname,xval,gg,nt)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    integer(ik4) , intent(in) :: nt
    integer(ik4) , dimension(:) , allocatable :: ival
    integer(ik4) :: i , j , ib , ivarid , mpierr
    if ( myid == iocpu ) then
      allocate(ival(numg))
      ivarid = searchvar(ncid,vname)
    end if
    call mpi_gatherv(xval,gg%gc(myid+1),mpi_integer4, &
                     ival,gg%gc,gg%gd,mpi_integer4,   &
                     iocpu,gg%icomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
    if ( myid /= iocpu ) return
    ib = 1
    do i = iout1 , iout2
      do j = jout1 , jout2
        if ( gg%gcmask(j,i) ) then
          ncid%i4buf(j,i) = ival(ib)
          ib = ib + 1
        end if
      end do
    end do
    deallocate(ival)
    istart(3) = nt
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = niout
    icount(1) = njout
    incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf,istart(1:3),icount(1:3))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_integer_2d_par_gg

  subroutine clm_writerec_logical_3d_par_gg(ncid,vname,xval,gg,nt)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    logical , dimension(:,:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    integer(ik4) , intent(in) :: nt
    logical , dimension(:) , allocatable :: lval
    integer(ik4) :: i , j , k , nk , ib , ivarid , mpierr
    nk = size(xval,2)
    if ( myid == iocpu ) then
      allocate(lval(numg))
      ivarid = searchvar(ncid,vname)
    end if
    do k = 1 , nk
      call mpi_gatherv(xval(:,k),gg%gc(myid+1),mpi_logical, &
                       lval,gg%gc,gg%gd,mpi_logical,   &
                       iocpu,gg%icomm,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
      if ( myid == iocpu ) then
        ib = 1
        do i = iout1 , iout2
          do j = jout1 , jout2
            if ( gg%gcmask(j,i) ) then
              if ( lval(ib) ) then
                ncid%i4buf(j,i) = 1
              else
                ncid%i4buf(j,i) = 0
              end if
              ib = ib + 1
            end if
          end do
        end do
        istart(4) = nt
        istart(3) = k
        istart(2) = 1
        istart(1) = 1
        icount(4) = 1
        icount(3) = 1
        icount(2) = niout
        icount(1) = njout
        incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf, &
                               istart(1:4),icount(1:4))
        call clm_checkncerr(__FILE__,__LINE__, &
          'Error read '//vname//' to file '//trim(ncid%fname))
      end if
    end do
    if ( myid /= iocpu ) return
    deallocate(lval)
  end subroutine clm_writerec_logical_3d_par_gg

  subroutine clm_writerec_integer_3d_par_gg(ncid,vname,xval,gg,nt)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    integer(ik4) , dimension(:,:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    integer(ik4) , intent(in) :: nt
    integer(ik4) , dimension(:) , allocatable :: ival
    integer(ik4) :: i , j , k , nk , ib , ivarid , mpierr
    nk = size(xval,2)
    if ( myid == iocpu ) then
      allocate(ival(numg))
      ivarid = searchvar(ncid,vname)
    end if
    do k = 1 , nk
      call mpi_gatherv(xval(:,k),gg%gc(myid+1),mpi_integer4, &
                       ival,gg%gc,gg%gd,mpi_integer4,   &
                       iocpu,gg%icomm,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
      if ( myid == iocpu ) then
        ib = 1
        do i = iout1 , iout2
          do j = jout1 , jout2
            if ( gg%gcmask(j,i) ) then
              ncid%i4buf(j,i) = ival(ib)
              ib = ib + 1
            end if
          end do
        end do
        istart(4) = nt
        istart(3) = k
        istart(2) = 1
        istart(1) = 1
        icount(4) = 1
        icount(3) = 1
        icount(2) = niout
        icount(1) = njout
        incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf, &
                               istart(1:4),icount(1:4))
        call clm_checkncerr(__FILE__,__LINE__, &
          'Error read '//vname//' to file '//trim(ncid%fname))
      end if
    end do
    if ( myid /= iocpu ) return
    deallocate(ival)
  end subroutine clm_writerec_integer_3d_par_gg

  subroutine clm_writerec_real4_2d_par_gg(ncid,vname,xval,gg,nt)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    integer(ik4) , intent(in) :: nt
    real(rk4) , dimension(:) , allocatable :: ival
    integer(ik4) :: i , j , ib , ivarid , mpierr
    if ( myid == iocpu ) then
      allocate(ival(numg))
      ivarid = searchvar(ncid,vname)
    end if
    call mpi_gatherv(xval,gg%gc(myid+1),mpi_real4, &
                     ival,gg%gc,gg%gd,mpi_real4,   &
                     iocpu,gg%icomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
    if ( myid /= iocpu ) return
    ib = 1
    do i = iout1 , iout2
      do j = jout1 , jout2
        if ( gg%gcmask(j,i) ) then
          ncid%i4buf(j,i) = ival(ib)
          ib = ib + 1
        end if
      end do
    end do
    deallocate(ival)
    istart(3) = nt
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = niout
    icount(1) = njout
    incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf,istart(1:3),icount(1:3))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real4_2d_par_gg

  subroutine clm_writerec_real4_3d_par_gg(ncid,vname,xval,gg,nt)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk4) , dimension(:,:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    integer(ik4) , intent(in) :: nt
    real(rk4) , dimension(:) , allocatable :: ival
    integer(ik4) :: i , j , k , nk , ib , ivarid , mpierr
    nk = size(xval,2)
    if ( myid == iocpu ) then
      allocate(ival(numg))
      ivarid = searchvar(ncid,vname)
    end if
    do k = 1 , nk
      call mpi_gatherv(xval(:,k),gg%gc(myid+1),mpi_real4, &
                       ival,gg%gc,gg%gd,mpi_real4,   &
                       iocpu,gg%icomm,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
      if ( myid == iocpu ) then
        ib = 1
        do i = iout1 , iout2
          do j = jout1 , jout2
            if ( gg%gcmask(j,i) ) then
              ncid%i4buf(j,i) = ival(ib)
              ib = ib + 1
            end if
          end do
        end do
        istart(4) = nt
        istart(3) = k
        istart(2) = 1
        istart(1) = 1
        icount(4) = 1
        icount(3) = 1
        icount(2) = niout
        icount(1) = njout
        incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf, &
                               istart(1:4),icount(1:4))
        call clm_checkncerr(__FILE__,__LINE__, &
          'Error read '//vname//' to file '//trim(ncid%fname))
      end if
    end do
    if ( myid /= iocpu ) return
    deallocate(ival)
  end subroutine clm_writerec_real4_3d_par_gg

  subroutine clm_writerec_real8_2d_par_gg(ncid,vname,xval,gg,nt)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    integer(ik4) , intent(in) :: nt
    real(rk8) , dimension(:) , allocatable :: ival
    integer(ik4) :: i , j , ib , ivarid , mpierr
    if ( myid == iocpu ) then
      allocate(ival(numg))
      ivarid = searchvar(ncid,vname)
    end if
    call mpi_gatherv(xval,gg%gc(myid+1),mpi_real8, &
                     ival,gg%gc,gg%gd,mpi_real8,   &
                     iocpu,gg%icomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
    if ( myid /= iocpu ) return
    ib = 1
    do i = iout1 , iout2
      do j = jout1 , jout2
        if ( gg%gcmask(j,i) ) then
          ncid%i4buf(j,i) = ival(ib)
          ib = ib + 1
        end if
      end do
    end do
    deallocate(ival)
    istart(3) = nt
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = niout
    icount(1) = njout
    incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf,istart(1:3),icount(1:3))
    call clm_checkncerr(__FILE__,__LINE__, &
      'Error read '//vname//' to file '//trim(ncid%fname))
  end subroutine clm_writerec_real8_2d_par_gg

  subroutine clm_writerec_real8_3d_par_gg(ncid,vname,xval,gg,nt)
    implicit none
    type(clm_filetype) , intent(inout) :: ncid
    character(len=*) , intent(in) :: vname
    real(rk8) , dimension(:,:) , intent(in) :: xval
    type(processor_type) , intent(in) :: gg
    integer(ik4) , intent(in) :: nt
    real(rk8) , dimension(:) , allocatable :: ival
    integer(ik4) :: i , j , k , nk , ib , ivarid , mpierr
    nk = size(xval,2)
    if ( myid == iocpu ) then
      allocate(ival(numg))
      ivarid = searchvar(ncid,vname)
    end if
    do k = 1 , nk
      call mpi_gatherv(xval(:,k),gg%gc(myid+1),mpi_real8, &
                       ival,gg%gc,gg%gd,mpi_real8,   &
                       iocpu,gg%icomm,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
      if ( myid == iocpu ) then
        ib = 1
        do i = iout1 , iout2
          do j = jout1 , jout2
            if ( gg%gcmask(j,i) ) then
              ncid%i4buf(j,i) = ival(ib)
              ib = ib + 1
            end if
          end do
        end do
        istart(4) = nt
        istart(3) = k
        istart(2) = 1
        istart(1) = 1
        icount(4) = 1
        icount(3) = 1
        icount(2) = niout
        icount(1) = njout
        incstat = nf90_put_var(ncid%ncid,ivarid,ncid%i4buf, &
                               istart(1:4),icount(1:4))
        call clm_checkncerr(__FILE__,__LINE__, &
          'Error read '//vname//' to file '//trim(ncid%fname))
      end if
    end do
    if ( myid /= iocpu ) return
    deallocate(ival)
  end subroutine clm_writerec_real8_3d_par_gg

  subroutine test_clmhelper
    implicit none
    type(clm_filetype) :: ncid
    logical , pointer , dimension(:) :: xval
    integer(ik4) , pointer , dimension(:) :: ival
    call clm_createfile('peppe.nc',ncid)
    call clm_adddim(ncid,'pippo',gcomm_pft%ns)
    call clm_adddim(ncid,'jx',njout)
    call clm_adddim(ncid,'iy',niout)
    call clm_addvar(clmvar_logical,ncid,'test',(/'pippo'/), &
            long_name='pillo',units='pallo')
    call clm_addvar(clmvar_integer,ncid,'mask',(/'jx','iy'/), &
            long_name='mask',units='1')
    call clm_enddef(ncid)
    allocate(xval(gcomm_pft%is:gcomm_pft%ie))
    allocate(ival(procinfo%ncells))
    ival = 1
    xval = .false.
    if ( myid == 1 ) xval = .true.
    call clm_writevar_par(ncid,'test',xval,gcomm_pft)
    call clm_writevar_par(ncid,'mask',ival,procinfo)
    deallocate(xval)
    deallocate(ival)
    call clm_closefile(ncid)
  end subroutine test_clmhelper

end module mod_clm_nchelper
