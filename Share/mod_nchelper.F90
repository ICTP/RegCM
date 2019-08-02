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

module mod_nchelper

  use netcdf
  use mod_realkinds
  use mod_stdio
  use mod_constants
  use mod_memutil
  use mod_dynparam
  use mod_message

  implicit none

  private

  public :: openfile_withname
  public :: createfile_withname
  public :: closefile
  public :: check_dims
  public :: check_var
  public :: ncd_inqdim
  public :: add_dimension
  public :: add_variable
  public :: add_attribute
  public :: get_attribute
  public :: check_dimlen
  public :: checkncerr

#ifdef SINGLE_PRECISION_REAL
  integer(ik4) , public :: regcm_vartype = nf90_real
#else
  integer(ik4) , public :: regcm_vartype = nf90_double
#endif

  interface read_var1d_static
    module procedure read_var1d_static_double_fix
    module procedure read_var1d_static_single_fix
    module procedure read_var1d_static_integer_fix
    module procedure read_var1d_static_double
    module procedure read_var1d_static_single
    module procedure read_var1d_static_integer
    module procedure read_var1d_static_text
  end interface read_var1d_static

  interface read_var2d_static
    module procedure read_var2d_static_double_fix
    module procedure read_var2d_static_single_fix
    module procedure read_var2d_static_integer_fix
    module procedure read_var2d_static_double
    module procedure read_var2d_static_single
    module procedure read_var2d_static_integer
  end interface read_var2d_static

  interface read_var3d_static
    module procedure read_var3d_static_double_fix
    module procedure read_var3d_static_single_fix
    module procedure read_var3d_static_integer_fix
    module procedure read_var3d_static_double
    module procedure read_var3d_static_single
    module procedure read_var3d_static_integer
  end interface read_var3d_static

  interface write_var1d_static
    module procedure write_var1d_static_double
    module procedure write_var1d_static_single
    module procedure write_var1d_static_integer
    module procedure write_var1d_static_text
  end interface write_var1d_static

  interface get_attribute
    module procedure get_attribute_char
    module procedure get_attribute_int4
    module procedure get_attribute_real4
    module procedure get_attribute_real8
  end interface get_attribute

  public :: read_var1d_static
  public :: read_var2d_static
  public :: read_var3d_static

  public :: write_var1d_static
  public :: write_var2d_static
  public :: write_var3d_static

  integer(ik4) :: incstat

  contains

  subroutine write_var1d_static_single(ncid,vnam,values,ipnt,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) :: vnam
    real(rk4) , dimension(:) , intent(in) :: values
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) , intent(in) , dimension(:) :: ivar
    incstat = nf90_put_var(ncid, ivar(ipnt), values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable '//vnam//' write')
    ipnt = ipnt + 1
    if ( debug_level > 2 ) then
      incstat = nf90_sync(ncid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error variable '//vnam//' sync')
    end if
  end subroutine write_var1d_static_single

  subroutine write_var1d_static_double(ncid,vnam,values,ipnt,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) :: vnam
    real(rk8) , dimension(:) , intent(in) :: values
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) , intent(in) , dimension(:) :: ivar
    incstat = nf90_put_var(ncid, ivar(ipnt), values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable '//vnam//' write')
    ipnt = ipnt + 1
    if ( debug_level > 2 ) then
      incstat = nf90_sync(ncid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error variable '//vnam//' sync')
    end if
  end subroutine write_var1d_static_double

  subroutine write_var1d_static_integer(ncid,vnam,values,ipnt,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) :: vnam
    integer(ik4) , dimension(:) , intent(in) :: values
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) , intent(in) , dimension(:) :: ivar
    incstat = nf90_put_var(ncid, ivar(ipnt), values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable '//vnam//' write')
    ipnt = ipnt + 1
    if ( debug_level > 2 ) then
      incstat = nf90_sync(ncid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error variable '//vnam//' sync')
    end if
  end subroutine write_var1d_static_integer

  subroutine write_var1d_static_text(ncid,vnam,values,ipnt,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) :: vnam
    character(len=*) , intent(in) :: values
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) , intent(in) , dimension(:) :: ivar
    incstat = nf90_put_var(ncid, ivar(ipnt), values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable '//vnam//' write')
    ipnt = ipnt + 1
    if ( debug_level > 2 ) then
      incstat = nf90_sync(ncid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error variable '//vnam//' sync')
    end if
  end subroutine write_var1d_static_text

  subroutine write_var2d_static(ncid,vnam,values,ipnt,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) :: vnam
    real(rk4) , dimension(:,:) , intent(in) :: values
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) , intent(in) , dimension(:) :: ivar
    incstat = nf90_put_var(ncid, ivar(ipnt), values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable '//vnam//' write')
    ipnt = ipnt + 1
    if ( debug_level > 2 ) then
      incstat = nf90_sync(ncid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error variable '//vnam//' sync')
    end if
  end subroutine write_var2d_static

  subroutine write_var3d_static(ncid,vnam,values,ipnt,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    real(rk4) , dimension(:,:,:) , intent(in) :: values
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) , intent(in) , dimension(:) :: ivar
    integer(ik4) , dimension(3) :: istart
    integer(ik4) , dimension(3) :: icount
    istart(1) = 1
    istart(2) = 1
    icount(1) = ubound(values,1)
    icount(2) = ubound(values,2)
    incstat = nf90_put_var(ncid,ivar(ipnt),values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable '//vnam//' write')
    if ( debug_level > 2 ) then
      incstat = nf90_sync(ncid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error variable '//vnam//' sync')
    end if
    ipnt = ipnt + 1
  end subroutine write_var3d_static

  subroutine read_var1d_static_text(ncid,vnam,values)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    character(len=*) , dimension(:) :: values
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid, vnam, ivarid)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error search '//vnam)
    incstat = nf90_get_var(ncid, ivarid, values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error read '//vnam)
  end subroutine read_var1d_static_text

  subroutine read_var1d_static_single(ncid,vnam,values,lerror)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    real(rk4) , pointer , dimension(:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid, vnam, ivarid)
    if ( incstat /= nf90_noerr ) then
      if ( present(lerror) ) then
        if ( lerror ) then
          lerror = .false.
          return
        end if
      end if
    end if
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error search '//vnam)
    incstat = nf90_get_var(ncid, ivarid, values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error read '//vnam)
  end subroutine read_var1d_static_single

  subroutine read_var1d_static_double(ncid,vnam,values,lerror)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    real(rk8) , pointer , dimension(:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid, vnam, ivarid)
    if ( incstat /= nf90_noerr ) then
      if ( present(lerror) ) then
        if ( lerror ) then
          lerror = .false.
          return
        end if
      end if
    end if
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error search '//vnam)
    incstat = nf90_get_var(ncid, ivarid, values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error read '//vnam)
  end subroutine read_var1d_static_double

  subroutine read_var1d_static_integer(ncid,vnam,values,lerror)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , pointer , dimension(:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid, vnam, ivarid)
    if ( incstat /= nf90_noerr ) then
      if ( present(lerror) ) then
        if ( lerror ) then
          lerror = .false.
          return
        end if
      end if
    end if
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error search '//vnam)
    incstat = nf90_get_var(ncid, ivarid, values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error read '//vnam)
  end subroutine read_var1d_static_integer

  subroutine read_var1d_static_double_fix(ncid,vnam,n,values,lerror)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , intent(in) :: n
    real(rk8) , dimension(n) :: values
    logical , optional , intent(inout) :: lerror
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid, vnam, ivarid)
    if ( incstat /= nf90_noerr ) then
      if ( present(lerror) ) then
        if ( lerror ) then
          lerror = .false.
          return
        end if
      end if
    end if
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error search '//vnam)
    incstat = nf90_get_var(ncid, ivarid, values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error read '//vnam)
  end subroutine read_var1d_static_double_fix

  subroutine read_var1d_static_single_fix(ncid,vnam,n,values,lerror)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , intent(in) :: n
    real(rk4) , dimension(n) :: values
    logical , optional , intent(inout) :: lerror
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid, vnam, ivarid)
    if ( incstat /= nf90_noerr ) then
      if ( present(lerror) ) then
        if ( lerror ) then
          lerror = .false.
          return
        end if
      end if
    end if
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error search '//vnam)
    incstat = nf90_get_var(ncid, ivarid, values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error read '//vnam)
  end subroutine read_var1d_static_single_fix

  subroutine read_var1d_static_integer_fix(ncid,vnam,n,values,lerror)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , intent(in) :: n
    integer(ik4) , dimension(n) :: values
    logical , optional , intent(inout) :: lerror
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid, vnam, ivarid)
    if ( incstat /= nf90_noerr ) then
      if ( present(lerror) ) then
        if ( lerror ) then
          lerror = .false.
          return
        end if
      end if
    end if
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error search '//vnam)
    incstat = nf90_get_var(ncid, ivarid, values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error read '//vnam)
  end subroutine read_var1d_static_integer_fix

  subroutine read_var2d_static_double(ncid,vnam,values,lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    real(rk8) , dimension(:,:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var2d_static_double

  subroutine read_var2d_static_single(ncid,vnam,values,lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    real(rk4) , dimension(:,:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var2d_static_single

  subroutine read_var2d_static_integer(ncid,vnam,values,lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , dimension(:,:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var2d_static_integer

  subroutine read_var2d_static_double_fix(ncid,vnam,n,m,values, &
                                          lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , intent(in) :: n , m
    real(rk8) , dimension(n,m) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var2d_static_double_fix

  subroutine read_var2d_static_single_fix(ncid,vnam,n,m,values, &
                                          lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , intent(in) :: n , m
    real(rk4) , dimension(n,m) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var2d_static_single_fix

  subroutine read_var2d_static_integer_fix(ncid,vnam,n,m,values, &
                                          lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , intent(in) :: n , m
    integer(ik4) , dimension(n,m) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var2d_static_integer_fix

  subroutine read_var3d_static_double(ncid,vnam,values,lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    real(rk8) , dimension(:,:,:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var3d_static_double

  subroutine read_var3d_static_single(ncid,vnam,values,lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    real(rk4) , dimension(:,:,:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var3d_static_single

  subroutine read_var3d_static_integer(ncid,vnam,values,lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , dimension(:,:,:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var3d_static_integer

  subroutine read_var3d_static_double_fix(ncid,vnam,n,m,l,values, &
                                          lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer , intent(in) :: n , m , l
    real(rk8) , dimension(n,m,l) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var3d_static_double_fix

  subroutine read_var3d_static_single_fix(ncid,vnam,n,m,l,values, &
                                          lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer , intent(in) :: n , m , l
    real(rk4) , dimension(n,m,l) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var3d_static_single_fix

  subroutine read_var3d_static_integer_fix(ncid,vnam,n,m,l,values, &
                                           lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer , intent(in) :: n , m , l
    integer(ik4) , dimension(n,m,l) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var3d_static_integer_fix

  subroutine add_dimension(ncid,dnam,nd,ipnt,idims)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: dnam
    integer(ik4) , intent(in) :: nd
    integer(ik4) , intent(inout) , dimension(:) :: idims
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) :: incstat
    incstat = nd
    if ( nd == -1 ) incstat = nf90_unlimited
    incstat = nf90_def_dim(ncid, dnam, incstat, idims(ipnt))
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding dimension '//dnam)
    ipnt = ipnt + 1
  end subroutine add_dimension

  subroutine createfile_withname(fname,ncid)
    implicit none
    character(len=*) , intent(in) :: fname
    integer(ik4) , intent(out) :: ncid
    incstat = nf90_create(fname, iomode, ncid)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error creating NetCDF output '//trim(fname))
  end subroutine createfile_withname

  subroutine openfile_withname(fname,ncid)
    implicit none
    character(len=*) , intent(in) :: fname
    integer(ik4) , intent(out) :: ncid
    incstat = nf90_open(fname, nf90_nowrite, ncid)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error open NetCDF input '//trim(fname))
  end subroutine openfile_withname

  subroutine ncd_inqdim(ncid,dname,dlen,lerror,lexist)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: dname
    integer(ik4) , intent(out) , optional :: dlen
    logical , intent(in) , optional :: lerror
    logical , intent(out) , optional :: lexist
    integer(ik4) :: istatus
    integer(ik4) :: idimid
    istatus = nf90_inq_dimid(ncid, dname, idimid)
    if ( istatus /= nf90_noerr ) then
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
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error search dimension '//dname)
    if ( present(dlen) ) then
      istatus = nf90_inquire_dimension(ncid, idimid, len=dlen)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read dimension '//dname)
    end if
  end subroutine ncd_inqdim

  subroutine add_attribute(ncid,aname,aval,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: aname
    character(len=*) , intent(in) :: aval
    integer(ik4) , intent(in) , optional :: ivar
    integer :: istat
    if ( present(ivar) ) then
      istat = nf90_put_att(ncid,ivar,aname,aval)
    else
      istat = nf90_put_att(ncid,nf90_global,aname,aval)
    end if
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error adding attribute '//aname)
  end subroutine add_attribute

  subroutine get_attribute_char(ncid,aname,aval,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: aname
    character(len=*) , intent(out) :: aval
    integer(ik4) , intent(in) , optional :: ivar
    integer :: istat
    if ( present(ivar) ) then
      istat = nf90_get_att(ncid,ivar,aname,aval)
    else
      istat = nf90_get_att(ncid,nf90_global,aname,aval)
    end if
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error reading attribute '//aname)
  end subroutine get_attribute_char

  subroutine get_attribute_real4(ncid,aname,aval,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: aname
    real(rk4) , intent(out) :: aval
    integer(ik4) , intent(in) , optional :: ivar
    integer :: istat
    if ( present(ivar) ) then
      istat = nf90_get_att(ncid,ivar,aname,aval)
    else
      istat = nf90_get_att(ncid,nf90_global,aname,aval)
    end if
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error reading attribute '//aname)
  end subroutine get_attribute_real4

  subroutine get_attribute_real8(ncid,aname,aval,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: aname
    real(rk8) , intent(out) :: aval
    integer(ik4) , intent(in) , optional :: ivar
    integer :: istat
    if ( present(ivar) ) then
      istat = nf90_get_att(ncid,ivar,aname,aval)
    else
      istat = nf90_get_att(ncid,nf90_global,aname,aval)
    end if
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error reading attribute '//aname)
  end subroutine get_attribute_real8

  subroutine get_attribute_int4(ncid,aname,aval,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: aname
    integer(ik4) , intent(out) :: aval
    integer(ik4) , intent(in) , optional :: ivar
    integer :: istat
    if ( present(ivar) ) then
      istat = nf90_get_att(ncid,ivar,aname,aval)
    else
      istat = nf90_get_att(ncid,nf90_global,aname,aval)
    end if
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error reading attribute '//aname)
  end subroutine get_attribute_int4

  logical function check_dimlen(ncid,dname,ival)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: dname
    integer(ik4) , intent(in) :: ival
    integer(ik4) :: idimid , dlen , istatus
    check_dimlen = .false.
    istatus = nf90_inq_dimid(ncid, dname, idimid)
    if ( istatus /= nf90_noerr ) return
    istatus = nf90_inquire_dimension(ncid, idimid, len=dlen)
    if ( istatus /= nf90_noerr ) return
    if ( dlen /= ival ) return
    check_dimlen = .true.
  end function check_dimlen

  subroutine check_dims(ncid)
    implicit none
    integer(ik4) , intent(in) :: ncid
    integer(ik4) :: istatus
    integer(ik4) :: idimid
    integer(ik4) :: iyy , jxx , kzz
    istatus = nf90_inq_dimid(ncid, 'jx', idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error search dimension JX')
    istatus = nf90_inquire_dimension(ncid, idimid, len=jxx)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read dimension JX')
    if ( jx /= jxx ) then
      write(stderr,*) 'DOMAIN FILE : ', jxx
      write(stderr,*) 'NAMELIST    : ', jx
      call die('Mismatch: JX in DOMAIN file /= JX in namelist')
    end if
    istatus = nf90_inq_dimid(ncid, 'iy', idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error search dimension IY')
    istatus = nf90_inquire_dimension(ncid, idimid, len=iyy)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read dimension IY')
    if ( iy /= iyy ) then
      write(stderr,*) 'DOMAIN FILE : ', iyy
      write(stderr,*) 'NAMELIST    : ', iy
      call die('Mismatch: IY in DOMAIN file /= IY in namelist')
    end if
    istatus = nf90_inq_dimid(ncid, 'kz', idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error search dimension KZ')
    istatus = nf90_inquire_dimension(ncid, idimid, len=kzz)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read dimension KZ')
    if ( kz /= kzz ) then
      write(stderr,*) 'DOMAIN FILE : ', kzz
      write(stderr,*) 'NAMELIST    : ', kz
      call die('Mismatch: KZ in DOMAIN file /= KZ in namelist')
    end if
  end subroutine check_dims

  subroutine add_variable(ncid,varname,long_name,units,idims,ipnt,ivars)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: varname
    character(len=*) , intent(in) :: long_name
    character(len=*) , intent(in) :: units
    integer(ik4) , dimension(:) , intent(in) :: idims
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) , dimension(:) , intent(inout) :: ivars
    integer(ik4) :: incstat
    incstat = nf90_def_var(ncid, varname, regcm_vartype, idims, ivars(ipnt))
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding variable '//varname)
#ifdef NETCDF4_HDF5
#if defined (NETCDF4_COMPRESS)
    incstat = nf90_def_var_deflate(ncid, ivars(ipnt), 1, 1, deflate_level)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error setting deflate on xlat')
#endif
#endif
    incstat = nf90_put_att(ncid, ivars(ipnt), 'long_name',long_name)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error long_name to '//varname)
    incstat = nf90_put_att(ncid, ivars(ipnt), 'units', units)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error units to '//varname)
    ipnt = ipnt + 1
  end subroutine add_variable

  subroutine check_var(ncid,vname,lerror)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , intent(inout) , optional :: lerror
    integer(ik4) :: istatus
    integer(ik4) :: ivarid
    istatus = nf90_inq_varid(ncid,vname,ivarid)
    if ( istatus /= nf90_noerr ) then
      if ( present(lerror) ) then
        if ( lerror ) then
          lerror = .false.
          return
        end if
      end if
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error search '//vname)
  end subroutine check_var

  subroutine closefile(ncid)
    implicit none
    integer(ik4) , intent(in) :: ncid
    integer(ik4) :: istatus
    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error closing file')
  end subroutine closefile

  subroutine checkncerr(ival,filename,line,arg)
    implicit none
    integer(ik4) , intent(in) :: ival , line
    character(len=8) :: cline
    character(*) , intent(in) :: filename , arg
    if ( ival /= nf90_noerr ) then
      write (cline,'(i8)') line
      write (stderr,*) nf90_strerror(ival)
      call die(filename,trim(cline)//':'//arg,ival)
    end if
  end subroutine checkncerr

end module mod_nchelper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
