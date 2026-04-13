
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

#ifdef OASIS

module mod_oasis_regcm_field

#ifndef TESTME
  use mod_oasis
  use mod_oasis_params
  use mod_memutil
#endif
  use mod_stdio
  use mod_intkinds
  use mod_realkinds
  use mod_oasis_regcm_grid

  implicit none

  private

  type :: regcm_oasis_field
    integer :: id
    character(len=16) :: fname
    logical :: factive
    type(regcm_oasis_grid), pointer :: grd => null( )
    class(regcm_oasis_field), pointer :: next => null( )
  end type regcm_oasis_field

  type, extends(regcm_oasis_field) :: regcm_oasis_field_2d
    real(rk8), pointer, contiguous, dimension(:,:) :: tp
  end type regcm_oasis_field_2d

  type, extends(regcm_oasis_field) :: regcm_oasis_field_3d
    integer :: nz
    real(rk8), pointer, contiguous, dimension(:,:,:) :: tp
  end type regcm_oasis_field_3d

  interface associate_storage
    module procedure associate_storage_field_2d
    module procedure associate_storage_field_3d
  end interface associate_storage

  interface get_associated_storage
    module procedure get_associated_storage_field_2d
    module procedure get_associated_storage_field_3d
  end interface get_associated_storage

  public :: regcm_oasis_field
  public :: regcm_oasis_field_2d
  public :: regcm_oasis_field_3d
  public :: populate_field_list
  public :: depopulate_field_list
  public :: add_field, lookup_field
  public :: print_field_list
  public :: field_is_active, activate_field
  public :: associate_storage, get_associated_storage
  public :: associate_grid_to_field
  public :: associate_grid_to_field_list
  public :: activate_regcm_export_flist
  public :: activate_regcm_import_flist

  integer, parameter :: oasis_2d_maxfields = 56
  integer, parameter :: oasis_3d_maxfields = 7

  type(regcm_oasis_field_2d), dimension(oasis_2d_maxfields) :: exch_2d = &
    [ regcm_oasis_field_2d(0,'RCM_SP',   .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_SST',  .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_WZ0',  .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_WUST', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_U10M', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_V10M', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_WSPD', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_SIT',  .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_WDIR', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_T2M',  .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_Q2M',  .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_SLP',  .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_TAUX', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_TAUY', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_Z0',   .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_USTR', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_EVAP', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_PREC', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_NUWA', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_ULHF', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_USHF', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_UWLW', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_DWLW', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_NULW', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_UWSW', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_DWSW', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_NDSW', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_RHOA', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_T10M', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_Q10M', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_TATM', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_UATM', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_VATM', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_QATM', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_ZATM', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_PATM', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_HGT',  .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_SWDIR',.false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_SWDIF',.false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_TGBB', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_RAM1', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_RAH1', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_SNCV', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_TLEF', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_SNOW', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_ALBEI',.false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_ALBED',.false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_FCH4', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_TLAI', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_TROF', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_SROF', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_SMELT',.false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_H2O10',.false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_ALBI2',.false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_ALBD2',.false.,null( ),null( ),null( )), &
      regcm_oasis_field_2d(0,'RCM_DEFLX',.false.,null( ),null( ),null( )) ]
  type(regcm_oasis_field_3d), dimension(oasis_3d_maxfields) :: exch_3d = &
    [ regcm_oasis_field_3d(0,'RCM_TSOI', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_3d(0,'RCM_H2OV', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_3d(0,'RCM_H2OL', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_3d(0,'RCM_H2OI', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_3d(0,'RCM_FVOC', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_3d(0,'RCM_FDST', .false.,null( ),null( ),null( )), &
      regcm_oasis_field_3d(0,'RCM_DDVEL',.false.,null( ),null( ),null( )) ]

#ifdef TESTME
  integer :: counter = 1
  integer, parameter :: oasis_in = 0
  integer, parameter :: oasis_out = 1
  character(len=*), parameter :: oasis_prefix = 'OASIS'
  interface getmem
    module procedure getmem2d
    module procedure getmem3d
  end interface getmem
#endif

  contains

  subroutine populate_field_list(flist)
    implicit none
    class(regcm_oasis_field), pointer, intent(inout) :: flist
    integer :: i
    call depopulate_field_list(flist)
    do i = 1, oasis_2d_maxfields
      call add_field(flist, exch_2d(i))
    end do
    do i = 1, oasis_3d_maxfields
      call add_field(flist, exch_3d(i))
    end do
  end subroutine populate_field_list

  subroutine get_associated_storage_field_2d(field,tp)
    implicit none
    class(regcm_oasis_field), pointer, intent(inout) :: field
    real(rk8), dimension(:,:), contiguous, intent(inout), pointer :: tp
    select type (field)
      type is ( regcm_oasis_field_2d )
        tp => field%tp
      class default
        continue
    end select
  end subroutine get_associated_storage_field_2d

  subroutine get_associated_storage_field_3d(field,tp)
    implicit none
    class(regcm_oasis_field), pointer, intent(inout) :: field
    real(rk8), dimension(:,:,:), contiguous, intent(inout), pointer :: tp
    select type (field)
      type is ( regcm_oasis_field_3d )
        tp => field%tp
      class default
        continue
    end select
  end subroutine get_associated_storage_field_3d

  subroutine associate_grid_to_field(field,grd)
    implicit none
    class(regcm_oasis_field), pointer, intent(inout) :: field
    type(regcm_oasis_grid) , pointer :: grd
    if ( associated(field%grd) ) then
      nullify(field%grd)
    end if
    field%grd => grd
  end subroutine associate_grid_to_field

  subroutine associate_grid_to_field_list(flist,grd)
    implicit none
    class(regcm_oasis_field), pointer, intent(inout) :: flist
    type(regcm_oasis_grid) , pointer :: grd
    class(regcm_oasis_field), pointer :: p
    p => flist
    do while(associated(p))
      p%grd => grd
      p => p%next
    end do
  end subroutine associate_grid_to_field_list

  subroutine associate_storage_field_2d(field,tp)
    implicit none
    class(regcm_oasis_field), pointer, intent(inout) :: field
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: tp
    select type (field)
      type is ( regcm_oasis_field_2d )
        field%tp => tp
      class default
        continue
    end select
  end subroutine associate_storage_field_2d

  subroutine associate_storage_field_3d(field,tp)
    implicit none
    class(regcm_oasis_field), pointer, intent(inout) :: field
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: tp
    select type (field)
      type is ( regcm_oasis_field_3d )
        field%tp => tp
        field%nz = size(tp,3)
      class default
        continue
    end select
  end subroutine associate_storage_field_3d

  subroutine activate_field(flist,fname)
    implicit none
    class(regcm_oasis_field), pointer, intent(inout) :: flist
    character(len=*), intent(in) :: fname
    class(regcm_oasis_field), pointer :: p
    logical :: lact
    p => lookup_field(flist,fname)
    if ( associated(p) ) then
      p%factive = .true.
    end if
  end subroutine activate_field

  function field_is_active(flist,fname) result(lact)
    implicit none
    class(regcm_oasis_field), pointer, intent(in) :: flist
    character(len=*), intent(in) :: fname
    class(regcm_oasis_field), pointer :: p
    logical :: lact
    lact = .false.
    p => lookup_field(flist,fname)
    if ( associated(p) ) then
      lact = p%factive
    end if
  end function field_is_active

  subroutine print_field_list(flist)
    implicit none
    class(regcm_oasis_field), pointer, intent(in) :: flist
    class(regcm_oasis_field), pointer :: p
    if ( .not. associated(flist) ) then
      write(6,*) 'Empty list'
      return
    end if
    write(6,*) 'FIELD_NAME     ACTIV GRID  STORAGE SIZE'
    write(6,*) '--------------------------------------------------------'
    p => flist
    do while (associated(p))
      select type (p)
        type is ( regcm_oasis_field_2d )
          write(6,'(1x,a,1l,4x,1l,1x)',advance='no') p%fname, p%factive, &
                        associated(p%grd)
          if ( associated(p%tp) ) then
            write(6,*) shape(p%tp)
          else
            write(6,*) ' Unassociated'
          end if
        type is ( regcm_oasis_field_3d )
          write(6,'(1x,a,1l,4x,1l,1x)',advance='no') p%fname, p%factive, &
                        associated(p%grd)
          if ( associated(p%tp) ) then
            write(6,*) shape(p%tp)
          else
            write(6,*) ' Unassociated'
          end if
        class default
          write(6,*) 'Unknown field type'
      end select
      p => p%next
    end do
  end subroutine print_field_list

  subroutine add_field(flist,f,nz)
    implicit none
    class(regcm_oasis_field), pointer, intent(inout) :: flist
    class(regcm_oasis_field), pointer, intent(in) :: f
    class(regcm_oasis_field), pointer :: p, tmp, fcopy
    type(regcm_oasis_field_2d), pointer :: fcopy2
    type(regcm_oasis_field_3d), pointer :: fcopy3
    integer, optional :: nz
    select type(f)
      class is ( regcm_oasis_field_2d )
        allocate(fcopy2)
        fcopy2%fname = f%fname
        fcopy2%factive = f%factive
        fcopy2%id = f%id
        fcopy2%grd => f%grd
        fcopy2%next => null( )
        fcopy2%tp => f%tp
        fcopy => fcopy2
      class is ( regcm_oasis_field_3d )
        allocate(fcopy3)
        fcopy3%fname = f%fname
        fcopy3%factive = f%factive
        fcopy3%id = f%id
        fcopy3%grd => f%grd
        fcopy3%next => null( )
        fcopy3%tp => f%tp
        if ( present(nz) ) then
          fcopy3%nz = nz
        else
          fcopy3%nz = 1
        end if
        fcopy => fcopy3
      class default
        return
    end select
    if ( .not. associated(flist) ) then
      flist => fcopy
      return
    end if
    p => flist
    tmp => flist
    do while(associated(p))
      ! name is key: check if already in list
      ! If in list, nothing to add ;)
      if ( p%fname == f%fname ) then
        return
      end if
      tmp => p
      p => p%next
    end do
    tmp%next => fcopy
  end subroutine add_field

  function lookup_field(flist,fname) result(p)
    implicit none
    class(regcm_oasis_field), pointer, intent(in) :: flist
    character(len=*), intent(in) :: fname
    class(regcm_oasis_field), pointer :: p
    p => flist
    do while(associated(p))
      if ( p%fname == fname ) then
        return
      end if
      p => p%next
    end do
    p => null( )
  end function lookup_field

  subroutine depopulate_field_list(flist)
    implicit none
    class(regcm_oasis_field), pointer, intent(inout) :: flist
    class(regcm_oasis_field), pointer :: p, tmp
    p => flist
    do while(associated(p))
      tmp => p
      p => p%next
      deallocate(tmp)
    end do 
    nullify(flist)
  end subroutine depopulate_field_list

  subroutine activate_regcm_export_flist(flist)
    implicit none
    class(regcm_oasis_field), pointer, intent(inout) :: flist
    class(regcm_oasis_field), pointer :: p
#ifdef DEBUG
    write(stdout,*) oasis_prefix, ' >> Activating field list'
    write(stdout,*) oasis_prefix, ' >> direction: export'
#endif
    p => flist
    do while(associated(p))
      call define_onefield(p,OASIS_Out)
      p => p%next
    end do 
#ifdef DEBUG
    write(stdout,*) oasis_prefix, ' >> Successful activation'
#endif
  end subroutine activate_regcm_export_flist

  subroutine activate_regcm_import_flist(flist)
    implicit none
    class(regcm_oasis_field), pointer, intent(inout) :: flist
    class(regcm_oasis_field), pointer :: p
#ifdef DEBUG
    write(stdout,*) oasis_prefix, ' >> Activating field list'
    write(stdout,*) oasis_prefix, ' >> direction: import'
#endif
    p => flist
    do while(associated(p))
      call define_onefield(p,OASIS_In)
      p => p%next
    end do 
  end subroutine activate_regcm_import_flist

  subroutine define_onefield(p, direction)
    implicit none
    class(regcm_oasis_field), pointer, intent(inout) :: p
    integer, intent(in) :: direction
    integer :: ierror
#ifdef DEBUG
    write(stdout,*) oasis_prefix, ' >>      name: ', p%fname
    write(stdout,*) oasis_prefix, ' >>      grid: ', &
            p%grd%naNM, '/', p%grd%naWM
#endif
#ifdef TESTME
    p%id = counter
    counter = counter + 1
    ierror = 0
    p%factive = .true.
#else
    call oasis_def_var(p%id, p%fname, p%grd%id, [1,1], &
            direction, [1,1], OASIS_Real, ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'oasis_def_var (', p%fname, ') abort compid ', comp_id
      call oasis_abort(comp_id, __FILE__, &
              'Problem in oasis_def_var call for field '//p%fname)
    end if
#endif
    select type(p)
      class is ( regcm_oasis_field_2d )
        if ( .not. associated(p%tp) ) then
          call getmem(p%tp,p%grd%j1l,p%grd%j2l,p%grd%i1l,p%grd%i2l, &
                  'OASIS: oasis_field_2d '//trim(p%fname))
        end if
      class is ( regcm_oasis_field_3d )
        if ( .not. associated(p%tp) ) then
          call getmem(p%tp,p%grd%j1l,p%grd%j2l,p%grd%i1l,p%grd%i2l, &
                  1, p%nz, 'OASIS: oasis_field_3d '//trim(p%fname))
        end if
      class default
        continue
    end select
#ifdef DEBUG
    write(stdout,"(' ',A,A,I3)") oasis_prefix, ' >>       id: ', p%id
#endif
  end subroutine define_onefield

#ifdef TESTME
  subroutine getmem2d(p,j1,j2,i1,i2,a)
    implicit none
    real(rk8), pointer, dimension(:,:), intent(inout) :: p
    integer, intent(in) :: j1,j2,i1,i2
    character(len=*), intent(in) :: a
    allocate(p(j1:j2,i1:i2))
    print *, 'ALLOCATE ', a, shape(p)
  end subroutine getmem2d
  subroutine getmem3d(p,j1,j2,i1,i2,k1,k2,a)
    implicit none
    real(rk8), pointer, dimension(:,:,:), intent(inout) :: p
    integer, intent(in) :: j1,j2,i1,i2,k1,k2
    character(len=*), intent(in) :: a
    allocate(p(j1:j2,i1:i2,k1:k2))
    print *, 'ALLOCATE ', a, shape(p)
  end subroutine getmem3d
#endif

end module mod_oasis_regcm_field

#ifdef TESTME
program testme
  use mod_realkinds
  use mod_oasis_regcm_grid
  use mod_oasis_regcm_field

  class(regcm_oasis_field), pointer :: flist => null( )
  class(regcm_oasis_field), pointer :: flist_regcm_export => null( )
  class(regcm_oasis_field), pointer :: flist_regcm_import => null( )
  class(regcm_oasis_field), pointer :: p
  real(rk8), dimension(:,:), contiguous, pointer :: f2d, f2d1
  real(rk8), dimension(:,:,:), contiguous, pointer :: f3d
  real(rk8), dimension(:,:,:), contiguous, pointer :: f3d_ref
  type(regcm_oasis_grid) , pointer :: pg

  allocate(f2d(2:10,2:10))
  allocate(f2d1(2:10,2:10))
  allocate(f3d(2:10,2:10,5), source=10.0_rk8)

  call populate_field_list(flist)

  ! jci1, jci2, 2, jx-2, ici1, ici2, 2, iy-2
  pg => new_oasis_regcm_grid('rcin','rcim', 2, 10, 2, 19, 2, 10, 2, 19)

  call add_field(flist_regcm_import, lookup_field(flist,'RCM_SST'))
  call add_field(flist_regcm_import, lookup_field(flist,'RCM_SIT'))
  call add_field(flist_regcm_import, lookup_field(flist,'RCM_TSOI'),5)

  call associate_grid_to_field_list(flist_regcm_import, pg)

  call associate_storage(lookup_field(flist_regcm_import,'RCM_SST'), f2d)
  call associate_storage(lookup_field(flist_regcm_import,'RCM_SIT'), f2d1)
  call associate_storage(lookup_field(flist_regcm_import,'RCM_TSOI'), f3d)

  call activate_regcm_import_flist(flist_regcm_import)

  call add_field(flist_regcm_export, lookup_field(flist,'RCM_U10M'))
  call add_field(flist_regcm_export, lookup_field(flist,'RCM_V10M'))
  call add_field(flist_regcm_export, lookup_field(flist,'RCM_PATM'))

  call associate_grid_to_field_list(flist_regcm_export, pg)
  call activate_regcm_export_flist(flist_regcm_export)

  call print_field_list(flist_regcm_import)
  call print_field_list(flist_regcm_export)

  if ( .not. field_is_active(flist_regcm_import,'RCM_SST') ) then
    print *, 'Error activating'
  end if

  p => lookup_field(flist_regcm_import,'RCM_SST')
  if ( associated(p) ) then
    print *, p%fname, p%factive
  end if
  p => lookup_field(flist_regcm_import,'RCM_T2M')
  if ( associated(p) ) then
    print *, p%fname, p%factive
  end if
  p => lookup_field(flist,'RCM_T2M')
  if ( associated(p) ) then
    print *, p%fname, p%factive
  end if
  call get_associated_storage( &
          lookup_field(flist_regcm_import,'RCM_TSOI'), f3d_ref)
  print *, maxval(f3d_ref)

  ! Cleanup storage
  call depopulate_field_list(flist)
  call depopulate_field_list(flist_regcm_export)
  call depopulate_field_list(flist_regcm_import)
  deallocate(f2d, f2d1, f3d)
  deallocate(pg)
end program testme
#endif

#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
