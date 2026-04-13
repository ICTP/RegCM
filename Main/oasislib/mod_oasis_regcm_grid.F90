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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               REMINDER FOR REGCM GRIDS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    jci1, jci2, 2, jx-2, ici1, ici2, 2, iy-2  ! Cross, no border
!    jce1, jce2, 1, jx-1, ice1, ice2, 1, iy-1  ! Cross, border
!    jdi1, jdi2, 2, jx-1, idi1, idi2, 2, iy-1  ! Dot, no border
!    jde1, jde2, 1, jx  , ide1, ide2, 1, iy    ! Dot, border
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mod_oasis_regcm_grid

#ifndef TESTME
  use mod_oasis
  use mod_oasis_params
#endif
  use mod_intkinds
  use mod_realkinds
  use mod_stdio

  implicit none

  private

  ! 3d fields are transered in 2d blocks per layer
  type :: regcm_oasis_grid
    integer :: id             ! partition identificator
    character(len=4) :: naNM  ! grid name (No Mask)
    character(len=4) :: naWM  ! grid name (With Mask)
    integer :: j1l, j2l       ! local indexes dim1 (lon)
    integer :: j1g, j2g       ! global indexes dim1 (lon)
    integer :: i1l, i2l       ! local indexes dim2 (lat)
    integer :: i1g, i2g       ! global indexes dim2 (lat)
    integer :: njl, nil       ! Number of points per each dim (local)
    integer :: njg, nig       ! Number of points per each dim (global)
    integer :: nc             ! number of corner per cell (4)
  end type regcm_oasis_grid

  interface new_oasis_regcm_grid
    module procedure new_grid_from_grid
    module procedure new_grid_from_args
  end interface new_oasis_regcm_grid

  public :: regcm_oasis_grid
  public :: new_oasis_regcm_grid

#ifdef TESTME
  integer :: counter = 1
  character(len=*), parameter :: oasis_prefix = 'OASIS'
#endif

  contains

  function check_grid(grd) result(lok)
    implicit none
    type(regcm_oasis_grid), pointer, intent(in) :: grd
    logical :: lok
    lok = .true.
    if ( grd%j1g <= 0 .or. grd%j2g <= 0 ) then
      lok = .false.
    end if
    if ( grd%j1l < grd%j1g .or. grd%j2l > grd%j2g ) then
      lok = .false.
    end if
    if ( grd%i1g <= 0 .or. grd%i2g <= 0 ) then
      lok = .false.
    end if
    if ( grd%i1l < grd%i1g .or. grd%i2l > grd%i2g ) then
      lok = .false.
    end if
  end function check_grid

  function new_grid_from_args(naNM,naWM, &
                  j1l,j2l,j1g,j2g,i1l,i2l,i1g,i2g) result(fg)
    implicit none
    character(len=4), intent(in) :: naNM, naWM
    integer, intent(in) :: j1l, j2l, j1g, j2g, i1l, i2l, i1g, i2g
    type(regcm_oasis_grid), pointer :: fg
    integer :: ioffset, ierror
    integer, dimension(5) :: ig_paral
    allocate(fg)
    fg%naNM = naNM
    fg%naWM = naWM
    fg%j1l = j1l
    fg%j2l = j2l
    fg%j1g = j1g
    fg%j2g = j2g
    fg%njl = j2l-j1l+1
    fg%njg = j2g-j1g+1
    fg%i1l = i1l
    fg%i2l = i2l
    fg%i1g = i1g
    fg%i2g = i2g
    fg%nil = i2l-i1l+1
    fg%nig = i2g-i1g+1
    fg%nc = 4
    if ( .not. check_grid(fg) ) then
      deallocate(fg)
      nullify(fg)
      return
    end if
#ifdef DEBUG
    write(stdout,*) oasis_prefix, ' >> Define ', fg%naNM, '/', &
                   fg%naWM, ' partitions'
#endif
#ifdef TESTME
    fg%id = counter
    counter = counter + 1
#else
    ! 2 indicates a box partition
    ioffset = (i1g-i1l)*fg%njg + (j1g-j1l)
    ig_paral = [ 2, ioffset, j1l, i1l, fg%njg ]
    call oasis_def_partition(fg%id, ig_paral, ierror)
    if (ierror /= 0) then
      write(stderr, *) 'oasis_def_partition (', naNM, '/', &
                naWM, ') abort compid ', comp_id
      call oasis_abort(comp_id, __FILE__, &
                'Problem in oasis_def_partition call for '// &
                naNM//'/'//naWM)
    end if
#endif
#ifdef DEBUG
    write(stdout,"(' ',A,A,I3)") oasis_prefix, &
            ' >> Box partition successfully defined -- id: ', fg%id
#endif
  end function new_grid_from_args

  function new_grid_from_grid(grd) result(fg)
    implicit none
    type(regcm_oasis_grid) :: grd
    type(regcm_oasis_grid), pointer :: fg
    integer :: ioffset, ierror
    integer, dimension(5) :: ig_paral
    allocate(fg)
    fg%naNM = grd%naNM
    fg%naWM = grd%naWM
    fg%j1l = grd%j1l
    fg%j2l = grd%j2l
    fg%j1g = grd%j1g
    fg%j2g = grd%j2g
    fg%njl = grd%j2l-grd%j1l+1
    fg%njg = grd%j2g-grd%j1g+1
    fg%i1l = grd%i1l
    fg%i2l = grd%i2l
    fg%i1g = grd%i1g
    fg%i2g = grd%i2g
    fg%nil = grd%i2l-grd%i1l+1
    fg%nig = grd%i2g-grd%i1g+1
    fg%nc = grd%nc
    if ( .not. check_grid(fg) ) then
      deallocate(fg)
      nullify(fg)
      return
    end if
#ifdef DEBUG
    write(stdout,*) oasis_prefix, ' >> Define ', fg%naNM, '/', &
                   fg%naWM, ' partitions'
#endif
#ifdef TESTME
    fg%id = counter
    counter = counter + 1
#else
    ! 2 indicates a box partition
    ioffset = (fg%i1g-fg%i1l)*fg%njg + (fg%j1g-fg%j1l)
    ig_paral = [ 2, ioffset, fg%j1l, fg%i1l, fg%njg ]
    call oasis_def_partition(fg%id, ig_paral, ierror)
    if (ierror /= 0) then
      write(stderr, *) 'oasis_def_partition (', fg%naNM, '/', &
                fg%naWM, ') abort compid ', comp_id
      call oasis_abort(comp_id, __FILE__, &
                'Problem in oasis_def_partition call for '// &
                fg%naNM//'/'//fg%naWM)
    end if
#endif
#ifdef DEBUG
    write(stdout,"(' ',A,A,I3)") oasis_prefix, &
            ' >> Box partition successfully defined -- id: ', fg%id
#endif
  end function new_grid_from_grid

end module mod_oasis_regcm_grid

#ifdef TESTMEPROG
program testme
  use mod_realkinds
  use mod_oasis_regcm_grid

  type(regcm_oasis_grid) , pointer :: pgde, pgdi
  type(regcm_oasis_grid) , pointer :: pgce, pgci

  !    jci1, jci2, 2, jx-2, ici1, ici2, 2, iy-2  ! Cross, no border (Internal)
  !    jce1, jce2, 1, jx-1, ice1, ice2, 1, iy-1  ! Cross, border    (External)
  !    jdi1, jdi2, 2, jx-1, idi1, idi2, 2, iy-1  ! Dot, no border   (Internal)
  !    jde1, jde2, 1, jx  , ide1, ide2, 1, iy    ! Dot, border      (External)

  ! Pointer to the Regcm grid dot  , external, nomask and mask
  pgde => new_oasis_regcm_grid('rden', 'rdem', 10, 21, 1, 21, 10, 21, 1, 21)
  ! Pointer to the Regcm grid dot  , internal, nomask and mask
  pgdi => new_oasis_regcm_grid('rdin', 'rdim', 2, 10, 2, 20, 2, 10, 2, 20)
  ! Pointer to the Regcm grid cross, external, nomask and mask
  pgce => new_oasis_regcm_grid('rcen', 'rcem', 10, 20, 1, 20, 10, 20, 1, 20)
  ! Pointer to the Regcm grid cross, internal, nomask and mask
  pgci => new_oasis_regcm_grid('rcin', 'rcim', 2, 10, 2, 19, 2, 10, 2, 19)
  
  deallocate(pgde)
  deallocate(pgdi)
  deallocate(pgce)
  deallocate(pgci)

end program testme
#endif

#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
