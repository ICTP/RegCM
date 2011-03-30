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
 
module mod_memutil

  use mod_constants
  use mod_dynparam
  use m_realkinds
  use m_die
  use m_stdio
  use m_mall

  private

  public :: memory_init , memory_destroy
  public :: getmem1d , relmem1d
  public :: getmem2d , relmem2d
  public :: getmem3d , relmem3d
  public :: getmem4d , relmem4d
  public :: getmem5d , relmem5d

  type pool1d
    real(dp) , allocatable , dimension(:) :: space
    type(pool1d) , pointer :: next
    type(pool1d) , pointer :: prev
  end type pool1d

  type pool2d
    real(dp) , allocatable , dimension(:,:) :: space
    type(pool2d) , pointer :: next
    type(pool2d) , pointer :: prev
  end type pool2d

  type pool3d
    real(dp) , allocatable , dimension(:,:,:) :: space
    type(pool3d) , pointer :: next
    type(pool3d) , pointer :: prev
  end type pool3d

  type pool4d
    real(dp) , allocatable , dimension(:,:,:,:) :: space
    type(pool4d) , pointer :: next
    type(pool4d) , pointer :: prev
  end type pool4d

  type pool5d
    real(dp) , allocatable , dimension(:,:,:,:,:) :: space
    type(pool5d) , pointer :: next
    type(pool5d) , pointer :: prev
  end type pool5d

  type(pool1d) , pointer :: root1d , curr1d , prev1d , next1d , last1d
  type(pool2d) , pointer :: root2d , curr2d , prev2d , next2d , last2d
  type(pool3d) , pointer :: root3d , curr3d , prev3d , next3d , last3d
  type(pool4d) , pointer :: root4d , curr4d , prev4d , next4d , last4d
  type(pool5d) , pointer :: root5d , curr5d , prev5d , next5d , last5d

  integer :: ials

  contains

  subroutine memory_init
    implicit none
    if ( debug_level > 2 ) then
      call mall_set()
    end if
    allocate(root1d, stat=ials)
    call checkalloc('memory_init: root1d')
    allocate(root2d, stat=ials)
    call checkalloc('memory_init: root2d')
    allocate(root3d, stat=ials)
    call checkalloc('memory_init: root3d')
    allocate(root4d, stat=ials)
    call checkalloc('memory_init: root4d')
    allocate(root5d, stat=ials)
    call checkalloc('memory_init: root5d')
    nullify(root1d%next)
    nullify(root1d%prev)
    nullify(root2d%next)
    nullify(root2d%prev)
    nullify(root3d%next)
    nullify(root3d%prev)
    nullify(root4d%next)
    nullify(root4d%prev)
    nullify(root5d%next)
    nullify(root5d%prev)
    last1d => root1d
    last2d => root2d
    last3d => root3d
    last4d => root4d
    last5d => root5d
  end subroutine memory_init

  subroutine getmem1d(a,nv,what)
    implicit none
    integer , intent (in) :: nv
    real(dp) , intent(out) , dimension(:) , pointer :: a
    character (len=*) , intent(in) :: what
    character (len=64) :: cspace
    write(cspace,'(i8)') nv
    allocate(last1d%space(nv),stat=ials)
    last1d%space(:) = d_zero
    call checkalloc('getmem1d : '//what//' : '//cspace)
    call mall_mci(last1d%space,'pool1d')
    a => last1d%space
    curr1d => last1d
    allocate(last1d%next)
    last1d => last1d%next
    last1d%prev => curr1d
    nullify(last1d%next)
  end subroutine getmem1d

  subroutine getmem2d(a,nv,nw,what)
    implicit none
    integer , intent (in) :: nv , nw
    real(dp) , intent(out) , dimension(:,:) , pointer :: a
    character (len=*) , intent(in) :: what
    character (len=64) :: cspace
    write(cspace,'(i8,a,i8)') nv, 'x', nw
    allocate(last2d%space(nv,nw),stat=ials)
    last2d%space(:,:) = d_zero
    call checkalloc('getmem2d : '//what//' : '//cspace)
    call mall_mci(last2d%space,'pool2d')
    a => last2d%space
    curr2d => last2d
    allocate(last2d%next)
    last2d => last2d%next
    last2d%prev => curr2d
    nullify(last2d%next)
  end subroutine getmem2d

  subroutine getmem3d(a,nv,nw,nx,what)
    implicit none
    integer , intent (in) :: nv , nw , nx
    real(dp) , intent(out) , dimension(:,:,:) , pointer :: a
    character (len=*) , intent(in) :: what
    character (len=64) :: cspace
    write(cspace,'(i8,a,i8,a,i8)') nv, 'x', nw, 'x', nx
    allocate(last3d%space(nv,nw,nx),stat=ials)
    last3d%space(:,:,:) = d_zero
    call checkalloc('getmem3d : '//what//' : '//cspace)
    call mall_mci(last3d%space,'pool3d')
    a => last3d%space
    curr3d => last3d
    allocate(last3d%next)
    last3d => last3d%next
    last3d%prev => curr3d
    nullify(last3d%next)
  end subroutine getmem3d

  subroutine getmem4d(a,nv,nw,nx,ny,what)
    implicit none
    integer , intent (in) :: nv , nw , nx , ny
    real(dp) , intent(out) , dimension(:,:,:,:) , pointer :: a
    character (len=*) , intent(in) :: what
    character (len=64) :: cspace
    integer :: i
    write(cspace,'(i8,a,i8,a,i8,a,i8)') nv, 'x', nw, 'x', nx, 'x', ny
    allocate(last4d%space(nv,nw,nx,ny),stat=ials)
    last4d%space(:,:,:,:) = d_zero
    call checkalloc('getmem4d : '//what//' : '//cspace)
    do i = 1 , nv
      call mall_mci(last4d%space(i,:,:,:),'pool4d')
    end do
    a => last4d%space
    curr4d => last4d
    allocate(last4d%next)
    last4d => last4d%next
    last4d%prev => curr4d
    nullify(last4d%next)
  end subroutine getmem4d

  subroutine getmem5d(a,nv,nw,nx,ny,nz,what)
    implicit none
    integer , intent (in) :: nv , nw , nx , ny , nz
    real(dp) , intent(out) , dimension(:,:,:,:,:) , pointer :: a
    character (len=*) , intent(in) :: what
    character (len=64) :: cspace
    integer :: i , j
    write(cspace,'(i8,a,i8,a,i8,a,i8,a,i8)') &
              nv, 'x', nw, 'x', nx, 'x', ny, 'x', nz
    allocate(last5d%space(nv,nw,nx,ny,nz),stat=ials)
    last5d%space(:,:,:,:,:) = d_zero
    call checkalloc('getmem5d : '//what//' : '//cspace)
    do i = 1 , nv
      do j = 1 , nw
        call mall_mci(last5d%space(i,j,:,:,:),'pool5d')
      end do
    end do
    a => last5d%space
    curr5d => last5d
    allocate(last5d%next)
    last5d => last5d%next
    last5d%prev => curr5d
    nullify(last5d%next)
  end subroutine getmem5d

  subroutine relmem1d(a)
    implicit none
    real(dp) , intent(in) , dimension(:) , pointer :: a
    if ( .not. associated(a) ) return
    curr1d => root1d
    do while ( associated(curr1d) )
      if ( associated(a,curr1d%space) ) then
        call mall_mco(curr1d%space,'pool1d')
        deallocate(curr1d%space)
        nullify(a)
        if ( associated(curr1d%prev) ) then
          prev1d => curr1d%prev
          if ( associated(curr1d%next) ) then
            next1d => curr1d%next
            next1d%prev => prev1d
            prev1d%next => next1d
          else
            nullify(prev1d%next)
            last1d => prev1d
          end if
          deallocate(curr1d)
        else
          if ( associated(curr1d%next) ) then
            next1d => curr1d%next
            nullify(next1d%prev)
            root1d => next1d
            deallocate(curr1d)
          end if
        end if
        return
      end if
      curr1d => curr1d%next
    end do
    write (stderr,*) 'This should not happen! 1D Array not found.'
    call die('mod_memutil','relmem1d',1)
  end subroutine relmem1d

  subroutine relmem2d(a)
    implicit none
    real(dp) , intent(in) , dimension(:,:) , pointer :: a
    if ( .not. associated(a) ) return
    curr2d => root2d
    do while ( associated(curr2d) )
      if ( associated(a,curr2d%space) ) then
        call mall_mco(curr2d%space,'pool2d')
        deallocate(curr2d%space)
        nullify(a)
        if ( associated(curr2d%prev) ) then
          prev2d => curr2d%prev
          if ( associated(curr2d%next) ) then
            next2d => curr2d%next
            next2d%prev => prev2d
            prev2d%next => next2d
          else
            nullify(prev2d%next)
            last2d => prev2d
          end if
          deallocate(curr2d)
        else
          if ( associated(curr2d%next) ) then
            next2d => curr2d%next
            nullify(next2d%prev)
            root2d => next2d
            deallocate(curr2d)
          end if
        end if
        return
      end if
      curr2d => curr2d%next
    end do
    write (stderr,*) 'This should not happen! 2D Array not found.'
    call die('mod_memutil','relmem2d',1)
  end subroutine relmem2d

  subroutine relmem3d(a)
    implicit none
    real(dp) , intent(in) , dimension(:,:,:) , pointer :: a
    if ( .not. associated(a) ) return
    curr3d => root3d
    do while ( associated(curr3d) )
      if ( associated(a,curr3d%space) ) then
        call mall_mco(curr3d%space,'pool3d')
        deallocate(curr3d%space)
        nullify(a)
        if ( associated(curr3d%prev) ) then
          prev3d => curr3d%prev
          if ( associated(curr3d%next) ) then
            next3d => curr3d%next
            next3d%prev => prev3d
            prev3d%next => next3d
          else
            nullify(prev3d%next)
            last3d => prev3d
          end if
          deallocate(curr3d)
        else
          if ( associated(curr3d%next) ) then
            next3d => curr3d%next
            nullify(next3d%prev)
            root3d => next3d
            deallocate(curr3d)
          end if
        end if
        return
      end if
      curr3d => curr3d%next
    end do
    write (stderr,*) 'This should not happen! 3D Array not found.'
    call die('mod_memutil','relmem3d',1)
  end subroutine relmem3d

  subroutine relmem4d(a)
    implicit none
    real(dp) , intent(in) , dimension(:,:,:,:) , pointer :: a
    integer :: i
    if ( .not. associated(a) ) return
    curr4d => root4d
    do while ( associated(curr4d) )
      if ( associated(a,curr4d%space) ) then
        do i = 1 , size(curr4d%space,1)
          call mall_mco(curr4d%space(i,:,:,:),'pool4d')
        end do
        deallocate(curr4d%space)
        nullify(a)
        if ( associated(curr4d%prev) ) then
          prev4d => curr4d%prev
          if ( associated(curr4d%next) ) then
            next4d => curr4d%next
            next4d%prev => prev4d
            prev4d%next => next4d
          else
            nullify(prev4d%next)
            last4d => prev4d
          end if
          deallocate(curr4d)
        else
          if ( associated(curr4d%next) ) then
            next4d => curr4d%next
            nullify(next4d%prev)
            root4d => next4d
            deallocate(curr4d)
          end if
        end if
        return
      end if
      curr4d => curr4d%next
    end do
    write (stderr,*) 'This should not happen! 4D Array not found.'
    call die('mod_memutil','relmem4d',1)
  end subroutine relmem4d

  subroutine relmem5d(a)
    implicit none
    real(dp) , intent(in) , dimension(:,:,:,:,:) , pointer :: a
    integer :: i , j
    if ( .not. associated(a) ) return
    curr5d => root5d
    do while ( associated(curr5d) )
      if ( associated(a,curr5d%space) ) then
        do i = 1 , size(curr5d%space,1)
          do j = 1 , size(curr5d%space,2)
            call mall_mco(curr5d%space(i,j,:,:,:),'pool5d')
          end do
        end do
        deallocate(curr5d%space)
        nullify(a)
        if ( associated(curr5d%prev) ) then
          prev5d => curr5d%prev
          if ( associated(curr5d%next) ) then
            next5d => curr5d%next
            next5d%prev => prev5d
            prev5d%next => next5d
          else
            nullify(prev5d%next)
            last5d => prev5d
          end if
          deallocate(curr5d)
        else
          if ( associated(curr5d%next) ) then
            next5d => curr5d%next
            nullify(next5d%prev)
            root5d => next5d
            deallocate(curr5d)
          end if
        end if
        return
      end if
      curr5d => curr5d%next
    end do
    write (stderr,*) 'This should not happen! 5D Array not found.'
    call die('mod_memutil','relmem5d',1)
  end subroutine relmem5d

  subroutine memory_destroy
    implicit none
    integer :: i , j
    curr1d => root1d
    do while ( associated(curr1d) )
      if ( allocated(curr1d%space) ) then
        call mall_mco(curr1d%space,'pool1d')
        deallocate(curr1d%space)
      end if
      next1d => curr1d%next
      deallocate(curr1d)
      curr1d => next1d
    end do
    curr2d => root2d
    do while ( associated(curr2d) )
      if ( allocated(curr2d%space) ) then
        call mall_mco(curr2d%space,'pool2d')
        deallocate(curr2d%space)
      end if
      next2d => curr2d%next
      deallocate(curr2d)
      curr2d => next2d
    end do
    curr3d => root3d
    do while ( associated(curr3d) )
      if ( allocated(curr3d%space) ) then
        call mall_mco(curr3d%space,'pool3d')
        deallocate(curr3d%space)
      end if
      next3d => curr3d%next
      deallocate(curr3d)
      curr3d => next3d
    end do
    curr4d => root4d
    do while ( associated(curr4d) )
      if ( allocated(curr4d%space) ) then
        do i = 1 , size(curr4d%space,1)
          call mall_mco(curr4d%space(i,:,:,:),'pool4d')
        end do
        deallocate(curr4d%space)
      end if
      next4d => curr4d%next
      deallocate(curr4d)
      curr4d => next4d
    end do
    curr5d => root5d
    do while ( associated(curr5d) )
      if ( allocated(curr5d%space) ) then
        do i = 1 , size(curr5d%space,1)
          do j = 1 , size(curr5d%space,2)
            call mall_mco(curr5d%space(i,j,:,:,:),'pool5d')
          end do
        end do
        deallocate(curr5d%space)
      end if
      next5d => curr5d%next
      deallocate(curr5d)
      curr5d => next5d
    end do
    if (debug_level > 2 ) then
      call mall_flush(stdout)
      call mall_set(.false.)
    end if
  end subroutine memory_destroy

  subroutine checkalloc(message)
    implicit none
    character (len=*) :: message
    if ( ials /= 0 ) then
      write (stderr,*) 'Memory error in allocation.'
      call die('mod_memutil',message,1)
    end if
  end subroutine checkalloc

end module mod_memutil
