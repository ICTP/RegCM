!
!(c) Matthew Kennel, Institute for Nonlinear Science (2004)
!
! The KDTREE2 software is licensed under the terms of the Academic Free
! Software License, listed herein.  In addition, users of this software
! must give appropriate citation in relevant technical documentation or
! journal paper to the author, Matthew B. Kennel, Institute For
! Nonlinear Science, preferably via a reference to the www.arxiv.org
! repository of this document, {\tt www.arxiv.org e-print:
! physics/0408067}.  This requirement will be deemed to be advisory and
! not mandatory as is necessary to permit the free inclusion of the
! present software with any software licensed under the terms of any
! version of the GNU General Public License, or GNU Library General
! Public License.
!
! Academic Free License
! Version 1.1
!
! This Academic Free License applies to any original work of authorship
! (the "Original Work") whose owner (the "Licensor") has placed the
! following notice immediately following the copyright notice for the
! Original Work: "Licensed under the Academic Free License version 1.1."
!
! Grant of License. Licensor hereby grants to any person obtaining a
! copy of the Original Work ("You") a world-wide, royalty-free,
! non-exclusive, perpetual, non-sublicenseable license (1) to use, copy,
! modify, merge, publish, perform, distribute and/or sell copies of the
! Original Work and derivative works thereof, and (2) under patent
! claims owned or controlled by the Licensor that are embodied in the
! Original Work as furnished by the Licensor, to make, use, sell and
! offer for sale the Original Work and derivative works thereof, subject
! to the following conditions.
!
! Right of Attribution. Redistributions of the Original Work must
! reproduce all copyright notices in the Original Work as furnished by
! the Licensor, both in the Original Work itself and in any
! documentation and/or other materials provided with the distribution of
! the Original Work in executable form.
!
! Exclusions from License Grant. Neither the names of Licensor, nor the
! names of any contributors to the Original Work, nor any of their
! trademarks or service marks, may be used to endorse or promote
! products derived from this Original Work without express prior written
! permission of the Licensor.
!
! WARRANTY AND DISCLAIMERS. LICENSOR WARRANTS THAT THE COPYRIGHT IN AND
! TO THE ORIGINAL WORK IS OWNED BY THE LICENSOR OR THAT THE ORIGINAL
! WORK IS DISTRIBUTED BY LICENSOR UNDER A VALID CURRENT LICENSE FROM THE
! COPYRIGHT OWNER. EXCEPT AS EXPRESSLY STATED IN THE IMMEDIATELY
! PRECEEDING SENTENCE, THE ORIGINAL WORK IS PROVIDED UNDER THIS LICENSE
! ON AN "AS IS" BASIS, WITHOUT WARRANTY, EITHER EXPRESS OR IMPLIED,
! INCLUDING, WITHOUT LIMITATION, THE WARRANTY OF NON-INFRINGEMENT AND
! WARRANTIES THAT THE ORIGINAL WORK IS MERCHANTABLE OR FIT FOR A
! PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY OF THE ORIGINAL
! WORK IS WITH YOU. THIS DISCLAIMER OF WARRANTY CONSTITUTES AN ESSENTIAL
! PART OF THIS LICENSE. NO LICENSE TO ORIGINAL WORK IS GRANTED HEREUNDER
! EXCEPT UNDER THIS DISCLAIMER.
!
! LIMITATION OF LIABILITY. UNDER NO CIRCUMSTANCES AND UNDER NO LEGAL
! THEORY, WHETHER TORT (INCLUDING NEGLIGENCE), CONTRACT, OR OTHERWISE,
! SHALL THE LICENSOR BE LIABLE TO ANY PERSON FOR ANY DIRECT, INDIRECT,
! SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES OF ANY CHARACTER ARISING
! AS A RESULT OF THIS LICENSE OR THE USE OF THE ORIGINAL WORK INCLUDING,
! WITHOUT LIMITATION, DAMAGES FOR LOSS OF GOODWILL, WORK STOPPAGE,
! COMPUTER FAILURE OR MALFUNCTION, OR ANY AND ALL OTHER COMMERCIAL
! DAMAGES OR LOSSES, EVEN IF SUCH PERSON SHALL HAVE BEEN INFORMED OF THE
! POSSIBILITY OF SUCH DAMAGES. THIS LIMITATION OF LIABILITY SHALL NOT
! APPLY TO LIABILITY FOR DEATH OR PERSONAL INJURY RESULTING FROM SUCH
! PARTY'S NEGLIGENCE TO THE EXTENT APPLICABLE LAW PROHIBITS SUCH
! LIMITATION. SOME JURISDICTIONS DO NOT ALLOW THE EXCLUSION OR
! LIMITATION OF INCIDENTAL OR CONSEQUENTIAL DAMAGES, SO THIS EXCLUSION
! AND LIMITATION MAY NOT APPLY TO YOU.
!
! License to Source Code. The term "Source Code" means the preferred
! form of the Original Work for making modifications to it and all
! available documentation describing how to access and modify the
! Original Work. Licensor hereby agrees to provide a machine-readable
! copy of the Source Code of the Original Work along with each copy of
! the Original Work that Licensor distributes. Licensor reserves the
! right to satisfy this obligation by placing a machine-readable copy of
! the Source Code in an information repository reasonably calculated to
! permit inexpensive and convenient access by You for as long as
! Licensor continues to distribute the Original Work, and by publishing
! the address of that information repository in a notice immediately
! following the copyright notice that applies to the Original Work.
!
! Mutual Termination for Patent Action. This License shall terminate
! automatically and You may no longer exercise any of the rights granted
! to You by this License if You file a lawsuit in any court alleging
! that any OSI Certified open source software that is licensed under any
! license containing this "Mutual Termination for Patent Action" clause
! infringes any patent claims that are essential to use that software.
!
! This license is Copyright (C) 2002 Lawrence E. Rosen. All rights
! reserved. Permission is hereby granted to copy and distribute this
! license without modification. This license may not be modified without
! the express written permission of its copyright owner.
!

module mod_kdtree2_priority

  use mod_realkinds
  use mod_message

  implicit none

  private

  integer , parameter :: sp = rk4
  integer , parameter :: dp = rk8
  integer , parameter :: kdkind = rk8

  !
  ! maintain a priority queue (PQ) of data, pairs of 'priority/payload',
  ! implemented with a binary heap.  This is the type, and the 'dis' field
  ! is the priority.
  !
  type kdtree2_result
    ! a pair of distances, indexes
    real(kdkind) :: dis  !=0.0
    integer :: idx       !=-1 Initializers cause some bugs in compilers.
  end type kdtree2_result
  !
  ! A heap-based priority queue lets one efficiently implement the following
  ! operations, each in log(N) time, as opposed to linear time.
  !
  ! 1)  add a datum (push a datum onto the queue, increasing its length)
  ! 2)  return the priority value of the maximum priority element
  ! 3)  pop-off (and delete) the element with the maximum priority, decreasing
  !     the size of the queue.
  ! 4)  replace the datum with the maximum priority with a supplied datum
  !     (of either higher or lower priority), maintaining the size of the
  !     queue.
  !
  ! In the k-d tree case, the 'priority' is the square distance of a point in
  ! the data set to a reference point.   The goal is to keep the smallest M
  ! distances to a reference point.  The tree algorithm searches terminal
  ! nodes to decide whether to add points under consideration.
  !
  ! A priority queue is useful here because it lets one quickly return the
  ! largest distance currently existing in the list.  If a new candidate
  ! distance is smaller than this, then the new candidate ought to replace
  ! the old candidate.  In priority queue terms, this means removing the
  ! highest priority element, and inserting the new one.
  !
  ! Algorithms based on Cormen, Leiserson, Rivest, _Introduction
  ! to Algorithms_, 1990, with further optimization by the author.
  !
  ! Originally informed by a C implementation by Sriranga Veeraraghavan.
  !
  ! This module is not written in the most clear way, but is implemented such
  ! for speed, as it its operations will be called many times during searches
  ! of large numbers of neighbors.
  !
  type pq
    !
    ! The priority queue consists of elements
    ! priority(1:heap_size), with associated payload(:).
    !
    ! There are heap_size active elements.
    ! Assumes the allocation is always sufficient.  Will NOT increase it
    ! to match.
    integer :: heap_size = 0
    type(kdtree2_result) , pointer , dimension(:) :: elems
  end type pq

  public :: kdtree2_result

  public :: pq
  public :: pq_create
  public :: pq_delete , pq_insert
  public :: pq_extract_max , pq_max , pq_replace_max , pq_maxpri

  contains

  type(pq) function pq_create(results_in) result(res)
    implicit none
    !
    ! Create a priority queue from ALREADY allocated
    ! array pointers for storage.  NOTE! It will NOT
    ! add any alements to the heap, i.e. any existing
    ! data in the input arrays will NOT be used and may
    ! be overwritten.
    !
    ! usage:
    !    real(kdkind), pointer :: x(:)
    !    integer, pointer :: k(:)
    !    allocate(x(1000),k(1000))
    !    pq => pq_create(x,k)
    !
    type(kdtree2_result) , target , dimension(:) :: results_in
    !
    !
    integer :: nalloc

    nalloc = size(results_in,1)
    if ( nalloc < 1 ) then
      call die('PQ_CREATE','Error, input arrays must be allocated.',1)
    end if
    res%elems => results_in
    res%heap_size = 0
    return
  end function pq_create

  !
  ! operations for getting parents and left + right children
  ! of elements in a binary heap.
  !

!
! These are written inline for speed.
!
!  integer function parent(i)
!    integer, intent(in) :: i
!    parent = (i/2)
!    return
!  end function parent

!  integer function left(i)
!    integer, intent(in) ::i
!    left = (2*i)
!    return
!  end function left

!  integer function right(i)
!    integer, intent(in) :: i
!    right = (2*i)+1
!    return
!  end function right

!  logical function compare_priority(p1,p2)
!    real(kdkind), intent(in) :: p1, p2
!
!    compare_priority = (p1 .gt. p2)
!    return
!  end function compare_priority

  subroutine heapify(a,i_in)
    implicit none
    !
    ! take a heap rooted at 'i' and force it to be in the
    ! heap canonical form.   This is performance critical
    ! and has been tweaked a little to reflect this.
    !
    type(pq) , pointer :: a
    integer , intent(in) :: i_in
    integer :: i, l, r, largest
    real(kdkind) :: pri_i , pri_l , pri_r , pri_largest
    type(kdtree2_result) :: temp

    i = i_in

    bigloop: &
    do
      l = 2*i ! left(i)
      r = l+1 ! right(i)
      !
      ! set 'largest' to the index of either i, l, r
      ! depending on whose priority is largest.
      !
      ! note that l or r can be larger than the heap size
      ! in which case they do not count.

      ! does left child have higher priority?
      if (l > a%heap_size) then
        ! we know that i is the largest as both l and r are invalid.
        exit
      else
        pri_i = a%elems(i)%dis
        pri_l = a%elems(l)%dis
        if (pri_l > pri_i) then
          largest = l
          pri_largest = pri_l
        else
          largest = i
          pri_largest = pri_i
        end if
        !
        ! between i and l we have a winner
        ! now choose between that and r.
        !
        if (r <= a%heap_size) then
          pri_r = a%elems(r)%dis
          if (pri_r > pri_largest) then
            largest = r
          end if
        end if
      end if

      if (largest /= i) then
        ! swap data in nodes largest and i, then heapify

        temp = a%elems(i)
        a%elems(i) = a%elems(largest)
        a%elems(largest) = temp
        !
        ! Canonical heapify() algorithm has tail-ecursive call:
        !
        !        call heapify(a,largest)
        ! we will simulate with cycle
        !
        i = largest
        cycle bigloop ! continue the loop
      else
        return   ! break from the loop
      end if
    end do bigloop
  end subroutine heapify

  subroutine pq_max(a,e)
    implicit none
    !
    ! return the priority and its payload of the maximum priority element
    ! on the queue, which should be the first one, if it is
    ! in heapified form.
    !
    type(pq) , pointer :: a
    type(kdtree2_result) , intent(out) :: e
    if ( a%heap_size > 0 ) then
      e = a%elems(1)
    else
      call die('PQ_MAX','ERROR, heap_size < 1',1)
    end if
  end subroutine pq_max

  real(kdkind) function pq_maxpri(a)
    implicit none
    type(pq) , pointer :: a
    if ( a%heap_size > 0 ) then
      pq_maxpri = a%elems(1)%dis
    else
      call die('PQ_MAX_PRI','ERROR, heapsize < 1',1)
    end if
  end function pq_maxpri

  subroutine pq_extract_max(a,e)
    implicit none
    !
    ! return the priority and payload of maximum priority
    ! element, and remove it from the queue.
    ! (equivalent to 'pop()' on a stack)
    !
    type(pq) , pointer :: a
    type(kdtree2_result) , intent(out) :: e

    if ( a%heap_size >= 1 ) then
      !
      ! return max as first element
      !
      e = a%elems(1)
      !
      ! move last element to first
      !
      a%elems(1) = a%elems(a%heap_size)
      a%heap_size = a%heap_size-1
      call heapify(a,1)
      return
    else
      call die('PQ_EXTRACT_MAX','error, attempted to pop non-positive PQ',1)
    end if
  end subroutine pq_extract_max

  real(kdkind) function pq_insert(a,dis,idx)
    implicit none
    !
    ! Insert a new element and return the new maximum priority,
    ! which may or may not be the same as the old maximum priority.
    !
    type(pq) , pointer :: a
    real(kdkind) , intent(in) :: dis
    integer , intent(in) :: idx
    ! type(kdtree2_result) , intent(in) :: e
    !
    integer :: i , isparent
    real(kdkind) :: parentdis
    !

    ! if (a%heap_size .ge. a%max_elems) then
    !   call die('PQ_INSERT', &
    !            'error, attempt made to insert element on full PQ',1)
    ! else
    a%heap_size = a%heap_size + 1
    i = a%heap_size

    do while ( i > 1 )
      isparent = int(i/2)
      parentdis = a%elems(isparent)%dis
      if ( dis > parentdis ) then
        ! move what was in i's parent into i.
        a%elems(i)%dis = parentdis
        a%elems(i)%idx = a%elems(isparent)%idx
        i = isparent
      else
        exit
      end if
    end do

    ! insert the element at the determined position
    a%elems(i)%dis = dis
    a%elems(i)%idx = idx

    pq_insert = a%elems(1)%dis
    ! end if
  end function pq_insert

!  subroutine pq_adjust_heap(a,i)
!    implicit none
!    type(pq) , pointer :: a
!    integer , intent(in) :: i
!    !
!    ! nominally arguments (a,i), but specialize for a=1
!    !
!    ! This routine assumes that the trees with roots 2 and 3 are
!    ! already heaps, i.e. the children of '1' are heaps.
!    !  When the procedure is completed, the tree rooted at 1 is a heap.
!    real(kdkind) :: prichild
!    integer :: parent , child , n
!
!    type(kdtree2_result) :: e
!
!    e = a%elems(i)
!
!    parent = i
!    child = 2*i
!    n = a%heap_size
!
!    do while ( child <= N )
!      if ( child < n ) then
!        if ( a%elems(child)%dis < a%elems(child+1)%dis ) then
!          child = child+1
!        end if
!      end if
!      prichild = a%elems(child)%dis
!      if ( e%dis >= prichild ) then
!        exit
!      else
!        ! move child into parent.
!        a%elems(parent) = a%elems(child)
!        parent = child
!        child = 2*parent
!      end if
!    end do
!    a%elems(parent) = e
!  end subroutine pq_adjust_heap

  real(kdkind) function pq_replace_max(a,dis,idx)
    implicit none
    !
    ! Replace the extant maximum priority element
    ! in the PQ with (dis,idx).  Return
    ! the new maximum priority, which may be larger
    ! or smaller than the old one.
    !
    type(pq),pointer         :: a
    real(kdkind), intent(in) :: dis
    integer, intent(in) :: idx
    ! type(kdtree2_result), intent(in) :: e
    ! not tested as well!

    integer :: parent , child , n
    real(kdkind) :: prichild , prichildp1

    type(kdtree2_result) :: etmp

    if ( .true. ) then
      n = a%heap_size
      if ( n >= 1 ) then
        parent = 1
        child = 2
        loop: do while ( child <= n )
          prichild = a%elems(child)%dis

          !
          ! posibly child+1 has higher priority, and if
          ! so, get it, and increment child.
          !

          if ( child < n ) then
            prichildp1 = a%elems(child+1)%dis
            if ( prichild < prichildp1 ) then
              child = child+1
              prichild = prichildp1
            end if
          end if

          if ( dis >= prichild ) then
            exit loop
            ! we have a proper place for our new element,
            ! bigger than either children's priority.
          else
            ! move child into parent.
            a%elems(parent) = a%elems(child)
            parent = child
            child = 2*parent
          end if
        end do loop
        a%elems(parent)%dis = dis
        a%elems(parent)%idx = idx
        pq_replace_max = a%elems(1)%dis
      else
        a%elems(1)%dis = dis
        a%elems(1)%idx = idx
        pq_replace_max = dis
      end if
    else
      !
      ! slower version using elementary pop and push operations.
      !
      call pq_extract_max(a,etmp)
      etmp%dis = dis
      etmp%idx = idx
      pq_replace_max = pq_insert(a,dis,idx)
    end if
  end function pq_replace_max

  subroutine pq_delete(a,i)
    implicit none
    !
    ! delete item with index 'i'
    !
    type(pq) , pointer :: a
    integer :: i

    if ( (i < 1) .or. (i > a%heap_size) ) then
      call die('PQ_DELETE','error, attempt to remove out of bounds element.',1)
    end if

    ! swap the item to be deleted with the last element
    ! and shorten heap by one.
    a%elems(i) = a%elems(a%heap_size)
    a%heap_size = a%heap_size - 1

    call heapify(a,i)

  end subroutine pq_delete

end module mod_kdtree2_priority

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
