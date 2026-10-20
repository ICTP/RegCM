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

! Reordering between subgrid (nnsg,j,i) and fine (jsg,isg) layouts.
!
! Template included by mod_mppparam.F90 once per real kind: the TK_*
! macros (type, MPI datatype, buffers and routine names) are set by
! mpp_tk_real8.inc or mpp_tk_real4.inc.
! It must be included AFTER the type/kind is included

  subroutine REORDER_ADD_SUBGRID_2D_TK(var3,var2,mask)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: var3
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: var2
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2(jj,ii) = var2(jj,ii) + var3((n2-1)*nsg+n1,j,i)
              else
                var2(jj,ii) = TK_MISSVAL
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2(jj,ii) = var2(jj,ii) + var3((n2-1)*nsg+n1,j,i)
            end do
          end do
        end do
      end do
    end if
  end subroutine REORDER_ADD_SUBGRID_2D_TK

  subroutine REORDER_SUBGRID_2D_TK(var3,var2,mask)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: var3
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: var2
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2(jj,ii) = var3((n2-1)*nsg+n1,j,i)
              else
                var2(jj,ii) = TK_MISSVAL
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2(jj,ii) = var3((n2-1)*nsg+n1,j,i)
            end do
          end do
        end do
      end do
    end if
  end subroutine REORDER_SUBGRID_2D_TK

  subroutine REORDER_ADD_SUBGRID_2D3D_TK(var3,var2_3,l,mask)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: var3
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: var2_3
    integer(ik4), optional, intent(in) :: l
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2, ll
    ll = 1
    if ( present(l) ) ll = l
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2_3(jj,ii,ll) = var2_3(jj,ii,ll) + var3((n2-1)*nsg+n1,j,i)
              else
                var2_3(jj,ii,ll) = TK_MISSVAL
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2_3(jj,ii,ll) = var2_3(jj,ii,ll) + var3((n2-1)*nsg+n1,j,i)
            end do
          end do
        end do
      end do
    end if
  end subroutine REORDER_ADD_SUBGRID_2D3D_TK

  subroutine REORDER_SUBGRID_2D3D_TK(var3,var2_3,l,mask)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(in) :: var3
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: var2_3
    integer(ik4), optional, intent(in) :: l
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2, ll
    ll = 1
    if ( present(l) ) ll = l
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2_3(jj,ii,ll) = var3((n2-1)*nsg+n1,j,i)
              else
                var2_3(jj,ii,ll) = TK_MISSVAL
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2_3(jj,ii,ll) = var3((n2-1)*nsg+n1,j,i)
            end do
          end do
        end do
      end do
    end if
  end subroutine REORDER_SUBGRID_2D3D_TK

  subroutine REORDER_ADD_SUBGRID_3D_TK(var4,var2,l,mask)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(in) :: var4
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: var2
    integer(ik4), intent(in) :: l
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2(jj,ii) = var2(jj,ii) + var4((n2-1)*nsg+n1,j,i,l)
              else
                var2(jj,ii) = TK_MISSVAL
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2(jj,ii) = var2(jj,ii) + var4((n2-1)*nsg+n1,j,i,l)
            end do
          end do
        end do
      end do
    end if
  end subroutine REORDER_ADD_SUBGRID_3D_TK

  subroutine REORDER_SUBGRID_3D_TK(var4,var2,l,mask)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(in) :: var4
    TK_TYPE, pointer, contiguous, dimension(:,:), intent(inout) :: var2
    integer(ik4), intent(in) :: l
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2(jj,ii) = var4((n2-1)*nsg+n1,j,i,l)
              else
                var2(jj,ii) = TK_MISSVAL
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2(jj,ii) = var4((n2-1)*nsg+n1,j,i,l)
            end do
          end do
        end do
      end do
    end if
  end subroutine REORDER_SUBGRID_3D_TK

  subroutine REORDER_SUBGRID_4D_TK(var4,var3,mask)
    implicit none
    TK_TYPE, pointer, contiguous, dimension(:,:,:,:), intent(in) :: var4
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: var3
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, l, ii, jj, n1, n2
    if ( present(mask) ) then
      do l = 1, size(var4,4)
        do i = ici1, ici2
          do j = jci1, jci2
            do n2 = 1, nsg
              ii = (i-1) * nsg + n2
              do n1 = 1, nsg
                jj = (j-1) * nsg + n1
                if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                  var3(jj,ii,l) = var4((n2-1)*nsg+n1,j,i,l)
                else
                  var3(jj,ii,l) = TK_MISSVAL
                end if
              end do
            end do
          end do
        end do
      end do
    else
      do l = 1, size(var4,4)
        do i = ici1, ici2
          do j = jci1, jci2
            do n2 = 1, nsg
              ii = (i-1) * nsg + n2
              do n1 = 1, nsg
                jj = (j-1) * nsg + n1
                var3(jj,ii,l) = var4((n2-1)*nsg+n1,j,i,l)
              end do
            end do
          end do
        end do
      end do
    end if
  end subroutine REORDER_SUBGRID_4D_TK

  subroutine INPUT_REORDER_TK(m1,m2,j1,j2,i1,i2)
    implicit none
    TK_TYPE, dimension(:,:), intent(in) :: m1
    TK_TYPE, pointer, contiguous, dimension(:,:,:), intent(inout) :: m2
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: i, j, ii, jj, n1, n2
    do i = i1, i2
      do j = j1, j2
        do n2 = 1, nsg
          ii = (i-1) * nsg + n2
          do n1 = 1, nsg
            jj = (j-1) * nsg + n1
            m2((n2-1)*nsg+n1,j,i) = m1(jj,ii)
          end do
        end do
      end do
    end do
  end subroutine INPUT_REORDER_TK

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
