!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine rdheadicbc(iy,jx,nx,ny,kz,clat,clon,ds,pt,sigf,sigh, &
                          & sighrev,xplat,xplon,f,xmap,dmap,xlat,xlon,  &
                          & zs,zssd,ls,iin,inhead,ibyte)
 
      implicit none
!
! Dummy arguments
!
      real(4) :: clat , clon , ds , pt , xplat , xplon
      integer :: ibyte , iin , nx , iy , ny , jx , kz
      character(70) :: inhead
      real(4) , dimension(nx,ny) :: dmap , f , ls , xlat , xlon , xmap ,&
                              &  zs , zssd
      real(4) , dimension(kz+1) :: sigf
      real(4) , dimension(kz) :: sigh , sighrev
      intent (in) ibyte , iin , inhead , nx , iy , ny , jx , kz
      intent (out) dmap , f , ls , sighrev , xlat , xlon , xmap , zs ,  &
                 & zssd
      intent (inout) clat , clon , ds , pt , sigf , sigh , xplat , xplon
!
! Local variables
!
      real(4) :: grdfac
      integer :: i , ibigend , ierr , igrads , j , k , kk , ni , nj , nk
      character(6) :: proj
      real(4) , dimension(iy,jx) :: tmp2d
!
      open (iin,file=inhead,status='old',form='unformatted',            &
          & recl=iy*jx*ibyte,access='direct')
      read (iin,rec=1,iostat=ierr) ni , nj , nk , ds , clat , clon ,    &
                                 & xplat , xplon , grdfac , proj ,      &
                                 & sigf , pt , igrads , ibigend
      print * , 'ni,nj,nk,ds='
      print * , ni , nj , nk , ds
      print * , 'sigf='
      print * , sigf
      print * , 'pt,clat,clon,xplat,xplon,proj='
      print * , pt , clat , clon , xplat , xplon , proj
      if ( ni/=jx .or. nj/=iy .or. kz/=nk ) then
        print * , 'Grid Dimensions DO NOT MATCH'
        print * , '  jx=' , jx , 'ix=' , iy , 'kx=' , kz
        print * , '  ni=' , ni , 'nj=' , nj , 'nk=' , nk
        print * , '  Also check ibyte in icbc.param: ibyte= ' , ibyte
        stop 'BAD DIMENSIONS (SUBROUTINE RDHEADICBC)'
      end if
!     print*,'ZS'
      read (iin,rec=2,iostat=ierr) tmp2d
      do j = 1 , ny
        do i = 1 , nx
          zs(i,j) = tmp2d(i+1,j+1)
        end do
      end do
!     print*,'ZSSD'
      read (iin,rec=3,iostat=ierr) tmp2d
      do j = 1 , ny
        do i = 1 , nx
          zssd(i,j) = tmp2d(i+1,j+1)
        end do
      end do
!     print*,'LU'
      read (iin,rec=4,iostat=ierr) tmp2d
      do j = 1 , ny
        do i = 1 , nx
          ls(i,j) = tmp2d(i+1,j+1)
        end do
      end do
!     print*,'XLAT'
      read (iin,rec=5,iostat=ierr) tmp2d
      do j = 1 , ny
        do i = 1 , nx
          xlat(i,j) = tmp2d(i+1,j+1)
        end do
      end do
!     print*,'XLON'
      read (iin,rec=6,iostat=ierr) tmp2d
      do j = 1 , ny
        do i = 1 , nx
          xlon(i,j) = tmp2d(i+1,j+1)
        end do
      end do
!     print*,'XMAP'
      read (iin,rec=9,iostat=ierr) tmp2d
      do j = 1 , ny
        do i = 1 , nx
          xmap(i,j) = tmp2d(i+1,j+1)
        end do
      end do
!     print*,'DMAP'
      read (iin,rec=10,iostat=ierr) tmp2d
      do j = 1 , ny
        do i = 1 , nx
          dmap(i,j) = tmp2d(i+1,j+1)
        end do
      end do
!     print*,'F'
      read (iin,rec=11,iostat=ierr) tmp2d
      do j = 1 , ny
        do i = 1 , nx
          f(i,j) = tmp2d(i+1,j+1)
        end do
      end do
 
      if ( ierr/=0 ) then
        print * , 'END OF FILE REACHED'
        print * , '  Check ibyte in postproc.param: ibyte= ' , ibyte
        stop 'EOF (SUBROUTINE RDHEADICBC)'
      end if
 
      do k = 1 , kz
        kk = kz - k + 1
        sigh(k) = (sigf(k)+sigf(k+1))/2.
        sighrev(kk) = sigh(k)
      end do
 
      end subroutine rdheadicbc
