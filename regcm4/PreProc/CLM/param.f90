      subroutine param(nx,ny,kz,ds,clat,clon,xplat,xplon,xlat,xlon,     &
                     & varmin,varmax,xlat1d,xlon1d,xlonmin,xlonmax,     &
                     & xlatmin,xlatmax,iadim,ndim)
 
      implicit none
!
! Dummy arguments
!
      real :: clat , clon , ds , xlatmax , xlatmin , xlonmax , xlonmin ,&
            & xplat , xplon
      integer :: kz , ndim , nx , ny
      integer , dimension(ndim) :: iadim
      real , dimension(ndim) :: varmax , varmin
      real , dimension(ny,nx) :: xlat , xlon
      real , dimension(ny) :: xlat1d
      real , dimension(nx) :: xlon1d
      intent (in) kz , ndim , nx , ny , xlat , xlon
      intent (out) iadim , varmax , varmin , xlat1d , xlatmax , xlatmin ,&
                 & xlon1d , xlonmax , xlonmin
!
! Local variables
!
      integer :: i , j
!
      varmin(1) = xlon(ny/2,1)
      varmin(2) = xlat(1,nx/2)
      varmin(3) = 1 !1050.
      varmax(1) = xlon(ny/2,nx)
      varmax(2) = xlat(ny,nx/2)
      varmax(3) = kz !1050.
      iadim(1) = nx
      iadim(2) = ny
      iadim(3) = kz
      do i = 1 , nx
        xlon1d(i) = xlon(ny/2,i)
      end do
      do j = 1 , ny
        xlat1d(j) = xlat(j,nx/2)
      end do
      xlonmin = minval(xlon)
      xlonmax = maxval(xlon)
      xlatmin = minval(xlat)
      xlatmax = maxval(xlat)
 
      end subroutine param
