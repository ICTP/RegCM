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

      module mod_interp
!
      use mod_constants
!
      implicit none
!
      private
!
      public :: interp , filter1plakes
!
      real(8) , dimension(4,4) :: c
      real(8) , dimension(16,16) :: wt
      integer , parameter :: maxbins = 20
      integer , dimension(maxbins) :: bincnt
      real(8) , dimension(maxbins) :: bmindist
      logical , dimension(2,maxbins) :: lndwt
      integer , dimension(2,maxbins) :: indwt
!
      data wt/1 , 0 , -3 , 2 , 4*0 , -3 , 0 , 9 , -6 , 2 , 0 , -6 , 4 , &
         & 8*0 , 3 , 0 , -9 , 6 , -2 , 0 , 6 , -4 , 10*0 , 9 , -6 ,     &
         & 2*0 , -6 , 4 , 2*0 , 3 , -2 , 6*0 , -9 , 6 , 2*0 , 6 , -4 ,  &
         & 4*0 , 1 , 0 , -3 , 2 , -2 , 0 , 6 , -4 , 1 , 0 , -3 , 2 ,    &
         & 8*0 , -1 , 0 , 3 , -2 , 1 , 0 , -3 , 2 , 10*0 , -3 , 2 ,     &
         & 2*0 , 3 , -2 , 6*0 , 3 , -2 , 2*0 , -6 , 4 , 2*0 , 3 , -2 ,  &
         & 0 , 1 , -2 , 1 , 5*0 , -3 , 6 , -3 , 0 , 2 , -4 , 2 , 9*0 ,  &
         & 3 , -6 , 3 , 0 , -2 , 4 , -2 , 10*0 , -3 , 3 , 2*0 , 2 , -2 ,&
         & 2*0 , -1 , 1 , 6*0 , 3 , -3 , 2*0 , -2 , 2 , 5*0 , 1 , -2 ,  &
         & 1 , 0 , -2 , 4 , -2 , 0 , 1 , -2 , 1 , 9*0 , -1 , 2 , -1 ,   &
         & 0 , 1 , -2 , 1 , 10*0 , 1 , -1 , 2*0 , -1 , 1 , 6*0 , -1 ,   &
         & 1 , 2*0 , 2 , -2 , 2*0 , -1 , 1/
!
       data lndwt /26*.false.,.true.,.true.,.true.,11*.false./
       data indwt /26*0,1,1,1,11*0/
!
      contains
!
      function inear(x,m,lwrap)
      implicit none
      integer :: inear
      real(8) , intent(in) :: x
      integer , intent(in) :: m
      logical , intent(in) :: lwrap
      if (.not. lwrap) then
        inear = min(max(nint(x),1),m)
      else
        inear = nint(x)
        if (inear < 1) then
          inear = m - inear
        end if
        if (inear > m) then
          inear = inear - m
        end if
      end if
      end function inear
!
      function jnear(y,n)
      implicit none
      integer :: jnear
      real(8) , intent(in) :: y
      integer , intent(in) :: n
      jnear = min(max(nint(y),1),n)
      end function jnear
!
      function ifloor(x,m,lwrap)
      implicit none
      integer :: ifloor
      real(8) , intent(in) :: x
      integer , intent(in) :: m
      logical , intent(in) :: lwrap
      if (.not. lwrap) then
        ifloor = min(max(floor(x),1),m)
      else
        ifloor = floor(x)
        if (ifloor < 1) then
          ifloor = m - ifloor
        end if
        if (ifloor > m) then
          ifloor = ifloor - m
        end if
      end if
      end function ifloor
!
      function jfloor(y,n)
      implicit none
      integer :: jfloor
      real(8) , intent(in) :: y
      integer , intent(in) :: n
      jfloor = min(max(floor(y),1),n)
      end function jfloor
!
      function nearpoint(x,y,m,n,grid,lwrap)
      implicit none
      real(8) :: nearpoint
      integer :: m , n
      real(8) :: x , y
      logical :: lwrap
      real(4) , dimension(m,n) :: grid
      intent (in) lwrap , m , n , grid , x , y
      nearpoint = grid(inear(x,m,lwrap),jnear(y,n))
      end function nearpoint
!
      function mostaround(x,y,m,n,grid,nbox,ibnty,h2opct,lwrap)
      implicit none
      real(8) :: mostaround
      integer , intent(in) :: m , n , nbox , ibnty
      real(8) , intent(in) :: x , y
      logical , intent(in) :: lwrap
      real(4) , intent(in) , dimension(m,n) :: grid
      real(4) , intent(in) :: h2opct
!
      real(8) , dimension(nbox*nbox) :: binval , bindist
      real(8) :: dist , rx , ry , wtp
      integer :: ii0 , jj0 , ii , jj
      integer :: totpoints , i , j , lastc , hbox

      hbox = nbox / 2
      totpoints = nbox*nbox
      ii0 = ifloor(x,m,lwrap)-hbox
      jj0 = jfloor(y,n)-hbox
      do i = 1 , nbox
        do j = 1 , nbox
          rx = ii0 + i - 1
          ry = jj0 + j - 1
          ii = ifloor(rx,m,lwrap)
          jj = jfloor(ry,n)
          binval((i-1)*nbox+j) = grid(ii,jj)
          bindist((i-1)*nbox+j) = sqrt((x-rx)**2+(y-ry)**2)
        end do
      end do
      bincnt = 0.0
      bmindist = 2*nbox
      do i = 1 , totpoints
        bincnt(int(binval(i))) = bincnt(int(binval(i))) + 1
        if (bindist(i) < bmindist(int(binval(i)))) &
          bmindist(int(binval(i))) = bindist(i)
      end do
!     Set point to land if less than fixed percent of water
      wtp = (sum(bincnt,mask=lndwt(ibnty,:))/totpoints)*100.0
      if (wtp > 0.0 .and. wtp < h2opct) then
        bincnt(indwt(ibnty,:)) = 0
      end if
      mostaround = -1
      lastc = -1
      do i = 1 , maxbins
        if (bincnt(i) > 0) then
          if (bincnt(i) > lastc) then
            mostaround = dble(i)
            dist = bmindist(i)
            lastc = bincnt(i)
          else if (bincnt(i) == lastc) then
            if (bmindist(i) < dist) then
              mostaround = dble(i)
              dist = bmindist(i)
              lastc = bincnt(i)
            end if
          end if
        end if
      end do
      end function mostaround
!
      function pctaround(x,y,m,n,grid,nbox,ival,lwrap)
      implicit none
      real(8) :: pctaround
      integer :: m , n , ival , nbox
      real(8) :: x , y , rx , ry
      logical :: lwrap
      real(4) , dimension(m,n) :: grid
      intent (in) lwrap , m , n , grid , x , y , ival
!
      integer :: ii0 , jj0 , ii , jj
      integer :: i , j
      real(8) :: pc

      pctaround = 0.0D0
      pc = nbox*nbox
      ii0 = ifloor(x,m,lwrap)
      jj0 = jfloor(y,n)
      do i = 1 , nbox
        do j = 1 , nbox
          rx = ii0 + i - nbox/2
          ry = jj0 + j - nbox/2
          ii = ifloor(rx,m,lwrap)
          jj = jfloor(ry,n)
          if (int(grid(ii,jj)) == ival) then
            pctaround = pctaround + 1
          end if
        end do
      end do
      pctaround = (pctaround / pc) * 100.0D0
      end function pctaround
!
      function bilinear(x,y,m,n,grid,lwrap)

      implicit none
!
      real(8) :: bilinear
      integer :: m , n
      real(8) :: x , y
      logical :: lwrap
      real(4) , dimension(m,n) :: grid
      intent (in) lwrap , m , n , grid , x , y
!
      real(8) :: dx, dy, p12, p03
      real(8) :: ii0, jj0, ii1, jj1, ii2, jj2, ii3, jj3
      integer :: i0, j0, i1, j1, i2, j2, i3, j3
!
!-----bilinear interpolation among four grid values
!
      if (.not. lwrap) then
        ii0 = min(max(floor(x),1),m)
        ii2 = min(max(ceiling(x),1),m)
        dx = (x-ii0)
      else
        ii0 = floor(x)
        ii2 = ceiling(x)
        dx = (x-ii0)
        if (ii0 < 1) then
          ii0 = m
        end if
        if (ii2 > m) then
          ii2 = 1
        end if
      end if
      ii1 = ii0
      ii3 = ii2
      jj0 = min(max(floor(y),1),n)
      dy = (y-jj0)
      jj3 = jj0
      jj1 = min(max(ceiling(y),1),n)
      jj2 = jj1

      i0 = int(ii0)
      j0 = int(jj0)
      i1 = int(ii1)
      j1 = int(jj1)
      i2 = int(ii2)
      j2 = int(jj2)
      i3 = int(ii3)
      j3 = int(jj3)

      p12 = dx*grid(i2,j2)+(1-dx)*grid(i1,j1)
      p03 = dx*grid(i3,j3)+(1-dx)*grid(i0,j0)
      bilinear = dy*p12+(1-dy)*p03

      end function bilinear
!
      function bicubic(x,y,m,n,grid,lwrap)
 
      implicit none
!
      real(8) :: bicubic
      integer :: m , n
      real(8) :: x , y
      logical :: lwrap
      real(4) , dimension(m,n) :: grid
      intent (in) grid , m , n , x , y , lwrap
!
      real(8) , dimension(4) :: f , f1 , f12 , f2
      real(8) :: xl , xu , yl , yu
      integer :: i , ii , j , mm , nn , im , imp1 , imn1

      mm = int(x)
      nn = int(y)
      if ( .not. lwrap) then
        mm = max(2, min(m-2,mm))
      end if
      nn = max(2, min(n-2,nn))

      xl = mm
      xu = mm + 1
      yl = nn
      yu = nn + 1
      do j = nn , nn + 1
        do i = mm , mm + 1
          ii = 1 + (i-mm) + 3*(j-nn)
          if ( ii==5 ) ii = 3
          if (lwrap) then
            im = i
            imp1 = i+1
            imn1 = i-1
            if ( i < 1 ) im = m-i
            if ( i > m ) im = i-m
            if ( imp1 < 1 ) imp1 = m-i+1
            if ( imp1 > m ) imp1 = i+1-m
            if ( imn1 < 1 ) imn1 = m-i-1
            if ( imn1 > m ) imn1 = i-1-m
          else
            im = i
            imp1 = i+1
            imn1 = i-1
          end if
          f(ii) = grid(im,j)
          f1(ii) = (grid(imp1,j)-grid(imn1,j))/(2D0)
          f2(ii) = (grid(im,j+1)-grid(im,j-1))/(2D0)
          f12(ii) = (grid(imp1,j+1)-grid(imp1,j-1)-grid(imn1,j+1)  &
                  & +grid(imn1,j-1))/(4D0)
        end do
      end do
 
      call bcuint(f,f1,f2,f12,xl,xu,yl,yu,x,y,bicubic)
 
      end function bicubic
!
      subroutine bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,a)
      implicit none
!
      real(8) :: a , x1 , x1l , x1u , x2 , x2l , x2u
      real(8) , dimension(4) :: y , y1 , y12 , y2
      intent (in) x1 , x1l , x1u , x2 , x2l , x2u , y , y1 , y12 , y2
      intent (out) a
!
      integer :: i
      real(8) :: t , u
!
      call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l)
      t = (x1-x1l)/(x1u-x1l)
      u = (x2-x2l)/(x2u-x2l)
      a = 0D0
      do i = 4 , 1 , -1
        a = t*a + ((c(i,4)*u+c(i,3))*u+c(i,2))*u + c(i,1)
      end do
      end subroutine bcuint
!
      subroutine bcucof(y,y1,y2,y12,d1,d2)
      implicit none
!
      real(8) :: d1 , d2
      real(8) , dimension(4) :: y , y1 , y12 , y2
      intent (in) d1 , d2 , y , y1 , y12 , y2
!
      real(8) , dimension(16) :: cl , x
      real(8) :: d1d2 , xx
      integer :: i , j , k , l

      d1d2 = d1*d2
      do i = 1 , 4
        x(i) = y(i)
        x(i+4) = y1(i)*d1
        x(i+8) = y2(i)*d2
        x(i+12) = y12(i)*d1d2
      end do
      do i = 1 , 16
        xx = 0D0
        do k = 1 , 16
          xx = xx + wt(i,k)*x(k)
        end do
        cl(i) = xx
      end do
      l = 0
      do i = 1 , 4
        do j = 1 , 4
          l = l + 1
          c(i,j) = cl(l)
        end do
      end do
      end subroutine bcucof
!
! Interpolates input regolar lat/lon grid on output model grid
!
      subroutine interp(iy,jx,xlat,xlon,omt,iniy,injx,milat,milon,imt, &
                        ntypec,itype,lwrap,lcross,ival,ibnty,h2opct)
 
      implicit none
!
      integer , intent(in) :: iy , jx , iniy , injx , ntypec , itype
      real(4) , intent(in) , dimension(iy, jx) :: xlat , xlon
      real(4) , intent(in) , dimension(injx, iniy) :: imt
      real(8) , intent(in) :: milat , milon
      logical , intent(in) :: lwrap , lcross
      integer , intent(in) , optional :: ival
      integer , intent(in) , optional :: ibnty
      real(4) , intent(in) , optional :: h2opct
      real(4) , intent(out) , dimension(iy, jx) :: omt
!
      integer :: nbox , ii , jj
      real(8) :: xx , yy , rinc
!
      rinc = 60.0D0/dble(ntypec)
!
      if (itype < 1 .or. itype > 5) then
        print *, 'Unknown interpolation type'
        stop
      end if
!
      do ii = 1 , iy
        do jj = 1 , jx
          yy = (dble(xlat(ii,jj))-milat)*rinc + 1.0D+00
          if (lcross) then
            xx = (mod((dble(xlon(ii,jj))+360.0D0),360.0D0)-milon) * &
                  rinc + 1.0D+00
          else
            xx = (dble(xlon(ii,jj))-milon)*rinc + 1.0D+00
          end if
 
!         yy and xx are the exact index values of a point i,j of the
!         mesoscale mesh when projected onto an earth-grid of lat_s
!         and lon_s for which terrain observations are available.  it
!         is assumed that the earth grid has equal spacing in both
!         latitude and longitude.
 
          select case (itype)
            case(1)
              omt(ii,jj) = bilinear(xx,yy,injx,iniy,imt,lwrap)
            case(2)
              omt(ii,jj) = bicubic(xx,yy,injx,iniy,imt,lwrap)
            case(3)
              omt(ii,jj) = nearpoint(xx,yy,injx,iniy,imt,lwrap)
            case(4,5)
              if (lwrap) then
                if (ii == 1 .or. ii == iy ) then 
                  omt(ii,jj) = nearpoint(xx,yy,injx,iniy,imt,lwrap)
                  cycle
                else
                  if (jj == jx) then
                    nbox = nint(max(abs(xlon(ii,jx-1)-(xlon(ii,1)+360))*&
                                rinc/2.0D0, 2.0D0))
                    nbox = min(nbox,8)
                  else if (jj == 1) then
                    nbox = nint(max(abs(xlon(ii,2)-(xlon(ii,jx)-360))*&
                                rinc/2.0D0, 2.0D0))
                    nbox = min(nbox,8)
                  else
                    nbox = nint(max(abs(xlon(ii,jj-1)-xlon(ii,jj+1))* &
                                rinc/2.0D0, 2.0D0))
                    nbox = min(nbox,8)
                  end if
                  nbox = nint(max(abs(xlat(ii-1,jj)-xlat(ii+1,jj))* &
                              rinc/2.0D0, dble(nbox)))
                end if
              else
                if (ii == 1 .or. jj == 1 .or. ii == iy .or. jj == jx) then
                  omt(ii,jj) = nearpoint(xx,yy,injx,iniy,imt,lwrap)
                  cycle
                else
                  nbox = nint(max(abs(xlon(ii,jj-1)-xlon(ii,jj+1))* &
                              rinc/2.0D0, 2.0D0))
                  nbox = nint(max(abs(xlat(ii-1,jj)-xlat(ii+1,jj))* &
                              rinc/2.0D0, dble(nbox)))
                end if
              end if
              nbox = nbox * abs(cos(xlat(ii,jj)*degrad)) + 1
              if (nbox < 2.0) then
                omt(ii,jj) = nearpoint(xx,yy,injx,iniy,imt,lwrap)
              else
                nbox = (nbox / 2) * 2
                if (itype == 4) then
                  omt(ii,jj) = mostaround(xx,yy,injx,iniy,imt,nbox, &
                                          ibnty,h2opct,lwrap)
                else
                  omt(ii,jj) = pctaround(xx,yy,injx,iniy,imt, &
                                         nbox,ival,lwrap)
                end if
              end if
          end select
        end do
      end do

      end subroutine interp

      subroutine filter1plakes(iy,jx,omt)
        implicit none
        integer , intent(in) :: iy , jx
        real(4) , intent(out) , dimension(iy, jx) :: omt
        integer , dimension(maxbins) :: cnt
        integer , dimension(9) :: around
        integer , parameter :: ilake = 14
        integer , parameter :: iocn = 15
        integer , parameter :: minlak = 2*ilake
        integer :: i , j , ii , jj , ip , jp , k , mpindex

        do i = 1 , iy
          do j = 1 , jx
            if (int(omt(i,j)) == ilake) then
              k = 1
              do ii = -1 , 1 , 1
                do jj = -1 , 1 , 1
                  ip = max(min(i+ii,iy),1)
                  jp = max(min(j+jj,jx),1)
                  around(k) = int(omt(ip,jp))
                  k = k + 1
                end do
              end do
              if (sum(around, around==ilake) .lt. minlak) then
                do k = 1 , maxbins
                  cnt(k) = sum(around/k,around==k)
                end do
                mpindex = 0
                do k = 1 , maxbins
                  if (k == ilake) cycle
                  if (k == iocn) cycle
                  if (cnt(k) > mpindex) mpindex = k
                end do
                omt(i,j) = mpindex
              end if
            end if
          end do
        end do
      end subroutine filter1plakes
!
      end module mod_interp
