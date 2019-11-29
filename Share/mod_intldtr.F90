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

module mod_intldtr

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_stdio
  use mod_message
  use mod_earth

  private

  public :: mxmnll , interp , filter1plakes

  real(rk8) , public :: grdlnmn , grdltmn , grdlnma , grdltma
  real(rk8) , public :: xmaxlat , xmaxlon , xminlat , xminlon
  integer(ik4) , public :: nlatin , nlonin
  logical , public :: lonwrap , lcrosstime

  real(rk8) , dimension(4,4) :: c
  real(rk8) , dimension(16,16) :: wt
  integer(ik4) , parameter :: maxbins = 20
  integer(ik4) , dimension(maxbins) :: bincnt
  real(rk8) , dimension(maxbins) :: bmindist
  logical , dimension(2,maxbins) :: lndwt

  data wt/1.0_rk8 , 0.0_rk8 , -3.0_rk8 , 2.0_rk8 , 4*0.0_rk8 , -3.0_rk8 ,   &
          0.0_rk8 , 9.0_rk8 , -6.0_rk8 , 2.0_rk8 , 0.0_rk8 , -6.0_rk8 ,     &
          4.0_rk8 , 8*0.0_rk8 , 3.0_rk8 , 0.0_rk8 , -9.0_rk8 , 6.0_rk8 ,    &
         -2.0_rk8 , 0.0_rk8 , 6.0_rk8 , -4.0_rk8 , 10*0.0_rk8 , 9.0_rk8 ,   &
         -6.0_rk8 , 2*0.0_rk8 , -6.0_rk8 , 4.0_rk8 , 2*0.0_rk8 , 3.0_rk8 ,  &
         -2.0_rk8 , 6*0.0_rk8 , -9.0_rk8 , 6.0_rk8 , 2*0.0_rk8 , 6.0_rk8 ,  &
         -4.0_rk8 , 4*0.0_rk8 , 1.0_rk8 , 0.0_rk8 , -3.0_rk8 , 2.0_rk8 ,    &
         -2.0_rk8 , 0.0_rk8 , 6.0_rk8 , -4.0_rk8 , 1.0_rk8 , 0.0_rk8 ,      &
         -3.0_rk8 , 2.0_rk8 , 8*0.0_rk8 , -1.0_rk8 , 0.0_rk8 , 3.0_rk8 ,    &
         -2.0_rk8 , 1.0_rk8 , 0.0_rk8 , -3.0_rk8 , 2.0_rk8 , 10*0.0_rk8 ,   &
         -3.0_rk8 , 2.0_rk8 , 2*0.0_rk8 , 3.0_rk8 , -2.0_rk8 , 6*0.0_rk8 ,  &
          3.0_rk8 , -2.0_rk8 , 2*0.0_rk8 , -6.0_rk8 , 4.0_rk8 , 2*0.0_rk8 , &
          3.0_rk8 , -2.0_rk8 , 0.0_rk8 , 1.0_rk8 , -2.0_rk8 , 1.0_rk8 ,     &
          5*0.0_rk8 , -3.0_rk8 , 6.0_rk8 , -3.0_rk8 , 0.0_rk8 , 2.0_rk8 ,   &
         -4.0_rk8 , 2.0_rk8 , 9*0.0_rk8 , 3.0_rk8 , -6.0_rk8 , 3.0_rk8 ,    &
          0.0_rk8 , -2.0_rk8 , 4.0_rk8 , -2.0_rk8 , 10*0.0_rk8 , -3.0_rk8 , &
          3.0_rk8 , 2*0.0_rk8 , 2.0_rk8 , -2.0_rk8 , 2*0.0_rk8 , -1.0_rk8 , &
          1.0_rk8 , 6*0.0_rk8 , 3.0_rk8 , -3.0_rk8 , 2*0.0_rk8 , -2.0_rk8 , &
          2.0_rk8 , 5*0.0_rk8 , 1.0_rk8 , -2.0_rk8 , 1.0_rk8 , 0.0_rk8 ,    &
         -2.0_rk8 , 4.0_rk8 , -2.0_rk8 , 0.0_rk8 , 1.0_rk8 , -2.0_rk8 ,     &
          1.0_rk8 , 9*0.0_rk8 , -1.0_rk8 , 2.0_rk8 , -1.0_rk8 , 0.0_rk8 ,   &
          1.0_rk8 , -2.0_rk8 , 1.0_rk8 , 10*0.0_rk8 , 1.0_rk8 , -1.0_rk8 ,  &
          2*0.0_rk8 , -1.0_rk8 , 1.0_rk8 , 6*0.0_rk8 , -1.0_rk8 , 1.0_rk8 , &
          2*0.0_rk8 , 2.0_rk8 , -2.0_rk8 , 2*0.0_rk8 , -1.0_rk8 , 1.0_rk8/

  data lndwt /26*.false.,.true.,.true.,.true.,11*.false./

  contains

  integer(ik4) function inear(x,m)
    implicit none
    real(rk8) , intent(in) :: x
    integer(ik4) , intent(in) :: m
    if (.not. lonwrap) then
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

  integer(ik4) function jnear(y,n)
    implicit none
    real(rk8) , intent(in) :: y
    integer(ik4) , intent(in) :: n
    jnear = min(max(nint(y),1),n)
  end function jnear

  integer(ik4) function ifloor(x,m)
    implicit none
    real(rk8) , intent(in) :: x
    integer(ik4) , intent(in) :: m
    if (.not. lonwrap) then
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

  integer(ik4) function jfloor(y,n)
    implicit none
    real(rk8) , intent(in) :: y
    integer(ik4) , intent(in) :: n
    jfloor = min(max(floor(y),1),n)
  end function jfloor

  real(rkx) function nearpoint(x,y,m,n,grid)
    implicit none
    integer(ik4) , intent(in) :: m , n
    real(rk8) , intent(in) :: x , y
    real(rkx) , intent(in) , dimension(m,n) :: grid
    nearpoint = grid(inear(x,m),jnear(y,n))
  end function nearpoint

  real(rkx) function mostaround(x,y,m,n,grid,nbox,ibnty,h2opct)
    implicit none
    integer(ik4) , intent(in) :: m , n , nbox , ibnty
    real(rk8) , intent(in) :: x , y
    real(rkx) , intent(in) , dimension(m,n) :: grid
    real(rkx) , intent(in) :: h2opct

    real(rk8) , dimension(nbox*nbox) :: binval , bindist
    real(rk8) :: dist , rx , ry , wtp
    integer(ik4) :: ii0 , jj0 , ii , jj
    integer(ik4) :: totpoints , i , j , lastc , hbox

    hbox = nbox / 2
    totpoints = nbox*nbox
    ii0 = ifloor(x,m)-hbox
    jj0 = jfloor(y,n)-hbox
    do i = 1 , nbox
      do j = 1 , nbox
        rx = real(ii0 + i - 1,rk8)
        ry = real(jj0 + j - 1,rk8)
        ii = ifloor(rx,m)
        jj = jfloor(ry,n)
        binval((i-1)*nbox+j) = grid(ii,jj)
        bindist((i-1)*nbox+j) = sqrt((x-rx)**2+(y-ry)**2)
      end do
    end do
    bincnt = 0
    bmindist = d_two*real(nbox,rk8)
    do i = 1 , totpoints
      bincnt(int(binval(i))) = bincnt(int(binval(i))) + 1
      if (bindist(i) < bmindist(int(binval(i)))) &
        bmindist(int(binval(i))) = bindist(i)
    end do
    ! Set point to land if less than fixed percent of water
    wtp = (real(sum(bincnt,mask=lndwt(ibnty,:)),rk8)/real(totpoints,rk8))*d_100
    if (wtp > d_zero .and. wtp < 100.0_rk8-h2opct) then
      if (ibnty == 1) then
        bincnt(14) = 0
        bincnt(15) = 0
      else
        bincnt(15) = 0
      end if
    end if
    mostaround = -d_one
    lastc = -1
    dist = maxval(bmindist) + d_one
    do i = 1 , maxbins
      if (bincnt(i) > 0) then
        if (bincnt(i) > lastc) then
          mostaround = real(i,rkx)
          dist = bmindist(i)
          lastc = bincnt(i)
        else if (bincnt(i) == lastc) then
          if (bmindist(i) < dist) then
            mostaround = real(i,rkx)
            dist = bmindist(i)
            lastc = bincnt(i)
          end if
        end if
      end if
    end do
  end function mostaround

  real(rkx) function pctaround(x,y,m,n,grid,nbox,ival)
    implicit none
    integer(ik4) , intent(in) :: m , n , ival , nbox
    real(rk8) , intent(in) :: x , y
    real(rkx) , intent(in) , dimension(m,n) :: grid
    integer(ik4) :: ii0 , jj0 , ii , jj
    integer(ik4) :: i , j
    real(rk8) :: rx , ry , pc , npc

    pc = d_zero
    npc = real(nbox*nbox,rk8)
    ii0 = ifloor(x,m)
    jj0 = jfloor(y,n)
    do i = 1 , nbox
      do j = 1 , nbox
        rx = real(ii0 + i - nbox/2,rk8)
        ry = real(jj0 + j - nbox/2,rk8)
        ii = ifloor(rx,m)
        jj = jfloor(ry,n)
        if (int(grid(ii,jj)) == ival) then
          pc = pc + 1
        end if
      end do
    end do
    pctaround = real(pc/npc*100_rk8,rkx)
  end function pctaround

  real(rkx) function bilinear(x,y,m,n,grid)
    implicit none
    integer(ik4) , intent(in) :: m , n
    real(rk8) , intent(in) :: x , y
    real(rkx) , intent(in) , dimension(m,n) :: grid

    real(rk8) :: dx, dy, p12, p03
    real(rk8) :: ii0, jj0, ii1, jj1, ii2, jj2, ii3, jj3
    integer(ik4) :: i0, j0, i1, j1, i2, j2, i3, j3

    !-----bilinear interpolation among four grid values

    if (.not. lonwrap) then
      ii0 = real(min(max(floor(x),1),m),rk8)
      ii2 = real(min(max(ceiling(x),1),m),rk8)
      dx = (x-ii0)
    else
      ii0 = floor(x)
      ii2 = ceiling(x)
      dx = (x-ii0)
      if (ii0 < 1) then
        ii0 = real(m,rk8)
      end if
      if (ii2 > m) then
        ii2 = d_one
      end if
    end if
    ii1 = ii0
    ii3 = ii2
    jj0 = real(min(max(floor(y),1),n),rk8)
    dy = (y-jj0)
    jj3 = jj0
    jj1 = real(min(max(ceiling(y),1),n),rk8)
    jj2 = jj1

    i0 = int(ii0)
    j0 = int(jj0)
    i1 = int(ii1)
    j1 = int(jj1)
    i2 = int(ii2)
    j2 = int(jj2)
    i3 = int(ii3)
    j3 = int(jj3)
    p12 = dx*grid(i2,j2)+(d_one-dx)*grid(i1,j1)
    p03 = dx*grid(i3,j3)+(d_one-dx)*grid(i0,j0)
    bilinear = real(dy*p12+(d_one-dy)*p03,rkx)
  end function bilinear

  real(rkx) function dwgt(x,y,m,n,grid,nb)
    implicit none
    integer(ik4) , intent(in) :: m , n , nb
    real(rk8) , intent(in) :: x , y
    real(rkx) , intent(in) , dimension(m,n) :: grid

    real(rk8) :: dx , pp , f
    real(rk8) :: ii0, jj0
    integer(ik4) :: i , j , ii , jj , iii , jjj , nh

    !-----distance weight with radius of influence

    nh = nb / 2
    if (.not. lonwrap) then
      ii0 = real(min(max(floor(x),1),m),rk8)
      jj0 = real(min(max(floor(y),1),n),rk8)
      ii = int(ii0)
      jj = int(jj0)
      pp = 0.0_rk8
      f = 0.0_rk8
      do j = jj - nh , jj + nh
        do i = ii - nh , ii + nh
          iii = max(1,min(m,i))
          jjj = max(1,min(n,j))
          dx = max(1.0e-10_rk8,sqrt((x-real(i,rk8))**2+(y-real(j,rk8))**2))
          pp = pp + d_one/dx
          f = f + grid(iii,jjj)/dx
        end do
      end do
    else
      ii0 = floor(x)
      ii = int(ii0)
      if ( ii < 1 ) ii = ii + m
      if ( ii > m ) ii = ii - m
      jj0 = real(min(max(floor(y),1),n),rk8)
      jj = int(jj0)
      pp = 0.0_rk8
      f = 0.0_rk8
      do j = jj - nh , jj + nh
        do i = ii - nh , ii + nh
          iii = i
          if ( i < 1 ) then
            iii = iii + m
          else if ( i > m ) then
            iii = iii - m
          end if
          jjj = max(1,min(n,j))
          dx = max(1.0e-10_rk8,sqrt((x-real(i,rk8))**2+(y-real(j,rk8))**2))
          pp = pp + d_one/dx
          f = f + grid(iii,jjj)/dx
        end do
      end do
    end if
    dwgt = real(f / pp,rkx)
  end function dwgt

  real(rkx) function bicubic(x,y,m,n,grid)
    implicit none
    integer(ik4) , intent(in) :: m , n
    real(rk8) , intent(in) :: x , y
    real(rkx) , intent(in) , dimension(m,n) :: grid
    real(rk8) , dimension(4) :: f , f1 , f12 , f2
    real(rk8) :: xl , xu , yl , yu
    integer(ik4) :: i , ii , j , mm , nn , im , imp1 , imn1

    mm = int(x)
    nn = int(y)
    if ( .not. lonwrap) then
      mm = max(2, min(m-2,mm))
    end if
    nn = max(2, min(n-2,nn))

    xl = real(mm,rk8)
    xu = real(mm + 1,rk8)
    yl = real(nn,rk8)
    yu = real(nn + 1,rk8)
    do j = nn , nn + 1
      do i = mm , mm + 1
        ii = 1 + (i-mm) + 3*(j-nn)
        if ( ii==5 ) ii = 3
        if (lonwrap) then
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
        f1(ii) = (grid(imp1,j)-grid(imn1,j))/(d_two)
        f2(ii) = (grid(im,j+1)-grid(im,j-1))/(d_two)
        f12(ii) = (grid(imp1,j+1)-grid(imp1,j-1)-&
                   grid(imn1,j+1)+grid(imn1,j-1))/(d_four)
      end do
    end do
    call bcuint(f,f1,f2,f12,xl,xu,yl,yu,x,y,bicubic)
  end function bicubic

  subroutine bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,a)
    implicit none
    real(rk8) , intent(in) :: x1 , x1l , x1u , x2 , x2l , x2u
    real(rkx) , intent(out) :: a
    real(rk8) , intent(in) , dimension(4) :: y , y1 , y12 , y2
    integer(ik4) :: i
    real(rk8) :: t , u

    call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l)
    t = (x1-x1l)/(x1u-x1l)
    u = (x2-x2l)/(x2u-x2l)
    a = d_zero
    do i = 4 , 1 , -1
      a = real(t*a + ((c(i,4)*u+c(i,3))*u+c(i,2))*u + c(i,1),rkx)
    end do
  end subroutine bcuint

  subroutine bcucof(y,y1,y2,y12,d1,d2)
    implicit none
    real(rk8) , intent(in) :: d1 , d2
    real(rk8) , intent(in) , dimension(4) :: y , y1 , y12 , y2
    real(rk8) , dimension(16) :: cl , x
    real(rk8) :: d1d2 , xx
    integer(ik4) :: i , j , k , l
    d1d2 = d1*d2
    do i = 1 , 4
      x(i) = y(i)
      x(i+4) = y1(i)*d1
      x(i+8) = y2(i)*d2
      x(i+12) = y12(i)*d1d2
    end do
    do i = 1 , 16
      xx = d_zero
      do k = 1 , 16
        xx = xx + wt(i,k)*x(k)
      end do
      cl(i) = xx
    end do
    l = 0
    do i = 1 , 4
      do j = 1 , 4
        l = l + 1
        c(j,i) = cl(l)
      end do
    end do
  end subroutine bcucof
  !
  ! Interpolates input regolar lat/lon grid on output model grid
  !
  subroutine interp(ds,jx,iy,xlat,xlon,omt, &
                    imt,itype,ival,ibnty,h2opct,rdem)
    implicit none
    integer(ik4) , intent(in) :: iy , jx , itype
    real(rkx) , intent(in) , dimension(jx,iy) :: xlat , xlon
    real(rkx) , intent(out) , dimension(jx,iy) :: omt
    real(rkx) , intent(in) , dimension(nlonin,nlatin) :: imt
    real(rkx) , intent(in) :: ds
    integer(ik4) , intent(in) , optional :: ival
    integer(ik4) , intent(in) , optional :: ibnty
    real(rkx) , intent(in) , optional :: h2opct , rdem

    integer(ik4) :: nbox , ii , jj
    integer , parameter :: nbmax = 100
    real(rk8) :: xx , yy , incx , incy , rincx , rincy , dd , drcm

    incy = (grdltma-grdltmn)/real(nlatin-1,rk8)
    if ( lcrosstime ) then
      if ( grdlnmn > grdlnma ) then
        incx = ((360.0_rk8+grdlnma)-grdlnmn)/real(nlonin-1,rk8)
      else
        incx = (grdlnma-grdlnmn)/real(nlonin-1,rk8)
      end if
    else if ( grdlnmn > grdlnma ) then
      incx = (grdlnma+(360.0_rk8-grdlnmn))/real(nlonin-1,rk8)
      grdlnmn = -(360.0_rk8-grdlnmn)
    else
      incx = (grdlnma-grdlnmn)/real(nlonin-1,rk8)
    end if
    rincx = d_one/incx
    rincy = d_one/incy
    if ( present(rdem) ) then
      drcm = rdem*ds*sqrt(2.0_rk8)
    end if

    ! yy and xx are the exact index values of a point j,i of the
    ! mesoscale mesh when projected onto an earth-grid of lat
    ! and lon for which terrain observations are available.

    select case (itype)
      case(1)
        do ii = 1 , iy
          do jj = 1 , jx
            yy = (xlat(jj,ii)-grdltmn)*rincy + d_one
            if (lcrosstime) then
              xx = (mod((xlon(jj,ii)+deg360),deg360)-grdlnmn) * &
                    rincx + d_one
            else
              xx = (xlon(jj,ii)-grdlnmn)*rincx + d_one
            end if
            omt(jj,ii) = bilinear(xx,yy,nlonin,nlatin,imt)
          end do
        end do
      case(2)
        do ii = 1 , iy
          do jj = 1 , jx
            yy = (xlat(jj,ii)-grdltmn)*rincy + d_one
            if (lcrosstime) then
              xx = (mod((xlon(jj,ii)+deg360),deg360)-grdlnmn) * &
                    rincx + d_one
            else
              xx = (xlon(jj,ii)-grdlnmn)*rincx + d_one
            end if
            omt(jj,ii) = bicubic(xx,yy,nlonin,nlatin,imt)
          end do
        end do
      case(3)
        do ii = 1 , iy
          do jj = 1 , jx
            yy = (xlat(jj,ii)-grdltmn)*rincy + d_one
            if (lcrosstime) then
              xx = (mod((xlon(jj,ii)+deg360),deg360)-grdlnmn) * &
                    rincx + d_one
            else
              xx = (xlon(jj,ii)-grdlnmn)*rincx + d_one
            end if
            omt(jj,ii) = nearpoint(xx,yy,nlonin,nlatin,imt)
          end do
        end do
      case(4,5,6)
        do ii = 1 , iy
          do jj = 1 , jx
            yy = (xlat(jj,ii)-grdltmn)*rincy + d_one
            if (lcrosstime) then
              xx = (mod((xlon(jj,ii)+deg360),deg360)-grdlnmn) * &
                    rincx + d_one
            else
              xx = (xlon(jj,ii)-grdlnmn)*rincx + d_one
            end if
            if (ii == 1 .or. ii == iy ) then
              omt(jj,ii) = nearpoint(xx,yy,nlonin,nlatin,imt)
              cycle
            end if
            if ( .not. lonwrap ) then
              if ( jj == 1 .or.  jj == jx ) then
                omt(jj,ii) = nearpoint(xx,yy,nlonin,nlatin,imt)
                cycle
              end if
            end if
            dd = gcdist_simple(real(grdltmn+incy*floor(yy),rkx),   &
                               real(grdlnmn+incx*floor(xx),rkx),   &
                               real(grdltmn+incy*ceiling(yy),rkx), &
                               real(grdlnmn+incx*ceiling(xx),rkx))
            nbox = min(max(nint(drcm/dd),2),nbmax)
            nbox = (nbox / 2) * 2
            if (itype == 4) then
              omt(jj,ii) = mostaround(xx,yy,nlonin,nlatin,imt,nbox,ibnty,h2opct)
            else if (itype == 5) then
              omt(jj,ii) = pctaround(xx,yy,nlonin,nlatin,imt,nbox,ival)
            else
              omt(jj,ii) = dwgt(xx,yy,nlonin,nlatin,imt,nbox)
            end if
          end do
        end do
      case default
        write(stderr,*) 'Unknown interpolation type'
        call die('interp')
    end select
  end subroutine interp

  subroutine filter1plakes(jx,iy,omt)
    implicit none
    integer(ik4) , intent(in) :: jx , iy
    real(rkx) , intent(inout) , dimension(jx,iy) :: omt
    integer(ik4) , dimension(maxbins) :: cnt
    integer(ik4) , dimension(9) :: around
    integer(ik4) , parameter :: ilake = 14
    integer(ik4) , parameter :: iocn = 15
    integer(ik4) , parameter :: minlak = 2*ilake
    integer(ik4) :: i , j , ii , jj , ip , jp , k , mpindex

    do i = 1 , iy
      do j = 1 , jx
        if (int(omt(j,i)) == ilake) then
          k = 1
          do ii = -1 , 1 , 1
            do jj = -1 , 1 , 1
              ip = max(min(i+ii,iy),1)
              jp = max(min(j+jj,jx),1)
              around(k) = int(omt(jp,ip))
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
            if (mpindex == 0) mpindex = iocn
            omt(j,i) = real(mpindex,rkx)
          end if
        end if
      end do
    end do
  end subroutine filter1plakes

  subroutine mxmnll(jx,iy,xlon,xlat,iband)
    implicit none
    integer(ik4) , intent(in) :: iy , jx , iband
    real(rkx) , dimension(jx,iy) , intent(in) :: xlat , xlon

    real(rk8) :: xtstlon1 , xtstlon2
    real(rk8) , parameter :: coord_eps = 0.001_rk8
    !
    ! PURPOSE : FINDS THE MAXIMUM AND MINIMUM LATITUDE AND LONGITUDE
    !
    xminlat = floor(minval(xlat))
    xmaxlat = ceiling(maxval(xlat))

    if ( iband.eq.1 ) then
      xminlon  = -deg180
      xmaxlon  =  deg180
      xtstlon1 = xminlon
      xtstlon2 = xmaxlon
    else if ( abs(xminlat+deg90) < coord_eps .or. &
              abs(xmaxlat-deg90) < coord_eps ) then
      xminlon  = -deg180
      xmaxlon  =  deg180
      xtstlon1 = xminlon
      xtstlon2 = xmaxlon
    else
      xminlon  = floor(minval(xlon(1,:)))
      xmaxlon  = ceiling(maxval(xlon(jx,:)))
      xtstlon1 = floor(maxval(xlon(1,:)))
      xtstlon2 = ceiling(minval(xlon(jx,:)))
    end if

    if ( abs(xminlon-xmaxlon) < coord_eps ) then
      xminlon  = -deg180
      xmaxlon  =  deg180
      xtstlon1 = xminlon
      xtstlon2 = xmaxlon
    end if

    lonwrap = .false.
    lcrosstime = .false.
    if ( (xmaxlon-xminlon) > (deg360-coord_eps) ) then
      lonwrap = .true.
      write(stdout,*) 'Special case for longitude wrapping'
    end if
    if ( abs(xminlon - xtstlon1) > deg180 .or.   &
         abs(xmaxlon - xtstlon2) > deg180 .or.   &
         xminlon > deg00 .and. xmaxlon < deg00 ) then
      lcrosstime = .true.
      if ( xminlon < deg00 .and. xtstlon1 > deg00 ) xminlon = xtstlon1
      if ( xmaxlon > deg00 .and. xtstlon2 < deg00 ) xmaxlon = xtstlon2
      write(stdout,*) 'Special case for timeline crossing'
    end if

    write(stdout,*) 'Calculated large extrema:'
    write(stdout,*) '         MINLAT = ', xminlat
    write(stdout,*) '         MAXLAT = ', xmaxlat
    write(stdout,*) '         MINLON = ', xminlon
    write(stdout,*) '         MAXLON = ', xmaxlon

  end subroutine mxmnll

end module mod_intldtr
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
