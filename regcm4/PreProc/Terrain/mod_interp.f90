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

      implicit none

      contains

      function oned(x,a,b,c,d)
      implicit none
!
! Dummy arguments
!
      real(8) :: a , b , c , d , x
      real(8) :: oned
      intent (in) a , b , c , d , x
!
      oned = 0.
      if ( x==0. ) oned = b
      if ( x==1. ) oned = c
      if ( b*c==0. ) return
      if ( a*d==0. ) then
        oned = b*(1.0-x) + c*x
        if ( a/=0.0 ) oned = b + x*(0.5*(c-a)+x*(0.5*(c+a)-b))
        if ( d/=0.0 ) oned = c + (1.0-x)                                &
                           & *(0.5*(b-d)+(1.0-x)*(0.5*(b+d)-c))
        return
      end if
      oned = (1.0-x)*(b+x*(0.5*(c-a)+x*(0.5*(c+a)-b)))                  &
           & + x*(c+(1.0-x)*(0.5*(b-d)+(1.0-x)*(0.5*(b+d)-c)))
      end function oned

      function bint(xx,yy,list,iii,jjj,flag)

      implicit none
!
! Dummy arguments
!
      logical :: flag
      integer :: iii , jjj
      real(8) :: xx , yy
      real(8) :: bint
      real(8) , dimension(iii,jjj) :: list
      intent (in) flag , iii , jjj , list , xx , yy
!
! Local variables
!
      real(8) :: a , b , c , d , e , f , g , h , x , y
      integer :: i , j , k , kk , knear , l , ll , lnear , n
      real(8) , dimension(4,4) :: stl
!
!-----bilinear interpolation among four grid values
!
      bint = 0.0
      n = 0
      i = int(xx)
      j = int(yy)
      x = xx - i
      y = yy - j
      if ( abs(x)>0.001 .and. abs(y)>0.001 ) then
        do k = 1 , 4
          kk = i + k - 2
          do l = 1 , 4
            stl(k,l) = 0.
            if ( .not.(flag .and. (l==1)) ) then
              if ( .not.(flag .and. (l==4)) ) then
                if ( .not.(flag .and. (k==1)) ) then
                  if ( .not.(flag .and. (k==4)) ) then
                    ll = j + l - 2
                    if ( kk>=1 .and. kk<=iii ) then
                      if ( ll<=jjj .and. ll>=1 ) then
                        stl(k,l) = list(kk,ll)
                        n = n + 1
                        if ( stl(k,l)<=0.0 ) stl(k,l) = -1.E-20
                      end if
                    end if
                  end if
                end if
              end if
            end if
          end do
        end do
!
!-----find index of closest point to xx,yy.
!
        knear = float(2) + x + 0.5
        lnear = float(2) + y + 0.5
        a = oned(x,stl(1,1),stl(2,1),stl(3,1),stl(4,1))
        b = oned(x,stl(1,2),stl(2,2),stl(3,2),stl(4,2))
        c = oned(x,stl(1,3),stl(2,3),stl(3,3),stl(4,3))
        d = oned(x,stl(1,4),stl(2,4),stl(3,4),stl(4,4))
        bint = oned(y,a,b,c,d)
!
!--------if closest point is ocean, automatically reset terrain to
!--------preserve coastline.
!
        if ( .not.flag .and. stl(knear,lnear)<=0.001 ) bint = -0.00001
        if ( n==16 ) return
        if ( flag .and. n==4 ) return
        e = oned(y,stl(1,1),stl(1,2),stl(1,3),stl(1,4))
        f = oned(y,stl(2,1),stl(2,2),stl(2,3),stl(2,4))
        g = oned(y,stl(3,1),stl(3,2),stl(3,3),stl(3,4))
        h = oned(y,stl(4,1),stl(4,2),stl(4,3),stl(4,4))
        bint = (bint+oned(x,e,f,g,h))/2.
        if ( .not.flag .and. stl(knear,lnear)<=0.001 ) bint = -0.00001
        return
      end if

      bint = list(i,j)

      end function bint

      subroutine interp(iy, jx, xlat , xlon , htg, htsdg, ntypec)
 
      use mod_block
      implicit none
!
! Dummy arguments
!
      integer :: iy , jx , ntypec
      real(4) , dimension(iy, jx) :: xlat , xlon , htg , htsdg
      intent(in) iy , jx , ntypec , xlat , xlon
      intent(out) htg , htsdg
!
! Local variables
!
      logical :: flag
      real(8) :: xx , yy , xd , yd , wc , wm , rinc , dsgrid
      integer :: ii , jj , iindex , jindex , mxi , mxj
      real(8) , dimension(:,:) , allocatable :: xin1 , xin2
      logical , dimension(:,:) , allocatable :: lmask
!
      rinc = nnc
      if (lcrosstime) then
        mxj = nint((mod((grdlnma+360.0),360.0)-grdlnmn)*rinc) + 1
      else
        mxj = nint((grdlnma-grdlnmn)*rinc) + 1
      end if
      mxi = nint((grdltma-grdltmn)*rinc) + 1

      if (lonwrap) then
        mxj = mxj + 4
      end if

      print *, 'Allocating 2x', mxi, 'x', mxj

      allocate(xin1(mxi,mxj))
      allocate(xin2(mxi,mxj))
      allocate(lmask(mxi,mxj))

      xin1 = 0.0
      xin2 = 0.0
      lmask = .false.
      flag = .false.
      dsgrid = float(ntypec)/60.

      do ii = 1 , nobs
        if (lcrosstime) then
          jindex = nint((mod((xobs(ii)+360.0),360.0)-grdlnmn)*rinc) + 1
        else
          jindex = nint((xobs(ii)-grdlnmn)*rinc) + 1
        end if
        iindex = nint((yobs(ii)-grdltmn)*rinc) + 1
        if ( iindex < 1 .or. iindex>mxi .or. &
             jindex < 1 .or. jindex>mxj ) then
          print 99001 , ii , xobs(ii) , nnc , yobs(ii) , iindex , jindex
          stop 400
        end if
        xin1(iindex,jindex) = ht(ii)/100.
        if (ht2(ii) <= 400000000.0) xin2(iindex,jindex) = ht2(ii)/100000.
        lmask(iindex,jindex) = .true.
      end do
!
!       Fill the matrix using nearest point to missing one.
!
      do ii = 1 , mxi-1
        do jj = 1 , mxj-1
          if (.not. lmask(ii,jj)) then
            yy = ((ii*dsgrid)-grdltmn)*rinc + 1.0D+00
            xx = ((jj*dsgrid)-grdlnmn)*rinc + 1.0D+00
            wc = 999.0
            wm = 999.0
            if (ii > 1) then
              if (lmask(ii-1,jj)) then
                yd = (((ii-1)*dsgrid)-grdltmn)*rinc + 1.0D+00
                wc = (yy-yd)
                xin1(ii,jj) = xin1(ii-1,jj)
                xin2(ii,jj) = xin2(ii-1,jj)
                wm = wc
              end if
            end if
            if (lmask(ii+1,jj)) then
              yd = (((ii+1)*dsgrid)-grdltmn)*rinc + 1.0D+00
              wc = (yy-yd)
              if (wc < wm) then
                xin1(ii,jj) = xin1(ii+1,jj)
                xin2(ii,jj) = xin2(ii+1,jj)
                wm = wc
              end if
            end if
            if (jj > 1) then
              if (lmask(jj > 1 .and. ii,jj-1)) then
                xd = (((jj-1)*dsgrid)-grdltmn)*rinc + 1.0D+00
                wc = (xx-xd)
                if (wc < wm) then
                  xin1(ii,jj) = xin1(ii,jj-1)
                  xin2(ii,jj) = xin2(ii,jj-1)
                  wm = wc
                end if
              end if
            end if
            if (lmask(ii,jj+1)) then
              xd = (((jj+1)*dsgrid)-grdltmn)*rinc + 1.0D+00
              wc = (xx-xd)
              if (wc < wm) then
                xin1(ii,jj) = xin1(ii,jj+1)
                xin2(ii,jj) = xin2(ii,jj+1)
                wm = wc
              end if
            end if
          end if
        end do
      end do
      deallocate(lmask)
 
      if (lonwrap) then
        xin1 = cshift(xin1,shift=-2,dim=2)
        xin1(:,1) = xin1(:,mxj-3)
        xin1(:,2) = xin1(:,mxj-4)
        xin1(:,mxj-2) = xin1(:,3)
        xin1(:,mxj-1) = xin1(:,4)
        xin2 = cshift(xin2,shift=-2,dim=2)
        xin2(:,1) = xin2(:,mxj-3)
        xin2(:,2) = xin2(:,mxj-4)
        xin2(:,mxj-2) = xin2(:,3)
        xin2(:,mxj-1) = xin2(:,4)
      end if

      do ii = 1 , iy
        do jj = 1 , jx
 
          yy = (xlat(ii,jj)-grdltmn)*rinc + 1.0D+00
          if (lcrosstime) then
            xx = (mod((xlon(ii,jj)+360.0),360.0)-grdlnmn)*rinc + &
                      1.0D+00
          else
            xx = (xlon(ii,jj)-grdlnmn)*rinc + 1.0D+00
          end if
 
!         yy and xx are the exact index values of a point i,j of the
!         mesoscale mesh when projected onto an earth-grid of lat_s
!         and lon_s for which terrain observations are available.  it
!         is assumed that the earth grid has equal spacing in both
!         latitude and longitude.
 
          htsdg(ii,jj) = bint(yy,xx,xin2,mxi,mxj,flag)
          htg(ii,jj) = bint(yy,xx,xin1,mxi,mxj,flag)
        end do
      end do

      deallocate(xin1)
      deallocate(xin2)
 
99001 format (1x,'ii = ',i6,' xobs(ii) = ',f10.4,' incr = ',i3,         &
             &'yobs(ii) = ',f10.4,' iindex = ',i10,'  jindex = ',i10)
 
      end subroutine interp

      subroutine anal2(iy, jx, htg, htsdg)

      use mod_block

      implicit none
!
! Dummy arguments
!
      integer :: iy , jx
      real(4) , dimension(iy,jx) :: htg , htsdg
      intent (in) iy , jx
      intent (out) htg , htsdg
!
! Local variables
!
      real(4) , dimension(2,iy,jx) :: cor
      real(4) , dimension(2,iy,jx) :: htsav
      real(4) , dimension(2,iy,jx) :: xsum
      real(4) , dimension(2,iy,jx) :: wtmax
      integer , dimension(2,iy,jx) :: ns
      real(4) :: riobs , ris , ris2 , rjobs , rsq , rx , ry , wt , x ,  &
               & xmaxj , xminj , y , ymaxi , ymini , h1 , h2
      integer :: i , j , kk , maxi , maxj , mini , minj
!
!     objective analysis to fill a grid based on observations
!     xobs and yobs are x and y positions on observations, not
!     necessarily grid points.
!
!-----grid lengths in x and y directions are unity.
!-----rin is radius of influence in grid units
!
      ris = sqrt(rin**2)
      ris2 = 2.*ris
!
      wtmax = 0.0
      ns = 0
      cor = 0.0
      xsum = 0.0
      htsav = 0.0
!
!-----begin to process the nobs observations
!
      do kk = 1 , nobs
!
!-----convert obs. location to grid increment values of i and j
!
        riobs = yobs(kk)
        rjobs = xobs(kk)
!
        ymaxi = riobs + rin
        maxi = int(ymaxi+0.99)
        maxi = min(maxi,iy)
!
        ymini = riobs - rin
        mini = int(ymini)
        mini = max(mini,1)
!
        xmaxj = rjobs + rin
        maxj = int(xmaxj+0.99)
        maxj = min(maxj,jx)
!
        xminj = rjobs - rin
        minj = int(xminj)
        minj = max(minj,1)
!
        do i = mini , maxi
          do j = minj , maxj
            x = j
            y = i
!-----compute distance of k_th station from i,jth grid point
            rx = x - rjobs
            ry = y - riobs
            rsq = (rx**2 + ry**2)
            if ( rsq<ris2 ) then
              wt = (ris-rsq)/(ris+rsq)
!
!-----save      max. weighting factor and terrain height to check if
!-----point     grid should be treated as a land or sea point.
!
              if ( wt>0.0 ) then

                h1 = ht(kk)/100.0
                h2 = ht2(kk)/100000.0

                wtmax(1,i,j) = max(wt,wtmax(1,i,j))
                ns(1,i,j) = ns(1,i,j) + 1
                cor(1,i,j) = cor(1,i,j) + wt*h1
                xsum(1,i,j) = xsum(1,i,j) + wt
                if ( abs(wt-wtmax(1,i,j))<0.00001 ) then
                  htsav(1,i,j) = h1
                end if
                if (h2 <= 400.0) then
                  wtmax(2,i,j) = max(wt,wtmax(2,i,j))
                  cor(2,i,j) = cor(2,i,j) + wt*h2
                  ns(2,i,j) = ns(2,i,j) + 1
                  xsum(2,i,j) = xsum(2,i,j) + wt
                  if ( abs(wt-wtmax(2,i,j))<0.00001 ) then
                    htsav(2,i,j) = h2
                  end if
                end if

              end if
            end if
          end do
        end do
      end do
!
!-----now apply summed weights and weighted observations to determine
!-----terrain value at i,j points
!--------if closest observation to i,j is ocean, override to
!--------preserve the coastline.
!
      do i = 1 , iy
        do j = 1 , jx
          if ( ns(1,i,j)/=0 ) then
            cor(1,i,j) = cor(1,i,j)/xsum(1,i,j)
            htg(i,j) = cor(1,i,j)
            if (htsav(1,i,j) .le. 0.001) htg(i,j) = htsav(1,i,j)
          end if
          if ( ns(2,i,j)/=0 ) then
            cor(2,i,j) = cor(2,i,j)/xsum(2,i,j)
            htsdg(i,j) = cor(2,i,j)
            if (htsav(2,i,j) .le. 0.001) htsdg(i,j) = htsav(2,i,j)
          end if
        end do
      end do

      end subroutine anal2

      subroutine surf(xlat,xlon,lnduse,iy,jx,incr,dsgrid,lndout,land,   &
                    & nrec,h2opct,nveg,aertyp,intext,texout,frac_tex,   &
                    & ntex)
      use mod_block
      implicit none
!
! Dummy arguments
!
      character(7) :: aertyp
      real(4) :: dsgrid , h2opct
      integer :: incr , iy , jx , nrec , ntex , nveg
      real(4) , dimension(iy,jx) :: lndout , texout , xlat , xlon
      real(4) , dimension(iy,jx,ntex) :: frac_tex
      integer , dimension(iy,jx) :: intext , lnduse
      real(4) , dimension(iy,jx,2) :: land
      intent (in) aertyp , dsgrid , h2opct , incr , &
                & iy , jx , nrec , ntex , nveg , xlat
      intent (out) frac_tex , intext , lnduse
      intent (inout) land , lndout , texout , xlon
!
! Local variables
!
      logical :: flag
      integer :: i , ii , iindex , ilev , j , jindex , jj , lrec ,  &
               & lengdo , nbase , mxi , mxj
      real(4) , dimension(iy,jx,2) :: itex
      real(8) :: xx , yy , rinc , xd , yd , wc , wm
      logical , dimension(:,:) , allocatable :: lmask
      real(8) , dimension(:,:) , allocatable :: lnd8
!
      flag = .true.
      lrec = 0
      land = 0.0
      itex = 0.0
      rinc = incr
!
      print *, 'Input data point MIN is at ', grdltmn , grdlnmn
      print *, 'Input data point MAX is at ', grdltma , grdlnma
      print *, 'Input data resolution is   ', dsgrid
      if (lcrosstime) then
        mxj = nint((mod((grdlnma+360.0),360.0)-grdlnmn)*rinc) + 1
      else
        mxj = nint((grdlnma-grdlnmn)*rinc) + 1
      end if
      mxi = nint((grdltma-grdltmn)*rinc) + 1

      if ( aertyp(7:7)=='1' ) then
        lengdo = nveg + ntex
      else
        lengdo = nveg
      end if

      if (lonwrap) then
        mxj = mxj + 4
      end if

      print *, 'Allocating ', mxi, 'x', mxj

      allocate(lmask(mxi,mxj))
      allocate(lnd8(mxi,mxj))

      do ilev = 1 , lengdo
        lmask(:,:) = .false.
        rewind (48)
        do lrec = 1 , nrec
          read (48) stores
          if (lcrosstime) then
            jindex = nint((mod((stores(2)+360.0),360.0)-grdlnmn)*       &
                           rinc) + 1
          else
            jindex = nint((stores(2)-grdlnmn)*rinc) + 1
          end if
          iindex = nint((stores(1)-grdltmn)*rinc) + 1
          if ( iindex < 1 .or. iindex>mxi .or. &
               jindex < 1 .or. jindex>mxj ) then
            print 99001 , iindex , jindex , lrec , stores(1) ,          &
                & stores(2) , grdltmn , grdlnmn
            stop 60
          end if
          lnd8(iindex,jindex) = stores(ilev+4)
          lmask(iindex,jindex) = .true.
        end do
!
!       Fill the matrix using nearest point to missing one.
!
        do ii = 1 , mxi-1
          do jj = 1 , mxj-1
            if (.not. lmask(ii,jj)) then
              yy = ((ii*dsgrid)-grdltmn)*rinc + 1.0D+00
              xx = ((jj*dsgrid)-grdlnmn)*rinc + 1.0D+00
              wc = 999.0
              wm = 999.0
              if (ii > 1) then
                if (lmask(ii-1,jj)) then
                  yd = (((ii-1)*dsgrid)-grdltmn)*rinc + 1.0D+00
                  wc = (yy-yd)
                  lnd8(ii,jj) = lnd8(ii-1,jj)
                  wm = wc
                end if
              end if
              if (lmask(ii+1,jj)) then
                yd = (((ii+1)*dsgrid)-grdltmn)*rinc + 1.0D+00
                wc = (yy-yd)
                if (wc < wm) then
                  lnd8(ii,jj) = lnd8(ii+1,jj)
                  wm = wc
                end if
              end if
              if (jj > 1) then
                if (lmask(jj > 1 .and. ii,jj-1)) then
                  xd = (((jj-1)*dsgrid)-grdltmn)*rinc + 1.0D+00
                  wc = (xx-xd)
                  if (wc < wm) then
                    lnd8(ii,jj) = lnd8(ii,jj-1)
                    wm = wc
                  end if
                end if
              end if
              if (lmask(ii,jj+1)) then
                xd = (((jj+1)*dsgrid)-grdltmn)*rinc + 1.0D+00
                wc = (xx-xd)
                if (wc < wm) then
                  lnd8(ii,jj) = lnd8(ii,jj+1)
                  wm = wc
                end if
              end if
            end if
          end do
        end do
!
!       Special case for wrapping of longitude
!
        if (lonwrap) then
          lnd8 = cshift(lnd8,shift=-2,dim=2)
          lnd8(:,1) = lnd8(:,mxj-3)
          lnd8(:,2) = lnd8(:,mxj-4)
          lnd8(:,mxj-2) = lnd8(:,3)
          lnd8(:,mxj-1) = lnd8(:,4)
        end if

        if ( ilev<=nveg ) then
          do ii = 1 , iy
            do jj = 1 , jx
              yy = (xlat(ii,jj)-grdltmn)*rinc + 1.0D+00
              if (lcrosstime) then
                xx = (mod((xlon(ii,jj)+360.0),360.0)-grdlnmn)*rinc + &
                      1.0D+00
              else
                xx = (xlon(ii,jj)-grdlnmn)*rinc + 1.0D+00
              end if
              lndout(ii,jj) = bint(yy,xx,lnd8,mxi,mxj,flag)
!
!             note: it is desirable to force grid boxes with less
!             than 75 percent water to be a land category,
!             even if water is the largest single category.
!
              if ( .not.((ilev==14 .or. ilev==15) .and.                 &
                     &    lndout(ii,jj)<h2opct) ) then
                if ( ilev/=25 .or. lndout(ii,jj) >=h2opct ) then
                  if ( lndout(ii,jj)>land(ii,jj,1) ) then
                    land(ii,jj,1) = lndout(ii,jj)
                    land(ii,jj,2) = ilev
                  end if
                end if
              end if
            end do
          end do
        else if ( ilev>nveg .and. aertyp(7:7) == '1' ) then
          nbase = nveg
          do ii = 1 , iy
            do jj = 1 , jx
              yy = (xlat(ii,jj)-grdltmn)*rinc + 1.0D+00
              if (lcrosstime) then
                xx = (mod((xlon(ii,jj)+360.0),360.0)-grdlnmn)*rinc + &
                     1.0D+00
              else
                xx = (xlon(ii,jj)-grdlnmn)*rinc + 1.0D+00
              end if
              texout(ii,jj) = bint(yy,xx,lnd8,mxi,mxj,flag)
              frac_tex(ii,jj,ilev-nbase) = texout(ii,jj)
!
!             note: it is desirable to force grid boxes with less
!             than 75 percent water to be a land category,
!             even if water is the largest single category.
!
              if ( ilev/=nveg+14 .or. texout(ii,jj) >=h2opct ) then
                if ( ilev/=nveg+18 .or. texout(ii,jj)>=h2opct ) then
                  if ( texout(ii,jj)>itex(ii,jj,1) ) then
                    itex(ii,jj,1) = texout(ii,jj)
                    itex(ii,jj,2) = ilev - nbase
                  end if
                end if
              end if
            end do
          end do
        end if
      end do
!
      deallocate(lmask)
      deallocate(lnd8)
!
      do i = 1 , iy
        do j = 1 , jx
          lndout(i,j) = land(i,j,2)
          lnduse(i,j) = int(land(i,j,2))
          if ( aertyp(7:7)=='1' ) then
            texout(i,j) = itex(i,j,2)
            intext(i,j) = int(itex(i,j,2))
          end if
        end do
      end do
!
99001 format (1x,'*** iindex = ',i4,'   jindex = ',i4,'   lrec = ',i5,  &
             &'   lat = ',f10.3,3x,'lon = ',f10.3,3x,'grdltmn = ',f10.3,&
            & 5x,'grdlnmn = ',f10.3)
!
      end subroutine surf

      end module mod_interp
