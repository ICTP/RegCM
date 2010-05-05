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

      subroutine interp(nx, ny, alat, alon, htg, htsdg, ntypec)
 
      use mod_block
      use mod_maps
      implicit none
!
! Dummy arguments
!
      integer :: nx , ny , ntypec
      real(4) , dimension(nx, ny) :: alat , alon , htg , htsdg
      intent(in) nx , ny , alat , ntypec
      intent(inout) alon
      intent(out) htg , htsdg
!
! Local variables
!
      real(4) :: dsgrid
      logical :: flag
      real(8) :: h1 , h2 , v21 , xx , yy
      integer :: i , ii , iindex , j , jindex
      real(8) , dimension(iter,jter) :: xin1 , xin2
!
      do i = 1 , iter
        do j = 1 , jter
          xin1(i,j) = 0.0
          xin2(i,j) = 0.0
        end do
      end do
      do ii = 1 , nobs
        jindex = (xobs(ii)-grdlnmn)*nnc + 1.1
        iindex = (yobs(ii)-grdltmn)*nnc + 1.1
        if ( iindex>iter .or. jindex>jter ) then
          print 99001 , ii , xobs(ii) , nnc , yobs(ii) , iindex , jindex
          stop 400
        end if
        h1 = max(ht(ii),0.0)
        xin1(iindex,jindex) = h1/100.
        h2 = max(ht2(ii),0.0)
        xin2(iindex,jindex) = h2/100000.
      end do
 
      flag = .false.
      dsgrid = float(ntypec)/60.
 
      do i = 1 , ny - 1
        do j = 1 , nx - 1
 
          yy = -(grdltmn-alat(i,j))/dsgrid + 1.0
          if ( grdlnmn<=-180.0 .and. alon(i,j)>0.0 ) alon(i,j)          &
             & = alon(i,j) - 360.
          xx = -(grdlnmn-alon(i,j))/dsgrid + 1.0
 
!         yy and xx are the exact index values of a point i,j of the
!         mesoscale mesh when projected onto an earth-grid of lat_s
!         and lon_s for which terrain observations are available.  it
!         is assumed that the earth grid has equal spacing in both
!         latitude and longitude.
 
          h1 = max(bint(yy,xx,xin1,iter,jter,flag),0.D0)*100.
          h2 = max(bint(yy,xx,xin2,iter,jter,flag),0.D0)*100000.
          htg(i,j) = h1
          v21 = h2 - h1**2
          htsdg(i,j) = sqrt(max(v21,0.D0))
 
        end do
      end do
 
99001 format (1x,'ii = ',i6,' xobs(ii) = ',f10.4,' incr = ',i3,         &
             &'yobs(ii) = ',f10.4,' iindex = ',i10,'  jindex = ',i10)
 
      end subroutine interp

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
      i = int(xx+0.00001)
      j = int(yy+0.00001)
      x = xx - i
      y = yy - j
      if ( abs(x)>0.0001 .or. abs(y)>0.0001 ) then
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
        go to 99999
      end if
      bint = list(i,j)
      return
99999 continue
      end function bint

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
        go to 99999
      end if
      oned = (1.0-x)*(b+x*(0.5*(c-a)+x*(0.5*(c+a)-b)))                  &
           & + x*(c+(1.0-x)*(0.5*(b-d)+(1.0-x)*(0.5*(b+d)-c)))
      return
99999 continue
      end function oned

      subroutine anal2(a2,asta,nsta,iy,jx,cor,xsum,ns,wtmax,htsav)

      use mod_block

      implicit none
!
! Dummy arguments
!
      integer :: iy , jx , nsta
      real(4) , dimension(iy,jx) :: a2 , cor , htsav , xsum , wtmax
      real(4) , dimension(nsta) :: asta
      integer , dimension(iy,jx) :: ns
      intent (in) asta , iy , jx , nsta
      intent (out) a2 , htsav
      intent (inout) cor , ns , xsum , wtmax
!
! Local variables
!
      real(4) :: deltas , dx , dxobs , dy , dyobs , riobs , ris , ris2 ,&
               & rjobs , rsq , rx , ry , wt , x , xcntr , xmaxj ,       &
               & xminj , y , ycntr , ymaxi , ymini
      integer :: i , ie , j , je , kk , maxi , maxj , mini , minj ,     &
               & nscan , nskip
!
!     objective analysis to fill a grid based on observations
!     xobs and yobs are x and y positions on observations, not
!     necessarily grid points.
!
      ie = iy - 1
      je = jx - 1
      nscan = 1
      deltas = dsinm
      print 99001 , rin , deltas
!
!-----grid lengths in x and y directions are unity.
!
      xcntr = jx/2.
      ycntr = iy/2.
      dy = deltas
      dx = deltas
      ris = rin**2
!-----rin is radius of influence in grid units
      ris2 = 2.*ris
!
      nskip = 1
      do j = 1 , jx
        do i = 1 , iy
          cor(i,j) = 0.0
          xsum(i,j) = 0.0
          ns(i,j) = 0
          wtmax(i,j) = 0.0
          htsav(i,j) = 0.0
        end do
      end do
!
!-----begin to process the nobs observations with loop 22
!
      do kk = 1 , nobs , nskip
        if ( asta(kk)<=400. ) then
!
!-----define max and min i and j values to limit the number of points
!-----must be considered.
!
!-----dind obs. location in terms of dx and dy
!
          dxobs = xobs(kk)/dx
          dyobs = yobs(kk)/dy
!-----convert obs. location to grid increment values of i and j
          rjobs = dxobs + xcntr - dxcen
          riobs = dyobs + ycntr - dycen
!
          ymaxi = riobs + rin
          maxi = int(ymaxi+0.99)
          maxi = min0(maxi,ie)
!
          ymini = riobs - rin
          mini = int(ymini)
          mini = max0(mini,1)
!
          xmaxj = rjobs + rin
          maxj = int(xmaxj+0.99)
          maxj = min0(maxj,je)
!
          xminj = rjobs - rin
          minj = int(xminj)
          minj = max0(minj,1)
!
          do i = mini , maxi
            do j = minj , maxj
              x = j - xcntr + dxcen
              y = i - ycntr + dycen
!-----compute distance of k_th station from i,jth grid point
              rx = x - xobs(kk)/dx
              ry = y - yobs(kk)/dy
              rsq = rx**2 + ry**2
              if ( rsq<ris2 ) then
                wt = (ris-rsq)/(ris+rsq)
!               60    continue
!
!-----save      max. weighting factor and terrain height to check if
!-----point     grid should be treated as a land or sea point.
!
                if ( wt>0.0 ) then
                  wtmax(i,j) = amax1(wt,wtmax(i,j))
                  if ( (wt-wtmax(i,j))==0.0 ) htsav(i,j) = asta(kk)
                  cor(i,j) = cor(i,j) + wt*asta(kk)
                  xsum(i,j) = xsum(i,j) + wt
                  ns(i,j) = ns(i,j) + 1
                end if
              end if
            end do
          end do
        end if
      end do
!
!-----now apply summed weights and weighted observations to determine
!-----terrain value at i,j points
!
      do i = 1 , ie
        do j = 1 , je
          if ( ns(i,j)/=0 ) then
            cor(i,j) = cor(i,j)/xsum(i,j)
            a2(i,j) = cor(i,j)
!--------if closest observation to i,j is ocean, override a2(i,j) to
!--------preserve the coastline.
!Sara
           if (htsav(i,j) .le. 0.001) a2(i,j) = htsav(i,j)
!Sara_
          end if
        end do
      end do

!-----may want to smooth final field a2 here

99001 format (1x,' rin,ds(m) =',2E12.3)
!     26 format(' no observations are within rin=',f7.2,
!     & ' grid lengths of i=',i3,' j=',i3)
 
      end subroutine anal2

      end module mod_interp
