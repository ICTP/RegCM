      subroutine anal2(a2,asta,nsta,iy,jx,cor,xsum,ns,wtmax,htsav)
      use mod_addstack
      use mod_block
      use mod_const
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
!           if (htsav(i,j) .le. 0.001) a2(i,j) = htsav(i,j)
!Sara_
          end if
        end do
      end do

!-----may want to smooth final field a2 here

99001 format (1x,' rin,ds(m) =',2E12.3)
!     26 format(' no observations are within rin=',f7.2,
!     & ' grid lengths of i=',i3,' j=',i3)
 
      end subroutine anal2
