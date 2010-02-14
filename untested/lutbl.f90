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

      subroutine lutbl(ptop)
!
      use mod_bmparam
      use mod_constants , only : cpd , rgas , rovcp
      implicit none
!
! PARAMETER definitions
!
      real(8) , parameter :: a2 = 17.2693882E0 , a3 = 273.16E0 ,        &
                           & a4 = 35.86E0 , pq0 = 379.90516E0 ,         &
                           & eps = 2.E-12 , eliwv = 2.72E6

!
! Dummy arguments
!
      real(8) :: ptop
      intent (in) ptop
!
! Local variables
!
      real(8) :: ape , dp , dqs , dth , dthe , p , ph , pt , qs , qs0k ,&
               & sqsk , sthek , th , the0k , thh
      real(8) , dimension(jtb) :: pnew , pold , &
                                & qsnew , qsold , thenew , theold ,     &
                                & tnew , told , y2p , y2t
      integer :: kp , kpm , kpm1 , kth , kthm , kthm1
!
!--------------coarse look-up table for saturation point----------------
!
      pt = ptop*1000.
!     ptop in pascal
 
      kthm = jtb
      kpm = itb
      kthm1 = kthm - 1
      kpm1 = kpm - 1
!
      thl = 210.
      thh = 385.
      pl = pt
      ph = 105000.
!
      dth = (thh-thl)/dble(kthm-1)
      dp = (ph-pl)/dble(kpm-1)
!
      rdth = 1.0/dth
      rdp = 1.0/dp
      rdq = kpm - 1
!
      th = thl - dth
 
!-----------------------------------------------------------------------
 
      do kth = 1 , kthm
        th = th + dth
        p = pl - dp
        do kp = 1 , kpm
          p = p + dp
          ape = (100000./p)**(rovcp)
          qsold(kp) = pq0/p*exp(a2*(th-a3*ape)/(th-a4*ape))
          pold(kp) = p
        end do
!
        qs0k = qsold(1)
        sqsk = qsold(kpm) - qsold(1)
        qsold(1) = 0.0
        qsold(kpm) = 1.0
!
        do kp = 2 , kpm1
          qsold(kp) = (qsold(kp)-qs0k)/sqsk
!wwwwwwwwwwwwww fix due to cyber half prec. limitation wwwwwwwwwwwwwwwww
          if ( (qsold(kp)-qsold(kp-1)).lt.eps ) qsold(kp) = qsold(kp-1) &
             & + eps
!wwwwwwwwwwwwww fix due to cyber half prec. limitation wwwwwwwwwwwwwwwww

        end do
!
        qs0(kth) = qs0k
        sqs(kth) = sqsk
!-----------------------------------------------------------------------
        qsnew(1) = 0.0
        qsnew(kpm) = 1.0
        dqs = 1.0/dble(kpm-1)
!
        do kp = 2 , kpm1
          qsnew(kp) = qsnew(kp-1) + dqs
        end do
!
        y2p(1) = 0.0
        y2p(kpm) = 0.0
!
        call spline(kpm,qsold,pold,y2p,kpm,qsnew,pnew)
!
        do kp = 1 , kpm
          ptbl(kp,kth) = pnew(kp)
        end do
!-----------------------------------------------------------------------
      end do
!-----------------------------------------------------------------------
 
!--------------coarse look-up table for t(p) from constant the----------
      p = pl - dp
      do kp = 1 , kpm
        p = p + dp
        th = thl - dth
        do kth = 1 , kthm
          th = th + dth
          ape = (100000./p)**(rovcp)
          qs = pq0/p*exp(a2*(th-a3*ape)/(th-a4*ape))
          told(kth) = th/ape
          theold(kth) = th*exp(eliwv*qs/(cpd*told(kth)))
        end do
!
        the0k = theold(1)
        sthek = theold(kthm) - theold(1)
        theold(1) = 0.0
        theold(kthm) = 1.0
!
        do kth = 2 , kthm1
          theold(kth) = (theold(kth)-the0k)/sthek
!wwwwwwwwwwwwww fix due to cyber half prec. limitation wwwwwwwwwwwwwwwww
          if ( (theold(kth)-theold(kth-1)).lt.eps ) theold(kth)         &
             & = theold(kth-1) + eps
!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
        end do
!
        the0(kp) = the0k
        sthe(kp) = sthek
!-----------------------------------------------------------------------
        thenew(1) = 0.0
        thenew(kthm) = 1.0
        dthe = 1.0/dble(kthm-1)
        rdthe = 1.0/dthe
!
        do kth = 2 , kthm1
          thenew(kth) = thenew(kth-1) + dthe
        end do
!
        y2t(1) = 0.0
        y2t(kthm) = 0.0
!
        call spline(kthm,theold,told,y2t,kthm,thenew,tnew)
!
        do kth = 1 , kthm
          ttbl(kth,kp) = tnew(kth)
        end do
!-----------------------------------------------------------------------
      end do
!-----------------------------------------------------------------------
!
      end subroutine lutbl
