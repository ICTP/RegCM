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
 
      subroutine cuparan(tten,qten,j)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      use mod_regcm_param
      use mod_param1 , only : dt , dtmin , nbatst
      use mod_param3 , only : ptop , a
      use mod_pmoist
      use mod_rad
      use mod_bats , only : pptc
      use mod_main
      use mod_constants , only : rgti
      use mod_date , only : jyear , jyear0 , ktau
      implicit none
!
! Dummy arguments
!
      integer :: j
      real(8) , dimension(iy,kz) :: qten , tten
      intent (inout) qten , tten
!
! Local variables
!
      real(8) :: aprdiv , calc , dtime , pkdcut , pkk , prainx , us , vs
      integer :: i , iconj , icut , iend , istart , k , kk
      integer , dimension(iy) :: kdet
      real(8) , dimension(iy,kz) :: outq , outt , p , po , q , qo , t , &
                                  & tn , vsp
      real(8) , dimension(iy) :: pret , psur , qcrit , ter11
!
!     zero out radiative clouds
!
      cldlwc = 0.0
      cldfra = 0.0
 
      icut = 0
      dtime = dt
      pkdcut = 75.
      istart = 2 + icut
      iend = iym2 - icut
!
!---  prepare input, erase output
!
      do i = istart , iend
        kdet(i) = 2
        qcrit(i) = 0.
        pret(i) = 0.
      end do

      do k = 1 , kz
        do i = 2 + icut , iym2 - icut
          kk = kz - k + 1
          us = ua(i,kk,j)/psb(i,j)
          vs = va(i,kk,j)/psb(i,j)
          us = 0.25*(ua(i,kk,j)/psb(i,j)+ua(i+1,kk,j)/psb(i+1,j)        &
             & +ua(i,kk,j+1)/psb(i,j+1)+ua(i+1,kk,j+1)/psb(i+1,j+1))
          vs = 0.25*(va(i,kk,j)/psb(i,j)+va(i+1,kk,j)/psb(i+1,j)        &
             & +va(i,kk,j+1)/psb(i,j+1)+va(i+1,kk,j+1)/psb(i+1,j+1))
          t(i,k) = tb(i,kk,j)/psb(i,j)
          q(i,k) = qvb(i,kk,j)/psb(i,j)
          if ( q(i,k).lt.1.E-08 ) q(i,k) = 1.E-08
          tn(i,k) = t(i,k) + (tten(i,kk))/psb(i,j)*dtime
          qo(i,k) = q(i,k) + (qten(i,kk))/psb(i,j)*dtime
          p(i,k) = 10.*psb(i,j)*a(kk) + 10.*ptop
          vsp(i,k) = dsqrt(us**2+vs**2)
          if ( qo(i,k).lt.1.E-08 ) qo(i,k) = 1.E-08
!
          po(i,k) = p(i,k)
          psur(i) = 10.*psb(i,j) + 10.*ptop
          outt(i,k) = 0.
          pkk = psur(i) - po(i,k)
          if ( pkk.le.pkdcut ) kdet(i) = kdet(i) + 1
          outq(i,k) = 0.
          ter11(i) = ht(i,j)*rgti
          if ( ter11(i).le.0. ) ter11(i) = 1.E-05
          qcrit(i) = qcrit(i) + qten(i,kk)
        end do
      end do
!
!---  call cumulus parameterization
!
      call cup(qcrit,t,q,ter11,tn,qo,po,pret,p,outt,outq,dtime,psur,vsp,&
             & istart,iend,kdet,j)
      do k = 1 , kz
        do i = istart , iend
          if ( pret(i).gt.0. ) then
            kk = kz - k + 1
            tten(i,kk) = psb(i,j)*outt(i,k) + tten(i,kk)
            qten(i,kk) = psb(i,j)*outq(i,k) + qten(i,kk)
          end if
        end do
      end do
!
!---  rain in cm.
!
      calc = .5
      iconj = 0
      do i = istart , iend
        if ( pret(i).gt.0. ) then
          rainc(i,j) = rainc(i,j) + pret(i)*calc*dt
!         print *,'rainc(',i,j,')=',rainc(i,j)
          iconj = iconj + 1
!.....................precipitation rate for bats (mm/s)
          aprdiv = dble(nbatst)
          if ( jyear.eq.jyear0 .and. ktau.eq.0 ) aprdiv = 1.
          prainx = pret(i)*calc*dt
          pptc(i,j) = pptc(i,j) + prainx/(dtmin*60.)/aprdiv
!.......................................................
        end if
      end do
      icon(j) = iconj
!
      end subroutine cuparan
