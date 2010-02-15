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
 
      subroutine nconvp(psa,psb,ta,tb,qva,qvb,qca,qcb)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the nonconvective precipitation, i.e.  c
!     removes the supersaturation as nonconvective precipitation.     c
!                                                                     c
!     ---the computation procedure uses the averaged value from       c
!        the two timesteps to check for supersaturation and removes   c
!        the supersaturation from both timesteps in order to conserve c
!        the total mass.                                              c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mod_regcm_param
      use mod_param3 , only : ptop , a
      use mod_pmoist
      use mod_constants , only : cpd , ep2 , wlhv , wlhvocp , svp1 ,    &
                             & svp2 , svp3 , tmelt
      implicit none
!
! Dummy arguments
!
#ifdef MPP1
      real(8) , dimension(ix,jxp) :: psa , psb
      real(8) , dimension(ix,kx,jxp) :: qca , qcb , qva , qvb , ta , tb
#else
      real(8) , dimension(ix,jx) :: psa , psb
      real(8) , dimension(ix,kx,jx) :: qca , qcb , qva , qvb , ta , tb
#endif
      intent (in) psa , psb
      intent (inout) qca , qcb , qva , qvb , ta , tb
!
! Local variables
!
      real(8) :: dqv , e1 , es1 , ex , exces , pavg , psx , qcas ,      &
               & qcavg , qcbs , qvas , qvavg , qvbs , r1 , tavg ,       &
               & tpavg , tta
      integer :: i , j , k
      real(8) , dimension(ix,kx) :: scr
!
!----------------------------------------------------------------------
#ifdef MPP1
      do j = 1 , jxp
#else
      do j = 1 , jxm1
#endif
        do k = 1 , kx
          do i = 1 , ixm1
            tavg = 0.5*(ta(i,k,j)/psa(i,j)+tb(i,k,j)/psb(i,j))
            tpavg = 0.5*(ta(i,k,j)+tb(i,k,j))
            qvavg = 0.5*(qva(i,k,j)/psa(i,j)+qvb(i,k,j)/psb(i,j))
            pavg = 0.5*(psa(i,j)+psb(i,j))
            psx = pavg*a(k) + ptop
            tta = tpavg/pavg
            if ( tta.gt.tmelt ) then
!             v8 svp formula
              e1 = svp1*dexp(svp2*(tta-tmelt)/(tta-svp3))
            else
              e1 = .611*dexp(22.514-6.15E3/tta)
            end if
            es1 = ep2*e1/(psx-e1)
            dqv = qvavg - es1*conf
            r1 = 1./(1.+wlhv*wlhv*es1/(rv*cpd*tavg*tavg))
            scr(i,k) = r1*dqv
          end do
        end do
!
        do k = 1 , kx
          do i = 1 , ixm1
            qcavg = dmin1(qca(i,k,j)/psa(i,j),qcb(i,k,j)/psb(i,j))
            exces = qcavg + scr(i,k)
            if ( exces.ge.0. ) then
              scr(i,k) = scr(i,k)
            else
              scr(i,k) = -qcavg
            end if
            qcas = qca(i,k,j) + scr(i,k)*psa(i,j)
            qcbs = qcb(i,k,j) + scr(i,k)*psb(i,j)
            qca(i,k,j) = dmax1(qcas,0.D0)
            qcb(i,k,j) = dmax1(qcbs,0.D0)
          end do
        end do
!
        do k = 1 , kx
          do i = 1 , ixm1
            qvas = qva(i,k,j) - scr(i,k)*psa(i,j)
            qvbs = qvb(i,k,j) - scr(i,k)*psb(i,j)
            qva(i,k,j) = dmax1(qvas,1.D-99)
            qvb(i,k,j) = dmax1(qvbs,1.D-99)
          end do
        end do
!
!.....compute nonconvective rainfall when cumulus parameterization
!       scheme is used (i.e. im = 1).
!       unit for precipitation is cm.
!
#ifdef MPP1
        if ( .not.((myid.eq.0 .and. j.eq.1) .or. (myid.eq.nproc-1 .and. &
           & j.eq.jxp-1)) ) then
#else
        if ( j.ne.1 .and. j.ne.jxm1 ) then
#endif
!
          do k = 1 , kx
            do i = 2 , ixm2
              ex = 1.
              ta(i,k,j) = ta(i,k,j) + wlhvocp*scr(i,k)*psa(i,j)*ex
              tb(i,k,j) = tb(i,k,j) + wlhvocp*scr(i,k)*psb(i,j)*ex
            end do
          end do
        end if
!
      end do     !end j=1,jxm1 loop
!
      end subroutine nconvp
