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
 
      subroutine condtq(j)
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the condensational or evaporational    c
!     heating term and the fallout term of precipitation from the     c
!     explicit moisture scheme.                                       c
!                                                                     c
!     ---the condensational or evaporational term are one step        c
!        adjustment based on asai (1965, j. meteo. soc. japan).       c
!                                                                     c
!     ---modified to include the effects of partial cloud cover       c
!        (see Pal et al 2000).  When partial clouds exist, the qvten  c
!        in/out of the clear and cloudy portions of the grid cell is  c
!        assumed to be at the same rate (i.e., if there is 80% cloud  c
!        cover, .2 of qvten goes to raising qv in the clear region    c
!        and .8 goes to condensation or evaporation of qc in the      c
!        cloudy portion).                                             c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      use mod_regcm_param
      use mod_param1
      use mod_param3
      use mod_main
      use mod_cvaria
      use mod_pmoist
      use mod_slice
      use mod_constants , only : cpd , ep2 , wlhv , wlhvocp , svp1 ,    &
                               & svp2 , svp3 , tmelt
      implicit none
!
! Dummy arguments
!
      integer :: j
      intent (in) j
!
! Local variables
!
! rhc    - Relative humidity at ktau+1
! rh0adj - Adjusted relative humidity threshold at ktau+1
! fccc   - Cloud fraction at ktau+1
!
      real(8) :: dqv , exces , fccc , pres , qvc_cld , qvs , r1 ,       &
               & rh0adj , rhc , satvp
      integer :: i , k
      real(8) , dimension(ix,kx) :: qccs , tmp1 , tmp2 , tmp3
!
 
!---------------------------------------------------------------------
!     1.  Compute t, qv, and qc at tau+1 without condensational term
!---------------------------------------------------------------------
      do k = 1 , kx
        do i = 2 , ixm2
          tmp3(i,k) = (tb(i,k,j)+dt*tten(i,k,j))/psc(i,j)
          qvcs(i,k) = dmax1((qvb(i,k,j)+dt*qvten(i,k,j))/psc(i,j),      &
                    & 1.D-30)
          qccs(i,k) = dmax1((qcb(i,k,j)+dt*qcten(i,k,j))/psc(i,j),      &
                    & 1.D-30)
        end do
      end do
 
!---------------------------------------------------------------------
!     2.  Compute the cloud condensation/evaporation term.
!---------------------------------------------------------------------
      do k = 1 , kx
        do i = 2 , ixm2
 
!         2a. Calculate the saturation mixing ratio and relative
!         humidity
          pres = (a(k)*psc(i,j)+ptop)*1000.
          if ( tmp3(i,k).gt.tmelt ) then
            satvp = svp1*1.E3*dexp(svp2*(tmp3(i,k)-tmelt)               &
                  & /(tmp3(i,k)-svp3))
          else
            satvp = .611*1.E3*dexp(22.514-6.15E3/tmp3(i,k))
          end if
          qvs = dmax1(ep2*satvp/(pres-satvp),1.D-30)
          rhc = dmax1(qvcs(i,k)/qvs,1.D-30)
 
          r1 = 1./(1.+wlhv*wlhv*qvs/(rv*cpd*tmp3(i,k)*tmp3(i,k)))
 
!         2b. Compute the relative humidity threshold at ktau+1
          if ( tmp3(i,k).gt.tc0 ) then
            rh0adj = rh0(i,j)
          else ! high cloud (less subgrid variability)
            rh0adj = rhmax - (rhmax-rh0(i,j))/(1.0+0.15*(tc0-tmp3(i,k)))
          end if
 
!         2c. Compute the water vapor in excess of saturation
          if ( rhc.ge.rhmax .or. rhc.lt.rh0adj ) then
                                                   ! Full or no cloud cover
            dqv = qvcs(i,k) - qvs*conf ! Water vapor in excess of sat
            tmp1(i,k) = r1*dqv
          else                                     ! Partial cloud cover
            fccc = 1. - sqrt(1.-(rhc-rh0adj)/(rhmax-rh0adj))
            fccc = dmin1(dmax1(fccc,0.01D0),1.0D0)
            qvc_cld = dmax1((qsb3d(i,k,j)+dt*qvten(i,k,j)/psc(i,j)),    &
                    & 0.0D0)
            dqv = qvc_cld - qvs*conf       ! qv diff between predicted qv_c
            tmp1(i,k) = r1*dqv*fccc        ! grid cell average
          end if
 
!         2d. Compute the new cloud water + old cloud water
          exces = qccs(i,k) + tmp1(i,k)
          if ( exces.ge.0. ) then
                              ! Some cloud is left
            tmp2(i,k) = tmp1(i,k)/dt
          else                ! The cloud evaporates
            tmp2(i,k) = -qccs(i,k)/dt
          end if
 
        end do
      end do
 
!---------------------------------------------------------------------
!     3.  Compute the tendencies.
!---------------------------------------------------------------------
      do k = 1 , kx
        do i = 2 , ixm2
          qvten(i,k,j) = qvten(i,k,j) - psc(i,j)*tmp2(i,k)
          qcten(i,k,j) = qcten(i,k,j) + psc(i,j)*tmp2(i,k)
          tten(i,k,j) = tten(i,k,j) + psc(i,j)*tmp2(i,k)*wlhvocp
        end do
      end do
 
      end subroutine condtq
