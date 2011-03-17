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
 
      module mod_condtq
!
! Heating term for explicit moisture
!
      use mod_runparams
      use mod_main
      use mod_cvaria
      use mod_pmoist
      use mod_slice
!
      private
      real(8) , parameter :: lowq = 1.0D-30
!
      public :: condtq
!
      contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     This subroutine computes the condensational or evaporational    c
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
! 
      subroutine condtq(j)
! 
      implicit none
!
      integer :: j
      intent (in) j
!
! rhc    - Relative humidity at ktau+1
! rh0adj - Adjusted relative humidity threshold at ktau+1
! fccc   - Cloud fraction at ktau+1
!
      real(8) :: dqv , exces , fccc , pres , qvc_cld , qvs , r1 ,       &
               & rh0adj , rhc , satvp
      integer :: i , k
      real(8) , dimension(iy,kz) :: qccs , tmp1 , tmp2 , tmp3
      real(8) , allocatable , save , dimension(:,:) :: qvcs

      if (.not.allocated(qvcs)) then
        allocate(qvcs(iy,kz))
      end if
!
!---------------------------------------------------------------------
!     1.  Compute t, qv, and qc at tau+1 without condensational term
!---------------------------------------------------------------------
      do k = 1 , kz
        do i = 2 , iym2
          tmp3(i,k) = (atm2%t(i,k,j)+dt*aten%t(i,k,j))/psc(i,j)
          qvcs(i,k) = dmax1((atm2%qv(i,k,j)+dt*aten%qv(i,k,j))/psc(i,j),&
                    & lowq)
          qccs(i,k) = dmax1((atm2%qc(i,k,j)+dt*aten%qc(i,k,j))/psc(i,j),&
                    & lowq)
        end do
      end do
 
!---------------------------------------------------------------------
!     2.  Compute the cloud condensation/evaporation term.
!---------------------------------------------------------------------
      do k = 1 , kz
        do i = 2 , iym2
 
!         2a. Calculate the saturation mixing ratio and relative
!         humidity
          pres = (a(k)*psc(i,j)+r8pt)*d_1000
          if ( tmp3(i,k) > tzero ) then
            satvp = svp1*d_1000*dexp(svp2*(tmp3(i,k)-tzero)               &
                  & /(tmp3(i,k)-svp3))
          else
            satvp = svp4*d_1000*dexp(svp5-svp6/tmp3(i,k))
          end if
          qvs = dmax1(ep2*satvp/(pres-satvp),lowq)
          rhc = dmax1(qvcs(i,k)/qvs,lowq)

          r1 = d_one/(d_one+wlhv*wlhv*qvs/ &
                      (rwat*cpd*tmp3(i,k)*tmp3(i,k)))
 
!         2b. Compute the relative humidity threshold at ktau+1
          if ( tmp3(i,k) > tc0 ) then
            rh0adj = rh0(i,j)
          else ! high cloud (less subgrid variability)
            rh0adj = rhmax - (rhmax-rh0(i,j))/ &
                      (d_one+0.15D0*(tc0-tmp3(i,k)))
          end if
 
!         2c. Compute the water vapor in excess of saturation
          if ( rhc >= rhmax .or. rhc < rh0adj ) then
                                                   ! Full or no cloud cover
            dqv = qvcs(i,k) - qvs*conf ! Water vapor in excess of sat
            tmp1(i,k) = r1*dqv
          else                                     ! Partial cloud cover
            fccc = d_one - dsqrt(d_one-(rhc-rh0adj)/(rhmax-rh0adj))
            fccc = dmin1(dmax1(fccc,0.01D0),d_one)
            qvc_cld = dmax1((qsb3d(i,k,j)+dt*aten%qv(i,k,j)/psc(i,j)),  &
                    & d_zero)
            dqv = qvc_cld - qvs*conf       ! qv diff between predicted qv_c
            tmp1(i,k) = r1*dqv*fccc        ! grid cell average
          end if
 
!         2d. Compute the new cloud water + old cloud water
          exces = qccs(i,k) + tmp1(i,k)
          if ( exces >= d_zero ) then
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
      do k = 1 , kz
        do i = 2 , iym2
          aten%qv(i,k,j) = aten%qv(i,k,j) - psc(i,j)*tmp2(i,k)
          aten%qc(i,k,j) = aten%qc(i,k,j) + psc(i,j)*tmp2(i,k)
          aten%t(i,k,j) = aten%t(i,k,j) + psc(i,j)*tmp2(i,k)*wlhvocp
        end do
      end do
 
      end subroutine condtq
!
      end module mod_condtq
