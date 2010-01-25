C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine condtq(j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the condensational or evaporational    c
c     heating term and the fallout term of precipitation from the     c
c     explicit moisture scheme.                                       c
c                                                                     c
c     ---the condensational or evaporational term are one step        c
c        adjustment based on asai (1965, j. meteo. soc. japan).       c
c                                                                     c
c     ---modified to include the effects of partial cloud cover       c
c        (see Pal et al 2000).  When partial clouds exist, the qvten  c
c        in/out of the clear and cloudy portions of the grid cell is  c
c        assumed to be at the same rate (i.e., if there is 80% cloud  c
c        cover, .2 of qvten goes to raising qv in the clear region    c
c        and .8 goes to condensation or evaporation of qc in the      c
c        cloudy portion).                                             c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer j
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/iunits.cb"
#     include "include/main.cb"
#     include "include/cvaria.cb"
#     include "include/pmoist.cb"
#     include "include/bats.cb2"
#     include "include/slice.cb"

      real(kind=8)  tmp1(ix,kx),tmp2(ix,kx),tmp3(ix,kx),qccs(ix,kx)
      integer i,k
      real(kind=8)  pres,satvp,qvs,dqv,r1,exces
      real(kind=8)  rhc     ! Relative humidity at ktau+1
      real(kind=8)  rh0adj  ! Adjusted relative humidity threshold at ktau+1
      real(kind=8)  fccc    ! Cloud fraction at ktau+1
      real(kind=8)  qvc_cld

c---------------------------------------------------------------------
c   1.  Compute t, qv, and qc at tau+1 without condensational term
c---------------------------------------------------------------------
      do k=1,kx
      do i=2,ilxm
        tmp3(i,k)=(tb(i,k,j)+dt*tten(i,k,j))/psc(i,j)
        qvcs(i,k)=dmax1((qvb(i,k,j)+dt*qvten(i,k,j))/psc(i,j),1.d-30)
        qccs(i,k)=dmax1((qcb(i,k,j)+dt*qcten(i,k,j))/psc(i,j),1.d-30)
      end do
      end do

c---------------------------------------------------------------------
c   2.  Compute the cloud condensation/evaporation term.
c---------------------------------------------------------------------
      do k=1,kx
      do i=2,ilxm

c   2a. Calculate the saturation mixing ratio and relative humidity
        pres=(a(k)*psc(i,j)+ptop)*1000.
        if(tmp3(i,k).gt.273.15)then
          satvp=svp1*1.e3*dexp(svp2*(tmp3(i,k)-273.15)/(tmp3(i,k)-svp3))
        else
          satvp=.611*1.e3*dexp(22.514-6.15e3/tmp3(i,k))
        end if
        qvs=dmax1(ep2*satvp/(pres-satvp),1.d-30)
        rhc=dmax1(qvcs(i,k)/qvs,1.d-30)

        r1=1./(1.+xlv*xlv*qvs/(rv*cp*tmp3(i,k)*tmp3(i,k)))

c   2b. Compute the relative humidity threshold at ktau+1
        if (tmp3(i,k).gt.tc0) then
          rh0adj = rh0(i,j)
        else   ! high cloud (less subgrid variability)
          rh0adj = rhmax - (rhmax-rh0(i,j))/(1.0+0.15*(tc0-tmp3(i,k)))
        end if

c   2c. Compute the water vapor in excess of saturation
        if (rhc.ge.rhmax .or. rhc.lt.rh0adj) then  ! Full or no cloud cover
          dqv = qvcs(i,k) - qvs*conf   ! Water vapor in excess of sat
          tmp1(i,k) = r1*dqv
        else                                       ! Partial cloud cover
          fccc = 1. - sqrt(1.-(rhc-rh0adj)/(rhmax-rh0adj))
          fccc = dmin1(dmax1(fccc,0.01D0),1.0d0)
          qvc_cld = dmax1((qsb3d(i,k,j)+dt*qvten(i,k,j)/psc(i,j)),0.0d0)
          dqv = qvc_cld - qvs*conf         ! qv diff between predicted qv_c
          tmp1(i,k) = r1*dqv*fccc          ! grid cell average 
        end if

c   2d. Compute the new cloud water + old cloud water
        exces = qccs(i,k) + tmp1(i,k)
        if (exces.ge.0.) then ! Some cloud is left
          tmp2(i,k) = tmp1(i,k)/dt
        else                  ! The cloud evaporates
          tmp2(i,k) = -qccs(i,k)/dt
        end if

      end do
      end do

c---------------------------------------------------------------------
c   3.  Compute the tendencies.
c---------------------------------------------------------------------
      do k=1,kx
      do i=2,ilxm
        qvten(i,k,j) = qvten(i,k,j) - psc(i,j)*tmp2(i,k)
        qcten(i,k,j) = qcten(i,k,j) + psc(i,j)*tmp2(i,k)
        tten(i,k,j)  = tten(i,k,j)  + psc(i,j)*tmp2(i,k)*xlvocp
      end do
      end do

      return
      end
