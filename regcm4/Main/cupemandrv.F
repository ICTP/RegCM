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

      subroutine cupemandrv(j)

c **********************************************
c **** Driver for Emanuel Convection Scheme ****
c **********************************************

      implicit none
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/main.cb"
#     include "include/pmoist.cb"
#     include "include/cvaria.cb"
#     include "include/slice.cb"
#     include "include/rad.cb"
#     include "include/bats2.cb"

      integer ntra, nl, i, j, k, kk, iconj
      parameter (ntra=0,nl=kx-1)

c **** Local definition of Input Variables to CUPEMAN
      real(kind=8) tcup(kx), qcup(kx), qscup(kx), ucup(kx), vcup(kx)
c    &     , tra(kx,ntra), pcup(kx), phcup(kxp1), dtime
     &     , tra(kx,1)   , pcup(kx), phcup(kxp1), dtime

c **** Local definition of Output Variables from CUPEMAN
      integer iflag
c     real(kind=8) ft(kx), fq(kx), fu(kx), fv(kx), ftra(kx,ntra)
      real(kind=8) ft(kx), fq(kx), fu(kx), fv(kx), ftra(kx,1)   
     &     , fppt, wd, tprime, qprime, cbmf

c **** Other local variable definition
      integer kbase, ktop, kclth
      real(kind=8) uconv, aprdiv, akclth


      dtime=dt
      uconv = 0.5*dt
      aprdiv=1.0/float(nbatst)
      if (jyear.eq.jyear0.and.ktau.eq.0) aprdiv=1.
      iconj=0
      do i=2,ilxm
        do k=1,kx
          kk = kxp1-k
          cldlwc(i,k)=0.            ! Zero out cloud water content
          cldfra(i,k)=0.            ! Zero out cloud fraction coverage
          tcup(k) = tb3d(i,kk,j)                          ! [k]
          qcup(k) = qvb3d(i,kk,j)/(1.+qvb3d(i,kk,j))      ! [kg/kg]
          qscup(k) = qsb3d(i,kk,j)/(1.+qsb3d(i,kk,j))     ! [kg/kg]
          ucup(k) = ubx3d(i,kk,j)                         ! [m/s]
          vcup(k) = vbx3d(i,kk,j)                         ! [m/s]
          pcup(k) = pb3d(i,kk,j)*10.                      ! [hPa]
        end do
        do k=1,kxp1
          kk = kxp1-k+1
          phcup(k) = (sigma(kk)*psb(i,j)+ptop)*10.        ! [hPa]
        end do
        cbmf = cbmf2d(i,j)                                ! [(kg/m**2)/s]

        call cupeman(tcup,qcup,qscup,ucup,vcup,tra,pcup,phcup
     &     , kx,kxp1,nl,ntra,dtime,iflag,ft,fq,fu
     &     , fv,ftra,fppt,wd,tprime,qprime,cbmf,kbase,ktop)

        cbmf2d(i,j) = cbmf

c iflag=0: No moist convection; atmosphere stable or surface temperature < 250K
c          or surface humidity is negative.
c iflag=1: Moist convection occurs.
c iflag=2: No moist convection: lifted condensation level above 200 mb.
c iflag=3: No moist convection: cloud base higher than level kx-2.
c iflag=4: Moist convection occurs, but CFL condition on the subsidence
c          warming is violated. (Does not terminate scheme.)
        if (iflag.eq.1 .or. iflag.eq.4) then  ! If moist convection

c         if (iflag.eq.4) then                ! If CFL violation
c           print*,'EMAN CFL VIOLATION: ',i,j,cbmf
c         end if

c **** Tendencies
          do k=1,kx
            kk = kxp1-k
            tten(i,kk,j) = ft(k)*psb(i,j) + tten(i,kk,j)
            qvten(i,kk,j) = fq(k)/(1.-fq(k))*psb(i,j) + qvten(i,kk,j)
c There is a bit of an inconsistency here...  The wind tendencies from
c convection are on cross points, but the model wants them on dot points.
            uten(i,kk,j) = fu(k)*psb(i,j) + uten(i,kk,j)
            vten(i,kk,j) = fv(k)*psb(i,j) + vten(i,kk,j)
          end do

c **** Cloud fraction and cloud water
          kclth=ktop-kbase+1
          akclth=1./float(kclth)
          do k=kbase,ktop
            kk = kxp1 - k
            cldlwc(i,kk)=cllwcv
            cldfra(i,kk)=1.-(1.-clfrcv)**akclth
          end do

c **** Precipitation
          if (fppt.gt.0.) then
            rainc(i,j) = rainc(i,j) + fppt*uconv ! mm
            pptc(i,j) = pptc(i,j) + fppt*aprdiv  ! mm/s
            iconj=iconj+1
          end if

        end if

      end do

      icon(j)=iconj


      return
      end subroutine cupemandrv
