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

      subroutine nconvp(psa,psb,ta,tb,qva,qvb,qca,qcb)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the nonconvective precipitation, i.e.  c
c     removes the supersaturation as nonconvective precipitation.     c
c                                                                     c
c     ---the computation procedure uses the averaged value from       c
c        the two timesteps to check for supersaturation and removes   c
c        the supersaturation from both timesteps in order to conserve c
c        the total mass.                                              c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/iunits.cb"
#     include "include/pmoist.cb"
#     include "include/bats2.cb"
#ifdef MPP1
#     include "mpiregcm.inc"
#     include "include/io.cb"
      real(kind=8)  psa(ix,jxp),psb(ix,jxp)
      real(kind=8)  ta(ix,kx,jxp), tb(ix,kx,jxp)
      real(kind=8)  qva(ix,kx,jxp),qvb(ix,kx,jxp)
      real(kind=8)  qca(ix,kx,jxp),qcb(ix,kx,jxp)
#else
      real(kind=8)  psa(ix,jx),psb(ix,jx)
      real(kind=8)  ta(ix,kx,jx),tb(ix,kx,jx)
      real(kind=8)  qva(ix,kx,jx),qvb(ix,kx,jx)
      real(kind=8)  qca(ix,kx,jx),qcb(ix,kx,jx)
#endif
c
      real(kind=8)  scr(ix,kx)
      integer i,j,k
      real(kind=8)  tavg,tpavg,qvavg,pavg,psx,tta,e1,es1,dqv,r1
      real(kind=8)  qcavg,exces,qcas,qcbs,qvas,qvbs,ex
c----------------------------------------------------------------------
#ifdef MPP1
      do 100 j = 1,jxp
#else
      do 100 j = 1,jx-1
#endif
c
       do 10 k=1,kx
       do 10 i = 1,ix-1
        tavg=0.5*(ta(i,k,j)/psa(i,j)+tb(i,k,j)/psb(i,j))
        tpavg=0.5*(ta(i,k,j)+tb(i,k,j))
        qvavg=0.5*(qva(i,k,j)/psa(i,j)+qvb(i,k,j)/psb(i,j))
        pavg=0.5*(psa(i,j)+psb(i,j))
        psx=pavg*a(k)+ptop
        tta=tpavg/pavg
        if(tta.gt.273.15)then
c     v8 svp formula
           e1=svp1*dexp(svp2*(tta-273.15)/(tta-svp3))
        else
           e1=.611*dexp(22.514-6.15e3/tta)
        end if
        es1=ep2*e1/(psx-e1)
        dqv=qvavg-es1*conf
        r1=1./(1.+xlv*xlv*es1/(rv*cp*tavg*tavg))
        scr(i,k)=r1*dqv
10     continue
c
       do 20 k=1,kx
       do 20 i = 1,ix-1
        qcavg=dmin1(qca(i,k,j)/psa(i,j),qcb(i,k,j)/psb(i,j))
        exces=qcavg+scr(i,k)
        if(exces.ge.0.)then
         scr(i,k)=scr(i,k)
        else
         scr(i,k)=-qcavg
        end if
        qcas=qca(i,k,j)+scr(i,k)*psa(i,j)
        qcbs=qcb(i,k,j)+scr(i,k)*psb(i,j)
        qca(i,k,j)=dmax1(qcas,0.d0)
        qcb(i,k,j)=dmax1(qcbs,0.d0)
20     continue
c
       do 60 k=1,kx
       do 60 i = 1,ix-1
        qvas=qva(i,k,j)-scr(i,k)*psa(i,j)
        qvbs=qvb(i,k,j)-scr(i,k)*psb(i,j)
        qva(i,k,j)=dmax1(qvas,1.d-99)
        qvb(i,k,j)=dmax1(qvbs,1.d-99)
60     continue
c
c.....compute nonconvective rainfall when cumulus parameterization
c     scheme is used (i.e. im = 1).
c     unit for precipitation is cm.
c
#ifdef MPP1
       if ((myid .eq. 0 .and. j .eq. 1) .or.
     &     (myid .eq. nproc-1 .and. j .eq.jxp-1)) go to 100
#else
       if (j .eq. 1 .or. j .eq. jx-1) go to 100
#endif
c
       do 90 k=1,kx
       do 90 i=2,ix-2
        ex = 1.
        ta(i,k,j)=ta(i,k,j)+xlvocp*scr(i,k)*psa(i,j)*ex
        tb(i,k,j)=tb(i,k,j)+xlvocp*scr(i,k)*psb(i,j)*ex
90     continue
c
100   continue   !end j=1,jx-1 loop
c
c
      return
      end subroutine nconvp
