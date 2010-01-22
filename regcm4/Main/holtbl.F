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
      subroutine holtbl
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      integer j
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'mpiregcm.h'
      INTEGER ierr
#endif
#     include "include/param1.cb"
#ifdef DIAG
#     include "include/diagnosis.cb"
#endif
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/iunits.cb"
#     include "include/main.cb"
#     include "include/mainchem.cb"
#     include "include/pbldim.cb"
#     include "include/cvaria.cb"
#     include "include/pmoist.cb"
#     include "include/bats.cb2"
#     include "include/slice.cb"
#     include "include/trachem.cb"
#ifdef MPP1
#     include "include/io.cb"
      common/kuv/ akxx1(ix-1,kx,0:jxp+1),akxx2(ix-1,kx,0:jxp+1)
#else
      common/kuv/ akxx1(ix-1,kx,jx-1), akxx2(ix-1,kx,jx-1)
#endif
      real(kind=8)  akxx1
      real(kind=8)  akxx2
Cchem
      real(kind=8)  vdep(ix-1,ntr)
Cchem_
c
#ifdef MPP1
      integer status(MPI_STATUS_SIZE)
      real(kind=8)  kvm(ix,kx,jxp),kvh(ix,kx,jxp)
      real(kind=8)  kvq(ix,kx,jxp),kvc(ix,kx,jxp),cgh(ix,kx,jxp)
      real(kind=8)  xhfx(ix,jxp),xqfx(ix,jxp),hfxv(ix,jxp)
      real(kind=8)  obklen(ix,jxp),th10(ix,jxp),ustr(ix,jxp)
#else
      real(kind=8)  kvm(ix,kx,jlx),kvh(ix,kx,jlx)
      real(kind=8)  kvq(ix,kx,jlx),kvc(ix,kx,jlx),cgh(ix,kx,jlx)
      real(kind=8)  xhfx(ix,jx),xqfx(ix,jx),hfxv(ix,jx)
      real(kind=8)  obklen(ix,jx),th10(ix,jx),ustr(ix,jx)
#endif
      COMMON /BLH_TMP/ kvm,kvh,kvq,kvc,cgh
     &       , xhfx,xqfx,obklen,hfxv,th10,ustr
C
#ifdef MPP1
      real(kind=8)  auxx(ilx,kx,jxp),avxx(ilx,kx,jxp),qcx(ilx,kx,jxp)
      real(kind=8)  dza(ilx,kx,jxp),rhohf(ix,kx,jxp)
      real(kind=8)  wksend(2*(ix-2)*kx),wkrecv(2*(ix-2)*kx)
      integer ii
#else
      real(kind=8)  auxx(ilx,kx,jlx),avxx(ilx,kx,jlx),qcx(ilx,kx,jlx)
      real(kind=8)  dza(ilx,kx,jlx),rhohf(ix,kx,jlx)
#endif
      real(kind=8)  govrth(ilx),ttnp(ilx,kx),rc(ilx,kx),kzm(ilx,kx),
     &        coef1(ix,kx),coef2(ix,kx),coef3(ix,kx),
     &        coefe(ix,kx),coeff1(ix,kx),coeff2(ix,kx),
     &        alphak(ix,kx),betak(ix,kx),tpred1(ix,kx),
     &        tpred2(ix,kx)

Cchem2
      real(kind=8) chix(ix,kx)
Cchem2_
      real(kind=8)  kzo
      real(kind=8)  szkm
      data    kzo/1./
      data    szkm/1600./
c
      integer i,k,idx,jdxm1,jdx,idxm1,itr
      real(kind=8)  dumr,tvcon,ps,ps2,kzmax,sh10
      real(kind=8)  ss,ri,drgdot,oblen,sf,uflxsf,vflxsf,uflxsfx,vflxsfx
c
c
c *********************************************************************
c
c   diagnostic on total evaporation
c
#ifdef DIAG
#ifdef MPP1
      call MPI_GATHER(qfx(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &             qfx_io(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &                              0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
        do j=2,mjx-2
          do i=2,ilxm
            tqeva=tqeva+qfx_io(i,j)*dx*dx*dtmin*60.
          enddo
        enddo
      endif
      call MPI_BCAST(tqeva,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#else
      do j=2,jx-2
        do i=2,ilxm
          tqeva=tqeva+qfx(i,j)*dx*dx*dtmin*60.
        enddo
      enddo
#endif
#endif
c----------------------------------------------------------------------
c-----some of the storage spaces for high-resolution pbl
c     are use to store the variables in this subroutine.
c     difft(i,k,j)   : temperature tendency (tten)
c     diffq(i,k,j)   : water vapor tendency (qvten)
c
c-----decouple flux-form variables to give u,v,t,theta,theta-vir.,
c     t-vir., qv, and qc at cross points and at ktau-1.
c
c     *** note ***
c         the boundary winds may not be adequately affected by friction,
c         so use only interior values of ubx3d and vbx3d to calculate
c         tendencies.
c
#ifdef MPP1
      call MPI_SENDRECV(psb(1,jxp),ix,MPI_DOUBLE_PRECISION,ieast,1,
     &                  psb(1,0),ix,MPI_DOUBLE_PRECISION,iwest,1,
     &                  MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(uvdrag(1,jxp),ix,MPI_DOUBLE_PRECISION,ieast,1,
     &                  uvdrag(1,0),ix,MPI_DOUBLE_PRECISION,iwest,1,
     &                  MPI_COMM_WORLD,status,ierr)
#endif
#ifdef MPP1
      do 1079 j=jbegin,jendx
#else
      do 1079 j=2,jlx
#endif
      do 20 k=1,kx
      do 20 i=2,ilx
       dumr = 4. / (psb(i,j)+psb(i,j-1)+psb(i-1,j)+psb(i-1,j-1))
       auxx(i,k,j)=ub(i,k,j) * dumr
       avxx(i,k,j)=vb(i,k,j) * dumr
20    continue

c
      do 40 k=1,kx
      do 40 i=2,ilx
       tvcon=(1.+ep1*qvb3d(i,k,j))
       thvx(i,k,j)=thx3d(i,k,j)*tvcon
40    continue
c
      do 50 k=1,kx
      do 50 i=2,ilx
       qcx(i,k,j)=qcb(i,k,j)/psb(i,j)
50    continue
c
c.....density at surface is stored in rhox2d(i,j), at half levels in
c     rhohf(i,k,j).
c
      do 130 k=1,kxm
      do 130 i=2,ilx
       dza(i,k,j)=za(i,k,j)-za(i,k+1,j)
       ps=(a(k)*psb(i,j)+ptop)*1000.
       ps2=(a(k+1)*psb(i,j)+ptop)*1000.
       rhohf(i,k,j)=(ps2-ps)/(g*dza(i,k,j))
130   continue
c
      do 140 i=2,ilx
       govrth(i)=g/thx3d(i,kx,j)
140   continue
c
c *********************************************************************
c
c-----compute the vertical diffusion term:
c
      do 160 k=2,kx
      do 160 i=2,ilx
       rc(i,k)=0.257*dzq(i,k,j)**0.175
160   continue
c
c-----compute the diffusion coefficient:
c
c  blackadar scheme above boundary layer top
c
      do 170 k=2,kx
      do 170 i=2,ilx
       kzmax=0.8*dza(i,k-1,j)*dzq(i,k,j)/dt
       ss=((ubx3d(i,k-1,j)-ubx3d(i,k,j))*(ubx3d(i,k-1,j)-ubx3d(i,k,j))
     &    +(vbx3d(i,k-1,j)-vbx3d(i,k,j))*(vbx3d(i,k-1,j)-vbx3d(i,k,j)))
     &   /(dza(i,k-1,j)*dza(i,k-1,j))+1.e-9
       ri=govrth(i)*(thvx(i,k-1,j)-thvx(i,k,j))/(ss*dza(i,k-1,j))
       if((ri-rc(i,k)).ge.0.)then
        kzm(i,k)=kzo
       else
        kzm(i,k)=kzo+dsqrt(ss)*(rc(i,k)-ri)*szkm/rc(i,k)
       endif
       kzm(i,k)=dmin1(kzm(i,k),kzmax)
170   continue
c
c *********************************************************************
c
c            holtslag pbl
c
c initialize bl diffusion coefficients and counter-gradient terms
c with free atmosphere values and make specific humidity
c
      do 180 k=2,kx
      do 180 i=2,ilx
c       eddy diffusivities for momentum, heat and moisture
        kvm(i,k,j) = kzm(i,k)
        kvh(i,k,j) = kzm(i,k)
        kvq(i,k,j) = kzm(i,k)
Cchem
        if(ichem.eq.1) kvc(i,k,j) = kzm(i,k)
Cchem_
c       counter gradient terms for heat and moisture
        cgh(i,k,j) = 0.0
180   continue

      do 190 i=2,ilx
c compute friction velocity
       idx = i
       idx = min0(idx,ilx)
       idxm1 = i-1
       idxm1 = max0(idxm1,2)
       jdx = j
#ifdef MPP1
       if (myid .eq. nproc-1) then
         jdx = min0(jdx,jendx)
       endif
#else
       jdx = min0(jdx,jlx)
#endif
       jdxm1 = j-1
#ifdef MPP1
       if (myid.eq.0) then
         jdxm1 = max0(jdxm1,2)
       endif
#else
       jdxm1 = max0(jdxm1,2)
#endif
       uflxsfx = uvdrag(idx,jdx)*ubx3d(i,kx,j)
       vflxsfx = uvdrag(idx,jdx)*vbx3d(i,kx,j)
       ustr(i,j) = dsqrt( dsqrt( uflxsfx*uflxsfx +
     &                      vflxsfx*vflxsfx ) / rhox2d(i,j) )

c convert surface fluxes to kinematic units
       xhfx(i,j) = hfx(i,j)/(cp*rhox2d(i,j))
       xqfx(i,j) = qfx(i,j)/rhox2d(i,j)
c compute virtual heat flux at surface
       hfxv(i,j) = xhfx(i,j) + 0.61*thx3d(i,kx,j)*xqfx(i,j)
190   continue
c
c estimate potential temperature at 10m via log temperature profile
c in the surface layer (brutsaert, p. 63).
c calculate mixing ratio at 10m by assuming a constant
c value from the surface to the lowest model level.
c

      do 200 i=2,ilx
       sh10 = qvb3d(i,kx,j)/(qvb3d(i,kx,j)+1)
c      th10(i,j) = ((thx3d(i,kx,j)+tgb(i,j))/2.0)*(1.0+0.61*sh10)
c      th10(i,j) = thvx(i,kx,j) + hfxv(i,j)/(karman*ustr(i,j))
c     1            *dlog(za(i,kx,j)/10.)

c "virtual" potential temperature
       if (hfxv(i,j).ge. 0.) then
        th10(i,j) = thvx(i,kx,j)
       else
c       th10(i,j) = (0.25*thx3d(i,kx,j)+0.75*tgb(i,j))*(1.0+0.61*sh10)
c----   first approximation for obhukov length
        oblen=-0.5*(thx3d(i,kx,j)
     &        +tgb(i,j))*(1.0+0.61*sh10)*ustr(i,j)**3/
     &      (g * karman * (hfxv(i,j) + dsign(1.d-10,hfxv(i,j))))
        if (oblen.ge.za(i,kx,j)) then
         th10(i,j) = thvx(i,kx,j) + hfxv(i,j)/(karman*ustr(i,j))
     &            *(dlog(za(i,kx,j)/10.)+5./oblen*(za(i,kx,j)-10.))
        elseif (oblen.lt.za(i,kx,j).and.oblen.gt.10.) then
         th10(i,j) = thvx(i,kx,j) + hfxv(i,j)/(karman*ustr(i,j))*
     &   (dlog(oblen/10.)+5./oblen*(oblen-10.)
     &                   +6*dlog(za(i,kx,j)/oblen))
        elseif (oblen.le.10.) then
         th10(i,j) = thvx(i,kx,j) + hfxv(i,j)/(karman*ustr(i,j))*
     &   6*dlog(za(i,kx,j)/10.)
        endif
        th10(i,j)=dmax1(th10(i,j),tgb(i,j))
       endif
cgtb   th10(i,j) = dmin1(th10(i,j),tgb(i,j))  ! gtb add to minimize obklen

c       compute obukhov length
       obklen(i,j) = - th10(i,j) * ustr(i,j)**3 /
     &              (g * karman * (hfxv(i,j) + dsign(1.d-10,hfxv(i,j))))
200   continue
c
c compute diffusivities and counter gradient terms
c
 1079 continue
      call blhnew
#ifdef MPP1
      do j=jbegin,jendx
        if((myid.ne.nproc-1).or.
     &     (myid.eq.nproc-1.and.j.lt.jendx)) then
#else
      do j=2,jlx-1
#endif
          do k=1,kx
          do i=2,ilx
        if(k.gt.1) akxx1(i,k,j)=rhohf(i,k-1,j)*kvm(i,k,j)/dza(i,k-1,j)
        akxx2(i,k,j)=g/(psb(i,j)*1000.)/dsigma(k)
          enddo
          enddo
#ifdef MPP1
        endif
#endif
      enddo
#ifdef MPP1
      ii=0
      do k=1,kx
      do i=2,ilx
        ii=ii+1
        wksend(ii) = akxx1(i,k,jxp)
      enddo
      enddo
      do k=1,kx
      do i=2,ilx
        ii=ii+1
        wksend(ii) = akxx2(i,k,jxp)
      enddo
      enddo
      call MPI_SENDRECV(wksend(1),(ix-2)*kx*2,MPI_DOUBLE_PRECISION,
     &          ieast,1,wkrecv(1),(ix-2)*kx*2,MPI_DOUBLE_PRECISION,
     &          iwest,1,MPI_COMM_WORLD,status,ierr)
      ii=0
      do k=1,kx
      do i=2,ilx
        ii=ii+1
        akxx1(i,k,0) = wkrecv(ii)
      enddo
      enddo
      do k=1,kx
      do i=2,ilx
        ii=ii+1
        akxx2(i,k,0) = wkrecv(ii)
      enddo
      enddo
#endif

#ifdef MPP1
      do 1080 j=jbegin,jendx
#else
      do 1080 j=2,jlx
#endif

c    calculate coefficients at dot points for u and v wind

#ifdef MPP1
      if ( myid .eq. 0.and. j .eq. 2) then
#else
      if(j.eq.2) then
#endif
       do k=1,kx
       do i=2,ilx
        idx = i
        idx = min0(idx,ilx-1)
        idxm1 = i-1
        idxm1 = max0(idxm1,2)
        if (k.gt.1) betak(i,k)=0.5*(akxx1(idx,k,j)+akxx1(idxm1,k,j))
        alphak(i,k)=0.5*(akxx2(idx,k,j)+akxx2(idxm1,k,j))
       enddo
       enddo
#ifdef MPP1
      else if(myid.eq.nproc-1 .and. j.eq.jendx) then
#else
      else if(j.eq.jlx) then
#endif
       do k=1,kx
       do i=2,ilx
        idx = i
        idx = min0(idx,ilx-1)
        idxm1 = i-1
        idxm1 = max0(idxm1,2)
        if (k.gt.1) betak(i,k)=0.5*(akxx1(idx,k,j-1)
     &                             +akxx1(idxm1,k,j-1))
        alphak(i,k)=0.5*(akxx2(idx,k,j-1)+akxx2(idxm1,k,j-1))
       enddo
       enddo
      else
       do k=1,kx
       do i=2,ilx
        idx = i
        idx = min0(idx,ilx-1)
        idxm1 = i-1
        idxm1 = max0(idxm1,2)
        if (k.gt.1)
     &  betak(i,k)=0.25*(akxx1(idx,k,j-1)+akxx1(idxm1,k,j-1)
     &                    +akxx1(idx,k,j)  +akxx1(idxm1,k,j))
        alphak(i,k)=0.25*(akxx2(idx,k,j-1)+akxx2(idxm1,k,j-1)
     &                     +akxx2(idx,k,j)  +akxx2(idxm1,k,j))
       enddo
       enddo
      endif
c
c *************************************************************************
c
c   start now procedure for implicit diffusion calculations
c   performed separately for wind (dot points)
c   and temperature and water vapor (cross points)
c   countergradient term is not included in the implicit diffusion scheme
c   its effect is included as in the old explicit scheme
c   calculations assume fluxes positive upward, so the sign in front
c   of uflxsf and vflxsf has been changed in the various terms
c

c
c   wind components
c

c
c    first compute coefficients of the tridiagonal matrix
c

      do 220 k=2,kx-1
      do 220 i=2,ilx
       coef1(i,k)=dt*alphak(i,k)*betak(i,k+1)
       coef2(i,k)=1.+dt*alphak(i,k)*(betak(i,k+1)+betak(i,k))
       coef3(i,k)=dt*alphak(i,k)*betak(i,k)
220   continue

      do 221 i=2,ilx
       coef1(i,1)=dt*alphak(i,1)*betak(i,2)
       coef2(i,1)=1.+dt*alphak(i,1)*betak(i,2)
       coef3(i,1)=0.
       coef1(i,kx)=0.
       coef2(i,kx)=1.+dt*alphak(i,kx)*betak(i,kx)
       coef3(i,kx)=dt*alphak(i,kx)*betak(i,kx)
221   continue

      do 224 i=2,ilx
       coefe(i,1)=coef1(i,1)/coef2(i,1)
       coeff1(i,1)=auxx(i,1,j)/coef2(i,1)
       coeff2(i,1)=avxx(i,1,j)/coef2(i,1)
224   continue

      do 225 k=2,kx-1
      do 225 i=2,ilx
       coefe(i,k)=coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
       coeff1(i,k)=(auxx(i,k,j)+coef3(i,k)*coeff1(i,k-1))/(coef2(i,k)-
     &           coef3(i,k)*coefe(i,k-1))
       coeff2(i,k)=(avxx(i,k,j)+coef3(i,k)*coeff2(i,k-1))/(coef2(i,k)-
     &           coef3(i,k)*coefe(i,k-1))
225   continue

      do 226 i=2,ilx
       idx = i
       idx = min0(idx,ilx)
       idxm1 = i-1
       idxm1 = max0(idxm1,2)
       jdx = j
#ifdef MPP1
       if (myid.eq.nproc-1) then
         jdx = min0(jdx,jendx)
       endif
#else
       jdx = min0(jdx,jlx)
#endif
       jdxm1 = j-1
#ifdef MPP1
       if (myid .eq. 0) then
         jdxm1 = max0(jdxm1,2)
       endif
#else
       jdxm1 = max0(jdxm1,2)
#endif
       drgdot = 0.25 * ( uvdrag(idxm1,jdxm1) + uvdrag(idxm1,jdx)
     &                + uvdrag(idx  ,jdxm1) + uvdrag(idx  ,jdx) )
       uflxsf = drgdot * auxx(i,kx,j)
       vflxsf = drgdot * avxx(i,kx,j)

       coefe(i,kx)=0.
       coeff1(i,kx)=(auxx(i,kx,j)-dt*alphak(i,kx)*uflxsf+coef3(i,kx)
     &     *coeff1(i,kx-1))/(coef2(i,kx)-coef3(i,kx)*coefe(i,kx-1))
       coeff2(i,kx)=(avxx(i,kx,j)-dt*alphak(i,kx)*vflxsf+coef3(i,kx)
     &     *coeff2(i,kx-1))/(coef2(i,kx)-coef3(i,kx)*coefe(i,kx-1))

226   continue
c
c   all coefficients have been computed, predict field and put it in
c   temporary work space tpred
c
      do 229 i=2,ilx
       tpred1(i,kx)=coeff1(i,kx)
       tpred2(i,kx)=coeff2(i,kx)
229   continue

      do 230 k=kx-1,1,-1
      do 230 i=2,ilx
       tpred1(i,k)=coefe(i,k)*tpred1(i,k+1)+coeff1(i,k)
       tpred2(i,k)=coefe(i,k)*tpred2(i,k+1)+coeff2(i,k)
230   continue

c
c   calculate tendency due to vertical diffusion using temporary
c   predicted field
c
      do 240 k=1,kx
      do 240 i=2,ilx
       dumr=0.25*(psb(i,j)+psb(i,j-1)+psb(i-1,j)+psb(i-1,j-1))
       uten(i,k,j)=uten(i,k,j)+(tpred1(i,k)-auxx(i,k,j))/dt*dumr
       vten(i,k,j)=vten(i,k,j)+(tpred2(i,k)-avxx(i,k,j))/dt*dumr
240   continue

c     temperature
c

c    calculate coefficients at cross points for temperature

      do 250 k=1,kx
      do 250 i=2,ilx
       if (k.gt.1) betak(i,k)=rhohf(i,k-1,j)*kvh(i,k,j)/dza(i,k-1,j)
       alphak(i,k)=g/(psb(i,j)*1000.)/dsigma(k)
250   continue

      do 260 k=2,kx-1
      do 260 i=2,ilx
       coef1(i,k)=dt*alphak(i,k)*betak(i,k+1)
       coef2(i,k)=1.+dt*alphak(i,k)*(betak(i,k+1)+betak(i,k))
       coef3(i,k)=dt*alphak(i,k)*betak(i,k)
260   continue

      do 261 i=2,ilx
       coef1(i,1)=dt*alphak(i,1)*betak(i,2)
       coef2(i,1)=1.+dt*alphak(i,1)*betak(i,2)
       coef3(i,1)=0.
       coef1(i,kx)=0.
       coef2(i,kx)=1.+dt*alphak(i,kx)*betak(i,kx)
       coef3(i,kx)=dt*alphak(i,kx)*betak(i,kx)
261   continue

      do 264 i=2,ilx
       coefe(i,1)=coef1(i,1)/coef2(i,1)
       coeff1(i,1)=thx3d(i,1,j)/coef2(i,1)
264   continue

      do 265 k=2,kx-1
      do 265 i=2,ilx
       coefe(i,k)=coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
       coeff1(i,k)=(thx3d(i,k,j)+coef3(i,k)*coeff1(i,k-1))/(coef2(i,k)-
     &           coef3(i,k)*coefe(i,k-1))
265   continue

      do 267 i=2,ilx
       coefe(i,kx)=0.
       coeff1(i,kx)=(thx3d(i,kx,j)
     &             +dt*alphak(i,kx)*hfx(i,j)/cp+coef3(i,kx)
     &    *coeff1(i,kx-1))/(coef2(i,kx)-coef3(i,kx)*coefe(i,kx-1))
267   continue

c
c   all coefficients have been computed, predict field and put it in
c   temporary work space tpred
c

      do 269 i=2,ilx
       tpred1(i,kx)=coeff1(i,kx)
269   continue

      do 270 k=kx-1,1,-1
      do 270 i=2,ilx
       tpred1(i,k)=coefe(i,k)*tpred1(i,k+1)+coeff1(i,k)
270   continue

c
c   calculate tendency due to vertical diffusion using temporary
c   predicted field
c
      do 280 k=1,kx
      do 280 i=2,ilx
       sf=tb(i,k,j)/thx3d(i,k,j)
       difft(i,k,j)=difft(i,k,j)+(tpred1(i,k)-thx3d(i,k,j))/dt*sf
280   continue
c
c   water vapor
c

c    calculate coefficients at cross points for water vapor

      do 290 k=1,kx
      do 290 i=2,ilx
       if (k.gt.1) betak(i,k)=rhohf(i,k-1,j)*kvq(i,k,j)/dza(i,k-1,j)
       alphak(i,k)=g/(psb(i,j)*1000.)/dsigma(k)
290   continue

      do 300 k=2,kx-1
      do 300 i=2,ilx
       coef1(i,k)=dt*alphak(i,k)*betak(i,k+1)
       coef2(i,k)=1.+dt*alphak(i,k)*(betak(i,k+1)+betak(i,k))
       coef3(i,k)=dt*alphak(i,k)*betak(i,k)
300   continue

      do 301 i=2,ilx
       coef1(i,1)=dt*alphak(i,1)*betak(i,2)
       coef2(i,1)=1.+dt*alphak(i,1)*betak(i,2)
       coef3(i,1)=0.
       coef1(i,kx)=0.
       coef2(i,kx)=1.+dt*alphak(i,kx)*betak(i,kx)
       coef3(i,kx)=dt*alphak(i,kx)*betak(i,kx)
301   continue

      do 304 i=2,ilx
       coefe(i,1)=coef1(i,1)/coef2(i,1)
       coeff1(i,1)=qvb3d(i,1,j)/coef2(i,1)
304   continue

      do 305 k=2,kx-1
      do 305 i=2,ilx
       coefe(i,k)=coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
       coeff1(i,k)=(qvb3d(i,k,j)+coef3(i,k)*coeff1(i,k-1))/(coef2(i,k)-
     &           coef3(i,k)*coefe(i,k-1))
305   continue

      do 306 i=2,ilx
       coefe(i,kx)=0.
       coeff1(i,kx)=(qvb3d(i,kx,j)+dt*alphak(i,kx)*qfx(i,j)+coef3(i,kx)*
     &     coeff1(i,kx-1))/(coef2(i,kx)-coef3(i,kx)*coefe(i,kx-1))
306   continue

c
c   all coefficients have been computed, predict field and put it in
c   temporary work space tpred

      do 309 i=2,ilx
       tpred1(i,kx)=coeff1(i,kx)
309   continue

      do 310 k=kx-1,1,-1
      do 310 i=2,ilx
       tpred1(i,k)=coefe(i,k)*tpred1(i,k+1)+coeff1(i,k)
310   continue

c
c   calculate tendency due to vertical diffusion using temporary
c   predicted field
c
      do 320 k=1,kx
      do 320 i=2,ilx
       diffq(i,k,j)=diffq(i,k,j)+(tpred1(i,k)-qvb(i,k,j)/psb(i,j))/
     &           dt*psb(i,j)
320   continue

c    calculate coefficients at cross points for cloud vater

       do 330 k=1,kx
       do 330 i=2,ilx
        if (k.gt.1) betak(i,k)=rhohf(i,k-1,j)*kvq(i,k,j)/dza(i,k-1,j)
        alphak(i,k)=g/(psb(i,j)*1000.)/dsigma(k)
330    continue

       do 340 k=2,kx-1
       do 340 i=2,ilx
        coef1(i,k)=dt*alphak(i,k)*betak(i,k+1)
        coef2(i,k)=1.+dt*alphak(i,k)*(betak(i,k+1)+betak(i,k))
        coef3(i,k)=dt*alphak(i,k)*betak(i,k)
340    continue

       do 341 i=2,ilx
        coef1(i,1)=dt*alphak(i,1)*betak(i,2)
        coef2(i,1)=1.+dt*alphak(i,1)*betak(i,2)
        coef3(i,1)=0.
        coef1(i,kx)=0.
        coef2(i,kx)=1.+dt*alphak(i,kx)*betak(i,kx)
        coef3(i,kx)=dt*alphak(i,kx)*betak(i,kx)
341    continue

       do 344 i=2,ilx
        coefe(i,1)=coef1(i,1)/coef2(i,1)
        coeff1(i,1)=qcx(i,1,j)/coef2(i,1)
344    continue

       do 345 k=2,kx-1
       do 345 i=2,ilx
        coefe(i,k)=coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
        coeff1(i,k)=(qcx(i,k,j)+coef3(i,k)*coeff1(i,k-1))/(coef2(i,k)-
     &           coef3(i,k)*coefe(i,k-1))
345    continue

       do 346 i=2,ilx
        coefe(i,kx)=0.
        coeff1(i,kx)=(qcx(i,kx,j)+coef3(i,kx)*
     &     coeff1(i,kx-1))/(coef2(i,kx)-coef3(i,kx)*coefe(i,kx-1))
346    continue

c
c   all coefficients have been computed, predict field and put it in
c   temporary work space tpred
c

       do 349 i=2,ilx
        tpred1(i,kx)=coeff1(i,kx)
349    continue

       do 350 k=kx-1,1,-1
       do 350 i=2,ilx
        tpred1(i,k)=coefe(i,k)*tpred1(i,k+1)+coeff1(i,k)
350    continue

c
c   calculate tendency due to vertical diffusion using temporary
c   predicted field
c
       do 360 k=1,kx
       do 360 i=2,ilx
        qcten(i,k,j)=qcten(i,k,j)+(tpred1(i,k)-qcb(i,k,j)/psb(i,j))/dt
     &  *psb(i,j)
360    continue

c
c ***********************************************************************
c
c   now add countergradient term to temperature and water vapor equation
c
Ctrapuv
      do i=2,ilx
        ttnp(i,1) =0.0d0
      enddo
Ctrapuv_
      do 380 k=2,kx
      do 380 i=2,ilx
       sf=tb(i,k,j)/(psb(i,j)*thx3d(i,k,j))
       ttnp(i,k)=sf*cp*rhohf(i,k-1,j)*kvh(i,k,j)*cgh(i,k,j)
380   continue
c
c-----compute the tendencies:
c
      do 400 i=2,ilx
       difft(i,kx,j)=difft(i,kx,j)-g*ttnp(i,kx)/(1000.*cp*dsigma(kx))
400   continue
c
      do 420 k=1,kxm
      do 420 i=2,ilx
       difft(i,k,j)=difft(i,k,j)+g*(ttnp(i,k+1)-ttnp(i,k))/(1000.*cp*
     &                    dsigma(k))
420   continue

c
Cchem2
      if(ichem.eq.1 .and. ichdrdepo .eq.1) then
c
c   coef1, coef2, coef3 and coefe are the same as for water vapor
c   and cloud water so they do not need to be recalculated

c
c    recalculation of coef1,2,3  with tracer diffusivity kvc 

       do 630 k=1,kx
       do 630 i=2,ilx
        if (k.gt.1) betak(i,k)=rhohf(i,k-1,j)*kvc(i,k,j)/dza(i,k-1,j)
        alphak(i,k)=g/(psb(i,j)*1000.)/dsigma(k)
630    continue

       do 640 k=2,kx-1
       do 640 i=2,ilx
        coef1(i,k)=dt*alphak(i,k)*betak(i,k+1)
        coef2(i,k)=1.+dt*alphak(i,k)*(betak(i,k+1)+betak(i,k))
        coef3(i,k)=dt*alphak(i,k)*betak(i,k)
640    continue

       do 641 i=2,ilx
        coef1(i,1)=dt*alphak(i,1)*betak(i,2)
        coef2(i,1)=1.+dt*alphak(i,1)*betak(i,2)
        coef3(i,1)=0.
        coef1(i,kx)=0.
        coef2(i,kx)=1.+dt*alphak(i,kx)*betak(i,kx)
        coef3(i,kx)=dt*alphak(i,kx)*betak(i,kx)
641    continue
c
c set the Vd for in case of prescribed deposition velocities
c

       do 430 itr=1,ntr
       do 430 i=2,ilx         
        if(veg2d(i,j).le.0.00001)then
         vdep(i,itr) = chtrdpv(itr,2)
        else
         vdep(i,itr) = chtrdpv(itr,1)
        endif
c provisoire test de la routine chdrydep pour les dust

         if ( chtrname(itr) .eq. 'DUST' ) then
           vdep(i,itr) = 0. 
         end if 
430    continue
c
       do 500 itr=1,ntr
c
       do 440 k=1,kx
       do 440 i=2,ilx
        chix(i,k)=chib(i,k,j,itr)/psb(i,j)
440    continue
c
       do 450 i=2,ilx
        coefe(i,1)=coef1(i,1)/coef2(i,1)
        coeff1(i,1)=chix(i,1)/coef2(i,1)
450    continue
c
       do 460 k=2,kx-1
       do 460 i=2,ilx
        coefe(i,k)=coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
        coeff1(i,k)=(chix(i,k)+coef3(i,k)*coeff1(i,k-1))/(coef2(i,k)-
     1           coef3(i,k)*coefe(i,k-1))
460    continue

       do 465 i=2,ilx
        coefe(i,kx)=0.

c add dry deposition option1
        coeff1(i,kx)=(chix(i,kx)-
     1     dt*alphak(i,kx)*chix(i,kx)*vdep(i,itr)*rhox2d(i,j)
     1     +coef3(i,kx)*
     1     coeff1(i,kx-1))/(coef2(i,kx)-coef3(i,kx)*coefe(i,kx-1))
465    continue
c
c   all coefficients have been computed, predict field and put it in
c   temporary work space tpred1
c
       do 470 i=2,ilx
        tpred1(i,kx)=coeff1(i,kx)
470    continue
c
       do 480 k=kx-1,1,-1
       do 480 i=2,ilx
        tpred1(i,k)=coefe(i,k)*tpred1(i,k+1)+coeff1(i,k)
480    continue
c
c
c   calculate tendency due to vertical diffusion using temporary
c   predicted field
c Dry deposition option 1 is included 

c
       do 490 k=1,kx
       do 490 i=2,ilx
cqian          chiten(i,k,j,itr)=chiten(i,k,j,itr)
CCGAFFE TEST diffusion/10
        chiten(i,k,j,itr)=chiten(i,k,j,itr)+(tpred1(i,k)-
     1  chix(i,k))/dt *psb(i,j)
c         chiten(i,k,j,itr)=chiten(i,k,j,itr)+0.1 *(tpred1(i,k)-
c     1  chix(i,k))/dt *psb(i,j)

490    continue
      do 502 i=2,ilx

      if ( chtrname (itr).ne. 'DUST') then
      remdrd(i,j,itr)=remdrd(i,j,itr)+chix(i,kx)*vdep(i,itr)
     1 * psb(i,j)*dt/2. *rhox2d(i,j)*g/ (psb(i,j)*1000.*dsigma(kx))
      end if

502   continue
500   continue
      endif
Cchem2_

 1080 continue
c
      return
      end
