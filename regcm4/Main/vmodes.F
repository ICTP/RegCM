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

      subroutine vmodes( lstand, sigmaf, kv1)
      implicit none
      integer kv1
      logical lprint, lhydro, lstand, lsigma
c
c  this subroutine determines the vertical modes of the psu/ncar meso-
c  scale model designated mm4.  it also computes associated transform
c  matrices used by the initialization software used with mm4.
c
c----------------------------------------------------------------------
c
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/split.cb"
#ifdef MPP1
#     include "mpiregcm.inc"
      integer ierr
#endif
c
c   the following are user set parameters:
c
c ix,jx  = dimension of horizontal grid (later called ni,nj).
c          as in mm4, ix is for n-s direction, jx for w-e direction.
c kx     = number of model data levels
c dt     = time step used to generate tendencies (later called delt)
c dx     = grid spacing at center in meters (later called delx).
c clat   = latitude of central point, used to determine coriolis param.
c ptop   = model top in units of cb (later called pt).
c
      integer nk1,nverts
      parameter (nk1=kx+1, nverts=7+8*kx+7*kx*kx+nk1*nk1)
c
c     programmed by ronald m. errico at ncar,  dec 1984.
c     revised by ronald errico and gary bates, nov 1987.
c     revised by ronald errico,                mar 1988.
c     for further info see: ncar tech note by errico and bates, 1988.
c
c  iunit is the output unit number for file of eigenvectors, etc.
c  lstand = .true. if standard atmosphere t to be used (ignore input
c            tbarh and ps in that case).  otherwise, ps and tbarh must
c            be defined on input.  note that in either case, pt must
c            also be defined on input (common block named cvert).
c
      real(kind=8)  rgas,xkappa
      parameter ( rgas=287., xkappa=.287 )
c
c  rgas and xkappa are atmospheric values for r and kappa.
c
      real(kind=8)  sigmaf(kv1)
      data lprint /.false./ ! true if all matrices to be printed
c
      real(kind=8)  work(kx)             !for suboutine invmtrx
      integer numerr,k,l,k1,mm,i,j,ier
      real(kind=8)  rm1,det(2),ps2,x
c
      numerr = 0
      lprint = .false.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c                            s  t  a  r  t
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c       set arrays describing vertical structure
c
c  set reference pressures
      if (lstand) ps=100.   ! standard ps in cb; otherwise ps set in tav
      pd=ps-pt
      r=rgas
c
c  read sigmaf (sigma at full (integral) index levels; kx+1 values as
c  in the mm4.   check that values are ordered properly.
c
#ifdef MPP1
      if (myid .eq. 0) then
#endif
      if (lstand) print 700
      if (.not.lstand) print 701
#ifdef MPP1
      end if
#endif
c
      lsigma=.false.
      if (sigmaf(1)  .ne.0.) lsigma=.true.
      if (sigmaf(nk1).ne.1.) lsigma=.true.
      do 1 k=1,kx
      if (sigmaf(k+1).le.sigmaf(k)) then
         lsigma=.true.
#ifdef MPP1
         if(myid .eq. 0) print 2,k,sigmaf(k+1),sigmaf(k)
#else
         print 2,k,sigmaf(k+1),sigmaf(k)
#endif
      end if
    1 continue
      if (lsigma) then
#ifdef MPP1
         if (myid .eq. 0) then
         write(*,*) 'sigma values in list vmode inappropriate'
         end if
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
         write(*,*) 'sigma values in list vmode inappropriate'
         stop
#endif
      end if
    2 format('0 for k=',i3,' sigmaf(k+1)=',f9.6,' .le. sigmaf(k)=',f9.6)
c
c  compute sigmah (sigma at half levels) and delta sigma
      do 10 k=1,kx
      sigmah(k) = 0.5 * (sigmaf(k)+sigmaf(k+1))
      dsigma(k) = sigmaf(k+1) - sigmaf(k)
   10 continue
      sigmah(nk1) = 1.0
c
c  set tbarh (temperature at half (data) levels: indexed k + 1/2)
      if (lstand) call vtlaps(tbarh,sigmah,r,pt,pd,kx)
      call vchekt(tbarh,sigmah,sigmaf,xkappa,pt,pd,kx,numerr)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c      determine thermodynamic matrix
c
c  compute thetah
      do 100 k=1,kx
      thetah(k)=tbarh(k)*((sigmah(k)+pt/pd)**(-xkappa))
  100 continue
c
c  compute tbarf and thetaf
      do 101 k=2,kx
      k1=k-1
      tbarf(k) = tbarh(k1)*(sigmah(k)-sigmaf(k))/(sigmah(k)-sigmah(k1))
     1         + tbarh(k)*(sigmaf(k)-sigmah(k1))/(sigmah(k)-sigmah(k1))
  101 continue
      tbarf(  1) = 0.
      tbarf(nk1) = 0.
c
      do 102 k=1,nk1
  102 thetaf(k) = tbarf(k) * ((sigmaf(k)+pt/pd)**(-xkappa))
c
c  define matrices for determination of thermodynamic matrix
      do 103 l=1,kx
      do 103 k=1,kx
      if (l.gt.k) e2(k,l)=0.
      if (l.le.k) e2(k,l)=1.
      e1(k,l)=1.
  103 continue
c
      do 104 l=1,kx
      do 104 k=1,kx
      a3(k,l)=0.
      d1(k,l)=0.
      d2(k,l)=0.
      s1(k,l)=0.
      s2(k,l)=0.
      x1(k,l)=0.
  104 continue
      do 105 k=1,kx
      a3(k,k)=-tbarh(k)
      d1(k,k)=sigmaf(k+1)-sigmaf(k)
      d2(k,k)=xkappa*tbarh(k)/(sigmah(k)+pt/pd)
      s1(k,k)=sigmaf(k)
      s2(k,k)=sigmah(k)
      x1(k,k)=1.
  105 continue
c
      do 107 k=1,kx
      do 106 l=1,kx
      e3(k,l)=0.
      g1(k,l)=0.
  106 continue
      e3(k,k)=1.
      if (k.gt.1)  g1(k,k)=tbarf(k)
      if (k.lt.kx) g1(k,k+1)=-tbarf(k+1)
      if (k.lt.kx) e3(k,k+1)=1.
  107 continue
c
c  compute g2 (i.e., the transform from divg. to sigma dot)
      call vsubtm(w1,e2,x1,kx)
      call vmultm(w2,w1,d1,kx)
      call vmultm(g2,e1,d1,kx)
      call vmultm(w1,s1,g2,kx)
      call vsubtm(g2,w1,w2,kx)
c
c  compute a1
      do 201 k=1,kx
      do 200 l=1,kx
  200 w2(k,l)=0.
  201 w2(k,k)=1./d1(k,k)
      call vmultm(w1,g1,g2,kx)
      call vmultm(a1,w2,w1,kx)
c
c  compute a2
      call vmultm(w1,e1,d1,kx)
      call vmultm(a2,s2,w1,kx)
      call vmultm(w2,e3,g2,kx)
      call vmultc(w2,w2,kx,0.5d0)
      call vsubtm(w1,w2,a2,kx)
      call vmultm(a2,d2,w1,kx)
c
c  compute a4
      call vmultm(w1,e1,d1,kx)
      call vmultm(a4,a3,w1,kx)
      call vmultc(a4,a4,kx,-1.d0)
c
c  compute a
      call vaddms(a,a1,a2,kx)
      call vaddms(a,a ,a3,kx)
      call vaddms(a,a ,a4,kx)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c       determine matrices for linearized determination of geopotential
c
c  compute delta log p
      do 300 k=2,kx
      w1(k,1) = dlog((sigmah(k)+pt/pd)/(sigmah(k-1)+pt/pd))
  300 continue
c
c  compute matrix which multiples t vector
      do 301 l=1,kx
      do 301 k=1,kx
      hydros(k,l)=0.
  301 continue
c
      do 302 k=1,kx-1
      do 302 l=k,kx-1
      hydros(k,l  ) = hydros(k,l  ) + w1(l+1,1)*dsigma(l  ) /
     a                (dsigma(l+1)+dsigma(l))
      hydros(k,l+1) = hydros(k,l+1) + w1(l+1,1)*dsigma(l+1) /
     b                (dsigma(l+1)+dsigma(l))
  302 continue
c
      do 303 k=1,kx
      hydros(k,kx) = hydros(k,kx)
     a             + dlog((1.+pt/pd)/(sigmah(kx)+pt/pd))
  303 continue
c
c  compute matirx which multiplies log(sigma*p+pt) vector
      do 305 l=1,nk1
      do 305 k=1,kx
  305 hydroc(k,l) = 0.
c
      tweigh(1) = 0.
      do 306 l=2,kx
  306 tweigh(l) = (tbarh(l)*dsigma(l)+tbarh(l-1)*dsigma(l-1))
     1        / (dsigma(l)+dsigma(l-1))
c
      do 307 l=2,kx-1
      do 307 k=1,l-1
  307 hydroc(k,l) = tweigh(l) - tweigh(l+1)
c
      do 308 l=1,kx-1
  308 hydroc(l,l) = tbarh(l) - tweigh(l+1)
c
      do 309 k=1,kx-1
  309 hydroc(k,kx) = tweigh(kx) - tbarh(kx)
c
      do 310 k=1,kx
  310 hydroc(k,nk1) = tbarh(kx)
c
c  test hydroc and hydros matrices (if correct, w1(k,1)=w1(k,2))
      lhydro=.false.
      do 313 k=1,kx
      w1(k,1) = 0.
      do 311 l=1,kx
      w1(k,1) = w1(k,1) + hydros(k,l) * tbarh(l)
  311 continue
      w1(k,2) = - tbarh(k) * dlog(sigmah(k)*pd+pt)
      do 312 l=1,nk1
      w1(k,2) = w1(k,2) + hydroc(k,l) * dlog(sigmah(l)*pd+pt)
  312 continue
      x=dabs(w1(k,1)-w1(k,2)) / (dabs(w1(k,1))+dabs(w1(k,2)))
      if (x.gt.1.e-8) lhydro=.true.
  313 continue
c
      if (lhydro) then
         numerr=numerr+1
         print 315
         call vprntv(w1(1,1),kx,'test1   ')
         call vprntv(w1(1,2),kx,'test2   ')
      end if
  315 format('0 problem with linearization of hydostatic equation')
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c       determine tau matrix
c
      do 400 l=1,kx
      do 400 k=1,nk1
      w3(k,l) = dsigma(l) / ( 1. + pt/(pd*sigmah(k)) )
  400 continue
c
      do 401 l=1,kx
      do 401 k=1,kx
      w2(k,l) = 0.
      do 401 mm=1,nk1
      w2(k,l) = w2(k,l) + hydroc(k,mm) * w3(mm,l)
  401 continue
c
      call vmultm(w1,hydros,a,kx)
      call vsubtm(tau,w1,w2,kx)
      rm1=-1.*r
      call vmultc(tau,tau,kx,rm1)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c       determine other matrices and vectors
c
c  compute eigenvalues and vectors for tau (rg calls eispack routines)
      call vmultc(w1,tau,kx,1.d0) ! copy tau since rg destroys input
      call rg(kx,kx,w1,hbar,w2,1,zmatx,iw2,w2(1,3),ier)
      call vcheki(ier,numerr,'zmatx   ')
      call vcheke(hbar,w2,kx,numerr,'tau     ')
      call vorder(zmatx,hbar,w1,w2,kx)
      call vnorml(zmatx,sigmaf,kx,nk1)
c
c  compute inverse of zmatx
      call invmtrx(zmatx,kx,zmatxr,kx,kx,det,iw2,ier,work)
      call vcheki(ier,numerr,'zmatxr  ')
c
c  compute inverse of hydros
      call invmtrx(hydros,kx,hydror,kx,kx,det,iw2,ier,work)
      call vcheki(ier,numerr,'hydror  ')
c
c  compute cpfac
      call invmtrx(tau,kx,w1,kx,kx,det,iw2,ier,work)
      call vcheki(ier,numerr,'taur    ')
c
      do 450 k=1,kx
      cpfac(k)=0.
      do 450 l=1,kx
      cpfac(k) = cpfac(k) + (sigmaf(l+1)-sigmaf(l)) * w1(l,k)
  450 continue
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c       determine arrays needed for daley's variational scheme
c             for determination of surface pressure changes
c
      do 500 i=1,kx      ! weight of t for different levels
      hweigh(i)=0.
  500 continue
      hweigh(kx)=1.      ! only lowest sigma level t considered
c
      do 502 j=1,kx
      do 502 i=1,kx      ! compute b(-1t) w/tbar**2 b(-1)
      w1(i,j)=0.
      do 502 k=1,kx
      w1(i,j)=hydror(k,i)*hydror(k,j)*hweigh(k)/(tbarh(k)**2) + w1(i,j)
  502 continue
c
      ps2=ps*ps
      do 504 j=1,nk1
      do 504 i=1,kx
      varpa1(i,j)=0.
      do 504 k=1,kx
      varpa1(i,j)=varpa1(i,j)+w1(i,k)*hydroc(k,j)*ps2
  504 continue
c
      do 506 j=1,nk1
      do 506 i=1,nk1
      varpa2(i,j)=0.
      do 506 k=1,kx
      varpa2(i,j)=varpa2(i,j)+hydroc(k,i)*varpa1(k,j)
  506 continue
c
      alpha1=hydros(kx,kx)*tbarh(kx)/ps
      alpha2=hweigh(kx)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c       output desired arrays
c
      call vprntv(sigmaf ,nk1,'sigmaf  ')
      call vprntv(tbarh  ,kx ,'t mean  ')
      call vprntv(ps     ,  1,'ps mean ')
      print 702,kx,numerr
c
c  printout if desired
      if (.not.lprint) return
      call vprntv(cpfac  ,kx ,'cpfac   ')
      call vprntv(dsigma ,kx ,'dsigma  ')
      call vprntv(hbar   ,kx ,'hbar    ')
      call vprntv(sigmah ,nk1,'sigmah  ')
      call vprntv(tbarf  ,nk1,'tbarf   ')
      call vprntv(thetah ,kx ,'thetah  ')
      call vprntv(thetaf ,nk1,'thetaf  ')
      call vprntv(hweigh ,kx ,'hweigh  ')
      print 703,alpha1,alpha2
      call vprntm(a     , kx, kx ,'a       ')
      call vprntm(hydros, kx, kx ,'hydros  ')
      call vprntm(hydror, kx, kx ,'hydror  ')
      call vprntm(hydroc, kx, nk1,'hydroc  ')
      call vprntm(tau   , kx, kx ,'tau     ')
      call vprntm(zmatx , kx, kx ,'zmatx   ')
      call vprntm(zmatxr, kx, kx ,'zmatxr  ')
      call vprntm(varpa1, kx, nk1,'varpa1  ')
      call vprntm(varpa2,nk1, nk1,'varpa2  ')
c
      return
  700 format('0 linearization about standard atmosphere (lstand=.t.)')
  701 format('0 linearization about horizontal mean of data ',
     a'(lstand=.f.)')
  702 format('0 vertical mode problem completed for kx=',i3,5x
     a,i1,' errors detected   (should be 0)')
  703 format('0alpha1 =',1p,1e16.5,'       alpha2 =',1p,1e16.5)
      end
