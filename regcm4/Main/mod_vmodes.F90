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
 
      module mod_vmodes

      private

      public :: vmodes

      contains

      subroutine vmodes(lstand,sigmaf,kv1)
!
      use mod_dynparam
      use mod_param1
      use mod_split
      use mod_constants , only : rgas , rovcp
      use mod_param3 , only : r8pt
      use mod_message
      implicit none
!
! Dummy arguments
!
      integer :: kv1
      logical :: lstand
      real(8) , dimension(kv1) :: sigmaf
      intent (in) kv1 , lstand
!
! Local variables
!
      real(8) , dimension(2) :: det
      integer :: ier , k , k1 , k2 , l , mm , numerr
      logical :: lhydro , lprint , lsigma
      real(8) :: ps2 , x
      real(8) , dimension(kz) :: work
      real(8) , dimension(1) :: pps
!
!  this subroutine determines the vertical modes of the psu/ncar meso-
!  scale model designated mm4.  it also computes associated transform
!  matrices used by the initialization software used with mm4.
!
!----------------------------------------------------------------------
!
!
!   the following are user set parameters:
!
! iy,jx  = dimension of horizontal grid (later called ni,nj).
!          as in mm4, iy is for n-s direction, jx for w-e direction.
! kz     = number of model data levels
! dt     = time step used to generate tendencies (later called delt)
! dx     = grid spacing at center in meters (later called delx).
! clat   = latitude of central point, used to determine coriolis param.
! r8pt   = model top in units of cb
!
!
!     programmed by ronald m. errico at ncar,  dec 1984.
!     revised by ronald errico and gary bates, nov 1987.
!     revised by ronald errico,                mar 1988.
!     for further info see: ncar tech note by errico and bates, 1988.
!
!  iunit is the output unit number for file of eigenvectors, etc.
!  lstand = .true. if standard atmosphere t to be used (ignore input
!            tbarh and ps in that case).  otherwise, ps and tbarh must
!            be defined on input.  note that in either case, r8pt must
!            also be defined on input (common block named cvert).
!
      data lprint/.false./  ! true if all matrices to be printed
!
      numerr = 0
      lprint = .false.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!                            s  t  a  r  t
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!       set arrays describing vertical structure
!
!  set reference pressures
      if ( lstand ) ps = 100.
                            ! standard ps in cb; otherwise ps set in tav
      pd = ps - r8pt
!
!  read sigmaf (sigma at full (integral) index levels; kzp1 values as
!  in the mm4.   check that values are ordered properly.
!
      write (aline,*) 'Calculating Vertical Modes'
      call say
      if ( lstand ) then
        write (aline,*) '- Linearization about standard atmosphere'
        call say
      else
        write (aline,*) '- Linearization about horizontal mean of data'
        call say
      end if
!
      lsigma = .false.
      if ( sigmaf(1).ne.0. ) lsigma = .true.
      if ( sigmaf(kzp1).ne.1. ) lsigma = .true.
      do k = 1 , kz
        if ( sigmaf(k+1).le.sigmaf(k) ) then
          lsigma = .true.
          write (aline,99001) k , sigmaf(k+1) , sigmaf(k)
99001     format ('0 for k=',i3,' sigmaf(k+1)=',f9.6,' .le. sigmaf(k)=',&
                & f9.6)
          call say
        end if
      end do
      if ( lsigma )                                                     &
      & call fatal(__FILE__,__LINE__,                                   &
      &          'Sigma values in list vmode inappropriate')
!
!  compute sigmah (sigma at half levels) and delta sigma
      do k = 1 , kz
        sigmah(k) = 0.5*(sigmaf(k)+sigmaf(k+1))
        dsigma(k) = sigmaf(k+1) - sigmaf(k)
      end do
      sigmah(kzp1) = 1.0
!
!  set tbarh (temperature at half (data) levels: indexed k + 1/2)
      if ( lstand ) call vtlaps(tbarh,sigmah,r8pt,pd,kz)
      call vchekt(tbarh,sigmah,sigmaf,rovcp,r8pt,pd,kz,numerr)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!      determine thermodynamic matrix
!
!  compute thetah
!
      thetah = tbarh*((sigmah+r8pt/pd)**(-rovcp))
!
!  compute tbarf and thetaf
!
      do k = 2 , kz
        k1 = k - 1
        tbarf(k) = tbarh(k1)*(sigmah(k)-sigmaf(k))                      &
                 & /(sigmah(k)-sigmah(k1)) + tbarh(k)                   &
                 & *(sigmaf(k)-sigmah(k1))/(sigmah(k)-sigmah(k1))
      end do
      tbarf(1) = 0.
      tbarf(kzp1) = 0.
!
      do k = 1 , kzp1
        if ( sigmaf(k).lt.1E-30 ) then
          thetaf(k) = tbarf(k)
        else
          thetaf(k) = tbarf(k)*((sigmaf(k)+r8pt/pd)**(-rovcp))
        end if
      end do
!
!  define matrices for determination of thermodynamic matrix
!
      do l = 1 , kz
        do k = 1 , kz
          if ( l.gt.k ) e2(k,l) = 0.
          if ( l.le.k ) e2(k,l) = 1.
          e1(k,l) = 1.
        end do
      end do
!
      a3 = 0.0
      d1 = 0.0
      d2 = 0.0
      s1 = 0.0
      s2 = 0.0
      x1 = 0.0
!
      do k = 1 , kz
        a3(k,k) = -tbarh(k)
        d1(k,k) = sigmaf(k+1) - sigmaf(k)
        d2(k,k) = rovcp*tbarh(k)/(sigmah(k)+r8pt/pd)
        s1(k,k) = sigmaf(k)
        s2(k,k) = sigmah(k)
        x1(k,k) = 1.
      end do
!
      do k = 1 , kz
        do l = 1 , kz
          e3(k,l) = 0.
          g1(k,l) = 0.
        end do
        e3(k,k) = 1.
        if ( k.gt.1 ) g1(k,k) = tbarf(k)
        if ( k.lt.kz ) g1(k,k+1) = -tbarf(k+1)
        if ( k.lt.kz ) e3(k,k+1) = 1.
      end do
!
!  compute g2 (i.e., the transform from divg. to sigma dot)
!
      w1 = e2 - x1
      w2 = matmul(w1,d1)
      g2 = matmul(e1,d1)
      w1 = matmul(s1,g2)
      g2 = w1 - w2
!
!  compute a1
!
      w2 = 0
      do k = 1 , kz
        w2(k,k) = 1.0/d1(k,k)
      end do
      w1 = matmul(g1,g2)
      a1 = matmul(w2,w1)
!
!  compute a2
!
      w1 = matmul(e1,d1)
      a2 = matmul(s2,w1)
      w2 = matmul(e3,g2)
      w2 = 0.5*w2
      w1 = w2-a2
      a2 = matmul(d2,w1)
!
!  compute a4
!
      w1 = matmul(e1,d1)
      a4 = matmul(a3,w1)
      a4 = -a4
!
!  compute a
!
      a = a1+a2+a3+a4
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!       determine matrices for linearized determination of geopotential
!
!  compute delta log p
      do k = 2 , kz
        w1(k,1) = dlog((sigmah(k)+r8pt/pd)/(sigmah(k-1)+r8pt/pd))
      end do
!
!  compute matrix which multiples t vector
!
      hydros = 0.0
!
      do k = 1 , kz - 1
        do l = k , kz - 1
          hydros(k,l) = hydros(k,l) + w1(l+1,1)*dsigma(l)               &
                      & /(dsigma(l+1)+dsigma(l))
          hydros(k,l+1) = hydros(k,l+1) + w1(l+1,1)*dsigma(l+1)         &
                        & /(dsigma(l+1)+dsigma(l))
        end do
      end do
!
      do k = 1 , kz
        hydros(k,kz) = hydros(k,kz)                                     &
                     & + dlog((1.+r8pt/pd)/(sigmah(kz)+r8pt/pd))
      end do
!
!  compute matirx which multiplies log(sigma*p+r8pt) vector
!
      hydroc = 0.0
!
      tweigh(1) = 0.
      do l = 2 , kz
        tweigh(l) = (tbarh(l)*dsigma(l)+tbarh(l-1)*dsigma(l-1))         &
                  & /(dsigma(l)+dsigma(l-1))
      end do
!
      do l = 2 , kz - 1
        do k = 1 , l - 1
          hydroc(k,l) = tweigh(l) - tweigh(l+1)
        end do
      end do
!
      do l = 1 , kz - 1
        hydroc(l,l) = tbarh(l) - tweigh(l+1)
      end do
!
      do k = 1 , kz - 1
        hydroc(k,kz) = tweigh(kz) - tbarh(kz)
      end do
!
      do k = 1 , kz
        hydroc(k,kzp1) = tbarh(kz)
      end do
!
!  test hydroc and hydros matrices (if correct, w1(k,1)=w1(k,2))
!
      lhydro = .false.
      do k = 1 , kz
        w1(k,1) = 0.
        do l = 1 , kz
          w1(k,1) = w1(k,1) + hydros(k,l)*tbarh(l)
        end do
        w1(k,2) = -tbarh(k)*dlog(sigmah(k)*pd+r8pt)
        do l = 1 , kzp1
          w1(k,2) = w1(k,2) + hydroc(k,l)*dlog(sigmah(l)*pd+r8pt)
        end do
        x = dabs(w1(k,1)-w1(k,2))/(dabs(w1(k,1))+dabs(w1(k,2)))
        if ( x.gt.1.E-8 ) lhydro = .true.
      end do
!
      if ( lhydro ) then
        numerr = numerr + 1
        print 99002
99002   format ('0 problem with linearization of hydostatic equation')
        call vprntv(w1(1,1),kz,'test1   ')
        call vprntv(w1(1,2),kz,'test2   ')
      end if
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!       determine tau matrix
!
      do l = 1 , kz
        do k = 1 , kzp1
          w3(k,l) = dsigma(l)/(1.+r8pt/(pd*sigmah(k)))
        end do
      end do
!
      do l = 1 , kz
        do k = 1 , kz
          w2(k,l) = 0.
          do mm = 1 , kzp1
            w2(k,l) = w2(k,l) + hydroc(k,mm)*w3(mm,l)
          end do
        end do
      end do
!
      w1 = matmul(hydros,a)
      tau = w1-w2
      tau = -rgas*tau
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!       determine other matrices and vectors
!
!  compute eigenvalues and vectors for tau (rg calls eispack routines)
!
      w1 = tau
      call rg(kz,w1,hbar,w2,1,zmatx,ier)
      call vcheki(ier,numerr,'zmatx   ')
      call vcheke(hbar,w2,kz,numerr,'tau     ')
      call vorder(zmatx,hbar,w1,w2,kz)
      call vnorml(zmatx,sigmaf,kz,kzp1)
!
!  compute inverse of zmatx
!
      call invmtrx(zmatx,kz,zmatxr,kz,kz,det,iw2,ier,work)
      call vcheki(ier,numerr,'zmatxr  ')
!
!  compute inverse of hydros
!
      call invmtrx(hydros,kz,hydror,kz,kz,det,iw2,ier,work)
      call vcheki(ier,numerr,'hydror  ')
!
!  compute cpfac
!
      call invmtrx(tau,kz,w1,kz,kz,det,iw2,ier,work)
      call vcheki(ier,numerr,'taur    ')
!
      do k = 1 , kz
        cpfac(k) = 0.
        do l = 1 , kz
          cpfac(k) = cpfac(k) + (sigmaf(l+1)-sigmaf(l))*w1(l,k)
        end do
      end do
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!       determine arrays needed for daley's variational scheme
!             for determination of surface pressure changes
!
      hweigh = 0.0
      hweigh(kz) = 1.    ! only lowest sigma level t considered
!
      do k1 = 1 , kz
        do k2 = 1 , kz    ! compute b(-1t) w/tbar**2 b(-1)
          w1(k2,k1) = 0.
          do k = 1 , kz
            w1(k2,k1) = hydror(k,k2)*hydror(k,k1)*hweigh(k)/            &
                    & (tbarh(k)**2)+w1(k2,k1)
          end do
        end do
      end do
!
      ps2 = ps*ps
      do k1 = 1 , kzp1
        do k2 = 1 , kz
          varpa1(k2,k1) = 0.
          do k = 1 , kz
            varpa1(k2,k1) = varpa1(k2,k1) + w1(k2,k)*hydroc(k,k1)*ps2
          end do
        end do
      end do
!
      do k1 = 1 , kzp1
        do k2 = 1 , kzp1
          varpa2(k2,k1) = 0.
          do k = 1 , kz
            varpa2(k2,k1) = varpa2(k2,k1) + hydroc(k,k2)*varpa1(k,k1)
          end do
        end do
      end do
!
      alpha1 = hydros(kz,kz)*tbarh(kz)/ps
      alpha2 = hweigh(kz)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!       output desired arrays
!
      call vprntv(sigmaf,kzp1,'sigmaf  ')
      call vprntv(tbarh,kz,'t mean  ')
      pps(1) = ps
      call vprntv(pps,1,'ps mean ')
      print 99003 , kz , numerr
99003 format ('0 vertical mode problem completed for kx=',i3,5x,i1,     &
             &' errors detected   (should be 0)')
!
!  printout if desired
      if ( .not.lprint ) then
        return
      end if
      call vprntv(cpfac,kz,'cpfac   ')
      call vprntv(dsigma,kz,'dsigma  ')
      call vprntv(hbar,kz,'hbar    ')
      call vprntv(sigmah,kzp1,'sigmah  ')
      call vprntv(tbarf,kzp1,'tbarf   ')
      call vprntv(thetah,kz,'thetah  ')
      call vprntv(thetaf,kzp1,'thetaf  ')
      call vprntv(hweigh,kz,'hweigh  ')
      print 99004 , alpha1 , alpha2
99004 format ('0alpha1 =',1p,1E16.5,'       alpha2 =',1p,1E16.5)
      call vprntm(a,kz,kz,'a       ')
      call vprntm(hydros,kz,kz,'hydros  ')
      call vprntm(hydror,kz,kz,'hydror  ')
      call vprntm(hydroc,kz,kzp1,'hydroc  ')
      call vprntm(tau,kz,kz,'tau     ')
      call vprntm(zmatx,kz,kz,'zmatx   ')
      call vprntm(zmatxr,kz,kz,'zmatxr  ')
      call vprntm(varpa1,kz,kzp1,'varpa1  ')
      call vprntm(varpa2,kzp1,kzp1,'varpa2  ')
!
      return
 
      end subroutine vmodes
!
!  Check that eigenvalues are real and positive valued.
!
      subroutine vcheke(er,ei,nk,numerr,aname)
      implicit none
!
! PARAMETER definitions
!
      real(8) , parameter :: tol = 1.E-9
!
! Dummy arguments
!
      character(8) :: aname
      integer :: nk , numerr
      real(8) , dimension(nk) :: ei , er
      intent (in) aname , ei , er , nk
      intent (inout) numerr
!
! Local variables
!
      real(8) :: emax
      integer :: n , nimag , numneg
!
      numneg = 0
      emax = 0.
      do n = 1 , nk
        if ( er(n).le.0. ) numneg = numneg + 1
        if ( er(n).gt.emax ) emax = er(n)
      end do
!
      nimag = 0
      do n = 1 , nk
        if ( ei(n)/emax.gt.tol ) nimag = nimag + 1
      end do
!
      if ( numneg+nimag.eq.0 ) then
        return
      end if
!
      numerr = numerr + 1
      print 99001 , aname , numneg , nimag
99001 format ('0 problem with equivalent depths determined from ',a8,/, &
            & 10x,i3,' depths are nonpositive valued',10x,i3,           &
             &' are not real')
!
      end subroutine vcheke
!
!  Flag of detected errors in linear algebra routines
!
      subroutine vcheki(ier,numerr,aname)
      implicit none
!
! Dummy arguments
!
      character(8) :: aname
      integer :: ier , numerr
      intent (in) aname , ier
      intent (inout) numerr
!
      if ( ier.ne.0 ) then
        numerr = numerr + 1
        print 99001 , aname , ier
99001   format ('0 error in determination of ',a8,                      &
               &' using library routine     ier=',i4)
      end if
!
      end subroutine vcheki
!
!  This routine normalizes the columns of z such that the component
!  with the largest absolute value is positive, and the sum of the
!  mass-weighted squares equals one.
!
      subroutine vnorml(z,s,nk,nk1)
      implicit none
!
! Dummy arguments
!
      integer :: nk , nk1
      real(8) , dimension(nk1) :: s
      real(8) , dimension(nk,nk) :: z
      intent (in) nk , nk1 , s
      intent (inout) z
!
! Local variables
!
      real(8) :: a , v , zmax
      integer :: k , kmax , l
!
      kmax = 1
      do l = 1 , nk
        zmax = -1.
        v = 0.
!
        do k = 1 , nk
          a = dabs(z(k,l))
          if ( a.gt.zmax ) then
            zmax = a
            kmax = k
          end if
          v = (s(k+1)-s(k))*a*a + v
        end do
!
        a = (z(kmax,l)/zmax)/dsqrt(v)
        do k = 1 , nk
          z(k,l) = a*z(k,l)
        end do
      end do
!
      end subroutine vnorml
!
!     Matrix inversion using linpack
!
      subroutine invmtrx(a,na,v,nv,n,d,ip,ier,work)
      use mod_message
      implicit none
      integer na , nv , n , ier , info , ip(n)
      real(kind=8) a(n,n) , v(n,n) , work(n) , d(2)
      integer i , j , job
!
!     08/23/91 version 1.0
!     12/10/92 updated to correct bugs
!     note: different from cray routine invmtx
!     uses subroutines sgefa/sgedi from library linpack
!     see dick valent, (SCD, consulting) if problems
!
      if ( n.ne.na .or. n.ne.nv ) call fatal(__FILE__,__LINE__,         &
          &'valent invmtx: equate n, na, nv')
!
      do j = 1 , n
        do i = 1 , n
          v(i,j) = a(i,j)
        end do
      end do
      call sgefa(v,n,n,ip,info)
      if ( info.ne.0 ) then
        write (aline,*) 'sgefa info = ' , info
        call say
        call fatal(__FILE__,__LINE__,'sgefa error')
      end if
      job = 11
      call sgedi(v,n,n,ip,d,work,job)
      ier = info
      end subroutine invmtrx
!
!  This routine orders the components of hbar so they are largest to
!  smallest valued.  The columns of z are reorded so that they
!  correspond to the same (but reordered) components of hbar.
!
      subroutine vorder(z,hbar,wz,wh,nk)
      implicit none
!
! Dummy arguments
!
      integer :: nk
      real(8) , dimension(nk) :: hbar
      real(8) , dimension(nk,2) :: wh
      real(8) , dimension(nk,nk) :: wz , z
      intent (in) nk
      intent (inout) hbar , wh , wz , z
!
! Local variables
!
      real(8) :: hmax
      integer :: k , kmax , l
!
      kmax = 1
      do k = 1 , nk
        wh(k,1) = hbar(k)
        wh(k,2) = 0.
        do l = 1 , nk
          wz(k,l) = z(k,l)
        end do
      end do
!
      do l = 1 , nk
        hmax = -1.D100
        do k = 1 , nk
          if ( (wh(k,2).eq.0.) .and. (wh(k,1).gt.hmax) ) then
            hmax = wh(k,1)
            kmax = k
          end if
        end do
!
        hbar(l) = hmax
        wh(kmax,2) = 1.
        do k = 1 , nk
          z(k,l) = wz(k,kmax)
        end do
      end do
!
      end subroutine vorder
!
!     Printout helpers
!
      subroutine vprntv(a,n,nam)
      implicit none
!
! Dummy arguments
!
      integer :: n
      character(8) :: nam
      real(8) , dimension(n) :: a
      intent (in) a , n , nam
 
      print 99001 , nam , a
99001 format ('0',a8,1x,1p,11G11.3,1x,/,9x,1p,11G11.3)
      end subroutine vprntv
!
!
! 
      subroutine vprntm(a,n1,n2,nam)
      implicit none
!
! Dummy arguments
!
      integer :: n1 , n2
      character(8) :: nam
      real(8) , dimension(n1,n2) :: a
      intent (in) a , n1 , n2 , nam
!
! Local variables
!
      integer :: k , l
 
      print 99001 , nam
      do k = 1 , n1
        print 99002 , k , (a(k,l),l=1,n2)
      end do
99001 format ('1',a8,/)
99002 format (1x,i3,5x,1p,11G11.3,1x,/,9x,1p,11G11.3)
      end subroutine vprntm
!
!  Check if tbar is stable.  This is not the actual stability condition
!  consistent with the model finite differences.
!
      subroutine vchekt(tbarh,sigmah,sigmaf,xkappa,pt,pd,nk,numerr)
      implicit none
!
! Dummy arguments
!
      integer :: nk , numerr
      real(8) :: pd , pt , xkappa
      real(8) , dimension(nk+1) :: sigmaf
      real(8) , dimension(nk) :: sigmah , tbarh
      intent (in) nk , pd , pt , sigmaf , sigmah , tbarh , xkappa
      intent (inout) numerr
!
! Local variables
!
      real(8) :: ds1 , ds2 , g1 , g2 , tb
      integer :: k
      logical :: lstab
!
      lstab = .true.
      do k = 1 , nk - 1
        ds1 = sigmaf(k+1) - sigmaf(k)
        ds2 = sigmaf(k+2) - sigmaf(k+1)
        tb = (ds1*tbarh(k)+ds2*tbarh(k+1))/(ds1+ds2)
        g1 = xkappa*tb/(sigmaf(k+1)+pt/pd)
        g2 = (tbarh(k+1)-tbarh(k))/(sigmah(k+1)-sigmah(k))
        if ( g1-g2.lt.0. ) lstab = .false.
      end do
      if ( .not.lstab ) then
        numerr = numerr + 1
        print 99001
99001   format ('0 indication that tbarh statically unstable')
      end if
!
      end subroutine vchekt
      end module mod_vmodes
