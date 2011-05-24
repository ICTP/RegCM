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

      use mod_runparams
      use mod_message
      use mod_service
      use linpack
      use eispack
      private

      public :: allocate_mod_vmodes , vmodes

      public :: a0 , a4
      public :: hbar , sigmah , tbarh
      public :: xps , pd
      public :: zmatx , zmatxr
      public :: tau , varpa1
      public :: hydroc , hydros

      real(8) :: xps , pd
      real(8) , allocatable , dimension(:,:) :: a0 , a4
      real(8) , allocatable , dimension(:) ::  hbar , sigmah , tbarh
      real(8) , allocatable , dimension(:,:) :: zmatx , zmatxr
      real(8) , allocatable , dimension(:,:) :: tau
      real(8) , allocatable , dimension(:,:) :: varpa1
      real(8) , allocatable , dimension(:,:) :: hydroc , hydros
!

      contains

      subroutine allocate_mod_vmodes
        implicit none
        allocate(a0(kz,kz))
        allocate(a4(kz,kz))
        allocate(hbar(kz))
        allocate(sigmah(kzp1))
        allocate(tbarh(kz))
        allocate(zmatx(kz,kz))
        allocate(zmatxr(kz,kz))
        allocate(tau(kz,kz))
        allocate(varpa1(kz,kzp1))
        allocate(hydroc(kz,kzp1))
        allocate(hydros(kz,kz))
        a0 = d_zero
        a4 = d_zero
        hbar = d_zero
        sigmah = d_zero
        tbarh = d_zero
        zmatx = d_zero
        zmatxr = d_zero
        tau = d_zero
        varpa1 = d_zero
        hydroc = d_zero
        hydros = d_zero

      end subroutine allocate_mod_vmodes

      subroutine vmodes(lstand,sigmaf,kv1)
!
      implicit none
!
      integer :: kv1
      logical :: lstand
      real(8) , dimension(kv1) :: sigmaf
      intent (in) kv1 , lstand
!
      real(8) , dimension(2) :: det
      integer :: ier , k , k1 , k2 , l , mm , numerr
      logical :: lhydro , lprint , lsigma
      real(8) :: ps2 , x
      real(8) , dimension(kz) :: work
      real(8) , dimension(1) :: pps
      real(8) , dimension(kz,kz) :: a1 , a2 , a3 , &
                & d1 , d2 , e1 , e2 , e3 , g1 , g2 , g3 , s1 , &
                & s2 , w1 , w2 , x1
      real(8) , dimension(kzp1,kz) :: w3
      integer , dimension(kz) :: iw2
      real(8) , dimension(kzp1) :: tbarf , thetaf
      real(8) , dimension(kz) :: thetah , tweigh
      real(8) :: alpha1 , alpha2
      real(8) , dimension(kz) :: cpfac , sdsigma , hweigh
      real(8) , dimension(kzp1,kzp1) :: varpa2
      real(8) , dimension(kz,kz) :: hydror
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
!            tbarh and xps in that case).  otherwise, xps and tbarh must
!            be defined on input.  note that in either case, r8pt must
!            also be defined on input (common block named cvert).
!
      data lprint/.false./  ! true if all matrices to be printed
!
      character (len=50) :: subroutine_name='vmodes'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
      a1 = d_zero
      a2 = d_zero
      a3 = d_zero
      d1 = d_zero
      d2 = d_zero
      e1 = d_zero
      e2 = d_zero
      e3 = d_zero
      g1 = d_zero
      g2 = d_zero
      g3 = d_zero
      s1 = d_zero
      s2 = d_zero
      w1 = d_zero
      w2 = d_zero
      x1 = d_zero
      iw2 = 0
      thetah = d_zero
      tweigh = d_zero
      tbarf = d_zero
      thetaf = d_zero
      w3 = d_zero
      cpfac = d_zero
      sdsigma = d_zero
      hweigh = d_zero
      hydror = d_zero
      varpa2 = d_zero
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
      if ( lstand ) xps = d_100
                            ! standard xps in cb; otherwise xps set in tav
      pd = xps - r8pt
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
      if ( dabs(sigmaf(1)) > dlowval ) lsigma = .true.
      if ( dabs(sigmaf(kzp1)-d_one) > dlowval ) lsigma = .true.
      do k = 1 , kz
        if ( sigmaf(k+1) < sigmaf(k) ) then
          lsigma = .true.
          write (aline,99001) k , sigmaf(k+1) , sigmaf(k)
          call say
        end if
      end do
      if ( lsigma )                                                     &
      & call fatal(__FILE__,__LINE__,                                   &
      &          'Sigma values in list vmode inappropriate')
!
!  compute sigmah (sigma at half levels) and delta sigma
      do k = 1 , kz
        sigmah(k) = (sigmaf(k)+sigmaf(k+1))*d_half
        sdsigma(k) = sigmaf(k+1) - sigmaf(k)
      end do
      sigmah(kzp1) = d_one
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
!  this array is never used: it is only computed and printed out 
!  the  following line causes a segfault on IBM SP6 when compiled with -q check 
!  because it uses F90 syntax arrays with different sizes  
!  S.C. 21/05/2010
!  I therefore decided to comment out this line 
!      thetah = tbarh*((sigmah+r8pt/pd)**(-rovcp))

!
!  compute tbarf and thetaf
!
      do k = 2 , kz
        k1 = k - 1
        tbarf(k) = tbarh(k1)*(sigmah(k)-sigmaf(k))                      &
                 & /(sigmah(k)-sigmah(k1)) + tbarh(k)                   &
                 & *(sigmaf(k)-sigmah(k1))/(sigmah(k)-sigmah(k1))
      end do
      tbarf(1) = d_zero
      tbarf(kzp1) = d_zero
!
      do k = 1 , kzp1
        if ( sigmaf(k) < dlowval ) then
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
          if ( l > k ) e2(k,l) = d_zero
          if ( l <= k ) e2(k,l) = d_one
          e1(k,l) = d_one
        end do
      end do
!
      do k = 1 , kz
        a3(k,k) = -tbarh(k)
        d1(k,k) = sigmaf(k+1) - sigmaf(k)
        d2(k,k) = rovcp*tbarh(k)/(sigmah(k)+r8pt/pd)
        s1(k,k) = sigmaf(k)
        s2(k,k) = sigmah(k)
        x1(k,k) = d_one
      end do
!
      do k = 1 , kz
        do l = 1 , kz
          e3(k,l) = d_zero
          g1(k,l) = d_zero
        end do
        e3(k,k) = d_one
        if ( k > 1 ) g1(k,k) = tbarf(k)
        if ( k < kz ) g1(k,k+1) = -tbarf(k+1)
        if ( k < kz ) e3(k,k+1) = d_one
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
      w2 = 0.0D0
      do k = 1 , kz
        w2(k,k) = d_one/d1(k,k)
      end do
      w1 = matmul(g1,g2)
      a1 = matmul(w2,w1)
!
!  compute a2
!
      w1 = matmul(e1,d1)
      a2 = matmul(s2,w1)
      w2 = matmul(e3,g2)
      w2 = w2*d_half
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
      a0 = a1+a2+a3+a4
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
      hydros = d_zero
!
      do k = 1 , kz - 1
        do l = k , kz - 1
          hydros(k,l) = hydros(k,l) + w1(l+1,1)*sdsigma(l)              &
                      & /(sdsigma(l+1)+sdsigma(l))
          hydros(k,l+1) = hydros(k,l+1) + w1(l+1,1)*sdsigma(l+1)        &
                        & /(sdsigma(l+1)+sdsigma(l))
        end do
      end do
!
      do k = 1 , kz
        hydros(k,kz) = hydros(k,kz)                                     &
                     & + dlog((d_one+r8pt/pd)/(sigmah(kz)+r8pt/pd))
      end do
!
!  compute matirx which multiplies log(sigma*p+r8pt) vector
!
      hydroc = d_zero
!
      tweigh(1) = d_zero
      do l = 2 , kz
        tweigh(l) = (tbarh(l)*sdsigma(l)+tbarh(l-1)*sdsigma(l-1))       &
                  & /(sdsigma(l)+sdsigma(l-1))
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
        w1(k,1) = d_zero
        do l = 1 , kz
          w1(k,1) = w1(k,1) + hydros(k,l)*tbarh(l)
        end do
        w1(k,2) = -tbarh(k)*dlog(sigmah(k)*pd+r8pt)
        do l = 1 , kzp1
          w1(k,2) = w1(k,2) + hydroc(k,l)*dlog(sigmah(l)*pd+r8pt)
        end do
        x = dabs(w1(k,1)-w1(k,2))/(dabs(w1(k,1))+dabs(w1(k,2)))
        if ( x > 1.0D-8 ) lhydro = .true.
      end do
!
      if ( lhydro ) then
        numerr = numerr + 1
        print 99002
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
          w3(k,l) = sdsigma(l)/(d_one+r8pt/(pd*sigmah(k)))
        end do
      end do
!
      do l = 1 , kz
        do k = 1 , kz
          w2(k,l) = d_zero
          do mm = 1 , kzp1
            w2(k,l) = w2(k,l) + hydroc(k,mm)*w3(mm,l)
          end do
        end do
      end do
!
      w1 = matmul(hydros,a0)
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
        cpfac(k) = d_zero
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
      hweigh = d_zero
      hweigh(kz) = d_one    ! only lowest sigma level t considered
!
      do k1 = 1 , kz
        do k2 = 1 , kz    ! compute b(-1t) w/tbar**2 b(-1)
          w1(k2,k1) = d_zero
          do k = 1 , kz
            w1(k2,k1) = hydror(k,k2)*hydror(k,k1)*hweigh(k)/            &
                    & (tbarh(k)**d_two)+w1(k2,k1)
          end do
        end do
      end do
!
      ps2 = xps*xps
      do k1 = 1 , kzp1
        do k2 = 1 , kz
          varpa1(k2,k1) = d_zero
          do k = 1 , kz
            varpa1(k2,k1) = varpa1(k2,k1) + w1(k2,k)*hydroc(k,k1)*ps2
          end do
        end do
      end do
!
      do k1 = 1 , kzp1
        do k2 = 1 , kzp1
          varpa2(k2,k1) = d_zero
          do k = 1 , kz
            varpa2(k2,k1) = varpa2(k2,k1) + hydroc(k,k2)*varpa1(k,k1)
          end do
        end do
      end do
!
      alpha1 = hydros(kz,kz)*tbarh(kz)/xps
      alpha2 = hweigh(kz)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!       output desired arrays
!
      if ( myid == 0 ) then
        call vprntv(sigmaf,kzp1,'sigmaf  ')
        call vprntv(tbarh,kz,'t mean  ')
        pps(1) = xps
        call vprntv(pps,1,'ps mean ')
        print 99003 , kz , numerr
      end if
!
!  printout if desired
      if ( .not.lprint ) then
        return
      end if
      call vprntv(cpfac,kz,'cpfac   ')
      call vprntv(sdsigma,kz,'sdsigma  ')
      call vprntv(hbar,kz,'hbar    ')
      call vprntv(sigmah,kzp1,'sigmah  ')
      call vprntv(tbarf,kzp1,'tbarf   ')
      call vprntv(thetah,kz,'thetah  ')
      call vprntv(thetaf,kzp1,'thetaf  ')
      call vprntv(hweigh,kz,'hweigh  ')
      print 99004 , alpha1 , alpha2
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

99001 format ('0 for k=',i3,' sigmaf(k+1)=',f9.6,' <= sigmaf(k)=',    &
             & f9.6)
99002 format ('0 problem with linearization of hydostatic equation')
99003 format ('0 vertical mode problem completed for kx=',i3,5x,i1,     &
             &' errors detected   (should be 0)')
99004 format ('0alpha1 =',1p,1E16.5,'       alpha2 =',1p,1E16.5)
 
      call time_end(subroutine_name,idindx)
      end subroutine vmodes
!
!  Check that eigenvalues are real and positive valued.
!
      subroutine vcheke(er,ei,nk,numerr,aname)
      implicit none
!
! PARAMETER definitions
!
      real(8) , parameter :: tol = 1.0D-9
!
      character(8) :: aname
      integer :: nk , numerr
      real(8) , dimension(nk) :: ei , er
      intent (in) aname , ei , er , nk
      intent (inout) numerr
!
      real(8) :: emax
      integer :: n , nimag , numneg
!
      numneg = 0
      emax = d_zero
      do n = 1 , nk
        if ( er(n) <= d_zero ) numneg = numneg + 1
        if ( er(n) > emax ) emax = er(n)
      end do
!
      nimag = 0
      do n = 1 , nk
        if ( ei(n)/emax > tol ) nimag = nimag + 1
      end do
!
      if ( numneg+nimag == 0 ) then
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
      character(8) :: aname
      integer :: ier , numerr
      intent (in) aname , ier
      intent (inout) numerr
!
      if ( ier /= 0 ) then
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
      integer :: nk , nk1
      real(8) , dimension(nk1) :: s
      real(8) , dimension(nk,nk) :: z
      intent (in) nk , nk1 , s
      intent (inout) z
!
      real(8) :: a , v , zmax
      integer :: k , kmax , l
!
      kmax = 1
      do l = 1 , nk
        zmax = -d_one
        v = d_zero
!
        do k = 1 , nk
          a = dabs(z(k,l))
          if ( a > zmax ) then
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
      implicit none
      integer :: na , nv , n , ier , info
      integer , dimension(n) :: ip
      real(8) :: a(n,n) , v(n,n) , work(n) , d(2)
      integer :: i , j , job
!
!     08/23/91 version 1.0
!     12/10/92 updated to correct bugs
!     note: different from cray routine invmtx
!     uses subroutines sgefa/sgedi from library linpack
!     see dick valent, (SCD, consulting) if problems
!
      character (len=50) :: subroutine_name='invmtrx'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
      if ( n /= na .or. n /= nv ) call fatal(__FILE__,__LINE__,         &
          &'valent invmtx: equate n, na, nv')
!
      do j = 1 , n
        do i = 1 , n
          v(i,j) = a(i,j)
        end do
      end do
      call sgefa(v,n,n,ip,info)
      if ( info /= 0 ) then
        write (aline,*) 'sgefa info = ' , info
        call say
        call fatal(__FILE__,__LINE__,'sgefa error')
      end if
      job = 11
      call sgedi(v,n,n,ip,d,work,job)
      ier = info
      call time_end(subroutine_name,idindx)  
      end subroutine invmtrx
!
!  This routine orders the components of hbar so they are largest to
!  smallest valued.  The columns of z are reorded so that they
!  correspond to the same (but reordered) components of hbar.
!
      subroutine vorder(z,hbar,wz,wh,nk)
      implicit none
!
      integer :: nk
      real(8) , dimension(nk) :: hbar
      real(8) , dimension(nk,2) :: wh
      real(8) , dimension(nk,nk) :: wz , z
      intent (in) nk
      intent (inout) hbar , wh , wz , z
!
      real(8) :: hmax
      integer :: k , kmax , l
!
      kmax = 1
      do k = 1 , nk
        wh(k,1) = hbar(k)
        wh(k,2) = d_zero
        do l = 1 , nk
          wz(k,l) = z(k,l)
        end do
      end do
!
      do l = 1 , nk
        hmax = -1.0D100
        do k = 1 , nk
          if ( (dabs(wh(k,2)) < dlowval) .and. (wh(k,1) > hmax) ) then
            hmax = wh(k,1)
            kmax = k
          end if
        end do
!
        hbar(l) = hmax
        wh(kmax,2) = d_one
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
      integer :: n1 , n2
      character(8) :: nam
      real(8) , dimension(n1,n2) :: a
      intent (in) a , n1 , n2 , nam
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
      integer :: nk , numerr
      real(8) :: pd , pt , xkappa
      real(8) , dimension(nk+1) :: sigmaf
      real(8) , dimension(nk) :: sigmah , tbarh
      intent (in) nk , pd , pt , sigmaf , sigmah , tbarh , xkappa
      intent (inout) numerr
!
      real(8) :: ds1 , ds2 , g1 , g2 , tb
      integer :: k
      logical :: lstab
!
      character (len=50) :: subroutine_name='vchekt'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
      lstab = .true.
      do k = 1 , nk - 1
        ds1 = sigmaf(k+1) - sigmaf(k)
        ds2 = sigmaf(k+2) - sigmaf(k+1)
        tb = (ds1*tbarh(k)+ds2*tbarh(k+1))/(ds1+ds2)
        g1 = xkappa*tb/(sigmaf(k+1)+pt/pd)
        g2 = (tbarh(k+1)-tbarh(k))/(sigmah(k+1)-sigmah(k))
        if ( g1-g2 < d_zero ) lstab = .false.
      end do
      if ( .not.lstab ) then
        numerr = numerr + 1
        print 99001
99001   format ('0 indication that tbarh statically unstable')
      end if
!
      call time_end(subroutine_name,idindx)
      end subroutine vchekt

!
      subroutine vtlaps(t,sigma,pt,pd,nk)

      implicit none
!
      real(8) , parameter :: tstrat = 218.15D0 , zstrat = 10769.0D0
      real(8) :: p0
!
      integer :: nk
      real(8) :: pd , pt
      real(8) , dimension(nk) :: sigma , t
      intent (in) nk , pd , pt , sigma
      intent (inout) t
!
      real(8) :: fac , p , z
      integer :: k
!
!  this routine computes the temperature corresponding to a u. s.
!  standard atmosphere (see text by hess). units of p are cb.
!
      p0 = stdp*d_r1000
      fac = rgas*lrate*regrav
      do k = 1 , nk
        p = sigma(k)*pd + pt
        t(k) = stdt*((p/p0)**fac)
        z = (stdt-t(k))/lrate
        if ( z > zstrat ) t(k) = tstrat
      end do
!
      end subroutine vtlaps

!
      end module mod_vmodes
