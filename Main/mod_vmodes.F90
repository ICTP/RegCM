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
  use mod_memutil
  use mod_mpmessage
  use mod_mppparam
  use mod_service
  use linpack
  use eispack

  private

  public :: allocate_mod_vmodes , vmodes

  public :: a0
  public :: sigmah , tbarh , hbar
  public :: xps , pd
  public :: zmatx , zmatxr
  public :: tau , varpa1
  public :: hydroc , hydros

  real(dp) :: xps , pd
  real(dp) , pointer , dimension(:,:) :: a0
  real(dp) , pointer , dimension(:) ::  sigmah , tbarh , hbar
  real(dp) , pointer , dimension(:,:) :: zmatx , zmatxr
  real(dp) , pointer , dimension(:,:) :: tau
  real(dp) , pointer , dimension(:,:) :: varpa1
  real(dp) , pointer , dimension(:,:) :: hydroc , hydros
!
  contains

  subroutine allocate_mod_vmodes
    implicit none
    call getmem2d(a0,1,kz,1,kz,'vmodes:a0')
    call getmem1d(hbar,1,kz,'vmodes:hbar')
    call getmem1d(sigmah,1,kzp1,'vmodes:sigmah')
    call getmem1d(tbarh,1,kz,'vmodes:tbarh')
    call getmem2d(zmatx,1,kz,1,kz,'vmodes:zmatx')
    call getmem2d(zmatxr,1,kz,1,kz,'vmodes:zmatxr')
    call getmem2d(tau,1,kz,1,kz,'vmodes:tau')
    call getmem2d(varpa1,1,kz,1,kzp1,'vmodes:varpa1')
    call getmem2d(hydroc,1,kz,1,kzp1,'vmodes:hydroc')
    call getmem2d(hydros,1,kz,1,kz,'vmodes:hydros')
  end subroutine allocate_mod_vmodes
!
!----------------------------------------------------------------------
!
!   This subroutine determines the vertical modes of the PSU/NCAR meso-
!   scale model designated MM4.  It also computes associated transform
!   matrices used by the initialization software used with MM4.
!   Adapted to be used in RegCM.
!
!----------------------------------------------------------------------
!
!
!   Programmed by Ronald M. Errico at NCAR,  Dec 1984.
!   Revised by Ronald Errico and Gary Bates, Nov 1987.
!   Revised by Ronald Errico,                Mar 1988.
!   For further info see: NCAR Tech Note by Errico and Bates, 1988.
!
!   lstand = .true. if standard atmosphere t to be used (ignore input
!             tbarh and xps in that case).  Otherwise, xps and tbarh must
!             be defined on input.
!
  subroutine vmodes(lstand)
!
    implicit none
!
    logical , intent(in) :: lstand
!
    real(dp) , dimension(2) :: det
    integer :: ier , k , k1 , k2 , l , mm , numerr
    logical :: lhydro , lprint , lsigma
    real(dp) :: ps2 , x
    real(dp) , dimension(kz) :: work
    real(dp) , dimension(1) :: pps
    real(dp) , dimension(kz,kz) :: a1 , a2 , a3 , a4 , d1 , d2 , &
                   e1 , e2 , e3 , g1 , g2 , g3 , s1 , s2 , w1 , w2 , x1
    real(dp) , dimension(kzp1,kz) :: w3
    integer , dimension(kz) :: iw2
    real(dp) , dimension(kzp1) :: tbarf , thetaf
    real(dp) , dimension(kz) :: thetah , tweigh
    real(dp) :: alpha1 , alpha2
    real(dp) , dimension(kz) :: cpfac , sdsigma , hweigh
    real(dp) , dimension(kzp1,kzp1) :: varpa2
    real(dp) , dimension(kz,kz) :: hydror
    data lprint/.false./  ! true if all matrices to be printed
!
    character (len=64) :: subroutine_name='vmodes'
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
    pd = xps - ptop
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
    if ( dabs(sigma(1)) > dlowval ) lsigma = .true.
    if ( dabs(sigma(kzp1)-d_one) > dlowval ) lsigma = .true.
    do k = 1 , kz
      if ( sigma(k+1) < sigma(k) ) then
        lsigma = .true.
        write (aline,99001) k , sigma(k+1) , sigma(k)
        call say
      end if
    end do
    if ( lsigma ) then
      call fatal(__FILE__,__LINE__,'Sigma values in vmode inappropriate')
    end if
!
!  compute sigmah (sigma at half levels) and delta sigma
    do k = 1 , kz
      sigmah(k) = (sigma(k)+sigma(k+1))*d_half
      sdsigma(k) = sigma(k+1) - sigma(k)
    end do
    sigmah(kzp1) = d_one
!
!  set tbarh (temperature at half (data) levels: indexed k + 1/2)
    if ( lstand ) call vtlaps
    call vchekt
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
!      thetah = tbarh*((sigmah+ptop/pd)**(-rovcp))

!
!  compute tbarf and thetaf
!
    do k = 2 , kz
      k1 = k - 1
      tbarf(k) = tbarh(k1)*(sigmah(k)-sigma(k))/(sigmah(k)-sigmah(k1)) + &
                 tbarh(k)*(sigma(k)-sigmah(k1))/(sigmah(k)-sigmah(k1))
    end do
    tbarf(1) = d_zero
    tbarf(kzp1) = d_zero
!
    do k = 1 , kzp1
      if ( sigma(k) < dlowval ) then
        thetaf(k) = tbarf(k)
      else
        thetaf(k) = tbarf(k)*((sigma(k)+ptop/pd)**(-rovcp))
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
      d1(k,k) = sigma(k+1) - sigma(k)
      d2(k,k) = rovcp*tbarh(k)/(sigmah(k)+ptop/pd)
      s1(k,k) = sigma(k)
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
      w1(k,1) = dlog((sigmah(k)+ptop/pd)/(sigmah(k-1)+ptop/pd))
    end do
!
!  compute matrix which multiples t vector
!
    hydros = d_zero
!
    do k = 1 , kz - 1
      do l = k , kz - 1
        hydros(k,l) = hydros(k,l) + &
                    w1(l+1,1)*sdsigma(l)/(sdsigma(l+1)+sdsigma(l))
        hydros(k,l+1) = hydros(k,l+1) + &
                    w1(l+1,1)*sdsigma(l+1)/(sdsigma(l+1)+sdsigma(l))
      end do
    end do
!
    do k = 1 , kz
      hydros(k,kz) = hydros(k,kz) + dlog((d_one+ptop/pd)/(sigmah(kz)+ptop/pd))
    end do
!
!  compute matirx which multiplies log(sigma*p+ptop) vector
!
    hydroc = d_zero
!
    tweigh(1) = d_zero
    do l = 2 , kz
      tweigh(l) = (tbarh(l)*sdsigma(l) + &
                   tbarh(l-1)*sdsigma(l-1))/(sdsigma(l)+sdsigma(l-1))
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
      w1(k,2) = -tbarh(k)*dlog(sigmah(k)*pd+ptop)
      do l = 1 , kzp1
        w1(k,2) = w1(k,2) + hydroc(k,l)*dlog(sigmah(l)*pd+ptop)
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
        w3(k,l) = sdsigma(l)/(d_one+ptop/(pd*sigmah(k)))
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
!   determine other matrices and vectors
!
!   compute eigenvalues and vectors for tau (rg calls eispack routines)
!
    w1 = tau
    call rg(kz,w1,hbar,w2,1,zmatx,ier)
    call vcheki(ier,numerr,'zmatx   ')
    call vcheke
    call vorder
    call vnorml
!
!   compute inverse of zmatx
!
    call invmtrx(zmatx,kz,zmatxr,kz,kz,det,iw2,ier,work)
    call vcheki(ier,numerr,'zmatxr  ')
!
!   compute inverse of hydros
!
    call invmtrx(hydros,kz,hydror,kz,kz,det,iw2,ier,work)
    call vcheki(ier,numerr,'hydror  ')
!
!   compute cpfac
!
    call invmtrx(tau,kz,w1,kz,kz,det,iw2,ier,work)
    call vcheki(ier,numerr,'taur    ')
!
    do k = 1 , kz
      cpfac(k) = d_zero
      do l = 1 , kz
        cpfac(k) = cpfac(k) + (sigma(l+1)-sigma(l))*w1(l,k)
      end do
    end do
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!   determine arrays needed for daley's variational scheme
!   for determination of surface pressure changes
!
    hweigh = d_zero
    hweigh(kz) = d_one  ! only lowest sigma level t considered
!
    do k1 = 1 , kz
      do k2 = 1 , kz    ! compute b(-1t) w/tbar**2 b(-1)
        w1(k2,k1) = d_zero
        do k = 1 , kz
          w1(k2,k1) = hydror(k,k2)*hydror(k,k1)*hweigh(k) / &
                    (tbarh(k)**d_two)+w1(k2,k1)
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
!   subract a4 from a for use in computing am.
!
    do k1 = 1 , kz
      do k2 = 1 , kz
        a0(k2,k1) = a0(k2,k1) - a4(k2,k1)
      end do
    end do
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!   output desired arrays
!
    if ( myid == iocpu ) then
      call vprntv(sigma,kzp1,'sigma   ')
      call vprntv(tbarh,kz,'t mean  ')
      pps(1) = xps
      call vprntv(pps,1,'ps mean ')
      print 99003 , kz , numerr
    end if
!
!   printout if desired
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
    call vprntm(hsigma,kz,kz,'a       ')
    call vprntm(hydros,kz,kz,'hydros  ')
    call vprntm(hydror,kz,kz,'hydror  ')
    call vprntm(hydroc,kz,kzp1,'hydroc  ')
    call vprntm(tau,kz,kz,'tau     ')
    call vprntm(zmatx,kz,kz,'zmatx   ')
    call vprntm(zmatxr,kz,kz,'zmatxr  ')
    call vprntm(varpa1,kz,kzp1,'varpa1  ')
    call vprntm(varpa2,kzp1,kzp1,'varpa2  ')
!

99001 format ('0 for k=',i3,' sigma(k+1)=',f9.6,' <= sigma(k)=',f9.6)
99002 format ('0 problem with linearization of hydostatic equation')
99003 format ('0 vertical mode problem completed for kx=',i3,5x,i1, &
              ' errors detected   (should be 0)')
99004 format ('0alpha1 =',1p,1E16.5,'       alpha2 =',1p,1E16.5)
 
    call time_end(subroutine_name,idindx)

    contains
!
!   Check if tbar is stable.  This is not the actual stability condition
!   consistent with the model finite differences.
!
    subroutine vchekt
      implicit none
!
      real(dp) :: ds1 , ds2 , g1 , g2 , tb
      integer :: k
      logical :: lstab
!
      character (len=64) :: subroutine_name='vchekt'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
      lstab = .true.
      do k = 1 , kz - 1
        ds1 = sigma(k+1) - sigma(k)
        ds2 = sigma(k+2) - sigma(k+1)
        tb = (ds1*tbarh(k)+ds2*tbarh(k+1))/(ds1+ds2)
        g1 = rovcp*tb/(sigma(k+1)+ptop/pd)
        g2 = (tbarh(k+1)-tbarh(k))/(sigmah(k+1)-sigmah(k))
        if ( g1-g2 < d_zero ) lstab = .false.
      end do
      if ( .not.lstab ) then
        numerr = numerr + 1
        print 99001
      end if
!
      call time_end(subroutine_name,idindx)
!
99001 format ('0 indication that tbarh statically unstable')
!
    end subroutine vchekt
!
!   This routine computes the temperature corresponding to a U.S.
!   standard atmosphere (see text by hess). Units of p are cb.
!
    subroutine vtlaps
      implicit none
!
      real(dp) , parameter :: tstrat = 218.15D0
      real(dp) , parameter :: zstrat = 10769.0D0
!
      real(dp) :: p0 , fac , p , z
      integer :: k
!
      p0 = stdp*d_r1000
      fac = rgas*lrate*regrav
      do k = 1 , kz
        p = sigmah(k)*pd + ptop
        tbarh(k) = stdt*((p/p0)**fac)
        z = (stdt-tbarh(k))/lrate
        if ( z > zstrat ) tbarh(k) = tstrat
      end do
!
    end subroutine vtlaps
!
!   This routine orders the components of hbar so they are largest to
!   smallest valued.  The columns of zmatx are reorded so that they
!   correspond to the same (but reordered) components of hbar.
!
    subroutine vorder
      implicit none
!
      real(dp) :: hmax
      integer :: k , kmax , l
!
      kmax = 1
      do k = 1 , kz
        w2(k,1) = hbar(k)
        w2(k,2) = d_zero
        do l = 1 , kz
          w1(k,l) = zmatx(k,l)
        end do
      end do
!
      do l = 1 , kz
        hmax = -1.0D100
        do k = 1 , kz
          if ( (dabs(w2(k,2)) < dlowval) .and. (w2(k,1) > hmax) ) then
            hmax = w2(k,1)
            kmax = k
          end if
        end do
!
        hbar(l) = hmax
        w2(kmax,2) = d_one
        do k = 1 , kz
          zmatx(k,l) = w1(k,kmax)
        end do
      end do
!
    end subroutine vorder
!
!   This routine normalizes the columns of z such that the component
!   with the largest absolute value is positive, and the sum of the
!   mass-weighted squares equals one.
!
    subroutine vnorml
      implicit none
!
      real(dp) :: a , v , zmax
      integer :: k , kmax , l
!
      kmax = 1
      do l = 1 , kz
        zmax = -d_one
        v = d_zero
!
        do k = 1 , kz
          a = dabs(zmatx(k,l))
          if ( a > zmax ) then
            zmax = a
            kmax = k
          end if
          v = (sigma(k+1)-sigma(k))*a*a + v
        end do
!
        a = (zmatx(kmax,l)/zmax)/dsqrt(v)
        do k = 1 , kz
          zmatx(k,l) = a*zmatx(k,l)
        end do
      end do
!
    end subroutine vnorml
!
!   Check that eigenvalues are real and positive valued.
!
    subroutine vcheke
      implicit none
!
      real(dp) , parameter :: tol = 1.0D-9
!
      real(dp) :: emax
      integer :: n , nimag , numneg
!
      numneg = 0
      emax = d_zero
      do n = 1 , kz
        if ( hbar(n) <= d_zero ) numneg = numneg + 1
        if ( hbar(n) > emax ) emax = hbar(n)
      end do
!
      nimag = 0
      do n = 1 , kz
        if ( w2(n,1)/emax > tol ) nimag = nimag + 1
      end do
!
      if ( numneg+nimag == 0 ) then
        return
      end if
!
      numerr = numerr + 1
      print 99001 , numneg , nimag
!
99001 format ('0 problem with equivalent depths determined from tau',/, &
         10x,i3,' depths are nonpositive valued',10x,i3,' are not real')
!
    end subroutine vcheke
!
  end subroutine vmodes
!
! Flag of detected errors in linear algebra routines
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
    end if
!
99001 format ('0 error in determination of ',a8, &
              ' using library routine ier=',i4)
  end subroutine vcheki
!
! Matrix inversion using linpack
!
  subroutine invmtrx(a,na,v,nv,n,d,ip,ier,work)
    implicit none
!
    integer :: na , nv , n , ier
    integer , dimension(n) :: ip
    real(dp) :: a(n,n) , v(n,n) , work(n) , d(2)
    integer :: i , j
!
!   08/23/91 Version 1.0
!   12/10/92 Updated to correct bugs
!   Note: different from cray routine invmtx
!   Uses subroutines sgefa/sgedi from library linpack
!
    character (len=64) :: subroutine_name='invmtrx'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
    if ( n /= na .or. n /= nv ) then
      call fatal(__FILE__,__LINE__,'invmtx: equate n, na, nv')
    end if
!
    do j = 1 , n
      do i = 1 , n
        v(i,j) = a(i,j)
      end do
    end do
    call sgefa(v,n,n,ip,ier)
    if ( ier /= 0 ) then
      write (aline,*) 'sgefa error info = ' , ier
      call say
      call fatal(__FILE__,__LINE__,'sgefa error')
    end if
    call sgedi(v,n,n,ip,d,work,11)
    call time_end(subroutine_name,idindx)  
  end subroutine invmtrx
!
! Printout helpers
!
  subroutine vprntv(a,n,nam)
    implicit none
!
    integer :: n
    character(8) :: nam
    real(dp) , dimension(n) :: a
    intent (in) a , n , nam
    integer :: k
    print * , nam
    do k = 1 , n
      print *, a(k)
    end do
  end subroutine vprntv
!
  subroutine vprntm(a,n1,n2,nam)
    implicit none
!
    integer :: n1 , n2
    character(8) :: nam
    real(dp) , dimension(n1,n2) :: a
    intent (in) a , n1 , n2 , nam
    integer :: k , l
    print * , nam
    do k = 1 , n1
      print * , k , (a(k,l),l=1,n2)
    end do
  end subroutine vprntm
!
end module mod_vmodes
