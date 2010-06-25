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
 
      module mod_outprt

      implicit none

      private

      public :: outprt

      contains

#ifdef MPP1
      subroutine outprt(iexec)

      implicit none
!
! Dummy arguments
!
      integer :: iexec
      intent (in) iexec
!
! Local variables
!
      integer dum
!
      dum = iexec
      return
      end subroutine outprt

#else

      subroutine outprt(iexec)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine creates the printer output for a quick check.   c
!     the complete analysis is handled by other program "dataflow"    c
!     later.                                                          c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_dynparam
      use mod_param1
      use mod_param2
      use mod_param3 , only : jxsex , kxout , r8pt , a
      use mod_main
      use mod_bdycod
      use mod_pmoist
      use mod_date
      use mod_constants , only : tzero , ep2 , svp1 , svp2 , svp3 ,     &
                      &          svp4 , svp5 , svp6 , gti
      implicit none
!
! Dummy arguments
!
      integer :: iexec
      intent (in) iexec
!
! Local variables
!
      character(40) :: c40nam
      real(8) :: const , pres , psab1 , psab2 , psabar , qvs , satvp ,  &
               & xth
      real(8) , dimension(jx,iy) :: fr , hscr1r , hscr2r , msfdr ,      &
                                  & psar , raincr , rainncr , tgar ,    &
                                  & xlatr , xlongr
      real(8) , dimension(iy,jx) :: hscr1 , hscr2
      integer :: i , iyn , j , jcx , jxn , k , kout , kxn
      real(8) , dimension(kz,iy) :: scrx1 , scrx2
      real(8) , dimension(iy,kz) :: scrx1r , scrx2r
!
      data iyn , jxn , kxn/2 , 2 , 1/
!
#ifndef BAND
!----------------------------------------------------------------------
!
      xth = xtime/60.
      jcx = jxsex
      kout = kxout
      print 99001 , ldatez + xtime/1440. , ktau , jyear
!
!-----north-south cross sections:
!
!.....u and v:
!
      if ( iexec.eq.1 ) then
!
!.....first time :
!
        if ( jcx.eq.1 ) then
!......for jcx = 1:
          do i = 2 , iym1
            do k = 1 , kz
              psabar = 0.5*(psa(i,jcx)+psa(i-1,jcx))
              scrx1(k,i) = ua(i,k,jcx)/psabar
              scrx2(k,i) = va(i,k,jcx)/psabar
            end do
          end do
          do k = 1 , kz
            scrx1(k,1) = ua(1,k,jcx)/psa(1,jcx)
            scrx1(k,iy) = ua(iy,k,jcx)/psa(iym1,jcx)
            scrx2(k,1) = va(1,k,jcx)/psa(1,jcx)
            scrx2(k,iy) = va(iy,k,jcx)/psa(iym1,jcx)
          end do
        else if ( jcx.eq.jx ) then
!......for jcx = jx:
          do i = 2 , iym1
            do k = 1 , kz
              psabar = 0.5*(psa(i,jxm1)+psa(i-1,jxm1))
              scrx1(k,i) = ua(i,k,jcx)/psabar
              scrx2(k,i) = va(i,k,jcx)/psabar
            end do
          end do
          do k = 1 , kz
            scrx1(k,1) = ua(1,k,jcx)/psa(1,jxm1)
            scrx1(k,iy) = ua(iy,k,jcx)/psa(iym1,jxm1)
            scrx2(k,1) = va(1,k,jcx)/psa(1,jxm1)
            scrx2(k,iy) = va(iy,k,jcx)/psa(iym1,jxm1)
          end do
        else
!......interior slice:
          do i = 2 , iym1
            do k = 1 , kz
              psabar = 0.25*(psa(i,jcx)+psa(i,jcx-1)+psa(i-1,jcx)       &
                     & +psa(i-1,jcx-1))
              scrx1(k,i) = ua(i,k,jcx)/psabar
              scrx2(k,i) = va(i,k,jcx)/psabar
            end do
          end do
          do k = 1 , kz
            psab1 = 0.5*(psa(1,jcx-1)+psa(1,jcx))
            psab2 = 0.5*(psa(iym1,jcx-1)+psa(iym1,jcx))
            scrx1(k,1) = ua(1,k,jcx)/psab1
            scrx1(k,iy) = ua(iy,k,jcx)/psab2
            scrx2(k,1) = va(1,k,jcx)/psab1
            scrx2(k,iy) = va(iy,k,jcx)/psab2
          end do
        end if
!
      else if ( iexec.gt.1 ) then
!
!.....subsequent calls:
!
        if ( jcx.eq.1 ) then
!......for jcx = 1:
          do i = 1 , iy
            do k = 1 , kz
              scrx1(k,i) = uj1(i,k)
              scrx2(k,i) = vj1(i,k)
            end do
          end do
        else if ( jcx.eq.jx ) then
!......for jcx = jx:
          do i = 1 , iy
            do k = 1 , kz
              scrx1(k,i) = ujl(i,k)
              scrx2(k,i) = vjl(i,k)
            end do
          end do
        else
!......interior slice:
          do i = 2 , iym1
            do k = 1 , kz
              psabar = 0.25*(psa(i,jcx)+psa(i,jcx-1)+psa(i-1,jcx)       &
                     & +psa(i-1,jcx-1))
              scrx1(k,i) = ua(i,k,jcx)/psabar
              scrx2(k,i) = va(i,k,jcx)/psabar
            end do
          end do
          do k = 1 , kz
            scrx1(k,1) = ui1(k,jcx)
            scrx1(k,iy) = uil(k,jcx)
            scrx2(k,1) = vi1(k,jcx)
            scrx2(k,iy) = vil(k,jcx)
          end do
        end if
!
      else        !end iexec test
      end if
!
      write (c40nam,99002) jcx
      call mapsmp(scrx1,scrx1r,kz,iy,1,kz,kxn,1,iy,iyn,0.D0,-1,c40nam,  &
                & xth)
      write (c40nam,99003) jcx
      call mapsmp(scrx2,scrx2r,kz,iy,1,kz,kxn,1,iy,iyn,0.D0,-1,c40nam,  &
                & xth)
!
!.....t:
!
      do i = 1 , iym1
        do k = 1 , kz
          scrx1(k,i) = ta(i,k,jcx)/psa(i,jcx)
        end do
      end do
      write (c40nam,99004) jcx
      call mapsmp(scrx1,scrx1r,kz,iy,1,kz,kxn,1,iym1,iyn,tzero,-1,      &
                & c40nam,xth)
!
!.....qv:
!
      do i = 1 , iym1
        do k = 1 , kz
          scrx1(k,i) = qva(i,k,jcx)/psa(i,jcx)
        end do
      end do
 
      write (c40nam,99005) jcx
      call mapsmp(scrx1,scrx1r,kz,iy,1,kz,kxn,1,iym1,iyn,0.D0,-1,c40nam,&
                & xth)
!
!....qc and qr:
!
      do i = 1 , iym1
        do k = 1 , kz
          scrx1(k,i) = qca(i,k,jcx)/psa(i,jcx)
        end do
      end do
      write (c40nam,99006) jcx
      call mapsmp(scrx1,scrx1r,kz,iy,1,kz,kxn,1,iym1,iyn,0.D0,-1,c40nam,&
                & xth)
!
!-----horizontal slices:
!
      if ( jyear.eq.jyear0 .and. ktau.eq.0 ) then
        c40nam = '       latitude at cross points        '
        call mapsmp(xlat,xlatr,iy,jx,1,iym1,5,1,jxm1,5,0.D0,1,c40nam,   &
                    & xth)
        c40nam = '       longitude at cross points       '
        call mapsmp(xlong,xlongr,iy,jx,1,iym1,5,1,jxm1,5,0.D0,1,c40nam, &
                  & xth)
        do j = 1 , jxm1
          do i = 1 , iym1
            hscr1(i,j) = ht(i,j)/gti
          end do
        end do
        c40nam = '     terrain height at cross points     '
        call mapsmp(hscr1,hscr1r,iy,jx,1,iym1,5,1,jxm1,5,0.D0,1,c40nam, &
                  & xth)
        c40nam = '    coriolis parameter at dot points    '
        call mapsmp(f,fr,iy,jx,1,iy,5,1,jx,5,0.D0,1,c40nam,xth)
        c40nam = '    map scale factor at dot points     '
        call mapsmp(msfd,msfdr,iy,jx,1,iy,5,1,jx,5,0.D0,1,c40nam,xth)
      end if
!
!....surface pressure
!
      const = -r8pt
      c40nam = '         surface pressure (cb)         '
      call mapsmp(psa,psar,iy,jx,1,iym1,iyn,1,jxm1,jxn,const,1,c40nam,  &
                  & xth)
!
!.....u and v:
!
      do j = 2 , jxm1
        do i = 2 , iym1
          psabar = 0.25*(psa(i,j)+psa(i,j-1)+psa(i-1,j)+psa(i-1,j-1))
          hscr1(i,j) = ua(i,kout,j)/psabar
          hscr2(i,j) = va(i,kout,j)/psabar
        end do
      end do
!
      if ( iexec.eq.1 ) then
!
!.....first time:
!
        do i = 2 , iym1
          psab1 = 0.5*(psa(i,1)+psa(i-1,1))
          psab2 = 0.5*(psa(i,jxm1)+psa(i-1,jxm1))
          hscr1(i,1) = ua(i,kout,1)/psab1
          hscr1(i,jx) = ua(i,kout,jx)/psab2
          hscr2(i,1) = va(i,kout,1)/psab1
          hscr2(i,jx) = va(i,kout,jx)/psab2
        end do
        do j = 2 , jxm1
          psab1 = 0.5*(psa(1,j)+psa(1,j-1))
          psab2 = 0.5*(psa(iym1,j)+psa(iym1,j-1))
          hscr1(1,j) = ua(1,kout,j)/psab1
          hscr1(iy,j) = ua(iy,kout,j)/psab2
          hscr2(1,j) = va(1,kout,j)/psab1
          hscr2(iy,j) = va(iy,kout,j)/psab2
        end do
        hscr1(1,1) = ua(1,kout,1)/psa(1,1)
        hscr1(1,jx) = ua(1,kout,jx)/psa(1,jxm1)
        hscr1(iy,1) = ua(iy,kout,1)/psa(iym1,1)
        hscr1(iy,jx) = ua(iy,kout,jx)/psa(iym1,jxm1)
        hscr2(1,1) = va(1,kout,1)/psa(1,1)
        hscr2(1,jx) = va(1,kout,jx)/psa(1,jxm1)
        hscr2(iy,1) = va(iy,kout,1)/psa(iym1,1)
        hscr2(iy,jx) = va(iy,kout,jx)/psa(iym1,jxm1)
!
      else if ( iexec.gt.1 ) then
!
        do i = 1 , iy
          hscr1(i,1) = uj1(i,kout)
          hscr1(i,jx) = ujl(i,kout)
          hscr2(i,1) = vj1(i,kout)
          hscr2(i,jx) = vjl(i,kout)
        end do
        do j = 2 , jxm1
          hscr1(1,j) = ui1(kout,j)
          hscr1(iy,j) = uil(kout,j)
          hscr2(1,j) = vi1(kout,j)
          hscr2(iy,j) = vil(kout,j)
        end do
!
      else        !end 2nd iexec test
      end if
 
      write (c40nam,99007) kout
      call mapsmp(hscr1,hscr1r,iy,jx,1,iy,iyn,1,jx,jxn,0.D0,1,c40nam,   &
                & xth)
      write (c40nam,99008) kout
      call mapsmp(hscr2,hscr2r,iy,jx,1,iy,iyn,1,jx,jxn,0.D0,1,c40nam,   &
                & xth)
!
!.....t:
!
      do j = 1 , jxm1
        do i = 1 , iym1
          hscr1(i,j) = ta(i,kout,j)/psa(i,j)
        end do
      end do
      write (c40nam,99009) kout
      call mapsmp(hscr1,hscr1r,iy,jx,1,iym1,iyn,1,jxm1,jxn,tzero,1,     &
                & c40nam,xth)
!
!.....relative humidity:
!
      do j = 1 , jxm1
        do i = 1 , iym1
          pres = a(kout)*psa(i,j) + r8pt
          if ( hscr1(i,j).gt.tzero ) then
!           v8 svp formula
            satvp = svp1*dexp(svp2*(hscr1(i,j)-tzero)/(hscr1(i,j)-svp3))
          else
            satvp = svp4*dexp(svp5-svp6/hscr1(i,j))
          end if
 
          qvs = ep2*satvp/(pres-satvp)
          hscr2(i,j) = (qva(i,kout,j)/psa(i,j))/qvs
        end do
      end do
      write (c40nam,99010) kout
      call mapsmp(hscr2,hscr2r,iy,jx,1,iym1,iyn,1,jxm1,jxn,0.D0,1,      &
                & c40nam,xth)
!
!.....qc and qr:
!
      do j = 1 , jxm1
        do i = 1 , iym1
          hscr1(i,j) = qca(i,kout,j)/psa(i,j)
        end do
      end do
      write (c40nam,99011) kout
      call mapsmp(hscr1,hscr1r,iy,jx,1,iym1,iyn,1,jxm1,jxn,0.D0,1,      &
                & c40nam,xth)
!
!.....surface temperature:
!
      if ( ibltyp.ne.0 ) then
        const = tzero
        c40nam = ' ground temperature (c)    '
        call mapsmp(tga,tgar,iy,jx,1,iym1,iyn,1,jxm1,jxn,const,1,c40nam,&
                  & xth)
      end if
!
!.....precipitation:
!
      c40nam = '   convective rainfall (mm)          '
      call mapsmp(rainc,raincr,iy,jx,1,iym1,iyn,1,jxm1,jxn,0.D0,1,      &
                  & c40nam,xth)
      c40nam = '  nonconvective rainfall (mm)      '
      call mapsmp(rainnc,rainncr,iy,jx,1,iym1,iyn,1,jxm1,jxn,0.D0,1,    &
                & c40nam,xth)
      do j = 1 , jxm1
        do i = 1 , iym1
          hscr1(i,j) = rainc(i,j) + rainnc(i,j)
        end do
      end do
      c40nam = '  total rainfall (mm)        '
      call mapsmp(hscr1,hscr1r,iy,jx,1,iym1,iyn,1,jxm1,jxn,0.D0,1,      &
                  & c40nam,xth)
99001 format (///1x,'--------------------------------------------------'&
            & ,/1x,'*****',4x,'large domain at t = ',f17.5,             &
             &' minutes, ktau = ',i7,' in year=',i4,3x,'*****'/1x,      &
         &'------------------------------------------------------------'&
        & ///)
99002 format ('  u  cross-section at j = ',i3,8x)
99003 format ('  v  cross-section at j = ',i3,8x)
99004 format ('  t  cross-section at j = ',i3,8x)
99005 format ('  qv  cross-section at j = ',i3,7x)
99006 format ('  qc  cross-section at j = ',i3,7x)
99007 format ('  u   at k = ',i3,24x)
99008 format ('  v   at k = ',i3,24x)
99009 format ('  t (c)  at k = ',i3,21x)
99010 format ('  relative humidity at k = ',i3,10x)
99011 format ('  qc  at k = ',i3,24x)
#endif
!
      end subroutine outprt

#endif
!
!
!
      subroutine mapsmp(fld,fldr,iyy,jxx,ia,ib,iny,ja,jb,jnx,const,     &
                      & ichos,c40nam,time)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!                                                                     c
!     this subroutine prints a sample of a two-dimensional data field c
!     on the line printer with 5 significant digits.                  c
!                                                                     c
!     *** note *** the values of fld(i,j) should be limited within    c
!                  1.e30 --- 1.e-30. if the value outside this        c
!                  range is desired, the program should be changed    c
!                  accordingly (in do loop 20).                       c
!                                                                     c
!                                                                     c
!     fld    : a two-dimensional array to hold the data field to be   c
!              sampled and printed. fld could be a horizontal slice,  c
!              fld(i,j), or a vertical slice fld(k,i) or fld(k,j).    c
!
!     fldr   : reverse array of fld; i.e., fldr(j,i)=fld(i,j)
!                                                                     c
!     iyy    : the first dimension of fld.                            c
!              for the horizontal slice, iyy is the dimension in the  c
!                                        y direction.                 c
!              for the vertical slice, iyy is the dimension in the    c
!                                      z direction.                   c
!                                                                     c
!     jxx    : the second dimension of fld.                           c
!              for the horizontal slice, jxx is the dimension in the  c
!                                        x direction.                 c
!              for the vertical slice, jxx is the dimension in either c
!                                      the x or y direction.          c
!                                                                     c
!     ia     : initial sampling point in the first dimension.         c
!                                                                     c
!     ib     : final sampling point in the first dimension.           c
!                                                                     c
!     iny    : sampling interval in the first dimension.              c
!                                                                     c
!     ja     : initial sampling point in the second dimension.        c
!                                                                     c
!     jb     : final sampling point in the second dimension.          c
!                                                                     c
!     jnx    : sampling interval in the second dimension.             c
!                                                                     c
!     const  : constant used to be subtracted from fldr.              c
!                                                                     c
!     ichos > 0 : for horizontal array fld(y,x)                       c
!           < 0 : for vertical cross section fld(z,y) or fld(z,x)     c
!                                                                     c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      implicit none
!
! Dummy arguments
!
      character(40) :: c40nam
      real(8) :: const , time
      integer :: ia , ib , ichos , iny , iyy , ja , jb , jnx , jxx
      real(8) , dimension(iyy,jxx) :: fld
      real(8) , dimension(jxx,iyy) :: fldr
      intent (in) c40nam , const , fld , ia , ib , ichos , iny , iyy ,  &
                & ja , jb , jnx , jxx , time
      intent (inout) fldr
!
! Local variables
!
      real(8) :: fldl , fldmax , fldmin , fldu
      integer :: i , i1 , i2 , iexp , ir , it , iy , j , j1 , j2 , j2n ,&
               & j3 , jj , jl , jn , jn1 , jt , jtn , k1 , k2 , k3 ,    &
               & k4 , ksigt , n , n1
      character(24) :: ifmt1 , ifmt2
      integer , dimension(100) :: jm
!
      data ksigt/5/
!
      do i = 1 , iyy
        do j = 1 , jxx
          fldr(j,i) = fld(i,j)
        end do
      end do
!
      n = 6
      k1 = ksigt + 2
      k2 = 124/k1
      k3 = ksigt/2
      k4 = ksigt - k3
!
      do i = ia , ib , iny
        do j = ja , jb , jnx
          fldr(j,i) = fldr(j,i) - const
        end do
      end do
!
      fldmax = fldr(ja,ia)
      fldmin = fldr(ja,ia)
      fldu = 10.**ksigt
      fldl = 10.**(ksigt-1)
      do i = ia , ib , iny
        do j = ja , jb , jnx
          if ( dabs(fldr(j,i)).le.1.E30 .and. dabs(fldr(j,i))           &
             & .ge.1.E-30 ) then
            if ( dabs(fldr(j,i)).gt.fldmax ) fldmax = dabs(fldr(j,i))
            if ( dabs(fldr(j,i)).lt.fldmin ) fldmin = dabs(fldr(j,i))
          end if
        end do
      end do
!
      if ( fldmax.ne.fldmin ) then
        iexp = 0
        do n1 = 1 , 500
          if ( fldmax.ge.fldu ) then
            fldmax = fldmax/10.
            iexp = iexp - 1
          else if ( fldmax.lt.fldl ) then
            fldmax = fldmax*10.
            iexp = iexp + 1
          else
            exit
          end if
        end do
        do i = ia , ib , iny
          do j = ja , jb , jnx
            fldr(j,i) = fldr(j,i)*10.**iexp
          end do
        end do
        iy = ib - ia + 1
        jn = k2*jnx
        jn1 = jn - 1
        write (n,99001) c40nam , iexp , time
        do j1 = ja , jb , jn
          jl = min0(j1+jn1,jb)
          jt = jl - j1 + 1
          jtn = (jt-1)/jnx + 1
          j2n = 0
          do j2 = 1 , jt , jnx
            j2n = j2n + 1
            jm(j2n) = j1 + j2 - 1
          end do
          write (ifmt1,99002) jtn , k4 , k3
          write (n,ifmt1) (jm(jj),jj=1,j2n)
          write (ifmt2,99003) jtn , k1
!110      format(1x,i2,1x,i2)
          it = (iy-1)/iny
          ir = iy - it*iny
          do i2 = ia , ib , iny
            i1 = ib + ia - i2 - ir + 1
            if ( ichos.lt.0 ) i1 = i2
            write (n,ifmt2) i1 , (fldr(j3,i1),j3=j1,jl,jnx) , i1
          end do
          write (n,ifmt1) (jm(jj),jj=1,j2n)
        end do
        do i = ia , ib , iny
          do j = ja , jb , jnx
            fldr(j,i) = fldr(j,i)/(10.**iexp) + const
          end do
        end do
      else
        do i = ia , ib , iny
          do j = ja , jb , jnx
            fldr(j,i) = fldr(j,i) + const
          end do
        end do
        write (n,99004) c40nam , fldmax , time
      end if      !end if(fldmax.ne.fldmin)test
99001 format (////' this is a list of  ',a40,'  ,scaled by  1.e',i3,5x, &
             &'at time = ',f10.3)
99002 format ('(/4x,',i2,'(',i2,'x,i2,',i2,'x)/)')
99003 format ('(1x,i2,',i2,'f',i2,'.0,2x,i2)')
99004 format (/'   all of the values of ',a40,' are equal to ',e15.5,5x,&
             &'at time = ',f10.3)
!
      end subroutine mapsmp
!
      end module mod_outprt
