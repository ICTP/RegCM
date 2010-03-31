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

      use mod_regcm_param
      use mod_param1
      use mod_param2
      use mod_param3 , only : jxsex , kxout , ptop , a
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
      const = -ptop
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
          pres = a(kout)*psa(i,j) + ptop
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
!
      end subroutine outprt

#endif
