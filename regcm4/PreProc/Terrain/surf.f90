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

      subroutine surf(xlat,xlon,lnduse,iy,jx,incr,dsgrid,     &
                    & lndout,land,nrec,h2opct,     &
                    & lsmtyp,sanda,sandb,claya,clayb,frac_lnd,nveg,     &
                    & aertyp,intext,texout,frac_tex,ntex)
      use mod_block
      use mod_interp , only : bint
      implicit none
!
! Dummy arguments
!
      character(7) :: aertyp
      real(4) :: dsgrid , h2opct
      integer :: incr , iy , jx , nrec , ntex , nveg
      character(4) :: lsmtyp
      real(4) , dimension(iy,jx) :: claya , clayb , lndout , sanda ,    &
                                  & sandb , texout , xlat , xlon
      real(4) , dimension(iy,jx,nveg) :: frac_lnd
      real(4) , dimension(iy,jx,ntex) :: frac_tex
      integer , dimension(iy,jx) :: intext , lnduse
      real(4) , dimension(iy,jx,2) :: land
      intent (in) aertyp , dsgrid , h2opct , incr , &
                & iy , jx , lsmtyp , nrec , ntex , nveg , xlat
      intent (out) claya , clayb , frac_lnd , frac_tex , intext ,       &
                 & lnduse , sanda , sandb
      intent (inout) land , lndout , texout , xlon
!
! Local variables
!
      logical :: flag
      integer :: i , ii , iindex , ilev , j , jindex , jj , k , lrec ,  &
               & lengdo , nbase
      real(4) , dimension(iy,jx,2) :: itex
      real(8) :: xx , yy
!
!---------------------------------------------------------------------
!
!     imx,jmx must correspond to iy,jx in the master input file;
!     otherwise the program will abort.
!
!-----------------------------------------------------------------------
!
      flag = .true.
      lrec = 0
      do k = 1 , 2
        do j = 1 , jx
          do i = 1 , iy
            land(i,j,k) = 0
            itex(i,j,k) = 0
          end do
        end do
      end do
!
      if ( aertyp(7:7)=='1' ) then
        if ( lsmtyp=='BATS' ) then
          lengdo = nveg + ntex
        else if ( lsmtyp=='USGS' ) then
          lengdo = nveg + 4 + ntex
        else
        end if
      else if ( lsmtyp=='BATS' ) then
        lengdo = nveg
      else if ( lsmtyp=='USGS' ) then
        lengdo = nveg + 4
      else
      end if
      do ilev = 1 , lengdo
        rewind (48)
        do lrec = 1 , nrec
          read (48) stores
          jindex = (stores(2)-grdlnmn)*incr + 1.1
          iindex = (stores(1)-grdltmn)*incr + 1.1
          if ( iindex>iter .or. jindex>jter ) then
            print 99001 , iindex , jindex , lrec , stores(1) , stores(2)&
                & , grdltmn , grdlnmn
            stop 60
          end if
          lnd8(iindex,jindex) = stores(ilev+4)
        end do
!
        if ( ilev<=nveg ) then
          do ii = 1 , iy
            do jj = 1 , jx
              yy = -(grdltmn-xlat(ii,jj))/dsgrid + 1.0
              if ( grdlnmn<=-180.0 .and. xlon(ii,jj)>0.0 ) xlon(ii,jj)  &
                 & = xlon(ii,jj) - 360.
              xx = -(grdlnmn-xlon(ii,jj))/dsgrid + 1.0
              lndout(ii,jj) = bint(yy,xx,lnd8,iter,jter,flag)
              frac_lnd(ii,jj,ilev) = lndout(ii,jj)
!
!             note: it is desirable to force grid boxes with less
!             than 75 percent water to be a land category,
!             even if water is the largest single category.
!
              if ( .not.(lsmtyp=='BATS' .and. (ilev==14 .or. ilev==15)  &
                 & .and. lndout(ii,jj)<h2opct) ) then
                if ( lsmtyp/='USGS' .or. ilev/=25 .or. lndout(ii,jj)    &
                   & >=h2opct ) then
 
                  if ( lndout(ii,jj)>land(ii,jj,1) ) then
                    land(ii,jj,1) = lndout(ii,jj)
                    land(ii,jj,2) = ilev
                  end if
                end if
              end if
            end do
          end do
        else if ( ((lsmtyp=='USGS' .and. ilev>nveg+4) .or.              &
                 &(lsmtyp=='BATS' .and. ilev>nveg)) .and. aertyp(7:7)   &
                 &=='1' ) then
          if ( lsmtyp=='BATS' ) then
            nbase = nveg
          else if ( lsmtyp=='USGS' ) then
            nbase = nveg + 4
          else
          end if
          do ii = 1 , iy
            do jj = 1 , jx
              yy = -(grdltmn-xlat(ii,jj))/dsgrid + 1.0
              if ( grdlnmn<=-180.0 .and. xlon(ii,jj)>0.0 ) xlon(ii,jj)  &
                 & = xlon(ii,jj) - 360.
              xx = -(grdlnmn-xlon(ii,jj))/dsgrid + 1.0
              texout(ii,jj) = bint(yy,xx,lnd8,iter,jter,flag)
              frac_tex(ii,jj,ilev-nbase) = texout(ii,jj)
!
!             note: it is desirable to force grid boxes with less
!             than 75 percent water to be a land category,
!             even if water is the largest single category.
!
              if ( lsmtyp/='BATS' .or. ilev/=nveg+14 .or. texout(ii,jj) &
                 & >=h2opct ) then
                if ( lsmtyp/='USGS' .or. ilev/=nveg+18 .or.             &
                   & texout(ii,jj)>=h2opct ) then
                  if ( texout(ii,jj)>itex(ii,jj,1) ) then
                    itex(ii,jj,1) = texout(ii,jj)
                    if ( lsmtyp=='BATS' ) then
                      itex(ii,jj,2) = ilev - nbase
                    else if ( lsmtyp=='USGS' ) then
                      itex(ii,jj,2) = ilev - nbase
                    else
                    end if
                  end if
                end if
              end if
            end do
          end do
        else if ( lsmtyp=='USGS' ) then
          do ii = 1 , iy
            do jj = 1 , jx
              yy = -(grdltmn-xlat(ii,jj))/dsgrid + 1.0
              if ( grdlnmn<=-180.0 .and. xlon(ii,jj)>0.0 ) xlon(ii,jj)  &
                 & = xlon(ii,jj) - 360.
              xx = -(grdlnmn-xlon(ii,jj))/dsgrid + 1.0
              if ( ilev==nveg+1 ) then
                sanda(ii,jj) = bint(yy,xx,lnd8,iter,jter,flag)
              else if ( ilev==nveg+2 ) then
                sandb(ii,jj) = bint(yy,xx,lnd8,iter,jter,flag)
              else if ( ilev==nveg+3 ) then
                claya(ii,jj) = bint(yy,xx,lnd8,iter,jter,flag)
              else if ( ilev==nveg+4 ) then
                clayb(ii,jj) = bint(yy,xx,lnd8,iter,jter,flag)
              else
              end if
            end do
          end do
        else
        end if
      end do
!
      do i = 1 , iy
        do j = 1 , jx
          lndout(i,j) = land(i,j,2)
          lnduse(i,j) = int(land(i,j,2))
          if ( aertyp(7:7)=='1' ) then
            texout(i,j) = itex(i,j,2)
            intext(i,j) = int(itex(i,j,2))
          end if
        end do
      end do
!
!
!-----grid the data.  grdltmn=minimum latitude of incoming data.
!-----grdlnmn = minimum longitude of incoming data.  point(1,1)
!-----is value at (grdltmn,grdlnmn)
!
99001 format (1x,'*** iindex = ',i3,'   jindex = ',i3,'   lrec = ',i5,  &
             &'   lat = ',f10.3,3x,'lon = ',f10.3,3x,'grdltmn = ',f10.3,&
            & 5x,'grdlnmn = ',f10.3)
!
      end subroutine surf
