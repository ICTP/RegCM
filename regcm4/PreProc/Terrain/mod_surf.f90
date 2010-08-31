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

      module mod_surf

      contains

      subroutine surf(xlat,xlon,lnduse,iy,jx,incr,dsgrid,lndout,land,   &
                    & nrec,h2opct,nveg,aertyp,intext,texout,frac_tex,   &
                    & ntex)
      use mod_block
      use mod_interp , only : bint
      implicit none
!
! Dummy arguments
!
      character(7) :: aertyp
      real(4) :: dsgrid , h2opct
      integer :: incr , iy , jx , nrec , ntex , nveg
      real(4) , dimension(iy,jx) :: lndout , texout , xlat , xlon
      real(4) , dimension(iy,jx,ntex) :: frac_tex
      integer , dimension(iy,jx) :: intext , lnduse
      real(4) , dimension(iy,jx,2) :: land
      intent (in) aertyp , dsgrid , h2opct , incr , &
                & iy , jx , nrec , ntex , nveg , xlat
      intent (out) frac_tex , intext , lnduse
      intent (inout) land , lndout , texout , xlon
!
! Local variables
!
      logical :: flag
      integer :: i , ii , iindex , ilev , j , jindex , jj , lrec ,  &
               & lengdo , nbase , mxi , mxj
      real(4) , dimension(iy,jx,2) :: itex
      real(8) :: xx , yy , rinc , xd , yd , wc , wm
      logical , dimension(:,:) , allocatable :: mask
      real(8) , dimension(:,:) , allocatable :: lnd8
!
      flag = .true.
      lrec = 0
      land = 0.0
      itex = 0.0
      rinc = incr
!
      print *, 'Input data point MIN is at ', grdltmn , grdlnmn
      print *, 'Input data point MAX is at ', grdltma , grdlnma
      print *, 'Input data resolution is   ', dsgrid
      if (lcrosstime) then
        mxj = nint((mod((grdlnma+360.0),360.0)-grdlnmn)*rinc) + 1
      else
        mxj = nint((grdlnma-grdlnmn)*rinc) + 1
      end if
      mxi = nint((grdltma-grdltmn)*rinc) + 1
      print *, 'Allocating ', mxi, 'x', mxj

      if ( aertyp(7:7)=='1' ) then
        lengdo = nveg + ntex
      else
        lengdo = nveg
      end if

      if (lonwrap) then
        mxj = mxj + 4
      end if

      allocate(mask(mxi,mxj))
      allocate(lnd8(mxi,mxj))

      do ilev = 1 , lengdo
        mask(:,:) = .false.
        rewind (48)
        do lrec = 1 , nrec
          read (48) stores
          if (lcrosstime) then
            jindex = nint((mod((stores(2)+360.0),360.0)-grdlnmn)*       &
                           rinc) + 1
          else
            jindex = nint((stores(2)-grdlnmn)*rinc) + 1
          end if
          iindex = nint((stores(1)-grdltmn)*rinc) + 1
          if ( iindex < 1 .or. iindex>iter .or. &
               jindex < 1 .or. jindex>jter ) then
            print 99001 , iindex , jindex , lrec , stores(1) ,          &
                & stores(2) , grdltmn , grdlnmn
            stop 60
          end if
          lnd8(iindex,jindex) = stores(ilev+4)
          mask(iindex,jindex) = .true.
        end do
!
!       Fill the matrix using nearest point to missing one.
!
        do ii = 1 , mxi-1
          do jj = 1 , mxj-1
            if (.not. mask(ii,jj)) then
              yy = ((ii*dsgrid)-grdltmn)*rinc + 1.0D+00
              xx = ((jj*dsgrid)-grdlnmn)*rinc + 1.0D+00
              wc = 999.0
              wm = 999.0
              if (ii > 1) then
                if (mask(ii-1,jj)) then
                  yd = (((ii-1)*dsgrid)-grdltmn)*rinc + 1.0D+00
                  wc = (yy-yd)
                  lnd8(ii,jj) = lnd8(ii-1,jj)
                  wm = wc
                end if
              end if
              if (mask(ii+1,jj)) then
                yd = (((ii+1)*dsgrid)-grdltmn)*rinc + 1.0D+00
                wc = (yy-yd)
                if (wc < wm) then
                  lnd8(ii,jj) = lnd8(ii+1,jj)
                  wm = wc
                end if
              end if
              if (jj > 1) then
                if (mask(jj > 1 .and. ii,jj-1)) then
                  xd = (((jj-1)*dsgrid)-grdltmn)*rinc + 1.0D+00
                  wc = (xx-xd)
                  if (wc < wm) then
                    lnd8(ii,jj) = lnd8(ii,jj-1)
                    wm = wc
                  end if
                end if
              end if
              if (mask(ii,jj+1)) then
                xd = (((jj+1)*dsgrid)-grdltmn)*rinc + 1.0D+00
                wc = (xx-xd)
                if (wc < wm) then
                  lnd8(ii,jj) = lnd8(ii,jj+1)
                  wm = wc
                end if
              end if
            end if
          end do
        end do
!
!       Special case for wrapping of longitude
!
        if (lonwrap) then
          lnd8 = cshift(lnd8,shift=-2,dim=2)
          lnd8(:,1) = lnd8(:,mxj-3)
          lnd8(:,2) = lnd8(:,mxj-4)
          lnd8(:,mxj-2) = lnd8(:,3)
          lnd8(:,mxj-1) = lnd8(:,4)
        end if

        if ( ilev<=nveg ) then
          do ii = 1 , iy
            do jj = 1 , jx
              yy = (xlat(ii,jj)-grdltmn)*rinc + 1.0D+00
              if (lcrosstime) then
                xx = (mod((xlon(ii,jj)+360.0),360.0)-grdlnmn)*rinc + &
                      1.0D+00
              else
                xx = (xlon(ii,jj)-grdlnmn)*rinc + 1.0D+00
              end if
              lndout(ii,jj) = bint(yy,xx,lnd8,mxi,mxj,flag)
!
!             note: it is desirable to force grid boxes with less
!             than 75 percent water to be a land category,
!             even if water is the largest single category.
!
              if ( .not.((ilev==14 .or. ilev==15) .and.                 &
                     &    lndout(ii,jj)<h2opct) ) then
                if ( ilev/=25 .or. lndout(ii,jj) >=h2opct ) then
                  if ( lndout(ii,jj)>land(ii,jj,1) ) then
                    land(ii,jj,1) = lndout(ii,jj)
                    land(ii,jj,2) = ilev
                  end if
                end if
              end if
            end do
          end do
        else if ( ilev>nveg .and. aertyp(7:7) == '1' ) then
          nbase = nveg
          do ii = 1 , iy
            do jj = 1 , jx
              yy = (xlat(ii,jj)-grdltmn)*rinc + 1.0D+00
              if (lcrosstime) then
                xx = (mod((xlon(ii,jj)+360.0),360.0)-grdlnmn)*rinc + &
                     1.0D+00
              else
                xx = (xlon(ii,jj)-grdlnmn)*rinc + 1.0D+00
              end if
              texout(ii,jj) = bint(yy,xx,lnd8,mxi,mxj,flag)
              frac_tex(ii,jj,ilev-nbase) = texout(ii,jj)
!
!             note: it is desirable to force grid boxes with less
!             than 75 percent water to be a land category,
!             even if water is the largest single category.
!
              if ( ilev/=nveg+14 .or. texout(ii,jj) >=h2opct ) then
                if ( ilev/=nveg+18 .or. texout(ii,jj)>=h2opct ) then
                  if ( texout(ii,jj)>itex(ii,jj,1) ) then
                    itex(ii,jj,1) = texout(ii,jj)
                    itex(ii,jj,2) = ilev - nbase
                  end if
                end if
              end if
            end do
          end do
        end if
      end do
!
      deallocate(mask)
      deallocate(lnd8)
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
99001 format (1x,'*** iindex = ',i4,'   jindex = ',i4,'   lrec = ',i5,  &
             &'   lat = ',f10.3,3x,'lon = ',f10.3,3x,'grdltmn = ',f10.3,&
            & 5x,'grdlnmn = ',f10.3)
!
      end subroutine surf

      end module mod_surf
