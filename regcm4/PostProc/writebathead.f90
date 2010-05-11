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

      subroutine writebathead(f,xmap,dmap,xlat,xlon,dlat,dlon,zs,ls,    &
                            & vvarmin,vvarmax,xlat1d,xlon1d,iadm,ndim,  &
                            & idout,xhro,iotyp)
 
      use mod_dynparam , only : jxm2 , iym2
      use mod_bats_param
      implicit none
!
! Dummy arguments
!
      integer :: idout , iotyp , ndim
      real(8) :: xhro
      real(4) , dimension(jxm2,iym2) :: dlat , dlon , dmap , f , ls ,   &
                             & xlat , xlon , xmap , zs
      integer , dimension(ndim) :: iadm
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(iym2) :: xlat1d
      real(4) , dimension(jxm2) :: xlon1d
      intent (in) dlat , dlon , dmap , f , ls , ndim , xlat , xlon ,    &
                & xmap , zs
!
! Local variables
!
      real(4) :: aaa , fact , misdat , offset , vmax , vmin , vmisdat , &
            & xmax , xmin
      integer :: i , j
      integer , dimension(jxm2,iym2) :: icol , ils , itex
      character(64) :: lname
      real(4) , dimension(jxm2,iym2) :: tmp2d
      character(64) :: units
      character(64) :: varnam
!
      iadm(3) = 1
      aaa = 2.**16. - 1.
!     CALL XTDOT(zs,zs,jxm2,iym2,1,jxm2-1,iym2-1)
!     CALL XTDOT(f,f,jxm2,iym2,1,jxm2-1,iym2-1)
      do j = 1 , iym2
        do i = 1 , jxm2
          ils(i,j) = nint(ls(i,j))
          itex(i,j) = iexsol(ils(i,j))
          icol(i,j) = kolsol(ils(i,j))
        end do
      end do
 
      varnam = 'ZB'
      print * , varnam
      lname = 'Terrain Elevation'
      units = 'm'
      xmax = 7000.
      xmin = -100.
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 3 , iym2 - 2
        do i = 3 , jxm2 - 2
          tmp2d(i-2,j-2) = zs(i,j)
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(zb)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(zb)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'LU'
      print * , varnam
      lname = 'Land Use Type'
      units = 'unitless'
      xmax = 21.
      xmin = -1.
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 3 , iym2 - 2
        do i = 3 , jxm2 - 2
          tmp2d(i-2,j-2) = ls(i,j)
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(lu)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(lu)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'F'
      print * , varnam
      lname = 'Coriolus'
      units = 'rad/sec'
      xmax = 0.001
      xmin = -0.001
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 3 , iym2 - 2
        do i = 3 , jxm2 - 2
          tmp2d(i-2,j-2) = f(i,j)
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(f)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(f)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'XMAP'
      print * , varnam
      lname = 'Cross-Grid Map Fact'
      units = 'unitless'
      xmax = 2.0
      xmin = 0.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 3 , iym2 - 2
        do i = 3 , jxm2 - 2
          tmp2d(i-2,j-2) = xmap(i,j)
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(xmap)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(xmap)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'DMAP'
      print * , varnam
      lname = 'Dot Grid Map Factor'
      units = 'degrees'
      xmax = 2.0
      xmin = 0.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 3 , iym2 - 2
        do i = 3 , jxm2 - 2
          tmp2d(i-2,j-2) = dmap(i,j)
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(dmap)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(dmap)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'XLAT'
      print * , varnam
      lname = 'Cross Grid Latitude'
      units = 'degrees'
      xmax = 100.0
      xmin = -100.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 3 , iym2 - 2
        do i = 3 , jxm2 - 2
          tmp2d(i-2,j-2) = xlat(i,j)
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(xlat)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(xlat)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'XLON'
      print * , varnam
      lname = 'Cross Grid Longitude'
      units = 'degrees'
      xmax = 200.0
      xmin = -200.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 3 , iym2 - 2
        do i = 3 , jxm2 - 2
          tmp2d(i-2,j-2) = xlon(i,j)
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(xlon)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(xlon)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'DLAT'
      print * , varnam
      lname = 'Dot Grid Latitude'
      units = 'degrees'
      xmax = 100.0
      xmin = -100.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 3 , iym2 - 2
        do i = 3 , jxm2 - 2
          tmp2d(i-2,j-2) = dlat(i,j)
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(dlat)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(dlat)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'DLON'
      print * , varnam
      lname = 'Dot Grid Longitude'
      units = 'degrees'
      xmax = 200.0
      xmin = -200.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 3 , iym2 - 2
        do i = 3 , jxm2 - 2
          tmp2d(i-2,j-2) = dlon(i,j)
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(dlon)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(dlon)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'land'
      print * , varnam
      lname = 'Land Mask'
      units = 'unitless'
      xmax = 1.5
      xmin = -0.5
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          if ( ils(i,j)==15 ) then
            tmp2d(i,j) = -999.
          else
            tmp2d(i,j) = 1.0
          end if
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(land)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(land)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'ROOTF'
      print * , varnam
      lname = 'Root Zone Root Frac'
      units = 'fraction'
      xmax = 1.5
      xmin = -0.5
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = rootf(ils(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(rootf)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(rootf)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'VEGC'
      print * , varnam
      lname = 'Max Vegetation Cover'
      units = 'fraction'
      xmax = 1.5
      xmin = -0.5
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = vegc(ils(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(vegc)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(vegc)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'SEASF'
      print * , varnam
      lname = 'VEGC(MAX)-VEGC @269K'
      units = 'fraction'
      xmax = 1.5
      xmin = -0.5
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = seasf(ils(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(seasf)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(seasf)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'ROUGH'
      print * , varnam
      lname = 'Roughness Length'
      units = 'm'
      xmax = 3.5
      xmin = -0.5
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = rough(ils(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(rough)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(rough)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'DISPLA'
      print * , varnam
      lname = 'Displacement Height'
      units = 'm'
      xmax = 30.0
      xmin = -2.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = displa(ils(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(displa)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(displa)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'RSMIN'
      print * , varnam
      lname = 'Min Stomatl Resist'
      units = 's/m'
      xmax = 300.0
      xmin = -10.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = rsmin(ils(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(rsmin)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(rsmin)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'XLA'
      print * , varnam
      lname = 'Max Leaf Area Index'
      units = 'unitless'
      xmax = 10.0
      xmin = -1.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = xla(ils(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(xla)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(xla)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'XLAI0'
      print * , varnam
      lname = 'Min Leaf Area Index'
      units = 'unitless'
      xmax = 10.0
      xmin = -1.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = xlai0(ils(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(xlai0)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(xlai0)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'SAI'
      print * , varnam
      lname = 'Stem Area Index'
      units = 'unitless'
      xmax = 10.0
      xmin = -1.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = sai(ils(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(sai)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(sai)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'SQRTDI'
      print * , varnam
      lname = 'Inv SQRT Leaf Dim'
      units = 'm-0.5'
      xmax = 50.0
      xmin = -10.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = sqrtdi(ils(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(sqrtdi)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(sqrtdi)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'FC'
      print * , varnam
      lname = 'Light Depend on RS'
      units = 'unitless'
      xmax = 50.0
      xmin = -10.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = fc(ils(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(fc)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(fc)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'DEPUV'
      print * , varnam
      lname = 'Top Soil Layer Depth'
      units = 'mm'
      xmax = 200.0
      xmin = -10.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = depuv(ils(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(depuv)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(depuv)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'DEPRV'
      print * , varnam
      lname = 'Root Soil Depth'
      units = 'mm'
      xmax = 4000.0
      xmin = -10.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = deprv(ils(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(deprv)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(deprv)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'DEPTV'
      print * , varnam
      lname = 'Total Soil Depth'
      units = 'mm'
      xmax = 8000.0
      xmin = -10.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = deptv(ils(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(deptv)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(deptv)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'IEXSOL'
      print * , varnam
      lname = 'Soil Texture'
      units = 'unitless'
      xmax = 20.0
      xmin = -1.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = iexsol(ils(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(iexsol)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(iexsol)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'KOLSOL'
      print * , varnam
      lname = 'Soil Color'
      units = 'unitless'
      xmax = 20.0
      xmin = -1.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = kolsol(ils(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(kolsol)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(kolsol)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'XMOPOR'
      print * , varnam
      lname = 'Soil Porosity'
      units = 'fraction'
      xmax = 1.5
      xmin = -1.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = xmopor(itex(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(xmopor)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(xmopor)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'XMOSUC'
      print * , varnam
      lname = 'Min Soil Suction'
      units = 'mm'
      xmax = 300.0
      xmin = -10.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = xmosuc(itex(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(xmosuc)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(xmosuc)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'XMOHYD'
      print * , varnam
      lname = 'Sat Soil Conductiv'
      units = 'mm/s'
      xmax = 1.0
      xmin = 0.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = xmohyd(itex(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(xmohyd)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(xmohyd)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'XMOWIL'
      print * , varnam
      lname = 'Soil Wilting Point'
      units = 'Fraction'
      xmax = 1.0
      xmin = 0.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = xmowil(itex(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(xmowil)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(xmowil)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'BEE'
      print * , varnam
      lname = 'Clapp Hornbereger'
      units = 'unitless'
      xmax = 20.0
      xmin = -1.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = bee(itex(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(bee)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(bee)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'SKRAT'
      print * , varnam
      lname = 'Thermal Conductivity'
      units = 'ratio'
      xmax = 20.0
      xmin = -1.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = skrat(itex(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(skrat)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(skrat)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'ALBVGS'
      print * , varnam
      lname = 'Veg Albedo < 0.7 um'
      units = 'fraction'
      xmax = 20.0
      xmin = -1.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = albvgs(itex(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(albvgs)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(albvgs)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      varnam = 'ALBVGL'
      print * , varnam
      lname = 'Veg Albedo > 0.7 um'
      units = 'fraction'
      xmax = 20.0
      xmin = -1.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      do j = 1 , iym2
        do i = 1 , jxm2
          tmp2d(i,j) = albvgl(itex(i,j))
        end do
      end do
      if ( iotyp==1 ) then
        call getminmax(tmp2d,jxm2,iym2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(albvgl)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(albvgl)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,tmp2d,jxm2,iym2,1,iadm,xhro,lname,     &
                & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,&
                & misdat,iotyp)
 
      end subroutine writebathead
