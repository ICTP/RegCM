!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine writehead(f,xmap,dmap,xlat,xlon,dlat,dlon,zs,zssd,ls,  &
                         & vvarmin,vvarmax,xlat1d,xlon1d,iadim,ndim,    &
                         & idout,xhro,iotyp)

      use mod_regcm_param , only : jxm2 , ixm2 

      implicit none
!
! Dummy arguments
!
      integer :: idout , iotyp , ndim
      real(8) :: xhro
      real(4) , dimension(jxm2,ixm2) :: dlat , dlon , dmap , f , ls ,   &
                               & xlat , xlon , xmap , zs , zssd
      integer , dimension(ndim) :: iadim
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(ixm2) :: xlat1d
      real(4) , dimension(jxm2) :: xlon1d
      intent (in) ndim
!
! Local variables
!
      real(4) :: aaa , fact , misdat , offset , vmax , vmin , vmisdat , &
            & xmax , xmin
      character(64) :: lname
      character(64) :: units
      character(64) :: varnam
!
      aaa = 2.**16. - 1.
!     CALL XTDOT(zs,zs,jxm2,ixm2,1,jxm2-1,ixm2-1)
!     CALL XTDOT(f,f,jxm2,ixm2,1,jxm2-1,ixm2-1)
 
      iadim(3) = 1
 
      varnam = 'HT'
      lname = 'Terrain Elevation'
      units = 'm'
      xmax = 7000.
      xmin = -100.
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      if ( iotyp==1 ) then
        call getminmax(zs,jxm2,ixm2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(zs)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(zs)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,zs,jxm2,ixm2,1,iadim,xhro,lname,units, &
                  & fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,    &
                  & misdat,iotyp)
 
      varnam = 'HTSD'
      lname = 'Elevation Std Dev'
      units = 'm'
      xmax = 5000.
      xmin = -100.
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      if ( iotyp==1 ) then
        call getminmax(zssd,jxm2,ixm2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(zssd)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(zssd)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,zssd,jxm2,ixm2,1,iadim,xhro,lname,     &
                  & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,&
                  & 0,misdat,iotyp)
 
      varnam = 'LU'
      lname = 'Land Use Type'
      units = 'unitless'
      xmax = 21.
      xmin = -1.
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      if ( iotyp==1 ) then
        call getminmax(ls,jxm2,ixm2,1,vmin,vmax,vmisdat)
        if ( vmin<xmin .or. vmax>xmax ) then
          print * , 'Values Out of Range:  FIELD=' , varnam
          print * , 'MINVAL(ls)=' , vmin , 'XMIN=' , xmin
          print * , 'MAXVAL(ls)=' , vmax , 'XMAX=' , xmax
          stop 999
        end if
        misdat = xmin
      else if ( iotyp==2 ) then
        misdat = vmisdat
      else
      end if
      call writecdf(idout,varnam,ls,jxm2,ixm2,1,iadim,xhro,lname,units, &
                  & fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,    &
                  & misdat,iotyp)
 
      varnam = 'F'
      lname = 'Coriolus'
      units = 'rad/sec'
      xmax = 0.001
      xmin = -0.001
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      if ( iotyp==1 ) then
        call getminmax(f,jxm2,ixm2,1,vmin,vmax,vmisdat)
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
      call writecdf(idout,varnam,f,jxm2,ixm2,1,iadim,xhro,lname,units,  &
                  & fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,0,    &
                  & misdat,iotyp)
 
      varnam = 'XMAP'
      lname = 'Cross-Grid Map Fact'
      units = 'unitless'
      xmax = 2.0
      xmin = 0.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      if ( iotyp==1 ) then
        call getminmax(xmap,jxm2,ixm2,1,vmin,vmax,vmisdat)
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
      call writecdf(idout,varnam,xmap,jxm2,ixm2,1,iadim,xhro,lname,     &
                  & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,&
                  & 0,misdat,iotyp)
 
      varnam = 'DMAP'
      lname = 'Dot Grid Map Factor'
      units = 'degrees'
      xmax = 2.0
      xmin = 0.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      if ( iotyp==1 ) then
        call getminmax(dmap,jxm2,ixm2,1,vmin,vmax,vmisdat)
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
      call writecdf(idout,varnam,dmap,jxm2,ixm2,1,iadim,xhro,lname,     &
                  & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,&
                  & 0,misdat,iotyp)
 
      varnam = 'XLAT'
      lname = 'Cross Grid Latitude'
      units = 'degrees'
      xmax = 100.0
      xmin = -100.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      if ( iotyp==1 ) then
        call getminmax(xlat,jxm2,ixm2,1,vmin,vmax,vmisdat)
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
      call writecdf(idout,varnam,xlat,jxm2,ixm2,1,iadim,xhro,lname,     &
                  & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,&
                  & 0,misdat,iotyp)
 
      varnam = 'XLON'
      lname = 'Cross Grid Longitude'
      units = 'degrees'
      xmax = 200.0
      xmin = -200.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      if ( iotyp==1 ) then
        call getminmax(xlon,jxm2,ixm2,1,vmin,vmax,vmisdat)
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
      call writecdf(idout,varnam,xlon,jxm2,ixm2,1,iadim,xhro,lname,     &
                  & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,&
                  & 0,misdat,iotyp)
 
      varnam = 'DLAT'
      lname = 'Dot Grid Latitude'
      units = 'degrees'
      xmax = 100.0
      xmin = -100.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      if ( iotyp==1 ) then
        call getminmax(dlat,jxm2,ixm2,1,vmin,vmax,vmisdat)
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
      call writecdf(idout,varnam,dlat,jxm2,ixm2,1,iadim,xhro,lname,     &
                  & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,&
                  & 0,misdat,iotyp)
 
      varnam = 'DLON'
      lname = 'Dot Grid Longitude'
      units = 'degrees'
      xmax = 200.0
      xmin = -200.0
      fact = (xmax-xmin)/aaa
      offset = (xmax+xmin)/2.
      if ( iotyp==1 ) then
        call getminmax(dlon,jxm2,ixm2,1,vmin,vmax,vmisdat)
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
      call writecdf(idout,varnam,dlon,jxm2,ixm2,1,iadim,xhro,lname,     &
                  & units,fact,offset,vvarmin,vvarmax,xlat1d,xlon1d,1.0,&
                  & 0,misdat,iotyp)
 
      end subroutine writehead
