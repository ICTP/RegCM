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

      subroutine xyobsll(iy,jx,iproj,clat,clon,plat,plon,truelath)
      use mod_block
      use mod_const
      use mod_constants , only : r2d => raddeg
      use mod_constants , only : d2r => degrad
      use mod_constants , only : pi => mathpi
      use mod_constants , only : a => erkm
      implicit none
!
! Dummy arguments
!
      real(4) :: clat , clon , plat , plon , truelath
      character(6) :: iproj
      integer :: iy , jx
      intent (in) clat , clon , iproj , iy , jx , truelath
!
! Local variables
!
      real(8) :: c1 , c2 , cell , cell2 , cntri , cntrj , flp , flpp ,  &
               & phi1 , phic , phir , phix , pole , psi1 , psx , r ,    &
               & xcntr , xlonx ,  xrot , ycntr
      real(8) :: xnr , ynr , xr , yr , pla , plo, cla , clo
      integer :: ie , ii , ilen , im , je
!
      pla = plat
      plo = plon
      cla = clat
      clo = clon
      ilen = iy*jx + 2
      ie = iy - 1
      je = jx - 1
      c1 = 1.454441E-4
      psi1 = 1.0E36
      pole = 90.
      if ( iproj=='LAMCON' ) psi1 = 90. - truelath
      if ( iproj=='POLSTR' ) psi1 = 30.0
      psi1 = psi1/r2d
!-----psi1 is colatitude of lat where cone or plane intersects earth

      if ( cla<0. ) then
        if ( truelath>0. ) then
          psi1 = -(90.-truelath)
        else
          psi1 = -(90.+truelath)
        end if
        pole = -90.
        psi1 = psi1/r2d
      end if
      if ( iproj=='LAMCON' .or. iproj=='POLSTR' ) then
        psx = (pole-cla)/r2d
        if ( iproj=='LAMCON' ) then
          cell = a*sin(psi1)/xn
          cell2 = (tan(psx/2.))/(tan(psi1/2.))
        else if ( iproj=='POLSTR' ) then
          cell = a*sin(psx)/xn
          cell2 = (1.+cos(psi1))/(1.+cos(psx))
        else
        end if
        r = cell*(cell2)**xn
        xcntr = 0.0
        ycntr = -r
      end if
      cntrj = float(je)/2.
      cntri = float(ie)/2.
!
!-----grid incoming data.  grdltmn=minimum latitude of incoming data.
!-----grdlnmn=minimum longitude of incoming data.
!
      do ii = 1 , nobs
        im = ii - 1
        if ( iproj=='LAMCON' .or. iproj=='POLSTR' ) then
          xrot = clo + 90./xn
          phix = yobs(ii)
          xlonx = xobs(ii)
          flpp = (xlonx-xrot)/r2d
          flp = xn*flpp
          psx = (pole-phix)/r2d
          if ( iproj=='LAMCON' ) then
            cell = a*sin(psi1)/xn
            cell2 = (tan(psx/2.))/(tan(psi1/2.))
          else if ( iproj=='POLSTR' ) then
            cell = a*sin(psx)/xn
            cell2 = (1.+cos(psi1))/(1.+cos(psx))
          else
          end if
          r = cell*(cell2)**xn
          xobs(ii) = (r*cos(flp)-xcntr)*1000.
          yobs(ii) = (r*sin(flp)-ycntr)*1000.
          if ( cla<0.0 ) xobs(ii) = -xobs(ii)
        end if
        if ( iproj=='NORMER' ) then
          phi1 = 0.0   ! plat/r2d
          phir = yobs(ii)/r2d
          phic = cla/r2d
          c2 = a*cos(phi1)
          cell = cos(phir)/(1.0+sin(phir))
          cell2 = cos(phic)/(1.0+sin(phic))
          ycntr = -c2*log(cell2)
          xobs(ii) = (c2*(xobs(ii)-clo)/r2d)*1000.
          yobs(ii) = (-c2*log(cell)-ycntr)*1000.
        end if
        if ( iproj=='ROTMER' ) then
          xcntr = plo - clo
          ycntr = pla - cla
          xnr = xobs(ii)
          ynr = yobs(ii)
          call nrot2rot(xnr,ynr,plo,pla,xr,yr)
          xobs(ii) = a*d2r*(xcntr+xr)*1000.
          yobs(ii) = a*d2r*(ycntr+yr)*1000.
        end if
        ht(ii) = ht(ii)/100.
        ht2(ii) = ht2(ii)/100000.
      end do

      end subroutine xyobsll
